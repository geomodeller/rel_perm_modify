"""
run_cmg_with_rpt_update.py

Velocity-adaptive RTYPE simulation loop for CMG GEM.

Workflow
--------
  Outer loop  (time steps 0 … N_STEPS-1, one calendar month each):
    Inner loop  (RPT convergence):
      1. Write rocktype.inc with current RTYPE assignment.
      2. Restore RST file to start-of-step state (restart from end of prev step).
      3. Generate datafile_step_{t}.dat from template (one month only).
      4. Run CMG GEM simulation.
      5. Read Darcy velocity magnitude from the last time step of the SR3 file.
      6. Map velocity → RPT (bins of 3.048 m/day: RPT 1 = [0, 3.048), RPT 2 = [3.048, 6.096), …).
      7. If any block changed RPT: update and repeat inner loop.
         If all blocks consistent: advance to next time step.

Restart mechanism
-----------------
  Step 0  : runs from scratch.  CMG writes datafile_step_0.rst.
  Step t>0: datafile_step_{t-1}.rst is copied to datafile_step_{t}.rst before each
             inner iteration so CMG restarts from the converged end-state of step t-1.

Usage
-----
  python script/run_cmg_with_rpt_update.py [--run-dir runs/run_11] [options]
"""

import argparse
import os
import re
import shutil
import sys
import platform

import numpy as np

# Make script/ importable regardless of cwd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmg_launcher import run_cmg_simulator
from sr3_reader import read_SR3, get_grid_properties


# ── Simulation constants ───────────────────────────────────────────────────────
NX, NY, NZ = 38, 42, 43
N_BLOCKS   = NX * NY * NZ   # 68628

N_RPT      = 10
BIN_SIZE   = 3.048           # m/day (= 10 ft/day per RPT bin)

MAX_INNER  = 50              # max RPT convergence iterations per time step

if platform.system() == "Windows":
    _CMG_EXE = r'"C:\Program Files\CMG\GEM\2024.20\Win_x64\EXE\gm202420.exe"'
else:
    _CMG_EXE = "/mnt/c/Program Files/CMG/GEM/2024.20/Win_x64/EXE/gm202420.exe"


# ── RTYPE helpers ─────────────────────────────────────────────────────────────

def velocity_to_rtype(vel_flat):
    """
    Map Darcy velocity magnitudes (m/day, length N_BLOCKS) to RPT integers (1-10).

    Bin boundaries: 0, 3.048, 6.096, …, 30.48 m/day.
    Velocities >= 10 * BIN_SIZE are clamped to RPT 10.
    """
    rpt = (vel_flat / BIN_SIZE).astype(int) + 1
    return np.clip(rpt, 1, N_RPT).astype(int)


def write_rocktype_inc(rtype_flat, run_dir):
    """
    Write *RTYPE *ALL block to <run_dir>/rocktype.inc.
    Values are in CMG I-J-K order (I varies fastest), 20 per line.
    """
    path = os.path.join(run_dir, 'rocktype.inc')
    lines = ['*RTYPE *ALL\n']
    for i in range(0, N_BLOCKS, 20):
        chunk = rtype_flat[i:i + 20]
        lines.append(' '.join(str(v) for v in chunk) + '\n')
    with open(path, 'w') as f:
        f.writelines(lines)


# ── Template parsing ──────────────────────────────────────────────────────────

def parse_template(template_path):
    """
    Split the template dat file into:

    header_lines : list[str]
        Every line from the start through the end of the GCONI block (inclusive).
        Contains $TIME_STEP$ placeholder in the FILENAME SR3-IN line.

    date_list : list[str]
        Active (uncommented) DATE lines found after the header, in order.
        date_list[0]  = 'DATE 2026 2 1.0'  (end of time step 0)
        date_list[23] = 'DATE 2028 1 1.0'  (end of time step 23)

    results_spec : list[str]
        Lines after the STOP keyword (RESULTS SPEC post-processing).
    """
    with open(template_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Locate end of GCONI block (GCONI line + its one data line)
    gconi_end = None
    date_re   = re.compile(r'^\s*DATE\s+', re.IGNORECASE)
    for i, line in enumerate(lines):
        if line.strip().startswith('GCONI'):
            j = i + 1
            while j < len(lines):
                s = lines[j].strip()
                if s and not s.startswith('**'):
                    gconi_end = j + 1   # one past the GCONI data line
                    break
                j += 1
            break

    if gconi_end is None:
        raise ValueError("GCONI block not found in template")

    # Collect active DATE lines from the recurrent section
    date_list = []
    stop_idx  = None
    for i, line in enumerate(lines[gconi_end:], start=gconi_end):
        s = line.strip()
        if s == 'STOP':
            stop_idx = i
            break
        if date_re.match(s) and not s.startswith('**'):
            date_list.append(s.rstrip())   # e.g. 'DATE 2026 2 1.0'

    results_spec = lines[stop_idx + 1:] if stop_idx is not None else []

    return lines[:gconi_end], date_list, results_spec


# ── Dat file generation ───────────────────────────────────────────────────────

def generate_dat_content(header_lines, target_date_str, results_spec, step):
    """
    Build a dat file that runs ONE calendar month.

    For step t: CMG restarts from the PREVIOUS step's SR3 file via REWIND 1.
    FILENAME SR3-IN is set to 'datafile_step_{step-1}.sr3' so CMG reads the
    end state of the previous step for restart and writes output there.

    Parameters
    ----------
    header_lines    : list[str]  Template header (up through GCONI).
    target_date_str : str        'DATE YYYY M D.0' for the end of this step.
    results_spec    : list[str]  Lines after STOP (post-processing metadata).
    step            : int        Outer time-step index (1-based; SR3-IN uses step-1).
    """
    sr3_in_step = step - 1  # FILENAME SR3-IN references the previous step's SR3
    out = [line.replace('$TIME_STEP$', str(sr3_in_step)) for line in header_lines]
    out.append(f'{target_date_str}\n')
    out.append('STOP\n')
    out.append('\n')
    out.extend(results_spec)
    return ''.join(out)


# ── SR3 velocity reader ───────────────────────────────────────────────────────

def read_last_velocity(sr3_path):
    """
    Read the Darcy velocity magnitude at the final time step of an SR3 file.

    Returns
    -------
    vel_flat : ndarray, shape (N_BLOCKS,)
        Velocity magnitude in m/day, flattened in CMG I-J-K order (I fastest).
        Inactive blocks receive zero.
    """
    sr3 = read_SR3(sr3_path)
    nt  = len(sr3.times['Days'])
    if nt == 0:
        raise RuntimeError(f"No time steps found in SR3: {sr3_path}")

    gp  = get_grid_properties(sr3, NX, NY, NZ, nt)

    zeros = np.zeros((nt, NZ, NY, NX))
    vx = gp.get('VELGXRC', zeros)
    vy = gp.get('VELGYRC', zeros)
    vz = gp.get('VELGZRC', zeros)

    # Magnitude at last time step, shape (NZ, NY, NX)
    vel_last = np.sqrt(vx[-1]**2 + vy[-1]**2 + vz[-1]**2)

    # Flatten in CMG order: i varies fastest → shape is already (NZ, NY, NX)
    # so C-order flatten gives last axis (NX=i) varying fastest.
    return vel_last.flatten()


# ── RST helpers ───────────────────────────────────────────────────────────────

def rst_path(run_dir, step):
    return os.path.join(run_dir, f'datafile_step_{step}.rst')


def start_rst_path(run_dir, step):
    return os.path.join(run_dir, f'_step_{step}_start.rst')


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--run-dir',  default=os.path.join('runs', 'run_11'),
                   help='Directory containing template.dat and *.inc files (default: runs/run_11)')
    p.add_argument('--template', default='template.dat',
                   help='Template dat filename inside --run-dir (default: template.dat)')
    p.add_argument('--n-steps',  type=int, default=24,
                   help='Number of time steps to run (default: 24)')
    p.add_argument('--start-step', type=int, default=1,
                   help='First time step index to run (default: 1, since step 0 already completed)')
    p.add_argument('--cmg-exe',  default=_CMG_EXE,
                   help='Path to CMG GEM executable')
    p.add_argument('--parasol',  type=int, default=4,
                   help='Number of parasol workers (default: 4)')
    p.add_argument('--max-inner', type=int, default=MAX_INNER,
                   help=f'Max RPT convergence iterations per step (default: {MAX_INNER})')
    args = p.parse_args()

    run_dir       = args.run_dir
    template_path = os.path.join(run_dir, args.template)

    if not os.path.isfile(template_path):
        sys.exit(f'ERROR: template not found: {template_path}')

    # ── Parse template ─────────────────────────────────────────────────────────
    header_lines, date_list, results_spec = parse_template(template_path)

    n_available = len(date_list)
    if n_available < args.n_steps:
        sys.exit(f'ERROR: template has {n_available} end-dates but --n-steps={args.n_steps}')

    start_step = args.start_step
    if start_step < 1:
        sys.exit('ERROR: --start-step must be >= 1 (step 0 must already be completed)')

    print(f'Template   : {template_path}')
    print(f'Run dir    : {run_dir}')
    print(f'Time steps : {start_step}–{args.n_steps-1}  ({date_list[start_step]} → {date_list[args.n_steps-1]})')

    # ── Initialise RTYPE from previous step SR3 (already completed) ───────────
    sr3_prev = os.path.join(run_dir, f'datafile_step_{start_step - 1}.sr3')
    if not os.path.isfile(sr3_prev):
        sys.exit(f'ERROR: SR3 for step {start_step - 1} not found: {sr3_prev}\n'
                 f'       Complete step {start_step - 1} before starting from step {start_step}.')

    print(f'Reading initial RTYPE from {sr3_prev} ...')
    vel0  = read_last_velocity(sr3_prev)
    rtype = velocity_to_rtype(vel0)
    write_rocktype_inc(rtype, run_dir)
    counts0 = [int(np.sum(rtype == r)) for r in range(1, N_RPT + 1)]
    active0 = [f'RPT{r}={c}' for r, c in enumerate(counts0, 1) if c > 0]
    print(f'Initial RPT distribution: {", ".join(active0)}\n')

    # ── Outer time-step loop ───────────────────────────────────────────────────
    for step in range(start_step, args.n_steps):
        target_date = date_list[step]          # 'DATE YYYY M D.0'
        dat_name    = f'datafile_step_{step}.dat'
        # FILENAME SR3-IN uses step-1, so CMG writes output into the previous
        # step's SR3 file (via REWIND 1 restart + appending new time steps).
        sr3_name    = f'datafile_step_{step - 1}.sr3'
        dat_path    = os.path.join(run_dir, dat_name)
        sr3_out     = os.path.join(run_dir, sr3_name)

        print(f'┌─ Step {step:2d}/{args.n_steps-1}  →  {target_date}')

        # Save the start-of-step RST so inner loop can re-run from the same state
        if step > 0:
            prev_rst  = rst_path(run_dir, step - 1)
            saved_rst = start_rst_path(run_dir, step)
            if not os.path.isfile(prev_rst):
                sys.exit(f'ERROR: RST from step {step-1} not found: {prev_rst}\n'
                         f'       Cannot restart step {step}.')
            shutil.copy(prev_rst, saved_rst)

        # ── Inner RPT convergence loop ─────────────────────────────────────────
        for inner in range(args.max_inner):

            # Restore start-of-step RST before each CMG run
            if step > 0:
                shutil.copy(start_rst_path(run_dir, step), rst_path(run_dir, step))

            # Write current RTYPE
            rtype[rtype>N_RPT] = N_RPT   # clamp to max RPT
            rtype[rtype<1]   = 1         # clamp to min RPT
            write_rocktype_inc(rtype, run_dir)

            # Generate dat file (one month only)
            dat_content = generate_dat_content(header_lines, target_date, results_spec, step)
            with open(dat_path, 'w', encoding='utf-8') as f:
                f.write(dat_content)

            # Run CMG GEM
            print(f'│  [inner {inner:2d}] Running CMG...', end=' ', flush=True)
            run_cmg_simulator(args.cmg_exe, run_dir, dat_name, args.parasol)

            if not os.path.isfile(sr3_out):
                sys.exit(f'\nERROR: SR3 not produced: {sr3_out}\n'
                         f'       Check {os.path.join(run_dir, dat_name.replace(".dat",".out"))}')

            # Read Darcy velocity and map to new RTYPE
            vel_flat  = read_last_velocity(sr3_out)
            new_rtype = velocity_to_rtype(vel_flat)

            n_changed = int(np.sum(new_rtype != rtype))
            print(f'changed={n_changed} blocks')

            if n_changed == 0:
                print(f'│  Converged in {inner+1} iteration(s)')
                rtype = new_rtype
                break

            rtype = new_rtype

        else:
            print(f'│  WARNING: did not converge within {args.max_inner} iterations')

        # Clean up temporary start RST
        sav = start_rst_path(run_dir, step)
        if os.path.isfile(sav):
            os.remove(sav)

        # Print current RPT distribution
        counts = [int(np.sum(rtype == r)) for r in range(1, N_RPT + 1)]
        active = [f'RPT{r}={c}' for r, c in enumerate(counts, 1) if c > 0]
        print(f'└─ RPT distribution: {", ".join(active)}\n')

    # ── Summary ────────────────────────────────────────────────────────────────
    print('=' * 60)
    print(f'Simulation complete.  Final RPT distribution:')
    for r in range(1, N_RPT + 1):
        n = int(np.sum(rtype == r))
        lo = (r - 1) * BIN_SIZE
        hi = r * BIN_SIZE
        print(f'  RPT {r:2d}  [{lo:6.3f}, {hi:6.3f}) m/day : {n:6d} blocks  ({n/N_BLOCKS*100:.1f}%)')


if __name__ == '__main__':
    main()
