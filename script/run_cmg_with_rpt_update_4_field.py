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
# GRID CORNER 38 42 43 → cartesian corner-point grid
NX, NY, NZ = 38, 42, 43
N_BLOCKS   = NX * NY * NZ   # 68628

N_RPT      = 30
BIN_SIZE   = 1.016           # m/day (30 bins over the same 0-30.48 m/day range)

INIT_RPT   = 16              # initial rock type (matches rocktype.inc *CON 16)

MAX_INNER  = 10              # max RPT convergence iterations per time step

if platform.system() == "Windows":
    _CMG_EXE = r'C:\Program Files\CMG\GEM\2024.20\Win_x64\EXE\gm202420.exe'
else:
    _CMG_EXE = "/mnt/c/Program Files/CMG/GEM/2024.20/Win_x64/EXE/gm202420.exe"


# ── RTYPE helpers ─────────────────────────────────────────────────────────────

def velocity_to_rtype(vel_flat):
    """
    Map Darcy velocity magnitudes (m/day, length N_BLOCKS) to RPT integers (1-30).

    Bin boundaries: 0, 1.016, 2.032, …, 30.48 m/day.
    Velocities >= 30 * BIN_SIZE are clamped to RPT 30.
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
def write_rocktype_inc_initial(run_dir):
    """
    Write *RTYPE *ALL block to <run_dir>/rocktype.inc.
    Values are in CMG I-J-K order (I varies fastest), 20 per line.
    """
    path = os.path.join(run_dir, 'rocktype.inc')
    lines = ['*RTYPE *ALL\n']
    rtype_flat = [INIT_RPT]*N_BLOCKS  # Initialize all blocks to the middle RPT
    for i in range(0, N_BLOCKS, 20):
        chunk = rtype_flat[i:i + 20]
        lines.append(' '.join(str(v) for v in chunk) + '\n')
    with open(path, 'w') as f:
        f.writelines(lines)


# ── Template parsing ──────────────────────────────────────────────────────────


# ── Dat file generation ───────────────────────────────────────────────────────

def generate_dat_content(template_path, step):
    """
    Build a dat file that runs ONE calendar month.

    The field template comments its restart directives with a '**$' prefix:
        **$FILENAME SR3-IN 'datafile_step_{$step$}.sr3'
        **$RESTART
        **$REWIND 1
        **$WRST TNEXT

    step 0  : keep them commented (fresh run from scratch).
    step t>0: strip the '**$' prefix to activate them and substitute
              {$step$} with step-1 so CMG restarts from the previous
              step's SR3 file.

    Parameters
    ----------
    template_path : str   Path to the template dat file.
    step          : int   Outer time-step index.
    """

    with open(template_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    if step == 0:
        return ''.join(lines)  # fresh run: restart directives stay commented
    else:
        out = []
        for line in lines:
            if line.startswith('**$'):
                line = line[3:].replace('{$step$}', str(step - 1))
            out.append(line)
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
    p.add_argument('--run-dir',  default=os.path.join('runs', 'field_porthos_with-update'),
                   help='Directory containing template.dat and *.inc files (default: runs/field_porthos_with-update)')
    p.add_argument('--template', default='template.dat',
                   help='Template dat filename inside --run-dir (default: template.dat)')
    p.add_argument('--start-step', type=int, default=0,
                   help='First time step index to run (default: 0)')
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

    start_step = args.start_step
    if start_step < 0:
        sys.exit('ERROR: --start-step must be >= 0')

    print(f'Template   : {template_path}')
    print(f'Run dir    : {run_dir}')

    # ── Read schedule template ─────────────────────────────────────────────────
    # template-schedule.inc holds DATE / WRST TNEXT pairs, one pair per step
    schedule_path = os.path.join(run_dir, 'schedule.inc')
    schedule_template_path = os.path.join(run_dir, 'template-schedule.inc')
    with open(schedule_template_path, 'r') as f:
        schedule_lines = [l for l in f.readlines() if l.strip()]
    n_steps = len(schedule_lines) // 2

    print(f'Time steps : {start_step}–{n_steps-1} ')
    # ── Outer time-step loop ───────────────────────────────────────────────────
    for step in range(start_step, n_steps):
        dat_name    = f'datafile_step_{step}.dat'
        sr3_name    = f'datafile_step_{step}.sr3'
        dat_path    = os.path.join(run_dir, dat_name)
        sr3_out     = os.path.join(run_dir, sr3_name)

        with open(schedule_path, 'w') as f:
            f.writelines(schedule_lines[:2*step+2])

        print(f'┌─ Step {step:2d} / {n_steps-1} ────────────────────────────────────────────────')

        # ── Inner RPT convergence loop ─────────────────────────────────────────
        for inner in range(args.max_inner):


            # Write current RTYPE
            if step == 0 and inner == 0:
                write_rocktype_inc_initial(run_dir)
                rtype = np.full(N_BLOCKS, INIT_RPT, dtype=int)
            else:
                write_rocktype_inc(rtype, run_dir)

            # Generate dat file (one month only)
            dat_content = generate_dat_content(template_path, step)
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
    print(f"CMG executable: {_CMG_EXE}")
    main()
