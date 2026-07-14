"""
Velocity-adaptive RTYPE loop for the Porthos field model WITH local grid
refinement (runs/field_porthos_with-update_lgr).

Same workflow as run_cmg_with_rpt_update_4_field.py, but LGR-aware:

  * The template contains REFINE 19:20 11:12 1:43 INTO 20 20 1 (5 m children).
  * Velocities are read for the fundamental grid AND for every refined child
    cell via sr3_reader.get_grid_properties_lgr().
  * rocktype.inc is written as
        *RTYPE *ALL            (68628 fundamental values)
        *RTYPE *RG i j k *ALL  (400 child values, one block per refined parent,
                                local-I-fastest order = SR3 child order)
  * RPT convergence is checked over fundamental + child cells together.
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
from sr3_reader import read_SR3, get_grid_properties_lgr


# ── Simulation constants ───────────────────────────────────────────────────────
# GRID CORNER 38 42 43 + REFINE 19:20 11:12 1:43 INTO 20 20 1
NX, NY, NZ = 38, 42, 43
N_BLOCKS   = NX * NY * NZ   # 68628 fundamental blocks

N_RPT      = 30
BIN_SIZE   = 1.016           # m/day (30 bins over the same 0-30.48 m/day range)

INIT_RPT   = 16              # initial rock type (all blocks, children inherit)

MAX_INNER  = 10              # max RPT convergence iterations per time step

if platform.system() == "Windows":
    _CMG_EXE = r'C:\Program Files\CMG\GEM\2024.20\Win_x64\EXE\gm202420.exe'
else:
    _CMG_EXE = "/mnt/c/Program Files/CMG/GEM/2024.20/Win_x64/EXE/gm202420.exe"


# ── RTYPE helpers ─────────────────────────────────────────────────────────────

def velocity_to_rtype(vel_flat):
    """
    Map Darcy velocity magnitudes (m/day) to RPT integers (1-30).

    Bin boundaries: 0, 1.016, 2.032, …, 30.48 m/day.
    Velocities >= 30 * BIN_SIZE are clamped to RPT 30.
    """
    rpt = (vel_flat / BIN_SIZE).astype(int) + 1
    return np.clip(rpt, 1, N_RPT).astype(int)


def _chunk_lines(values, per_line=20):
    return [' '.join(str(v) for v in values[i:i + per_line]) + '\n'
            for i in range(0, len(values), per_line)]


def write_rocktype_inc(rtype_fund, rtype_child, lgr_table, run_dir):
    """
    Write rocktype.inc with fundamental + refined-grid rock types.

    rtype_fund : (N_BLOCKS,) ints — fundamental grid, CMG natural order
                 (values of refined parents are ignored by CMG).
    rtype_child: (n_child_cells,) ints — refined cells in SR3 child order.
    lgr_table  : from sr3_reader.get_lgr_table(); gives each refined parent's
                 (i,j,k) and its slice in rtype_child.  Child values are
                 local-I-fastest, matching what *RG ... *ALL expects.
    """
    path = os.path.join(run_dir, 'rocktype.inc')
    lines = ['*RTYPE *ALL\n']
    lines += _chunk_lines(rtype_fund)
    for rec in lgr_table:
        i, j, k = rec['parent_ijk']
        lines.append(f'*RTYPE *RG {i} {j} {k} *ALL\n')
        lines += _chunk_lines(rtype_child[rec['start']:rec['stop']])
    with open(path, 'w') as f:
        f.writelines(lines)


def write_rocktype_inc_initial(run_dir):
    """Uniform initial rock type; refined children inherit from their parent."""
    path = os.path.join(run_dir, 'rocktype.inc')
    with open(path, 'w') as f:
        f.write(f'*RTYPE *CON {INIT_RPT}\n')


# ── Dat file generation ───────────────────────────────────────────────────────

def generate_dat_content(template_path, step):
    """
    Build a dat file for one time step (same convention as the field script).

    The template comments its restart directives with a '**$' prefix:
        **$FILENAME SR3-IN 'datafile_step_{$step$}.sr3'
        **$RESTART
        **$REWIND 1
        **$WRST TNEXT

    step 0  : keep them commented (fresh run from scratch).
    step t>0: strip the '**$' prefix and substitute {$step$} with step-1.
    """
    with open(template_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    if step == 0:
        return ''.join(lines)
    out = []
    for line in lines:
        if line.startswith('**$'):
            line = line[3:].replace('{$step$}', str(step - 1))
        out.append(line)
    return ''.join(out)


# ── SR3 velocity reader (LGR-aware) ──────────────────────────────────────────

def read_last_velocity(sr3_path):
    """
    Read Darcy velocity magnitudes at the final time step of an SR3 file.

    Returns
    -------
    vel_fund  : (N_BLOCKS,) magnitudes on the fundamental grid, CMG natural
                order (refined parents and inactive blocks are zero).
    vel_child : (n_child_cells,) magnitudes on the refined cells, SR3 child
                order (see lgr_table).
    lgr_table : refined-grid table (parent addresses + slices), needed to
                write the *RG rocktype blocks.
    """
    sr3 = read_SR3(sr3_path)
    nt  = len(sr3.times['Days'])
    if nt == 0:
        raise RuntimeError(f"No time steps found in SR3: {sr3_path}")

    fund, children, lgr_table = get_grid_properties_lgr(sr3, NX, NY, NZ, nt)

    zeros_f = np.zeros((nt, NZ, NY, NX))
    vxf = fund.get('VELGXRC', zeros_f)
    vyf = fund.get('VELGYRC', zeros_f)
    vzf = fund.get('VELGZRC', zeros_f)
    vel_fund = np.sqrt(vxf[-1]**2 + vyf[-1]**2 + vzf[-1]**2).flatten()

    n_child = lgr_table[-1]['stop'] if lgr_table else 0
    zeros_c = np.zeros((nt, n_child))
    vxc = children.get('VELGXRC', zeros_c)
    vyc = children.get('VELGYRC', zeros_c)
    vzc = children.get('VELGZRC', zeros_c)
    vel_child = np.sqrt(vxc[-1]**2 + vyc[-1]**2 + vzc[-1]**2)

    return vel_fund, vel_child, lgr_table


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--run-dir',  default=os.path.join('runs', 'field_porthos_with-update_lgr'),
                   help='Directory containing template.dat and *.inc files '
                        '(default: runs/field_porthos_with-update_lgr)')
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

    # ── Read schedule template (DATE / WRST TNEXT pairs, one per step) ────────
    schedule_path = os.path.join(run_dir, 'schedule.inc')
    schedule_template_path = os.path.join(run_dir, 'template-schedule.inc')
    with open(schedule_template_path, 'r') as f:
        schedule_lines = [l for l in f.readlines() if l.strip()]
    n_steps = len(schedule_lines) // 2

    print(f'Time steps : {start_step}–{n_steps-1} ')

    # ── State: current RTYPE (fundamental + children) ─────────────────────────
    rtype_fund  = np.full(N_BLOCKS, INIT_RPT, dtype=int)
    rtype_child = None            # unknown until the first SR3 is read
    lgr_table   = None

    if start_step > 0:
        # Resume: initialise RTYPE from the previous step's SR3
        sr3_prev = os.path.join(run_dir, f'datafile_step_{start_step - 1}.sr3')
        if not os.path.isfile(sr3_prev):
            sys.exit(f'ERROR: cannot resume from step {start_step}: {sr3_prev} not found')
        print(f'Initialising RTYPE from {sr3_prev} ...')
        vel_fund, vel_child, lgr_table = read_last_velocity(sr3_prev)
        rtype_fund  = velocity_to_rtype(vel_fund)
        rtype_child = velocity_to_rtype(vel_child)

    # ── Outer time-step loop ───────────────────────────────────────────────────
    for step in range(start_step, n_steps):
        dat_name = f'datafile_step_{step}.dat'
        sr3_name = f'datafile_step_{step}.sr3'
        dat_path = os.path.join(run_dir, dat_name)
        sr3_out  = os.path.join(run_dir, sr3_name)

        with open(schedule_path, 'w') as f:
            f.writelines(schedule_lines[:2*step+2])

        print(f'┌─ Step {step:2d} / {n_steps-1} ────────────────────────────────────────────────')

        # ── Inner RPT convergence loop ─────────────────────────────────────────
        for inner in range(args.max_inner):

            # Write current RTYPE
            if rtype_child is None:
                write_rocktype_inc_initial(run_dir)     # uniform, children inherit
            else:
                write_rocktype_inc(rtype_fund, rtype_child, lgr_table, run_dir)

            # Generate dat file
            dat_content = generate_dat_content(template_path, step)
            with open(dat_path, 'w', encoding='utf-8') as f:
                f.write(dat_content)

            # Run CMG GEM
            print(f'│  [inner {inner:2d}] Running CMG...', end=' ', flush=True)
            run_cmg_simulator(args.cmg_exe, run_dir, dat_name, args.parasol)

            if not os.path.isfile(sr3_out):
                sys.exit(f'\nERROR: SR3 not produced: {sr3_out}\n'
                         f'       Check {os.path.join(run_dir, dat_name.replace(".dat",".out"))}')

            # Read velocities (fundamental + children) and map to new RTYPE
            vel_fund, vel_child, lgr_table = read_last_velocity(sr3_out)
            new_fund  = velocity_to_rtype(vel_fund)
            new_child = velocity_to_rtype(vel_child)

            if rtype_child is None:
                rtype_child = np.full(len(new_child), INIT_RPT, dtype=int)

            n_changed = int(np.sum(new_fund != rtype_fund)) \
                      + int(np.sum(new_child != rtype_child))
            print(f'changed={n_changed} blocks '
                  f'(fund={int(np.sum(new_fund != rtype_fund))}, '
                  f'lgr={int(np.sum(new_child != rtype_child))})')

            rtype_fund, rtype_child = new_fund, new_child

            if n_changed == 0:
                print(f'│  Converged in {inner+1} iteration(s)')
                break

        else:
            print(f'│  WARNING: did not converge within {args.max_inner} iterations')

        # Print current RPT distribution (fundamental + children together)
        all_rpt = np.concatenate([rtype_fund, rtype_child])
        counts = [int(np.sum(all_rpt == r)) for r in range(1, N_RPT + 1)]
        active = [f'RPT{r}={c}' for r, c in enumerate(counts, 1) if c > 0]
        print(f'└─ RPT distribution: {", ".join(active)}\n')

    # ── Summary ────────────────────────────────────────────────────────────────
    all_rpt = np.concatenate([rtype_fund, rtype_child])
    n_all   = len(all_rpt)
    print('=' * 60)
    print(f'Simulation complete.  Final RPT distribution (fund + LGR children):')
    for r in range(1, N_RPT + 1):
        n = int(np.sum(all_rpt == r))
        lo = (r - 1) * BIN_SIZE
        hi = r * BIN_SIZE
        print(f'  RPT {r:2d}  [{lo:6.3f}, {hi:6.3f}) m/day : {n:6d} blocks  ({n/n_all*100:.1f}%)')


if __name__ == '__main__':
    print(f"CMG executable: {_CMG_EXE}")
    main()
