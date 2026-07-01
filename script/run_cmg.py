"""
Run a CMG GEM ensemble simulation.

Usage
-----
  python script/run_cmg.py --template <cmg_data_folder> --ensemble <ensemble_folder> [options]

Examples
--------
  # Single deterministic run (no parameter variation)
  python script/run_cmg.py --template cmg_datafile --ensemble runs/base

  # 20-realization ensemble with xvar sampling
  python script/run_cmg.py --template cmg_datafile --ensemble runs/ens_01 \\
      --num-reals 20 --xvars xvars.json

  # Parallel ensemble (8 workers)
  python script/run_cmg.py --template cmg_datafile --ensemble runs/ens_01 \\
      --num-reals 20 --xvars xvars.json --parallel --workers 8

  # Overwrite an existing ensemble folder
  python script/run_cmg.py --template cmg_datafile --ensemble runs/ens_01 --overwrite
"""

import argparse
import glob
import json
import os
import shutil
import sys
import platform

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cmg_launcher import (
    generate_ensemble_folder,
    run_cmg_for_ensemble,
    run_cmg_for_ensemble_parallel,
)

# ── Default CMG executable ────────────────────────────────────────────────────

if platform.system() == "Windows":
    _DEFAULT_CMG_EXE = r'"C:\Program Files\CMG\GEM\2024.20\Win_x64\EXE\gm202420.exe"'
else:
    _DEFAULT_CMG_EXE = "/mnt/c/Program Files/CMG/GEM/2024.20/Win_x64/EXE/gm202420.exe"  # Adjust for other operating systems if needed


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Generate a CMG ensemble folder and run GEM simulations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--template", type=str, default="cmg_datafile_v2", metavar="PATH",
        help="CMG data template folder (contains one .dat file and *.inc files).",
    )
    p.add_argument(
        "--ensemble", type=str, default="runs/run_11", metavar="PATH",
        help="Output ensemble folder path (e.g. runs/ens_01).",
    )
    p.add_argument(
        "--num-reals", type=int, default=1, metavar="N",
        help="Number of realizations to generate (default: 1).",
    )
    p.add_argument(
        "--xvars", default=None, metavar="FILE",
        help="Path to xvars JSON file for parameter sampling. "
             "If omitted, the .dat file is copied verbatim.",
    )
    p.add_argument(
        "--parallel", type=bool, default=False,
        help="Run realizations in parallel using multiple processes.",
    )
    p.add_argument(
        "--workers", type=int, default=4, metavar="N",
        help="Number of parallel workers (default: 4; only used with --parallel).",
    )
    p.add_argument(
        "--cmg-exe", default=_DEFAULT_CMG_EXE, metavar="PATH",
        help=f"Path to the CMG GEM executable (default: {_DEFAULT_CMG_EXE}).",
    )
    p.add_argument(
        "--overwrite", type=bool, default=True,
        help="Delete and recreate the ensemble folder if it already exists.",
    )
    p.add_argument(
        "--no-run", action="store_true",
        help="Only generate the ensemble folder; do not launch simulations.",
    )
    p.add_argument(
        "--parasol", type=int, default=4, metavar="N",
        help="Number of parasol workers (default: 4; only used with --parallel).",
    )
    return p.parse_args()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    # ── Validate template folder ───────────────────────────────────────────────
    if not os.path.isdir(args.template):
        sys.exit(f"ERROR: template folder not found: {args.template}")

    dat_files = glob.glob(os.path.join(args.template, "*.dat"))
    if not dat_files:
        sys.exit(f"ERROR: no .dat file found in {args.template}")
    cmg_dat_name = os.path.basename(dat_files[0])

    # ── Load xvars ─────────────────────────────────────────────────────────────
    xvars = None
    if args.xvars:
        if not os.path.isfile(args.xvars):
            sys.exit(f"ERROR: xvars file not found: {args.xvars}")
        with open(args.xvars) as f:
            xvars = json.load(f)
        print(f"Loaded xvars from {args.xvars}  ({len(xvars)} variables)")

    # ── Handle existing ensemble folder ────────────────────────────────────────
    if os.path.exists(args.ensemble):
        if args.overwrite:
            print(f"Removing existing ensemble folder: {args.ensemble}")
            shutil.rmtree(args.ensemble)
        else:
            sys.exit(
                f"ERROR: ensemble folder already exists: {args.ensemble}\n"
                "Use --overwrite to replace it."
            )

    # ── Generate ensemble ──────────────────────────────────────────────────────
    print(f"Generating {args.num_reals} realization(s) in: {args.ensemble}")
    print(f"  template : {args.template}")
    print(f"  dat file : {cmg_dat_name}")
    print(f"  xvars    : {'none (deterministic)' if xvars is None else args.xvars}")

    generate_ensemble_folder(
        ensemble_dir=args.ensemble,
        cmg_template=args.template,
        num_of_reals=args.num_reals,
        xvars=xvars,
    )
    print(f"Ensemble folder ready: {args.ensemble}")

    # ── Run simulations ────────────────────────────────────────────────────────
    if args.no_run:
        print("--no-run flag set; skipping simulation launch.")
        return

    print(f"\nLaunching CMG GEM  ({args.cmg_exe})")
    if args.parallel and args.num_reals > 1:
        print(f"Mode: parallel  ({args.workers} workers)")
        run_cmg_for_ensemble_parallel(
            cmg_exe=args.cmg_exe,
            ensemble_dir=args.ensemble,
            cmg_data_file=cmg_dat_name,
            num_parasol_workers=args.parasol,
            max_workers=args.workers,
        )
    else:
        print("Mode: sequential")
        run_cmg_for_ensemble(
            cmg_exe=args.cmg_exe,
            ensemble_dir=args.ensemble,
            cmg_data_file=cmg_dat_name,
            num_parasol_workers=args.parasol,
        )

    print("\nAll simulations complete.")


if __name__ == "__main__":
    main()
