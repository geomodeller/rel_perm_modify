"""
CMG GEM simulation launcher utilities.

Functions
---------
generate_ensemble_folder   : copy template into real_0 … real_N folders, optionally sampling xvars
run_cmg_simulator          : run a single CMG simulation in its folder
run_cmg_for_ensemble       : run all realizations sequentially
run_cmg_for_ensemble_parallel : run all realizations in parallel
generate_realization       : sample one set of xvar values
"""

import os, subprocess
import glob
import random
import shutil
import concurrent.futures
import numpy as np


# ── CMG execution ─────────────────────────────────────────────────────────────

def run_cmg_simulator(cmg_exe, real_dir, cmg_data_file, num_parasol_workers):
    """
    Run a single CMG simulation.

    Parameters
    ----------
    cmg_exe       : str  Full path to the CMG GEM executable (may include quotes).
    real_dir      : str  Realization directory to cd into before running.
    cmg_data_file : str  Data file name (basename only, relative to real_dir).
    num_parasol_workers : int  Number of parasol workers to use.
    """
    current_dir = os.getcwd()
    os.chdir(real_dir)
    subprocess.run([cmg_exe, "-f", cmg_data_file, "-parasol", str(num_parasol_workers)],
                    stdout=subprocess.DEVNULL,   # Throw away normal output
                    stderr=subprocess.DEVNULL)   # Throw away errors
                    
    os.chdir(current_dir)


def run_cmg_for_ensemble(cmg_exe, ensemble_dir, cmg_data_file, num_parasol_workers):
    """Run all real_* folders in ensemble_dir sequentially."""
    real_dirs = sorted(glob.glob(os.path.join(ensemble_dir, 'real*')))
    for real_dir in real_dirs:
        run_cmg_simulator(cmg_exe, real_dir, cmg_data_file, num_parasol_workers)


def run_cmg_for_ensemble_parallel(cmg_exe, ensemble_dir, cmg_data_file, num_parasol_workers, max_workers=4):
    """Run all real_* folders in ensemble_dir in parallel."""
    real_dirs = sorted(glob.glob(os.path.join(ensemble_dir, 'real*')))
    n = len(real_dirs)
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(run_cmg_simulator,
                     [cmg_exe] * n,
                     real_dirs,
                     [cmg_data_file] * n,
                     [num_parasol_workers] * n)


# ── Ensemble folder generation ────────────────────────────────────────────────

def generate_ensemble_folder(ensemble_dir, cmg_template, num_of_reals, xvars=None):
    """
    Build ensemble_dir/real_0 … real_{N-1} from the template folder.

    Each realization folder receives:
      - all .inc files (copied verbatim)
      - the .dat file (copied verbatim when xvars is None;
                       written with substitutions when xvars is provided)

    Parameters
    ----------
    ensemble_dir  : str   Output top-level directory (must not already exist).
    cmg_template  : str   Template folder containing one .dat and *.inc files.
    num_of_reals  : int   Number of realizations to create.
    xvars         : dict  Sampling configuration (see generate_realization).
                          None → copy .dat without substitution.
    """
    dat_files = glob.glob(os.path.join(cmg_template, '*.dat'))
    if not dat_files:
        raise FileNotFoundError(f"No .dat file found in {cmg_template}")
    cmg_dat_src = dat_files[0]
    cmg_dat_name = os.path.basename(cmg_dat_src)

    inc_files = glob.glob(os.path.join(cmg_template, '*.inc'))

    if os.path.exists(ensemble_dir):
        raise FileExistsError(f"Ensemble folder already exists: {ensemble_dir}")
    os.makedirs(ensemble_dir)

    for i in range(num_of_reals):
        real_dir = os.path.join(ensemble_dir, f'real_{i}')
        os.makedirs(real_dir)

        # Always copy .inc files verbatim
        copy_selected_files(inc_files, real_dir)

        # Copy or substitute .dat file
        dat_dst = os.path.join(real_dir, cmg_dat_name)
        if xvars is not None:
            xvar_real = generate_realization(xvars)
            replace_words_in_file(cmg_dat_src, dat_dst, xvar_real)
        else:
            shutil.copy2(cmg_dat_src, dat_dst)


# ── xvar sampling ─────────────────────────────────────────────────────────────

def generate_realization(xvars):
    """
    Sample one realization from the xvars specification.

    Accepted formats per variable
    ------------------------------
    [min, max]                           → uniform in [min, max]
    {"type"|"distribution": "uniform",  "range": [min, max]}
    {"type"|"distribution": "gaussian", "mean": m, "std": s, "range": [lo, hi]}
    {"type"|"distribution": "categorical", "choices"|"class": [...]}

    Returns
    -------
    dict mapping each key to a string value ready for text substitution.
    """
    xvars_real = {}

    for key, config in xvars.items():
        # Shorthand: bare [min, max] list
        if isinstance(config, list):
            xvars_real[key] = str(round(random.uniform(config[0], config[1]), 6))
            continue

        # Accept "distribution" as alias for "type"
        var_type = config.get("type") or config.get("distribution", "uniform")

        if var_type == "uniform":
            lo, hi = config.get("range", [0, 1])
            xvars_real[key] = str(round(random.uniform(lo, hi), 6))

        elif var_type == "gaussian":
            mean = config.get("mean", 0)
            std  = config.get("std", 1)
            lo, hi = config.get("range", [-np.inf, np.inf])
            value = np.clip(np.random.normal(mean, std), lo, hi)
            xvars_real[key] = str(round(float(value), 6))

        elif var_type == "categorical":
            # Accept "class" as alias for "choices"
            choices = config.get("choices") or config.get("class", ["default"])
            xvars_real[key] = str(random.choice(choices))

        else:
            raise ValueError(f"Unknown variable type '{var_type}' for key '{key}'")

    return xvars_real


# ── File utilities ────────────────────────────────────────────────────────────

def copy_selected_files(file_list, destination_dir, verbose=False):
    """Copy each file in file_list into destination_dir."""
    os.makedirs(destination_dir, exist_ok=True)
    for src in file_list:
        if os.path.isfile(src):
            dst = os.path.join(destination_dir, os.path.basename(src))
            shutil.copy2(src, dst)
            if verbose:
                print(f"Copied: {os.path.basename(src)}")
        else:
            print(f"File not found: {src}")
    if verbose:
        print("Done.")


def replace_words_in_file(input_file, output_file, xvar_real, verbose=False):
    """
    Write output_file as a copy of input_file with all xvar_real keys
    replaced by their corresponding values.
    """
    with open(input_file, 'r', encoding='utf-8') as f:
        content = f.read()
    for key, value in xvar_real.items():
        content = content.replace(key, value)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(content)
    if verbose:
        print(f"Written: {output_file}")
