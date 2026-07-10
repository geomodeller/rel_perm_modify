"""
Generate cmg_datafile_v2/rock_fluid_properties_30bins.inc
with 30 rock types (RPT 1 - RPT 30) using the same Corey power-law model as
cmg_datafile_v2/rock_fluid_properties.inc (10 rock types).

Same TOTAL parameter change as the 10-bin file, split into 30 bins
(29 steps) instead of 10 (9 steps):

  Varying across RPT 1 -> RPT 30 (linear interpolation)
    S_wc    (residual/connate water saturation) : 0.55 -> 0.40
    krg_max (endpoint CO2 rel perm)             : 0.25 -> 0.45

  Fixed (all RPTs)
    krw_max = 1.0,  S_gr (residual CO2) = 0.0,  Pc = 0

  Model (Corey power-law, matching the original file)
    krw  = Se_w^4,          Se_w = (Sw - S_wc) / (1 - S_wc)
    krow = (1 - Se_w)^2
    krg  = krg_max * Se_g^2, Se_g = Sg / (1 - S_wc)
    krog = 0

The script also verifies its formulas by regenerating RPT 1 and RPT 10 of the
original 10-bin file and comparing every table value.

Usage
-----
  python3 script/generate_30_binned_relperm_curve.py
"""

import os
import re
import sys

import numpy as np

# ── Constants ─────────────────────────────────────────────────────────────────

N_RPT = 30          # number of rock types (bins) in the new file
N_SAT = 46          # points in the mobile saturation range (same as 10-bin file)

# Total parameter change identical to the 10-bin file:
S_WC_START,  S_WC_END  = 0.55, 0.40   # residual water saturation
KRG_MAX_START, KRG_MAX_END = 0.25, 0.45  # endpoint CO2 (gas) rel perm

KRW_MAX = 1.0       # fixed endpoint water rel perm
N_W, N_OW, N_G = 4, 2, 2   # Corey exponents: krw, krow, krg

ORIG_INC = os.path.join('cmg_datafile_v2', 'rock_fluid_properties.inc')
OUT_INC  = os.path.join('cmg_datafile_v2', 'rock_fluid_properties_30bins.inc')


# ── Corey table builders ──────────────────────────────────────────────────────

def build_swt(S_wc):
    """Return (Sw, krw, krow, Pc) arrays for the water-gas SWT table."""
    Sw   = np.concatenate([[0.0], np.linspace(S_wc, 1.0, N_SAT)])
    Se_w = np.clip((Sw - S_wc) / (1.0 - S_wc), 0.0, 1.0)
    krw  = KRW_MAX * Se_w ** N_W
    krow = (1.0 - Se_w) ** N_OW
    Pc   = np.zeros_like(Sw)
    return Sw, krw, krow, Pc


def build_sgt(S_wc, krg_max):
    """Return (Sg, krg, krog, Pc) arrays for the gas SGT table."""
    Sg   = np.linspace(0.0, 1.0 - S_wc, N_SAT)
    Se_g = Sg / (1.0 - S_wc)
    krg  = krg_max * Se_g ** N_G
    krog = np.zeros_like(Sg)
    Pc   = np.zeros_like(Sg)
    return Sg, krg, krog, Pc


# ── Parameter variation across RPTs ──────────────────────────────────────────

def make_rpt_params(n_rpt):
    """Linear interpolation of (S_wc, krg_max) from RPT 1 to RPT n_rpt."""
    params = []
    for n in range(n_rpt):
        frac    = n / (n_rpt - 1)   # 0.0 at RPT 1,  1.0 at RPT n_rpt
        S_wc    = S_WC_START  + frac * (S_WC_END  - S_WC_START)
        krg_max = KRG_MAX_START + frac * (KRG_MAX_END - KRG_MAX_START)
        params.append((S_wc, krg_max))
    return params


# ── CMG file formatter (same conventions as the original 10-bin file) ────────

def _fmt(v):
    """Format a single float value for a CMG table column."""
    if v == 0.0:
        return f"{'0':>18}"
    elif 0 < abs(v) < 1e-4:
        return f"{v:>18.6E}"
    else:
        return f"{v:>18.10g}"


def _write_rpt(fh, rpt_num, S_wc, krg_max):
    fh.write(f"\nRPT {rpt_num}\n")
    fh.write(f"** S_wc={S_wc:.4f}  krg_max={krg_max:.4f}  krw_max={KRW_MAX:.1f}  S_gr=0.0\n")

    fh.write("**        Sw             krw            krow            Pc\n")
    fh.write("SWT\n")
    for row in zip(*build_swt(S_wc)):
        fh.write(''.join(_fmt(v) for v in row) + '\n')

    fh.write("\n")

    fh.write("**        Sg             krg            krog            Pc\n")
    fh.write("SGT\n")
    for row in zip(*build_sgt(S_wc, krg_max)):
        fh.write(''.join(_fmt(v) for v in row) + '\n')


# ── Verification against the original 10-bin file ────────────────────────────

def _parse_rpt_tables(path, rpt_num):
    """Extract (swt_rows, sgt_rows) as float arrays for one RPT of an inc file."""
    with open(path) as f:
        lines = f.readlines()

    swt, sgt, mode = [], [], None
    in_rpt = False
    for line in lines:
        s = line.strip()
        m = re.match(r'^RPT\s+(\d+)$', s)
        if m:
            in_rpt = int(m.group(1)) == rpt_num
            mode = None
            continue
        if not in_rpt or not s or s.startswith('**'):
            continue
        if s == 'SWT':
            mode = swt
            continue
        if s == 'SGT':
            mode = sgt
            continue
        if mode is not None:
            mode.append([float(v) for v in s.split()])
    return np.array(swt), np.array(sgt)


def verify_against_original():
    """RPT 1/10 of the original must equal our formulas at S_wc/krg_max endpoints."""
    checks = [(1, S_WC_START, KRG_MAX_START), (10, S_WC_END, KRG_MAX_END)]
    for rpt, S_wc, krg_max in checks:
        swt_orig, sgt_orig = _parse_rpt_tables(ORIG_INC, rpt)
        swt_new = np.column_stack(build_swt(S_wc))
        sgt_new = np.column_stack(build_sgt(S_wc, krg_max))
        if not (np.allclose(swt_orig, swt_new, atol=1e-9) and
                np.allclose(sgt_orig, sgt_new, atol=1e-9)):
            sys.exit(f'ERROR: regenerated tables do not match original RPT {rpt} '
                     f'in {ORIG_INC} — model formulas are wrong, aborting.')
    print(f'Verified: formulas reproduce RPT 1 and RPT 10 of {ORIG_INC}')


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    verify_against_original()

    rpt_params = make_rpt_params(N_RPT)

    print(f"\n{'RPT':>4}  {'S_wc':>8}  {'krg_max':>8}")
    print("-" * 26)
    for i, (S_wc, krg_max) in enumerate(rpt_params):
        print(f"{i+1:>4}  {S_wc:>8.4f}  {krg_max:>8.4f}")

    with open(OUT_INC, "w") as fh:
        fh.write("** ============================================================\n")
        fh.write("** Rock-Fluid Properties - 2D CO2-Water Flow\n")
        fh.write("** System  : CO2 (gas/supercritical) + Brine, 2-phase\n")
        fh.write("** Units   : SI  (Pc in kPa, saturations as fractions)\n")
        fh.write("**\n")
        fh.write("** Model   : Corey power-law\n")
        fh.write("**   krw  = Se_w^4,   Se_w = (Sw - S_wc) / (1 - S_wc)\n")
        fh.write("**   krg  = krg_max * Se_g^2,  Se_g = Sg / (1 - S_wc)\n")
        fh.write("**\n")
        fh.write("** Fixed (all RPTs)\n")
        fh.write("**   krw_max = 1.0,  S_gr (residual CO2) = 0.0\n")
        fh.write("**\n")
        fh.write(f"** Varying across RPT 1 -> RPT {N_RPT} (linear interpolation)\n")
        fh.write(f"**   S_wc    : {S_WC_START:.2f} -> {S_WC_END:.2f}\n")
        fh.write(f"**   krg_max : {KRG_MAX_START:.2f} -> {KRG_MAX_END:.2f}\n")
        fh.write("**\n")
        fh.write(f"** Same total change as the 10-bin file, split into {N_RPT} bins\n")
        fh.write("** Generated by script/generate_30_binned_relperm_curve.py\n")
        fh.write("** ============================================================\n")

        for i, (S_wc, krg_max) in enumerate(rpt_params):
            _write_rpt(fh, i + 1, S_wc, krg_max)

    print(f"\nWrote: {OUT_INC}")
