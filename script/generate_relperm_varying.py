"""
Generate cmg_datafile/rock_fluid_properties_varying.inc
with 10 rock types (RPT 1 - RPT 10) using Brooks-Corey model.

RPT 1 : original BC-fitted parameters
RPT 2-10 : linear variation of two parameters across the 9 steps:
    S_wc / S_lr   (connate water / liquid sat) : 0.13 -> 0.03  (-0.10 total)
    k_ep          (endpoint kr: k_rw_max = k_rog_max) : 0.35 -> 0.60  (+0.25 total)

All other BC parameters are held fixed at the RPT-1 fitted values:
    S_or = S_gr = 0.20
    k_row_max = k_rg_max = 1.0
    n_w = n_l = 3.5 ,  n_o = n_g = 1.5
    P_e = 0.3837 ,  lam = 1.5447

Gas-liquid table uses SGT (indexed by Sg, increasing from 0 to Sg_max = 1 - S_lr).
"""

import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_relperm import parse_relperm
from fit_brooks_corey import fit_swt, fit_sgt

# ── Constants ─────────────────────────────────────────────────────────────────

N_RPT    = 10       # total number of rock types
N_MOBILE = 46       # linspace points within mobile range (gives 26 rows total)
PC_MAX   = 110  # Pc cap at zero wetting-phase saturation


# ── Saturation grid ───────────────────────────────────────────────────────────

def _sat_grid(S_wr, S_nwr):
    """
    Build 26-point saturation array:
      [0.0, S_wr, ...22 inner..., 1-S_nwr, 1.0]
    """
    mobile = np.linspace(S_wr, 1.0 - S_nwr, N_MOBILE)
    return np.concatenate([[0.0], mobile, [1.0]])


# ── Brooks-Corey evaluation ───────────────────────────────────────────────────

def _eval_bc(S, S_wr, S_nwr, k_rw, k_rnw, n_w, n_nw, P_e, lam):
    """
    Evaluate kr_wet, kr_nwet, Pc over saturation array S.

    Wetting phase:     kr_wet  = k_rw  * Se ^ n_w
    Non-wetting phase: kr_nwet = k_rnw * (1-Se) ^ n_nw
    Capillary pressure: Pc     = P_e   * Se ^ (-1/lam)

    Edge cases handled:
      S <= 0      : kr_wet=0,    kr_nwet=k_rnw,  Pc=PC_MAX
      S = S_wr    : kr_wet=0,    kr_nwet=k_rnw,  Pc at Se_min (first inner step)
      S >= 1-S_nwr: kr_nwet=0,  Pc=P_e
      S >= 1.0    : kr_wet=1.0  (CMG convention outside two-phase region)
    """
    S    = np.asarray(S, dtype=float)
    span = 1.0 - S_wr - S_nwr
    Se   = np.clip((S - S_wr) / span, 0.0, 1.0)

    # ── relative permeabilities ───────────────────────────────────────────────
    kr_wet  = k_rw  * Se ** n_w
    kr_nwet = k_rnw * (1.0 - Se) ** n_nw

    kr_wet[S  < S_wr]         = 0.0    # wetting phase immobile below S_wr
    kr_nwet[S > 1.0 - S_nwr]  = 0.0   # non-wetting immobile above 1-S_nwr

    # Extend for CMG table endpoints outside two-phase region
    kr_wet[S  >= 1.0] = 1.0    # fully liquid-saturated: kr=1
    kr_nwet[S <= 0.0] = k_rnw  # fully gas-saturated:   kr=k_rnw_max

    # ── capillary pressure ────────────────────────────────────────────────────
    # At Se=0 (S=S_wr), Pc diverges; use Pc at first inner mobile step instead.
    Se_min = 1.0 / (N_MOBILE - 1)          # = 1/23, same for every RPT
    Se_pc  = np.maximum(Se, Se_min)

    Pc                      = P_e * Se_pc ** (-1.0 / lam)
    Pc[S <= 0.0]            = PC_MAX        # hard cap at zero saturation
    Pc[S >= 1.0 - S_nwr]   = P_e           # upper boundary = entry pressure

    return kr_wet, kr_nwet, Pc


# ── Table builders for SWT and SLT ───────────────────────────────────────────

def build_swt(p):
    """Return (Sw, krw, krow, Pc) arrays from SWT BC params dict."""
    Sw = _sat_grid(p["S_wc"], p["S_or"])
    krw, krow, Pc = _eval_bc(
        Sw,
        p["S_wc"], p["S_or"],
        p["k_rw_max"],  p["k_row_max"],
        p["n_w"],       p["n_o"],
        p["P_e_wo"],    p["lam_wo"],
    )
    return Sw, krw, krow, Pc


def build_sgt(p):
    """Return (Sg, krg, krog, Pc) arrays from SGT BC params dict.

    SGT is indexed by gas saturation Sg (increasing from 0 to Sg_max = 1 - S_lr).
    Builds the saturation grid in the Sg domain, evaluates via _eval_bc in the
    equivalent Sl = 1 - Sg domain, then converts back.
    """
    Sg_max = 1.0 - p["S_lr"]
    # Mobile range: S_gr to Sg_max.  If S_gr=0, the linspace already starts at 0.
    mobile = np.linspace(p["S_gr"], Sg_max, N_MOBILE)
    if p["S_gr"] > 1e-9:
        # Add a Sg=0 sentinel below the mobile range (mirrors SWT's Sw=0 sentinel)
        Sg = np.concatenate([[0.0], mobile])
    else:
        Sg = mobile  # mobile starts at 0.0, no separate sentinel needed

    Sl = 1.0 - Sg
    krog, krg, Pc = _eval_bc(
        Sl,
        p["S_lr"],      p["S_gr"],
        p["k_rog_max"], p["k_rg_max"],
        p["n_l"],       p["n_g"],
        p["P_e_gl"],    p["lam_gl"],
    )
    # At Sg=0 there is no gas-oil interface: Pcog must be 0 by convention
    Pc[0] = 0.0
    return Sg, krg, krog, Pc


# ── Parameter variation across RPTs ──────────────────────────────────────────

def make_rpt_params(p_swt_base, p_sgt_base):
    """
    Generate parameter dicts for all N_RPT rock types.

    Linear interpolation from RPT 1 (n=0) to RPT 10 (n=9):
      S_wc / S_lr  : 0.13  -> 0.03   (frac * 0.10 decrease)
      k_ep         : 0.35  -> 0.60   (frac * 0.25 increase)
    """
    rpts = []
    for n in range(N_RPT):
        frac  = n / (N_RPT - 1)    # 0.0 at RPT 1,  1.0 at RPT 10
        S_wc_n = 0.13 - frac * 0.10
        k_ep_n = 0.35 + frac * 0.25

        p_swt_n = dict(p_swt_base)
        p_swt_n["S_wc"]     = S_wc_n
        p_swt_n["k_rw_max"] = k_ep_n

        p_sgt_n = dict(p_sgt_base)
        p_sgt_n["S_lr"]      = S_wc_n
        p_sgt_n["k_rog_max"] = k_ep_n

        rpts.append((p_swt_n, p_sgt_n))
    return rpts


# ── CMG file formatter ────────────────────────────────────────────────────────

def _fmt(v):
    """Format a single float value for CMG table column."""
    if v == 0.0:
        return f"{'0':>18}"
    elif 0 < abs(v) < 1e-4:
        return f"{v:>18.6E}"
    else:
        return f"{v:>18.10g}"


def _write_rpt(fh, rpt_num, p_swt, p_sgt):
    fh.write(f"\nRPT {rpt_num}\n")

    # SWT
    fh.write("**        Sw             krw            krow            Pc\n")
    fh.write("SWT\n")
    Sw, krw, krow, Pc_wo = build_swt(p_swt)
    for s, a, b, c in zip(Sw, krw, krow, Pc_wo):
        fh.write(f"{_fmt(s)}{_fmt(a)}{_fmt(b)}{_fmt(c)}\n")

    fh.write("\n")

    # SGT
    fh.write("**        Sg             krg            krog            Pc\n")
    fh.write("SGT\n")
    Sg, krg, krog, Pc_gl = build_sgt(p_sgt)
    for s, a, b, c in zip(Sg, krg, krog, Pc_gl):
        fh.write(f"{_fmt(s)}{_fmt(a)}{_fmt(b)}{_fmt(c)}\n")


# ── Main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Load original table and fit BC parameters
    swt, sgt = parse_relperm("cmg_datafile/rock_fluid_properties.inc")
    p_swt_base = fit_swt(swt)
    p_sgt_base = fit_sgt(sgt)

    # Build parameter set for all 10 RPTs
    rpt_params = make_rpt_params(p_swt_base, p_sgt_base)

    # Print summary table
    print(f"{'RPT':>4}  {'S_wc/S_lr':>12}  {'k_ep (k_rw_max=k_rog_max)':>26}")
    print("-" * 48)
    for i, (ps, pg) in enumerate(rpt_params):
        print(f"{i+1:>4}  {ps['S_wc']:>12.4f}  {ps['k_rw_max']:>26.4f}")

    # Write CMG .inc file
    out_path = os.path.join("cmg_datafile", "rock_fluid_properties_varying.inc")
    with open(out_path, "w") as fh:
        fh.write("** Rock-fluid properties with varying endpoint kr and residual saturation\n")
        fh.write("** RPT 1 : original BC-fitted curve\n")
        fh.write("** RPT 2-10 : S_wc/S_lr 0.13->0.03, k_rw_max/k_rog_max 0.35->0.60\n")
        fh.write(f"** Generated by generate_relperm_varying.py\n")

        for i, (ps, pg) in enumerate(rpt_params):
            _write_rpt(fh, i + 1, ps, pg)

    print(f"\nWrote: {out_path}")
