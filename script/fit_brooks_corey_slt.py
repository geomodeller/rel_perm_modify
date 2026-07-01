"""
Brooks-Corey (BC) parameterization for SWT and SLT tables (CMG format).

════════════════════════════════════════════════════════════════════════════════
 TABLE LAYOUTS
════════════════════════════════════════════════════════════════════════════════
SWT  (water-oil system)
  swt[:, 0]  Sw    water saturation
  swt[:, 1]  krw   water relative permeability        (wetting phase)
  swt[:, 2]  krow  oil relative permeability           (non-wetting phase)
  swt[:, 3]  Pc_wo water-oil capillary pressure  Pc = Po - Pw

SLT  (gas-liquid system)
  slt[:, 0]  Sl    liquid saturation  (Sl = Sw + So; in CCS with no oil, Sl = Sw)
  slt[:, 1]  krg   gas relative permeability           (non-wetting phase, e.g. scCO2)
  slt[:, 2]  krog  liquid relative permeability        (wetting phase)
  slt[:, 3]  Pc_gl gas-liquid capillary pressure  Pc = Pg - Pl

════════════════════════════════════════════════════════════════════════════════
 BROOKS-COREY MODEL  (same structure for both tables)
════════════════════════════════════════════════════════════════════════════════
Effective (normalized) saturation of the wetting phase:
    Se = (S_w - S_wr) / (1 - S_wr - S_nwr)      Se ∈ [0, 1]

Wetting-phase kr:
    kr_w(Se)  = k_rw_max  * Se ^ n_w

Non-wetting-phase kr:
    kr_nw(Se) = k_rnw_max * (1 - Se) ^ n_nw

Capillary pressure:
    Pc(Se)    = P_e * Se ^ (-1 / lam)

════════════════════════════════════════════════════════════════════════════════
 ALL PARAMETERS
════════════════════════════════════════════════════════════════════════════════
SWT parameters (dict key → meaning):
  S_wc      connate (irreducible) water saturation   — read from table
  S_or      residual oil saturation                  — read from table
  k_rw_max  endpoint water kr  at Sw = 1 - S_or      — read from table
  k_row_max endpoint oil kr    at Sw = S_wc           — read from table
  n_w       water  Corey exponent                    — fitted
  n_o       oil    Corey exponent                    — fitted
  P_e_wo    water-oil capillary entry pressure        — fitted
  lam_wo    pore size distribution index (water-oil)  — fitted

SLT parameters (dict key → meaning):
  S_lr      residual liquid saturation               — read from table
  S_gr      residual gas saturation                  — read from table
  k_rog_max endpoint liquid kr at Sl = 1 - S_gr      — read from table
  k_rg_max  endpoint gas kr   at Sl = S_lr           — read from table
  n_l       liquid Corey exponent                    — fitted
  n_g       gas    Corey exponent                    — fitted
  P_e_gl    gas-liquid capillary entry pressure       — fitted
  lam_gl    pore size distribution index (gas-liquid) — fitted

════════════════════════════════════════════════════════════════════════════════
 HOW TO GET WATER (BRINE) AND CO2 kr AND Pc IN A CCS CONTEXT
════════════════════════════════════════════════════════════════════════════════
In CCS (CO2 storage into a saline aquifer) there are only two phases:
  liquid = brine (water)    →  Sw + So = Sl = Sw  (So = 0)
  gas    = scCO2

  Brine kr    : krw(Sw)          from SWT table, column 1
  CO2   kr    : krg(Sl = Sw)     from SLT table, column 1
                  (Sl equals Sw because there is no oil)
  CO2-brine Pc: Pc_gl(Sl = Sw)  from SLT table, column 3
                  (= Pco2 − Pbrine)

The SWT Pc_wo (= Po − Pw) is the water-oil capillary pressure; in a
two-phase CO2-brine system it is NOT used directly — only Pc_gl matters.

In a three-phase system (water + oil + CO2):
  krw        : krw(Sw)            from SWT  (unchanged)
  krg        : krg(Sl = Sw + So)  from SLT  (liquid includes both water and oil)
  kro        : Stone's model combining krow(Sw) from SWT and krog(Sl) from SLT
  Pc_wo      : Po − Pw            from SWT
  Pc_gl      : Pg − Pl = Pg − Po from SLT  (gas-oil capillary pressure)
  Pc_gw      : Pg − Pw = Pc_gl + Pc_wo

Summary table:
  Property          Two-phase (CO2-brine)         Three-phase
  ─────────────────────────────────────────────────────────────
  Brine kr          krw(Sw)    from SWT            krw(Sw) from SWT
  CO2 kr            krg(Sl=Sw) from SLT            krg(Sl=Sw+So) from SLT
  Oil kr            —                              Stone: f(krow, krog)
  CO2-brine Pc      Pc_gl(Sl)  from SLT            Pc_gw = Pc_gl + Pc_wo
"""

import sys
import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from read_relperm import parse_relperm


# ── Brooks-Corey functional forms ─────────────────────────────────────────────

def _bc_kr_wet(Se, n):
    return Se ** n

def _bc_kr_nwet(Se, n):
    return (1.0 - Se) ** n

def _bc_pc(Se, P_e, lam):
    return P_e * Se ** (-1.0 / lam)


# ── Generic two-phase BC fitter ───────────────────────────────────────────────

def _fit_bc_2phase(S, kr_wet, kr_nwet, Pc):
    """
    Generic Brooks-Corey fit for any two-phase table.

    S       : saturation of the wetting phase (Sw for SWT, Sl for SLT)
    kr_wet  : wetting-phase kr  (krw for SWT, krog for SLT)
    kr_nwet : non-wetting kr    (krow for SWT, krg  for SLT)
    Pc      : capillary pressure data

    Returns dict with generic keys:
      S_wr, S_nwr, k_rw_max, k_rnw_max, n_w, n_nw, P_e, lam
    """
    # ── Residual saturations from the table ───────────────────────────────────
    # S_wr : highest wetting-phase S where kr_wet is still zero
    S_wr = float(S[kr_wet == 0.0].max())

    # S_nwr : 1 - (lowest S where kr_nwet is zero, on the high-S side)
    nwet_zero_high = S[(kr_nwet == 0.0) & (S > 0.5)]
    S_nwr = float(1.0 - nwet_zero_high.min())

    # ── Endpoint kr values ────────────────────────────────────────────────────
    idx_wr  = np.argmin(np.abs(S - S_wr))
    idx_nwr = np.argmin(np.abs(S - (1.0 - S_nwr)))
    k_rnw_max = float(kr_nwet[idx_wr])   # non-wetting kr at S = S_wr
    k_rw_max  = float(kr_wet[idx_nwr])   # wetting kr   at S = 1 - S_nwr

    # ── Normalized saturation ─────────────────────────────────────────────────
    Se = np.clip((S - S_wr) / (1.0 - S_wr - S_nwr), 0.0, 1.0)

    # ── Fit wetting kr (n_w) ──────────────────────────────────────────────────
    mask_w = (kr_wet > 0.0) & (S <= 1.0 - S_nwr)
    (n_w,), _ = curve_fit(
        _bc_kr_wet, Se[mask_w], kr_wet[mask_w] / k_rw_max,
        p0=[2.0], bounds=(0.1, 20.0),
    )

    # ── Fit non-wetting kr (n_nw) ─────────────────────────────────────────────
    mask_nw = (kr_nwet > 0.0) & (S <= 1.0 - S_nwr)
    (n_nw,), _ = curve_fit(
        _bc_kr_nwet, Se[mask_nw], kr_nwet[mask_nw] / k_rnw_max,
        p0=[2.0], bounds=(0.1, 20.0),
    )

    # ── Fit Pc (P_e, lam) — mobile range, Se > 0 ─────────────────────────────
    mask_pc = (S > S_wr) & (S <= 1.0 - S_nwr)
    (P_e, lam), _ = curve_fit(
        _bc_pc, Se[mask_pc], Pc[mask_pc],
        p0=[0.5, 0.5], bounds=([1e-6, 0.01], [1e4, 20.0]),
    )

    return {
        "S_wr":     S_wr,
        "S_nwr":    S_nwr,
        "k_rw_max": k_rw_max,
        "k_rnw_max":k_rnw_max,
        "n_w":      float(n_w),
        "n_nw":     float(n_nw),
        "P_e":      float(P_e),
        "lam":      float(lam),
    }


# ── Table-specific fitting wrappers ──────────────────────────────────────────

def fit_swt(swt):
    """
    Fit BC to SWT table.  Returns SWT-specific parameter dict.
    SWT columns: Sw, krw, krow, Pc_wo
    """
    Sw   = swt[:, 0]
    krw  = swt[:, 1]   # wetting phase
    krow = swt[:, 2]   # non-wetting phase
    Pc   = swt[:, 3]

    g = _fit_bc_2phase(Sw, krw, krow, Pc)

    return {
        "S_wc":     g["S_wr"],      # connate water saturation
        "S_or":     g["S_nwr"],     # residual oil saturation
        "k_rw_max": g["k_rw_max"],  # max water kr at Sw = 1 - S_or
        "k_row_max":g["k_rnw_max"], # max oil kr   at Sw = S_wc
        "n_w":      g["n_w"],       # water Corey exponent
        "n_o":      g["n_nw"],      # oil   Corey exponent
        "P_e_wo":   g["P_e"],       # water-oil entry pressure
        "lam_wo":   g["lam"],       # pore size distribution index (water-oil)
    }


def fit_slt(slt):
    """
    Fit BC to SLT table.  Returns SLT-specific parameter dict.
    SLT columns: Sl, krg, krog, Pc_gl
    """
    Sl   = slt[:, 0]
    krg  = slt[:, 1]   # non-wetting phase (gas/CO2)
    krog = slt[:, 2]   # wetting phase (liquid)
    Pc   = slt[:, 3]

    g = _fit_bc_2phase(Sl, krog, krg, Pc)   # note: swap so kr_wet=krog, kr_nwet=krg

    return {
        "S_lr":     g["S_wr"],      # residual liquid saturation
        "S_gr":     g["S_nwr"],     # residual gas saturation
        "k_rog_max":g["k_rw_max"],  # max liquid kr at Sl = 1 - S_gr
        "k_rg_max": g["k_rnw_max"], # max gas kr   at Sl = S_lr
        "n_l":      g["n_w"],       # liquid Corey exponent
        "n_g":      g["n_nw"],      # gas    Corey exponent
        "P_e_gl":   g["P_e"],       # gas-liquid entry pressure
        "lam_gl":   g["lam"],       # pore size distribution index (gas-liquid)
    }


def fit_sgt(sgt):
    """
    Fit BC to SGT table.  Returns same parameter dict keys as fit_slt.
    SGT columns: Sg, krg, krog, Pc_og  (Sg increasing from 0 to Sg_max)

    Converts to liquid-saturation domain (Sl = 1 - Sg, reversed arrays) so
    _fit_bc_2phase can be reused with wetting=liquid(krog), non-wetting=gas(krg).
    """
    Sg   = sgt[:, 0]
    krg  = sgt[:, 1]
    krog = sgt[:, 2]
    Pc   = sgt[:, 3]

    Sl   = (1.0 - Sg)[::-1].copy()
    krog = krog[::-1].copy()
    krg  = krg[::-1].copy()
    Pc   = Pc[::-1].copy()

    g = _fit_bc_2phase(Sl, krog, krg, Pc)

    return {
        "S_lr":     g["S_wr"],
        "S_gr":     g["S_nwr"],
        "k_rog_max":g["k_rw_max"],
        "k_rg_max": g["k_rnw_max"],
        "n_l":      g["n_w"],
        "n_g":      g["n_nw"],
        "P_e_gl":   g["P_e"],
        "lam_gl":   g["lam"],
    }


# ── Prediction ────────────────────────────────────────────────────────────────

def predict_swt(Sw, p):
    """Evaluate BC kr and Pc for SWT.  Returns (krw, krow, Pc_wo)."""
    Se = np.clip((Sw - p["S_wc"]) / (1.0 - p["S_wc"] - p["S_or"]), 0.0, 1.0)

    krw  = np.where(Sw >= p["S_wc"],         p["k_rw_max"]  * Se              ** p["n_w"], 0.0)
    krow = np.where(Sw <= 1.0 - p["S_or"],   p["k_row_max"] * (1.0 - Se)      ** p["n_o"], 0.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        Pc_wo = np.where(
            (Sw > p["S_wc"]) & (Sw <= 1.0 - p["S_or"]),
            p["P_e_wo"] * Se ** (-1.0 / p["lam_wo"]),
            np.nan,
        )
    return krw, krow, Pc_wo


def predict_slt(Sl, p):
    """Evaluate BC kr and Pc for SLT.  Returns (krg, krog, Pc_gl)."""
    Se = np.clip((Sl - p["S_lr"]) / (1.0 - p["S_lr"] - p["S_gr"]), 0.0, 1.0)

    krog = np.where(Sl >= p["S_lr"],         p["k_rog_max"] * Se              ** p["n_l"], 0.0)
    krg  = np.where(Sl <= 1.0 - p["S_gr"],   p["k_rg_max"]  * (1.0 - Se)     ** p["n_g"], 0.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        Pc_gl = np.where(
            (Sl > p["S_lr"]) & (Sl <= 1.0 - p["S_gr"]),
            p["P_e_gl"] * Se ** (-1.0 / p["lam_gl"]),
            np.nan,
        )
    return krg, krog, Pc_gl


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_all(swt, slt, p_swt, p_slt, save_path="brooks_corey_fit.png"):
    S_fine = np.linspace(0.0, 1.0, 400)

    krw_bc, krow_bc, Pcwo_bc = predict_swt(S_fine, p_swt)
    krg_bc, krog_bc, Pcgl_bc = predict_slt(S_fine, p_slt)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Brooks-Corey fit  —  SWT (top) and SLT (bottom)", fontsize=13)

    # ── Row 0: SWT ────────────────────────────────────────────────────────────
    ax = axes[0, 0]
    ax.plot(swt[:, 0], swt[:, 1], "o", color="steelblue", ms=5, label="Data $k_{rw}$")
    ax.plot(S_fine, krw_bc, "-", color="steelblue",
            label=f"B-C  $n_w={p_swt['n_w']:.3f}$")
    ax.set_title("SWT — Water kr"); ax.set_xlabel("$S_w$"); ax.set_ylabel("$k_{rw}$")
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(swt[:, 0], swt[:, 2], "s", color="tomato", ms=5, label="Data $k_{row}$")
    ax.plot(S_fine, krow_bc, "-", color="tomato",
            label=f"B-C  $n_o={p_swt['n_o']:.3f}$")
    ax.set_title("SWT — Oil kr"); ax.set_xlabel("$S_w$"); ax.set_ylabel("$k_{row}$")
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    m = (swt[:, 0] > p_swt["S_wc"]) & (swt[:, 0] <= 1.0 - p_swt["S_or"])
    ax.plot(swt[m, 0], swt[m, 3], "^", color="seagreen", ms=5, label="Data $P_{c,wo}$")
    ax.plot(S_fine, Pcwo_bc, "-", color="seagreen",
            label=f"B-C  $P_e={p_swt['P_e_wo']:.3f}$,  $\\lambda={p_swt['lam_wo']:.3f}$")
    ax.set_title("SWT — Water-oil $P_c$"); ax.set_xlabel("$S_w$"); ax.set_ylabel("$P_c$")
    ax.set_xlim(0, 1); ax.legend(); ax.grid(True, alpha=0.3)

    # ── Row 1: SLT ────────────────────────────────────────────────────────────
    ax = axes[1, 0]
    ax.plot(slt[:, 0], slt[:, 1], "o", color="darkorange", ms=5, label="Data $k_{rg}$")
    ax.plot(S_fine, krg_bc, "-", color="darkorange",
            label=f"B-C  $n_g={p_slt['n_g']:.3f}$")
    ax.set_title("SLT — CO$_2$ (gas) kr"); ax.set_xlabel("$S_l$"); ax.set_ylabel("$k_{rg}$")
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(slt[:, 0], slt[:, 2], "s", color="mediumpurple", ms=5, label="Data $k_{rog}$")
    ax.plot(S_fine, krog_bc, "-", color="mediumpurple",
            label=f"B-C  $n_l={p_slt['n_l']:.3f}$")
    ax.set_title("SLT — Liquid kr"); ax.set_xlabel("$S_l$"); ax.set_ylabel("$k_{rog}$")
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    m = (slt[:, 0] > p_slt["S_lr"]) & (slt[:, 0] <= 1.0 - p_slt["S_gr"])
    ax.plot(slt[m, 0], slt[m, 3], "^", color="teal", ms=5, label="Data $P_{c,gl}$")
    ax.plot(S_fine, Pcgl_bc, "-", color="teal",
            label=f"B-C  $P_e={p_slt['P_e_gl']:.3f}$,  $\\lambda={p_slt['lam_gl']:.3f}$")
    ax.set_title("SLT — Gas-liquid $P_c$"); ax.set_xlabel("$S_l$"); ax.set_ylabel("$P_c$")
    ax.set_xlim(0, 1); ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.show()
    print(f"Plot saved: {save_path}")


# ── Parameter printer ─────────────────────────────────────────────────────────

def print_params(p_swt, p_slt):
    print("=" * 62)
    print(" Brooks-Corey Parameters - SWT (water-oil)")
    print("=" * 62)
    print("  Saturation endpoints  (read from table):")
    print(f"    S_wc      connate water saturation          = {p_swt['S_wc']:.4f}")
    print(f"    S_or      residual oil saturation           = {p_swt['S_or']:.4f}")
    print("  Kr endpoints  (read from table):")
    print(f"    k_rw_max  max water kr  at Sw = 1 - S_or   = {p_swt['k_rw_max']:.4f}")
    print(f"    k_row_max max oil kr    at Sw = S_wc        = {p_swt['k_row_max']:.4f}")
    print("  Corey exponents  (fitted):")
    print(f"    n_w       water Corey exponent              = {p_swt['n_w']:.4f}")
    print(f"    n_o       oil   Corey exponent              = {p_swt['n_o']:.4f}")
    print("  Capillary pressure  (fitted):")
    print(f"    P_e_wo    water-oil entry pressure          = {p_swt['P_e_wo']:.4f}")
    print(f"    lam_wo    pore size distribution index      = {p_swt['lam_wo']:.4f}")

    print()
    print("=" * 62)
    print(" Brooks-Corey Parameters - SLT (gas-liquid)")
    print("=" * 62)
    print("  Saturation endpoints  (read from table):")
    print(f"    S_lr      residual liquid saturation        = {p_slt['S_lr']:.4f}")
    print(f"    S_gr      residual gas saturation           = {p_slt['S_gr']:.4f}")
    print("  Kr endpoints  (read from table):")
    print(f"    k_rog_max max liquid kr at Sl = 1 - S_gr   = {p_slt['k_rog_max']:.4f}")
    print(f"    k_rg_max  max gas kr    at Sl = S_lr        = {p_slt['k_rg_max']:.4f}")
    print("  Corey exponents  (fitted):")
    print(f"    n_l       liquid Corey exponent             = {p_slt['n_l']:.4f}")
    print(f"    n_g       gas    Corey exponent             = {p_slt['n_g']:.4f}")
    print("  Capillary pressure  (fitted):")
    print(f"    P_e_gl    gas-liquid entry pressure         = {p_slt['P_e_gl']:.4f}")
    print(f"    lam_gl    pore size distribution index      = {p_slt['lam_gl']:.4f}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    filepath = "cmg_datafile/rock_fluid_properties.inc"
    swt, sgt = parse_relperm(filepath)

    p_swt = fit_swt(swt)
    p_sgt = fit_sgt(sgt)

    print_params(p_swt, p_sgt)
    # plot_all expects SLT (Sl-indexed) data; convert sgt for plotting
    slt_equiv = np.column_stack([
        1.0 - sgt[:, 0],   # Sl = 1 - Sg
        sgt[:, 1],          # krg
        sgt[:, 2],          # krog
        sgt[:, 3],          # Pc
    ])[::-1]
    plot_all(swt, slt_equiv, p_swt, p_sgt)
