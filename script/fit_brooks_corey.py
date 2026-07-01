"""
Brooks-Corey (BC) parameterization for SWT and SGT tables (CMG format).

════════════════════════════════════════════════════════════════════════════════
 TABLE LAYOUTS
════════════════════════════════════════════════════════════════════════════════
SWT  (water-oil system)
  swt[:, 0]  Sw    water saturation
  swt[:, 1]  krw   water relative permeability        (wetting phase)
  swt[:, 2]  krow  oil relative permeability           (non-wetting phase)
  swt[:, 3]  Pc_wo water-oil capillary pressure  Pc = Po - Pw

SGT  (gas-liquid system)
  sgt[:, 0]  Sg    gas saturation   (increasing from 0 to Sg_max = 1 - S_lr)
  sgt[:, 1]  krg   gas relative permeability           (non-wetting phase)
  sgt[:, 2]  krog  liquid relative permeability        (wetting phase)
  sgt[:, 3]  Pc_og gas-liquid capillary pressure  Pc = Pg - Po
               Pc_og increases with Sg; equals 0 at Sg = 0

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

For SGT, Se is computed from the equivalent liquid saturation Sl = 1 - Sg.

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

SGT parameters (dict key → meaning):
  S_lr      connate liquid saturation  (= 1 - Sg_max) — read from table
  S_gr      critical gas saturation                  — read from table
  k_rog_max endpoint liquid kr at Sg = S_gr           — read from table
  k_rg_max  endpoint gas kr   at Sg = Sg_max          — read from table
  n_l       liquid Corey exponent                    — fitted
  n_g       gas    Corey exponent                    — fitted
  P_e_gl    gas-liquid capillary entry pressure       — fitted
  lam_gl    pore size distribution index (gas-liquid) — fitted

════════════════════════════════════════════════════════════════════════════════
 HOW TO GET BRINE AND CO2/CH4 kr AND Pc IN A CO2 + CH4 + BRINE CCS CONTEXT
════════════════════════════════════════════════════════════════════════════════
GEM assigns kr to phases, not components:
  liquid phase = brine (water + dissolved CO2)
  gas phase    = CH4-rich (+ trace CO2)
  oil phase    = CO2-rich dense phase

Two-phase (hydrocarbon single-phase + brine):
  Brine  kr    : krw(Sw)          from SWT
  Gas/CO2 kr   : krg(Sg = 1-Sw)  from SGT   (Sg = 1 - Sl = 1 - Sw when no oil)
  Pc_og        : Pc_og(Sg)        from SGT

Three-phase (CO2-dense oil + CH4 gas + brine):
  Brine  kr    : krw(Sw)                   from SWT
  CH4 gas kr   : krg(Sg)                   from SGT
  CO2 oil kr   : kro via Stone's 2nd model  using krow(Sw) + krog(Sg)
  Pc_wo        : Po - Pw                   from SWT
  Pc_og        : Pg - Po                   from SGT
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

    S       : wetting-phase saturation  (Sw for SWT; Sl = 1-Sg reversed for SGT)
    kr_wet  : wetting-phase kr          (krw for SWT; krog for SGT)
    kr_nwet : non-wetting kr            (krow for SWT; krg for SGT)
    Pc      : capillary pressure data

    Returns dict with generic keys:
      S_wr, S_nwr, k_rw_max, k_rnw_max, n_w, n_nw, P_e, lam
    """
    # ── Residual saturations from the table ───────────────────────────────────
    S_wr = float(S[kr_wet == 0.0].max())

    nwet_zero_high = S[(kr_nwet == 0.0) & (S > 0.5)]
    S_nwr = float(1.0 - nwet_zero_high.min())

    # ── Endpoint kr values ────────────────────────────────────────────────────
    idx_wr  = np.argmin(np.abs(S - S_wr))
    idx_nwr = np.argmin(np.abs(S - (1.0 - S_nwr)))
    k_rnw_max = float(kr_nwet[idx_wr])
    k_rw_max  = float(kr_wet[idx_nwr])

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

    # ── Fit Pc (P_e, lam) — mobile range only, exclude Pc=0 boundary rows ────
    # Pc=0 occurs at the wetting-phase sentinel (e.g. Sg=0 for SGT) where the
    # BC power-law diverges rather than returning 0; excluding avoids bias.
    mask_pc = (S > S_wr) & (S <= 1.0 - S_nwr) & (Pc > 0.0)
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
        "S_wc":     g["S_wr"],
        "S_or":     g["S_nwr"],
        "k_rw_max": g["k_rw_max"],
        "k_row_max":g["k_rnw_max"],
        "n_w":      g["n_w"],
        "n_o":      g["n_nw"],
        "P_e_wo":   g["P_e"],
        "lam_wo":   g["lam"],
    }


def fit_sgt(sgt):
    """
    Fit BC to SGT table.  Returns SGT-specific parameter dict.
    SGT columns: Sg, krg, krog, Pc_og  (Sg increasing from 0 to Sg_max = 1-S_lr)

    Internally converts to liquid-saturation domain (Sl = 1-Sg, reversed) so
    _fit_bc_2phase can be reused: wetting = liquid (krog), non-wetting = gas (krg).
    """
    Sg   = sgt[:, 0]
    krg  = sgt[:, 1]
    krog = sgt[:, 2]
    Pc   = sgt[:, 3]

    # Reverse to Sl-domain (increasing Sl = decreasing Sg)
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

    krw  = np.where(Sw >= p["S_wc"],       p["k_rw_max"]  * Se         ** p["n_w"], 0.0)
    krow = np.where(Sw <= 1.0 - p["S_or"], p["k_row_max"] * (1.0 - Se) ** p["n_o"], 0.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        Pc_wo = np.where(
            (Sw > p["S_wc"]) & (Sw <= 1.0 - p["S_or"]),
            p["P_e_wo"] * Se ** (-1.0 / p["lam_wo"]),
            np.nan,
        )
    return krw, krow, Pc_wo


def predict_sgt(Sg, p):
    """Evaluate BC kr and Pc for SGT.  Returns (krg, krog, Pc_og).

    Converts internally to Sl = 1 - Sg for the BC formulas, then maps back.
    Pc_og is set to 0 at Sg = 0 (no gas-oil interface).
    """
    Sl = 1.0 - Sg
    Se = np.clip((Sl - p["S_lr"]) / (1.0 - p["S_lr"] - p["S_gr"]), 0.0, 1.0)

    krog = np.where(Sl >= p["S_lr"],       p["k_rog_max"] * Se         ** p["n_l"], 0.0)
    krg  = np.where(Sl <= 1.0 - p["S_gr"], p["k_rg_max"]  * (1.0 - Se) ** p["n_g"], 0.0)

    with np.errstate(divide="ignore", invalid="ignore"):
        Pc_og = np.where(
            (Sl > p["S_lr"]) & (Sl <= 1.0 - p["S_gr"]) & (Se > 0),
            p["P_e_gl"] * Se ** (-1.0 / p["lam_gl"]),
            np.nan,
        )
    Pc_og = np.where(Sg <= 0.0, 0.0, Pc_og)  # Pcog = 0 at Sg = 0
    return krg, krog, Pc_og


# ── Plotting ──────────────────────────────────────────────────────────────────

def plot_all(swt, sgt, p_swt, p_sgt, save_path="brooks_corey_fit.png"):
    Sw_fine = np.linspace(0.0, 1.0, 400)
    Sg_fine = np.linspace(0.0, 1.0 - p_sgt["S_lr"], 400)

    krw_bc, krow_bc, Pcwo_bc   = predict_swt(Sw_fine, p_swt)
    krg_bc, krog_bc, Pcog_bc   = predict_sgt(Sg_fine, p_sgt)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Brooks-Corey fit  —  SWT (top) and SGT (bottom)", fontsize=13)

    # ── Row 0: SWT ────────────────────────────────────────────────────────────
    ax = axes[0, 0]
    ax.plot(swt[:, 0], swt[:, 1], "o", color="steelblue", ms=5, label="Data $k_{rw}$")
    ax.plot(Sw_fine, krw_bc, "-", color="steelblue",
            label=f"B-C  $n_w={p_swt['n_w']:.3f}$")
    ax.set_title("SWT — Water kr"); ax.set_xlabel("$S_w$"); ax.set_ylabel("$k_{rw}$")
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(swt[:, 0], swt[:, 2], "s", color="tomato", ms=5, label="Data $k_{row}$")
    ax.plot(Sw_fine, krow_bc, "-", color="tomato",
            label=f"B-C  $n_o={p_swt['n_o']:.3f}$")
    ax.set_title("SWT — Oil kr"); ax.set_xlabel("$S_w$"); ax.set_ylabel("$k_{row}$")
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[0, 2]
    m = (swt[:, 0] > p_swt["S_wc"]) & (swt[:, 0] <= 1.0 - p_swt["S_or"])
    ax.plot(swt[m, 0], swt[m, 3], "^", color="seagreen", ms=5, label="Data $P_{c,wo}$")
    ax.plot(Sw_fine, Pcwo_bc, "-", color="seagreen",
            label=f"B-C  $P_e={p_swt['P_e_wo']:.3f}$,  $\\lambda={p_swt['lam_wo']:.3f}$")
    ax.set_title("SWT — Water-oil $P_c$"); ax.set_xlabel("$S_w$"); ax.set_ylabel("$P_c$ (kPa)")
    ax.set_xlim(0, 1); ax.legend(); ax.grid(True, alpha=0.3)

    # ── Row 1: SGT ────────────────────────────────────────────────────────────
    Sg_max = 1.0 - p_sgt["S_lr"]

    ax = axes[1, 0]
    ax.plot(sgt[:, 0], sgt[:, 1], "o", color="darkorange", ms=5, label="Data $k_{rg}$")
    ax.plot(Sg_fine, krg_bc, "-", color="darkorange",
            label=f"B-C  $n_g={p_sgt['n_g']:.3f}$")
    ax.set_title("SGT — Gas kr"); ax.set_xlabel("$S_g$"); ax.set_ylabel("$k_{rg}$")
    ax.set_xlim(0, Sg_max); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    ax.plot(sgt[:, 0], sgt[:, 2], "s", color="mediumpurple", ms=5, label="Data $k_{rog}$")
    ax.plot(Sg_fine, krog_bc, "-", color="mediumpurple",
            label=f"B-C  $n_l={p_sgt['n_l']:.3f}$")
    ax.set_title("SGT — Liquid kr"); ax.set_xlabel("$S_g$"); ax.set_ylabel("$k_{rog}$")
    ax.set_xlim(0, Sg_max); ax.set_ylim(-0.02, 1.05); ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1, 2]
    m = sgt[:, 0] > 0
    ax.plot(sgt[m, 0], sgt[m, 3], "^", color="teal", ms=5, label="Data $P_{c,og}$")
    ax.plot(Sg_fine, Pcog_bc, "-", color="teal",
            label=f"B-C  $P_e={p_sgt['P_e_gl']:.3f}$,  $\\lambda={p_sgt['lam_gl']:.3f}$")
    ax.set_title("SGT — Gas-liquid $P_c$"); ax.set_xlabel("$S_g$"); ax.set_ylabel("$P_c$ (kPa)")
    ax.set_xlim(0, Sg_max); ax.legend(); ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.show()
    print(f"Plot saved: {save_path}")


# ── Parameter printer ─────────────────────────────────────────────────────────

def print_params(p_swt, p_sgt):
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
    print(" Brooks-Corey Parameters - SGT (gas-liquid)")
    print("=" * 62)
    print("  Saturation endpoints  (read from table):")
    print(f"    S_lr      connate liquid saturation         = {p_sgt['S_lr']:.4f}")
    print(f"    S_gr      critical gas saturation           = {p_sgt['S_gr']:.4f}")
    print("  Kr endpoints  (read from table):")
    print(f"    k_rog_max max liquid kr at Sg = S_gr        = {p_sgt['k_rog_max']:.4f}")
    print(f"    k_rg_max  max gas kr    at Sg = Sg_max      = {p_sgt['k_rg_max']:.4f}")
    print("  Corey exponents  (fitted):")
    print(f"    n_l       liquid Corey exponent             = {p_sgt['n_l']:.4f}")
    print(f"    n_g       gas    Corey exponent             = {p_sgt['n_g']:.4f}")
    print("  Capillary pressure  (fitted):")
    print(f"    P_e_gl    gas-liquid entry pressure         = {p_sgt['P_e_gl']:.4f}")
    print(f"    lam_gl    pore size distribution index      = {p_sgt['lam_gl']:.4f}")


# ── Entry point ───────────────────────────────────────────────────────────────

if __name__ == "__main__":
    filepath = "cmg_datafile/rock_fluid_properties.inc"
    swt, sgt = parse_relperm(filepath)

    p_swt = fit_swt(swt)
    p_sgt = fit_sgt(sgt)

    print_params(p_swt, p_sgt)
    plot_all(swt, sgt, p_swt, p_sgt)
