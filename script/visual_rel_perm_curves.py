"""
Visualize relative permeability curves from a CMG GEM .inc file.

Usage
-----
  python script/visual_rel_perm_curves.py [inc_file] [--out PNG]

  inc_file  Path to the .inc file (default: cmg_2D_co2_flow/rock_fluid_properties.inc)
  --out     Output PNG path            (default: <inc_file_dir>/rel_perm_curves.png)

The script parses all RPT blocks in the file and produces a single merged figure:
  krw (from SWT) and krg (from SGT, converted via Sw = 1 - Sg) vs Sw.
  X-axis spans from the minimum residual water saturation across all RPTs to 1.0.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm


# ── Parser ────────────────────────────────────────────────────────────────────

def _is_data_row(line):
    """True if the line contains only numbers (a table data row)."""
    stripped = line.strip()
    if not stripped or stripped.startswith("**"):
        return False
    try:
        [float(x) for x in stripped.split()]
        return True
    except ValueError:
        return False


def parse_inc(filepath):
    """
    Parse a CMG rock-fluid .inc file.

    Returns
    -------
    list of dict with keys:
        rpt   : int
        swt   : ndarray (N, 4)  columns: Sw, krw, krow, Pc
        sgt   : ndarray (M, 4)  columns: Sg, krg, krog, Pc
    """
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    rpts = []
    current_rpt = None
    current_table = None   # "SWT" or "SGT"
    collecting = False

    for line in lines:
        stripped = line.strip()

        # New RPT block
        if stripped.startswith("RPT ") and stripped[4:].strip().isdigit():
            if current_rpt is not None:
                rpts.append(current_rpt)
            current_rpt = {"rpt": int(stripped[4:].strip()), "swt": [], "sgt": []}
            current_table = None
            collecting = False
            continue

        if current_rpt is None:
            continue

        # Table keywords
        if stripped == "SWT":
            current_table = "swt"
            collecting = True
            continue
        if stripped == "SGT":
            current_table = "sgt"
            collecting = True
            continue

        # Stop collecting on a new keyword (non-comment, non-numeric, non-blank)
        if collecting and stripped and not stripped.startswith("**"):
            if not _is_data_row(line):
                collecting = False
                current_table = None

        # Collect data rows
        if collecting and current_table and _is_data_row(line):
            row = [float(x) for x in stripped.split()]
            current_rpt[current_table].append(row)

    if current_rpt is not None:
        rpts.append(current_rpt)

    for r in rpts:
        r["swt"] = np.array(r["swt"]) if r["swt"] else np.empty((0, 4))
        r["sgt"] = np.array(r["sgt"]) if r["sgt"] else np.empty((0, 4))

    return rpts


# ── Plot ──────────────────────────────────────────────────────────────────────

def plot_relperm(rpts, out_png):
    n = len(rpts)
    colors = cm.viridis(np.linspace(0, 1, n))

    fig, ax = plt.subplots(figsize=(8, 5.5))

    # Determine x-axis lower bound: minimum S_wc across all RPTs
    swc_values = []
    for r in rpts:
        swt = r["swt"]
        if swt.size:
            # S_wc = first Sw where krw is still 0 (last zero-krw row)
            zero_mask = swt[:, 1] == 0.0
            if zero_mask.any():
                swc_values.append(swt[zero_mask, 0].max())
    x_min = min(swc_values) if swc_values else 0.0

    for i, r in enumerate(rpts):
        col   = colors[i]
        label = f"RPT {r['rpt']}" if i in (0, n - 1) else None
        lw    = 2.0 if i in (0, n - 1) else 1.0

        # krw vs Sw — only mobile region (Sw >= S_wc)
        swt = r["swt"]
        if swt.size:
            mask = swt[:, 0] >= (swt[swt[:, 1] == 0.0, 0].max() if np.any(swt[:, 1] == 0.0) else 0.0)
            Sw  = swt[mask, 0]
            krw = swt[mask, 1]
            ax.plot(Sw, krw, color=col, lw=lw, ls="-", label=label)

        # krg vs Sw — convert SGT via Sw = 1 - Sg, then restrict to Sw >= S_wc
        sgt = r["sgt"]
        if sgt.size:
            Sw_g = 1.0 - sgt[:, 0]   # reverse mapping
            krg  = sgt[:, 1]
            # Sort by Sw ascending (Sg is ascending, so Sw_g is descending — reverse)
            order = np.argsort(Sw_g)
            Sw_g, krg = Sw_g[order], krg[order]
            swc = swt[swt[:, 1] == 0.0, 0].max() if (swt.size and np.any(swt[:, 1] == 0.0)) else 0.0
            mask = Sw_g >= swc
            ax.plot(Sw_g[mask], krg[mask], color=col, lw=lw, ls="--")

    # ── Axis formatting ────────────────────────────────────────────────────────
    ax.set_xlabel("Water Saturation  $S_w$", fontsize=12)
    ax.set_ylabel("Relative Permeability", fontsize=12)
    ax.set_title("$k_{rw}$ and $k_{rg}$ vs. $S_w$  (all RPTs)", fontsize=13)
    ax.set_xlim(x_min, 1.0)
    ax.set_ylim(-0.01, 1.05)
    ax.grid(True, alpha=0.3)

    # ── Legend: line style + RPT endpoints ────────────────────────────────────
    style_handles = [
        plt.Line2D([0], [0], color="grey", lw=2, ls="-",  label="$k_{rw}$  (solid)"),
        plt.Line2D([0], [0], color="grey", lw=2, ls="--", label="$k_{rg}$  (dashed)"),
    ]
    rpt_handles = [
        plt.Line2D([0], [0], color=colors[0],   lw=2, label=f"RPT {rpts[0]['rpt']}"),
        plt.Line2D([0], [0], color=colors[n-1], lw=2, label=f"RPT {rpts[-1]['rpt']}"),
    ]
    ax.legend(handles=style_handles + rpt_handles, fontsize=10, loc="center left")

    # ── Colorbar ──────────────────────────────────────────────────────────────
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(1, n))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("RPT number", fontsize=11)
    cbar.set_ticks(np.linspace(1, n, min(n, 10), dtype=int))

    # Pull S_wc and krg_max from first/last RPT for subtitle
    swc_1  = rpts[0]["swt"][rpts[0]["swt"][:, 1] == 0.0, 0].max()
    swc_10 = rpts[-1]["swt"][rpts[-1]["swt"][:, 1] == 0.0, 0].max()
    krg_1  = rpts[0]["sgt"][-1, 1]
    krg_10 = rpts[-1]["sgt"][-1, 1]
    fig.suptitle(
        f"Relative Permeability  ({n} RPTs)\n"
        f"$S_{{wc}}$: {swc_1:.2f} → {swc_10:.2f}   |   "
        f"$k_{{rg,max}}$: {krg_1:.3f} → {krg_10:.3f}",
        fontsize=12,
    )
    plt.tight_layout()
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_png}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "inc_file", nargs="?",
        default="cmg_2D_co2_flow/rock_fluid_properties.inc",
        help="Path to the CMG .inc file (default: cmg_2D_co2_flow/rock_fluid_properties.inc)",
    )
    p.add_argument(
        "--out", default=None, metavar="PNG",
        help="Output PNG path (default: <inc_dir>/rel_perm_curves.png)",
    )
    args = p.parse_args()

    if not os.path.isfile(args.inc_file):
        sys.exit(f"ERROR: file not found: {args.inc_file}")

    out_png = args.out or os.path.join(
        os.path.dirname(os.path.abspath(args.inc_file)), "rel_perm_curves.png"
    )

    print(f"Parsing: {args.inc_file}")
    rpts = parse_inc(args.inc_file)
    if not rpts:
        sys.exit("ERROR: no RPT blocks found in the file.")
    print(f"Found {len(rpts)} RPT block(s)")
    plot_relperm(rpts, out_png)


if __name__ == "__main__":
    main()
