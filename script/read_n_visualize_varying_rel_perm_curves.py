import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

INC_FILE = "../cmg_datafile/rock_fluid_properties_varying.inc"


def parse_inc(filepath):
    """Parse RPT blocks from CMG .inc file, returning dict keyed by RPT number."""
    data = {}
    current_rpt = None
    current_table = None

    with open(filepath) as f:
        lines = f.readlines()

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("**"):
            continue

        rpt_match = re.match(r"^RPT\s+(\d+)\s*$", stripped)
        if rpt_match:
            current_rpt = int(rpt_match.group(1))
            data[current_rpt] = {"swt": [], "sgt": []}
            current_table = None
            continue

        if stripped == "SWT":
            current_table = "swt"
            continue
        if stripped == "SGT":
            current_table = "sgt"
            continue

        if current_rpt is not None and current_table is not None:
            values = stripped.split()
            try:
                row = [float(v) for v in values]
                if len(row) >= 3:
                    data[current_rpt][current_table].append(row)
            except ValueError:
                pass

    # Convert to numpy arrays; deduplicate trailing repeated rows
    for rpt in data:
        for tbl in ("swt", "sgt"):
            arr = np.array(data[rpt][tbl])
            # Remove duplicate last row (CMG pads with a repeated endpoint)
            if len(arr) > 1 and np.allclose(arr[-1], arr[-2]):
                arr = arr[:-1]
            data[rpt][tbl] = arr

    return data


def main():
    data = parse_inc(INC_FILE)
    n_rpt = len(data)

    colors = cm.tab10(np.linspace(0, 1, n_rpt))

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    ax_krw, ax_krow = axes[0]
    ax_krg, ax_krog = axes[1]

    for idx, rpt in enumerate(sorted(data)):
        swt = data[rpt]["swt"]  # cols: Sw, krw, krow, Pc
        sgt = data[rpt]["sgt"]  # cols: Sg, krg, krog, Pc
        label = f"RPT {rpt}"
        c = colors[idx]

        ax_krw.plot(swt[:, 0], swt[:, 1], color=c, label=label)
        ax_krow.plot(swt[:, 0], swt[:, 2], color=c, label=label)
        ax_krg.plot(sgt[:, 0], sgt[:, 1], color=c, label=label)
        ax_krog.plot(sgt[:, 0], sgt[:, 2], color=c, label=label)

    ax_krw.set_xlabel("Water Saturation Sw")
    ax_krw.set_ylabel("krw")
    ax_krw.set_title("Water Relative Permeability (krw)")
    ax_krw.legend(fontsize=7, ncol=2)
    ax_krw.grid(True, alpha=0.3)

    ax_krow.set_xlabel("Water Saturation Sw")
    ax_krow.set_ylabel("krow")
    ax_krow.set_title("Oil Relative Permeability vs Water Sat. (krow)")
    ax_krow.legend(fontsize=7, ncol=2)
    ax_krow.grid(True, alpha=0.3)

    ax_krg.set_xlabel("Gas Saturation Sg")
    ax_krg.set_ylabel("krg")
    ax_krg.set_title("Gas Relative Permeability (krg)")
    ax_krg.legend(fontsize=7, ncol=2)
    ax_krg.grid(True, alpha=0.3)

    ax_krog.set_xlabel("Gas Saturation Sg")
    ax_krog.set_ylabel("krog")
    ax_krog.set_title("Oil Relative Permeability vs Gas Sat. (krog)")
    ax_krog.legend(fontsize=7, ncol=2)
    ax_krog.grid(True, alpha=0.3)

    fig.suptitle("Varying Relative Permeability Curves (RPT 1–10)", fontsize=14)
    plt.tight_layout()
    plt.savefig("varying_rel_perm_curves.png", dpi=150, bbox_inches="tight")
    print("Saved: varying_rel_perm_curves.png")
    plt.show()


if __name__ == "__main__":
    main()
