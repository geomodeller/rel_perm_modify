import numpy as np

def parse_relperm(filepath, verbose=False):
    ## TODO: slt should be renamed to sgt for clarity
    swt_rows = []
    sgt_rows = []
    current_table = None

    with open(filepath, "r") as f:
        for line in f:
            stripped = line.strip()

            # Skip empty lines and comments
            if not stripped or stripped.startswith("**"):
                continue

            if stripped.startswith("RPT"):
                continue

            if stripped == "SWT":
                current_table = "SWT"
                continue

            if stripped == "SGT":
                current_table = "SGT"
                continue

            # Try to parse a data row
            try:
                values = [float(v) for v in stripped.split()]
                if len(values) == 4:
                    if current_table == "SWT":
                        if verbose: print("SWT row:", values)  # Debugging line
                        swt_rows.append(values)
                    elif current_table == "SGT":
                        if verbose: print("SGT row:", values)  # Debugging line
                        sgt_rows.append(values)
            except ValueError:
                continue

    swt = np.array(swt_rows)  # columns: Sw, krw, krow, Pc
    sgt = np.array(sgt_rows)  # columns: Sg, krg, krog, Pc

    assert not np.array_equal(swt, sgt), "SWT and SGT data should not be identical"
    return swt, sgt


if __name__ == "__main__":
    filepath = "cmg_datafile/rock_fluid_properties.inc"
    swt, sgt = parse_relperm(filepath)

    print("=== SWT (water-oil) ===")
    print(f"Shape: {swt.shape}")
    print(f"{'Sw':>10} {'krw':>14} {'krow':>14} {'Pc':>14}")
    for row in swt:
        print(f"{row[0]:>10.4f} {row[1]:>14.8f} {row[2]:>14.8f} {row[3]:>14.6f}")

    print("\n=== SGT (gas-oil) ===")
    print(f"Shape: {sgt.shape}")
    print(f"{'Sg':>10} {'krg':>14} {'krog':>14} {'Pc':>14}")
    for row in sgt:
        print(f"{row[0]:>10.4f} {row[1]:>14.8f} {row[2]:>14.8f} {row[3]:>14.6f}")

    # Export as .npy files
    np.save("swt.npy", swt)
    np.save("sgt.npy", sgt)
    print("\nSaved: swt.npy, sgt.npy")

    # Export as CSV files
    header_swt = "Sw,krw,krow,Pc"
    header_sgt = "Sg,krg,krog,Pc"
    np.savetxt("swt.csv", swt, delimiter=",", header=header_swt, comments="")
    np.savetxt("sgt.csv", sgt, delimiter=",", header=header_sgt, comments="")
    print("Saved: swt.csv, sgt.csv")
