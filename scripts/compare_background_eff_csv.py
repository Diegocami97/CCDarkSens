#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: compare_background_eff_csv.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: compare_background_eff_csv.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_background_eff_csv.py -- Compare
#  outputs/scan_pattern/Background_efficiencies.csv to
#  data/Background_efficiencies.csv.
# ============================================================================

"""Compare outputs/scan_pattern/Background_efficiencies.csv to data/Background_efficiencies.csv."""
import pandas as pd
import sys

# ----------------------------------------------------------------------------
# load_csv
#   Read a Background_efficiencies.csv into a DataFrame indexed by the identified pattern (iden_pat).
# ----------------------------------------------------------------------------
def load_csv(path):
    df = pd.read_csv(path)
    df = df.set_index("iden_pat")
    # iden_pat may be int or str; normalize to string for column matching
    df.index = df.index.astype(str)
    return df

# ----------------------------------------------------------------------------
# main
#   Compare the background-efficiency CSV produced by the app with the collaboration reference: structure (rows and columns), the diagonal, the 15 largest cell differences and the row sums.
# ----------------------------------------------------------------------------
def main():
    ref_path = "data/Background_efficiencies.csv"
    out_path = "outputs/scan_pattern/Background_efficiencies.csv"
    if len(sys.argv) >= 3:
        ref_path, out_path = sys.argv[1], sys.argv[2]
    ref = load_csv(ref_path)
    out = load_csv(out_path)

    print("=== Structure ===")
    print(f"Reference: {len(ref)} rows, columns: {list(ref.columns)[:5]}...")
    print(f"Output:    {len(out)} rows, columns: {list(out.columns)[:5]}...")
    ref_rows = set(ref.index)
    out_rows = set(out.index)
    print(f"Reference row ids (sorted): {sorted(ref_rows, key=lambda x: (len(x), x))}")
    print(f"Output row ids (sorted):   {sorted(out_rows, key=lambda x: (len(x), x))}")
    if ref_rows != out_rows:
        print(f"Only in reference: {ref_rows - out_rows}")
        print(f"Only in output:   {out_rows - ref_rows}")
    assert list(ref.columns) == list(out.columns), "Column order differs"

    print("\n=== Diagonal (correct identification) ===")
    for pid in sorted(ref.index, key=lambda x: (len(x), x)):
        col = f"eff_{pid}"
        if col not in ref.columns or pid not in ref.index:
            continue
        r = ref.loc[pid, col] if pid in ref.index else 0.0
        o = out.loc[pid, col] if pid in out.index else 0.0
        diff = o - r
        print(f"  iden_pat={pid:>3}  ref={r:.6f}  out={o:.6f}  diff={diff:+.6f}")

    print("\n=== Largest absolute differences (by cell) ===")
    common_idx = ref.index.intersection(out.index)
    diffs = []
    for i in common_idx:
        for c in ref.columns:
            rv = ref.loc[i, c]
            ov = out.loc[i, c] if i in out.index else 0.0
            diffs.append((abs(ov - rv), i, c, rv, ov))
    for abs_diff, i, c, rv, ov in sorted(diffs, reverse=True)[:15]:
        print(f"  {i} / {c}:  ref={rv:.4e}  out={ov:.4e}  |diff|={abs_diff:.4e}")

    print("\n=== Row sum (should be ~1.0 per row) ===")
    for pid in sorted(ref.index, key=lambda x: (len(x), x))[:8]:
        rs = ref.loc[pid].sum()
        os = out.loc[pid].sum() if pid in out.index else 0.0
        print(f"  iden_pat={pid:>3}  ref sum={rs:.6f}  out sum={os:.6f}")

if __name__ == "__main__":
    main()
