#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: verify_p100K_scaled.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  verify_p100K_scaled.py -- Sanity checks for scaled p100K tables vs
#  original reference (Si band-gap QA)
# ============================================================================
"""Sanity checks for scaled p100K tables vs original reference."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_p100K_scaled import (  # noqa: E402
    REF_EH_EV,
    REF_GAP_EV,
    build_scaled_table,
    extend_energy_grid,
    load_p100k_csv,
    map_E_prime,
)

REF = Path("data/p100K_table.csv")
MANIFEST = Path("configs/band_gap_pheno_scenarios.json")


# ----------------------------------------------------------------------------
# row_at
#   P(n | E = e0) for every n by linear interpolation (zero outside the grid).
# ----------------------------------------------------------------------------
def row_at(E: np.ndarray, P: np.ndarray, e0: float) -> np.ndarray:
    out = np.zeros(P.shape[0])
    for n in range(P.shape[0]):
        out[n] = np.interp(e0, E, P[n], left=0.0, right=0.0)
    return out


# ----------------------------------------------------------------------------
# main
#   Verify the p100K tables: an identity rebuild of the reference must reproduce it, and every scaled table of the manifest must have the right grid, be zero below its gap and have consistent probabilities above it; prints OK or FAIL for each check.
# ----------------------------------------------------------------------------
def main() -> None:
    E_raw, P_raw = load_p100k_csv(REF)
    print("=" * 70)
    print("p100K table verification")
    print("=" * 70)

    print(f"\n[original] {REF}")
    print(f"  E: [{E_raw[0]:.4g}, {E_raw[-1]:.4g}] eV  n_rows={E_raw.size}  n_ne={P_raw.shape[0]}")
    print(f"  step dE ≈ {E_raw[1]-E_raw[0]:.4g} eV")

    E_ext = extend_energy_grid(E_raw, E_min=0.05)
    P_id = build_scaled_table(E_ext, P_raw, E_raw, REF_GAP_EV, REF_EH_EV)
    # Identity map should reproduce ref at E >= 1.1
    mask = E_ext >= E_raw[0]
    max_diff = 0.0
    for i in np.where(mask)[0]:
        d = np.max(np.abs(P_id[:, i] - np.array([np.interp(E_ext[i], E_raw, P_raw[n]) for n in range(P_raw.shape[0])])))
        max_diff = max(max_diff, d)
    print(f"\n[identity rebuild] extended grid + gap={REF_GAP_EV} eh={REF_EH_EV}")
    print(f"  E: [{E_ext[0]:.4g}, {E_ext[-1]:.4g}] eV  n_rows={E_ext.size}")
    print(f"  max |P_rebuilt - P_ref_interp| for E>={E_raw[0]:.4g} eV: {max_diff:.2e}  "
          f"{'OK' if max_diff < 1e-10 else 'FAIL'}")

    with open(MANIFEST, encoding="utf-8") as f:
        man = json.load(f)
    entries = [man["reference"]] + man.get("scenarios", [])

    print("\n" + "-" * 70)
    print("Scaled tables from manifest")
    print("-" * 70)

    all_ok = True
    for ent in entries:
        sc = ent.get("scenario", "?")
        csv = Path(ent["ionization_csv"])
        gap = float(ent["band_gap_eV"])
        eh = float(ent["eh_pair_eV"])
        if not csv.exists():
            print(f"\n[{sc}] MISSING {csv}")
            all_ok = False
            continue

        E, P = load_p100k_csv(csv)
        ok_e = E[0] <= 0.05 + 1e-9 and E[-1] >= 50.0 - 1e-9
        ok_n = P.shape == P_raw.shape and E.size == E_ext.size

        print(f"\n[{sc}] {csv.name}")
        print(f"  label: {ent.get('label', '')}")
        print(f"  gap={gap} eV  eh={eh} eV  scale={REF_EH_EV/eh:.4g}")
        print(f"  E: [{E[0]:.4g}, {E[-1]:.4g}] eV  rows={E.size}  "
              f"grid_end={'OK' if ok_e else 'FAIL'}")

        # Below gap: all zero
        below = E < gap - 1e-9
        max_below = float(P[:, below].max()) if below.any() else 0.0
        print(f"  max P below gap (E<{gap} eV): {max_below:.2e}  "
              f"{'OK' if max_below < 1e-12 else 'FAIL'}")

        # Row sums ~ 1 where ionization on
        sums = P.sum(axis=0)
        above = E >= gap
        s_above = sums[above]
        bad_sum = np.sum(np.abs(s_above - 1.0) > 0.02)
        print(f"  rows with |sum P - 1| > 0.02 (E>=gap): {bad_sum}/{above.sum()}  "
              f"{'OK' if bad_sum == 0 else 'warn'}")

        # Anchor at gap
        if np.any(np.isclose(E, gap)):
            p1_gap = float(P[0, np.isclose(E, gap)][0])
        else:
            p1_gap = float(row_at(E, P, gap)[0])
        p1_ref = float(row_at(E_raw, P_raw, REF_GAP_EV)[0])
        ok_anchor = abs(p1_gap - p1_ref) < 0.02
        print(f"  anchor P(n=1) @ E={gap}: new={p1_gap:.4f}  ref@1.2={p1_ref:.4f}  "
              f"{'OK' if ok_anchor else 'CHECK'}")

        # High-E tail: compare to ref at same E' for a few points
        test_E = [10.0, 20.0, 50.0]
        print(f"  high-E check (E, E', max|P_new-P_ref(E')|):")
        for e0 in test_E:
            if e0 < gap:
                continue
            ep = map_E_prime(e0, gap, REF_GAP_EV, eh, REF_EH_EV)
            pn = row_at(E, P, e0)
            pr = row_at(E_raw, P_raw, ep)
            d = float(np.max(np.abs(pn - pr)))
            print(f"    E={e0:5.1f}  E'={ep:6.2f}  max_diff={d:.2e}")

        if sc == "ref":
            # ref table on disk should match raw (no extension in file)
            if E.size != E_raw.size:
                print(f"  note: ref CSV on disk has {E.size} rows (extended tables have {E_ext.size})")

        if not (ok_e and max_below < 1e-12):
            all_ok = False

    print("\n" + "=" * 70)
    print("Overall:", "PASS" if all_ok else "see FAIL/CHECK above")
    print("=" * 70)


if __name__ == "__main__":
    main()
