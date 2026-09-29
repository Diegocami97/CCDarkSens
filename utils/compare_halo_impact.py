#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: compare_halo_impact.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  compare_halo_impact.py -- Quick script comparing integrated QEDark rates
#  for default vs updated SHM halo velocity parameters at several masses.
# ============================================================================

"""Compare rate impact of halo (220,232) vs (238,263) km/s. Quick estimate."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import numpy as np
from ccdarkphys.qedark.entry import compute_dRdE

# ----------------------------------------------------------------------------
# integrated_rate_first_bin
#   Rate summed over 1.2 <= E < 5 eV for a heavy-mediator silicon point at sigma_e = 1e-36 cm^2 with the given halo, in events per kg per year.
# ----------------------------------------------------------------------------
def integrated_rate_first_bin(mchi_MeV, halo, dE=0.1):
    res = compute_dRdE(
        material="Si", mediator="heavy",
        mchi_eV=mchi_MeV * 1e6, sigma_e_cm2=1e-36,
        halo=halo, binsize_eV=dE,
    )
    E_lo, E_hi = 1.2, 5.0
    mask = (res["E_eV"] >= E_lo) & (res["E_eV"] < E_hi)
    return np.sum(res["dRdE_kg_year_eV"][mask] * dE)

# ----------------------------------------------------------------------------
# main
#   Print the rate ratio between the halo (238, 263 km/s) and the old halo (220, 232 km/s) for several masses, and the total rate for 1 MeV.
# ----------------------------------------------------------------------------
def main():
    halo_old = {"v0_kms": 220.0, "vE_kms": 232.0, "vesc_kms": 544.0}
    halo_new = {"v0_kms": 238.0, "vE_kms": 263.0, "vesc_kms": 544.0}

    masses = [0.5, 1.0, 5.0, 10.0, 100.0]  # MeV
    print("Halo impact: (220,232) vs (238,263) km/s, first 3.8 eV bin, sigma_e=1e-36")
    print("m_chi [MeV]   R_old [1/(kg·y)]   R_new [1/(kg·y)]   R_new/R_old")
    print("-" * 65)
    for m in masses:
        r_old = integrated_rate_first_bin(m, halo_old)
        r_new = integrated_rate_first_bin(m, halo_new)
        ratio = r_new / r_old if r_old > 0 else float("nan")
        print(f"  {m:6.1f}       {r_old:.4e}      {r_new:.4e}      {ratio:.3f}")
    print()
    # One more: total rate over full range
    print("Total rate (sum over all E bins) for m_chi=1 MeV:")
    for name, halo in [("old 220,232", halo_old), ("new 238,263", halo_new)]:
        res = compute_dRdE("Si", "heavy", 1e6, 1e-36, halo)
        total = np.sum(res["dRdE_kg_year_eV"] * 0.1)
        print(f"  {name}: {total:.4e} events/(kg·year)")
    r_old = np.sum(compute_dRdE("Si", "heavy", 1e6, 1e-36, halo_old)["dRdE_kg_year_eV"] * 0.1)
    r_new = np.sum(compute_dRdE("Si", "heavy", 1e6, 1e-36, halo_new)["dRdE_kg_year_eV"] * 0.1)
    print(f"  Ratio (new/old): {r_new/r_old:.3f}")

if __name__ == "__main__":
    main()
