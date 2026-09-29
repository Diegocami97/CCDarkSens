#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: verify_dRdE_units.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  verify_dRdE_units.py -- Self-consistency check that ccdarkphys QEDark
#  outputs are in events/(kg·year·eV) and integrates correctly over energy
#  bins.
# ============================================================================

"""
Verify whether our dRdE output should be (notebook_return / dE) or (notebook_return as-is).

The C++ pipeline does:  counts = rate * exposure_kg_year * dE
So the CSV must contain rate in events/(kg·year·eV). We need to know:
  (A) Notebook dRdE at one Ee returns "integrated rate in bin dE" [events/(kg·year)] → we divide by dE → CSV correct.
  (B) Notebook dRdE at one Ee returns "differential dR/dE" [events/(kg·year·eV)] → we must NOT divide by dE.

Check: For a fixed energy bin [E_lo, E_hi], the INTEGRATED rate should be
  integral_{E_lo}^{E_hi} (dR/dE) dE  [events/(kg·year)].

We compute with our code (with division by dE) and form:
  integrated_A = sum over E in [E_lo, E_hi] of (dRdE_kg_year_eV * dE)  = sum(dRdE_kg_year)  [if we use divided]
  So integrated_A = sum of our "per bin" values = what notebook would sum in dRdnearray.

If the notebook's dRdE returns "per dE bin" (events/(kg·year)), then for one big bin
  notebook_bin_rate = sum of dRdE at each 0.1 eV point in that bin (no extra dE factor).
So we can't call the notebook from here; we can only check self-consistency:

  (1) Our dRdE_kg_year_eV = dRdE_kg_year / dE.
  (2) Integrated rate in [E1, E2] = sum(dRdE_kg_year_eV[i] * dE) for E_eV[i] in [E1,E2] = sum(dRdE_kg_year) in that range.
  (3) So the same physics either way; the only question is what we WRITE to CSV. The C++ expects events/(kg·year·eV), so we must write dRdE_kg_year_eV = dRdE_kg_year / dE.

Run this script and optionally compare the printed "integrated rate in first 3.8 eV bin" with the first non-zero value from the notebook's dRdnearray for the same (mX, sigma_e, Ebin=3.8). If they match, the convention is correct.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

import numpy as np
from ccdarkphys.qedark.entry import compute_dRdE

# ----------------------------------------------------------------------------
# main
#   Check the units of the QEDark dR/dE for one point (0.5 MeV, sigma_e = 1e-36 cm^2): integrate dR/dE over the first 3.8 eV bin (1.2-5.0 eV) in two ways and print the values so they can be compared with the QEDark notebook.
# ----------------------------------------------------------------------------
def main():
    mchi_MeV = 0.5
    sigma_e_cm2 = 1e-36
    halo = {"v0_kms": 220.0, "vE_kms": 232.0, "vesc_kms": 544.0}
    dE = 0.1

    res = compute_dRdE(
        material="Si",
        mediator="heavy",
        mchi_eV=mchi_MeV * 1e6,
        sigma_e_cm2=sigma_e_cm2,
        halo=halo,
        binsize_eV=dE,
    )
    E_eV = res["E_eV"]
    dRdE_kg_year_eV = res["dRdE_kg_year_eV"]

    # Reconstruct "per bin" rate (what notebook returns at each point)
    dRdE_kg_year_per_bin = dRdE_kg_year_eV * dE

    # First bin of width 3.8 eV (like notebook Si epsilon): E from 1.2 to 1.2+3.8 = 5.0 eV
    E_lo, E_hi = 1.2, 5.0
    mask = (E_eV >= E_lo) & (E_eV < E_hi)
    n_bins = int(np.round((E_hi - E_lo) / dE))

    # Integrated rate in [E_lo, E_hi] using our dR/dE (events/(kg·year·eV))
    integrated_from_differential = np.sum(dRdE_kg_year_eV[mask] * dE)

    # Same using "per bin" interpretation (sum of notebook-like returns)
    integrated_from_per_bin = np.sum(dRdE_kg_year_per_bin[mask])

    print("Verify dRdE units (mchi=0.5 MeV, sigma_e=1e-36 cm^2)")
    print("  dE =", dE, "eV")
    print("  Bin [1.2, 5.0] eV (first 3.8 eV bin):")
    print("    Integrated rate (sum of dRdE_kg_year_eV * dE) =", integrated_from_differential, "events/(kg·year)")
    print("    Sum of (dRdE_kg_year_eV * dE) same as sum of per-bin values:", np.allclose(integrated_from_differential, integrated_from_per_bin))
    print("  First few E_eV and dRdE_kg_year_eV (should be events/(kg·year·eV)):")
    for i in range(min(5, len(E_eV))):
        if dRdE_kg_year_eV[i] > 0:
            print(f"    E={E_eV[i]:.2f} eV  dRdE_kg_year_eV={dRdE_kg_year_eV[i]:.6e}")
    print()
    print("To confirm: run the QEdark notebook dRdnearray(Si, mX=0.5e6, Ebin=3.8, ...) with sigma_e=1")
    print("  and same halo. First non-zero bin value should be ~ integrated rate in that bin [events/(kg·year)].")
    print("  Our integrated_rate above should match that first bin value (within sigma_e and halo).")

if __name__ == "__main__":
    main()
