#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: validate_wimp_nucleon_rate_parity.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  validate_wimp_nucleon_rate_parity.py -- Bit-for-bit parity check:
#  ccdarkphys.wimp_nucleon.rate vs WIMPyCCD's own printed reference number.
# ============================================================================

"""
Reproduces, using the ported ccdarkphys.wimp_nucleon.rate module, the exact
number WIMPyCCD's own rate_with_NR.ipynb notebook printed for its own
default point:

    Notebook cell (id "30260b12"):
        params = DMParams(m_D_GeV=5, sigma_pb=0.1, exposure_kgdays=11,
                           E_th_KeV=0.63)
        integrate.quad(lambda E: get_rate_ben(E, dm_params=params), 0, 10)[0]
        -> 316.80360741504927   (events / kg / day, integrated 0-10 keV_nr)

DMParams defaults not overridden in that call (from WIMPyCCD's
analysis/dm_params.py): rho_GeV_cm3=0.3, v0_kms=238, vE_kms=263,
vesc_kms=544, mass_number_A=28. sigma_pb -> sigma_n_cm2 via
DMParams.sigma_m2's own convention (sigma_pb * 1e-40 m^2 = sigma_pb * 1e-36
cm^2), so sigma_pb=0.1 -> sigma_n_cm2=1e-37.

This is a pure numerical-parity check of the PORT (does our code reproduce
WIMPyCCD's own printed number), not a validation that the underlying
physics convention is correct -- see ccdarkphys/wimp_nucleon/rate.py's
module docstring for the two caveats found during the post-parity sanity
pass, BOTH now confirmed as bugs and fixed in the production
dRdE_nr_kg_day_keV: item 1 (form factor identically 1.0, missing a
momentum->wavenumber conversion) and item 2 (DM-nucleus reduced mass used
in the rate prefactor where the standard per-nucleon-sigma convention
needs DM-nucleon -- a mass-dependent effect reaching ~99% at 10 GeV for
Si). This test deliberately calls dRdE_nr_kg_day_keV_wimpyccd_literal,
which reproduces the ORIGINAL WIMPyCCD formula with both bugs intact, so
it keeps checking the historical port rather than silently start
comparing against a moving target -- and so this repo can reproduce
exactly what the original WIMPyCCD code would have computed, if ever
needed for comparison or provenance.

Usage:
    PYTHONPATH=python python3 utils/validate_wimp_nucleon_rate_parity.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "python"))

import scipy.integrate as integrate

from ccdarkphys.wimp_nucleon.rate import dRdE_nr_kg_day_keV_wimpyccd_literal

WIMPYCCD_REFERENCE_RATE_KG_DAY = 316.80360741504927
WIMPYCCD_REFERENCE_INPUTS = dict(
    mchi_GeV=5.0,
    sigma_n_cm2=1.0e-37,  # sigma_pb=0.1 -> *1e-36 cm^2/pb
    A=28,
    rho_chi_gev_cm3=0.3,
    v0_kms=238.0,
    vE_kms=263.0,   # WIMPyCCD's own default, NOT the Baxter/Migdal-consistent 253.7 used elsewhere
    vesc_kms=544.0,
)
INTEGRATION_BOUNDS_KEV = (0.0, 10.0)
RTOL = 1.0e-6


# ----------------------------------------------------------------------------
# main
#   Integrate the ported WIMP-nucleus rate over the reference bounds and compare it with the WIMPyCCD reference number (events/kg/day); PASS if the relative difference is within RTOL (exit code 0), otherwise FAIL (1).
# ----------------------------------------------------------------------------
def main() -> int:
    integrand = lambda E: dRdE_nr_kg_day_keV_wimpyccd_literal(E, **WIMPYCCD_REFERENCE_INPUTS)
    ported_rate, abserr = integrate.quad(integrand, *INTEGRATION_BOUNDS_KEV)

    rel_diff = abs(ported_rate - WIMPYCCD_REFERENCE_RATE_KG_DAY) / WIMPYCCD_REFERENCE_RATE_KG_DAY

    print(f"[parity] WIMPyCCD reference : {WIMPYCCD_REFERENCE_RATE_KG_DAY!r} events/kg/day")
    print(f"[parity] ported rate        : {ported_rate!r} events/kg/day (quad abserr={abserr:.3e})")
    print(f"[parity] relative difference: {rel_diff:.3e}  (tolerance {RTOL:.0e})")

    if rel_diff > RTOL:
        print("[parity] FAIL — ported rate does not match WIMPyCCD's reference number.")
        return 1

    print("[parity] PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
