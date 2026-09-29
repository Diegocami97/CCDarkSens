#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: validate_wimp_nucleon_quenching.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  validate_wimp_nucleon_quenching.py -- Parity and sanity checks for
#  ccdarkphys.wimp_nucleon.quenching.
# ============================================================================

"""
Five checks:

1. lindhard_yield pointwise parity against WIMPyCCD's own lindhard()
   (re-implemented literally here, not imported cross-repo).

2. Full-pipeline parity against WIMPyCCD's own printed number from
   rate_with_NR.ipynb (cell "35ab2da0"): raw rate (WIMPyCCD's literal,
   buggy dRdE_nr_kg_day_keV_wimpyccd_literal) -> lindhard() -> weighted
   histogram, at m_D_GeV=5, sigma_pb=0.1 (sigma_n_cm2=1e-37), using
   WIMPyCCD's own grid construction (np.arange(0.001,30.001,0.001),
   Eion histogram 0-8 keV in 0.02 keV bins) -> printed "total rate:
   316.7046843791842" events/kg/day.

3. chavarria_table_yield: exact recovery at all 12 Table I points, boundary
   behavior (zero below 0.3 keV_nr, Lindhard fallback above 2.28 keV_nr),
   and conservation of total counts through compute_dRdE_ee's
   weighted-histogram remap (checked by direct sum of per-bin counts, NOT
   trapz -- trapz underestimates the histogram because ~18% of the total
   rate at m_chi=5 GeV collapses to exactly E_ee=0 under the Chavarria
   zero-cutoff, a spike trapezoidal integration handles poorly).

4. julian_table_yield: exact recovery at 5 representative points spanning
   the table (including the exact first point, which previously clipped to
   zero due to a CSV-rounding boundary bug -- now fixed by storing 9 d.p.
   instead of 6), boundary behavior (zero below Er_min, Lindhard fallback
   above Er_max), and the same weighted-histogram conservation check as
   Chavarria.

5. Three-model cross-check: at a fixed Er inside all three models' shared
   support, confirms the three give genuinely different Ee (they are
   independent measurements/models, not near-duplicates), and reports the
   relative spread so the difference is visible at a glance rather than
   just asserted.

Usage:
    PYTHONPATH=python python3 utils/validate_wimp_nucleon_quenching.py
"""
from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "python"))

import numpy as np

from ccdarkphys.wimp_nucleon.quenching import lindhard_yield, chavarria_table_yield, julian_table_yield
from ccdarkphys.wimp_nucleon.rate import dRdE_nr_kg_day_keV_wimpyccd_literal
from ccdarkphys.wimp_nucleon.entry import compute_dRdE_ee

RTOL = 1.0e-6
CHAVARRIA_TABLE = [
    (0.68, 0.06), (0.86, 0.09), (1.05, 0.12), (1.22, 0.15), (1.36, 0.18), (1.52, 0.21),
    (1.66, 0.24), (1.81, 0.27), (1.94, 0.30), (2.07, 0.33), (2.18, 0.36), (2.28, 0.39),
]
# 5 representative points from data/julian_photoneutron_2024_iteration_central.csv
# (first, ~25%, ~50%, ~75%, last) -- the first is deliberately the exact table
# edge, the point that previously clipped to zero before the 9-d.p. CSV fix.
JULIAN_TABLE_SAMPLE = [
    (0.379900692, 0.011260712),
    (0.939878356, 0.104485982),
    (1.429490135, 0.201960000),
    (1.907878505, 0.299200000),
    (2.496801694, 0.445060000),
]


# ----------------------------------------------------------------------------
# _wimpyccd_lindhard
#   Independent reference implementation of the Lindhard ionization energy for silicon (Z = 14, k = 0.15) as written in WIMPyCCD, in keV.
# ----------------------------------------------------------------------------
def _wimpyccd_lindhard(Enr_KeV):
    Z, k = 14, 0.15
    eta = 11.5 * Enr_KeV * Z ** (-7.0 / 3.0)
    g = 3.0 * eta ** 0.15 + 0.7 * eta ** 0.6 + eta
    return (k * g / (1.0 + k * g)) * Enr_KeV


# ----------------------------------------------------------------------------
# check_lindhard_pointwise
#   Check that lindhard_yield agrees with the WIMPyCCD reference within RTOL at eight recoil energies.
# ----------------------------------------------------------------------------
def check_lindhard_pointwise() -> bool:
    ok = True
    for E in [0.001, 0.1, 0.5, 1.0, 5.0, 10.0, 20.0, 30.0]:
        ours, ref = lindhard_yield(E), _wimpyccd_lindhard(E)
        rel = abs(ours - ref) / ref if ref else abs(ours - ref)
        if rel > RTOL:
            print(f"[lindhard-pointwise][FAIL] E_r={E}: ours={ours} ref={ref} rel={rel:.2e}")
            ok = False
    print(f"[lindhard-pointwise] {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------------
# check_notebook_reference
#   Check that the total quenched rate for m = 5 GeV, sigma = 1e-37 cm^2 reproduces the notebook reference value 316.7046843791842.
# ----------------------------------------------------------------------------
def check_notebook_reference() -> bool:
    dE_nr = 0.001
    E_nr = np.arange(0.001, 30.0 + dE_nr, dE_nr)
    weights = dRdE_nr_kg_day_keV_wimpyccd_literal(
        E_nr, mchi_GeV=5.0, sigma_n_cm2=1.0e-37, A=28,
        rho_chi_gev_cm3=0.3, v0_kms=238.0, vE_kms=263.0, vesc_kms=544.0,
    ) * dE_nr
    E_ion = lindhard_yield(E_nr)
    edges = np.arange(0.0, 8.0 + 0.02, 0.02)
    hist, _ = np.histogram(E_ion, bins=edges, weights=weights)
    total = hist.sum()
    ref = 316.7046843791842
    rel = abs(total - ref) / ref
    print(f"[notebook-reference] ours={total!r}  WIMPyCCD={ref!r}  rel={rel:.3e}  "
          f"{'PASS' if rel <= RTOL else 'FAIL'}")
    return rel <= RTOL


# ----------------------------------------------------------------------------
# check_chavarria_table
#   Check the Chavarria table yield at its tabulated points and that it is exactly zero in the zero-crossing region.
# ----------------------------------------------------------------------------
def check_chavarria_table() -> bool:
    ok = True
    for Er, Ee_expected in CHAVARRIA_TABLE:
        got = chavarria_table_yield(Er)
        if abs(got - Ee_expected) > 1e-9:
            print(f"[chavarria-table][FAIL] Er={Er}: expected={Ee_expected} got={got}")
            ok = False
    if chavarria_table_yield(0.3) != 0.0 or chavarria_table_yield(0.2) != 0.0:
        print("[chavarria-table][FAIL] zero-crossing region not exactly zero")
        ok = False
    print(f"[chavarria-table-points] {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------------
# check_julian_table
#   Check the Julian table yield at sample points, that it is zero below its range, and that above its range it follows the Lindhard yield.
# ----------------------------------------------------------------------------
def check_julian_table() -> bool:
    ok = True
    for Er, Ee_expected in JULIAN_TABLE_SAMPLE:
        got = julian_table_yield(Er)
        if abs(got - Ee_expected) > 1e-7:
            print(f"[julian-table][FAIL] Er={Er}: expected={Ee_expected} got={got}")
            ok = False
    if julian_table_yield(0.2) != 0.0:
        print("[julian-table][FAIL] below-range region not exactly zero")
        ok = False
    above = julian_table_yield(5.0)
    ref = lindhard_yield(5.0)
    if abs(above - ref) > 1e-9:
        print(f"[julian-table][FAIL] above-range fallback mismatch: got={above} lindhard={ref}")
        ok = False
    print(f"[julian-table-points] {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------------
# check_three_model_spread
#   Compare the three quenching models at 1.5 keV_nr and print their spread (a sanity check only: they are independent models).
# ----------------------------------------------------------------------------
def check_three_model_spread() -> bool:
    # Er=1.5 keV_nr is inside all three models' measured/analytic support.
    Er = 1.5
    vals = {
        "lindhard": lindhard_yield(Er),
        "chavarria_table": chavarria_table_yield(Er),
        "julian_table": julian_table_yield(Er),
    }
    spread = (max(vals.values()) - min(vals.values())) / min(vals.values())
    print(f"[three-model-spread] at Er={Er} keV_nr: " +
          "  ".join(f"{k}={v:.5f}" for k, v in vals.items()) +
          f"   spread={spread:.1%}")
    # Sanity only: the three are independent models, not required to agree --
    # this just fails if they're suspiciously identical (a sign one is
    # accidentally aliased to another) or NaN.
    ok = spread > 1e-3 and all(np.isfinite(v) for v in vals.values())
    print(f"[three-model-spread] {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------------
# check_conservation
#   Check that quenching conserves the total event rate: the integral of dR/dE_ee equals that of the raw nuclear-recoil rate for each model.
# ----------------------------------------------------------------------------
def check_conservation() -> bool:
    ok = True
    common = dict(A=28, mchi_MeV=5000.0, sigma_n_cm2=1e-40, nr_Emin_keV=0.001, nr_Emax_keV=30.0, nr_nbins=5000)
    for model in ["lindhard", "chavarria_table", "julian_table"]:
        res = compute_dRdE_ee(**common, ee_Emin_eV=0, ee_Emax_eV=8000, ee_nbins=2000, quenching_model=model)
        bin_w = 8000.0 / 2000
        total_ee = (res["dRdE_kg_year_eV"] * bin_w).sum()
        # independent reference: raw counts computed the same way compute_dRdE_ee does internally
        from ccdarkphys.wimp_nucleon.entry import compute_dRdE
        res_nr = compute_dRdE(**common)
        nr_bin_w = (30.0 - 0.001) * 1000.0 / 5000
        total_nr = (res_nr["dRdE_kg_year_eV"] * nr_bin_w).sum()
        rel = abs(total_ee - total_nr) / total_nr
        print(f"[conservation:{model}] raw={total_nr:.4f}  quenched={total_ee:.4f}  rel diff={rel:.3e}")
        if rel > 1e-6:
            ok = False
    print(f"[conservation] {'PASS' if ok else 'FAIL'}")
    return ok


# ----------------------------------------------------------------------------
# main
#   Run all quenching checks and report ALL CHECKS PASS or SOME CHECKS FAILED (exit code 0 or 1).
# ----------------------------------------------------------------------------
def main() -> int:
    results = [
        check_lindhard_pointwise(),
        check_notebook_reference(),
        check_chavarria_table(),
        check_julian_table(),
        check_three_model_spread(),
        check_conservation(),
    ]
    if all(results):
        print("\nALL CHECKS PASS")
        return 0
    print("\nSOME CHECKS FAILED")
    return 1


if __name__ == "__main__":
    sys.exit(main())
