# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: quenching.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  quenching.py -- Nuclear recoil ionization efficiency (quenching / yield)
#  models: E_ee = Gamma(E_nr).
# ============================================================================

"""
Nuclear recoil ionization efficiency (quenching / yield) models: E_ee = Gamma(E_nr).

Four independent models, all taking/returning keV:

1. lindhard_yield -- analytic Lindhard theory (Lindhard, Nielsen, Scharff &
   Thomsen 1963), the standard textbook model. Ported from WIMPyCCD's
   analysis/yield_functions.py::lindhard() (k=0.15, Z=14 hardcoded for Si --
   matches our target).

2. chavarria_table_yield -- the actual MEASURED Si nuclear-recoil ionization
   efficiency, Chavarria et al. 2016 (Phys. Rev. D 94, 082007), Table I:
   12 points spanning Er = 0.68-2.28 keV_nr, measured with a 124Sb-9Be
   photoneutron source down to 60 eV_ee. The paper's headline result is that
   the measured efficiency deviates SIGNIFICANTLY from the Lindhard
   extrapolation in this range -- this is real calibration data, not a
   refinement of Lindhard's functional form.

   Extrapolation outside the measured range follows the same convention
   WIMPyCCD's own (unfinished) damicm_no_fano/damicm_fano functions were
   already structured for, and matches PhysRevD.94.082006's own stated
   threshold:
     - below 0.68 keV_nr: linear to zero at Er = 0.3 keV_nr (082006's own
       stated value, presumably derived from this exact measurement)
     - above 2.28 keV_nr: falls back to lindhard_yield(). No measurement
       exists there. NOTE: this produces a real discontinuity at the
       Er=2.28 boundary, since the whole point of this measurement is that
       Lindhard does NOT match the data in the measured range. Accepted as
       a known simplification (the same one WIMPyCCD's own code makes),
       not smoothed over -- flagged here rather than hidden.

3. chavarria_izraelevitch_table_yield -- Chavarria et al. 2016 below 2.28
   keV_nr, PLUS Izraelevitch et al. 2017 (JINST 12, P06014, arXiv:1702.00873,
   "antonella1" in PhysRevD.94.082006's own citation) from 2.28 up to 20.67
   keV_nr -- real measured data covering that whole range, matching what
   PhysRevD.94.082006 itself actually does (Sec. IV), rather than
   chavarria_table_yield's fallback to unconstrained Lindhard theory above
   2.28 keV_nr. See data/izraelevitch2017_table1.csv for provenance.
   Lindhard fallback only above 20.67 keV_nr (no measurement from either
   dataset there).

4. julian_table_yield -- a separate, more recent, UNPUBLISHED photo-neutron
   calibration measurement (J. Cuevas-Zepeda, "PhotoNeutronAnalysis-2024",
   provided directly by the author, not the Chavarria paper reprocessed).
   117 points, Er = 0.38-2.50 keV_nr -- lower threshold and finer sampling
   than Chavarria, and reports systematically LOWER ionization yield than
   Chavarria across most of the shared range. See
   data/julian_photoneutron_2024_iteration_central.csv for full provenance
   (iteration method, gamma-normalized central value; the author's own
   analysis also has a +-15% systematic band and an independent "integral
   method" cross-check, neither reproduced here).

   Because this dataset's own lower boundary (0.38 keV_nr) is already well
   below Chavarria's, and no independently-stated physical zero-crossing
   energy is available for it (unlike Chavarria/082006's stated 0.3 keV_nr),
   this function does NOT invent a taper-to-zero below the measured range --
   it clamps to 0 instead. Above the range it falls back to lindhard_yield,
   same convention as chavarria_table_yield.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d

from ccdarkphys.common import io as CIO

_LINDHARD_Z_DEFAULT = 14  # Si
_LINDHARD_K_DEFAULT = 0.15

_CHAVARRIA_ER_MIN_KEV = 0.68
_CHAVARRIA_ER_MAX_KEV = 2.28
_CHAVARRIA_ZERO_CROSSING_KEV = 0.3  # PhysRevD.94.082006's stated threshold


def lindhard_yield(E_r_keV, Z: int = _LINDHARD_Z_DEFAULT, k: float = _LINDHARD_K_DEFAULT):
    """Standard Lindhard theory quenching, exact port of WIMPyCCD's lindhard()."""
    E_r_keV = np.asarray(E_r_keV, dtype=float)
    eta = 11.5 * E_r_keV * Z ** (-7.0 / 3.0)
    g = 3.0 * eta ** 0.15 + 0.7 * eta ** 0.6 + eta
    nr_yield = k * g / (1.0 + k * g)
    return nr_yield * E_r_keV


_chavarria_interp = None  # lazily built on first use


# ----------------------------------------------------------------------------
# _load_chavarria_table
#   Linear interpolator E_ee(E_R) from the tabulated Chavarria 2016 quenching data (NaN outside the table).
# ----------------------------------------------------------------------------
def _load_chavarria_table():
    path = CIO.data_path(__file__, "data", "chavarria2016_table1.csv")
    data = np.loadtxt(path, delimiter=",", comments="#")
    Er, Ee = data[:, 0], data[:, 1]
    return interp1d(Er, Ee, kind="linear", bounds_error=False, fill_value=np.nan)


def chavarria_table_yield(E_r_keV):
    """
    Measured Si nuclear recoil ionization efficiency, Chavarria et al. 2016
    (Phys. Rev. D 94, 082007), Table I. See module docstring for the
    extrapolation convention outside the measured range (0.68-2.28 keV_nr).
    """
    global _chavarria_interp
    if _chavarria_interp is None:
        _chavarria_interp = _load_chavarria_table()

    scalar_input = np.ndim(E_r_keV) == 0
    E_r = np.atleast_1d(np.asarray(E_r_keV, dtype=float))
    E_ee = np.empty_like(E_r)

    below = E_r < _CHAVARRIA_ER_MIN_KEV
    above = E_r > _CHAVARRIA_ER_MAX_KEV
    within = ~below & ~above

    E_ee[within] = _chavarria_interp(E_r[within])

    # Below measured range: linear to zero at the 082006-stated threshold.
    Ee_at_Ermin = float(_chavarria_interp(_CHAVARRIA_ER_MIN_KEV))
    slope = Ee_at_Ermin / (_CHAVARRIA_ER_MIN_KEV - _CHAVARRIA_ZERO_CROSSING_KEV)
    below_taper = below & (E_r >= _CHAVARRIA_ZERO_CROSSING_KEV)
    below_zero = below & (E_r < _CHAVARRIA_ZERO_CROSSING_KEV)
    E_ee[below_taper] = slope * (E_r[below_taper] - _CHAVARRIA_ZERO_CROSSING_KEV)
    E_ee[below_zero] = 0.0

    # Above measured range: Lindhard fallback (known discontinuity -- see module docstring).
    E_ee[above] = lindhard_yield(E_r[above])

    return float(E_ee[0]) if scalar_input else E_ee


_IZRAELEVITCH_ER_MIN_KEV = 1.79
_IZRAELEVITCH_ER_MAX_KEV = 20.67

_izraelevitch_interp = None  # lazily built on first use


# ----------------------------------------------------------------------------
# _load_izraelevitch_table
#   Linear interpolator E_ee(E_R) from the tabulated Izraelevitch 2017 quenching data (NaN outside the table).
# ----------------------------------------------------------------------------
def _load_izraelevitch_table():
    path = CIO.data_path(__file__, "data", "izraelevitch2017_table1.csv")
    data = np.loadtxt(path, delimiter=",", comments="#")
    Er, Ee = data[:, 0], data[:, 1]
    return interp1d(Er, Ee, kind="linear", bounds_error=False, fill_value=np.nan)


def chavarria_izraelevitch_table_yield(E_r_keV):
    """
    Combined measured Si nuclear-recoil ionization efficiency: Chavarria et
    al. 2016 below 2.28 keV_nr, Izraelevitch et al. 2017 (JINST 12, P06014,
    arXiv:1702.00873) from 2.28 up to 20.67 keV_nr -- matching
    PhysRevD.94.082006's own stated approach (Sec. IV: "We adopt new results
    [antonella1, Chavarria:2016xsi] ... covering most of the energy range
    relevant for low-mass WIMP searches"), where "antonella1" is the
    Izraelevitch measurement. chavarria_table_yield alone instead falls back
    to unconstrained Lindhard theory above 2.28 keV_nr -- a real gap, since
    that fallback region's share of the total signal rate grows from ~6% at
    4 GeV to ~54% at 10 GeV WIMP mass (the higher the WIMP mass, the more of
    the recoil spectrum extends into it). Below 0.68 keV_nr: same
    linear-to-zero taper as chavarria_table_yield. Above 20.67 keV_nr (no
    measurement from either dataset): Lindhard fallback, same convention as
    chavarria_table_yield -- but note that region is a small and
    ever-shrinking fraction of the total rate (0.09% at 10 GeV per the
    same check above), unlike the now-closed 2.28-20.67 keV_nr gap.
    """
    global _chavarria_interp, _izraelevitch_interp
    if _chavarria_interp is None:
        _chavarria_interp = _load_chavarria_table()
    if _izraelevitch_interp is None:
        _izraelevitch_interp = _load_izraelevitch_table()

    scalar_input = np.ndim(E_r_keV) == 0
    E_r = np.atleast_1d(np.asarray(E_r_keV, dtype=float))
    E_ee = np.empty_like(E_r)

    below = E_r < _CHAVARRIA_ER_MIN_KEV
    chavarria_range = (E_r >= _CHAVARRIA_ER_MIN_KEV) & (E_r < _CHAVARRIA_ER_MAX_KEV)
    izraelevitch_range = (E_r >= _CHAVARRIA_ER_MAX_KEV) & (E_r <= _IZRAELEVITCH_ER_MAX_KEV)
    above = E_r > _IZRAELEVITCH_ER_MAX_KEV

    E_ee[chavarria_range] = _chavarria_interp(E_r[chavarria_range])
    E_ee[izraelevitch_range] = _izraelevitch_interp(E_r[izraelevitch_range])

    # Below measured range: linear to zero at the 082006-stated threshold
    # (same convention as chavarria_table_yield).
    Ee_at_Ermin = float(_chavarria_interp(_CHAVARRIA_ER_MIN_KEV))
    slope = Ee_at_Ermin / (_CHAVARRIA_ER_MIN_KEV - _CHAVARRIA_ZERO_CROSSING_KEV)
    below_taper = below & (E_r >= _CHAVARRIA_ZERO_CROSSING_KEV)
    below_zero = below & (E_r < _CHAVARRIA_ZERO_CROSSING_KEV)
    E_ee[below_taper] = slope * (E_r[below_taper] - _CHAVARRIA_ZERO_CROSSING_KEV)
    E_ee[below_zero] = 0.0

    # Above both measured ranges: Lindhard fallback (small and shrinking
    # share of the total rate -- see docstring).
    E_ee[above] = lindhard_yield(E_r[above])

    return float(E_ee[0]) if scalar_input else E_ee


_julian_interp = None  # lazily built on first use
_julian_er_min_kev = None  # exact table endpoints, read from the data itself
_julian_er_max_kev = None  # (avoids boundary mismatches from rounded constants)


# ----------------------------------------------------------------------------
# _load_julian_table
#   Interpolator E_ee(E_R) from the central photoneutron-iteration table (2024) plus the E_R range it covers.
# ----------------------------------------------------------------------------
def _load_julian_table():
    path = CIO.data_path(__file__, "data", "julian_photoneutron_2024_iteration_central.csv")
    data = np.loadtxt(path, delimiter=",", comments="#")
    Er, Ee = data[:, 0], data[:, 1]
    interp = interp1d(Er, Ee, kind="linear", bounds_error=False, fill_value=np.nan)
    return interp, float(Er.min()), float(Er.max())


def julian_table_yield(E_r_keV):
    """
    Measured Si nuclear recoil ionization efficiency, unpublished DAMIC-M
    photo-neutron calibration (J. Cuevas-Zepeda, 2024/2025). See module
    docstring and data/julian_photoneutron_2024_iteration_central.csv for
    provenance and the extrapolation convention (clamped to 0 below the
    measured range, Lindhard fallback above it).
    """
    global _julian_interp, _julian_er_min_kev, _julian_er_max_kev
    if _julian_interp is None:
        _julian_interp, _julian_er_min_kev, _julian_er_max_kev = _load_julian_table()

    scalar_input = np.ndim(E_r_keV) == 0
    E_r = np.atleast_1d(np.asarray(E_r_keV, dtype=float))
    E_ee = np.empty_like(E_r)

    below = E_r < _julian_er_min_kev
    above = E_r > _julian_er_max_kev
    within = ~below & ~above

    E_ee[within] = _julian_interp(E_r[within])
    E_ee[below] = 0.0
    E_ee[above] = lindhard_yield(E_r[above])

    return float(E_ee[0]) if scalar_input else E_ee


QUENCHING_MODELS = {
    "lindhard": lindhard_yield,
    "chavarria_table": chavarria_table_yield,
    "chavarria_izraelevitch_table": chavarria_izraelevitch_table_yield,
    "julian_table": julian_table_yield,
}
