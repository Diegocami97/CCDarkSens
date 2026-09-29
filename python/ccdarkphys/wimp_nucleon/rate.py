# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  rate.py -- Diego Venegas-Vargas DAMIC-M collaboration CCDarkSens Framework
#  rate.py -- Spin-independent WIMP-nucleus elastic scattering rate, dR/dE_R.
# ============================================================================

"""
Spin-independent WIMP-nucleus elastic scattering rate, dR/dE_R.

Ported from WIMPyCCD's analysis/rate_models.py (nuclear_form_factor,
get_rate_ben; "Ben Loer's thesis sec. 1.2.2"). Kept numerically faithful to
WIMPyCCD -- see utils/validate_wimp_nucleon_rate_parity.py, which reproduces
WIMPyCCD's own printed reference number bit-for-bit -- rather than rederived
from a textbook formula, per the decision to match WIMPyCCD exactly for this
first pass.

Two items were found during the post-parity sanity pass; both are now
resolved (fixed below):

1. RESOLVED -- form factor bug: WIMPyCCD's own nuclear_form_factor carries
   an unresolved author comment ("check why it's always so close to !???
   compare with F from Levin"). Root cause: q = sqrt(2*m_T*E_R) is a
   MOMENTUM in SI units (kg*m/s), but it's multiplied directly by r_n
   (in meters) as if it were a wavenumber (1/m) -- missing a division by
   hbar. The resulting q*r_n is always ~1e-35, deep in the spherical
   Bessel function's small-argument limit j1(x)/x -> 1/3, so F evaluates
   to EXACTLY 1.000000 for every E_R and every A -- not "close to 1",
   identically 1, always. That's exactly the bug the author's own comment
   was flagging. Confirmed numerically (see conversation/plan record):
   for Si (A=28) the corrected form factor is 0.9998 at E_R=0.1 keV down
   to 0.9525 at E_R=30 keV -- physically sensible, mild coherent-limit
   suppression -- versus the buggy version's flat 1.0 everywhere. Since
   the rate scales as F^2, the bug overestimates the rate by up to ~9%
   at 30 keV (negligible near threshold, where a 1-10 GeV Si search's
   rate concentrates anyway).

   nuclear_form_factor() below is now the CORRECTED version (q properly
   converted to a wavenumber via q/hbar before the Bessel call) and is
   the production default. The literal WIMPyCCD behavior (F==1 always)
   is kept as nuclear_form_factor_wimpyccd_literal() purely so
   utils/validate_wimp_nucleon_rate_parity.py still has an exact,
   unambiguous target to check the port against -- it should never be
   used for physics.

2. RESOLVED -- reduced-mass bug: WIMPyCCD's DMParams.mu_n is built from
   m_T_kg (the full nucleus mass, A*m_proton), i.e. it is the DM-NUCLEUS
   reduced mass (mu_chiN), not DM-nucleon (mu_chin), despite the "_n" name,
   and the ORIGINAL code uses that same mu_chiN in TWO places: (a) v_min
   inside the halo integral, and (b) the explicit 1/(2 A mu^2 m_chi) rate
   prefactor.

   Re-deriving the standard SI rate from scratch (sigma_n as a per-nucleon
   cross section, isospin-conserving coupling): dsigma/dE_R =
   (m_T/(2 mu_chiN^2 v^2)) * sigma_0 * F^2, and sigma_0 = sigma_n *
   (mu_chiN/mu_chin)^2 * A^2 converts the nucleus-level cross section to
   the per-nucleon convention. Substituting, mu_chiN CANCELS OUT of the
   rate prefactor entirely -- mu_chin is what's left explicit. mu_chiN
   only legitimately survives through v_min's kinematics (a real fact
   about the 2-body WIMP-nucleus collision), inside eta. So usage (a) was
   always correct; usage (b) was the bug.

   Magnitude, Si (A=28), across this project's 1-10 GeV range: using
   mu_chiN instead of mu_chin in the prefactor makes WIMPyCCD's rate too
   LOW by 75% at 1 GeV, 88% at 2 GeV, 97% at 5 GeV, 99% at 10 GeV (rate
   scales as 1/mu^2, and mu_chin/mu_chiN shrinks fast as m_chi approaches
   and exceeds the nucleon mass) -- a mass-dependent, order-of-magnitude
   effect, not a small correction like item 1.

   dRdE_nr_kg_day_keV() below is now the CORRECTED production version:
   mu_chi_nucleus_kg (mu_chiN) is still used for v_min/eta (that usage was
   always right), but a separate mu_chi_nucleon_kg (mu_chin, built from a
   fixed nucleon mass rather than the full nucleus mass) is used in the
   explicit rate prefactor. The full original WIMPyCCD formula -- BOTH
   bugs (item 1's F==1 and item 2's mu_chiN-in-the-prefactor) reproduced
   exactly -- is preserved as dRdE_nr_kg_day_keV_wimpyccd_literal(), kept
   ONLY so this module can reproduce what the original WIMPyCCD repo
   actually computes (for provenance, comparison, or re-deriving a past
   result) and so utils/validate_wimp_nucleon_rate_parity.py keeps an
   exact, unambiguous target. Never use the _wimpyccd_literal functions
   for physics.
"""
from __future__ import annotations

import numpy as np
from scipy.special import spherical_jn

from ccdarkphys.wimp_nucleon.halo import (
    RHO_CHI_GEV_CM3_DEFAULT,
    V0_KMS_DEFAULT,
    VE_KMS_DEFAULT,
    VESC_KMS_DEFAULT,
    mean_inverse_speed,
)

# WIMPyCCD's DMParams constants, ported exactly (including the two slightly
# different keV->J constants used in different places -- see module docstring
# in halo.py for the velo_int one; nuclear_form_factor and get_rate_ben each
# use their own literal below, kept distinct from halo.py's 1.602e-16 for
# bit-for-bit parity with WIMPyCCD).
_N0_PER_KG = 6.02e26            # Avogadro's number, "per kg at A=1" convention
_PROTON_MASS_GEV = 0.938
_GEV_TO_KG = 1.78e-27
_RHO_GEV_CM3_TO_SI = 1.783e-21  # GeV/cm^3 -> kg/m^3
_KEV_TO_J_FORMFACTOR = 1.6e-16  # nuclear_form_factor's own constant
_KEV_TO_J_RATE = 1.6e-16        # get_rate_ben's final unit-conversion constant
_SECONDS_PER_DAY = 24 * 3600
_HBAR_JS = 1.054571817e-34      # J*s -- needed to convert a SI momentum into a wavenumber (1/m)


# ----------------------------------------------------------------------------
# _nucleus_mass_kg
#   Nucleus mass in kg, approximated as A times the proton mass.
# ----------------------------------------------------------------------------
def _nucleus_mass_kg(A: float) -> float:
    return A * _PROTON_MASS_GEV * _GEV_TO_KG


# ----------------------------------------------------------------------------
# _nucleon_mass_kg
#   Nucleon (proton) mass in kg.
# ----------------------------------------------------------------------------
def _nucleon_mass_kg() -> float:
    return _PROTON_MASS_GEV * _GEV_TO_KG


# ----------------------------------------------------------------------------
# _reduced_mass_kg
#   Reduced mass m1*m2/(m1 + m2) in kg.
# ----------------------------------------------------------------------------
def _reduced_mass_kg(m1_kg: float, m2_kg: float) -> float:
    return m1_kg * m2_kg / (m1_kg + m2_kg)


def _q_momentum_kg_m_s(E_R_keV, A: float):
    """sqrt(2 m_T E_R): a MOMENTUM (kg*m/s), not yet a wavenumber."""
    E_R_J = np.asarray(E_R_keV, dtype=float) * _KEV_TO_J_FORMFACTOR
    m_T_kg = _nucleus_mass_kg(A)
    return np.sqrt(2.0 * m_T_kg * E_R_J)


def nuclear_form_factor(E_R_keV, A: float):
    """
    Simplified Helm form factor -- CORRECTED (production default). Same
    functional form and same rn/s parameters as WIMPyCCD's
    nuclear_form_factor (rn = 1.14*A^(1/3) fm used directly as the
    spherical-Bessel argument, no Lewin-Smith "R1=sqrt(rn^2-5s^2)"
    correction, s = 0.9 fm), but q is properly converted from a momentum
    to a wavenumber (q/hbar) before multiplying by r_n, fixing the bug
    documented in the module docstring's item 1.
    """
    q_momentum = _q_momentum_kg_m_s(E_R_keV, A)
    q_wavenumber = q_momentum / _HBAR_JS  # 1/m
    r_n = 1.14 * A ** (1.0 / 3.0) * 1e-15  # m
    s = 0.9 * 1e-15  # m
    x = q_wavenumber * r_n
    return 3.0 * spherical_jn(1, x) / x * np.exp(-((q_wavenumber * s) ** 2) / 2.0)


def nuclear_form_factor_wimpyccd_literal(E_R_keV, A: float):
    """
    Literal, BUGGY port of WIMPyCCD's nuclear_form_factor -- q (a momentum,
    kg*m/s) multiplied directly by r_n (m) without converting to a
    wavenumber first, so the spherical-Bessel argument is never
    dimensionless and always sits in the x->0 limit: F is identically
    1.000000 for every E_R and A. Kept ONLY so
    utils/validate_wimp_nucleon_rate_parity.py has an exact target to
    check the port against -- never use this for physics. See module
    docstring, item 1.
    """
    q_momentum = _q_momentum_kg_m_s(E_R_keV, A)
    r_n = 1.14 * A ** (1.0 / 3.0) * 1e-15  # m
    s = 0.9 * 1e-15  # m
    x = q_momentum * r_n  # NOT dimensionless -- this is the bug, preserved for parity
    return 3.0 * spherical_jn(1, x) / x * np.exp(-((q_momentum * s) ** 2) / 2.0)


def dRdE_nr_kg_day_keV(
    E_R_keV,
    mchi_GeV: float,
    sigma_n_cm2: float,
    A: float,
    *,
    rho_chi_gev_cm3: float = RHO_CHI_GEV_CM3_DEFAULT,
    v0_kms: float = V0_KMS_DEFAULT,
    vE_kms: float = VE_KMS_DEFAULT,
    vesc_kms: float = VESC_KMS_DEFAULT,
):
    """
    dR/dE_R for SI WIMP-nucleus elastic scattering, in events/kg/day/keV
    (nuclear recoil energy axis). PRODUCTION version: corrected form
    factor (module docstring item 1) and correct DM-nucleon reduced mass
    in the rate prefactor (item 2), while still using the DM-nucleus
    reduced mass for v_min/eta -- that usage was always physically
    correct. sigma_n_cm2 here is a genuine per-nucleon cross section, as
    the naming implies.

    For a byte-for-byte reproduction of the ORIGINAL WIMPyCCD repo (both
    known bugs included), see dRdE_nr_kg_day_keV_wimpyccd_literal below.
    """
    m_T_kg = _nucleus_mass_kg(A)
    m_chi_kg = mchi_GeV * _GEV_TO_KG
    mu_chi_nucleus_kg = _reduced_mass_kg(m_chi_kg, m_T_kg)          # for v_min / eta kinematics -- always correct
    mu_chi_nucleon_kg = _reduced_mass_kg(m_chi_kg, _nucleon_mass_kg())  # for the rate prefactor -- the fix
    rho_SI = rho_chi_gev_cm3 * _RHO_GEV_CM3_TO_SI
    sigma_m2 = sigma_n_cm2 * 1.0e-4  # cm^2 -> m^2

    F2 = nuclear_form_factor(E_R_keV, A) ** 2
    eta = mean_inverse_speed(
        E_R_keV, m_T_kg, mu_chi_nucleus_kg, v0_kms=v0_kms, vE_kms=vE_kms, vesc_kms=vesc_kms
    )

    rate_SI = (
        (_N0_PER_KG * sigma_m2 * rho_SI * m_T_kg) / (2.0 * A * mu_chi_nucleon_kg ** 2 * m_chi_kg)
        * A ** 2 * F2 * eta
    )
    return rate_SI * _SECONDS_PER_DAY * _KEV_TO_J_RATE


def dRdE_nr_kg_day_keV_wimpyccd_literal(
    E_R_keV,
    mchi_GeV: float,
    sigma_n_cm2: float,
    A: float,
    *,
    rho_chi_gev_cm3: float = RHO_CHI_GEV_CM3_DEFAULT,
    v0_kms: float = V0_KMS_DEFAULT,
    vE_kms: float = VE_KMS_DEFAULT,
    vesc_kms: float = VESC_KMS_DEFAULT,
):
    """
    Byte-for-byte reproduction of WIMPyCCD's get_rate_ben, BOTH known bugs
    included: nuclear_form_factor_wimpyccd_literal (F==1 always, item 1)
    and the DM-NUCLEUS reduced mass used in the explicit rate prefactor
    where the standard convention needs DM-nucleon (item 2). This is the
    function utils/validate_wimp_nucleon_rate_parity.py checks against
    WIMPyCCD's own printed reference number, and the one to reach for if
    you ever need to reproduce a number the original WIMPyCCD repo would
    have produced. NEVER use this for physics -- see module docstring.
    """
    m_T_kg = _nucleus_mass_kg(A)
    m_chi_kg = mchi_GeV * _GEV_TO_KG
    mu_chi_nucleus_kg = _reduced_mass_kg(m_chi_kg, m_T_kg)  # used in BOTH places below -- this is the bug
    rho_SI = rho_chi_gev_cm3 * _RHO_GEV_CM3_TO_SI
    sigma_m2 = sigma_n_cm2 * 1.0e-4  # cm^2 -> m^2

    F2 = nuclear_form_factor_wimpyccd_literal(E_R_keV, A) ** 2
    eta = mean_inverse_speed(
        E_R_keV, m_T_kg, mu_chi_nucleus_kg, v0_kms=v0_kms, vE_kms=vE_kms, vesc_kms=vesc_kms
    )

    rate_SI = (
        (_N0_PER_KG * sigma_m2 * rho_SI * m_T_kg) / (2.0 * A * mu_chi_nucleus_kg ** 2 * m_chi_kg)
        * A ** 2 * F2 * eta
    )
    return rate_SI * _SECONDS_PER_DAY * _KEV_TO_J_RATE
