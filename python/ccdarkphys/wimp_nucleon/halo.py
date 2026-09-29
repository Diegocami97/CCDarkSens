# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  halo.py -- Standard Halo Model velocity integral (mean inverse speed
#  eta(v_min)) for spin-independent WIMP-nucleus elastic scattering.
# ============================================================================

"""
Standard Halo Model velocity integral (mean inverse speed eta(v_min)) for
spin-independent WIMP-nucleus elastic scattering.

Ported from WIMPyCCD's analysis/rate_models.py::velo_int ("Ben Loer's thesis
sec. 1.2.2"), generalized to take explicit halo parameters instead of a
DMParams object so it can be driven from CCDarkSens JSON configs. The
piecewise truncated-Maxwell-Boltzmann form is standard (Lewin & Smith 1996 /
Savage, Freese & Gondolo 2006) -- only the parametrization was changed here,
not the physics.

Default halo parameters follow Baxter et al. 2021 (EPJC 81, 907, Table 1),
matching the values already used by ccdarkphys.migdal.entry for darkelf
calls: v0 = 238 km/s, vE = 253.7 km/s (annual-averaged |v0 + v_sun|),
vesc = 544 km/s, rho_chi = 0.3 GeV/cm^3. WIMPyCCD's own defaults (vE=263,
rho=0.3) differ slightly -- kept overridable per-call so the parity check in
utils/validate_wimp_nucleon_rate_parity.py can reproduce WIMPyCCD's own
printed numbers using WIMPyCCD's own parameters, while production configs
use the Baxter/Migdal-consistent defaults below.
"""
from __future__ import annotations

import numpy as np
from scipy.special import erf

V0_KMS_DEFAULT = 238.0
VE_KMS_DEFAULT = 253.7
VESC_KMS_DEFAULT = 544.0
RHO_CHI_GEV_CM3_DEFAULT = 0.3

# WIMPyCCD's velo_int uses 1.602e-16 J/keV for its own internal E_R -> SI
# conversion (nuclear_form_factor uses a slightly different 1.6e-16 -- see
# rate.py). Kept exactly as WIMPyCCD has it for bit-for-bit parity.
_KEV_TO_J = 1.602e-16


def v_min_mps(E_R_keV, m_T_kg: float, mu_kg: float):
    """
    Minimum WIMP speed able to produce nuclear recoil energy E_R.

    mu_kg is the DM-nucleus reduced mass (see note in rate.py about
    WIMPyCCD naming this `mu_n` despite it being nucleus-, not nucleon-,
    level) -- kept consistent with how WIMPyCCD's velo_int uses it.
    """
    E_R_J = np.asarray(E_R_keV, dtype=float) * _KEV_TO_J
    return np.sqrt(m_T_kg * E_R_J / 2.0) / mu_kg


def mean_inverse_speed(
    E_R_keV,
    m_T_kg: float,
    mu_kg: float,
    *,
    v0_kms: float = V0_KMS_DEFAULT,
    vE_kms: float = VE_KMS_DEFAULT,
    vesc_kms: float = VESC_KMS_DEFAULT,
):
    """
    eta(v_min) = <1/v> over the truncated Maxwell-Boltzmann distribution,
    boosted into the Earth frame. Exact port of WIMPyCCD's velo_int (same
    three-case piecewise structure; the v_min > vesc + vE case is implicitly
    zero via the np.zeros_like initialization, same as the original).

    Units: SI in (km/s halo params get converted internally), SI out (s/m).
    """
    v0 = v0_kms * 1.0e3
    vE = vE_kms * 1.0e3
    vesc = vesc_kms * 1.0e3

    inverse_k = (np.pi ** 1.5) * v0 ** 3 * (
        erf(vesc / v0) - (2.0 * vesc / (np.sqrt(np.pi) * v0)) * np.exp(-(vesc ** 2) / v0 ** 2)
    )
    k = 1.0 / inverse_k

    scalar_input = np.ndim(E_R_keV) == 0
    vmin = np.atleast_1d(v_min_mps(E_R_keV, m_T_kg, mu_kg)).astype(float)
    eta = np.zeros_like(vmin)

    prefactor = (np.pi ** 1.5) * v0 ** 3 * k / (2.0 * vE)

    # Case 1: v_min <= v_esc - v_E
    mask1 = vmin <= (vesc - vE)
    eta[mask1] = (
        erf((vE - vmin[mask1]) / v0)
        + erf((vE + vmin[mask1]) / v0)
        - (4.0 * vE / (np.sqrt(np.pi) * v0)) * np.exp(-(vesc ** 2) / v0 ** 2)
    )

    # Case 2: v_esc - v_E < v_min <= v_esc + v_E
    mask2 = (vmin > (vesc - vE)) & (vmin <= (vesc + vE))
    eta[mask2] = (
        erf((vE - vmin[mask2]) / v0)
        + erf(vesc / v0)
        - (2.0 * (vE + vesc - vmin[mask2]) / (np.sqrt(np.pi) * v0)) * np.exp(-(vesc ** 2) / v0 ** 2)
    )

    # Case 3: v_min > v_esc + v_E -> eta stays 0

    result = prefactor * eta
    return float(result[0]) if scalar_input else result
