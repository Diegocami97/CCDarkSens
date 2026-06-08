"""
DM–electron differential rates for crystal scattering (QCDark-style kernel).

Adapted from the physics in the reference QCDark ``dark_matter_rates`` module
(eta_MB, momentum integral, prefactor chain). Implemented in-tree so CCDarkSens
does not depend on ``collab_frameworks`` or external checkout paths at import time.

Halo **parameters** (``v0``, ``vEarth``, ``vEscape`` in km/s; ``rhoX`` in eV/cm³)
are passed through ``astro_model`` from ``qcdark.entry`` — use the same numerical
values as in your QEDark JSON if you want matched astrophysical inputs; the SHM
η factor here remains the reference QCDark ``eta_MB`` formulation (distinct from
``qedark.entry``'s ``eta_shm_numeric``).
"""

from __future__ import annotations

import numpy as np
from scipy import special
from scipy.integrate import simpson

from ccdarkphys.common import constants as QEC

# Reference QCDark conventions (km/s for halo inputs in ``astro_model``)
_LIGHT_SPEED_KM_S = 299792.458
_PI = np.pi

# Same as reference: 1 cm in seconds, 1 s in years (for unit conversion of d_rate)
_CM2SEC = 1.0 / _LIGHT_SPEED_KM_S * 1e-5
_SEC2YR = 1.0 / (60.0 * 60.0 * 24.0 * 365.25)


def _reduced_mass_mXe(m_x_eV: float) -> float:
    me = float(QEC.me_eV)
    return (m_x_eV * me) / (m_x_eV + me)


def eta_mb(q_arr: np.ndarray, e_ev: float, m_x_eV: float, astro_model: dict) -> np.ndarray:
    """
    η for the truncated Maxwellian (same piecewise structure as reference QCDark).

    ``astro_model`` keys: v0, vEarth, vEscape in **km/s**; rhoX in eV/cm³; sigma_e in cm².
    q_arr: momentum-transfer array in **eV** (same convention as reference integrand).
    Returns η in (cm/s)⁻¹-compatible scaling as in the reference (uses c=1 velocity algebra).
    """
    v_esc = astro_model["vEscape"] / _LIGHT_SPEED_KM_S
    v_e = astro_model["vEarth"] / _LIGHT_SPEED_KM_S
    v_0 = astro_model["v0"] / _LIGHT_SPEED_KM_S

    q_flat = np.asarray(q_arr, dtype=float).ravel()
    val = np.zeros_like(q_flat, dtype=float)
    for i, q in enumerate(q_flat):
        v_min = q / (2.0 * m_x_eV) + e_ev / q
        if v_min < v_esc - v_e:
            val[i] = (
                -4.0 * v_e * np.exp(-((v_esc / v_0) ** 2))
                + np.sqrt(_PI)
                * v_0
                * (
                    special.erf((v_min + v_e) / v_0)
                    - special.erf((v_min - v_e) / v_0)
                )
            )
        elif v_min < v_esc + v_e:
            val[i] = (
                -2.0 * (v_e + v_esc - v_min) * np.exp(-((v_esc / v_0) ** 2))
                + np.sqrt(_PI)
                * v_0
                * (special.erf(v_esc / v_0) - special.erf((v_min - v_e) / v_0))
            )
        else:
            val[i] = 0.0

    kk = (v_0**3) * (
        -2.0 * _PI * (v_esc / v_0) * np.exp(-((v_esc / v_0) ** 2))
        + (_PI**1.5) * special.erf(v_esc / v_0)
    )
    return ((v_0**2) * _PI / (2.0 * v_e * kk)) * val


def f_dm(q_arr: np.ndarray, fdm_exp: int) -> np.ndarray:
    """|F_DM(q)| as in reference (heavy: exponent 0 → 1)."""
    q = np.maximum(np.asarray(q_arr, dtype=float), 1e-30)
    return (QEC.alpha * QEC.me_eV / q) ** fdm_exp


def tf_screening(q_arr: np.ndarray, e_ev: float, screening: dict) -> np.ndarray:
    """Thomas–Fermi-style screening factor (1 if disabled)."""
    if not screening.get("DoScreen", False):
        return np.ones_like(np.asarray(q_arr, dtype=float), dtype=float)
    eps0 = screening["eps0"]
    alpha_s = screening["alphaS"]
    q_tf = screening["qTF"]
    omega_p = screening["omegaP"]
    me = float(QEC.me_eV)
    q = np.asarray(q_arr, dtype=float)
    val = (
        1.0 / (eps0 - 1.0)
        + alpha_s * ((q / q_tf) ** 2)
        + q**4 / (4.0 * (me**2) * (omega_p**2))
        - (e_ev / omega_p) ** 2
    )
    return 1.0 / (1.0 + 1.0 / val)


def _lindhard_f(u: np.ndarray, z: np.ndarray) -> np.ndarray:
    return 0.5 + (1.0 / (8.0 * z)) * (
        (1.0 - (z - u) ** 2) * np.log((z - u + 1.0) / (z - u - 1.0))
        + (1.0 - (z + u) ** 2) * np.log((z + u + 1.0) / (z + u - 1.0))
    )


def lindhard_screening(q_arr: np.ndarray, e_ev: float, _screening: dict) -> np.ndarray:
    """
    Lindhard dielectric response, matching the pydme QCDark branch:
      epsilon = 1 + epsilon_pre * f(u, z)
      Screening factor in rate integrand is |1/epsilon|^2.
    """
    omega_p = 16.6  # eV, pydme default for Si
    alpha_s = float(QEC.alpha)
    me = float(QEC.me_eV)
    v_f = ((3.0 * _PI * omega_p**2) / (4.0 * alpha_s * me**2)) ** (1.0 / 3.0)
    q = np.asarray(q_arr, dtype=float)
    z = q / (2.0 * me * v_f)
    # small imaginary regulator to avoid branch singularities (same idea as pydme)
    u = (float(e_ev) + 1j * 1e-8) / (q * v_f)
    epsilon_pre = (3.0 * omega_p**2) / ((q * v_f) ** 2)
    epsilon = 1.0 + epsilon_pre * _lindhard_f(u, z)
    return 1.0 / epsilon


def screening_factor(q_arr: np.ndarray, e_ev: float, screening: dict) -> np.ndarray:
    if not screening.get("DoScreen", False):
        return np.ones_like(np.asarray(q_arr, dtype=float), dtype=float)
    method = str(screening.get("method", "ThomasFermi")).lower()
    if method == "lindhard":
        return lindhard_screening(q_arr, e_ev, screening)
    return tf_screening(q_arr, e_ev, screening)


def _momentum_integrand(
    dq: float,
    d_e: float,
    q_index: np.ndarray,
    e_index: int,
    m_x_eV: float,
    f_crystal2: np.ndarray,
    fdm_exp: int,
    screening: dict,
    astro_model: dict,
) -> np.ndarray:
    q_arr = dq * q_index + dq / 2.0
    e_val = d_e * e_index + d_e / 2.0
    ff_slice = f_crystal2[q_index, e_index]
    eta = eta_mb(q_arr, e_val, m_x_eV, astro_model)
    fd = f_dm(q_arr, fdm_exp)
    scr = screening_factor(q_arr, e_val, screening)
    return (
        e_val / (q_arr**2) * eta * (fd**2) * ff_slice * (np.abs(scr) ** 2)
    )


def _d_rate_fixed_e(
    dq: float,
    d_e: float,
    m_cell: float,
    e_index: int,
    m_x_eV: float,
    f_crystal2: np.ndarray,
    fdm_exp: int,
    screening: dict,
    astro_model: dict,
) -> float:
    rho_x = astro_model["rhoX"]
    sigma_e = astro_model["sigma_e"]
    prefactor = (
        (rho_x / m_x_eV)
        * (5.609588e35 / m_cell)
        * sigma_e
        * QEC.alpha
        * ((QEC.me_eV / _reduced_mass_mXe(m_x_eV)) ** 2)
    )
    qi = np.arange(np.shape(f_crystal2)[0], dtype=int)
    y = _momentum_integrand(
        dq, d_e, qi, e_index, m_x_eV, f_crystal2, fdm_exp, screening, astro_model
    )
    x_q = dq * qi.astype(float) + dq / 2.0
    return float(prefactor * simpson(y, x=x_q))


def differential_rate_spectrum(
    m_x_eV: float,
    dq: float,
    d_e: float,
    m_cell: float,
    f_crystal2: np.ndarray,
    fdm_exp: int,
    screening: dict,
    astro_model: dict,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (E_bin_centers_eV, dR/dE in events / kg / year / eV).

    Uses the same normalization chain as reference ``d_rate``:
    ``vals / (cm2sec * sec2yr * E)``.
    """
    num_e = int(np.shape(f_crystal2)[1])
    vals = np.empty(num_e, dtype=float)
    for j in range(num_e):
        vals[j] = _d_rate_fixed_e(
            dq,
            d_e,
            m_cell,
            j,
            m_x_eV,
            f_crystal2,
            fdm_exp,
            screening,
            astro_model,
        )
    e_arr = np.arange(num_e, dtype=float) * d_e + d_e / 2.0
    d_rd_e = vals / (_CM2SEC * _SEC2YR * e_arr)
    d_rd_e[~np.isfinite(d_rd_e)] = 0.0
    return e_arr, d_rd_e
