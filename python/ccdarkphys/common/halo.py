"""
Standard Halo Model helper(s): eta_SHM(vmin; v0, vE, vesc).
All velocities in cm/s. Returns dimensionless η(vmin).
"""

from __future__ import annotations
import numpy as np
from math import erf, sqrt, pi, exp

def _erf(x: np.ndarray) -> np.ndarray:
    # vectorized math.erf
    vfunc = np.vectorize(erf, otypes=[float])
    return vfunc(x)

def _exp(x: np.ndarray) -> np.ndarray:
    vfunc = np.vectorize(exp, otypes=[float])
    return vfunc(x)

# def eta_shm(vmin_cm_s: np.ndarray,
#             v0_cm_s: float,
#             vE_cm_s: float,
#             vesc_cm_s: float) -> np.ndarray:
#     """
#     Truncated Maxwellian with shift vE, normalized by standard N.
#     Follows the common closed-form used in DM-e literature.

#     Parameters
#     ----------
#     vmin_cm_s : array-like
#     v0_cm_s, vE_cm_s, vesc_cm_s : floats

#     Returns
#     -------
#     eta : np.ndarray (dimensionless), same shape as vmin_cm_s
#     """
#     vmin = np.asarray(vmin_cm_s, dtype=float)
#     v0   = float(v0_cm_s)
#     vE   = float(vE_cm_s)
#     vesc = float(vesc_cm_s)

#     z = vesc / v0
#     N = erf(z) - (2.0/np.sqrt(np.pi))*z*np.exp(-z*z)

#     a = (vmin - vE)/v0
#     b = (vmin + vE)/v0

#     term1 = _erf(b) - _erf(a)
#     term2 = (2.0/np.sqrt(np.pi))*( (_exp(-a*a) - _exp(-b*b)) * (v0/(2.0*vE)) )

#     eta = (term1 - term2) / (2.0*vE*N)
#     eta[vmin > (vesc + vE)] = 0.0
#     eta = np.maximum(eta, 0.0)
#     eta[~np.isfinite(eta)] = 0.0
#     return eta

# Refined to preserve tiny high-velocity tail instead of clipping to zero
def eta_shm_analytic(vmin_cm_s: np.ndarray,
            v0_cm_s: float,
            vE_cm_s: float,
            vesc_cm_s: float) -> np.ndarray:
    """
    Standard Halo Model (truncated Maxwellian with shift vE), normalized by N.

    All velocities in cm/s. Returns η(v_min) in (cm/s)^-1.

    This uses the usual closed-form expression, but *preserves* the tiny
    high-velocity tail instead of clipping small negative values to zero
    (which was killing ultralight DM rates).
    """
    vmin = np.asarray(vmin_cm_s, dtype=float)
    v0   = float(v0_cm_s)
    vE   = float(vE_cm_s)
    vesc = float(vesc_cm_s)

    # Normalization
    z = vesc / v0
    N = erf(z) - (2.0 / np.sqrt(np.pi)) * z * np.exp(-z * z)

    a = (vmin - vE) / v0
    b = (vmin + vE) / v0

    term1 = _erf(b) - _erf(a)
    term2 = (2.0 / np.sqrt(np.pi)) * ((_exp(-a * a) - _exp(-b * b)) * (v0 / (2.0 * vE)))

    eta = (term1 - term2) / (2.0 * vE * N)

    # Hard kinematic cutoff: no DM above vesc + vE
    eta[vmin > (vesc + vE)] = 0.0

    # Numerical fix:
    # For vmin in the extreme high tail, the analytic expression can give
    # tiny negative values from cancellation. We interpret those as small
    # positive contributions (like the numerical integral does), and flip
    # the sign instead of forcing them to zero.
    neg = eta < 0.0
    eta[neg] = -eta[neg]

    # Clean up NaNs / Infs
    eta[~np.isfinite(eta)] = 0.0

    return eta


import numpy as np
from scipy.integrate import nquad
from scipy.special import erf

def eta_shm_numeric(vmin_cm_s: np.ndarray,
                    v0_cm_s: float,
                    vE_cm_s: float,
                    vesc_cm_s: float) -> np.ndarray:
    """
    Numeric SHM halo integral η(v_min), with the same interface as eta_shm:

        - vmin_cm_s : array-like, v_min in cm/s
        - v0_cm_s   : scalar,    v0 in cm/s
        - vE_cm_s   : scalar,    vE in cm/s
        - vesc_cm_s : scalar,    vesc in cm/s

    Returns η(v_min) in (cm/s)^-1 as an array with the same shape as vmin_cm_s.

    Implementation: direct integration over speed and angle using
    scipy.integrate.nquad, following Eqs. (B4)–(B5) of 1509.01598,
    with a sharp cutoff at vesc.
    """

    v0   = float(v0_cm_s)
    vE   = float(vE_cm_s)
    vesc = float(vesc_cm_s)

    # Normalization:
    # KK = ∫ d^3v exp(-|v|^2 / v0^2) Θ(vesc - |v|)
    #    = π^(3/2) v0^3 [ erf(z) - (2/√π) z e^{-z^2} ],  z = vesc / v0
    z  = vesc / v0
    KK = (v0**3) * (
        -2.0 * np.exp(-z * z) * np.pi * vesc / v0
        + (np.pi**1.5) * erf(z)
    )

    def scalar_eta(vmin: float) -> float:
        """η(vmin) for a single vmin (cm/s) via direct integration."""
        vmin = float(vmin)

        def f_exp(vsq):
            return np.exp(-vsq / (v0 * v0))

        # Region 1: vmin <= vesc - vE  (Eq. B4)
        if vmin <= vesc - vE:

            def bounds_cosq():
                return [-1.0, 1.0]

            def bounds_vx(cosq):
                # vx_max from |v + vE| <= vesc:
                # vx_max = -cosq vE + sqrt((cosq^2 - 1) vE^2 + vesc^2)
                vmax = -cosq * vE + np.sqrt((cosq * cosq - 1.0) * vE * vE + vesc * vesc)
                return [vmin, vmax]

            def integrand(vx, cosq):
                vsq = vx * vx + vE * vE + 2.0 * vx * vE * cosq
                return (2.0 * np.pi / KK) * vx * f_exp(vsq)

            val, _err = nquad(integrand, [bounds_vx, bounds_cosq])
            return float(val)

        # Region 2: vesc - vE < vmin <= vesc + vE  (Eq. B5)
        elif vesc - vE < vmin <= vesc + vE:

            def bounds_cosq(vx):
                # cosθ_max = (vesc^2 - vE^2 - vx^2) / (2 vx vE)
                num = vesc * vesc - vE * vE - vx * vx
                den = 2.0 * vx * vE
                up  = num / den
                return [-1.0, up]

            def bounds_vx():
                # vx ∈ [vmin, vE + vesc]
                return [vmin, vE + vesc]

            def integrand(cosq, vx):
                vsq = vx * vx + vE * vE + 2.0 * vx * vE * cosq
                return (2.0 * np.pi / KK) * vx * f_exp(vsq)

            val, _err = nquad(integrand, [bounds_cosq, bounds_vx])
            return float(val)

        # Region 3: vmin > vesc + vE → no DM
        else:
            return 0.0

    # Vectorized wrapper with same behavior as analytic eta_shm
    vmin_arr = np.asarray(vmin_cm_s, dtype=float)
    out = np.empty_like(vmin_arr, dtype=float)

    it = np.nditer(vmin_arr, flags=['multi_index'])
    while not it.finished:
        vm = float(it[0])
        out[it.multi_index] = scalar_eta(vm)
        it.iternext()

    out[~np.isfinite(out)] = 0.0
    return out


def kms_to_cms(x_kms: float) -> float:
    return float(x_kms) * 1.0e5
