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

def eta_shm(vmin_cm_s: np.ndarray,
            v0_cm_s: float,
            vE_cm_s: float,
            vesc_cm_s: float) -> np.ndarray:
    """
    Truncated Maxwellian with shift vE, normalized by standard N.
    Follows the common closed-form used in DM-e literature.

    Parameters
    ----------
    vmin_cm_s : array-like
    v0_cm_s, vE_cm_s, vesc_cm_s : floats

    Returns
    -------
    eta : np.ndarray (dimensionless), same shape as vmin_cm_s
    """
    vmin = np.asarray(vmin_cm_s, dtype=float)
    v0   = float(v0_cm_s)
    vE   = float(vE_cm_s)
    vesc = float(vesc_cm_s)

    z = vesc / v0
    N = erf(z) - (2.0/np.sqrt(np.pi))*z*np.exp(-z*z)

    a = (vmin - vE)/v0
    b = (vmin + vE)/v0

    term1 = _erf(b) - _erf(a)
    term2 = (2.0/np.sqrt(np.pi))*( (_exp(-a*a) - _exp(-b*b)) * (v0/(2.0*vE)) )

    eta = (term1 - term2) / (2.0*vE*N)
    eta[vmin > (vesc + vE)] = 0.0
    eta = np.maximum(eta, 0.0)
    eta[~np.isfinite(eta)] = 0.0
    return eta

def kms_to_cms(x_kms: float) -> float:
    return float(x_kms) * 1.0e5
