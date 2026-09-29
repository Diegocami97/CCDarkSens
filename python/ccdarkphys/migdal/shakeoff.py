# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: shakeoff.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  shakeoff.py -- Port of the DAMIC-M collaboration's own Migdal-effect rate
#  calculation (collab_frameworks/dim, scripts/calculate_rates/migdal.py,
#  git commit c790012), used to produce the published DAMIC-M PRL 2025
#  heavy-mediator Migdal exclusion curve. This is a *different* calculation
#  from ccdarkphys.migdal.entry's darkelf.dRdomega_migdal(method="Lindhard")
#  path: it computes the electron shake-off probability directly from a
#  tabulated dielectric function (Eq. 18 of Phys. Rev. D 105, 015014 (2022),
#  Eq. 5.7 of JHEP 01 (2023) 023) rather than going through DarkELF's own
#  Migdal machinery, and only the heavy-mediator (F_DM=1) case is implemented
#  -- the reference script's light-mediator path is an incomplete, commented-
#  out sketch, so I have not ported it. Use the darkelf-Lindhard path in
#  entry.py for the light mediator.
#
#  I ported this essentially verbatim, including a couple of quirks in the
#  original that I chose to preserve rather than "fix", since the point is
#  to match the reference curve as closely as possible, not to write my own
#  improved version of it:
#    - dPdomega's kmin/kmax handling: the caller always ends up with kmin=1
#      (Python's `not 0` is True, so passing kmin=0 explicitly still falls
#      through to the kmin=1 default) and kmax=22000 fixed by hand, which is
#      larger than the tabulated dielectric function's actual k range
#      (~10^4). scipy's RectBivariateSpline extrapolates rather than raising
#      past that edge; I keep the same 22000 value for consistency with
#      whatever curve was actually produced and published.
#    - mNSi = 28.0855e9 eV: this reads as the Si atomic weight in amu with
#      "e9" tacked on as an approximate amu->eV conversion (the precise
#      factor is 9.31494e8 eV/amu, about 7% smaller) -- an actual imprecision
#      in the original, kept as-is for the same reason.
# ============================================================================

"""
Port of the DAMIC-M collaboration's own Migdal-effect (heavy mediator) rate
calculation from collab_frameworks/dim/scripts/calculate_rates/migdal.py.

Public:
    compute_dRdE(material, mchi_eV, sigma_n_cm2, Emin_eV, Emax_eV, binsize_eV,
                 darkelf_dir=None)
"""

from __future__ import annotations

import os
from math import pi, sqrt, exp

import numpy as np
import pandas as pd
from scipy.special import erf
from scipy.integrate import quad
from scipy.interpolate import interp1d, RectBivariateSpline


# speed of light [cm/s]
_C_CMS = 3e10
# reference nuclear recoil energy [eV] the shake-off probability table is
# tabulated at; dRdomegadv rescales to the physical Er via (Er^2)/(2*Ennorm)
_E_N_NORMALIZED = 100.0

_CACHE: dict[str, "_ShakeoffModel"] = {}


# ----------------------------------------------------------------------------
# ShakeOffProbability
#   Electron shake-off probability dP/domega from a tabulated dielectric
#   function and Zion(k), following Eq. 18 of Phys. Rev. D 105, 015014
#   (2022) / Eq. 0 of JHEP 01 (2023) 023.
# ----------------------------------------------------------------------------
class ShakeOffProbability:
    def __init__(self, dielectric_file: str, zion_file: str):
        Re_eps, Im_eps, omega_arr, k_arr = self._read_dielectric(dielectric_file)
        self.Re_eps = Re_eps
        self.Im_eps = Im_eps
        self.kmin, self.kmax = k_arr.min(), k_arr.max()
        self.omega_min, self.omega_max = omega_arr.min(), omega_arr.max()
        self.zion = self._read_zion(zion_file)

    def ELF(self, omega, k):
        """Energy loss function Im[-1/eps] = Im[eps] / |eps|^2."""
        re = self.Re_eps(omega, k, grid=False)
        im = self.Im_eps(omega, k, grid=False)
        return im / (re ** 2 + im ** 2)

    @staticmethod
    def _read_dielectric(dielectric_file: str):
        data = pd.read_csv(dielectric_file, sep=r"\s+", header=None, skiprows=1,
                            names=["omega", "k", "eps1", "eps2"])
        df_R = data.pivot(index="omega", columns="k", values="eps1")
        df_I = data.pivot(index="omega", columns="k", values="eps2")
        R_arr = np.array(df_R)
        I_arr = np.array(df_I)
        omega_arr = df_R.index.values
        k_arr = df_R.columns.values
        Re_epsilon = RectBivariateSpline(omega_arr, k_arr, R_arr)
        Im_epsilon = RectBivariateSpline(omega_arr, k_arr, I_arr)
        return Re_epsilon, Im_epsilon, omega_arr, k_arr

    @staticmethod
    def _read_zion(zion_file: str):
        zion_data = np.loadtxt(zion_file, skiprows=1).T
        k_arr = zion_data[0]
        zion_arr = zion_data[1]
        return interp1d(k_arr, zion_arr, fill_value=(zion_arr[0], zion_arr[-1]), bounds_error=False)

    def dPdomegadk(self, omega, k, E_N):
        """dP/(domega dk) shake-off probability density [1/eV^2]."""
        m_N = 2.632e10  # mass of Si atom [eV]
        alpha = 1.0 / 137.0
        numerator = 4 * alpha * E_N
        denominator = 3 * pi ** 2 * omega ** 4 * m_N
        integrand = self.zion(k) ** 2 * k ** 2 * self.ELF(omega, k)
        return numerator / denominator * integrand

    def dPdomega(self, omega, E_N, kmin=None, kmax=None):
        """dP/domega [1/eV], integrated over k in log space."""
        if not kmax:
            kmax = self.kmax
        if not kmin:
            kmin = 1

        def log_dRdomegadk(logk):
            k = 10 ** logk
            prefactor = k * np.log(10)
            result = prefactor * self.dPdomegadk(omega, k, E_N)
            return max(result, 0)

        result, _ = quad(log_dRdomegadk, kmin, np.log10(kmax), limit=100)
        return result


# ----------------------------------------------------------------------------
# Migdal
#   Heavy-mediator (F_DM=1) Migdal-effect differential rate dR/domega,
#   following Eq. 5.7 of JHEP 01 (2023) 023.
# ----------------------------------------------------------------------------
class Migdal:
    def __init__(self, dPdomega_table, rho_X=0.3e9, v0=238e5, vE=253.7e5, vesc=544e5):
        self.dPdomega = dPdomega_table
        self.rhoX = rho_X
        self.v0 = v0 / _C_CMS
        self.vE = vE / _C_CMS
        self.vesc = vesc / _C_CMS
        self.vmax = self.vesc + self.vE

        self.sigman = 1.0  # cm^2, rate is linear in sigma_n so this is a unit reference
        self.ASi = 28
        self.mNSi = 28.0855e9  # eV (see module docstring: an approximate amu->eV factor)
        self.NTSi = 1.0 / ((self.mNSi * 1e-9) * 1.66e-27)  # Si atoms/kg

    def get_Er(self, v, omega, mX, mN):
        """Both analytic solutions of the recoil-energy quadratic."""
        denominator = mN ** 2 + 2 * mN * mX + mX ** 2
        sqrt_argument = mN * mX ** 3 * (mN * mX * v ** 2 - 2 * mN * omega - 2 * mX * omega)
        sqrt_term = np.sqrt(np.maximum(sqrt_argument, 0))
        Er1 = (mN * mX ** 2 * v ** 2 - mN * mX * omega - mX ** 2 * omega - v * sqrt_term) / denominator
        Er2 = (mN * mX ** 2 * v ** 2 - mN * mX * omega - mX ** 2 * omega + v * sqrt_term) / denominator
        return Er1, Er2

    def get_vminr(self, omega, mX, mN):
        """Minimum DM velocity for a given omega (from Er1 == Er2)."""
        return np.sqrt(2 * (mN + mX) * omega / (mN * mX))

    def fXint(self, v, vmins):
        """Integrated Standard Halo Model (Maxwell-Boltzmann) velocity distribution."""
        ve = self.vE
        vesc = self.vesc
        v0 = self.v0
        N0 = pi ** 1.5 * v0 ** 2 * (v0 * erf(vesc / v0) - 2 * (vesc / sqrt(pi)) * exp(-(vesc ** 2) / v0 ** 2))
        return v0 ** 2 / (2 * N0 * ve * v) * (
            np.exp(-((v - ve) ** 2) / v0 ** 2) * np.heaviside(vesc - (v - ve), 1)
            - np.exp(-((v + ve) ** 2) / v0 ** 2) * np.heaviside(vesc - (v + ve), 1)
        )

    def dRdomegadv(self, mX, omega, v):
        """Rate integrand before the velocity integral (Eq. 5.7, JHEP 01 (2023) 023)."""
        vminr = self.get_vminr(omega, mX, self.mNSi)
        mu_xnSi = 1 / (1 / mX + 1 / (self.mNSi / self.ASi))
        prefactor = ((self.mNSi / (2 * _E_N_NORMALIZED)) * self.NTSi * (self.rhoX / mX)
                     * (self.sigman / 2 / mu_xnSi ** 2) * _C_CMS * (60 * 60 * 24 * 366) * self.ASi ** 2 * 2 * pi)
        fxi = self.fXint(v, vminr)
        Er1, Er2 = self.get_Er(v, omega, mX, self.mNSi)
        return prefactor * fxi * v * (Er2 ** 2 - min(Er1 ** 2, Er2 ** 2)) * self.dPdomega(omega)

    def dRdomega(self, mX, omega):
        """dR/domega [events/kg/year/eV] at reference sigma_n = 1 cm^2."""
        vminr = self.get_vminr(omega, mX, self.mNSi)

        def integrand(v):
            return max(self.dRdomegadv(mX, omega, v), 0)

        rate, _ = quad(integrand, vminr, self.vmax, limit=100)
        return rate


# ----------------------------------------------------------------------------
# _ShakeoffModel
#   Lazily-built, process-cached (ShakeOffProbability, dPdomega table)
#   pair -- the expensive dielectric-file load and dP/domega tabulation
#   are mass-independent, so I build them once per process and reuse
#   across every mass point in a generation run.
# ----------------------------------------------------------------------------
class _ShakeoffModel:
    def __init__(self, dielectric_file: str, zion_file: str, om_min: float, om_max: float):
        self.sop = ShakeOffProbability(dielectric_file, zion_file)
        om_list = np.arange(max(om_min, 0.1), om_max + 0.2, 0.1)
        dPdom = np.array([self.sop.dPdomega(om, _E_N_NORMALIZED, 0, 22000) for om in om_list])
        self.dPdom_table = interp1d(om_list, dPdom, bounds_error=False, fill_value=0.0)


def _resolve_data_files(darkelf_dir: str) -> tuple[str, str]:
    dielectric_file = os.path.join(darkelf_dir, "data", "Si", "Si_gpaw_withLFE.dat")
    zion_file = os.path.join(darkelf_dir, "data", "Si", "Si_Zion.dat")
    if not os.path.isfile(dielectric_file):
        raise FileNotFoundError(dielectric_file)
    if not os.path.isfile(zion_file):
        raise FileNotFoundError(zion_file)
    return dielectric_file, zion_file


def _get_model(darkelf_dir: str, om_min: float, om_max: float) -> _ShakeoffModel:
    key = f"{darkelf_dir}:{om_min}:{om_max}"
    if key not in _CACHE:
        dielectric_file, zion_file = _resolve_data_files(darkelf_dir)
        _CACHE[key] = _ShakeoffModel(dielectric_file, zion_file, om_min, om_max)
    return _CACHE[key]


# ----------------------------------------------------------------------------
# compute_dRdE
#   Same outward contract as ccdarkphys.migdal.entry.compute_dRdE, restricted
#   to the heavy mediator (F_DM=1) -- the reference implementation this is
#   ported from has no working light-mediator path.
# ----------------------------------------------------------------------------
def compute_dRdE(
    material: str,
    mchi_eV: float,
    sigma_n_cm2: float,
    Emin_eV: float = 0.0,
    Emax_eV: float = 20.0,
    binsize_eV: float = 0.1,
    *,
    darkelf_dir: str | None = None,
) -> dict:
    if material.lower() not in ("si", "silicon"):
        raise NotImplementedError("Only Silicon is supported by the shake-off Migdal port.")
    if not darkelf_dir:
        darkelf_dir = os.environ.get("CCDARK_SENS_DARKELF_DIR")
    if not darkelf_dir:
        raise FileNotFoundError("darkelf_dir is required (for the Si_gpaw_withLFE.dat / Si_Zion.dat data files).")

    n_bins = max(1, int(round((Emax_eV - Emin_eV) / binsize_eV)))
    e_edges = np.linspace(Emin_eV, Emax_eV, n_bins + 1)
    e_centers = 0.5 * (e_edges[:-1] + e_edges[1:])

    model = _get_model(darkelf_dir, float(e_centers.min()), float(e_centers.max()))
    migdal = Migdal(model.dPdom_table, rho_X=0.3e9, v0=238e5, vE=253.7e5, vesc=544e5)

    # dR/domega scales linearly in sigma_n; Migdal.dRdomega is at sigma_n = 1 cm^2.
    drde = np.array([migdal.dRdomega(float(mchi_eV), float(om)) for om in e_centers]) * float(sigma_n_cm2)

    meta = {
        "material": material,
        "mediator": "heavy",
        "darkelf_mediator": "n/a (shakeoff, not darkelf)",
        "rate_method": "shakeoff",
        "reference": "collab_frameworks/dim scripts/calculate_rates/migdal.py (JHEP 01 (2023) 023 Eq. 5.7)",
        "mchi_eV": float(mchi_eV),
        "sigma_n_cm2": float(sigma_n_cm2),
        "Emin_eV": float(Emin_eV),
        "Emax_eV": float(Emax_eV),
        "binsize_eV": float(binsize_eV),
        "v0_kms": 238.0,
        "vE_kms": 253.7,
        "vesc_kms": 544.0,
        "darkelf_dir": darkelf_dir,
    }
    return {
        "E_eV": np.asarray(e_centers, dtype=float),
        "dRdE_kg_year_eV": np.asarray(drde, dtype=float),
        "meta": meta,
    }
