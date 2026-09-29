# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: entry.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  entry.py -- QCDark entry point for DM–electron scattering (Silicon crystal
#  |F|^2 from HDF5).
# ============================================================================

"""
QCDark entry point for DM–electron scattering (Silicon crystal |F|^2 from HDF5).

Uses an in-tree rate kernel compatible with the reference QCDark ``dark_matter_rates``
normalization (no imports from ``collab_frameworks``).

**Crystal table:** place the HDF5 file under ``<repo>/data/qcdark/`` (e.g.
``Si_f2_qcdark.h5`` or ``Si_final.hdf5`` copied from standalone QCDark), or set
``CCDARK_SENS_QCDARK_FORM_FACTOR``, or pass ``form_factor_h5=`` to ``compute_dRdE``.

**Contract:** Same outward API as ``ccdarkphys.qedark.entry.compute_dRdE`` — return
keys, halo dict (``v*_kms`` or ``v*_cm_s``), mediator labels, CSV metadata via
``common.io.write_csv(..., entry="QCDark")``.

Use the **same numerical halo inputs** as in QEDark configs (``v0_kms``, ``vE_kms``,
``vesc_kms``); local density is ``constants.rho_X_eVcm3``. Physics remains QCDark-native
(kernel, η, crystal HDF5). Downstream ROOT/JSON workflows only need ``rates_dir`` pointing
at these CSVs (same two-column ``E``, ``dRdE`` …/kg/year/eV format as QEDark).

Output grid follows the HDF5 binning (bin centers ``(i+1/2)*dE``); ``binsize_eV``
must equal the file's ``dE`` (validated).

Public:
    compute_dRdE(..., form_factor_h5=None)
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

from ccdarkphys.common import constants as QEC
from ccdarkphys.common import halo as HALO
from ccdarkphys.common import io as CIO
from ccdarkphys.common.mediator_map import MEDIATOR_TO_FDM_INDEX
from ccdarkphys.qcdark.form_factor import CrystalFormFactor
from ccdarkphys.qcdark import kernel as QCK

_ENV_FORM_FACTOR = "CCDARK_SENS_QCDARK_FORM_FACTOR"
_DEFAULT_SCREENING_SI = {
    "DoScreen": True,
    "method": "Lindhard",
    "eps0": 11.3,
    "qTF": 4.13e3,  # eV
    "omegaP": 16.6,  # eV
    "alphaS": 1.563,
}


def repo_qcdark_data_dir() -> Path:
    """``<repository_root>/data/qcdark`` — canonical place for crystal HDF5 files."""
    return Path(__file__).resolve().parents[3] / "data" / "qcdark"


# ----------------------------------------------------------------------------
# _resolve_form_factor_path
#   Find the crystal HDF5 file: the explicit argument first, then the CCDARK_SENS_QCDARK_FORM_FACTOR environment
#   variable, then the usual names under data/qcdark/, then the package data directory. Raises FileNotFoundError
#   with instructions if none exists.
# ----------------------------------------------------------------------------
def _resolve_form_factor_path(explicit: str | None) -> str:
    if explicit:
        return os.path.abspath(os.path.expanduser(explicit))
    env = os.environ.get(_ENV_FORM_FACTOR)
    if env:
        return os.path.abspath(os.path.expanduser(env.strip()))
    data_qcdark = repo_qcdark_data_dir()
    for name in ("Si_f2_qcdark.h5", "Si_final.hdf5", "demo_Si_f2_qcdark.h5"):
        cand = data_qcdark / name
        if cand.is_file():
            return str(cand.resolve())
    try:
        return CIO.data_path(__file__, "data", "Si_f2_qcdark.h5")
    except FileNotFoundError as exc:
        dq = repo_qcdark_data_dir()
        raise FileNotFoundError(
            "QCDark crystal HDF5 not found. Copy your standalone-QCDark output into "
            f"{dq}/ (e.g. Si_f2_qcdark.h5 or Si_final.hdf5), run ``python -m ccdarkphys.qcdark.demo_run`` "
            "for a synthetic demo_Si_f2_qcdark.h5 there, or set "
            f"{_ENV_FORM_FACTOR}. Layout: reference QCDark ``form_factor`` / ``results/f2``."
        ) from exc


# ----------------------------------------------------------------------------
# _halo_kms
#   Halo velocities (v0, vE, vesc) in km/s from either the *_cm_s or the *_kms keys of the halo dict.
# ----------------------------------------------------------------------------
def _halo_kms(halo: dict) -> tuple[float, float, float]:
    if "v0_cm_s" in halo:
        return (
            float(halo["v0_cm_s"]) / 1.0e5,
            float(halo["vE_cm_s"]) / 1.0e5,
            float(halo["vesc_cm_s"]) / 1.0e5,
        )
    return (
        float(halo["v0_kms"]),
        float(halo["vE_kms"]),
        float(halo["vesc_kms"]),
    )


def compute_dRdE(
    material: str,
    mediator: str,
    mchi_eV: float,
    sigma_e_cm2: float,
    halo: dict,
    band_gap_eV: float = 1.2,
    eh_pair_eV: float = 3.8,
    binsize_eV: float = 0.1,
    *,
    form_factor_h5: str | None = None,
) -> dict:
    """
    Same signature as ``qedark.entry.compute_dRdE``, plus optional keyword-only
    ``form_factor_h5`` for the crystal HDF5 path.

    Returns:
      dict with ``E_eV``, ``dRdE_kg_year_eV`` (events / kg / year / eV), ``meta``.
    """
    if material.lower() not in ("si", "silicon"):
        raise NotImplementedError("Only Silicon is supported (material='Si').")
    n_fdm = MEDIATOR_TO_FDM_INDEX.get(mediator)
    if n_fdm is None:
        raise ValueError("mediator: heavy/massive/0 or light/massless/2")

    path = _resolve_form_factor_path(form_factor_h5)
    ff = CrystalFormFactor(path)

    tol = max(1e-12, abs(ff.dE) * 1e-6)
    if abs(float(binsize_eV) - ff.dE) > tol:
        raise ValueError(
            f"binsize_eV ({binsize_eV}) must match crystal table dE={ff.dE} eV "
            "(set detector.binsize_eV in JSON to this value)."
        )

    v0_kms, v_e_kms, v_esc_kms = _halo_kms(halo)

    if "v0_cm_s" in halo:
        v0_cm_s = float(halo["v0_cm_s"])
        v_e_cm_s = float(halo["vE_cm_s"])
        v_esc_cm_s = float(halo["vesc_cm_s"])
    else:
        v0_cm_s = HALO.kms_to_cms(halo["v0_kms"])
        v_e_cm_s = HALO.kms_to_cms(halo["vE_kms"])
        v_esc_cm_s = HALO.kms_to_cms(halo["vesc_kms"])

    astro_model = {
        "v0": v0_kms,
        "vEarth": v_e_kms,
        "vEscape": v_esc_kms,
        "rhoX": float(QEC.rho_X_eVcm3),
        "sigma_e": float(sigma_e_cm2),
    }

    screening = dict(_DEFAULT_SCREENING_SI)

    e_ev, d_rd_e = QCK.differential_rate_spectrum(
        float(mchi_eV),
        ff.dq,
        ff.dE,
        ff.mCell,
        ff.ff,
        int(n_fdm),
        screening,
        astro_model,
    )

    meta = {
        "material": material,
        "mediator": mediator,
        "table_path": ff._source_path,
        "table_sha1": CIO.sha1sum(ff._source_path),
        "v0_cm_s": v0_cm_s,
        "vE_cm_s": v_e_cm_s,
        "vesc_cm_s": v_esc_cm_s,
        "mchi_eV": float(mchi_eV),
        "sigma_e_cm2": float(sigma_e_cm2),
        "eh_pair_eV": float(eh_pair_eV),
        "band_gap_eV": float(band_gap_eV),
        "band_gap_file_eV": float(ff.band_gap),
        "screening": str(screening["method"]),
        "screening_DoScreen": bool(screening["DoScreen"]),
        "screening_eps0": float(screening["eps0"]),
        "screening_qTF_eV": float(screening["qTF"]),
        "screening_omegaP_eV": float(screening["omegaP"]),
        "screening_alphaS": float(screening["alphaS"]),
    }

    return {
        "E_eV": e_ev,
        "dRdE_kg_year_eV": np.asarray(d_rd_e, dtype=float),
        "meta": meta,
    }


# ----------------------------------------------------------------------------
# _cli
#   Command-line front end: compute one rate table for the given mediator, mass and cross section and write it to --out_csv.
# ----------------------------------------------------------------------------
def _cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--material", default="Si")
    ap.add_argument(
        "--mediator", required=True, choices=["heavy", "massive", "light", "massless"]
    )
    ap.add_argument("--mchi_MeV", type=float, required=True)
    ap.add_argument("--sigma_e_cm2", type=float, required=True)
    ap.add_argument("--v0_kms", type=float, default=220.0)
    ap.add_argument("--vE_kms", type=float, default=232.0)
    ap.add_argument("--vesc_kms", type=float, default=544.0)
    ap.add_argument("--band_gap_eV", type=float, default=1.2)
    ap.add_argument("--eh_pair_eV", type=float, default=3.8)
    ap.add_argument("--binsize_eV", type=float, default=0.1)
    ap.add_argument(
        "--form_factor_h5",
        default=None,
        help=f"Override HDF5 path (else env {_ENV_FORM_FACTOR}, then data/qcdark/*.h5).",
    )
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    res = compute_dRdE(
        material=args.material,
        mediator=args.mediator,
        mchi_eV=args.mchi_MeV * 1.0e6,
        sigma_e_cm2=args.sigma_e_cm2,
        halo={"v0_kms": args.v0_kms, "vE_kms": args.vE_kms, "vesc_kms": args.vesc_kms},
        band_gap_eV=args.band_gap_eV,
        eh_pair_eV=args.eh_pair_eV,
        binsize_eV=args.binsize_eV,
        form_factor_h5=args.form_factor_h5,
    )
    e, r, meta = res["E_eV"], res["dRdE_kg_year_eV"], res["meta"]
    CIO.write_csv(args.out_csv, e, r, meta, entry="QCDark")
    print(f"[qcdark] wrote {args.out_csv}  (N={len(e)})")


if __name__ == "__main__":
    _cli()
