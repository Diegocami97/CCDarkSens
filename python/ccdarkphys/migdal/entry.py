"""
DarkELF entry point for the Migdal effect.

Supports two material classes:
  - Single-nucleus targets (e.g. "Si"): standard darkelf path, YAML + ELF loaded normally.
  - Compound targets (e.g. "srcd2sb2"): multi-nucleus summation with Lindhard electron
    response. Each nuclear species is handled by a separate darkelf object initialised
    from the Si YAML (ELF file unused with method="Lindhard"), with A, mN, and Zion
    overridden per species. The compound rate is the mass-fraction-weighted sum:
        dR/dω = Σ_ν  f_ν × dR_ν/dω  [events / kg_compound / yr / eV]
    where f_ν = N_ν·A_ν / M_cell.

Public:
    compute_dRdE(material, mchi_eV, sigma_n_cm2, mediator="heavy",
                 Emin_eV=0.0, Emax_eV=20.0, binsize_eV=0.1,
                 darkelf_dir=None, darkelf_kwargs=None)
      -> dict(E_eV, dRdE_kg_year_eV, meta)

CLI example:
    PYTHONPATH=python python3 -m ccdarkphys.migdal.entry \
      --material Si --mchi_MeV 1000 --sigma_n_cm2 1e-38 --mediator heavy \
      --darkelf_dir /path/to/DarkELF \
      --out_csv data/migdal_rates/Si/heavy/dRdE_Si28_heavy_m1000.000000_s1e-38.csv
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

from ccdarkphys.common import io as CIO
from ccdarkphys.common.mediator_map import MEDIATOR_TO_FDM_INDEX

_ENV_DARKELF_DIR = "CCDARK_SENS_DARKELF_DIR"

# darkelf source-file mapping per single-nucleus material.
_DARKELF_FILES = {
    "si": dict(target="Si", filename="Si_mermin.dat", phonon_filename="Si_epsphonon_data6K.dat"),
}

# Compound materials: list of nuclear species with their properties.
#   A    — average atomic mass (amu); determines nuclear recoil kinematics and A² coherence factor.
#   Zion — free electrons per atom for the Lindhard electron response (valence-only convention).
#           Uncertainty bracket for Cd: Zion=2 (valence) vs Zion=12 (4d¹⁰5s²).
#   N    — stoichiometric count per formula unit.
# The Si YAML is used to initialise darkelf (ELF file not read with method="Lindhard");
# A, mN, and Zion are overridden after construction — Si object is never modified.
_COMPOUND_NUCLEI: dict[str, list[dict]] = {
    "srcd2sb2": [
        {"name": "Sr", "A": 87.62,  "Zion": 2,  "N": 1},
        {"name": "Cd", "A": 112.41, "Zion": 2,  "N": 2},
        {"name": "Sb", "A": 121.76, "Zion": 5,  "N": 2},
    ],
}


def _resolve_darkelf_dir(explicit: str | None) -> str:
    if explicit:
        return os.path.abspath(os.path.expanduser(explicit))
    env = os.environ.get(_ENV_DARKELF_DIR)
    if env:
        return os.path.abspath(os.path.expanduser(env.strip()))
    raise FileNotFoundError(
        f"darkelf source directory is required. Pass darkelf_dir=... or set {_ENV_DARKELF_DIR}."
    )


def _load_darkelf(darkelf_dir: str):
    if darkelf_dir not in sys.path:
        sys.path.insert(0, darkelf_dir)
    try:
        from darkelf import darkelf as darkelf_cls
    except ImportError as exc:
        raise ImportError(
            f"darkelf package not importable from {darkelf_dir!r}. "
            f"Pass darkelf_dir=... pointing at the darkelf repo root, or set {_ENV_DARKELF_DIR}."
        ) from exc
    return darkelf_cls


def _make_energy_grid(Emin_eV: float, Emax_eV: float, binsize_eV: float):
    n_bins = max(1, int(round((Emax_eV - Emin_eV) / binsize_eV)))
    e_edges = np.linspace(Emin_eV, Emax_eV, n_bins + 1)
    return 0.5 * (e_edges[:-1] + e_edges[1:])


def _call_migdal(obj, e_centers, sigma_n_cm2: float) -> np.ndarray:
    # Enth=0: include all nuclear recoil energies. The default Enth=4*ombar≈0.12 eV
    # kills all phase space for mchi < 15 MeV on Si (E_R,max < 0.12 eV at those masses).
    # method="Lindhard": analytic Thomas-Fermi approximation — avoids k-grid boundary
    # artifacts at low mchi and requires only Zion (no ELF data file).
    return obj.dRdomega_migdal(e_centers, sigma_n=float(sigma_n_cm2),
                               Enth=0.0, method="Lindhard")


def compute_dRdE(
    material: str,
    mchi_eV: float,
    sigma_n_cm2: float,
    mediator: str = "heavy",
    Emin_eV: float = 0.0,
    Emax_eV: float = 20.0,
    binsize_eV: float = 0.1,
    *,
    darkelf_dir: str | None = None,
    darkelf_kwargs: dict | None = None,
) -> dict:
    """
    Same outward contract as other CCDarkSens backends.

    Returns:
      dict with ``E_eV``, ``dRdE_kg_year_eV`` (events / kg / year / eV), ``meta``.
    """
    key = material.lower()
    if mediator not in MEDIATOR_TO_FDM_INDEX:
        raise ValueError("mediator: heavy/massive/0 or light/massless/2")
    darkelf_mediator = "massive" if MEDIATOR_TO_FDM_INDEX[mediator] == 0 else "massless"

    darkelf_path = _resolve_darkelf_dir(darkelf_dir)
    darkelf_cls = _load_darkelf(darkelf_path)
    e_centers = _make_energy_grid(Emin_eV, Emax_eV, binsize_eV)

    # ── Compound material path ────────────────────────────────────────────────
    if key in _COMPOUND_NUCLEI:
        nuclei = _COMPOUND_NUCLEI[key]
        M_cell = sum(nuc["N"] * nuc["A"] for nuc in nuclei)
        si_kwargs = dict(_DARKELF_FILES["si"])  # Si YAML satisfies constructor; ELF unused
        drde_total = np.zeros_like(e_centers)
        for nuc in nuclei:
            f_nu = nuc["N"] * nuc["A"] / M_cell
            obj = darkelf_cls(mX=float(mchi_eV), v0kms=238.0, vekms=253.7,
                              vesckms=544.0, **si_kwargs)
            obj.rhoX = 0.3e9
            # Override nuclear parameters BEFORE update_params so muxN is correct.
            obj.A    = nuc["A"]
            obj.mN   = nuc["A"] * obj.mp
            obj.Zion = nuc["Zion"]
            obj.update_params(mX=float(mchi_eV), mediator=darkelf_mediator)
            drde_total += f_nu * _call_migdal(obj, e_centers, sigma_n_cm2)
        meta = {
            "material": material,
            "compound_nuclei": [n["name"] for n in nuclei],
            "mass_fractions": {n["name"]: round(n["N"] * n["A"] / M_cell, 4) for n in nuclei},
            "Zion_per_nucleus": {n["name"]: n["Zion"] for n in nuclei},
            "mediator": mediator,
            "darkelf_mediator": darkelf_mediator,
            "mchi_eV": float(mchi_eV),
            "sigma_n_cm2": float(sigma_n_cm2),
            "Emin_eV": float(Emin_eV),
            "Emax_eV": float(Emax_eV),
            "binsize_eV": float(binsize_eV),
            "darkelf_dir": darkelf_path,
        }
        return {
            "E_eV": np.asarray(e_centers, dtype=float),
            "dRdE_kg_year_eV": np.asarray(drde_total, dtype=float),
            "meta": meta,
        }

    # ── Single-nucleus path (Si, Ge, …) — unchanged ──────────────────────────
    if key not in _DARKELF_FILES:
        raise NotImplementedError(
            f"No darkelf file mapping for material={material!r}. "
            f"Known single-nucleus: {list(_DARKELF_FILES)}; "
            f"known compounds: {list(_COMPOUND_NUCLEI)}."
        )
    file_kwargs = dict(_DARKELF_FILES[key])
    if darkelf_kwargs:
        file_kwargs.update(darkelf_kwargs)

    # Halo parameters from Baxter et al. 2021 (EPJC 81, 907)
    obj = darkelf_cls(mX=float(mchi_eV), v0kms=238.0, vekms=253.7,
                      vesckms=544.0, **file_kwargs)
    obj.rhoX = 0.3e9
    obj.update_params(mX=float(mchi_eV), mediator=darkelf_mediator)
    drde = _call_migdal(obj, e_centers, sigma_n_cm2)

    meta = {
        "material": material,
        "mediator": mediator,
        "darkelf_mediator": darkelf_mediator,
        "mchi_eV": float(mchi_eV),
        "sigma_n_cm2": float(sigma_n_cm2),
        "Emin_eV": float(Emin_eV),
        "Emax_eV": float(Emax_eV),
        "binsize_eV": float(binsize_eV),
        "E_gap_darkelf_eV": float(getattr(obj, "E_gap", float("nan"))),
        "darkelf_dir": darkelf_path,
    }
    return {
        "E_eV": np.asarray(e_centers, dtype=float),
        "dRdE_kg_year_eV": np.asarray(drde, dtype=float),
        "meta": meta,
    }


def _header_lines(meta: dict) -> list:
    lines = [
        "# Differential Rates computed with CCDarkSens (DarkELF-Migdal entry)",
        f"# material = {meta['material']}, mediator = {meta['mediator']} "
        f"(darkelf: {meta['darkelf_mediator']})",
        f"# mchi (eV) = {meta['mchi_eV']}",
        f"# sigma_n (cm^2, DM-nucleon) = {meta['sigma_n_cm2']}",
        f"# Emin_eV = {meta['Emin_eV']}, Emax_eV = {meta['Emax_eV']}, binsize_eV = {meta['binsize_eV']}",
    ]
    if "compound_nuclei" in meta:
        lines.append(f"# compound nuclei = {meta['compound_nuclei']}, "
                     f"mass fractions = {meta['mass_fractions']}, "
                     f"Zion = {meta['Zion_per_nucleus']}")
        lines.append("# rate = sum_nu f_nu * dR_nu/domega  [events/kg_compound/yr/eV]  (Lindhard, Enth=0)")
    else:
        lines.append(f"# darkelf E_gap (eV) = {meta.get('E_gap_darkelf_eV', 'n/a')}")
    lines += [
        f"# darkelf_dir = {meta['darkelf_dir']}",
        "# Output units: dR/dE_e in events / kg / year / eV (post-Migdal ionization spectrum)",
        "# Columns: E (eV), dRdE (events/kg/year/eV)",
    ]
    return lines


def _cli():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--material", default="Si")
    ap.add_argument("--mchi_MeV", type=float, required=True)
    ap.add_argument("--sigma_n_cm2", type=float, required=True)
    ap.add_argument("--mediator", default="heavy", choices=["heavy", "massive", "light", "massless"])
    ap.add_argument("--Emin_eV", type=float, default=0.0)
    ap.add_argument("--Emax_eV", type=float, default=20.0)
    ap.add_argument("--binsize_eV", type=float, default=0.1)
    ap.add_argument("--darkelf_dir", default=None)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    res = compute_dRdE(
        material=args.material,
        mchi_eV=args.mchi_MeV * 1.0e6,
        sigma_n_cm2=args.sigma_n_cm2,
        mediator=args.mediator,
        Emin_eV=args.Emin_eV,
        Emax_eV=args.Emax_eV,
        binsize_eV=args.binsize_eV,
        darkelf_dir=args.darkelf_dir,
    )
    CIO.write_csv_generic(args.out_csv, res["E_eV"], res["dRdE_kg_year_eV"], _header_lines(res["meta"]))
    print(f"[migdal] wrote {args.out_csv}  (N={len(res['E_eV'])})")


if __name__ == "__main__":
    _cli()
