#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — migdal_generate_grid
#  Grid driver that expands a (mchi, sigma_n) JSON grid and writes dR/dE_e
#  CSV rate tables for the Migdal effect via darkelf.
#
#  Author: Diego Venegas-Vargas
# ============================================================================

"""
Expand a (mchi_MeV, sigma_n_cm2) grid into Migdal rate CSV files.

Sigma_n fast path (linear rescaling)
--------------------------------------
For fixed mchi, the Migdal ionization probability I(omega) is independent of
sigma_n.  The rate is simply:

    R(sigma_n) = sigma_n / sigma_n_ref * R(sigma_n_ref)

so we call darkelf once per mchi (at sigma_n_ref = sigma_n_list[0] that still
needs a file) and rescale all other sigma_n values analytically.  This is
exact because sigma_n enters only as an overall prefactor in dRdomega_migdal.

JSON keys
---------
material, mediator ("heavy" | "light"), darkelf_dir (or CCDARK_SENS_DARKELF_DIR),
rates_dir, filename_template,
grid: {"mchi_MeV": ..., "sigma_n_cm2": ...},
detector: {"band_gap_eV": null, "Emin_eV": 0.0, "Emax_eV": 20.0, "binsize_eV": 0.1},
options: {"skip_existing": true, "overwrite": false, "parallel": 0,
          "format": {"mchi": ".6f", "sigma": ".1e"}}.

Usage::

    python3 utils/migdal_generate_grid.py configs/migdal_generate_si_heavy.json
"""
from __future__ import annotations

import contextlib
import json
import os
import sys
from pathlib import Path
from multiprocessing import Pool

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qedark_generate_grid import _expand_axis  # shared grid-axis parser
from ccdarkphys.common import io as CIO
from ccdarkphys.migdal.entry import compute_dRdE, _header_lines

# Map material name → nuclear target label used in filenames.
# For compounds the label is the material name itself (no single nucleus).
_MATERIAL_TO_TARGET_NUCLEUS = {
    "si": "Si28",
    "ge": "Ge74",
    "srcd2sb2": "SrCd2Sb2",
}


def _build_out_path(base_dir: Path, filename_template: str, material: str,
                    mediator: str, mchi_str: str, sigma_str: str) -> Path:
    target_nucleus = _MATERIAL_TO_TARGET_NUCLEUS.get(material.lower(), material)
    fname = filename_template.format(
        material=material,
        target_nucleus=target_nucleus,
        mediator=mediator,
        mchi_MeV=mchi_str,
        sigma_n_cm2=sigma_str,
    )
    base_dir.mkdir(parents=True, exist_ok=True)
    return base_dir / fname


def _exists_any(out_path: Path) -> bool:
    return out_path.exists() or Path(str(out_path) + ".gz").exists()


def _one_mass_task(task: dict) -> tuple[bool, str]:
    """
    Process all sigma_n values for one mchi point.

    Calls darkelf once at sigma_ref = sigma_items[0], then rescales remaining
    sigma values linearly — exact because R ∝ sigma_n.
    """
    material = task["material"]
    mediator = task["mediator"]
    mchi_eV = task["mchi_eV"]
    mchi_str = task["mchi_str"]
    sigma_items = task["sigma_items"]  # list of (sigma_val, sigma_str, out_path)
    Emin_eV = task["Emin_eV"]
    Emax_eV = task["Emax_eV"]
    binsize_eV = task["binsize_eV"]
    darkelf_dir = task["darkelf_dir"]
    idx = task.get("index", 0)
    n_tasks = task.get("n_tasks", 0)
    prefix = f"[{idx}/{n_tasks}] " if idx and n_tasks else ""

    sigma_ref_val, sigma_ref_str, _ = sigma_items[0]

    print(
        f"[migdal-grid] {prefix}mchi={mchi_str} MeV — "
        f"init at sigma_ref={sigma_ref_str}, then {len(sigma_items)} linear rescales",
        flush=True,
    )

    errors = []

    # One darkelf call at the reference sigma_n
    try:
        with open(os.devnull, "w") as devnull:
            with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                res_ref = compute_dRdE(
                    material=material,
                    mchi_eV=mchi_eV,
                    sigma_n_cm2=float(sigma_ref_val),
                    mediator=mediator,
                    Emin_eV=Emin_eV,
                    Emax_eV=Emax_eV,
                    binsize_eV=binsize_eV,
                    darkelf_dir=darkelf_dir,
                )
    except Exception as e:
        return False, f"[migdal-grid][ERROR] {prefix}mchi={mchi_str}: darkelf call failed: {e}"

    E_ref = res_ref["E_eV"]
    drde_ref = res_ref["dRdE_kg_year_eV"]
    meta_ref = res_ref["meta"]

    for sigma_val, sigma_str, out_path in sigma_items:
        try:
            scale = float(sigma_val) / float(sigma_ref_val)
            drde = drde_ref * scale

            meta = dict(meta_ref)
            meta["sigma_n_cm2"] = float(sigma_val)

            CIO.write_csv_generic(str(out_path), E_ref, drde, _header_lines(meta))
        except Exception as e:
            errors.append(
                f"[migdal-grid][ERROR] {prefix}mchi={mchi_str} sigma={sigma_str}: {e}"
            )

    if errors:
        return False, "\n".join(errors)
    return True, ""


def run_from_config(cfg: dict) -> None:
    material = cfg.get("material", "Si")
    mediator = cfg.get("mediator", "heavy")
    darkelf_dir = cfg.get("darkelf_dir") or os.environ.get("CCDARK_SENS_DARKELF_DIR")
    if not darkelf_dir:
        raise FileNotFoundError(
            "darkelf_dir is required (set in config or CCDARK_SENS_DARKELF_DIR env var)."
        )

    det = cfg.get("detector", {})
    Emin_eV = float(det.get("Emin_eV", 0.0))
    Emax_eV = float(det.get("Emax_eV", 20.0))
    binsize_eV = float(det.get("binsize_eV", 0.1))

    rates_dir = Path(cfg["rates_dir"])
    templ = cfg["filename_template"]

    opts = cfg.get("options", {})
    skip_existing = bool(opts.get("skip_existing", True))
    overwrite = bool(opts.get("overwrite", False))
    parallel = int(opts.get("parallel", 0))
    fmt = opts.get("format", {"mchi": ".6f", "sigma": ".1e"})
    fmt_mchi = fmt.get("mchi", ".6f")
    fmt_sigma = fmt.get("sigma", ".1e")

    gspec = cfg["grid"]
    mchi_list_MeV = _expand_axis(gspec["mchi_MeV"], "mchi_MeV")
    sigma_list = _expand_axis(gspec["sigma_n_cm2"], "sigma_n_cm2")

    tasks = []
    for mchi_MeV in mchi_list_MeV:
        mchi_eV = mchi_MeV * 1.0e6
        mchi_str = format(mchi_MeV, fmt_mchi)
        sigma_items = []
        for sigma in sigma_list:
            sigma_str = format(sigma, fmt_sigma)
            out_path = _build_out_path(rates_dir, templ, material, mediator, mchi_str, sigma_str)
            if (skip_existing or not overwrite) and _exists_any(out_path):
                continue
            sigma_items.append((sigma, sigma_str, out_path))
        if sigma_items:
            tasks.append({
                "material": material,
                "mediator": mediator,
                "mchi_eV": mchi_eV,
                "mchi_str": mchi_str,
                "sigma_items": sigma_items,
                "Emin_eV": Emin_eV,
                "Emax_eV": Emax_eV,
                "binsize_eV": binsize_eV,
                "darkelf_dir": darkelf_dir,
            })

    if not tasks:
        print("[migdal-grid] all files already exist, nothing to do.")
        return

    n_tasks = len(tasks)
    total_sigma = sum(len(t["sigma_items"]) for t in tasks)
    print(
        f"[migdal-grid] {n_tasks} mass points × up to {len(sigma_list)} sigma_n values "
        f"= {total_sigma} files to write  (1 darkelf init per mass)"
    )
    for i, task in enumerate(tasks, 1):
        task["index"] = i
        task["n_tasks"] = n_tasks

    if parallel and parallel > 1:
        with Pool(parallel) as pool:
            results = pool.map(_one_mass_task, tasks)
    else:
        results = [_one_mass_task(t) for t in tasks]

    for ok, msg in results:
        if not ok and msg:
            print(msg, file=sys.stderr, flush=True)


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: python3 utils/migdal_generate_grid.py <config.json>")
        sys.exit(1)
    cfg_path = Path(sys.argv[1]).expanduser()
    with open(cfg_path) as f:
        cfg = json.load(f)
    run_from_config(cfg)


if __name__ == "__main__":
    main()
