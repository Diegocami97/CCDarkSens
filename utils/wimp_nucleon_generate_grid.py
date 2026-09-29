#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: wimp_nucleon_generate_grid.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: wimp_nucleon_generate_grid.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  wimp_nucleon_generate_grid.py -- Grid driver that expands a (mchi,
#  sigma_n) JSON grid and writes dR/dE_R CSV rate tables for spin-independent
#  WIMP-nucleus elastic scattering.
# ============================================================================

"""
Expand a (mchi_MeV, sigma_n_cm2) grid into WIMP-nucleon SI rate CSV files.
Structure mirrors migdal_generate_grid.py.

Sigma_n fast path (linear rescaling)
--------------------------------------
dR/dE_R is exactly linear in sigma_n (it enters only as an overall
prefactor -- see ccdarkphys/wimp_nucleon/rate.py), so, as in the Migdal
generator, we compute once per mchi at sigma_ref = sigma_list[0] and rescale
all other sigma_n values analytically:

    R(sigma_n) = sigma_n / sigma_n_ref * R(sigma_n_ref)

JSON keys
---------
target_nucleus, A, mediator (informational only -- SI is a contact
interaction, no light/heavy distinction),
recoil: {"Emin_keV":..., "Emax_keV":..., "nbins":...},
halo: {"rho_chi_gev_cm3":..., "v0_kms":..., "vE_kms":..., "vesc_kms":...}   (optional; Baxter/Migdal defaults if omitted)
quenching: {"model": "lindhard" | "chavarria_table",
            "ee_Emin_eV":..., "ee_Emax_eV":..., "ee_nbins":...}
           (OPTIONAL -- omit entirely to get the Phase 1 behavior: raw dR/dE_R
           CSVs, unquenched. If present, output is dR/dE_ee instead, quenched
           via ccdarkphys.wimp_nucleon.quenching.)
rates_dir, filename_template,
grid: {"mchi_MeV": ..., "sigma_n_cm2": ...},
options: {"skip_existing": true, "overwrite": false, "parallel": 0,
          "format": {"mchi": ".6f", "sigma": ".1e"}}.

Usage::

    python3 utils/wimp_nucleon_generate_grid.py configs/wimp_nucleon_generate_si_heavy.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from multiprocessing import Pool

_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT / "python"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from qedark_generate_grid import _expand_axis  # shared grid-axis parser
from ccdarkphys.common import io as CIO
from ccdarkphys.wimp_nucleon.entry import compute_dRdE, compute_dRdE_ee, _header_lines, _header_lines_ee


# ----------------------------------------------------------------------------
# _build_out_path
#   Output file path of one (mass, cross-section) point from the filename template; creates the directory.
# ----------------------------------------------------------------------------
def _build_out_path(base_dir: Path, filename_template: str, target_nucleus: str,
                    mediator: str, mchi_str: str, sigma_str: str) -> Path:
    fname = filename_template.format(
        target_nucleus=target_nucleus,
        mediator=mediator,
        mchi_MeV=mchi_str,
        sigma_n_cm2=sigma_str,
    )
    base_dir.mkdir(parents=True, exist_ok=True)
    return base_dir / fname


# ----------------------------------------------------------------------------
# _exists_any
#   True if the output file exists, either plain or gzip-compressed.
# ----------------------------------------------------------------------------
def _exists_any(out_path: Path) -> bool:
    return out_path.exists() or Path(str(out_path) + ".gz").exists()


def _one_mass_task(task: dict) -> tuple[bool, str]:
    """Process all sigma_n values for one mchi point: one rate call at
    sigma_ref, then exact linear rescaling for the rest."""
    target_nucleus = task["target_nucleus"]
    A = task["A"]
    mediator = task["mediator"]
    mchi_MeV = task["mchi_MeV"]
    mchi_str = task["mchi_str"]
    sigma_items = task["sigma_items"]  # list of (sigma_val, sigma_str, out_path)
    recoil = task["recoil"]
    halo = task["halo"]
    quenching = task["quenching"]  # None -> raw E_R output (Phase 1 behavior)
    idx = task.get("index", 0)
    n_tasks = task.get("n_tasks", 0)
    prefix = f"[{idx}/{n_tasks}] " if idx and n_tasks else ""

    sigma_ref_val, sigma_ref_str, _ = sigma_items[0]

    print(
        f"[wimp_nucleon-grid] {prefix}mchi={mchi_str} MeV — "
        f"init at sigma_ref={sigma_ref_str}, then {len(sigma_items)} linear rescales"
        + (f"  (quenching={quenching['model']})" if quenching else ""),
        flush=True,
    )

    common_kwargs = dict(
        A=A,
        mchi_MeV=mchi_MeV,
        sigma_n_cm2=float(sigma_ref_val),
        target_nucleus=target_nucleus,
        mediator=mediator,
        nr_Emin_keV=recoil["Emin_keV"],
        nr_Emax_keV=recoil["Emax_keV"],
        nr_nbins=recoil["nbins"],
        rho_chi_gev_cm3=halo.get("rho_chi_gev_cm3"),
        v0_kms=halo.get("v0_kms"),
        vE_kms=halo.get("vE_kms"),
        vesc_kms=halo.get("vesc_kms"),
    )

    try:
        if quenching:
            res_ref = compute_dRdE_ee(
                **common_kwargs,
                ee_Emin_eV=quenching["ee_Emin_eV"],
                ee_Emax_eV=quenching["ee_Emax_eV"],
                ee_nbins=quenching["ee_nbins"],
                quenching_model=quenching["model"],
            )
        else:
            res_ref = compute_dRdE(**common_kwargs)
    except Exception as e:
        return False, f"[wimp_nucleon-grid][ERROR] {prefix}mchi={mchi_str}: rate call failed: {e}"

    E_ref = res_ref["E_eV"]
    drde_ref = res_ref["dRdE_kg_year_eV"]
    meta_ref = res_ref["meta"]
    header_fn = _header_lines_ee if quenching else _header_lines

    errors = []
    for sigma_val, sigma_str, out_path in sigma_items:
        try:
            scale = float(sigma_val) / float(sigma_ref_val)
            drde = drde_ref * scale

            meta = dict(meta_ref)
            meta["sigma_n_cm2"] = float(sigma_val)

            CIO.write_csv_generic(str(out_path), E_ref, drde, header_fn(meta))
        except Exception as e:
            errors.append(
                f"[wimp_nucleon-grid][ERROR] {prefix}mchi={mchi_str} sigma={sigma_str}: {e}"
            )

    if errors:
        return False, "\n".join(errors)
    return True, ""


# ----------------------------------------------------------------------------
# run_from_config
#   Generate the WIMP-nucleon rate grid described by the config: for every mass and cross section compute the spin-independent rate (quenched to E_ee if a quenching model is given) and write the CSV.
# ----------------------------------------------------------------------------
def run_from_config(cfg: dict) -> None:
    target_nucleus = cfg.get("target_nucleus", "Si28")
    A = float(cfg.get("A", 28))
    mediator = cfg.get("mediator", "heavy")

    recoil = cfg.get("recoil", {})
    recoil = {
        "Emin_keV": float(recoil.get("Emin_keV", 0.001)),
        "Emax_keV": float(recoil.get("Emax_keV", 30.0)),
        "nbins": int(recoil.get("nbins", 3000)),
    }
    halo = cfg.get("halo", {})  # left as-is; compute_dRdE applies Baxter/Migdal defaults for any missing key

    quenching_cfg = cfg.get("quenching")
    quenching = None
    if quenching_cfg:
        quenching = {
            "model": quenching_cfg["model"],
            "ee_Emin_eV": float(quenching_cfg.get("ee_Emin_eV", 0.0)),
            "ee_Emax_eV": float(quenching_cfg.get("ee_Emax_eV", 8000.0)),
            "ee_nbins": int(quenching_cfg.get("ee_nbins", 800)),
        }

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
        mchi_str = format(mchi_MeV, fmt_mchi)
        sigma_items = []
        for sigma in sigma_list:
            sigma_str = format(sigma, fmt_sigma)
            out_path = _build_out_path(rates_dir, templ, target_nucleus, mediator, mchi_str, sigma_str)
            if (skip_existing or not overwrite) and _exists_any(out_path):
                continue
            sigma_items.append((sigma, sigma_str, out_path))
        if sigma_items:
            tasks.append({
                "target_nucleus": target_nucleus,
                "A": A,
                "mediator": mediator,
                "mchi_MeV": mchi_MeV,
                "mchi_str": mchi_str,
                "sigma_items": sigma_items,
                "recoil": recoil,
                "halo": halo,
                "quenching": quenching,
            })

    if not tasks:
        print("[wimp_nucleon-grid] all files already exist, nothing to do.")
        return

    n_tasks = len(tasks)
    total_sigma = sum(len(t["sigma_items"]) for t in tasks)
    print(
        f"[wimp_nucleon-grid] {n_tasks} mass points × up to {len(sigma_list)} sigma_n values "
        f"= {total_sigma} files to write  (1 rate eval per mass)"
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


# ----------------------------------------------------------------------------
# main
#   Command line: python3 utils/wimp_nucleon_generate_grid.py <config.json>; reads the config and calls run_from_config().
# ----------------------------------------------------------------------------
def main() -> None:
    if len(sys.argv) != 2:
        print("usage: python3 utils/wimp_nucleon_generate_grid.py <config.json>")
        sys.exit(1)
    cfg_path = Path(sys.argv[1]).expanduser()
    with open(cfg_path) as f:
        cfg = json.load(f)
    run_from_config(cfg)


if __name__ == "__main__":
    main()
