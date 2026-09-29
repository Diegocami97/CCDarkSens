#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: darkphoton_generate_grid.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  darkphoton_generate_grid.py -- Grid driver that expands a (mA', epsilon)
#  JSON grid and writes dR/dE CSV rate tables for hidden-photon absorption
#  via darkelf.
# ============================================================================

"""
Expand a (mA', epsilon) grid into dark-photon-absorption rate CSV files.

NOT built on the qedark_generate_grid.py shared driver: that driver's fast
path rescales linearly in sigma, which is wrong for absorption (R ∝ ε²).

Epsilon fast path (ε² rescaling)
---------------------------------
For fixed m_A', the ELF is fixed, and darkelf gives

    R(ε) = ε² × foo × ELF(m_A', m_A')

so R(ε) = R(ε_ref) × (ε/ε_ref)².  This is exact — ε only enters as an
overall prefactor.  The script therefore calls darkelf once per m_A'
(at the first ε in the list that needs a file), then rescales all remaining
ε values analytically.  For an 80×30 grid this reduces DarkELF inits from
2 400 → 80.

The "one init per sigma" linear shortcut used by qedark_generate_grid is
still forbidden: R ∝ σ (DM-electron) vs R ∝ ε² (dark-photon) differ in
the power, so the two strategies are not interchangeable.

JSON keys
---------
material, darkelf_dir (or CCDARK_SENS_DARKELF_DIR), rates_dir,
filename_template, grid ({"mA_eV": ..., "epsilon": ...}),
detector ({"band_gap_eV": null, "binsize_eV": 0.1}), options.

Usage::

    python3 utils/darkphoton_generate_grid.py configs/darkphoton_generate_si_demo.json
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
from ccdarkphys.darkphoton.entry import compute_dRdE, _header_lines


# ----------------------------------------------------------------------------
# _build_out_path
#   Output file path for one (mass, epsilon) point from the filename template; creates the directory.
# ----------------------------------------------------------------------------
def _build_out_path(base_dir: Path, filename_template: str, material: str,
                    mA_str: str, eps_str: str) -> Path:
    fname = filename_template.format(material=material, mA_eV=mA_str, epsilon=eps_str)
    base_dir.mkdir(parents=True, exist_ok=True)
    return base_dir / fname


# ----------------------------------------------------------------------------
# _exists_any
#   True if the output file exists, either plain or gzip-compressed.
# ----------------------------------------------------------------------------
def _exists_any(out_path: Path) -> bool:
    return out_path.exists() or Path(str(out_path) + ".gz").exists()


def _boxcar(mA_eV: float, rate_kg_yr: float, binsize_eV: float):
    """Build the 4-point boxcar representing a monochromatic line at mA_eV."""
    half = 0.5 * binsize_eV
    edge = 1e-3 * binsize_eV
    height = rate_kg_yr / binsize_eV
    E = np.array([mA_eV - half - edge, mA_eV - half, mA_eV + half, mA_eV + half + edge])
    dRdE = np.array([0.0, height, height, 0.0])
    return E, dRdE


def _one_mass_task(task: dict) -> tuple[bool, str]:
    """
    Process all epsilon values for one m_A' point.

    Calls darkelf once at eps_list[0] (the reference), then rescales
    remaining epsilon values with (eps/eps_ref)^2 — exact because R ∝ ε².
    """
    material = task["material"]
    mA_eV = task["mA_eV"]
    mA_str = task["mA_str"]
    eps_items = task["eps_items"]   # list of (eps_val, eps_str, out_path)
    band_gap_eV = task["band_gap_eV"]
    binsize_eV = task["binsize_eV"]
    darkelf_dir = task["darkelf_dir"]
    idx = task.get("index", 0)
    n_tasks = task.get("n_tasks", 0)
    prefix = f"[{idx}/{n_tasks}] " if idx and n_tasks else ""

    errors = []
    eps_csv = task.get("eps_csv")
    density_g_cm3 = task.get("density_g_cm3")
    eps_ref_val, eps_ref_str, _ = eps_items[0]

    print(
        f"[darkphoton-grid] {prefix}mA'={mA_str} eV — "
        f"init at eps_ref={eps_ref_str}, then {len(eps_items)} ε² rescales",
        flush=True,
    )

    # One darkelf init + R_absorption call at the reference epsilon
    try:
        with open(os.devnull, "w") as devnull:
            with contextlib.redirect_stdout(devnull), contextlib.redirect_stderr(devnull):
                res_ref = compute_dRdE(
                    material=material,
                    mA_eV=mA_eV,
                    epsilon=eps_ref_val,
                    band_gap_eV=band_gap_eV,
                    binsize_eV=binsize_eV,
                    eps_csv=eps_csv,
                    density_g_cm3=density_g_cm3,
                    darkelf_dir=darkelf_dir,
                )
    except Exception as e:
        return False, f"[darkphoton-grid][ERROR] {prefix}mA'={mA_str}: darkelf call failed: {e}"

    rate_ref = res_ref["meta"]["rate_total_kg_yr"]
    meta_ref = res_ref["meta"]

    # Write all epsilon CSVs using ε² rescaling
    for eps_val, eps_str, out_path in eps_items:
        try:
            scale = (eps_val / eps_ref_val) ** 2
            rate = rate_ref * scale
            E, dRdE = _boxcar(mA_eV, rate, binsize_eV)

            meta = dict(meta_ref)
            meta["epsilon"] = float(eps_val)
            meta["rate_total_kg_yr"] = rate

            CIO.write_csv_generic(str(out_path), E, dRdE, _header_lines(meta))
        except Exception as e:
            errors.append(
                f"[darkphoton-grid][ERROR] {prefix}mA'={mA_str} eps={eps_str}: {e}"
            )

    if errors:
        return False, "\n".join(errors)
    return True, ""


# ----------------------------------------------------------------------------
# run_from_config
#   Generate the dark-photon absorption rate grid described by the config: for every mass and epsilon compute the DarkELF rate and write the CSV, skipping existing files unless overwrite is requested.
# ----------------------------------------------------------------------------
def run_from_config(cfg: dict) -> None:
    material = cfg["material"]
    darkelf_dir = cfg.get("darkelf_dir") or os.environ.get("CCDARK_SENS_DARKELF_DIR")
    if not darkelf_dir:
        raise FileNotFoundError(
            "darkelf_dir is required (set in config or CCDARK_SENS_DARKELF_DIR env var)."
        )

    det = cfg.get("detector", {})
    band_gap_eV = det.get("band_gap_eV")
    binsize_eV = det.get("binsize_eV", 0.1)
    eps_csv = cfg.get("eps_csv")
    density_g_cm3 = cfg.get("density_g_cm3")

    rates_dir = Path(cfg["rates_dir"])
    templ = cfg["filename_template"]

    opts = cfg.get("options", {})
    skip_existing = bool(opts.get("skip_existing", True))
    overwrite = bool(opts.get("overwrite", False))
    parallel = int(opts.get("parallel", 0))
    fmt = opts.get("format", {"mA": ".6f", "epsilon": ".3e"})
    fmt_mA = fmt.get("mA", ".6f")
    fmt_eps = fmt.get("epsilon", ".3e")

    gspec = cfg["grid"]
    mA_list = _expand_axis(gspec["mA_eV"], "mA_eV")
    eps_list = _expand_axis(gspec["epsilon"], "epsilon")

    # Build per-mass tasks; each task carries all epsilon values that need files
    tasks = []
    for mA in mA_list:
        mA_str = format(mA, fmt_mA)
        eps_items = []
        for eps in eps_list:
            eps_str = format(eps, fmt_eps)
            out_path = _build_out_path(rates_dir, templ, material, mA_str, eps_str)
            if (skip_existing or not overwrite) and _exists_any(out_path):
                continue
            eps_items.append((eps, eps_str, out_path))
        if eps_items:
            tasks.append({
                "material": material,
                "mA_eV": mA,
                "mA_str": mA_str,
                "eps_items": eps_items,
                "band_gap_eV": band_gap_eV,
                "binsize_eV": binsize_eV,
                "darkelf_dir": darkelf_dir,
                "eps_csv": eps_csv,
                "density_g_cm3": density_g_cm3,
            })

    if not tasks:
        print("[darkphoton-grid] all files already exist, nothing to do.")
        return

    n_tasks = len(tasks)
    total_eps = sum(len(t["eps_items"]) for t in tasks)
    print(
        f"[darkphoton-grid] {n_tasks} mass points × up to {len(eps_list)} ε values "
        f"= {total_eps} files to write  (1 darkelf init per mass)"
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
#   Command line: python3 utils/darkphoton_generate_grid.py <config.json>.
# ----------------------------------------------------------------------------
def main() -> None:
    if len(sys.argv) != 2:
        print("usage: python3 utils/darkphoton_generate_grid.py <config.json>")
        sys.exit(1)
    cfg_path = Path(sys.argv[1]).expanduser()
    with open(cfg_path) as f:
        cfg = json.load(f)
    run_from_config(cfg)


if __name__ == "__main__":
    main()
