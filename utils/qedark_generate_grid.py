#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: qedark_generate_grid.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  qedark_generate_grid.py -- Multiprocess grid driver that expands (mχ, σe)
#  JSON grids and writes dR/dE CSVs via QEDark, QCDark, or QCDark2 backends.
# ============================================================================

"""
Expand (mχ, σ_e) grids into differential-rate CSV files.

JSON keys: ``material``, ``mediator``, ``halo``, ``detector``, ``rates_dir``,
``filename_template``, ``grid``, ``options``. Optional ``backend``: ``qedark``
(default), ``qcdark``, or ``qcdark2``. For QCDark, set ``form_factor_h5`` to
your crystal HDF5 (same layout as reference QCDark ``results/f2``). For QCDark2,
set ``epsilon_h5`` to a dielectric-function HDF5. In both cases, set
``detector.binsize_eV`` equal to the file ``dE``.

Convenience wrapper for QCDark defaults: ``utils/qcdark_generate_grid.py``.
"""
import os, sys, json, math, gzip, shutil
from pathlib import Path
from itertools import product
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, Dict, Any

# ensure local python package importable when run from repo root
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ccdarkphys.common.io import write_csv


def _dispatch_compute(backend: str):
    """Lazy import so multiprocessing workers only load the chosen backend."""
    b = (backend or "qedark").lower().strip()
    if b == "qcdark":
        from ccdarkphys.qcdark.entry import compute_dRdE

        return compute_dRdE, "QCDark"
    if b == "qcdark2":
        from ccdarkphys.qcdark2.entry import compute_dRdE

        return compute_dRdE, "QCDark2"
    if b == "qedark":
        from ccdarkphys.qedark.entry import compute_dRdE

        return compute_dRdE, "QEDark"
    raise ValueError(f"Unknown backend {backend!r}; use 'qedark', 'qcdark', or 'qcdark2'.")

# ---------------------------
# Grid expansion helpers
# ---------------------------

# ----------------------------------------------------------------------------
# _uniq_preserve
#   The items of xs without duplicates, in their original order.
# ----------------------------------------------------------------------------
def _uniq_preserve(xs):
    seen = set(); out = []
    for x in xs:
        if x in seen: continue
        seen.add(x); out.append(x)
    return out

def _expand_axis(spec: Any, kind: str) -> List[float]:
    """
    spec can be:
      - {"values":[...]}
      - {"linspace":{"start":..., "stop":..., "num":..., "endpoint":true}}
      - {"logspace":{"start_exp":..., "stop_exp":..., "num":..., "endpoint":true}}
      - OR a plain list (backward compatible)
    kind is only used for nicer error messages.
    Returns a list of floats.
    """
    import numpy as np

    if isinstance(spec, dict):
        vals = []
        if "values" in spec:
            vals += [float(v) for v in spec["values"]]
        if "linspace" in spec:
            p = spec["linspace"]
            start, stop = float(p["start"]), float(p["stop"])
            num = int(p["num"])
            endpoint = bool(p.get("endpoint", True))
            vals += list(np.linspace(start, stop, num=num, endpoint=endpoint, dtype=float))
        if "logspace" in spec:
            p = spec["logspace"]
            # Support two equivalent conventions:
            #   1) logspace.start_exp/stop_exp (log10 mass bounds)
            #   2) logspace.start/stop (literal mass bounds in MeV) → convert to log10 internally
            if "start_exp" in p or "stop_exp" in p:
                start_exp, stop_exp = float(p["start_exp"]), float(p["stop_exp"])
            else:
                start, stop = float(p["start"]), float(p["stop"])
                if start <= 0 or stop <= 0:
                    raise ValueError(f"logspace({kind}) requires positive start/stop when using start/stop mode")
                start_exp, stop_exp = math.log10(start), math.log10(stop)
            num = int(p["num"])
            endpoint = bool(p.get("endpoint", True))
            exps = np.linspace(start_exp, stop_exp, num=num, endpoint=endpoint, dtype=float)
            vals += list(np.power(10.0, exps))
        if not vals:
            raise ValueError(f"Grid spec for {kind!r} is empty or malformed: {spec}")
        return _uniq_preserve(vals)
    elif isinstance(spec, list):
        return [float(v) for v in spec]
    else:
        raise TypeError(f"Grid spec for {kind!r} must be dict or list, got {type(spec)}")

# ---------------------------
# Worker
# ---------------------------

# ----------------------------------------------------------------------------
# _build_out_path
#   Output path of one grid point from the file-name and sub-directory templates (with a .gz suffix when compressing); creates the directory.
# ----------------------------------------------------------------------------
def _build_out_path(base_dir: Path,
                    filename_template: str,
                    subdir_template: str,
                    material: str,
                    mediator: str,
                    mchi_MeV_str: str,
                    sigma_str: str,
                    compress: bool) -> Path:
    subdir = subdir_template.format(
        material=material, mediator=mediator,
        mchi_MeV=mchi_MeV_str, sigma_e_cm2=sigma_str
    ) if subdir_template else ""
    out_dir = (base_dir / subdir) if subdir else base_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    fname = filename_template.format(
        material=material, mediator=mediator,
        mchi_MeV=mchi_MeV_str, sigma_e_cm2=sigma_str
    )
    if compress and not fname.endswith(".gz"):
        fname += ".gz"
    return out_dir / fname

def _exists_any(out_path: Path) -> bool:
    """Check for either .csv or .csv.gz existing counterpart."""
    if out_path.exists():
        return True
    # Also consider the alternate extension (csv <-> csv.gz)
    if out_path.suffix == ".gz":
        alt = out_path.with_suffix("")  # drop .gz
        return alt.exists()
    else:
        gz = Path(str(out_path) + ".gz")
        return gz.exists()

# ----------------------------------------------------------------------------
# _write_csv_maybe_gz
#   Write the rate CSV, gzip-compressing it (through a temporary file) when compress is set.
# ----------------------------------------------------------------------------
def _write_csv_maybe_gz(
    out_path: Path, E, R, meta, compress: bool, *, entry: str = "QEDark"
):
    if not compress:
        write_csv(str(out_path), E, R, meta, entry=entry)
        return
    # write to tmp .csv then gzip
    tmp_csv = out_path.with_suffix(".tmp.csv")
    write_csv(str(tmp_csv), E, R, meta, entry=entry)
    with open(tmp_csv, "rb") as f_in, gzip.open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    tmp_csv.unlink(missing_ok=True)

def _one_task(task: Dict[str, Any]) -> Tuple[bool, str]:
    """
    Single (mass, sigma) task: compute one spectrum, write one CSV.
    Returns (ok, message).
    """
    try:
        compute_dRdE, _ = _dispatch_compute(task["backend"])
        kw = dict(
            material=task["material"],
            mediator=task["mediator"],
            mchi_eV=task["mchi_MeV"] * 1.0e6,
            sigma_e_cm2=task["sigma_val"],
            halo=task["halo"],
            band_gap_eV=task["detector"]["band_gap_eV"],
            eh_pair_eV=task["detector"]["eh_pair_eV"],
            binsize_eV=task["detector"]["binsize_eV"],
        )
        ff = task.get("form_factor_h5")
        if ff:
            kw["form_factor_h5"] = ff
        eps = task.get("epsilon_h5")
        if eps:
            kw["epsilon_h5"] = eps
        res = compute_dRdE(**kw)
        E, R, meta = res["E_eV"], res["dRdE_kg_year_eV"], res["meta"]
        _write_csv_maybe_gz(
            task["out_path"],
            E,
            R,
            meta,
            task["compress"],
            entry=task["csv_entry"],
        )
        return True, f"[grid] wrote {task['out_path']}"
    except Exception as e:
        return False, f"[grid][ERROR] {task['out_path']}: {e}"

def _one_task_mass(task: Dict[str, Any]) -> List[Tuple[bool, str]]:
    """
    One spectrum per mass, then scale and write one CSV per sigma (fast path).
    Rate is linear in sigma_e: R(sigma) = R(sigma_ref) * (sigma / sigma_ref).
    Returns list of (ok, message) for each file written.
    """
    sigma_ref = task["sigma_ref"]
    outputs = task["outputs"]  # list of (sigma_str, sigma_val, out_path)
    if not outputs:
        return []
    try:
        compute_dRdE, _ = _dispatch_compute(task["backend"])
        kw = dict(
            material=task["material"],
            mediator=task["mediator"],
            mchi_eV=task["mchi_MeV"] * 1.0e6,
            sigma_e_cm2=sigma_ref,
            halo=task["halo"],
            band_gap_eV=task["detector"]["band_gap_eV"],
            eh_pair_eV=task["detector"]["eh_pair_eV"],
            binsize_eV=task["detector"]["binsize_eV"],
        )
        ff = task.get("form_factor_h5")
        if ff:
            kw["form_factor_h5"] = ff
        eps = task.get("epsilon_h5")
        if eps:
            kw["epsilon_h5"] = eps
        res = compute_dRdE(**kw)
        E = res["E_eV"]
        R_ref = res["dRdE_kg_year_eV"]
        meta_base = res["meta"]
        compress = task["compress"]
        messages = []
        for sigma_str, sigma_val, out_path in outputs:
            R = R_ref * (float(sigma_val) / sigma_ref)
            meta = {**meta_base, "sigma_e_cm2": float(sigma_val)}
            _write_csv_maybe_gz(out_path, E, R, meta, compress, entry=task["csv_entry"])
            messages.append((True, f"[grid] wrote {out_path}"))
        return messages
    except Exception as e:
        return [(False, f"[grid][ERROR] m={task['mchi_MeV']}: {e}")]

# ---------------------------
# Main
# ---------------------------

# ----------------------------------------------------------------------------
# run_from_config
#   Generate the whole (mass, cross-section) rate grid described by the config for the chosen backend (qedark, qcdark or qcdark2): compute every point with the backend's entry point and write the CSVs.
# ----------------------------------------------------------------------------
def run_from_config(cfg: Dict[str, Any]) -> None:
    backend = str(cfg.get("backend", "qedark")).lower().strip()
    _, csv_entry = _dispatch_compute(backend)
    form_factor_h5 = cfg.get("form_factor_h5")  # optional; QCDark HDF5 path
    epsilon_h5 = cfg.get("epsilon_h5")          # optional; QCDark2 epsilon HDF5 path

    material = cfg["material"]                   # "Si"
    mediator = cfg["mediator"]                   # "heavy" | "massless" | ...
    halo     = cfg["halo"]                       # speeds in km/s or cm/s (entry handles both)
    det      = cfg.get("detector", {})
    det.setdefault("band_gap_eV", 1.2)
    det.setdefault("eh_pair_eV", 3.8)
    det.setdefault("binsize_eV", 0.1)

    rates_dir = Path(cfg["rates_dir"])
    templ     = cfg["filename_template"]         # "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv"
    subtempl  = cfg.get("subdir_template", "")   # e.g. "m={mchi_MeV}"

    opts      = cfg.get("options", {})
    skip_existing = bool(opts.get("skip_existing", True))
    overwrite     = bool(opts.get("overwrite", False))
    compress      = bool(opts.get("compress", False))
    parallel      = int(opts.get("parallel", 0))  # 0/1 => serial
    progress      = bool(opts.get("progress", True))
    # Fast path: one spectrum per mass, scale by sigma when writing (~N_sigma fewer compute_dRdE calls)
    use_sigma_scaling = bool(opts.get("use_sigma_scaling", True))
    sigma_ref         = float(opts.get("sigma_ref", 1e-36))
    fmt           = opts.get("format", {"mchi": ".6f", "sigma": ".1e"})
    fmt_mchi      = fmt.get("mchi", ".6f")
    fmt_sigma     = fmt.get("sigma", ".1e")

    print(f"[grid] backend={backend}  csv_tag={csv_entry}")

    # Expand grids
    gspec = cfg["grid"]
    mchi_list = _expand_axis(gspec["mchi_MeV"], "mchi_MeV")
    sigma_list_raw = _expand_axis(gspec["sigma_e_cm2"], "sigma_e_cm2")  # floats
    sigma_pairs = [(format(s, fmt_sigma), float(s)) for s in sigma_list_raw]

    if use_sigma_scaling:
        # One task per mass: compute once at sigma_ref, write one CSV per sigma (scaled)
        tasks = []
        for m in mchi_list:
            m_str = format(m, fmt_mchi)
            outputs = []
            for s_str, s_val in sigma_pairs:
                out_path = _build_out_path(
                    rates_dir, templ, subtempl, material, mediator, m_str, s_str, compress
                )
                if skip_existing and _exists_any(out_path):
                    if progress:
                        print(f"[grid] skip (exists): {out_path}")
                    continue
                if overwrite is False and _exists_any(out_path):
                    if progress:
                        print(f"[grid] skip (overwrite disabled): {out_path}")
                    continue
                outputs.append((s_str, s_val, out_path))
            if outputs:
                t = {
                    "backend": backend,
                    "csv_entry": csv_entry,
                    "material": material,
                    "mediator": mediator,
                    "halo": halo,
                    "detector": det,
                    "mchi_MeV": float(m),
                    "sigma_ref": sigma_ref,
                    "outputs": outputs,
                    "compress": compress,
                }
                if form_factor_h5:
                    t["form_factor_h5"] = form_factor_h5
                if epsilon_h5:
                    t["epsilon_h5"] = epsilon_h5
                tasks.append(t)
        if not tasks:
            print("[grid] nothing to do.")
            return
        n_spectra = len(tasks)
        n_files = sum(len(t["outputs"]) for t in tasks)
        if progress:
            print(f"[grid] fast path: {n_spectra} spectra → {n_files} files (sigma_ref={sigma_ref:.2e})")
        if parallel and parallel > 1:
            nproc = min(parallel, cpu_count())
            if progress:
                print(f"[grid] launching {n_spectra} mass jobs with {nproc} workers...")
            with Pool(processes=nproc) as pool:
                for messages in pool.imap_unordered(_one_task_mass, tasks):
                    for ok, msg in messages:
                        print(msg)
        else:
            if progress:
                print(f"[grid] running {n_spectra} mass jobs serially...")
            for i, t in enumerate(tasks, 1):
                for ok, msg in _one_task_mass(t):
                    if progress:
                        print(f"[{i}/{n_spectra}] {msg}")
    else:
        # Legacy: one task per (mass, sigma)
        tasks = []
        for m in mchi_list:
            m_str = format(m, fmt_mchi)
            for s_str, s_val in sigma_pairs:
                out_path = _build_out_path(
                    rates_dir, templ, subtempl, material, mediator, m_str, s_str, compress
                )
                if skip_existing and _exists_any(out_path):
                    if progress:
                        print(f"[grid] skip (exists): {out_path}")
                    continue
                if overwrite is False and _exists_any(out_path):
                    if progress:
                        print(f"[grid] skip (overwrite disabled): {out_path}")
                    continue
                tt = {
                    "backend": backend,
                    "csv_entry": csv_entry,
                    "material": material,
                    "mediator": mediator,
                    "halo": halo,
                    "detector": det,
                    "mchi_MeV": float(m),
                    "sigma_val": float(s_val),
                    "out_path": out_path,
                    "compress": compress,
                }
                if form_factor_h5:
                    tt["form_factor_h5"] = form_factor_h5
                if epsilon_h5:
                    tt["epsilon_h5"] = epsilon_h5
                tasks.append(tt)
        if not tasks:
            print("[grid] nothing to do.")
            return
        if parallel and parallel > 1:
            nproc = min(parallel, cpu_count())
            if progress:
                print(f"[grid] launching {len(tasks)} jobs with {nproc} workers...")
            with Pool(processes=nproc) as pool:
                for ok, msg in pool.imap_unordered(_one_task, tasks):
                    print(msg)
        else:
            if progress:
                print(f"[grid] running {len(tasks)} jobs serially...")
            for i, t in enumerate(tasks, 1):
                ok, msg = _one_task(t)
                if progress:
                    print(f"[{i}/{len(tasks)}] {msg}")


# ----------------------------------------------------------------------------
# main
#   Read the JSON config and run run_from_config().
# ----------------------------------------------------------------------------
def main(cfg_path: str) -> None:
    with open(cfg_path, "r") as f:
        cfg = json.load(f)
    run_from_config(cfg)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python3 utils/qedark_generate_grid.py <config.json>")
        sys.exit(1)
    main(sys.argv[1])
