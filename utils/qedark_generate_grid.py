#!/usr/bin/env python3
import os, sys, json, math, gzip, shutil
from pathlib import Path
from itertools import product
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, Dict, Any

# ensure local python package importable when run from repo root
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

from ccdarkphys.qedark.entry import compute_dRdE
from ccdarkphys.common.io import write_csv

# ---------------------------
# Grid expansion helpers
# ---------------------------

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
            start_exp, stop_exp = float(p["start_exp"]), float(p["stop_exp"])
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

def _write_csv_maybe_gz(out_path: Path, E, R, meta, compress: bool):
    if not compress:
        write_csv(str(out_path), E, R, meta)
        return
    # write to tmp .csv then gzip
    tmp_csv = out_path.with_suffix(".tmp.csv")
    write_csv(str(tmp_csv), E, R, meta)
    with open(tmp_csv, "rb") as f_in, gzip.open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    tmp_csv.unlink(missing_ok=True)

def _one_task(task: Dict[str, Any]) -> Tuple[bool, str]:
    """
    Returns (ok, message)
    """
    try:
        res = compute_dRdE(
            material=task["material"],
            mediator=task["mediator"],
            mchi_eV=task["mchi_MeV"] * 1.0e6,
            sigma_e_cm2=task["sigma_val"],
            halo=task["halo"],
            band_gap_eV=task["detector"]["band_gap_eV"],
            eh_pair_eV=task["detector"]["eh_pair_eV"],
            binsize_eV=task["detector"]["binsize_eV"],
        )
        E, R, meta = res["E_eV"], res["dRdE_kg_year_eV"], res["meta"]
        _write_csv_maybe_gz(task["out_path"], E, R, meta, task["compress"])
        return True, f"[grid] wrote {task['out_path']}"
    except Exception as e:
        return False, f"[grid][ERROR] {task['out_path']}: {e}"

# ---------------------------
# Main
# ---------------------------

def main(cfg_path: str):
    with open(cfg_path, "r") as f:
        cfg = json.load(f)

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
    fmt           = opts.get("format", {"mchi": ".6f", "sigma": ".1e"})
    fmt_mchi      = fmt.get("mchi", ".6f")
    fmt_sigma     = fmt.get("sigma", ".1e")

    # Expand grids
    gspec = cfg["grid"]
    mchi_list = _expand_axis(gspec["mchi_MeV"], "mchi_MeV")
    sigma_list_raw = _expand_axis(gspec["sigma_e_cm2"], "sigma_e_cm2")  # floats
    # Preserve both pretty string & float value
    sigma_pairs = [(format(s, fmt_sigma), float(s)) for s in sigma_list_raw]

    # Build task list
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
            tasks.append({
                "material": material,
                "mediator": mediator,
                "halo": halo,
                "detector": det,
                "mchi_MeV": float(m),
                "sigma_val": float(s_val),
                "out_path": out_path,
                "compress": compress,
            })

    if not tasks:
        print("[grid] nothing to do.")
        return

    # Run
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

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python3 utils/qedark_generate_grid.py <config.json>")
        sys.exit(1)
    main(sys.argv[1])
