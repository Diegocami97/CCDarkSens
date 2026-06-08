"""
Runnable QCDark demo (standalone CCDarkSens; no collab_frameworks).

Creates a **tiny synthetic** crystal HDF5 (valid layout for ``CrystalFormFactor``),
runs ``compute_dRdE`` once, and optionally writes a CSV next to the fixture.

Requires: numpy, scipy, h5py

Writes the synthetic fixture under ``<repo>/data/qcdark/demo_Si_f2_qcdark.h5``.

Usage from the repository root::

    PYTHONPATH=python python3 -m ccdarkphys.qcdark.demo_run

Or with a custom output CSV::

    PYTHONPATH=python python3 -m ccdarkphys.qcdark.demo_run --out_csv data/qcdark_rates/demo/demo_qcdark.csv
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np

from ccdarkphys.common import constants as QEC
from ccdarkphys.common import io as CIO
from ccdarkphys.qcdark.entry import compute_dRdE, repo_qcdark_data_dir


def _demo_h5_path() -> Path:
    return repo_qcdark_data_dir() / "demo_Si_f2_qcdark.h5"


def write_demo_fixture(path: Path | str, *, seed: int = 0) -> None:
    """
    Write a minimal HDF5 file compatible with ``CrystalFormFactor``.

    The |F|^2 values are **not** physical; they only exercise the pipeline.
    """
    try:
        import h5py
    except ImportError as exc:
        raise ImportError("Demo needs h5py: pip install h5py") from exc

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    nq, n_e = 31, 40
    d_e = 0.1
    dq_attr = 0.02  # dimensionless; multiplied by α m_e in the loader

    rng = np.random.default_rng(seed)
    ff = np.abs(rng.standard_normal((nq, n_e))) * 1e-30 + 1e-40

    m_cell_si = 2.0 * 28.0855 * QEC.amu_kg

    with h5py.File(path, "w") as f:
        rs = f.create_group("run_settings")
        rs.create_dataset("a", data=np.eye(3, dtype=np.float64) * 5.431)
        rs.attrs["atom"] = "Si"
        rs.attrs["basis"] = "TZP"
        rs.attrs["ecp"] = "None"
        rs.attrs["numcon"] = "all"
        rs.attrs["numval"] = "all"
        rs.attrs["xc"] = "pbe"
        rs.attrs["df"] = "MDF"
        rs.create_dataset("rcut", data=np.array(4.0, dtype=np.float64))
        rs.create_dataset("precision", data=np.array(1e-8, dtype=np.float64))
        rs.create_dataset("kpts", data=np.zeros((1, 3), dtype=np.float64))

        res = f.create_group("results")
        res.attrs["VCell"] = np.float64(160.12)
        res.attrs["mCell"] = np.float64(m_cell_si)
        res.attrs["dq"] = np.float64(dq_attr)
        res.attrs["dE"] = np.float64(d_e)
        res.attrs["bandgap"] = np.float64(1.11)
        res.attrs["scissor"] = np.bool_(True)
        res.create_dataset("f2", data=ff.astype(np.float64))


def main() -> None:
    ap = argparse.ArgumentParser(description="CCDarkSens QCDark demo (synthetic HDF5 + one rate).")
    ap.add_argument(
        "--out_csv",
        default="",
        help="If set, write this CSV (default: print summary only).",
    )
    ap.add_argument(
        "--fixture",
        default="",
        help="HDF5 path (default: data/qcdark/demo_Si_f2_qcdark.h5 under repo root).",
    )
    args = ap.parse_args()

    fixture = Path(args.fixture) if args.fixture else _demo_h5_path()
    if not fixture.is_file():
        print(f"[demo] writing synthetic fixture → {fixture}")
        write_demo_fixture(fixture)
    else:
        print(f"[demo] using existing fixture {fixture}")

    halo = {"v0_kms": 238.0, "vE_kms": 263.0, "vesc_kms": 544.0}
    res = compute_dRdE(
        material="Si",
        mediator="heavy",
        mchi_eV=10e6,
        sigma_e_cm2=1e-37,
        halo=halo,
        band_gap_eV=1.2,
        eh_pair_eV=3.8,
        binsize_eV=0.1,
        form_factor_h5=str(fixture),
    )
    e = res["E_eV"]
    r = res["dRdE_kg_year_eV"]
    meta = res["meta"]

    print("[demo] OK — same API as qedark.entry.compute_dRdE")
    print(f"       N_bins={len(e)}, E[0..2]={e[:3]}, max(dRdE)={float(np.max(r)):.6e}")
    print(f"       table_sha1={meta['table_sha1'][:12]}…")

    if args.out_csv:
        CIO.write_csv(args.out_csv, e, r, meta, entry="QCDark")
        print(f"[demo] wrote {args.out_csv}")

    print()
    print("Same grid driver as QEDark; QCDark convenience entry point:")
    print(f"  python3 utils/qcdark_generate_grid.py configs/qcdark_generate_demo.json")
    print("(Run this once before grids if the demo HDF5 is missing:)")
    print(f"  PYTHONPATH=python python3 -m ccdarkphys.qcdark.demo_run")


if __name__ == "__main__":
    main()
