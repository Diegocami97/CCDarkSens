#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — qcdark2_regenerate_epsilon
#  Regenerate QCDark2 dielectric-function HDF5 files with scissor band-gap control
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Regenerate QCDark2 dielectric-function HDF5 files (proper band-gap control via scissor).

Runs ``python -m qcdark2.dielectric_pyscf`` from a local QCDark2 install, then packages
``<name>_resources/epsilon.hdf5`` into a slim file usable by ``ccdarkphys.qcdark2``.

Band gap: set ``scissor_bandgap`` in the .in file (eV). DFT is cached under
``save_path/DFT_resources``; changing only scissor reuses DFT and recomputes ε(ω,q).

Examples:
  # Smoke test (4x4x4 k-grid, q_max=1)
  python3 utils/qcdark2_regenerate_epsilon.py --input configs/qcdark2/Si_smoke.in

  # Production LFE segment, gaps 1.1 and 1.2 eV
  python3 utils/qcdark2_regenerate_epsilon.py \\
    --template configs/qcdark2/Si_lfe8q.in --scissor 1.1 1.2

  # Package only (after a manual QCDark2 run)
  python3 utils/qcdark2_regenerate_epsilon.py --package-only \\
    --resources-dir data/qcdark2_epsilon/Si/Si_lfe8q_gap1.20_resources \\
    --out data/qcdark2_epsilon/Si/Si_gap1.20_lfe8q.h5
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
# Set CCDARK_QCDARK2_DIR to your local QCDark2 checkout, or pass --qcdark2 on the CLI.
_qcdark2_env = os.environ.get("CCDARK_QCDARK2_DIR", "")
DEFAULT_QCDARK2 = Path(_qcdark2_env) if _qcdark2_env else Path("QCDark2")
DEFAULT_VENV_PY = DEFAULT_QCDARK2 / ".venv" / "bin" / "python3"


def _resolve_python(qcdark2_root: Path, python_exe: str | None) -> Path:
    if python_exe:
        return Path(python_exe).expanduser().resolve()
    cand = qcdark2_root / ".venv" / "bin" / "python3"
    if cand.is_file():
        return cand
    return Path(sys.executable)


def _render_template(template_path: Path, scissor_eV: float) -> str:
    text = template_path.read_text()
    if "{SCISSOR}" not in text:
        raise ValueError(f"Template {template_path} has no {{SCISSOR}} placeholder.")
    gap = f"{scissor_eV:.4f}".rstrip("0").rstrip(".")
    return text.replace("{SCISSOR}", gap).replace("{SCISSOR_TAG}", gap.replace(".", "p"))


def _run_qcdark2(python: Path, qcdark2_root: Path, input_path: Path) -> None:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(qcdark2_root) + os.pathsep + env.get("PYTHONPATH", "")
    cmd = [str(python), "-m", "qcdark2.dielectric_pyscf", str(input_path)]
    print(f"[qcdark2] {' '.join(cmd)}")
    subprocess.run(cmd, cwd=str(qcdark2_root), env=env, check=True)


def _infer_name_from_input(inp: Path) -> str:
    for line in inp.read_text().splitlines():
        line = line.split("#", 1)[0].strip().replace(" ", "")
        if line.startswith("name="):
            return line.split("=", 1)[1]
    raise ValueError(f"Could not find name= in {inp}")


def _infer_save_path_from_input(inp: Path) -> Path:
    for line in inp.read_text().splitlines():
        line = line.split("#", 1)[0].strip().replace(" ", "")
        if line.startswith("save_path="):
            return Path(line.split("=", 1)[1])
    raise ValueError(f"Could not find save_path= in {inp}")


def package_epsilon(
    resources_dir: Path,
    out_path: Path,
    *,
    scissor_eV: float | None = None,
) -> None:
    """Copy epsilon.hdf5 to CCDarkSens layout (epsilon, q, E, M_cell, V_cell, dE)."""
    src = resources_dir / "epsilon.hdf5"
    if not src.is_file():
        raise FileNotFoundError(f"Missing {src}")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(src, "r") as fin, h5py.File(out_path, "w") as fout:
        for key in ("epsilon", "q", "E"):
            if key not in fin:
                raise KeyError(f"{src} missing dataset {key!r}")
            fout.create_dataset(key, data=fin[key][:], compression="gzip")
        for attr in ("M_cell", "V_cell", "dE"):
            if attr in fin.attrs:
                fout.attrs[attr] = fin.attrs[attr]
        if scissor_eV is not None:
            fout.attrs["scissor_bandgap_eV"] = float(scissor_eV)
        raw = fin.attrs.get("scissor_bandgap")
        if raw is not None and scissor_eV is None:
            try:
                fout.attrs["scissor_bandgap_eV"] = float(raw)
            except (TypeError, ValueError):
                pass

    with h5py.File(out_path, "r") as h:
        q, e = h["q"][:], h["E"][:]
    print(
        f"[package] wrote {out_path}  "
        f"q=[{float(q.min()):.4g}, {float(q.max()):.4g}]  "
        f"E=[{float(e.min()):.4g}, {float(e.max()):.4g}]"
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qcdark2-root", type=Path, default=DEFAULT_QCDARK2)
    ap.add_argument("--python", default=None, help="Python with qcdark2+pyscf (default: QCDark2/.venv)")
    ap.add_argument("--input", type=Path, help="Ready-to-run .in file")
    ap.add_argument("--template", type=Path, help=".in template with {SCISSOR} placeholders")
    ap.add_argument("--scissor", type=float, nargs="+", metavar="eV", help="Scissor band gaps (eV)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--package-only", action="store_true")
    ap.add_argument("--resources-dir", type=Path)
    ap.add_argument("--out", type=Path, help="Output .h5 for --package-only")
    args = ap.parse_args()

    qcdark2_root = args.qcdark2_root.expanduser().resolve()
    python = _resolve_python(qcdark2_root, args.python)

    if args.package_only:
        if not args.resources_dir or not args.out:
            ap.error("--package-only requires --resources-dir and --out")
        package_epsilon(args.resources_dir.resolve(), args.out.resolve())
        return

    if args.input:
        inputs = [(args.input.resolve(), None)]
    elif args.template and args.scissor:
        inputs = []
        for g in args.scissor:
            inputs.append((args.template.resolve(), float(g)))
    else:
        ap.error("Provide --input FILE, or --template FILE --scissor eV [eV ...]")

    for template_or_inp, scissor in inputs:
        if scissor is not None:
            body = _render_template(template_or_inp, scissor)
            tag = f"{scissor:.4f}".rstrip("0").rstrip(".").replace(".", "p")
            with tempfile.NamedTemporaryFile(
                mode="w",
                suffix=".in",
                prefix=f"qcdark2_gap{tag}_",
                delete=False,
            ) as tf:
                tf.write(body)
                inp_path = Path(tf.name)
        else:
            inp_path = template_or_inp

        name = _infer_name_from_input(inp_path)
        save_path = _infer_save_path_from_input(inp_path)
        save_path.mkdir(parents=True, exist_ok=True)
        resources = save_path / f"{name}_resources"

        if args.dry_run:
            print(f"[dry-run] would run with {inp_path} -> {resources}/epsilon.hdf5")
            continue

        try:
            _run_qcdark2(python, qcdark2_root, inp_path)
        finally:
            if scissor is not None:
                inp_path.unlink(missing_ok=True)

        gap = scissor
        if gap is None:
            m = re.search(r"gap([\dp]+)", name)
            if m:
                gap = float(m.group(1).replace("p", "."))

        out_name = name
        if gap is not None and "gap" not in out_name:
            out_name = f"{name}_scissor{gap:.2f}eV"
        out_h5 = save_path / f"{out_name}.h5"
        package_epsilon(resources, out_h5, scissor_eV=gap)


if __name__ == "__main__":
    main()
