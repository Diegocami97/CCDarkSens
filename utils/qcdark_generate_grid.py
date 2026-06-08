#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — qcdark_generate_grid
#  CLI wrapper that sets backend=qcdark and invokes qedark_generate_grid.py with resolved crystal form-factor HDF5 paths.
#
#  Author: Diego Venegas-Vargas
# ============================================================================

"""
Generate QCDark differential-rate CSV grids — same JSON schema and driver as
``qedark_generate_grid.py``, but defaults ``backend`` to ``qcdark``.

Typical workflow:
  1. Produce or point to a crystal |F|^2 HDF5 (standalone QCDark workflow), **or**
     run ``python -m ccdarkphys.qcdark.demo_run`` once to create a synthetic demo table.
  2. Set ``form_factor_h5`` in the JSON (absolute path, or repo-relative when run from repo root).
  3. Set ``detector.binsize_eV`` equal to ``dE`` inside that HDF5.

Usage (from repository root)::

    python3 utils/qcdark_generate_grid.py configs/qcdark_generate_demo.json

Equivalent to ``qedark_generate_grid.py`` with ``\"backend\": \"qcdark\"`` in the config.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_driver():
    path = Path(__file__).resolve().parent / "qedark_generate_grid.py"
    spec = importlib.util.spec_from_file_location("qedark_generate_grid", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: python3 utils/qcdark_generate_grid.py <config.json>")
        sys.exit(1)

    cfg_path = Path(sys.argv[1]).expanduser()
    with open(cfg_path, "r") as f:
        cfg = json.load(f)

    cfg.setdefault("backend", "qcdark")

    # Resolve optional HDF5 path: relative paths are from cwd (usually repo root).
    ff = cfg.get("form_factor_h5") or os.environ.get("CCDARK_SENS_QCDARK_FORM_FACTOR")
    if ff:
        p = Path(ff).expanduser()
        if not p.is_absolute():
            p = Path.cwd() / p
        cfg["form_factor_h5"] = str(p.resolve())

    sys.path.insert(0, str(_REPO_ROOT / "python"))
    driver = _load_driver()
    driver.run_from_config(cfg)


if __name__ == "__main__":
    main()
