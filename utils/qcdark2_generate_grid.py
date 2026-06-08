#!/usr/bin/env python3
# ============================================================================
#  CCDarkSens — qcdark2_generate_grid
#  CLI wrapper that sets backend=qcdark2 and calls the shared grid driver to emit QCDark2 dR/dE CSV rate tables from JSON.
#
#  Author: Diego Venegas-Vargas
# ============================================================================

"""
Generate QCDark2 differential-rate CSV grids using the same JSON schema as
``qedark_generate_grid.py``, but defaulting ``backend`` to ``qcdark2``.

Usage (from repository root)::

    python3 utils/qcdark2_generate_grid.py configs/qcdark2_generate_si_comp_demo.json
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
        print("usage: python3 utils/qcdark2_generate_grid.py <config.json>")
        sys.exit(1)

    cfg_path = Path(sys.argv[1]).expanduser()
    with open(cfg_path, "r") as f:
        cfg = json.load(f)

    cfg.setdefault("backend", "qcdark2")

    eps = cfg.get("epsilon_h5") or os.environ.get("CCDARK_SENS_QCDARK2_EPSILON")
    if eps:
        p = Path(eps).expanduser()
        if not p.is_absolute():
            p = Path.cwd() / p
        cfg["epsilon_h5"] = str(p.resolve())

    sys.path.insert(0, str(_REPO_ROOT / "python"))
    driver = _load_driver()
    driver.run_from_config(cfg)


if __name__ == "__main__":
    main()
