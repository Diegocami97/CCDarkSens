"""
Small I/O helpers for data discovery and CSV writing.
"""

from __future__ import annotations
import os
import hashlib

def data_path(pkg_file: str, subdir: str, filename: str) -> str:
    """
    Resolve a data file shipped inside the package.
    pkg_file: usually __file__ of the calling module
    """
    base = os.path.dirname(os.path.abspath(pkg_file))
    cand = os.path.join(base, subdir, filename)
    if not os.path.isfile(cand):
        raise FileNotFoundError(f"Data file not found: {cand}")
    return cand

def sha1sum(path: str) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

def _csv_header_lines(entry: str) -> list:
    return [
        f"# Differential Rates computed with CCDarkSens ({entry} entry)",
        "# material = {material}, mediator = {mediator}, table = {table_path}",
        "# table_sha1 = {table_sha1}",
        "# halo (cm/s): v0={v0_cm_s}, vE={vE_cm_s}, vesc={vesc_cm_s}",
        "# mX (eV) = {mchi_eV}",
        "# sigma_e (cm^2) = {sigma_e_cm2}",
        "# Output units: dR/dE in events / kg / year / eV",
        "# Columns: E (eV), dRdE (events/kg/year/eV)",
    ]

# Kept for callers that expect a static name; first line is the same as _csv_header_lines("QEDark")[0]
CSV_HEADER = _csv_header_lines("QEDark")

def write_csv(out_path: str, E, R, meta: dict, *, entry: str = "QEDark") -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        for line in _csv_header_lines(entry):
            f.write(line.format(**meta) + "\n")
        f.write("E,dRdE\n")
        for e, v in zip(E, R):
            f.write(f"{e:.8g},{v:.10g}\n")
