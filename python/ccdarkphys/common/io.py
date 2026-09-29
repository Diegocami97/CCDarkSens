# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  io.py -- Diego Venegas-Vargas DAMIC-M collaboration CCDarkSens Framework
#  io.py -- Small I/O helpers for data discovery and CSV writing.
# ============================================================================

"""
Small I/O helpers for data discovery and CSV writing.
"""

from __future__ import annotations
import os
import hashlib

# ----------------------------------------------------------------------------
# data_path
#   Full path of a data file shipped inside the package; raises FileNotFoundError if it does not exist.
# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# sha1sum
#   SHA-1 hex digest of a file (read in 64 kB chunks); I record it in the CSV headers so a rate table can be traced to its input table.
# ----------------------------------------------------------------------------
def sha1sum(path: str) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()

# ----------------------------------------------------------------------------
# _csv_header_lines
#   Template of the comment header of a DM-electron rate CSV (material, mediator, table hash, halo, mass, cross section, units).
# ----------------------------------------------------------------------------
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
CSV_HEADER = _csv_header_lines("QEDark")  # static copy of the QEDark header for callers that expect this name

# ----------------------------------------------------------------------------
# write_csv
#   Write a DM-electron rate table: the header template filled from meta, then the columns E [eV] and dRdE [events/kg/year/eV].
# ----------------------------------------------------------------------------
def write_csv(out_path: str, E, R, meta: dict, *, entry: str = "QEDark") -> None:
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        for line in _csv_header_lines(entry):
            f.write(line.format(**meta) + "\n")
        f.write("E,dRdE\n")
        for e, v in zip(E, R):
            f.write(f"{e:.8g},{v:.10g}\n")

def write_csv_generic(out_path: str, E, R, header_lines: list) -> None:
    """Same E,dRdE body format as write_csv, but with caller-supplied header
    text instead of the DM-electron-specific template (halo velocities,
    sigma_e) — for backends (dark photon, Migdal) whose metadata doesn't fit
    that template."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        for line in header_lines:
            f.write(line + "\n")
        f.write("E,dRdE\n")
        for e, v in zip(E, R):
            f.write(f"{e:.8g},{v:.10g}\n")
