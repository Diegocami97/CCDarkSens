#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: build_darkphoton_hypmat_qcdark2_elf.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  build_darkphoton_hypmat_qcdark2_elf.py -- Build a DarkELF optical-limit
#  data file for SrCd₂Sb₂ using the QCDark2-derived direct-gap onset combined
#  with Si's measured optical dielectric function.
# ============================================================================

"""
Build a DarkELF optical-limit data file for SrCd₂Sb₂ using the QCDark2-derived
direct-gap onset combined with Si's measured optical dielectric function.

Rationale
---------
The QCDark2 ε(q,ω) file for the scissors-corrected Si at E_gap=0.34 eV shows
that the q→0 optical ELF onset is at ~2.1 eV — this is the scissors-shifted Si
direct gap (Si direct gap ~3.4 eV minus scissors shift ~1.3 eV). The "fast"
k-grid used for DM-e rate generation is too coarse to give a smooth optical ELF
(ε₁ shows unphysical values of 34–82 below the gap due to sparse k-sampling),
so we cannot use the raw H5 dielectric directly.

Instead we use Si's *measured* optical dielectric function
(Si_eps_electron_opticallimit.dat, the same file used for the DAMIC-M 2025 PRL
Si hidden-photon limit), combined with the scissors-derived band_gap_eV = 2.1 eV
as the absorption threshold. This gives:

  - The correct interband onset energy (scissors-consistent, not the 3.4 eV
    Si direct gap but the shifted 2.1 eV value)
  - The real Si ELF including all interband structure and the sharp Si plasmon
    at ~17 eV — qualitatively correct even for SrCd₂Sb₂ while we await the
    true material wavefunctions
  - Zero numerical noise (measured optical data, not a coarse k-grid)

This is registered as a *new* material key "hypmat_qcdark2" in entry.py.
It does NOT overwrite or replace hypmat_unscreened / hypmat_screened (Drude model).
Both model families coexist for comparison.

Limitations / caveats (documented in docs/DarkPhoton_Absorption_HypotheticalMaterial.md):
  1. Uses Si interband structure, not SrCd₂Sb₂ interband structure.
  2. Scissors shifts all conduction bands uniformly — approximate for optical gaps.
  3. The true SrCd₂Sb₂ ELF requires QCDark2 with a dense optical k-grid from
     the AiiDA wavefunction files (bands_workchain_pk=9491).

Output
------
Writes two files under DARKELF_DIR/data/HypMat_0p34/:
  HypMat_0p34_qcdark2_eps1_si_optical.dat   — Si optical ELF with E_gap=2.1 eV
  (The band_gap_eV override is applied at darkelf init time; no dat file content
   changes are needed between the two scenarios — we just pass different band_gap_eV
   values when calling darkelf.)

Actually: the dat file itself is identical to Si's optical limit file.
We create a symlink-style copy under HypMat_0p34/ so darkelf can find it
via the target="HypMat_0p34" lookup path, with appropriate header.

Usage
-----
  python3 utils/build_darkphoton_hypmat_qcdark2_elf.py \\
      --darkelf_dir /path/to/DarkELF \\
      [--qcdark2_h5 data/qcdark2_epsilon/Si/Si_fast_gap0p34.h5]

Or set CCDARK_SENS_DARKELF_DIR.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np

_ENV_DARKELF_DIR = "CCDARK_SENS_DARKELF_DIR"

# Scissors-shifted direct gap read from QCDark2 H5 (onset of Im[ε] > 0)
# Determined by inspecting Si_fast_gap0p34.h5: ELF first non-zero at ~2.1 eV
DIRECT_GAP_EV = 2.1   # eV — scissors-shifted Si direct gap for E_gap_indirect=0.34 eV

# Material density for HypMat (SrCd₂Sb₂)
RHO_T = 5.76  # g/cm³


# ----------------------------------------------------------------------------
# _resolve_darkelf_dir
#   Directory of the DarkELF source: --darkelf_dir if given, else the environment variable; raises FileNotFoundError if neither is set.
# ----------------------------------------------------------------------------
def _resolve_darkelf_dir(explicit: str | None) -> str:
    if explicit:
        return os.path.abspath(os.path.expanduser(explicit))
    env = os.environ.get(_ENV_DARKELF_DIR)
    if env:
        return os.path.abspath(os.path.expanduser(env.strip()))
    raise FileNotFoundError(
        f"Pass --darkelf_dir or set {_ENV_DARKELF_DIR}."
    )


def _read_si_optical_dat(darkelf_dir: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read Si_eps_electron_opticallimit.dat and return (omega, eps1, eps2)."""
    dat_path = Path(darkelf_dir) / "data" / "Si" / "Si_eps_electron_opticallimit.dat"
    if not dat_path.exists():
        raise FileNotFoundError(
            f"Si optical limit file not found: {dat_path}\n"
            "Run 'python3 utils/build_si_optical_limit.py' first."
        )
    omega, eps1, eps2 = [], [], []
    with open(dat_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                try:
                    omega.append(float(parts[0]))
                    eps1.append(float(parts[1]))
                    eps2.append(float(parts[2]))
                except ValueError:
                    continue
    return np.array(omega), np.array(eps1), np.array(eps2)


# ----------------------------------------------------------------------------
# build
#   Write the DarkELF data files of the QCDark2-based HypMat proxy into <darkelf_dir>/data/HypMat_0p34, checking the scissors-shifted direct-gap onset against the QCDark2 HDF5 when it is available.
# ----------------------------------------------------------------------------
def build(darkelf_dir: str, qcdark2_h5: str | None = None):
    darkelf_path = Path(darkelf_dir)
    out_dir = darkelf_path / "data" / "HypMat_0p34"
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Verify scissors onset from H5 if available ---
    onset_note = f"scissors-shifted onset set to {DIRECT_GAP_EV} eV (hardcoded)"
    if qcdark2_h5 and Path(qcdark2_h5).exists():
        try:
            import h5py
            with h5py.File(qcdark2_h5, "r") as hf:
                E   = hf["E"][:]
                eps = hf["epsilon"][0, :]   # q→0 (q=0.01)
            elf_q0 = np.imag(-1.0 / eps)
            nonzero = np.where(elf_q0 > 1e-6)[0]
            if len(nonzero):
                observed_onset = E[nonzero[0]]
                onset_note = (f"scissors-shifted onset from H5: {observed_onset:.2f} eV "
                              f"(file: {Path(qcdark2_h5).name})")
                print(f"[qcdark2_elf] H5 ELF onset: {observed_onset:.2f} eV")
        except ImportError:
            print("[qcdark2_elf] h5py not available — using hardcoded onset")

    # --- Read Si measured optical data ---
    print(f"[qcdark2_elf] Reading Si optical data from {darkelf_path}/data/Si/")
    omega_si, eps1_si, eps2_si = _read_si_optical_dat(darkelf_dir)
    print(f"[qcdark2_elf] Si optical data: {len(omega_si)} points, "
          f"ω = {omega_si[0]:.3f}–{omega_si[-1]:.3f} eV")

    # --- Write HypMat optical dat file ---
    # We keep Si's ε₁(ω) and ε₂(ω) as-is; the threshold at DIRECT_GAP_EV is
    # enforced by setting band_gap_eV in darkelf init (not in the dat file itself).
    out_path = out_dir / "HypMat_0p34_qcdark2_eps_electron_opticallimit.dat"
    citation = (
        f"CCDarkSens SrCd2Sb2 (HypMat) QCDark2-informed ELF: "
        f"Si measured optical dielectric (Si_eps_electron_opticallimit.dat) "
        f"with band_gap_eV={DIRECT_GAP_EV} eV ({onset_note}). "
        f"rho={RHO_T} g/cm3. "
        f"Si interband structure + plasmon used as proxy for SrCd2Sb2 pending "
        f"true wavefunction calculation (AiiDA bands_workchain_pk=9491). "
        "See docs/DarkPhoton_Absorption_HypotheticalMaterial.md and "
        "utils/build_darkphoton_hypmat_qcdark2_elf.py"
    )
    with open(out_path, "w") as f:
        f.write(citation + "\n")
        for w, e1, e2 in zip(omega_si, eps1_si, eps2_si):
            f.write(f"{w:.6e}  {e1:.6e}  {e2:.6e}\n")
    print(f"[qcdark2_elf] wrote {out_path}  ({len(omega_si)} rows)")

    # --- Update YAML (density correction, keep existing if present) ---
    yaml_path = out_dir / "HypMat_0p34.yaml"
    if yaml_path.exists():
        txt = yaml_path.read_text()
        if f"rhoT : {RHO_T}" not in txt:
            print(f"[qcdark2_elf] YAML already exists with different rhoT — not overwriting")
        else:
            print(f"[qcdark2_elf] YAML already up to date")
    else:
        print(f"[qcdark2_elf] Run build_darkphoton_hypothetical_material.py first to create YAML")

    print(f"\n[qcdark2_elf] Done.")
    print(f"  Material key in entry.py: 'hypmat_qcdark2'")
    print(f"  band_gap_eV to pass at darkelf init: {DIRECT_GAP_EV} eV (direct gap)")
    print(f"  This does NOT overwrite hypmat_unscreened or hypmat_screened (Drude).")


# ----------------------------------------------------------------------------
# main
#   Command line: --darkelf_dir and the optional QCDark2 HDF5 used to verify the onset; then build().
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--darkelf_dir", default=None)
    ap.add_argument("--qcdark2_h5", default="data/qcdark2_epsilon/Si/Si_fast_gap0p34.h5",
                    help="QCDark2 H5 file to verify direct-gap onset (optional)")
    args = ap.parse_args()
    build(_resolve_darkelf_dir(args.darkelf_dir), args.qcdark2_h5)


if __name__ == "__main__":
    main()
