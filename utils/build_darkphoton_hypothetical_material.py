#!/usr/bin/env python3
"""
Generate synthetic DarkELF data files for SrCd₂Sb₂ (codename HypMat) with:
  E_gap    = 0.34 eV   (R2SCAN no-SOC indirect gap; direct gap ~0.695 eV)
  rho_T    = 5.76 g/cm³ (from lattice params a=4.44 Å, c=28.196 Å, MW=555.96 g/mol, Z=3)
  sigma_DC = 300 Ω⁻¹cm⁻¹  (DC optical conductivity from screening spreadsheet)
  Z_val    = 16  (Sr: 2s, 2×Cd: 2s each, 2×Sb: 5 valence each — outer shell only)

Two bracketing scenarios for the unknown high-frequency dielectric constant ε∞:
  - eps1_1  : ε∞ = 1   (unscreened; upper bound on absorption rate)
  - eps1_12 : ε∞ = 12  (Si-like screening; lower bound on absorption rate)

Dielectric model — Drude (replaces the previous step-function σ_DC model):
---------------------------------------------------------------------------
The previous model set ε₂(ω) = σ_DC/(ε₀ω) for ω ≥ E_gap (constant conductivity).
This gives a monotonically falling ELF with no plasmon feature. The Drude model
is physically better motivated and produces a sharp plasmon peak in the ELF:

  ε(ω) = ε∞ - ωp²/(ω² + iγω)

  ε₁(ω) = ε∞ - ωp²/(ω² + γ²)
  ε₂(ω) = ωp²γ / (ω(ω² + γ²))

Parameters derived from known material properties:
  ωp  = sqrt(n_e e²/(ε₀ mₑ))  — free-carrier plasma frequency
      ≈ 11.7 eV  (from rho=5.76 g/cm³, MW=555.96 g/mol, Z_val=16)
  γ   = ε₀ ωp² / σ_DC          — Drude damping from DC conductivity
      ≈ 0.062 eV  (from σ_DC=300 Ω⁻¹cm⁻¹; narrow resonance, sharp peak)

Screened plasmon position (where ELF peaks, Re(ε)=0):
  ωs = sqrt(ωp²/ε∞ - γ²) ≈ ωp/sqrt(ε∞)
  ε∞=1  → ωs ≈ 11.7 eV   (unscreened — peak near 12 eV)
  ε∞=12 → ωs ≈  3.4 eV   (screened — peak near Si gap region)

Both ε₁ and ε₂ are now frequency-dependent (Drude), producing a realistic
plasmon resonance instead of the featureless 1/ω falloff. The ε∞ bracket
controls the plasmon position, not just the overall scale.

Note: the ε₂ is identical for both brackets; only ε₁ (= ε∞ + Drude real part)
differs. Both brackets share the same ωp and γ — ε∞ is the only free parameter.

Outputs (under DARKELF_DIR/data/HypMat_0p34/):
  HypMat_0p34.yaml
  HypMat_0p34_eps_electron_opticallimit_eps1_1.dat
  HypMat_0p34_eps_electron_opticallimit_eps1_12.dat

Usage:
  python3 utils/build_darkphoton_hypothetical_material.py \\
      --darkelf_dir /path/to/DarkELF

Or set CCDARK_SENS_DARKELF_DIR in your environment and omit --darkelf_dir.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

_ENV_DARKELF_DIR = "CCDARK_SENS_DARKELF_DIR"

# --- Material parameters (SrCd₂Sb₂) ---
E_GAP_EV     = 0.34      # eV — R2SCAN indirect band gap (no SOC)
RHO_T        = 5.76      # g/cm³ — from lattice params (corrected from 8.0 placeholder)
SIGMA_DC_CGS = 300.0     # Ω⁻¹cm⁻¹ — DC optical conductivity
MW_G_MOL     = 555.96    # g/mol — molecular weight of SrCd₂Sb₂
Z_VAL        = 16        # valence electrons per formula unit (Sr:2 + 2×Cd:4 + 2×Sb:10)

# --- Physical constants (SI) ---
E_CHARGE     = 1.602e-19  # C
EPS0_SI      = 8.854e-12  # F/m
M_E_KG       = 9.109e-31  # kg
HBAR_EV_S    = 6.582e-16  # eV·s
N_AVO        = 6.022e23   # mol⁻¹

# --- Drude parameters derived from material properties ---
SIGMA_DC_SI = SIGMA_DC_CGS * 100.0  # S/m

# Valence electron number density (m⁻³)
RHO_KG_M3   = RHO_T * 1e3             # kg/m³
MW_KG_MOL   = MW_G_MOL * 1e-3         # kg/mol
N_E_M3      = (RHO_KG_M3 / MW_KG_MOL) * N_AVO * Z_VAL  # electrons/m³

# Plasma frequency (rad/s then eV)
OMEGAP_RAD_S = np.sqrt(N_E_M3 * E_CHARGE**2 / (EPS0_SI * M_E_KG))
OMEGAP_EV    = OMEGAP_RAD_S * HBAR_EV_S  # ≈ 11.7 eV

# Drude damping from σ_DC = ε₀ ωp² / γ  →  γ = ε₀ ωp² / σ_DC
GAMMA_RAD_S  = EPS0_SI * OMEGAP_RAD_S**2 / SIGMA_DC_SI
GAMMA_EV     = GAMMA_RAD_S * HBAR_EV_S   # ≈ 0.062 eV

print(f"[HypMat] ωp  = {OMEGAP_EV:.3f} eV")
print(f"[HypMat] γ   = {GAMMA_EV:.4f} eV")
print(f"[HypMat] Screened plasmon (ε∞=1):  ωs ≈ {OMEGAP_EV:.2f} eV")
print(f"[HypMat] Screened plasmon (ε∞=12): ωs ≈ {OMEGAP_EV/np.sqrt(12):.2f} eV")

# --- Placeholder material constants (needed by darkelf init, not absorption) ---
LATTICE_SPACING_ANG  = 6.5     # Å
CLA_KMS              = 3.0     # LA sound speed km/s
CTA_KMS              = 1.5     # TA sound speed km/s
A_PLACEHOLDER        = MW_G_MOL / 5.0  # avg atomic mass (5 atoms per f.u.)

# --- ω grid (eV) ---
OMEGA_BELOW  = np.array([0.05, 0.10, 0.20, 0.33])
OMEGA_AT_GAP = np.array([E_GAP_EV])
OMEGA_ABOVE  = np.concatenate([
    np.linspace(E_GAP_EV + 0.01, 1.0, 50),
    np.linspace(1.05, 20.0, 100),
    np.linspace(20.5, 100.0, 30),
])
OMEGA_ALL = np.concatenate([OMEGA_BELOW, OMEGA_AT_GAP, OMEGA_ABOVE])


def drude_eps(omega: np.ndarray, eps_inf: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Drude dielectric function components.

    ε₁(ω) = ε∞ - ωp²/(ω² + γ²)
    ε₂(ω) = ωp²γ / (ω(ω² + γ²))    [zero for ω < E_gap]

    Returns (ε₁, ε₂) arrays.
    """
    wp2  = OMEGAP_EV**2
    gam  = GAMMA_EV
    denom = omega**2 + gam**2

    eps1 = np.full_like(omega, eps_inf) - wp2 / denom

    eps2 = np.zeros_like(omega)
    mask = omega >= E_GAP_EV
    eps2[mask] = wp2 * gam / (omega[mask] * denom[mask])

    return eps1, eps2


def _resolve_darkelf_dir(explicit: str | None) -> str:
    if explicit:
        return os.path.abspath(os.path.expanduser(explicit))
    env = os.environ.get(_ENV_DARKELF_DIR)
    if env:
        return os.path.abspath(os.path.expanduser(env.strip()))
    raise FileNotFoundError(
        f"Pass --darkelf_dir or set {_ENV_DARKELF_DIR}."
    )


def _write_optical_dat(path: str, omega: np.ndarray, eps_inf: float, citation: str):
    eps1_vals, eps2_vals = drude_eps(omega, eps_inf)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(citation + "\n")
        for w, e1, e2 in zip(omega, eps1_vals, eps2_vals):
            f.write(f"{w:.6e}  {e1:.6e}  {e2:.6e}\n")
    print(f"  wrote {path}  ({len(omega)} rows)")


def _write_yaml(path: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    lines = [
        f"rhoT : {RHO_T}  # g/cm³ — from lattice params a=4.44A c=28.196A MW=555.96 Z=3",
        f"E_gap : {E_GAP_EV}  # eV — R2SCAN indirect gap (no SOC)",
        "e0 : 1.6  # eV — charge-pair creation energy",
        f"A : {A_PLACEHOLDER:.2f}  # avg atomic mass per atom (MW/5 atoms per f.u.)",
        f"omegap : {OMEGAP_EV:.3f}  # eV — Drude plasma frequency from valence electron density",
        f"lattice_spacing : {LATTICE_SPACING_ANG}  # Å placeholder",
        f"cLAkms : {CLA_KMS}  # km/s placeholder",
        f"cTAkms : {CTA_KMS}  # km/s placeholder",
        "ombar : 0.03  # eV — average phonon energy placeholder",
        "Enl_list : []",
        "LOvec : [0.03]",
        "atoms : ['HypMat']",
        "unitcell : {'HypMat': {'A': 150.0, 'mult': 1}}",
        f"# SrCd2Sb2 (codename HypMat): E_gap={E_GAP_EV} eV, rho={RHO_T} g/cm3,",
        f"# sigma_DC={SIGMA_DC_CGS} Ohm^-1 cm^-1, Drude model: omegap={OMEGAP_EV:.3f} eV, gamma={GAMMA_EV:.4f} eV",
        "# See utils/build_darkphoton_hypothetical_material.py",
    ]
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  wrote {path}")


def build(darkelf_dir: str):
    data_dir = os.path.join(darkelf_dir, "data", "HypMat_0p34")
    print(f"Building HypMat_0p34 in {data_dir}")

    _write_yaml(os.path.join(data_dir, "HypMat_0p34.yaml"))

    citation_base = (
        f"CCDarkSens SrCd2Sb2 (HypMat) Drude model: "
        f"rho={RHO_T} g/cm3, E_gap={E_GAP_EV} eV, sigma_DC={SIGMA_DC_CGS} Ohm^-1cm^-1, "
        f"omegap={OMEGAP_EV:.3f} eV, gamma={GAMMA_EV:.4f} eV. "
        "See utils/build_darkphoton_hypothetical_material.py"
    )

    _write_optical_dat(
        os.path.join(data_dir, "HypMat_0p34_eps_electron_opticallimit_eps1_1.dat"),
        OMEGA_ALL, eps_inf=1.0,
        citation=citation_base + " | eps_inf=1 (unscreened upper bound, plasmon ~11.7 eV)",
    )
    _write_optical_dat(
        os.path.join(data_dir, "HypMat_0p34_eps_electron_opticallimit_eps1_12.dat"),
        OMEGA_ALL, eps_inf=12.0,
        citation=citation_base + " | eps_inf=12 (Si-like screened lower bound, plasmon ~3.4 eV)",
    )
    print("Done.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--darkelf_dir", default=None,
                    help="Path to DarkELF repo root (default: $CCDARK_SENS_DARKELF_DIR)")
    args = ap.parse_args()
    build(_resolve_darkelf_dir(args.darkelf_dir))


if __name__ == "__main__":
    main()
