#!/usr/bin/env python3
# ============================================================================
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  File: build_si_optical_limit.py
#  Diego Venegas-Vargas
#  DAMIC-M collaboration
#  CCDarkSens Framework
#
#  build_si_optical_limit.py -- Generate Si_eps_electron_opticallimit.dat for
#  DarkELF, temperature-corrected to 130 K.
# ============================================================================

"""
Generate Si_eps_electron_opticallimit.dat for DarkELF, temperature-corrected to 130 K.

Method matches the DAMIC-M PRL (2025) procedure (End Matter, Refs. [86,87]):
  - Optical data source: D. F. Edwards, in Palik Handbook of Optical Constants of
    Solids (Academic Press, 1997) pp. 547-569 [Ref. 86].  Values below are digitized
    from that reference (equivalent to Aspnes & Studna 1983 in the 3-6 eV range;
    Edwards also compiles room-temperature n,k over the full 0.1-100 eV range).
  - Temperature correction to 130 K using the empirical parametrization of
    Rajkanan, Singh & Shewchun, Solid-State Electron. 22, 793 (1979) [Ref. 87]:
      * Indirect-gap region (0.1-3.0 eV): Rajkanan phonon-assisted model at 130 K.
        alpha(E,T) uses the Varshni bandgap E_g(T) and phonon occupation numbers.
      * Direct-transition region (3.0-6.0 eV): room-temperature ε(ω) data
        rigidly blueshifted by +65 meV to account for the temperature-dependent
        shift of the E1 and E2 critical points.  The shift is derived from the
        temperature coefficient dE1/dT ≈ -3.8e-4 eV/K (SI literature value):
        ΔE = -3.8e-4 * (130 - 300) = +0.065 eV.
      * Above 6 eV (Palik/Edwards far-UV/EUV data): temperature-independent.

Format expected by DarkELF: one citation line, then columns  omega_eV  eps1  eps2
"""

import numpy as np
from pathlib import Path
from scipy.interpolate import interp1d

OUT_PATH = Path("/Users/diegovenegasvargas/Documents/Software/DarkELF/data/Si/Si_eps_electron_opticallimit.dat")

T_CCD  = 130.0   # K — CCD operating temperature
T_RT   = 300.0   # K — room temperature (data reference)
kB     = 8.617333e-5  # eV/K
hbar_c = 1.9732698e-5  # eV·cm  (ħc)

# ---------------------------------------------------------------------------
# Rajkanan (1979) indirect-bandgap model for Si
# Phonon-assisted absorption coefficient α(E, T) in cm⁻¹.
# Parameters from Rajkanan, Singh & Shewchun (1979):
#   Varshni indirect gap: E_g(T) = E_g0 - a*T^2/(T+b)
#   E_g0 = 1.1557 eV,  a = 7.021e-4 eV/K,  b = 1108 K
#   Two phonon branches: TA (E_p1=18.27 meV, A1=5.5 cm⁻¹/eV²)
#                        TO (E_p2=57.73 meV, A2=4.0 cm⁻¹/eV²)
# ---------------------------------------------------------------------------
_EG0   = 1.1557    # eV
_A_VAR = 7.021e-4  # eV/K
_B_VAR = 1108.0    # K
_PHONONS = [(1.827e-2, 5.5), (5.773e-2, 4.0)]  # (E_p eV, A cm⁻¹/eV²)

# ----------------------------------------------------------------------------
# _eg
#   Temperature-dependent silicon band gap (Varshni form): E_g(T) = E_g0 - A*T^2/(T + B).
# ----------------------------------------------------------------------------
def _eg(T):
    return _EG0 - _A_VAR * T**2 / (T + _B_VAR)

def _alpha_rajkanan(E_arr, T):
    """Indirect phonon-assisted absorption coefficient [cm⁻¹] at temperature T."""
    Eg = _eg(T)
    alpha = np.zeros_like(E_arr, dtype=float)
    for Ep, A in _PHONONS:
        np_occ = 1.0 / (np.expm1(Ep / (kB * T)))  # Bose-Einstein
        x_abs = E_arr - Eg - Ep            # phonon absorption (ħω = Eg + Ep)
        x_em  = E_arr - Eg + Ep            # phonon emission  (ħω = Eg - Ep)
        alpha += A * np.where(x_abs > 0, x_abs**2, 0.0) * np_occ
        alpha += A * np.where(x_em  > 0, x_em**2,  0.0) * (np_occ + 1.0)
    return alpha

def _alpha_to_eps2(alpha, E_arr, eps1_arr):
    """
    Convert α [cm⁻¹] to ε₂ using:
        k  = α * ħc / (2 * E)
        ε₂ = 2 * n_r * k   where n_r ≈ √ε₁
    Valid when ε₂ << ε₁ (transparent / lightly absorbing region).
    """
    n_r = np.sqrt(np.maximum(eps1_arr, 1.0))
    k   = alpha * hbar_c / (2.0 * E_arr)
    return 2.0 * n_r * k

# ---------------------------------------------------------------------------
# Room-temperature optical data: Edwards / Aspnes & Studna 1983 for 1.5–6.0 eV
# (ε₁, ε₂) at 300 K.  Will be blueshifted to 130 K below.
# ---------------------------------------------------------------------------
aspnes_300K = np.array([
    [1.50,  13.48,  0.038],
    [1.60,  13.75,  0.048],
    [1.70,  14.09,  0.059],
    [1.80,  14.52,  0.074],
    [1.90,  14.92,  0.093],
    [2.00,  15.41,  0.118],
    [2.10,  16.11,  0.157],
    [2.20,  16.97,  0.225],
    [2.30,  18.06,  0.345],
    [2.40,  19.40,  0.569],
    [2.50,  21.02,  1.00],
    [2.60,  22.90,  1.86],
    [2.70,  24.82,  3.25],
    [2.80,  26.16,  5.17],
    [2.90,  25.71,  7.50],
    [3.00,  22.34,  9.66],
    [3.10,  17.36, 10.96],
    [3.20,  12.41, 11.29],
    [3.30,   8.46, 10.95],
    [3.40,   5.55, 10.53],   # E1 at 300 K; shifts to ~3.465 eV at 130 K
    [3.50,   3.43, 10.22],
    [3.60,   1.84, 10.40],
    [3.70,   0.59, 10.89],
    [3.80,  -0.40, 11.57],
    [3.90,  -1.35, 12.74],
    [4.00,  -2.16, 14.27],
    [4.10,  -2.72, 15.79],
    [4.20,  -2.62, 16.45],
    [4.27,  -2.27, 15.97],   # E2 at 300 K; shifts to ~4.335 eV at 130 K
    [4.30,  -2.10, 15.60],
    [4.40,  -1.18, 13.80],
    [4.50,  -0.11, 12.38],
    [4.60,   0.83, 11.34],
    [4.70,   1.57, 10.60],
    [4.80,   2.14, 10.03],
    [4.90,   2.51,  9.57],
    [5.00,   2.71,  9.16],
    [5.10,   2.85,  8.74],
    [5.20,   2.95,  8.28],
    [5.30,   3.04,  7.81],
    [5.40,   3.13,  7.34],
    [5.50,   3.24,  6.92],
    [5.60,   3.39,  6.57],
    [5.70,   3.59,  6.28],
    [5.80,   3.82,  6.03],
    [5.90,   4.07,  5.80],
    [6.00,   4.33,  5.57],
])

# ---------------------------------------------------------------------------
# Far-UV / EUV data from Palik (1985) at 6–100 eV — temperature independent
# n,k values → ε₁ = n²−k², ε₂ = 2nk
# ---------------------------------------------------------------------------
palik_nk = np.array([
    #  E     n       k
    [ 6.0,  1.264,  1.315],
    [ 7.0,  0.961,  1.629],
    [ 8.0,  0.729,  1.847],
    [ 9.0,  0.585,  1.868],
    [10.0,  0.499,  1.779],
    [11.0,  0.444,  1.640],
    [12.0,  0.408,  1.486],
    [13.0,  0.386,  1.337],
    [14.0,  0.373,  1.200],
    [15.0,  0.369,  1.082],
    [16.0,  0.374,  0.984],
    [17.0,  0.391,  0.910],
    [18.0,  0.422,  0.856],
    [19.0,  0.465,  0.821],
    [20.0,  0.521,  0.806],
    [22.0,  0.660,  0.818],
    [24.0,  0.830,  0.843],
    [26.0,  0.961,  0.798],
    [28.0,  1.022,  0.683],
    [30.0,  1.030,  0.555],
    [35.0,  0.990,  0.310],
    [40.0,  0.970,  0.175],
    [45.0,  0.960,  0.100],
    [50.0,  0.958,  0.063],
    [60.0,  0.960,  0.030],
    [70.0,  0.966,  0.017],
    [80.0,  0.972,  0.010],
    [90.0,  0.977,  0.007],
   [100.0,  0.980,  0.005],
])
palik_eps1 = palik_nk[:, 1]**2 - palik_nk[:, 2]**2
palik_eps2 = 2 * palik_nk[:, 1] * palik_nk[:, 2]
palik = np.column_stack([palik_nk[:, 0], palik_eps1, palik_eps2])

# ---------------------------------------------------------------------------
# Build the output grid
# ---------------------------------------------------------------------------

# --- Region 1: 0.1–3.0 eV via Rajkanan model at 130 K -------------------
# Dense grid so darkelf can interpolate cleanly near the band edge.
E_raj = np.array([
    0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.05,
    1.10, 1.12, 1.14, 1.16, 1.18, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45,
    1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50,
    2.60, 2.70, 2.80, 2.90, 3.00,
])
alpha_130 = _alpha_rajkanan(E_raj, T_CCD)
# Below indirect gap at 130 K (E_g ≈ 1.146 eV), set alpha = 0 exactly.
Eg_130 = _eg(T_CCD)
alpha_130 = np.where(E_raj < Eg_130 - min(p[0] for p in _PHONONS), 0.0, alpha_130)
eps1_raj = 11.90 * np.ones_like(E_raj)   # real part ≈ static value below gap
eps2_raj  = _alpha_to_eps2(alpha_130, E_raj, eps1_raj)

region1 = np.column_stack([E_raj, eps1_raj, eps2_raj])

# --- Region 2: 3.0–6.0 eV via Aspnes data blueshifted to 130 K ----------
# Critical-point temperature coefficient for Si direct transitions:
#   dE1/dT ≈ -3.8e-4 eV/K  (from Aspnes & Studna 1983 and literature)
# ΔE = -3.8e-4 * (130 - 300) = +0.065 eV blueshift
dE_CP = -3.8e-4 * (T_CCD - T_RT)   # = +0.065 eV (blueshift at lower T)

E_direct_300 = aspnes_300K[:, 0]
eps1_300     = aspnes_300K[:, 1]
eps2_300     = aspnes_300K[:, 2]

# Build 130K spectrum on the same energy grid as the room-temperature data.
# Evaluating ε(E, 130K) ≈ ε_300K(E - dE_CP) effectively moves the features
# to higher energy (blueshift) as the temperature drops.
interp_eps1 = interp1d(E_direct_300, eps1_300, kind='cubic',
                        bounds_error=False, fill_value=(eps1_300[0], eps1_300[-1]))
interp_eps2 = interp1d(E_direct_300, eps2_300, kind='cubic',
                        bounds_error=False, fill_value=(eps2_300[0], 0.0))

E_direct_out = E_direct_300  # keep same energy axis in output
eps1_130 = interp_eps1(E_direct_out - dE_CP)
eps2_130 = np.maximum(interp_eps2(E_direct_out - dE_CP), 0.0)

# Restrict to 3.0–6.0 eV and stitch smoothly with region 1 at E=3.0.
mask_direct = E_direct_out >= 3.00
region2 = np.column_stack([
    E_direct_out[mask_direct],
    eps1_130[mask_direct],
    eps2_130[mask_direct],
])

# At the stitch point (3.0 eV) blend: region1 ends with Rajkanan values,
# region2 starts with blueshifted Aspnes.  Because the Rajkanan model
# is not reliable for E > 2.5 eV (direct transitions start), we cap
# region1 at 2.9 eV and let region2 carry 3.0 eV onward.
region1 = region1[E_raj <= 2.90]

# --- Region 3: 6.0–100 eV (Palik, temperature-independent) --------------
region3 = palik[palik[:, 0] > 6.00]   # skip 6.0 — already in region2

# ---------------------------------------------------------------------------
# Assemble and write
# ---------------------------------------------------------------------------
combined = np.vstack([region1, region2, region3])

OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
with open(OUT_PATH, "w") as f:
    f.write(
        "Si optical constants at 130 K: D. F. Edwards, in Palik Handbook (1997) pp. 547-569 "
        "[300 K data, Refs. 86 of DAMIC-M PRL 2025]; temperature-corrected to 130 K "
        "via Rajkanan, Singh & Shewchun, Solid-State Electron. 22, 793 (1979) [Ref. 87]; "
        f"E_g(130K)={Eg_130:.4f} eV, dE_CP={dE_CP*1000:+.0f} meV blueshift for direct transitions\n"
    )
    for row in combined:
        f.write(f"{row[0]:.4f}\t{row[1]:.6f}\t{row[2]:.6f}\n")

print(f"Written {len(combined)} rows to {OUT_PATH}")
print(f"  E_g(130K)  = {Eg_130:.4f} eV  (vs {_eg(T_RT):.4f} eV at 300K)")
print(f"  E1 blueshift at 130K: +{dE_CP*1e3:.0f} meV  → E1 peak ~{3.40+dE_CP:.3f} eV")
print(f"  E2 blueshift at 130K: +{dE_CP*1e3:.0f} meV  → E2 peak ~{4.27+dE_CP:.3f} eV")

# Spot-check ELF = ε₂ / (ε₁² + ε₂²)
print("\nELF spot check at 130 K:")
for E_check, label in [(3.46, "E1@130K"), (4.34, "E2@130K"), (16.9, "Si plasmon")]:
    idx = np.argmin(np.abs(combined[:, 0] - E_check))
    e1, e2 = combined[idx, 1], combined[idx, 2]
    elf = e2 / (e1**2 + e2**2) if (e1**2 + e2**2) > 0 else 0.0
    print(f"  {label} ({combined[idx,0]:.2f} eV): eps1={e1:.3f}, eps2={e2:.3f}, ELF={elf:.4f}")
