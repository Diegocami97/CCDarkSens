# Dark Photon Absorption — SrCd₂Sb₂ (HypMat)
## Physics Framework, Formal Calculations, and Assumptions

**Document:** `docs/DarkPhoton_Absorption_HypotheticalMaterial.md`  
**Author:** Diego Venegas-Vargas  
**Date:** 2026-06-23  
**Status:** Phenomenological study — bracket calculation pending full material characterization

---

## Changelog

| Date | Change |
|---|---|
| 2026-07-03 | Material identified as **SrCd₂Sb₂** (from collaboration AM₂X₂ screening, R-3m structure) |
| 2026-07-03 | **Density corrected**: 8.0 → **5.76 g/cm³** (from lattice params a=4.44 Å, c=28.196 Å, MW=555.96 g/mol, Z=3) |
| 2026-07-03 | **Optical model upgraded**: step-function σ_DC replaced by **Drude model** (see §A1-updated below) |
| 2026-07-03 | Band structure received (R2SCAN no-SOC): indirect gap 0.560 eV, direct gap 0.695 eV — using 0.34 eV pending SOC correction |
| 2026-07-03 | **QCDark2-informed ELF** added as third model variant (`hypmat_qcdark2`) — see §QCDark2-ELF below |

---

## A1-updated — Drude dielectric model (replaces step-function assumption)

### Why the step-function model was inadequate

The original Assumption A1 set σ₁(ω) = σ_DC = constant for ω ≥ E_gap, giving:

$$\varepsilon_2(\omega) = \frac{\sigma_\text{DC}}{\varepsilon_0\,\omega} \approx \frac{2.23}{\omega_\text{eV}}$$

with ε₁ = constant. The resulting ELF = Im[−1/ε] is **monotonically decreasing** with
no plasmon feature. Real materials have a sharp ELF peak at the screened plasmon
frequency, producing a characteristic sharp minimum in the sensitivity curve — visible
in Si (plasmon ~16 eV) but absent in the step-function model.

### Drude replacement

The Drude dielectric function with sub-gap cutoff:

$$\varepsilon_1(\omega) = \varepsilon_\infty - \frac{\omega_p^2}{\omega^2 + \gamma^2}$$

$$\varepsilon_2(\omega) = \begin{cases} \dfrac{\omega_p^2\,\gamma}{\omega(\omega^2 + \gamma^2)} & \omega \geq E_\text{gap} \\ 0 & \omega < E_\text{gap} \end{cases}$$

Both ε₁ and ε₂ are now frequency-dependent. The ELF peaks sharply at the screened
plasmon ωs ≈ ωp/√ε∞, producing the physically expected sharp minimum in sensitivity.

### Drude parameters for SrCd₂Sb₂

All parameters derived from known material properties — **no additional DFT input needed**:

**Plasma frequency ωp** from valence electron density (Z_val = 16 per formula unit:
Sr 2s, 2×Cd 2s each, 2×Sb 5 valence each):

$$n_e = \frac{\rho\,N_A}{M_W}\,Z_\text{val} = \frac{5.76 \times 6.022\times10^{23}}{555.96}\times 16 \approx 9.98\times10^{22}\ \text{cm}^{-3}$$

$$\omega_p = \sqrt{\frac{n_e\,e^2}{\varepsilon_0\,m_e}} \approx 11.7\ \text{eV}$$

**Drude damping γ** from DC conductivity via σ_DC = ε₀ωp²/γ:

$$\gamma = \frac{\varepsilon_0\,\omega_p^2}{\sigma_\text{DC}} \approx 0.062\ \text{eV} \ll \omega_p$$

The narrow damping (γ/ωp ≈ 0.005) means a sharp, well-resolved ELF peak.

**Screened plasmon position** (where sensitivity minimum appears):

| Bracket | ε∞ | ωs = ωp/√ε∞ |
|---|---|---|
| Unscreened | 1 | ~11.7 eV |
| Screened (Si-like) | 12 | ~3.4 eV |

### What this model still misses

The Drude model is a significant improvement over the step-function but remains approximate:
- **Interband transitions** above the gap are not captured (Drude is a free-carrier model)
- **k-dependence**: ε(q,ω) at finite q for DM-electron rates requires QCDark2 output
- **ε∞**: bracketed at 1 and 12; the true value is unknown without DFT

The definitive calculation requires the full ε(q,ω) from QCDark2 using AiiDA
wavefunction files (bands_workchain_pk=9491). The Drude model is the best
physics-motivated placeholder available from the screening spreadsheet inputs alone.

---

## QCDark2-ELF variant (`hypmat_qcdark2`)

A third model variant sits between the featureless Drude and the true SrCd₂Sb₂ ELF.

### Motivation

Inspecting the q→0 slice of `Si_fast_gap0p34.h5` (the QCDark2 dielectric computed
for the scissors-shifted Si used in DM-e rate generation) reveals:

- ELF onset at **~2.1 eV** — the scissors-shifted Si *direct* gap
  (Si direct gap ~3.4 eV minus scissors correction ~1.3 eV)
- Si-like interband structure from 2–15 eV
- **Sharp Si plasmon peak at ~19.5 eV** (ELF ≈ 34.5)

The "fast" k-grid is too coarse to use ε(q,ω) directly as an optical ELF (ε₁ shows
unphysical values of 34–82 below the gap due to sparse k-sampling). Instead we use:

**Si's measured optical dielectric function** (`Si_eps_electron_opticallimit.dat`,
the same file used for the DAMIC-M 2025 PRL Si hidden-photon limit), with
`band_gap_eV = 2.1 eV` set at darkelf init time to enforce the scissors-derived onset.

### What this gives

- Correct absorption onset at 2.1 eV (direct gap, consistent with scissors result)
- Real Si interband structure and sharp Si plasmon at ~17 eV
- A **sharp minimum in the dark photon sensitivity** near 17 eV — qualitatively more
  realistic than the featureless Drude curves
- rho = 5.76 g/cm³ (correct SrCd₂Sb₂ density)

### Caveats

- The interband structure belongs to Si, not SrCd₂Sb₂
- The plasmon position reflects Si's electron density, not SrCd₂Sb₂'s
- The scissors correction is approximate for optical transitions
- The true ELF will differ quantitatively once AiiDA wavefunctions are available

### Three-model summary

| Key | ELF model | Plasmon | Status |
|---|---|---|---|
| `hypmat_unscreened` | Drude (ε∞=1), σ_DC=300 Ω⁻¹cm⁻¹ | None (overdamped) | Current default |
| `hypmat_screened` | Drude (ε∞=12), σ_DC=300 Ω⁻¹cm⁻¹ | None (overdamped) | Conservative bracket |
| `hypmat_qcdark2` | Si measured optical + scissors onset | Si plasmon ~17 eV | Intermediate — use while awaiting AiiDA wavefunctions |

### File generation

```bash
# Step 1 — build the .dat file (copies Si optical data with HypMat header)
python3 utils/build_darkphoton_hypmat_qcdark2_elf.py \
    --darkelf_dir /path/to/DarkELF \
    --qcdark2_h5 data/qcdark2_epsilon/Si/Si_fast_gap0p34.h5

# Step 2 — generate rate tables
python3 utils/darkphoton_generate_grid.py configs/darkphoton_generate_hypmat_qcdark2.json

# Step 3 — run scan
./build/ccdarksens_scan_dmelectron_pattern configs/darkphoton_scan_hypmat_qcdark2_ne.json
```

---

## 1. Motivation

The CCDarkSens band-gap phenomenology study (heavy and light mediator,
`docs/band_gap_pheno_ionization.md`) explores how DM-electron scattering
sensitivity changes when the semiconductor target has a lower band gap
than silicon. The same physics question applies to a complementary channel:
**hidden-photon (dark photon) absorption**, where the signal is a single
ionization event at energy E = m_A' rather than a continuous recoil spectrum.

A hypothetical target material with the following properties has been identified:

| Property | Symbol | Value | Units |
|---|---|---|---|
| Band gap | E_gap | 0.34 | eV |
| Mean e-h pair energy | ε_h | 1.6 | eV |
| DC optical conductivity | σ₁(DC) | 300 | Ω⁻¹cm⁻¹ |

This document records the complete formal framework for computing the dark
photon absorption rate for this material, explicitly labeling all
approximations made when only these three numbers are available.

---

## 2. Dark Photon Absorption: Physical Process

### 2.1 The signal

A kinetically-mixed hidden photon (dark photon) A' of mass m_A' and
kinetic-mixing parameter κ can be absorbed by an electron in the target
material if m_A' ≥ E_gap. The full rest energy m_A' is deposited as a
single ionization event — unlike DM-electron scattering, which produces a
continuous recoil spectrum. The signal is therefore **monochromatic** at
E_deposit = m_A'.

### 2.2 Coupling to matter

The dark photon A' couples to ordinary photons through kinetic mixing:

$$
\mathcal{L} \supset -\frac{\kappa}{2}\, F_{\mu\nu}^{A'} F^{\mu\nu}_\gamma
$$

In the presence of a medium, the dark photon acquires a medium-induced
self-energy whose imaginary part determines the absorption rate. In the
non-relativistic, long-wavelength limit (k → 0, q → 0) appropriate
for m_A' in the eV range, this self-energy is directly related to the
complex dielectric function of the target material ε(ω).

---

## 3. Absorption Rate: Formal Expression

### 3.1 General formula

The dark photon absorption rate per unit target mass is (Knapen, Kozaczuk
& Lin 2021, Eq. 2.6; Hochberg et al. 2021):

$$
R_\text{abs}(m_{A'}) = \kappa^2 \cdot \frac{\rho_\chi}{\rho_T\, m_{A'}}
\cdot \text{Im}\!\left[-\frac{1}{\varepsilon(m_{A'}, 0)}\right]
\cdot \frac{1}{\text{eV}} \times (\text{kg·yr})^{-1}
$$

where:
- κ — kinetic-mixing parameter (dimensionless)
- ρ_χ — local DM density = 0.3 GeV/cm³ = 0.3×10⁹ eV/cm³
- ρ_T — target material mass density [g/cm³]
- m_A' — dark photon mass [eV]
- ε(ω, k=0) — complex dielectric function in the optical (k→0) limit

### 3.2 Energy-loss function (ELF)

The quantity Im[−1/ε] is the **energy-loss function** (ELF):

$$
\text{ELF}(\omega) \equiv \text{Im}\!\left[-\frac{1}{\varepsilon(\omega)}\right]
= \frac{\varepsilon_2(\omega)}{\varepsilon_1^2(\omega) + \varepsilon_2^2(\omega)}
$$

where ε₁(ω) and ε₂(ω) are the real and imaginary parts of the dielectric
function. The ELF is dimensionless and experimentally accessible via
electron energy-loss spectroscopy (EELS) or optical reflectance measurements.

### 3.3 Weak-absorption (semiconductor) limit

For a semiconductor with modest conductivity (ε₂ « ε₁ at photon energies
near E_gap), the ELF simplifies to:

$$
\text{ELF}(\omega) \approx \frac{\varepsilon_2(\omega)}{\varepsilon_1^2(\omega)}
\quad \text{when } \varepsilon_2 \ll \varepsilon_1
$$

**This limit will be checked a posteriori** once ε₁ is known (see §5.3).
For the present hypothetical material it is adopted as **Assumption A3**
(see §5).

---

## 4. Connecting Optical Conductivity to the Dielectric Function

### 4.1 Formal relation

The complex dielectric function and complex optical conductivity σ(ω) =
σ₁(ω) + iσ₂(ω) are related by (standard electrodynamics, SI units):

$$
\varepsilon(\omega) = \varepsilon_\infty + \frac{i\,\sigma(\omega)}{\varepsilon_0\,\omega}
$$

where ε_∞ is the high-frequency (above-band-gap) dielectric constant. Taking
imaginary parts:

$$
\boxed{\varepsilon_2(\omega) = \frac{\sigma_1(\omega)}{\varepsilon_0\,\omega}}
\quad \text{(SI units)}
$$

and taking real parts:

$$
\varepsilon_1(\omega) = \varepsilon_\infty - \frac{\sigma_2(\omega)}{\varepsilon_0\,\omega}
$$

### 4.2 Kramers-Kronig relation (what we would do with full spectral data)

Given the full frequency-dependent conductivity spectrum σ₁(ω'), ε₁(ω) can be
obtained rigorously via the Kramers-Kronig dispersion relation:

$$
\varepsilon_1(\omega) - 1 = \frac{2}{\pi}\,\mathcal{P}
\int_0^\infty \frac{\omega'\,\varepsilon_2(\omega')}{\omega'^2 - \omega^2}\,d\omega'
$$

This integral requires knowledge of ε₂(ω') = σ₁(ω')/ε₀ω' over all frequencies,
which in turn requires a measurement or ab initio calculation of σ₁(ω) as a
function of ω — not available for this hypothetical material.

### 4.3 Penn model estimate of ε₁(0) (what we would do with electron density)

In the absence of a full Kramers-Kronig integration, the static dielectric
constant can be estimated via the Penn model:

$$
\varepsilon_1(0) \approx 1 + \left(\frac{\hbar\omega_p}{E_\text{gap}}\right)^2
$$

where ω_p = √(n_e e²/ε₀ m_e) is the plasma frequency and n_e is the valence
electron density. This requires knowledge of n_e (electrons per unit volume),
which is not currently available for this material.

---

## 5. Assumptions Made in the Present Calculation

With only (E_gap = 0.34 eV, ε_h = 1.6 eV, σ₁(DC) = 300 Ω⁻¹cm⁻¹) available,
the following assumptions are made. Each is explicitly labeled and its impact
on the result is characterized.

---

### Assumption A1 — Frequency-independent conductivity

**Statement:** σ₁(ω) is treated as constant for ω ≥ E_gap:

$$
\sigma_1(\omega) = \begin{cases}
0 & \omega < E_\text{gap} \\
300\;\Omega^{-1}\text{cm}^{-1} & \omega \geq E_\text{gap}
\end{cases}
$$

**Physical justification:** The DC conductivity is a direct measure of σ₁ in
the ω → 0 limit. In a simple Drude metal or doped semiconductor at energies
well below the next higher interband transition, σ₁ varies slowly. For m_A'
values scanned near E_gap, this is plausible.

**When it breaks down:** At photon energies well above E_gap (say, several
× E_gap), real semiconductor conductivities typically rise substantially due
to interband transitions. The step-function model becomes an increasingly
poor approximation at higher m_A'. The direction of the bias is known:
the step-function underestimates σ₁(ω) at high ω (where σ₁ would rise) →
**underestimates the absorption rate at high m_A'**.

**Impact on result:** Primarily affects the shape of R(m_A') at m_A' ≫ E_gap.
Near threshold (m_A' ≈ E_gap), the impact is minimal.

---

### Assumption A2 — Static dielectric constant ε₁ = 1 (upper bracket) or ε₁ = 12 (lower bracket)

**Statement:** In the absence of a Kramers-Kronig calculation or Penn model
estimate (both of which require unavailable spectral data), ε₁(ω) is treated
as a constant. Two bracketing values are used:

| Scenario | ε₁ | Physical interpretation |
|---|---|---|
| **No screening** | 1.0 | Vacuum dielectric constant — vacuum coupling, zero medium polarization |
| **Si-like screening** | 12.0 | Typical semiconductor, ε₁ ∼ n² where n ≈ 3.5 |

The ELF scales as ε₂/ε₁² in the weak-absorption limit (A3), so the
absorption rate scales as **1/ε₁²**. A factor of 12 in ε₁ therefore shifts
the limit curve by a factor of **144** on the rate, or equivalently a factor
of **12** on κ — a very significant effect.

**What would pin this down:** A single measurement of the static refractive
index n (via ellipsometry or reflectance) immediately gives ε₁ = n², which
is the dominant uncertainty in this calculation.

---

### Assumption A3 — Weak-absorption limit (ε₂ « ε₁)

**Statement:** ELF ≈ ε₂/ε₁² is used rather than the exact expression.

**Check:** With σ₁ = 300 Ω⁻¹cm⁻¹ and ω = E_gap = 0.34 eV:

$$
\varepsilon_2(E_\text{gap}) = \frac{\sigma_1}{\varepsilon_0\,\omega}
= \frac{3\times10^4\;\text{S/m}}{8.854\times10^{-12}\;\text{F/m}
\;\times\; 0.34\;\text{eV}\;\times\; 1.519\times10^{15}\;\text{rad/s/eV}}
\approx \frac{3\times10^4}{4.57\times10^3} \approx 6.6
$$

For the Si-like screening scenario (ε₁ = 12), ε₂ ≈ 6.6 « ε₁ = 12 — marginal
but still in the weak-absorption regime. For the no-screening scenario
(ε₁ = 1), ε₂ = 6.6 > ε₁ → the weak-absorption approximation **breaks down**.
The exact ELF formula must be used in that case:

$$
\text{ELF}(E_\text{gap}) = \frac{6.6}{1 + 6.6^2} \approx \frac{6.6}{44.6} \approx 0.15
$$

This is significantly different from ε₂/1 = 6.6, so the ε₁=1 bracket
requires the exact ELF formula. The code handles this exactly (`darkelf`
uses the full ELF = ε₂/(ε₁²+ε₂²) always — see `darkelf/epsilon.py:339`).

---

### Assumption A4 — Target density ρ_T = 8 g/cm³

**Statement:** The mass density is taken to be 8 g/cm³, consistent with
heavy narrow-gap semiconductors in the same class (PbTe: 8.16 g/cm³,
HgTe: 8.09 g/cm³).

**Impact:** R_abs scales as 1/ρ_T linearly. A factor of 2 error in ρ_T
shifts the entire rate (and hence the κ limit curve) by a factor of √2 on
κ. This is the **smallest uncertainty** in the present calculation.

---

### Assumption A5 — ε_h enters the ionization yield only, not the absorption rate

**Statement:** The mean electron-hole pair energy ε_h = 1.6 eV determines
how the deposited energy m_A' maps to an integer number of electron-hole pairs
n_e = max(m_A' − E_gap, 0) / ε_h. It does **not** enter the absorption rate
formula.

**Where it enters the pipeline:** In `response.charge_ionization.eh_pair_eV`
in the CCDarkSens JSON scan config. The existing `ChargeIonization::FoldToNe`
infrastructure handles this exactly as it does for the Si band-gap pheno study,
with a rescaled p100K table built using `(E_gap, ε_h) = (0.34, 1.6)` eV.

---

## 6. Complete Calculation with Available Numbers

### 6.1 ε₂(ω) from σ₁(DC)

Unit conversion (SI → eV-based numerical value):

$$
\sigma_1 = 300\;\Omega^{-1}\text{cm}^{-1} = 3\times10^4\;\text{S/m}
$$

$$
\varepsilon_2(\omega) = \frac{3\times10^4}{8.854\times10^{-12}
\;\times\; \omega[\text{eV}]\;\times\; 1.5193\times10^{15}}
= \frac{3\times10^4}{1.345\times10^4\;\times\;\omega[\text{eV}]}
= \frac{2.23}{\omega[\text{eV}]}
$$

Evaluated at representative energies:

| m_A' [eV] | ε₂ (step-function model) |
|---|---|
| 0.34 (threshold) | 6.56 |
| 1.0 | 2.23 |
| 5.0 | 0.45 |
| 10.0 | 0.22 |
| 50.0 | 0.045 |

### 6.2 ELF(ω) for the two bracketing scenarios

$$
\text{ELF}(\omega) = \frac{\varepsilon_2(\omega)}{\varepsilon_1^2 + \varepsilon_2^2(\omega)}
$$

| m_A' [eV] | ELF (ε₁=1) | ELF (ε₁=12) |
|---|---|---|
| 0.34 | 0.148 | 0.046 |
| 1.0 | 0.910 | 0.015 |
| 5.0 | 0.41 | 0.0031 |
| 10.0 | 0.22 | 0.0015 |

Note: at m_A' = 1 eV, ε₂ = 2.23 > 1, so the ε₁=1 curve peaks near ε₂=1 and
the weak-absorption approximation is poor there. The full ELF formula is
required and is used in the code.

### 6.3 Absorption rate scaling

The rate scales as:

$$
R_\text{abs}(\kappa, m_{A'}) = \kappa^2 \cdot C \cdot \frac{\text{ELF}(m_{A'})}{m_{A'}}
$$

where the prefactor C depends only on (ρ_χ, ρ_T, unit conversions) — the
`darkelf`-internal prefactor `foo` in `absorption.py`:

$$
C = \frac{\rho_\chi}{\rho_T} \times \frac{c^3}{\hbar^2}
\times (\text{conversion to events/kg/yr})
$$

For ρ_T = 8 g/cm³ and ρ_χ = 0.3 GeV/cm³, `darkelf` evaluates C internally
using its unit conversions (`eVtoInvYr`, `c0cms`). The κ upper limit from
an exposure-limited measurement is then:

$$
\kappa^{90\%}(m_{A'}) = \sqrt{\frac{N_{90\%}}{R_\text{abs}(\kappa=1,\,m_{A'})
\;\times\; \text{exposure}\;[\text{kg·yr}]}}
$$

where N_90% = 2.44 events (Poisson 90% CL, zero observed background).

---

## 7. What the Complete Calculation Would Require

For the record, a fully rigorous calculation (no approximations) would need:

| Input | What it determines | How to get it |
|---|---|---|
| σ₁(ω) spectrum [full ω range] | ε₂(ω) via ε₂=σ₁/ε₀ω | Optical reflectance / ellipsometry / ab initio DFPT |
| ε₁(ω) spectrum | ELF numerator screening | Kramers-Kronig from σ₁(ω), or ellipsometry |
| ρ_T [g/cm³] | Rate normalization (linear) | Crystal density measurement or DFT lattice |
| Band structure / phonon spectrum | Migdal probability, phonon absorption at ω<E_gap | Full DFT+DFPT calculation |

The present calculation uses:
- σ₁(ω) → step function anchored at σ₁(DC) = 300 Ω⁻¹cm⁻¹ (A1)
- ε₁(ω) → two bracketing constants (A2): ε₁=1 and ε₁=12
- ρ_T → 8 g/cm³ (A4)
- Phonon absorption → not included (ω < E_gap region excluded)

---

## 8. Implementation

### 8.1 Synthetic data files

Two files are generated programmatically by
`utils/build_darkphoton_hypothetical_material.py` and placed under
`DarkELF/data/HypMat_0p34/`:

**`HypMat_0p34.yaml`** — minimal material YAML required by darkelf:
```yaml
rhoT:    8.0       # g/cm³ (A4)
E_gap:   0.34      # eV
e0:      1.6       # eV (ε_h, used only downstream in CCDarkSens)
unitcell: {'HypMat': {'A': 50.0, 'mult': 1}}  # placeholder A=50
atoms: ['HypMat']
```

**`HypMat_0p34_eps_electron_opticallimit.dat`** — 3-column optical file:
```
Hypothetical material: E_gap=0.34eV, sigma_DC=300/Ohm/cm, eps1=<bracket>
# omega[eV]   eps1   eps2
0.01          1.0    0.0        <- below gap: eps2 = 0
...
0.34          1.0    6.56       <- threshold
0.50          1.0    4.46
1.00          1.0    2.23
...
100.0         1.0    0.022
```

The eps₁ column is filled with the bracketing constant (1.0 or 12.0).
Two files are written — one per ε₁ assumption — with both registered
in `_DARKELF_FILES` as `"hypmat_unscreened"` and `"hypmat_screened"`.

### 8.2 Pipeline integration

After file generation, `darkphoton/entry.py::compute_dRdE` is called
identically to the Si case:

```python
res = compute_dRdE(
    material="hypmat_unscreened",  # or "hypmat_screened"
    mA_eV=mA_eV,
    epsilon=epsilon,
    darkelf_dir="/path/to/DarkELF",
)
```

darkelf routes through `load_eps_electron_opticallimit` → `R_absorption`
via the `eps_electron_opticallimit_loaded` path, handling all prefactor
and unit conversion internally. The output is the same boxcar-format CSV
that feeds directly into `DarkPhotonModel::Configure` → `ChargeIonization::FoldToNe`
→ PLR scan with no C++ changes.

### 8.3 Ionization yield

For the C++ scan config, the hypothetical material's ionization parameters
are set in `response.charge_ionization`:

```json
"charge_ionization": {
  "table_csv": "data/p100K_gap0p34_eh1p6.csv",
  "band_gap_eV": 0.34,
  "eh_pair_eV": 1.6,
  "scenario": "HypMat-0p34"
}
```

The p100K table is built via the existing `utils/build_p100K_scaled.py`
infrastructure — the same rescaling formula already used for the DM-electron
band-gap pheno study:

$$
E'(E) = E_\text{gap}^\text{ref} + \bigl(E - 0.34\bigr) \times \frac{3.8}{1.6}
\quad\text{for } E \geq 0.34\;\text{eV}
$$

---

## 9. Result Interpretation and Caveats

The output of this calculation is a **bracket on the κ vs m_A' exclusion
curve** for this hypothetical material, bounded by two ε₁ assumptions:

- **ε₁ = 1 (no screening):** upper bound on absorption rate → **strongest
  (most optimistic) limit** on κ
- **ε₁ = 12 (Si-like screening):** suppressed rate → **weaker (conservative)
  limit** on κ, closer to real material expectation

The width of this bracket (factor of ε₁ = 12 on κ) represents the dominant
systematic uncertainty. The frequency-independence assumption (A1) primarily
affects the shape at m_A' ≫ E_gap = 0.34 eV — likely an underestimate
there. The density uncertainty (A4) is subdominant.

**The result should be reported as:** "projected sensitivity for a hypothetical
material with E_gap=0.34 eV, σ_DC=300 Ω⁻¹cm⁻¹, assuming frequency-independent
conductivity and static dielectric constant in the range ε₁=1–12."

This is the same spirit as the DM-electron band-gap pheno study (B-thresh /
D-equal brackets), and carries the same methodological caveat: not a prediction
for any specific real material, but a demonstration of how sensitivity scales
with material parameters.

---

## 10. Connection to the Band-Gap Pheno Study

| Feature | DM-electron band-gap pheno | Dark photon absorption (this doc) |
|---|---|---|
| Rate calculation | QCDark2 scissor-shifted DFT | darkelf optical-limit absorption |
| Material input (rates) | Full ε(ω,k) from DFT | Synthetic ε₁,ε₂(ω) from σ_DC |
| Material input (ionization) | Rescaled p100K table | Same rescaled p100K table |
| Known unknowns | eh_pair_eV has no fundamental relation to E_gap | ε₁(ω) unknown; σ(ω) extrapolated from DC |
| Bracketing strategy | B-thresh vs D-equal | ε₁=1 vs ε₁=12 |
| What is claimed | Sensitivity bracket, not real material prediction | Same |

Both studies are explicit bracketing exercises under transparent, documented
assumptions. Neither claims to predict the sensitivity of a specific real
material without further experimental input.
