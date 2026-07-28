# Hypothetical Low-Gap Material — Unified Signal Model Reference
## Dark Photon Absorption and DM-Electron Scattering: Full Derivations, Assumptions, and Pipeline

**Document:** `docs/HypMat_Signal_Models_Unified.md`
**Author:** Diego Venegas-Vargas
**Date:** 2026-06-30
**Status:** Phenomenological study — two independent signal channels, explicit assumptions labeled

---

## 0. Scope

This document is the single reference for **both** signal models studied for a
hypothetical low-gap target material:

1. **Dark photon (hidden photon) absorption** — a kinetically-mixed dark photon
   A' of mass m_A' is absorbed by the material, depositing energy m_A' as a
   monochromatic ionization event.
2. **DM-electron scattering** — a sub-GeV dark matter particle χ scatters off
   a bound electron, depositing a continuous recoil-energy spectrum dR/dE.

Both studies target the same hypothetical material:

| Property | Symbol | Value | Units |
|---|---|---|---|
| Band gap | E_gap | 0.34 | eV |
| Mean e-h pair energy | ε_h | 1.6 | eV |
| DC optical conductivity | σ₁(DC) | 300 | Ω⁻¹cm⁻¹ |
| Mass density | ρ_T | 8 | g/cm³ |
| Target mass | M | 1 | kg |
| Livetime | T | 1 | yr |
| Dark current | DC | 10⁻⁵ | e⁻/pix/day |

The document gives full step-by-step derivations of each signal rate, then
explicitly catalogs what each calculation assumes and what it misses.

---

## Part I — Dark Photon Absorption

### I.1 Physical process

A kinetically-mixed hidden photon A' of mass m_A' and kinetic-mixing parameter
κ couples to ordinary photons via:

$$
\mathcal{L} \supset -\frac{\kappa}{2}\, F_{\mu\nu}^{A'} F^{\mu\nu}_\gamma
$$

When a dark photon passes through a material it can be resonantly absorbed,
depositing its full rest energy as a single electron-hole pair creation event.
This is the condensed-matter analogue of photoelectric absorption: a photon of
energy ω = m_A' is absorbed by an electron in the valence band, promoting it
to the conduction band. The signal is therefore **monochromatic** at

$$
E_\text{deposit} = m_{A'} \qquad \text{(eV)}
$$

This distinguishes dark photon absorption from DM-electron scattering (Part II),
which produces a continuous differential rate dR/dE.

### I.2 In-medium self-energy

In a material, the photon propagator is modified by the medium's dielectric
response. The dark photon acquires an in-medium effective mass through the
transverse photon self-energy Π_T(ω, k):

$$
\Pi_T(\omega, k) = \omega^2 \left[1 - \varepsilon(\omega, k)\right]
$$

In the non-relativistic, long-wavelength limit appropriate for m_A' in the
eV–keV range (where k ≈ m_A' v_DM → 0 for v_DM ~ 10⁻³ c), the dielectric
function approaches its **optical limit** k → 0:

$$
\varepsilon(\omega) \equiv \varepsilon(\omega, k=0) = \varepsilon_1(\omega) + i\,\varepsilon_2(\omega)
$$

### I.3 Absorption rate: step-by-step derivation

**Starting point:** the imaginary part of the dark photon self-energy determines
its decay width into the medium. The total absorption rate per unit target mass
is obtained by integrating the dark photon flux against the medium's response
(Hochberg et al. 2017; Knapen, Kozaczuk & Lin 2021):

**Step 1 — Dark photon flux.**
The local DM number density is n_χ = ρ_χ / m_A', with ρ_χ = 0.3 GeV/cm³.
In the rest frame of the detector, dark photons arrive with velocity v_DM ~ 10⁻³
and flux Φ = n_χ v_DM.

**Step 2 — Coupling of A' to photon.**
After kinetic mixing diagonalization, the A' couples to charged matter with
effective coupling κ e. In the medium, the A' propagator carries a self-energy
whose imaginary part (from Im[ε]) gives the probability per unit time of
absorption.

**Step 3 — Rate per target electron.**
The absorption rate per target electron at photon energy ω = m_A' is:

$$
\Gamma_\text{abs}(\omega) = \frac{\kappa^2\,\omega}{n_e}\,\text{Im}\!\left[-\frac{1}{\varepsilon(\omega)}\right]
\cdot \frac{\rho_\chi}{m_{A'}}
$$

where n_e is the electron number density in the material.

**Step 4 — Rate per unit target mass.**
Dividing by the target mass density ρ_T and integrating the monochromatic signal
(the δ-function collapses the ω integral):

$$
\boxed{
R_\text{abs}(m_{A'}) = \frac{\kappa^2\,\rho_\chi}{\rho_T\,m_{A'}}
\cdot \text{Im}\!\left[-\frac{1}{\varepsilon(m_{A'})}\right]
\cdot \mathcal{N}
\quad \left[\text{events}\,\text{kg}^{-1}\,\text{yr}^{-1}\right]
}
$$

where 𝒩 converts natural units (eV) to events/kg/yr using ℏ, c, and the
electron charge.

**Step 5 — Energy-loss function (ELF).**
The imaginary part of the inverse dielectric function is called the
**energy-loss function** (ELF):

$$
\text{ELF}(\omega) \equiv \text{Im}\!\left[-\frac{1}{\varepsilon(\omega)}\right]
= \frac{\varepsilon_2(\omega)}{\varepsilon_1^2(\omega) + \varepsilon_2^2(\omega)}
$$

This quantity is dimensionless and physically measurable via electron energy-loss
spectroscopy (EELS) or optical ellipsometry. The absorption rate is therefore:

$$
R_\text{abs}(m_{A'}) \propto \frac{\kappa^2}{m_{A'}} \cdot \text{ELF}(m_{A'})
$$

**Step 6 — Upper limit on κ.**
Given an observed event count N_obs and a Poisson 90% CL upper limit N_90% (e.g.
N_90% = 2.44 for N_obs = 0), the exclusion limit is:

$$
\kappa^{90\%}(m_{A'}) = \sqrt{\frac{N_{90\%}}{R_\text{abs}(\kappa=1,\,m_{A'})
\;\times\; M \cdot T}}
$$

where M·T is the exposure in kg·yr.

### I.4 Connecting σ₁ to the dielectric function

For a material characterized only by its DC optical conductivity σ₁, we
construct ε(ω) from the fundamental relation (SI units):

$$
\varepsilon(\omega) = \varepsilon_\infty + \frac{i\,\sigma(\omega)}{\varepsilon_0\,\omega}
$$

Taking real and imaginary parts:

$$
\varepsilon_2(\omega) = \frac{\sigma_1(\omega)}{\varepsilon_0\,\omega}
\qquad
\varepsilon_1(\omega) = \varepsilon_\infty - \frac{\sigma_2(\omega)}{\varepsilon_0\,\omega}
$$

**Numerical evaluation** with σ₁ = 300 Ω⁻¹cm⁻¹ = 3×10⁴ S/m:

$$
\varepsilon_2(\omega) = \frac{3\times10^4}{8.854\times10^{-12} \times \omega[\text{eV}]
\times 1.5193\times10^{15}} = \frac{2.23}{\omega[\text{eV}]}
$$

| m_A' [eV] | ε₂ | ELF (ε₁=1) | ELF (ε₁=12) |
|---|---|---|---|
| 0.34 | 6.56 | 0.148 | 0.046 |
| 1.0  | 2.23 | 0.910 | 0.015 |
| 5.0  | 0.45 | 0.41  | 0.0031 |
| 10.0 | 0.22 | 0.22  | 0.0015 |

Note: at m_A' near 1 eV, ε₂ ≈ 2.23 > ε₁ = 1, so the weak-absorption
approximation ELF ≈ ε₂/ε₁² fails there. The exact formula is always used
in the code.

### I.5 Ionization yield: from E_deposit to n_e

Once the deposited energy E_dep = m_A' is known, the mean number of
electron-hole pairs created is:

$$
\langle n_e \rangle = \frac{E_\text{dep} - E_\text{gap}}{\varepsilon_h}
\qquad \text{for } E_\text{dep} \geq E_\text{gap}
$$

The actual n_e is not sharp but follows a Fano-broadened distribution P(n_e|E).
This is encoded in the ionization table `p100K_gap0p34_eh1p6.csv`, derived via
the anchored energy map from the Si reference table (see §III.3).

### I.6 Bracketing strategy for ε₁

Since ε₁(ω) is unknown for this hypothetical material, two extreme scenarios
bracket the true answer:

| Scenario | ε₁ | Interpretation | Rate scaling |
|---|---|---|---|
| Unscreened | 1 | No dielectric screening (vacuum) | Maximum rate (optimistic κ limit) |
| Si-like screened | 12 | Typical semiconductor n ≈ 3.5 | Minimum rate (conservative κ limit) |

The rate scales as 1/ε₁² in the weak-absorption limit, so the spread in κ
between the two brackets is a factor of ε₁ = 12 — the dominant systematic
uncertainty in this calculation.

---

## Part II — DM-Electron Scattering

### II.1 Physical process

A sub-GeV DM particle χ scatters off a bound electron in the target via a
mediator. We consider two standard mediator scenarios:

- **Light mediator** — massless or ultra-light dark photon: mediator propagator
  ~ 1/q², rate ∝ F_DM(q)² = (α m_e / q²)².
- **Heavy mediator** — mediator mass ≫ q_typical: propagator ~ 1/m_med², rate
  independent of q.

The DM-electron interaction Lagrangian (benchmark: dark photon mediator) is:

$$
\mathcal{L} \supset -g_\chi\,\bar{\chi}\gamma^\mu\chi\,A'_\mu
- e\,Q_e\,\bar{e}\gamma^\mu e\,A_\mu
- \frac{\kappa}{2}\,F^{A'}_{\mu\nu}F^{\mu\nu}
$$

with coupling constant σ_e (the DM-electron cross section at fixed reference
momentum q_0 = α m_e) as the scan parameter.

### II.2 Reference cross section

The DM-free-electron cross section at reference momentum q₀ = α m_e is
defined as:

$$
\bar{\sigma}_e = \frac{\mu_{\chi e}^2\, g_\chi^2\, g_e^2}{\pi\, m_{A'}^4}
\bigg|_{q=q_0}
$$

where μ_χe = m_χ m_e / (m_χ + m_e) is the DM-electron reduced mass.

The full cross section at momentum transfer q and DM mass m_χ is:

$$
\frac{d\sigma}{d\ln q} = \bar{\sigma}_e\, F_\text{DM}^2(q)\, \left(\frac{q_0}{q}\right)^2
$$

with:
- Heavy mediator: F_DM(q) = 1
- Light mediator: F_DM(q) = (q₀/q)²

### II.3 Rate derivation: step-by-step

**Step 1 — Transition matrix element.**
The probability that DM of mass m_χ, incoming velocity v, deposits energy ω
and momentum **q** in the crystal is determined by the matrix element:

$$
|\mathcal{M}|^2 \propto \bar{\sigma}_e\, F_\text{DM}^2(q)
\cdot |f_\text{crystal}(\mathbf{q}, \omega)|^2
$$

where f_crystal(q, ω) is the crystal form factor encoding the electronic
structure of the target.

**Step 2 — Crystal form factor.**
The crystal form factor is related to the imaginary part of the inverse
dielectric function at **finite momentum transfer** q:

$$
|f_\text{crystal}(\mathbf{q}, \omega)|^2 \propto \frac{1}{V_\text{cell}}
\sum_{i \to f} |\langle f | e^{i\mathbf{q}\cdot\mathbf{r}} | i \rangle|^2
\delta(E_f - E_i - \omega)
$$

This is exactly Im[-1/ε(q, ω)] in the RPA, which is what QCDark2 computes
via the dielectric function HDF5.

**Step 3 — DM velocity average.**
The DM velocity distribution in the detector frame is the Standard Halo
Model (SHM), a Maxwell-Boltzmann truncated at v_esc:

$$
f(\mathbf{v}) = \frac{1}{N_\text{esc}} \exp\!\left(-\frac{|\mathbf{v}+\mathbf{v}_E|^2}{v_0^2}\right)
\Theta(v_\text{esc} - |\mathbf{v}+\mathbf{v}_E|)
$$

with v₀ = 238 km/s, v_E = 263 km/s (Earth velocity), v_esc = 544 km/s.

**Step 4 — Kinematic constraint.**
For DM of mass m_χ and speed v, the minimum DM speed to deposit (q, ω) is:

$$
v_\text{min}(q, \omega) = \frac{\omega}{q} + \frac{q}{2m_\chi}
$$

This sets the lower bound of the velocity integral.

**Step 5 — Differential rate.**
Combining steps 1–4, the differential scattering rate per unit mass per unit
recoil energy is:

$$
\frac{dR}{dE}\bigg|_E = \frac{\rho_\chi}{\rho_T\, m_\chi}\,\bar{\sigma}_e
\int \frac{d^3q}{(2\pi)^3}\,
F_\text{DM}^2(q)\,\frac{1}{q}\,\frac{2\pi^2}{\omega}
\,\text{Im}\!\left[-\frac{1}{\varepsilon(q,\omega)}\right]
\cdot g(v_\text{min})
$$

where:

$$
g(v_\text{min}) = \int_{v > v_\text{min}} \frac{f(\mathbf{v})}{v}\,d^3v
$$

is the mean inverse speed integrated above threshold.

**Step 6 — Numerical evaluation.**
QCDark2 evaluates this integral numerically using the precomputed ε(q, ω)
HDF5 on a grid of (q, ω) values, summing contributions from all allowed
DM velocities. The output is dR/dE in units of events/(kg·yr·eV), stored
as a 2-column CSV indexed by (m_χ, σ_e).

**Step 7 — Upper limit on σ_e.**
Given observed counts and background, the profile likelihood ratio (PLR) scan
over σ_e at each mass m_χ yields:

$$
\sigma_e^{90\%}(m_\chi) = \text{PLR}^{-1}_{90\%}\!\left[
\frac{dR}{dE}(\bar{\sigma}_e=1),\; M\cdot T,\; \text{background},\; \text{efficiency}
\right]
$$

### II.4 The scissors approximation

Since a full DFT calculation of ε(q, ω) for the hypothetical material is not
available, we use a **scissors-corrected Si dielectric function** as a proxy.

**What the scissors operator does (step by step):**

1. Run DFT on bulk Si (fixed lattice, PBE functional, k-grid 4×4×4).
   This gives occupied/virtual molecular orbital (MO) energies and coefficients.

2. Set HOMO = 0 eV. Measure the DFT gap:
   Δ_DFT = ε_LUMO − ε_HOMO (typically 0.6–1.0 eV for PBE/Si).

3. Apply scissor: shift **all conduction-band energies** by a rigid correction:
   $$
   \varepsilon_n^{(\text{scissor})} = \varepsilon_n^{(\text{DFT})} + (\text{target gap} - \Delta_\text{DFT})
   \qquad \forall n \in \text{conduction}
   $$

4. Recompute Im ε(ω, q) in the RPA using the shifted MO energies but
   **unchanged MO coefficients (wavefunctions)**. This modifies the transition
   energies (where peaks appear in ε₂) without changing the transition matrix
   elements (how strong each peak is).

5. Package the resulting ε(q, ω) into `data/qcdark2_epsilon/Si/Si_fast_gap0p34.h5`.

**What changes vs what does not:**

| Quantity | With scissor gap = 0.34 eV |
|---|---|
| Lattice vectors, atom positions | Unchanged |
| DFT Hamiltonian, SCF density | Unchanged (reused from cache) |
| MO **coefficients** (wavefunctions) | Unchanged |
| MO **energies** (conduction bands) | Rigidly shifted down by ~(1.1 − 0.34) eV |
| ε(ω, q), DM rates, limits | Changed — threshold at 0.34 eV, more phase space |

### II.5 From ε to DM rate: key physical consequence of the scissors

With the conduction band threshold at 0.34 eV:
- DM masses as low as m_χ ≈ m_e × v_max / Δv_min ≈ a few MeV can now
  deposit ω ≥ E_gap.
- The velocity integral g(v_min) opens up for lower m_χ compared to the
  Si reference (gap = 1.2 eV).
- The spectral shape of dR/dE near threshold reflects the phase space opening,
  which in turn reflects the Si conduction band density of states — a Si-like
  band structure, not the true hypmat band structure.

---

## Part III — Shared Downstream Chain

Both signal models produce a dR/dE spectrum in CSV format. From that point the
**identical** downstream pipeline is used:

```
dR/dE(E)   [events / kg / yr / eV]
    ↓
× exposure (M·T) × dE  →  events per energy bin
    ↓
× P(n_e | E)  →  spread into n_e integer bins     [ChargeIonization::FoldToNe]
    ↓
pattern / diffusion / readout response             [EfficiencyMC]
    ↓
profile likelihood ratio scan over coupling        [PLR]
    ↓
upper limit on κ (dark photon) or σ_e (DM-e)
```

### III.1 Monochromatic vs continuous input

| Model | Input to pipeline |
|---|---|
| Dark photon | Boxcar dR/dE: zero except in [m_A' ± δ/2] where height = R_abs/δ |
| DM-electron | Continuous dR/dE(E) from QCDark2 rate table |

The boxcar representation ensures that ChargeIonization::FoldToNe integrates
the correct total rate when folded against P(n_e|E). The bin width δ (set by
`binsize_eV` in the config) must match downstream energy binning.

### III.2 Observable space

Both studies use **n_e bins** (not pattern bins) for the hypothetical material:

```json
"experiment": { "observable_bins": "ne", "roi_bins": [1,2,3,4,5,6,7] },
"response":   { "analysis_space":  "ne" }
```

This covers ionization events up to 7 e⁻, corresponding to recoil energies up
to approximately 7 × ε_h + E_gap = 7 × 1.6 + 0.34 ≈ 11.5 eV.

### III.3 Ionization table: the anchored pheno map

The function P(n_e | E) is not computed from first principles for the
hypothetical material. Instead it is rescaled from the Si reference table
`data/p100K_table.csv` (Si: E_gap = 1.2 eV, ε_h = 3.8 eV) via the
**anchored energy map**:

$$
E'(E) = E_\text{gap}^\text{ref} + \bigl(E - E_\text{gap}^\text{new}\bigr)
\times \frac{\varepsilon_h^\text{ref}}{\varepsilon_h^\text{new}}
\qquad \text{for } E \geq E_\text{gap}^\text{new}
$$

Numerically for (E_gap_new = 0.34 eV, ε_h_new = 1.6 eV):

$$
E'(E) = 1.2 + (E - 0.34) \times \frac{3.8}{1.6} = 1.2 + (E - 0.34) \times 2.375
$$

**Physical interpretation:** the map asks "at what energy does bulk Si have the
same fractional excitation above threshold as the new material at energy E?"
The Si P(n_e|E') table is then used directly, anchoring the ionization
probability to the Si reference at the same relative excitation level.

The resulting table `data/p100K_gap0p34_eh1p6.csv` has:
- Zero ionization probability for E < 0.34 eV
- Mean n_e = 1 at E ≈ 0.34 + 1.6 = 1.94 eV
- Mean n_e = k at E ≈ 0.34 + k × 1.6 eV

### III.4 Efficiency and background

Both studies use the same efficiency CSV derived from the Si LBC measurement:

```
data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv
```

Background model: `dc_flat_migration` with dark current scaled to:

```
λ_e = DC × 365.25 = 10⁻⁵ × 365.25 = 3.65×10⁻³  e⁻/pix/yr
```

This is purely a rate scaling of the Si DC model. It does not account for
material-specific sources of background, carrier generation rates, or
ionization yield from radioactive contaminants in the hypothetical material.

---

## Part IV — Complete Assumption Catalog

### IV.1 Dark photon absorption assumptions

| Label | Assumption | Impact | How to fix |
|---|---|---|---|
| **A1** | σ₁(ω) = σ₁(DC) = const for ω ≥ E_gap | Underestimates rate at m_A' ≫ E_gap (real σ₁ rises with interband transitions) | Measure σ₁(ω) spectrum via ellipsometry |
| **A2** | ε₁(ω) bracketed by {1, 12} | Factor of ε₁² = 144 uncertainty on rate; factor of 12 on κ | Measure static n (refractive index); gives ε₁ = n² |
| **A3** | Exact ELF formula used (not weak-absorption approx) | None — full formula is correct; approximation is only noted to flag where ε₂ > ε₁ | N/A |
| **A4** | ρ_T = 8 g/cm³ | Linear in rate → √2 on κ limit if ρ_T off by 2× | Crystal density measurement or DFT lattice |
| **A5** | ε_h enters only ionization, not rate | Correct by construction — ε_h is a charge yield parameter, not relevant to EM absorption | N/A |
| **A6** | Efficiency borrowed from Si | Unknown error in signal acceptance | Re-simulate diffusion for hypmat geometry |

### IV.2 DM-electron scattering assumptions

| Label | Assumption | Impact | How to fix |
|---|---|---|---|
| **B1** | Band topology is Si (scissors only shifts gap) | Wrong spectral shape of dR/dE — MO coefficients (matrix elements) are Si's, not hypmat's | Full DFT+Wannier for hypmat |
| **B2** | Scissors shifts conduction bands rigidly | Band dispersion (effective masses), Fermi velocity, screening at finite q remain Si-like | DFT of actual material |
| **B3** | k-grid is 4×4×4 (fast) | ε(q, ω) less accurate than 8×8×8 production; limit shape noisier | Rerun with `Si_lfe8q.in` at 8×8×8 |
| **B4** | No LFE at this level (Si_fast_scissor.in) | Local field effects missing at small q; affects low-mass DM most | Use LFE template for production |
| **B5** | Density ρ_T = 8 g/cm³ vs Si 2.33 g/cm³ | Linear prefactor in rate — need to verify it's applied in config, not absorbed into QCDark2 Si normalization | Check `density_g_cm3` config field |
| **B6** | Valence electron count = Si (4 per cell) | Rate normalization: if hypmat has different Z_val per cell, number of targets per kg differs | Specify unit cell composition |
| **B7** | Ionization table rescaled from Si | Shape of P(n_e|E) is Si-like (Fano factor, sub-gap phonon tail) | DFT phonon + ionization model for hypmat |
| **B8** | Efficiency borrowed from Si | Same as A6 | Same |
| **B9** | SHM halo model | Standard assumption, not material-specific | Use alternative halo models for systematics |

### IV.3 Assumptions common to both

| Assumption | Impact |
|---|---|
| E_gap = 0.34 eV — threshold for both rate and ionization | If wrong, shifts minimum accessible DM mass |
| ε_h = 1.6 eV — mean charge yield scale | If wrong, shifts n_e distribution: fewer/more electrons per event |
| DC = 10⁻⁵ e/pix/day — background level | If higher, weakens limit (more background events in n_e=1 bin) |
| Exposure M·T = 1 kg·yr — linear in rate limit | N/A — purely a design parameter |

---

## Part V — Comparative Summary

### V.1 Side-by-side

| Feature | Dark photon absorption | DM-e scattering |
|---|---|---|
| Signal shape | Monochromatic line at E = m_A' | Continuous dR/dE(E) |
| Rate formula | ∝ κ² × ELF(m_A') / m_A' | Integral over ε(q,ω), DM velocity |
| Material input (rate) | Synthetic ε₁, ε₂(ω) from σ_DC — two-bracket approach | QCDark2 scissors ε(q,ω) with Si wavefunctions |
| Material input (ionization) | Rescaled p100K table (same for both) | Rescaled p100K table (same for both) |
| q-dependence | None (optical limit k→0 exact) | Full finite-q ε(q,ω) from QCDark2 |
| Dominant unknown | ε₁(ω) — factor 12 on κ | Band structure shape (B1) — spectral shape of dR/dE |
| Bracketing strategy | ε₁ = {1, 12} (unscreened / Si-like) | 4×4×4 fast ↔ future 8×8×8 production |
| Pheno validity | Strong: optical limit is exact; only σ₁(ω) extrapolated | Moderate: scissors gap is right; wavefunction character is wrong |
| What is claimed | Sensitivity bracket for material with these σ_DC, E_gap, ε_h | Sensitivity projection with Si-proxy form factor at correct threshold |

### V.2 Which channel is more reliable?

For this hypothetical material:

- **Dark photon is the cleaner case.** The rate formula is exact given ε(ω).
  The only inputs are σ₁(DC) (measured or specified) and ε₁ (bracketed).
  The two-bracket strategy gives a rigorous upper and lower bound on sensitivity
  that can only be narrowed with a refractive index measurement.

- **DM-electron is more approximate.** The scissors ε(q, ω) is physically
  motivated and standard in the literature, but the spectral shape of dR/dE
  reflects Si band structure. The result is a reasonable estimate of the
  sensitivity trend (how it scales with mass, how it compares to Si), but not
  a quantitative prediction for a specific real material.

Both channels use the same ionization table, efficiency model, and PLR
infrastructure — systematic differences between them are therefore confined
to the rate calculation only.

### V.3 Honest reporting language

**Dark photon:**
> "Projected exclusion bracket on kinetic mixing κ for a material with
> E_gap = 0.34 eV, σ_DC = 300 Ω⁻¹cm⁻¹, assuming frequency-independent
> conductivity and static dielectric constant in the range ε₁ = 1–12."

**DM-electron:**
> "Projected sensitivity on DM-electron cross section σ_e for a material
> with ionization threshold E_gap = 0.34 eV and mean pair energy ε_h = 1.6 eV,
> using a scissors-corrected Si dielectric function (QCDark2) as a proxy for
> the target crystal form factor."

---

## Part VI — Implementation Reference

### VI.1 Dark photon: files and commands

| Step | Command |
|---|---|
| Build optical file | `python3 utils/build_darkphoton_hypothetical_material.py --darkelf_dir /path/to/DarkELF` |
| Generate rate CSVs (unscreened) | `python3 utils/darkphoton_generate_grid.py configs/darkphoton_generate_hypmat_unscreened.json` |
| Generate rate CSVs (screened) | `python3 utils/darkphoton_generate_grid.py configs/darkphoton_generate_hypmat_screened.json` |
| Run scan (unscreened) | `./build/ccdarksens_scan_dmelectron_pattern configs/darkphoton_scan_hypmat_unscreened.json` |
| Run scan (screened) | `./build/ccdarksens_scan_dmelectron_pattern configs/darkphoton_scan_hypmat_screened.json` |
| Plot | `./build/ccdarksens_plot_dmelectron_limit <ROOT files> --dark-photon --from-qhist --batch` |

Key config files: [`configs/darkphoton_generate_hypmat_unscreened.json`](../configs/darkphoton_generate_hypmat_unscreened.json), [`configs/darkphoton_scan_hypmat_unscreened.json`](../configs/darkphoton_scan_hypmat_unscreened.json)

### VI.2 DM-electron: files and commands

| Step | Command |
|---|---|
| Build ionization table | `python3 utils/build_p100K_scaled.py --band-gap-eV 0.34 --eh-pair-eV 1.6 --scenario "D-equal" --E-min 0.05` |
| Build ε HDF5 (scissors 0.34 eV) | `python3 utils/qcdark2_regenerate_epsilon.py --template configs/qcdark2/Si_fast_scissor.in --scissor 0.34` |
| Generate rate CSVs | `python3 utils/qcdark2_generate_grid.py configs/qcdark2_generate_hypmat_gap0p34_light.json` |
| Run scan | `./build/ccdarksens_scan_dmelectron_pattern configs/scan_hypmat_gap0p34_light.json` |
| Plot | `./build/ccdarksens_plot_dmelectron_limit <ROOT files> --from-qhist --batch` |

Key output files:
- `data/p100K_gap0p34_eh1p6.csv` — ionization table
- `data/qcdark2_epsilon/Si/Si_fast_gap0p34.h5` — dielectric function HDF5
- `data/qcdark2_rates/Si/light/Si_fast_gap0p34/` — rate CSV grid (to be generated)

### VI.3 Config parameters for both channels

```json
// Shared detector block
"detector": {
  "mass_kg": 1.0,
  "density_g_cm3": 8.0,
  "livetime_days": 365.25
},

// Shared ionization block
"charge_ionization": {
  "table_csv":    "data/p100K_gap0p34_eh1p6.csv",
  "band_gap_eV":  0.34,
  "eh_pair_eV":   1.6,
  "scenario":     "D-equal"
},

// Shared background block
"run": {
  "background_source": "dc_flat_migration",
  "background_model":  "scale"
},
"experiment": {
  "lambda_e_per_pix_per_year": 3.65e-3,
  "observable_bins": "ne",
  "roi_bins": [1,2,3,4,5,6,7]
}
```

---

## Part VII — Connection to the Band-Gap Pheno Study

The hypothetical material studies documented here are an extension of the
broader band-gap phenomenological study (`docs/band_gap_pheno_ionization.md`,
`docs/qcdark2_dielectric_workflow.md`), which scans the sensitivity of
DM-electron scattering across multiple (E_gap, ε_h) combinations using the
same Si scissors approach.

The key distinction:

| Band-gap pheno study | Hypmat study (this doc) |
|---|---|
| Scans E_gap ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 1.2} eV | Fixed E_gap = 0.34 eV |
| Both heavy and light mediator | Both channels (DM-e + dark photon) |
| Multiple ε_h scenarios (B-thresh, D-equal, Klein) | Single (E_gap, ε_h) = (0.34, 1.6) eV |
| No dark photon channel | Dark photon with two ε₁ brackets |
| Si efficiency/background throughout | Same, extended to n_e = 7 |
| Purpose: scan how sensitivity scales with material properties | Purpose: full signal model for a specific hypothetical target |

The hypmat study therefore represents a **fixed point** within the band-gap
pheno parameter space, enriched with the additional dark photon absorption
channel and an explicit two-bracket systematic for ε₁.
