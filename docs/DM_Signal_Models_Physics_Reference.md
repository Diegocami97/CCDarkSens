# Dark Matter Signal Models — Physics Reference
## QEDark/QCDark2 DM–Electron Scattering, Dark Photon Absorption, and the Migdal Effect
### Standard Silicon and SrCd₂Sb₂ Narrow-Gap Projections

**Document:** `docs/DM_Signal_Models_Physics_Reference.md`
**Author:** Diego Venegas-Vargas
**Status:** Consolidated physics reference — derivations, assumptions, and pipeline connections for all three signal channels implemented in CCDarkSens

---

## 0. Scope and how this document fits with the rest of `docs/`

CCDarkSens is a **general-purpose** dark matter sensitivity framework for
DAMIC-M: it computes profile-likelihood-ratio (PLR) upper limits over a
shared detector-response pipeline (n_e or pattern bins), and that pipeline
is agnostic to which physical process supplies the input deposited-energy
spectrum dR/dE. The framework's intended scope includes both
**electron-recoil-like** searches (observable ionization produced directly
by the interaction, no separate nuclear-recoil energy scale) and
**nuclear-recoil-like** searches for canonical WIMP dark matter (elastic
DM–nucleus scattering, converted to ionization via a material quenching
factor). Each physical process plugs into the pipeline as its own
rate/signal module.

**This document covers only the electron-recoil-like channels currently
implemented** — three independent signal physics channels that share the
pipeline described in Part V:

| Channel | DM couples to | Rate engine | y-axis of the limit |
|---|---|---|---|
| DM–electron scattering | bound electron | QEDark / QCDark / QCDark2 | σ_e [cm²] |
| Dark photon absorption | bound electron (EM) | DarkELF `absorption` module | κ (kinetic mixing) |
| Migdal effect | nucleus | DarkELF `Migdal` module | σ_n [cm²] (per nucleon) |

**Explicitly out of scope:** elastic WIMP–nucleus scattering (canonical
spin-independent/spin-dependent nuclear recoil searches). Unlike the Migdal
effect — where the nuclear recoil is an intermediate step that still
terminates in an *electronic* shake-off signal, and which is therefore
covered here — a direct WIMP search would use the primary nuclear recoil
energy itself, converted to ionization via a quenching-factor model
distinct from Part IV's charge-ionization treatment. This is a planned
extension of CCDarkSens, not yet implemented; it would enter the same PLR
pipeline as its own rate module once added, and will be documented
separately when it exists.

This document derives the rate formulas for the three implemented channels
from first principles, states the standard **silicon** case for each, and
then states the **SrCd₂Sb₂** (internal codename HypMat / Sr2Cb2Sd)
narrow-gap projection case for each, with explicit, itemized assumptions.
It also derives the **charge-ionization model** (Klein's formula, the
P(n_e|E) table, and the "D-equal" scenario convention) that converts any of
the three deposited-energy spectra into the observable n_e / pattern
spectrum, and shows how all three channels feed the identical downstream
PLR chain.

This document **consolidates and extends** several channel-specific documents
already in this repository rather than replacing them:

| Existing document | Relationship to this one |
|---|---|
| [`HypMat_Signal_Models_Unified.md`](HypMat_Signal_Models_Unified.md) | Full step-by-step derivation of dark photon + DM-e for HypMat with **old** parameters (ε_h=1.6 eV, ρ=8 g/cm³, ε₁ step-function). This document supersedes those parameter values with the corrected ones (ε_h=1.7778 eV via Klein, ρ=5.76 g/cm³, Drude ε(ω)). |
| [`DarkPhoton_Absorption_HypotheticalMaterial.md`](DarkPhoton_Absorption_HypotheticalMaterial.md) | Assumption-by-assumption catalog for the dark photon channel; the Drude upgrade (§A1-updated) is the model used here. |
| [`Migdal_Integration_Plan.md`](Migdal_Integration_Plan.md) | Full Migdal derivation and CCDarkSens integration status (Si only; SrCd₂Sb₂ explicitly out of scope — see §3.5 below for why). |
| [`band_gap_pheno_ionization.md`](band_gap_pheno_ionization.md) | Defines the ionization-table rescaling scenarios (B-thresh, D-equal, A-ratio) used in §4. |
| [`LBC_QEDark_Reproduction_Guide.md`](LBC_QEDark_Reproduction_Guide.md) | Operational guide for reproducing the DAMIC-M PRL 2025 QEDark limit; §1 below gives the underlying physics. |

---

## 1. The common detector observable

Regardless of channel, the DM interaction deposits an energy ω (electronic
excitation energy, eV) in the target. This document derives dR/dω (or
dR/dE — the same quantity) for each channel; from that point the identical
chain applies (derived in full in Part V):

```
dR/dω(ω)  [events / kg / yr / eV]
    ↓  × exposure (M·T) × dω
events per energy bin
    ↓  × P(n_e | ω)                         [ChargeIonization::FoldToNe]
n_e spectrum
    ↓  pattern / diffusion / readout response [EfficiencyMC, PatternRates]
    ↓  profile likelihood ratio scan          [PLR]
upper limit on σ_e, κ, or σ_n
```

---

## Part I — DM–Electron Scattering (QEDark / QCDark / QCDark2)

### 1.1 Physical process and Lagrangian

A sub-GeV DM particle χ scatters off a bound valence electron in the target
crystal, promoting it to the conduction band and depositing a continuous
recoil-energy spectrum. The benchmark interaction (dark-photon mediator,
Essig, Mardon & Volansky 2012; Essig et al. 2016 "QEDark") is:

$$
\mathcal{L} \supset -g_\chi\,\bar{\chi}\gamma^\mu\chi\,A'_\mu
- e\,Q_e\,\bar{e}\gamma^\mu e\,A_\mu
- \frac{\kappa}{2}\,F^{A'}_{\mu\nu}F^{\mu\nu}
$$

After kinetic-mixing diagonalization this generates an effective
DM–electron vertex with coupling g_χ g_e κ / m²_med, where m_med is the
mediator mass. Two limits define the two "mediator cases" used throughout
CCDarkSens:

- **Heavy mediator** (m_med ≫ typical momentum transfer q ~ αm_e): the
  propagator is momentum-independent, F_DM(q) = 1.
- **Light (ultralight/massless) mediator** (m_med ≪ q): the propagator
  ∝ 1/q², giving F_DM(q) = (α m_e / q)².

### 1.2 Reference cross section and mediator form factor

The free-electron DM cross section, evaluated at the reference momentum
transfer q₀ = α m_e (the typical momentum transfer in atomic-scale
scattering), is:

$$
\bar{\sigma}_e = \frac{\mu_{\chi e}^2\, g_\chi^2\, g_e^2}{\pi\, m_{A'}^4}\bigg|_{q=q_0}
$$

where μ_χe = m_χ m_e/(m_χ + m_e) is the DM–electron reduced mass. This
single number, σ_e, is the scan parameter reported on the y-axis of the
final limit. The momentum-transfer dependence factors out into the
dimensionless mediator form factor:

$$
F_\text{DM}(q) = \begin{cases}
1 & \text{heavy } (n_\text{FDM}=0) \\
(\alpha m_e / q)^2 & \text{light } (n_\text{FDM}=2)
\end{cases}
$$

CCDarkSens's Python rate modules (`ccdarkphys.common.mediator_map.MEDIATOR_TO_FDM_INDEX`)
encode this exactly as the QEDark convention F_DM(q) ∝ q^(−n_FDM), n_FDM ∈ {0, 2}.

### 1.3 The crystal form factor

The electronic structure of the target enters through the **crystal form
factor** f_crystal(q, ω), which encodes the probability of a momentum-q,
energy-ω transition from an occupied valence Bloch state to an empty
conduction state:

$$
|f_\text{crystal}(\mathbf{q},\omega)|^2 \propto \frac{1}{V_\text{cell}}
\sum_{i\to f} |\langle f | e^{i\mathbf{q}\cdot\mathbf{r}} | i\rangle|^2\,
\delta(E_f - E_i - \omega)
$$

In the random-phase approximation (RPA) this is exactly proportional to the
imaginary part of the inverse dielectric function at finite momentum:

$$
|f_\text{crystal}(\mathbf{q},\omega)|^2 \;\propto\; \text{Im}\!\left[-\frac{1}{\varepsilon(\mathbf{q},\omega)}\right]
$$

CCDarkSens has **two independent implementations** of this quantity, both
consumed downstream by the identical C++ `RateTable`/`ccdarksens_scan_dmelectron_pattern`
pipeline:

| Engine | Source of |f_crystal|² or ε(q,ω) | Files |
|---|---|---|---|
| **QEDark** | Precomputed |f_crystal|² lookup table (`Si_f2.txt`, from Essig et al. 2016 DFT calculation on bulk Si) | `python/ccdarkphys/qedark/entry.py` |
| **QCDark2** | Full ε(q,ω) HDF5 computed by DFT+RPA (this collaboration's in-house dielectric-function code), evaluated on a (q, E) grid with local-field effects | `python/ccdarkphys/qcdark2/entry.py`, `data/qcdark2_epsilon/<material>/*.h5` |

QCDark2 is the more general engine: it accepts an arbitrary DFT-derived
ε(q,ω) HDF5 (attrs `M_cell`, `V_cell`, `dE`), so the same rate-generation
code path handles standard Si and any other material for which a dielectric
function file exists (§1.7).

### 1.4 DM velocity distribution and kinematics

The local DM population follows the truncated Standard Halo Model (SHM):

$$
f(\mathbf{v}) = \frac{1}{N_\text{esc}}\exp\!\left(-\frac{|\mathbf{v}+\mathbf{v}_E|^2}{v_0^2}\right)\Theta(v_\text{esc}-|\mathbf{v}+\mathbf{v}_E|)
$$

CCDarkSens's production configs (`configs/examples/qedark_generate_*.json`,
`qcdark2_generate_*.json`) use:

$$
v_0 = 238\ \text{km/s}, \qquad v_E = 263\ \text{km/s}, \qquad v_\text{esc} = 544\ \text{km/s}
$$

For a deposited (q, ω), the minimum DM speed able to supply that energy and
momentum transfer is:

$$
v_\text{min}(q,\omega) = \frac{\omega}{q} + \frac{q}{2m_\chi}
$$

and the halo integral reduces to the mean inverse speed above threshold:

$$
g(v_\text{min}) = \int_{v>v_\text{min}} \frac{f(\mathbf{v})}{v}\,d^3v
$$

### 1.5 The dR/dE master formula

Combining §1.2–1.4, the differential DM–electron scattering rate per unit
target mass per unit deposited energy is:

$$
\boxed{
\frac{dR}{dE}\bigg|_\omega = \frac{\rho_\chi}{\rho_T\, m_\chi}\,\bar{\sigma}_e
\int \frac{d^3q}{(2\pi)^3}\, F_\text{DM}^2(q)\,\frac{1}{q}\,\frac{2\pi^2}{\omega}\,
\text{Im}\!\left[-\frac{1}{\varepsilon(q,\omega)}\right]\, g(v_\text{min})
}
$$

with ρ_χ = 0.3 GeV/cm³ the local DM density and ρ_T the target mass density.
This is evaluated numerically by both QEDark and QCDark2, producing a CSV of
dR/dE [events/(kg·yr·eV)] versus E for each (m_χ, σ_e, mediator) grid point,
in exactly the two-column format read by `RateTable::LoadCSV` and rebinned
onto a uniform ROOT histogram by `RateTable::MakeTH1D` (`src/io/RateTable.cc`),
which **linearly interpolates the tabulated rate at each output bin's
center**. For DM-electron and Migdal spectra — smooth and continuous over
many bins — bin-center sampling is an excellent approximation. It is *not*
safe for the dark photon channel's monochromatic boxcar signal; see §2.6
for why, and for the numerical fix (`nbins`) required there.

**Threshold behavior:** ε(q,ω) (and hence dR/dE) is zero for ω below the
material's band gap E_gap — this is where the target's band structure sets
the absolute low-mass reach of the search.

### 1.6 Standard silicon case

| Parameter | Value | Note |
|---|---|---|
| Band gap E_gap | 1.2 eV (QEDark/scissor convention) or 1.11 eV (darkelf optical, 130 K) | Both values appear in this codebase for different purposes — see §4.1 caveat |
| Mean e-h pair energy ε_h | 3.8 eV | Measured, 100 K Si, PRD 102, 063026 |
| Target density ρ_T | 2.329 g/cm³ | Bulk Si |
| Valence electron count | Z_val = 4 (Si, 3s²3p²) | Sets the number of scattering targets per unit cell |
| Crystal form factor source | `Si_f2.txt` (QEDark) or `Si_comp.h5` (QCDark2, full DFT+RPA) | Both reproduce the same physics; QCDark2 includes local-field effects |
| σ_e reference q₀ | α m_e ≈ 3.7 keV | Standard convention |
| SHM (v₀, v_E, v_esc) | (238, 263, 544) km/s | This codebase's production default |
| Mediator cases | Heavy (n_FDM=0), light (n_FDM=2) | Both scanned in the DAMIC-M PRL |

This reproduces the DAMIC-M PRL 2025 result: freeze-in excluded for an
ultralight mediator over 3.5–490 MeV/c², freeze-out excluded for a heavy
mediator over 2.9–21.5 MeV/c² (see `docs/LBC_QEDark_Reproduction_Guide.md`
for the operational reproduction steps).

### 1.7 SrCd₂Sb₂ (HypMat) case

No first-principles DFT wavefunctions for SrCd₂Sb₂ are yet available
(AiiDA `bands_workchain_pk=9491`, pending). The DM-electron rate for this
material is therefore computed with the **scissors approximation**: a
scissor-shifted Si dielectric function stands in as a phenomenological proxy
for the true SrCd₂Sb₂ ε(q,ω).

**Procedure (scissors operator), step by step:**

1. Run DFT on bulk Si (PBE functional, 4×4×4 k-grid — the "fast" template
   `Si_fast_scissor.in`), obtaining Kohn–Sham orbital energies and
   coefficients.
2. Measure the DFT gap Δ_DFT = ε_LUMO − ε_HOMO (typically 0.6–1.0 eV
   for PBE/Si, an underestimate of the true 1.12 eV Si gap — the well-known
   DFT band-gap problem).
3. Rigidly shift all **conduction-band** energies:
   $$
   \varepsilon_n^{(\text{scissor})} = \varepsilon_n^{(\text{DFT})} +
   (\text{target gap} - \Delta_\text{DFT}) \qquad \forall\, n \in \text{conduction}
   $$
   with target gap = 0.34 eV (SrCd₂Sb₂ indirect gap, working value).
4. Recompute Im ε(q,ω) in the RPA using the shifted energies but
   **unchanged** wavefunction coefficients — this moves *where* the
   spectral weight sits (transition energies) without changing *how strong*
   each transition is (matrix elements, which remain Si's).
5. Package the result as `data/qcdark2_epsilon/Si/Si_fast_gap0p34.h5`, used
   identically to a native SrCd₂Sb₂ ε(q,ω) file by the rate generator.

| Parameter | Value | Note |
|---|---|---|
| Band gap (used for DM-e threshold) | 0.34 eV (indirect, working value; R2SCAN no-SOC gives indirect 0.560 eV / direct 0.695 eV) | Scissors target |
| Mean e-h pair energy ε_h | **1.7778 eV** (Klein's formula, §4.1) | Supersedes the earlier 1.6 eV placeholder |
| Target density ρ_T | **5.76 g/cm³** (from lattice a=4.44 Å, c=28.196 Å, MW=555.96 g/mol, Z=3) | Supersedes the earlier 8.0 g/cm³ placeholder |
| Valence electron count | Z_val = 16 per formula unit (Sr 2s + 2×Cd 2s + 2×Sb 5 valence) | Used for Drude ω_p in §2.5; DM-e rate keeps Si's 4/cell (assumption B6, unresolved) |
| Crystal form factor source | Scissor-shifted `Si_fast_gap0p34.h5` (Si wavefunctions, shifted energies) | **Band topology and matrix elements are Si's**, not SrCd₂Sb₂'s |
| k-grid | 4×4×4 ("fast") | Coarser than the 8×8×8 production template; noisier ε(q,ω) |
| Local field effects | Not included at this level | Most relevant at small q, i.e. lowest-mass DM |

**Explicit residual assumptions (full catalog in `HypMat_Signal_Models_Unified.md` §IV.2, labels B1–B9):**
the dominant one is **B1** — the spectral *shape* of dR/dE (where the peaks
in ε₂ sit, how fast the phase space opens above threshold) reflects Si's
band structure, not SrCd₂Sb₂'s. The scissors correction is trustworthy for
the **threshold location**; it is not trustworthy for the detailed shape.
This is a genuine sensitivity **projection using a proxy form factor**, not
a first-principles prediction, and must be reported as such.

---

## Part II — Dark Photon (Hidden Photon) Absorption

### 2.1 Physical process

A kinetically-mixed dark photon A′ of mass m_A′ couples to ordinary
electromagnetism via:

$$
\mathcal{L} \supset -\frac{\kappa}{2}\,F_{\mu\nu}^{A'}F^{\mu\nu}_\gamma
$$

If m_A′ ≥ E_gap, the dark photon can be **absorbed**, promoting a valence
electron to the conduction band and depositing its **full rest energy**
m_A′ as a single ionization event. Unlike DM–electron scattering, the
signal is **monochromatic**:

$$
E_\text{deposit} = m_{A'}
$$

### 2.2 In-medium self-energy and the energy-loss function

The dark photon's absorption rate in a medium is governed by the imaginary
part of its transverse self-energy, which in the non-relativistic,
long-wavelength (k → 0) limit relevant for m_A′ in the eV range reduces to
the material's optical-limit complex dielectric function ε(ω) = ε₁(ω) +
iε₂(ω). Following Hochberg et al. (2017) and Knapen, Kozaczuk & Lin (2021),
the absorption rate per unit target mass is:

$$
\boxed{
R_\text{abs}(m_{A'}) = \kappa^2\,\frac{\rho_\chi}{\rho_T\, m_{A'}}\,
\text{Im}\!\left[-\frac{1}{\varepsilon(m_{A'})}\right] \times \mathcal{N}
\quad\left[\text{events kg}^{-1}\,\text{yr}^{-1}\right]
}
$$

where 𝒩 converts natural units to events/kg/yr. The quantity

$$
\text{ELF}(\omega) \equiv \text{Im}\!\left[-\frac{1}{\varepsilon(\omega)}\right]
= \frac{\varepsilon_2(\omega)}{\varepsilon_1^2(\omega)+\varepsilon_2^2(\omega)}
$$

is the **energy-loss function** — dimensionless, measurable via EELS or
optical ellipsometry, and the sole material input to the rate once ρ_T is
fixed. The rate is therefore a **boxcar in energy**: zero except in a bin
of width δ centered at m_A′, with height R_abs/δ, which is exactly how
CCDarkSens represents it in the shared CSV format consumed by
`ChargeIonization::FoldToNe` (Part V).

The 90% CL upper limit on κ from an exposure M·T with N_90% = 2.44
(Poisson, zero observed background) is:

$$
\kappa^{90\%}(m_{A'}) = \sqrt{\dfrac{N_{90\%}}{R_\text{abs}(\kappa=1,\,m_{A'})\times M\!\cdot\! T}}
$$

### 2.3 Standard silicon case

The DAMIC-M PRL 2025 hidden-photon limit uses **Si's measured optical
dielectric function** (`Si_eps_electron_opticallimit.dat`) — real
ellipsometry/EELS data, not a model. Key features:

- Absorption onset at the direct gap ≈ 3.4 eV (indirect gap 1.11 eV does
  not contribute to k → 0 dipole absorption)
- A sharp Si plasmon peak at ω ≈ 16–17 eV, ELF ≈ 30, producing a
  characteristic sharp *minimum* in the exclusion curve (a plasmon is a
  resonant enhancement of the absorption rate, which strengthens — not
  weakens — the limit at that mass)
- Reported as **most stringent hidden-photon limit for m_A′ = 2.5–24 eV**
  in the PRL

This is the "textbook" case: no bracketing, no scissors, no synthetic
model — the material's real measured ε(ω) is used directly.

### 2.4 SrCd₂Sb₂ (HypMat) case

No optical measurement of SrCd₂Sb₂ exists. Three model variants bracket the
unknown ε(ω), from crudest to most structured:

#### 2.4.1 Drude model (primary bracketing model)

The free-carrier Drude dielectric function with a sub-gap cutoff (no
absorption below the gap, since interband transitions are forbidden there):

$$
\varepsilon_1(\omega) = \varepsilon_\infty - \frac{\omega_p^2}{\omega^2+\gamma^2}
\qquad
\varepsilon_2(\omega) = \begin{cases}
\dfrac{\omega_p^2\,\gamma}{\omega(\omega^2+\gamma^2)} & \omega \geq E_\text{gap} \\
0 & \omega < E_\text{gap}
\end{cases}
$$

Both parameters follow from known material properties with **no additional
DFT input**:

**Plasma frequency** ω_p from the valence electron density (Z_val = 16 per
formula unit — Sr 2s, 2×Cd 2s, 2×Sb 5 valence electrons):

$$
n_e = \frac{\rho\,N_A}{M_W}\,Z_\text{val}
= \frac{5.76\times6.022\times10^{23}}{555.96}\times 16 \approx 9.98\times10^{22}\ \text{cm}^{-3}
$$

$$
\omega_p = \sqrt{\frac{n_e\,e^2}{\varepsilon_0\,m_e}} \approx 11.7\ \text{eV}
$$

**Damping** γ from the DC conductivity σ_DC = 300 Ω⁻¹cm⁻¹ via the Drude
relation σ_DC = ε₀ ω_p² / γ:

$$
\gamma = \frac{\varepsilon_0\,\omega_p^2}{\sigma_\text{DC}} \approx 0.062\ \text{eV} \ll \omega_p
$$

This narrow damping (γ/ω_p ≈ 0.005) gives a *sharp, well-resolved* ELF
peak at the screened plasmon frequency ω_s = ω_p/√ε_∞:

| Bracket | ε_∞ | ω_s |
|---|---|---|
| Unscreened | 1 | ≈ 11.7 eV |
| Si-like screened | 12 | ≈ 3.4 eV |

Both brackets are run as independent scan configs (`hypmat_unscreened`,
`hypmat_screened`), giving an upper/lower rate bracket that differs by
≈ (ε_∞)² in the weak-absorption regime.

**What the Drude model still misses:** it is a free-carrier model and does
not capture interband transitions above the gap; it has no q-dependence
(only relevant for the optical/dark-photon channel, which is exactly k→0);
ε_∞ itself remains bracketed rather than measured.

#### 2.4.2 QCDark2-ELF variant (intermediate model)

A third variant uses Si's **measured** optical dielectric function
(the same file as §2.3) with `band_gap_eV = 2.1 eV` enforced at darkelf
init time — the scissors-shifted *direct* gap consistent with the DM-e
scissors calculation (§1.7). This trades the featureless Drude curve for
one with real (if Si-specific) interband structure and a sharp Si plasmon
near 17–19.5 eV, at the cost of the interband structure and plasmon
position belonging to Si rather than SrCd₂Sb₂.

| Model key | ELF source | Plasmon | Status |
|---|---|---|---|
| `hypmat_unscreened` | Drude, ε_∞=1, σ_DC=300 Ω⁻¹cm⁻¹ | None (overdamped in the original step-function version; sharp in the corrected Drude version, §2.4.1) | Primary optimistic bracket |
| `hypmat_screened` | Drude, ε_∞=12, σ_DC=300 Ω⁻¹cm⁻¹ | Sharp at ≈3.4 eV | Primary conservative bracket |
| `hypmat_qcdark2` | Si measured optical + 2.1 eV scissors onset | Si plasmon ≈17–19.5 eV | Intermediate — awaiting AiiDA wavefunctions |

**Correct density is now used in all three variants:** ρ_T = 5.76 g/cm³
(previously 8.0 g/cm³ placeholder — the 8.0 → 5.76 correction lowers the
rate by 5.76/8.0 = 0.72×, weakening the projected κ limit by
√(8.0/5.76) ≈ 1.18×, since R_abs ∝ 1/ρ_T linearly).

The Drude bracket (ε_∞ = 1 vs 12) remains the dominant systematic — a
factor of 12 on κ — and can only be narrowed by a refractive-index
measurement (ε_∞ = n²).

### 2.6 Numerical binning requirement for the monochromatic signal

The boxcar representation of §2.2 is exact at the CSV level (rate density
R_abs/δ inside a window of physical width δ = `binsize_eV`, zero outside).
Turning that CSV into the histogram consumed by `ChargeIonization::FoldToNe`
requires `RateTable::MakeTH1D` (§1.5) to resample it onto `nbins` output
bins — and because that resampling is **bin-center interpolation**, the
output bin width must be fine enough to actually land inside the δ-wide
boxcar at essentially every mass point on the scan grid. This is a
correctness requirement specific to the dark photon channel: DM-electron
and Migdal spectra are smooth over eV-scale features, so bin-center
sampling at nbins≈200 introduces no significant error there.

**Two artifacts from under-resolving the boxcar (nbins=200, bin-center):**

- **Missing low-mass coverage.** With `Emin_eV=0.34`, `Emax_eV≈100`, and
  nbins=200, the first bin center sits at ≈0.589 eV. Every scanned mass
  m_A′ between the physical threshold (0.34 eV) and ≈0.589 eV falls between
  bin centers and is never sampled inside its own boxcar → zero signal is
  assigned to masses that are, physically, fully accessible. Limit coverage
  therefore starts well above the true band-gap threshold.
- **~10× signal overcounting at "lucky" bins.** When a bin center *does*
  happen to fall inside the δ=0.05 eV-wide boxcar, `MakeTH1D` assigns that
  bin the boxcar's rate density R_abs/δ as its content. The histogram's
  implicit physical interpretation is that this density applies across the
  *entire* bin width (≈0.498 eV at nbins=200), so downstream integration
  effectively counts R_abs × (0.498/0.05) ≈ 10× the true rate. Because the
  90% CL limit scales as κ ∝ (rate)^(−1/2), a 10× rate inflation makes the
  sensitivity at those ~23 "lucky" mass points appear ≈3× better in κ than
  physically justified — while the surrounding masses (missed entirely) show
  zero.

The net effect at nbins=200 was a **staircase artifact**: sharp, jagged,
alternating over/under-estimates of sensitivity as a function of m_A′,
rather than the smooth curve the physics predicts. An earlier attempt fixed
this by switching `RateTable::MakeTH1D` from bin-center sampling to
trapezoid integration over each bin (correctly captures the boxcar
regardless of where it falls within a bin), but this produced a *different*
artifact: `ChargeIonization::FoldToNe` still evaluates the ionization table
P(n_e|E) at the bin center, which — at nbins=200 (bin half-width ≈0.25 eV)
— can be up to 0.25 eV away from the true mass m_A′. Near an ionization
threshold (E_gap + n×ε_h), that mismatch flips events between adjacent n_e
bins, producing exaggerated, misaligned staircase steps. Trapezoid
integration was therefore reverted, leaving `MakeTH1D` at bin-center
sampling (current state of `src/io/RateTable.cc`).

**Correct fix: increase `nbins` to resolve the boxcar directly.** Setting
`nbins=4000` over `Emin_eV=0.34`–`Emax_eV=100` gives a bin width of
≈0.025 eV — matched to the boxcar half-width (δ/2 = 0.025 eV for
δ=0.05 eV) — so that:

- every scanned m_A′ has a bin center within ≈0.012 eV of the true mass →
  P(n_e|E) is evaluated at essentially the correct energy, and the
  staircase steps that remain are the **physical** ones, landing precisely
  at the ionization thresholds E_gap + n×ε_h = 0.34 + n×1.7778 eV
  (§4.1, Klein), not numerical artifacts;
- coverage extends down to the true band-gap threshold rather than the
  first surviving bin center of a coarse grid;
- the ~10× overcounting bias vanishes because the bin width now
  approximately matches the physical boxcar width, so bin-center sampling
  and the true integrated rate agree to good approximation — reproducing
  the same physics the (reverted) trapezoid-integration fix targeted,
  without its ionization-table misalignment side effect.

Six dark-current-sweep HypMat-unscreened scan configs were updated to
`nbins=4000` accordingly:
`darkphoton_scan_hypmat_unscreened_ne_dc1e{1,2,3,4,6}.json` plus
`darkphoton_scan_hypmat_unscreened_ne_1kgyr_fresh.json` (the sixth DC point
in the sweep), along with a new Si-Mermin reference projection,
`configs/darkphoton_scan_si_mermin_ne_proj_nbins4000.json` — a sibling of
the existing `darkphoton_scan_si_mermin_ne_proj.json` (still `nbins=200`).
Other dark photon configs not yet migrated to `nbins=4000` (most
`srcd2sb2_*`, `eu5in2sb6_*`, the `*_dc1e2/3_{1,10,100}gyr` variants, and
the `old_grid_*` files) still carry `nbins=200` and inherit the same
coverage/overcounting caveats described above until updated — any figure
drawing from those configs should be treated as provisional.

**Physics finding enabled by the fix — DC sweep is signal-limited above
~2 eV.** With correct normalization, the six DC levels (10⁻¹–10⁻⁶ e⁻/pix/day)
produce **nearly identical** κ(m_A′) curves for m_A′ ≳ 2 eV. Dark-current
background enters the n_e≥2 bins only through pileup-like combinatorics
(∝ DC²), which is negligible at every DC level considered once n_e≥2 is
populated — i.e. sensitivity there is set entirely by signal statistics and
exposure, not background. DC only matters on the **low-mass plateau**
(0.34–1 eV), where only the n_e=1 bin is populated and single-pixel dark
current (∝ DC) is the dominant background. This is a genuine physics result
of the corrected calculation, not an artifact — and it means DC-level
comparisons for this channel should be read as differing only near
threshold, not across the full mass range.

---

## Part III — The Migdal Effect

### 3.1 Physical process: why it is needed

For elastic DM–nucleus scattering, the maximum nuclear recoil energy is
E_n^max = 2μ²_χN v²/m_N. For a Si nucleus (m_N ≃ 26 GeV) and v ~ 10⁻³c:

| m_χ | E_n^max | Detectable in a CCD? |
|---|---|---|
| 1 GeV | ~100 eV | Yes |
| 100 MeV | ~1 eV | Marginal |
| 10 MeV | ~0.01 eV | No |

Below ~100 MeV, the nuclear recoil alone is invisible. The **Migdal
effect** converts part of the recoil momentum into an *electronic*
excitation via a sudden-perturbation (shake-off) process, which *can* be
detected as an ionization signal — extending sensitivity down to ~1 MeV.

### 3.2 The sudden approximation and factorization

Because the nuclear scattering timescale (~1/q ~ 10⁻²² s) is far shorter
than the electronic relaxation timescale (~1/ω ~ 10⁻¹⁷–10⁻¹⁵ s), the
process factorizes into two sequential steps:

```
Step 1:  χ + nucleus_A  →  χ' + nucleus_A*     (nuclear scattering)
Step 2:  nucleus_A*  →  nucleus_A + e⁻ + ...    (electronic shake-off)
```

The differential rate for an electronic excitation of energy ω is:

$$
\frac{dR}{d\omega} = \frac{\rho_\chi}{m_N m_\chi}\, I(\omega)
\int_{v_\text{min}}^{v_\text{esc}+v_E} v\,f(v)\,J(v,\omega)\, dv
$$

**I(ω) — shake-off probability (material physics):**

$$
I(\omega) = \frac{1}{E_n}\frac{dP}{d\omega}
= \frac{2\alpha_\text{EM}\,Z_\text{ion}^2}{3\pi^2\,\omega^4}
\int_0^{k_\text{max}} k^2\,\text{Im}\!\left[-\frac{1}{\varepsilon(\omega,k)}\right] dk
$$

driven by the **same ELF**, Im[−1/ε(ω,k)], that governs DM–electron
scattering (§1.3) and dark photon absorption (§2.2) — the Migdal effect is
literally the electron energy-loss function of the material, probed by a
recoiling nucleus instead of an external probe.

**J(v,ω) — nuclear kinematics (DM/mediator physics), free-nucleus
approximation:**

$$
J(v,\omega) \propto \frac{A^2\,\sigma_n}{v\,\mu_{\chi N}^2}\,F_\text{DM}^2(q)
$$

with A² = 784 for Si (A=28) the coherent nuclear enhancement, and σ_n the
DM-nucleon reference cross section at q₀ = αm_e.

**Kinematic threshold:**

$$
v_\text{min}(\omega) = \sqrt{2\omega/\mu_{\chi N}}
$$

For ω ≃ E_gap = 1.11 eV and m_χ = 10 MeV, v_min ≈ 6×10⁻³c — already
strongly suppressed by the halo tail, setting the effective low-mass floor
of the Migdal search.

### 3.3 Mediator form factor (DM–nucleus vertex)

- **Heavy mediator**: F_DM²(q) = (q₀²+m_med²)²/(q²+m_med²)² → 1 (thermal
  freeze-out benchmark, q-independent)
- **Light (massless) mediator**: F_DM²(q) = (q₀/q)⁴ (freeze-in benchmark,
  strongly enhanced at low q)

This is a **distinct namespace** from the DM–electron mediator of Part I —
same "heavy/light" language, different vertex (DM-nucleus vs DM-electron).

### 3.4 Standard silicon case

| Parameter | Value | Source |
|---|---|---|
| E_gap | 1.11 eV | darkelf Si optical data, 130 K |
| ombar (phonon frequency) | 0.03 eV | Si |
| E_nth (nuclear recoil cutoff) | 0.12 eV = 4×ombar | free-nucleus approximation validity floor |
| Z_ion(0) | 4.0 | effective ionic charge (`Si_Zion.dat`, Brown et al. 2006) |
| A | 28 | Si mass number → A²=784 coherent enhancement |
| ELF source | `Si_mermin.dat` (Mermin extension of optical data) | darkelf |
| Approximation | free-nucleus, `method="grid"` | impulse approximation and Ibe atomic Migdal deferred |
| SHM | v₀=238, v_E=253.7, v_esc=544 km/s | darkelf default |

**Reference rate (heavy mediator, σ_n=10⁻³⁶ cm², integrated ω∈[E_gap,20 eV]):**

| m_χ | Rate [events/kg/yr] |
|---|---|
| 100 MeV | ≈3.5×10⁴ |
| 1 GeV | ≈1.1×10⁶ |

Rate scales **exactly linearly** in σ_n (it enters only as an overall
prefactor of J), which the grid generator exploits as a fast path
(§4.2 of `Migdal_Integration_Plan.md`). The DAMIC-M PRL 2025 reports the
Migdal channel (using darkelf) as the **most stringent limit for
m_χ = 1–35 MeV/c²**, with the explicit caveat that "the Migdal effect is not
yet experimentally calibrated."

The detector pipeline from `ChargeIonization` onward is **100% identical**
to DM–electron scattering: darkelf's `dRdomega_migdal(ω)` already outputs
the spectrum in electronic excitation energy ω, the same physical quantity
as E_e in Part I.

### 3.5 SrCd₂Sb₂ case: explicitly out of scope

Unlike the DM-electron and dark-photon channels, the Migdal effect for
SrCd₂Sb₂ is **not currently computed** in this framework — and should not
be reported as a projection. The reason is structural, not a missing
config:

- I(ω) (§3.2) requires the material's electron ELF Im[−1/ε(ω,k)] over a
  **finite-k grid**, exactly the same input DarkELF needs for absorption.
  DarkELF ships this for Si (`Si_mermin.dat`) but has **no equivalent file
  for SrCd₂Sb₂**.
- Unlike the dark-photon channel (§2.4), there is no simple Drude
  substitute here: I(ω) is a **k-integrated** quantity (0 → k_max), so a
  frequency-only model (Drude ε(ω) at k=0) cannot stand in — the Migdal
  rate genuinely needs the momentum-dependent ELF, which is precisely the
  finite-q information the Drude bracketing strategy was built to avoid
  needing.
- The scissors trick used for DM–electron scattering (§1.7) shifts *Si's*
  band structure; it does not, by itself, produce Zion(k) or a
  finite-k ELF file in the format DarkELF's Migdal module expects.

**What would be required to add this channel:** a DarkELF-format ELF grid
Im[−1/ε(ω,k)] for SrCd₂Sb₂ over the relevant (ω, k) range, plus an
effective ionic charge Z_ion(k) for the constituent atoms (Sr, Cd, Sb).
Both require either full ab-initio calculation or, at minimum, a
material-specific momentum-dependent model — beyond what the Drude/σ_DC
bracketing approach can supply. This is recorded as explicitly out of
scope in `Migdal_Integration_Plan.md` §4.8 and remains so here.

---

## Part IV — Charge Ionization Model (shared by all three channels)

### 4.1 Klein's formula: mean electron–hole pair creation energy

Once a deposited energy ω ≥ E_gap is fixed, the **mean** number of
electron–hole pairs created is set by the mean pair-creation energy ε_h
(historically called the "Fano" or "ionization" energy):

$$
\langle n_e\rangle = \frac{\omega - E_\text{gap}}{\varepsilon_h}
$$

For real silicon, ε_h = 3.8 eV is **directly measured** (100 K, PRD 102,
063026) — not derived from a formula. For a hypothetical material with no
measurement, Klein (1968) proposed an empirical linear relation between
ε_h and the band gap, calibrated across several elemental and compound
semiconductors (Si, Ge, GaAs, CdS, ...):

$$
\varepsilon_h \approx a\,E_\text{gap} + b
$$

CCDarkSens uses this class of relation with the coefficients

$$
\boxed{\varepsilon_h = 2.67\,E_\text{gap} + 0.87\ \text{eV}}
$$

which for SrCd₂Sb₂ (E_gap = 0.34 eV) gives:

$$
\varepsilon_h = 2.67\times0.34 + 0.87 = 0.9078 + 0.87 = \mathbf{1.7778\ eV}
$$

This value (`configs/darkphoton_one_point_eu5in2sb6.json` and the Sr2Cb2Sd
scan configs) **supersedes** the earlier ε_h = 1.6 eV placeholder used in
the original HypMat unified-model document. Note the repository also
contains a second, coarser Klein-type fit — the widely cited Alig–Bloom
form ε_h = 2.8 E_gap + 0.5 eV — used only in exploratory heatmap plots
(`utils/plot_band_gap_2d_heatmap.py`); it is **not** the parameterization
used for the production SrCd₂Sb₂ configs and is noted here only to avoid
confusion between the two curves if both appear in a figure.

The actual n_e for a given ω is **not sharp** — it follows a
Fano-broadened distribution P(n_e|ω), captured in a lookup table
(§4.3), not by rounding ⟨n_e⟩.

### 4.2 Ionization scenarios: the "D-equal" convention and its siblings

Because there is **no known fundamental relation** ε_h(E_gap) valid for
arbitrary low-gap materials (Klein's fit is itself only a phenomenological
average over a handful of real semiconductors), CCDarkSens's band-gap
phenomenology framework treats (E_gap, ε_h) as two **independently
specifiable** knobs and labels each combination with a scenario ID:

| Scenario ID | E_gap | ε_h | Interpretation |
|---|---|---|---|
| **ref** | 1.2 eV | 3.8 eV | Reference Si, unscaled `p100K_table.csv` |
| **B-thresh** | E_gap (new) | 3.8 eV (fixed, Si-like) | Lower the threshold only, keep Si's pair scale |
| **D-equal** | E_gap (new) | **= E_gap (new)** | Aggressive equal-scale pheno; not real bulk Si |
| **A-ratio** | E_gap (new) | 3.8 × E_gap/1.2 | Preserve Si's ε_h/E_gap ratio |
| **Klein** | E_gap (new) | 2.67×E_gap + 0.87 | This document's §4.1 relation; used for SrCd₂Sb₂ |
| **C-grid** | user | user | Free 2D grid, no rule |

The SrCd₂Sb₂ production configs use the **Klein** scenario specifically
(`ion_scenario: "Klein"` in `configs/Sr2Cb2Sd/write_configs.py`), i.e. ε_h
is *not* pinned equal to E_gap (that would be D-equal) nor fixed at Si's
3.8 eV (B-thresh) — it follows the Klein formula above.

### 4.3 The P(n_e|E) table and the anchored energy map

The full ionization model is a table P(n_e|E): CSV columns
`Er_eV, P1, P2, ..., Pk`, loaded by `ChargeIonization::LoadCSV_` and
consumed by `ChargeIonization::FoldToNe`, which convolves a dR/dE histogram
against this table (weighted by exposure) to produce the observable n_e
spectrum — see `include/ccdarksens/response/ChargeIonization.hh`.

For materials without a dedicated Monte Carlo P(n_e|E) calculation, the
table is **rescaled from the Si reference table** `data/p100K_table.csv`
(E_gap=1.2 eV, ε_h=3.8 eV, 100 K) via the anchored energy map:

$$
E'(E) = E_\text{gap}^\text{ref} + \bigl(E - E_\text{gap}^\text{new}\bigr)\times
\frac{\varepsilon_h^\text{ref}}{\varepsilon_h^\text{new}} \qquad \text{for } E \geq E_\text{gap}^\text{new}
$$

$$
P_\text{new}(n_e\mid E) = \begin{cases}
0 & E < E_\text{gap}^\text{new} \\
P_\text{ref}\bigl(n_e \mid E'(E)\bigr) & \text{otherwise}
\end{cases}
$$

(linear interpolation on the reference table). **Physical interpretation:**
the map asks "at what energy does bulk Si have the same *fractional*
excitation above threshold as the new material at energy E?" — it
anchors the ionization *shape* (Fano broadening, multiplicity structure) to
Si's measured behavior while relocating the *threshold and scale* to the
new material's (E_gap, ε_h).

For SrCd₂Sb₂ (E_gap=0.34 eV, ε_h=1.7778 eV via Klein):

$$
E'(E) = 1.2 + (E-0.34)\times\frac{3.8}{1.7778} = 1.2 + (E-0.34)\times 2.138
$$

built by `utils/build_p100K_scaled.py --band-gap-eV 0.34 --eh-pair-eV 1.7778 --scenario Klein`.

**What this does and does not claim:** the turn-on energy (0.34 eV) and the
mean-n_e scale (via ε_h) are physically motivated (Klein's relation,
however approximate); the detailed *shape* of P(n_e|E) — Fano factor,
sub-gap phonon tail — is inherited from Si and is **not** a first-principles
SrCd₂Sb₂ calculation.

---

## Part V — The Unified PLR Sensitivity Pipeline

All three channels (Parts I–III) and the ionization model (Part IV) feed
the **same** downstream chain, implemented once in C++ and dispatched by
`model.type` in the JSON config (`dm_electron`, `dark_photon`, `migdal`):

```
Rate CSV [events/kg/yr/eV]           ← Part I, II, or III
  (continuous dR/dE, or boxcar for dark photon)
        │
        ▼
RateTable::MakeTH1D                   bin-center interpolation onto uniform E histogram
                                       (needs fine nbins for narrow boxcar features — §2.6)
        │
        ▼
ChargeIonization::FoldToNe            × exposure (M·T), fold against P(n_e|E)  ← Part IV
        │
        ▼
n_e spectrum  S(n_e)
        │
        ▼
EfficiencyMC / PatternRates           diffusion, readout noise, pattern classification
        │  (observable_bins: "ne" or "pattern")
        ▼
Profile likelihood ratio scan          background model (dc_flat_migration or Bp_theta_Br)
        │
        ▼
Upper limit at 90% CL:
   σ_e(m_χ)   — DM-electron scattering
   κ(m_A')    — dark photon absorption
   σ_n(m_χ)   — Migdal (per-nucleon)
```

### 5.1 What differs between channels, and what does not

| Stage | DM-electron | Dark photon | Migdal |
|---|---|---|---|
| Rate engine | QEDark / QCDark2 | DarkELF `absorption` | DarkELF `Migdal` |
| Input spectrum shape | Continuous dR/dE(E) | Monochromatic boxcar at m_A' | Continuous dR/dω, peaked near E_gap |
| Scan axis (x) | m_χ | m_A' | m_χ |
| Scan axis (y, linear in) | σ_e | κ² | σ_n |
| `ChargeIonization::FoldToNe` | identical call | identical call | identical call |
| `EfficiencyMC` / `PatternRates` | identical | identical | identical |
| PLR machinery | identical | identical | identical |
| Background model | identical (`dc_flat_migration` or `Bp_theta_Br`) | identical | identical |

Because the y-axis enters each rate **linearly** (σ_e in Part I, κ² in
Part II, σ_n in Part III), all three channels reuse the same
fast-rescaling trick when regenerating grids: compute once at a reference
coupling, then scale analytically for every other grid point — this is
exact by construction, not an approximation.

### 5.2 Config field summary (`model.type`)

```json
// DM-electron
"model": { "type": "dm_electron", "mediator": "heavy|light",
           "rates_dir": "data/qcdark2_rates/...", ... }

// Dark photon
"model": { "type": "dark_photon", "rates_dir": "data/darkphoton_rates/...", ... }

// Migdal
"model": { "type": "migdal", "mediator": "heavy|light",
           "target_nucleus": "Si28", "nuclear_A": 28, "nuclear_Z": 14,
           "rates_dir": "data/migdal_rates/...", ... }
```

```json
// Shared ionization block (Part IV) — used identically by all three
"response": {
  "charge_ionization": {
    "table_csv":   "data/p100K_gap0p34_eh1p7778.csv",
    "band_gap_eV": 0.34,
    "eh_pair_eV":  1.7778,
    "scenario":    "Klein"
  }
}
```

---

## Part VI — Parameter Summary: Si vs SrCd₂Sb₂ Across All Channels

| Quantity | Standard Si | SrCd₂Sb₂ (HypMat) | Status for SrCd₂Sb₂ |
|---|---|---|---|
| Band gap (DM-e / ionization) | 1.2 eV (QEDark) / 1.11 eV (darkelf) | 0.34 eV (indirect, working value) | Scissors proxy (§1.7) |
| Band gap (dark photon, direct) | ≈3.4 eV (measured) | 2.1 eV (scissors direct) | Model-dependent per variant (§2.4) |
| ε_h (mean e-h pair energy) | 3.8 eV (measured) | 1.7778 eV (Klein: 2.67×0.34+0.87) | Phenomenological (§4.1) |
| Density ρ_T | 2.329 g/cm³ | 5.76 g/cm³ (corrected from 8.0 placeholder) | Crystallographic (§1.7) |
| DM-e crystal form factor | `Si_f2.txt` / `Si_comp.h5` (native Si DFT) | Scissor-shifted Si (`Si_fast_gap0p34.h5`) | Si-proxy — shape not native |
| Dark photon ε(ω) | Measured optical data | Drude (ε_∞={1,12}) or QCDark2-ELF (Si proxy) | Bracketed / proxy |
| Migdal ELF Im[−1/ε(ω,k)] | `Si_mermin.dat` | **Not available** | **Out of scope** (§3.5) |
| Valence electrons/cell (DM-e) | 4 (native) | 4 (inherited from Si proxy, B6 unresolved) | Open assumption |
| Valence electrons (Drude ω_p) | n/a | Z_val=16/formula unit (Sr+2Cd+2Sb) | Used only in dark photon channel |
| SHM (v₀, v_E, v_esc) | (238, 263, 544) km/s | same | Not material-dependent |

**Reading this table honestly:** for SrCd₂Sb₂, the dark-photon channel is
the most defensible (its rate formula is exact given ε(ω); only ε(ω)
itself is bracketed/modeled). The DM-electron channel is a genuine
sensitivity trend, not a quantitative prediction, because the crystal form
factor's shape is inherited from Si. The Migdal channel does not currently
exist for this material and would require new ELF input data, not a config
change.

---

## References

- Klein, C. A., *Bandgap Dependence and Related Features of Radiation
  Ionization Energies in Semiconductors*, J. Appl. Phys. **39**, 2029
  (1968) — origin of the linear ε_h(E_gap) empirical relation used in §4.1.
- Essig, R., Mardon, J., Volansky, T., *Direct Detection of Sub-GeV Dark
  Matter*, Phys. Rev. D 85, 076007 (2012).
- Essig, R. et al., *Direct Detection of sub-GeV Dark Matter with
  Semiconductor Targets*, JHEP 05 (2016) 046 (the "QEDark" paper).
- Dreyer, C., Essig, R., Fernandez-Serra, M., Singal, A., Zhen, C.,
  *A New Direction for Dark Matter Direct Detection: Silicon Devices with
  Sub-eV Threshold*, Phys. Rev. D 109, 115008 (2024) — QCDark, used for the
  DM–electron channel in the DAMIC-M PRL.
- Hochberg, Y. et al., *Absorption of light dark matter in semiconductors*,
  Phys. Rev. D 95, 023013 (2017).
- Knapen, S., Kozaczuk, J., Lin, T., *DarkELF: A Python package for dark
  matter scattering in dielectric targets*, Phys. Rev. D 104, 015031
  (2021), arXiv:2104.12785 — the DarkELF package used for both the dark
  photon and Migdal channels.
- Ibe, M., Nakayama, W., Shigeki, Y., Yanagida, T., arXiv:1707.07258 —
  atomic Migdal calculation (not used by this framework; ELF-based method
  preferred for semiconductors).
- Essig, R., Pradler, J., Sholapurkar, M., Yu, T.-T., *Relation between the
  Migdal Effect and Dark Matter–Electron Scattering in Isolated Atoms and
  Semiconductors*, Phys. Rev. Lett. 124, 021801 (2020), arXiv:1908.10881.
- DAMIC-M Collaboration, *Probing Benchmark Models of Hidden-Sector Dark
  Matter with DAMIC-M*, Phys. Rev. Lett. **135**, 071002 (2025),
  DOI: 10.1103/2tcc-bqck. Reports DM-electron (heavy and ultralight
  mediator), hidden-photon absorption, and Migdal-effect limits using
  exactly the QEDark/QCDark, DarkELF-absorption, and DarkELF-Migdal engines
  documented here. Limits repository:
  https://github.com/DAMIC-M/LBC_2025_HSBenchmark_Limits
- Si charge-yield reference (100 K, ε_h=3.8 eV): Phys. Rev. D 102, 063026.
