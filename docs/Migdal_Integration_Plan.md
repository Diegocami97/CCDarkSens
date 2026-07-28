# Migdal Effect — Physics, Detector Response, and CCDarkSens Integration Plan

---

## 1. Physics of the Migdal effect

### 1.1 Motivation: why standard nuclear recoils fail at low mass

A dark matter particle of mass m_χ scattering elastically off a nucleus of
mass m_N at velocity v transfers a nuclear recoil energy:

    E_n = q² / (2 m_N)   where   q ∈ [q_min, q_max] = μ_χN v [1 ∓ cos θ_cm]

The maximum nuclear recoil energy is:

    E_n^max = 2 μ_χN² v² / m_N

For a Si nucleus (m_N ≃ 26 GeV) and v ~ 10⁻³ c (SHM peak):

| m_χ    | E_n^max   | Detectable? |
|--------|-----------|-------------|
| 1 GeV  | ~100 eV   | Yes (ionization) |
| 100 MeV| ~1 eV     | Marginal (1 phonon) |
| 10 MeV | ~0.01 eV  | No (sub-phonon) |
| 1 MeV  | ~0.1 meV  | Completely invisible |

For DM masses below ~100 MeV, nuclear recoils alone produce no detectable
signal in a Si CCD detector.  The **Migdal effect** extends sensitivity by
converting part of the nuclear recoil momentum into an electronic excitation
that *can* be detected as an ionization signal.

### 1.2 The Migdal mechanism

When a nucleus recoils suddenly (faster than the electron relaxation time), the
electron cloud cannot follow instantaneously.  In the center-of-mass frame of
the recoiling nucleus the electron cloud acquires a collective velocity −v_N,
where v_N = q_N / m_N is the nuclear recoil velocity.  This sudden boost drives
electronic transitions from the ground state to excited states — the
**sudden-perturbation (shake-off) mechanism**.

The process can be visualized in two steps:

```
Step 1:  χ + nucleus_A  →  χ' + nucleus_A*     (nuclear scattering, time ~ 1/q)
Step 2:  nucleus_A*  →  nucleus_A + e⁻ + ...    (electronic de-excitation, time ~ 1/ω)
```

Because the nuclear scattering time (~1/q ~ 10⁻²² s) is much shorter than the
electron relaxation time (~1/ω ~ 10⁻¹⁷–10⁻¹⁵ s), the sudden approximation
applies: the nuclear scattering and the electronic emission factorize.

The differential rate for producing an electronic excitation of energy ω is:

    dR/dω = (ρ_χ / m_N m_χ) × I(ω) × ∫_{v_min}^{v_esc+v_E} v f(v) J(v,ω) dv

where each factor has a clear physical interpretation:

**I(ω) — shake-off probability (material physics):**

    I(ω) = (1/E_n) dP/dω

This is the ionization probability per unit nuclear recoil energy,
integrated over the momentum k transferred to the electronic system:

    I(ω) = (2 α_EM Z_ion²) / (3π² ω⁴)  ×  ∫_0^{k_max} k² × Im[-1/ε(ω,k)] dk

The integrand contains:
- `α_EM` — fine structure constant
- `Z_ion(k)` — effective ionic charge (momentum-dependent, loaded from
  `Si_Zion.dat` in darkelf; for Si, Z_ion(0) = 4.0)
- `Im[-1/ε(ω,k)]` — the electron energy-loss function (ELF) of the material,
  evaluated at excitation energy ω and momentum k.  This encodes the material's
  electronic response.  For Si, darkelf uses the Mermin extension
  (`Si_mermin.dat`) derived from optical data.

Note that I(ω) scales as E_n, so dR/dω ∝ I(ω) ∝ E_n × (nuclear recoil
factor) — meaning larger nuclear recoils (higher m_χ) produce more
electronic emission.

**J(v,ω) — nuclear kinematics (DM and mediator physics):**

In the free-nucleus approximation (used by default in darkelf):

    J(v,ω) = (2π² A² σ_n) / (v μ_χN²) × F_DM²(q) × ...  [integrated over E_n]

where:
- `A²` — coherent nuclear enhancement (A=28 for Si → factor of 784 over
  single-nucleon cross section)
- `σ_n` — DM-nucleon reference cross section at q₀ = α m_e ≃ 3.7 keV
- `F_DM²(q)` — mediator form factor (see §1.4)
- the integral runs over nuclear recoil energies E_n consistent with DM
  velocity v and electronic excitation ω, subject to E_n > E_nth (the
  threshold below which the free approximation breaks down)

**f(v) — velocity distribution:**

darkelf uses the Standard Halo Model (SHM):
- v₀ = 238 km/s (most probable speed)
- v_E = 253.7 km/s (Earth velocity)
- v_esc = 544 km/s (escape velocity)

**v_min — kinematic lower bound on DM velocity:**

    v_min(ω) = sqrt(2ω / μ_χN)

This sets the low-mass kinematic cutoff: for ω = E_gap ≃ 1.11 eV and
m_χ = 10 MeV, v_min ≃ 6×10⁻³ c — comparable to v_esc and already strongly
suppressed by the halo velocity integral.

### 1.3 Reference rate values (darkelf, Si, heavy mediator)

At σ_n = 10⁻³⁶ cm² and integrating dR/dω over ω ∈ [E_gap, 20 eV]:

| m_χ    | Total rate (events/kg/yr) |
|--------|--------------------------|
| 100 MeV| ~3.5 × 10⁴             |
| 1 GeV  | ~1.1 × 10⁶             |

The rate scales linearly with σ_n (exact, by construction — σ_n only enters
as an overall prefactor in J).  This allows the sigma fast path in the
generator (§4.2).

### 1.4 Mediator form factor (DM-nucleus vertex)

The DM-nucleus interaction is mediated by a boson with mass m_med.  The
momentum-transfer form factor is:

- **Heavy mediator** (m_med >> q_typ):
  F_DM²(q) = (q₀² + m_med²)² / (q² + m_med²)²  →  1   (momentum-independent)
  Benchmark: thermal freeze-out.  Cross section is constant in q.

- **Light (massless) mediator** (m_med → 0):
  F_DM²(q) = (q₀ / q)⁴   (strong enhancement at low q)
  Benchmark: freeze-in, ultralight mediator.  Cross section grows steeply
  at low momentum transfer — relevant for the lightest accessible DM masses.

In darkelf, the mediator is selected via:
- `mMed = 1e30` (or any value >> v m_χ) → heavy
- `mMed = -1` → massless/light

This is a distinct namespace from the DM-electron mediator in DMElectronModel —
the same "heavy/light" language describes a different physical vertex (DM-nucleus
vs DM-electron).  `DMNucleonConfig::mediator` carries this field.

### 1.5 Approximations and validity

darkelf's Migdal implementation uses the **soft/free-nucleus approximation**:

- **Soft limit**: the photon emitted by the Migdal process has wavelength much
  larger than the atomic radius, so the nucleus appears as a point charge.
  This is valid for ω << q²/(2 m_e) ≃ 10 eV for typical Migdal momenta.
  Since we integrate to ω_max ~ 20 eV, this is marginally satisfied.

- **Free nucleus**: the nucleus is treated as unbound (ignoring lattice binding).
  The impulse approximation alternative accounts for phonon binding and is
  more accurate at low E_n.  darkelf supports `approximation="impulse"` but it
  requires the phonon spectrum `ombar` in the YAML file.  For Si, ombar = 0.03 eV
  and E_nth = 4×ombar = 0.12 eV.  We use the free approximation (default) with
  E_nth = 0.12 eV as the nuclear recoil lower cutoff.

- **ELF method**: `method="grid"` uses the precomputed dielectric function
  grid from Mermin/GPAW for the momentum integral in I(ω).  This is the most
  accurate option for Si.

The `Si_Migdal_FAC.dat` file (Ibe et al. atomic Migdal) is also present in
the darkelf data directory.  We do not use it (`method="Ibe"`) — the ELF-based
calculation is more reliable for a semiconductor target.

### 1.6 Comparison with DM-electron scattering

| Feature | DM-electron | Migdal |
|---|---|---|
| DM couples to | electron | nucleus |
| y-axis | σ_e [cm²] | σ_n [cm²] (per nucleon) |
| Coherent factor | 1 | A² = 784 |
| Form factor F_DM | same q-dependence | same q-dependence |
| Observable | ω → n_e pairs | ω → n_e pairs (**identical**) |
| Rate generator | QEDark/QCDark/QCDark2 | darkelf.dRdomega_migdal |
| Relevant mass range | ~1 MeV – 1 GeV | ~10 MeV – 10 GeV |
| Band-gap sensitivity | yes (lower gap → more signal) | yes (same) |
| Material model needed | crystal form factor | electron ELF |

**The detector pipeline from ChargeIonization onward is 100% identical.**
darkelf.dRdomega_migdal(ω) already outputs the spectrum in the electronic
excitation energy ω — the same quantity as E_e in DM-electron scattering.
The ChargeIonization table p(Q|ω) is applied identically to convert
deposited energy → number of electron-hole pairs Q → n_e spectrum.

### 1.7 Why DAMIC-M is competitive for Migdal

1. **Low threshold**: DAMIC-M's 2e⁻ threshold (1–2 e-h pairs) means sensitivity
   to ω near the band gap (1.11 eV for Si).  Most WIMP-search experiments have
   much higher thresholds.

2. **Large exposure**: 1 kg-yr projected exposure and low dark-current background
   give competitive rates even for the suppressed Migdal cross sections.

3. **Pattern analysis**: the multi-pixel pattern read-out suppresses backgrounds
   by selecting only single-pixel (1e⁻, 2e⁻) events consistent with
   low-energy ionization — the expected signature of Migdal events.

4. **A² enhancement**: the coherent A² = 784 factor over single-nucleon
   cross section gives Migdal limits on σ_n that are competitive with direct
   DM-electron limits on σ_e despite the additional nuclear recoil suppression.

### 1.8 DIM framework status

The collaboration DIM framework (`collab_frameworks/dim/`) has **no Migdal
implementation**.

- `calculate_rates/DarkElf_DMe.py` — the only rate calculator present; calls
  `Si.dRdomega_electron()` for DM-electron scattering only
- `DIM/DetectorResponse/include/dRdEHandler.h` — the string "Migdal" appears
  once in a comment listing example mediator label strings
  (`/// Model (i.e. massive, massless, Migdal, etc.)`); it is never
  instantiated or dispatched anywhere in the codebase
- `DIM/DetectorResponse/src/dRdEHandler.cpp` — generic CSV-table consumer;
  would accept Migdal rate tables but none are generated by DIM

The Migdal integration is therefore entirely new work in CCDarkSens.

---

## 2. Charge ionization step — is it the same?

**Yes, completely.** This is the key architectural insight.

`darkelf.dRdomega_migdal(omega)` integrates out all nuclear physics and
outputs dR/dω where **ω is the electronic excitation energy** [eV] —
the energy deposited into the electron cloud after the nucleus recoils.
This is physically and dimensionally identical to E_e in DM-electron
scattering: energy available to create electron-hole pairs in the detector.

The pipeline is:

```
DM-electron:
  dR/dE_e  [CSV, from QEDark/QCDark2]
     → ChargeIonization::FoldToNe(E_e)  →  p(Q|E_e)  →  n_e spectrum
     → EfficiencyMC → PatternRates → PLR

Migdal:
  dR/dω    [CSV, from darkelf via migdal/entry.py]
     → ChargeIonization::FoldToNe(ω)    →  p(Q|ω)    →  n_e spectrum  ← identical
     → EfficiencyMC → PatternRates → PLR                               ← identical
```

Two practical differences in the CSV compared to DM-electron tables:

1. **Energy range**: Migdal ω is kinematically bounded by ~m_χ v²/2.
   For m_χ = 100 MeV and v = v₀, ω_max ~ a few eV.  For m_χ = 1 GeV,
   ω_max can reach ~20 eV.  We set `Emax_eV = 20.0` in the rate tables
   (covers the full interesting range) and `Emin_eV = 0.0` (darkelf handles
   the band-gap threshold internally via `Enth`).

2. **Units**: same as DM-electron tables — events / kg / year / eV.
   The generator rescales from darkelf's 1/kg/yr/eV output directly.

---

## 3. darkelf Si parameters (confirmed from loaded instance)

| Parameter | Value | Source |
|---|---|---|
| `E_gap` | 1.11 eV | Si at 130 K (optical data) |
| `ombar` | 0.03 eV | Si phonon frequency |
| `Enth` (default) | 0.12 eV = 4×ombar | nuclear recoil lower cutoff |
| `ommax` | 99.3 eV | max ω in ELF grid |
| `Zion` | 4.0 | effective ionic charge for Migdal |
| `A` | 28.0 | Si mass number |
| ELF file | `Si_mermin.dat` | Mermin extension |
| Zion(k) file | `Si_Zion.dat` | Brown et al. 2006 |
| `electron_ELF_loaded` | True | — |

---

## 4. CCDarkSens integration plan

### 4.1 Status of existing components

All C++ infrastructure and the Python entry point are complete from prior work:

| Component | File | Status |
|---|---|---|
| Python rate entry point | `python/ccdarkphys/migdal/entry.py` | **Done** |
| DM-nucleon config struct | `include/ccdarksens/model/DMNucleonConfig.hh` | **Done** |
| C++ model header | `include/ccdarksens/model/MigdalModel.hh` | **Done** |
| C++ model implementation | `src/model/MigdalModel.cc` | **Done** |
| Scan app dispatch | `apps/ccdarksens_scan_dmelectron_pattern.cc` ~line 171 | **Done** |
| CMakeLists.txt entry | `CMakeLists.txt` | **Done** |
| Grid generator script | `utils/migdal_generate_grid.py` | **Done** |
| darkelf installation | `/Users/.../Software/DarkELF`, Python 3.11.9 | **Done** |

### 4.2 The sigma_n fast path in the generator

R ∝ σ_n exactly (σ_n enters only as an overall prefactor in J(v,ω)).
Therefore the generator (`migdal_generate_grid.py`) initializes darkelf once
per mass point (at sigma_ref = sigma_list[0]) and rescales all other
sigma_n values analytically:

    dR/dω(sigma_n) = (sigma_n / sigma_ref) × dR/dω(sigma_ref)

For a 50×30 (mass × coupling) grid this reduces darkelf inits from 1500 → 50.

### 4.3 Filename convention

The filename template follows `DMNucleonConfig::filename_template`:

```
dRdE_{target_nucleus}_{mediator}_m{mchi_MeV}_s{sigma_n_cm2}.csv
```

Example: `dRdE_Si28_heavy_m100.000000_s1.0e-36.csv`

Token substitution is performed by `MigdalModel::ResolvePath_()` in C++.
The generator uses the same format strings to produce matching filenames.

Formatting:
- `mchi_MeV`: `".6f"` — six decimal places (e.g. `100.000000`)
- `sigma_n_cm2`: `".1e"` — one significant digit in exponent (e.g. `1.0e-36`)

### 4.4 JSON config field reference for `model.type = "migdal"`

In a scan or one-point config, the `model` block uses these fields:

```json
"model": {
  "type": "migdal",
  "material": "Si",               // kept for documentation; not used in filename
  "target_nucleus": "Si28",       // fills {target_nucleus} in filename template
  "nuclear_A": 28,                // mass number (default 28 for Si)
  "nuclear_Z": 14,                // atomic number (default 14 for Si)
  "mediator": "heavy",            // "heavy" | "light" — DM-nucleus mediator
  "rates_dir": "data/migdal_rates/Si/heavy",
  "filename_template": "dRdE_{target_nucleus}_{mediator}_m{mchi_MeV}_s{sigma_n_cm2}.csv",
  "Emin_eV": 0.0,
  "Emax_eV": 20.0,
  "nbins": 200,
  "example_point": {
    "mchi_MeV": 100.0,
    "sigma_n_cm2": 1e-36
  },
  "grid": {
    "mchi_MeV": { "logspace": {"start": 10.0, "stop": 1000.0, "num": 200} },
    "sigma_n_cm2": { "logspace": {"start_exp": -42, "stop_exp": -30, "num": 200} },
    "format": { "mchi": ".6f", "sigma": ".1e" }
  }
}
```

Differences from `dm_electron`:
- `target_nucleus`, `nuclear_A`, `nuclear_Z` — new fields, parsed by ConfigManager
- `sigma_n_cm2` axis replaces `sigma_e_cm2`; the scan app resolves this via
  `resolve_axis({"sigma_n_cm2", "sigma_e_cm2"}, ...)` at ~line 1281

### 4.5 Remaining steps (to be implemented)

**Step 1 — Rate generation configs**

```
configs/migdal_generate_si_heavy.json
configs/migdal_generate_si_light.json
```

Tell the generator: mass grid, sigma_n grid, energy binning, output directory,
filename template, and darkelf_dir.

**Step 2 — Generate rate tables**

```bash
python3 utils/migdal_generate_grid.py configs/migdal_generate_si_heavy.json
python3 utils/migdal_generate_grid.py configs/migdal_generate_si_light.json
```

Output: `data/migdal_rates/Si/heavy/dRdE_Si28_heavy_m{...}_s{...}.csv`
         `data/migdal_rates/Si/light/dRdE_Si28_light_m{...}_s{...}.csv`

**Step 3 — One-point scan configs**

```
configs/migdal_one_point_si_heavy.json
configs/migdal_one_point_si_light.json
```

Point at a single pre-generated rate CSV.  Run with
`ccdarksens_example_one_point_pattern` to verify the full pipeline produces
a finite limit.

**Step 4 — Full grid scan configs**

```
configs/migdal_scan_si_heavy.json
configs/migdal_scan_si_light.json
```

Same structure as the DM-electron scan configs but with `model.type="migdal"`,
`sigma_n_cm2` axis, and pointing at the Migdal rates directory.

**Step 5 — Plot app: Migdal mode** (future)

`apps/ccdarksens_plot_dmelectron_limit.cc` currently labels the y-axis σ_e.
For Migdal:
- `--migdal` CLI flag to enter Migdal mode
- Y-axis label: `σ_n [cm²]` (DM-nucleon, per nucleon)
- Literature curves: CRESST-III, CDEX-1B, SuperCDMS Migdal limits
- Theory targets on the σ_n axis (freeze-out for heavy, freeze-in for light)
This is deferred until the scan results exist.

### 4.6 Mass range rationale

For heavy mediator Migdal, the signal is suppressed at low DM mass by the
kinematic factor v_min(ω)/v_esc.  Sensitivity typically begins around
m_χ ~ 10–50 MeV and peaks near 100 MeV – 1 GeV, where the coherent A²
enhancement compensates the reduced Migdal probability.  The grid spans:

- `mchi_MeV`: 10 – 1000 MeV (log-spaced, 200 points)
- `sigma_n_cm2`: 1e-42 – 1e-30 cm² (log-spaced, 200 points)

### 4.7 Verification checklist

1. **Unit check**: generate one CSV at (100 MeV, 1e-36 cm², heavy), confirm
   total integrated rate ≃ 3.5×10⁴ events/kg/yr (matches §1.3 reference value).

2. **Shape check**: dR/dω should peak at ω ~ E_gap (1.11 eV) and fall off at
   high ω; the spectrum should be continuous (no boxcar structure, unlike
   dark photon absorption).

3. **One-point pipeline**: run `ccdarksens_example_one_point_pattern` on
   `migdal_one_point_si_heavy.json` → confirm `S_true_ne` histogram is
   non-zero for n_e = 1, 2.

4. **Sigma linearity**: confirm that the (1e-36 cm²) and (1e-37 cm²) CSVs
   have rates differing by exactly 10× (fast-path correctness).

5. **Regression**: confirm `model.type = "dm_electron"` configs produce
   identical results before and after (no code path interference).

6. **Sensitivity order-of-magnitude**: the expected 90% CL limit on σ_n for
   1 kg-yr exposure should be O(10⁻³⁶) cm² at m_χ ~ 100 MeV — consistent
   with published Migdal sensitivities for comparable Si detectors.

### 4.8 Scope boundaries

**In scope:**
- Si target, heavy and light mediator
- Free-nucleus approximation (`approximation="free"`)
- ELF-based Migdal (`method="grid"`, using `Si_mermin.dat`)
- n_e observable (same detector pipeline as DM-electron)

**Out of scope (explicitly deferred):**
- Impulse approximation (`approximation="impulse"`)
- Ibe et al. atomic Migdal (`method="Ibe"`)
- Elastic nuclear recoil channel (`dRdEn_nuclear`) — different observable (E_n)
- SrCd₂Sb₂ target — darkelf has no ELF data for this material
- Non-Si targets
- Plot app Migdal mode — deferred until scan results exist

---

## 5. File inventory

### New files to be created

| File | Purpose |
|---|---|
| `configs/migdal_generate_si_heavy.json` | Generator config for heavy mediator |
| `configs/migdal_generate_si_light.json` | Generator config for light mediator |
| `configs/migdal_one_point_si_heavy.json` | One-point scan config, heavy |
| `configs/migdal_one_point_si_light.json` | One-point scan config, light |
| `configs/migdal_scan_si_heavy.json` | Full grid scan config, heavy |
| `configs/migdal_scan_si_light.json` | Full grid scan config, light |
| `data/migdal_rates/Si/heavy/dRdE_Si28_heavy_m*.csv` | Rate tables (generated) |
| `data/migdal_rates/Si/light/dRdE_Si28_light_m*.csv` | Rate tables (generated) |

### Existing files (already complete, no changes)

| File | Role |
|---|---|
| `python/ccdarkphys/migdal/entry.py` | Python rate entry point |
| `include/ccdarksens/model/DMNucleonConfig.hh` | DM-nucleon config struct |
| `include/ccdarksens/model/MigdalModel.hh` | C++ model header |
| `src/model/MigdalModel.cc` | C++ model implementation |
| `utils/migdal_generate_grid.py` | Grid generator driver |
| `apps/ccdarksens_scan_dmelectron_pattern.cc` | Scan app (dispatch done) |
| `CMakeLists.txt` | Build system (MigdalModel.cc included) |

---

## 6. References

- **DAMIC-M PRL 2025**: "Probing Benchmark Models of Hidden-Sector Dark Matter with DAMIC-M",
  *Physical Review Letters* **135**, 071002 (2025). DOI: 10.1103/2tcc-bqck.
  Published 13 August 2025; corrected 19 December 2025.
  Limits repository: https://github.com/DAMIC-M/LBC_2025_HSBenchmark_Limits
  → Reports Migdal limits (End Matter): most stringent for DM masses **1–35 MeV/c²** using DarkELF.
  → Caveat: "Migdal effect not yet experimentally calibrated."

- Ibe, Nakayama, Shigeki, Yanagida (2017) arXiv:1707.07258 — atomic Migdal calculation
- Knapen, Kozaczuk, Lin (2021) arXiv:2104.12785 — darkelf paper, ELF-based Migdal
- Essig, Pradler, Sholapurkar, Yu (2019) arXiv:1911.09955 — Migdal effect for direct detection review
- Dreyer, Essig, Fernandez-Serra, Singal, Zhen, PRD 109, 115008 (2024) — QCDark (used for DM-e⁻ in PRL 2025)
- darkelf source: `/Users/diegovenegasvargas/Documents/Software/DarkELF`
  (installed as editable package on Python 3.11.9)
