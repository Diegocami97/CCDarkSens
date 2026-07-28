# CCDarkSens framework: slide organization + PyDME rate generators
*(Living document / slide script for collaboration talks.)*

## Goals and audience
- **Audience**: collaborators who know DM–electron searches but not the CCDarkSens codebase.
- **Goal**: explain what the framework computes, how it’s structured, and how to run/extend it.
- **Core takeaway**: a reproducible, configuration-driven pipeline from **physics inputs → detector response → likelihood → limits**.

---

## Suggested slide deck (titles, bullets, speaker notes)

### Slide 1 — Title & one-line summary
- **Title**: *CCDarkSens: a configurable framework for DM–electron sensitivity and limits*
- **One-liner**: “From microphysics to projected limits via modular detector-response + statistics pipeline.”
- **What will be covered**: architecture, components, workflows, validation, outputs.

**Speaker notes**
- Frame the framework as: “standardized assumptions + reproducible scans + modular swap-in pieces.”

---

### Slide 2 — Physics question and deliverables
- **Inputs** (examples): \(m_\chi\), \(\sigma_e\) / coupling, mediator model (heavy/light), halo parameters, target material, exposure, thresholds/ROI, efficiencies, backgrounds, nuisance priors.
- **Outputs**:
  - differential rates / binned yields
  - profile-likelihood curves (diagnostics)
  - expected limits / discovery reach (Asimov and/or toys)
  - plots and artifacts (PDF/ROOT/CSV)

**Speaker notes**
- Clarify terminology you use internally: “sensitivity” vs “upper limit” vs “discovery.”

---

### Slide 3 — End-to-end pipeline (one diagram slide)
- **Pipeline**:
  1. **DM–electron rate model** (microphysics + astrophysics + material)
  2. **Detector response** (ionization/transport/diffusion/pattern/cluster, etc.)
  3. **Background model(s)**
  4. **Binning / observable space** (e.g. \(n_e\), pattern space)
  5. **Likelihood + test statistic**
  6. **Scan** over parameter grid → **limit curve(s)**

**Speaker notes**
- Emphasize: each step is modular and mostly controlled by JSON config.

---

### Slide 4 — Repository map (where things live)
- **Core library (C++)**: detector/response/statistics machinery.
- **Apps (C++)**: entry points for scans and plotting.
- **Configs (JSON)**: single source of truth for a run.
- **Python utilities**: rate generation (QEdark-style), helpers, plotting.
- **Data**: tables (rates, efficiencies, probabilities), reference curves.

**Speaker notes**
- Stress that analyses should be reproducible from a config + commit hash + environment.

---

### Slide 5 — Config-driven runs (how analyses are specified)
- **Principle**: “config describes the analysis; code is the engine.”
- **Typical config sections**:
  - `run`: minimizer/stat options, UL settings, output directory
  - `detector`, `experiment`: geometry/exposure/binning/observable choice
  - `response`: response-mode choice + response parameters
  - `backgrounds`: background parameterization and sources
  - `model`: signal model type + rate-table inputs (CSV dirs/templates) + scan grid

**Speaker notes**
- If helpful, name one canonical config you usually show (e.g. `configs/scan_dmelectron_ne_pydme_minuit.json`).

---

### Slide 6 — Rate engine (DM–electron)
- **Astrophysics**: SHM halo integral \( \eta(v_{\min}) \) and conventions/units.
- **Material response**: crystal form factors / probabilities (e.g. Si tables).
- **Rate tables**: precomputed/interpolated inputs to keep scans fast.
- **Cross-checks**: ensure rates match reference frameworks (e.g. PyDME/QEdark).

**Speaker notes**
- Keep equations minimal here; push details to appendix if needed.

---

### Slide 7 — Detector response modules (high level)
- **Charge/ionization**: convert energy deposition to \(n_e\) (probability tables).
- **Transport + diffusion**: maps “true charge” to reconstructed observables.
- **Classification spaces**:
  - \(n_e\)-space (single-pixel charge)
  - pattern space (multi-pixel topology categories)
- **Efficiencies**: trigger/reconstruction/selection/pattern efficiency.

**Speaker notes**
- Present response as a pipeline where you can choose “fast scan mode” vs “more detailed mode.”

---

### Slide 8 — Background modeling
- **Background sources**: dark current + flat components (and/or templates).
- **Template approaches**: pydme-style \(B = B_p + \theta B_r\) where applicable.
- **Nuisances**: prior constraints and profiling options.

**Speaker notes**
- Clarify which parameters are floated vs fixed in the typical scan.

---

### Slide 9 — Statistics model and test statistic
- **Likelihood**: Poisson per bin/category (pattern bins or \(n_e\) bins).
- **Profiling**: nuisance minimization (Brent vs Minuit2 vs “PyDME-style” 2D).
- **Test statistic**: PLR (profile likelihood ratio) with Asimov median (and toys if used).

**Speaker notes**
- Mention that there is explicit effort to match PyDME’s likelihood/minimization conventions in specific modes.

---

### Slide 10 — Scan mechanics and performance
- **Grid scan**: over \(m_\chi\) and \(\sigma_e\) (or coupling).
- **Caching**: reuse signal templates across \(\sigma\) values when possible.
- **Artifacts**: per-mass diagnostic profile-likelihood plots + summary limit curves.

**Speaker notes**
- Give a “typical runtime” example if you have it (even rough).

---

### Slide 11 — Example workflow (end-to-end)
- Choose config (e.g. “NE Asimov + PyDME-style profiling”).
- Run scan app → produce outputs (q-grid, UL curve, diagnostics).
- Plot limit curve and inspect profile-likelihood PDFs at representative masses.

**Speaker notes**
- Keep this concrete: show exact file names and what to look at in outputs.

---

### Slide 12 — Validation and cross-checks
- Reproduce reference curves (previous limits; PyDME comparisons).
- Stability checks: binning, interpolation density, nuisance bounds, boundary minima.
- Variation studies: halo parameters, efficiency variants, background assumptions.

**Speaker notes**
- Show at least one “sanity plot” you trust (e.g. profile-likelihood curve at \(m_\chi=10\) MeV).

---

### Slide 13 — Extending the framework (developer-facing)
- Add a new background component (new parameterization + config hook).
- Add/modify a response stage (pattern MC, diffusion, classifier, etc.).
- Add a new rate-table source or material model.
- Add a new observable space or binning scheme.

**Speaker notes**
- Emphasize keeping changes config-driven and auditable.

---

### Slide 14 — Summary and next steps
- **Summary**:
  - modular pipeline: physics → detector → stats
  - config-driven reproducibility
  - supports fast scans + detailed studies
- **Next steps**: planned features, pending validations, contributions requested.

---

## PyDME framework: rate-generator types available (from vendored source)
This list is based on the PyDME copy vendored in this repo at `collab_frameworks/pydme/`.

### 1) QEdark-based DM–electron scattering rates: `QEDark4DM`
- **Where**: `collab_frameworks/pydme/pydme/mqedark/qedark4dm.py`
- **Class**: `QEDark4DM(mX, xsec_e, FDMn, ...)`
- **What it generates**:
  - differential energy rates \(dR/dE\) in **events/kg/year** internally
  - convenience outputs in **events/g/day** for saved CSVs
  - conversions from energy-space to:
    - **ionized-electron yields** (per \(n_e\)) via `get_ne_from_energy_rates(...)`
    - **pattern-space rates** via `get_pattern_rates(...)`
    - **post-diffusion observed-charge rates** via `get_diffused_rates(...)`
- **DM form factor choices** (via `FDMn`):
  - `0`: heavy mediator (\(F_\mathrm{DM}=1\))
  - `1`: dipole-like (\(\propto 1/q\)) (supported by code path)
  - `2`: light mediator (\(\propto 1/q^2\))
- **Halo / flux input modes** (how \(\eta\) is computed):
  - **Analytic SHM**: `get_rates_SHM(...)` requires `self.halo == "SHM"`
  - **Verne fv-files (daily modulation / directional info)**: `get_rates_from_fvfiles(...)` reading “fv” files (gamma-dependent speed distributions)
  - **SRDM / DaMaSCUS(-SUN) style flux files**: `get_rates_from_fvfiles(..., SRDMFiles=True, ...)` reading SRDM files via `get_data_from_SRDMFile(...)`

### 2) QCDark-based DM–electron scattering rates: `QCDark4DM`
- **Where**: `collab_frameworks/pydme/pydme/mqedark/qcdark4dm.py`
- **Class**: `QCDark4DM(mX, xsec_e, FDMn, material="Si"/"Ge", ...)`
- **What it generates**: same “rate product” family as `QEDark4DM` (energy rates → \(n_e\) → pattern/diffused rates), but using a QCDark-style crystal form factor stored in HDF5.
- **Targets**:
  - `material="Si"` (supported)
  - `material="Ge"` (supported)
- **Screening / dielectric models (QCDark-only)**
  - **None**: `DoScreen=False` (no screening factor)
  - **Lindhard**: `screening_method='Lindhard'` (default if screening enabled)
  - **Thomas–Fermi**: `screening_method='ThomasFermi'`
  - **Mermin**: `screening_method='Mermin'`
  These enter as a multiplicative factor \(|1/\epsilon(E,q)|^2\) in the integrand.
- **Halo / flux input modes**: same pattern as QEDark4DM
  - analytic SHM (`get_rates_SHM`)
  - Verne fv-files (`get_rates_from_fvfiles`)
  - SRDM flux files (`get_rates_from_fvfiles(..., SRDMFiles=True, ...)`)

### 3) Halo-model “rate driver” options used inside PyDME QEdark code paths
*(This is not a separate scattering generator, but it controls the halo integral used by the generators.)*

- **Where**:
  - `collab_frameworks/pydme/pydme/mqedark/qedark.py` (rate kernels)
  - `collab_frameworks/pydme/pydme/mqedark/DM_halo_dist.py` (eta implementations + SHM/SRDM switch)
- **Halo model keywords present in code**:
  - `shm`: standard halo model (via `etaSHM`)
  - `tsa`: Tsallis (via `etaTsa`)
  - `dpl`: double power law (via `etaDPL`)
  - `msw`: Mao–Strigari–Wechsler empirical model (via `etaMSW`)
  - `DailyModulation`: uses gamma-dependent interpolation (`etagamma` / `etagamma_vectorized_smart`)
  - In `DM_halo_dist.etagamma_vectorized_smart(...)`, the `halo` switch includes:
    - `'SHM'` (analytic/numeric SHM integral driver)
    - `'SRDM'` (ROOT-accelerated interpolation/integration over SRDM flux tables)

---

## Appendix (optional slides)
- **Parameter definitions**: \(\sigma_e\) convention, form factor index `FDMn`, halo parameters, ROI definitions.
- **Likelihood details**: explicit NLL and constraint terms in “PyDME-matching” mode.
- **Config glossary**: a table mapping config keys → what they change in the pipeline.
- **Troubleshooting**: common failure modes (bad boundary fits, missing rate tables, efficiency CSV mismatches).

