# Band-gap study: architecture (QCDark2, CCDarkSens, ionization)

> **Companion docs:** [band_gap_pheno_ionization.md](band_gap_pheno_ionization.md) (procedure + p100K scaling), [band_gap_study_roadmap.md](band_gap_study_roadmap.md) (what to run next), [band_gap_pheno_p100K_scaling_explained.md](band_gap_pheno_p100K_scaling_explained.md), [qcdark2_dielectric_workflow.md](qcdark2_dielectric_workflow.md).

This document records **where each parameter lives in code**, what QCDark2 actually does (verified against `/Users/diegovenegasvargas/Documents/Software/QCDark2/qcdark2/`), and how that maps to the CCDarkSens pheno study.

---

## 1. Two independent physics layers

| Layer | Physics | Knob in study | Affects |
|-------|---------|---------------|---------|
| **A — Bulk / bands** | RPA dielectric ε(ω, **q**) from DFT + scissor | **`scissor_bandgap`** when building HDF5 | **`dR/dE(E)`** |
| **B — Ionization** | How deposited energy becomes \(n_e\) | **Scaled p100K table** (`E_gap`, `eh` in builder) | **`S(n_e)`**, limits |

Changing **only** layer B at fixed `rates_dir` isolates ionization. Changing layer A requires a **new ε HDF5 + rate grid** per scissor gap.

**There is no `eh_pair_eV` in QCDark2 ε or `get_dR_dE`.** The name `eh_pair_eV` in CCDarkSens JSON is bookkeeping + input to **`build_p100K_scaled.py`**.

---

## 2. QCDark2 source (external package)

Path on your machine: `/Users/diegovenegasvargas/Documents/Software/QCDark2/qcdark2/`

### 2.1 Dielectric generation — band gap only

| Item | Location |
|------|----------|
| Input keyword | `scissor_bandgap` in `.in` files (`materials/Si.in`, `configs/qcdark2/Si_fast_scissor.in`) |
| Parsed | `dielectric_pyscf/input_parameters.py` |
| Applied | `dielectric_pyscf/dft_routines.py` → `convert_to_eV_and_scissor()` |
| Effect | Shifts **all conduction bands** so fundamental gap = `scissor_bandgap` |
| Output | `<name>_resources/epsilon.hdf5` → packaged as `Si_fast_gap*.h5` |

**Not used in ε:** `eh_pair_eV`, `band_gap_eV` (CCDarkSens JSON names).

DFT under `save_path/DFT_resources/` is **reused** when only `scissor_bandgap` and `name` change (same lattice, basis, k-grid).

### 2.2 DM rates — `dR/dE` only

| Item | Location |
|------|----------|
| Entry | `dark_matter_rates.py` → `get_dR_dE(epsilon, m_X, mediator, astro_model, screening, velocity_dist)` |
| Inputs | Loaded HDF5 (`epsilon`, `q`, `E`, `M_cell`, `V_cell`, `dE`), \(m_\chi\), mediator, halo, screening |
| Output | `dR_dE` [events / kg / year / eV], energy array `E` |

**Not passed to `get_dR_dE`:** any ionization parameter, `eh_pair_eV`, `band_gap_eV`.

### 2.3 Ionization inside QCDark2 (optional, **not** used by CCDarkSens scans)

| Item | Location |
|------|----------|
| Function | `dark_matter_rates.recoil_spectrum(dR, ionization_file=...)` |
| Default file | `../secondary_ionization/p100K.dat` (path to table file) |
| Mechanism | Reads **columns** `E` + `pair_creation_prob[n]`; **no** scalar `eh_pair_eV` |
| Used in | `examples/DM_calculations.ipynb` |

CCDarkSens uses **`ChargeIonization` + CSV** (`p100K_table.csv` or scaled `p100K_gap*_eh*.csv`) instead of `recoil_spectrum`.

### 2.4 Summary table: does QCDark2 use `eh_pair_eV`?

| Code path | `scissor_bandgap` / gap | `eh_pair_eV` |
|-----------|-------------------------|--------------|
| `dielectric_pyscf` | ✅ | ❌ |
| `get_dR_dE` | ❌ (in ε already) | ❌ |
| `recoil_spectrum` | ❌ | ❌ (uses full table **file**) |
| CCDarkSens `qcdark2/entry.py` | ❌ (metadata only) | ❌ (metadata only) |

---

## 3. CCDarkSens wiring

### 3.1 Rate generation

`python/ccdarkphys/qcdark2/entry.py` calls `get_dR_dE` only.  
`band_gap_eV` and `eh_pair_eV` from JSON → **`meta`** in rate CSV headers only.

### 3.2 Scans

`response.charge_ionization` in scan JSON (`ConfigManager` + `ccdarksens_scan_dmelectron_pattern.cc`):

```json
"charge_ionization": {
  "table_csv": "data/p100K_gap0p1_eh0p1.csv",
  "band_gap_eV": 0.1,
  "eh_pair_eV": 0.1,
  "scenario": "D-equal"
}
```

Only **`table_csv`** changes physics; other fields are labels / consistency checks.

### 3.3 p100K rescaling (layer B)

`utils/build_p100K_scaled.py` — anchored map ([band_gap_pheno_p100K_scaling_explained.md](band_gap_pheno_p100K_scaling_explained.md)):

- Extended grid **0.05–50 eV** (0.05 eV steps)
- `E'(E) = E_gap_ref + (E - E_gap_new) × (eh_ref / eh_new)`
- Below `E_gap_new`: \(P = 0\)

---

## 4. What to regenerate per scenario

| Change | New ε HDF5? | New rate grid? | New p100K CSV? |
|--------|-------------|----------------|----------------|
| Scissor gap 0.1 → 0.7 eV | **Yes** | **Yes** | Optional (match gap in builder) |
| Same gap, B-thresh → D-equal | No | No | **Yes** |
| JSON `eh_pair_eV` only | No | No | No (unless you rebuild table) |

**Do not** run 16 dielectric jobs for 4 gaps × 4 eh values. Run **one ε + one rate grid per gap**, then **several ionization tables per gap**.

---

## 5. Pheno scenario IDs (ionization)

| ID | `eh_pair_eV` | Role |
|----|--------------|------|
| **ref** | 3.8 | `p100K_table.csv` (Si reference) |
| **B-thresh** | 3.8 | Low gap, Si-like multiplicity |
| **D-equal** | = gap | Equal-scale pheno (aggressive) |
| **A-ratio** | \(3.8 \times \mathrm{gap}/1.2\) | Keeps Si ratio \( \varepsilon_h / E_\mathrm{gap} \) |

---

## 6. Plotting & validation tools

| Script | Purpose |
|--------|---------|
| `utils/build_p100K_scaled.py` | Build scaled CSVs; `--from-manifest`, `--force`, `--E-min 0.05` |
| `utils/eval_p100K_scaling.py` | Console QA + `eval_*.pdf` |
| `utils/plot_p100K_scaling_compare.py` | Overlays + **`Pne_fan/`** (all \(P(n_e\mid E_r)\), PRD-style) |
| `utils/verify_p100K_scaled.py` | Grid / anchor / high-E checks |
| `utils/plot_band_gap_pheno_dRdE.py` | Compare `dR/dE` across gaps |

---

## 7. Methods wording (short)

- **Rates:** phenomenological bulk Si DFT + **scissor** conduction-band shift to effective gap \(E_\mathrm{gap}^\mathrm{scissor}\).
- **Ionization:** PRD 100 K Si \(P(n_e\mid E)\) **rescaled** with independent \((E_\mathrm{gap}^\mathrm{ion}, \varepsilon_h)\); not new ab initio ionization MC.
- **Do not claim** QCDark2 used `eh_pair_eV` in ε or `dR/dE`.

---

## 8. File manifest

`configs/band_gap_pheno_scenarios.json` — scenarios + `gap_sweep` list (gaps, paths, status).

See [band_gap_study_roadmap.md](band_gap_study_roadmap.md) for phased execution including **0.7 eV** and **0.9 eV** scissor points.
