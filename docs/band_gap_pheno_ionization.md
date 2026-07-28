# Band-gap pheno study: `eh_pair_eV`, `band_gap_eV`, and `p100K_table.csv`

> **Related:** [band_gap_study_architecture.md](band_gap_study_architecture.md) (QCDark2 vs CCDarkSens — **read this first**), [band_gap_study_roadmap.md](band_gap_study_roadmap.md) (execution plan incl. **0.7 & 0.9 eV**), [qcdark2_dielectric_workflow.md](qcdark2_dielectric_workflow.md), [band_gap_pheno_p100K_scaling_explained.md](band_gap_pheno_p100K_scaling_explained.md).

This document defines the **full pheno procedure**: couple scissor-shifted QCDark2 rates with **configurable** ionization tables, document that **\(E_\mathrm{gap}\) and \(\varepsilon_h\) are independent choices** (including **equal-scale** pheno, e.g. 0.1 eV for both), and lay out **step-by-step work with a figure at every stage** so you can see what changes.

**Status:** Tools and `charge_ionization` wiring implemented. Gap sweep extended to **0.1–1.2 eV** including **0.7, 0.9** (see roadmap). Full scans: in progress.

---

## 1. Problem statement

The scissor-only band-gap study changes **only** `dR/dE`. Ionization is frozen to `data/p100K_table.csv` (Si, turn-on ~**1.2 eV**, pair scale ~**3.8 eV**):

| Layer | Knob | Current behavior |
|-------|------|------------------|
| QCDark2 ε | `scissor_bandgap` in HDF5 | ✅ varies per `Si_fast_gap*.h5` |
| Rate JSON | `band_gap_eV`, `eh_pair_eV` | ❌ not used by QCDark2 `get_dR_dE` (metadata only in CSV) |
| Scan ionization | `response.charge_ionization.table_csv` | ✅ per-scenario scaled p100K |

**Goal:** for each pheno scenario, align (1) `epsilon.h5`, (2) `rates_dir`, (3) `charge_ionization.table_csv`, and (4) documented `(band_gap_eV, eh_pair_eV)` — and **plot after each stage**.

---

## 2. Physics: three quantities (no fixed link between gap and \(\varepsilon_h\))

| Quantity | Role | Reference Si |
|----------|------|----------------|
| **`band_gap_eV`** | Threshold for producing charge / cutting `dR/dE` (QEdark) or turn-on in p100K | ~1.2 eV |
| **`eh_pair_eV`** | Mean energy per e–h pair; sets \(n_e\) multiplicity scale in p100K rescaling | ~3.8 eV |
| **`p100K_table.csv`** | Tabular \(P(n_e \mid E_r)\), PRD 102.063026 (100 K Si) | encodes ~1.2 / ~3.8 |

**Important:** there is **no known fundamental relation** \(\varepsilon_h = f(E_\mathrm{gap})\) for arbitrary low-gap pheno. Real Si has \(\varepsilon_h \gg E_\mathrm{gap}\). You **choose** \((E_\mathrm{gap}, \varepsilon_h)\) per scenario and document it.

---

## 3. Pheno scenarios: how to choose \((E_\mathrm{gap}, \varepsilon_h)\)

Use a **scenario ID** in configs, plot labels, and output paths.

| ID | `band_gap_eV` (scissor) | `eh_pair_eV` (ionization table) | Interpretation |
|----|-------------------------|----------------------------------|----------------|
| **ref** | 1.2 | 3.8 | Reference Si (`p100K_table.csv`, unscaled) |
| **B-thresh** | 0.1 (or 0.3, 0.5) | **3.8** | Lower **threshold only**; Si-like pair scale |
| **D-equal** | 0.1 (or 0.3, 0.5) | **same as gap** (e.g. 0.1 & 0.1) | **Equal-scale pheno** — aggressive, not bulk Si |
| **A-ratio** | 0.1, 0.3, 0.5 | \(3.8 \times E_\mathrm{gap}/1.2\) | Keeps Si **ratio** \(\varepsilon_h/E_\mathrm{gap}\) (optional default) |
| **C-grid** | user | user | Independent 2D grid |

**Recommended first comparison at scissor 0.1 eV** (same `Si_fast_gap0p1` rates, different ionization only):

| Scenario | `ionization_csv` | What you learn |
|----------|------------------|----------------|
| `0p1_eh3p8` | `p100K_gap0p1_eh3p8.csv` | Effect of opening threshold with Si-like \(\varepsilon_h\) |
| `0p1_eh0p1` | `p100K_gap0p1_eh0p1.csv` | Effect of **also** compressing pair scale to 0.1 eV |

Then repeat **D-equal** (or **A-ratio**) for 0.3, 0.5, 1.2 eV scissor points as needed.

---

## 4. End-to-end pipeline

```mermaid
flowchart LR
  S0[Step 0 Layout]
  S1[Step 1 dRdE vs E]
  S2[Step 2 p100K Pne]
  S3[Step 3 S_ne fold]
  S4[Step 4 Scan]
  S5[Step 5 Limits]
  S0 --> S1 --> S2 --> S3 --> S4 --> S5
```

**Figure output root (create once):**

```text
outplots/band_gap_pheno/
  step1_dRdE/
  step2_p100K/
  step3_Sne/
  step4_scan_diag/
  step5_limits/
```

---

## 5. Step-by-step procedure (with figures)

Work from repo root (`CCDarkSens/`). Each step produces at least one PDF/PNG under `outplots/band_gap_pheno/`.

### Step 0 — Layout and scenario table

**Do**

1. Create `outplots/band_gap_pheno/` subdirs (above).
2. Fill `configs/band_gap_pheno_scenarios.json` (see §11) with scenarios you will run.
3. Confirm ε and rates exist:

```bash
ls data/qcdark2_epsilon/Si/Si_fast_gap*.h5
ls data/qcdark2_rates/Si/heavy/Si_fast_gap0p1/ | head
```

**Figure:** none (inventory only).

**Check:** HDF5 attrs `scissor_bandgap_eV` match intended gap (0.1, 0.3, 0.5, 1.2).

---

### Step 1 — Dielectric / differential rates `dR/dE(E)`

**What changes:** scissor gap shifts where `dR/dE` is nonzero and its shape at low \(E\).

**Do**

1. Pick one reference point: e.g. `mchi = 0.5` MeV, `sigma = 1e-26` cm² (file exists in each `Si_fast_gap*` folder).
2. Overlay `dR/dE` vs \(E\) for gaps 0.1, 0.3, 0.5, 1.2 eV.

**Commands (after `utils/plot_band_gap_pheno_dRdE.py` exists — TODO)**

```bash
python3 utils/plot_band_gap_pheno_dRdE.py \
  --mass-mev 0.5 --sigma 1e-26 \
  --rates-dirs \
    data/qcdark2_rates/Si/heavy/Si_fast_gap0p1 \
    data/qcdark2_rates/Si/heavy/Si_fast_gap0p3 \
    data/qcdark2_rates/Si/heavy/Si_fast_gap0p5 \
    data/qcdark2_rates/Si/heavy/Si_fast_gap1p2 \
  --labels "gap 0.1" "gap 0.3" "gap 0.5" "gap 1.2" \
  --out outplots/band_gap_pheno/step1_dRdE/dRdE_m0p5_s1e-26.pdf
```

**Interim (no new script):** adapt `utils/plot_qedark_rates.py` logic or a short Python one-off reading one CSV per gap.

**Figure to expect**

- Log–log or semilog \(E\) [eV] vs `dR/dE` [events/(kg·year·eV)].
- **Look for:** lowest \(E\) where each curve rises; integral at low \(E\); ratio of curves vs ref (1.2 eV).

**Save also:** integrated rate below 2 eV, 2–10 eV (table in log or CSV next to PDF).

---

### Step 2 — Ionization tables `P(n_e | E)`

**What changes:** turn-on energy and how fast \(P(n_e{>}1|E)\) grows with \(E\) (controlled by `eh_pair_eV` in rescaling).

**Do**

1. Build tables (TODO: `utils/build_p100K_scaled.py`):

```bash
# Reference (copy or identity)
cp data/p100K_table.csv data/p100K_gap1p2_eh3p8.csv

# Equal-scale examples
python3 utils/build_p100K_scaled.py \
  --base data/p100K_table.csv \
  --band-gap-eV 0.1 --eh-pair-eV 0.1 \
  --out data/p100K_gap0p1_eh0p1.csv \
  --plot outplots/band_gap_pheno/step2_p100K/Pne_gap0p1_eh0p1.pdf

# Threshold-only at 0.1 eV gap, Si-like eh
python3 utils/build_p100K_scaled.py \
  --base data/p100K_table.csv \
  --band-gap-eV 0.1 --eh-pair-eV 3.8 \
  --out data/p100K_gap0p1_eh3p8.csv \
  --plot outplots/band_gap_pheno/step2_p100K/Pne_gap0p1_eh3p8.pdf
```

Repeat for 0.3, 0.5 eV as needed (equal and/or 3.8 eh).

2. Overlay plots for \(P(n_e{=}1\mid E)\) and optionally \(\sum_{n\ge1} P(n\mid E)\).

**Commands (plot-only — TODO or use `--plot` on builder)**

```bash
python3 utils/plot_band_gap_pheno_p100K.py \
  --tables \
    data/p100K_table.csv \
    data/p100K_gap0p1_eh3p8.csv \
    data/p100K_gap0p1_eh0p1.csv \
  --labels "Si ref" "0.1/3.8 thresh" "0.1/0.1 equal" \
  --out outplots/band_gap_pheno/step2_p100K/Pne_compare_0p1_scenarios.pdf
```

**Figures to expect**

| Plot | What to look for |
|------|------------------|
| \(P(n_e{=}1\mid E)\) vs \(E\) | Turn-on at new `band_gap_eV`; ref turns on ~1.2 eV |
| \(P(n_e{=}2\mid E)\) vs \(E\) | **Equal eh** → multi-\(n_e\) structure compressed to lower \(E\) |
| Heatmap \(P(n_e\mid E)\) | Optional: 2D view per scenario |

**Key check:** for `gap0p1_eh0p1`, confirm \(P>0\) for \(E \gtrsim 0.1\) eV. For **ref p100K**, \(P=0\) below 1.2 eV.

---

### Step 3 — Folded signal `S(n_e)` at one DM point

**What changes:** same `dR/dE`, different ionization → different \(S_{n_e}\). Isolates ionization from rate shape.

**Do**

1. Use `ccdarksens_example_one_point_pattern` (or minimal scan at one \((m_\chi,\sigma)\)) with:
   - **Same** `rates_dir` = `Si_fast_gap0p1`
   - **Different** `response.charge_ionization.table_csv` per scenario
2. Dump `S_pat_validation` / `S_ne` from ROOT or app stdout.

**Configs (TODO):** e.g. `configs/example_band_gap_pheno_0p1_eh0p1.json` pointing at one rate point + `p100K_gap0p1_eh0p1.csv`.

**Figure**

- Bar chart or line: \(S_{n_e}\) or \(S_\mathrm{pat}\) per bin for scenarios `ref`, `0p1_eh3p8`, `0p1_eh0p1`.
- Save: `outplots/band_gap_pheno/step3_Sne/Sne_m0p5_gap0p1_scenarios.pdf`

**Look for**

- With **only** scissor rates but **ref** p100K: low bins may still be **zero** (p100K blocks \(E<1.2\) eV).
- With `0p1_eh3p8` or `0p1_eh0p1`: **nonzero** low-\(n_e\) signal if `dR/dE` has weight below 1.2 eV.

---

### Step 4 — Full grid scan (limits input)

**What changes:** profile likelihood uses \(S_{n_e}\) + background across full \((m_\chi,\sigma)\) grid.

**Do**

1. Wire `response.charge_ionization` in scan JSON (TODO in C++).
2. One scan per **(scissor gap, ionization scenario)** you care about, e.g.:
   - `outputs/scan_band_gap_0p1_eh0p1/`
   - `outputs/scan_band_gap_0p1_eh3p8/`
3. Keep exposure, background, ROI **identical** to `scan_dmelectron_band_gap_study_*.json`.

```bash
build/ccdarksens_scan_dmelectron_pattern \
  configs/scan_band_gap_pheno_0p1_eh0p1.json
```

**Diagnostic figure (mid-scan)**

- For one mass: `q(m_\chi, \sigma)` heatmap or slice at fixed \(m_\chi\) — compare ref vs new ionization.
- Save: `outplots/band_gap_pheno/step4_scan_diag/q_slice_m10_gap0p1.pdf` (from ROOT `q_mchi_sigma` histogram).

---

### Step 5 — Limit curves \(\sigma_\mathrm{UL}(m_\chi)\)

**What changes:** final observable for collaboration plots.

**Do**

```bash
build/ccdarksens_plot_dmelectron_limit \
  outputs/scan_band_gap_0p1_eh0p1/scan_dmelectron_pattern.root "0.1 eV gap, 0.1 eh equal" \
  outputs/scan_band_gap_0p1_eh3p8/scan_dmelectron_pattern.root "0.1 eV gap, 3.8 eh" \
  outputs/scan_dmelectron_band_gap_study_1p2/scan_dmelectron_pattern.root "1.2 eV ref" \
  heavy --from-qhist --batch \
  --title "Band-gap pheno: ionization scenarios" \
  --out-pdf outplots/band_gap_pheno/step5_limits/limit_compare_0p1_ionization.pdf \
  --out-csv outplots/band_gap_pheno/step5_limits/limit_compare_0p1_ionization.csv
```

**Figures**

| Output | Purpose |
|--------|---------|
| Overlay PDF | Limit curves |
| CSV | Contours for further plotting |
| Ratio plot (TODO script) | \(\sigma_\mathrm{UL}(m_\chi) / \sigma_\mathrm{UL}^\mathrm{ref}(m_\chi)\) vs \(m_\chi\) |

```bash
python3 utils/plot_band_gap_pheno_ratio.py \
  --csv outplots/band_gap_pheno/step5_limits/limit_compare_0p1_ionization.csv \
  --reference "1.2 eV ref" \
  --out outplots/band_gap_pheno/step5_limits/ratio_vs_ref.pdf
```

---

### Step 6 — Full gap sweep (optional)

Repeat Steps 1–5 for scissor gaps **0.1, 0.3, 0.5, 1.2** with your chosen ionization rule (**D-equal** or **A-ratio** or **B-thresh**).

**Summary figure:** grid of limit curves or ratio-to-ref for all gaps — `outplots/band_gap_pheno/step5_limits/limit_all_gaps_equal_eh.pdf`.

---

## 6. Building scaled `p100K` tables

**Full explanation with examples:** [band_gap_pheno_p100K_scaling_explained.md](band_gap_pheno_p100K_scaling_explained.md).

### Rescaling formula (Tier 1 pheno — anchored map)

Reference: \(E_\mathrm{gap}^\mathrm{ref}=1.2\) eV, \(\varepsilon_h^\mathrm{ref}=3.8\) eV.

\[
E'(E) = E_\mathrm{gap}^\mathrm{ref} + \bigl(E - E_\mathrm{gap}^\mathrm{new}\bigr) \times \frac{\varepsilon_h^\mathrm{ref}}{\varepsilon_h^\mathrm{new}}
\]

\[
P_\mathrm{new}(n_e \mid E) =
\begin{cases}
0 & E < E_\mathrm{gap}^\mathrm{new} \\
P_\mathrm{ref}\bigl(n_e \mid E'(E)\bigr) & \text{otherwise}
\end{cases}
\]

(linear interpolation on the reference table)

**Equal-scale case:** set `band-gap-eV` and `eh-pair-eV` to the **same** value (e.g. 0.1 and 0.1). No special code path — only the two numbers passed to the builder.

### Naming

```text
data/p100K_gap{gap}_eh{eh}.csv
# Examples:
#   p100K_gap0p1_eh0p1.csv   ← equal 0.1 / 0.1
#   p100K_gap0p1_eh3p8.csv  ← threshold 0.1, Si eh
#   p100K_gap0p3_eh0p3.csv   ← equal 0.3 / 0.3
```

---

## 7. Code and config changes

### 7.1 `utils/build_p100K_scaled.py` (done)

Output energy grid is extended to **0.05 eV** in 0.05 eV steps (21 rows below the reference table’s 1.1 eV start) so pheno turn-ons at 0.1 / 0.3 eV sit on explicit table rows.

```bash
# One scenario
python3 utils/build_p100K_scaled.py --band-gap-eV 0.1 --eh-pair-eV 0.1 --scenario D-equal --plot

# All tables in manifest (use --force to overwrite)
python3 utils/build_p100K_scaled.py --from-manifest configs/band_gap_pheno_scenarios.json --force
```

Optional: `--E-min 0.05` (default), `--force`.

### 7.2 Plot helpers (done)

| Script | Step |
|--------|------|
| `utils/plot_band_gap_pheno_dRdE.py` | 1 |
| `utils/plot_band_gap_pheno_p100K.py` | 2 |
| `utils/plot_band_gap_pheno_ratio.py` | 5 |

### 7.3 `response.charge_ionization` in JSON + `ConfigManager` (done)

```json
"response": {
  "charge_ionization": {
    "table_csv": "data/p100K_gap0p1_eh0p1.csv",
    "band_gap_eV": 0.1,
    "eh_pair_eV": 0.1,
    "scenario": "D-equal"
  }
}
```

Default: `data/p100K_table.csv` if omitted.

### 7.4 Scan configs (done)

| Config | Ionization | Rates |
|--------|------------|-------|
| `configs/scan_band_gap_pheno_0p1_eh0p1.json` | `p100K_gap0p1_eh0p1.csv` (D-equal) | `Si_fast_gap0p1` |
| `configs/scan_band_gap_pheno_0p1_eh3p8.json` | `p100K_gap0p1_eh3p8.csv` (B-thresh) | `Si_fast_gap0p1` |
| `configs/scan_band_gap_pheno_0p3_eh0p3.json` | `p100K_gap0p3_eh0p3.csv` | `Si_fast_gap0p3` |
| `configs/scan_band_gap_pheno_0p3_eh0p95.json` | `p100K_gap0p3_eh0p95.csv` (A-ratio) | `Si_fast_gap0p3` |

---

## 8. What to say in methods (wording)

| Do say | Do not say |
|--------|------------|
| Pheno study with scissor-shifted ε and **rescaled** p100K ionization; \((E_\mathrm{gap}, \varepsilon_h)\) chosen per scenario (including equal-scale) | “Fundamental \(\varepsilon_h(E_\mathrm{gap})\) relation” without evidence |
| Show **B-thresh** and **D-equal** at same scissor gap to bracket ionization uncertainty | “True 0.1 eV band gap silicon” |

---

## 9. What this does **not** claim

- Not real low-gap bulk Si; MO character fixed.
- Rescaled p100K ≠ new PRD MC.
- Pattern efficiencies, diffusion, readout unchanged.

---

## 10. Relation to scissor-only work

Already done: `Si_fast_gap*`, rate grids, `outputs/scan_dmelectron_band_gap_study_*`, `outputs/band_gap_study_compare.pdf` (modest separation — expected until Step 2–3 fix ionization).

---

## 11. Example scenario manifest

**`configs/band_gap_pheno_scenarios.json`** (template — create when implementing):

```json
{
  "reference": {
    "label": "Si ref",
    "band_gap_eV": 1.2,
    "eh_pair_eV": 3.8,
    "epsilon_h5": "data/qcdark2_epsilon/Si/Si_fast_gap1p2.h5",
    "rates_dir": "data/qcdark2_rates/Si/heavy/Si_fast_gap1p2",
    "ionization_csv": "data/p100K_table.csv",
    "scenario": "ref"
  },
  "scenarios": [
    {
      "label": "0.1 eV equal eh",
      "band_gap_eV": 0.1,
      "eh_pair_eV": 0.1,
      "epsilon_h5": "data/qcdark2_epsilon/Si/Si_fast_gap0p1.h5",
      "rates_dir": "data/qcdark2_rates/Si/heavy/Si_fast_gap0p1",
      "ionization_csv": "data/p100K_gap0p1_eh0p1.csv",
      "scenario": "D-equal",
      "scan_config": "configs/scan_band_gap_pheno_0p1_eh0p1.json",
      "outdir": "outputs/scan_band_gap_0p1_eh0p1"
    },
    {
      "label": "0.1 eV gap, Si eh",
      "band_gap_eV": 0.1,
      "eh_pair_eV": 3.8,
      "epsilon_h5": "data/qcdark2_epsilon/Si/Si_fast_gap0p1.h5",
      "rates_dir": "data/qcdark2_rates/Si/heavy/Si_fast_gap0p1",
      "ionization_csv": "data/p100K_gap0p1_eh3p8.csv",
      "scenario": "B-thresh",
      "scan_config": "configs/scan_band_gap_pheno_0p1_eh3p8.json",
      "outdir": "outputs/scan_band_gap_0p1_eh3p8"
    }
  ]
}
```

Same `rates_dir` for both 0.1 eV rows — **only ionization differs** in Step 3–5.

---

## 12. Implementation checklist

### Documentation

- [x] [band_gap_study_architecture.md](band_gap_study_architecture.md) — QCDark2 source map, `eh_pair_eV` clarification
- [x] [band_gap_study_roadmap.md](band_gap_study_roadmap.md) — phased plan + **0.7 / 0.9 eV**
- [x] Equal-scale (**D-equal**) and independent \((E_\mathrm{gap}, \varepsilon_h)\) choices
- [x] Step-by-step procedure with figure paths

### Python tools

- [x] `utils/build_p100K_scaled.py`, `eval_p100K_scaling.py`, `plot_p100K_scaling_compare.py`, `verify_p100K_scaled.py`
- [x] `utils/plot_band_gap_pheno_dRdE.py`, `plot_band_gap_pheno_ratio.py`
- [x] `configs/band_gap_pheno_scenarios.json`

### C++ / configs

- [x] `response.charge_ionization` in `ConfigManager` + scan app
- [x] `configs/scan_band_gap_pheno_0p1_*`, `0p3_*`
- [x] `configs/scan_dmelectron_band_gap_study_0p7.json`, `_0p9.json` (scissor-only baseline)
- [ ] `configs/scan_band_gap_pheno_0p7_*`, `_0p9_*` (coupled pheno — copy from 0p1 template)

### QCDark2 (per gap)

- [x] ε + rates: 0.1, 0.3, 0.5, 1.2 eV
- [ ] ε + rates: **0.7, 0.9 eV** ← next (Phase A in roadmap)

### Runs & figures

- [x] Step 2 PDF + `Pne_fan/`: `outplots/band_gap_pheno/step2_p100K/`
- [ ] Step 1 / 3 / 5 for full gap sweep after new rates exist

---

## 13. Quick command reference (full chain, equal-scale 0.1 eV example)

```bash
mkdir -p outplots/band_gap_pheno/{step1_dRdE,step2_p100K,step3_Sne,step4_scan_diag,step5_limits}

# Step 2 — ionization table (once script exists)
python3 utils/build_p100K_scaled.py \
  --band-gap-eV 0.1 --eh-pair-eV 0.1 \
  --out data/p100K_gap0p1_eh0p1.csv \
  --plot outplots/band_gap_pheno/step2_p100K/Pne_gap0p1_eh0p1.pdf

# Step 4 — scan (once config + C++ wired)
build/ccdarksens_scan_dmelectron_pattern configs/scan_band_gap_pheno_0p1_eh0p1.json

# Step 5 — limits
build/ccdarksens_plot_dmelectron_limit \
  outputs/scan_band_gap_0p1_eh0p1/scan_dmelectron_pattern.root "0.1/0.1 equal" \
  outputs/scan_dmelectron_band_gap_study_1p2/scan_dmelectron_pattern.root "ref 1.2" \
  heavy --from-qhist --batch \
  --out-pdf outplots/band_gap_pheno/step5_limits/limit_0p1_equal_vs_ref.pdf \
  --out-csv outplots/band_gap_pheno/step5_limits/limit_0p1_equal_vs_ref.csv
```

See [qcdark2_dielectric_workflow.md §14](qcdark2_dielectric_workflow.md#14-phenomenological-low-band-gap-dm-study-using-si) for dielectric generation; this file covers **ionization + validation plots**.
