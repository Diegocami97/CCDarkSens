# How pattern efficiencies flow into the pipeline — example app walkthrough

This document walks through how **signal** and **DC background** pattern efficiencies are obtained, stored, and used in the example app `ccdarksens_example_one_point_pattern` when `observable_bins = "pattern"`.

---

## 0. S_true → S_rec → S_obs → S_pat (pattern-space analysis)

In **pattern mode**, the chain is:

| Symbol | Meaning | Where it lives |
|--------|--------|----------------|
| **S_true(n_e)** | Expected signal count per **true** electron number n_e (from ionization). dR/dE × exposure folded by ionization probabilities. | Inside `DetectorResponsePipeline::Apply`: step 1 (ionization). |
| **S_rec(n_e_obs)** | Expected count per **reconstructed** n_e (after charge smearing). S_true folded by P(n_obs \| n_true) from PCD: diffusion, readout noise, DC, threshold → P(q \| n_e) → kernel. | Same pipeline: step 2 (PCD kernel fold). **No extra diffusion** is applied here — the PCD kernel already encodes it. |
| **S_obs(n_e_obs)** | Observable spectrum in n_e space. If pattern efficiency is applied in the pipeline: S_obs = ε(n_e_obs) × S_rec (bin-by-bin). If `skip_pattern_efficiency` is true (our pattern case): pipeline returns S_rec and we treat it as S_obs; efficiency is applied only in the pattern fold below. | Pipeline step 3; returned as `TH1D` to the app. |
| **S_pat** | Signal rate **per pattern bin**. Must use **S_true** (not S_obs) with ε(pattern \| n_true): S_pat[p] = Σ_{n_true} S_true(n_true) × ε(p \| n_true). | In the app: `FoldNeToPatternRates(S_true, ..., pattern_roi, pattern_eff_map)`. |

**Why S_true for the pattern fold:** The efficiency table ε(pattern \| n_e) is **P(classified as pattern \| true n_e)** — PatternMC injects **true** n_e and runs full detector (diffusion, etc.). If we folded **S_obs** (reconstructed spectrum, already smeared by PCD/diffusion) with ε(pattern \| n_true), we would apply the detector response twice: once in the PCD kernel (S_true→S_rec) and again in ε. So we **fold S_true with ε(pattern \| n_true)**; then diffusion and full detector are applied only once (inside ε). **Background** must be in true n_e as well: the flat-background component is built with **ion->FoldToNe(dRdE_flat, ...)** in pattern mode (not pipe.Apply), so B_tot is true n_e and B_pat = FoldNeToPatternRates(B_tot, ...) is correct and avoids double-counting.

**Diffusion in the pipeline:** When PCD machinery is present, the pipeline uses the PCD kernel (built from full MC that already includes diffusion in ChargeTransport/PixelSimulator). The optional `Diffusion` object set by the app is only used in the **else** branch when PCD is *not* available; in that case histogram-level diffusion is applied once. So either (PCD path → no `diff_` call) or (no PCD → one `diff_->Apply()`).

---

## 1. High-level flow

```
Config (pattern_roi, roi_bins, optional efficiency_csv)
        │
        ▼
┌─────────────────────────────────────────────────────────────────┐
│  pattern_eff_map:  (pattern_id, n_e) → ε(pattern | n_e)          │
│  Filled from: CSV file OR PatternMC::BuildPatternTable (2D/1D)   │
└─────────────────────────────────────────────────────────────────┘
        │
        ├──────────────────────────────────┬──────────────────────────────────┐
        ▼                                  ▼                                  ▼
┌───────────────────┐            ┌───────────────────┐            ┌───────────────────┐
│  B_tot(n_e)       │            │  S_true(n_e)      │            │  efficiency_per_   │
│  (background      │            │  (signal true     │            │  pattern.csv       │
│   spectrum,       │            │   n_e spectrum;  │            │  (output)          │
│   true n_e)       │            │   used for fold)  │            │                   │
└─────────┬─────────┘            └─────────┬─────────┘            └───────────────────┘
          │                                │
          │  FoldNeToPatternRates          │  FoldNeToPatternRates
          │  (same pattern_eff_map         │  (ε(pattern|n_true); avoids
          │   ε(pattern|n_true))           │   double-counting)
          ▼                                ▼
┌───────────────────┐            ┌───────────────────┐
│  B_pat            │            │  S_pat            │
│  (background      │            │  (signal rate     │
│   rate per        │            │   per pattern     │
│   pattern bin)   │            │   bin)            │
└─────────┬─────────┘            └─────────┬─────────┘
          │                                │
          └──────────────┬─────────────────┘
                         ▼
              PLR test statistic: q(B_pat, S_pat + B_pat, B_pat)
```

**Important:** The **same** efficiency table **ε(pattern | n_e)** is used for both signal and background. There is no separate “background efficiency table” in this pipeline; background is modeled as a spectrum in n_e (B_tot) and folded with the same ε. When **background_efficiency_csv** is set, DC (and full background) are folded through the background migration matrix P(identified | true pattern) from that CSV instead.

---

## 2. Step-by-step in the example app

### 2.1 Config and mode

- **observable_bins** = `"pattern"` → we work in pattern space; the app sets `use_pattern_bins = true`.
- **pattern_roi** = list of pattern IDs (e.g. 1, 2, 3, 4, 5, 11, 21, …, 311) that define the observable bins.
- **roi_bins** = list of n_e values (e.g. 1..5) used when writing the efficiency CSV and for ROI diagnostics.
- **ne_min_bkg** = `min(ne_min, 0)` so n_e = 0 (dark current) is included when folding background.

### 2.2 Filling `pattern_eff_map` (signal efficiency)

**pattern_eff_map** is a `std::map<std::pair<int,int>, double>`: key = (pattern_id, n_e), value = **ε(pattern | n_e)** = P(event classified as that pattern | n_e true electrons).

Two sources (in order):

1. **From CSV (if configured)**  
   If `response.pattern_mc.efficiency_csv` is set (e.g. to `outputs/scan_pattern/efficiency_per_pattern.csv`):
   - The app reads the file (columns: `pattern`, `ne`, `Efficiency`).
   - Each row → `pattern_eff_map[{pattern_code, ne_val}] = eff_val`.
   - This is the **signal** efficiency table (same format as produced by the example/scan app or the notebook).

2. **From PatternMC (if no CSV)**  
   If no CSV is given:
   - PatternMC is built with the classifier and (optionally) the 2D image generator.
   - **BuildPatternTable(ne_min_bkg, ne_max, Ee_ref_eV)** is triggered (via PrecomputeEpsilon or when the table is first needed).
   - For each n_e, MC generates events (2D image or 1D row), runs the classifier, and fills **pattern_table_[n_e][label]** = P(label | n_e).
   - The app then copies the table into **pattern_eff_map**: for each (label, n_e), `pattern_id = encode(label)`, `pattern_eff_map[{pattern_id, n_e}] = P`.

So:

- **Signal efficiency** = **ε(pattern | n_e)** from either:
  - Precomputed CSV (e.g. from a previous run or the notebook), or
  - PatternMC (2D or 1D) in this run.

### 2.3 Pipeline epsilon (when not using pattern bins)

Even in pattern mode, the app builds a **DetectorResponsePipeline** and a **PatternEfficiency** object backed by a histogram **h_eps_ne** (ε_total(n_e) for n_e):

- If **pattern_eff_map** was already filled (e.g. from CSV), **PrecomputeEpsilonWithPatternEff** is used to build h_eps_ne from that map (e.g. sum over pattern_roi per n_e).
- Otherwise **PrecomputeEpsilon** uses PatternMC’s table to fill h_eps_ne.

Then:

- **pipe.SetPatternEfficiency(pe)** and **pipe.SetSkipPatternEfficiency(use_pattern_bins)**.
- When **use_pattern_bins** is true, the pipeline **skips** applying this ε in its Apply(); the pattern-space fold is done explicitly later with **pattern_eff_map**.

So in the example app, when we are in pattern mode, the pipeline does **not** apply ε(n_e) inside Apply; all pattern weighting is done via **pattern_eff_map** and **FoldNeToPatternRates**.

### 2.4 Background spectrum in n_e: B_tot(n_e)

- **B_dc_ne** = **BackgroundBuilder::BuildBkgAsimov()** → dark-current spectrum in n_e (and optionally other components).
- **B_flat_ne** = pipeline applied to a flat dRdE (if flat background is configured).
- **B_tot** = B_dc_ne + B_flat_ne (and any other configured background in n_e).

So **B_tot** is the **total background rate per n_e** (including n_e = 0 for DC). No pattern efficiency is applied inside the background builder when `use_pattern_bins` is true; pattern folding is done in the next step.

### 2.5 Folding B_tot into pattern space → B_pat

When **use_pattern_bins** is true, **B_pat** is computed in one of two ways:

- **If `background_efficiency_csv` is set:**  
  The **DC background (and full background)** are folded through the **background efficiency** (migration matrix) from that CSV: P(identified pattern | true pattern). B_tot(n_e) is mapped to true single-pixel pattern rates (n_e = 0..5 → pattern codes 0..5), then **B_pat = FoldBackgroundWithMigration(B_true_pat, pattern_roi, migration_map)**. This applies the background-specific identification/migration, not the signal ε(pattern|n_e).

- **If `background_efficiency_csv` is empty:**  
  Background is folded with the same **signal** efficiency ε(pattern | n_e):

```cpp
B_pat = FoldNeToPatternRates(*B_tot, ne_min_bkg, ne_max, summary.pattern_roi, pattern_eff_map);
```

**FoldNeToPatternRates** (in `PatternRates.cc`) does, for each pattern_id in pattern_roi:

```
rate(pattern_id) = Σ_{n_e = ne_min}^{ne_max}  B_tot(n_e) × ε(pattern_id | n_e)
```

So:

- With no background CSV: **same ε(pattern | n_e)** is used for background as for signal.
- With background CSV: **DC and full background** are folded through the **background** migration matrix.
- In both cases, background in pattern space = **B_pat** = one rate per pattern bin in **pattern_roi**.

### 2.6 Signal spectrum: S_true(n_e) and S_obs(n_e)

- **dRdE_sig** = DM spectrum dR/dE for the chosen (m_chi, sigma_e).
- **S_true** = **ion->FoldToNe(dRdE_sig, exposure, ne_min, ne_max)** → signal in **true** n_e (ionization only). Used for the pattern fold so ε(pattern | n_true) is applied once.
- **S_obs** = **pipe.Apply(dRdE_sig, ...)** → signal in **reconstructed** n_e (for n_e-space analysis or debug).

### 2.7 Folding into pattern space → S_pat

```cpp
S_pat = FoldNeToPatternRates(*S_true, ne_min, ne_max, summary.pattern_roi, pattern_eff_map);
```

Same formula as for B_pat, but with **S_true** and the same **pattern_eff_map** ε(pattern | n_true):

```
S_pat[bin] = Σ_{n_true}  S_true(n_true) × ε(pattern_roi[bin] | n_true).
```

### 2.8 Test statistic and output

- **model_test** = S_pat + B_pat (per pattern bin).
- **q_ts** = PLR test statistic: **EvaluateRatio(B_pat, model_test, B_pat)** (e.g. profile likelihood ratio with background as null).
- The app writes **efficiency_per_pattern.csv** (pattern_roi × roi_bins) from **pattern_eff_map**, and optionally ROOT histograms and a plot (e.g. **example_one_point_pattern_rates_per_pattern.png**) of S_pat vs B_pat per pattern bin.

---

## 3. Where the two efficiency products are used

| Product | Meaning | Where it is used in the example app |
|--------|---------|-------------------------------------|
| **efficiency_per_pattern.csv** (or in-memory **pattern_eff_map**) | **Signal efficiency** ε(pattern \| n_e): P(classified as pattern \| n_e electrons). | Used for **both** B_pat and S_pat via **FoldNeToPatternRates**. So it is the **only** pattern efficiency in the pipeline. |
| **Background_efficiencies.csv** (from `ccdarksens_pattern_background_eff`) | **Migration matrix** P(identified pattern \| true pattern) for ideal clusters (e.g. DC/noise). | **Not** used in the example app. It is for analyses that model background as a mixture of true patterns and then apply this matrix. In the current pipeline, background is B_tot(n_e) folded with the same ε(pattern \| n_e). |

So:

- **Signal and DC background** in the example app both use the **same** ε(pattern | n_e) (from CSV or PatternMC).
- **Background_efficiencies.csv** is an optional, more detailed description of pattern migration (e.g. for a Poisson-mixture background model) and can be used by separate analysis code, not by this pipeline.

---

## 4. Summary diagram (example app, pattern mode)

```
                    ┌─────────────────────────────────────┐
                    │  pattern_eff_map                     │
                    │  (pattern_id, n_e) → ε(pattern|n_e)  │
                    │  from: CSV or PatternMC table        │
                    └─────────────────┬───────────────────┘
                                       │
         ┌─────────────────────────────┼─────────────────────────────┐
         │                             │                             │
         ▼                             ▼                             ▼
  B_tot(n_e)                    S_true(n_e)                efficiency_per_
  (DC + flat bkg)               (pipe.Apply(signal))       pattern.csv
         │                             │                    (written out)
         │ FoldNeToPatternRates        │ FoldNeToPatternRates
         ▼                             ▼
  B_pat[]                        S_pat[]
  (one per pattern_roi)          (one per pattern_roi)
         │                             │
         └─────────────┬───────────────┘
                       ▼
              model_test = S_pat + B_pat
                       ▼
              q = PLR(B_pat, model_test, B_pat)
```

---

## 5. Config snippets that matter

- **experiment.observable_bins** = `"pattern"` → use pattern bins and the flow above.
- **experiment.pattern_roi** → list of pattern IDs (defines bins).
- **experiment.roi_bins** → n_e values for CSV output and ROI (e.g. [1,2,3,4,5]).
- **response.pattern_mc.efficiency_csv** (optional) → path to (pattern, ne, Efficiency) CSV; if set, **pattern_eff_map** is loaded from it and PatternMC table is not used for the fold (BuildPatternTable may still run for PrecomputeEpsilonWithPatternEff / h_eps_ne if needed).
- **response.pattern_mc.use_2d_image_efficiency** → if true and no CSV, PatternMC uses 2D image + isolation to build the table (notebook-style).

This is how signal and DC background efficiencies enter the pipeline and how the example app uses them end-to-end.

---

## 6. Should we use Background_efficiencies instead?

**Short answer:** Use **Background_efficiencies** for the **background** fold when you model dark current (or similar) as a **rate per true (ideal) pattern** and want the classifier’s migration (e.g. 1→0, 2→3) to be the one measured with ideal clusters. Keep **ε(pattern | n_e)** for **signal**. The current “one ε for both” is simpler and consistent when background is already a spectrum in n_e and you assume the same classifier response.

### Two ways to fold background into pattern space

| Approach | Background model | Fold | When it’s appropriate |
|----------|------------------|------|-------------------------|
| **Current** | B_tot(n_e): rate per n_e (from DC, etc.) | B_pat = Σ_n_e B_tot(n_e) × **ε(pattern \| n_e)** | Same classifier response for signal and background; background naturally given as spectrum in n_e. |
| **Migration matrix** | B_true_pat[true_pattern]: rate per ideal pattern | B_pat[identified] = Σ_{true_pat} B_true_pat[true_pat] × **P(identified \| true_pat)** from **Background_efficiencies.csv** | DC modeled as ideal clusters that then get misclassified; you want the migration measured from ideal-cluster simulations. |

### When to use Background_efficiencies for background

- **Use the migration matrix** when:
  - You think of DC (or surface) background as **ideal clusters** (true patterns 0, 1, 2, …, 311) with rates B_true_pat, and
  - You want **P(identified \| true pattern)** from the background-efficiency simulation (noise + isolation, same as notebook) rather than ε(pattern | n_e) from signal MC.
- **Keep the current approach** when:
  - Background is already a histogram in n_e (e.g. B_tot from a DC model in electron space), and
  - You’re okay assuming the same ε(pattern | n_e) for background as for signal (one table, simpler).

### What you need to use Background_efficiencies in the pipeline

1. **Background in true-pattern space**  
   You need **B_true_pat[true_pattern]** (one rate per ideal pattern). For DC often only single-pixel patterns (0)–(5) are used:
   - B_true_pat[(0)] = B_tot(0), B_true_pat[(1)] = B_tot(1), …, B_true_pat[(5)] = B_tot(5), and B_true_pat[multi-pixel] = 0 (or from a separate model).
2. **Load the migration matrix**  
   Read **Background_efficiencies.csv**: row = true pattern, columns = P(identified | true). Build a (true_pat, iden_pat) → probability map.
3. **Fold**  
   For each observed pattern bin in pattern_roi:  
   `B_pat[iden_pat] = Σ_{true_pat} B_true_pat[true_pat] × P(iden_pat | true_pat)`.
4. **Signal unchanged**  
   Keep S_pat = FoldNeToPatternRates(S_true, …, **pattern_eff_map**); only the **background** fold switches to the migration matrix.

So: **yes, you should use Background_efficiencies for the background fold** if your physics and background model are “rate per true pattern” + migration. The example app doesn’t do this yet; the app uses **background_efficiency_csv** in the `backgrounds` config: when set, it loads the migration matrix, builds B_true_pat from B_tot (single-pixel n_e 0..5 → patterns (0)–(5); multi-pixel true rates 0), and folds with FoldBackgroundWithMigration. Omit or leave empty to keep B_pat = FoldNeToPatternRates(B_tot, …, pattern_eff_map). To add more background types with misidentification, extend B_true_pat (e.g. from other components) and use the same or component-specific migration CSV.
