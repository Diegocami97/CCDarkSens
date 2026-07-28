# Exposure usage audit

**Goal:** Confirm `exposure_kg_year` is applied **once** when converting rates to counts, and is not double-applied anywhere.

**Definition:** `exposure_kg_year = livetime_days * duty_cycle * mass_kg / 365.25` (from config; set in `ExperimentSetup::prepare_summary()`).

---

## 1. Where exposure is applied (rate → counts)

Exposure is used **only** when converting a **rate** (events per kg·year per eV) into **counts**:

| Location | What | Formula | Applied |
|----------|------|---------|---------|
| `ChargeIonization::FoldToNe(dRdE, exposure_kg_year, ...)` | dR/dE → S(n_e) counts | `counts = rate * exposure_kg_year * dE` per bin | **Once** per call |
| `DetectorResponsePipeline::Apply(dRdE, exposure_kg_year, ...)` | Pattern/PCD path | Calls `ion_->FoldToNe(dRdE, exposure_kg_year, ...)` **once**, then kernel/pattern folds (no exposure) | **Once** |
| `DetectorResponsePipeline::ApplyEDependent(...)` | E-dependent path | `weight = dRdE_val * exposure_kg_year * dE` per bin | **Once** per bin |

So every use is of the form: **rate × exposure_kg_year → counts**. There is no second multiplication by exposure downstream.

---

## 2. Scan app (`ccdarksens_scan_dmelectron_pattern.cc`) — usage list

| Line (approx) | Code | Role |
|---------------|------|------|
| 723 | `B_flat_ne = ion->FoldToNe(*dRdE_flat, summary.exposure_kg_year, ...)` | Flat background: rate → counts. **Once.** |
| 726 | `pipe.Apply(*dRdE_flat, ..., summary.exposure_kg_year, ...)` | Flat bkg in non-pattern branch. Apply() uses exposure once inside. **Once.** |
| 1119 | `S_true = ion->FoldToNe(*dRdE_sig, summary.exposure_kg_year, ...)` | S_grid for pydme/minuit2d. **Once** per σ. |
| 1212 | `S_true = ion->FoldToNe(*dRdE_sig, summary.exposure_kg_year, ...)` | Signal at each (mχ, σ) in main loop. **Once.** |
| 1216 | `pipe.Apply(*dRdE_sig, ..., summary.exposure_kg_year, ...)` | Signal in non-pattern branch. **Once** inside Apply. |

**Downstream of the above (no exposure):**

- `FoldNeToPatternRates(S_true, ...)` → takes **counts** in n_e, returns pattern-space **counts** (sum over n_e with ε). No exposure parameter; does not multiply by exposure.
- `B_pat` when `background_source == "bp_br_template"` → from config (Bp+Br). No exposure.
- `B_dc_ne` → from `BackgroundBuilder::BuildBkgAsimov()`. Uses **livetime_days**, **n_exposures**, **n_pixels** (and λ per pix per year). Does **not** use `exposure_kg_year`. So DC and signal/flat use different exposure concepts (pixel-exposures vs kg·year), by design; no double-count.

---

## 3. Other components (no double application)

- **ProfileLikelihood:** Uses vectors of **counts** (S_pat, B_pat, data). No exposure parameter; NLL is in counts.
- **Rate tables (dR/dE):** Stored as rate per kg·year per eV. Exposure is applied only in `FoldToNe` / `Apply` as above.
- **BackgroundBuilder (DC):** Uses livetime and number of exposures (and pixels), not `exposure_kg_year`. So no overlap with the kg·year exposure used for signal/flat.

---

## 4. Conclusion

- **exposure_kg_year** is applied **exactly once** per spectrum when going from dR/dE (rate) to counts: in `ChargeIonization::FoldToNe` or inside `DetectorResponsePipeline::Apply` / `ApplyEDependent`.
- **FoldNeToPatternRates** and the likelihood operate on **counts**; they do not take or multiply by exposure.
- **B_pat** in the template path is fixed from config; in the DC+flat path it is built from B_dc (livetime-based) and B_flat (one use of exposure_kg_year in FoldToNe). No second application of exposure anywhere.

**No double application of exposure was found.**

---

## 5. Signal rate and counts units (for limit scan)

| Quantity | Units | Where |
|----------|--------|--------|
| **Rate table CSV** (dR/dE) | events/(kg·year·eV) | RateTable loads column as R_kg_year_eV; CSV header must match (see e.g. `# Output units: dR/dE in events / kg / year / eV`). |
| **exposure_kg_year** | kg·year | `ExperimentSetup::prepare_summary()`: livetime_days × duty_cycle × mass_kg / 365.25. |
| **FoldToNe** | counts = rate × exposure_kg_year × dE | ChargeIonization.cc. Input dRdE bin content [events/(kg·year·eV)]; output TH1D in n_e has **counts**. |
| **FoldNeToPatternRates** | input: counts per n_e; output: **counts** per pattern | PatternRates.cc. Sum over n_e of (S_ne × ε(pattern, n_e)); no exposure; result is counts. |
| **S_pat, B_pat, data** in ProfileLikelihood | **counts** (per pattern bin) | Same units for signal, background, and observed data. |

**Check:** Rate tables must be normalized in events/(kg·year·eV). If the reference uses different units (e.g. events/(kg·day·eV)), convert before loading or scale exposure accordingly (e.g. use exposure in kg·day and ensure rate is per kg·day·eV so product is counts).

---

## 6. Rate generation (where dR/dE CSVs are produced)

The pipeline expects CSV columns **E (eV), dRdE** with **dRdE in events/(kg·year·eV)**. RateTable does not convert units; it loads the second column as R_kg_year_eV.

| Generator | Units produced | Used by |
|-----------|----------------|--------|
| **ccdarkphys/qedark/entry.py** `compute_dRdE()` | events/(kg·year·eV) | Internal: prefactor (kg·year)^−1, then ÷ dE. Returns `dRdE_kg_year_eV`. |
| **ccdarkphys/common/io.py** `write_csv()` | Writes R as given; header states "events / kg / year / eV". | Callers must pass `dRdE_kg_year_eV`. |
| **utils/qedark_generate_grid.py** | Correct. Uses `compute_dRdE` → `dRdE_kg_year_eV`, writes via common/io. | Use this (or entry CLI) for pipeline rate tables. |
| **utils/qedark_generate.py** (bridge) | Previously wrote **events/(g·day·eV)** from the bridge. Now converts to kg·year·eV before writing so output matches pipeline. | Single-file generation; output is now pipeline-ready. |

**Summary:** Use **qedark_generate_grid.py** (or `python -m ccdarkphys.qedark.entry` with common/io) to produce rate tables. If you use **qedark_generate.py**, it now converts bridge output to events/(kg·year·eV) so the written CSV is correct for RateTable.
