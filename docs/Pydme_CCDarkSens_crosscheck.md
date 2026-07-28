# pydme vs CCDarkSens cross-check

Reference: `collab_frameworks/pydme/`  
Config: `configs/scan_dmelectron_pattern_pydme_minuit.json`

---

## 0. pydme framework findings (for agreement beyond 1 MeV)

| Item | pydme value | Location |
|------|-------------|----------|
| **Exposure** | `exposure = np.sum(Npix * texp * Mpix)` in **[g·day]** | `lbc_dmanalysis_upperlimits.py:64` |
| **Exposure (commented override)** | `1250` g·day | `lbc_dmanalysis_upperlimits.py:67` |
| **Exposure default in NLLAnalysis** | `250` (unused when set from data) | `dme_nLL_exclusion.py:326` |
| **Conversion** | 1250 g·day = 1.25 kg·day ≈ **0.00342 kg·year** | 1250/1000/365.25 |
| **CCDarkSens exposure** | 1.3 kg·day = 85.356 × 0.01523 / 365.25 ≈ **0.00356 kg·year** | Config |
| **Efficiency for signal generation** | `Efficiencies_patterns_Nsims100000_DCTrue_alpha1.csv` (100k) | `srdm_rates_with_pydme.py:45` |
| **theta_limits** | `[[0.5, 10]]` (overrides default 1e-6–0.1) | `lbc_dmanalysis_upperlimits.py:84` |
| **Halo (SRDM analysis)** | `SRDMModulation` | `lbc_dmanalysis_upperlimits.py:43` |
| **Halo (SHM)** | Standard halo, signal ∝ σ_e | Used when loading SHM signal files |
| **Signal units** | events/g/day (in pattern signal files) | `data/__init__.py:68` |
| **Mpix (binning=100)** | V=0.0015×0.0015×0.0669×100 cm³ × 2.33 g/cm³ ≈ 3.5e-5 g | `get_mpix_grams()` |

**Note:** The DAMIC-M 2025 reference (solid line) may come from a different analysis (e.g. SHM vs SRDM, different exposure). The pydme SRDM analysis uses SRDMModulation; CCDarkSens uses SHM (standard halo) with QEDark rates.

---

## 1. Background model (Bp + θ·Br)

| Item | pydme (`background_pattern.py`) | CCDarkSens (`scan_dmelectron_pattern_pydme_minuit.json`) |
|------|---------------------------------|----------------------------------------------------------|
| **Bp** | `[141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]` | `[141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]` ✓ |
| **Br** | `[0.039, 0.039, 0.016, 0.052, 0.011, 0.035]` | `[0.039, 0.039, 0.016, 0.052, 0.011, 0.035]` ✓ |
| **Pattern order** | `[11, 21, 111, 31, 22, 211]` | `[11, 21, 111, 31, 22, 211]` ✓ |
| **Constrain** | `-θ·Br + 98·ln(θ·Br)` per bin | Same (constrain_prior_strength=98) ✓ |
| **B normalization** | pydme divides Bp, Br by `len(gamma)` (number of data rows). For 6 patterns × 1 row each, total B = Bp + θ·Br. | CCDarkSens uses B_pat = Bp + Br directly (theta=1 nominal). Same total. ✓ |

---

## 2. Theta bounds

| Item | pydme (`lbc_dmanalysis_upperlimits.py`) | CCDarkSens |
|------|----------------------------------------|------------|
| **theta_limits** | `[[0.5, 10]]` | `theta_lo: 0.5`, `theta_hi: 10` ✓ |

---

## 3. Efficiency file — mismatch

| Item | pydme (reference) | CCDarkSens `scan_dmelectron_pattern_pydme_minuit.json` |
|------|-------------------|--------------------------------------------------------|
| **efficiency_csv** | Uses `Efficiencies_patterns_Nsims100000_DCTrue_alpha1.csv` (srdm_rates_with_pydme) or similar | **`data/efficiencies_paolo.csv`** |
| **efficiency_csv_reference** | — | **`data/efficiencies_paolo.csv`** |

**Difference:** Paolo uses Nsims=100k, lambda=0.00015; 1M reference uses Nsims=1M, lambda=0.00041. Several bins differ (e.g. 22@4: Paolo 0.133 vs 1M 0.114; 31@5: Paolo 0.49 vs 1M 0.065).

**To match pydme / DAMIC-M 2025 reference:** Use
```json
"efficiency_csv": "data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv",
"efficiency_csv_reference": "data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv"
```

---

## 4. Exposure

| Item | pydme | CCDarkSens |
|------|-------|------------|
| **Formula** | `exposure = sum(Npix * texp * Mpix)` (gram-days) | `exposure_kg_year = livetime_days * duty_cycle * mass_kg / 365.25` |
| **Config** | From data (texp, Npix per row) | `livetime_days: 85.467`, `mass_kg: 0.01523`, `duty_cycle: 1.0` |
| **Result** | Depends on dataset | ~0.00356 kg·year |

**Check:** Confirm the reference exposure (kg·year) and align `livetime_days` and `mass_kg` if needed.

---

## 5. Signal rates

| Item | pydme | CCDarkSens |
|------|-------|------------|
| **Source** | `signal_rates` CSV, columns S11;S21;S111;S31;S22;S211 | `model.rates_dir` + `filename_template`, QEDark dR/dE tables |
| **Path** | `analysis/SRDM/upper_limit/signal_files/` | `data/qedark_rates/Si/heavy/after_eta_fix/` |
| **Format** | Pre-folded pattern rates | dR/dE(E) → FoldToNe → FoldNeToPatternRates |

**Check:** Same mediator (heavy), same material (Si), same rate tables for a fair comparison.

---

## 6. Summary of actions to match reference

1. **Efficiency:** Switch to `Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv` in the config.
2. **Exposure:** Verify reference exposure (kg·year) and adjust `livetime_days` or `mass_kg` if different.
3. **Rates:** Ensure `rates_dir` and rate tables match the reference (e.g. heavy mediator, after_eta_fix).
4. **theta_hi:** If using pattern 22@4 with higher efficiency (e.g. 0.112), set `theta_hi: 100` to avoid boundary issues.

---

## 7. Agreement beyond 1 MeV — checklist

If your dotted curve (CCDarkSens) disagrees with the solid reference (DAMIC-M 2025) for m_χ > 1 MeV:

1. **Exposure:** Match kg·year. CCDarkSens: `livetime_days * duty_cycle * mass_kg / 365.25`. pydme uses g·day from data; 1250 g·day ≈ 0.00342 kg·year.
2. **Efficiency:** Use the same file as the reference. pydme’s rate generator uses `Efficiencies_patterns_Nsims100000` (100k); the 1M file may differ in several bins.
3. **Halo model:** Confirm reference uses SHM (standard halo). pydme SRDM uses SRDMModulation (different signal shape).
4. **Rate tables:** Same mediator (heavy), material (Si), and dR/dE source. pydme signal files are pre-folded with efficiency; CCDarkSens folds at runtime.
5. **Bp, Br, theta bounds:** Already aligned (see §1, §2).
