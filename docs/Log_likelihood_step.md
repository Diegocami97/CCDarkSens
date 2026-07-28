# Log-likelihood step (current implementation)

This document describes how the binned log-likelihood and test statistic are currently computed in the pattern pipeline and scan.

---

## 1. Interface and implementation

- **Interface:** `ccdarksens::stats::ITestStatistic` (`include/ccdarksens/stats/ITestStatistic.hh`)
  - **EvaluateNLL(data, model):** \(-2\ln L\) contribution from a Poisson binned likelihood (up to constants).
  - **EvaluateRatio(data, model_test, model_null):**  
    \(q = -2\ln\bigl[ L(\text{data}\,|\,\text{model\_test}) / L(\text{data}\,|\,\text{model\_null}) \bigr]\).

- **Implementation:** `PoissonAsimovPLR` (`src/stats/PoissonAsimovPLR.cc`)

  - **NLL (per bin):**  
    For bin \(i\) with observed count \(n_i\) and expected rate \(\mu_i\):  
    \[
    -\ln L_i \;=\; \mu_i - n_i \ln(\mu_i)
    \]  
    (plus constant \(\ln(n_i!)\) omitted).  
    - \(\mu_i \leq 0\) with \(n_i > 0\): large penalty (1e9).  
    - \(\mu_i = n_i = 0\): no contribution.

  - **Likelihood ratio:**  
    \[
    q \;=\; 2\,\bigl[ \text{NLL}(\text{data}\,|\,\text{model\_test}) - \text{NLL}(\text{data}\,|\,\text{model\_null}) \bigr].
    \]  
    Returned value is \(\max(q, 0)\) for numerical safety.

  - **Asimov (background-only) special case:**  
    When \(\text{data} = \text{model\_null} = B\) and \(\text{model\_test} = B + S\):  
    \[
    q \;=\; 2\sum_i \Bigl[ S_i - B_i \ln\Bigl(1 + \frac{S_i}{B_i}\Bigr) \Bigr].
    \]  
    (With the convention that \(B_i=0,\,S_i>0 \Rightarrow\) term \(= 2S_i\).)

- **Factory:** `MakeTestStatistic(StatisticsConfig)` in `TestStatisticFactory.cc` builds a test statistic from config; currently only **"PLR"** is used, which returns `PoissonAsimovPLR`.

---

## 2. How the example app uses it (pattern space)

In `ccdarksens_example_one_point_pattern.cc` (pattern bins):

1. **Bins:** One bin per element of `pattern_roi` (pattern IDs).

2. **Vectors passed to PLR:**
   - **data** = `B_pat` (expected counts under background-only).
   - **model_null** = `B_pat` (background-only hypothesis).
   - **model_test** = `S_pat + B_pat` (signal + background).

3. **Call:**  
   `q_ts = ts->EvaluateRatio(B_pat, model_test, B_pat);`

So this is the **Asimov** setup: “data” is set to the background expectation; the test statistic measures how much the likelihood degrades when moving from the best fit under signal+background to the background-only fit (which is perfect for this Asimov data). So \(q\) is the usual **profile likelihood ratio for discovery** (null = background, test = signal+background), evaluated at the Asimov background-only “data”.

---

## 3. How the scan app uses it

In `ccdarksens_scan_dmelectron_pattern.cc`:

- **Pattern bins** (`use_pattern_bins == true`):
  - **data** = `B_pat`
  - **model_null** = `B_pat`
  - **model_test** = `S_pat + B_pat` (with `S_pat` from `FoldNeToPatternRates(S_obs, ..., pattern_eff_map)` and the same `B_pat` as in the example app).
  - Same Asimov convention as the example app.

- **n_e ROI bins** (`use_pattern_bins == false`):
  - Bins = `roi_bins` (e.g. n_e = 1, 2, 3, 4, 5).
  - **data** = `B_tot(ne)` for each ROI n_e.
  - **model_null** = same (background).
  - **model_test** = `S_obs(ne) + B_tot(ne)`.
  - Again Asimov: “data” = background expectation.

In both cases the actual call is:

```cpp
double q_ts = ts->EvaluateRatio(data, model_test, model_null);
```

So the **log-likelihood step** is: build the three vectors (data, model_test, model_null) from S and B in the chosen space (pattern or n_e), then compute \(q\) via `PoissonAsimovPLR::EvaluateRatio`. No profiling over nuisance parameters or multiple likelihood terms—just one Poisson binned likelihood in the chosen observable space.

---

## 4. Summary

| Item | Current implementation |
|------|------------------------|
| **Likelihood** | Poisson, one count per bin (pattern or n_e). |
| **NLL** | \(\sum_i \bigl[\mu_i - n_i\ln(\mu_i)\bigr]\), with safe handling of \(\mu_i\leq 0\) and \(n_i=0\). |
| **Test statistic** | \(q = 2[\text{NLL}(\text{data}|\text{test}) - \text{NLL}(\text{data}|\text{null})]\). |
| **Usage** | Asimov: data = B, null = B, test = S+B ⇒ \(q = 2\sum_i [S_i - B_i\ln(1+S_i/B_i)]\). |
| **Bins** | Pattern space: one bin per `pattern_roi`; n_e space: one bin per `roi_bins`. |
| **Nuisance parameters** | None; rates are fixed (no profiling or systematics in the likelihood). |

For the full scan, the same PLR is evaluated at each grid point. See **docs/Log_likelihood_pydme_comparison.md** for comparison with pydme and options to proceed.

For the full scan, the same PLR is evaluated at each (m_χ, σ_e) grid point using that point’s S_pat (or S_obs) and the same B_pat (or B_tot).
