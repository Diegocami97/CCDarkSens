# Profile-likelihood minimization flow (detail)

This document describes step-by-step how the minimization over the nuisance parameter θ is done in the pattern scan, including inputs, methods, and data flow.

---

## 1. Config inputs (run section)

From the JSON config, the scan reads:

| Key | Meaning | Used as |
|-----|---------|--------|
| `use_profile_likelihood` | Enable profile likelihood | If true, create `ProfileLikelihood` and minimize over θ at each (mχ, σ). |
| `profile_minimizer` | `"brent"`, `"minuit"`, `"minuit2d"`, or `"pydme"` | **brent** / **minuit**: 1D over θ at each grid σ; UL = first grid σ with q_μ ≥ target_q. **minuit2d**: one 2D Minuit over (log₁₀(σ_e), θ) per mχ for nll_min; q_μ and UL still on grid. **pydme**: same 2D fit as minuit2d; **UL from bracketing + bisection in log₁₀(σ)** (match pydme `compute_upper_limit_fast`); q_μ curve still filled on grid for plotting. |
| `theta_lo`, `theta_hi` | Bounds for θ | Passed as `profile_param_lo`, `profile_param_hi` into every minimization call. |
| `background_model` | `"Bp_theta_Br"` or `"scale"` | If `Bp_theta_Br`, background is B_i = Bp_i + θ·Br_i; θ is the parameter we minimize over. |
| `background_Bp`, `background_Br` | Vectors (per pattern) | Set on `ProfileLikelihood` via `SetBpBr(Bp, Br)`. |
| `constrain_prior_strength` | e.g. 98 | Weight of the pydme-style prior term in NLL (0 = no prior). |
| `data_path` | Path to data or empty | If empty, Asimov: data = B_pat; otherwise data = observed counts per pattern. |
| `single_bin_likelihood` | Collapse to one bin | If true, S and data are summed to one number; otherwise one bin per pattern. |

---

## 2. One-time setup (before the grid loop)

**File:** `apps/ccdarksens_scan_dmelectron_pattern.cc`

1. **Create `ProfileLikelihood`** (around line 942)  
   `profile_pl = std::make_unique<ccdarksens::stats::ProfileLikelihood>();`

2. **Set data**  
   - Asimov (no `data_path`): `profile_pl->SetData(B_pat)` — data = expected background per pattern.  
   - With data file: `profile_pl->SetData(data)` — data = observed counts per pattern (same length as `B_pat`).

3. **Set background model**  
   - **Bp_theta_Br:**  
     `profile_pl->SetBpBr(run.background_Bp, run.background_Br)`  
     Then: `profile_pl->SetConstrainPriorStrength(run.constrain_prior_strength)`.  
   - **Scale mode:**  
     `profile_pl->SetBTemplate(B_pat)` (and optional Gaussian constrain if configured).

4. **Null hypothesis (background only)**  
   The null case is **background only**: no signal. We set  
   `profile_S_null` = vector of zeros (same length as pattern bins, or `{0.0}` in single-bin).  
   So the **model** for the null is μ = S + B(θ) = **0 + B(θ) = B(θ)**.  
   When we compute `nll_null`, we minimize NLL(data | μ = B(θ)) over θ — i.e. we fit the background model B(θ) to the data with no signal. So the null is exactly “background only”.

5. **Minimization bounds**  
   `profile_param_lo = run.theta_lo`, `profile_param_hi = run.theta_hi` (when `UseBpBr()`).

---

## 3. Per–(mχ, σ) flow: what gets minimized

For each mass `mchi` and cross-section `sigma_val` in the grid:

1. **Build signal**  
   Detector pipeline (rates → pattern fold, etc.) produces **S_pat**: expected signal counts per pattern bin for this (mχ, σ).

2. **Input to minimization**  
   - **S_for_pl** = `S_pat` (or, in single-bin mode, a 1-element vector with the sum of S_pat).  
   - **Bounds** = `(profile_param_lo, profile_param_hi)` = `(theta_lo, theta_hi)` from config.

3. **Minimization call**  
   - **Brent:**  
     `nll = profile_pl->MinimizeOverScale(S_for_pl, profile_param_lo, profile_param_hi).second`  
   - **Minuit:**  
     `nll_m = profile_pl->MinimizeOverScaleMinuit(S_for_pl, profile_param_lo, profile_param_hi).second`  
     `nll_b = profile_pl->MinimizeOverScale(S_for_pl, profile_param_lo, profile_param_hi).second`  
     `nll = std::min(nll_m, nll_b)`  
   So we always store one NLL value per (mχ, σ).

4. **After the σ-loop for this mχ**  
   - **Grid / minuit:** `nll_min` = minimum of all `nll_values` over σ.  
   - **minuit2d:** Before the σ-loop, precompute S at all grid σ; build S(log₁₀(σ)) interpolator; run one 2D Minuit over (log₁₀(σ_e), θ) → nll_min; at each grid σ store nll(σ) from 1D min over θ. UL = first grid σ with q_μ ≥ target_q.  
   - **pydme:** Same 2D fit and grid q_μ curve as minuit2d. **Upper limit** = **bracketing + bisection in log₁₀(σ)** (as in pydme `compute_upper_limit_fast`): start at σ̂ from the 2D fit, step right until q_μ(σ) ≥ target_q, then bisect to find the crossing; UL = 10^(that log₁₀(σ)). If bracketing never reaches target_q, fall back to grid UL.  
   - `q_μ(σ)` = 2·(nll(σ) − nll_min).  
   - Upper limit: **grid / minuit2d** = smallest grid σ with q_μ ≥ target_q; **pydme** = bisection crossing in log₁₀(σ).

So the **minimization** is: for a **fixed** S (and fixed data, Bp, Br, prior strength), find the θ in [theta_lo, theta_hi] that minimizes NLL(S, θ), and return that minimum NLL (and, inside the minimizer, the best-fit θ). With **minuit2d**, the **global** minimum is found in one 2D run (so nll_min can be slightly lower than the minimum over the grid), matching pydme’s single 2D fit per mass.

---

## 4. NLL (objective function)

**File:** `src/stats/ProfileLikelihood.cc`  
**Method:** `double ProfileLikelihood::NLL(const std::vector<double>& S, double param) const`

- **Inputs:**  
  - `S` = signal expectation per bin (length = number of pattern bins, or 1 in single-bin).  
  - `param` = θ (nuisance parameter).

- **Internal state (set at setup):**  
  - `data_` = observed or Asimov counts per bin.  
  - `Bp_`, `Br_` (Bp_theta_Br mode) or `B_template_` (scale mode).  
  - `constrain_prior_strength_`, optional `constrain_` (e.g. Gaussian on scale).

- **Bp_theta_Br mode (your use case):**  
  For each bin i:  
  - B_i = Bp_i + param * Br_i  
  - μ_i = S[i] + B_i  
  - NLL += μ_i − data_[i]*ln(μ_i)  
  If `constrain_prior_strength_ > 0`:  
  - NLL += Σ_i ( −θ·Br_i + constrain_prior_strength_ * ln(θ·Br_i) )  
  If `constrain_` is set: NLL += constrain_(param).

- **Output:** One number = NLL(S, θ). The minimizer varies θ and calls this repeatedly.

---

## 5. Brent minimization

**File:** `src/stats/ProfileLikelihood.cc`  
**Method:** `std::pair<double, double> ProfileLikelihood::MinimizeOverScale(const std::vector<double>& S, double scale_lo, double scale_hi, double tol) const`

- **Inputs:**  
  - `S` = S_for_pl (signal per bin for this (mχ, σ)).  
  - `scale_lo`, `scale_hi` = theta_lo, theta_hi.  
  - `tol` = 1e-6 (default).

- **Algorithm:**  
  - 1D Brent in [scale_lo, scale_hi].  
  - At each step calls `NLL(S, x)` for a trial x; no fixed grid of θ values.  
  - Stops when bracket width is small relative to `tol` or after 200 iterations.

- **Output:**  
  - `(theta_hat, nll_min)` = best θ and minimum NLL.  
  - The scan uses only `.second` (nll_min).

---

## 6. Minuit2 (Simplex) minimization

**File:** `src/stats/ProfileLikelihood.cc`  
**Method:** `std::pair<double, double> ProfileLikelihood::MinimizeOverScaleMinuit(const std::vector<double>& S, double scale_lo, double scale_hi) const`

- **Inputs:** Same S, scale_lo, scale_hi (no tol; Minuit has its own convergence).

- **Steps:**  
  1. Create minimizer: `ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex")`.  
  2. Wrap NLL: `NLLFunctor` implements `DoEval(const double* x)` → `pl->NLL(S, x[0])`.  
  3. Set one variable: name `"theta"`, initial value `0.5*(scale_lo+scale_hi)`, step `0.01*(scale_hi-scale_lo)`, limits `[scale_lo, scale_hi]`.  
  4. `SetTolerance(1e-8)`, `SetPrintLevel(0)`.  
  5. Call `min->Minimize()`.  
  6. If it fails: return `MinimizeOverScale(S, scale_lo, scale_hi)` (Brent).  
  7. If best-fit θ is within 1% of lo or hi: same Brent fallback.  
  8. Otherwise return `(min->X()[0], min->MinValue())`.

- **Output:** Same as Brent: (theta_hat, nll_min); scan uses only nll_min.

---

## 7. Summary diagram

```
Config (theta_lo, theta_hi, profile_minimizer, Bp, Br, data, …)
    ↓
Setup: ProfileLikelihood::SetData, SetBpBr, SetConstrainPriorStrength
    ↓
For each mχ:
    nll_null = minimize over θ in [theta_lo, theta_hi] with S = 0  (one call; null = background only, μ = B(θ))
    For each σ:
        S_pat = signal for (mχ, σ)
        S_for_pl = S_pat (or sum in single-bin)
        nll = minimize over θ in [theta_lo, theta_hi] with S = S_for_pl
              → Brent: MinimizeOverScale(S_for_pl, lo, hi)
              → Minuit: min(MinimizeOverScaleMinuit(...), MinimizeOverScale(...))
        nll_values.push_back(nll)
    nll_min = min(nll_values)
    q_μ(σ) = 2*(nll(σ) - nll_min)
    upper_limit = smallest σ with q_μ ≥ target_q
```

Each “minimize over θ” call uses only the **bounds** (theta_lo, theta_hi); there is no fixed number of θ steps—Brent and Simplex both search the interval until convergence.
