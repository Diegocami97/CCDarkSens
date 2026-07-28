# Log-likelihood: comparison with pydme and how to proceed

This document compares the current CCDarkSens log-likelihood step with the **pydme** framework and outlines options to align or extend our implementation.

---

## 0. Current alignment (pattern scan with profile likelihood)

When **`use_profile_likelihood`** is true and **`background_model`** is **`Bp_theta_Br`**, the scan does the following, matching pydme:

| Step | pydme | CCDarkSens |
|------|--------|------------|
| **NLL** | Poisson: `sum_i [ (S_i+B_i) - D_i*ln(S_i+B_i) ]` + constrain | Same: `mu_i - n_i*ln(mu_i)` per bin (ProfileLikelihood::NLL). |
| **Background** | B = Bp + theta\*Br (Background_pattern) | B_i = Bp_i + theta\*Br_i from config (SetBpBr). |
| **Constrain** | L = sum_i ( -theta\*Br_i + **98**\*ln(theta\*Br_i) ) | Same: `constrain_prior_strength_` (default **98**) → -theta\*Br_i + strength\*ln(theta\*Br_i). |
| **Minimize** | Minuit over (xsec_e, theta); for UL fix xsec_e, minimize over theta. | For each (mchi, sigma): **Brent** 1D over theta only (MinimizeOverScale(S, theta_lo, theta_hi)). Sigma is on a grid (we do not minimize over sigma). |
| **nll_null** | NLL at sigma=0, theta profiled | MinimizeOverScale(profile_S_null, …) with S_null=0. |
| **nll(sigma)** | NLL at fixed sigma_e, theta profiled | MinimizeOverScale(S_for_pl, …) for each grid sigma. |
| **nll_min** | min over sigma of nll(sigma) | min over sigma_list of nll_values. |
| **q_mu** | 2\*(nll(sigma) - nll_min) | 2\*(nll_values[k] - nll_min). |
| **target_q** | 2\*delta_nll, delta_nll = 0.5\*(norm.ppf(CL))² → target_q = (norm.ppf(CL))² | target_q = NormQuantile(CL)² (e.g. 2.71 for 90% CL). |
| **Upper limit** | Smallest sigma_e where q_mu ≥ target_q (bracketing + bisection in log(sigma)) | Smallest **grid** sigma where q_mu ≥ target_q (no interpolation). |

**Theta bounds (important for Bp+Br from data):**

- With **Bp+Br** derived from data, the nominal background is **B = Bp + Br**, i.e. **θ = 1**. The minimizer must be allowed to reach θ ≈ 1, so **theta_lo and theta_hi must bracket 1**.
- **Ranges that work:** e.g. **`theta_lo: 0.5`, `theta_hi: 10`** or **`0.1` / `10`**, **`0.3` / `5`**, **`0.01` / `100`**. Any interval that contains 1 is valid.
- **Range that breaks the limit curve:** pydme’s generic default **θ ∈ [1e-6, 0.1]** does *not* contain 1. Using it forces the fit to the boundary (θ = 0.1), flattens q_μ and can produce a flat/wrong exclusion. For this Bp+Br setup, use bounds that include 1 (e.g. **0.5 and 10**).

**Other parameters:**

- **constrain_prior_strength:** 98 matches pydme’s Background_pattern prior; can be changed if your prior differs.
- **profile_minimizer:** `"brent"` (default) or `"minuit"`. Both minimize NLL over θ in [theta_lo, theta_hi]. **Brent** is built-in; **minuit** uses ROOT Minuit2 when enabled at build time (same family as pydme's iminuit). If built without Minuit2, `"minuit"` falls back to Brent. To enable Minuit2: compile with `-DCCDARKSENS_USE_MINUIT2` and link `ROOT::Minuit2`.

---

## 1. Where pydme does the likelihood

- **Main class:** `NLLAnalysis` in **`collab_frameworks/pydme/pydme/dme_nLL_exclusion.py`**
- **Usage:** e.g. **`pydme/analysis/SRDM/upper_limit/lbc_dmanalysis_upperlimits.py`** (pattern analysis, upper limits)
- **Minimizer:** iminuit (Minuit) for NLL minimization

---

## 2. Poisson NLL (same core formula)

- **pydme** (`nLL_Poisson`):  
  `-ln L = sum_i [ (S_i + B_i) - D_i * ln(S_i + B_i) ] + constrain`  
  So per bin: model = S+B, data = D; same Poisson term `mu - n*ln(mu)`. The **constrain** term depends on background parameters theta (e.g. from `Background_pattern`: prior/penalty on theta).

- **CCDarkSens** (`PoissonAsimovPLR::EvaluateNLL`):  
  Same per-bin term: `mu_i - n_i*ln(mu_i)`; no constrain.

So the **Poisson likelihood term** is the same; pydme adds a constrain term and **profiles over theta**.

---

## 3. Parameters and minimization

| Aspect | CCDarkSens (current) | pydme |
|--------|----------------------|--------|
| **Parameters** | None; S and B are fixed (precomputed). | **xsec_e** (log10) and **theta** (background params, e.g. one per pattern). |
| **Minimization** | None. | **Minuit** minimizes NLL over (xsec_e, theta). |
| **Background** | Fixed B_pat (or B_tot). | B = B(theta); theta is fitted (optionally constrained). |
| **Signal** | S_pat(sigma_e) at each grid point. | S(xsec_e) from model; xsec_e is the parameter of interest. |

---

## 4. Test statistic and upper limit

| Aspect | CCDarkSens (current) | pydme |
|--------|----------------------|--------|
| **Discovery q** | q = 2[NLL(B \| S+B) - NLL(B \| B)] with **Asimov data = B**. | **q0** = 2[NLL(sigma=0, theta_hat_hat) - NLL(sigma_hat, theta_hat)] (profile over theta). |
| **Upper-limit q_mu** | Not implemented. | **q_mu** = 2[NLL(sigma_e, theta_hat_hat) - NLL(sigma_hat, theta_hat)]; find sigma_e where q_mu = CL threshold. |
| **Upper limit** | Not implemented (only q at one sigma_e). | **Bracketing + bisection** (fast) or **scan** (slow) over sigma_e. |

pydme does **full profile likelihood**: minimize over theta, then define q_mu and q0, and compute upper limits by inverting the test statistic.

---

## 5. Binning and single-bin option

- **pydme**: Can run **binned** (one bin per pattern or per time/gamma) or **single-bin** (`do_single_bin_likelihood = True`): sum all rates and data into one Poisson count. The SRDM upper-limit script uses single-bin.
- **CCDarkSens**: Always binned (pattern_roi or roi_bins); no single-bin option.

---

## 6. Summary of differences

1. **Nuisance parameters**: pydme profiles over background parameters (theta); we use fixed B.
2. **Test statistic**: We use a single Asimov ratio; pydme uses profile-likelihood q_mu and q0.
3. **Upper limit**: We only evaluate q at given (m_chi, sigma_e); pydme finds sigma_e at which q_mu reaches the CL threshold.
4. **Constrain term**: pydme adds a term in theta to the NLL (e.g. from Background_pattern); we do not.
5. **Single vs multi-bin**: pydme can collapse to one bin; we stay multi-bin.

---

## 7. How to proceed

### Option A: Keep current method (simplest)

- **Pros:** Already implemented; fast; no minimizer; easy to scan (m_chi, sigma_e) and plot q.
- **Cons:** No upper limit curve; no profiling over background; not directly comparable to pydme limits.
- **Use when:** You want a quick discovery-style q map or relative sensitivity; fixed B is acceptable.

### Option B: Match pydme formula only (no profiling)

- Use the **same Poisson NLL** as now (we already do).
- Add an optional **constrain** term in the NLL if we later have a parametric B(theta).
- Still **no minimization**: B fixed. We only evaluate q at given (m_chi, sigma_e).
- **Result:** Numerical agreement with pydme when pydme is run with fixed theta and no constrain.

### Option C: Add profile likelihood (match pydme fully)

1. **Parametric background:** B_pat = B(theta) (e.g. one theta per pattern or a small set of DC parameters) and optionally constrain(theta).
2. **Minimizer:** Use a C++ minimizer (e.g. Minuit via ROOT/IMinuit) to minimize NLL(data, S(sigma_e) + B(theta)) over theta.
3. **Test statistics:**  
   - **q0:** Fix sigma_e = 0, minimize over theta → NLL_0; q0 = 2(NLL_0 - NLL_min).  
   - **q_mu:** Fix sigma_e = mu, minimize over theta → NLL_mu; q_mu = 2(NLL_mu - NLL_min).
4. **Upper limit:** For each m_chi, find sigma_e such that q_mu(sigma_e) = CL threshold (bracketing + bisection, as in pydme).
- **Result:** Same logic as pydme; direct comparison of limits and discovery q0.

---

## 8. How to do profiling in the Asimov case

In the **Asimov** setup we set the observed counts to the background expectation: **data = B** (one value per bin). Profiling means we let the background be parametrized by **theta** and minimize the NLL over theta for each hypothesis.

### 8.0 Does this match the implementation in pydme?

**Structure: yes.** The procedure (parametric B(theta), minimize NLL over theta for each hypothesis, then form q = 2*(NLL_test - NLL_null) or q_mu = 2*(NLL(sigma_e) - NLL_min)) is the same as in pydme. So the **formulas and profiling steps** match.

**Data: one important difference.** In pydme, **data** is the **real observed counts** `D` (from the experiment, e.g. pattern counts from a CSV). So by default pydme does **not** use Asimov data. The SRDM script passes `D = data['Count_Candidate_{pat}'].values` into `NLLAnalysis`; that is real data.

- If we implement **Asimov profiling** (data = B, then profile over theta), we are doing the **same statistical procedure** as pydme, but with **data = B** instead of real D. So we would match pydme **numerically** only if we ran pydme with its input data set to the expected background B (e.g. replace `D` by `B(theta_ref)` for some reference theta, or by our fixed B_pat).
- If we implement **profile likelihood with real data** (data = D, profile over theta), then we match pydme **exactly** (same data, same NLL, same minimizer logic).

**Summary:** The Asimov profiling recipe in Section 8 matches pydme’s **implementation pattern** (minimize over theta, form q and q_mu the same way). It does **not** use the same **data** as pydme’s default (they use real D). For a numerical match in the Asimov case, run pydme with data = B, or implement the same B(theta), S(sigma_e), and constrain(theta) in C++ and use data = B.

### 8.1 Parametric background

- Define **B(theta)** so that each bin’s expected background is a function of a small set of parameters (e.g. one overall scale, or one scale per pattern).
- Example: **B_i(theta) = theta_0 * b_i** with fixed template **b_i** (e.g. from a DC model). Then theta_0 is the only nuisance parameter.
- Or: **B_i(theta) = theta_i** (one parameter per bin). Then theta is just the vector of expected backgrounds.

### 8.2 Asimov profile likelihood ratio (discovery)

1. **Asimov data:**  
   **data = B** (vector of expected background counts per bin). So we pretend we observed exactly the background.

2. **Null (background only):**  
   Model: **mu = B(theta)**.  
   - Minimize **NLL(B | B(theta))** over theta → get **theta_hat_null**.  
   - With data = B, the best fit is when B(theta_hat_null) = B (if that’s in the parameter space).  
   - **NLL_null = NLL(B | B(theta_hat_null))** (often 0 or minimal if B is in the model).

3. **Test (signal + background):**  
   Model: **mu = S + B(theta)**.  
   - Fix S (from the chosen m_chi, sigma_e).  
   - Minimize **NLL(B | S + B(theta))** over theta → **theta_hat_hat**.  
   - So we fit S+B(theta) to “data” B; the fit pulls B(theta) toward B − S.  
   - **NLL_test = NLL(B | S + B(theta_hat_hat))** (≥ NLL_null).

4. **Profile likelihood ratio (Asimov):**  
   **q = 2 * (NLL_test − NLL_null)**  
   i.e. **q = 2 * (NLL(B | S + B(theta_hat_hat)) − NLL(B | B(theta_hat_null)))**.

So in the Asimov case, “profiling” means: for the null, minimize over theta so that B(theta) fits B; for the test, minimize over theta so that S + B(theta) fits B. Then q is the usual likelihood-ratio statistic with those profiled values.

### 8.3 Implementation sketch (Asimov + profiling)

1. **Define B(theta):**  
   - Either **B(theta) = theta** (one parameter per bin; then theta_hat_null = B exactly, NLL_null = 0).  
   - Or a low-dimensional model, e.g. **B_i(theta) = theta_0 * b_i**, with **b** a fixed template; then minimize over theta_0 (and any other params).

2. **Null fit:**  
   - Minimizer: find theta that minimizes **NLL(data=B | B(theta))** (and optionally add **constrain(theta)**).  
   - Store **theta_hat_null** and **NLL_null**.

3. **Test fit (for a given S):**  
   - Minimizer: find theta that minimizes **NLL(data=B | S + B(theta))** (and optionally constrain(theta)).  
   - Store **NLL_test**.

4. **q:**  
   - **q = 2 * (NLL_test − NLL_null)**.

5. **Optional – Upper limit (still Asimov):**  
   - For each sigma_e (or a grid), S = S(sigma_e). Do the test fit for that S, get NLL_test(sigma_e).  
   - **q_mu(sigma_e) = 2 * (NLL_test(sigma_e) − NLL_min)**, where NLL_min is the global minimum over (sigma_e, theta) (or use NLL_null as reference for the “no signal” case).  
   - Find sigma_e such that q_mu(sigma_e) = CL threshold (e.g. 90% CL).

### 8.4 Relation to current (no profiling) formula

- If **B is not parametrized** (fixed B), then:  
  - **theta_hat_null** is irrelevant; **NLL_null = NLL(B | B)** is the minimum (0 if we define NLL so that data = model gives 0).  
  - There is no theta in the test either, so **NLL_test = NLL(B | S + B)**.  
  - So **q = 2 * NLL(B | S + B)** (with NLL_null = 0), which is exactly the current Asimov formula **q = 2 * sum_i [ S_i − B_i ln(1 + S_i/B_i) ]**.

- Adding profiling changes q only when B is parametrized and the best-fit B(theta) under the two hypotheses differs (e.g. when the test fit can “absorb” some of the signal by adjusting theta).

---

## 9. Implemented profile likelihood (CCDarkSens)

A **profile likelihood** path is implemented to match the pydme procedure (Section 8).

- **Class:** `ccdarksens::stats::ProfileLikelihood` (`include/ccdarksens/stats/ProfileLikelihood.hh`, `src/stats/ProfileLikelihood.cc`).
- **Background model (two modes):**
  - **scale:** B = scale × B_template. `SetBTemplate(B_pat)`; bounds default 0.01–10.
  - **Bp_theta_Br (pydme):** B_i = Bp_i + θ×Br_i. `SetBpBr(Bp, Br)`. Bounds `theta_lo`, `theta_hi` (default 0.5–10). If config does not set Bp/Br, Bp=0 and Br=B_pat (so B = θ×B_pat).
- **Constrain:** Optional `SetConstrain(scale/θ)` (e.g. Gaussian prior). For Bp_theta_Br, **pydme L(θ):** `SetConstrainPriorStrength(N)` adds L = Σ_i (−θ·Br_i + N·ln(θ·Br_i)); use N=98 to match pydme.
- **NLL:** Poisson term + constrain (and pydme L when prior strength > 0).
- **Data:** **Asimov:** `SetData(B_pat)`. **Real:** `SetData(loaded_counts)` from CSV.
- **Config:** `run.use_profile_likelihood`, `run.data_path`, `run.background_model` ("scale" | "Bp_theta_Br"), `run.background_Bp`, `run.background_Br` (arrays, same order as pattern_roi), `run.constrain_prior_strength` (0=off, 98=pydme), `run.theta_lo`, `run.theta_hi`.
- **Apps:** Example and scan use profile q when `use_profile_likelihood` is true; scan also outputs q0, q_mu, upper limit at grid resolution.

---

## 10. How close are we to pydme?

### What already matches

| Aspect | pydme | CCDarkSens (with profile) |
|--------|--------|---------------------------|
| **Poisson NLL** | `model − data·ln(model)` per bin | Same formula. |
| **Data** | Real D or Asimov | Asimov (data = B) or real (CSV via `data_path`). |
| **Profiling** | Minimize NLL over θ for null and test | Minimize over **scale** for null (S=0) and test (S=S_pat). |
| **Discovery q** | q0 = 2(NLL(σ=0, θ̂̂) − NLL_min) | q = 2(NLL_test − NLL_null) with null = S=0, test = S_pat (both profiled). Same idea. |
| **Binned pattern space** | One bin per pattern | Same (pattern_roi). |

So for **discovery-style q** (null vs signal+background, with one nuisance), the **procedure** matches: we profile over a background parameter and form the likelihood ratio. Numerically, results will agree only if the background model and data are equivalent (see gaps below).

### What is different or missing

| Gap | pydme | CCDarkSens |
|-----|--------|------------|
| **q_mu and upper limit** | Free fit over (xsec_e, θ) → NLL_min; q_mu(σ) = 2(NLL(σ) − NLL_min); find σ where q_mu = CL threshold (bracketing + bisection). | Not implemented. We only evaluate q at fixed (m_χ, σ_e) grid points; no “free fit”, no q_mu curve, no upper-limit solving. |
| **Parameter of interest in minimizer** | **xsec_e** (log10) is a Minuit parameter; minimizer varies both xsec_e and θ. | **σ_e** is not a minimizer parameter; S is fixed at each grid point. So we get “q at this σ_e”, not “q_mu(σ_e)” relative to a global NLL_min. |
| **Background parametrization** | B = B(θ) with e.g. one θ per pattern (e.g. B = Bp + θ·Br) and optional **constrain** term L(θ) in NLL. | Single **scale**: B = scale × B_template (one global nuisance). No per-pattern θ, no constrain term used by default. |
| **Constrain term** | NLL includes L (e.g. from `Background_pattern`: prior/penalty on θ). | `SetConstrain(scale)` exists but is optional and not wired in config; not the same form as pydme. |
| **Minimizer** | Minuit (iminuit) over multiple parameters. | 1D Brent over scale only. |
| **Single-bin option** | Can collapse to one Poisson count (`do_single_bin_likelihood`). | Always multi-bin (pattern bins). |

### Summary

- **Discovery q (with profiling):** We are **conceptually aligned** and **implementationally close**: same NLL, same idea of profiling over a background nuisance, Asimov or real data. Remaining differences: we use one global scale instead of pydme’s per-pattern θ and constrain; no q_mu / upper-limit machinery.
- **Upper limits:** We are **not** close: pydme’s upper limit comes from a **free fit** (minimize over xsec_e and θ), then **q_mu(σ)** and **bisection** to find σ at the CL threshold. We would need: (1) a way to minimize over (σ_e, scale) or to define NLL_min from a scan, (2) q_mu(σ_e) = 2(NLL(σ_e) − NLL_min), (3) bracketing + bisection (or similar) to get the exclusion σ_e.

To **numerically** match pydme’s discovery q for a given setup, use the same data (Asimov or real), the same B template, and a comparable background model (e.g. one scale; or extend to per-pattern θ and constrain if needed). To match **limits**, implement the q_mu + upper-limit step above.

---

---

## 11. Checklist: get as close as possible to pydme (keep precomputed grid)

Everything below keeps **(m_χ, σ_e)** as a **precomputed grid** (no σ_e in the minimizer). We emulate pydme’s statistics by defining NLL_min and q_mu from the grid.

### 11.1 Test statistics on the grid

| Item | What to do |
|------|------------|
| **NLL at each grid point** | For each (m_χ, σ_e) we already have S_pat. Minimize NLL over θ (scale or full θ) → store **NLL(σ_e)** and optionally **θ̂(σ_e)**. |
| **NLL_min per mass** | For each m_χ, **NLL_min = min over σ_e (in grid)** of NLL(σ_e). Store **σ_e_at_min** (grid point where minimum is reached). |
| **q_mu(σ_e)** | At each grid point: **q_mu = 2 × (NLL(σ_e) − NLL_min)**. Write this into the same 2D histogram as now (or a dedicated “q_mu” map). |
| **q0** | **q0 = 2 × (NLL(S=0) − NLL_min)**. For each m_χ: one “null” fit (S=0, minimize over θ) → NLL_0; then q0 = 2×(NLL_0 − NLL_min). |

**Where:** Scan app (and optionally example app). After the grid loop over (m_χ, σ_e), add a per–m_χ pass: (1) find NLL_min and σ_e_at_min over the σ_e grid; (2) compute one null fit for q0; (3) when writing q, either store q_mu instead of “discovery q” or store both.

### 11.2 Upper limit (grid-only or with S interpolation)

| Item | What to do |
|------|------------|
| **Target q** | **target_q = (norm_ppf(CL))²** (e.g. 90% CL → ≈ 2.71). Same as pydme’s 2×ΔNLL threshold. |
| **Limit at grid resolution** | For each m_χ: **upper limit = smallest grid σ_e** such that q_mu(σ_e) ≥ target_q. If no such point, use largest grid σ_e or mark “no limit”. Resolution = grid step. |
| **Finer limit (optional)** | Build **interpolant S(m_χ, σ_e)** from the existing grid. For limit finding: **bracketing + bisection** in σ_e; at each trial σ_e evaluate S via interpolant, minimize NLL over θ, get q_mu. No extra pipeline runs; limit resolution set by bisection tolerance. |

### 11.3 Background and nuisance parameters

| Item | What to do |
|------|------------|
| **Keep one global scale** | Current **B = scale × B_template** is enough for “one nuisance”. No code change if you don’t need per-pattern. |
| **Per-pattern θ (pydme-like)** | **B_i(θ) = θ_i** (one θ per pattern bin) or **B_i = Bp_i + θ_i×Br_i** (template + per-pattern scale). Requires **vector θ** and a **multi-D minimizer** (e.g. ROOT Minuit2). |
| **Constrain term L(θ)** | Add **L(θ)** to NLL (e.g. -ln prior(θ) or pydme’s Background_pattern formula). API: e.g. `SetConstrain(std::function<double(const std::vector<double>& theta)>)`. Config: optional path to a formula or coefficients (e.g. Br_i, prior strength). |
| **Minimizer** | 1D (scale only): keep **Brent**. Multi-D (per-pattern θ): add **Minuit2** (ROOT) or equivalent; same NLL, just minimize over (θ_1,…,θ_n). |

### 11.4 Single-bin option

| Item | What to do |
|------|------------|
| **Config** | e.g. **run.single_bin_likelihood** (bool, default false). |
| **Logic** | If true: sum **data** over pattern bins → one count; sum **model** (S + B(θ)) over pattern bins → one expectation; **one Poisson term** in NLL. Still minimize over θ (scale or full θ). q_mu, q0, and limit unchanged in definition, just with a single bin. |

### 11.5 Config and output

| Item | What to do |
|------|------------|
| **Confidence level** | **run.cl** (e.g. 0.9) already exists; use it to compute **target_q** for the upper limit. |
| **Output** | Per m_χ: **q0**, **NLL_min**, **σ_e_at_min**, **upper limit** (grid or bisection). Optional: **q_mu(m_χ, σ_e)** 2D map (or overwrite current q map with q_mu). |

### 11.6 Summary order of implementation

1. **NLL_min and q_mu on the grid** (11.1) — **Done.** Scan app: per m_χ, null fit + collect NLL(σ_e), then NLL_min, q_mu, q0.
2. **Upper limit at grid resolution** (11.2) — **Done.** target_q = (NormQuantile(CL))²; smallest grid σ_e with q_mu ≥ target_q; written to `upper_limit_sigma_e_mchi`.
3. **Constrain term** (11.3) — **Done.** Config `constrain_scale_prior_mean` and `constrain_scale_prior_sigma` (optional); Gaussian prior on scale.
4. **Single-bin option** (11.4) — **Done.** Config `single_bin_likelihood`; collapse data and B (and S) to one bin in scan and example apps.
5. **Per-pattern θ + Minuit2** (11.3) — Not implemented; would require multi-D minimizer for B_i(θ).
6. **Finer upper limit** (11.2, optional) — Not implemented; would require S(σ_e) interpolant and bisection.

This gets you as close as possible to pydme **while keeping the precomputed (m_χ, σ_e) grid**: same NLL form, same q_mu/q0/limit definitions, optional per-pattern θ and constrain, optional single-bin; only difference is NLL_min and the limit are derived from the grid (and optionally from an interpolated S(σ_e)) instead of from a minimizer that varies σ_e.

---

### Recommended next steps

1. **Short term:** Keep **Option A** for the full scan: produce a **q(m_chi, sigma_e)** map (and optionally a simple sensitivity curve by defining a q threshold and reading off sigma_e at each m_chi). Document that this is Asimov, fixed B, no upper limit.
2. **Comparison:** Run pydme with **fixed theta** (and, if possible, same B and S) and compare **q** from both frameworks to check that the Poisson term and ratio match when there is no profiling.
3. **Medium term:** If you need upper limits and/or comparison with pydme limits, implement **Option C** (or Option B first, then Option C): parametric B(theta), minimizer, q_mu, and upper-limit finding.
