# How pydme does the profile-likelihood minimization

This document describes the minimization flow in **pydme** (pattern analysis) so it can be compared with CCDarkSens. The main code is in `collab_frameworks/pydme/pydme/dme_nLL_exclusion.py` and `collab_frameworks/pydme/pydme/detector/background_models/background_pattern.py`.

---

## 1. Parameters and bounds

**Minimizer:** iminuit (Minuit).

**Parameters (pattern analysis):**

| Parameter | Meaning | Default bounds (theta_limits) |
|-----------|---------|------------------------------|
| **x0** (xsec_e) | log₁₀(cross-section) in cm² | `xsec_e_lims` (e.g. [-35, -29.5]) |
| **x1** (theta_0) | Single nuisance for background | **[[1e-6, 0.1]]** (hardcoded in `NLLAnalysis.__init__`) |

So pydme minimizes over **two** parameters at once: **log₁₀(σ_e)** and **θ**. The internal variable for the cross-section is **log₁₀(xsec_e)**; the background uses **θ** (one number for pattern analysis).

---

## 2. Background model (pattern)

**File:** `pydme/detector/background_models/background_pattern.py`  
**Function:** `Background_pattern(theta, pattern, gamma, ...)`

- **theta:** list; for the standard pattern analysis **theta[0]** is the single global θ.
- **Bp, Br:** hardcoded per pattern (same 6 patterns and values as in our config):
  - Bp = [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]
  - Br = [0.039, 0.039, 0.016, 0.052, 0.011, 0.035]
- **Per-pattern:** B = Bp[row] + **theta[0]** * Br[row] (then normalized by len(gamma)).
- **Constrain (prior):** L = -theta[0]*Br[row] + **98** * log(theta[0]*Br[row]) — same functional form as our `constrain_prior_strength=98`.

So the background model is **B = Bp + θ·Br** (one θ for all patterns), and the prior term uses 98.

---

## 3. Model and NLL

**Model:** For pattern analysis, `Model_loader` calls `Model_pattern`, which:

1. **Signal S:** If xsec_e &lt; 0 (interpreted as log₁₀), xsec_e → 10^xsec_e and signal is loaded (SHM, pattern rates × exposure). If xsec_e ≥ 0, treated as “no signal” (S = 0) for the q_μ=0 case.
2. **Background:** For each pattern, `Background_pattern(theta, pattern, N_pix)` returns B and L (constrain).
3. **Expected rate:** R = s + b (per pattern), then NLL uses **model = R**, **data** = observed (or Asimov) counts.

**NLL (Poisson):** `nLL_Poisson(self, **params)` in `dme_nLL_exclusion.py` (lines 539–589):

- Unpacks **xsec_e** and **theta** from `params` (Minuit passes these).
- For each dataset i, builds **theta** array and **pattern** array for that dataset.
- Calls **Model_loader(mX, xsec_e, ..., theta, ...)** → returns **Nim** (model prediction per bin) and **l** (constrain per bin).
- **NLL = Σ ( model - data*ln(model) ) + Σ constrain**  
  Same as our formula: Poisson term plus pydme prior term.

So the **objective** is NLL(xsec_e, θ); Minuit minimizes it over **both** parameters.

---

## 4. Flow: full fit and nll_min

**Method:** `minimize_nll(self, theta)` (lines 590–646)

1. Build **param_dict**: `xsec_e` = self.xsec_e (initial log₁₀(xsec)), `theta_0` = self.theta[0], ...
2. **Minuit(wrapped_nLL, *param_dict.values())** with **errordef = Minuit.LIKELIHOOD**, **strategy**, **print_level**, **tol = 0.1**.
3. **Limits:**  
   `migrad.limits["x0"]` = xsec_e_lims  
   `migrad.limits["x1"]` = theta limits = **[[1e-6, 0.1]]** (one theta).
4. **migrad.migrad()** → minimizes NLL over **xsec_e** and **theta** together.
5. **nll_min** = migrad.fval, **xsec_at_nll_min** = migrad.values['x0'] (log₁₀), best-fit theta from migrad.values.

So the **global minimum** is found in one Minuit run over (log₁₀(σ_e), θ). There is no separate “null” run for nll_min; the null (S=0) is used only for **q0** and for **q_μ** at fixed σ.

---

## 5. Flow: q0 (null = background only)

**In compute_upper_limit_fast (lines 669–696):**

1. After step 1 (minimize_nll → nll_min, xsec_at_nll_min):
2. **Fix x0 to 0:** `self.migrad.fixto('x0', 0.0)`  
   So log₁₀(xsec) = 0 → xsec_e = 10^0 = 1 in linear space, but in pydme the convention for “no signal” may be a special value; in Model_pattern, **xsec_e ≥ 0** gives S = 0 (see line 129–132). So fixing x0 to **0** gives **no signal** (null = background only).
3. **migrad.migrad()** → minimizes over **theta only** (x0 fixed).
4. **nll_0** = migrad.fval.
5. **tmu0 (q0)** = 2 * (nll_0 - nll_min).

So the **null case is background only**: fix signal (x0=0 → S=0), minimize over θ, then q0 = 2*(nll_0 - nll_min).

---

## 6. Flow: q_μ at fixed σ (for upper limit)

**Method:** `_qmu_at(self, x, strategy=0)` (lines 648–667)

- **x** = value of **x0** = log₁₀(xsec_e) at which to evaluate q_μ.
1. **Fix x0 to x:** `self.migrad.fixto('x0', x)`.
2. **migrad.migrad()** → minimizes over **theta only** (x0 fixed at x).
3. **nll** = migrad.fval.
4. **q_μ** = 2 * (nll - self.nll_min).
5. Return (q_μ, nll).

So at each **fixed** cross-section (log₁₀), pydme **profiles over θ only**, then forms q_μ = 2*(nll - nll_min). Same idea as us: at fixed S(σ), minimize NLL over θ, then compare to global nll_min.

---

## 7. Upper limit (fast mode)

**Method:** `compute_upper_limit_fast` (lines 669–775)

1. **minimize_nll** → get **nll_min**, **xsec_at_nll_min** (both parameters free).
2. **target_q** = 2 * (0.5 * norm.ppf(CL)^2) = norm.ppf(CL)^2 (e.g. 2.71 for 90% CL).
3. **q0:** fix x0=0, migrad → nll_0, tmu0 = 2*(nll_0 - nll_min).
4. **Bracketing:** Start at lo = xsec_at_nll_min. Step **hi** to the right in log₁₀(xsec) (step 0.3, then doubled) until **_qmu_at(hi)** ≥ target_q.
5. **Bisection:** Bisect in log₁₀(xsec) until |q_mu - target_q| &lt; q_tol or interval small; each iteration calls **_qmu_at(mid)** (fix x0, minimize over θ).
6. **Upper limit** = the xsec_e (in log₁₀) where q_μ = target_q.

So the UL is found by **bracketing + bisection in log₁₀(σ_e)**; at each evaluation, **θ is profiled** (Minuit with x0 fixed).

---

## 8. Summary: pydme vs CCDarkSens

| Aspect | pydme | CCDarkSens (grid: brent/minuit) | CCDarkSens (minuit2d) |
|--------|--------|---------------------------------|------------------------|
| **Parameters in minimizer** | **Two:** log₁₀(xsec_e) and θ. One Minuit, both free. | **One:** θ only. σ_e on a **grid**; no minimization over σ. | **Two:** one 2D Minuit over (log₁₀(σ_e), θ) per mχ; then grid for q_μ. |
| **nll_min** | One **migrad()** (both free) → global minimum. | min over **grid** of (min_θ NLL(S(σ), θ)). | min(2D fit, grid min); 2D can sit between grid points. |
| **Null (S=0)** | **fixto('x0', 0)** then migrad() over θ → nll_0. | S = 0 vector, minimize over θ → nll_null. | Same: S = 0, min over θ. |
| **q_μ at fixed σ** | **fixto('x0', x)** then **migrad()** over θ → nll; q_μ = 2*(nll − nll_min). | At each grid σ: 1D min over θ (Brent/Minuit) → nll(σ); q_μ = 2*(nll(σ) − nll_min). | Same: at each grid σ, 1D min over θ; q_μ = 2*(nll(σ) − nll_min). |
| **Upper limit** | **Bracketing + bisection** in log₁₀(σ); at each trial σ call _qmu_at(x) (fix x0, migrad). UL = crossing with target_q. | Smallest **grid** σ (no refinement). | **minuit2d**: first grid σ. **pydme**: **bracketing + bisection** in log₁₀(σ), matching pydme. |
| **θ bounds** | **[[1e-6, 0.1]]** in code; scripts often override to e.g. [0.5, 10]. | **theta_lo, theta_hi** from config. | Same config. |
| **Algorithm** | Minuit **Migrad** (gradient). | Brent or Minuit2 **Simplex** (1D over θ). | Minuit2 **Simplex** for 2D and 1D profile. |

So: **pydme** does one 2D fit, then **profiles by fixing x0 and running migrad**; UL is from **bracketing + bisection**, not a fixed grid. **We** use a **σ grid** and 1D minimization over θ at each point; in **minuit2d** we add a 2D fit only to define nll_min (and optionally use it), but we still take the UL as the **first grid point** above the threshold (no bisection).

---

## 9. Where the θ-bound issue can come from

- In pydme, **theta is in [1e-6, 0.1]** by default. So the **nominal** B = Bp + Br (θ=1) is **outside** that range. For Asimov data = B, the best-fit θ would be 1, so pydme’s default bounds are **wrong** for that nominal. (They may be intended for a different parametrization or a different analysis.)
- In CCDarkSens, we use **configurable** theta_lo, theta_hi. For Bp+Br with nominal θ=1 we need **[0.5, 10]** or **[0.1, 100]** etc. so that the minimizer can reach θ ≈ 1. If we used [1e-6, 0.1] we would hit the same problem as pydme’s default: best fit at boundary, flat q_μ.
- When you add efficiency for pattern 22 @ ne=4, the best-fit θ can **move** (e.g. above 10); then **theta_hi = 10** is too low and we hit the boundary again. So we need **wider bounds** (e.g. theta_hi = 100) when including that efficiency.

---

## 10. pydme flow diagram

```
NLLAnalysis.__init__:
  theta_limits = { 'xsec_e': xsec_e_lims, 'theta': [[1e-6, 0.1]] }
  nLL = nLL_Poisson (model = S(xsec_e) + B(theta), data, constrain = 98 term)

minimize_nll():
  Minuit(nLL), params = xsec_e (log10), theta_0
  limits: x0 in xsec_e_lims, x1 in [1e-6, 0.1]
  migrad() → nll_min, xsec_at_nll_min, theta_hat

q0:
  fixto('x0', 0)  [S=0]
  migrad() → nll_0
  q0 = 2*(nll_0 - nll_min)

q_μ(x) at fixed log10(xsec)=x:
  fixto('x0', x)
  migrad() → nll
  q_μ = 2*(nll - nll_min)

Upper limit (fast):
  Bracketing: find hi such that q_μ(hi) >= target_q
  Bisection: refine until q_μ(mid) ≈ target_q
  UL = that log10(xsec) value
```

---

## 11. Matching the reference (DAMIC-M Pattern) limit curve

If your CCDarkSens limit curve (dashed) sits **above** the target (solid) by a roughly constant factor (e.g. ~1.5–2×) for m_χ ≳ 1 MeV, the limit is less stringent than the reference. These knobs affect the curve and should be aligned with the reference setup:

| Knob | Effect | What to check |
|------|--------|----------------|
| **Exposure** | Signal ∝ exposure; UL on σ_e ∝ 1/exposure. **Higher exposure → lower (stricter) limit.** | `experiment.livetime_days`, `detector.mass_kg`. Exposure in code: `exposure_kg_year = livetime_days * duty_cycle * mass_kg / 365.25`. If the reference used more kg·year, increase livetime or mass to match. |
| **CL** | `run.cl` (e.g. 0.9) → target_q = NormQuantile(CL)². Same CL as reference for a fair comparison. | Usually 90% CL (target_q ≈ 1.64). |
| **Background (Bp, Br, θ)** | Different B or θ bounds change nll_min and the shape of q_μ. | Same Bp/Br and θ bounds as pydme/reference (e.g. theta_lo=0.5, theta_hi=10 when nominal θ=1). |
| **Pattern ROI & efficiency** | Different patterns or ε(pattern\|n_e) change S and B. | Same `pattern_roi` and efficiency table as the reference. |
| **Rates (dR/dE)** | Same rate tables and σ grid as reference. | `model.rates_dir`, `model.filename_template`, grid `sigma_e_cm2`. |

**Practical step:** Get the **exposure in kg·year** used for the solid curve. Set your config so that `livetime_days * duty_cycle * mass_kg / 365.25` equals that value (or scale one of livetime/mass). Re-run the scan; if the offset was mainly from exposure, your curve should move down toward the solid one.
