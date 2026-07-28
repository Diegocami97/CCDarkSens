# Upper limit computation: pydme vs CCDarkSens

Comparison of how the 90% CL upper limit on σ_e is computed in pydme (`compute_upper_limit_fast`) and in CCDarkSens (pydme-mode scan).

---

## pydme (`dme_nLL_exclusion.compute_upper_limit_fast`)

1. **Free fit**  
   Full minimization over **x0** (= log10(xsec_e)) and **theta** (and any other params). Gives `nll_min` and `xsec_at_nll_min` (MLE in log10 space).

2. **target_q**  
   `target_q = 2 * delta_nll_threshold` with `delta_nll_threshold = 0.5 * (NormQuantile(CL)^2)` (same as CCDarkSens, 90% CL → target_q ≈ 2.71).

3. **Bracketing (right of MLE)**  
   - `lo = max(xsec_at_nll_min, theta_limits['xsec_e'][0])`  
   - Step **right** in log10(xsec): `hi = min(hi + ul_brack_step, ub)` with `ul_brack_step` (0.3), doubled each iteration.  
   - At each `hi` call `_qmu_at(hi)`: fix x0=hi, minimize over theta (Minuit migrad), return `q_mu = 2*(nll - nll_min)`.  
   - Stop when `q_hi >= target_q` or hit upper bound.

4. **Bisection**  
   - Interval `[left, right]` with `q_mu(left) < target_q <= q_mu(right)`.  
   - `mid = (left+right)/2`, `q_mid, nll_mid = _qmu_at(mid)`.  
   - If `q_mid >= target_q` → `right = mid`, else `left = mid`.  
   - Repeat until `|right - left| < tol_x` or `|q_mid - target_q| < ul_q_tol`.  
   - **Upper limit (log10)** = `right`; in cm² they use `10**right`.

5. **_qmu_at(x)**  
   Fix x0 (log10(xsec)) = x, run Minuit to minimize over theta only, return `q_mu = 2*(nll - nll_min)`. So at any x they get **exact** S(x) from their model (rate tables / interpolation) and profile over theta.

---

## CCDarkSens (pydme mode, after the fix)

1. **No global free fit in log10(σ)**  
   We have a **sigma grid** and (optionally) a 2D Minuit over (log10_sigma, theta). If the 2D fit is accepted we get `log10_sigma_hat_pydme`; if not we only have the grid.

2. **target_q**  
   Same: `target_q = NormQuantile(CL)^2` (so 2 * 0.5 * …), 90% CL.

3. **Bracketing**  
   - Start: `lo = use_2d_nll_min ? max(log10_sigma_hat_pydme, log10_lo, log10_sigma_best) : max(log10_lo, log10_sigma_best)` (grid point with smallest NLL).  
   - Step right: `hi = min(hi + step, log10_hi)`, step 0.3, doubled each time.  
   - At each point we call `q_mu_at(log10_sig)`: **interpolate** S from `S_grid` at that log10_sigma, then `MinimizeOverScaleMinuit(S, theta_lo, theta_hi)` (profile over theta), return `2*(nll - nll_min)`.

4. **Bisection**  
   Same idea: bisect in log10(σ), at each `mid` interpolate S(mid) from S_grid, profile over theta, get q_mu; move left/right until converged.  
   **Upper limit** = `10^right` (continuous value in cm²).

5. **nll_min**  
   We use **nll_min_grid** = min over the sigma grid of the profile NLL (not a free 2D minimum). So the denominator of q_mu is the best NLL on the grid, to avoid spurious 2D minima flattening the curve.

---

## Summary: same idea, small differences

| Aspect | pydme | CCDarkSens |
|--------|--------|------------|
| **Algorithm** | Bracketing right of MLE + bisection in log10(xsec) | Same: bracket right of start + bisection in log10(σ) |
| **target_q** | 2 × (0.5 × NormQuantile(CL)²) | Same |
| **At each log10(σ)** | Fix x0, minimize over theta (Minuit), q_mu = 2×(NLL−nll_min) | Interpolate S from grid, minimize over theta (Minuit/Brent), same q_mu |
| **Start for bracket** | MLE from free fit (xsec_at_nll_min) | 2D fit result or grid best (log10_sigma_best) |
| **S at arbitrary σ** | Exact from model (x0 → S inside NLL) | Linear interpolation in log10(σ) between S_grid points |
| **nll_min** | Global minimum from free fit | Minimum over sigma grid of profile NLL |

So we do the **same** bracketing + bisection and the **same** profiling over theta at each σ. The differences are:

- **Signal**: pydme evaluates S exactly at any log10(xsec); we use **interpolated** S from the grid (negligible if the grid is dense).
- **Starting point**: we use the grid best when the 2D fit is rejected (e.g. theta at boundary); pydme always has a free fit so they always have xsec_at_nll_min.
- **Denominator of q_mu**: we use the grid minimum on purpose so the curve is not flattened by a bad 2D minimum; pydme uses the global minimum. So we are slightly more conservative.

Overall the method is the same; the limit curve should match closely up to interpolation and the nll_min choice.

---

## Matching pydme exactly: `pydme_style_ul`

Set **`"pydme_style_ul": true`** in the run config (e.g. `scan_dmelectron_pattern_pydme_minuit.json`) to:

1. **Accept 2D fit at boundary**  
   The 2D (log10_sigma, theta) minimizer is no longer rejected when the minimum lies on theta_lo/theta_hi (or sigma bounds). So we get a valid MLE and bracket start even when theta is at 0.5, like pydme.

2. **Use 2D nll_min for the UL**  
   When the 2D fit is used, the denominator of q_mu and the UL is **nll_min_2d** (global minimum from the 2D fit) instead of the grid minimum. So the same formula as pydme: q_mu = 2*(NLL − nll_min_global).

With `pydme_style_ul: true`, bracket start and nll_min match pydme. The only remaining difference is S(σ) at bisection points (we interpolate from the grid; pydme evaluates the model). For a dense sigma grid this is negligible.
