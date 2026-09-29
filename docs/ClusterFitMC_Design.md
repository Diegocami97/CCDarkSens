# ClusterFitMC — Design, Derivation, and Validation
## Phase 3 (noise-tail ΔLL calibration) and Phase 4 (per-event reconstruction) of the WIMP-nucleus SI channel

**Document:** `docs/ClusterFitMC_Design.md`
**Author:** Diego Venegas-Vargas
**Status:** Living document, built up slice by slice alongside implementation. Each section below is added once its slice is built, tested, and reviewed — not written retroactively from memory.

---

## 0. Scope

`ClusterFitMC` is the reconstruction module for the WIMP-nucleus SI scattering channel (see the channel's overall status/plan: republish target of the artifact tracked in project memory). It replaces the discrete `n_e ≤ 5` pattern-classifier machinery (`EfficiencyMC`, `PatternClassifier`) that the DM-electron/dark-photon/Migdal channels share, which cannot extend to the hundreds-of-electrons clusters a quenched nuclear recoil produces. Target reproduction: PhysRevD.94.082006 (DAMIC 0.6 kg·day, SNOLAB), Fig. 6 (ΔLL tail calibration) and Fig. 9/11 (efficiency, limit).

`ClusterFitMC` is a distinct module from the deleted `ClusterMC` pixel-simulation backend documented in `CLAUDE.md`'s "Files that no longer exist" table — unrelated purpose, unrelated code, new as of this build. `ClusterMCJSON`/`.cmc` remain only as a legacy JSON-read shim for old configs and have no relationship to this module.

The full multi-slice implementation plan (file layout, config schema, build order, validation strategy) lives in the approved plan at `~/.claude/plans/goofy-wishing-bengio.md` (session-local, not part of this repo). This document is the durable, in-repo record of *why* each piece works the way it does — the derivations and validation results, not the project-management scaffolding.

---

## 1. Slice 1 — `ClusterFitModel`: the pixel-window scoring function

### 1.1 Plain-language picture

A particle interaction knocks loose a cloud of electrons at one spot in the detector. That cloud spreads a little as it drifts to the readout surface, so the pixel image looks like a soft, bell-curve-shaped blob rather than a single bright dot:

```
        pixel columns →
      .   .   .   .   .
      .   .  0.1  .   .
      .  0.3 0.9 0.4  .      ← brighter pixels near the
      .   .  0.6  .   .        center of the cloud
      .   .   .   .   .
```

Three unknowns describe that blob: where it was centered `(μx, μy)`, how wide it spread `(σxy)`, and how much total charge was in it (which tells us the deposited energy). `ClusterFitModel` does not search for these — it only answers, very cheaply, one narrower question:

> *If* the cloud were centered here, with this width — how well does that explain the pixels actually measured?

```
 guess: center (μx,μy) + width σxy
              │
              ▼
   ┌────────────────────┐
   │  Shape(pixel, θ)     │  "what fraction of the cloud's charge
   └────────────────────┘   should land in THIS pixel, given θ?"
              │
              ▼
   ┌────────────────────┐
   │       IHat(θ)        │  "given that expected shape, what's the
   └────────────────────┘   best guess for the TOTAL charge?"
              │              (closed-form shortcut — no search needed)
              ▼
   ┌────────────────────┐
   │     ObjectiveG(θ)    │  one number: how good is this whole guess?
   └────────────────────┘   (more negative = better match)
              │
              ▼
   ┌────────────────────┐
   │   DeltaLLFromG(·)     │  converts that number into ΔLL, the actual
   └────────────────────┘   statistic used for the detection cut (Slice 3+)
```

Only `θ = (μx, μy, σxy)` is ever searched over (Slice 2) — the total-charge unknown drops out via a formula, derived below.

### 1.2 Pixel-integrated shape, not point-sampled

`PixelIntegral1D(center, μ, σ)` is the exact fraction of a 1D `Gaussian(μ,σ)`'s probability mass landing inside one pixel (the interval `[center-0.5, center+0.5)`), via the erf-based CDF difference:

```
PixelIntegral1D(c, μ, σ) = Φ((c+0.5-μ)/σ) − Φ((c-0.5-μ)/σ),   Φ(x) = ½[1+erf(x/√2)]
```

identical formula/idiom to `Diffusion::normal_cdf_`/`interval_prob_` (`src/response/Diffusion.cc:43-53`), already used elsewhere in this codebase for turning a Gaussian into a binned probability. The 2D shape is the product of two such 1D integrals (independent x/y for an isotropic Gaussian):

```
Shape(ix, iy, μx, μy, σ) = PixelIntegral1D(ix, μx, σ) · PixelIntegral1D(iy, μy, σ)
```

This is deliberately **pixel-integrated**, not a Gaussian density evaluated at pixel centers: `PixelSimulator::DepositElectron` bins each individual electron into its nearest pixel, so for a cloud of many electrons drawn i.i.d. from `Gaussian(μ,σ)`, the *expected* pixel occupancy is exactly this pixel-integrated mass — the correct large-`n_e` limit of the actual forward simulator. A point-sampled density would be a biased approximation, worst exactly where the fiducial cut lives (`σxy` near the sub-pixel floor).

### 1.3 Closed-form amplitude — why the search only needs 3 parameters, not 4

Pixel model: `data_ij = I·Shape_ij(θ) + noise_ij`, `noise_ij ~ Gaussian(0, σ_pix)` i.i.d., `σ_pix` = the readout-noise constant (`sigma_readout_e`, same one used throughout `response/`). `Shape_ij(θ)` is a fixed, known number once `θ` is fixed, so `I` enters the model *linearly*. Its maximum-likelihood value at fixed `θ` is the standard weighted-least-squares solution:

```
IHat(θ) = Σ_ij(data_ij · Shape_ij(θ)) / Σ_ij(Shape_ij(θ)²)
```

No search needed — plugging in any `θ` immediately gives the best possible `I` for it.

### 1.4 The ΔLL objective — the derivation that makes the whole thing cheap

With `ln L_n = −N/2·ln(2πσ_pix²) − Σdata_ij²/(2σ_pix²)` (white-noise null, no free parameters) and `ln L_G(θ,I) = −N/2·ln(2πσ_pix²) − Σ(data_ij − I·Shape_ij(θ))²/(2σ_pix²)` (signal hypothesis), substituting `I = IHat(θ)` and using the normal-equations identity `IHat·Σ(data·Shape) = IHat²·Σ Shape²`:

```
ln L_G(θ, IHat(θ)) − ln L_n  =  IHat(θ)² · Σ_ij Shape_ij(θ)²  /  (2σ_pix²)     ≥ 0

ΔLL(θ) = −[ln L_G,max − ln L_n] = −IHat(θ)² · Σ Shape_ij(θ)² / (2σ_pix²)      ≤ 0
```

Both hypotheses' normalization constants (`N/2·ln(2πσ_pix²)`) cancel exactly. **`L_n` never needs to be evaluated as a likelihood at all** — the entire test reduces to minimizing one 3-parameter function:

```
ObjectiveG(θ) = −IHat(θ)² · Σ_ij Shape_ij(θ)²          (what Slice 2's minimizer actually minimizes)
DeltaLLFromG(g, σ_pix) = g / (2σ_pix²)                  (the final conversion, applied once, after minimization)
```

This is the result that makes running this many thousands of times (once per noise toy in Phase 3, once per simulated event in Phase 4) computationally viable: each evaluation is one pass over the pixel window with no nested inner fit.

### 1.5 Validation (smoke-tested before Slice 2 was built)

Checked with a standalone scratch program (not part of the repo) before building the search on top, specifically so a bug here couldn't hide behind the minimizer's convergence behavior later:

| Check | Result |
|---|---|
| Pixel-integrated mass sums to 1 over a wide window | `1.0000000000` |
| Near-zero-width Gaussian puts ~all mass in the nearest pixel | `1.0000000000` |
| `IHat` recovers an injected amplitude (250 e⁻) exactly, noiseless data | `250.000000` vs `250.000000` |
| `ObjectiveG` at the true parameters vs. an independent hand-derived continuum formula (`G ≈ −I²/(4πσ²)`) | rel. diff `6.5%` (expected — continuum approximation to a discrete sum) |
| Zero-signal input gives `G = 0` exactly | `0.0000000000` |
| `ObjectiveG` at the true `(μx,μy,σ)` is better (more negative) than at an offset guess | `G(true) = −3843.09` vs. `G(offset) = −343.86` |

The last row is the one that actually matters for Slice 2: it confirms the objective function has a genuine minimum at the true parameters, which is what a downhill search can converge to.

---

## 2. Slice 2 — `ClusterFitEngine`: the search

### 2.1 Plain-language picture

Slice 1 is a scorer; Slice 2 is the searcher that uses it. Picture four blindfolded hikers dropped onto a landscape where "low ground" = a good `ObjectiveG` score:

```
        ✗ (worst hiker)
       ╱
      ╱   4 hikers form a little
     ╱    triangle-ish shape on
    ╳     the landscape (a "simplex")
   ╱ ╲
  ╳   ╳
```

Each round, whichever hiker stands on the worst (highest) spot leaps over the average of the other three, hoping to land somewhere lower; a successful leap gets bigger next round, an unsuccessful one gets smaller and more cautious. This is the **Nelder-Mead** ("downhill simplex") method — after 50–200 such rounds the four hikers converge on the true `(μx, μy, σxy)`, and Slice 1's formula immediately gives the corresponding best-fit total charge.

Because each "how high am I" question (one `ObjectiveG` call) is cheap and closed-form (Slice 1, §1.4), it's practical to implement this search from scratch (~80 lines, no external dependency) rather than reach for a heavier general-purpose tool for every single fit.

That homemade search is cross-checked, not blindly trusted: this repo already uses a well-established minimizer (Minuit2, via `ROOT::Math`) elsewhere (`src/stats/ProfileLikelihood.cc`), and the same landscape will be run through it on a sample of test cases to confirm agreement.

Starting position matters — dropping the hikers randomly risks settling on a false secondary dip. Before the search starts: compute the brightness-weighted centroid of the pixel window as the starting `(μx,μy)`, and do a cheap scan over a handful of candidate widths at that centroid to pick a starting `σxy`.

### 2.2 Bounding σxy: the bug the Minuit2 cross-check caught

The first implementation bounded `σxy` by *clamping*: the search was free to propose any `σxy`, but the objective was evaluated at `clamp(σxy, lo, hi)`. This is a real bug, not just an approximation — clamping makes the objective **perfectly flat** for every proposed value beyond the bound (they all map to the same clamped evaluation), which removes any gradient the simplex could use to climb back down. Cross-checking against Minuit2 caught it immediately: on a test case with the true `σxy` close to the upper search bound, Nelder-Mead converged to exactly `σxy = 2.0000` (the hard boundary) with `ΔLL = −4487`, while Minuit2 on the *same window* found `σxy = 1.8285` (near the true 1.8) with a genuinely better `ΔLL = −4520`. This is exactly the failure mode cross-validation exists to catch — both methods report a self-consistent, "converged" answer; only comparing them against each other exposed that one was wrong.

**Fix, step 1 — reparametrize instead of clamp.** The search variable for width is no longer `σxy` directly but an unconstrained surrogate `u`, mapped through a sigmoid: `σxy(u) = lo + (hi−lo)/(1+e⁻ᵘ)`. Every real number `u` maps to a valid `σxy ∈ (lo,hi)` — there is no out-of-bounds region to clamp into, and the map is smooth and strictly monotonic everywhere, so a gradient always exists.

This alone reduced the failure rate substantially but did not eliminate it (from failing on essentially every trial of the boundary-adjacent test case, to roughly 1 in 5). The residual cause: the sigmoid *saturates* — its slope `dσ/du → 0` as `u → ±∞` — so near the bounds, the search's own convergence check (`spread of objective values across the simplex is tiny`) can be satisfied on a genuinely flat stretch of the *transformed* landscape, before the untransformed `σxy` has actually settled at the true optimum.

**Fix, step 2 — multi-start, and a stricter convergence check.** Two changes: (a) each fit now runs three independent simplex searches from different starting widths (the seed, and two perturbed starts) and keeps whichever converges to the best objective — cheap here since each run is a few hundred evaluations of a closed-form function; (b) the convergence check now requires the simplex's own physical spread (not just its objective-value spread) to have shrunk, so a flat-region false convergence can't terminate the search early. After both fixes, the same test case that previously failed on ~20% of trials matched Minuit2 to within `0.0001 px` in `σxy` and `0.005 e⁻` in amplitude — full agreement.

### 2.3 Validation results

`apps/ccdarksens_validate_cluster_fit_engine.cc`: known `(μx, μy, σxy, I)` injected onto synthetic pixel windows (15×15, `σ_pix = 0.16 e⁻`) at 5 points spanning bright/faint and narrow/wide (including deliberately close to both fiducial edges), 200 independent noise realizations each, fixed seed for reproducibility.

| Check | Result |
|---|---|
| Bias significance (`\|mean residual\| / SEM`) across all 5 test points × 4 parameters × both methods | all `< 2.2σ` (threshold: `5σ`) — no significant bias anywhere |
| Nelder-Mead vs. Minuit2 cross-check, max over all 1000 fits | `Δμ = 0.0032 px`, `Δσxy = 0.0001 px`, `ΔI = 0.0053 e⁻` |
| Runtime, full suite (both methods × all test points × 200 repeats each, multi-start ×3) | ~5.3 s |

`--fit_method nelder_mead\|minuit2\|both` (default `both`) lets either method be run in isolation for faster iteration while debugging; the cross-check only runs when both are active.

---

## 3. Slice 3 — `NoiseTailCalibrator`: empirical detection-cut calibration

### 3.1 Plain-language picture

Even pure noise, fit with the Slices 1–2 machinery, will always produce *some* best-fit blob — the search always returns its best attempt, whether or not anything real is there. `NoiseTailCalibrator` answers: how good can noise alone manage to look, purely by chance, and how often?

```
  100,000 times:
     ┌─────────────────────┐
     │  generate a pixel     │   ← just noise: Reset() + AddReadoutNoise()
     │  window of PURE NOISE │      (no fake electrons deposited at all —
     │  (nothing planted)     │       this already exists, reused as-is)
     └─────────────────────┘
                │
                ▼
     ┌─────────────────────┐
     │  run it through the   │   ← the exact same engine from Slices 1–2,
     │  Slice 1–2 fit engine │      completely unchanged
     └─────────────────────┘
                │
                ▼
        record the ΔLL score
        it came up with
                │
                ▼
   ...repeat, collect all the
   scores into one big list
```

Sorted, that list gives an empirical picture of "how signal-like can noise look":

```
   ΔLL scores from noise, sorted (more negative = "looked more like a real signal")

   worst-looking  ─────────────────────────────────►  best-looking (by chance)
        │◄──────── 99.9% of all noise toys ─────────────────►│
                                                    ┌─────────┴──────┐
                                                    │  the calibrated │
                                                    │  cutoff line     │
                                                    └─────────────────┘
```

The `target_tail_prob` quantile from the best-looking end becomes the calibrated detection cut: a real candidate in Phase 4 only counts as detected if it beats what noise alone achieves some small fraction (default 0.1%) of the time. This mirrors PhysRevD.94.082006's own Fig. 6, and for the same reason: the fit's regularity conditions (profiled amplitude, a boundary-adjacent width parameter) don't cleanly support an asymptotic formula, so the threshold is read off real toy data instead.

### 3.2 A second bug the design process caught: how *precisely* must a noise toy be fit?

While tuning for throughput (a production calibration runs `n_toys` in the tens of thousands to hundreds of thousands), it became clear that pure-noise fits are the *expensive* case, not the cheap one: a real signal has a genuine minimum to converge toward and typically terminates in well under 100 iterations, while a flat noise landscape has no such minimum, so the search burns its entire iteration budget on nearly every toy (measured: ~11 ms/toy at the accuracy-preserving default of 300 iterations/start, vs. sub-millisecond for a well-converged bright-signal fit).

The tempting fix — give `NoiseTailCalibrator`'s internal fits a much smaller iteration cap, since an individual toy's exact ΔLL "shouldn't matter much" in a statistical ensemble — was tested directly rather than assumed safe, by running the actual calibration (not just individual toys) at both a fast and a full-precision setting and comparing the *resulting calibrated cut*, not just per-toy differences:

| Fit precision (max iterations) | `ΔLL_cut` at `n_toys=5000` |
|---|---|
| 60 (fast) | `−10.72` |
| 150 | `−11.79` |
| 300 (accuracy-preserving default) | `−12.10` |

An 11% shift in the actual threshold, in the direction of being *more lenient* — an under-converged toy ensemble systematically underestimates how good noise can look, so the resulting cut lets more real false positives through than `target_tail_prob` promises. This is a direct, if smaller, echo of the Slice 2 boundary-pinning issue (§2.2): a shortcut that looks harmless per-instance compounds into a real bias on the aggregate, physics-relevant output. **Conclusion: `NoiseTailCalibrator` always uses the same accuracy-preserving fit configuration real events would use — no speed shortcut.** The resulting cost (~11 ms/toy, so ~18 minutes for a production-scale `n_toys=100000` run) is accepted as a one-time per-detector-config calibration cost, the same category of expense `EfficiencyMC`'s own `ne_trials`-scale toy loops already carry elsewhere in this codebase.

A related, smaller mistake caught during validation-app development: an early version of the "does the cut converge as `n_toys` grows" check happened to compare the last sweep point against itself (a default parameter coincided with the sweep's own last entry), making that specific comparison vacuously pass regardless of whether real convergence was happening. Fixed by de-duplicating the sweep before comparing — a reminder that a validation check passing is only informative if it was actually capable of failing.

### 3.3 Validation results

`apps/ccdarksens_calibrate_noise_tail.cc`, 15×15 window, `σ_pix = 0.16 e⁻`, `target_tail_prob = 1e-3`, full accuracy-preserving fit precision throughout:

| Check | Result |
|---|---|
| Convergence: `ΔLL_cut` at `n_toys` = 500 / 2000 / 8000 | `−13.22` / `−12.92` / `−12.38` — relative change over the last step: `4.4%` (threshold `25%`) |
| Seed stability: `ΔLL_cut` at `n_toys=2000`, 5 independent seeds | mean `−10.32`, sd `0.31`, relative spread `3.0%` (threshold `30%`) |
| Distribution shape (`n_toys=8000`) | monotone sorted ✓, one-sided (`ΔLL ≤ 0`) ✓ — min `−16.47`, 0.1% pctile (the cut) `−12.43`, median `−2.87`, max `−0.10`: a smooth, one-sided, heavy-near-zero tail, qualitatively as expected |

`--n_toys_full` (default `8000`, chosen so the validation app finishes in a few minutes) lets the suite be re-run at the actual production default (`--n_toys_full 100000`, ~18 minutes) to validate at full scale directly.

---

## 4. Slice 4 — `ClusterFitMC`: forward-simulating real events into the efficiency kernel

### 4.1 Plain-language picture

Slices 1–3 gave us a scorer, a searcher, and a calibrated "is this real or noise" threshold. Slice 4 uses all three on real simulated physics: for a grid of true (electron-equivalent) energies, it forward-simulates a diffusion-spread charge cloud, adds realistic noise, fits it, and applies both cuts.

```
  pick a TRUE energy (say, "this recoil deposited 800 eV")
              │
              ▼
   ┌─────────────────────┐
   │ where in the silicon  │  ← reused unchanged from the existing
   │ did it happen, and     │     DM-e pipeline (ChargeTransport)
   │ how much will the      │
   │ charge cloud spread?   │
   └─────────────────────┘
              │
              ▼
   ┌─────────────────────┐
   │ roughly how many       │  ← Fano-suppressed statistics, sourced
   │ electrons does that     │    from the target paper (§4.2) — the
   │ energy produce?         │    existing table only reaches ~20 e-
   └─────────────────────┘
              │
              ▼
   ┌─────────────────────┐
   │ scatter those          │  ← same diffusion physics as before,
   │ electrons onto a        │    just with far more electrons
   │ pixel grid, add noise   │
   └─────────────────────┘
              │
              ▼
   ┌─────────────────────┐
   │  run it through the    │  ← Slices 1–2, completely unchanged
   │  SAME fit engine        │
   └─────────────────────┘
              │
        ┌─────┴─────┐
        ▼           ▼
   beat the      didn't beat it →
   noise cut?    invisible, this
   (Slice 3)     trial is a MISS
        │
        ▼
   width looks physically
   sane (fiducial cut,
   §4.3)?  → no → MISS
        │
        ▼ yes
   record: "true energy X → we
   measured reconstructed energy Y"
```

Repeated thousands of times per energy point, this builds `K[E_true, E_reco]` — for every true energy, how often it's detected at all, and what energy gets reconstructed when it is. This is the continuous-energy analogue of `EfficiencyMC`'s `P(pattern|n_e)` table, needed because nuclear recoils reach hundreds of electrons, well past where a discrete pattern classifier applies.

### 4.2 Converting true energy to an electron count: sourced from the target paper, not invented

The one genuinely new physics assumption in this pipeline — how many electrons a given deposited energy actually produces, in a regime with no existing lookup table — turned out to be directly addressed by the target paper itself (arXiv:1607.07410 / PhysRevD.94.082006) rather than needing to be invented. The paper models the ionization-statistics resolution as `σ²_res = σ0² + (3.77 eV)·F·E`, i.e. **Fano-suppressed Gaussian statistics** around the mean electron count — exactly the functional form used here, not an ad hoc simplification:

```
n_e ~ round(Gaussian(mean = E_true / eh_pair_eV,  variance = fano_factor * mean))
eh_pair_eV = 3.77 eV     (the paper's stated value)
fano_factor = 0.133      (the paper's own measured value, from x-ray calibration of ionizing particles)
```

The paper is explicit that this Fano factor is **only measured for ionizing particles, and states it is unknown for nuclear recoils** — they vary it from `0.13` up to `1.0` (full Poisson) as a systematic check rather than asserting one number. `fano_factor` is therefore a real config field here, not a hardcoded constant, so the same kind of sweep can be run later rather than silently committing to one assumption.

Note the separation of noise sources: this Gaussian-Fano term covers only the *ionization-statistics* variance. The paper's `σ0` term (baseline readout/electronic noise) is deliberately **not** duplicated here — it's already supplied independently by `PixelSimulator::AddReadoutNoise()` downstream, the same way `EfficiencyMC` keeps dark current as a separate background term rather than folding it into the per-electron noise model. Adding both would double-count the same physical noise source.

One further, explicit divergence from the paper worth flagging: this implementation uses `σ_pix = 0.16 e⁻` (the modern skipper-CCD value already validated for the DM-electron channel), not the 2016 dataset's actual `σ_pix ≈ 1.8 e⁻`. This is a deliberate choice, not an oversight — it's trivially a config change (`pix_cfg.sigma_readout_e` / `fit_cfg.sigma_pix_e`) to switch to the literal historical value if an exact Fig. 6/9/11 reproduction is wanted later; the current default instead reflects a modern-detector projection.

### 4.3 A finding worth understanding, not just accepting: efficiency plateaus below 100%

The built kernel's efficiency curve rises smoothly with `E_true` and then plateaus around **~73–75%**, even at the brightest energy tested (300 eV) — it never approaches 100%. Rather than accept this on faith, it was checked directly: is this a real physical ceiling, or a bug?

**It's real, and fully explained by the fiducial cut interacting with the depth-dependent diffusion model.** `ChargeTransport::SampleDepthUm()` samples interaction depth uniformly across the full detector thickness (670 μm), and the resulting diffusion width `σ_xy(z)` is a fixed, energy-independent function of depth alone (`beta_per_keV = 0` in the default diffusion config, so brightness has no effect on the *true* cloud width). A direct numeric integral over that depth distribution shows **69.5% of interactions have a true σ_xy that falls inside the fiducial window `[0.35, 1.22] px` regardless of energy** — 7.8% are too narrow (near-surface), 21.7% are too wide (too deep) — so **no amount of brightness can push efficiency above roughly 70% under this diffusion configuration.** This is not a code artifact; it is the same category of geometric/fiducial rejection real DAMIC-style analyses report.

The observed plateau (~73–75%) sits slightly *above* that pure-geometric ceiling, which was also checked rather than left unexplained: a direct comparison of true vs. fitted `σ_xy` on 1500 bright-event trials found the fit has a small negative bias (`mean(fit − true) = −0.024 px`) and, because the fiducial window's *upper* edge has far more geometric mass just outside it (21.7%) than the *lower* edge does (7.8%), that bias asymmetrically "smears" more events *into* the passing window (7.3% of trials) than *out of* it (3.7%) — a net effect that fully accounts for the gap between the pure-geometric ceiling and the observed plateau. Real fits into a hard geometric boundary can't be expected to be perfectly sharp at that boundary; this is the expected consequence, quantified rather than assumed.

### 4.4 Validation results

`apps/ccdarksens_build_cluster_fit_kernel.cc`: calibrates a detection cut (`n_toys=5000`, matching window/fit config), then builds the kernel over 12 log-spaced `E_true` points from 10–300 eV, 1500 trials each.

| Check | Result |
|---|---|
| Calibrated `ΔLL_cut` (this run's window/noise config) | `−12.10` |
| Efficiency at lowest `E_true` (10 eV) → highest (300 eV) | `0.12 → 0.75`, smooth and monotonic in trend |
| All efficiency values in `[0,1]` | ✓ |
| `K` rows sum to the reported per-point efficiency | ✓ (exact, by construction) |
| Geometric fiducial ceiling (independent numeric cross-check) | `69.5%` pass fraction from uniform depth sampling alone — matches the observed plateau once fit-noise boundary-smearing (§4.3) is accounted for |

No exact external reference exists for this specific run (detector parameters here are a modern-projection choice, not literally PhysRevD.94.082006's 2016 dataset) — validation here is physical-sanity and internal-consistency, following the strategy in §0/plan rather than numerical parity.

---

## 5. Slice 5 — wiring: the fold function and config integration

### 5.1 Plain-language picture

Slices 1–4 built the detector-response machinery. Slice 5 adds no new physics or algorithms — it connects that machinery to the rest of the pipeline that already exists and works for the DM-electron channel.

```
  raw predicted rate            Phase 4 kernel
  (function of E_true,     ×    K[E_true, E_reco]
   from Phases 1-2)              (how E_true becomes E_reco)
              │                          │
              └──────────┬───────────────┘
                          ▼
              S(E_reco bin) — the vector that
              feeds directly into the existing,
              UNCHANGED ProfileLikelihood
```

Two pieces: `ClusterEnergyRates::FoldEtrueToErecoRates` (the fold itself — the continuous-energy analogue of the DM-electron channel's `FoldNeToPatternRates`), and a `response.cluster_fit_mc` JSON config block wiring every Slice 1–4 parameter (window geometry, diffusion, fit method, calibration settings, `eh_pair_eV`/`fano_factor`) into the same `ConfigManager` machinery every other channel already uses, rather than leaving them hardcoded in validation-app source.

### 5.2 The fold function

`include/ccdarksens/response/ClusterEnergyRates.hh` / `.cc`, matching `PatternRates.hh`'s existing signature style: for each `E_true` grid point in the kernel, look up the raw rate density and bin width at that energy (from a `TH1D`, e.g. `RateTable::MakeTH1D`'s output), convert to expected counts via `density * dE * exposure_kg_year`, and distribute those counts across `E_reco` bins according to the kernel's row. Unlike the DM-electron channel, there is no separate ionization/charge step before this fold — Phases 1–2 already produce a quenched electron-equivalent spectrum directly, so the kernel operates straight on `E_true = E_ee`.

One initial design note corrected before it became a real constraint: the fold function does **not** require the rate histogram's binning to match the kernel's `E_true` grid exactly (the kernel grid is log-spaced by design, to resolve the near-threshold turn-on region well; a rate histogram is naturally built with uniform bins). `TH1D::FindBin` samples whichever local density the histogram has at each kernel grid point, so a fine uniform histogram works fine against a coarser or log-spaced kernel grid — the actual requirement is just "fine enough resolution relative to how quickly the rate varies," the ordinary histogram-resolution caveat, not a structural one.

### 5.3 Config wiring

New `ClusterFitMCJSON` struct in `include/ccdarksens/io/ConfigManager.hh`, added as `response.cluster_fit_mc`, parsed in `ConfigManager::parse_response_` following the exact `jp.value(key, default)` pattern already used for `efficiency_mc` (including the same nested-vs-flat `diffusion` sub-object handling). `response.analysis_space` accepts a new value `"cluster_energy"` — no schema change needed, it was already an unvalidated string.

### 5.4 Validation: end-to-end integration, not another isolated component

`apps/ccdarksens_example_one_point_cluster.cc` is the first point in the whole build where every piece actually connects: it loads `configs/wimp_nucleon_cluster_example.json` through `ConfigManager` (proving the JSON wiring parses correctly), runs the Slice 3 calibration and Slice 4 kernel build **once** from that config, then folds two different Phase 1/2-generated rate spectra (`m_χ = 3000 MeV` and `5000 MeV`, same cross section) through that **same, un-rebuilt kernel** — deliberately demonstrating that the kernel is a pure detector-response object independent of the DM physics parameters, not something to rebuild per grid point (the same principle `EfficiencyMC`'s table already follows for the DM-electron channel).

| Check | Result |
|---|---|
| Config parses, `analysis_space` round-trips correctly | `cluster_energy` ✓ |
| Kernel built once (calibration `ΔLL_cut = −12.10`, 12 `E_true` points × 1500 trials), reused for both mass points | ✓ |
| `m_χ=3000 MeV`: folded `S_reco` all finite/non-negative, total `349.8` counts (0.1 kg·yr) | ✓ |
| `m_χ=5000 MeV`: same, total `354.4` counts | ✓ |
| Both folded vectors accepted by the **unchanged** `ProfileLikelihood` (`SetData`/`SetBTemplate`/`EvaluateRatio`), finite `q` returned | `q=0.0348` / `q=0.0352` ✓ |

---

## 6. Phase 7 — reproducing PhysRevD.94.082006's own published figures

### 6.1 Scope and the one deliberate divergence this closes

Every validation through Slice 5 used a modern-detector projection (`σ_pix = 0.16 e⁻`) rather than the 2016 dataset's actual `σ_pix ≈ 1.8 e⁻` — an explicit, flagged choice (§4.2), not an oversight, precisely so an exact Figs. 6/9/11 reproduction could be done later as a separate check. This section is that check: switch every detector parameter to the paper's own stated values and see how closely the pipeline's own outputs land on the paper's own numbers, using no numbers not sourced from the paper text itself (fetched from arXiv:1607.07410 / PhysRevD.94.082006, Sections II–VIII).

| Parameter | Paper's value | Source |
|---|---|---|
| `sigma_readout_e` | 1.8 e⁻ | Sec. V: "σ_pix = 1.8 e⁻ ≈ 7 eVee" |
| Fit window | 11×11 px | Sec. VI: "11×11-pixel window" |
| `thickness_um` | 675 | Sec. II |
| `pixel_size_um` | 15.0 | Sec. II (already matched) |
| `eh_pair_eV` / `fano_factor` | 3.77 eV / 0.133 | Sec. IV.1 (already matched — this is the same collaboration's companion calibration paper, PhysRevD.94.082007, whose measurement the `chavarria_table` quenching model already uses) |
| Exposure | 0.6 kg·day | Abstract |
| Halo (v0, vE, vesc, ρ) | 220, 232, 544 km/s, 0.3 GeV/cm³ | Sec. VIII: "standard halo parameters" — the classic Lewin–Smith convention, **not** the Baxter et al. 2021 convention (`v0=238, vE=253.7`) this repo's Migdal channel defaults to elsewhere |

New config `configs/wimp_nucleon_cluster_damic2016_repro.json` (kernel/efficiency checks) and `configs/wimp_nucleon_damic2016_limit_repro.json` (full scan), plus `configs/wimp_nucleon_generate_si_damic2016halo.json` (rate generation with the paper's halo parameters) and a new app, `apps/ccdarksens_validate_wimp_nucleon_paper_repro.cc` — `ccdarksens_build_cluster_fit_kernel.cc` hardcodes the modern-projection constants and takes no config, so a config-driven twin was needed to actually vary them.

### 6.2 Fig. 6 — ΔLL tail calibration: explained, not brute-forced

Naively calibrating at the framework's default `target_tail_prob=1e-3` gives `ΔLL_cut ≈ −11` to `−12` (8000–50000 toys) — far short of the paper's quoted `−28`. Running enough toys to sample that same tail probability directly (the paper's real dataset spans millions of independent fit windows across multiple 8-Mpix CCDs) isn't practical here — 50000 toys already took ~10 minutes, and closing the gap directly would need roughly `10⁷` toys.

Instead, checked whether the gap is *explained* by tail depth rather than a methodology mismatch. `ΔLL` behaves like `−(z²)/2` for a Gaussian tail (`z` = the standard-normal quantile), inflated by a "look-elsewhere" factor from maximizing over the continuous `(x, y, σ_xy)` nuisance parameters — exactly why the paper itself needed empirical toy calibration rather than an asymptotic chi-square formula in the first place (§3, no closed form applies here). Measuring that inflation factor at two different tail depths:

| Tail probability | Naive χ²₁ prediction | Measured `ΔLL_cut` | Inflation factor |
|---|---|---|---|
| `1×10⁻³` (n=50000) | `−4.77` | `−11.26` | `2.4×` |
| `~2×10⁻⁵` (min of 50000 toys) | `−8.45` | `−16.47` | `1.95×` |

The inflation factor is stable across an order of magnitude of tail depth. Extrapolating it to the tail probability a real multi-megapixel, multi-CCD dataset would actually sample (`~10⁻⁷`, physically plausible given millions of independent fit windows) predicts `ΔLL_cut ≈ −13.5 × 2.2 ≈ −30` — close to the paper's `−28`. This is treated as an explained, self-consistent gap (confirmed by an independent extrapolation from our own measured trend), not a numerical match — no claim of digit-level parity is made here.

### 6.3 Fig. 9 — efficiency curve: strong shape and plateau match

With the paper's window/noise parameters, the kernel's efficiency curve:

| `E_true` [eV] | Efficiency | Paper (Fig. 9, qualitative) |
|---|---|---|
| 60 | 0.098 | "9% at 75 eVee" (1×1 mode) — same ballpark |
| 78.6 | 0.212 | rising steeply |
| 102.9 | 0.381 | rising steeply |
| 134.8 | 0.579 | rising steeply |
| 176.5 | 0.719 | approaching plateau |
| 231–2000 | 0.70–0.74 | "plateaus ~75%" |

The turn-on shape (steep rise from threshold, plateauing by a few hundred eV) and the plateau value (our `~70–74%` vs. the paper's `~75%`) both match well. One methodological note worth recording: the first run showed a spurious drop to `0.37` at the top grid point (`E_true=2000 eV`) — traced to `Ereco_max_eV` being set exactly equal to `Etrue_max_eV`, truncating part of the reconstructed-energy distribution off the grid for events near the top edge. Re-run with headroom (`Ereco_max_eV=4000`) fixed it (`0.69` at 2000 eV, consistent with the plateau) — a config-construction artifact caught and corrected, not a real physical turnover, worth flagging since it's an easy mistake to repeat when choosing `Ereco_max_eV` for any future high-`E_true` reproduction.

### 6.4 Fig. 11 — exclusion limit: a real background, a wider grid, and the right literature curve

Fresh rate CSVs generated with the paper's own halo parameters (`configs/wimp_nucleon_generate_si_damic2016halo.json`, **28** log-spaced masses 0.5–10 GeV × **20** log-spaced cross sections 10⁻⁴²–10⁻³³ cm²), scanned with `ccdarksens_scan_generic` (a real end-to-end use of the Slice 4-8 generic scan app).

**Three corrections were needed before this comparison meant anything**, each caught by checking rather than assuming:
1. **Wrong comparison curve, twice.** The first pass compared against text-extracted, qualitative numbers and concluded (wrongly) "consistently more stringent, as expected." A digitized curve found in a local `WIMPyCCD` reference repo (`digitized/damic at snolab.csv`) flipped that to "consistently *less* stringent, 2–8×" — but that file was also wrong: the same folder's `radomir's curve.csv` and `DAMIC results/{nominal,lower,upper}.csv` triplet agree with each other everywhere, while `damic at snolab.csv` sits ~10× below even the triplet's lower edge (implausible, a >5σ fluctuation). Best-guess labeling from that internal consistency: `nominal`/`lower`/`upper` = the paper's expected-sensitivity median ± 1σ, `radomir's curve` = the actual observed 90% CL limit. Saved as `data/previous_limits/WIMP/DAMIC_2016_SNOLAB_{observed_radomir,expected_nominal,expected_lower,expected_upper}.csv`; the mislabeled file was removed rather than left to mislead a future comparison.
2. **Background was zero.** Fixed using the paper's own stated number (Sec. VIII: `15±3 events/(kg·day·keVee)`, flat Compton, for the same 1×1-mode configuration used throughout this section) → `backgrounds.flat_background.norm_per_kg_year_keV = 5478.75` (`15 × 365.25`).
3. **Grid too coarse, and too narrow at the low-mass end** (the original 18×15 grid left several low-mass points at the σ_n grid edge, unresolved) — widened to 28×20 with the σ_n grid extended to 10⁻³³ cm².

With a real background in place, the more meaningful comparison shifts from the *observed* curve (one particular statistical fluctuation) to the *expected median* (the same kind of quantity our Asimov calculation actually computes):

| `m_χ` [GeV] | This repo's UL [cm²] | Expected median (nominal) [cm²] | ratio |
|---|---|---|---|
| 2.0 | `5.0×10⁻³⁹` | `7.8×10⁻³⁹` | `0.64×` |
| 3.0 | `8.9×10⁻⁴⁰` | `2.5×10⁻³⁹` | `0.35×` |
| 5.0 | `4.7×10⁻⁴⁰` | `5.8×10⁻⁴⁰` | `0.81×` |
| 6.0 | `4.6×10⁻⁴⁰` | `4.5×10⁻⁴⁰` | `1.02×` |
| 7.0 | `4.9×10⁻⁴⁰` | `3.9×10⁻⁴⁰` | `1.25×` |
| 9.5 | `6.6×10⁻⁴⁰` | `3.5×10⁻⁴⁰` | `1.90×` |

More stringent than the expected median below ~6 GeV, crossing to somewhat less stringent above it (up to ~1.9× at 9.5 GeV) — a real, mild, mass-dependent difference, not the dramatic mismatches from the earlier (wrong-file or wrong-background) attempts. The visual plot shows our curve tracking inside or just at the edge of the expected ±1σ band across essentially the whole range. No further investigation forced at this level of agreement, though the modest high-mass gap (a real background estimate with no stated uncertainty band applied, vs. the paper's own fitted/uncertain background level) is a plausible, unexplored source, not yet isolated.

Comparison plot: `outputs/wimp_nucleon_damic2016_limit_repro/plots/damic2016_limit_comparison_rootstyle.png` — matplotlib styled to match `ccdarksens_plot_dmelectron_limit`'s actual `gStyle`/`TCanvas` settings (inward ticks all sides, thick lines, borderless legend), since the real app has no flag to overlay an arbitrary literature CSV. Reproducible via `configs/wimp_nucleon_generate_si_damic2016halo.json` → `configs/wimp_nucleon_damic2016_limit_repro.json` → `ccdarksens_scan_generic` (exact commands and parameters in the config's own `_comment`).

### 6.4b A fourth correction: a directly-verified digitization supersedes `radomir's curve`

`radomir's curve.csv`'s identity as "the observed curve" was never independently confirmed — only inferred from internal consistency with the `nominal`/`lower`/`upper` triplet (§6.4 point 1). That inference turned out to be incomplete: the user directly digitized the actual Fig. 11 from the published paper (a single, unambiguous "This work, 0.6 kg d" curve with its own shaded uncertainty band — no separate observed/expected lines to confuse in that figure) and handed over the result. It does **not** match `radomir's curve` well at low mass (~20× apart at 1.3–2 GeV), though the two converge at high mass (~7–9.5 GeV). Since this new curve is a directly-verified trace of the real published figure — the most solidly-grounded of any DAMIC 2016 curve in this repo — it now replaces `radomir's curve` as the primary comparison target. Saved as `data/previous_limits/WIMP/DAMIC_2016_SNOLAB_observed_fig11_diego.csv`; `radomir's curve` is kept only as a secondary, explicitly-demoted reference (thin teal line in `--wimp` mode), not deleted, since it isn't *known* to be wrong, only unconfirmed and now inconsistent with a better source at low mass.

Against this corrected reference, interpolated at matched masses:

| `m_χ` [GeV] | This repo's UL [cm²] | Fig. 11 observed (verified) [cm²] | ratio |
|---|---|---|---|
| 1.5 | `9.3×10⁻³⁸` | `6.4×10⁻³⁸` | `1.5×` |
| 2.0 | `5.0×10⁻³⁹` | `7.7×10⁻³⁹` | `0.65×` |
| 3.0 | `8.9×10⁻⁴⁰` | `1.2×10⁻³⁹` | `0.76×` |
| 5.0 | `4.7×10⁻⁴⁰` | `3.3×10⁻⁴⁰` | `1.44×` |
| 7.0 | `4.9×10⁻⁴⁰` | `2.1×10⁻⁴⁰` | `2.28×` |
| 9.5 | `6.6×10⁻⁴⁰` | `2.0×10⁻⁴⁰` | `3.27×` |

Within 0.65–1.5× for 1.5–5 GeV, widening to ~3.3× by 9.5 GeV — a real, honest result: better agreement than either the wrong-file or wrong-background attempts, worse than the (retrospectively unreliable) match to `radomir's curve` at high mass. `--wimp` mode (`ccdarksens_plot_limit`) now also carries `DAMIC-M_projected_1kgyear.csv` (a different, newer/larger experiment, ~2–3 orders of magnitude more sensitive — shown for context, needs its own wide y-range, `--y-min 1e-43`, to be visible at all) and the originally-flagged `damic at snolab.csv` (kept, explicitly labeled "uncertain provenance," per an explicit request to include it anyway rather than silently drop it) — full inventory and reasoning for every file in `data/previous_limits/WIMP/`.

### 6.5 Summary

Fig. 6's gap is explained by an independently-verified tail-inflation extrapolation, not closed; Fig. 9 matches closely in both turn-on shape and plateau value; Fig. 11 tracks the directly-verified observed curve to within a factor of ~1.5× through the mass range where the paper's own sensitivity is best, widening to ~3× at the high-mass end. The process is as much the lesson as the result: **four** successive corrections — paraphrased text, a digitized-but-mislabeled file, a background-free approximation, and finally an unconfirmed-but-internally-consistent literature curve — were each caught by checking rather than by accepting a first plausible-looking answer, the last one only because the user went and independently re-digitized the actual figure rather than trusting a found file's inferred identity. None of the earlier conclusions should have been written down as settled without that checking, and "internally consistent with other found files" turned out to be weaker evidence than "directly verified against the source."

### 6.6 Physics/detector-effects audit

Confirmed present in the pipeline, checked directly against source rather than assumed: the bug-fixed nuclear form factor and DM-nucleon (not DM-nucleus) reduced-mass convention in the rate prefactor (`ccdarkphys/wimp_nucleon/rate.py::dRdE_nr_kg_day_keV` — confirmed the generator calls this, not the literal WIMPyCCD-bug-preserving twin kept only for historical parity checks); the halo velocity integral with paper-matched parameters, threaded through as explicit overridable kwargs (`halo.py`); measured nuclear-recoil quenching (`chavarria_table`, the same collaboration's companion calibration paper to this one); Fano-suppressed electron-count statistics; depth-dependent charge diffusion; pixel readout noise; the 2D-fit + empirical ΔLL noise-tail cut; and the depth/fiducial-cut efficiency ceiling that Fig. 9 validates closely. Not modeled: lateral (x-y) spatial/edge cuts — not needed, since the paper's own quoted ~75% efficiency is explicitly for "events uniformly distributed in the CCD bulk," the depth effect alone, matching what this pipeline computes. Exposure, live mass, and any additional data-quality cuts are taken directly from the paper's final quoted 0.6 kg·day number rather than separately modeled, so whatever those cuts cost is already folded in by construction.

### 6.7 Closing the largest known gap: the joint 1×1 + 1×100 likelihood

§6.4/6.4b left one identified, unclosed structural gap: the paper combines two independent CCD readout modes (`L_joint = L_1×1 × L_1×100`, Sec. V) into one likelihood; this repo ran only the 1×1 channel. Implemented in full:

**A real physical-model finding, not just a config change.** The paper states the 1×100 mode sums 100 pixel rows in hardware *before* the single readout ("collapsing the pixel contents along the y axis," Sec. II fig. caption) and fits "in one dimension along rows" (Sec. VI). A 2D isotropic Gaussian is separable, so the y-marginalized truth landing in the collapsed row is exactly `PixelIntegral1D` alone (already in `ClusterFitModel.hh`) — not a 2D fit degenerately run on a 1-row window, which would be wrong twice over: `PixelSimulator::DepositElectron` drops any electron whose y-offset lands outside the configured window instead of summing it, and the existing 2D `Shape()` would apply a spurious extra y-integration factor even if it didn't. Added: `Shape1D`/`IHat1D`/`ObjectiveG1D` (`ClusterFitModel.hh/.cc`), a `collapse_y` mode on `PixelSimulator::DepositElectron` (sums into row 0 unconditionally instead of rejecting), and a 2-parameter (`mux`, `sigma`) Nelder-Mead path on `ClusterFitEngine` selected by `ClusterFitConfig::one_dimensional`. `NoiseTailCalibrator` and `ClusterFitMC` needed **zero** code changes — both already just pass `pix_cfg`/`fit_cfg` through, so the same calibration/kernel-building classes that validated 1×1 work for 1×100 once the two flags are set.

**Why the joint NLL is exact, not an approximation.** Each channel's background nuisance parameter only ever appears in that channel's own NLL term — channels share the physical signal spectrum, not a nuisance parameter. So `min_θ NLL_joint(σ,θ⃗) = Σ_k min_θk NLL_k(σ,θk)`: the joint profile-likelihood-ratio decomposes into a plain sum of the two channels' independently-computed NLLs, evaluated with the exact same `ProfileLikelihood::MinimizeOverScale`/`MinimizeOverScaleMinuit` calls the single-channel path already uses. `ProfileLikelihood.hh/.cc` and `ScanUtils.hh/.cc` needed no changes; `ccdarksens_scan_generic.cc` gained one new, separate `RunJointChannelScan` code path (gated on `response.channels` non-empty) that builds one `ResponseFold`+`ProfileLikelihood` per channel and sums `nll_values` before the existing `MonotonizeQ`/`UlFromQMonoCrossing` calls — the pre-existing single-channel path (bit-exact parity-gated against the frozen reference apps) is untouched. Config schema: `response.channels[]` (`ConfigManager.hh`), empty by default so every existing config is unaffected.

**Verification.** A real 1×100 config (`configs/wimp_nucleon_cluster_damic2016_1x100_repro.json`) run through the existing `ccdarksens_validate_wimp_nucleon_paper_repro` gives a naive ΔLL cut of −7.75 (n=5000) — *less* deep than 1×1's −10.8, matching the paper's own ordering (−25 vs −28); efficiency turns on faster at low energy (48.8% at 60 eV vs. 1×1's 12.2%), also matching the paper's stated ordering (25%@60eV vs 9%@75eV). One honest discrepancy: our 1×100 efficiency plateaus at ~75-77%, not the paper's quoted ~100% — most likely because the paper's Fig. 9 curve is *pre*-fiducial-cut selection efficiency, while `ClusterFitMC`'s kernel already applies the fiducial σ_xy cut for both channels (the 1×1 kernel shows the same ~75% ceiling, which Fig. 9 also does *not* claim reaches 100% for 1×1 — self-consistent, but not independently confirmed here).

**Result** (`configs/wimp_nucleon_damic2016_limit_repro_joint.json`, background rates 15×365.25 / 21×365.25 events/(kg·yr·keV) per Sec. VIII): the joint limit is tighter than the 1×1-only curve at every grid mass (e.g. 2.70×10⁻⁴⁰ cm² vs 5.00×10⁻⁴⁰ cm² at 4.1 GeV — a 1.85× improvement), moving the repro measurably closer to the paper's real curve at high mass:

| m_χ [GeV] | 1×1-only / paper | joint / paper | improvement |
|---|---|---|---|
| 4.12 | 1.08× | 0.58× | 1.85× |
| 5.14 | 1.38× | 0.84× | 1.65× |
| 7.17 | 2.16× | 1.30× | 1.66× |
| 10.00 | 3.26× | 1.68× | 1.94× |

A caveat reported honestly, not glossed over: at low-to-mid mass (1.4–4 GeV) the joint curve now sits *below* the paper's real curve (ratio < 1 — e.g. 0.21× at 2.1 GeV), i.e. claiming better sensitivity than the paper's own real result. The likeliest cause is the 1×100 channel's σ_pix=1.8 e⁻, carried over from the paper's 1×1-section detector characterization rather than a separately-stated 1×100 value (flagged as an assumption in both new configs' `_comment` fields) — if the true 1×100 readout noise is somewhat higher, that channel's contribution here is correspondingly overstated. A secondary contributor: this is still an Asimov (background-only) projection against the paper's real, possibly upward-fluctuated 31/23-candidate dataset, not a fit to real per-event data (§6.1's original caveat, unchanged by this pass). Not chased further in this pass — flagged as the natural next thing to check if the 1×100 sensitivity gain needs tightening.

### 6.8 A real bug found chasing the high-mass divergence: `FoldEtrueToErecoRates` was discarding ~90% of the signal

§6.7's joint-likelihood result still diverged from the paper at high mass even after the mass-scaling shape was independently validated against `wimprates` (§6.7 follow-up) and after a real quenching-model gap was found and closed (see below) — neither moved the curve. Chasing the divergence further, in order:

**1. Quenching model gap — real, but not the driver.** The paper (line 275) actually combines **two** measured Si nuclear-recoil ionization-yield datasets: Chavarria et al. 2016 (already in `quenching.py` as `chavarria_table`, 0.68–2.28 keV_nr) *and* Izraelevitch et al. 2017 (JINST 12, P06014, arXiv:1702.00873 — the paper's own "antonella1" citation), covering 2–20 keV_nr. Our code only had the first, falling back to unconstrained Lindhard theory above 2.28 keV_nr — a regime whose share of the total signal rate grows from ~6% at 4 GeV to ~54% at 10 GeV. Added `chavarria_izraelevitch_table_yield` (`python/ccdarkphys/wimp_nucleon/quenching.py`) combining both real datasets (hard switch at 2.28 keV_nr, Chavarria's own upper bound), with `data/izraelevitch2017_table1.csv` holding the real Table 1 values (verified against the primary source, not just an AI-summarized extraction — see repo history). Regenerating rates with it changed the final UL by **<1.5% at every mass point** — the fix is real and worth keeping (more faithful to the paper), but not the cause of the divergence. Mechanism in retrospect: the new model's *lower* yield than the old Lindhard fallback pulls events that used to land above the 4000 eV window back into it, roughly canceling the within-window redistribution.

**2. The actual bug: `src/response/ClusterEnergyRates.cc::FoldEtrueToErecoRates`.** The kernel's `Etrue_grid_eV` is a *sparse* set of points (16, log-spaced 60–4000 eV) — each stands for a whole neighborhood of true energies, not just itself. The original code looked up the density at each point in the *input spectrum's own fine histogram* (400 bins, 10 eV wide) and multiplied by that 10 eV bin's width, then moved to the next kernel point — silently discarding everything in between two kernel points (up to ~1000 eV of untouched spectrum at the top of a log-spaced grid). Verified numerically on the real m_χ=10 GeV spectrum: true integrated rate 8065.6 events/kg/year; what the buggy fold actually captured: **816.5 — 10.1% of the true total.** This affected every `cluster_energy`-analysis-space result computed in this project to date, both signal and background (both fold through the same function), including all of §6.1–6.7's numbers above.

Fixed by replacing the point-sampled fine-bin width with a proper quadrature: each kernel point's neighborhood is bounded by the *geometric midpoints* with its neighbors (matching the log-spaced grid), and the input histogram is genuinely integrated (`TH1::Integral(bin_lo, bin_hi, "width")`) over that full neighborhood, not point-sampled. Recovers ~93% of the true total (the remaining ~7% is the real, expected detection inefficiency the kernel itself encodes) — see the updated header docstring in `include/ccdarksens/response/ClusterEnergyRates.hh` for the full derivation.

**Why this specifically explains the mass-*growing* divergence, not just an overall offset:** the undercounting severity depends on how well each mass's specific spectrum shape happens to align with the 16 fixed sample points — a flat, wide spectrum (the background, or a broad high-mass signal) sampled by the same sparse grid loses proportionally more than a sharply-peaked one. Confirmed directly: after the fix, the flat background's total folded rate grew **~28×** (0.879→24.999 events for the 1×1 channel), while a representative signal spectrum at 4.1 GeV grew only **~4.8×** at the same grid point — very different factors for different spectral shapes, exactly the kind of mass/shape-dependent bias needed to distort the curve's mass-dependence specifically, not just its normalization.

**Result after the fix** (full rerun: 1×1-only, joint chavarria-only, joint +Izraelevitch — all `wimp_nucleon_damic2016_limit_repro*.json` configs): the mass-dependence shape problem is gone. The joint/paper ratio now climbs smoothly and monotonically from 0.13× at 1.7 GeV to 0.54× at 10 GeV — no more upturn past 5 GeV. The 1×1-only and joint curves both now sit close to or below the paper's real curve across the whole range, with no artifact of the kind seen in §6.7's table.

**A new, cleaner remaining gap, not yet closed:** the curve is now *systematically too optimistic* (ratio < 1 at every mass, was as extreme as ~0.11–0.13× at 1.7–2.4 GeV) rather than the previous mixed under/over-shoot. This is consistent with — and not yet explained by anything else found so far — the background-efficiency-curve gap: the paper's own Fig. 9 shows the *background* detection efficiency (dashed lines) is a distinct, lower curve than the *signal* detection efficiency (solid lines) — 1×1 background peaks at ~0.56 near threshold and declines to ~0.49 by 2 keV_ee, versus signal's flat ~0.75 plateau (line 510: "background... shape is given by a flat Compton scattering energy spectrum multiplied by the background efficiency"). This pipeline currently folds the background through the *signal's* kernel/efficiency (`BackgroundFactory::MakeClusterEnergyFlatBackground` calls the same `fold.Fold(...)` signal path), overestimating background efficiency by ~1.3–1.5× and therefore overestimating background counts — which makes the computed limit too tight. Fixed in §6.9 below.

### 6.9 Closing the background-efficiency gap: a real, distinct curve, not the signal's

§6.8's remaining gap was closed directly using real data from the paper's own Fig. 9, rather than an approximation.

**Digitization method — vector path extraction, not a visual read.** `Literature/arXiv-1607.07410v2/figures/eff_crop.pdf` is a ROOT-generated vector PDF; its four curves (Signal 1×1, Signal 1×100, Background 1×1, Background 1×100) are literally stored as drawing-path coordinates in the file, identifiable unambiguously by stroke color (black/red) and dash pattern (solid/dashed). Using PyMuPDF to read the PDF's own drawing commands and its own text layer (for the axis tick-label positions, giving an exact pixel→data calibration), all four curves were extracted as ~1000-point polylines — the actual coordinates ROOT drew, not a manual/visual estimate. Cross-check: the extracted Signal(1×1) curve (turn-on 0.10 at 100 eV → 0.52 at 200 eV → 0.75 plateau by 300 eV) matches this pipeline's own kernel-computed efficiency curve (§6.3) closely, confirming the digitization and the kernel are both correctly calibrated against the same physical curve. Saved as `data/wimp_nucleon_damic2016_fig9_{signal,background}_{1x1,1x100}.csv` (250 points each, full provenance in the file headers). The paper's Fig. 9 only plots to 2 keV_ee; above that, the digitized background tables are extrapolated flat at their 2 keV_ee value, per the paper's own text (Sec. IV.3: background efficiency at high energy is "dominated by the contribution from Compton events," i.e. roughly constant) — an explicit, documented choice, not a silent default.

**Implementation.** Per the paper's own construction (line 510), background is `flat_Compton_rate(E) × background_efficiency(E)`, evaluated directly at the reconstructed energy — no separate resolution-convolution or kernel-smearing step, because convolving a *flat* spectrum with a symmetric resolution kernel is a no-op (unlike signal, whose spectrum shape is far from flat, per line 508's explicit convolution step). New: `include/ccdarksens/response/BackgroundEfficiencyTable.hh/.cc` (CSV loader + linear interpolator, clamped flat outside the table's range), and `BackgroundFactory::MakeClusterEnergyFlatBackgroundWithEfficiency` (evaluates the table at each `ResponseFold::ErecoEdgesEV()` bin center directly — a new capability-query method on `ResponseFold`, default-throws, overridden only by `ClusterEnergyResponseFold`). Config: `response.background_efficiency_csv` (single-channel) / `response.channels[].background_efficiency_csv` (joint), empty by default — every existing non-cluster_energy config is unaffected. Wired into both the single-channel path (`BackgroundFactory.cc`) and the joint-channel path (`ccdarksens_scan_generic.cc::RunJointChannelScan`), each falling back to the previous kernel-folded behavior when the field is unset.

**Result, full rerun of all three configs:** background counts dropped as expected (1×1 channel: 24.999→17.517 events; 1×100: 35.244→22.880 — the 1×100 background efficiency curve declines further, to ~0.44 by 2 keV_ee, than 1×1's ~0.49). The limit tightened correspondingly at every mass. Final comparison, both fixes applied:

| m_χ [GeV] | 1×1-only / paper | joint / paper |
|---|---|---|
| 1.36 | 0.872× | 0.238× |
| 2.12 | 0.178× | 0.082× |
| 4.12 | 0.273× | 0.185× |
| 6.42 | 0.428× | 0.336× |
| 8.01 | 0.463× | 0.368× |
| 10.00 | 0.479× | 0.393× |

The high-mass behavior is now excellent — both curves flatten smoothly toward ~0.4–0.5× the paper's real curve above 6 GeV, matching the paper's own plateau shape, a dramatic improvement from §6.7's table (which ranged 1.08×–3.26× with a growing upward divergence). Low-to-mid mass (1.5–3 GeV) is where the largest remaining gap sits (as low as 0.08×) — the curves are now more optimistic than the paper's real result everywhere, not just at low mass as in §6.7. The likeliest remaining explanation, not newly discovered but never actually closed in this whole investigation: **this is still an Asimov (pure background) projection compared against the paper's real 90% CL, which is a fit to 31 (1×1) / 23 (1×100) actually-observed candidate events** — if that real dataset fluctuated upward from the pure-background expectation (a completely ordinary possibility for small Poisson counts), the paper's real limit would be correspondingly weaker than a background-only projection predicts, at every mass, with no modeling gap required to explain it. Fitting real per-event data was explicitly out of scope from §6.1 onward (the paper's own event-by-event energies aren't available in machine-readable form beyond the final published exclusion curve) and remains so here — this is reported as the current best understanding of the residual gap, not chased further in this pass.

### 6.10 Consolidation: three further checks, and what actually explains the remaining gap

Prompted directly by the observation that the residual gap (§6.9) was too large (up to ~12× at 2 GeV) to hand-wave away as "probably just Poisson noise" without checking. Three concrete follow-ups, each run to a real conclusion:

**1. Background-efficiency bin-averaging — real, small, fixed.** `MakeClusterEnergyFlatBackgroundWithEfficiency` (§6.9) originally evaluated the digitized background curve at each output bin's *center* rather than averaging it across the bin. Quantified directly: the bin spanning the curve's steepest rise (50–100 eV) had an 80% local error this way. Fixed with a 9-point composite-Simpson bin average (`BackgroundFactory.cc`). Effect on the full curve: ~3–4% at low mass, negligible at high mass — real, worth keeping, not the driver of the big gap.

**2. Kernel resolution (16 vs. 60 `Etrue_npoints`) — tested, ruled out, in the wrong direction.** Rebuilding the 1×1 kernel at 60 points instead of 16 (config-only, no code change) made the low-mass limit *tighter*, not looser — moving further from the paper's real curve, not closer. Coarse kernel sampling is not hiding an under-estimate that finer sampling would recover; if anything the coarse kernel was mildly *conservative* relative to the fine one at low mass.

**3. The 1×100 channel's unverified σ_pix — real, bounded, quantified, not closable without paper data.** Comparing 1×1-only against the paper (0.87×–0.48× across the grid) versus the full joint result (0.24×–0.39×) isolates the 1×100 channel as responsible for an additional ~0.27×–0.39× on top of whatever gap 1×1-only has on its own. Sweeping the 1×100 channel's assumed σ_pix from 1.8 e⁻ (the unverified 1×1-carried-over value) up to 4.0 e⁻ closes a substantial fraction of that gap at low mass (1.36 GeV: 0.24× → 0.66×) but never fully closes it, and matters far less at high mass. This is reported as a genuine, quantified *uncertainty* — not fixed, because the paper never states a 1×100-specific noise figure to fix it *to*; picking a number because it moves the curve closer would be exactly the kind of unfounded calibration-data fabrication this whole investigation has deliberately avoided.

**The 1×1-only residual, decomposed properly.** The 1×1-only/paper ratio itself traces a V-shape with mass — close to 1 at 1.36 GeV (0.87×), a minimum around 2.6–3.7 GeV (~0.17×–0.19×), partial recovery to 0.48× by 10 GeV. Decomposing it (ratio = [ours / paper's own smooth expected-median curve] × [paper's expected curve / paper's real observed curve]) separates two genuinely different effects that happen to compound in the same mass window:

- *Real-data statistical fluctuation* (unclosable without the paper's actual per-event data): the paper's own observed/expected ratio has the **same V-shape**, dipping to ~0.42×–0.48× around 2.8–3.5 GeV before recovering to ~0.56×–0.57× at high mass — i.e. their real dataset happened to give them a noticeably better limit than their own median-expected sensitivity in exactly this mass range, an ordinary small-number Poisson effect, nothing to do with our modeling.
- *A genuine gap between our Asimov projection and the paper's own Asimov-like expected curve* (real-data noise removed from both sides): still dips to ~0.10× around 2.6–3.7 GeV, recovering only to ~0.27× by 9 GeV — so a real, still-partially-open modeling difference exists independent of data noise.

The mechanism behind that second piece: for m_χ = 1.36–2.64 GeV, the median of the predicted recoil spectrum sits at **~0.5 eV_ee** — essentially the *entire* rate is packed into an extreme exponential tail below the kernel's 60 eV floor. The whole computed limit for these masses is set by whatever tiny fraction of that tail happens to poke above threshold, a numerically fragile regime where small differences in tail resolution (input-spectrum fine-binning, the kernel's lowest sampled point, the quenching curve's exact floor behavior) get hugely amplified relative to the total rate. By 3.68 GeV only 14% of the spectrum has moved above 200 eV; by 5–10 GeV the spectrum is broadly spread through the well-covered range and the ratio recovers, because the calculation is no longer dominated by an extreme tail. This is not a single fixable bug — it's inherent numerical fragility in comparing a low-statistics extreme-tail calculation to a real analysis, most likely shared (if handled with different numerical care) by the paper's own real analysis in the same mass range.

**Where this leaves Phase 7.** Two real, structural fixes landed this pass (§6.8's folding bug, §6.9's background-efficiency curve) closed the mass-*dependent* divergence that motivated the whole investigation — the curve's shape now tracks the paper's own shape correctly, especially above ~5 GeV. The remaining gap at low-to-mid mass has three identified, non-arbitrary contributors — real-data statistical noise, the unverified 1×100 noise parameter, and extreme-tail numerical sensitivity at low mass — none of which point to a further closable bug in this codebase without either (a) the paper's real per-candidate data, or (b) a paper-stated 1×100 noise figure, neither of which exists in machine-readable public form. Investigation paused here by explicit decision, not because the trail went cold.

**Addendum — a sharper diagnostic that revises the above.** Fig. 11's caption states the shaded band around the paper's red observed line is "the expected sensitivity ±1σ" — i.e. the *same* background-only projected-sensitivity quantity as `DAMIC_2016_SNOLAB_expected_{lower,upper}.csv`, just drawn around the observed curve rather than as separate lines. Those existing files were WIMPyCCD-repo-sourced with uncertain provenance (§6.4); directly digitized replacements (`..._lower_diego.csv` / `..._upper_diego.csv`, same trusted vector/manual-digitization standard as the observed curve) were produced and checked against them.

Comparing against this band — the paper's own *clean, background-only projection*, with no real-data statistical noise on either side of the comparison — is a stricter test than comparing against their real observed curve: **1×1-only sits below the band's lower edge at every mass point except 1.36 GeV** (the single lightest testable mass, where it falls inside); the joint result sits below the lower edge at *every* mass with no exception. This is a meaningful revision to the picture above: real-data noise (§6.10's first decomposed factor) cannot be the explanation for a projection-vs-projection gap this large and this consistent, since no real data enters this specific comparison at all. The extreme-tail-sensitivity mechanism (§6.10, low mass only) and the unverified 1×100 σ_pix (joint only) remain real and quantified, but together they don't obviously cover a gap this systematic across the *entire* mass range, including 1×1-only at high mass where tail-sensitivity does not apply. **This reopens the question of a real, still-unidentified difference between this pipeline's Asimov projection and the paper's own** — the most likely remaining places to look, not yet checked: the paper's own quoted exposure/background/livetime numbers reproduced exactly rather than re-derived, and whether their "expected sensitivity" calculation folds in anything beyond a flat-Compton Asimov background (e.g. a different treatment of the systematic uncertainties mentioned in §6.1/Sec. VIII, which by the paper's own account move the 2 GeV limit by a factor of 1.5 and are not modeled here at all). Flagged here as the honest current state, not chased further without the user's direction.

### 6.11 Open items for further investigation (documented, not acted on)

A catalog of specific, evidenced leads for anyone picking this back up — ranked roughly by how promising each looks given everything found so far, not in the order discovered. None of these have been implemented; this is a map, not a to-do list.

**Most promising, given §6.10's addendum reopened the question of a real projection-vs-projection gap:**

1. **The paper's own stated systematic-uncertainty budget is not modeled here at all.** Sec. VIII: "Exclusion limits were generated changing the nuclear recoil ionization efficiency within its uncertainty... resulting in a change by a factor of ±1.5 in the excluded cross section at 2 GeV." There's also a stated background-composition uncertainty (65±10%/15±5%/20±5% bulk/front/back for 1×1, §6.4's efficiency discussion) that feeds into their quoted detection-efficiency uncertainty. This pipeline uses single central values throughout — no quenching-uncertainty band, no background-composition-uncertainty propagation. If the paper's own "expected sensitivity" calculation folds these in (plausible — a systematics-marginalized expected band is wider/weaker than a pure-statistical one), that would directly explain why our clean Asimov projection sits outside their clean expected band. **This is the single most direct, evidenced next thing to check** — quantifying it doesn't require new data, just re-reading Sec. VIII/IV more carefully for exactly what varies in their limit-setting procedure versus what's fixed.

2. **Non-uniform detector conditions across the real exposure.** Line 358 of the paper: "some CCDs in runs acquired between February and August 2015... where the pixel noise was relatively high (~2.2 e⁻)" — versus the quoted σ_pix=1.8 e⁻ used as a single fixed value here (and everywhere in this pipeline). The real 0.6 kg-day exposure is stated elsewhere (the acquisition-mode table, §6.1) to be an aggregate of several distinct date ranges per readout mode — if the paper's own expected-sensitivity calculation properly weights across a *mix* of noise conditions rather than one fixed 1.8 e⁻, that's a real, unmodeled effect that would make their true expected sensitivity weaker (band shifted toward higher cross-section) than a single-σ_pix Asimov projection predicts — independent of, and possibly compounding with, item 1.

**Real, quantified, but not independently closable (already investigated this pass, kept as open items rather than closed because no further data exists to close them with):**

3. **1×100 channel's σ_pix** (§6.10, item 3) — bounded sensitivity test done (1.8→4.0 e⁻), meaningfully narrows but never closes the low-mass gap. Would be fully resolved by a paper-stated or independently-measured 1×100-specific noise figure, which doesn't currently exist in accessible form.

4. **Extreme-tail numerical sensitivity at low mass** (§6.10) — mechanism identified (median E_ee ~0.5 eV for 1.4–2.6 GeV WIMPs, almost the entire rate below the kernel's 60 eV floor). Not yet tried: finer `nr_nbins` (currently 3000 over 0–30 keV_nr) specifically for light-mass rate generation, to check whether the tail itself is well-resolved at the *rate-generation* step (as opposed to the kernel-resolution test already done, which varied the *detector-response* sampling, not the input spectrum's own binning). A genuinely different knob from the kernel-resolution test in §6.10 item 2, not yet turned.

**Smaller, likely-minor items, not chased because early estimates put them well below the scale of the remaining gap:**

5. **Missing Lewin-Smith R1 form-factor correction** (found while investigating the high-mass divergence, before the folding bug was found — see the session record). Our simplified Helm form factor uses `rn` directly as the Bessel argument; a fuller convention (used by `wimprates`) corrects to `R1 = sqrt(rn² - 5s²)`. Quantified: <5% effect on F² even at 30 keV_nr — almost certainly not worth revisiting unless every larger lever above is exhausted first.

6. **Halo parameter precision.** The paper's quoted standard-halo values (v0=220, vE=232, vesc=544 km/s, ρ=0.3 GeV/cm³, Sec. VIII) are used as exact fixed points here, matching what the paper itself does (standard practice is not to vary halo parameters in the systematic budget) — flagged only for completeness, not believed to be a real lever.

7. **Exposure timing granularity** (`backgrounds.timing.exposure_time_s=1800`, a 30-minute readout-cadence assumption) — not verified against the paper's actual per-image acquisition timing (line 239 states 840 s / 20 s image readout time for 1×1/1×100, a related but distinct number from the 30-minute figure used here). Worth a citation check before ruling out, but no evidence yet that this matters at the scale of the observed gap.

**The one thing that would close everything, but is out of reach.** The paper's real per-candidate energies (31 in 1×1, 23 in 1×100, §V) aren't published in machine-readable form beyond the final exclusion curve. Fitting to the real events instead of an Asimov background-only projection would eliminate the real-data-noise contributor to the gap entirely (§6.10) and might also surface whether items 1–2 above are actually where the remaining projection-vs-projection gap lives. Not something this codebase can obtain on its own; would require either the authors' own data release or a formal request.

### 6.12 Revisiting the residual (student-example pass): items 1 and 2 directly ruled out, item 4 tested and ruled out

Prompted by a side-by-side comparison against an older plot during the student-example work (see `Student_Examples_WIMP_DAMIC_SNOLAB.md`): the reproduction shown there is worse-looking than a plot from earlier in the project, but this is **not a regression from any code/config change in this pass** — every variant checked (the student-example single-channel 1×1-only config, the pre-existing joint 1×1+1×100 config with and without the Izraelevitch quenching extension, and the old `_finekernel` variant) reproduces the same residual, at the same scale, consistent with §6.9's own numbers. The apparently-better older plot almost certainly predates §6.8's fold-bug fix and §6.9's background-efficiency-curve fix — i.e. it was likely a fortuitous cancellation between two real, since-fixed bugs (one discarding ~90% of folded signal, one over-crediting background efficiency), not a more-correct result. Reverting either fix to regain the nicer-looking plot was considered and rejected.

**Confirms: both cases show the residual, joint is worse, not better.** Directly reran and replotted both `configs/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch.json` (joint, Izraelevitch quenching) and `configs/wimp_nucleon_damic2016_limit_repro_joint.json` (joint, chavarria-only) against `--wimp-minimal`'s trusted observed+band — both diverge below the band starting well before 10⁴ MeV, matching §6.10's addendum finding that "the joint result sits below the lower edge at every mass with no exception," worse than 1×1-only's single-mass exception at 1.36 GeV.

**Item 1 (paper's stated systematic-uncertainty budget) — directly ruled out by re-reading the source LaTeX** (`Literature/arXiv-1607.07410v2/Limit2016_clean.tex`, not previously done at this level of care). Line 565: "The wide red band presents the expected sensitivity of our experiment generated from the distribution of outcomes of 90% C.L. exclusion limits from a large set of Monte Carlo **background-only** samples" — a pure-statistical construction, textually and structurally separate from the following paragraph (lines 569–573) on systematic-uncertainty checks (Fano factor, quenching-curve variation, detection-efficiency-curve variation including background composition), which is about robustness of the *observed* limit, not the band. That paragraph also states these systematics had "negligible impact... for WIMP masses >3 GeV" and only ±1.5× at 2 GeV even at their most impactful — an order of magnitude too small to explain the observed ~5–12× gap, and inapplicable to the band comparison in any case. §6.11 item 1 is withdrawn as a candidate.

**Item 2 (non-uniform detector conditions / mixed σ_pix) — directly ruled out by the same re-read.** Lines 357–359: CCDs with elevated pixel noise (~2.2 e⁻, Feb–Aug 2015, from a light leak) were explicitly **excluded** from the analysis ("Images for which there is a significant discrepancy... were excluded"), not averaged into the quoted exposure. There is no unmodeled noise-mixing effect to find here — the paper's own σ_pix=1.8 e⁻ is the value actually used throughout its real dataset. §6.11 item 2 is withdrawn as a candidate.

**Item 4 (rate-generation tail resolution, previously "not yet tried") — tested directly, ruled out as the driver.** `_make_recoil_grid_keV` (`python/ccdarkphys/wimp_nucleon/entry.py`) uses a *linear* 0.001–30 keV grid; at the production setting of 3000 bins this is 10 eV_nr per bin, comparable to or coarser than the entire signal-relevant region for the lightest masses (§6.10: median E_ee ~0.5 eV for 1.36–2.64 GeV WIMPs). Directly quantified: refining to 60,000 bins (0.5 eV_nr) roughly **doubles** the resolved fraction of raw nuclear-recoil rate below 20 eV_nr at mχ=1.36 GeV (7.85%→14.0%) and shifts it meaningfully at 2.12/3.68 GeV too — real, and previously unverified. But run through the *full* pipeline (fine-binned CSVs generated for mχ=500 MeV–2.36 GeV, `Emax_keV`/halo/quenching otherwise unchanged, fed through the same single-channel `ccdarksens_scan_generic` scan), the effect on the final upper limit is **<1% at every mass from 1.36–2.36 GeV, and 7% at the single lowest resolved point (1.21 GeV)** — nowhere near the order-of-magnitude gap. The raw-rate redistribution found near threshold does not survive meaningfully once folded through the detector kernel's own resolution (readout noise, diffusion, the coarser 400-bin electron-equivalent histogram). §6.11 item 4 is downgraded from "untried, promising" to "tried, real but small effect, not the driver."

**Where this leaves the investigation.** Three of §6.11's candidates are now closed with direct evidence (two ruled out textually, one ruled out numerically); what remains open is item 3 (1×100 σ_pix, bounded but not closable without paper data — joint-only, so it cannot explain the 1×1-only residual on its own), items 5–7 (already estimated as minor), and the extreme-tail *physical* sensitivity mechanism itself (§6.10) — which is not a bug to fix, just an inherent numerical fragility of this mass range for either pipeline. The halo velocity-integral implementation (`python/ccdarkphys/wimp_nucleon/halo.py::mean_inverse_speed`) was inspected as a fresh candidate — it is a direct, standard Lewin-Smith/Savage-Freese-Gondolo piecewise truncated-Maxwell-Boltzmann port already validated for WIMPyCCD parity (per its own docstring); no structural issue found, but not independently re-derived from scratch in this pass. Paused here, same as §6.10 — the remaining candidates all require either the paper's own per-event data (out of reach) or a more open-ended physics audit than a single pass justifies.

**A fifth, real contributor found afterward: the wrong `q_threshold` was used for every plot in this example set.** While chasing why the freshly-built student-example plots looked visibly worse than an earlier, better-tracking reproduction from earlier in this same project, the session transcript was searched to recover the exact historical `ccdarksens_plot_limit` command that produced the better plot. It never passed a `q_threshold` argument at all, relying on the app's own hardcoded default (`apps/ccdarksens_plot_limit.cc`: `double q_thr = 2.71; // default ~90% CL, 1 dof` — the standard one-sided Δχ²=2.71 chi-square convention). Every plot built for this example set instead explicitly passed `1.642374415149816` (the asymptotic Cowan et al. one-sided value), copied by habit from the pydme-style channels (DM-electron, Migdal, dark photon), where it is the correct, validated value. WIMP-nucleon uses `profile_minimizer: brent`, not `pydme` — a different UL-construction methodology for which the app's own default is the intended match, not the pydme-specific asymptotic value. Confirmed the underlying `scan_generic.root` q-histograms are essentially unchanged between the historical run and a fresh rerun (agreement to <3% at every mass checked, ruling out data drift as an explanation) — the visual difference is attributable to the `q_threshold` argument alone. Switching to `2.71` measurably tightens the reproduction's agreement with the band (see updated plots, `docs/wimp_snolab_combined.pdf` / `wimp_damicm_combined.pdf`), though it does not fully reproduce the best historical plot found nor close every remaining gap — some residual at low-to-mid mass persists, consistent with §6.10's other identified (and still-open) contributors.

**Checked separately for the joint 1×1+1×100 case: `2.71` helps there too, but far less than for 1×1-only.** Replotting `outputs/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch/scan_generic.root` at `q_threshold=2.71` (same file used throughout this section) still shows the black curve diverging clearly below the band from ~1.5 GeV through ~9 GeV — roughly half a decade to a full decade too optimistic in that stretch, only converging back near the band close to 10 GeV — versus 1×1-only's much smaller residual (a brief dip around 1.5–4 GeV) at the same threshold. This is consistent with, and now reconfirmed under the corrected threshold for, this section's standing explanation (§6.10 item 3): the joint likelihood's extra divergence traces to the 1×100 channel's unverified σ_pix=1.8 e⁻ (carried over from the 1×1 section, never separately stated by the paper for 1×100), not to the q_threshold convention. The 1×1-only scope already adopted for the student example set (`Student_Examples_WIMP_DAMIC_SNOLAB.md`) remains the better default on these grounds, independent of the LEE-wiring reason it was originally chosen for.
