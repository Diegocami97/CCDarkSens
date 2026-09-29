<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_WIMP_DAMIC_SNOLAB.md -- WIMP-nucleon SI, DAMIC at SNOLAB (2016) reproduction, baseline + low-energy-excess variant.
-->

# Student Example — WIMP-Nucleon SI, DAMIC at SNOLAB (2016) Reproduction

Reproduction of the published DAMIC at SNOLAB (2016) WIMP-nucleon spin-independent limit (PhysRevD.94.082006), with a second variant that turns on the low-energy-excess (LEE) background. Companion: [`Student_Examples_WIMP_DAMIC_M_Projection.md`](Student_Examples_WIMP_DAMIC_M_Projection.md) (a hypothetical DAMIC-M projection, not a reproduction).

---

## 1. What this example does

Two configs, same detector/exposure/background, differing only in whether the LEE background term is enabled:

| | Baseline | +LEE |
|---|---|---|
| Config | [`wimp_damic_snolab_2016.json`](../configs/examples/wimp_damic_snolab_2016.json) | [`wimp_damic_snolab_2016_lee.json`](../configs/examples/wimp_damic_snolab_2016_lee.json) |
| LEE background | off (`backgrounds.lee.enabled` absent/false) | **on**, `preset: "damic_snolab_2023_skipper"` (arXiv:2306.01717) |

| Item | Value |
|---|---|
| Signal model | WIMP-nucleon, spin-independent, elastic |
| Observable | `cluster_energy` — **1×1 channel only** (scope reduced from an earlier joint 1×1+1×100 attempt — see the boxed note below) |
| Exposure | DAMIC at SNOLAB (2016): 0.01 kg × 60 days |
| Background | Flat Compton rate 5478.75 events/(kg·yr·keV) + Fig. 9-digitized 1×1 background efficiency curve (paper values) |
| Data | Asimov |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σ_SI(mχ) |
| True-energy floor | 60 eV (DAMIC 2016's own value) |

> **Why 1×1-only, and a real bug this comparison surfaced.** These two configs were originally built as a joint 1×1+1×100 likelihood (matching the paper's own Fig. 11 construction), the same shape as `configs/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch.json`. Running that joint pair produced **bit-identical upper limits, baseline vs. +LEE, at all 28 mass points** — a dead giveaway that the LEE term wasn't actually being applied, not that its effect was merely small. Traced to `RunJointChannelScan` (`apps/ccdarksens_scan_generic.cc`): it builds each channel's background directly (`MakeClusterEnergyFlatBackground`/`...WithEfficiency`) and never checks `cfg.backgrounds().has_lee_bkg` or calls `MakeClusterEnergyLEEBackground` — that call only exists in the single-channel path's `MakeBackground()` (`BackgroundFactory.cc`). `backgrounds.lee` is parsed correctly into the config; it's just never read by the joint path — a silent gap, since `ConfigManager` has no unknown/unused-field validation anywhere. Per explicit instruction, scope was reduced to 1×1-only for now (which already correctly applies LEE via the working single-channel path) rather than patching the joint path immediately — confirmed the two configs now produce genuinely different upper limits at every mass point (e.g. at mχ=10 GeV: 9.76×10⁻⁴¹ baseline vs. 1.22×10⁻⁴⁰ with LEE, ~25% weaker as expected from added background). A joint-channel LEE fix is a known, deferred follow-up.

These two configs are still the **first working on/off comparison** of the LEE feature in this repo — `MakeClusterEnergyLEEBackground` existed in the codebase already but had never been exercised end-to-end in a committed config before this example set.

---

## 2. Prerequisites

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_generic ccdarksens_plot_limit
```

Needs the WIMP-nucleon rate/kernel machinery: `NoiseTailCalibrator`, `ClusterFitMC`, `ClusterFitEngine` (all built as part of `ccdarksens_scan_generic`).

---

## 3. Pipeline overview

```
Step 1   Calibrate the pure-noise ΔLL tail (per channel)             -- part of the scan run
           |
Step 2   Build the fit-recovery kernel K[E_true, E_reco]             -- part of the scan run
           |
Step 3   Fold signal + background through the kernel, run the PLR    -- part of the scan run
           |
Step 4   Plot the limit curve against the paper's own digitized curves
```

**What the pipeline is physically doing:** for the WIMP-nucleon channel, the observable is a per-event reconstructed cluster energy rather than an n_e or pattern bin. `NoiseTailCalibrator` first characterizes the pure-noise fluctuation tail of the Nelder-Mead cluster fit (needed to set the analysis threshold), then `ClusterFitMC` builds a response kernel mapping true nuclear-recoil energy to reconstructed energy, folding in diffusion, readout noise, and (for these configs) the paper's own 60 eV true-energy floor, for the single 1×1-channel readout (see the boxed note in §1 for why this is 1×1-only rather than the paper's joint 1×1+1×100 construction). The LEE variant adds an extra exponential background term (`backgrounds.lee`) on top of the flat rate, using the preset decay constant from the cited skipper-CCD low-energy-excess measurement.

---

## 4. Step 1 — Rate/kernel generation (informational — happens inside the scan run)

Unlike the other channels, there is no separate offline rate-generation step here — `ccdarksens_scan_generic` runs the noise-tail calibration and kernel build internally at scan time for this analysis space. If you want to inspect these steps in isolation first, use `ccdarksens_calibrate_noise_tail` and `ccdarksens_build_cluster_fit_kernel` directly (see the app table in the top-level project documentation).

---

## 5. Step 2 — Run the scans

```bash
build/ccdarksens_scan_generic configs/examples/wimp_damic_snolab_2016.json
build/ccdarksens_scan_generic configs/examples/wimp_damic_snolab_2016_lee.json
```

Outputs: `outputs/wimp_damic_snolab_2016/scan_generic.root`, `outputs/wimp_damic_snolab_2016_lee/scan_generic.root`.

**Runtime:** the noise-tail calibration and kernel build both run Monte Carlo trials per grid point — this is one of the heavier scans in the repo. For a quick check, copy the config and reduce the mass grid.

---

## 6. Step 3 — Plot the limit curve

**The recommended view puts both curves on one canvas against the trusted reference and its band:**

```bash
build/ccdarksens_plot_limit --batch --wimp-minimal \
  --out-pdf outplots/wimp_snolab_combined.pdf \
  outputs/wimp_damic_snolab_2016/scan_generic.root "DAMIC SNOLAB 2016 (reproduced)" \
  outputs/wimp_damic_snolab_2016_lee/scan_generic.root "DAMIC SNOLAB 2016 +LEE" \
  2.71 heavy
```

**Use `q_threshold = 2.71`, not the `1.642374415149816` value used for the pydme-style channels (DM-electron, Migdal, dark photon).** `2.71` is `ccdarksens_plot_limit`'s own hardcoded default (`double q_thr = 2.71; // default ~90% CL, 1 dof`, the standard Δχ²=2.71 one-sided chi-square convention) — the right one for this channel because WIMP-nucleon uses `profile_minimizer: brent`, not `pydme`; `1.642374415149816` is the asymptotic Cowan et al. one-sided value specifically validated for the pydme-style channels' own UL convention, not this one. This was found by tracing an earlier, better-tracking version of this reproduction back through the session history to the exact command that produced it, which never overrode the default. Using 2.71 measurably tightens the reproduction's agreement with the band; it does not fully close every remaining gap (see §7 below and `docs/ClusterFitMC_Design.md` §6.12).

`--wimp-minimal` always reads the upper limit from the q(mχ, σ) histogram and loads only the paper's own trusted digitization: the observed curve (Fig. 11) and its expected ±1σ band — it deliberately drops the plain `--wimp` mode's extra secondary/uncertain-provenance curves (an alternate "radomir" observed digitization, an uncertain-provenance curve, the expected median, and a literature DAMIC-M projected-1kg-year curve — all four are still loaded internally, just not drawn or added to the legend in minimal mode) so the comparison isn't cluttered by curves this example doesn't need. A rendered copy is checked in at [`wimp_snolab_combined.pdf`](wimp_snolab_combined.pdf).

If you want each curve in isolation, or the full literature context (all four extra curves), swap `--wimp-minimal` for `--wimp` and/or drop one of the two ROOT/label pairs.

---

## 7. Known, documented caveat — read before drawing conclusions

This reproduction tracks the paper's real curve within roughly 1.5–3× through most of the mass range at `q_threshold=1.642374415149816`, and sits outside the paper's own published expected ±1σ band at nearly every mass — a real, unexplained residual that survived two genuine bug fixes made while building this pipeline (a folding-quadrature error that had discarded ~90% of signal, and a background that had used the signal's own efficiency curve instead of its own). Using the correct `q_threshold=2.71` for this channel (see §6 above) measurably tightens the agreement — the reproduction now tracks much closer to the band through most of the range — but does not close every residual: some gap remains at low-to-mid mass. This affects **both** the 1×1-only config used here and the joint 1×1+1×100 likelihood — the joint case is worse even at the corrected `q_threshold=2.71` (reverified directly: it still diverges clearly below the band from ~1.5 GeV through ~9 GeV), versus 1×1-only's much smaller residual (a brief dip around 1.5–4 GeV) at the same threshold — consistent with the joint likelihood's extra gap tracing to the 1×100 channel's unverified σ_pix (§6.11 item 3), not to the threshold convention. See `docs/ClusterFitMC_Design.md` §6.10–§6.12 for the full write-up: §6.12 (added during this example set's development) directly ruled out the two previously-top-ranked candidate explanations by re-reading the paper's own LaTeX source (the systematic-uncertainty budget only applies to the observed curve, not the expected-sensitivity band it's compared against; the noisy ~2.2 e⁻ CCD runs were explicitly excluded from the analysis, not averaged in), tested/ruled out a third (finer rate-generation binning near the kinematic threshold has a real but small, <7%, effect — not the order-of-magnitude driver), and identified a fourth, real contributor: this whole example set had been plotted with the wrong `q_threshold` for a non-pydme (`brent`-minimized) channel. Use this example to learn the WIMP-nucleon pipeline's mechanics, not as a fully validated reproduction of the paper's absolute exclusion curve.

The LEE variant has not been independently cross-checked against any external LEE-inclusive analysis — it demonstrates that the feature runs and shifts the limit, not that the shift is quantitatively correct.
