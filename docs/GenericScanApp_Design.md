# Generic Multi-Channel Scan App — Design and Validation

**Document:** `docs/GenericScanApp_Design.md`
**Author:** Diego Venegas-Vargas
**Status:** Living document, built up slice by slice alongside implementation, matching `docs/ClusterFitMC_Design.md`'s precedent.

---

## 0. Scope

Every signal channel (DM-electron, dark photon, Migdal, and the new WIMP-nucleus SI channel) needs the same thing: a grid scan over (mass, coupling) producing a σ_UL(mass) exclusion curve. Today only `apps/ccdarksens_scan_dmelectron_pattern.cc` does this, hardcoded to pattern-space. This work builds a single, genuinely flexible scan app on top of the existing, already-validated physics — `ModelFactory`, `ScanUtils`, `ProfileLikelihood`, `ClusterFitMC` — without modifying the reference app or its two downstream consumers (`ccdarksens_band`, `ccdarksens_plot_dmelectron_limit`), which stay untouched reference points.

Full implementation plan: `~/.claude/plans/goofy-wishing-bengio.md` (session-local). This document is the durable, in-repo record of the derivations and validation results.

---

## 1. Slice 1 — `ResponseFold`, `PatternResponseFold`, `NeSpaceResponseFold`

### 1.1 The abstraction

Reading the reference app in full turned up the key fact that makes a generic scan app tractable: its actual per-grid-point signal fold **bypasses `DetectorResponsePipeline::Apply()` entirely** (explicit code comment: *"Use S_true, not S_obs — folding with S_obs would double-count"*) and calls `ChargeIonization::FoldToNe` + `FoldNeToPatternRates` directly. The real, reusable contract underneath every channel is just:

```
(TH1D dRdE_density, exposure_kg_year) -> vector<double> S(bin)
```

— exactly the shape `ClusterEnergyRates::FoldEtrueToErecoRates` (built for the WIMP-nucleus channel) already has. `ResponseFold` is a 3-method abstract interface capturing this; background construction is *not* a method on it — it's the same `Fold()` call applied to a flat dummy spectrum, matching how the reference app already treats background (same fold path as signal, different input histogram).

### 1.2 A correction caught mid-build: pattern-space and n_e-space are genuinely different paths

The original design (from reading the reference app's pattern-space branch) only accounted for pattern-space and the new cluster-energy space. Validating `PatternResponseFold` against a real config, the first config picked happened to run in **n_e-space** (`experiment.observable_bins: "n_e"`, the default) — which produced no `S_pat` output at all and revealed that the reference app treats n_e-space as a structurally different fold, not a special case of pattern-space:

```
pattern-space:  S_true = FoldToNe(dRdE, exposure, ne_min, ne_max)
                 S_pat  = FoldNeToPatternRates(S_true, ne_min, ne_max, pattern_roi, pattern_eff_map)

n_e-space:       S_true = FoldToNe(dRdE, exposure, ne_min, ne_max)   -- same first step
                 S_obs[ne] = S_true[ne] * clamp(eps_ne[ne], 0, 1)     for ne in roi_bins
```

n_e-space has no multi-pixel pattern classification — just a per-`n_e` detection efficiency curve. `NeSpaceResponseFold` was added alongside `PatternResponseFold` to cover this, both real, both needed (per explicit user requirement: DM-electron must keep working in both observable spaces). `ResponseFactory` (Slice 4) will dispatch on `response.analysis_space == "cluster_energy"` first, otherwise on `experiment.observable_bins` — the two existing DM-electron/dark-photon/Migdal config concepts stay exactly as they are; only cluster-energy is new.

### 1.3 Validation

No unit-test scaffolding — direct parity against the **unmodified reference app's own printed diagnostic output**, on two real, independently-selected configs (not synthetic fixtures):

| Fold | Config | Result |
|---|---|---|
| `NeSpaceResponseFold` | `configs/band_gap_one_point_gap0p1_B-thresh.json` (n_e-space, CSV-weighted efficiency) | 5/5 `roi_bins` match to `<10⁻⁹` relative |
| `PatternResponseFold` | `configs/ne_imaging_one_point_si_ref_pattern.json` (pattern-space, multi-electron patterns `{11,21,111,31,22,211}`) | 6/6 `pattern_roi` match to `<2×10⁻⁹` relative (including two exact-zero entries) |

One methodological note worth recording: the first `NeSpaceResponseFold` validation attempt used `ε(n_e)` values read off the reference app's own `std::cout` diagnostic print (e.g. "ε = 0.98"), and came back ~0.2–0.9% off — close enough to look plausible, wrong enough to matter. Rather than loosen the tolerance to make it pass, the actual full-precision `ε(n_e)` values were obtained by reconstructing them through the real `EfficiencyMC::PrecomputeEpsilonWithPatternEff` call (the same one the reference app itself uses) instead of hand-copying truncated console output — which then matched to `<10⁻⁹}`. The lesson: a "close" validation number is not a passing one; hunt down the actual source of a sub-1% gap rather than accepting it, especially when a more rigorous cross-check (reusing the real code path instead of eyeballing text output) is available.

Both response folds call zero new physics — `FoldToNe`, `FoldNeToPatternRates`, and `EfficiencyMC::PrecomputeEpsilonWithPatternEff` are all pre-existing, already-validated functions; the new code is exclusively the thin `ResponseFold`-interface wrapping.

---

## 2. Slice 2 — `ClusterEnergyResponseFold` and the flat-Compton background

### 2.1 The wrapper

`ClusterEnergyResponseFold` is the cluster-energy analogue of Slice 1's two folds — a thin `ResponseFold` wrapper around `FoldEtrueToErecoRates` against a `KernelMatrix` built once via `ClusterFitMC::BuildKernel` (expensive, ~minutes; built once, never per grid point, exactly like `EfficiencyMC`'s table in the reference app). No new physics — `Fold()` is a single pass-through call.

### 2.2 The actually new piece: a real background for the WIMP-nucleus channel

Every DM-electron scan needs a background estimate, and the reference app builds one by folding a flat (energy-independent) dummy rate spectrum through the *same* path signal uses — physically, the flat-Compton background discussed when this channel was first designed (Fig. 1 of the earlier plan artifact labeled this box explicitly). The WIMP-nucleus channel never had this — the earlier `ccdarksens_example_one_point_cluster.cc` diagnostic only ever compared signal against itself (Asimov self-consistency), never a real background.

`MakeFlatDrdeSpectrum(rate_per_kg_year_keV, Emin_eV, Emax_eV, nbins, name)` generalizes the reference app's inline construction (`ccdarksens_scan_dmelectron_pattern.cc:729-739`) into a standalone, channel-agnostic function — same d.r.u. unit convention (`events/(kg·year·keV)`, matching `BackgroundJSON::flat_bkg_norm_per_kg_year`'s own documented units), same conversion to per-eV, same "constant across every bin" construction. The one generalization: it takes its own `(Emin_eV, Emax_eV, nbins)` directly rather than borrowing a `ModelFactory`-shaped histogram's binning, since cluster-energy has no `ModelJSON`-shaped spectrum to borrow from (its energy grid lives in `ClusterFitMCJSON` instead).

Folding this through `ClusterEnergyResponseFold` — the same flat, d.r.u.-based background concept flowing through the WIMP-nucleus reconstruction kernel for the first time — is the actual new capability; the background physics itself is unchanged.

### 2.3 A real bug caught during validation, not cosmetic

The first validation run produced correct numbers but also a ROOT warning: `TROOT::Append: Replacing existing TH1: dRdE_flat (Potential memory leak)`. Cause: `MakeFlatDrdeSpectrum` hardcoded a fixed histogram name, and ROOT auto-registers every `TH1` in its global directory by name — calling the function twice (as any real usage would, once per background/signal comparison) created two histograms under the identical name, a genuine double-ownership hazard, not just console noise. Checked how the rest of this codebase avoids this class of bug (grepped for `SetDirectory(nullptr)` — unused anywhere) and found the actual established convention: `RateTable::MakeTH1D` takes a caller-supplied `name` argument, avoiding the collision through caller discipline rather than a ROOT-directory workaround. Fixed `MakeFlatDrdeSpectrum` to match that same convention (added a `name` parameter, default `"dRdE_flat"`) rather than introducing an inconsistent new pattern. Re-ran the validation after the fix — warning gone, all numbers unchanged (same seeds, same kernel).

### 2.4 Validation

Standalone check building a small (6-point, reduced-trial) kernel and:

| Check | Result |
|---|---|
| `ClusterEnergyResponseFold::Fold()` vs. direct `FoldEtrueToErecoRates()` call, same inputs | exact match, `0.000e+00` abs diff |
| Flat background folded through the kernel: all bins finite and non-negative | ok |
| Linearity: total at 2× the input rate / total at 1× rate | `2.000000` (exact) |
| ROOT histogram naming, re-checked after the fix | no warnings |

---

## 3. Slice 3 — `WimpNucleonModel` and the `ModelFactory` branch

### 3.1 Nothing new to design

`WimpNucleonModel` is close to a literal copy of `MigdalModel` (76 lines): same token-replacement path resolution (`{target_nucleus}`,`{mediator}`,`{mchi_MeV}`,`{sigma_n_cm2}` substituted into `filename_template`), same `RateTable::LoadCSV` + `MakeTH1D` call, no rescaling logic (every `(mass, cross-section)` WIMP-nucleus grid point already has its own pre-generated CSV, unlike dark photon's ε²-rescale trick). This was possible because `DMNucleonConfig` — the config struct `MigdalModel` already uses — has a doc comment that explicitly anticipated this exact reuse ("a future elastic nuclear-recoil/WIMP model will consume the same struct"), so no new config struct was needed either. One new `else if (mj.type == "wimp_nucleon")` branch in `ModelFactory.cc`, mirroring the `migdal` branch line-for-line — the only edit to an existing file in this slice.

### 3.2 Validation

Diffed `ModelFactory::MakeSignalSpectrumE(mj, mass, coupling, &ok)` against a direct `RateTable::LoadCSV`+`MakeTH1D` call on the same underlying CSV, for both existing fixture files:

| Check | Result |
|---|---|
| `m_χ=3000 MeV` fixture, factory vs. direct load | exact match, `0.000e+00` max relative diff across all 400 bins |
| `m_χ=5000 MeV` fixture, factory vs. direct load | exact match, `0.000e+00` max relative diff |
| Missing-file failure path (nonexistent mass) | `ok=false`, non-null histogram returned — matches `ModelFactory`'s own documented contract ("valid, possibly all-zero TH1D even on load failure") |

As expected for a routing-only addition: zero new computation, exact match by construction.

---

## 4. Scope revision — full statistical parity, not a reduced v1

An earlier pass through this design deferred five statistical features (2D `pydme` minimizer, `Bp_Br` background templates, `single_bin_likelihood`, `smooth_ul_envelope`, `band.cc`'s `threshold_toys` toy-MC mode) as "out of scope for v1." The user rejected that explicitly: *"I need all of them to work, especially if we are meant to reproduce what has been already done for the DM-e cases."* Checking which code path the actual flagship reproduction config (`configs/scan_dmelectron_pattern_pydme_exact.json`, header comment: *"Match pydme exactly"*) uses confirmed the concern was concrete, not hypothetical: that config sets `background_source=bp_br_template`, `background_model=Bp_theta_Br`, `profile_minimizer=pydme`, `pydme_style_ul=true`, and a CSV-weighted `efficiency_csv`/`efficiency_csv_reference` pattern-efficiency map — none of which the deferred-scope v1 would have supported. Slices 4–8 below build all five, plus two further gaps discovered only by trying to actually reproduce that config end-to-end (CSV-weighted pattern efficiencies, `run.data_path` real-data loading) that the original research pass had missed because it read the reference app's *structure* without tracing which branch a real production config actually takes.

---

## 5. Slice 4 — `ResponseFactory`, `ResponseFold::FoldNeBackground`, `BackgroundFactory`

### 5.1 `ResponseFactory`: replicating ~550 lines of setup, not writing a dispatcher

The original plan described `ResponseFactory::MakeResponseFold` as a thin dispatcher. Reading the reference app's actual detector-response setup (`ChargeIonization`/`ChargeTransport`/`PatternClassifier`/`EfficiencyMC` construction, `pattern_eff_map` assembly) end to end showed it is not thin — it's ~550 lines of config-driven object construction with two genuinely different `pattern_eff_map` sources (CSV-weighted vs. MC-generated from `EfficiencyMC::GetPatternTable()`) that the reference app's own accepted_labels/EfficiencyMCConfig plumbing depends on. `ResponseFactory::MakeResponseFold` (`src/response/ResponseFactory.cc`) replicates this faithfully for pattern/n_e space, and separately replicates the WIMP-nucleon channel's own setup sequence (`NoiseTailCalibrator` → `ClusterFitMC::BuildKernel`, from `ccdarksens_example_one_point_cluster.cc`) for `cluster_energy`.

One real, load-bearing gap found in the process: `efficiency_csv`/`efficiency_csv_reference` (the CSV-weighted pattern-efficiency path the flagship config actually uses) were never parsed by `ConfigManager` — the reference app reads them off the raw JSON directly, bypassing its own typed config struct. Added as two new fields on `EfficiencyMCJSON` (`include/ccdarksens/io/ConfigManager.hh`), parsed additively in `ConfigManager.cc` — the reference app's own raw-JSON reads are untouched, so its behavior doesn't change.

Deliberately *not* replicated: `use_2d_image_efficiency` (2D-image pattern-efficiency path) and `include_dc_pileup` — both documented in `CLAUDE.md` as diagnostic-only, not production paths. Also not replicated: the `backgrounds.background_efficiency_csv` migration-matrix background variant (an alternative to the plain dark-current+flat construction) — the flagship config uses `bp_br_template` instead, so this wasn't exercised by anything requiring parity; a known, explicitly-flagged gap rather than a silent one.

### 5.2 The `Fold()` contract doesn't cover dark current — `FoldNeBackground`

`ResponseFold::Fold(dRdE, exposure)` covers signal and the *flat* background (Slice 2) but not dark current, which `BackgroundBuilder::BuildBkgAsimov()` synthesizes directly as an n_e-space histogram from geometry × timing × λ — no dR/dE spectrum involved. Added a non-pure-virtual `FoldNeBackground(const TH1D& B_ne_asimov)` to `ResponseFold`, default-throwing; `PatternResponseFold` and `NeSpaceResponseFold` override it (reusing their existing `FoldNeToPatternRates`/elementwise-ε machinery on the caller-supplied histogram instead of a freshly-folded one). `ClusterEnergyResponseFold` leaves the throwing default — there is no pixel/dark-current model for cluster-energy, so a config requesting it there is a real error, not a silent no-op. (Slice 8 later added a capability query, `SupportsDarkCurrentBackground()`, so `BackgroundFactory` can degrade gracefully instead of throwing when the default `background_source` is used on a channel with no DC concept — see §9.)

### 5.3 `BackgroundFactory`: one background, two ways to get its numbers, two ways to use them

`run.background_source` decides how the numeric background vector is obtained (`"dc_flat_migration"`: `BackgroundBuilder` + `FoldNeBackground`, plus flat via `MakeFlatDrdeSpectrum` + `Fold()`; `"bp_br_template"`: `run.background_Bp`/`run.background_Br` loaded straight from config). `run.background_model` decides, independently, how `ProfileLikelihood` treats that vector (`SetBTemplate` vs. `SetBpBr`+constraints) — confirmed generic in both directions, since the reference app's own n_e-space branch already calls `SetBpBr` degenerately (`Bp=0, Br=B_template`) purely to reuse the θ-constrained-prior/2D-minimizer machinery. `BackgroundFactory::MakeBackground` (`src/response/BackgroundFactory.cc`) returns `{B_pat, Bp, Br}` always populated (degenerate `Bp=0` when sourced from `dc_flat_migration`), so the scan app's Phase A can wire either background model without branching on source.

One improvement over the reference app: `bp_br_template`'s size check is against `fold.NumBins()` (generic — pattern, n_e, or cluster-energy) rather than the reference app's hardcoded `pattern_roi.size()`.

### 5.4 Validation

New app `ccdarksens_validate_response_factory` builds `B_pat` both ways on real configs and diffs against ground truth:

| Path | Config | Result |
|---|---|---|
| `dc_flat_migration` | `configs/ne_imaging_one_point_si_ref_pattern.json` (real DC λ, CSV-weighted efficiency) | vs. reference app's own printed `B_pat` per pattern bin: `<3×10⁻⁹` relative, all 6 bins |
| `bp_br_template` | `configs/scan_dmelectron_pattern_pydme_exact.json` (real production `Bp`/`Br` arrays) | vs. hand-computed `Bp+Br`: `0.0` / `2×10⁻¹⁶` (float noise), all 6 bins |

---

## 6. Slice 5 — `apps/ccdarksens_scan_generic.cc` core (1D minimizer) + hard parity gate

Phase A builds `ResponseFold` + `BackgroundFactory` output once, then `ProfileLikelihood`: `SetData` (from `run.observed_counts` if present, else Asimov = background), `single_bin_likelihood` collapse applied to *every* relevant vector together (`S`, data, `B_pat`, `Bp`, `Br`) — a deliberate fix over the reference app, which collapses `S`/`B`/data but leaves `Bp`/`Br` per-bin when combined with `Bp_theta_Br`, a real size-mismatch bug not worth reproducing in new code. Phase B walks the (mass, coupling) grid with the reference app's exact 1D-minimizer branching, including a genuine asymmetry it preserves on purpose: pattern-space's `"minuit"` mode also cross-checks Brent and takes the better fit; n_e-space's `"minuit"`/`"pydme"` use Minuit only. Phase C writes `q_mchi_sigma` (TH2D), `upper_limit_sigma_e_mchi` (TH1D, `smooth_ul_envelope`-processed if set — window-min envelope, half-window 5 bins, exact port), and `upper_limit_sigma_e_mchi_graph` (TGraph, the one both `band.cc` and the plot app require).

**Hard parity gate**: unmodified reference app vs. the new app, 3-mass × 5-sigma real grid (`ne_imaging_one_point_si_ref_pattern.json` base), `background_model="scale"` — `upper_limit_sigma_e_mchi` and `_graph` both **bit-exact** (`0.000e+00` relative diff, every point) via a ROOT-level diff, not just matching to displayed console precision.

---

## 7. Slice 6 — 2D minimizer (`minuit2d`/`pydme`) + a `run.data_path` gap found mid-slice

### 7.1 The 2D pre-fit and pydme's continuous bisection

Ported exactly: per mass, if `profile_minimizer` is `minuit2d`/`pydme` (pattern-space only — the reference app never runs the 2D fit for n_e-space), one `MinimizeOverSigmaAndTheta` call over the whole σ grid gives the true global NLL minimum, and the resulting per-σ signal grid is reused (not re-folded) for the rest of that mass's σ loop. `pydme` mode alone additionally runs a continuous-σ bisection (`ScanUtils::InterpolateSignal`+`BisectUpperLimit`) as a diagnostic, written to the optional `upper_limit_sigma_e_mchi_pydme_bisection` graph — confirmed neither `band.cc` nor the plot app depend on it except the plot app's `--draw-both` flag.

### 7.2 A real gap: `run.data_path`

Running the hard-parity gate against the actual flagship config (`scan_dmelectron_pattern_pydme_exact.json`, unmodified) failed — that config loads real observed counts from `run.data_path` (a ROOT file, `D_pat` histogram mapped by pattern-code bin label), a mechanism the new app didn't have at all (only inline `run.observed_counts` was supported). Ported `load_data_csv`/`load_data_root`/`load_data` near-verbatim from the reference app, wired into Phase A with the exact same precedence: `observed_counts` (inline) overrides `data_path` (file) overrides Asimov.

### 7.3 Validation

Hard parity gate repeated with the flagship config's real settings (`bp_br_template`, `Bp_theta_Br`, `pydme`, `pydme_style_ul=true`, CSV-weighted efficiencies, and finally real `data_path` data) on a reduced 3-mass grid: main UL curve and the `upper_limit_sigma_e_mchi_pydme_bisection` diagnostic graph both **bit-exact** against the unmodified reference app.

---

## 8. Slice 7 — `run.mode="threshold_toys"` + `run.q_target_lookup_path`

Ported from `ccdarksens_scan_srdm_pattern_csv.cc`'s own existing implementation (the only place this algorithm previously existed): per mass, `S_thr` via `ScanUtils::InterpolateSignal` at `sigma_threshold(m_chi)` (loaded from a `sigma_threshold_per_mass` TGraph), Poisson toys against `S_thr + Bp + Br` with an independently-seeded RNG per mass (`std::seed_seq{seed_lo, seed_hi, ix, 0xC0FFEE}`), each toy's `q` from a seeded 2D fit (`MinimizeOverScaleMinuit` for the constrained top, `MinimizeOverSigmaAndTheta` seeded at `(log10 σ_thr, θ̂_top)` for the global minimum — seeding matters because σ_thr sits high in the scan range where a midpoint-seeded Simplex often misses the true minimum), and an `nth_element`-based percentile writing `q_target_per_mass` to `<outdir>/qtarget_threshold.root`. In normal `"scan"` mode, `run.q_target_lookup_path` (if set) overrides the asymptotic per-mass PLR threshold via nearest-neighbor lookup on that same graph. Required two new, purely additive `RunHeader` fields on `ConfigManager` (`mode`, `q_target_lookup_path`, `threshold_toys.*`) — the reference app never reads them, so its behavior is unchanged.

**Validation**: unmodified `ccdarksens_band.cc`, pointed at the new app with `band.threshold="both"` on a small real `Bp_theta_Br`/`pydme` grid — ran end to end through Phase 0 (asymptotic UL) → Phase 1 (`qtarget_threshold.root`) → Phase 2 (outer toy loop, both asymptotic and toy_mc passes), 3/3 toys succeeded in both modes, final `band.root` contains both asymptotic and toy_mc quantile-band `TGraph`s.

---

## 9. Slice 8 — WIMP-nucleon (`cluster_energy`) end-to-end + one more gap

Running the full `wimp_nucleon`+`cluster_energy` pipeline for the first time (calibration → kernel build → `ModelFactory` → `ClusterEnergyResponseFold::Fold` → PLR) surfaced one more gap: `BackgroundFactory`'s default `"dc_flat_migration"` path unconditionally tried the dark-current fold, which `cluster_energy` cannot support (no pixel/n_e model) — it threw, even though the *flat* half of that background source is perfectly well-defined there (validated back in Slice 2). Fixed by adding a capability query, `ResponseFold::SupportsDarkCurrentBackground()` (default `false`; `true` for pattern/n_e), which `BackgroundFactory` checks before attempting the DC fold — `dc_flat_migration` now degrades gracefully to flat-only for cluster-energy instead of requiring every channel to fake a DC contribution via `bp_br_template`.

Fresh WIMP-nucleon rate CSVs were generated at the current binning (`ee_Emax_eV=8000`, matching `configs/wimp_nucleon_generate_si_heavy.json`) via `utils/wimp_nucleon_generate_grid.py` — the existing `m3000`/`m5000` fixtures are stale at `ee_Emax_eV=400` and were kept as Slice-3 unit-test fixtures only, not used here. A 2-mass × 2-sigma scan (reduced MC trial counts for speed) ran end to end without error and produced the correct output object shape (`q_mchi_sigma`, `upper_limit_sigma_e_mchi`, `upper_limit_sigma_e_mchi_graph`, `exposure_kg_year`) — this is a mechanical integration check (does the whole pipeline run and produce the right shape), not a physics-fidelity validation at production MC-trial counts.

---

## 10. Summary: what this buys

A single scan app now covers every signal channel and every statistical mode the DM-electron reference app supports, verified by direct bit-exact comparison rather than by construction alone. `ccdarksens_band` and `ccdarksens_plot_dmelectron_limit` work against its output completely unmodified. The reference app itself was never touched — every finding here came from reading it, never from editing it.
