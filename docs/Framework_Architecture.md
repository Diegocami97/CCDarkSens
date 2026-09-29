<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Framework_Architecture.md -- Standalone architecture reference for the whole framework: how a config file turns into an exclusion curve, module by module.
-->

# CCDarkSens Framework Architecture

**Document:** `docs/Framework_Architecture.md`
**Author:** Diego Venegas-Vargas
**Status:** Standalone reference, grounded directly in the current source (not the older per-topic docs, some of which predate recent refactors). Written independently of `Beginners_Guide.md`, `DM_Signal_Models_Physics_Reference.md`, and `GenericScanApp_Design.md` — some overlap with those is intentional rather than cross-referenced away.

---

## 0. Scope and how to read this document

This document answers one question: **given a JSON config file, how does the code actually turn it into a σ_UL(mχ) exclusion curve?** It covers the framework generically — all three analysis spaces (`ne`, `pattern`, `cluster_energy`) and all four signal models (`dm_electron`, `dark_photon`, `migdal`, `wimp_nucleon`) — with `ccdarksens_scan_generic` called out throughout as the primary, recommended entry point (`docs/GenericScanApp_Design.md` covers its validation history in more depth; this document covers its architecture).

The pipeline has six layers, each with its own section below, read top to bottom in the order data actually flows through them:

```
JSON config
   |
1. ConfigManager            -- parses the file into typed config structs
   |
2. Detector + ExperimentSetup -- geometry/mass, exposure_kg_year, ROI/binning
   |
3. ModelFactory              -- loads a rate CSV for one (mass, coupling) point into a dR/dE histogram
   |
4. ResponseFactory/ResponseFold -- folds dR/dE through detector physics into observable-space counts S(bin)
   |
5. BackgroundFactory         -- builds the expected background B(bin) (or Bp/Br split) in the same bins
   |
6. ProfileLikelihood + ScanUtils -- likelihood, nuisance-parameter minimization, q-statistic, UL crossing
   |
sigma_UL(mchi) exclusion curve, written to a ROOT file
```

A seventh section covers the orchestration app (`ccdarksens_scan_generic.cc`) that calls all of this once per grid point, and an eighth covers a genuinely important cross-cutting subtlety (the `q_threshold` convention) that caused real confusion earlier in this project because it lives partly in the scan app and partly in the separate plotting app.

---

## 1. Config layer — `ConfigManager`

Files: `include/ccdarksens/io/ConfigManager.hh`, `src/io/ConfigManager.cc`.

`ConfigManager(path)` stores the path; `parse()` reads the file as one `nlohmann::json` object and dispatches each top-level block to its own parser:

| JSON section | Required? | Parser | Result accessor |
|---|---|---|---|
| `run` | **Yes** (`.at("run")`, throws if absent) | `parse_run_` | `cfg.run()` → `RunHeader` |
| `detector` | No | `parse_detector_` | `cfg.detector()` → `Detector` (see §2) |
| `experiment` | No | `parse_experiment_` | `cfg.experiment_cfg()` → `ExperimentConfig` |
| `backgrounds` | No | `parse_backgrounds_` (also fills `TimingJSON` from the nested `backgrounds.timing` sub-object — there is no top-level `"timing"` key) | `cfg.backgrounds()`, `cfg.timing()` |
| `response` | No | `parse_response_` | `cfg.response()` → `ResponseJSON` |
| `model` | No | `parse_model_` | `cfg.model()` → `ModelJSON` |

**There is no schema validation and no rejection of unknown fields anywhere in `ConfigManager`.** Every leaf is read with `nlohmann::json::value(key, default)`, so a misspelled key is silently ignored — no warning, no error. This is a real, load-bearing fact about the framework: several bugs found during the student-example work (a truncated `Emax_eV`, an unwired `backgrounds.lee` block in the joint-channel scan path) were invisible precisely because nothing checks that a config field is actually consumed by the code path it's meant to affect. The only things that *do* throw are a handful of `.at()` calls for fields genuinely required to construct something (`detector.rows/cols/pixel_size_um/thickness_mm/target_element/density_g_cm3`; `experiment.mode/livetime_days/binning.ne_min/binning.ne_max`), plus `experiment.mode` being rejected if it's not exactly `"observed"`, `"asimov"`, or `"toys"`.

One accessor is genuinely unsafe: `cfg.detector()` dereferences a `std::unique_ptr<Detector>` that is only constructed if the JSON has a `"detector"` key at all. A config that omits `detector` entirely makes any call to `cfg.detector()` undefined behavior — there is no null check. Every other accessor (`experiment_cfg()`, `backgrounds()`, `timing()`, `response()`, `model()`) returns a value member with an in-struct default, so it's always safe to call even if the section was absent from the JSON (though downstream code, e.g. `ExperimentSetup`'s constructor, will then throw on the default `livetime_days = 0.0`).

### What each section actually controls

- **`run`** — statistics/run-control, not physics: `cl` (confidence level, default 0.90 — feeds `target_q`, see §6), `test_stat`, `profile_minimizer` (`brent`/`minuit`/`minuit2d`/`pydme`), `background_source` (`dc_flat_migration`/`bp_br_template`), `background_model` (`scale`/`Bp_theta_Br`), `single_bin_likelihood`, `pydme_style_ul`, `smooth_ul_envelope`, `mode` (`scan`/`threshold_toys`), plus the `Bp_theta_Br`-specific constraint/prior knobs (`constrain_prior_strength`, `theta_lo`/`theta_hi`, etc.). Note `run.exposure_kg_year` is a distinct, legacy field from the one the physics pipeline actually uses (`ExperimentSummary::exposure_kg_year`, computed fresh in §2) — setting the `run`-level field has no effect on the scan.
- **`detector`** — geometry and target material; see §2.
- **`experiment`** — exposure inputs (`livetime_days`, `duty_cycle`, `binning.ne_min/ne_max`) and observable definition (`observable_bins`, `roi_bins`, `pattern_roi`); see §2.
- **`backgrounds`** — dark current, flat/d.r.u. background, the low-energy-excess (LEE) term, and pattern-migration efficiency; see §5.
- **`response`** — `analysis_space` (the single most important field in this section — selects which of the three response backends gets built, see §4), plus each backend's own sub-config (`efficiency_mc`/`pattern_mc`, `cluster_fit_mc`, `charge_ionization`, `pattern_classifier`, `channels[]` for the joint-likelihood case).
- **`model`** — `type` (`dm_electron`/`dark_photon`/`migdal`/`wimp_nucleon`), `rates_dir`/`filename_template` (where to find rate CSVs and how to name them), `Emin_eV`/`Emax_eV`/`nbins` (the signal-spectrum histogram binning — see §3 for why getting this wrong silently drops signal), and the mass/coupling `grid` block.

### `model.type` and `response.analysis_space` are independent knobs — nothing enforces a valid pairing

This is worth stating plainly because it is not obvious from reading either factory in isolation. `ModelFactory::MakeSignalSpectrumE` dispatches purely on `model.type`; it has no awareness of `response.analysis_space`. `ResponseFactory::MakeResponseFold` dispatches purely on `response.analysis_space`; it has no awareness of `model.type`. Nothing in `ConfigManager`, `ExperimentSetup`, `ModelFactory`, or `ResponseFactory` validates that, e.g., `model.type=="wimp_nucleon"` is paired with `response.analysis_space=="cluster_energy"` (the intended pairing, per scattered comments elsewhere) rather than left at the default `"pattern"`. A mismatched config would not be rejected up front — it would either fail later for an unrelated reason or silently produce a signal folded through the wrong response path. This has not been observed to actually happen in practice (every config in the repo pairs them correctly by convention), but it is a real gap, not a hypothetical one, and worth knowing before writing a new config from scratch.

---

## 2. Detector and experiment setup

Files: `include/ccdarksens/detector/Detector.hh`/`.cc`, `include/ccdarksens/experiment/ExperimentSetup.hh`/`.cc`.

**`Detector`** validates its geometry/material at construction (throws `std::invalid_argument` for non-positive rows/cols/pixel size/thickness/density, out-of-range `active_fraction`, or an empty target-element string) and exposes exactly one derived quantity beyond plain accessors: `mass_kg()`. If the config gives an explicit `detector.mass_kg`, that value wins; otherwise `mass_kg()` computes it from geometry: `area_cm2 = (cols × pitch_cm) × (rows × pitch_cm)`, `volume_cm3 = area_cm2 × thickness_cm × active_fraction`, `mass_g = density_g_cm3 × volume_cm3`.

**`ExperimentSetup`** is where `exposure_kg_year` — the single number that scales every folded signal and background spectrum — actually gets computed, and it happens in exactly one place:

```
exposure_kg_year = livetime_days * duty_cycle * detector_mass_kg / 365.25
```

The constructor throws if `livetime_days <= 0`, `duty_cycle` is outside `(0, 1]`, or `binning.ne_max <= binning.ne_min`. `prepare_summary()` returns an `ExperimentSummary` — the struct that gets threaded through the rest of the pipeline (`ResponseFactory`, `BackgroundFactory`, the scan loop all take this, not the raw `ExperimentConfig`):

```cpp
struct ExperimentSummary {
  double exposure_kg_year;
  BinningNE binning;
  std::vector<int> roi_bins;      // sorted, de-duplicated n_e ROI
  std::vector<int> pattern_roi;   // pattern IDs, pattern space only
  std::string observable_bins;    // "n_e" or "pattern"
  uint64_t rng_seed_used;
  std::string mode_string;
};
```

`roi_bins` gets normalized here (if the config didn't supply one, it's built as the full inclusive `[ne_min, ne_max]` range; if it did, it's sorted and de-duplicated). `pattern_roi` and `observable_bins` are passed through **verbatim, with no validation** — in particular, `ExperimentSetup` does not check that `pattern_roi` is non-empty when `observable_bins=="pattern"`. That one real validation (`std::runtime_error` if `observable_bins=="pattern"` but `experiment.pattern_roi` is empty) lives downstream, in `ResponseFactory::MakePatternOrNeSpaceFold`.

---

## 3. Model layer — `ModelFactory`

Files: `include/ccdarksens/model/ModelFactory.hh`, `src/model/ModelFactory.cc`, one model class per `model.type` (`DMElectronModel`, `DarkPhotonModel`, `MigdalModel`, `WimpNucleonModel`), `include/ccdarksens/io/RateTable.hh`/`.cc`.

`MakeSignalSpectrumE(mj, mass_val, coupling_str, &ok)` is a short, pure dispatcher (~80 lines) that copies the relevant `ModelJSON` fields into a model-specific config struct, calls that model's `Configure()` + `MakeSpectrum_E()`, and returns a `TH1D`. The four model classes are thin, near-identical wrappers around the same two-step recipe — they are **not** where the physics lives (that's in the Python rate-generation scripts under `python/ccdarkphys/`, which produce the CSVs on disk ahead of time); their only job is filename resolution and CSV loading:

1. **`ResolvePath_()`** — token-substitutes `{material}`, `{mediator}`, `{mchi_MeV}` (formatted to exactly 6 decimals), `{sigma_e_cm2}`/whatever coupling field name applies into `model.filename_template`, prepended with `model.rates_dir`.
2. **`RateTable::LoadCSV(path)`** — reads a two-column `(E_eV, dR/dE)` CSV (comma or whitespace separated, `#`-comments and blank lines skipped, first data line treated as a header and discarded). Returns `false` (not an exception) if the file can't be opened or the two columns end up different lengths — `MakeSignalSpectrumE`'s `ok` output parameter is exactly this boolean, and the scan loop treats a load failure as "use an all-zero signal" rather than aborting (with a warning at `verbosity>=1`).
3. **`RateTable::MakeTH1D(name, Emin_eV, Emax_eV, nbins)`** — builds a `TH1D(name, name, nbins, Emin_eV, Emax_eV)` and, **for each of the `nbins` bins of that fixed range**, evaluates the loaded table by linear interpolation at the bin *center* (not integrated over the bin width — the code's own comment flags this as bin-center sampling, "will be replaced with proper trapezoid integration later"). Below the table's first tabulated energy, every bin gets the table's first value; above the table's last tabulated energy, the bin gets exactly `0`.

**This is the exact mechanism behind a real bug found this session**: if `model.Emax_eV` is set smaller than where a rate CSV's actual physical content sits (a dark-photon absorption line at `mA' = 40 eV` with `model.Emax_eV = 20`, say), the histogram's bins simply never sample the table at 40 eV — that whole slice of the mass grid contributes a genuinely all-zero signal, no warning, no error, no different in any log line from a rate table that was correctly zero there. `model.Emin_eV`/`Emax_eV`/`nbins` are not descriptive metadata; they define the literal energy range the fold ever sees.

The `model.grid` block (mass/coupling axis definitions) is expanded by `ccdarksens::utils::ExpandAxis` (`include/ccdarksens/utils/AppUtils.hh`), which the scan apps call once per axis via `ExpandAxisFromGrid`, itself just trying candidate key names in order (e.g. `{"mA_eV", "mchi_MeV"}` for `dark_photon`, `{"mchi_MeV"}` otherwise) until one exists in the JSON. `ExpandAxis` accepts a plain array, `{"values":[...]}`, `{"linspace":{...}}`, or `{"logspace":{...}}`; the expanded list is **not sorted or de-duplicated** — an out-of-order or duplicate `values` list will scan exactly as given. `model.grid.format.mchi`/`format.sigma` (e.g. `.6f`, `.1e`) are the printf-style precision hints used to turn a numeric grid value into the exact string substituted into the filename template — get this wrong (or out of sync with whatever generated the CSVs) and every file lookup for that axis fails silently in the same "load returns false, fold with zero" way described above.

---

## 4. Detector response layer — `ResponseFactory` / `ResponseFold`

Files: `include/ccdarksens/response/ResponseFold.hh` (interface), `include/ccdarksens/response/ResponseFactory.hh`/`.cc` (dispatcher), `NeSpaceResponseFold.hh`/`.cc`, `PatternResponseFold.hh`/`.cc`, `ClusterEnergyResponseFold.hh`/`.cc` (the three concrete backends), plus the physics engines each one leans on: `ChargeIonization`, `EfficiencyMC`, `PatternClassifier`, `PatternRates.hh` for pattern/n_e space, and `NoiseTailCalibrator`/`ClusterFitMC` for cluster-energy space.

This is the most complex layer in the framework, and the one where a genuinely important design decision is easy to miss: **the fold classes do not reuse the older `DetectorResponsePipeline` class** (still compiled, still referenced by several of the pre-`scan_generic` single-purpose apps) **even though it looks like it should do the same job.** The reason, stated explicitly in `ResponseFold.hh`'s own design comment, is that the reference app's actual signal path already bypassed `DetectorResponsePipeline::Apply()` — calling `ChargeIonization::FoldToNe` and `FoldNeToPatternRates` directly instead, with a comment reading *"Use S_true, not S_obs — folding with S_obs would double-count."* `ResponseFold` generalizes what the reference app *actually does*, not what the older pipeline class implements. Treat `DetectorResponsePipeline` as legacy: still used by older PCD/diffusion cross-check apps, not part of the current production fold path, and not something new work should build on.

### 4.1 The interface (`ResponseFold`)

A four-method abstract base:

```cpp
virtual std::vector<double> Fold(const TH1D& dRdE_density, double exposure_kg_year) const = 0;
virtual std::size_t NumBins() const = 0;
virtual std::vector<std::string> BinLabels() const = 0;
virtual std::vector<double> FoldNeBackground(const TH1D& B_ne_asimov) const;      // default: throws
virtual bool SupportsDarkCurrentBackground() const { return false; }              // default: false
virtual std::vector<double> ErecoEdgesEV() const;                                 // default: throws
```

`Fold()` is the one contract every backend must implement: take a `dR/dE` density histogram [events/(kg·year·eV)] and an exposure, return expected counts per output bin. Background is deliberately *not* a separate abstraction — it's the same `Fold()` call applied to a flat dummy spectrum in the two spaces that have an energy axis. `FoldNeBackground`/`SupportsDarkCurrentBackground` exist because dark current is synthesized directly in n_e space (no `dR/dE` spectrum involved) and only makes sense for `NeSpaceResponseFold`/`PatternResponseFold` — both override the default; `ClusterEnergyResponseFold` does not, so a `cluster_energy` config with `background_source=dc_flat_migration` degrades gracefully to flat-only rather than throwing (`BackgroundFactory` checks `SupportsDarkCurrentBackground()` before ever calling `FoldNeBackground`). `ErecoEdgesEV()` exists purely so the background layer can evaluate a background-efficiency curve or the LEE exponential against real reconstructed-energy bin edges — only `ClusterEnergyResponseFold` has bin edges that mean anything, so it's the only override.

### 4.2 The dispatcher (`ResponseFactory::MakeResponseFold`)

A ~14-line function, purely a string switch on `response.analysis_space`:

- `"cluster_energy"` → `MakeClusterEnergyFoldFromJSON(response.cluster_fit_mc)`
- `"pcd"` → throws unconditionally, not implemented in this path ("the reference app's signal path bypasses its own PCD wiring too")
- anything else (including the default, `"pattern"`) → `MakePatternOrNeSpaceFold(cfg, summary, config_path)`, which internally branches again on `summary.observable_bins` (`"pattern"` vs `"n_e"`)

The return value, `ResponseFactoryResult`, bundles the constructed `fold`, the `ChargeIonization` instance used to build it (`ion` — needed later so the background layer can fold a flat spectrum through the *same* ionization table the signal used; null for `cluster_energy`, which has no such table), and `ne_min_bkg`/`ne_max` (pattern space includes n_e=0 for dark-current purposes; n_e space does not; both are unused for cluster_energy).

### 4.3 `n_e` space — `NeSpaceResponseFold` + `ChargeIonization`

`ChargeIonization` loads a table-driven P(n_e | E) CSV (columns `E_eV, P1, P2, ..., Pk`, linearly interpolated and clamped to `[0,1]`) and provides `FoldToNe(dRdE, exposure_kg_year, ne_min, ne_max)`: for each `dR/dE` bin, `counts = rate × exposure × bin_width_eV`, distributed across n_e = 1..k by the table's interpolated probability at that bin's center energy. This is the true, un-thresholded n_e spectrum `S_true(ne)`.

`NeSpaceResponseFold::Fold()` is then just `S_true(ne) × clamp(eps_ne[ne], 0, 1)` for each `ne` in the ROI — and `eps_ne` is **not** a pure charge-quantization efficiency; it comes from `EfficiencyMC`, the same Monte-Carlo/CSV pattern-efficiency machinery pattern space uses, just configured with single-pixel `{ne}` labels instead of multi-pixel pattern codes. So n_e-space and pattern-space share the entire efficiency-construction machinery described in §4.4 below — they differ only in what they do with the resulting per-(pattern-or-ne) efficiency map, not in how it's built.

### 4.4 Pattern space — `PatternResponseFold` + `EfficiencyMC` + `PatternClassifier`

`PatternResponseFold::Fold()` does the same `FoldToNe` first step, then `FoldNeToPatternRates(S_true, ne_min, ne_max, pattern_roi, pattern_eff_map)`:

```
for pattern_id in pattern_roi:
  rate = sum over ne in [ne_min, ne_max] of  S_true(ne) * pattern_eff_map[{pattern_id, ne}]   // missing entry = 0
```

(This is exactly the function whose `n_e >= 10` special case was the source of a real bug fixed this session — see `docs/ClusterFitMC_Design.md` and the absorption/DM-electron-intermediate student docs for the full story; the current code has no special case at all, a missing `(pattern, ne)` table entry is always zero acceptance regardless of `ne`.)

The `pattern_eff_map` — a `(pattern_code, ne) → efficiency` lookup — is built one of two ways, and which one happens is decided entirely by whether `response.efficiency_mc.efficiency_csv` is set:

- **CSV supplied** → `EfficiencyMC::PrecomputeEpsilonWithPatternEff` takes the **fast path**: it reads the CSV values directly and never runs any Monte Carlo. `PatternClassifier` is constructed (unconditionally, by `ResponseFactory`) but its scanning methods are never called in this path.
- **No CSV** → `EfficiencyMC::PrecomputeEpsilon` calls `BuildPatternTable`, which actually runs the pattern-efficiency Monte Carlo: simulate charge clouds (`PixelSimulator` + `ChargeTransport::SampleCloudXY`) at each n_e, classify the resulting pixel pattern (`PatternClassifier::ScanRow`, or `ScanAllImage2DWithIsolation` for the 2D-image diagnostic path gated by `use_2d_image_efficiency`), and tabulate `P(label | ne)` over many trials (`efficiency_mc.ne_trials`, default 50,000 per n_e — genuinely expensive, which is why every production config in this repo ships a pre-computed `efficiency_csv` instead of regenerating it at scan time).

**`PatternClassifier` therefore only ever matters when *building* an efficiency table from scratch — it has no role in the per-grid-point `Fold()` call during an actual scan.** This is worth knowing before trying to debug a scan-time pattern issue by looking at `PatternClassifier`: if the config has an `efficiency_csv`, that class is never exercised at all during the scan.

### 4.5 Cluster-energy space — `ClusterEnergyResponseFold`

Built by `MakeClusterEnergyFoldFromJSON(cluster_fit_mc_json)`: calibrate a pure-noise ΔLL cut (`NoiseTailCalibrator::Run()`), build a `ChargeTransport` and a `ClusterFitMC`, call `ClusterFitMC::BuildKernel(Etrue_grid, Ereco_edges)` to get a `KernelMatrix` (this step is the expensive one — Monte Carlo per point on a sparse `Etrue` grid, minutes per channel), and wrap the kernel in a `ClusterEnergyResponseFold` whose `Fold()` is just `FoldEtrueToErecoRates(dRdE, exposure, kernel)`. The kernel-building math (the sparse-grid quadrature bug found and fixed, the background-efficiency-curve distinction, the noise-tail calibration itself) is covered in depth in `docs/ClusterFitMC_Design.md` and deliberately not repeated here.

The one thing worth stating precisely in *this* document: **the single-channel vs. joint-channel (`response.channels[]`) split lives almost entirely in the orchestration app, not in the response classes.** `ResponseFactory`/`ClusterEnergyResponseFold` know nothing about "joint" — they only ever build one fold from one `ClusterFitMCJSON` block. `ResponseFactory.hh` exposes `MakeClusterEnergyFoldFromJSON` as a public function precisely so `apps/ccdarksens_scan_generic.cc::RunJointChannelScan` can call it once per entry in `response.channels[]`, building one independent `ResponseFold` *and one independent `ProfileLikelihood`* per channel, folding the same physical signal spectrum through each, profiling each channel's NLL separately, and summing the per-channel NLLs before the usual `MonotonizeQ`/`UlFromQMonoCrossing` crossing search (exact, not approximate, because each channel's background nuisance parameter appears only in its own NLL term — see §6). This is also where the real, previously-hidden bug lived that the WIMP-nucleon joint-channel LEE background was never wired in: `RunJointChannelScan` builds each channel's background directly from `MakeClusterEnergyFlatBackground`/`...WithEfficiency` and never checks `backgrounds.lee` at all, unlike the single-channel path (§5) which does.

---

## 5. Background layer — `BackgroundFactory`

Files: `include/ccdarksens/response/BackgroundFactory.hh`, `src/response/BackgroundFactory.cc`. (A separate, older `include/ccdarksens/backgrounds/BackgroundBuilder.hh`/`.cc` exists specifically for the Poisson dark-current piece — `BackgroundFactory` constructs and calls into it, it is not a competing/duplicate path.)

`MakeBackground(cfg, summary, fold, response)` is the single-channel public entry point, dispatching on `run.background_source`:

- **`"bp_br_template"`** — `MakeBpBrTemplateBackground`: takes `run.background_Bp`/`background_Br` directly from the config (throws if their length doesn't match the fold's bin count), sets `B_pat = Bp + Br`. This is the path used to reproduce a real published analysis's own pre-computed background template exactly, bypassing any physics reconstruction of it.
- **`"dc_flat_migration"`** (the default) — `MakeDcFlatMigrationBackground`, which sums up to three independently-constructed pieces, each already folded into the response's own output bins:
  1. **Dark current** — only if `fold.SupportsDarkCurrentBackground()` (pattern/n_e space only; degrades to zero for cluster_energy rather than throwing). `BackgroundBuilder` builds a Poisson dark-current spectrum in raw n_e space from `backgrounds.dark_current.lambda_e_per_pix_per_year`/`norm_scale` plus the detector geometry and timing, and `fold.FoldNeBackground(...)` folds it through the same pattern/n_e efficiency the signal used.
  2. **Flat background** — a flat d.r.u. spectrum over `model.Emin_eV`/`Emax_eV`/`nbins` (note: it reuses the **signal model's** binning fields, not `backgrounds.flat_background`'s own `Emin_eV`/`Emax_eV`/`nbins`, which exist in the parsed config but are not what this path actually uses), folded through `fold.Fold(...)` exactly like a signal spectrum — **unless** `response.background_efficiency_csv` is set, in which case `MakeClusterEnergyFlatBackgroundWithEfficiency` is used instead, weighting the flat rate by a background-specific detection-efficiency curve (bin-averaged via 9-point Simpson's rule, not point-sampled — the digitized curve used for the WIMP-nucleon reproduction rises by tens of percent within a single bin width near threshold, so point-sampling was a real, since-fixed error).
  3. **Low-energy excess (LEE)** — only if `backgrounds.lee.enabled`. `MakeClusterEnergyLEEBackground` integrates the exponential `dR/dE = (1/ε)·exp(−E/ε)` exactly (closed form, no quadrature) over each reconstructed-energy bin. This path calls `fold.ErecoEdgesEV()`, so it is meaningful for `cluster_energy` only (per its own header comment, "WIMP-nucleon channel only") — and, per §4.5, it is only ever reached from the *single-channel* code path; the joint-channel path never checks `has_lee_bkg` at all.

Every piece returns a `BackgroundResult { B_pat, Bp, Br }`; for every `dc_flat_migration` sub-path, **`Bp` is always all-zero and `Br` is always the full total** — meaning `background_model="Bp_theta_Br"` on top of `background_source="dc_flat_migration"` reduces to a single global scale on the entire background (`B = 0 + theta × B_pat`), not a genuine two-component split. A real `Bp`/`Br` split only exists when `background_source="bp_br_template"` supplies both vectors directly from the config. This distinction — `BackgroundFactory` only ever *produces* `Bp`/`Br`; the stats layer (§6) is what actually treats `scale` vs `Bp_theta_Br` differently — is worth being explicit about, since it's easy to assume incorrectly that `dc_flat_migration` has some inherent primary/residual background split baked in.

---

## 6. Stats layer — `ProfileLikelihood` + `ScanUtils`

Files: `include/ccdarksens/stats/ProfileLikelihood.hh`/`src/stats/ProfileLikelihood.cc`, `include/ccdarksens/scan/ScanUtils.hh`/`src/scan/ScanUtils.cc`. This is the layer with the thinnest pre-existing documentation, so it's covered in the most mechanical detail here.

### 6.1 The likelihood

`ProfileLikelihood` holds observed data (`SetData`) and a background in one of two modes:

- **Scale mode** (`SetBTemplate(B_template)`): the model is `mu_i = S_i + scale * B_template_i`.
- **Bp/theta/Br mode** (`SetBpBr(Bp, Br)`, "pydme mode"): the model is `mu_i = S_i + Bp_i + theta * Br_i`, with `theta` the single global nuisance parameter (not `scale`).

`NLL(S, param)` computes the Poisson negative log-likelihood, dropping the `n_i!`-dependent (parameter-independent) term:

```
NLL = sum_i [ mu_i - n_i * safe_log(mu_i) ]
```

(`safe_log` — `include/ccdarksens/stats/StatsUtils.hh` — guards against `log(0)`; a `mu_i <= 0` with `n_i > 0` is penalized with a large constant, `1e9`, rather than producing `NaN` or `inf`.) In `Bp_theta_Br` mode, an optional prior/constraint term on `theta` can be added on top (`constrain_prior_strength`, with pydme-matching sign/tau-weighting/multi-bin options) — this is what lets a `Bp_theta_Br` config reproduce a real published analysis's own constrained-background fit, not just a scale factor.

Two minimizers are available for finding `param_hat` at fixed `S`: `MinimizeOverScale` (Brent, bounded) and `MinimizeOverScaleMinuit` (ROOT Minuit2, falling back to Brent if unavailable). `run.profile_minimizer` selects which one the scan loop actually calls (`"brent"` → `MinimizeOverScale`; `"minuit"`/`"minuit2d"`/`"pydme"` → `MinimizeOverScaleMinuit`, with pattern-space `"minuit"` additionally cross-checking against Brent and keeping whichever gives the lower NLL). `"pydme"`/`"minuit2d"` additionally trigger a genuine 2D joint fit over `(log10(sigma), theta)` in pattern space (`MinimizeOverSigmaAndTheta`) — this is not a separate code path elsewhere; it's the same `ProfileLikelihood` class, called differently by the scan app depending on `run.profile_minimizer` and `use_pattern_bins`.

### 6.2 From NLL to an upper limit — `ScanUtils`

Given `nll_values[k]` (one NLL per coupling grid point, at fixed mass) and the global minimum `nll_min` across that row:

- **`MonotonizeQ(nll_values, nll_min)`**: `q[k] = max(0, 2*(nll_values[k] - nll_min))`, then a running maximum so `q(sigma)` never decreases. This is necessary because the raw `q(sigma)` curve from independent per-point minimizations is not guaranteed monotonic (numerical minimizer noise, or — as found this session for dark-photon absorption — genuine but unphysical non-monotonicity from a since-fixed detector-response bug); a one-sided upper limit search requires monotonicity to have a well-defined single crossing.
- **`UlFromQMonoCrossing(sigma_list, q_mono, target_q)`**: finds the first grid interval where `q_mono` crosses `target_q` and interpolates **log-linearly in sigma** (linear in `log10(sigma)`, linear in `q`) between the two bracketing grid points. If `q_mono` never reaches `target_q` anywhere (too-good sensitivity, grid doesn't extend far enough), returns the largest grid sigma; if `q_mono` is already above target everywhere, returns the smallest.
- **`BisectUpperLimit`** (used only in `pydme_style_ul` mode): a genuine bracket-and-bisect search directly in `q_mu(log10 sigma)` space rather than interpolating a fixed grid — used to get a continuous, grid-resolution-independent crossing for the pydme-style curve.

### 6.3 The `q_threshold` — computed correctly by the scan, hardcoded differently by the plotter

`target_q` — the value `q_mono` must cross to define the upper limit — is computed **the same way, every time, in both the single-channel and joint-channel scan paths**, directly from `run.cl`:

```cpp
const double target_q = std::pow(TMath::NormQuantile(run.cl), 2);
```

This is the correct asymptotic one-sided formula (Cowan et al. 2011): for `cl=0.90`, `NormQuantile(0.90) = 1.2816`, giving `target_q ≈ 1.6424` — the `1.642374415149816` constant that shows up throughout the student-example docs and configs. **This is computed dynamically from `run.cl`, not hardcoded, and it is identical regardless of `profile_minimizer` or analysis space.** The scan's own stored `upper_limit_sigma_e_mchi` graph (written into the output ROOT file, see §7) always uses this value.

The separate plotting app, `apps/ccdarksens_plot_limit.cc`, reimplements its own crossing search directly on a saved `q(mchi,sigma)` 2D histogram (`--from-qhist` mode; essentially the same log-linear-interpolation logic as `UlFromQMonoCrossing`, but operating on a `TH2D` read back from disk rather than a live `nll_values` vector) — and its own `q_thr` **defaults to a hardcoded `2.71`** (`apps/ccdarksens_plot_limit.cc`: `double q_thr = 2.71; // default ~90% CL, 1 dof` — the classic two-sided-derived Δχ²=2.71 chi-square convention, a *different* number from the scan's own asymptotic one-sided value) unless an explicit `q_threshold` argument overrides it. Every student-example plot command in this repo passes `1.642374415149816` explicitly to match the pydme-style channels' own convention — **except** the WIMP-nucleon channel, where the correct value to use is actually the plotting app's own default of `2.71` (that channel uses `profile_minimizer: brent`, not `pydme`, and a real, extensively-traced investigation this session confirmed the historical, better-tracking WIMP-nucleon reproduction plots were made without overriding the default at all). This is a genuine, easy-to-miss subtlety: **the "right" `q_threshold` argument to pass to the plotting app depends on which channel's convention you're reproducing, and it is not the same number the scan app used internally to compute its own stored upper-limit graph.**

---

## 7. Orchestration — `ccdarksens_scan_generic.cc`

File: `apps/ccdarksens_scan_generic.cc` (952 lines). This is the app every layer above feeds into, and the one meant for all new work (the older, per-channel/per-space apps like `ccdarksens_scan_dmelectron_pattern.cc` are kept as frozen, unmodified parity references, not as the recommended path). Its own header comment describes it precisely: built on the same already-validated primitives above, generalized over `ResponseFold` so one loop works for every model × analysis-space combination.

`main()`:

1. **Setup** — parse the config; if `response.channels` is non-empty, hand off entirely to `RunJointChannelScan` (§4.5/§7's joint variant) and return from there instead.
2. **Phase A** (built once, before the grid loop) — `ExperimentSetup::prepare_summary()` → `MakeResponseFold` (§4) → `MakeBackground` (§5) → resolve observed data with precedence `run.observed_counts` (inline array) > `run.data_path` (CSV or ROOT `D_pat` histogram) > Asimov (`data = background`) → construct the `ProfileLikelihood` in scale or `Bp_theta_Br` mode (§6.1).
3. **Grid** — expand the mass/coupling axes from `model.grid` (§3), with the axis-name pair chosen by `model.type`; compute `target_q` (§6.3); optionally load a per-mass `q_target_lookup_path` override (used by `ccdarksens_band`'s toy-MC calibration).
4. **Threshold-toys early exit** — if `run.mode=="threshold_toys"`, run a toy-MC calibration of the per-mass `q` threshold instead of the normal scan, write `qtarget_threshold.root`, and return — this is a distinct mode from the normal scan, not a post-processing step on top of it.
5. **Phase B** (the grid loop) — for every mass: get the null-hypothesis NLL (`S=0`); optionally run the pattern-space-only 2D `(log10 sigma, theta)` fit (`minuit2d`/`pydme`); for every coupling value, build the signal spectrum (`MakeSignalSpectrumE`, §3), fold it (`fold.Fold(...)`, §4), and minimize the NLL (§6.1); take the global minimum NLL across the row (or the 2D fit's minimum, if `pydme_style_ul` requests it); `MonotonizeQ` + `UlFromQMonoCrossing` (§6.2) to get that mass's upper limit; optionally also run the continuous `BisectUpperLimit` pydme-style crossing.
6. **Phase C** — optional `smooth_ul_envelope` post-processing (a ±5-point running-minimum envelope over the mass axis — note from §6.3 that this has **no effect at all** on any plot made with `--from-qhist`, since that mode rebuilds the curve directly from the raw `q(mchi,sigma)` histogram, bypassing the stored, possibly-smoothed `upper_limit_sigma_e_mchi` graph entirely); write the output ROOT file.

**Output ROOT file** (`<outdir>/scan_generic.root`), object names chosen to match the frozen reference apps exactly so `ccdarksens_band` and `ccdarksens_plot_dmelectron_limit` work against either app unmodified:

| Object | Type | Contents |
|---|---|---|
| `q_mchi_sigma` | `TH2D` | the raw, monotonized `q` value at every (mass, coupling) grid point — this is what `--from-qhist` rebuilds a curve from |
| `upper_limit_sigma_e_mchi` | `TH1D` | the scan's own stored upper limit per mass, using the scan's own `target_q` (§6.3) |
| `upper_limit_sigma_e_mchi_graph` | `TGraph` | same data as the histogram above, as a graph (required by both downstream consumer apps) |
| `upper_limit_sigma_e_mchi_pydme_bisection` | `TGraph` (optional) | the continuous `BisectUpperLimit` curve, only written if `pydme_mode` produced any points |
| `exposure_kg_year` | `TParameter<double>` | the exposure computed in §2 |

`RunJointChannelScan` (the branch taken when `response.channels` is non-empty) mirrors this same Phase A/B/C structure per channel independently, then sums NLLs across channels before the same `MonotonizeQ`/crossing step (§4.5) — it currently supports only `background_model=scale` (not `Bp_theta_Br`), only `profile_minimizer` `brent`/`minuit`, Asimov data only, and does not support `run.mode=threshold_toys`; these are explicit, checked-and-thrown scope limits in the code, not silent gaps.

---

## 8. What to actually read, by task

- **Writing a new example config for an existing channel/space** — you mostly need §1 (which fields exist and are silently ignored if misspelled) and §3 (the `Emin_eV`/`Emax_eV` truncation trap).
- **Debugging "my scan gives zero signal everywhere"** — check §3 first (rate file not found, or found but outside the declared energy range), then §4.4 if pattern-space (missing `pattern_eff_map` entries).
- **Adding a new signal model** — §3 only; `ModelFactory` is deliberately the single place to add a new `model.type` branch.
- **Adding a new analysis space / detector-response backend** — §4; implement the four-method `ResponseFold` interface, add one branch to `ResponseFactory::MakeResponseFold`.
- **Understanding why two plots of the "same" scan look different** — §6.3, almost always. Check what `q_threshold` argument was actually passed to `ccdarksens_plot_limit`, and whether `--from-qhist` was used (which ignores `run.smooth_ul_envelope` entirely, per §7).
- **The WIMP-nucleon `cluster_energy` kernel-building math itself** (noise-tail calibration, the sparse-grid quadrature fix, the background-efficiency curve) — not here; go to `docs/ClusterFitMC_Design.md`, which covers it in the depth §4.5 deliberately doesn't repeat.
