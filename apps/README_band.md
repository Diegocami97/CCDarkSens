<!--
  CCDarkSens — README_band
  User-facing quick-start for ccdarksens_band: build, phased runs, config keys, and how band ROOT output feeds the limit plotter.

  Author: Diego Venegas-Vargas
-->

# `ccdarksens_band` — generic sensitivity-band orchestrator

A standalone C++ executable (`build/ccdarksens_band`) that drives any
CCDarkSens scan binary as a subprocess to produce a frequentist *expected*
sensitivity envelope (the "Brazilian band"): median expected and
± 1σ / ± 2σ quantiles of the σ_UL curve under the background-only
hypothesis.

The full design and rationale (toy-sampling source, threshold modes,
contracts between the orchestrator and the scan binaries) lives in
[`docs/Plan_sensitivity_band_app.md`](../docs/Plan_sensitivity_band_app.md).
This README is the **practical** quick-start.

---

## What it does

For each scan-binary run the orchestrator orchestrates these phases:

| Phase | When | Run | Purpose |
|------:|:------|:----|:--------|
| 0 | `threshold ∈ {toy_mc, both}` | scan binary, **Asimov** data (`observed_counts` stripped) | Get `σ_threshold(m_χ)` (= asymptotic UL on Asimov) used as the test point for Phase 1. |
| 1 | `threshold ∈ {toy_mc, both}` | scan binary in `run.mode="threshold_toys"` | Sample N sub-toys per mass under H_μ at `(m_χ, σ_threshold(m_χ))`, compute `q_μ`, store the configurable percentile (default 90 %) as `q_target_per_mass`. |
| 2 | always | scan binary K times in parallel, fed Poisson(Bp+Br) toy data | Harvest one σ_UL curve per toy. After all toys: pivot to per-mass arrays and take quantiles. |

Outer toys are drawn from `Poisson(Bp_i + Br_i)` (B-only, θ = 1 nominal).
Observed counts **never** influence the band — they are plotted on top of it
by `ccdarksens_plot_dmelectron_limit` from a separate nominal scan run.

---

## Build

```bash
cmake -S . -B build
cmake --build build -j
```

This produces `build/ccdarksens_band` alongside the scan binaries.

---

## Run

The orchestrator reads the **same JSON config** as the target scan binary,
plus a top-level `band` block:

```bash
build/ccdarksens_band configs/scan_srdm_pattern_csv_per_mass_pydme_match.json
```

**Phases only (C++ orchestrator — same subprocesses as a full run, for timing or
incremental work):** all minimization is inside `scan_binary`; the orchestrator
only writes JSON, calls `std::system`, and aggregates ROOT. Use:

```bash
build/ccdarksens_band configs/scan_srdm_pattern_csv_per_mass_pydme_match.json --stop-after phase0
build/ccdarksens_band configs/scan_srdm_pattern_csv_per_mass_pydme_match.json --stop-after phase1
```

`phase0` runs the Asimov scan and writes `_toys/sigma_threshold.root`.
`phase0` / `phase1` require `band.threshold` ∈ {`toy_mc`, `both`}.
`phase1` runs Phase 0 **and** Phase 1, then exits (writes `_toys/phase1/qtarget_threshold.root`).
Omit the flag (or pass `--stop-after phase2`) for the full outer toy loop and `band.root`.

Output:

```
outputs/band_srdm_pattern_csv_per_mass_pydme_match/band.root
outputs/band_srdm_pattern_csv_per_mass_pydme_match/_toys/phase0/...   (kept)
outputs/band_srdm_pattern_csv_per_mass_pydme_match/_toys/phase1/...   (kept)
outputs/band_srdm_pattern_csv_per_mass_pydme_match/_toys/toy_*/       (deleted unless keep_per_toy_outputs=true)
```

---

## Config: the `band` block

Add this top-level block alongside `run`, `detector`, etc. Required keys are
in **bold**:

```jsonc
"band": {
  "scan_binary": "build/ccdarksens_scan_srdm_pattern_csv",   // REQUIRED
  "n_toys": 200,                                             // REQUIRED
  "outdir": "outputs/band_<label>",                          // REQUIRED

  "rng_seed": 7777,
  "n_workers": 0,                  // 0 → hardware_concurrency()
  "threshold": "both",             // "asymptotic" | "toy_mc" | "both"
  "n_threshold_toys": 5000,
  "threshold_percentile": 0.90,
  "threshold_seed": 8888,
  "background_bp_key": "run.background_Bp",
  "background_br_key": "run.background_Br",
  "ul_graph_name":     "upper_limit_sigma_e_mchi_graph",
  "save_per_toy_curves": false,
  "keep_per_toy_outputs": false,
  "abort_on_failure": false,
  "min_success_fraction": 0.95,
  "tmp_dir": "<outdir>/_toys",     // override scratch directory
  "quantiles": { "median": 0.50, "low_1sigma": 0.16, "high_1sigma": 0.84, "low_2sigma": 0.025, "high_2sigma": 0.975 }
}
```

### Parameter reference

Every key under the top-level `band` block, in the order they appear above:

#### Required

| Key | Type | Meaning |
|:----|:-----|:--------|
| `scan_binary` | string | Path (relative to the working directory or absolute) to the scan executable that satisfies the §2 contract in `Plan_sensitivity_band_app.md` — i.e. it accepts the same JSON config you pass to `ccdarksens_band` and writes a TGraph named `ul_graph_name` (see below) into a `*.root` in `run.outdir`. The orchestrator will invoke this binary three times (Phase 0, Phase 1, K outer toys) for `threshold ∈ {toy_mc, both}`, or just K times for `threshold = asymptotic`. |
| `n_toys` | int ≥ 1 | K, the number of outer-loop toy pseudo-experiments used to build the Brazilian band. Each toy is a single Poisson realization of the background-only data and a full scan-binary call. Statistical noise on the band quantiles scales like `1/sqrt(K)`. Typical values: 20 for a smoke test, 200 for a publishable plot, 1000+ for tail-quantile precision. |
| `outdir` | string | Where `band.root` is written. Created if missing. Toy scratch files (Phase 0/1 outputs and per-toy directories unless `keep_per_toy_outputs=true`) live in `outdir/_toys` by default — see `tmp_dir`. |

#### RNG / parallelism

| Key | Type | Default | Meaning |
|:----|:-----|:-------:|:--------|
| `rng_seed` | uint64 | `12345` | Master seed for the outer (Phase 2) Poisson sampler. Combined with the toy index via `std::seed_seq` to give every toy an independent `mt19937_64` stream — re-running with the same seed and same `n_toys` reproduces every toy's input data byte-for-byte. |
| `n_workers` | int ≥ 0 | `0` (= `std::thread::hardware_concurrency()`) | Number of outer toys executed in parallel. Each worker spawns one scan-binary subprocess at a time, which is itself multi-threaded for Minuit/Minimization, so the practical optimum is often `n_cores / scan_internal_threads`. Set to `1` to debug a single-toy reproduction. Capped internally at `n_toys`. |

#### Threshold mode (rejection threshold for q_μ)

| Key | Type | Default | Meaning |
|:----|:-----|:-------:|:--------|
| `threshold` | `"asymptotic"` \| `"toy_mc"` \| `"both"` | `"asymptotic"` | How `q_target` (the rejection threshold for the test statistic q_μ that the scan binary uses to find σ_UL at each mass) is determined. `asymptotic` uses the χ² shortcut (≈ 1.642 for one-sided 90% CL); `toy_mc` runs Phase 0 + Phase 1 to derive a per-mass `q_target` from the empirical q_μ distribution under H_μ; `both` runs both and writes two complete band sets to the same `band.root` for direct comparison. |
| `n_threshold_toys` | int ≥ 1 | `10000` | Phase 1 only. Number of sub-toys per mass drawn from `Poisson(μ·s + b)` at `(m_χ, σ_threshold(m_χ))` to build the q_μ histogram from which the percentile is taken. 5 000 is enough for a stable 90 % quantile; tail percentiles (e.g. 95 %) need more. Cost: `n_threshold_toys × N_mass` full minimizations per band run. |
| `threshold_percentile` | float ∈ (0, 1) | `0.90` | Percentile of the q_μ distribution to use as `q_target`. Must match the CL the scan binary expects (90 % CL → 0.90). Changing this changes only Phase 1; the asymptotic threshold is hard-wired to `q_target = (Φ⁻¹(CL))²` in the scan binary. |
| `threshold_seed` | uint64 | `23456` | Independent seed used by the scan binary inside Phase-1 `threshold_toys` mode. Disjoint from `rng_seed` so changing the outer-toy seed does not invalidate the Phase-1 `q_target_per_mass` lookup (and vice versa). |

#### Where the orchestrator finds the background rates

| Key | Type | Default | Meaning |
|:----|:-----|:-------:|:--------|
| `background_bp_key` | dotted string | `"run.background_Bp"` | JSON dotted path (e.g. `"run.background_Bp"` resolves to `cfg["run"]["background_Bp"]`) to the per-pattern primary background rate vector `Bp_i`. Outer-toy means are `λ_i = Bp_i + Br_i`. Change this only if your config stores the rates under a different key. |
| `background_br_key` | dotted string | `"run.background_Br"` | Same as `background_bp_key` but for the residual / spurious-charge component `Br_i`. Must have the same length as the resolved `Bp` vector. |

#### Output / scratch / I/O contract

| Key | Type | Default | Meaning |
|:----|:-----|:-------:|:--------|
| `ul_graph_name` | string | `"upper_limit_sigma_e_mchi_graph"` | Name of the TGraph the scan binary writes to its output ROOT file containing `(m_χ, σ_UL)` points. The orchestrator opens every `*.root` in the toy's `outdir` and uses the first one that contains a TGraph by this name, so changing this is only needed if you swap to a non-standard scan binary. |
| `tmp_dir` | string | `"<outdir>/_toys"` | Scratch root for Phase 0/1 (`tmp_dir/phase0`, `tmp_dir/phase1`) and per-toy dirs (`tmp_dir/toy_<t>`). Phase 0/1 dirs are always kept (for debugging or re-use); per-toy dirs are removed unless `keep_per_toy_outputs=true`. Put this on a fast local disk (not network home) for big runs. |
| `save_per_toy_curves` | bool | `false` | If true, write a `band_per_toy` TTree to `band.root` with one row per (toy_idx, mchi) carrying both `sigma_UL_asy` and `sigma_UL_toy`. Useful for debugging individual outliers; disabled by default to keep `band.root` small. |
| `keep_per_toy_outputs` | bool | `false` | If true, the per-toy directory `tmp_dir/toy_<t>/{asymptotic,toy_mc}/` (containing the per-toy config, scan log, and ROOT output) is kept after the orchestrator harvests its σ_UL curve. Off by default — N can easily exceed 1000 toys × 2 passes × ~MB of artefacts each. Phase 0/1 outputs are *always* kept regardless of this flag. |

#### Failure-handling

| Key | Type | Default | Meaning |
|:----|:-----|:-------:|:--------|
| `abort_on_failure` | bool | `false` | If true, the orchestrator exits non-zero as soon as the surviving-toy fraction (per pass) drops below `min_success_fraction`. If false (default), it prints a `WARNING` and proceeds to write the band from whatever toys did survive. |
| `min_success_fraction` | float ∈ [0, 1] | `0.95` | Floor on the per-pass survival fraction (separately enforced for `_asymptotic` and `_toy_mc`). Failed toys = scan binary returns non-zero, doesn't write a ROOT, or doesn't write the expected TGraph. Below this fraction either a warning or an error fires (see `abort_on_failure`). Drop to `0.5` for noisy debug runs; raise to `0.99` for production. |

#### Quantiles (the band itself)

| Key | Type | Default | Meaning |
|:----|:-----|:-------:|:--------|
| `quantiles.median` | float ∈ [0, 1] | `0.50` | Quantile drawn as the central "median expected" curve (`median_expected_sigma_e_mchi_*`). |
| `quantiles.low_1sigma` | float ∈ [0, 1] | `0.16` | Lower edge of the 1σ Brazilian band (`band_1sigma_low_sigma_e_mchi_*`). 0.16 ≈ 50 % − 34 %. |
| `quantiles.high_1sigma` | float ∈ [0, 1] | `0.84` | Upper edge of the 1σ band. 0.84 ≈ 50 % + 34 %. |
| `quantiles.low_2sigma` | float ∈ [0, 1] | `0.025` | Lower edge of the 2σ band. 0.025 ≈ 50 % − 47.5 %. |
| `quantiles.high_2sigma` | float ∈ [0, 1] | `0.975` | Upper edge of the 2σ band. 0.975 ≈ 50 % + 47.5 %. |

The defaults reproduce the standard HEP "Brazilian band" convention. Override them only if you want a different summary (e.g. set `low_1sigma=0.05`, `high_1sigma=0.95` for a 90 % central interval) — the scheme keeps producing 5 graphs per threshold variant regardless.

### Threshold modes

| Mode | Phase 0 | Phase 1 | Phase 2 passes per toy | Use when |
|:-----|:-------:|:-------:|:----------------------:|:---------|
| `asymptotic` | skip | skip | 1 (no `q_target_lookup_path`) | Quick comparison; relies on χ² coverage of `q_μ`. |
| `toy_mc` | run | run | 1 (with `q_target_lookup_path`) | Robust band, uses toy-calibrated `q_target` per mass. |
| `both` | run | run | 2 (asymptotic *and* toy_mc) | Validate the asymptotic shortcut against the toy-MC truth. |

### Background sampling source

Outer toys for Phase 2 are sampled from
`λ_i = Bp_i + Br_i` (default, B-only). The keys
`background_bp_key` / `background_br_key` are dotted paths into the JSON
(default `run.background_Bp` / `run.background_Br`); change them if your
config stores the background rates elsewhere (e.g. a different analysis).

---

## Output: `band.root`

| Object | Type | Meaning |
|:-------|:-----|:--------|
| `median_expected_sigma_e_mchi_<asy|toy>` | TGraph | 50 % quantile of σ_UL across surviving toys. |
| `band_1sigma_low_sigma_e_mchi_<asy|toy>` | TGraph | 16 % quantile. |
| `band_1sigma_high_sigma_e_mchi_<asy|toy>` | TGraph | 84 % quantile. |
| `band_2sigma_low_sigma_e_mchi_<asy|toy>` | TGraph | 2.5 % quantile. |
| `band_2sigma_high_sigma_e_mchi_<asy|toy>` | TGraph | 97.5 % quantile. |
| `sigma_threshold_per_mass` | TGraph | (Asimov asymptotic UL used as Phase-1 test point.) |
| `q_target_per_mass` | TGraph | (Phase-1 toy-MC `q_μ` percentile used as rejection threshold.) |
| `band_per_toy` | TTree | Per-toy σ_UL points, only if `save_per_toy_curves=true`. |
| `n_toys`, `rng_seed`, `n_workers`, `n_ok_asy`, `n_ok_toy`, `threshold_*` | TParameter\<double\> | Provenance. |

---

## Plot it

The plotter has two new flags:

```bash
build/ccdarksens_plot_dmelectron_limit \
  outputs/scan_srdm_pattern_csv_per_mass_pydme_match/scan_srdm_pattern_csv_per_mass_pydme_match.root \
  "Observed (per-mass, pydme match)" \
  2.71 heavy \
  --band outputs/band_srdm_pattern_csv_per_mass_pydme_match/band.root \
  --band-mode both \
  --batch
```

`--band-mode` accepts `asymptotic`, `toy_mc`, or `both`. Without it the
plotter auto-detects which variants are present in `band.root` and draws
all available ones. The legend changes labels accordingly.

`--batch` (optional) puts ROOT in non-interactive mode: no GUI window is
opened and the binary returns as soon as the PDF/ROOT outputs are written.
Use it when running headless (CI, remote shells, scripted plotting). Omit
it for interactive use where you want to inspect the canvas live.

Drawing order (back to front): 2σ band → 1σ band → external fills →
external contours → median expected (dashed black for toy_mc, dotted gray
for asymptotic when both are shown) → observed limit (your scan).

---

## Performance & failure modes

* The outer toy loop is parallel: K ≤ `n_workers` scan-binary subprocesses
  run simultaneously. Each subprocess runs the same multi-threaded
  Minuit/Migrad path the nominal scan uses, so total CPU usage easily
  saturates a workstation. Tune `n_workers` if needed.
* Per-toy ROOT files are deleted as soon as the orchestrator harvests
  σ_UL — set `keep_per_toy_outputs=true` to inspect them.
* Toys can fail (Minuit non-convergence, NaNs, etc.). Failed toys are
  excluded from the per-mass quantile. By default a survival fraction
  below `min_success_fraction=0.95` triggers a warning; set
  `abort_on_failure=true` to make it an error.
* Phase-1 toy throughput is bottlenecked by `n_threshold_toys × N_mass`
  full likelihood minimisations. 5 000 sub-toys × 6 masses ≈ a few minutes
  on the SRDM example. Reduce for quick tests.

---

## Reproducibility

* Outer-toy RNG: `mt19937_64` seeded by `seed_seq{rng_seed, t_idx, …}`,
  giving an independent stream per toy index. Re-running with the same
  `rng_seed` and `n_toys` reproduces the band exactly **modulo
  scan-binary nondeterminism**.
* Phase-1 RNG: same scheme, seeded by `threshold_seed` inside the scan
  binary's `threshold_toys` mode.
* Provenance: all relevant seeds, toy counts, and survival counts are
  written to `band.root` as TParameter objects.
