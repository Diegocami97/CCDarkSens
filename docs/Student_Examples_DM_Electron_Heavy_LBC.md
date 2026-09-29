<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_DM_Electron_Heavy_LBC.md -- DM-electron, heavy mediator, real LBC pattern-space background/counts reproduction.
-->

# Student Example — DM-Electron, Heavy Mediator, LBC Reproduction

A worked example using the **real LBC pattern-space background template and real observed pattern counts** — this reproduces an actual published-style analysis, not a projection. Companion case: [`Student_Examples_DM_Electron_Heavy_Projection.md`](Student_Examples_DM_Electron_Heavy_Projection.md).

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | DM-electron scattering, Si, **heavy** mediator (QEDark) |
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Exposure | **1.3 kg-day** (verified from the config: `mass_kg=0.01523 × livetime_days=85.356 = 1.29997 kg-day`) |
| Background | `B = Bp + θ·Br`, the **real** LBC values `Bp=[141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]`, `Br=[0.039, 0.039, 0.016, 0.052, 0.011, 0.035]` |
| Data | **Real observed pattern counts**: `[144, 0, 0, 1, 0, 0]` |
| Statistic | Profile likelihood ratio (pydme-style), 90% CL upper limit on σₑ(mχ) |

**Config:** [`configs/examples/dm_electron_heavy_lbc_1p3kgday.json`](../configs/examples/dm_electron_heavy_lbc_1p3kgday.json)

**Grid:** 800 masses (0.2–1000 MeV) × 300 cross sections (10⁻⁴⁶–10⁻²⁶ cm²) — verified complete (895,363 files on disk; a config-level bug fix was applied here, see §7).

---

## 2. Prerequisites

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_generic ccdarksens_plot_limit
```

**Required input:** `build/data_pattern.root` (observed pattern counts + exposure metadata from the pydme LBC export) must exist. It is present in this checkout.

---

## 3. Pipeline overview

```
Step 1   Generate QEDark rate CSVs (Python)   -- already done in this repo checkout
           |
Step 2   Run the pattern-space scan (C++)
           |
Step 3   Plot the limit curve with the pydme reference overlay (C++)
```

**What the pipeline is physically doing:** at each grid point, a QEDark dR/dE(E) spectrum is folded to S_true(n_e) via the ionization table, then to S(pattern) via the pattern efficiency table (`data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`). The background per pattern bin is `Bp + θ·Br` with the real digitized values above. θ is profiled with a tau-weighted Gamma prior (strength 98) at each trial σₑ, giving q(σ) = 2ΔlnL against the real observed counts, and the 90% CL upper limit is where q crosses 1.642374415149816.

---

## 4. Step 1 — Rate generation (informational — already done)

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_heavy_dense.json
```

Output: `data/qedark_rates/Si/heavy/long_scan/` — already populated, 895,363 files (more than the declared 240,000-point grid; harmless extras from earlier iterations).

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday.json
```

Output: `outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root`.

**Runtime:** ~90 seconds for the full 800×300 grid with profiling on a normal laptop (measured directly; earlier guidance in this doc overstated this as "expect many hours" — it does not take that long). For a quick check anyway, copy the config and reduce `model.grid.mchi_MeV.logspace.num` to e.g. 5.

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this reproduction on the same canvas as the 1 kg-year projection and the real published DAMIC-M curve**, so you can see at a glance both how well this reproduces the paper and how much a 1 kg-year exposure would gain over it. Run the projection companion's scan first ([`Student_Examples_DM_Electron_Heavy_Projection.md`](Student_Examples_DM_Electron_Heavy_Projection.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --show-damic \
  --out-pdf outplots/dm_electron_heavy_combined.pdf \
  outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root "DM-e heavy LBC (reproduced)" \
  outputs/dm_electron_heavy_projection_1kgyear/scan_generic.root "DM-e heavy, 1 kg-yr projection" \
  1.642374415149816 heavy
```

`--from-qhist` forces both curves to be rebuilt directly from each file's own q(mχ,σ) histogram rather than any pre-stored UL graph, so the two are computed identically and are directly comparable to each other and to `--show-damic`'s canonical paper-export DAMIC-M (2025) curve. This is the primary, recommended way to look at this example's output.

If you only want this reproduction's own curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --draw-both --show-damic \
  --out-pdf outplots/dm_electron_heavy_lbc_1p3kgday.pdf \
  outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root \
  "DM-e heavy LBC" 1.642374415149816 heavy
```

`--draw-both` draws the pydme-bisection diagnostic curve alongside the stored UL curve, if present.

---

## 7. A bug this example surfaced (fixed here, not yet fixed upstream)

While building this example, running the *original*, already-committed `configs/examples/lbc_qedark_heavy_mediator.json` (this config's ancestor) failed with:

```
ERROR: ResponseFactory: failed to open pattern efficiency CSV: .../configs/data/Efficiencies_patterns_...csv
```

This is a real bug in `ResolveConfigRelativePath` (`src/response/ResponseFactory.cc`): it resolves a config's relative paths against `config_dir/..`, which is correct for a config living directly under `configs/`, but is one directory too shallow for anything under `configs/examples/` — it looks in `configs/data/...` instead of the repo-root `data/...`. It reproduces identically with the untouched reference app `ccdarksens_scan_dmelectron_pattern`, so it is not new or scan_generic-specific.

**No C++ was changed.** The fix applied in this example's config (and every other new example config that uses `efficiency_mc.efficiency_csv`) is a one-line value change: `"data/Efficiencies_..."` → `"../data/Efficiencies_..."`. The pre-existing `lbc_qedark_heavy_mediator.json` / `lbc_qedark_light_mediator.json` / `migdal_scan_si_heavy.json` / `migdal_scan_si_light.json` still have the old, broken value and will fail with the same error until someone applies the same one-line fix (or the C++ resolver is corrected).
