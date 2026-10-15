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
| Signal model | DM-electron scattering, Si, **heavy** mediator (QEDark), standard halo of Baxter et al.: v₀ = 238, **v_E = 253.7**, v_esc = 544 km/s |
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Exposure | **1.3 kg-day** (verified from the config: `mass_kg=0.01523 × livetime_days=85.356 = 1.29997 kg-day`) |
| Background | `B = Bp + θ·Br`, the **real** LBC values `Bp=[141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]`, `Br=[0.039, 0.039, 0.016, 0.052, 0.011, 0.035]` |
| Data | **Real observed pattern counts**: `[144, 0, 0, 1, 0, 0]` |
| Statistic | Profile likelihood ratio (pydme-style), one-sided q with a robust reference minimum, 90% CL upper limit on σₑ(mχ); asymptotic threshold 1.642, checked with toy Monte Carlo (§5) |

**Config:** [`configs/examples/dm_electron_heavy_lbc_1p3kgday.json`](../configs/examples/dm_electron_heavy_lbc_1p3kgday.json)

**Grid:** 800 masses (0.2–1000 MeV) × 300 cross sections (10⁻⁴⁶–10⁻²⁶ cm²), 240,000 rate files in `data/qedark_rates/Si/heavy/long_scan_vE253p7` (a config-level bug fix was applied here, see §7).

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
Step 1   Generate QEDark rate CSVs (Python, v_E = 253.7 km/s)
           |
Step 2   Run the pattern-space scan (C++): the result
           |         optional: check the threshold with toys and scan again
Step 3   Plot the limit curve with the published DAMIC-M overlay (C++)
```

**What the pipeline is physically doing:** at each grid point, a QEDark dR/dE(E) spectrum is folded to S_true(n_e) via the ionization table, then to S(pattern) via the pattern efficiency table (`data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`). The background per pattern bin is `Bp + θ·Br` with the real digitized values above. θ is profiled with a tau-weighted Gamma prior (strength 98) at each trial σₑ, giving q(σ) = 2ΔlnL against the real observed counts, and the 90% CL upper limit is where q crosses the threshold: 1.642 for the asymptotic one-sided test (step 1 of §5, the result). A per-mass threshold from toy Monte Carlo (step 3 of §5) is a check; it is higher where the limit is set by bins that expect about 0.1 events and saw none.

---

## 4. Step 1 — Rate generation

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_heavy_dense.json
```

Output: `data/qedark_rates/Si/heavy/long_scan_vE253p7/`, 240,000 files, about 8 GB, a few minutes on 8 processes. The halo is the one of Baxter et al. (arXiv:2105.00599, Table 1): v_E = |v₀ + v_⊙ + v_⊕(March 9)| = 253.7 km/s. The older `long_scan` folder was made with v_E = 263 km/s, which makes the limit 30–40% too strong below 10 MeV, and is no longer read by this example.

---

## 5. Step 2 — Run the scan

Step 1 is the result. Steps 2 and 3 are an optional check of its threshold; the config of each step is written for the next:

```bash
# 1. the result: asymptotic scan (q threshold 1.642), one-sided q, robust minimum
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday.json
# 2. toy Monte Carlo: the 90th percentile of q at the cross section of the step-1 limit, at every mass
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday_toys.json
# 3. the same scan with that per-mass threshold (the check)
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday_mc.json
```

Outputs: `outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root` (step 1), `outputs/dm_electron_heavy_lbc_1p3kgday_toys/qtarget_threshold.root` (step 2) and `outputs/dm_electron_heavy_lbc_1p3kgday_mc/scan_generic.root` (step 3, the check of step 1).

**Why the check.** The asymptotic threshold is accurate where one bin with many counts decides the limit (`{11}`, above 10 MeV and below 3 MeV) and too strong where `{21}` and `{111}` decide it, between about 3 and 7 MeV: those bins expect about 0.1 events, saw none, and the true 90% threshold needs about 2.3 signal events rather than the 0.8 that q = 1.642 corresponds to. The toy step measures it, and step 3 shows how much the limit moves (a median 1.44 between 3 and 10 MeV, nothing elsewhere). It uses the one-sided statistic (`run.threshold_toys.one_sided: true`): a toy whose best fit lies above the tested cross section contributes q = 0. The two-sided statistic recommended by Baxter et al. gave a limit 20–35% weaker than the published curve above 1 MeV in a test. Against the published DAMIC-M curve the step-1 limit is within 4% at 1–3 MeV and above 12 MeV, 27% below it at 3–6 MeV and within 8% at 6–12 MeV; the Monte Carlo check is within 7% at 3–6 MeV and 17–34% above it between 6 and 12 MeV. Neither reproduces the published curve at every mass; the walkthrough (§9.1) gives the table and the two readings.

**Runtime:** about 2 minutes for step 1 and for step 3, about 16 minutes for step 2 (800 masses × 10,000 toys), on 8 cores. For a quick check, copy the three configs and reduce `model.grid.mchi_MeV.logspace.num` to e.g. 5.

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this reproduction on the same canvas as the 1 kg-year projection and the real published DAMIC-M curve**, so you can see at a glance both how well this reproduces the paper and how much a 1 kg-year exposure would gain over it. Run the projection companion's scan first ([`Student_Examples_DM_Electron_Heavy_Projection.md`](Student_Examples_DM_Electron_Heavy_Projection.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --heavy --from-qhist --show-damic \
  --out-pdf outplots/dm_electron_heavy_combined.pdf \
  outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root "DM-e heavy LBC (reproduced)" \
  outputs/dm_electron_heavy_projection_1kgyear/scan_generic.root "DM-e heavy, 1 kg-yr projection"
```

`--from-qhist` rebuilds both curves from each file's own q(mχ,σ) at 90% C.L., so they are computed the same way and are comparable with each other and with the published curve. The projection scan is an Asimov scan with its own QCDark2 rates and halo; it is not part of this threshold check.

To see what the Monte Carlo check does to this reproduction, draw the stored limits of step 1 (asymptotic 1.642, the result) and step 3 (Monte Carlo) together and print the ratio (`--use-stored-limit`, because `--from-qhist` would rebuild step 3 at 1.642 as well):

```bash
build/ccdarksens_plot_limit --batch --heavy --use-stored-limit --show-damic --ratio-to 1 \
  --out-pdf outplots/dm_electron_heavy_lbc_threshold.pdf \
  outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root "asymptotic 1.642" \
  outputs/dm_electron_heavy_lbc_1p3kgday_mc/scan_generic.root "toy Monte Carlo"
```

---

## 7. A bug this example surfaced (fixed here, not yet fixed upstream)

While building this example, running the *original*, already-committed `configs/examples/lbc_qedark_heavy_mediator.json` (this config's ancestor) failed with:

```
ERROR: ResponseFactory: failed to open pattern efficiency CSV: .../configs/data/Efficiencies_patterns_...csv
```

This is a real bug in `ResolveConfigRelativePath` (`src/response/ResponseFactory.cc`): it resolves a config's relative paths against `config_dir/..`, which is correct for a config living directly under `configs/`, but is one directory too shallow for anything under `configs/examples/` — it looks in `configs/data/...` instead of the repo-root `data/...`. It reproduces identically with the untouched reference app `ccdarksens_scan_dmelectron_pattern`, so it is not new or scan_generic-specific.

**No C++ was changed.** The fix applied in this example's config (and every other new example config that uses `efficiency_mc.efficiency_csv`) is a one-line value change: `"data/Efficiencies_..."` → `"../data/Efficiencies_..."`. The pre-existing `lbc_qedark_heavy_mediator.json` / `lbc_qedark_light_mediator.json` / `migdal_scan_si_heavy.json` / `migdal_scan_si_light.json` still have the old, broken value and will fail with the same error until someone applies the same one-line fix (or the C++ resolver is corrected).
