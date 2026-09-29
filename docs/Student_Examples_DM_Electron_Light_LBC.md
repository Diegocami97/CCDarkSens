<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_DM_Electron_Light_LBC.md -- DM-electron, ultralight mediator, real LBC pattern-space background/counts reproduction.
-->

# Student Example — DM-Electron, Light (Ultralight) Mediator, LBC Reproduction

Real LBC background and observed counts, ultralight/massless mediator case. Companion: [`Student_Examples_DM_Electron_Light_Projection.md`](Student_Examples_DM_Electron_Light_Projection.md).

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | DM-electron scattering, Si, **ultralight/massless** mediator (QEDark) |
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Exposure | 1.3 kg-day (same detector/livetime as the heavy-mediator LBC case) |
| Background | Real LBC `Bp`/`Br` (same values as the heavy case — the background template describes the detector/data, not the signal model) |
| Data | Real observed counts `[144, 0, 0, 1, 0, 0]` (same dataset as the heavy case) |
| Statistic | Profile likelihood ratio (pydme-style), 90% CL upper limit on σₑ(mχ) |

**Config:** [`configs/examples/dm_electron_light_lbc_1p3kgday.json`](../configs/examples/dm_electron_light_lbc_1p3kgday.json)

**Grid:** 800 masses × 300 cross sections; rates in `data/qedark_rates/Si/ultralight/long_scan/` (240,000 files, complete).

---

## 2–3. Prerequisites and pipeline

Identical to [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) §2–3 — same `build/data_pattern.root` requirement, same fold → pattern-efficiency → Bp+θBr → profile-likelihood chain. Only the rate table differs.

---

## 4. Step 1 — Rate generation (informational — already done)

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_light_dense.json
```

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/dm_electron_light_lbc_1p3kgday.json
```

Output: `outputs/dm_electron_light_lbc_1p3kgday/scan_generic.root`.

**Runtime:** ~90 seconds for the full 800×300 grid with profiling on a normal laptop (measured directly).

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this reproduction on the same canvas as the 1 kg-year projection and the real published DAMIC-M curve**, so you can see at a glance both how well this reproduces the paper and how much a 1 kg-year exposure would gain over it. Run the projection companion's scan first ([`Student_Examples_DM_Electron_Light_Projection.md`](Student_Examples_DM_Electron_Light_Projection.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --show-damic \
  --out-pdf outplots/dm_electron_light_combined.pdf \
  outputs/dm_electron_light_lbc_1p3kgday/scan_generic.root "DM-e light LBC (reproduced)" \
  outputs/dm_electron_light_projection_1kgyear/scan_generic.root "DM-e light, 1 kg-yr projection" \
  1.642374415149816 light
```

`--from-qhist` forces both curves to be rebuilt directly from each file's own q(mχ,σ) histogram rather than any pre-stored UL graph, so the two are computed identically and are directly comparable to each other and to `--show-damic`'s canonical paper-export DAMIC-M (2025) curve. This is the primary, recommended way to look at this example's output.

If you only want this reproduction's own curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --draw-both --show-damic \
  --out-pdf outplots/dm_electron_light_lbc_1p3kgday.pdf \
  outputs/dm_electron_light_lbc_1p3kgday/scan_generic.root \
  "DM-e light LBC" 1.642374415149816 light
```

`--draw-both` draws the pydme-bisection diagnostic curve alongside the stored UL curve, if present.

To compare heavy and light mediator LBC curves on one plot:

```bash
build/ccdarksens_plot_limit --batch --show-damic \
  --out-pdf outplots/dm_electron_lbc_both_mediators.pdf \
  outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root "heavy" \
  outputs/dm_electron_light_lbc_1p3kgday/scan_generic.root "light" \
  1.642374415149816 heavy
```

---

## 7. Note

Same `efficiency_csv` path fix as the heavy case applies (already in this config) — see the heavy-LBC doc §7 for the underlying bug.
