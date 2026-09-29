<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_DM_Electron_Light_Projection.md -- DM-electron, ultralight mediator, hypothetical 1 kg-year sensitivity projection.
-->

# Student Example — DM-Electron, Light (Ultralight) Mediator, Projection

Hypothetical 1 kg-year sensitivity projection for the ultralight (massless) mediator case. Companion: [`Student_Examples_DM_Electron_Light_LBC.md`](Student_Examples_DM_Electron_Light_LBC.md).

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | DM-electron scattering, Si, **ultralight/massless** mediator (QEDark) |
| Observable | n_e bins 1–5 |
| Exposure | 1 kg · 1 year (hypothetical) |
| Background | Flat dark current (1×10⁻⁵ e⁻/pixel/day) + flat d.r.u., θ-scaled |
| Data | Asimov (S = B) |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σₑ(mχ) |

**Config:** [`configs/examples/dm_electron_light_projection_1kgyear.json`](../configs/examples/dm_electron_light_projection_1kgyear.json)

**Note on this config's construction:** unlike the heavy-mediator projection (which uses a QCDark2 rate table), no QCDark2 light-mediator dense grid exists yet in this repo. This example reuses the already-complete **QEDark ultralight** rate table (`data/qedark_rates/Si/ultralight/long_scan/`, 240,000 files, same table the LBC case uses) instead, just analyzed in n_e space with a flat-DC background rather than pattern space with the real LBC background. If a QCDark2 light-mediator grid is generated later, swap `model.rates_dir`.

**Grid:** 800 masses × 300 cross sections, same range as the heavy case.

---

## 2. Prerequisites

Same as the heavy-mediator projection: build `ccdarksens_scan_generic` and `ccdarksens_plot_limit`.

---

## 3. Pipeline overview

Identical to the heavy-mediator projection (§3 of that doc) — same ionization table, same flat-DC/flat-d.r.u. background construction, same profile likelihood — only the rate table and the resulting σₑ vs. mχ curve's shape differ (ultralight mediator has an extra 1/q⁴ enhancement at low momentum transfer).

---

## 4. Step 1 — Rate generation (informational — already done)

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_light_dense.json
```

Output: `data/qedark_rates/Si/ultralight/long_scan/` — already populated.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/dm_electron_light_projection_1kgyear.json
```

Output: `outputs/dm_electron_light_projection_1kgyear/scan_generic.root`.

**Runtime:** ~90 seconds for the full 800×300 grid on a normal laptop (measured directly).

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this projection on the same canvas as the reproduced LBC result and the real published DAMIC-M curve**, so you can see at a glance both how well the LBC reproduction tracks the paper and how much a 1 kg-year exposure would gain over it. Run the LBC companion's scan first ([`Student_Examples_DM_Electron_Light_LBC.md`](Student_Examples_DM_Electron_Light_LBC.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --show-damic \
  --out-pdf outplots/dm_electron_light_combined.pdf \
  outputs/dm_electron_light_lbc_1p3kgday/scan_generic.root "DM-e light LBC (reproduced)" \
  outputs/dm_electron_light_projection_1kgyear/scan_generic.root "DM-e light, 1 kg-yr projection" \
  1.642374415149816 light
```

`--from-qhist` forces both curves to be rebuilt directly from each file's own q(mχ,σ) histogram rather than any pre-stored UL graph, so the two are computed identically and are directly comparable to each other and to `--show-damic`'s canonical paper-export DAMIC-M (2025) curve. This is the primary, recommended way to look at this example's output.

If you only want this projection's own curve in isolation:

```bash
build/ccdarksens_plot_limit --batch \
  --out-pdf outplots/dm_electron_light_projection_1kgyear.pdf \
  outputs/dm_electron_light_projection_1kgyear/scan_generic.root \
  "DM-e light, 1 kg-yr projection" 1.642374415149816 light
```

The trailing `light` selects the light-mediator literature overlay set.

---

## 7. Note

The same `efficiency_csv` path bug described in the heavy-LBC doc's §7 applies here too; this config already has the one-line fix (`"../data/Efficiencies_..."`) applied.
