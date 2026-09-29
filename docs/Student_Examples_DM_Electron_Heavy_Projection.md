<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_DM_Electron_Heavy_Projection.md -- DM-electron, heavy mediator, hypothetical 1 kg-year sensitivity projection.
-->

# Student Example — DM-Electron, Heavy Mediator, Projection

A worked example for a **hypothetical / future-exposure sensitivity projection**: no real background measurement is used, only a flat dark-current + flat d.r.u. background model, at a 1 kg-year exposure. This is the workflow to copy when you want to ask "what would DAMIC-M see for model X at exposure Y?" for a model that doesn't have real background data yet.

Companion case: [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) — same physics model, but the real LBC background and observed counts instead of this flat placeholder.

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | DM-electron scattering, Si, **heavy** mediator (QCDark2, Si_comp composite dielectric) |
| Observable | n_e bins 1–5 |
| Exposure | 1 kg · 1 year (hypothetical — we do not have 1 kg-year of real DAMIC-M data yet) |
| Background | Flat dark current (1×10⁻⁵ e⁻/pixel/day) + flat d.r.u., folded to B(n_e); nuisance θ scales it |
| Data | **Asimov** (S = B at every grid point — there is no real dataset here) |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σₑ(mχ) |

**Config:** [`configs/examples/dm_electron_heavy_projection_1kgyear.json`](../configs/examples/dm_electron_heavy_projection_1kgyear.json)

**Grid:** 800 masses (0.2–1000 MeV, log) × 300 cross sections (10⁻⁴⁶–10⁻²⁶ cm², log) — 240,000 rate points, already generated in this repo.

This is the same physics case documented in [`QCDark2_NE_Exposure_Projections_Guide.md`](QCDark2_NE_Exposure_Projections_Guide.md) at 1 kg-year; this doc is the same workflow, repackaged and renamed for a clean student-facing set. If you need 0.5 or 2.0 kg-year instead, that guide's sibling configs cover those.

---

## 2. Prerequisites

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_generic ccdarksens_plot_limit
```

Needs ROOT (with Minuit2), nlohmann/json, and for rate generation: Python 3.8+, NumPy, and the `ccdarkphys` package (`PYTHONPATH=python`).

---

## 3. Pipeline overview

```
Step 1   Generate QCDark2 rate CSVs (Python)   -- already done in this repo checkout
           |
Step 2   Run the n_e-space Asimov scan (C++)
           |
Step 3   Plot the limit curve (C++)
```

**What the pipeline is physically doing:** a QCDark2 dR/dE(E) spectrum for silicon (using the composite Si_comp dielectric function) is folded through the ionization table P(n_e|E) to get S(n_e), a flat-DC + flat-d.r.u. Asimov background B(n_e) is built from the detector geometry and exposure, and a profile-likelihood ratio scans a nuisance parameter θ that rescales B to get the 90% CL upper limit on σₑ at each mass point. No real detector data enters anywhere in this path — everything is a projection.

---

## 4. Step 1 — Rate generation (informational — already done)

```bash
python3 utils/qcdark2_generate_grid.py configs/examples/qcdark2_generate_si_comp_dense.json
```

Output: `data/qcdark2_rates/Si/heavy/Si_comp_long_scan/` (240,000 files, verified complete: 800 masses × 300 cross sections, exact match). This directory already exists in the repo — you do not need to regenerate it to run this example. Regenerate only if you change the halo, gap, or grid.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_projection_1kgyear.json
```

Output: `outputs/dm_electron_heavy_projection_1kgyear/scan_generic.root`, containing `q_mchi_sigma` (TH2D), `upper_limit_sigma_e_mchi_graph` (TGraph, the primary result), and `exposure_kg_year`.

To change the exposure, copy this config and edit `detector.mass_kg` (exposure_kg_year = mass_kg when `livetime_days=365.25`, `duty_cycle=1.0`) and `run.outdir`/`run.label`. To change the dark-current assumption, edit `backgrounds.dark_current.lambda_e_per_pix_per_year` (1×10⁻⁵ e⁻/pixel/day = 0.00365 e⁻/pixel/year).

**Runtime:** ~90 seconds for the full 800×300 grid on a normal laptop (measured directly; earlier guidance in this doc overstated this as "tens of minutes to a few hours" — it does not take that long). For a quick check anyway, copy the config and reduce `model.grid.mchi_MeV.logspace.num` to e.g. 10.

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this projection on the same canvas as the reproduced LBC result and the real published DAMIC-M curve**, so you can see at a glance both how well the LBC reproduction tracks the paper and how much a 1 kg-year exposure would gain over it. Run the LBC companion's scan first ([`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --show-damic \
  --out-pdf outplots/dm_electron_heavy_combined.pdf \
  outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root "DM-e heavy LBC (reproduced)" \
  outputs/dm_electron_heavy_projection_1kgyear/scan_generic.root "DM-e heavy, 1 kg-yr projection" \
  1.642374415149816 heavy
```

`--from-qhist` forces both curves to be rebuilt directly from each file's own q(mχ,σ) histogram rather than any pre-stored UL graph, so the two are computed identically and are directly comparable to each other and to `--show-damic`'s canonical paper-export DAMIC-M (2025) curve. This is the primary, recommended way to look at this example's output.

If you only want this projection's own curve in isolation:

```bash
build/ccdarksens_plot_limit --batch \
  --out-pdf outplots/dm_electron_heavy_projection_1kgyear.pdf \
  outputs/dm_electron_heavy_projection_1kgyear/scan_generic.root \
  "DM-e heavy, 1 kg-yr projection" 1.642374415149816 heavy
```

`1.642374415149816` is the 90% CL, 1-dof threshold used for pydme-style parity elsewhere in this repo — keep it for consistency, or use ROOT's own `TMath::NormQuantile(0.90)^2` if you prefer to derive it yourself. The trailing `heavy` selects the heavy-mediator literature overlay set.
