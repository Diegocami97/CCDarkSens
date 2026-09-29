<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_Migdal_Projection.md -- Migdal effect, heavy + light mediator, hypothetical 1 kg-year sensitivity projection.
-->

# Student Example — Migdal Effect, Projection (Heavy + Light Mediator)

Hypothetical 1 kg-year sensitivity projection for the Migdal-effect channel, silicon target, both mediator cases. Companion: [`Student_Examples_Migdal_LBC.md`](Student_Examples_Migdal_LBC.md).

> **KNOWN LIMITATION — read this before trusting any number from this example.** DarkELF (the Migdal rate calculator in `python/ccdarkphys/migdal/entry.py`) reconstructs its internal dielectric/form-factor objects from scratch on every single `compute_dRdE` call — there is no caching across grid points. At ~1-2 seconds of overhead per point regardless of the physics itself, the declared 200×200-*cross-section*-point dense grid would take a long time at the mass counts used elsewhere in this example set. For this example I generated a reduced **150-mass, 200-cross-section "demo" grid** instead (same mass *range*, somewhat coarser mass spacing than the full 200-point declared grid). Confirmed complete via filename parsing. The resulting curve is coarser in mass than what the full grid would show, but not wrong — treat it as good enough to learn the workflow and see the right qualitative shape, not as the production-resolution DAMIC-M projection.

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | Migdal effect (DM-nucleus scattering with electron ionization), Si, both **heavy** and **light** mediator |
| Observable | n_e bins |
| Exposure | 1 kg · 1 year (hypothetical) |
| Background | Flat dark current + flat d.r.u., θ-scaled |
| Data | Asimov (S = B) |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σₙ(mχ) |

**Configs:**
- [`configs/examples/migdal_heavy_projection_1kgyear.json`](../configs/examples/migdal_heavy_projection_1kgyear.json) — rates in `data/migdal_rates/Si/heavy_demo150x200/` (30,000 files, complete)
- [`configs/examples/migdal_light_projection_1kgyear.json`](../configs/examples/migdal_light_projection_1kgyear.json) — rates in `data/migdal_rates/Si/light_demo150x200/` (30,000 files, complete)

**Mass range: 0.1–1000 MeV** (150 log-spaced points), matching the LBC companion rather than an earlier, narrower 10–1000 MeV range this example originally inherited from `configs/migdal_scan_si_heavy.json`'s declared grid. That 10 MeV floor was never a physics or DarkELF limitation — a direct check (`compute_dRdE` down to 0.5 MeV) confirms DarkELF's Lindhard-method Migdal calculation works fine at low mass — it was just copied from a config that happened to start there. It mattered because DAMIC-M's own PRL reports the Migdal channel as **its most stringent limit for m_χ = 1–35 MeV** (`docs/DM_Signal_Models_Physics_Reference.md`), a range this example was previously cutting off entirely.

---

## 2. Prerequisites

Build `ccdarksens_scan_generic` and `ccdarksens_plot_limit` as usual. Regenerating rates (only needed if you want the full-resolution grid) additionally needs DarkELF installed and importable from `python/ccdarkphys/migdal/entry.py`.

---

## 3. Pipeline overview

```
Step 1   Generate Migdal dR/dE(E) rate CSVs via DarkELF (Python)   -- already done (demo grid)
           |
Step 2   Run the n_e-space Asimov scan (C++)
           |
Step 3   Plot the limit curve (C++)
```

**What the pipeline is physically doing:** DarkELF computes the Migdal-effect differential rate dR/dE for a DM-nucleus collision that ionizes a bound electron, for silicon, at each (mχ, σₙ) grid point. This spectrum is folded through the same ionization table used elsewhere in the repo to get S(n_e), then combined with a flat-DC/flat-d.r.u. Asimov background exactly as in the DM-electron projection examples — the statistical machinery downstream of the rate table is identical across all DM-electron and Migdal projections.

---

## 4. Step 1 — Rate generation (informational — the demo grid is already built)

```bash
python3 utils/migdal_generate_grid.py configs/examples/migdal_generate_si_heavy_demo150x200.json
python3 utils/migdal_generate_grid.py configs/examples/migdal_generate_si_light_demo150x200.json
```

150×200 reduced-grid configs (mass range extended to 0.1–1000 MeV at 150 points instead of the originally-inherited 10–1000 MeV/40 points — see the box above; cross section bumped from an earlier 40 to 200 — see below), `parallel: 8`. They only reduce the grid density relative to a full-scale config — no physics settings differ.

To regenerate the full 200×200 grid, use the full-resolution generate config and budget several hours, or first add caching to `ccdarkphys.migdal.entry` (a real performance bug, flagged but not fixed as part of this task since it requires a source-code change).

**Why 200 cross-section points, not 40 — see the companion LBC doc §4 for the full story.** The short version: `--from-qhist` linearly interpolates the upper limit between the two σ grid points bracketing the q = q_threshold crossing, so a coarse σ grid produces a visibly faceted *and numerically biased* curve, not just a "less smooth-looking" one — this was originally caught on the LBC pair (which has a real published curve to compare against) and then applied here too for consistency, even though there's no external reference for the projection case to validate against directly. σ resolution is cheap to increase (the physics is computed once per mass, all σ values come from an exact linear rescale), so there's no reason to stay at the old, coarser grid.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/migdal_heavy_projection_1kgyear.json
build/ccdarksens_scan_generic configs/examples/migdal_light_projection_1kgyear.json
```

Outputs: `outputs/migdal_heavy_projection_1kgyear/scan_generic.root`, `outputs/migdal_light_projection_1kgyear/scan_generic.root`.

**Runtime:** a few seconds each (8,000 points, measured directly).

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this projection on the same canvas as the reproduced LBC result and the published literature curves.** Run the LBC companion's scans first ([`Student_Examples_Migdal_LBC.md`](Student_Examples_Migdal_LBC.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --migdal --heavy \
  --out-pdf outplots/migdal_heavy_combined.pdf \
  outputs/migdal_heavy_lbc_1p3kgday/scan_generic.root "Migdal heavy LBC (reproduced)" \
  outputs/migdal_heavy_projection_1kgyear/scan_generic.root "Migdal heavy, 1 kg-yr projection" \
  1.642374415149816 heavy

build/ccdarksens_plot_limit --batch --from-qhist --migdal --light \
  --out-pdf outplots/migdal_light_combined.pdf \
  outputs/migdal_light_lbc_1p3kgday/scan_generic.root "Migdal light LBC (reproduced)" \
  outputs/migdal_light_projection_1kgyear/scan_generic.root "Migdal light, 1 kg-yr projection" \
  1.642374415149816 light
```

`--migdal` switches the plot to σₙ (DM-nucleon) axes and automatically loads the published Migdal literature curves from `data/previous_limits/Migdal_heavy/` (XENON1T, PandaX-4T, DarkSide-50, SENSEI — plus DAMIC-M 2025 itself in `--heavy` mode only, since that digitization is heavy-mediator-specific and `data/previous_limits/Migdal_light/` has no light-mediator equivalent yet), building a shaded lower-envelope exclusion band from all of them — no `--show-damic` needed for this mode. `--from-qhist` forces both curves to be rebuilt directly from each file's own q(mχ,σ) histogram, so they're computed identically and directly comparable. This is the primary, recommended way to look at this example's output. Rendered copies are checked in at [`migdal_heavy_combined.pdf`](migdal_heavy_combined.pdf) and [`migdal_light_combined.pdf`](migdal_light_combined.pdf).

If you only want one mediator's projection curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --migdal --heavy \
  --out-pdf outplots/migdal_heavy_projection_1kgyear.pdf \
  outputs/migdal_heavy_projection_1kgyear/scan_generic.root \
  "Migdal heavy, 1 kg-yr projection" 1.642374415149816 heavy
```

(swap `--light`/light paths for that mediator.)
