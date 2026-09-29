<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_DM_Electron_Intermediate_Projection.md -- DM-electron, intermediate-mass mediator (mA'=5,10 keV), hypothetical 1 kg-year sensitivity projection.
-->

# Student Example — DM-Electron, Intermediate-Mass Mediator, Projection

Hypothetical 1 kg-year sensitivity projection for two intermediate-mass-mediator slices (mA' = 5 keV and 10 keV). Companion: [`Student_Examples_DM_Electron_Intermediate_LBC.md`](Student_Examples_DM_Electron_Intermediate_LBC.md).

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | DM-electron scattering, Si, **intermediate-mass mediator**, mA' = **5 keV and 10 keV** (QCDark2, composite Si_comp.h5 dielectric function — same DFT treatment as the heavy-mediator case) |
| Observable | n_e bins 1–5 |
| Exposure | 1 kg · 1 year (hypothetical) |
| Background | Flat dark current (1×10⁻⁵ e⁻/pixel/day) + flat d.r.u., θ-scaled |
| Data | Asimov (S = B) |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σ̄ₑ(mχ) |

**Configs:**
- [`configs/examples/dm_electron_intermediate_mA5keV_projection_1kgyear.json`](../configs/examples/dm_electron_intermediate_mA5keV_projection_1kgyear.json)
- [`configs/examples/dm_electron_intermediate_mA10keV_projection_1kgyear.json`](../configs/examples/dm_electron_intermediate_mA10keV_projection_1kgyear.json)

Both are identical except for `model.rates_dir` (`mA_5keV_comp/` vs `mA_10keV_comp/`) and `run.label`/`run.outdir`.

**Grid:** 800 masses (0.2–1000 MeV, log) × 8 cross sections (10⁻³⁵–10⁻⁴² cm²) = 6,400 rate points per slice — the same dense mass grid used by the heavy/light-mediator examples, generated fresh for this example set (verified complete, zero errors).

This config's `model.type` is `dm_electron` with a numeric `mediator` (mA' in eV, e.g. `5000.0`) — QCDark2 natively supports an arbitrary intermediate mediator mass via `F_DM = (m_A² + (αm_e)²) / (m_A² + (αm_e·q)²)` (`qcdark2/dark_matter_rates.py::get_F_DM`), the standard interpolating form factor between the light (`m_A→0`) and heavy (`m_A→∞`) limits. This is a real, first-class QCDark2 capability, not a special-cased hack.

---

## 2. Prerequisites

Same build as the other DM-electron examples (`ccdarksens_scan_generic`, `ccdarksens_plot_limit`). Rate generation, if you ever need to redo it, additionally needs a working `qcdark2` Python install (see [`Student_Examples_DM_Electron_Heavy_Projection.md`](Student_Examples_DM_Electron_Heavy_Projection.md) §2) and the composite `Si_comp.h5` dielectric-function file from your own QCDark2 checkout.

---

## 3. Pipeline overview

Identical fold → flat-DC-background → profile-likelihood chain as the heavy/light-mediator projections (see [`Student_Examples_DM_Electron_Heavy_Projection.md`](Student_Examples_DM_Electron_Heavy_Projection.md) §3). The only difference from the heavy-mediator case is the `mediator` value passed into QCDark2 — a float (mA' in eV) instead of the string `"heavy"` — which selects the intermediate form factor instead of `F_DM=1`.

---

## 4. Step 1 — Rate generation (informational — already done)

```bash
python3 utils/qcdark2_generate_grid.py configs/qcdark2_generate_Si_intermediate_mA5keV_comp.json
python3 utils/qcdark2_generate_grid.py configs/qcdark2_generate_Si_intermediate_mA10keV_comp.json
```

Each of these configs sets `"mediator": "5000.0"` / `"10000.0"` (mA' in eV) and `"epsilon_h5"` pointing at your local QCDark2 checkout's `dielectric_functions/composite/Si_comp.h5` (the checked-in config has a `path/to/QCDark2/...` placeholder — substitute your own path or set `CCDARK_SENS_QCDARK2_EPSILON`, same as the heavy-mediator generate config).

Output: `data/qcdark2_rates/Si/intermediate/mA_5keV_comp/` and `mA_10keV_comp/`, 6,400 files each (800 masses × 8 cross sections), already generated and complete in this repo checkout. Took ~95–105 seconds per slice on a normal laptop (measured directly).

**A note on this being a first-time-working combination:** getting this to run required fixing a stale `qcdark2` package install (the environment's `site-packages` copy predated the intermediate-mediator support added to the QCDark2 checkout source, and a leftover orphaned copy from an older non-editable install was shadowing the corrected editable install even after reinstalling). If you hit `ValueError: mediator must be set to "light" or "heavy" to determine form of F_DM` when generating your own rates, check that `python3 -c "import qcdark2.dark_matter_rates as dm; print(dm.__file__)"` resolves to your actual QCDark2 checkout, not a stale `site-packages` copy — reinstall with `pip install -e .` from the QCDark2 repo root and delete any leftover non-editable `qcdark2/` directory under `site-packages` if the editable install doesn't take effect immediately.

---

## 5. Step 2 — Run the scans

```bash
build/ccdarksens_scan_generic configs/examples/dm_electron_intermediate_mA5keV_projection_1kgyear.json
build/ccdarksens_scan_generic configs/examples/dm_electron_intermediate_mA10keV_projection_1kgyear.json
```

Output: `outputs/dm_electron_intermediate_mA5keV_projection_1kgyear/scan_generic.root` and `outputs/dm_electron_intermediate_mA10keV_projection_1kgyear/scan_generic.root`.

**Runtime:** a few seconds each for the full 800×8 grid (measured directly).

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts all four curves — both mass slices, both exposures — on one canvas** (there is no published literature curve for the intermediate mediator to add, unlike the heavy/light cases). Run both LBC companion scans first ([`Student_Examples_DM_Electron_Intermediate_LBC.md`](Student_Examples_DM_Electron_Intermediate_LBC.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --intermediate \
  --out-pdf outplots/dm_electron_intermediate_combined.pdf \
  outputs/dm_electron_intermediate_mA5keV_lbc_1p3kgday/scan_generic.root "mA'=5 keV, LBC" \
  outputs/dm_electron_intermediate_mA5keV_projection_1kgyear/scan_generic.root "mA'=5 keV, 1 kg-yr" \
  outputs/dm_electron_intermediate_mA10keV_lbc_1p3kgday/scan_generic.root "mA'=10 keV, LBC" \
  outputs/dm_electron_intermediate_mA10keV_projection_1kgyear/scan_generic.root "mA'=10 keV, 1 kg-yr" \
  1.642374415149816 intermediate
```

`--intermediate` suppresses the (inapplicable) heavy/light literature overlays. `--from-qhist` forces all four curves to be rebuilt directly from each file's own q(mχ,σ) histogram, so they're computed identically and directly comparable to each other. This is the primary, recommended way to look at this example's output. A rendered copy is checked in at [`dm_electron_intermediate_combined.pdf`](dm_electron_intermediate_combined.pdf) — smooth curves across the full mass range, with no kink or truncation.

**Known cosmetic issue:** the plot's `F_DM` legend annotation is currently hard-coded to only two states in `apps/ccdarksens_plot_limit.cc` (`"heavy"` → `F_DM=1`, anything else → the light-mediator `1/q²`-type formula) — for `intermediate` mode it incorrectly shows the light-mediator formula. This is a display-only bug (does not affect the computed curve), not yet fixed since it requires a source change.

If you only want one slice's projection curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --intermediate \
  --out-pdf outplots/dm_electron_intermediate_mA5keV_projection_1kgyear.pdf \
  outputs/dm_electron_intermediate_mA5keV_projection_1kgyear/scan_generic.root \
  "DM-e intermediate mA'=5 keV" 1.642374415149816 intermediate
```

(swap in the mA10keV path/label for that slice.)

---

## 7. To use a different mA' value

Copy `configs/qcdark2_generate_Si_intermediate_mA5keV_comp.json`, change `"mediator"` to your desired mA' in eV, and pick a new `rates_dir`. QCDark2 accepts any positive float here — you are not limited to a fixed set of slices the way the earlier EXCEED-DM-based approach was (see §4's note on why this now works for arbitrary mA').

---

## 8. Superseded approach (for context, not for use)

An earlier version of this example used `python/ccdarkphys/exdm/entry.py` to convert EXCEED-DM `binned_scatter_rate` HDF5 output into rate CSVs (`data/exdm_rates/Si/intermediate/mA_{X}keV/`, still present on disk for slices 1/3/5/10/30/50/100 keV). That path had a known limitation — the available HDF5 was computed with EXCEED-DM's tutorial-scale 2×2×2 k-mesh, giving a numerically wrong absolute reach and an artificial kink above mχ ≳ 100 MeV — and only covered 7 fixed mA' values on a coarse 30-mass grid. QCDark2's native intermediate-mediator support (this doc) supersedes it: no k-mesh caveat, any mA' value, and the full 800-mass dense grid used everywhere else in this example set. The EXDM-based configs and data are left in place but are no longer used by the student examples.
