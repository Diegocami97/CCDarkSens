# Example configs for collaboration workflows

## QEDark — LBC pydme reproduction (dense grid)

Guide: [`docs/LBC_QEDark_Reproduction_Guide.md`](../../docs/LBC_QEDark_Reproduction_Guide.md)

| Step | Heavy | Light (massless) |
|------|-------|------------------|
| Rates | `qedark_generate_heavy_dense.json` | `qedark_generate_light_dense.json` |
| Scan | `lbc_qedark_heavy_mediator.json` | `lbc_qedark_light_mediator.json` |

## QCDark2 — pattern-count sensitivity (dense grid)

Guide: [`docs/QCDark2_Pattern_Counts_Guide.md`](../../docs/QCDark2_Pattern_Counts_Guide.md)

| Step | Heavy | Light |
|------|-------|-------|
| Rates | `qcdark2_generate_heavy_dense.json` | `qcdark2_generate_light_dense.json` |
| Scan | `qcdark2_pattern_counts_heavy_mediator.json` | `qcdark2_pattern_counts_light_mediator.json` |

**Vary observed data:** edit `run.observed_counts` in the scan JSON (see `_comment_change_counts` at top of file). One entry per pattern in `experiment.pattern_roi` order `[11, 21, 111, 31, 22, 211]`.

Both workflow families use the same dense grid: **800 masses × 300 cross sections**.

## QCDark2 — n_e exposure projections (flat DC + flat d.r.u)

Guide: [`docs/QCDark2_NE_Exposure_Projections_Guide.md`](../../docs/QCDark2_NE_Exposure_Projections_Guide.md)

These are **Asimov** (sensitivity) projections in **n_e space** using the flat DC and flat d.r.u background model.

| Step | Config |
|------|--------|
| Rates (Si_comp, if missing) | `qcdark2_generate_si_comp_dense.json` |

| kg·year | Scan JSON |
|---------|-----------|
| 0.5     | `qcdark2_lbc_ne_flatbkg_proj_0p5kgy.json` |
| 1.0     | `qcdark2_lbc_ne_flatbkg_proj_1p0kgy.json` |
| 2.0     | `qcdark2_lbc_ne_flatbkg_proj_2p0kgy.json` |

**Vary exposure:** edit `detector.mass_kg` (or `experiment.livetime_days` / `duty_cycle`) and `run.outdir` — see `_comment_change_exposure` in the scan JSONs.
