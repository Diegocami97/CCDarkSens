<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_DM_Electron_Intermediate_LBC.md -- DM-electron, intermediate-mass mediator (mA'=5,10 keV), real LBC background/counts reproduction.
-->

# Student Example — DM-Electron, Intermediate-Mass Mediator, LBC Reproduction

Real LBC background and observed counts, intermediate-mass-mediator case (mA' = 5 keV and 10 keV). Companion: [`Student_Examples_DM_Electron_Intermediate_Projection.md`](Student_Examples_DM_Electron_Intermediate_Projection.md).

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | DM-electron scattering, Si, **intermediate-mass mediator**, mA' = **5 keV and 10 keV** (QCDark2, composite Si_comp.h5 dielectric function) |
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Exposure | 1.3 kg-day (same LBC dataset as the heavy/light-mediator LBC cases) |
| Background | Real LBC `Bp`/`Br` (same values as the heavy/light cases) |
| Data | Real observed counts `[144, 0, 0, 1, 0, 0]` |
| Statistic | Profile likelihood ratio (pydme-style), 90% CL upper limit on σ̄ₑ(mχ) |

**Configs:**
- [`configs/examples/dm_electron_intermediate_mA5keV_lbc_1p3kgday.json`](../configs/examples/dm_electron_intermediate_mA5keV_lbc_1p3kgday.json)
- [`configs/examples/dm_electron_intermediate_mA10keV_lbc_1p3kgday.json`](../configs/examples/dm_electron_intermediate_mA10keV_lbc_1p3kgday.json)

Both were built by taking the heavy-mediator LBC config's background/data/pattern-efficiency block and swapping in the intermediate-mediator model block (`model.rates_dir` pointed at `data/qcdark2_rates/Si/intermediate/mA_5keV_comp/` or `mA_10keV_comp/`) — the LBC background template describes the detector and dataset, not the signal model, so it is identical across all mediator cases and both mass slices.

**Grid:** 800 masses × 8 cross sections = 6,400 rate points per slice (same dense grid as the projection case — see [`Student_Examples_DM_Electron_Intermediate_Projection.md`](Student_Examples_DM_Electron_Intermediate_Projection.md) §1 for the QCDark2/Si_comp.h5 rate-generation details).

---

## 2–3. Prerequisites and pipeline

Identical to [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) §2–3 — same `build/data_pattern.root` requirement, same fold → pattern-efficiency → Bp+θBr → profile-likelihood chain. Only the rate table differs (QCDark2 with a numeric intermediate mediator mass instead of `"heavy"`).

---

## 4. Step 1 — Rate generation (informational — already done, see companion doc)

No action needed; the `data/qcdark2_rates/Si/intermediate/mA_5keV_comp/` and `mA_10keV_comp/` CSVs are reused as-is from the projection example. See [`Student_Examples_DM_Electron_Intermediate_Projection.md`](Student_Examples_DM_Electron_Intermediate_Projection.md) §4 for how they were generated (including a note on a stale-`qcdark2`-install issue you may hit if regenerating).

---

## 5. Step 2 — Run the scans

```bash
build/ccdarksens_scan_generic configs/examples/dm_electron_intermediate_mA5keV_lbc_1p3kgday.json
build/ccdarksens_scan_generic configs/examples/dm_electron_intermediate_mA10keV_lbc_1p3kgday.json
```

Output: `outputs/dm_electron_intermediate_mA5keV_lbc_1p3kgday/scan_generic.root` and `outputs/dm_electron_intermediate_mA10keV_lbc_1p3kgday/scan_generic.root`.

**Runtime:** a few seconds each for the full 800×8 grid (measured directly).

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts all four curves — both mass slices, both exposures — on one canvas** (there is no published literature curve for the intermediate mediator to add). Run both projection companion scans first ([`Student_Examples_DM_Electron_Intermediate_Projection.md`](Student_Examples_DM_Electron_Intermediate_Projection.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --intermediate \
  --out-pdf outplots/dm_electron_intermediate_combined.pdf \
  outputs/dm_electron_intermediate_mA5keV_lbc_1p3kgday/scan_generic.root "mA'=5 keV, LBC" \
  outputs/dm_electron_intermediate_mA5keV_projection_1kgyear/scan_generic.root "mA'=5 keV, 1 kg-yr" \
  outputs/dm_electron_intermediate_mA10keV_lbc_1p3kgday/scan_generic.root "mA'=10 keV, LBC" \
  outputs/dm_electron_intermediate_mA10keV_projection_1kgyear/scan_generic.root "mA'=10 keV, 1 kg-yr" \
  1.642374415149816 intermediate
```

`--intermediate` suppresses the heavy/light literature overlays (not applicable here). `--from-qhist` forces all four curves to be rebuilt directly from each file's own q(mχ,σ) histogram, so they're computed identically and directly comparable to each other. This is the primary, recommended way to look at this example's output. A rendered copy is checked in at [`dm_electron_intermediate_combined.pdf`](dm_electron_intermediate_combined.pdf) — smooth curves across the full mass range, with no kink or truncation. See the projection doc's Step 6 for a known cosmetic `F_DM` legend-label issue in `--intermediate` mode.

If you only want one slice's LBC curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --intermediate \
  --out-pdf outplots/dm_electron_intermediate_mA5keV_lbc_1p3kgday.pdf \
  outputs/dm_electron_intermediate_mA5keV_lbc_1p3kgday/scan_generic.root \
  "DM-e intermediate mA'=5 keV LBC" 1.642374415149816 intermediate
```

(swap in the mA10keV path/label for that slice.)

---

## 7. To use a different mA' value

Same as the projection case: QCDark2 accepts any positive mA' value in eV — see [`Student_Examples_DM_Electron_Intermediate_Projection.md`](Student_Examples_DM_Electron_Intermediate_Projection.md) §7.

---

## 8. Superseded approach (for context, not for use)

This example previously used EXCEED-DM-converted rates (`data/exdm_rates/Si/intermediate/mA_{X}keV/`) with a known tutorial-k-mesh limitation. See the projection doc's §8 for details — QCDark2's native intermediate-mediator support supersedes it.
