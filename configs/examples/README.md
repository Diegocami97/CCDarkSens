<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

README.md -- Index of the example configs in this folder: for every physics case (QEDark and QCDark2 DM-electron scans, Migdal, dark photon, WIMP-nucleon) it lists which JSON generates the rate tables and which JSON runs the scan, with a pointer to the matching guide in docs/.
-->

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

## WIMP-Nucleon SI — DAMIC 2016 joint 1×1+1×100 reproduction

Guide: [`docs/WIMP_Nucleon_Reproduction_Guide.md`](../../docs/WIMP_Nucleon_Reproduction_Guide.md)

| Step | Config |
|------|--------|
| Rates | `wimp_nucleon_generate_si_damic2016halo_izraelevitch.json` |
| Kernel sanity check (1×1 only, optional) | `wimp_nucleon_cluster_damic2016_repro.json` |
| Joint scan | `wimp_nucleon_damic2016_limit_repro_joint_izraelevitch.json` |

Cluster-energy analysis space, not pattern/n_e — this is the WIMP-nucleon channel's own detector-response pipeline (per-event cluster fit → noise-tail-calibrated kernel), combining the 1×1 and 1×100 readout channels via a joint likelihood.

## Migdal Effect — Si heavy/light mediator projection

Guide: [`docs/Migdal_Reproduction_Guide.md`](../../docs/Migdal_Reproduction_Guide.md)

| Step | Heavy | Light |
|------|-------|-------|
| Rates | `migdal_generate_si_heavy.json` | `migdal_generate_si_light.json` |
| Scan | `migdal_scan_si_heavy.json` | `migdal_scan_si_light.json` |

Both use the same grid: **200 masses × 200 cross sections**. n_e-space, 1 kg·yr Asimov projection, DC=10⁻⁵ e⁻/pix/day.

## Dark Photon Absorption — SrCd₂Sb₂ (HypMat) unscreened projection

Guide: [`docs/DarkPhoton_Reproduction_Guide.md`](../../docs/DarkPhoton_Reproduction_Guide.md)

| Step | Config |
|------|--------|
| Rates | `darkphoton_generate_hypmat_unscreened.json` |
| Scan | `darkphoton_scan_hypmat_unscreened_ne.json` |

n_e-space, 1 kg·yr Asimov projection, DC=10⁻⁵ e⁻/pix/day. One of three coexisting HypMat model families (Drude unscreened/screened, QCDark2-ELF proxy) — see the guide for the other two.

## Student examples — full channel set (Projection + LBC per channel)

Fifteen self-contained configs handed off to the group's students, one pair (Projection = hypothetical flat-DC background, LBC = real DAMIC-M pattern-space background and observed counts) per DM-electron mediator case, plus a Migdal pair and an absorption pair, plus two WIMP-nucleon cases. Each has its own step-by-step doc under `docs/Student_Examples_*.md` covering rate generation, the scan step, and the plotting step with literature overlays.

| Channel | Rate-generation config(s) | Projection config | Projection guide | LBC config | LBC guide |
|---|---|---|---|---|---|
| DM-e, heavy mediator | `qedark_generate_heavy_dense.json` | `dm_electron_heavy_projection_1kgyear.json` | [`Student_Examples_DM_Electron_Heavy_Projection.md`](../../docs/Student_Examples_DM_Electron_Heavy_Projection.md) | `dm_electron_heavy_lbc_1p3kgday.json` | [`Student_Examples_DM_Electron_Heavy_LBC.md`](../../docs/Student_Examples_DM_Electron_Heavy_LBC.md) |
| DM-e, light (ultralight) mediator | `qedark_generate_light_dense.json` | `dm_electron_light_projection_1kgyear.json` | [`Student_Examples_DM_Electron_Light_Projection.md`](../../docs/Student_Examples_DM_Electron_Light_Projection.md) | `dm_electron_light_lbc_1p3kgday.json` | [`Student_Examples_DM_Electron_Light_LBC.md`](../../docs/Student_Examples_DM_Electron_Light_LBC.md) |
| DM-e, intermediate mediator (mA'=5 and 10 keV) | `../qcdark2_generate_Si_intermediate_mA5keV_comp.json`, `../qcdark2_generate_Si_intermediate_mA10keV_comp.json` (in `configs/`, not `configs/examples/`) | `dm_electron_intermediate_mA5keV_projection_1kgyear.json`, `dm_electron_intermediate_mA10keV_projection_1kgyear.json` | [`Student_Examples_DM_Electron_Intermediate_Projection.md`](../../docs/Student_Examples_DM_Electron_Intermediate_Projection.md) | `dm_electron_intermediate_mA5keV_lbc_1p3kgday.json`, `dm_electron_intermediate_mA10keV_lbc_1p3kgday.json` | [`Student_Examples_DM_Electron_Intermediate_LBC.md`](../../docs/Student_Examples_DM_Electron_Intermediate_LBC.md) |
| Migdal (heavy + light, one doc per case) | LBC: `migdal_generate_si_heavy_150x300.json`/`_light_150x300.json` (in `configs/`, not `configs/examples/`). Projection (reduced demo grid): `migdal_generate_si_heavy_demo150x200.json`, `migdal_generate_si_light_demo150x200.json` | `migdal_heavy_projection_1kgyear.json`, `migdal_light_projection_1kgyear.json` | [`Student_Examples_Migdal_Projection.md`](../../docs/Student_Examples_Migdal_Projection.md) | `migdal_heavy_lbc_1p3kgday.json`, `migdal_light_lbc_1p3kgday.json` | [`Student_Examples_Migdal_LBC.md`](../../docs/Student_Examples_Migdal_LBC.md) |
| Dark photon absorption | `darkphoton_generate_si_lbc_full.json` (also reused, mA'≤50 eV subset, by the projection) | `absorption_si_projection_1kgyear.json` (Si, reuses the LBC rate table, no separate generate config) | [`Student_Examples_Absorption_Projection.md`](../../docs/Student_Examples_Absorption_Projection.md) | `absorption_si_lbc_1p3kgday.json` (Si) | [`Student_Examples_Absorption_LBC.md`](../../docs/Student_Examples_Absorption_LBC.md) |

| WIMP-nucleon case | Config(s) | Guide |
|---|---|---|
| DAMIC at SNOLAB (2016) reproduction, baseline + LEE (1×1-only, see limitations below) | `wimp_damic_snolab_2016.json`, `wimp_damic_snolab_2016_lee.json` | [`Student_Examples_WIMP_DAMIC_SNOLAB.md`](../../docs/Student_Examples_WIMP_DAMIC_SNOLAB.md) |
| DAMIC-M projection, floor60, 1.3 kg-day (LBC-exposure) + 1 kg-year | `wimp_damicm_lbc_1p3kgday.json`, `wimp_damicm_projection_floor60_1kgyear.json` | [`Student_Examples_WIMP_DAMIC_M_Projection.md`](../../docs/Student_Examples_WIMP_DAMIC_M_Projection.md) |

**LBC exposure** across all pattern-space configs above: 1.3 kg-day (`mass_kg=0.01523 × livetime_days=85.356`), same detector/dataset for every mediator/channel — the LBC background template describes the detector and dataset, not the signal model.

**Known limitations carried by specific configs in this set** (see each config's own doc for the full explanation — not repeated here):
- Migdal projection rates use a reduced 150-mass-point demo grid (`*_demo150x200/`) instead of a full-resolution grid — DarkELF has no caching, so a denser mass grid was impractically slow to generate; the cross-section axis is *not* reduced (200 points) for the reason below. Mass range is 0.1–1000 MeV, covering the 1–35 MeV band DAMIC-M's own PRL reports as this channel's most stringent-limit region (an earlier iteration inherited a 10 MeV floor that cut off part of that range — fixed). The Migdal LBC pair uses the real, complete 150×300 grids.
- All Migdal grids (both projection and LBC) use at least 200-300 cross-section points, not the smaller counts an earlier iteration used (30-40) — `--from-qhist` linearly interpolates the upper limit between adjacent cross-section grid points, so a coarse cross-section grid produces both a visibly faceted **and a numerically biased** curve (confirmed against the published DAMIC-M heavy-mediator curve: a 30-point cross-section grid averaged ~20% too strong vs. the reference; 300 points brought that to ~8% with a much tighter spread). Cross-section resolution is nearly free to increase (rate is computed once per mass, then rescaled exactly for every cross-section value), so there's no reason to use a coarse cross-section grid for any Migdal channel.
- The dark-photon LBC rates were regenerated fresh into `data/darkphoton_rates/Si_lbc_full/` after the pre-existing `Si/` directory was found to hold wrong-coupling-range files; the old directory and the pre-existing `darkphoton_scan_si_lbc.json` were left untouched.
- The dark-photon LBC config's `model.Emax_eV` was found set to 20 eV while its mass grid runs to 100 eV — any mA' above 20 eV had its rate CSV's energy content fall entirely outside the declared histogram range and got silently dropped (zero signal, no warning), pinning the UL to the coupling-grid edge for that whole region. Fixed by raising `Emax_eV` to 105 eV (`nbins` scaled to match). The Si projection config, which reuses the LBC rate table but restricts itself to mA' ≤ 50 eV, sets its own `Emax_eV` to 55 eV for the same reason.
- A shared-code bug in `src/response/PatternRates.cc` (`FoldNeToPatternRates`) forced full (100%) detection efficiency into *every* declared pattern bin for any event depositing `n_e >= 10` electrons, instead of correctly treating an untabulated, out-of-range cluster shape as zero acceptance. This spuriously improved the dark-photon-absorption LBC upper limit in the mA'≈24-35 eV range (where the Fano-broadened n_e distribution's tail first reaches n_e=10); a duplicate of the same hack in `src/response/ResponseFactory.cc` was also removed, though a direct A/B check confirmed it had no effect on the n_e-space projection examples. DM-electron heavy/light/Migdal LBC examples never populate n_e>=10 and are numerically unaffected (verified bit-identical before/after); the DM-electron **intermediate**-mediator LBC examples do populate it and were measurably affected (mA'=5 keV: ~0.3% shift; mA'=10 keV: ~11% shift) — both were rerun and their checked-in PDF updated. See [`Student_Examples_Absorption_LBC.md`](../../docs/Student_Examples_Absorption_LBC.md) §7.3 for the full writeup.
- The WIMP-nucleon SNOLAB baseline/+LEE configs were originally a joint 1×1+1×100 likelihood, but `RunJointChannelScan` (`apps/ccdarksens_scan_generic.cc`) never wires the `backgrounds.lee` low-energy-excess term into either channel's background construction (only the single-channel `MakeBackground()` path does) — the +LEE joint config produced bit-identical upper limits to its baseline sibling at all 28 mass points, the tell that the term wasn't applied at all. Scope reduced to 1×1-only for now (which already applies LEE correctly); a joint-channel fix is a known, deferred follow-up. See [`Student_Examples_WIMP_DAMIC_SNOLAB.md`](../../docs/Student_Examples_WIMP_DAMIC_SNOLAB.md) for the full writeup.
- All configs under this set that use `efficiency_mc.efficiency_csv` apply a `"../data/..."` prefix to work around a pre-existing path-resolution bug in `ResolveConfigRelativePath` (`src/response/ResponseFactory.cc`) that only manifests for configs one directory deeper than `configs/`. The pre-existing `lbc_qedark_heavy_mediator.json`, `lbc_qedark_light_mediator.json`, `migdal_scan_si_heavy.json`, and `migdal_scan_si_light.json` above still have the old, unfixed path and will fail until the same one-line fix is applied to them (or the resolver itself is fixed in a separate, approved change).
