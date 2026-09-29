<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_Absorption_Projection.md -- Dark-photon absorption, real Si target, hypothetical 1 kg-year sensitivity projection.
-->

# Student Example — Dark Photon Absorption, Projection (Si)

Hypothetical 1 kg-year sensitivity projection for dark-photon absorption on real silicon, meant to be viewed alongside the LBC reproduction. Companion: [`Student_Examples_Absorption_LBC.md`](Student_Examples_Absorption_LBC.md) (plain Si, real LBC data).

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | Dark photon absorption, **Si** |
| Observable | n_e bins (extended to n_e ≤ 20) |
| Exposure | 1 kg · 1 year (hypothetical) |
| Background | Flat dark current (1e-5 e⁻/pixel/day) only, `norm_per_kg_year_keV: 0.0` |
| Data | Asimov (S = B) |
| Statistic | Profile likelihood ratio, 90% CL upper limit on ε(mA') |

**Config:** [`configs/examples/absorption_si_projection_1kgyear.json`](../configs/examples/absorption_si_projection_1kgyear.json)

**Grid:** 170 masses (mA' = 1–50 eV, explicit list) × 40 couplings (ε = 10⁻²⁰–10⁻¹⁰, log); reuses the LBC companion's rate table, `data/darkphoton_rates/Si_lbc_full/`.

This projection exists so the channel has the same "real material, hypothetical exposure, plotted next to the real-data LBC reproduction" pairing that every other channel in this example set has (DM-electron, Migdal). It deliberately reuses the LBC config's own rate table rather than generating a fresh one — both configs already share the same mass range, and the p100K ionization table caps how high in mass an n_e-space Si config can usefully go anyway (see §4).

---

## 2. Prerequisites

Build `ccdarksens_scan_generic` and `ccdarksens_plot_limit`. No new rate generation is needed — this config reuses the LBC companion's rate table (see §4).

---

## 3. Pipeline overview

```
Step 1   (none — reuses the LBC companion's rate CSVs)
           |
Step 2   Run the n_e-space Asimov scan (C++)
           |
Step 3   Plot the limit curve (C++)
```

**What the pipeline is physically doing:** for dark-photon absorption, the DM particle itself is absorbed (not scattered), depositing its full rest mass mA' as energy — the rate depends on the material's dielectric function Im[-1/ε(ω)] rather than a momentum-transfer form factor. The resulting dR/dE(E) is folded through the ionization table exactly as in the scattering channels, then combined with the same flat-DC Asimov background and profiled the same way; only the physical origin of the rate table and the final-state variable (mA', ε) instead of (mχ, σ) differ.

---

## 4. Step 1 — Rate generation (none needed — reuses the LBC table)

`model.rates_dir` points straight at `data/darkphoton_rates/Si_lbc_full/`, the same 200-mass × 40-ε table the LBC companion uses (see that doc §7.1). This config's own mass grid is restricted to the 170 values ≤ 50 eV (an explicit `grid.mA_eV.values` list, not a logspace), because the n_e-space response here uses the real-Si charge-ionization table `data/p100K_gap1p2_eh3p8.csv`, which only tabulates n_e-conversion probabilities up to E = 50 eV (20 n_e columns) — any higher mA' would deposit energy the table has no ionization probabilities for, which the code would (correctly, but silently) treat as no signal, the same failure mode as the `Emax_eV` bug in the LBC doc §7.2. `model.Emax_eV` is set to 55 eV (`nbins`=550) with the same 5 eV margin logic as that fix.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/absorption_si_projection_1kgyear.json
```

Output: `outputs/absorption_si_projection_1kgyear/scan_generic.root`.

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this projection on the same canvas as the LBC reproduction** (see the LBC doc §6 for the full command):

```bash
build/ccdarksens_plot_limit --batch --from-qhist --dark-photon \
  --out-pdf outplots/absorption_si_combined.pdf \
  outputs/absorption_si_lbc_1p3kgday/scan_generic.root "Si LBC (reproduced)" \
  outputs/absorption_si_projection_1kgyear/scan_generic.root "Si, 1 kg-yr projection" \
  1.642374415149816 heavy
```

`--dark-photon` switches the plot to (mA', ε) axes and loads the published dark-photon literature curves (XENON1T bracket, stellar cooling limits) from `data/previous_limits/dark_photon/` for a qualitative comparison. `--from-qhist` forces both curves to be rebuilt directly from each file's own q(mA',ε) histogram rather than the scan's own internally-stored UL graph, used consistently across every example in this set. A rendered copy is checked in at [`absorption_si_combined.pdf`](absorption_si_combined.pdf).

If you only want this projection's own curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --dark-photon \
  --out-pdf outplots/absorption_si_projection_1kgyear.pdf \
  outputs/absorption_si_projection_1kgyear/scan_generic.root \
  "Dark photon absorption, Si, 1 kg-yr projection" 1.642374415149816 heavy
```

**Why the curve has visible step/staircase structure, and why that's real, not a grid artifact.** This config's background is dark-current-only (`flat_background.norm_per_kg_year_keV: 0.0`) with a monochromatic (single-energy, single-n_e-bin) signal — a combination that makes the UL sensitive to exactly which discrete n_e bin the deposit lands in and to the pattern-classification thresholds (`thr_M`/`thr_MN`/`thr_MNL`, `Qmin_e`) near the low-energy edge. This was checked directly rather than assumed: an A/B test regenerating the same mass grid's rates at 40 vs. 200 ε points gave **pixel-identical** step positions, ruling out ε-grid resolution as the cause (unlike the genuine σ-grid-resolution bias found for Migdal, see that doc's LBC companion §4). The steps track real n_e-bin transitions in the ionization table (`E = band_gap_eV + n × eh_pair_eV`) and threshold effects near the ROI edge — structure, not noise.

---

## 7. Note

This example's `efficiency_csv` field already has the `"../data/..."` path fix described in [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) §7.
