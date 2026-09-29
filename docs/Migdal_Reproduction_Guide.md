# Migdal Effect Sensitivity Guide

A step-by-step walkthrough for collaboration members who want to project a **DAMIC-M Migdal-effect** sensitivity curve — DM-nucleus scattering off Si with the recoiling nucleus's Migdal-effect ionization providing the observable, extending sensitivity below the nuclear-recoil threshold where elastic scattering alone is invisible.

This guide assumes you have read the overview in [`Beginners_Guide.md`](Beginners_Guide.md). For the full physics scope, darkelf integration details, and validation checklist, see [`Migdal_Integration_Plan.md`](Migdal_Integration_Plan.md).

---

## 1. What you will reproduce

| Item | Value |
|------|-------|
| Observable | Six pattern bins folded from n_e (same detector pipeline as DM-electron) |
| Signal model | Migdal effect, Si target, free-nucleus approximation, ELF-based (`Si_mermin.dat` via darkelf) |
| Mediators | **Heavy** and **light** — same grid, different rate tables |
| Exposure | 1 kg · 1 year (Asimov projection) |
| Background | Flat dark-current migration, λ_DC = 3.65×10⁻³ e⁻/pix/img (DC = 10⁻⁵ e⁻/pix/day) |
| Data | Asimov (background-only pseudo-data) |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σₙ(mχ) |

**Example configs** (ready to run):

| Step | Heavy | Light |
|------|-------|-------|
| Rates | [`configs/examples/migdal_generate_si_heavy.json`](../configs/examples/migdal_generate_si_heavy.json) | [`configs/examples/migdal_generate_si_light.json`](../configs/examples/migdal_generate_si_light.json) |
| Scan | [`configs/examples/migdal_scan_si_heavy.json`](../configs/examples/migdal_scan_si_heavy.json) | [`configs/examples/migdal_scan_si_light.json`](../configs/examples/migdal_scan_si_light.json) |

Both use the same grid: **200 masses** (10–1000 MeV, log-spaced) × **200 cross sections** (10⁻⁴² – 10⁻³⁰ cm², log-spaced).

---

## 2. Prerequisites

### Software

- C++17 compiler, CMake ≥ 3.18
- [ROOT](https://root.cern/) (Core, Hist, RIO, Minuit2)
- [nlohmann/json](https://github.com/nlohmann/json)
- Python 3.8+ with NumPy, [darkelf](https://github.com/tongyanlin/DarkELF) (editable install) for rate generation

### Build CCDarkSens

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_generic ccdarksens_plot_limit
```

### Input data files

| File | Purpose |
|------|---------|
| `data/p100K_gap1p2_eh3p8.csv` | Charge-ionization yield table, Si at 100K, E_gap=1.2 eV, eh_pair=3.8 eV |
| `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv` | Precomputed P(pattern \| nₑ) efficiencies |

---

## 3. Pipeline overview

```
Step 1   Generate Migdal rate CSVs     (Python + darkelf, one-time per mediator/grid)
           ↓
Step 2   Run n_e-space scan            (C++, finds UL at each mχ)
           ↓
Step 3   Plot limit curve              (C++, σₙ axis)
```

---

## 4. Step 1 — Generate Migdal rate tables

```bash
python3 utils/migdal_generate_grid.py configs/examples/migdal_generate_si_heavy.json
python3 utils/migdal_generate_grid.py configs/examples/migdal_generate_si_light.json
```

Output directories: `data/migdal_rates/Si/heavy/`, `data/migdal_rates/Si/light/`

Filename pattern: `dRdE_Si28_{mediator}_m{mchi}_s{sigma}.csv`

### Notes

- **Grid size:** 200 × 200 = **40 000** CSV files per mediator. This is a full production grid — expect a long run. For a quick sanity check, temporarily reduce `grid.mchi_MeV.logspace.num` / `grid.sigma_n_cm2.logspace.num` to e.g. 5 in a copy of the config.
- **Verify:**

  ```bash
  ls data/migdal_rates/Si/heavy/dRdE_Si28_heavy_m100.000000_s1.0e-36.csv
  ```

  Confirm linear σₙ scaling as a fast correctness check (rates at 1e-36 and 1e-37 cm² should differ by exactly 10×).

The scan config's `model.rates_dir` and `model.grid` must match the rate-generation JSON — the example configs are already aligned.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/migdal_scan_si_heavy.json
build/ccdarksens_scan_generic configs/examples/migdal_scan_si_light.json
```

`ccdarksens_scan_generic` dispatches on `model.type: "migdal"` the same way it does for `dm_electron`/`dark_photon`/`wimp_nucleon` — the rest of the pipeline (charge ionization → pattern efficiency → flat DC-migration background → profile likelihood) is identical to a DM-electron n_e-space scan; only the rate CSVs and the y-axis quantity (σₙ instead of σₑ) differ.

### Output

| Path | Content |
|------|---------|
| `outputs/migdal/si_heavy_dc1e5/scan_generic.root` | Heavy-mediator results |
| `outputs/migdal/si_light_dc1e5/scan_generic.root` | Light-mediator results |

### Runtime expectation

Comparable to a DM-electron n_e-space projection scan at this grid size — expect on the order of tens of minutes to a few hours depending on machine, dominated by the 200-mass × 200-σₙ profile likelihood evaluation.

---

## 6. Step 3 — Plot the limit curve

```bash
build/ccdarksens_plot_limit \
  outputs/migdal/si_heavy_dc1e5/scan_generic.root "Migdal heavy" \
  outputs/migdal/si_light_dc1e5/scan_generic.root "Migdal light" \
  --migdal --batch --from-qhist --plain-legend \
  --out-pdf outplots/migdal_si_heavy_light.pdf
```

`--migdal` selects mχ vs. σₙ (DM-nucleon) axes. `--from-qhist` reads the q-histogram format from the scan ROOT output; `--plain-legend` suppresses the auto-generated mass/sigma annotation.

---

## 7. Scope boundaries

This example covers the **in-scope** baseline documented in `Migdal_Integration_Plan.md` §4.8: Si target, free-nucleus approximation, ELF-based Migdal via darkelf's `Si_mermin.dat`, n_e observable. The impulse approximation, Ibe et al. atomic Migdal, the elastic nuclear-recoil channel, and non-Si targets (including the SrCd₂Sb₂ configs elsewhere in `configs/`) are later, separately-scoped extensions — not part of this reproducible baseline.
