# Dark Photon Absorption Reproduction Guide

A step-by-step walkthrough for collaboration members who want to project a **dark photon absorption** sensitivity curve for the hypothetical material **SrCd₂Sb₂ (HypMat)**, using a Drude-model energy-loss function.

This guide assumes you have read the overview in [`Beginners_Guide.md`](Beginners_Guide.md). For the full physics derivation (absorption rate formula, ELF construction, the Drude-model assumptions and their brackets), see [`DarkPhoton_Absorption_HypotheticalMaterial.md`](DarkPhoton_Absorption_HypotheticalMaterial.md).

---

## 1. What you will reproduce

| Item | Value |
|------|-------|
| Observable | n_e electron-count spectrum (40 bins) |
| Signal model | Dark photon absorption, SrCd₂Sb₂ target, Drude ELF (ε∞=1, "unscreened" bracket) |
| Material | Density 5.76 g/cm³, direct gap 2.1 eV (scissors-shifted), charge-ionization gap 0.34 eV, eh_pair=1.7778 eV |
| Exposure | 1 kg · 1 year (Asimov projection) |
| Background | Flat dark-current migration, DC = 10⁻⁵ e⁻/pix/day |
| Data | Asimov (background-only pseudo-data) |
| Statistic | Profile likelihood ratio (pydme-style), 90% CL upper limit on ε(m_A') |

**Example configs** (ready to run):

| Step | Config |
|------|--------|
| Rates | [`configs/examples/darkphoton_generate_hypmat_unscreened.json`](../configs/examples/darkphoton_generate_hypmat_unscreened.json) |
| Scan | [`configs/examples/darkphoton_scan_hypmat_unscreened_ne.json`](../configs/examples/darkphoton_scan_hypmat_unscreened_ne.json) |

This is one of three coexisting HypMat model families (Drude unscreened, Drude screened ε∞=12, QCDark2-ELF proxy) — this guide covers the unscreened bracket only. See `configs/darkphoton_generate_hypmat_screened.json` / `configs/darkphoton_generate_hypmat_qcdark2.json` for the other two if you need them; they follow the same pattern.

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
| `data/p100K_gap0p34_eh1p78.csv` | Charge-ionization yield table for the 0.34 eV gap material |
| `data/darkphoton_rates/HypMat_0p34_unscreened/` | Drude ELF synthetic data + generated rate CSVs |

If the Drude ELF input files are missing, regenerate them first:

```bash
python3 utils/build_darkphoton_hypothetical_material.py --darkelf_dir <path-to-DarkELF>
```

---

## 3. Pipeline overview

```
Step 1   Generate dark-photon rate CSVs     (Python + darkelf, one-time per model/grid)
           ↓
Step 2   Run n_e-space scan                 (C++, finds UL at each m_A')
           ↓
Step 3   Plot limit curve                   (C++, ε vs. m_A' axes)
```

---

## 4. Step 1 — Generate dark photon rate tables

```bash
python3 utils/darkphoton_generate_grid.py configs/examples/darkphoton_generate_hypmat_unscreened.json
```

Output directory: `data/darkphoton_rates/HypMat_0p34_unscreened/`

Filename pattern: `dRdE_hypmat_unscreened_absorption_m{mA_eV}_e{epsilon}.csv`

### Notes

- **Grid size:** 300 masses (0.1–100 eV, log-spaced) × 40 couplings (10⁻²⁰ – 10⁻⁸, log-spaced) = **12 000** CSV files.
- **Physics note:** the Drude model here is overdamped (γ=61.7 eV ≫ ωp=11.7 eV, from σ_DC=300 Ω⁻¹cm⁻¹), so this bracket gives a featureless curve with no sharp plasmon resonance — that's expected, not a bug. The QCDark2-ELF proxy family (`darkphoton_generate_hypmat_qcdark2.json`) uses Si's measured optical dielectric as a stand-in and shows a sharp plasmon feature instead, for comparison.
- **Verify:**

  ```bash
  ls data/darkphoton_rates/HypMat_0p34_unscreened/dRdE_hypmat_unscreened_absorption_m1.000000_e1.0e-14.csv
  ```

The scan config's `model.rates_dir` and `model.grid` must match the rate-generation JSON — the example configs are already aligned.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/darkphoton_scan_hypmat_unscreened_ne.json
```

`ccdarksens_scan_generic` dispatches on `model.type: "dark_photon"` the same way it does for `dm_electron`/`migdal`/`wimp_nucleon` — charge ionization uses the material's own gap/eh_pair (0.34 eV / 1.7778 eV), not silicon's.

### Output

| Path | Content |
|------|---------|
| `outputs/darkphoton/hypmat_unscreened_ne/scan_generic.root` | Scan results |

### Runtime expectation

Comparable to a DM-electron n_e-space projection scan at this grid size (300×40 points) — expect on the order of tens of minutes.

---

## 6. Step 3 — Plot the limit curve

```bash
build/ccdarksens_plot_limit \
  outputs/darkphoton/hypmat_unscreened_ne/scan_generic.root "SrCd2Sb2 unscreened (DC=1e-5)" \
  --dark-photon --plain-legend --batch --from-qhist \
  --out-pdf outplots/darkphoton_hypmat_unscreened.pdf
```

`--dark-photon` selects ε vs. m_A' axes and loads the stellar-cooling/direct-detection literature overlay curves. `--from-qhist` reads the q-histogram format from the scan ROOT output; `--plain-legend` suppresses the auto-generated mass/coupling annotation.

To compare several dark-current levels on one plot, pass multiple `(root, label)` pairs — see `configs/darkphoton_scan_hypmat_unscreened_ne_dc1e2.json` / `_dc1e3.json` for the other DC-level configs referenced elsewhere in `configs/`.

---

## 7. Open items

The true SrCd₂Sb₂ energy-loss function awaits ab-initio wavefunction data; both Drude brackets (ε∞=1, ε∞=12) and the QCDark2-ELF proxy are stand-ins bracketing the true answer, not a measured result. See `DarkPhoton_Absorption_HypotheticalMaterial.md` §9 for the full caveat discussion.
