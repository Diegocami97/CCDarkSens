# WIMP-Nucleon SI Reproduction Guide

A step-by-step walkthrough for collaboration members who want to reproduce the **DAMIC 2016 (0.6 kg-day, SNOLAB) WIMP-nucleon spin-independent exclusion limit** (PhysRevD.94.082006 / arXiv:1607.07410) using CCDarkSens's cluster-fit response chain and the joint 1×1 + 1×100 readout-channel likelihood.

This guide assumes you have read the overview in [`Beginners_Guide.md`](Beginners_Guide.md) and focuses on one concrete analysis path: **cluster-energy analysis space + joint 1×1/1×100 likelihood + `chavarria_izraelevitch_table` quenching**. For the full physics derivation and validation history, see [`ClusterFitMC_Design.md`](ClusterFitMC_Design.md).

---

## 1. What you will reproduce

| Item | Value |
|------|-------|
| Observable | Reconstructed cluster energy `E_reco` (per-channel Etrue→Ereco kernel) |
| Signal model | WIMP-nucleon SI recoil on Si, DAMIC 2016's own halo (v0=220, vE=232, vesc=544 km/s, ρ=0.3 GeV/cm³) |
| Quenching | `chavarria_izraelevitch_table` — measured Si nuclear-recoil yield (Chavarria 2016 below 2.28 keV_nr, Izraelevitch 2017 2.28–20.67 keV_nr, Lindhard fallback above) |
| Channels | Joint 1×1 (11×11 px window) + 1×100 (row-summed, 1D fit) readout, combined via `L_joint = L_1×1 × L_1×100` |
| Background | Flat Compton rate per channel (15 / 21 events·kg⁻¹·day⁻¹·keV⁻¹), each with its own digitized Fig. 9 detection-efficiency curve |
| Data | Asimov (background-only pseudo-data) — **not** a fit to the paper's real 31+23 candidate events |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σₙ(mχ) |

**Example configs** (ready to run):

| Step | Config |
|------|--------|
| Rates | [`configs/examples/wimp_nucleon_generate_si_damic2016halo_izraelevitch.json`](../configs/examples/wimp_nucleon_generate_si_damic2016halo_izraelevitch.json) |
| Kernel sanity check (1×1 only, optional) | [`configs/examples/wimp_nucleon_cluster_damic2016_repro.json`](../configs/examples/wimp_nucleon_cluster_damic2016_repro.json) |
| Joint scan | [`configs/examples/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch.json`](../configs/examples/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch.json) |

---

## 2. Prerequisites

### Software

- C++17 compiler, CMake ≥ 3.18
- [ROOT](https://root.cern/) (Core, Hist, RIO, Minuit2)
- [nlohmann/json](https://github.com/nlohmann/json)
- Python 3.8+ with NumPy (for rate generation)

### Build CCDarkSens

```bash
cmake -B build -S .
cmake --build build -j8 --target ccdarksens_scan_generic ccdarksens_validate_wimp_nucleon_paper_repro ccdarksens_plot_limit
```

### Input data files

These must already exist (they do, in a repo checkout) — no separate download needed:

| File | Purpose |
|------|---------|
| `data/wimp_nucleon_damic2016_fig9_background_1x1.csv` | Digitized background detection-efficiency curve, 1×1 channel (paper Fig. 9) |
| `data/wimp_nucleon_damic2016_fig9_background_1x100.csv` | Same, 1×100 channel |

---

## 3. Pipeline overview

```
Step 1   Generate WIMP-nucleon rate CSVs        (Python, one-time per halo/grid)
           ↓
Step 2   (Optional) Sanity-check the 1×1 kernel  (C++, ΔLL cut + K[Etrue,Ereco])
           ↓
Step 3   Run the joint 1×1+1×100 scan            (C++, finds UL at each mχ)
           ↓
Step 4   Plot the limit curve                    (C++, WIMP-nucleon axes)
```

---

## 4. Step 1 — Generate WIMP-nucleon rate tables

```bash
python3 utils/wimp_nucleon_generate_grid.py configs/examples/wimp_nucleon_generate_si_damic2016halo_izraelevitch.json
```

Output directory: `data/wimp_nucleon_rates_ee/Si/heavy_damic2016halo_izraelevitch/`

Filename pattern: `dRdE_ee_Si28_heavy_m{mchi}_s{sigma}.csv`

These CSVs already exist in a repo checkout, so this step is normally a no-op re-run (the generator skips files already on disk) — useful mainly if you're regenerating after a quenching-model or halo change.

**Verify:**

```bash
ls data/wimp_nucleon_rates_ee/Si/heavy_damic2016halo_izraelevitch/dRdE_ee_Si28_heavy_m10000.000000_s1.0e-33.csv
```

---

## 5. Step 2 — (Optional) Sanity-check the 1×1 kernel in isolation

Before running the full joint scan, you can inspect the 1×1 channel's noise-tail calibration and Etrue→Ereco kernel on their own — useful for debugging a detector-parameter change without paying for the full joint minimization:

```bash
build/ccdarksens_validate_wimp_nucleon_paper_repro configs/examples/wimp_nucleon_cluster_damic2016_repro.json \
  --dump-toys=outputs/wimp_nucleon_cluster_damic2016_repro/toys.csv \
  --dump-kernel=outputs/wimp_nucleon_cluster_damic2016_repro/kernel.csv
```

This prints the calibrated ΔLL cut (should extrapolate to ≈ −28.2, close to the paper's −28) and the marginal detection efficiency (should turn on around 11% at 60 eV and plateau near 75–79% above ~300 eV, matching paper Fig. 9's shape). `--dump-toys`/`--dump-kernel` are optional and only needed if you want to plot the raw ΔLL distribution or kernel heatmap yourself.

---

## 6. Step 3 — Run the joint 1×1 + 1×100 scan

```bash
build/ccdarksens_scan_generic configs/examples/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch.json
```

### What the scan does at each mass point

1. Load the shared `dR/dE_ee` CSV (both channels see the same true-energy spectrum).
2. Fold through each channel's own kernel → `S_1×1(E_reco)`, `S_1×100(E_reco)`.
3. Build each channel's own flat background × its own digitized efficiency curve → `B_1×1`, `B_1×100`.
4. Minimize each channel's NLL independently over its own background scale θ (the two channels don't share a nuisance parameter, only the signal strength) — this is exact, not an approximation, for `background_model: "scale"` (see `ClusterFitMC_Design.md` §6.7 for the derivation).
5. Sum the two channels' NLLs → `NLL_joint`, compute `q(σ) = 2ΔNLL_joint`.
6. Find σ_UL where q crosses **1.642** (90% CL, 1 dof).

### Output

| Path | Content |
|------|---------|
| `outputs/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch/scan_generic.root` | `TH2D h_q`, `TH1D h_upper_limit`, `TGraph g_upper_limit`, exposure metadata |

### Runtime expectation

Much lighter than the DM-electron pattern-space scans — the grid here is a handful of masses × a moderate σₙ range, and each point's kernel is pre-built once before the grid loop. Expect this to finish in minutes, not hours.

---

## 7. Step 4 — Plot the limit curve

```bash
build/ccdarksens_plot_limit \
  outputs/wimp_nucleon_damic2016_limit_repro_joint_izraelevitch/scan_generic.root "CCDarkSens joint repro" \
  --wimp --batch \
  --x-min 400 --x-max 12000 --y-min 1e-41 --y-max 1e-35 \
  --out-pdf outplots/wimp_nucleon_damic2016_joint_repro.pdf
```

`--wimp` selects mχ vs. σₙ axes and the WIMP-nucleon literature/paper overlay curves. `--x-min`/`--x-max`/`--y-min`/`--y-max` set the plot range in MeV / cm².

---

## 8. Known, documented caveats

This reproduction tracks the paper's real curve within ~1.5–3× through most of the mass range, and sits outside the paper's own published expected ±1σ band at nearly every mass — a real, unexplained residual that survived two genuine bug fixes made while building this (a folding-quadrature error that discarded ~90% of signal, and a background using the signal's own efficiency instead of its own). See `ClusterFitMC_Design.md` §6.10–§6.11 for the full consolidation and the ranked list of what to check next (quenching-uncertainty systematics is the top candidate).
