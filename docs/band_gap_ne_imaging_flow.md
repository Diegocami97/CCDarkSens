# Band-gap n_e imaging figures: flow & debugging guide

> Scope: the three "n_e imaging" figures for the band-gap pheno study and the
> exact data pipeline behind them. Written so the flow can be debugged
> end-to-end. Companion: [band_gap_pheno_ionization.md](band_gap_pheno_ionization.md),
> [band_gap_study_architecture.md](band_gap_study_architecture.md).

All outputs live under `outplots/band_gap_pheno/ne_imaging/`.
Reference DM point for the data figures: **m_χ = 1.007781 MeV** (closest grid
mass to 1 MeV), **σ̄_e = 1.1e-35 cm²**, **heavy** mediator, exposure **0.5 kg·yr**
(`mass_kg=0.5`, `livetime_days=365.25`).

---

## 0. TL;DR repro

```bash
cd /Users/diegovenegasvargas/Documents/CCDarkSens

# Figure 1 (analytic, no data needed)
python3 utils/plot_ne_vs_Er_v2.py

# Figure 3 (analytic, uses diffusion params only, no data needed)
python3 utils/plot_cluster_visualization.py

# Figure 2 data (one n_e-space scan run per scenario; ~25 s each, 3 runs total):
# Each dump contains S_true_ne__<tag>, S_obs_ne__<tag> (= S_true·ε(n_e)), and B_tot_ne.
for c in si_ref gap0p5_eh1p0 gap0p1_eh0p5; do
  build/ccdarksens_scan_dmelectron_pattern configs/ne_imaging_one_point_${c}.json
done

# Figure 2 (reads the ROOT files above)
python3 utils/plot_ne_signal_spectrum_compare.py
```

---

## 1. The three figures

| Fig | Script | Needs ROOT data? | Output |
|-----|--------|------------------|--------|
| 1 | `utils/plot_ne_vs_Er_v2.py` | No (pure analytic) | `ne_vs_Er_conventions.pdf` |
| 2 | `utils/plot_ne_signal_spectrum_compare.py` | **Yes** (uproot) | `ne_signal_spectrum_compare.pdf` |
| 3 | `utils/plot_cluster_visualization.py` | No (analytic + diffusion params) | `cluster_visualization.pdf` |

Each script also writes a `.png` alongside the `.pdf`.

### Fig 1 — ⟨n_e⟩ vs E_r (two conventions)
Pure analytic. Two panels:
- Convention 1: `⟨n_e⟩ = (E_r − E_gap)/ε_h` for `E_r ≥ E_gap`, else 0.
- Convention 2: `⟨n_e⟩ = E_r/ε_h` for `E_r ≥ E_gap`, else 0.

Scenarios are hard-coded in `SCENARIOS` (Si 1.2/3.8, 0.5/1.0, 0.1/0.5, 0.1/0.1).
No framework dependency — safe to edit freely.

### Fig 3 — pixel cluster visualization
Pure analytic, but uses the **real** CCDarkSens diffusion model so the cluster
shapes match the framework. Formula mirrors
`include/ccdarksens/response/DiffusionPhysics.hh::ComputeSigmaXYUm`:

```
σ_xy(z,E) = sqrt(-A · log(1 - b·z)) · (α + β·E_keV)     [µm]
```

Parameters (hard-coded, match the configs / C++ defaults):
`A=803.25 µm²`, `b=0.00065 µm⁻¹`, `α=1.0`, `β=0`, `σ_ro=0.16 e⁻`, pitch `15 µm`,
thickness `0.67 mm`. For the 4 eV event: σ_xy = 5.15 / 14.09 / 20.36 µm at
z = 50 / 337 / 620 µm.

Pipeline per panel: gaussian charge cloud (⟨n_e⟩ electrons) → integrate over
pixels (erf, separable) → add `N(0, σ_ro)` readout noise → round to integers.
RNG is seeded (`RNG_SEED=20260609`) so panels are reproducible.

### Fig 2 — S(n_e) signal spectrum compare
Reads existing ROOT outputs (see §2–§4). Two panels:
- Left: `S_true(n_e)` — pure ionization (`ChargeIonization::FoldToNe`).
- Right: `S_obs(n_e) = S_true(n_e) · ε(n_e)` — after the **single-pixel** efficiency
  ε(n_e) (the correct ROI for the n_e analysis space).

Both panels overlay `B_tot(n_e)` (dark current + flat, true n_e) from the same
run. Bars are log-scale; values below `LOG_FLOOR=1e-3` are not drawn. All three
spectra come from **one** scan dump (see §2).

---

## 2. Data pipeline for Figure 2

A **single** app run per scenario produces all three spectra in one
self-consistent n_e-space file. (Earlier this used a two-app workaround; that is
no longer needed after the efficiency/ROI fix — see §5.)

```
                 configs/ne_imaging_one_point_<id>.json          (observable_bins = "ne")
                          │  dump_point_spectra_root = true
                          ▼
   build/ccdarksens_scan_dmelectron_pattern
                          │
                          ▼
   outputs/band_gap_ne_imaging/<id>/scan_dmelectron_pattern.root
       ├─ dRdE__<tag>
       ├─ S_true_ne__<tag>      ◄── Fig 2 LEFT panel (ionization signal)
       ├─ S_obs_ne__<tag>       ◄── Fig 2 RIGHT panel (= S_true · ε(n_e), single-pixel eff)
       └─ B_tot_ne              ◄── Fig 2 background overlay (true n_e: DC + flat)
```

`<id>` ∈ {`si_ref`, `gap0p5_eh1p0`, `gap0p1_eh0p5`}.
`<tag>` = `mchi_1p007781__sigma_1p1em35` (sanitized point tag).

---

## 3. Scenarios and inputs

| id | (E_gap, ε_h) eV | rates_dir | ionization table |
|----|------------------|-----------|------------------|
| `si_ref` | (1.2, 3.8) | `data/qcdark2_rates/Si/heavy/Si_fast_gap1p2` | `data/p100K_gap1p2_eh3p8.csv` |
| `gap0p5_eh1p0` | (0.5, 1.0) | `.../Si_fast_gap0p5` | `data/p100K_gap0p5_eh1p0.csv` |
| `gap0p1_eh0p5` | (0.1, 0.5) | `.../Si_fast_gap0p1` | `data/p100K_gap0p1_eh0p5.csv` |

Rate file resolved per point: `dRdE_Si_heavy_m1.007781_s1.1e-35.csv`
(`m` formatted `%.6f` via `DMElectronModel::format_mchi_6f`; `s` via
`grid.format.sigma = ".1e"`).

Config differences are only: `run.label`, `run.outdir`, `model.rates_dir`,
`response.charge_ionization.{table_csv,band_gap_eV,eh_pair_eV,scenario}`, and
`experiment.observable_bins` (`ne` for scan, `pattern` for example).

---

## 4. Histogram reference (what to inspect when debugging)

Read with `uproot`. Map bins → n_e by **axis centers** (do NOT assume
`values()[0]` is n_e=1; the histograms include an n_e=0 bin):

```python
import uproot, numpy as np
def by_ne(h, ne_lo=1, ne_hi=5):
    c = h.axis().centers(); v = h.values()
    return [float(v[int(np.argmin(np.abs(c-ne)))]) for ne in range(ne_lo, ne_hi+1)]
```

Expected values at the reference point (sanity baseline):

| id | S_true[1..5] | S_obs[1..5] (= S_true·ε) | B_tot[1..5] |
|----|--------------|--------------------------|-------------|
| si_ref | 1863, 1.7, 0, 0, 0 | 1762, 0.6, 0, 0, 0 | 2.77e5, 0.10, 2.8e-4, 1.5e-4, 7.3e-5 |
| gap0p5_eh1p0 | 5051, 3197, 1423, 140, 2.3 | 4776, 1082, 220, 12.3, 0.14 | 2.77e5, 0.10, 7.5e-5, 4.4e-5, 2.9e-5 |
| gap0p1_eh0p5 | 2336, 3513, 5933, 5144, 3482 | 2209, 1190, 919, 452, 205 | 2.77e5, 0.10, 3.8e-5, 2.2e-5, 1.5e-5 |

The `S_obs/S_true` ratio reproduces the **single-pixel diagonal efficiency**
`ε(n_e) ≈ [0.95, 0.34, 0.15, 0.09, 0.06]` (= `P(single-pixel{n_e}|n_e) ·
Eff_csv(n_e,n_e)`), decreasing with n_e because more electrons spread out of one
pixel.

Physics check: the **n_e=1 background is ~2.77e5 (dark current)**; it falls to
~0.1 at n_e=2 and ~1e-4 at n_e≥3. Even after the (decreasing) single-pixel
efficiency, the low-gap signal still populates n_e=2–5 where the background is
negligible → that is the discrimination message.

> The huge `1.4e11` you may see at `values()[0]` is the **n_e=0** dark-current
> bin (outside ROI), not n_e=1. Use `by_ne` to avoid this trap.

---

## 5. Gotchas (the non-obvious stuff that caused detours)

1. **Efficiency ROI must follow the analysis space (fixed bug).**
   `EfficiencyMC` builds `ε(n_e)` by summing over `emc_cfg.accepted_labels`.
   Previously both apps built `accepted_labels` **only** from
   `experiment.pattern_roi = [11,21,111,31,22,211]` regardless of
   `observable_bins`. In `ne` mode that folded the n_e observable against the
   *multi-pixel* ROI, which excludes the single-pixel patterns the n_e space
   actually is → `ε(n_e) ≈ 0` and `S_obs(n_e) ≈ 0` (wrong).

   Fix (in `ccdarksens_scan_dmelectron_pattern.cc` and
   `ccdarksens_example_one_point_pattern.cc`): `accepted_labels` now follows the
   space —
   - `pattern` mode → labels from `pattern_roi` (multi-pixel), unchanged.
   - `ne` mode → labels from `roi_bins` as **single-pixel** patterns
     (`lab.q = {ne}`, one pixel holding n_e electrons).

   Result: `ε(n_e) = P(single-pixel{n_e}|n_e) · Eff_csv(n_e,n_e)` → the sensible
   decreasing diagonal `[~0.95, 0.34, 0.15, 0.09, 0.06]`. This affects **all
   `observable_bins="ne"` results** (including n_e-space limit curves), which
   were previously suppressed. Pattern-mode results are unchanged.

   > The single-pixel efficiency itself was always fine (~0.98 at n_e=1 in
   > `data/efficiencies_paolo.csv` and
   > `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`). The earlier
   > suppression was the ROI mismatch, not a tiny efficiency.

2. **Removed hard-coded `h_eps_ne` override in the example app.**
   `ccdarksens_example_one_point_pattern.cc` used to overwrite the computed
   `h_eps_ne` with placeholder values `[1.0, 0.38, 0.65, 0.79, 0.86]` right after
   computing it. That masked bug #1 (and was physically backwards — efficiency
   increasing with n_e). It has been removed; the app now uses the real computed
   single-pixel efficiency. Fig 2 no longer depends on this app.

3. **`S_true` is the same in both apps** (it is `ion->FoldToNe(dRdE, exposure)`),
   so Fig 2 takes it from the scan dump. If scan `S_true` and an example-app
   `S_true` ever disagree, the rate file or exposure differs.

4. **Background units / n_e=0 bin.** `B_tot` is dominated by the n_e=0 and n_e=1
   dark-current bins. Always select by axis center, not raw index.

5. **Mass formatting.** The point must resolve to an existing rate file. `m` is
   `%.6f`; passing `1.0` on the CLI yields `m1.000000` (no such file). Use the
   exact grid mass `1.007781`.

6. **Matplotlib under sandbox.** Running the plot scripts inside a restricted
   sandbox can segfault on font-cache writes. Run them normally (outside the
   sandbox) if that happens.

---

## 6. Files created for this test

Configs (new; existing scan configs untouched). The `_pattern.json` variants are
no longer needed for Fig 2 (kept for reference only):
- `configs/ne_imaging_one_point_si_ref.json`
- `configs/ne_imaging_one_point_gap0p5_eh1p0.json`
- `configs/ne_imaging_one_point_gap0p1_eh0p5.json`

Scripts:
- `utils/plot_ne_vs_Er_v2.py`
- `utils/plot_ne_signal_spectrum_compare.py`
- `utils/plot_cluster_visualization.py`

C++ changed (efficiency/ROI fix, see §5):
- `apps/ccdarksens_scan_dmelectron_pattern.cc`
- `apps/ccdarksens_example_one_point_pattern.cc`

Generated data (one file per scenario, contains all three spectra):
- `outputs/band_gap_ne_imaging/<id>/scan_dmelectron_pattern.root`

Figures:
- `outplots/band_gap_pheno/ne_imaging/ne_vs_Er_conventions.{pdf,png}`
- `outplots/band_gap_pheno/ne_imaging/ne_signal_spectrum_compare.{pdf,png}`
- `outplots/band_gap_pheno/ne_imaging/cluster_visualization.{pdf,png}`

---

## 7. Debugging checklist

- [ ] Rate file exists: `ls data/qcdark2_rates/Si/heavy/Si_fast_gap0p5/dRdE_Si_heavy_m1.007781_s1.1e-35.csv`
- [ ] Ionization table exists: `ls data/p100K_gap0p5_eh1p0.csv`
- [ ] Binary built: `ls build/ccdarksens_scan_dmelectron_pattern`
- [ ] Scan ROOT has `S_true_ne__*`, `S_obs_ne__*`, `B_tot_ne`: `python3 -c "import uproot;print(uproot.open('outputs/band_gap_ne_imaging/si_ref/scan_dmelectron_pattern.root').keys())"`
- [ ] `S_obs/S_true` ≈ single-pixel diagonal `[0.95, 0.34, 0.15, 0.09, 0.06]`
- [ ] Values match §4 baseline (use `by_ne`)
- [ ] If a bar is missing in Fig 2: it is < `LOG_FLOOR` (1e-3), expected for empty bins
