<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_Migdal_LBC.md -- Migdal effect, heavy + light mediator, real LBC pattern-space background/counts reproduction.
-->

# Student Example — Migdal Effect, LBC Reproduction (Heavy + Light Mediator)

Real LBC pattern-space background and observed counts, Migdal-effect channel, both mediator cases. Companion: [`Student_Examples_Migdal_Projection.md`](Student_Examples_Migdal_Projection.md).

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | Migdal effect, Si, both **heavy** and **light** mediator |
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Exposure | 1.3 kg-day (same LBC dataset as the DM-electron LBC examples) |
| Background | Real LBC `Bp`/`Br` |
| Data | Real observed counts `[144, 0, 0, 1, 0, 0]` |
| Statistic | Profile likelihood ratio (pydme-style), 90% CL upper limit on σₙ(mχ) |

**Configs:**
- [`configs/examples/migdal_heavy_lbc_1p3kgday.json`](../configs/examples/migdal_heavy_lbc_1p3kgday.json) — rates: `data/migdal_rates/Si/heavy_150x300/` (150 masses × 300 cross sections)
- [`configs/examples/migdal_light_lbc_1p3kgday.json`](../configs/examples/migdal_light_lbc_1p3kgday.json) — rates: `data/migdal_rates/Si/light_150x300/` (150 masses × 300 cross sections)

Generated from [`configs/migdal_generate_si_heavy_150x300.json`](../configs/migdal_generate_si_heavy_150x300.json) / [`configs/migdal_generate_si_light_150x300.json`](../configs/migdal_generate_si_light_150x300.json).

---

## 2–3. Prerequisites and pipeline

Identical structure to [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) §2–3 — same `build/data_pattern.root` requirement, same fold → pattern-efficiency → Bp+θBr → profile-likelihood chain, only the rate table (Migdal instead of QEDark/QCDark2) and the resulting cross-section axis (σₙ instead of σₑ) differ.

---

## 4. Step 1 — Rate generation (informational — already done)

```bash
python3 utils/migdal_generate_grid.py configs/migdal_generate_si_heavy_150x300.json
python3 utils/migdal_generate_grid.py configs/migdal_generate_si_light_150x300.json
```

Output: `data/migdal_rates/Si/heavy_150x300/` and `light_150x300/`, 45,000 files each (150 masses × 300 cross sections), already generated and complete. Took ~5 minutes each (measured directly) — cheap despite the large file count, because the rate only needs to be computed once per mass; all 300 cross-section values at fixed mass are obtained by linear rescaling (`R(σ) = σ/σ_ref × R(σ_ref)`, exact since the rate is linear in σₙ), so this is not 45,000× the per-point cost.

**Why 300 sigma points, not 30 — this matters far more than it sounds.** An earlier version of this example used only 80 masses × 30 cross sections, later bumped to 150 masses × 30, and the reproduced LBC curve was visibly "faceted" (not smooth) at either resolution — see git history / the session that built this if you want the full story. I initially assumed this was a mass-grid-density issue; increasing masses alone (80→150) barely helped. The actual cause: `--from-qhist` reconstructs the upper limit by linearly interpolating, in log(σ), between the two σ grid points that bracket the q(mχ,σ) = q_threshold crossing (`apps/ccdarksens_plot_limit.cc`, `build_and_push_q`). With only 30 σ points spanning 20 decades, there are only ~30 possible interpolation segments across the whole curve, and each transition between segments introduces a small kink — *and*, because q(σ) is not actually linear in log(σ) between such widely-spaced points, the coarse grid doesn't just look jagged, it's **numerically biased**: comparing against the real published DAMIC-M heavy-mediator curve, the 30-σ-point LBC reproduction averaged ~20% too strong (mean ratio 0.80 to the reference, std 0.12); the 300-σ-point version averages **8% off with a much tighter spread** (mean ratio 1.08, std 0.04) — close to the DM-electron-level agreement. Bumping σ resolution costs essentially nothing (see above), so there's no reason not to use 300 points for any Migdal grid you build yourself. The DM-electron LBC examples were never affected by this because they'd already used 300 σ points from the start.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/migdal_heavy_lbc_1p3kgday.json
build/ccdarksens_scan_generic configs/examples/migdal_light_lbc_1p3kgday.json
```

Outputs: `outputs/migdal_heavy_lbc_1p3kgday/scan_generic.root`, `outputs/migdal_light_lbc_1p3kgday/scan_generic.root`.

**Runtime:** ~10 seconds each (measured directly — the pattern-space profile likelihood scales with total grid points, 45,000 here vs 2,400 at the old 80×30 resolution).

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this reproduction on the same canvas as the 1 kg-year projection and the published literature curves.** Run the projection companion's scans first ([`Student_Examples_Migdal_Projection.md`](Student_Examples_Migdal_Projection.md) Step 2), then:

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

`--migdal` loads the published XENON1T/PandaX-4T/DarkSide-50/SENSEI Migdal curves automatically, plus DAMIC-M 2025 itself when `--heavy` (see the projection doc §6 for why it's excluded from `--light`). `--from-qhist` forces both curves to be rebuilt directly from each file's own q(mχ,σ) histogram, so they're computed identically and directly comparable. This is the primary, recommended way to look at this example's output. Rendered copies are checked in at [`migdal_heavy_combined.pdf`](migdal_heavy_combined.pdf) and [`migdal_light_combined.pdf`](migdal_light_combined.pdf).

If you only want one mediator's LBC curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --migdal --heavy \
  --out-pdf outplots/migdal_heavy_lbc_1p3kgday.pdf \
  outputs/migdal_heavy_lbc_1p3kgday/scan_generic.root \
  "Migdal heavy LBC" 1.642374415149816 heavy
```

(swap `--light`/light paths for that mediator.)

---

## 7. Note

Like the DM-electron LBC configs, these use `efficiency_mc.efficiency_csv` with the `"../data/..."` path fix described in [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) §7. The pre-existing `migdal_scan_si_heavy.json` / `migdal_scan_si_light.json` under `configs/examples/` still have the old, broken path value.
