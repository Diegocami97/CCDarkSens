<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_Absorption_LBC.md -- Dark-photon absorption, real Si target, real LBC pattern-space background/counts reproduction.
-->

# Student Example — Dark Photon Absorption, LBC Reproduction (Si)

Real LBC pattern-space background and observed counts, dark-photon absorption channel, plain silicon target. Companion: [`Student_Examples_Absorption_Projection.md`](Student_Examples_Absorption_Projection.md) — a same-material **Si** 1 kg-year projection meant to be viewed alongside this LBC reproduction.

---

## 1. What this example does

| Item | Value |
|---|---|
| Signal model | Dark photon absorption, **Si** |
| Observable | Six pattern bins: `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` |
| Exposure | 1.3 kg-day (same LBC dataset as the other pattern-space examples) |
| Background | Real LBC `Bp`/`Br` |
| Data | Real observed counts `[144, 0, 0, 1, 0, 0]` |
| Statistic | Profile likelihood ratio (pydme-style), 90% CL upper limit on ε(mA') |

**Config:** [`configs/examples/absorption_si_lbc_1p3kgday.json`](../configs/examples/absorption_si_lbc_1p3kgday.json), built from `configs/darkphoton_scan_si_lbc.json`.

**Grid:** 200 masses × 40 couplings = 8,000 rate points, in `data/darkphoton_rates/Si_lbc_full/` — regenerated fresh for this example (see §7).

---

## 2–3. Prerequisites and pipeline

Same structure as [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) §2–3 — same `build/data_pattern.root` requirement, same fold → pattern-efficiency → Bp+θBr → profile-likelihood chain, only the rate table (dark-photon absorption instead of scattering) and cross-section axis (ε(mA') instead of σₑ(mχ)) differ. See the projection companion doc §3 for the absorption-specific physics description.

---

## 4. Step 1 — Rate generation (informational — regenerated for this example, see §7)

```bash
python3 utils/darkphoton_generate_grid.py configs/examples/darkphoton_generate_si_lbc_full.json
```

Output: `data/darkphoton_rates/Si_lbc_full/`, 8,000 files, verified as an exact match to the declared 200×40 grid.

---

## 5. Step 2 — Run the scan

```bash
build/ccdarksens_scan_generic configs/examples/absorption_si_lbc_1p3kgday.json
```

Output: `outputs/absorption_si_lbc_1p3kgday/scan_generic.root`.

---

## 6. Step 3 — Plot the limit curve

**The most useful view puts this reproduction on the same canvas as the Si 1 kg-year projection.** Run the projection companion's scan first ([`Student_Examples_Absorption_Projection.md`](Student_Examples_Absorption_Projection.md) Step 2), then:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --dark-photon \
  --out-pdf outplots/absorption_si_combined.pdf \
  outputs/absorption_si_lbc_1p3kgday/scan_generic.root "Si LBC (reproduced)" \
  outputs/absorption_si_projection_1kgyear/scan_generic.root "Si, 1 kg-yr projection" \
  1.642374415149816 heavy
```

`--dark-photon` loads the published XENON1T-bracket and stellar-cooling literature curves from `data/previous_limits/dark_photon/` automatically. `--from-qhist` rebuilds each curve directly from its own q(mA',ε) histogram rather than the scan's own stored UL graph — used consistently across every example in this set. A rendered copy is checked in at [`absorption_si_combined.pdf`](absorption_si_combined.pdf).

If you only want this reproduction's own curve in isolation:

```bash
build/ccdarksens_plot_limit --batch --from-qhist --dark-photon \
  --out-pdf outplots/absorption_si_lbc_1p3kgday.pdf \
  outputs/absorption_si_lbc_1p3kgday/scan_generic.root \
  "Dark photon absorption, Si, LBC" 1.642374415149816 heavy
```

---

## 7. Three real bugs this example surfaced and fixed

**7.1 — Wrong-range rate files (data problem).** The pre-existing `data/darkphoton_rates/Si/` directory (used by the already-committed `configs/darkphoton_scan_si_lbc.json`) turned out to contain rate files with ε values in the wrong range: filenames parsed to ε ≈ 10⁻²⁶–10⁻⁴⁶, while the config's declared grid is ε = 10⁻²⁰–10⁻¹⁰. This is not a partial/incomplete-grid situation — the files simply describe a different coupling range entirely, so they cannot be used as-is for this config.

**Fix:** rather than touch the existing `data/darkphoton_rates/Si/` directory (which the pre-existing, already-committed config still points at), the full 200×40 grid was regenerated fresh into a new directory, `data/darkphoton_rates/Si_lbc_full/`, and this example's config (`absorption_si_lbc_1p3kgday.json`) points `model.rates_dir` there instead. No existing files were modified or deleted. The pre-existing `configs/darkphoton_scan_si_lbc.json` still points at the old, wrong-range `Si/` directory and will need the same fix (either regenerate its rates or repoint it at `Si_lbc_full/`) before it can be trusted.

**7.2 — `model.Emax_eV` silently truncating the signal above 20 eV (a real code-path bug, config-level fix applied).** A direct point-by-point comparison of this reproduction against the real published DAMIC-M curve (`data/previous_limits/dark_photon/DAMIC-M_2025_DAMICmodel_HP.txt`) showed excellent agreement (mean ratio ≈ 0.9–1.0) for mA' below ~20 eV, but for every mass above that the upper limit was pinned to exactly the ε-grid's weak edge (1e-10) — a dead giveaway of *no signal at all* reaching the likelihood, not a genuine physics limit.

Root cause: `model.Emin_eV`/`Emax_eV`/`nbins` in the scan config set the actual histogram binning that `ModelFactory::MakeSignalSpectrumE` fills the rate CSV into (`src/model/ModelFactory.cc`) — this is *not* just descriptive metadata. This config declared `Emax_eV: 20`, but the mass grid runs up to 100 eV, and each absorption rate CSV is a narrow "boxcar" line centered exactly at its own mA'. For any mA' > 20 eV, the boxcar's energy content falls entirely outside the declared histogram range and is silently dropped — `MakeSignalSpectrumE` still returns a valid all-zero-signal TH1D (by design, so the scan doesn't crash), but with zero events, the likelihood no longer prefers any of the tested ε values and the UL search returns the grid boundary as a "no exclusion found here" sentinel.

**Fix applied (config-only, no C++ change):** `model.Emax_eV` raised from 20 to 105 eV (`nbins` scaled to 1050 to keep the same 0.1 eV bin width the rate CSVs themselves use), covering the full mass grid with margin. After the fix, the reproduction tracks the published curve closely across the *entire* range, including its sharp post-minimum rise near mA'≈20–25 eV, which was completely absent before. This is the same class of bug as the DM-electron heavy-LBC path-resolution issue (§7 of that doc) — a pre-existing config value that happened to work for the narrower range it was originally built for, but silently breaks once the grid is extended.

This example's `efficiency_csv` field also has the same `"../data/..."` path fix described in the DM-electron heavy-LBC doc §7.

**7.3 — Shared-code bug: `n_e ≥ 10` spuriously granted full efficiency into every declared pattern (fixed in `src/response/PatternRates.cc`).** After fixing §7.2, a further zoomed-in comparison against the reference curve showed the reproduction tracking it closely up to mA'≈24 eV, then turning over and dropping back down instead of continuing to rise with the reference all the way to its last tabulated point at 30 eV — the opposite of the physically-expected trend (rising mass should mean *falling* acceptance once the deposited charge exceeds what the six declared patterns describe, hence a *rising*, not falling, upper limit).

Root cause, in `FoldNeToPatternRates` (`src/response/PatternRates.cc`, used by every `analysis_space: "pattern"` config): the function summed `h_ne(n_e) × efficiency(pattern, n_e)` over `n_e` for each pattern in the ROI, with a hardcoded shortcut —

```cpp
constexpr int kNeFullEfficiency = 10;
...
double eff = 1.0;
if (ne < kNeFullEfficiency) {
  auto it = pattern_eff_map.find({pattern_id, ne});
  eff = (it != pattern_eff_map.end()) ? it->second : 0.0;
}
```

— that forced `eff = 1.0` for every `n_e ≥ 10`, for *every* pattern in the ROI simultaneously. A real cluster with 10+ ionized electrons has a large, bright shape that isn't any of the six small declared patterns (`{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}` — total charge 2–4 e⁻ each); the original comment's intent ("a bright cluster is certainly detected") was reasonable, but the implementation credited that "certain detection" to *all six* small patterns at once instead of routing it to whatever large-cluster classification a real 10+ e⁻ event would actually get (which isn't one of these six, so it should contribute zero to this ROI). The efficiency table used here (`Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`) only tabulates `n_e = 1..5` in the first place, so this fallback covered a large, entirely untabulated range.

For dark-photon absorption specifically, the charge-ionization table's Fano-broadened P(n_e|E) puts non-negligible tail probability at n_e≥10 once the central n_e is around 6–9 — which for this material's band gap/eh-pair values corresponds to mA' ≈ 24–35 eV, exactly the turnover region found above. This channel is unusually exposed to the bug because its signal is monochromatic and its mass grid deliberately extends past that threshold; the DM-electron and Migdal LBC examples never populate n_e ≥ 10 in practice, so their results are completely unaffected (verified directly: identical upper limits at every checked mass point before and after the fix, e.g. DM-electron heavy at mχ=1000 MeV: 1.35941e-37 both ways; Migdal heavy at mχ=1000 MeV: 1.13658e-37 both ways).

**Fix applied (shared C++ source, approved before editing):** removed the `n_e ≥ 10` special case entirely — a missing `(pattern, n_e)` entry now always counts as 0 acceptance, regardless of how large `n_e` is, consistent with how the function already treated missing entries below 10. A duplicate of the same `kNeFullEfficiency = 10` hack existed in `src/response/ResponseFactory.cc` (feeding the n_e-space fold's combined efficiency histogram) and was removed too; a direct A/B check confirmed it had **zero effect** on the Si n_e-space projection example (bit-identical upper limits before/after at every mass point checked in 35–50 eV), so the projection curve's step structure documented in the companion doc is unaffected and still attributed correctly to real threshold/n_e-bin-transition physics, not this bug. After the fix, the LBC reproduction's red curve continues rising past mA'≈24 eV right alongside the reference black curve instead of dipping back down.
