<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_WIMP_DAMIC_M_Projection.md -- WIMP-nucleon SI, DAMIC-M hypothetical sensitivity projections (1.3 kg-day and 1 kg-year, floor60).
-->

# Student Example — WIMP-Nucleon SI, DAMIC-M Projection (floor60)

Two hypothetical DAMIC-M sensitivity projections for the WIMP-nucleon spin-independent channel, differing only in exposure: a 1.3 kg-day case matched to the real LBC dataset's exposure duration (used elsewhere in this repo), and a 1 kg-year case. Companion: [`Student_Examples_WIMP_DAMIC_SNOLAB.md`](Student_Examples_WIMP_DAMIC_SNOLAB.md) (the real DAMIC 2016 paper reproduction these projections are built from).

> **Two caveats baked into this config's name and content — read both before using this for anything beyond learning the workflow.**
> 1. **Background is a DAMIC-2016 placeholder, not a DAMIC-M measurement.** The flat 5478.75 events/(kg·yr·keV) rate and the Fig. 9-digitized 1×1 background efficiency curve are reused as-is from the DAMIC 2016 SNOLAB reproduction. DAMIC-M does not have its own WIMP-nucleon background model yet.
> 2. **"floor60" — the true-energy floor is DAMIC 2016's own 60 eV, carried over unchanged**, not rescaled for DAMIC-M's ~11× lower readout noise (0.16 e⁻ vs. 1.8 e⁻). This is a controlled, conservative baseline that **understates** DAMIC-M's real low-mass reach. A "floorscaled" variant exists conceptually (scaling the floor down to ~5 eV using the noise ratio) but is not validated and not used here — per project decision, floor60 is the one to use "until we have an actual measurement of the quenching at that level."

---

## 1. What this example does

| Item | 1.3 kg-day (LBC-exposure) | 1 kg-year |
|---|---|---|
| Signal model | WIMP-nucleon, spin-independent, elastic | same |
| Observable | `cluster_energy` — 1×1 channel only | same |
| Exposure | DAMIC-M: 1 kg × 1.3 days (hypothetical — mass_kg=1.0, livetime_days=1.3, matching the real LBC dataset's exposure duration used elsewhere in this repo) | DAMIC-M: 1 kg × 1 year (mass_kg=1.0, livetime_days=365.25) |
| Background | DAMIC 2016's flat rate + Fig. 9 background-efficiency curve (placeholder, see caveat above) | same |
| Data | Asimov | Asimov |
| Statistic | Profile likelihood ratio, 90% CL upper limit on σ_SI(mχ) | same |
| True-energy floor | 60 eV, unscaled (floor60 variant) | same |
| Diffusion / noise / halo | DAMIC-M values sourced from pydme + this repo's existing defaults only (not from the excluded `DamicMSignal-main` reference repo) | same |
| Config | [`configs/examples/wimp_damicm_lbc_1p3kgday.json`](../configs/examples/wimp_damicm_lbc_1p3kgday.json) | [`configs/examples/wimp_damicm_projection_floor60_1kgyear.json`](../configs/examples/wimp_damicm_projection_floor60_1kgyear.json) |

Both configs are otherwise identical, built from `configs/wimp_nucleon_damicm_1x1_floor60_1p3kgday.json`/`_1kgyear.json` — only `livetime_days` differs. The σ grid was extended down to 10⁻⁴⁵ cm² (50 points) in both — the original 10⁻⁴² lower bound left the true crossing unresolved near the grid edge for a projection this sensitive. The 1.3 kg-day case exists so this channel has the same "near-term LBC-scale exposure vs. eventual 1 kg-year" pairing every other channel in this example set has — it is still a hypothetical Asimov projection, not a real-data reproduction (there is no real WIMP-nucleon LBC dataset/analysis to reproduce).

---

## 2. Prerequisites

Same as the SNOLAB reproduction — build `ccdarksens_scan_generic` and `ccdarksens_plot_limit` (the WIMP-nucleon cluster-fit machinery is included in both).

---

## 3. Pipeline overview

Identical mechanics to the SNOLAB reproduction (see [`Student_Examples_WIMP_DAMIC_SNOLAB.md`](Student_Examples_WIMP_DAMIC_SNOLAB.md) §3): noise-tail calibration → kernel build → fold → profile likelihood, all internal to the scan run. Only the detector exposure (1 kg-year vs. 0.01 kg × 60 days), readout noise (0.16 e⁻ vs. 1.8 e⁻), and channel count (1×1 only, no 1×100) differ — the background model and true-energy floor are deliberately carried over unchanged from the SNOLAB case as the current best-available placeholder.

---

## 4. Step 1 — Rate/kernel generation (informational — happens inside the scan run)

As with the SNOLAB example, there is no separate offline step; `ccdarksens_scan_generic` runs calibration and kernel-building internally. To inspect the DAMIC-M-specific halo/diffusion defaults this config draws on, see `configs/wimp_nucleon_generate_si_damicm_halo.json`'s own `_comment` for their exact provenance.

---

## 5. Step 2 — Run the scans

```bash
build/ccdarksens_scan_generic configs/examples/wimp_damicm_lbc_1p3kgday.json
build/ccdarksens_scan_generic configs/examples/wimp_damicm_projection_floor60_1kgyear.json
```

Outputs: `outputs/wimp_damicm_lbc_1p3kgday/scan_generic.root`, `outputs/wimp_damicm_projection_floor60_1kgyear/scan_generic.root`.

**Runtime:** ~50 seconds each (measured directly) — heavier than an n_e/pattern-space scan (noise-tail + kernel MC trials per point), but not extreme at this grid size.

---

## 6. Step 3 — Plot the limit curve

**The recommended view puts both exposures on one canvas against the trusted SNOLAB reference and its band:**

```bash
build/ccdarksens_plot_limit --batch --wimp-minimal \
  --out-pdf outplots/wimp_damicm_combined.pdf \
  outputs/wimp_damicm_lbc_1p3kgday/scan_generic.root "DAMIC-M, 1.3 kg-day (LBC-exposure)" \
  outputs/wimp_damicm_projection_floor60_1kgyear/scan_generic.root "DAMIC-M, 1 kg-yr projection" \
  2.71 heavy
```

**`q_threshold = 2.71`, not `1.642374415149816`** — same reasoning as the SNOLAB companion doc §6: `2.71` is `ccdarksens_plot_limit`'s own default (the standard Δχ²=2.71 one-sided convention), correct for this `profile_minimizer: brent` channel; `1.642374415149816` is the asymptotic Cowan et al. value validated for the pydme-style channels only.

`--wimp-minimal` loads only the DAMIC 2016 paper's trusted observed curve and its expected ±1σ band (dropping `--wimp`'s extra secondary curves — median, two uncertain-provenance digitizations, and a literature DAMIC-M-projected-1kg-year curve — kept out here for a cleaner comparison; use plain `--wimp` if you want those back). A rendered copy is checked in at [`wimp_damicm_combined.pdf`](wimp_damicm_combined.pdf). Both DAMIC-M curves sit far below the DAMIC 2016 band, as expected for a much lower-noise, hypothetical future detector.

---

## 7. Note

This example does not have a "floorscaled" sibling in this set — floor60 was the version the group chose to standardize on. If a real quenching/threshold measurement becomes available at low true energy, that would be the trigger to build and validate a floorscaled variant.
