<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_Row_Cut_Recovered_Events.md -- From an event recovered by a modified row cut to its effect on the LBC limit: check the pattern by hand, account for every count that changed, follow what a mask change does to the background and exposure, rerun the scan, and read the result; with every assumption and every input still needed from the student stated explicitly.
-->

# Student Example — An Event Recovered by a Modified Row Cut: What Does It Do to the Limit?

This example answers the last slide of Salma Douieb's "Multiplicity cut on LBC science run data" ("How does this affect our limits?"). It starts from one event and ends at a limit curve. Every step can be run, and every step says what it teaches, what it assumes, and what it still needs from Salma.

Companions: [`Student_Examples_Background_Pattern_Counts.md`](Student_Examples_Background_Pattern_Counts.md) (the background formula), [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md) (the scan this one modifies), and [`Pipeline_Raw_Data_to_Limit.md`](Pipeline_Raw_Data_to_Limit.md) (the whole chain from raw images).

---

## 0. How to read this document

Four labelled boxes appear in every step. They are the point of the document.

> **WHAT YOU LEARN.** The idea the step teaches, in one or two sentences.

> **ASSUMPTION (A1, A2, …).** Something I took as true without being able to check it. If an assumption is wrong, the numbers below it change. All of them are collected in Appendix C.

> **INPUT NEEDED FROM SALMA (S1, S2, …).** Something only her analysis of the data can supply. The example runs without it, but the result is then an estimate, and the step says so. All of them are collected in §11.

> **NEEDED FROM THE COLLABORATION (C1, C2, …).** Something neither Salma nor this repository has: it has to come from the owners of the mask code, the simulation, or the published analysis.

### Vocabulary

| Word | Meaning |
|---|---|
| e⁻ | One electron of charge in a pixel. A pixel reads a number like 1.85 e⁻ because the readout has noise. |
| σ (read noise) | The spread of the readout noise of one CCD, about 0.17 e⁻ for CCD D. A pixel with true charge 1 e⁻ is read as 1.00 ± 0.17. |
| 1 e⁻ bound | The charge above which a pixel is called "at least 1 e⁻": 3.75σ, which is 0.64 e⁻ for CCD D. The 2 e⁻ bound is 1.64 e⁻. |
| Pattern | Neighbouring pixels in a row that together look like one small deposit: `{11}` is two 1 e⁻ pixels side by side, `{21}` a 2 and a 1, `{31}` a 3 and a 1, and so on. The search bins are `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}`. |
| Mask | A rule that removes pixels, rows or columns that are known to be bad before patterns are searched. |
| Row cut | The mask that removes a row when it has more than one pixel at or above 2 e⁻, unless those pixels are next to each other. |
| N_pix | The number of pixels that survive all masks. The published value is 1,853,807,441. |
| D_p | The observed count in pattern bin p. Published: (144, 0, 0, 1, 0, 0) for `{11}, {21}, {111}, {31}, {22}, {211}`. |
| B^rc_p (Bp) | The background in bin p from *random coincidences* of dark-current electrons. It is computed, not measured. |
| Br | The radiogenic background in bin p, per unit of its nuisance parameter θ. It comes from a Geant4 simulation. |
| θ | The nuisance parameter that scales Br. The fit may move it from 1. |
| Limit | The largest cross section σ_e that the data still allow at 90% confidence. A weaker limit is a larger number. |
| q | The test statistic of the profile likelihood. A cross section is excluded when q is above 1.642. |

---

## 1. What this example does

| Item | Value |
|---|---|
| Starting point | One (3.17, 1.85) e⁻ pair in row 1574 of CCD D (module 104), columns 6083–6084, removed by the second multiplicity-cut clause. A modified row cut gives counts (158, 0, 0, 2, 0, 0) instead of (144, 0, 0, 1, 0, 0). |
| Question | How do the counts, the background and the limit change? |
| Tools | `build/ccdarksens_classify_pair`, `build/ccdarksens_compute_brc`, `build/ccdarksens_scan_generic`, `build/ccdarksens_plot_limit` (each one is taught in [`Student_Examples_Pattern_Pipeline_Walkthrough.md`](Student_Examples_Pattern_Pipeline_Walkthrough.md)) |
| Scans | Five configs `configs/examples/dm_electron_heavy_lbc_rowcut_{A,B,C,D,E}.json`, all copies of the published LBC config (Step 7) |
| Result | With the published background held fixed the limit weakens by 1.4–1.7× at heavy mass and about 1.9× below 3 MeV. If pixels were returned, the background and exposure rise and the weakening is smaller (1.4× below 3 MeV and 1.4× at heavy mass in the most generous case, f = 6.2%). |

Three questions, kept separate on purpose:

1. **Is the event a pattern?** (Steps 2 and 3. A hand calculation.)
2. **Why did the counts change?** (Step 4. Every change must be attributed to an event.)
3. **What does the change do to the limit?** (Steps 5 to 8. Counts, background and exposure move together.)

---

## 2. Prerequisites

```bash
cmake --build build -j8 --target \
  ccdarksens_classify_pair ccdarksens_compute_brc ccdarksens_scan_generic ccdarksens_plot_limit
```

The heavy-mediator rate tables must be in `data/qedark_rates/Si/heavy/long_scan/` (240,000 files, not in the repository). They are generated by `python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_heavy_dense.json`, which Stage R of the walkthrough explains; at the measured speed this is of the order of ten minutes. Each scan then takes about two minutes on a laptop; run them in parallel.

---

## 3. Step 1 — Counts, background and exposure come from the same pixels

The published analysis takes calibrated images and removes pixels with masks. The pattern search runs on what is left, and each bin's count is compared with a background that was computed **from the same surviving pixels**:

```
 raw images ──► masks ──► surviving pixels N_pix ──► pattern search ──► counts D_p
                              │
                              └──► rates of isolated 1, 2, 3 e- pixels
                                          │
                                          └──► B^rc_p = N_pix · Σ_q Π Poisson · B[p|q]
```

A change to a mask moves three things at once: the counts `D_p`; `N_pix` (and so `B^rc_p`, which is proportional to it); and the exposure that multiplies the signal. Replacing `observed_counts` in a config and rerunning changes only the first of the three.

> **WHAT YOU LEARN.** A limit is a comparison between *counts*, *expected background* and *expected signal*. All three depend on which pixels were kept, so none of them can be edited alone.

---

## 4. Step 2 — Is the event a pattern? Do it by hand

The data-selection note defines, for two consecutive pixels with charges `q_j`, `q_j+1`,

```
p_mn = -log[ cdf(q_j ; m) · cdf(q_j+1 ; n) ] ,     cdf(q ; m) = Φ((q − m)/σ)
```

where `Φ` is the standard normal CDF. In words: `cdf(q; m)` is the chance that a pixel that really holds `m` electrons reads *at most* `q`. `p_mn` is small when both pixels read at least as much as the hypothesis `(m, n)` needs.

The pair is a 1+1 candidate if `p_11 < 4`, promoted to 2+1 if `min(p_21, p_12) < 4`, and to 3+1 if `min(p_31, p_13) < 4`. The total charge must stay below 5.6 e⁻.

Salma's plot gives the 1 e⁻ bound of CCD D as 0.64 e⁻, and the bound is 3.75σ, so σ = 0.64/3.75 = 0.171 e⁻.

```bash
build/ccdarksens_classify_pair 3.17 1.85 --bound 0.64
```

```
sigma = 0.1707 e-   Qmin = 0.6400 e-   Qmax = 5.6000 e-   pixels = 3.17 1.85   total = 5.0200 e-

p_mn = -log[cdf(q1;m) cdf(q2;n)] for the first pair (3.17, 1.85)
   (1,1)      0.00
   (2,1)      0.00
   (1,2)      1.66
   (3,1)      0.17
   (1,3)     25.55
   (2,2)      1.66
   (3,2)      1.84
   (2,3)     25.55
chain: p_11 = 0.00  p_21 = min(p_21,p_12) = 0.00  p_31 = min(p_31,p_13) = 0.17  p_22 = 1.66   (pairs pass below 4.00, triples below 5.50)

selection chain : {3,1} (start pixel 0, total 5.02 e-)
selection tree  : {3,2} (start pixel 0, total 5.02 e-)
```

How to read it:

- **The CDF is a lower tail.** A pixel far *above* the hypothesis costs nothing (`cdf = 1.0000`). A pixel far *below* it is expensive (`cdf = 0.0000`, `p = 25.55` in the (1,3) row). The test only asks "is this pixel too small to hold `m` electrons?".
- **So the chain climbs until it fails.** 1+1 passes, 2+1 passes, 3+1 passes (`p_31 = 0.17`). A 4+1 hypothesis fails because 3.17 e⁻ is too small to hold four electrons (`p_41 = 14`, far above 4), so the chain stops at 3+1.
- **There is no 3+2 bin in the search.** The pair takes the label of the highest class that passes, 3+1, and falls in `{31}`. That is also why a 5 e⁻ deposit is part of the `{31}` signal efficiency at n_e = 5.
- **The label depends on the selection rule.** The framework's `chain` selection (the published one) returns `(3,1)` for this pair. The notebook `tree` selection of `efficiencies.py` has its own (3,2) class and returns `(3,2)`, which is not a search bin. Same pixels, different bin. Always state which selection a table was made with.

Try other pairs to see the other branches. A clean 1+1+1 (the third pixel must pass `p_111 < 5.5`):

```bash
build/ccdarksens_classify_pair 1.05 0.80 --sigma 0.16 --left 0.95
```

> **WHAT YOU LEARN.** How a pattern label comes out of a few lines of arithmetic, why a label is a property of the *selection rule* and not of the pixels, and how to check any candidate yourself.

> **ASSUMPTION (A7).** σ = 0.171 e⁻ for CCD D comes from reading "1 e⁻ bound 0.64" off her plot and dividing by 3.75. The charges 3.17 and 1.85 are taken as calibrated electrons.

> **INPUT NEEDED FROM SALMA (S8).** Confirm: image ID, module and CCD names, the calibrated charges of the two pixels, and the σ (or 1 e⁻ bound) of that CCD in that image.

---

## 5. Step 3 — Read the rules exactly

The published pattern definition (PRL 135, 071002 and the note) adds three isolation clauses:

1. the pixels in the rows above and below the pattern are < 1 e⁻;
2. the pattern pixels in the other CCDs are < 1 e⁻;
3. no other pattern or pixel ≥ 2 e⁻ is in the columns of the candidate pattern.

Separately, the mask step removes rows and columns "with N_{q≥2e⁻} > 1, **unless these pixels are consecutive along a row**". This is the clause that fired.

Check the event against each, using her slides:

| Check | Evidence on the slides | Result |
|---|---|---|
| Rows above and below < 1 e⁻ | CCD D panel: all neighbours below 0.5 e⁻ | pass |
| Same pixels in CCDs A and B < 1 e⁻ | 0.09, 0.01 in A; 0.03, 0.14 in B | pass |
| No other ≥ 2 e⁻ pixel in columns 6083–6084 | occupancy histograms; no other hits ≥ 2 e⁻ | pass |
| Row cut: pixels ≥ 2 e⁻ in row 1574 | 3.17 and 1.85 e⁻ only (the 2 e⁻ bound is 1.64); they are adjacent | **should be exempt** |

On her slides the row carries mask code −701 anyway. I cannot read that code from the slides. Three explanations are possible and the example does not choose between them:

- the exemption is implemented differently from the sentence quoted above (for example it needs more than two pixels, or a different charge threshold);
- the exemption exists, but a different clause with the same code removed the row;
- the exemption is not applied in the code that produced the published masks.

> **WHAT YOU LEARN.** How to turn a published sentence into a checklist and test one event against it. And that a rule written in a paper is not the same as the rule in the code.

> **ASSUMPTION (A8).** The event is a genuine ionisation event. Her three checks (not a hot column, no crosstalk, not a CTI trail) argue for it, and I do not re-verify them.

> **INPUT NEEDED FROM SALMA (S7).** (a) Where in the mask code the −701 flag is set, and the line that tests "consecutive". (b) Whether any *other* row in the dataset has two adjacent pixels ≥ 2 e⁻, and what happened to it. She found one event; the question is whether the exemption works at all.

> **NEEDED FROM THE COLLABORATION (C1).** The owner of the mask code says which of the three explanations holds.

**Why it matters beyond one event.** The exemption exists to protect patterns with two adjacent pixels ≥ 2 e⁻, which are exactly `{22}` and the (3,2) members of `{31}`. If it does not work in general, the signal-efficiency tables (which were simulated without this mask) overstate the efficiency of those shapes.

> **NEEDED FROM THE COLLABORATION (C3).** If the exemption was not applied, whether the `{22}` and (3,2) signal efficiencies of the published analysis need a correction.

---

## 6. Step 4 — Account for every count that changed

Slide 4 compares the published counts with a "row cut modified" column:

| Bin | Published | Modified |
|---|---|---|
| {11} | 144 | 158 |
| {31} | 1 | 2 |
| {21}, {22}, {111}, {211} | 0 | 0 |

The tracked per-image table gives the published side directly:

```bash
awk -F, 'NR>1{a+=$6; b+=$7; n+=$5} END{printf "Count_Candidate_11 %d\nCount_Candidate_31 %d\nNusedpix %d\n", a, b, n}' \
  data/Final_Combined_Image_Data.csv
awk -F, 'NR>1 && $7>0{print $1}' data/Final_Combined_Image_Data.csv       # images with a {31} candidate
```

```
Count_Candidate_11 144
Count_Candidate_31 1
Nusedpix 1853807441
4493.0
```

Four things to take from this:

- **The one recovered event moves only `{31}` by +1.** It cannot account for +14 in `{11}`.
- **It is a second `{31}` event, not the published one.** The published candidate is the only image with `Count_Candidate_31 = 1` (image 4493) and lies in column 2103 of CCD 1-D according to the PRL; hers is in CCD D of module 104, column 6083.
- **The table is missing a column.** Her replication with the *original* cut is not shown. Without it, `158 − 144` mixes two effects: the modified cut, and whatever differs between her replication and the published selection. The two must be separated before the 158 means anything.
- **Attribute event by event.** For every bin and every image where the count differs, list the event (CCD, row, columns, charges) and the reason: "recovered by the row cut", "lost by the row cut", "classified differently". A change with no event attached is a difference in the replication, not in the data.

The PRL also gives a control that the modified selection can be run through directly: the same pattern selection applied to *vertical* and *diagonal* pixel pairs, where no ionisation signal is expected, found **139** and **130** `{11}` patterns against 141.4 expected from random coincidences. If the modified cut produces 158 `{11}` along the row but still about 140 vertically and diagonally, the extra counts are specific to the row direction and need an explanation. If the vertical and diagonal counts rise too, the change is a change in the pixels kept and `N_pix` must follow (Step 5).

> **WHAT YOU LEARN.** A count is only meaningful with a list of the events behind it. A number that changed without a named event is a warning about the *procedure*, not a result about the data.

> **ASSUMPTION (A2).** The modified cut gives (158, 0, 0, 2, 0, 0) exactly as on slide 4. I do not know how it is defined.

> **INPUT NEEDED FROM SALMA (S1).** The exact definition of "row cut modified": the changed rule in words or the code diff.
>
> **INPUT NEEDED FROM SALMA (S2).** The baseline column: pattern counts of her replication with the *original* row cut, per bin and per image.
>
> **INPUT NEEDED FROM SALMA (S3).** An event list for every candidate that differs between the published and the modified counts: image, CCD, row, columns, charges, and which cut removed or kept it.
>
> **INPUT NEEDED FROM SALMA (S5).** The vertical and diagonal `{11}` control counts under the modified cut (the published values are 139 and 130).

---

## 7. Step 5 — What a mask change does to the background

`B^rc_p` is proportional to the number of pixels that survived the masks. If the modified cut returns a fraction `f` of the pixels, then to first order, at fixed rates,

```
B^rc_p  →  B^rc_p · (1 + f)
```

`ccdarksens_compute_brc` reproduces this proportionality, here for `{11}` with uniform rates (θ = 2.9, 2, 2):

```bash
for f in 0 0.001 0.01 0.0622; do
  n=$(awk -v f=$f 'BEGIN{printf "%.0f", 1853807441*(1+f)}')
  echo "f=$f"
  build/ccdarksens_compute_brc --config configs/examples/background_pattern_rates.json \
    --n-pixels $n | grep -E "^11 " | head -1
done
```

(The first number on each line is `B^rc_{11}` in scientific notation: `1.4684e+02` is 146.84.)

| f (fraction of pixels returned) | {11} from the app (uniform rates) | {11} from published Bp = 141.4 × (1+f) |
|---|---|---|
| 0 | 146.8 | 141.4 |
| 0.1% | 147.0 | 141.5 |
| 1% | 148.3 | 142.8 |
| 6.2% | 156.0 | 150.2 |

The two columns differ by a constant 4%: the script includes the lone-pixel misidentification term and uses one rate for all images, while the published 141.4 is the diagonal term. The scaling is what matters. **When recomputing Bp, always compute the old and the new value with the same procedure and take the ratio; never put the script's 146.8 next to the published 141.4.**

**An upper bound on f.** The masks keep 94% of the pixels (Table II of the note gives 0.9414 and 0.9474 for the two datasets). Even if *every* masked pixel came back, `f = 1/0.9414 − 1 = 6.2%`, and the expectation in `{11}` would rise only to 150.2. Her 158 is still 7.8 above that. A single row-cut clause returns far fewer pixels, so an expectation near 142–143 is more realistic and the excess is closer to +15 than +8.

Two consequences:

- the background must be recomputed with the new `N_pix`, not left at the published numbers (Step 6);
- the exposure (and the signal efficiency of the mask) rises by the same fraction `f`, which makes the limit slightly *stronger*. This goes the other way from the count increase.

> **WHAT YOU LEARN.** Background is not a number written in a config; it is the product of a pixel count and a rate. Returning pixels returns background with them, so extra events from returned pixels are partly expected, not a signal.

> **ASSUMPTION (A4).** The returned pixels have the same dark-current rates as the surviving ones, so `B^rc` scales with `1 + f`. Rows removed by a cut may be more active than average, in which case this underestimates the background.
>
> **ASSUMPTION (A5).** The signal exposure scales with `1 + f` as well, because the surviving fraction of every image rises by `f`.

> **INPUT NEEDED FROM SALMA (S4).** The pixels returned by the modified cut: total, and per image and CCD. This is `f`. Until she supplies it, the example uses f = 1% and f = 6.2% as the plausible range.

---

## 8. Step 6 — Recompute Bp

What the repository can do, and what has to come from the data:

| Quantity | Repository | From her pass over the data |
|---|---|---|
| `B[p\|q]` matrix | `data/Background_efficiencies.csv`; regenerate with `ccdarksens_pattern_background_eff` | nothing new |
| Signal efficiencies | `data/Efficiencies_patterns_*.csv`, `data/efficiencies_paolo.csv` | nothing, unless the exemption question (C3) changes the `{22}` / (3,2) efficiency |
| `N_pix` after the modified cut | the published total 1,853,807,441 | the returned pixels per image (S4) |
| Isolated 1, 2, 3 e⁻ rates | per-image 1 e⁻ counts only | per image **and per 10-column group** for 1, 2 and 3 e⁻, after the modified cut (S6) |
| `Br` (radiogenic) | published vector, from Geant4 | not changed by a mask unless the radiogenic control sample was affected (C2) |

With only the pixel total, scale the published `Bp` by `1 + f` (one line):

```bash
awk -v f=0.01 'BEGIN{split("141.4 0.111 0.042 0.019 2.5e-5 5.8e-5", b, " ");
                     for (i = 1; i <= 6; i++) printf "%.6g ", b[i]*(1+f); print ""}'
```

With per-segment rates, use `ccdarksens_compute_brc --segments table.csv` (columns `n_pixels, lambda_1, lambda_2, lambda_3plus`, one row per image, CCD and column group; Stage 3 of the walkthrough). That reproduces the structure of Eq. (4) of the paper.

> **WHAT YOU LEARN.** Which of the ingredients of Bp are properties of the detector (the matrix, fixed), which are properties of the data (rates and pixel counts, her job), and which are not hers at all (Br).

> **ASSUMPTION (A1).** `B[p|q]` does not depend on the row cut: it is a noise-only simulation of 3×5 pixels and does not include masks.

> **INPUT NEEDED FROM SALMA (S6).** Isolated 1, 2 and 3 e⁻ pixel counts per image and per 10-column group after the modified cut, **and the same numbers for the original cut**, so the effect of the cut can be isolated. If she cannot produce them, Bp stays a scaling estimate.

> **NEEDED FROM THE COLLABORATION (C2).** Whether the radiogenic control sample (98 clusters at 2.5–7.5 keV) passes through the row cut, and the Geant4 τ_p needed to give any new bin a Br.

State clearly in anything written up which of the two methods was used. The pixel-total scaling is an estimate; the segment sum is the calculation.

---

## 9. Step 7 — Run the scans

Five configs keep everything of the published LBC example except the observed counts and, in D and E, the background and exposure:

| Config | `observed_counts` | Bp and exposure | What it isolates |
|---|---|---|---|
| (published) `dm_electron_heavy_lbc_1p3kgday.json` | 144, 0, 0, 1, 0, 0 | published | reference |
| `…_rowcut_A.json` | 144, 0, 0, **2**, 0, 0 | published | the recovered event alone |
| `…_rowcut_C.json` | **158**, 0, 0, 1, 0, 0 | published | the `{11}` change alone |
| `…_rowcut_B.json` | **158**, 0, 0, **2**, 0, 0 | published | both count changes, background not updated |
| `…_rowcut_D.json` | **158**, 0, 0, **2**, 0, 0 | × 1.01, exposure × 1.01 | B with f = 1% |
| `…_rowcut_E.json` | **158**, 0, 0, **2**, 0, 0 | × 1.0622, exposure × 1.0622 | B with f = 6.2% |

A, C and B are a decomposition: A is the `{31}` effect, C the `{11}` effect, B both. D and E show how much of B is undone once the background is allowed to follow the pixels.

```bash
for k in A B C D E; do
  build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_rowcut_$k.json &
done
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday.json
wait
```

### 9.1 Do not read the stored limit when the data show an excess

Open the log of run A:

```
[scan-generic] mchi=50.0162 -> upper_limit=1e-46
```

`1e-46 cm²` is the first point of the cross-section grid. It is not a limit. The stored limit is the first crossing of `q(σ)` through the threshold (1.642). With two events in `{31}` against 0.07 expected, `q_0` at `σ → 0` is **already above the threshold** for many masses, there is no crossing, and the code returns the lowest grid point. Read it as "the data prefer a signal at this mass".

| Run | Mass points (of 800) with the limit on the grid floor |
|---|---|
| published | 0 |
| A ({31} = 2) | 271 |
| C ({11} = 158) | 103 |
| B | 466 |
| D (f = 1%) | 356 |
| E (f = 6.2%) | 294 |

The scan also stores a continuous limit that starts at the best-fit σ and steps upward (`upper_limit_sigma_e_mchi_pydme_bisection`). That one is a proper one-sided limit when there is an excess. For the published counts it agrees with the grid crossing to 3% at 0.7–3 MeV and 14–21% above (median over each window), with spikes where a mass sits between n_e thresholds. Use it for the comparison and quote **window medians**, not single mass points.

### 9.2 Compare

```bash
mkdir -p outputs/row_cut_recovered_events
build/ccdarksens_plot_limit --batch --heavy --draw-both --ratio-to published \
  --out-pdf outputs/row_cut_recovered_events/limit_compare.pdf \
  --out-root outputs/row_cut_recovered_events/limit_compare.root \
  --out-csv outputs/row_cut_recovered_events/limit_compare.csv \
  outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root published \
  outputs/dm_electron_heavy_lbc_rowcut_A/scan_generic.root A \
  outputs/dm_electron_heavy_lbc_rowcut_B/scan_generic.root B \
  outputs/dm_electron_heavy_lbc_rowcut_C/scan_generic.root C \
  outputs/dm_electron_heavy_lbc_rowcut_D/scan_generic.root D \
  outputs/dm_electron_heavy_lbc_rowcut_E/scan_generic.root E \
  1.642374415149816 heavy
```

`--draw-both` makes the app compare the continuous bisection limit of each input. To count the points of the *stored* limit that sit on the sigma floor (the table of §9.1), run the same command without `--draw-both`: the ratio table then reports "271 of 800 points on the sigma floor" for A, and so on.

Measured with this checkout. Each entry is the ratio of the limit to the published one, median and 16–84% range over the mass window; bigger means a weaker limit.

| m_χ window (MeV) | A: {31} only | C: {11} only | B: both, Bp fixed | D: B with f = 1% | E: B with f = 6.2% |
|---|---|---|---|---|---|
| 0.7 – 3 | 1.00 | 1.90 | 1.87 | 1.74 | 1.38 |
| 3 – 10 | 1.58 | 1.41 | 2.44 | 2.25 | 1.79 |
| 10 – 30 | 1.71 | 1.13 | 1.87 | 1.93 | 1.66 |
| 30 – 100 | 1.42 | 1.05 | 1.66 | 1.66 | 1.42 |
| 100 – 1000 | 1.42 | 1.05 | 1.66 | 1.66 | 1.42 |

(The app prints the 16–84% range of each window as well. The limit has a staircase at the n_e thresholds, so differences between two columns smaller than about 10% should be treated as noise.)

### 9.3 Hand checks

Neither number below uses the framework, and both land close to the table.

- **`{31}`, classical Poisson 90% upper limit on the signal count** with 0.07 expected background: 3.8 events for 1 observed, 5.3 for 2 observed, ratio 1.37. The table has 1.42 in A at 30–1000 MeV, where `{31}` dominates the signal.
- **`{11}`, Gaussian approximation** `s_up ≈ max(0, n − b) + 1.64 √b`:

  | Case | expected b | s_up | ratio to published |
  |---|---|---|---|
  | published (144 observed) | 141.4 | 22.1 | 1 |
  | 158 observed, Bp fixed | 141.4 | 36.1 | 1.63 |
  | 158 observed, f = 1% | 142.8 | 34.8 | 1.57 |
  | 158 observed, f = 6.2% | 150.2 | 27.9 | 1.26 |

  The framework gives 1.87, 1.74 and 1.38 below 3 MeV in B, D and E. It profiles the radiogenic nuisance, so it is not expected to match to better than that, but it moves in the same direction and by a similar amount.

> **WHAT YOU LEARN.** How to build a controlled comparison: change one thing per run, keep a run that changes only one input, and check each result with an independent hand estimate. And that the number a program prints is not always a limit (§9.1).

> **ASSUMPTION (A6).** The 90% threshold used is the asymptotic one-sided value 1.642, as in the published example. For small counts in `{31}` the toy-Monte-Carlo threshold is higher (2.7–3.0 in the note for low masses), which would make a limit *weaker*, not stronger.
>
> **ASSUMPTION (A3).** The pattern selection is the published `chain`. Signal efficiencies are from `Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`, as in the published example, and they are *not* recomputed for the modified cut.

> **INPUT NEEDED FROM SALMA (S1, S4).** Scans A–E are scenarios because her cut definition (S1) and returned pixel fraction (S4) are not known. With them, only the one matching scenario needs to be run, with the real Bp.

---

## 10. Step 8 — Read the result

Four questions, in order.

**Which bin drives which mass?** Below 3 MeV the signal is one or two electrons, so only `{11}` matters: A has no effect and C carries all of B. Above 30 MeV deposits reach 4–5 e⁻ and `{31}` dominates: A carries B and C changes the limit by 5%. In between both matter and the ratios multiply roughly (1.58 × 1.41 ≈ 2.2 against the measured 2.44 for B).

**How much does letting the background follow the pixels matter?** Mostly where `{11}` drives the limit. Going from B to E, the ratio below 3 MeV falls from 1.87 to 1.38. At 30–1000 MeV it falls from 1.66 to 1.42, through the 6% larger exposure and the higher background in the other bins; I have not separated the two. The honest result, until S4 is supplied, is the range B to E.

**How surprising are the counts?** With `Bp + θ·Br ≈ 0.019 + 0.052 = 0.07` expected in `{31}`:

```
P(n ≥ 1) = 1 − e^(−0.07)               = 6.9%     (published count, one event)
P(n ≥ 2) = 1 − e^(−0.07)(1 + 0.07)     = 0.24%    (two events, about 2.8σ one-sided)
```

Multiplying by the six bins gives about 1.4%, so this is not a discovery claim. It is a mild excess in one bin.

**How much of that to believe?** Two reasons not to treat 0.07 as exact:

- the reproduction audit found that this framework's `B^rc` for `{31}` is **28% below** the published value when built from the available inputs, and `{111}` 10% above (see `Pipeline_Raw_Data_to_Limit.md`). The per-band 3 e⁻ populations are the missing input (S6);
- `Br` for `{31}` (0.052) comes from a Geant4 simulation that is not in the repository (C2).

> **WHAT YOU LEARN.** How to state a result honestly when inputs are missing: as a range with the missing input named, and with the significance of the counts shown next to it.

> **ASSUMPTION (A9).** `Bp` and `Br` for `{31}` are the published ones. The "2.8σ" is a counting statement with those numbers; it is not a statement about the physical origin.

---

## 11. Step 9 — Everything that is still needed

### From Salma

| # | What | Why | Used in |
|---|---|---|---|
| S1 | Exact definition of "row cut modified" (rule in words or code diff) | tells which scenario is real and whether a replication differs from the published selection | Steps 4, 7 |
| S2 | Baseline counts of her replication with the **original** cut, per bin and per image | separates the effect of the cut from the effect of the replication | Step 4 |
| S3 | Event list for every candidate that differs (image, CCD, row, columns, charges, cut that removed or kept it) | attributes each count change to an event | Step 4 |
| S4 | Pixels returned by the modified cut, total and per image and CCD (= f) | sets how much Bp and exposure rise | Steps 5, 6, 7 |
| S5 | Vertical and diagonal `{11}` control counts with the modified cut (published 139 and 130) | tells whether extra counts are specific to the row direction | Step 4 |
| S6 | Isolated 1, 2, 3 e⁻ counts per image and per 10-column group, original and modified cut | allows Bp to be recomputed with Eq. (4) rather than scaled | Step 6 |
| S7 | Location of the −701 flag and the "consecutive" test in the mask code; whether other rows with two adjacent ≥ 2 e⁻ pixels exist and what happened to them | decides which of the three explanations holds | Step 3 |
| S8 | Image ID, module/CCD names, calibrated charges and σ for the event | fixes the inputs of Step 2 | Step 2 |

### From the collaboration

| # | What | Who | Used in |
|---|---|---|---|
| C1 | Which of the three explanations holds for the row-cut exemption | owner of the mask code | Step 3 |
| C2 | Geant4 τ_p for any new bin, and whether the radiogenic control sample passes through the row cut | owner of the radiogenic background | Steps 6, 8 |
| C3 | Whether the published `{22}` and (3,2) signal efficiencies need a correction if the exemption was not applied | owners of the signal-efficiency table | Step 3 |

With S1–S4 the counts and `f` are settled and the limit comparison can be redone with a consistent background. With S5 the row-direction question is settled. With S7 and C1 the efficiency question is settled. With S6 `Bp` can be recomputed with the structure of Eq. (4) instead of a pixel-total scaling.

---

## 12. Limits of this example

- A, B and C use the published `Bp` and `Br`, so they bound the effect of the count change. D and E add the background and exposure scaling but, until S4 is known, they are scenarios, not a result for the modified analysis.
- The bisection limit is the quantity compared; its median differs from the grid-crossing limit by up to 21% for the published counts. Ratios of window medians are robust to that, single points are not.
- Sensitivity bands (Asimov, toys) are not recomputed: the modified counts are data, so the observed limit is what moves.
- Whether the event is signal, background or a detector artefact is not decided here. `{31}` has a non-zero background from random coincidence and radiogenics that this example cannot refine (S6, C2).

---

## Appendix A — Sketch of the event-selection pass

The framework does not yet contain the pass from calibrated images to pattern candidates. This is the outline of what Salma's replication does and where each rule of this example sits. Each stage is a place to count pixels, which is what Steps 4 and 6 need.

```
for each image, for each CCD:
    q = calibrated charge in electrons, rows x columns (100x1 binned)
    keep = all pixels

    # stage 1: masks (note, section "Masks")
    hot columns        : column median/MAD > threshold        -> mask column +/- 5
    clusters/crosstalk : pixels near > 100 e- or in other CCD -> mask
    trailing           : pixels after > 100 e-                -> mask
    noise/variance     : row/col noise or correlation failing -> mask
    row/col multiplicity: count N(q >= 2e-) per row and col;
                          mask if N > 1 UNLESS the pixels are consecutive in a row
    isolated columns   : columns left alone between masked ones -> mask
    N_pix[image, ccd, column_group] = count(keep)

    # stage 2: rates for the background
    n1, n2, n3 = isolated pixels with 1, 2, >=3 electrons   (per image, per group)

    # stage 3: patterns (note, Eqs. 9-13)
    for each pair (j, j+1) with both pixels >= Qmin and both kept:
        p11 < 4 ? -> p21 < 4 ? -> p31 < 4 ?          (promotion chain)
        extend 11 / 21 with a third 1 e- pixel if p111 or p211 < 5.5
        keep the highest total charge, require total < 5.6 e-
    require: rows above/below < 1 e-, same pixels in other CCDs < 1 e-,
             no other pixel >= 2 e- in the candidate columns
    count into {11,21,111,31,22,211}
```

Checks this pass must satisfy, with the numbers it has to hit:

| Check | Target |
|---|---|
| Candidate totals | (144, 0, 0, 1, 0, 0) for the published selection |
| Used pixels | 1,853,807,441 (Table II retention 0.9414 and 0.9474) |
| Isolated-pixel totals (D2, PRL) | 184 at 2 e⁻, 17 at 3 e⁻, 1 at 4 e⁻ |
| Control selection on vertical and diagonal pairs (PRL) | 139 and 130 `{11}` patterns, against 141.4 expected |

---

## Appendix B — Files

| File | Role |
|---|---|
| `apps/ccdarksens_classify_pair.cc` | hand check of `p_mn` and of the label for the chain and tree selections |
| `apps/ccdarksens_compute_brc.cc` | `B^rc_p` from rates and pixel counts |
| `apps/ccdarksens_plot_limit.cc` | curves, CSV export and ratios of limit curves per mass window |
| `configs/examples/dm_electron_heavy_lbc_rowcut_{A,B,C,D,E}.json` | scans with the modified observed counts (D and E also scale Bp and exposure) |
| `data/Final_Combined_Image_Data.csv` | per-image used pixels and published candidate counts |

---

## Appendix C — Assumptions in one place

| # | Assumption | If it is wrong |
|---|---|---|
| A1 | `B[p\|q]` is unchanged by the modified cut | Bp changes by the change of the matrix, not only by `f` |
| A2 | The modified cut gives (158, 0, 0, 2, 0, 0) as on slide 4 | every scan except the published one must be redone with the real counts |
| A3 | The published `chain` selection and signal-efficiency table are used for the modified cut | signal rows for `{22}` and `{31}` change (C3) |
| A4 | Returned pixels have the same rates as surviving ones | Bp is under- or overestimated; S6 settles it |
| A5 | Exposure scales with `1 + f` | D and E are off by the difference |
| A6 | The 90% threshold is the asymptotic 1.642 | with the toy threshold the limits are weaker |
| A7 | σ = 0.171 e⁻ for CCD D and the charges are calibrated electrons | the p-values of Step 2 shift; the label is robust unless a p_mn is near 4 |
| A8 | The event is a genuine ionisation event | the event should be treated as background, not data |
| A9 | `Bp` and `Br` for `{31}` are the published ones | the significance of the `{31}` excess changes (the audit says Bp({31}) may be 28% low) |
