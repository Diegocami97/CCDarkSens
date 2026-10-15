<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_Pattern_Pipeline_Walkthrough.md -- The pattern-space pipeline one command at a time: classify a pair, generate the signal efficiency table and the background identification matrix, check both against the references, compute Bp, wire them into a scan config, run the scan and plot the limit; what every app does, what goes in and out, the number to check, and what fails.
-->

# Student Example — The Pattern Pipeline, One Command at a Time

This is the walkthrough of every tool that stands between "I have a pattern selection and a signal model" and "I have a limit curve", starting from the generation of the signal rates. It is a full pipeline test in the sense of the other student examples, now for the pattern space. Nothing is hidden behind a driver: each command is run by hand, its output is read, and the number to check is stated. Every command and every number below was run on the checkout this document sits in.

The example that needs it is [`Student_Examples_Row_Cut_Recovered_Events.md`](Student_Examples_Row_Cut_Recovered_Events.md) (what a recovered event does to the limit). The physics behind each stage is in [`Pipeline_Raw_Data_to_Limit.md`](Pipeline_Raw_Data_to_Limit.md). The scan this walkthrough ends with is the one of [`Student_Examples_DM_Electron_Heavy_LBC.md`](Student_Examples_DM_Electron_Heavy_LBC.md).

---

## Start here: the order

Do the sections in this order. The times are what the commands took on the checkout this was written on; the reading is the long part.

| # | Section | What you do | Time |
|---|---|---|---|
| 1 | §0–§1 | read the map and the two points in it that are easy to miss | 15 min |
| 2 | §2 | build the programs; check the files that must be present | 5 min |
| 3 | §3 (Stage 0) | classify Salma's pair by hand: `{3,1}` for the chain, `{3,2}` for the tree | 1 min |
| 4 | §4 (Stage 1) | generate the signal tables and check them against the references | 15 s of running |
| 5 | §5 (Stage 2) | generate the background matrices and check them | 1 min |
| 6 | §6 (Stage 3) | compute Bp from the rates | instant |
| 7 | §7 (Stage R) | generate the rate files (605 files). **They are needed from the next step on.** | 1 s |
| 8 | §8–§10 (Stages 4–6) | read the three scan configs, run the three scans, compare them | 2 min |
| 9 | §11 | run the whole sequence in one block and compare with the CHECK values | 1 min |
| 10 | §14 (optional) | confirm four of the results with independent Python implementations | 2 min |

After step 8 every tool has been run and every number read. The next document is [`Student_Examples_Row_Cut_Recovered_Events.md`](Student_Examples_Row_Cut_Recovered_Events.md), which uses these tools on Salma's event.

---

## 0. How to read this document

Each stage has the same parts: what the program does inside, the usage line with every argument, what the input and output files look like, the check, and what goes wrong. Boxes mark the things worth stopping at.

> **WHAT YOU LEARN.** The idea the stage teaches.

> **CHECK.** A number that tells you the stage worked. If you do not get it, stop here.

> **WATCH OUT.** A way this stage fails, or fails quietly.

> **ASSUMPTION (A1, A2, …).** Something taken as true that this walkthrough cannot check.

> **INPUT NEEDED FROM SALMA (S4, S6).** Something only her analysis of the data can supply. The numbering is the one of the row-cut example.

> **NEEDED FROM THE COLLABORATION (C2).** Something neither she nor this repository has.

### Words used in this document

| Word | Meaning |
|---|---|
| Signal table P(pattern \| n_e) | For an n_e-electron deposit, the probability that the selection finds each pattern. One row per (pattern, n_e). |
| Matrix B[p \| q] | The probability that an *ideal* pattern q of dark-current electrons plus readout noise is *identified* as pattern p. Rows are identified patterns, columns are injected patterns. |
| Bp (B^rc_p) | The expected count of random coincidences in bin p: N_pix × Σ_q P(q) × B[p \| q]. |
| Selection | The rule that turns pixel charges into a pattern label: `tree` (the notebook decision tree used for the pydme tables) or `chain` (the published analysis). |
| Trials N | The number of simulated images (or clusters) per row. The statistical error of an entry is about √(p(1−p)/N). |
| Pull | (new − reference) / σ, with σ from a pooled binomial test of the two tables. A faithful reproduction has an RMS pull close to 1. |

---

## 1. The map

```
 rates.json ──► [R] qedark_generate_grid.py ──► dRdE_Si_heavy_m<mass>_s<sigma>.csv  (one file per grid point)
   (halo, mass and σ axes)                                  │
                                                            │
 config.json ──► [1] ccdarksens_pattern_signal_eff ──► Efficiencies_patterns.csv ──┐
     (selection, α, σ_ro, pile-up)                       P(pattern | n_e)          │
                                                                                   │
 config.json ──► [2] ccdarksens_pattern_background_eff ──► Background_efficiencies.csv
                                                            B[p | q]               │
                                                                 │                 │
 rates.json ──► [3] ccdarksens_compute_brc ◄──────────────────────┘                 │
   (λ1, λ2, λ3, N_pix)        │                                                     │
                              ▼                                                     │
                       background_Bp  (six numbers) ──────────┐                     │
                                                              ▼                     ▼
 observed counts  ──────────────────────────────►  [4] scan config.json  (efficiency_csv, background_Bp, ...)
                                                              │
                                                              ▼
              rate files ──►  [5] ccdarksens_scan_generic ──► scan_generic.root
                                                              │
                                                              ▼
                                         [6] ccdarksens_plot_limit ──► curves, CSV, ratios

 [0] ccdarksens_classify_pair   hand check of one pattern, outside the chain
 [c] ccdarksens_validate_pattern_efficiency / ccdarksens_validate_background_efficiency
                                check [1] and [2] against a reference table
```

Two things in this map are easy to miss:

- The matrix of stage 2 reaches the scan **only through Bp**. In the pattern-space template mode of the scan (`background_source = "bp_br_template"`), the six numbers of `run.background_Bp` are pasted into the config; the scan never reads the matrix. Changing the matrix without recomputing and pasting Bp changes nothing.
- The signal table of stage 1 reaches the scan through one config key, `response.efficiency_mc.efficiency_csv`; the rate files of stage R through another, `model.rates_dir`.

| Stage | Program | What it makes | Time measured on this checkout |
|---|---|---|---|
| 0 | `ccdarksens_classify_pair` | the label and p_mn of one pair | instant |
| R | `utils/qedark_generate_grid.py` (Python) | the rate files | 605 files in about 1 s; 240,000 files in a few minutes |
| 1 | `ccdarksens_pattern_signal_eff` | the signal table | 7 s for 10⁵ trials per n_e (n_e = 1..5) |
| 2 | `ccdarksens_pattern_background_eff` | the background matrix | 3 s for 10⁵ trials, 30 s for 10⁶ |
| c | `ccdarksens_validate_*` | a pass/fail and a table of pulls | instant for `--new` |
| 3 | `ccdarksens_compute_brc` | the six Bp values | instant |
| 4–5 | `ccdarksens_scan_generic` | the limit | 1 s for 5 masses, about 2 min for the full 800-mass grid; the toy calibration of §9.1 about 20 min for 800 masses at 10,000 toys |
| 6 | `ccdarksens_plot_limit` | curves, CSV, ratios | seconds |

---

## 2. Before you start

```bash
cd /path/to/CCDarkSens                      # every command below is run from the repository root
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j8 --target \
  ccdarksens_classify_pair ccdarksens_pattern_signal_eff ccdarksens_pattern_background_eff \
  ccdarksens_validate_pattern_efficiency ccdarksens_validate_background_efficiency \
  ccdarksens_compute_brc ccdarksens_scan_generic ccdarksens_plot_limit
mkdir -p outputs/walkthrough
```

The build type matters for time: a build with no type is not optimised and the pattern Monte Carlo runs about 14 times slower (99 s against 7 s for stage 1). The `CMakeLists.txt` now defaults to `Release`; the flag above makes it explicit. All times in this document are for `Release`.

What must be present:

| Needed by | What | Where it comes from |
|---|---|---|
| all stages | the build above | `cmake` |
| stages 1–3, validators | `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv`, `data/efficiencies_paolo.csv`, `data/Background_efficiencies.csv`, `data/Final_Combined_Image_Data.csv` | in the repository |
| stage R | Python 3 with numpy, and `python/ccdarkphys/qedark/data/Si_f2.txt` | numpy from your Python; the data file is in the repository |
| stages 4–6 | the rate files of stage R in `data/qedark_rates/Si/heavy/walkthrough_vE253p7/` | **you generate them in stage R**; they are not in the repository (`data/**` is ignored by git) |

Stages 0–3 need only the tracked files. Do stage R before the scans.

Two path rules decide whether a run works:

- **The generators and validators take paths relative to the working directory.** Run from the repository root and `configs/examples/...`, `data/...` and `outputs/...` all work. `ccdarksens_plot_limit` is the exception: it finds its literature files from the repository root whatever the working directory.
- **Paths inside a config are relative to the folder above the config file.** A config in `configs/examples/` that says `"../data/Efficiencies_patterns_….csv"` means `<root>/data/…`, because `configs/examples/` plus `..` is `configs/`, and `configs/../data` is `<root>/data`. This is the rule of `ResolveConfigRelativePath` in `ResponseFactory.cc`. A path written as `data/…` instead looks in `configs/data/…` and fails.

> **WATCH OUT.** The second rule is the usual cause of "failed to open pattern efficiency CSV". Count the `..` from the config file, not from where you are standing.

### 2.1 The configuration files of this walkthrough

Every program that needs a JSON is listed here with the file it reads. A program that takes none (stages 0 and 6) is configured by its command-line arguments only; the usage lines in the stages say which.

| Stage | Program | JSON file | What the file is | Blocks it reads |
|---|---|---|---|---|
| R | `python3 utils/qedark_generate_grid.py` | `configs/examples/qedark_generate_walkthrough.json` | five masses × 121 cross sections, heavy mediator, QEDark silicon | `material`, `mediator`, `halo`, `detector`, `rates_dir`, `filename_template`, `grid`, `options` |
| 1 | `ccdarksens_pattern_signal_eff` | `configs/examples/pattern_signal_eff_tree_pydme.json` | signal table, `tree` selection, α = 1, pile-up 4.1×10⁻⁴: the pydme reference | `detector`, `response.efficiency_mc`, `response.pattern_image`, `response.pattern_classifier` |
| 1 | same | `configs/examples/pattern_signal_eff_chain.json` | signal table, published `chain` selection, α = 0.83, pile-up 1.5×10⁻⁴ | same |
| 2 | `ccdarksens_pattern_background_eff` | `configs/examples/pattern_background_eff.json` | background matrix, `tree` selection, 3 × 5 clusters | same |
| 2 | same | `configs/examples/pattern_background_eff_chain.json` | background matrix, `chain` selection | same |
| c | `ccdarksens_validate_*` with `--config` | the four files above | the validator regenerates the table with the config of the generator | same |
| 3 | `ccdarksens_compute_brc --config` | `configs/examples/background_pattern_rates.json` | rates λ₁, λ₂, λ₃, pixel count, matrix path, patterns | `rates`, `pixels`, `migration_matrix_csv`, `patterns`, `breakdown_pattern`, `output_json` |
| 4–5 | `ccdarksens_scan_generic` | `configs/examples/dm_electron_heavy_lbc_1p3kgday.json` (+ `_toys.json`, `_mc.json`) | the published LBC scan: pydme table, published Bp and Br, observed counts, 800 × 300 grid; step 1 asymptotic, step 2 toy calibration, step 3 Monte Carlo threshold (§9.1) | `run`, `detector`, `experiment`, `response`, `model` |
| 4–5 | same | `configs/examples/walkthrough_quick_scan_reference.json` | five masses, the shipped pydme table, published Bp: the reference | same |
| 4–5 | same | `configs/examples/walkthrough_quick_scan_reference_toys.json`, `…_reference_mc.json` | the reference scan with the Monte Carlo threshold (steps 2 and 3 of §9.1) | `run.threshold_toys`, `run.q_target_lookup_path` |
| 4–5 | same | `configs/examples/walkthrough_quick_scan.json` | five masses, your chain signal table, published Bp | same |
| 4–5 | same | `configs/examples/walkthrough_quick_scan_own_tables.json` | five masses, your chain signal table **and your Bp** | same |
| 0, 6 | `ccdarksens_classify_pair`, `ccdarksens_plot_limit` | none | arguments on the command line | – |

The keys that matter in each file are explained in Appendix A. Two facts about reading them:

- **Unknown keys are ignored without a message.** A misspelled key leaves the default in place.
- **Not every key present in a file is read by every program.** The generators read only the blocks in the last column; `run.*` and most of `detector` are present in their configs because the same parser reads them, and are not used (Appendix A.1 lists exactly what is).

---

## 3. Stage 0 — Classify one pair by hand

**What it does.** Builds a three-row image whose middle row holds the pixels you give it, runs the framework's own `PatternClassifier` with both selections, and prints the statistics that decide the label. It is the same code the generators run, so there is no second implementation to disagree with.

**Usage.**

```bash
build/ccdarksens_classify_pair q1 q2 [q3 ...] (--sigma S | --bound B) [--qmax Q]
                               [--left L] [--right R] [--above a1 a2 ..] [--below b1 b2 ..]
```

| Argument | Meaning |
|---|---|
| `q1 q2 …` | calibrated charges (e⁻) of consecutive pixels of one row |
| `--sigma S` | readout noise of the CCD (default 0.16 e⁻) |
| `--bound B` | the 1 e⁻ bound of the CCD, which is 3.75σ; an alternative to `--sigma` |
| `--qmax Q` | maximum total charge of a pattern (default 5.6 e⁻) |
| `--left L`, `--right R` | charge of the pixel to the left or right of the list |
| `--above …`, `--below …` | charges of the rows above and below the pixels (default empty) |

**Run it on Salma's pair.**

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

**How to read it.**

- `cdf(q ; m)` is the chance that a pixel that really holds m electrons reads *at most* q. A pixel far above the hypothesis costs nothing; a pixel far below it is expensive (`(1,3)` and `(2,3)` give 25.55, because 1.85 e⁻ is far too small to hold 3 electrons).
- The chain climbs while the statistic stays below 4: 1+1 passes, 2+1 passes, 3+1 passes (0.17). It stops at 3+1.
- The published selection (`chain`) calls the pair `{3,1}`. The notebook tree has a separate (3,2) class, which is not one of the six search bins. **The label is a property of the selection, not of the pixels.**

A clean 1+1+1, where the third pixel is checked against 5.5:

```bash
build/ccdarksens_classify_pair 1.05 0.80 --sigma 0.16 --left 0.95
```

gives `{1,1,1}` for both selections.

> **WHAT YOU LEARN.** How a label comes out of a few lines of arithmetic, and that a table made with one selection cannot be used with the other.

> **CHECK.** `3.17 1.85 --bound 0.64` gives chain `{3,1}` and tree `{3,2}`, with p_31 = 0.17.

> **INPUT NEEDED FROM SALMA (S8).** The image ID, CCD names, calibrated charges and the σ (or 1 e⁻ bound) of the CCD for the event. The bound of 0.64 e⁻ was read off her plot.

> **CROSS-CHECK.** `utils/pattern_pvalues.py` re-implements the formulas of the note independently of the C++ classifier (§14.1).

---

## 4. Stage 1 — The signal efficiency table P(pattern | n_e)

**What the program does inside.** For each n_e = 1..5 it generates N images. Each image is built like this:

1. `n_e` electrons are placed at a random depth z in the 669 µm CCD; each one diffuses to a Gaussian spread σ_xy(z) = √(−A ln(1 − b z)) × (α + β·E), with A = 803.25 µm², b = 6.5×10⁻⁴ µm⁻¹ (E = 0 here, so only α matters);
2. the cloud centre is drawn at random inside the central columns of a 3 × 50 image (the cloud is binned over 100 rows, so three rows are three binned rows);
3. Poisson dark-current pile-up with mean `lambda_dc` per binned pixel and Gaussian readout noise with σ_ro = 0.16 e⁻ are added;
4. the selection scans the middle row and counts every pattern it finds (an image can hold more than one, for example two lone electrons, so a row of probabilities can add to more than 1).

The result is the fraction of images in which each (pattern, n_e) was found.

**Usage.**

```bash
build/ccdarksens_pattern_signal_eff CONFIG.json [Nsims=100000] [outdir=.] [ne_max=5]
```

**The two configs.** They differ in exactly the fields below; everything else is the same.

| Field | `pattern_signal_eff_tree_pydme.json` | `pattern_signal_eff_chain.json` |
|---|---|---|
| `response.efficiency_mc.alpha` | 1.0 | 0.83 |
| `response.pattern_image.lambda_dc` | 0.00041 | 0.00015 |
| `response.pattern_classifier.selection` | `tree` | `chain` |

The first reproduces the collaboration's pydme table; the second is the published selection with the diffusion and pile-up of the analysis note.

**Run both.**

```bash
mkdir -p outputs/walkthrough/signal_tree outputs/walkthrough/signal_chain
build/ccdarksens_pattern_signal_eff configs/examples/pattern_signal_eff_tree_pydme.json  100000 outputs/walkthrough/signal_tree
build/ccdarksens_pattern_signal_eff configs/examples/pattern_signal_eff_chain.json       100000 outputs/walkthrough/signal_chain
```

Output file (`Efficiencies_patterns.csv`):

```
# Simulation Parameters: Nsims=100000, alpha=1, beta=0, A=803.25,b=0.00065,lambda=0.00041, selection=tree
pattern,ne,Efficiency
1,1,0.98652
2,1,0.00041
...
```

The tree run writes 15 patterns, the chain run writes 11, because the chain has no 32, 41, 311 or 221 class.

**Check 1 — the tree table against the pydme reference.**

```bash
build/ccdarksens_validate_pattern_efficiency data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv \
  --new outputs/walkthrough/signal_tree/Efficiencies_patterns.csv
```

```
signal: 80 entries, 49 non-zero in either table
  RMS pull (non-zero entries): 1.02   chi2/n: 1.03
  entries beyond 3 sigma: 0 (expected 0.1)   beyond 2 sigma: 3 (expected 2.2)   max |pull|: 2.60
  per n_e:  sum(ref)  sum(new)  RMS pull
     1     0.9874   0.9877      0.43
     2     1.2151   1.2150      1.58
     3     1.2788   1.2763      1.28
     4     1.2343   1.2330      0.80
     5     1.1385   1.1374      0.78
  RESULT: PASS (RMS pull <= 1.5 and no entry beyond 5 sigma)
```

**Check 2 — the chain table against the published efficiencies.**

```bash
build/ccdarksens_validate_pattern_efficiency data/efficiencies_paolo.csv \
  --new outputs/walkthrough/signal_chain/Efficiencies_patterns.csv --patterns 11 21 22 31 111 211
```

```
  per n_e:  sum(ref)  sum(new)  RMS pull
     2     0.3837   0.3882      1.10
     3     0.6499   0.6564      1.15
     4     0.7872   0.7932      2.79
     5     0.8630   0.8632     11.12
  largest pulls:
    (111, 5): ref 0.027699  new 0.01557  pull -18.64
    (31, 5): ref 0.4893  new 0.52471  pull +15.84
```

Read this in two parts, the statistical one and the physical one.

**Statistically** the n_e = 5 row is inconsistent: the pull test says so with large numbers, because the reference has 10⁵ trials and an entry near 0.5 has an error of only about 0.0016. The sums over the six search patterns, which are what the PRL quotes ((38, 65, 79, 86)%), agree: at n_e = 5, 0.8632 against 0.8630. The n_e = 2 and 3 rows agree entry by entry. At n_e = 5 the charge is distributed differently among the patterns: `{111}` is 0.016 against 0.028, `{31}` is 0.525 against 0.489, `{211}` 0.261 against 0.277, `{11}` 0.011 against 0.015. The row `P({11} | 1 e⁻)` is exactly 0 in the published file and 0.0003 here. The exact published rules beyond the note are not known.

**Physically** the difference is small. The largest entry difference is 0.036 in absolute terms (`{31}`). Running the same five-bin scan on 25 masses from 1 to 1000 MeV with the published table instead of the chain table, with the published Bp, Br and counts, moves the limit by **0.7% at the median and 1.2% at most** (the chain-table limit is 0.7% lower at the high masses and 1% lower below 3 MeV). For comparison, using the pydme tree table instead moves it by up to +7.6% at high mass, and a 0.25 dex cross-section grid by up to 7%. So the validator's FAIL is a statement about statistical consistency of a simulation with a published table, not about the size of the effect on the result.

**What the validator does with the flags.**

| Form | What it does |
|---|---|
| `validate_pattern_efficiency [N] [ref.csv]` | runs its own built-in simulation (the constants fixed in the code) and compares |
| `validate_pattern_efficiency [N] [ref.csv] --config CONFIG.json` | generates the table with the same code and parameters as stage 1 and compares, so no separate numbers can drift |
| `validate_pattern_efficiency ref.csv --new t1.csv [t2.csv …]` | does not simulate; compares existing tables (averaged if several) |
| `--n-ref N`, `--n-new N` | trials behind each entry; default is the `Nsims=` value in the CSV header |
| `--patterns c1 c2 …` | restrict to these pattern codes |
| `--out file.csv` | write the per-entry comparison |

> **WHAT YOU LEARN.** What the signal table is made of (diffusion, random centre, pile-up, noise, selection), and that "reproduces the reference" has a precise meaning: RMS pull near 1.

> **CHECK.** Tree table vs pydme reference: RMS pull ≈ 1 (1.02 here). Chain table vs published: six-pattern sums within 1.5% for n_e = 2..5.

> **WATCH OUT.** (1) The statistical error of an entry shrinks as 1/√N. At 10⁵ trials a 0.003 entry has a ±0.0002 error; do not read differences below that. (2) A table is only valid with the selection that made it. (3) The `Efficiency` values are not normalised to 1: they are the number of times a pattern was found per image.

> **ASSUMPTION (A3).** The chain implementation follows the data-selection note and the PRL. The n_e = 5 differences show it is not identical to the collaboration's code.

> **CROSS-CHECK.** `utils/compare_tables.py` repeats the table comparison independently (§14.2).

---

## 5. Stage 2 — The background identification matrix B[p | q]

**What the program does inside.** For each of the 26 ideal patterns q (the empty pattern, 1..5 e⁻ in one pixel, and every two- or three-pixel combination with total charge ≤ 5) it simulates N noisy 3 × 5 clusters: the ideal charges plus Gaussian readout noise, no diffusion, no dark current. It classifies each cluster and counts what comes out. A cluster where nothing is found counts as the empty pattern `{0}`. Dividing by N gives B[p | q].

**Usage.**

```bash
build/ccdarksens_pattern_background_eff CONFIG.json [Nsims=100000] [outdir=.]
```

The configs are `pattern_background_eff.json` (tree) and `pattern_background_eff_chain.json` (chain). The program warns if the noise or thresholds differ from the values the shipped table was validated with: σ_ro = 0.16, Qmin = 0.6, Qmax = 5.6, thresholds 3.5 / 4.0 / 5.5.

**Run it.**

```bash
mkdir -p outputs/walkthrough/background_tree
build/ccdarksens_pattern_background_eff configs/examples/pattern_background_eff.json 100000 outputs/walkthrough/background_tree
```

Output file (`Background_efficiencies.csv`) has one row per *identified* pattern (string-sorted: `0, 1, 11, 111, 112, …`) and one column `eff_<injected>`:

```
iden_pat,eff_0,eff_1,eff_2,...,eff_311
0,9.999700e-01,3.058000e-02,3.700000e-04,...
```

**Check — against the shipped matrix.**

```bash
build/ccdarksens_validate_background_efficiency data/Background_efficiencies.csv \
  --new outputs/walkthrough/background_tree/Background_efficiencies.csv --n-new 100000
```

```
background: 676 entries, 180 non-zero in either table
  RMS pull (non-zero entries): 1.00   chi2/n: 1.00
  entries beyond 3 sigma: 0 (expected 0.5)   beyond 2 sigma: 4 (expected 8.2)   max |pull|: 2.72
```

The background matrix has no `Nsims=` in its header, so `--n-new` is needed; the reference is taken as 10⁵ trials per entry.

The chain matrix is *not* expected to reproduce the shipped one:

```bash
build/ccdarksens_pattern_background_eff configs/examples/pattern_background_eff_chain.json 100000 outputs/walkthrough/background_chain
```

(RMS pull 135 against the shipped matrix.) The shipped matrix has rows 32, 41, 311 and 221; the chain has no such classes, so a (3,2) cluster is identified as `{31}` in the chain matrix and as `{32}` in the shipped one. Whichever selection the signal table used, the matrix must come from the same one.

> **WHAT YOU LEARN.** The matrix is a property of the selection and the noise, not of the dark-current rate. The same matrix is used for any set of rates.

> **CHECK.** `B[{11} | {11}] ≈ 0.908` and RMS pull ≈ 1 against the shipped matrix (the tree run gives 0.9095 at 10⁵ trials, 0.9083 at 10⁶).

> **WATCH OUT.** Rare entries need many trials, and they matter. `B[11 | 1]`, the probability that a *lone* 1 e⁻ pixel is called `{11}`, is about 10⁻⁵. There are 5×10⁵ lone pixels, so it adds 5 to 7 counts to the 141 of `{11}`. At 10⁵ trials that entry is one or zero counts. Use at least 10⁶ trials for anything that goes into a limit (30 s), and do not interpret a single rare entry.

---

## 6. Stage 3 — Bp from the rates

**What the program does inside.** For each of the six bins p it sums over the 26 injected patterns q:

```
B^rc_p = N_pix × Σ_q  P(q) × B[p | q]          P(q) = Π_i Poisson(q_i | λ_{q_i})
```

where `λ_k` is the mean number of pixels with k electrons per pixel readout (k = 1, 2, and 3 or more share one rate for 3, 4, 5), the empty pattern has P = 1, and `N_pix` is the number of pixel readouts. It prints the total, the diagonal (q = p) part and the misidentified (q ≠ p) part.

**Usage.**

```bash
build/ccdarksens_compute_brc [theta1] [theta2] [theta3] [N_total] [beff_csv] [image_csv]
build/ccdarksens_compute_brc --config rates.json [--lambda1 L] [--lambda2 L] [--lambda3 L]
                             [--n-pixels N] [--segments table.csv] [--matrix B.csv] [--out f.json]
```

`theta_k = λ_k × 10⁴`, so `theta = 3` means λ = 3×10⁻⁴ per pixel. A flag overrides the matching entry of the config.

**Run it with the shipped matrix** (this is the notebook check):

```bash
build/ccdarksens_compute_brc 3 2 2
```

```
pattern           B^rc_p    diagonal q=p   misidentified
11            1.5694e+02      1.5138e+02      5.5631e+00
21            2.3512e-02      1.0083e-02      1.3429e-02
111           4.7022e-02      4.3687e-02      3.3352e-03
31            1.5683e-06      6.7328e-07      8.9500e-07
22            1.0073e-06      6.7381e-07      3.3353e-07
211           1.1137e-05      3.0275e-06      8.1100e-06
  {"background_Bp":[156.942,0.0235117,0.0470225,1.56828e-06,1.00734e-06,1.11374e-05]}
```

These are the numbers the `DataSimulator.ipynb` notebook of pydme prints for θ = (3, 2, 2).

**Run it from a config, with the matrix you generated:**

```bash
build/ccdarksens_compute_brc --config configs/examples/background_pattern_rates.json \
  --matrix outputs/walkthrough/background_tree/Background_efficiencies.csv \
  --out outputs/walkthrough/background_Bp_tree.json
```

The config is `configs/examples/background_pattern_rates.json` with blocks `rates`, `pixels`, `migration_matrix_csv`, `patterns`, `breakdown_pattern`, `output_json`. With θ = (2.9, 2, 2) the six values are

| Matrix | {11} | {21} | {111} | {31} | {22} | {211} |
|---|---|---|---|---|---|---|
| shipped | 146.84 | 0.02262 | 0.04258 | 1.516e-6 | 9.96e-7 | 1.038e-5 |
| generated, tree, 10⁵ trials | 147.09 | 0.02143 | 0.04735 | 1.432e-6 | 6.72e-7 | 1.187e-5 |
| generated, tree, 10⁶ trials | 148.52 | 0.02148 | 0.04266 | 1.457e-6 | 7.92e-7 | 1.236e-5 |

The `{11}` diagonal is 141.5 to 141.7 in all three (this is the published 141.4); the differences are in the misidentified column, which depends on a handful of counts (Stage 2, "watch out").

**Use per-segment rates.** The paper sums over every image and every 10-column group. `--segments table.csv` takes a CSV with columns `n_pixels`, `lambda_1` and optionally `lambda_2`, `lambda_3plus`, one row per segment, and sums over the rows with their own rates. With a single segment this is the uniform formula above.

> **WHAT YOU LEARN.** Bp is a product of a pixel count, Poisson probabilities from measured rates, and a matrix. Change a rate and it moves as a power: `{11}` as λ₁², `{111}` as λ₁³.

> **CHECK.** θ = (3, 2, 2) with the shipped matrix gives 156.942, 0.0235117, 0.0470225, 1.56828e-06, 1.00734e-06, 1.11374e-05.

> **WATCH OUT.** Use the matrix of the same selection as the signal table. With the chain signal table the Bp must come from the chain matrix (rows 32, 41, 311, 221 are empty there and their mass goes to `{31}` and `{211}`); with the pydme tree table, from the tree matrix.

> **WATCH OUT.** The published Bp is (141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5). This app does not return it: `{11}` matches the diagonal term only, and the other five bins depend on per-band populations and non-Poisson 2 e⁻ and 3 e⁻ pixel counts (184 and 17 in D2) that a uniform rate cannot represent. Use the published numbers for any result meant to reproduce the paper, and this app to understand and vary the background.

> **CROSS-CHECK.** `utils/background_pattern_counts.py` computes the same sum independently (§14.3).

> **INPUT NEEDED FROM SALMA (S4, S6).** `N_pix` after her modified cut (S4), and the isolated 1, 2, 3 e⁻ counts per image and per 10-column group, for the original and the modified cut (S6). Without them, scale Bp by (1 + f) as in the row-cut example.

---

## 7. Stage R — The signal rate tables

**What they are.** The scan folds a differential rate dR/dE for every (m_χ, σ_e) point. These rates are the physics input of the whole chain: DM–electron scattering in silicon with the QEDark form factor (`python/ccdarkphys/qedark/`, data file `Si_f2.txt`), the standard halo of Baxter et al. (v₀ = 238, v_E = 253.7, v_esc = 544 km/s; see A10), a heavy mediator, a band gap of 1.2 eV, 3.8 eV per electron–hole pair, and 0.1 eV energy bins. One CSV is written per (m_χ, σ_e) point.

**This is the one Python program of the pipeline.** The rates come from a Python physics code (a port of QEDark), not from a C++ app. Nothing else in this walkthrough needs Python.

**Usage.**

```bash
python3 utils/qedark_generate_grid.py CONFIG.json          # needs Python 3 and numpy; run from the repository root
```

**The config** is `configs/examples/qedark_generate_walkthrough.json` (Appendix A.4 explains every key). It is a copy of the dense configuration of the published scan with three changes: only the five masses of the quick scans, a cross-section axis reduced to 10⁻⁴² to 10⁻³⁴ cm² with 121 points (0.067 dex steps, the step of the published grid), and the output folder `data/qedark_rates/Si/heavy/walkthrough_vE253p7`.

**Run it.**

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_walkthrough.json
```

```
[grid] wrote data/qedark_rates/Si/heavy/walkthrough_vE253p7/dRdE_Si_heavy_m199.961135_s1.0e-34.csv
```

It writes 605 files (5 masses × 121 cross sections) in about 1 s on 4 processes. One file looks like this:

```
# Differential Rates computed with CCDarkSens (QEDark entry)
# material = Si, mediator = heavy, table = .../python/ccdarkphys/qedark/data/Si_f2.txt
# halo (cm/s): v0=23800000.0, vE=25370000.0, vesc=54400000.0
# mX (eV) = 5995706.0
# sigma_e (cm^2) = 1e-38
# Output units: dR/dE in events / kg / year / eV
# Columns: E (eV), dRdE (events/kg/year/eV)
E,dRdE
1.2,116.5229741
1.3,430.3625329
...
```

**The file name is the key.** The scan does not search the folder: it builds the name `dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv` from the config, with the mass to six decimals (`m5.995706`) and the cross section to two significant digits (`s1.0e-38`). A mass or cross-section value in the scan config that does not reproduce an existing file name finds no rate.

> **WHAT YOU LEARN.** The physics input is a table of files named by (m_χ, σ_e). The scan config's grid decides which names it asks for, so the two grids must be generated and read with the same numbers and the same format.

> **CHECK.** `ls data/qedark_rates/Si/heavy/walkthrough_vE253p7 | wc -l` gives 605; `head` of a file shows `# Columns: E (eV), dRdE (events/kg/year/eV)`.

> **WATCH OUT.** (1) `options.skip_existing` is true: running the generator again does not rewrite an existing file. If you change the halo or the form factor, set `overwrite` to true or use a new `rates_dir`. (2) A coarse cross-section axis changes the limit. With 33 points (0.25 dex steps) the reference scan below came out up to 7% away from the published one; with 121 points it is within 0.5%. Keep the step at 0.07 dex or finer. (3) The full published scan uses 800 masses × 300 cross sections = 240,000 files, `configs/examples/qedark_generate_heavy_dense.json`, written to `data/qedark_rates/Si/heavy/long_scan_vE253p7` (about 8 GB; it took a few minutes on 8 processes). A folder `long_scan` made with the older v_E = 263 km/s may sit next to it; the scan configs of this walkthrough read `long_scan_vE253p7` only.

> **ASSUMPTION (A10).** The halo is the one recommended by Baxter et al. (arXiv:2105.00599, Table 1 and Eqs. 11–12), which the PRL cites: ρ = 0.3 GeV/cm³, v₀ = 238 km/s, v_esc = 544 km/s, and the speed of the laboratory in the galactic frame, v_E = |v₀ + v_⊙ + v_⊕(March 9)| = |(0, 238, 0) + (11.1, 12.2, 7.3) + (29.2, −0.1, 5.9)| = |(40.3, 250.1, 13.2)| = 253.7 km/s. Rate tables made with v_E = 263 km/s (the default of older pydme configs) raise the rate in the low-mass tail and make the limit up to 30–40% too strong below 10 MeV. The band gap (1.2 eV), 3.8 eV per pair and 0.1 eV bins are the ones of the published analysis; I did not re-derive the rates here.

---

## 8. Stage 4 — Wire the results into a scan config

The scan reads a config. Three entries connect it to what you just made; everything else is the published LBC configuration.

| Config key | Takes | Written by |
|---|---|---|
| `model.rates_dir` | the folder of rate files | stage R |
| `response.efficiency_mc.efficiency_csv` | the signal table, path relative to the folder above the config | stage 1 |
| `run.background_Bp` | six numbers in the order of `experiment.pattern_roi` (11, 21, 111, 31, 22, 211) | stage 3 |
| `run.observed_counts` | six counts in the same order | the selection run on data |

Three ready-made configs are provided. All are copies of `dm_electron_heavy_lbc_1p3kgday.json` with these changes:

- `model.rates_dir` is `data/qedark_rates/Si/heavy/walkthrough_vE253p7` (stage R);
- `run.one_sided_q` is `true` (q is 0 at and below the best-fit cross section; §9) and `run.robust_minimum` is `true` (the reference NLL minimum cannot be above the grid minimum; §9);
- `model.grid.mchi_MeV` is the five masses `{"values": [1.99987, 5.995706, 19.997409, 59.953172, 199.961135]}` and `model.grid.sigma_e_cm2` is `logspace(-42, -34, 121)`, the same numbers as the generator config;
- `run.label` and `run.outdir` are set so the output does not overwrite the published scan.

| Config | `efficiency_csv` | `run.background_Bp` | Meaning |
|---|---|---|---|
| `walkthrough_quick_scan_reference.json` | the shipped pydme table | published (141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5) | the five-mass version of the published scan; the reference for the other two |
| `walkthrough_quick_scan.json` | **your** chain table (stage 1) | published | only the signal table is yours |
| `walkthrough_quick_scan_own_tables.json` | **your** chain table | **yours**: 148.539, 0.0214786, 0.0472488, 1.4576e-06, 7.9154e-07, 9.99519e-06 (the chain matrix of stage 2 at 10⁶ trials, θ = (2.9, 2, 2), through stage 3) | signal table **and** Bp |

Open the three and diff them. To use your own Bp in the second file, replace `run.background_Bp` by the line printed by `ccdarksens_compute_brc`.

All three keep the published `background_Br`, `observed_counts` = (144, 0, 0, 1, 0, 0) and the published exposure, so only what you change moves.

The limits at the five masses (cm²), grid-crossing values from the scan log:

| m_χ (MeV) | published scan (800 masses) | reference quick scan | chain table | chain table and own Bp |
|---|---|---|---|---|
| 1.99987 | 1.690e-36 | 1.683e-36 | 1.647e-36 | 1.070e-36 |
| 5.99571 | 6.716e-39 | 6.692e-39 | 6.858e-39 | 6.881e-39 |
| 19.9974 | 4.908e-39 | 4.891e-39 | 4.708e-39 | 4.706e-39 |
| 59.9532 | 1.014e-38 | 1.016e-38 | 9.503e-39 | 9.518e-39 |
| 199.961 | 2.998e-38 | 3.001e-38 | 2.789e-38 | 2.795e-38 |

How to read it:

- **The reference reproduces the published scan** to 0.4% or better at all five masses (these are asymptotic limits with q = 1.642; §9 adds the Monte Carlo threshold). This is the check that stage R and the quick-scan configs are right.
- **The chain table** moves the limit by −7% to +3% against the pydme table. The signal efficiencies of the two selections differ most at n_e = 4 and 5.
- **Your Bp** matters at 2 MeV, where `{11}` drives the limit: your Bp for `{11}` is 148.5, above the 144 observed, so the data sit below the background and the limit becomes 36% stronger than the reference. The published Bp is 141.4, below the 144 observed. At the other masses your Bp changes the limit by 1% or less, because the `{11}` background is not what limits there.

> **WATCH OUT.** (1) **The masses must be points of the rate grid.** The rate files are looked up by mass with six decimals. A mass that is not in the folder gives no signal, and the scan reports `upper_limit=1e-26`, the top of the cross-section grid, with no error message. The published grid is `logspace(0.2, 1000, 800)`. (2) A misspelled config key is ignored silently; the scan runs with the default. (3) `observed_counts` and `background_Bp` must have as many entries as `pattern_roi`. (4) **A single mass is not a smooth number.** The crossing at 200 MeV for the chain table moved from 3.1e-38 to 2.6e-38 when the cross-section grid was shifted by half a step, a 15% change with the same physics. The limit is a staircase in n_e and has features in σ; quote window medians over many masses (stage 6 on the full scan), not five points.

---

## 9. Stage 5 — Run the scan

**What the program does inside.** For each (m_χ, σ_e) grid point it loads the rate file, folds it through the charge response into n_e bins and then into pattern bins with the signal table, adds the background `Bp + θ·Br`, and minimises the profile likelihood over the radiogenic nuisance θ. The statistic q(σ) is built over the cross-section grid and the 90% limit is where q crosses the threshold, 1.642 for the asymptotic one-sided test.

**Usage.**

```bash
build/ccdarksens_scan_generic CONFIG.json
```

**Run the three quick scans.**

```bash
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_reference.json
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan.json
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_own_tables.json
```

Each takes about a second. The reference prints:

```
[scan-generic] analysis_space=pattern observable_bins=pattern exposure=0.00355913 kg*year
[scan-generic] mchi=1.99987 -> upper_limit=1.68302e-36
[scan-generic] mchi=5.99571 -> upper_limit=6.69171e-39
[scan-generic] mchi=19.9974 -> upper_limit=4.89089e-39
[scan-generic] mchi=59.9532 -> upper_limit=1.01551e-38
[scan-generic] mchi=199.961 -> upper_limit=3.00072e-38
[scan-generic] Done. Output written to outputs/walkthrough/scan_quick_reference/scan_generic.root
```

**The output file** `scan_generic.root` holds `q_mchi_sigma` (the statistic on the grid), `upper_limit_sigma_e_mchi` and `upper_limit_sigma_e_mchi_graph` (the limit from the grid crossing), `upper_limit_sigma_e_mchi_pydme_bisection` (the continuous limit; present when the profile minimiser is `pydme`) and `exposure_kg_year`.

**Run the full published scan** (800 masses × 300 cross sections, about 2 min). It needs the dense rate grid of stage R, 240,000 files:

```bash
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_heavy_dense.json
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday.json
```

### 9.1 Checking the threshold with toys (the Poisson regime)

The limit of the scans above is the result of this walkthrough. It crosses q(σ) with the asymptotic value 1.642, which is accurate when the limit is set by a bin with many counts (`{11}`, below 3 MeV and above 10 MeV). It is too strong when the limit is set by bins that expect about 0.1 events and saw none (`{21}`, `{111}`): the 90% threshold for a counting experiment with zero observed needs about 2.3 signal events, and q = 1.642 corresponds to about 0.8. Between 3 and 7 MeV the signal of n_e = 3 lands in exactly those bins. The published analysis built the distribution of the test statistic with toy Monte Carlo (PRL, "Results"), and Baxter et al. (arXiv:2105.00599, §2.3) recommend that asymptotic results be checked with toys. The app can do that check in two more runs, with the same config except for three keys; you do not need it for the other stages.

1. **Step 1** (what you ran above, `…reference.json`): the scan, with `run.one_sided_q = true` and `run.robust_minimum = true`. This is the headline limit.
2. **Step 2** (`…reference_toys.json`, `run.mode = "threshold_toys"`): at each mass, 10,000 Poisson toys of the expected counts at the step-1 limit σ_thr, each fitted with σ and θ free; q of a toy is 2(NLL at σ_thr − NLL at the best fit), set to 0 when the best fit lies above σ_thr (one-sided). The 90th percentile of those q values is written to `qtarget_threshold.root`.
3. **Step 3** (`…reference_mc.json`, `run.q_target_lookup_path`): the scan again, with that value as the threshold at each mass. This is the check of step 1, not a replacement for it.

```bash
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_reference_toys.json
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_reference_mc.json
```

```
[scan-generic] threshold_toys mchi=1.99987 sigma_thr=1.68302e-36 q_target(p=0.9)=1.7553
[scan-generic] threshold_toys mchi=5.99571 sigma_thr=6.69171e-39 q_target(p=0.9)=2.86719
[scan-generic] threshold_toys mchi=19.9974 sigma_thr=4.89089e-39 q_target(p=0.9)=1.63944
[scan-generic] threshold_toys mchi=59.9532 sigma_thr=1.01551e-38 q_target(p=0.9)=1.62159
[scan-generic] threshold_toys mchi=199.961 sigma_thr=3.00072e-38 q_target(p=0.9)=1.59292
...
[scan-generic] mchi=1.99987 -> upper_limit=1.73615e-36
[scan-generic] mchi=5.99571 -> upper_limit=1.0146e-38
[scan-generic] mchi=19.9974 -> upper_limit=4.88595e-39
[scan-generic] mchi=59.9532 -> upper_limit=1.00886e-38
[scan-generic] mchi=199.961 -> upper_limit=2.95393e-38
```

At 20, 60 and 200 MeV (`{11}` with many counts) the toy threshold is 1.59–1.64, the asymptotic value, and the limit moves by at most 3%. At 6 MeV it is 2.87 and the limit rises by 52%. At 2 MeV, where `{11}` decides, it is 1.76 (+3%).

**The published scan and its check** (800 masses; about 2 min, 16 min and 2 min on 8 cores):

```bash
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday.json        # step 1
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday_toys.json   # step 2
build/ccdarksens_scan_generic configs/examples/dm_electron_heavy_lbc_1p3kgday_mc.json     # step 3
```

The headline limit is `outputs/dm_electron_heavy_lbc_1p3kgday/scan_generic.root` (step 1); `outputs/dm_electron_heavy_lbc_1p3kgday_mc/scan_generic.root` is the same scan with the Monte Carlo threshold.

**How close this is to the published curve.** Ratio of the limit of this scan to the published DAMIC-M curve (the collaboration's paper-export file `DAMIC-M_2025_QEDark_DMe_heavymediator.txt`, 25 masses), median per window. The config's exposure is 1.274 kg-day against 1.257 kg-day in the PRL, which makes the limit 1.3% lower on its own.

| m_χ window (MeV) | masses | asymptotic (step 1, headline) | Monte Carlo (step 3, check) |
|---|---|---|---|
| 0.7 – 1.1 | 5 | 0.80 | 0.80 |
| 1.1 – 3 | 4 | 0.96 | 0.96 |
| 3 – 6 | 3 | 0.73 | 0.93 |
| 6 – 12 | 3 | 0.92 | 1.32 |
| 12 – 100 | 5 | 0.97 | 0.96 |
| 100 – 1100 | 5 | 0.97 | 0.94 |

On the dense 800-mass grid the Monte Carlo limit divided by the asymptotic one has a median of 1.00 below 3 MeV, **1.44 between 3 and 10 MeV** (16–84%: 1.24–1.49) and 0.99–1.01 above 10 MeV.

How to read it:

- **Above 12 MeV and between 1 and 3 MeV** both thresholds reproduce the published curve to 4–6%; the Monte Carlo step does nothing there, as it should.
- **The 0.75 MeV point** (ratio 0.37) is on the steep edge of the curve, where the limit changes by two orders of magnitude over 10% in mass. It is the same for both thresholds and is not a threshold effect.
- **3–6 MeV:** the asymptotic limit is 27% below the published one and the Monte Carlo limit 7% below. This is the effect the toys correct.
- **6–12 MeV:** the Monte Carlo limit is 17–34% *above* the published curve (1.12 at 5.5 MeV, 1.32 at 7.2 MeV, 1.34 at 9 MeV, 1.17 at 11.8 MeV) where the asymptotic limit was within 8%. So the published curve lies between the two thresholds there. I do not know which choice in the published toy procedure accounts for it; candidates not tested are the handling of the discreteness of q for few counts (the toy threshold jumps between 2.25 and 2.95 at neighbouring masses near 9–10 MeV) and how the nuisance θ is fixed in the toys. Read it this way: the step-1 limit is the result; between 3 and 10 MeV the asymptotic threshold understates the limit by a median 1.44 (up to 1.5), and the toy result is the better-motivated number there. Quote both if the window matters, and say which threshold a number uses.

> **WATCH OUT.** (1) **Plot step 3 with `--use-stored-limit`.** `--from-qhist` (the default) rebuilds every curve from the stored q-map at the asymptotic 1.642 and would draw the step-1 limit again; the stored limit of the `_mc` file is the one with the per-mass threshold. (2) The toy cross section is the limit of step 1, which is lower than the final one by up to 30% at a few masses. Repeating step 2 with the step-3 file as `sigma_threshold_graph_path` changed the limits by 0.5% in the median in a test at 25 masses; one pass is enough for a table, not for a single mass. (3) Masses with no signal (below about 0.7 MeV) have `sigma_thr = 1e-26` and `q_target = 0`; a zero threshold is replaced by the asymptotic value in the scan, and the limit there stays at `1e-26`, no limit. (4) The **two-sided** statistic (`run.threshold_toys.one_sided = false`, which is the default of the key) is the one Baxter et al. recommend; in a 25-mass test it gave limits 20–35% weaker than the published curve above 1 MeV, so it is not what the published curve shows.

> **WHAT YOU LEARN.** A limit is a crossing of q(σ) through a threshold. Everything upstream (rates, signal table, Bp, counts) enters through the expected counts in the six bins.

> **CHECK.** The reference quick scan at m_χ = 5.99571 MeV gives `upper_limit=6.69171e-39`, within 0.4% of the published scan (`6.71562e-39`); at all five masses within 0.5%.

> **WATCH OUT.** `upper_limit=1e-46` is the *bottom* of the cross-section grid, not a limit. It appears when q(σ → 0), the improvement of the best fit over the background-only fit, is already above the threshold, so there is no crossing. With the default q construction (`run.one_sided_q` false) the scan forces q to be non-decreasing, so a fit that prefers a small signal removes the limit at that mass; the row-cut example shows this on scan A of the older configs, where the largest q(σ → 0) is only 1.81 against 1.642. The configs of this walkthrough set `run.one_sided_q` to true, which puts q = 0 at and below the best-fit cross section, so the limit is the first crossing above it and the vertical lines do not appear. They also set `run.robust_minimum` to true: the reference minimum is then the lowest of the 2D fit, the grid minimum and the background-only fit. The 2D Minuit fit sometimes stops 0.1–0.2 in NLL above the true minimum, which moves q for the whole mass and makes the limit jump by up to 35% between neighbouring masses (seen in the {31} = 2 scans of the row-cut example above 100 MeV); the robust minimum removes that and only ever lowers a limit. Do not read the stored value `1e-46` as a limit if it does appear. `upper_limit=1e-26` is the *top* of the grid: no signal was found, usually because the mass or the rates folder is wrong (stage R). At masses below about 0.7 MeV the same value is correct: the deposit cannot make two electrons, and the paper's curve ends there too.

> **ASSUMPTION (A6).** The headline limit uses the asymptotic threshold 1.642. It is accurate where the limit is set by a bin with many counts and too strong at 3–10 MeV, where it is set by bins that expect about 0.1 events and saw none; §9.1 checks it with a Monte Carlo threshold, which reproduces the published curve at 3–6 MeV and overshoots it at 6–12 MeV.

---

## 10. Stage 6 — Plot and compare

**Usage.**

```bash
build/ccdarksens_plot_limit --batch --heavy [--cl 0.90] [--show-damic]
    [--ratio-to N|LABEL] [--sigma-floor X] [--x-min X --x-max X --y-min Y --y-max Y]
    [--out-pdf F] [--out-root F] [--out-csv F]
    scan1.root "label1" [scan2.root "label2" ...]
```

| Flag | Effect |
|---|---|
| `--batch` | no window; write the files and exit |
| `--heavy`, `--light`, `--intermediate` | which mediator's literature curves to overlay |
| `--cl X` | **the confidence level of every plotted limit; default 0.90.** One-sided asymptotic, q = NormQuantile(X)² (0.90 gives 1.6424, the value of the scan). It is printed in the log, written in the CSV header and drawn on the plot as "90% C.L." |
| `--from-qhist` | the default for DM-electron plots: every curve is rebuilt from the q(m_χ, σ) histogram stored in its own scan file, at `--cl`, so all curves are computed the same way. Where q at σ → 0 is already above the threshold there is no crossing and the curve is drawn on the floor of the σ grid |
| `--use-stored-limit` | use the limit stored in the scan file instead; it is at the C.L. the scan was run with (`run.cl`) |
| `--q-threshold Q` | a threshold other than the one of `--cl` (the old positional number does the same). Accepted but flagged: the log warns and the plot is labelled "q_thr = Q (…% C.L. one-sided)", never "90%" |
| `--x-min`, `--x-max`, `--y-min`, `--y-max` | axis ranges; use them to zoom, or to keep floor values off the plot |
| `--draw-both` | additionally draw the continuous bisection limit (used for the cross-check of §14.4; not needed for the plots) |
| `--show-damic` | overlay the published DAMIC-M curve |
| `--ratio-to N` or a label | print, per mass window, the median (16–84%) ratio of every other curve to input N (1-based) or to the input with that label |
| `--sigma-floor X` | the lowest cross section of the scan grid (default 1e-46); points at it are reported as "no crossing" |
| `--out-csv F` | write every curve as `label,mchi_MeV,sigma_e_cm2` |
| `--legacy-migdal-canvas` | also write the old Migdal projection canvas (off by default; it writes `dme_migdal_limit_curve.*` in the working directory) |

**One confidence level for every plot.** All limits are drawn at 90% C.L. unless `--cl` says otherwise. Older commands ended with `1.642374415149816 heavy`; the number is the 90% value and is now derived from `--cl`, and `--heavy` already gives the mediator, so both can be left out. A different number (for example `2.71`, which is 95.0% one-sided) is still accepted, but the log prints `WARNING … not 90%` and the label on the plot shows the threshold, so a plot can no longer carry "90%" with another threshold by accident.

**Compare the three quick scans.**

```bash
build/ccdarksens_plot_limit --batch --heavy --ratio-to reference \
  --out-pdf outputs/walkthrough/limit_compare.pdf \
  --out-root outputs/walkthrough/limit_compare.root \
  --out-csv outputs/walkthrough/limit_compare.csv \
  outputs/walkthrough/scan_quick_reference/scan_generic.root "reference" \
  outputs/walkthrough/scan_quick/scan_generic.root "chain table" \
  outputs/walkthrough/scan_quick_own_tables/scan_generic.root "own Bp"
```

The log starts with `[limit] confidence level = 90% one-sided asymptotic, q_threshold = 1.64237 (standard: --cl 0.900000)`; the same level is in the CSV header (`# cl_one_sided=0.9`) and in the top-left corner of the plot.

```
[limit] Forcing curve build from q(mchi,sigma) histogram.
[limit] ratio of each curve to input 1 ("reference"); median and 16-84% range per mass window, bigger = weaker limit
  "chain table" (g_dme_limit_1_q, 0 of 5 points on the sigma floor)
         0.7 -      3 MeV  n=   1 of    1   0.98 (0.98 - 0.98)
         3 -     10 MeV  n=   1 of    1   1.02 (1.02 - 1.02)
        10 -     30 MeV  n=   1 of    1   0.96 (0.96 - 0.96)
        30 -  1e+02 MeV  n=   1 of    1   0.94 (0.94 - 0.94)
     1e+02 -  1e+03 MeV  n=   1 of    1   0.93 (0.93 - 0.93)
  "own Bp" (g_dme_limit_2_q, 0 of 5 points on the sigma floor)
         0.7 -      3 MeV  n=   1 of    1   0.64 (0.64 - 0.64)
         3 -     10 MeV  n=   1 of    1   1.03 (1.03 - 1.03)
        10 -     30 MeV  n=   1 of    1   0.96 (0.96 - 0.96)
        30 -  1e+02 MeV  n=   1 of    1   0.94 (0.94 - 0.94)
     1e+02 -  1e+03 MeV  n=   1 of    1   0.93 (0.93 - 0.93)
```

(`n = 1 of 1` per window because the quick scans have five masses; with the full 800-mass scans the same table gives window medians. "n of N" counts the masses that have a crossing out of the masses of the window. The ratio is computed at the masses of the curve being compared, by interpolating the reference.)

**What the plots should look like** (`--out-pdf outputs/walkthrough/limit_compare.pdf`, labelled "90% C.L." in the top-left corner): the three curves are short segments between 2 and 200 MeV that lie on the published DAMIC-M line, which is the check that the rates, the signal table and the configs are right. Zoom with `--x-min 1.5 --x-max 300 --y-min 2e-39 --y-max 5e-36` to see them: at 2–4 MeV the own-Bp curve is below the other two, and above 20 MeV the chain-table curves are a few percent below the reference. The curves are straight between points because there are only five masses.

> **WATCH OUT.** When the data have an excess the q-map has no crossing at some masses (q at σ → 0 is above the threshold). The curve then drops to the bottom of the σ grid there, which looks like a vertical line on a plot; it means "no limit from this method at this mass", and the data prefer a signal. The table prints how many points of each curve are on the floor and leaves them out of the ratio (the `n of N` column). Use `--y-min` to keep the floor off the plot.

> **CROSS-CHECK.** `utils/compare_limit_curves.py` prints the same ratios from the ROOT files without the plot app (§14.4).

> **WATCH OUT.** With `--show-damic` or the literature overlays the program needs the files in `data/previous_limits/`; it finds them from the repository root, so the working directory no longer matters, but a file that is missing from the repository is reported as `ERROR: cannot open CSV file`. The DAMIC-M reference curve is first looked for in `collab_frameworks/` (not part of a plain clone); if it is not there the program prints a `note:` and uses the copy in `data/previous_limits/<mediator>_mediator/`.

---

## 11. Everything in one place

The commands of stages 0–6 and R, in order, as one copy-paste block. Nothing is skipped; this is the sequence a driver would run.

```bash
# 0. check a pattern by hand
build/ccdarksens_classify_pair 3.17 1.85 --bound 0.64

# 1. signal table (tree = pydme reference selection, chain = published selection)
#    expected: the tree check PASSES; the chain check prints FAIL (statistically different at n_e = 5; its effect on the limit is at most 1.2%, see Stage 1)
mkdir -p outputs/walkthrough/signal_tree outputs/walkthrough/signal_chain
build/ccdarksens_pattern_signal_eff configs/examples/pattern_signal_eff_tree_pydme.json 100000 outputs/walkthrough/signal_tree
build/ccdarksens_pattern_signal_eff configs/examples/pattern_signal_eff_chain.json      100000 outputs/walkthrough/signal_chain
build/ccdarksens_validate_pattern_efficiency data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv \
  --new outputs/walkthrough/signal_tree/Efficiencies_patterns.csv
build/ccdarksens_validate_pattern_efficiency data/efficiencies_paolo.csv \
  --new outputs/walkthrough/signal_chain/Efficiencies_patterns.csv --patterns 11 21 22 31 111 211

# 2. background matrices (use 1000000 trials for anything that goes into a limit)
mkdir -p outputs/walkthrough/background_tree outputs/walkthrough/background_chain
build/ccdarksens_pattern_background_eff configs/examples/pattern_background_eff.json       1000000 outputs/walkthrough/background_tree
build/ccdarksens_pattern_background_eff configs/examples/pattern_background_eff_chain.json 1000000 outputs/walkthrough/background_chain
build/ccdarksens_validate_background_efficiency data/Background_efficiencies.csv \
  --new outputs/walkthrough/background_tree/Background_efficiencies.csv --n-new 1000000   # tree reproduces the shipped matrix

# 3. Bp, from the matrix of the same selection as the signal table (chain)
build/ccdarksens_compute_brc --config configs/examples/background_pattern_rates.json \
  --matrix outputs/walkthrough/background_chain/Background_efficiencies.csv \
  --out outputs/walkthrough/background_Bp_chain.json

# R. rate tables for the five masses (605 files)
python3 utils/qedark_generate_grid.py configs/examples/qedark_generate_walkthrough.json
ls data/qedark_rates/Si/heavy/walkthrough_vE253p7 | wc -l      # 605

# 4-5. wire and scan: reference (shipped table), then your signal table, then your signal table and Bp
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_reference.json
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan.json
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_own_tables.json

# 5b. Monte Carlo threshold for the reference: toys at the step-1 limits, then the scan again
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_reference_toys.json
build/ccdarksens_scan_generic configs/examples/walkthrough_quick_scan_reference_mc.json

# 6. compare with the reference
build/ccdarksens_plot_limit --batch --heavy --ratio-to reference \
  --out-pdf outputs/walkthrough/limit_compare.pdf --out-root outputs/walkthrough/limit_compare.root \
  --out-csv outputs/walkthrough/limit_compare.csv \
  outputs/walkthrough/scan_quick_reference/scan_generic.root "reference" \
  outputs/walkthrough/scan_quick/scan_generic.root "chain table" \
  outputs/walkthrough/scan_quick_own_tables/scan_generic.root "own Bp"
```

The numbers each stage must give are the boxes marked CHECK above.

---

## 12. When it does not work

| Symptom | Cause | Fix |
|---|---|---|
| `failed to open pattern efficiency CSV: …/configs/data/…` | `efficiency_csv` path written from the repository root instead of from the folder above the config | write `../data/…` for a config in `configs/examples/` |
| every `upper_limit=1e-26` | no rate file was found: the masses or cross sections of the scan config are not the ones the generator wrote, or `model.rates_dir` is wrong | `ls` the folder; use the same mass values (six decimals) and σ axis in both configs |
| `ModuleNotFoundError: numpy` in stage R | Python without numpy | install numpy; it is the only Python dependency |
| the generator ran and nothing changed | `options.skip_existing` is true and the files exist | set `overwrite` to true, or use a new `rates_dir` |
| `upper_limit=1e-46` at many masses | q(σ→0) above threshold (excess in the data); no crossing | read it as "the data prefer a signal here"; in the `--from-qhist` plot it is a vertical line; the continuous limit (`--draw-both`, §14.4) gives a one-sided limit above the best fit |
| a changed matrix does not move the limit | the scan uses the pasted `background_Bp`, not the matrix | recompute Bp with `ccdarksens_compute_brc` and paste it |
| RMS pull ≫ 1 against a reference | the table was made with a different selection, α, pile-up or σ_ro | compare the `# Simulation Parameters` header; use the matching reference |
| validator says "give --n-ref and --n-new" | the CSV header has no `Nsims=` value | pass the trial counts explicitly |
| Bp for `{11}` jumps by 1–2 counts between matrices | rare entries such as B[11\|1] depend on few counts | use 10⁶ trials or more; compare the diagonal term, which is stable |
| `[bg-eff] WARNING: … differs from the validated value` | σ_ro, Qmin, Qmax or thresholds differ from the validated ones | intended for a study; otherwise restore the validated values |
| a misspelled config key does nothing | `ConfigManager` does not reject unknown keys | check the key against a shipped config |

---

## 13. Adding a new pattern bin

The recipe is in `Pipeline_Raw_Data_to_Limit.md` §15. In terms of the stages above:

1. Check that the selection can output the pattern (stage 0 with a pair or triple that has that shape; `chain_max_pair` is 3 by default, 4 adds `{41}`).
2. Add the new bin to the signal table: the generator already writes a row for every pattern it finds (stage 1).
3. Regenerate the matrix with the same selection and read the new row B[p_new | q] (stage 2).
4. Compute `B^rc_new` with `ccdarksens_compute_brc --config …` after adding the code to the `patterns` list of the rates config and setting `breakdown_pattern` to it (stage 3). The leading term is the product of the single-pixel populations of its pixels, so it depends on measured 2 e⁻, 3 e⁻ or 4 e⁻ pixel counts, not on the 1 e⁻ rate alone (S6).
5. Add the code to `experiment.pattern_roi`, a number to `run.background_Bp` and `run.background_Br`, and a count to `run.observed_counts` (stage 4).

> **NEEDED FROM THE COLLABORATION (C2).** `Br` for a new bin comes from the Geant4 simulation (τ_p), which is not in the repository. Until it is supplied the limit has to be given with `Br_new = 0` and with an upper bound, and both stated.

---

## 14. Independent cross-checks in Python

Four small Python scripts in `utils/` repeat one part of the pipeline with a separate implementation. They are not part of the workflow above; they exist so that a number from an app can be confirmed by code written separately. If an app and its script disagree, one of them has a bug, and finding which is the point.

| Script | Cross-checks | Stage | Needs |
|---|---|---|---|
| `utils/pattern_pvalues.py` | `ccdarksens_classify_pair` | 0 | standard library |
| `utils/compare_tables.py` | `ccdarksens_validate_*` with `--new` | 1, 2 | `pandas` |
| `utils/background_pattern_counts.py` | `ccdarksens_compute_brc --config` | 3 | standard library |
| `utils/compare_limit_curves.py` | `ccdarksens_plot_limit --draw-both --ratio-to` | 6 | `numpy`, `uproot` (`matplotlib` for `--plot`) |

Run them from the repository root after the §11 block has produced its outputs. In each case the numbers below are the ones that must agree with the stage.

### 14.1 `pattern_pvalues.py` against `ccdarksens_classify_pair`

```bash
python3 utils/pattern_pvalues.py 3.17 1.85 --bound 0.64
```

```
  p_11 = 0.00   p_21 = min(p_21,p_12) = 0.00   p_31 = min(p_31,p_13) = 0.17   p_22 = 1.66
  p_11 < 4.0: candidate 1+1
  p_21 < 4.0: promoted to 2+1
  p_31 < 4.0: promoted to 3+1

  result: {31}   (total charge 5.02 < 5.6)
```

Same p_mn table and the same `{31}` as Stage 0. Options are the same as the app: `--sigma` or `--bound`, `--left`, `--right`, `--qmax`. The script knows only the chain of the note; the app also runs the tree. A disagreement in a p_mn would mean the C++ classifier differs from the formulas of the note.

### 14.2 `compare_tables.py` against `ccdarksens_validate_*  --new`

```bash
python3 utils/compare_tables.py signal \
  --ref data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv \
  --new outputs/walkthrough/signal_tree/Efficiencies_patterns.csv \
  --n-ref 1000000 --n-new 100000

python3 utils/compare_tables.py background \
  --ref data/Background_efficiencies.csv \
  --new outputs/walkthrough/background_tree/Background_efficiencies.csv \
  --n-ref 100000 --n-new 1000000
```

```
signal: 80 entries, 49 non-zero in either table
  RMS pull (non-zero entries): 1.02   chi2/n: 1.03
  entries beyond 3 sigma: 0 (expected 0.1)   beyond 2 sigma: 3 (expected 2.2)   max |pull|: 2.60
...
background: 676 entries, 199 non-zero in either table
  RMS pull (non-zero entries): 0.98   chi2/n: 0.96
  entries beyond 3 sigma: 2 (expected 0.5)   beyond 2 sigma: 7 (expected 9.1)   max |pull|: 3.16
```

These are the lines of Stages 1 and 2, digit for digit. The script always needs `--n-ref` and `--n-new` (the app reads `Nsims=` from the header); `--patterns` restricts the signal comparison and `--out file.csv` writes the per-entry table. Several `--new` files are averaged.

### 14.3 `background_pattern_counts.py` against `ccdarksens_compute_brc`

```bash
python3 utils/background_pattern_counts.py configs/examples/background_pattern_rates.json \
  --lambda1 3e-4 --lambda2 2e-4 --lambda3 2e-4
```

```
  {11}        1.5694e+02       1.5138e+02       5.5631e+00
  {21}        2.3512e-02       1.0083e-02       1.3429e-02
  ...
  {"background_Bp": [156.942, 0.0235117, 0.0470225, 1.56828e-06, 1.00734e-06, 1.11374e-05]}
```

Same six numbers as `build/ccdarksens_compute_brc 3 2 2` in Stage 3. The script takes the same config file and the same overrides (`--lambda1`, `--lambda2`, `--lambda3`, `--n-pixels`); it has no `--matrix` flag, so the matrix is the one named in `migration_matrix_csv` of the config (copy the config and edit that line to test another). A table of segments goes in `pixels.segments_csv` as for the app.

### 14.4 `compare_limit_curves.py` against `ccdarksens_plot_limit --draw-both --ratio-to`

```bash
python3 utils/compare_limit_curves.py \
  outputs/walkthrough/scan_quick_reference/scan_generic.root \
  outputs/walkthrough/scan_quick/scan_generic.root \
  outputs/walkthrough/scan_quick_own_tables/scan_generic.root \
  --labels reference "chain table" "own Bp"
```

```
    window (MeV)     n           chain table                own Bp
    0.7-    3.0     1      1.02 (1.02-1.02)      0.65 (0.65-0.65)
    3.0-   10.0     1      1.00 (1.00-1.00)      0.96 (0.96-0.96)
   10.0-   30.0     1      0.94 (0.94-0.94)      0.97 (0.97-0.97)
   30.0-  100.0     1      0.94 (0.94-0.94)      0.94 (0.94-0.94)
  100.0- 1000.0     1      0.95 (0.95-0.95)      0.94 (0.94-0.94)
```

The script compares the **continuous bisection limit**, so the numbers match the app only with `--draw-both`. With `--from-qhist`, which is what the plots and Stage 6 use, the app compares the limit built from the q-map, which is what Stage 6 prints (0.98, 1.02, 0.96, 0.94, 0.93). At a single mass the two limits differ by up to 7% (here 1.03 against 0.96 at 6 MeV for the own-Bp curve); for the published counts the bisection limit over the stored one has window medians over the 800-mass scan of 1.03 (0.7–3 MeV) to 1.20 (3–10 and 100–1000 MeV). Add `--plot out.png` for a figure and `--sigma-floor` as in the app.

> **WHAT YOU LEARN.** A number that two separately written programs agree on is much more trustworthy than one a single program gives. When they differ, the difference is itself a finding (here: the bisection limit and the stored limit are not the same quantity at a single mass).

> **WATCH OUT.** The scripts are independent re-implementations, not tests that run automatically. If the apps are changed on purpose, the scripts must be changed with them.

---

## Appendix A — Key reference of the JSON files

Keys not listed here are not read by the program named.

### A.1 Generator configs (`ccdarksens_pattern_signal_eff`, `ccdarksens_pattern_background_eff`, validators with `--config`)

| Key | What it controls | Value in the signal configs / background configs |
|---|---|---|
| `detector.thickness_mm` | CCD thickness; the electron depth z is drawn uniformly in it | 0.669 (chain), 0.67 (background) |
| `response.efficiency_mc.A_um2` | diffusion constant A of σ_xy² = −A ln(1 − b z) (α + β E)² | 803.25 µm² |
| `response.efficiency_mc.b_umInv` | diffusion constant b | 6.5×10⁻⁴ µm⁻¹ |
| `response.efficiency_mc.alpha` | the α of (α + β E); the main diffusion knob | 1.0 (tree/pydme), 0.83 (chain) |
| `response.efficiency_mc.beta_per_keV` | energy dependence of the width; no effect at E = 0 | 0 |
| `response.efficiency_mc.rng_seed` | seed of the charge transport | 987654321 (signal), 123456789 (background) |
| `response.pattern_image.nrows_binned`, `ncols` | size of the binned image | 3 × 50 (signal), 3 × 5 (background cluster) |
| `response.pattern_image.row_binning`, `col_binning` | how many raw rows / columns make one binned pixel | 100 × 1 (signal), 1 × 1 (background) |
| `response.pattern_image.pixel_size_um` | pixel size | 15 µm |
| `response.pattern_image.sigma_readout_e` | Gaussian readout noise added to every pixel | 0.16 e⁻ |
| `response.pattern_image.lambda_dc` | Poisson dark-current pile-up per binned pixel; 0 switches it off | 4.1×10⁻⁴ (tree), 1.5×10⁻⁴ (chain), 0 (background) |
| `response.pattern_image.rng_seed` | seed of the image noise | as above |
| `response.pattern_classifier.Qmin_e` | seed threshold: a pixel below it starts no pattern | 0.6 e⁻ (= 3.75σ) |
| `response.pattern_classifier.neighbor_Qmax_e` | isolation threshold for neighbouring pixels and rows | 0.6 e⁻ |
| `response.pattern_classifier.Qmax_e` | maximum total charge of a pattern | 5.6 e⁻ |
| `response.pattern_classifier.sigma_res_e` | σ used in the CDFs of p_mn | 0.16 e⁻ |
| `response.pattern_classifier.thr_M`, `thr_MN`, `thr_MNL` | pass thresholds for one, two and three pixels | 3.5, 4.0, 5.5 |
| `response.pattern_classifier.max_e_per_pixel` | highest electron count tried per pixel | 5 |
| `response.pattern_classifier.enable_MN`, `enable_MNL` | allow two- and three-pixel patterns | true |
| `response.pattern_classifier.selection` | `tree` or `chain` | `tree` (default), `chain` |
| `response.pattern_classifier.chain_max_pair` | highest pair class of the chain: 3 gives up to 3+1, 4 adds 4+1 | 3 |
| `response.pattern_classifier.chain_allow_22` | label an un-extended 2+1 as 2+2 when p_22 < thr_MN | true |
| `response.pattern_classifier.chain_horizontal_isolation` | also require empty horizontal neighbours (the note does not) | false |
| `response.pattern_classifier.chain_first_only` | after the first multi-pixel pattern drop later ones (single pixels still counted) | true |

Present in the files but **not read** by the generators: all of `run`, `detector.rows`, `cols`, `pixel_size_um`, `target_element`, `density_g_cm3`, `mass_kg` (in the background configs `rows` and `cols` are used only if `pattern_image` gives no size), `response.mode`, `response.analysis_space`, and `response.efficiency_mc.sigma_readout_e`, `Qmin_e`, `Qmax_e`, `enable_MN`, `enable_MNL` (the classifier block is the one read).

### A.2 Rates config (`ccdarksens_compute_brc --config`)

| Key | Meaning |
|---|---|
| `rates.lambda_1`, `lambda_2`, `lambda_3plus` | mean number of pixels with 1, 2, 3-or-more electrons per pixel readout (θ × 10⁻⁴) |
| `pixels.n_pixel_readouts` | pixel readouts N_pix; `null` = sum of `Nusedpix` of `image_table_csv` (1,853,807,441) |
| `pixels.image_table_csv` | the per-image table with `Nusedpix` |
| `pixels.segments_csv` | optional table of segments with their own rates; replaces the uniform rates |
| `migration_matrix_csv` | the matrix B[p\|q] to use |
| `patterns` | the pattern codes to report, in the order of the output (default 11, 21, 111, 31, 22, 211) |
| `breakdown_pattern` | the pattern whose injected-pattern contributions are listed |
| `output_json` | where to write the `background_Bp` list |

### A.3 Scan config (`ccdarksens_scan_generic`)

| Key | Meaning | Value here |
|---|---|---|
| `run.label`, `run.outdir` | name and output folder of the run | set per run |
| `run.cl` | confidence level; sets the q threshold (0.90 gives 1.642) | 0.9 |
| `run.one_sided_q` | q = 0 at and below the best-fit cross section on the grid; the limit is the first crossing above it. Default false (q is made non-decreasing from the lowest σ) | true |
| `run.robust_minimum` | the reference NLL minimum is min(2D fit, grid minimum, background-only). Default false | true |
| `run.mode` | `scan` (default) or `threshold_toys` (toy calibration of the threshold, then exit) | `scan`; `threshold_toys` in the `_toys` configs |
| `run.threshold_toys.sigma_threshold_graph_path` | ROOT file with the cross section per mass to calibrate at: a `ccdarksens_band` output or the `scan_generic.root` of step 1 | step-1 output |
| `run.threshold_toys.n_threshold_toys`, `percentile`, `rng_seed` | toys per mass, the quantile of q that becomes the threshold (0.90), the seed | 10000, 0.90, 23456 |
| `run.threshold_toys.one_sided` | true: a toy whose best fit lies above the tested σ contributes q = 0 (the published curve); false: it contributes t > 0 (two-sided) | true |
| `run.q_target_lookup_path` | `qtarget_threshold.root` of the toy step; the scan then uses that per-mass threshold instead of 1.642 | set in the `_mc` configs |
| `run.background_source` | `bp_br_template` uses the pasted Bp and Br; `dc_flat_migration` computes the background from the dark current | `bp_br_template` |
| `run.background_model` | `Bp_theta_Br`: expected = S + Bp + θ·Br | `Bp_theta_Br` |
| `run.background_Bp`, `background_Br` | six numbers each, in the order of `experiment.pattern_roi` | published, or yours |
| `run.observed_counts` | the observed counts per bin; takes precedence over `run.data_path` and over the Asimov data | 144, 0, 0, 1, 0, 0 |
| `run.constrain_*`, `theta_lo`, `theta_hi` | the constraint on θ from the radiogenic control sample (98 clusters) and its bounds | tau-weighted, 0.5 to 10 |
| `run.profile_minimizer`, `pydme_style_ul` | the 2D `pydme` fit and the continuous bisection limit | `pydme`, true |
| `detector.mass_kg`, `experiment.livetime_days`, `duty_cycle` | exposure = mass × livetime × duty / 365.25 | 0.01523 kg, 85.356 d, 1 → 0.003559 kg·yr |
| `experiment.pattern_roi`, `observable_bins` | the pattern bins and the observable space | 11, 21, 111, 31, 22, 211; `pattern` |
| `response.analysis_space` | selects the response backend | `pattern` |
| `response.efficiency_mc.efficiency_csv` | the signal table, path relative to the folder above the config | pydme table, or yours |
| `model.type`, `rates_dir`, `filename_template` | the signal model and where its rate files are | `dm_electron`; `data/qedark_rates/Si/heavy/long_scan_vE253p7` (published) or `…/walkthrough_vE253p7` (quick scans) |
| `model.grid.mchi_MeV`, `model.grid.sigma_e_cm2` | the mass axis and the cross-section axis; they must equal the generator's axes | published: 800 logspace 0.2–1000 MeV, 300 logspace 10⁻⁴⁶–10⁻²⁶ cm²; quick scans: five masses, 121 logspace 10⁻⁴²–10⁻³⁴ cm² |

### A.4 Rate generator config (`utils/qedark_generate_grid.py`)

| Key | Meaning | Value here |
|---|---|---|
| `material`, `mediator` | the target and the mediator class: `heavy` (also `massive`) or `light` (also `massless`) | `Si`, `heavy` |
| `halo.v0_kms`, `vE_kms`, `vesc_kms` | standard halo model: circular speed, Earth speed, escape speed | 238, 253.7, 544 km/s |
| `detector.band_gap_eV` | threshold below which no electron is produced | 1.2 eV |
| `detector.eh_pair_eV` | energy per electron–hole pair above the gap | 3.8 eV |
| `detector.binsize_eV` | energy bin of the output spectrum | 0.1 eV |
| `rates_dir` | folder the files are written to; the scan config's `model.rates_dir` must be this | `data/qedark_rates/Si/heavy/walkthrough_vE253p7` |
| `filename_template` | file name pattern; the scan builds the same name | `dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv` |
| `grid.mchi_MeV` | mass axis: `values`, `linspace` or `logspace` | five values |
| `grid.sigma_e_cm2` | cross-section axis | `logspace` 10⁻⁴²–10⁻³⁴, 121 points |
| `grid.format` | number formats used in the file names | mass `.6f`, σ `.1e` |
| `options.parallel` | number of worker processes | 4 |
| `options.skip_existing`, `overwrite`, `compress` | do not rewrite existing files; replace them; write `.gz` | true, false, false |

---

## Appendix B — Files

| File | Role |
|---|---|
| `apps/ccdarksens_classify_pair.cc` | stage 0 |
| `apps/ccdarksens_pattern_signal_eff.cc`, `apps/ccdarksens_pattern_background_eff.cc` | stages 1 and 2 (thin wrappers) |
| `src/response/PatternTables.cc` | the generation, reading, writing and comparison shared by the generators and the validators |
| `apps/ccdarksens_validate_pattern_efficiency.cc`, `apps/ccdarksens_validate_background_efficiency.cc` | the checks |
| `apps/ccdarksens_compute_brc.cc` | stage 3 |
| `utils/qedark_generate_grid.py`, `python/ccdarkphys/qedark/` | stage R: the grid driver and the QEDark rate code |
| `apps/ccdarksens_scan_generic.cc` | stages 4–5 |
| `apps/ccdarksens_plot_limit.cc` | stage 6 |
| `utils/pattern_pvalues.py`, `utils/compare_tables.py`, `utils/background_pattern_counts.py`, `utils/compare_limit_curves.py` | independent Python cross-checks (§14) |
| `configs/examples/pattern_signal_eff_{tree_pydme,chain}.json` | stage 1 inputs |
| `configs/examples/pattern_background_eff{,_chain}.json` | stage 2 inputs |
| `configs/examples/background_pattern_rates.json` | stage 3 input |
| `configs/examples/qedark_generate_walkthrough.json` | stage R input |
| `configs/examples/dm_electron_heavy_lbc_1p3kgday{,_toys,_mc}.json` | the published scan in its three threshold steps (§9.1) |
| `configs/examples/walkthrough_quick_scan_reference_{toys,mc}.json` | steps 2 and 3 for the five-mass reference |
| `configs/examples/walkthrough_quick_scan_reference.json`, `walkthrough_quick_scan.json`, `walkthrough_quick_scan_own_tables.json` | stages 4–5: shipped table / your signal table / your table and Bp |
| `data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv` | pydme signal reference |
| `data/efficiencies_paolo.csv` | published signal efficiencies |
| `data/Background_efficiencies.csv` | shipped background matrix |
| `data/Final_Combined_Image_Data.csv` | per-image used pixels and candidate counts |
