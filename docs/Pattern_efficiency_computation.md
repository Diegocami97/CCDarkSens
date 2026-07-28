# How pattern efficiencies are computed — inputs and outputs

Pattern efficiencies are **P(pattern | n_e)**: the probability that an event with **n_e** true electrons is classified as a given **pattern**. They are computed by Monte Carlo in **PatternMC::BuildPatternTable**, then used to fold spectra from n_e to pattern space.

**In this C++ code:** We use a **single** map **ε(pattern | n_e)** for both signal and background. Only the **n_e range** and the **input spectrum** differ (background: B_tot(n_e), n_e = 0…ne_max; signal: S_obs(n_e), n_e = ne_min…ne_max). The notebook instead defines **different** efficiencies for signal vs background; see §10 below.

---

## 1. Where it happens

- **PatternMC::BuildPatternTable(ne_min, ne_max, Ee_ref_eV)**  
  Fills the internal table **pattern_table_[n_e][label] = P(label | n_e)**.
- The app then copies this into **pattern_eff_map**: `(pattern_id, n_e) → P(pattern | n_e)` (pattern_id = decimal encoding of label, e.g. 111 for (1,1,1)).

---

## 2. Inputs to the efficiency computation

### 2a. Parameters (config / constructor)

| Input | Meaning | Source |
|-------|--------|--------|
| **ne_min, ne_max** | Range of true electron number | Experiment binning (e.g. 0–20) |
| **Ee_ref_eV** | Reference energy (eV) for diffusion | Fixed (e.g. 50 eV) |
| **ne_trials** | MC trials per n_e | `pattern_mc.n_events_per_ne` (e.g. 20000) |

### 2b. Physics / detector (2D path)

| Input | Meaning | Source |
|-------|--------|--------|
| **ChargeTransport** | Depth z, σ_xy(z,E), cloud sampling | Detector thickness, pattern_mc (A_um2, b_umInv, α, β) |
| **PatternImageGenerator** | 2D binned image (n_e, Ee_eV) → image | pattern_image (nrows_binned, ncols, row/col_binning, pixel_size_um, sigma_readout_e, lambda_dc) + detector raw size if used |
| **PatternClassifier** | Hit threshold, M/MN/MNL thresholds | pattern_classifier (Qmin_e, thr_M, thr_MN, thr_MNL, enable_MN, enable_MNL, max_e_per_pixel) |

### 2c. Physics / detector (1D row path, when no 2D generator)

| Input | Meaning | Source |
|-------|--------|--------|
| **ChargeTransport** | As above | Same |
| **PixelSimulator** | 1D row of pixel charges | pattern_mc (row_length, pixel_size_um, sigma_readout_e, lambda_dc) |
| **PatternClassifier** | As above | Same |

### 2d. Optional external input

| Input | Meaning | Source |
|-------|--------|--------|
| **efficiency_csv** | Precomputed (pattern_id, n_e, efficiency) | If set, **replaces** MC-built table for filling **ε(n_e)** and for **pattern_eff_map** when loading from file instead of PatternMC. |

---

## 3. Algorithm (how P(pattern | n_e) is computed)

### 3a. 2D image path (when `use_2d_image_efficiency` and PatternImageGenerator is set)

For each **n_e** in [ne_min, ne_max]:

1. For **trial = 1 … ne_trials**:
   - **Generate one 2D image**: `image_2d = img_gen->GenerateImage(n_e, Ee_ref_eV)`.
     - Sample depth z, get σ_xy(z, E); sample n_e cloud positions; fill raw pixel grid; apply row/column binning; add readout noise (and optional dark current).
   - Take **rows 0, 1, 2** as row_above, row_middle, row_below.
   - **Classify**: `best = classifier->ScanImage2DWithIsolation(row_above, row_middle, row_below)`.
     - Scan 5-pixel windows on the middle row; require **isolation** (see below); for each valid window try M, MN, MNL and pick best (e.g. by total charge or statistic).
   - If `best.valid`, increment **counts[best.label]**.
2. **P(label | n_e) = counts[label] / n_trials** for each label that appeared.
3. Store **pattern_table_[n_e] = { label → P(label|n_e) }**.

**What “isolated” means:** For the 2D image we use three rows (above, middle, below). A pattern on the **middle row** is counted only if it is **isolated**: in the same column range (the 5-pixel window), the **row above** and **row below** must have charge **below** the isolation threshold (e.g. `Qmin_e`). So “isolated” = no significant charge in the pixels directly above or below the pattern. That rejects events where a second cluster or dark current in a neighboring row could contaminate the pattern, and matches the notebook definition (only patterns with empty neighbors in the perpendicular direction).

### 3b. 1D row path (when no 2D generator)

For each **n_e** in [ne_min, ne_max]:

1. For **trial = 1 … ne_trials**:
   - Sample z, σ_xy(z, E); sample n_e positions; deposit on 1D row (PixelSimulator); add readout noise (and optional DC).
   - **Classify**: `ScanRow(row)` → pick best pattern (e.g. by total charge).
   - If best valid, increment **counts[best.label]**.
2. **P(label | n_e) = counts[label] / n_trials**.
3. Store **pattern_table_[n_e]**.

---

## 4. Outputs

### 4a. From BuildPatternTable

- **pattern_table_** (internal to PatternMC):  
  `pattern_table_[n_e]` = map **PatternLabel → double**  
  So for each n_e you get a probability distribution over **labels** (e.g. (1), (1,1,1)) that sum to ≤ 1 (can be < 1 if some events are invalid/unclassified).

### 4b. pattern_eff_map (used by the apps)

- Apps build **pattern_eff_map: (pattern_id, n_e) → efficiency**:
  - **pattern_id** = decimal encoding of label (e.g. 1 for (1), 11 for (1,1), 111 for (1,1,1)).
  - **efficiency** = P(pattern | n_e) from the table (or from CSV if provided).
- So the **output** of the efficiency computation, as used later, is:
  - **pattern_eff_map**: for each (pattern_id, n_e) in the table (or CSV), the value **ε(pattern_id, n_e) = P(event classified as that pattern | n_e)**.

### 4c. Optional file output

- **efficiency_per_pattern.csv** (example app): columns `pattern_id, n_e, efficiency` for each pattern in **pattern_roi** and each n_e (efficiency = 0 if not in table).

---

## 5. How the efficiencies are used (pattern folding)

- **FoldNeToPatternRates(h_ne, ne_min, ne_max, pattern_roi, pattern_eff_map)**:
  - For each **pattern_id** in **pattern_roi**:
    - **rate(pattern_id) = Σ_{n_e} h_ne(n_e) × ε(pattern_id, n_e)**.
  - So **ε(pattern_id, n_e)** is the “efficiency” to go from n_e to that pattern bin.
- Applied **separately** to:
  - **B_tot(n_e)** → **B_pat** (background),
  - **S_obs(n_e)** → **S_pat** (signal).

---

## 6. Summary table

| Stage | Inputs | Outputs |
|-------|--------|--------|
| **BuildPatternTable** | n_e range, Ee_ref_eV, ne_trials; ChargeTransport; PatternImageGenerator (2D) or PixelSimulator (1D); PatternClassifier | **pattern_table_[n_e][label] = P(label \| n_e)** |
| **Fill pattern_eff_map** | pattern_table_ (or efficiency_csv); pattern_roi (for CSV writing) | **pattern_eff_map: (pattern_id, n_e) → ε** |
| **FoldNeToPatternRates** | h_ne (e.g. B_tot or S_obs), ne_min, ne_max, pattern_roi, pattern_eff_map | **Vector of rates per pattern bin** (one per pattern_roi entry) |

So: **efficiencies per pattern** are **ε(p, n_e) = P(pattern p | n_e)**, with inputs = (n_e, E_ref, detector/response config, classifier config, number of trials), and outputs = table **pattern_table_** → **pattern_eff_map** used to fold **B_tot** and **S_obs** into **B_pat** and **S_pat**.

---

## 7. Step-by-step: how P(pattern | n_e) is computed

1. **For each n_e** (e.g. 0, 1, 2, …, 20) we run **ne_trials** independent MC events (e.g. 20k or 50k).

2. **Each trial:**
   - **2D path:** Generate one 2D binned image: sample depth z, diffusion σ_xy(z,E), place n_e electrons in the cloud, fill raw pixels, apply row/column binning, add readout noise (and optional dark current).  
   - **1D path:** Generate one 1D row: same physics, but deposit onto a single row (PixelSimulator).

3. **Classify that event:**
   - **2D:** Use only **rows 0, 1, 2** as (above, middle, below). Call **ScanImage2DWithIsolation**: scan 5-pixel windows on the **middle row (row 1)**; for each window require **isolation** (rows 0 and 2 below threshold in that column range); if a valid pattern is found, take the best one (e.g. by total charge). If **best.valid**, we count **counts[best.label]++**.  
   - **1D:** Call **ScanRow(row)**, pick best pattern; if found, count it.

4. **No filtering by pattern_roi in the MC:** We count **every** label the classifier returns. The table stores **all** labels that appeared (e.g. (1), (2), (1,1), (2,1), …). So the MC does not know about “pattern_roi”; it just records the empirical distribution over **PatternLabel**.

5. **Normalize:** P(label | n_e) = counts[label] / n_trials for each label that appeared. Labels that never appeared have P = 0 (implicit).

6. **pattern_eff_map:** When building the map from the table, we convert each **PatternLabel** to an integer **pattern_id** (e.g. (1,1) → 11, (2,1,1) → 211) and store (pattern_id, n_e) → P. When we **fold** or **report** for **pattern_roi**, we look up (pattern_id, n_e) in this map; if the MC never produced that pattern for that n_e, the efficiency is **0**.

---

## 8. Why can some (or many) patterns have 0 efficiency?

- **We only count when the classifier returns a valid pattern.** If in the 2D path we use **only the middle row (row 1)** and require **isolation**, then:
  - If most of the charge lands in **row 0 or row 2** (e.g. diffusion + row binning), the middle row may have no above-threshold seed → **no pattern found** → no count for any label. So the **total** efficiency (sum of P over all patterns) can be **&lt; 1**.
  - If the middle row has charge but **row above or below** often has charge above the isolation threshold, we **reject** the event (not isolated) → again no count.
  - So many trials can yield **best.valid = false**. Among the trials that do pass, only **certain** labels may appear (e.g. single-pixel patterns (1), (2), (3) on row 1 with empty neighbors). Multi-pixel patterns (11, 21, 31, 22, 211, 111) may **rarely or never** appear on that single row with isolation satisfied → **0 efficiency** for those pattern_roi entries.

- **So “0 efficiency” does not mean the detector is wrong.** It means: in our MC, for that n_e, we **never** (or almost never) got that pattern when we only look at the **middle row** and require **isolation**. That can happen if:
  - The 2D image has only 3 binned rows and we fix “middle” = row 1, so we are effectively restricting to events that land in that one row with empty neighbors.
  - Diffusion and binning make that a small or skewed subset of events, so the empirical P(pattern|n_e) is 0 for many pattern IDs.
  - **Previously:** C++ used the valid pattern with **maximum total charge** (the window on the cluster peak), while the notebook uses the **first** valid pattern in scan order. That made C++ almost always pick single-pixel patterns and gave 0 for 11, 21, 31, 22, 211, 111. **Now** C++ uses the first valid pattern, matching the notebook, so all those patterns get non-zero efficiencies when the MC produces them.

- **If you see 0 for *all* patterns** for some n_e: then almost no trials passed the 2D isolation + middle-row pattern finding (or the 1D row rarely had a valid pattern). That can be fixed by: increasing trials, relaxing isolation, or (in the 2D case) revisiting row assignment / image size so the “middle” row more often contains the cluster.

---

## 9. Parity with the notebook (How_to_calculate_efficiencies_for_patterns.ipynb)

Yes. The C++ 2D path matches the notebook’s efficiency calculation:

| Step | Notebook (`efficiencies.py`) | C++ (PatternMC + PatternImageGenerator + PatternClassifier) |
|------|-------------------------------|-------------------------------------------------------------|
| Image | `generate_image_E(ne, Ee)` → (3, 50) binned; adds **dark current** when `poison=True` (Poisson(λ) per pixel) | `PatternImageGenerator::GenerateImage(n_e, Ee_eV)` → same; adds dark current when `pattern_image.lambda_dc > 0` (`include_dark_current`) |
| Rows | Uses **row 1** as middle; rows 0, 2 as up/down (`image[0]`, `image[1]`, `image[2]`) | Uses **row 1** as middle; rows 0, 2 as above/below in `ScanImage2DWithIsolation` |
| Isolation | Require `image_ext_up[...] < qmin` and `image_ext_down[...] < qmin` in the same column range | Require `row_above[k] < thr` and `row_below[k] < thr` for the 5-pixel window (thr = Qmin_e) |
| Pattern ID | 5-pixel window → `identify_pattern(segment)` (Pm, Pmn, Pmnl) | 5-pixel window → `classify_from_seed` (M, MN, MNL, same statistics) |
| Which pattern per event | **First** valid isolated pattern in scan order (left to right) | **First** valid isolated pattern in scan order (same as notebook) |
| Efficiency | `pattern_simulation(Nsims, ne, Ee)` → count patterns, **efficiency = count / Nsims** | `BuildPatternTable` → for each n_e run `ne_trials`, count labels → **P(label\|n_e) = count / n_trials** |

So the notebook and C++ both: generate a 2D binned image per event, take the middle row (row 1), scan 5-pixel windows left-to-right, require isolation (above and below below threshold), and count the **first** valid pattern found (not the one with maximum charge). Taking “first” matches the notebook and yields non-zero efficiencies for multi-pixel patterns (11, 21, 31, 22, 211, 111); taking “max charge” would bias toward single-pixel patterns at the cluster peak and give zeros for many pattern_roi. For **signal** the procedure is the same in both; the notebook additionally defines **different** efficiencies for **background** (see §10).

---

## 10. Notebook: different efficiencies for signal vs background

The notebook defines **two different** efficiency frameworks:

### Signal (same as C++)

- **pattern_simulation(ne, Ee)**: generate **n_e electrons**, diffuse, bin to (3,50), add readout noise and **dark current** (when `poison=True`: Poisson(λ) per pixel), then **scan_image** → one pattern per event. So the signal path includes dark current; e.g. (1,1) from 1 electron gets an important contribution from DC landing next to the signal.
- C++ does the same: when **pattern_image.lambda_dc > 0**, `PatternImageGenerator` adds Poisson(lambda_dc) per pixel (`include_dark_current`), so the built ε(pattern | n_e) includes DC in the signal path.
- Output: **ε(pattern | n_e)** = P(identified pattern | n_e electrons). Same as our **BuildPatternTable** / **pattern_eff_map**.

### Background (implemented in C++)

- **pattern_simulation_background()** (notebook): for each **true pattern** (0), (1), (2), (1,1), (2,1), … (all combinations up to 5 electrons), call **simulate_cluster(m, n, l)**.
  - **simulate_cluster** builds a **3×5 ideal cluster**: middle row has exactly m, n, l electrons in three pixels, **no diffusion**. Only **readout noise** is added.
  - Then **scan_image_background** (and **identify_pattern_background**, which also allows pattern (0)) → identified pattern.
- Output: **P(identified pattern | true pattern)** = a **migration matrix** between patterns (including (0) for empty). So background uses "true pattern → readout only → identified pattern".
- In the formula **B_mnl = P(m)·P(n)·P(l)·ε_mnl + M**: P(m)P(n)P(l) are Poisson(dark current) probabilities for true pattern (m,n,l), ε_mnl is the probability of **correctly** identifying it, and **M** is the rate from **misidentification** (other true patterns that get identified as mnl).

| | Signal | Background (notebook) |
|---|--------|------------------------|
| **Input** | n_e (number of electrons) | True pattern (m,n,l) or (0) |
| **Physics** | Diffusion + readout (+ optional DC) | **Readout only** (no diffusion) |
| **Output** | ε(pattern \| n_e) | P(identified pattern \| true pattern) = migration matrix |

### What C++ does for background

- **Same map ε(pattern | n_e) for B_tot:** The scan/example apps fold B_tot(n_e) with ε(pattern | n_e) as before (no separate background efficiency table).
- **Background efficiency matrix (notebook parity):** The app **ccdarksens_pattern_background_eff** implements **pattern_simulation_background**:
  - Enumerates ideal patterns: **(0)** first, then (1)..(5), (1,1), (1,2), …, (1,1,1), … with sum ≤ 5 (same as notebook `all_pats`).
  - For each ideal, runs **Nsims** trials: **SimulateCluster(b, c, d)** → 3×5 image (readout noise only), then **ScanImage2DWithIsolation** with a classifier that has **allow_pattern_zero = true** (notebook **identify_pattern_background**: single-pixel branch can return pattern (0) when charge is below 1-e threshold).
  - Counts identified pattern per trial; when no valid pattern is found, counts as **(0)**.
  - Writes **Background_efficiencies.csv**: columns **iden_pat**, **eff_0**, **eff_1**, **eff_11**, …; each row = one ideal pattern, each cell = P(identified = column pattern | ideal = row pattern). This file can be used by analysis code (e.g. **Background_pattern_** in the repo) to build the background model with misidentification.
- **Config:** **pattern_classifier.allow_pattern_zero** (default false) is only set true when building the classifier inside the background app; signal/scan apps keep it false.
