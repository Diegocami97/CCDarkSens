# Background efficiency: CCDarkSens vs pydme

This note explains how background is applied in the two frameworks and why results can differ.

## 1. Pydme (NLL exclusion)

- **File:** `collab_frameworks/pydme/pydme/dme_nLL_exclusion.py` uses `Background_pattern(theta, pattern, Npix)` from `background_pattern.py`.
- **Background_pattern (used in NLL):** Does **not** read `Background_efficiencies.csv`. It uses a **hardcoded** parametric model:
  - **Patterns:** 6 only: `[11, 21, 111, 31, 22, 211]`.
  - **Formula:** `B(pattern) = Bp[pattern] + theta[0]*Br[pattern]`, with `Bp` and `Br` fixed arrays (and divided by number of time bins `len(gamma)`).
  - **Typical values in code:**  
    `Bp = [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]`,  
    `Br = [0.039, 0.039, 0.016, 0.052, 0.011, 0.035]`.
- So in pydme the “background efficiency” in the exclusion NLL is **not** a migration matrix; it is this **Bp + θ·Br** template. No CSV is used there.

The other function **Background_pattern_** in the same file **does** read the CSV and builds B from a Poisson DC model and the migration matrix; it is **not** used by the main exclusion code (`dme_nLL_exclusion` calls `Background_pattern`, not `Background_pattern_`).

## 2. CCDarkSens (pattern scan / example)

- When **`background_efficiency_csv`** is **set**:
  - **B_tot(n_e)** is built from DC + flat (rates in n_e).
  - **B_true_pat:** only single-pixel true patterns are filled: `B_true_pat[0..5] = B_tot(0..5)` (n_e → pattern codes 0–5); multi-pixel true rates are 0.
  - **Migration:** `Background_efficiencies.csv` is read as **row = true pattern, column = identified pattern**, i.e. `value = P(identified | true)`.
  - **B_pat** is computed as:  
    `B_pat[iden] = Σ_{true} B_true_pat[true] × P(iden | true)`  
    (`FoldBackgroundWithMigration`).
- When **`background_efficiency_csv`** is **empty**:
  - **B_pat** is obtained by folding **B_tot(n_e)** with the **signal** efficiency **ε(pattern | n_e)** (same as signal), no migration matrix.

So in CCDarkSens, “background efficiency” in the **data configs** means: **apply the migration matrix from the CSV** (true pattern → identified pattern). In pydme’s exclusion, there is **no** such step; they use Bp/Br directly.

## 3. Main differences

| Aspect | Pydme (exclusion NLL) | CCDarkSens (with `background_efficiency_csv`) |
|--------|------------------------|-----------------------------------------------|
| **Background model in NLL** | B = **Bp + θ·Br** (hardcoded Bp, Br) | B from **DC + flat** → **B_tot(n_e)** → **B_true_pat(0..5)** → **migration CSV** → **B_pat** |
| **Use of Background_efficiencies.csv** | **Not used** in main NLL | **Used** as migration P(identified \| true) |
| **Pattern set** | 6 patterns: 11, 21, 111, 31, 22, 211 | Configurable **pattern_roi** (e.g. 11 patterns) |
| **B source** | Fixed Bp, Br per pattern | Physics: DC rate, exposure, then migration |

So:

- **Pydme:** B is a **parametric template** (Bp, Br); no migration applied in the NLL.
- **CCDarkSens:** B is **physics-based** (DC + flat in n_e, then migration matrix to pattern space).

Numerical agreement is only expected if:

- You compare using the **same** pattern set (e.g. restrict to pydme’s 6 patterns), and  
- Either:
  - You feed pydme with a template that was produced by the same DC + migration as in CCDarkSens, or  
  - You turn off migration in CCDarkSens and use a Bp/Br-style template (if we add that option).

## 4. CSV convention (for migration)

- **CCDarkSens writer** (`ccdarksens_pattern_background_eff`):  
  **Row = true (ideal) pattern**, **column = identified pattern**  
  → `value = P(identified | true)`.
- **CCDarkSens reader** (scan/example):  
  Reads `(first_col, eff_X)` as `(true_pat, iden_pat)` and stores  
  `migration_map[(true_pat, iden_pat)] = P(iden | true)`.  
  So the formula `B_pat[iden] = Σ_true B_true[true] × P(iden|true)` is correct; **no transpose**.
- **Pydme Background_pattern_** (not used in exclusion):  
  Uses `row = np.where(iden_pat == pattern)[0][0]` and then columns; it effectively expects a different indexing (row = identified). So if the **same** CSV file (row = true, col = iden) were used there, that function would be using the matrix in a different way. The main pydme exclusion does not use that function.

## 5. What to do if results don’t match

1. **Pattern set**  
   Use the same patterns in both. Pydme uses 6: `[11, 21, 111, 31, 22, 211]`. If you use more or different patterns in CCDarkSens, B_pat and the likelihood will differ.

2. **Source of B**  
   - In pydme, B is **Bp + θ·Br** (no migration in NLL).  
   - In CCDarkSens, B is **DC + flat → migration**.  
   So B values will match only if pydme’s Bp/Br were generated from the same DC and migration (or equivalent).

3. **To align with pydme’s B model**  
   - Set **pattern_roi** to pydme’s 6 patterns.  
   - In config, set **background_Bp** and **background_Br** to pydme’s arrays (same order as pattern_roi).  
   - Use **background_model: "Bp_theta_Br"**.  
   - Then the **profile likelihood** in CCDarkSens uses the same parametric form as pydme.  
   - To avoid mixing two B models: either **clear `background_efficiency_csv`** and feed B only via Bp/Br (would require adding a “B from template only” path), or accept that our **central** B is from DC+migration and only the **profile** (Bp + θ·Br) matches pydme’s form.

4. **To keep CCDarkSens physics-based**  
   Keep **background_efficiency_csv** and DC/flat as now. Then B is uniquely defined by DC, exposure, and migration. Agreement with pydme would require using the same migration matrix and the same DC/pattern assumptions when generating pydme’s Bp/Br elsewhere.

## 6. Summary

- **Pydme exclusion NLL:** B = Bp + θ·Br (hardcoded); **Background_efficiencies.csv is not used**.
- **CCDarkSens (with background_efficiency_csv):** B = DC + flat → B_tot(n_e) → B_true_pat(0..5) → **migration from CSV** → B_pat.
- Differences in results are expected unless pattern set and B source (and migration, if any) are aligned. Use this doc to decide whether to match pydme’s template (Bp/Br, 6 patterns) or to keep the current DC+migration pipeline and only compare when B is generated consistently.
