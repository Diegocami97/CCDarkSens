# Efficiency CSV comparison: three files

## Files

| File | Path | Parameters |
|------|------|------------|
| **Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv** | `data/` or pydme `.../pattern_efficiency/` | Nsims=1000000, alpha=1, beta=0, A=803.25, b=0.00065, **lambda=0.00041**, DC=True |
| **efficiency_per_pattern.csv** | `outputs/scan_pattern/efficiency_per_pattern.csv` | CCDarkSens output (no param header) |
| **efficiencies_paolo.csv** | `data/efficiencies_paolo.csv` | Nsims=100000, alpha=1, beta=0, A=803.25, b=0.00065, **lambda=0.00015**, DC=True |

Columns: **pattern, ne, Efficiency**. ROI patterns: **[11, 21, 111, 31, 22, 211]**.

---

## Three-way comparison (all ROI patterns, ne = 1..5)

| pattern | ne | Efficiencies_Nsims1000000_DCTrue_alpha1 | efficiency_per_pattern (CCDarkSens) | efficiencies_paolo |
|---------|----|----------------------------------------|-------------------------------------|-------------------|
| 11 | 1 | 0.00073 | 7.20e-04 | 0 |
| 11 | 2 | 0.37418 | 3.93e-01 | 0.38327 |
| 11 | 3 | 0.13724 | 4.89e-01 | 0.10816 |
| 11 | 4 | 0.05675 | 4.62e-01 | 0.03285 |
| 11 | 5 | 0.03035 | 4.09e-01 | 0.01526 |
| 21 | 1 | 0 | 0 | 0 |
| 21 | 2 | 0.00057 | 5.70e-04 | 0.00023 |
| 21 | 3 | 0.36556 | 3.66e-01 | 0.41192 |
| 21 | 4 | 0.12272 | 1.23e-01 | 0.10566 |
| 21 | 5 | 0.04917 | 4.90e-02 | 0.03008 |
| 111 | 1 | 0 | 0 | 0 |
| 111 | 2 | 0.00036 | 3.50e-04 | 0.00016 |
| 111 | 3 | 0.13529 | 1.35e-01 | 0.12945 |
| 111 | 4 | 0.03363 | 3.40e-02 | 0.06120 |
| 111 | 5 | 0.009 | 9.00e-03 | 0.02770 |
| 31 | 1 | 0 | 0 | 0 |
| 31 | 2 | 0 | 0 | 0 |
| 31 | 3 | 0.00023 | 2.30e-04 | 0.00012 |
| 31 | 4 | 0.19667 | 1.97e-01 | 0.23724 |
| 31 | 5 | 0.06525 | 6.50e-02 | 0.48930 |
| 22 | 1 | 0 | 0 | 0 |
| 22 | 2 | 0 | 0 | 1e-07 |
| 22 | 3 | 0.00016 | 1.60e-04 | 0.00007 |
| 22 | 4 | 0.11409 | 0 | 0.13332 |
| 22 | 5 | 0.02996 | 0 | 0.02361 |
| 211 | 1 | 0 | 0 | 0 |
| 211 | 2 | 0 | 0 | 0 |
| 211 | 3 | 0.00053 | 5.30e-04 | 0.00022 |
| 211 | 4 | 0.22307 | 2.23e-01 | 0.21692 |
| 211 | 5 | 0.06319 | 6.30e-02 | 0.27704 |

---

## Summary

1. **Pattern 11**  
   CCDarkSens is **higher** at ne=2,3,4,5 than both the 1M reference and Paolo (e.g. at ne=3: 0.49 vs 0.14 vs 0.11). Paolo is close to the 1M reference at ne=2; at ne=3–5 Paolo is lower than the 1M ref and much lower than CCDarkSens.

2. **Patterns 21, 31, 22, 211**  
   **CCDarkSens:** all **zero** for these multi-pixel patterns at ne≥3.  
   **1M reference and Paolo:** both have **substantial** efficiency (21, 31, 22, 211 at ne=3–5). Paolo is similar to or slightly larger than the 1M ref in several bins (e.g. 21@3: 0.41 vs 0.37; 31@5: **0.49** vs 0.065 — Paolo has a very large 31 at ne=5; 211@5: 0.28 vs 0.06).

3. **Pattern 111**  
   Reference: 0.135, 0.034, 0.009 at ne=3,4,5. CCDarkSens: much lower (0.005–0.006). Paolo: **0.13, 0.06, 0.03** — close to 1M ref at ne=3, **higher** at ne=4,5.

**Paolo vs 1M reference:** Same qualitative shape (21, 31, 22, 211 and 111 all non-zero at ne≥3). Main differences: Paolo uses **lambda=0.00015** and **Nsims=100k** (vs 0.00041 and 1M); 31@5 is much larger in Paolo (0.49 vs 0.065); 211@5 larger in Paolo (0.28 vs 0.06); 11@3–5 lower in Paolo.

**Conclusion:** CCDarkSens misses all multi-pixel ROI efficiency (21, 31, 22, 211 and most of 111) → **worse limit**. Both the 1M reference and **efficiencies_paolo.csv** restore that and would strengthen the limit; Paolo’s table is a valid alternative reference (different λ and stats).

---

## Recommendation: use reference overlay

To match the reference limit (same exposure) **without** changing the classifier/MC, use the **efficiency reference overlay**:

1. In `configs/scan_dmelectron_pattern_pydme.json` (or your pattern scan config):
   - Set **`efficiency_csv`** to **`""`** (empty) so efficiencies are generated from PatternMC, **or** keep it pointing to an existing CSV.
   - Add **`efficiency_csv_reference`**: **`"data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv"`** or **`"data/efficiencies_paolo.csv"`** (Paolo: 100k sims, lambda=0.00015; similar shape to 1M ref, some bins higher).

2. The scan will:
   - If `efficiency_csv` is non-empty: load that CSV, then **overlay** (pattern, ne) values from the reference file.
   - If `efficiency_csv` is empty: fill from PatternMC, then **overlay** with the reference file.
   - Write **efficiency_per_pattern.csv** (with reference values for any (pattern, ne) in the reference) and use it for S_pat/B_pat.

So the limit curve uses the reference efficiencies and should align with the solid line. Paths are resolved relative to the config file, so `data/...` works from any working directory.

Long-term, the CCDarkSens pattern MC/classifier can be aligned with the reference (thresholds, 2D scan, DC) so that the generated table matches without overlay.
