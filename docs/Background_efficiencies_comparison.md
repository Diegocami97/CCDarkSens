# Comparison: outputs/scan_pattern/Background_efficiencies.csv vs data/Background_efficiencies.csv

## C++ changes to match the notebook (efficiencies.py)

The following were aligned with the notebook so that re-running the notebook and the C++ app should give consistent results:

1. **Single-pixel (identify_M)**  
   The notebook uses a **sequential threshold check**: return [m] when Pm(q,1)<thr,…,Pm(q,m)<thr and Pm(q,m+1)≥thr (i.e. **largest m** with Pm(q,m) < thr). The C++ no longer uses “argmin Pm”; it now uses this same sequential logic (and checks Pm(q,6) to decide [5] vs none/[0]).

2. **SimulateCluster**  
   The notebook does `np.round(cl + noise, 5)`. The C++ now adds N(0, σ) and rounds to 5 decimals the same way.

3. **MNL (3-pixel)**  
   The notebook uses Pmnl = **min over all 6 permutations** of (m,n,l). The C++ `identify_MNL` now uses the minimum of `pattern_stat_MNL` over the six permutations of (m,n,l) for each candidate label.

4. **MN (2-pixel)**  
   The notebook uses min over the 2 permutations of (m,n); the C++ already did that in `pattern_stat_MN`.

5. **Scan order**  
   MN → MNL → M (2-pixel, then 3-pixel, then 1-pixel) is unchanged and matches the notebook.

If **data/Background_efficiencies.csv** was produced by an older notebook version or different logic (e.g. different single-pixel rule), its numbers can still differ from the current C++ output. Re-run the notebook’s `pattern_simulation_background` with the same Nsims and parameters to compare.

---

## Structure

| | Reference (data/) | C++ output (outputs/) |
|---|------------------|------------------------|
| **Header** | Same: iden_pat, eff_0, eff_1, …, eff_311 | Same |
| **Row order** | String-sorted iden_pat: 0, 1, 11, 111, 112, 113, 12, 121, 122, 13, 131, 14, 2, 21, 211, 212, 22, 221, 23, 3, 31, 311, 32, 4, 41, 5 | Numeric/size order: 0, 1, 2, 3, 4, 5, 11, 12, …, 311 |
| **Rows** | 26 data rows | 26 data rows (same set of ideal patterns) |

So the **column set and row set match**; only the **order of rows** differs. Downstream code that indexes by `iden_pat` (e.g. `eff[eff['iden_pat']==pattern]`) will see the same rows in both files.

---

## Main numerical differences

### 1. Ideal pattern (0) – empty / noise

| | Reference | C++ output |
|---|-----------|------------|
| **eff_0** | ~0.99996 | **1.0** |
| **Other** | Non-zero leakage to eff_1 (~3.1%), eff_11 (~8%), eff_21 (~1.3%), etc. | All 0 |

- **Reference:** When the true pattern is (0), noise is often classified as 1, 11, 21, etc.
- **C++:** When the true pattern is (0), we always classify as (0). So we have no “noise → pattern” misidentification in this row.

Possible causes: different readout noise model, different handling of “no pattern” (e.g. when no seed is found), or different Nsims/RNG.

---

### 2. Ideal pattern (1) – one electron

| | Reference | C++ output |
|---|-----------|------------|
| **eff_1** | ~0.9691 | ~0.9709 |
| **eff_0** | ~4e-5 | ~0.0291 |
| **eff_2** | ~0.0301 | 0 |

- **Reference:** Small misidentification 1→2; almost no 1→0.
- **C++:** Noticeable 1→0 (~2.9%); no 1→2.

So we push more 1-e events to (0) and never to (2); the reference does the opposite.

---

### 3. Ideal patterns (2), (3), (4), (5) – single-pixel 2–5 electrons

| iden_pat | Reference diagonal (eff_X for X=iden_pat) | C++ diagonal |
|----------|-------------------------------------------|--------------|
| 2 | eff_2 ≈ 0.969 | **eff_2 = 0**, eff_1 ≈ 0.9994 |
| 3 | eff_3 ≈ 0.969 | **eff_3 = 0**, eff_1 ≈ 0.9994 |
| 4 | eff_4 ≈ 0.969 | **eff_4 = 0**, eff_1 ≈ 0.990 (+ small eff_11) |
| 5 | eff_5 ≈ 0.969 | **eff_5 = 0**, eff_1 ≈ 0.9987 |

- **Reference:** Strong diagonal: true (2) → identified (2), etc., with ~3% leakage to adjacent (e.g. 2→3).
- **C++:** Almost all 2-, 3-, 4-, 5-e single-pixel events are classified as **(1)**. So our single-pixel branch is effectively returning (1) for all these cases instead of (2)–(5).

This points to a difference in how the **single-pixel (M) branch** is applied (e.g. charge scale, thresholds, or which statistic is used) or how **SimulateCluster** places charge (e.g. one pixel with b electrons vs reference implementation).

---

### 4. Two-pixel patterns (e.g. 11, 21)

| iden_pat | Reference main confusions | C++ main confusions |
|----------|---------------------------|---------------------|
| 11 | eff_11≈0.91, eff_21≈0.074, eff_12≈0.074 | eff_11≈0.91, eff_0≈0.028, eff_1≈0.062 |
| 21 | eff_21≈0.91, eff_11≈0.073, eff_12≈0.073, eff_22≈0.09 | eff_11≈0.98 (identified as 11), eff_0≈5e-4, eff_1≈0.018 |

- **Reference:** Strong permutation confusion: 11↔21, 11↔12, 21↔22, etc.
- **C++:** Little permutation confusion; more “downgrade” to (0) or (1). For ideal 21 we mostly get eff_11 (identified as 11), not eff_21.

So the reference keeps **pattern shape** (two pixels) but swaps digits; we often **collapse to one pixel or empty**.

---

### 5. Three-pixel patterns (111, 211, 221, 311)

- **Reference:** 111→111 (~0.87), 211→211 (~0.87), 221→221 (~0.86), 311→311 (~0.86), with permutation confusion (e.g. 211→111, 211→121, 221→212).
- **C++:** 111→111 (~0.91), 111→11 (~0.09); 211→11 (~0.98); 221→11 (~1.0); 311→11 (~0.98). So we mostly **identify 3-pixel ideals as 11** (two-pixel), not as the correct 3-pixel pattern.

So for 3-pixel ideals, the reference keeps 3-pixel identification with permutation mix; we strongly downgrade to 2-pixel (11).

---

## Summary

| Aspect | Reference (data/) | C++ (outputs/) |
|--------|-------------------|----------------|
| **Row order** | String-sorted iden_pat | Size then lexicographic |
| **Ideal (0)** | Some leakage to 1, 11, 21, … | Always (0) |
| **Single-pixel (1)–(5)** | Strong diagonal; 2↔3 etc. | (2),(3),(4),(5) almost all → (1) |
| **Two-pixel** | Permutation confusion (11↔21, 12) | More 11→0, 11→1; 21→11 |
| **Three-pixel** | Correct 3-pixel with permutations | Mostly → 11 |

So the **structure** (columns, set of rows) matches, but the **entries differ** in ways that suggest:

1. **Single-pixel (M) branch:** Our implementation often assigns (1) where the reference assigns (2)–(5). Worth checking: charge scale from `SimulateCluster`, use of `pattern_stat_M`, and thresholds (thr_M, Qmin_e).
2. **Permutation confusion:** Reference has 11↔21, 12↔21, etc.; we have little of that. That could be from a different **order of checking** (e.g. which of 11 vs 21 is tried first) or different 2-pixel statistics.
3. **3-pixel vs 2-pixel:** We often return 11 for 111/211/221/311. So either we are not finding 3-pixel patterns (e.g. window/segment layout) or we prefer 2-pixel over 3-pixel when both pass.
4. **Ideal (0):** Different noise or “no pattern” handling can explain 0→0 only in C++ vs leakage in the reference.

For **drop-in use** of the C++ CSV in place of `data/Background_efficiencies.csv`, the **row order** may need to be made consistent (e.g. output string-sorted iden_pat like the reference) if the analysis assumes a fixed row order. The **numerical differences** will change background rates in pattern space until the classifier and/or cluster simulation are aligned with the notebook/reference.
