# Band-gap pheno: full process and p100K scaling formula

> **Related:** [band_gap_pheno_ionization.md](band_gap_pheno_ionization.md) (implementation plan), [qcdark2_dielectric_workflow.md](qcdark2_dielectric_workflow.md) (scissor ε).

This note explains what each piece is, how **p100K scaling** works (with **\(E_\mathrm{gap}^\mathrm{new}\)** in the formula), and how it fits into the band-gap phenomenological study.

---

## 1. What you are trying to study

You want to ask: *“If the material behaved as if it had a lower band gap (and maybe a different e–h pair scale), how would DM-electron **limits** change?”*

You are **not** simulating a real 0.1 eV crystal from scratch. You use **Si DFT + tricks**:

1. **Scissor** in QCDark2 → shifts conduction bands → changes **\(dR/dE(E)\)** (how many DM events per eV of recoil energy).
2. **Rescaled p100K table** → changes **\(P(n_e \mid E)\)** (how recoil energy becomes observable electron count).
3. **Same** backgrounds, pattern efficiencies, exposure (for a clean comparison).

Two **independent** knobs in the ionization step:

| Name | Meaning | Typical Si |
|------|---------|------------|
| **\(E_\mathrm{gap}\)** | Minimum energy before ionization can start | ~1.2 eV |
| **\(\varepsilon_h\)** (`eh_pair_eV`) | Scale for “how much energy per e–h pair” (multiplicity spread) | ~3.8 eV |

There is **no** known law \(\varepsilon_h = f(E_\mathrm{gap})\). You can set both to **0.1 eV** (equal-scale pheno), or gap 0.1 eV and \(\varepsilon_h = 3.8\) eV (threshold only), etc.

---

## 2. The full chain (what happens in order)

```text
DM physics (QCDark2 + scissor ε)
        ↓
   dR/dE(E)     events / (kg·year·eV)  at each recoil energy E
        ↓
   × exposure × dE  →  number of events in each E bin
        ↓
   × P(n_e | E)     →  spread into n_e bins  (p100K table)
        ↓
   pattern / diffusion / readout  →  pattern or n_e counts
        ↓
   likelihood  →  limit curve σ_UL(m_χ)
```

**Today’s problem:** step 1 can use scissor gap 0.1 eV, but step 3 still uses **fixed** `p100K_table.csv` (Si: no \(n_e\) below ~1.2 eV). So low-\(E\) signal in \(dR/dE\) is **thrown away** in ionization.

**Fix:** build a **new** `p100K_....csv` for each \((E_\mathrm{gap}^\mathrm{new}, \varepsilon_h^\mathrm{new})\) and point the scan at it.

---

## 3. What `p100K_table.csv` is

It is a **lookup table**, not a formula in code.

- **Rows:** recoil energy \(E\) (eV), from ~1.1 to 50 eV in small steps.
- **Columns:** \(P(n_e{=}1\mid E), P(n_e{=}2\mid E), \ldots\) — probability that **exactly** \(n\) electrons are produced given deposited energy \(E\).

Built from PRD 102.063026 (100 K Si). Reference behavior:

- Below ~**1.2 eV:** essentially all \(P = 0\).
- At ~1.2 eV: usually \(P(n_e{=}1) \approx 1\).
- Higher \(E\): probability spreads to \(n_e = 2, 3, \ldots\) on a scale set by ~**3.8 eV** pair physics.

`ChargeIonization` in CCDarkSens reads this CSV and does:

```text
for each bin of dR/dE at energy E:
    counts = rate × exposure × dE
    for each n_e:
        S(n_e) += counts × P(n_e | E)   # interpolate P from table
```

So **whatever** is in the table at low \(E\) directly controls whether you get signal there.

---

## 4. Why we do not re-run PRD MC for every gap

Regenerating \(P(n_e|E)\) from scratch for arbitrary gap/eh would need the full ionization model from the paper. For a **fast pheno study**, we **derive** new tables from the Si reference by an **energy remap**.

You choose:

- \(E_\mathrm{gap}^\mathrm{ref} = 1.2\) eV, \(\varepsilon_h^\mathrm{ref} = 3.8\) eV (the reference table),
- \(E_\mathrm{gap}^\mathrm{new}\), \(\varepsilon_h^\mathrm{new}\) (your scenario).

---

## 5. The scaling formula (anchored version)

### Idea in words

For each **output** energy \(E\) (on the same grid as today):

1. Decide which **reference** energy \(E'\) has “the same ionization physics” we want at \(E\).
2. When \(E = E_\mathrm{gap}^\mathrm{new}\), that should match reference at **its** gap: \(E' = E_\mathrm{gap}^\mathrm{ref}\).
3. Energies **above** the new gap are stretched/compressed by the ratio of pair scales \(\varepsilon_h^\mathrm{ref}/\varepsilon_h^\mathrm{new}\).
4. Below the new gap: force \(P = 0\).

### Math

**Map output \(E\) → reference \(E'\):**

\[
\boxed{
E'(E) = E_\mathrm{gap}^\mathrm{ref}
+ \bigl(E - E_\mathrm{gap}^\mathrm{new}\bigr)
\times \frac{\varepsilon_h^\mathrm{ref}}{\varepsilon_h^\mathrm{new}}
}
\]

**Probabilities:**

\[
P_\mathrm{new}(n_e \mid E) =
\begin{cases}
0 & \text{if } E < E_\mathrm{gap}^\mathrm{new} \\[6pt]
P_\mathrm{ref}\bigl(n_e \mid E'(E)\bigr) & \text{otherwise}
\end{cases}
\]

\(P_\mathrm{ref}\) is read from `p100K_table.csv` with **linear interpolation** in \(E'\) (same as C++ `ChargeIonization`).

### What each symbol does

| Symbol | Role |
|--------|------|
| \(E\) | Energy on the **new** table’s axis (same grid as reference) |
| \(E'\) | Where you **look up** on the **Si reference** table |
| \(E_\mathrm{gap}^\mathrm{ref}\) | Reference turn-on (~1.2 eV) |
| \(E_\mathrm{gap}^\mathrm{new}\) | **Your** pheno gap — anchors the map |
| \(\varepsilon_h^\mathrm{ref}\) | Reference pair scale (~3.8 eV) |
| \(\varepsilon_h^\mathrm{new}\) | **Your** pheno pair scale |

The factor \(\varepsilon_h^\mathrm{ref}/\varepsilon_h^\mathrm{new}\) is “how many times faster ionization structure changes with \(E\)” compared to Si:

- \(\varepsilon_h^\mathrm{new} = 3.8\) → factor **1** → same spread as Si above threshold (only gap position changes if \(E_\mathrm{gap}^\mathrm{new} \neq 1.2\)).
- \(\varepsilon_h^\mathrm{new} = 0.1\) → factor **38** → very fast change in \(P(n_e|E)\) as \(E\) increases (compressed pheno).

### Older (weaker) map — do not use for equal-scale pheno

Some drafts used only:

\[
E'(E) = E_\mathrm{gap}^\mathrm{ref} + (E - E_\mathrm{gap}^\mathrm{ref}) \times \frac{\varepsilon_h^\mathrm{ref}}{\varepsilon_h^\mathrm{new}}
\]

with \(P=0\) for \(E < E_\mathrm{gap}^\mathrm{new}\). That map does **not** put \(E_\mathrm{gap}^\mathrm{new}\) inside the affine part, so at \(E = E_\mathrm{gap}^\mathrm{new}\) you do **not** automatically sample the reference at its gap. Use the **anchored** formula above instead.

---

## 6. Two examples (same scissor rates, different tables)

Assume **same** `dR/dE` from `Si_fast_gap0p1.h5` (scissor 0.1 eV). Only the p100K table changes.

### A) Threshold only: \(E_\mathrm{gap}^\mathrm{new} = 0.1\) eV, \(\varepsilon_h^\mathrm{new} = 3.8\) eV

Factor \(= 3.8/3.8 = 1\).

| Output \(E\) | \(E'(E)\) | Effect |
|--------------|-----------|--------|
| &lt; 0.1 eV | — | \(P = 0\) |
| 0.1 eV | \(1.2 + (0.1-0.1)\times1 = 1.2\) eV | Same as ref at gap → ionization turns on |
| 1.2 eV | \(1.2 + (1.2-0.1)\times1 = 2.3\) eV | Samples ref at 2.3 eV |
| 5 eV | \(1.2 + 4.9 = 6.1\) eV | Samples ref at 6.1 eV |

**Interpretation:** “Low threshold like 0.1 eV, but **one e–h pair every ~3.8 eV** still (Si-like multiplicity).”

### B) Equal-scale: \(E_\mathrm{gap}^\mathrm{new} = 0.1\) eV, \(\varepsilon_h^\mathrm{new} = 0.1\) eV

Factor \(= 3.8/0.1 = 38\).

| Output \(E\) | \(E'(E)\) | Effect |
|--------------|-----------|--------|
| &lt; 0.1 eV | — | \(P = 0\) |
| 0.1 eV | **1.2 eV** | Turn-on at ref gap |
| 1.2 eV | \(1.2 + 1.1\times38 \approx 43\) eV | Samples ref at **high** \(E\) → lots of multi-\(n_e\) already at modest output \(E\) |

**Interpretation:** “Both scales 0.1 eV” in pheno means **threshold and pair scale are both small** — ionization is very “efficient” per eV above 0.1 eV. This is **aggressive**; not what real Si does (\(\varepsilon_h \gg E_\mathrm{gap}\)).

### C) Reference Si (no scaling)

Use `p100K_table.csv` as-is; scissor 1.2 eV HDF5 and matching rates. \(E_\mathrm{gap}^\mathrm{new} = 1.2\), \(\varepsilon_h^\mathrm{new} = 3.8\) with the anchored formula gives \(E'(E)=E\) (identity map).

---

## 7. The whole process step by step (with what to plot)

### Step A — Dielectric / scissor (already done)

- Run QCDark2 with `scissor_bandgap = 0.1, 0.3, 0.5, 1.2` eV.
- Output: `Si_fast_gap0p1.h5`, etc.
- **Plot:** `dR/dE(E)` for one \((m_\chi, \sigma)\) — curves separate at **low** \(E\).
- **QCDark2 ignores** JSON `band_gap_eV` / `eh_pair_eV`; only the HDF5 matters.

### Step B — Rate grids

- `qcdark2_generate_grid.py` reads each HDF5 → folder of CSVs `dR/dE_...`.
- JSON `detector.band_gap_eV` / `eh_pair_eV` are **labels in file headers** only.

### Step C — Build scaled p100K tables

For each scenario (e.g. `p100K_gap0p1_eh0p1.csv`, `p100K_gap0p1_eh3p8.csv`):

1. Loop over each \(E\) on the grid.
2. If \(E < E_\mathrm{gap}^\mathrm{new}\): write zeros.
3. Else: compute \(E'(E)\), interpolate \(P_\mathrm{ref}(n_e|E')\), write row.
4. **Plot:** \(P(n_e{=}1|E)\) vs \(E\) for ref vs new tables — you **see** turn-on and multiplicity change.

### Step D — Fold to \(S(n_e)\) (one DM point)

- Same `rates_dir` (e.g. gap 0.1 scissor).
- Different `charge_ionization.table_csv`.
- **Plot:** bar chart of \(S(n_e)\) — shows whether low bins go from **zero** (old p100K) to **nonzero** (new table + low-\(E\) rates).

### Step E — Full scan + limits

- Scan with profile likelihood, same background/exposure.
- **Plot:** \(\sigma_\mathrm{UL}(m_\chi)\) overlay; optional ratio to reference.

Figures go under `outplots/band_gap_pheno/` (see [band_gap_pheno_ionization.md](band_gap_pheno_ionization.md) §5).

---

## 8. How scenarios fit together

| Scenario | Scissor HDF5 | \(E_\mathrm{gap}^\mathrm{new}\) | \(\varepsilon_h^\mathrm{new}\) | Ionization file |
|----------|--------------|----------------------------------|--------------------------------|-----------------|
| Reference | `Si_fast_gap1p2` | 1.2 | 3.8 | `p100K_table.csv` |
| Low gap, Si pairs | `Si_fast_gap0p1` | 0.1 | 3.8 | `p100K_gap0p1_eh3p8.csv` |
| Equal 0.1 / 0.1 | `Si_fast_gap0p1` | 0.1 | 0.1 | `p100K_gap0p1_eh0p1.csv` |

Same **`rates_dir`** for the last two rows → limit differences come from **ionization only**.

Manifest template: `configs/band_gap_pheno_scenarios.json`.

---

## 9. What is *not* changed (limitations)

- Si lattice, MOs, scissor only shifts CB energies.
- Pattern efficiencies, diffusion, readout noise.
- Rescaled p100K is a **pheno proxy**, not new ab initio ionization MC.

### How to describe results

| Do say | Do not say |
|--------|------------|
| Pheno study with scissor-shifted ε and rescaled p100K; \((E_\mathrm{gap}, \varepsilon_h)\) chosen per scenario | “Fundamental \(\varepsilon_h(E_\mathrm{gap})\) relation” |
| Show threshold-only and equal-scale at same scissor to bracket ionization | “True 0.1 eV band gap silicon” |

---

## 10. Short recap

1. **Scissor** changes **where** DM energy deposition appears in \(dR/dE(E)\).
2. **p100K scaling** changes **how** that energy becomes \(n_e\), using a map that:
   - aligns **new gap** → **reference gap** (\(E_\mathrm{gap}^\mathrm{new} \to E_\mathrm{gap}^\mathrm{ref}\)),
   - rescales energy above the gap by \(\varepsilon_h^\mathrm{ref}/\varepsilon_h^\mathrm{new}\),
   - zeros everything below \(E_\mathrm{gap}^\mathrm{new}\).
3. **\(E_\mathrm{gap}^\mathrm{new}\)** belongs **inside** \((E - E_\mathrm{gap}^\mathrm{new})\), not only as a separate cutoff — that is what makes “0.1 eV for both” mean “turn on at 0.1 eV like Si turns on at 1.2 eV.”

Implementation: `utils/build_p100K_scaled.py` (TODO) should implement the **anchored** formula; wire via `response.charge_ionization` in scan configs.
