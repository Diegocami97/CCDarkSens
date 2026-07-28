# Matching the reference limit curve (dme_limit_curve.pdf)

**Setup:** Both the dotted line (our scan) and the solid line (reference) are for **the same exposure: 1.3 kg-day**. So the ~order-of-magnitude gap in sensitivity is **not** from exposure; it comes from efficiency, background, or other analysis choices.

---

## 1. Pattern efficiency

If our **ε(pattern|n_e)** is lower than in the reference, S_pat is smaller → we need **higher σ** to reach the same q → **worse limit**.

- **`response.pattern_mc.efficiency_csv`** — Must match the efficiency used for the reference (same detector/analysis). If the reference uses a different efficiency table or no extra factor, we should align.
- **`backgrounds.pattern_efficiency`** — e.g. flat ε = 0.95. If this is applied on top of the CSV in our chain but not in the reference, we are over-suppressing signal.
- **Double suppression** — Check we are not applying diffusion/readout/pattern cuts twice (e.g. both in the CSV and again in the pipeline).

---

## 2. Background (Bp + Br)

With **`background_source: "bp_br_template"`**, B_pat = Bp + θ·Br. For the **same 1.3 kg-day exposure**:

- **Higher B_pat** than the reference → we need higher σ for the same q → **worse limit**.
- Bp/Br in the config must be the **expected rates per pattern at 1.3 kg-day**. If they were taken from a different exposure or from “per kg-year” and not scaled to 1.3 kg-day, B_pat will be wrong.
- Confirm with the reference: same 6 patterns, same Bp/Br source (e.g. same fit or same Asimov at 1.3 kg-day).

---

## 3. Rates and σ convention

- **Heavy mediator, F_DM = 1**, same **σ̄_e** definition (cm²).
- **dR/dE** normalization: same as reference (e.g. events/kg/year from the same QEDark/output). Check **`model.rates_dir`** and **`filename_template`** and that we are not applying an extra normalization or form factor.

---

## 4. ROI and analysis

- **Pattern set:** Same 6 patterns as reference: `[11, 21, 111, 31, 22, 211]` (already in pydme config).
- **Profile likelihood:** Same treatment of θ (Bp + θ·Br), same bounds (theta_lo, theta_hi), same prior strength if applicable.

---

## 5. CL / q threshold

90% CL → target q = 2.71. Unlikely to explain a factor of ~10.

---

## Recommended checks (same 1.3 kg-day)

1. **Pattern efficiency** — Compare ε(pattern|n_e) or integrated S_pat at one (m_χ, σ) with the reference. If we are ~√10 lower in efficiency, that can give ~10× worse limit.
2. **B_pat** — Compare B_pat (and Bp, Br) with the reference at 1.3 kg-day. If our B_pat is larger, limits will be weaker.
3. **Signal rates** — At a fixed (m_χ, σ), compare S_pat or dR/dE with the reference to catch normalization or convention differences.
