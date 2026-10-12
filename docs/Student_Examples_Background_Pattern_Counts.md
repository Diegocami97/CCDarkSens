<!--
Diego Venegas-Vargas
DAMIC-M collaboration
CCDarkSens Framework

Student_Examples_Background_Pattern_Counts.md -- Pointer to the walkthrough for the random-coincidence background B^rc_p in the six pattern bins, with the account of why the published Bp cannot be reproduced exactly from uniform rates.
-->

# Student Example — Background Pattern Counts from Dark-Current Rates

The step-by-step version of this example is part of
[`Student_Examples_Pattern_Pipeline_Walkthrough.md`](Student_Examples_Pattern_Pipeline_Walkthrough.md), Stage 3. This page keeps the short version and the account of what the tool can and cannot reproduce.

---

## 1. What it does

Given the mean number of 1-, 2- and 3-or-more-electron pixels per pixel readout (`lambda_1`, `lambda_2`, `lambda_3plus`) and the number of pixel readouts, `ccdarksens_compute_brc` returns the random-coincidence background B^rc_p in the six bins `{11}`, `{21}`, `{111}`, `{31}`, `{22}`, `{211}`:

```
B^rc_p = N_pix × Σ_q  P(q) × B[p | q]          P(q) = Π_i Poisson(q_i | λ_{q_i})
```

with `B[p | q]` the matrix of `data/Background_efficiencies.csv` (the probability that an ideal pattern q is identified as p), and the empty pattern weighted 1. The same formula is `Background_pattern_` of pydme.

```bash
build/ccdarksens_compute_brc --config configs/examples/background_pattern_rates.json
build/ccdarksens_compute_brc 3 2 2            # theta = lambda x 10^4
```

The output has two columns that matter: the **diagonal** term (q = p, the ideal pattern identified correctly) and the **misidentified** term (other patterns identified as p; for `{11}` dominated by lone 1 e⁻ pixels, with `B[11|1]` of order 10⁻⁵).

## 2. Check

For theta = (3, 2, 2) the six values are `156.94, 0.023512, 0.047022, 1.5683e-6, 1.0073e-6, 1.1137e-5`, identical to the numbers `DataSimulator.ipynb` of pydme prints. For theta = (2.9, 2, 2) the `{11}` diagonal term is 141.46, the published 141.4.

## 3. Can this reproduce the published Bp?

The published values, hard-coded in pydme's `Background_pattern` and printed in Eq. (5) of the DAMIC-M solar-reflected dark matter paper (arXiv:2607.09352), are

```
B^rc = (141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5)
```

**This tool does not reproduce them exactly, and the reason is known.**

- pydme does not compute these numbers anywhere: `Background_pattern`, the function the limit code calls, contains them as literals. The formula function `Background_pattern_` is only used in a simulation notebook, and this tool matches that to five digits.
- The paper's Eq. (4) is not one uniform rate times N_pix. It sums over every image and every group of 10 columns, with the isolated-pixel rates measured in that group and the number of available adjacent positions after the masks. Non-uniform rates raise the sum, and the effect grows with the number of pixels in the pattern.
- The 2 e⁻ and 3 e⁻ pixel populations are not Poissonian given the 1 e⁻ rate (D2: 184 isolated 2 e⁻ pixels against 70 expected, 17 isolated 3 e⁻ pixels against 7×10⁻³). They have to be measured numbers, which is why no uniform rate set can return all six values.
- `--segments table.csv` (columns `n_pixels`, `lambda_1`, optionally `lambda_2`, `lambda_3plus`, one row per image, CCD or column group) supports the structure of Eq. (4). The per-band inputs are not in this repository.

What is missing to reproduce the published numbers: the isolated-pixel occupancies for 1, 2 and 3 e⁻ per image **and per 10-column group**; the masked number of available positions per pattern and group; the exact version of `B[p|q]` used. Use this tool to understand and vary the background, and quote the published `Bp` for any result meant to reproduce the paper. The full discussion is in [`Pipeline_Raw_Data_to_Limit.md`](Pipeline_Raw_Data_to_Limit.md) §6 and §8.
