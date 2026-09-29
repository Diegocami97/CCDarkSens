# WIMP-Nucleon Phase 7 — Slide Bullets

First-person bullet points summarizing the joint 1×1+1×100 likelihood implementation and the DAMIC 2016 (0.6 kg-day) reproduction investigation. Source: `docs/ClusterFitMC_Design.md` §6.7–6.11.

## Slide: Joint 1×1 + 1×100 Likelihood

- I implemented the DAMIC 1×100 readout channel from scratch, including a genuine 1D detector-response model (not just a config change) — the 1×100 mode sums 100 pixel rows in hardware, so I built a proper 1D Gaussian fit path alongside the existing 2D one
- I showed the joint likelihood needs no new statistical machinery: I proved each channel's background nuisance parameter only enters its own likelihood term, so the joint fit is exactly a sum of two independently-computed likelihoods
- I verified the joint channel end-to-end: ΔLL cut and efficiency turn-on both matched the paper's stated ordering between 1×1 and 1×100

## Slide: Two Real Bugs Found and Fixed

- I found that the code folding signal/background through the detector-response kernel was discarding ~90% of the true rate — verified numerically (8065.6 → 816.5 events/kg/year on a real spectrum) and fixed with proper quadrature
- I found the background was using the *signal's* detection efficiency instead of its own — I digitized the paper's real background-efficiency curve directly from the figure's vector PDF data and wired it in
- Together, these fixes closed the mass-dependent divergence I'd originally been asked to investigate — the curve's shape now tracks the paper's own shape correctly

## Slide: Validating Against the Published Result

- I directly extracted the paper's observed curve and expected ±1σ band data from Figure 11 for a rigorous, quantitative comparison
- I compared our projection against the paper's own *expected* sensitivity band (not just their real result) to remove real-data noise from the comparison entirely
- I found our result sits outside their expected band at nearly every mass — meaning there's still a real, unexplained difference between our projection and theirs, independent of any statistical fluctuation in their actual dataset

## Slide: Open Questions for Further Work

- I haven't yet modeled the paper's own stated systematic uncertainties (quenching uncertainty alone shifts their limit by ±1.5× at 2 GeV) — I consider this the most promising next lever
- I flagged that the paper mentions non-uniform detector noise across their real exposure (some runs at ~2.2 e⁻ vs. the 1.8 e⁻ I use as one fixed value)
- I bounded the 1×100 channel's noise assumption via a sensitivity sweep, but can't close it without a number the paper never states
- I identified that the real fix — fitting to the paper's actual candidate events instead of a background-only projection — isn't possible without their per-event data, which isn't public
