# Plan: generic sensitivity-band tool

Plan for adding a **Brazilian-band** (median expected limit ± 1σ / ± 2σ green/yellow band)
generator that works with **any** scan binary in `apps/`, not just the SRDM pattern-CSV
one. This is the canonical "expected sensitivity" representation used in HEP DM-search
papers (CMS, ATLAS, DAMIC-M, etc.).

---

## 1. What the band represents and how toys are sampled

### Conceptual primer: what the band is, why two loops, and what the quantiles are

This subsection is the orientation document for someone reading the plan
cold. The rest of §1 then makes the same statements precise; §2 / §3
implement them.

#### The big picture

The Brazilian band answers a single question:

> *"If our universe had no signal, what range of upper limits would we have
> set, given the statistical fluctuations of the background?"*

To answer it, we need three things:

1. A **procedure** that maps a dataset to a `σ_UL(m_χ)` curve — that is the
   existing scan binary.
2. A way to **simulate many possible datasets** under the no-signal
   hypothesis — pseudo-experiments, a.k.a. toys.
3. A way to **summarize the spread** of `σ_UL` curves across those toys —
   the band itself.

Each toy is one realization of "what would our data have looked like if no
new physics is present?". We sample background-only data many times, run
the scan binary on each, then look at the *distribution* of `σ_UL` at every
mass. The spread of that distribution **is** the band.

#### Why two loops?

The non-trivial part is step 1: how does the scan binary actually extract
`σ_UL` from a single dataset? The PLR recipe is:

For each mass `m_χ` and candidate `σ`, compute the test statistic

```
q_μ(σ) = -2 [ ln L(μ = σ ; θ̂̂) − ln L(μ̂, θ̂) ]
```

and define

```
σ_UL(m_χ) = min { σ : q_μ(σ) ≥ q_target }.
```

That `q_target` is the **rejection threshold for `H_μ` at the chosen
confidence level**. By the Neyman construction it is the upper-tail
percentile of the `q_μ` distribution **under the signal hypothesis `H_μ`**:

```
P( q_μ ≥ q_target  |  H_μ ) = 1 − CL.
```

There are two ways to obtain `q_target`:

- **Asymptotic shortcut** (Wilks/Wald). `q_μ → χ²_1` in the large-N limit,
  so `q_target = (Φ⁻¹(CL))² ≈ 1.642` for one-sided 90 % CL. One number, no
  toys. Cheap. Can fail when statistics are low, when a nuisance hits a
  boundary, or when `q_μ` is non-smooth (signal-pattern degeneracies).
- **Toy-MC** (proper Neyman). Don't trust Wilks. *Measure* the `q_μ`
  distribution under `H_μ` at the relevant `(m_χ, σ)` and read the 90 %
  quantile off it directly. This is what `n_threshold_toys` controls.

The two loops in the band tool correspond to these two completely different
statistical questions:

```
INNER (Phase 1, n_threshold_toys)
   "Given H_μ at (m_χ, σ_threshold), what value does q_μ exceed only
    1−CL of the time?"
   → answer = q_target(m_χ)               (one number per mass)

OUTER (Phase 2, n_toys)
   "Given B-only data, what σ_UL does the procedure yield?"
   → answer = σ_UL distribution           (one full curve per toy)

Aggregate of OUTER → quantiles per mass → the Brazilian band.
```

The inner loop characterizes the **test** itself (it runs once per
analysis). The outer loop characterizes our **expected sensitivity** under
no signal (this is the band you actually plot).

When `band.threshold = "asymptotic"`, the inner loop is skipped entirely:
`q_target` is hard-wired to `(Φ⁻¹(CL))²` and only the outer loop runs.
When `band.threshold = "toy_mc"` or `"both"`, both loops run and the outer
loop's scan binary uses the per-mass `q_target` lookup produced by the
inner loop.

#### Why the inner loop is run at `σ_threshold` and not at every σ

Strictly, the PLR test has a different `q_μ` distribution at every `σ`,
but in the regime where the experiment is sensitive (Wilks/Wald) the
distribution depends primarily on `m_χ`. Cowan-Cranmer-Gross-Vitells (2011)
showed that picking `σ` close to where the limit will fall — i.e. the
asymptotic UL on Asimov data, our `σ_threshold(m_χ)` from Phase 0 — gives
the small-correction toy threshold without scanning the entire `(m_χ, σ)`
grid in toys. This is why **Phase 0** runs the asymptotic scan on Asimov
data once, only to pick that calibration point per mass.

#### Why outer toys are sampled under `H_0` (B-only)

Because the Brazilian band represents *expected sensitivity in the absence
of signal*. We want the limits we would set on average when the universe
is empty — the headline frequentist quantity in any limit-setting paper.
The observed limit (your real data) is then plotted on top of the band so
the reader can judge: did we get lucky, unlucky, or perfectly typical
compared to no-signal expectations?

Sampling outer toys under `H_μ` for some fixed `μ` would instead produce an
*exclusion-power band* — a different and rarely shown plot.

#### Quantiles: what they actually are

A **quantile** `q ∈ [0, 1]` of a distribution is the value below which a
fraction `q` of the probability lies. For a sample of `N` values
`{x_1, …, x_N}`, sort them and the `q`-quantile is approximately
`x_{⌈qN⌉}`, with linear interpolation between neighbors when `qN` is not
an integer (numpy `"linear"` convention; this is what `quantile_inplace`
in `apps/ccdarksens_band.cc` implements).

Examples on a sorted sample of 200 outer-toy `σ_UL` values at one mass:

| Quantile | Index in sorted list (≈) | Meaning |
|:--:|:--:|:--|
| 0.50  | 100–101 | **median** |
| 0.16  | 32      | lower edge of the 1σ-equivalent central interval |
| 0.84  | 168     | upper edge of the 1σ-equivalent central interval |
| 0.025 | 5       | lower edge of the 2σ-equivalent central interval |
| 0.975 | 195     | upper edge of the 2σ-equivalent central interval |

#### Why these five specific quantiles?

For a Gaussian distribution, the central intervals corresponding to "± k
standard deviations" are:

| Interval | Probability content | Lower percentile | Upper percentile |
|:--------:|:--:|:--:|:--:|
| ±1σ | 68.27 % | 0.5 − 0.3413 = **0.1587** | 0.5 + 0.3413 = **0.8413** |
| ±2σ | 95.45 % | 0.5 − 0.4772 = **0.0228** | 0.5 + 0.4772 = **0.9772** |

Rounded to **0.16 / 0.84 / 0.025 / 0.975**, plus the **0.50** median.
That is the entire HEP convention encoded in `band.quantiles`.

Interpretation:

- 68 % of B-only outer toys produce a `σ_UL` inside the **green** band.
- 95 % of B-only outer toys produce a `σ_UL` inside the **yellow** band.
- 50 % above and 50 % below the dashed median curve.

Crucially these are **percentile-based intervals**, not "mean ± standard
deviation". If the distribution of `σ_UL` across toys is skewed — and at
low statistics it usually is (Poisson asymmetry, σ-grid discreteness, q_μ
plateaus) — the band is visibly asymmetric. The "σ" label is a Gaussian-
equivalent shorthand for *probability content*, not a real statement about
the second moment.

#### How the five percentiles become the band's five TGraphs

For each mass `m_χ`, after `K = n_toys` outer toys, you have `K` values of
`σ_UL`. Sort them, take the five quantiles, and you get five points that
get connected across masses to make the five band-edge `TGraph`s:

```
σ_UL distribution at m_χ (over K toys):

  (low tail) ─── 2.5% ─── 16% ─── 50% ─── 84% ─── 97.5% ─── (high tail)
                  │         │       │       │         │
             band_2sigma_   │   median_     │    band_2sigma_
                 low        │   expected    │        high
                      band_1sigma_   band_1sigma_
                          low            high
```

When the plotter draws filled regions between `(low_2σ, high_2σ)` and
`(low_1σ, high_1σ)`, you get the two-tone Brazilian-band ribbon that the
observed limit overlays on top.

That's all the conceptual machinery: the two loops do genuinely different
jobs (one calibrates the test, one builds the band), and the five
quantiles are nothing more than a percentile prescription with
Gaussian-σ-flavored labels. The remaining sections (1 cont., 2, 3, 4) make
each of these statements precise and turn them into a single C++
orchestrator binary `apps/ccdarksens_band.cc`.

### Standard practice: toys are sampled from the **background model**, not from observed counts

The Brazilian band is by convention an **expected sensitivity band under the null
hypothesis (background-only)**. Each toy dataset is drawn as

```
D_i^{(t)} ~ Poisson(B_p^i + B_r^i),   i = 1..N_pattern,  t = 1..K
```

with the nuisance `θ` held at its nominal value of 1.
(`DataSimulator.ipynb` cell 21 does exactly this, splitting `Poisson(Bp_i) + Poisson(Br_i)`,
which is mathematically `Poisson(Bp_i + Br_i)`.) The full limit pipeline is then run on
each toy → one `σ_UL(m_χ)` curve per toy → the **16/50/84 %** and **2.5/97.5 %**
quantiles are taken across toys at each mass to get the median expected limit and the
±1σ (green) / ±2σ (yellow) bands. The **observed limit** from the real
`D = [144, 0, 0, 1, 0, 0]` is then drawn on top of the band as a separate curve.

### Why background-only and not observed?

In roughly increasing importance:

1. **Pre-unblinding / experiment-design meaning.** The band answers
   "what limit would this experiment expect to set if there is no signal?". That is
   defined entirely by the background model, completely independent of what the data
   turned out to be. You can quote it before you ever look at the data.

2. **Frequentist coverage.** `σ_UL` is a function of the data; the band traces the
   sampling distribution of that function under the null ensemble. Quantile properties
   (90 % CL, etc.) are well-defined in this construction. There is no equivalent
   well-defined ensemble if you bootstrap from a single observed dataset.

3. **No circularity / no centering trick.** If you sampled
   `D^{(t)} ~ Poisson(D_obs)` per bin, the band would by construction be centered on
   whatever limit `D_obs` happens to give, so the observed curve would always sit at
   (or very near) the median. The whole point of the band is to *measure* how far the
   observed sits relative to expectation — sampling from observed defeats that.

4. **Communicates the right physics message.** In our specific case
   `D = [144, ...]` is `+2.5` events above `B = [141.5, ...]` in pattern-11. With
   B-only toys the band is centered on the median-expected limit (which corresponds
   to `D ≈ 141`), and the observed curve sits *slightly above* (i.e. a bit looser
   than) the median. That is the message you want: *our observed limit is mildly
   weaker than median-expected because pattern-11 had a small upward fluctuation in
   this dataset*. With observed-counts toys you'd lose that signal entirely.

### Common variants you'll see in HEP papers

For completeness, and so we know what the future extension points are:

| Variant | Toy mean per bin | When it's used |
|---|---|---|
| **B-only nominal** (this plan) | `Bp + 1·Br` | Pre-unblinding sensitivity band, what we want. |
| **B-only post-fit** | `Bp + θ̂·Br` (θ̂ from null fit to data) | "Expected with constrained nuisance" — folds in what data tells you about θ. Subtler interpretation; tends to *narrow* the band. |
| **Asimov / asymptotic** | no toys, analytic | Same band quantiles via Cowan et al. 2011 PLR formalism. Much cheaper for many masses; valid when Poisson counts aren't tiny. |
| **B+S injection** | `Bp + Br + μ·S` | Discovery sensitivity / power study, *not* a 90 %-CL exclusion band. |
| **Bootstrap from observed** | `Poisson(D_obs)` | Almost never used; gives an "uncertainty on the observed result", not an expected sensitivity. |

These extensions are deliberately compatible with the architecture in §2 — they all
reduce to "swap the toy sampler" or "skip the toy loop".

### How `q_target` is set: asymptotic vs toy-MC

A second axis, **independent** of toy sampling, is how each per-toy `σ_UL` is
extracted from the PLR `q_μ(D)`. The upper limit is the σ where
`q_μ(D) = q_target(m_χ, σ)`. The choice of `q_target` is what distinguishes
"Carlos-style" from "Antoine-style" limit setting (cf. the conversation
comparing the two reference notebooks):

| Mode | `q_target` | Cost | When valid | Comment |
|---|---|---|---|---|
| **`asymptotic`** | `q_target = (Φ⁻¹(CL))² ≈ 1.642` for 90 % CL one-sided, **constant** in `(m_χ, σ)` | cheap (sub-second) | Wilks/Wald regime: large counts in the bin(s) that drive the limit | What our existing scan binary does today; matches Carlos / pydme |
| **`toy_mc`** | `q_target(m_χ) = (1 − α)`-quantile of the `q_μ` distribution under `H_μ` at `(m_χ, σ_threshold(m_χ))`, estimated from `N_threshold ~ 10⁴` sub-toys | order(N_threshold) sub-toys per mass, **once** per band run (not per outer toy) | always — this is the proper Neyman construction | What Antoine's notebook does for the test statistic; the gold standard |
| **`both`** (default for first runs) | both the above, in one band run | sum of the two | development & validation | Produces two bands in the same output ROOT; lets you see how much the asymptotic approximation under/over-covers |

**Crucially, all toy data — outer band toys *and* inner threshold sub-toys —
are sampled from the background model alone. Observed counts never enter the
band at any phase.** The toy_mc threshold sub-toys are drawn from
`Poisson(s_i + b_i)` at the *trial* signal hypothesis `H_μ`, which is itself a
pure model prediction; no `D_obs` involved. The σ at which we anchor those
sub-toys (`σ_threshold(m_χ)`) is the **expected (Asimov) limit** computed by
running the scan binary once on Asimov data `D_i = Bp_i + Br_i` — again, the
background model only. The observed limit is drawn on the final plot as a
separate curve, but it does not feed into the band construction at any point.

**Why we want both available, with `both` as the development default.** In our
SRDM regime the (1,1) pattern bin (`B ≈ 141`) dominates the limit, so Wilks
holds well there and `asymptotic` should be excellent (consistent with the
~8 % numerical agreement we already have with pydme). But for low-mass DM
where the signal lives in higher-multiplicity bins (e.g. (3,1), (2,1,1),
(1,1,1) where `B ~ O(0.01)`), the asymptotic approximation can deviate.
Running both side-by-side **measures that deviation directly** rather than
guessing it; once we trust the result we can ship just the better one.

**Where the threshold depends on σ.** The proper Neyman quantile depends on
the full hypothesis `(m_χ, σ)`, but for an upper limit only the value at the
crossing matters. The band tool computes `q_target` at one σ per mass —
specifically at `σ_threshold(m_χ) = σ_UL_asimov_asy(m_χ)`, the **expected**
(Asimov B-only) asymptotic crossing, which is essentially free to compute. A
2-D `(m_χ, σ)` threshold grid is a future refinement, useful only if
`q_target(m_χ, σ)` varies meaningfully in σ near the crossing (it usually
doesn't in our regime).

**Why this couples cleanly to the architecture in §2.** The toy generation
step is unchanged across modes — same B-only Poisson at pattern level. What
changes is the per-toy `σ_UL` extraction inside the scan binary, which we
control with one new optional field `run.q_target_lookup_path`. When absent
(default), the scan binary uses the asymptotic constant `q_target ≈ 1.642`;
when present, it uses the per-mass lookup. Toy generation, signal model,
constraint term, and minimization are all untouched.

---

## 2. Architecture: one config, one binary, reuse the existing scan

Two non-negotiable requirements drive the design:

1. **Reuse the existing minimization code unchanged.** Whatever runs today when you
   call the scan binary is exactly what runs per toy. We do not reimplement the
   likelihood, the constraint term, the bisection, or anything else — we *call*
   the scan binary.
2. **Single JSON config.** The user provides one JSON file that already specifies
   the signal models, the background to sample from, and the run options; we just
   add a small optional `band` block to that same file with the band-specific
   parameters (number of toys, seed, etc.). No second config to keep in sync.

### Goal

The tool must work with **any** of the existing scan binaries — currently:

- `ccdarksens_scan_srdm_pattern_csv`
- `ccdarksens_scan_dmelectron_grid`
- `ccdarksens_scan_dmelectron_pattern`
- `ccdarksens_scan_dmelectron_grid_dualspace`
- `ccdarksens_scan_dmelectron_pattern_Edep`
- `ccdarksens_scan_dmelectron_toy`

— and any future scan binary, **without** the band tool knowing anything about the
underlying physics, observable space, or signal model.

### Design choice: subprocess orchestrator that reads the same config the scan reads

Three architectures were considered:

1. **Pure C++ in-process app, SRDM-specific** (the previous draft of this document).
   Fast per toy, no subprocess overhead. Rejected because it bakes in SRDM CSV loading
   and is not generic, and because it duplicates the minimization code that already
   lives in the scan binary.

2. **Generic in-process via abstract `LimitFitter` C++ interface.** Refactor every scan
   binary so its core is a `LimitFitter` subclass with `SetData()` / `ComputeUL(cl)`;
   the band binary takes a fitter and runs toys in-memory. Fast, but requires
   substantial refactor of six existing binaries — violates "reuse the existing
   minimization code unchanged" and turns into a maintenance commitment forever.

3. **Generic subprocess orchestrator, reading the same JSON the scan reads.** A small
   new binary (`ccdarksens_band`) that knows nothing about any specific analysis. The
   user adds a `band` block to their existing scan config; running
   `ccdarksens_band <config>` generates toys, patches the per-toy data into the same
   config (with `band` stripped and `observed_counts` / `outdir` injected), invokes
   the scan binary as a subprocess, and harvests the resulting
   `upper_limit_sigma_e_mchi_graph` from the output ROOT. Works with all six binaries
   today and any future one, with **zero refactor** of the scans themselves.

We pick architecture 3. It satisfies both requirements: the existing scan binary's
minimization is reused exactly (one subprocess invocation per toy), and the user
sees one JSON file with one new block in it.

### User experience

The user takes their existing scan config, e.g.
`configs/scan_srdm_pattern_csv_per_mass_pydme_match.json`, and adds an optional
`band` block at the top level (see Step 2 below). Then:

```bash
# Nominal observed limit (as today, unchanged):
build/ccdarksens_scan_srdm_pattern_csv \
    configs/scan_srdm_pattern_csv_per_mass_pydme_match.json

# Expected sensitivity band (new). The band block in the same config selects
# threshold: "asymptotic" (default, fast), "toy_mc" (proper Neyman, slower),
# or "both" (recommended for first runs — produces both bands in one output
# ROOT for direct comparison).
build/ccdarksens_band \
    configs/scan_srdm_pattern_csv_per_mass_pydme_match.json
```

Same config, two commands, one for the observed limit and one for the band. The
scan binary ignores the `band` block; the band binary internally strips it before
passing per-toy configs to the scan binary.

### Contract that scan binaries must satisfy

For a scan binary to be drivable by `ccdarksens_band`, it must:

1. **Take a single JSON config path as its CLI argument** — already true for all scans.
2. **Honor `run.observed_counts`** as an optional override of the data vector
   (Asimov default if omitted) — already true for `ccdarksens_scan_srdm_pattern_csv`;
   need to verify / port to the other five binaries (small, uniform change).
3. **Honor `run.outdir`** — already true for all scans.
4. **Write a `TGraph` of `(m_χ, σ_UL)` named `upper_limit_sigma_e_mchi_graph`** to the
   output ROOT — already true for `ccdarksens_scan_srdm_pattern_csv`; need to verify /
   port to the other five (also small).
5. **Tolerate (ignore) unknown top-level keys** like `band` — `nlohmann::json` parsing
   already does this; just don't add strict-key validation.
6. **Honor an optional `run.q_target_lookup_path`** (string). When present,
   point at a ROOT file containing a `TGraph` named `q_target_per_mass` of
   `(m_χ, q_target)` pairs (one entry per mass in the scan); the binary uses
   **that** value as the PLR crossing threshold instead of the asymptotic
   `(Φ⁻¹(CL))²`. When absent, behavior is unchanged from today (asymptotic
   threshold). One-time patch per scan binary, ~10 lines.
7. **Support an optional `run.mode = "threshold_toys"`** (default `"scan"`).
   In `threshold_toys` mode the binary does *not* compute `σ_UL`; instead it
   reads `run.threshold_toys` (see below), runs `N_threshold` Poisson(`s_i +
   b_i`) sub-toys per mass at the configured `(m_χ, σ_threshold(m_χ))` point,
   computes `q_μ` for each, and writes a `TGraph` named `q_target_per_mass`
   of `(m_χ, c_μ)` where `c_μ` is the configured percentile. Reuses the same
   `ProfileLikelihood` minimization that scan mode uses — no new likelihood
   code, only a different driver loop.

Items 1–5 are the **base contract** required for the asymptotic band; items 6
and 7 are the **toy-MC threshold extension** required for the toy-MC band. A
scan binary that only implements 1–5 still gets a usable asymptotic band; the
band tool detects this and disables `toy_mc` / `both` modes for that analysis
with a clear error.

Anything else (Spat templates, signal models, observable space, nuisance treatment)
is the scan binary's private business.

### The phases the band tool drives

```
                    ┌──────────────────────────────────────────────┐
                    │             ccdarksens_band                  │
                    │  reads cfg, owns toy generation & quantiles  │
                    └─────────────────┬────────────────────────────┘
                                      │ subprocess invocations
                                      ▼
   Phase 0 (only in toy_mc / both): ONE call to scan_binary in scan mode, with
   run.observed_counts STRIPPED so the binary falls back to its Asimov default
   D_i = Bp_i + Br_i.  Result: σ_UL_asimov_asy(m_χ), saved as
   sigma_threshold_per_mass.

   Phase 1 (only in toy_mc / both): ONE call to scan_binary in
   threshold_toys mode.  Reads run.threshold_toys (the σ_threshold graph from
   Phase 0, N_threshold, percentile); writes q_target_per_mass TGraph.

   Phase 2 (always): K_outer calls to scan_binary in scan mode.  Each toy
   passes its own D^{(t)} via run.observed_counts and, in toy_mc / both,
   the q_target lookup from Phase 1 via run.q_target_lookup_path.
```

**No phase ever uses `D_obs`.** Phase 0 strips `run.observed_counts` from the
config it passes to the scan binary, so the underlying scan defaults to
Asimov B-only data. Phase 1's inner sub-toys are drawn from
`Poisson(s_i + b_i)` under `H_μ`. Phase 2's outer toys are drawn from
`Poisson(Bp_i + Br_i)` under the null. The user's observed `D_obs = [144, 0,
0, 1, 0, 0]` is plotted on top of the band as a separate curve produced by a
*separate* nominal scan run — it is never an input to the band itself.

In `both` mode the band tool runs Phase 0 + Phase 1 once and then Phase 2
*twice* per toy (asymptotic pass with no `q_target_lookup_path`, then toy-MC
pass with the lookup attached). The outer-toy `D^{(t)}` realizations are
seeded from the same RNG stream in both passes, so the two bands are computed
on **identical Poisson realizations**, making the comparison per-toy.

### Why subprocess overhead is acceptable here

Per-toy overhead = (process spawn + ROOT init + signal-model load + Spat rebuild +
ROOT-file write/read). For the SRDM scan that is ~0.5–1 s of fixed cost on top of
~0.5–2 s of actual fit work. Parallelism over `n_workers` threads makes the wall
clock manageable:

| K | serial wall time | 8-way parallel |
|---|---|---|
| 200 | ~5–10 min | ~40–80 s |
| 1000 | ~25–50 min | ~3–6 min |
| 10000 | ~4–8 hr | ~30–60 min |

For typical Brazilian-band publications K = 1000 is plenty; the 8-way parallel cost
is small and a one-time investment per analysis.

If a specific analysis ever becomes the bottleneck, the scan binary can be optimized
internally (e.g. cache pre-loaded templates in a `mmap`'d file, reduce ROOT I/O) — and
those optimizations also benefit the nominal scan, not just the band tool. The band
tool itself stays unchanged.

---

## 3. Implementation plan

Five concrete pieces, in order. Each is committable on its own.

### Step 1 — Verify and uniformize the contract on all scan binaries

Two-tier porting effort: **base contract first** (gives every analysis the
asymptotic band immediately), then **threshold extension** (adds toy-MC mode
where it's wanted).

#### 1a. Base contract (gates `band.threshold = "asymptotic"`)

For each of the six scan binaries:

- **Read of `run.observed_counts`.** Verify the binary accepts an optional
  `observed_counts` array and uses it as the data vector when present (matching
  `ccdarksens_scan_srdm_pattern_csv`'s behavior). If absent, port the same
  ~15 lines of parsing + `SetData(observed)` logic.
- **Output of `upper_limit_sigma_e_mchi_graph`.** Verify the binary writes a `TGraph`
  with this exact name. If it currently only writes a `TH1D`, also build the `TGraph`
  with the *exact* `mchi` values (we already did this fix for the SRDM binary).

These changes are uniform and small (~30 lines per binary) and independently
useful — `observed_counts` benefits every analysis even without bands.

#### 1b. Threshold extension (gates `band.threshold = "toy_mc" | "both"`)

Two additions, each ~30 lines, both only activated by new optional config
keys (zero behavior change when those keys are absent):

- **`run.q_target_lookup_path` (scan mode).** Currently the bisection in the
  scan binary uses `q_target = (Φ⁻¹(cl))²`. Add: if `run.q_target_lookup_path`
  is non-empty, open that ROOT file, read `TGraph q_target_per_mass`, and use
  its value at the current mass instead of the asymptotic constant. Per-mass
  lookup is exact (mass grid matches by construction); no interpolation
  needed.

- **`run.mode = "threshold_toys"`.** New driver loop, *reusing the same
  `ProfileLikelihood` instance* the scan mode uses:

  ```
  for each m_χ in the scan grid:
      σ_thr = run.threshold_toys.sigma_threshold_per_mass[m_χ]    // see step 2
      s_i = signal expected counts at (m_χ, σ_thr)                // existing helper
      b_i = Bp_i + 1 · Br_i                                       // θ=1 nominal
      λ_i = s_i + b_i
      for k = 1..run.threshold_toys.n_threshold_toys:
          n_i^(k) ~ Poisson(λ_i)                                  // independent RNG stream
          plnll.SetData(n_i^(k))
          nll_top  = plnll.MinimizeOverScaleMinuit(σ_thr, ...).second   // numerator
          nll_glob = plnll.MinimizeOverSigmaAndTheta(...).first         // denominator
          q^(k)    = max(0, 2·(nll_top − nll_glob))
      c_μ(m_χ) = quantile(q[], run.threshold_toys.percentile)     // default 0.9
  write TGraph q_target_per_mass (m_χ, c_μ) to <outdir>/qtarget_threshold.root
  ```

  No new likelihood code: we call the existing `MinimizeOverScaleMinuit` and
  `MinimizeOverSigmaAndTheta` from `ProfileLikelihood`. The only new things
  are the Poisson sampler, the inner loop, and the `std::nth_element`
  quantile.

  The `run.threshold_toys` block (consumed only in this mode) is:

  ```json
  "threshold_toys": {
    "sigma_threshold_graph_path": "tmp/sigma_threshold.root", // TGraph "sigma_threshold_per_mass"
    "n_threshold_toys": 10000,
    "percentile": 0.90,
    "rng_seed": 23456
  }
  ```

  The band tool writes `sigma_threshold.root` itself in Phase 0 (one TGraph:
  `(m_χ, σ_UL_asimov_asy)` from the Asimov-data scan run) before invoking
  threshold mode — the scan binary just reads it.

If a scan binary doesn't fit naturally (e.g. it's not a per-mass UL setter),
it just won't be drivable by the band tool, which is fine — bands only make
sense for limit scans. A binary with only the base contract (1a) but not the
threshold extension (1b) gets `asymptotic` bands; the band tool errors out
with a clear message if `toy_mc` / `both` is requested.

### Step 2 — New executable `ccdarksens_band` (reads the same JSON the scan reads)

**File:** `apps/ccdarksens_band.cc`  
**Build:** add `add_executable(ccdarksens_band ...)` to `CMakeLists.txt`.

The user's existing scan config gets one new optional top-level block, `band`. Every
other key in the file is exactly what `ccdarksens_scan_*` already consumes — we don't
duplicate anything. Concretely, for the SRDM-pattern-CSV case
(`configs/scan_srdm_pattern_csv_per_mass_pydme_match.json`) the user adds:

```json
{
  "experiment": { ... },          // unchanged: signal models live here
  "models":     { ... },          //   (or wherever the scan reads them today)
  "run": {
    "background_Bp": [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5],
    "background_Br": [0.039, 0.039, 0.016, 0.052, 0.011, 0.035],
    "observed_counts": [144, 0, 0, 1, 0, 0],
    "outdir": "outputs/scan_srdm_pattern_csv_per_mass_pydme_match",
    ...
  },

  "band": {
    "_comment": "Optional. Read by ccdarksens_band; ignored by the scan binary.",
    "scan_binary": "build/ccdarksens_scan_srdm_pattern_csv",
    "n_toys": 200,
    "rng_seed": 12345,
    "n_workers": 8,
    "outdir": "outputs/band/srdm_pattern_csv_per_mass_pydme_match",
    "threshold": "both",
    "n_threshold_toys": 10000
  }
}
```

Three required `band` keys (`scan_binary`, `n_toys`, `outdir`) and a couple of
optional ones with sensible defaults. Everything else is auto-discovered:

| Auto-discovered field | How |
|---|---|
| Models / signal grids | Already specified in the rest of the config; the scan binary picks them up. |
| Background means for Poisson sampling | Read directly from `run.background_Bp` and `run.background_Br` in the same config. |
| Output ROOT filename per toy | Glob `*.root` inside the toy's `outdir`; pick the one containing `upper_limit_sigma_e_mchi_graph`. |
| UL graph name in the output ROOT | Hard-coded to `upper_limit_sigma_e_mchi_graph` (the contract from §2). |

Full optional `band` keys (with defaults):

```json
"band": {
  "scan_binary": "...",                          // required
  "n_toys": 200,                                 // required
  "outdir": "outputs/band/...",                  // required
  "rng_seed": 12345,                             // default 12345
  "n_workers": 0,                                // 0 = std::thread::hardware_concurrency()
  "save_per_toy_curves": false,                  // default false
  "keep_per_toy_outputs": false,                 // default false (cleaned after success)
  "abort_on_failure": false,                     // default false (skip failing toys)
  "min_success_fraction": 0.95,                  // error out if survival < this
  "tmp_dir": "<outdir>/_toys",                   // default
  "background_bp_key": "run.background_Bp",      // dotted path, default
  "background_br_key": "run.background_Br",      // dotted path, default
  "ul_graph_name": "upper_limit_sigma_e_mchi_graph",
  "quantiles": {
    "median": 0.50,
    "low_1sigma": 0.16,  "high_1sigma": 0.84,
    "low_2sigma": 0.025, "high_2sigma": 0.975
  },

  // q_target threshold (see "How q_target is set" in §1)
  "threshold": "asymptotic",                     // "asymptotic" | "toy_mc" | "both"
  "n_threshold_toys": 10000,                     // sub-toys per mass for toy_mc
  "threshold_percentile": 0.90,                  // 1 − α for one-sided UL at CL
  "threshold_seed": 23456                        // independent RNG stream from outer toys
}
```

The dotted-path keys (`background_bp_key`, `background_br_key`) let the same band
tool work with scans that name their B-only sampling means slightly differently —
e.g. a future scan that uses `experiment.background.expected_counts` instead of
`run.background_Bp`. Default values match `ccdarksens_scan_srdm_pattern_csv`.

#### Algorithm

The orchestration has up to **three phases** depending on `band.threshold`:

```
            ┌──────────────────────────────────────────────────────────────┐
            │ Phase 0 — Asimov asymptotic UL (only if toy_mc / both)       │
            │ One subprocess: scan_binary with observed_counts STRIPPED    │
            │ (scan falls back to its Asimov default D_i = Bp_i + Br_i)    │
            │ → reads σ_UL_asimov_asy(m_χ); writes sigma_threshold.root   │
            └─────────────┬────────────────────────────────────────────────┘
                          ▼
            ┌──────────────────────────────────────────────────────────────┐
            │ Phase 1 — q_target threshold toys (only if toy_mc / both)    │
            │ One subprocess: scan_binary --mode=threshold_toys            │
            │ N_threshold sub-toys × N_masses → qtarget_threshold.root     │
            └─────────────┬────────────────────────────────────────────────┘
                          ▼
            ┌──────────────────────────────────────────────────────────────┐
            │ Phase 2 — outer toy loop (always)                            │
            │ K_outer subprocesses (parallel): scan_binary on toy data     │
            │ D^{(t)} ~ Poisson(Bp + Br), invoked once with asymptotic     │
            │ q_target (asymptotic mode) AND/OR once with q_target lookup  │
            │ (toy_mc mode). Same D^{(t)} for both passes.                 │
            └──────────────────────────────────────────────────────────────┘
```

Detailed steps:

1. **Parse & validate.** Read input config; extract the `band` block; verify
   required keys (`scan_binary`, `n_toys`, `outdir`). Validate
   `band.threshold ∈ {asymptotic, toy_mc, both}`. If `toy_mc` / `both` and the
   scan binary doesn't advertise threshold support (see §2 contract item 7),
   error out with a clear message.

2. **Resolve sampling means.** Resolve `Bp` and `Br` from the same config via
   `band.background_bp_key` / `band.background_br_key`. Compute outer-toy
   Poisson means `λ_i = Bp_i + Br_i`; verify non-negative, same length.

3. **Phase 0 — Asimov asymptotic UL** *(only if `threshold ∈ {toy_mc, both}`)*.
   Build a Phase-0 config: stripped input config (no `band` block), with
   **`run.observed_counts` removed**, and `run.outdir = tmp_dir/phase0`.
   Removing `observed_counts` makes the scan binary fall back to its Asimov
   default (`D_i = Bp_i + Br_i`) — *the same B-only model that drives toy
   sampling*. Run the scan binary once; read `upper_limit_sigma_e_mchi_graph`
   from its output (this is `σ_UL_asimov_asy(m_χ)`). Write a one-graph ROOT
   file `tmp_dir/sigma_threshold.root` with TGraph `sigma_threshold_per_mass`.
   Cost: one scan invocation (~0.5–2 s for SRDM). The observed counts in the
   user's config are *only* read by the band tool to assert their length
   matches `Bp / Br`; they are stripped from every per-phase config the band
   tool builds and never reach a phase that contributes to the band.

4. **Phase 1 — `q_target` threshold toys** *(only if `threshold ∈ {toy_mc,
   both}`)*. Build a Phase-1 config: stripped input config + injected
   `run.mode = "threshold_toys"`, `run.threshold_toys` block (pointing at the
   `sigma_threshold.root` from Phase 0, with `n_threshold_toys`, `percentile`,
   `rng_seed = band.threshold_seed`), and `run.outdir = tmp_dir/phase1`.
   `run.observed_counts` is also stripped here (threshold mode does not need
   data). Run the scan binary once. On completion read the output
   `qtarget_threshold.root` containing `q_target_per_mass`. Cost:
   `n_threshold_toys × N_masses` sub-toy minimizations (~30 s for SRDM at
   `n_threshold_toys = 10⁴`).

5. **Phase 2 — outer toy loop** *(always)*. Spawn a thread pool of
   `n_workers` threads. For `t = 0..K_outer−1`:
   - Seed `std::mt19937_64` with `rng_seed + t` (independent stream per toy).
   - Sample `D^{(t)}_i ~ Poisson(λ_i)` (background-only nominal).
   - Build the per-toy config: stripped input config + `run.observed_counts =
     D^{(t)}` (the *toy* counts, not `D_obs`) + `run.outdir =
     tmp_dir/toy_<t>/asy` (or `/toymc` for the toy_mc pass).
   - **Asymptotic pass** *(if `threshold ∈ {asymptotic, both}`)*: invoke the
     scan binary with the per-toy config as written (no
     `q_target_lookup_path`). Harvest `σ_UL_asy^{(t)}(m_χ)`.
   - **Toy-MC pass** *(if `threshold ∈ {toy_mc, both}`)*: invoke the scan
     binary with `run.q_target_lookup_path` = path to
     `qtarget_threshold.root` from Phase 1. Harvest `σ_UL_toy^{(t)}(m_χ)`.
   - Both passes use the **same `D^{(t)}`** so the two bands are computed on
     identical Poisson realizations, making per-toy comparison meaningful.
   - If `keep_per_toy_outputs == false`, delete the toy dir.

6. **Quantiles.** After all outer toys complete, for each threshold mode that
   ran, pivot the K curves into a per-mass matrix and compute the five
   quantiles via `std::nth_element`.

7. **Write `<band.outdir>/band.root`.** Per threshold mode that ran, write:

   | Object | Asymptotic suffix | Toy-MC suffix |
   |---|---|---|
   | `median_expected_sigma_e_mchi` | `_asymptotic` | `_toy_mc` |
   | `band_1sigma_sigma_e_mchi`     | `_asymptotic` | `_toy_mc` |
   | `band_2sigma_sigma_e_mchi`     | `_asymptotic` | `_toy_mc` |

   In single-mode runs (only `asymptotic` or only `toy_mc`) the suffix is
   still written so plotting code is uniform; the plotter accepts either
   flavor. Plus one `meta` TTree with full provenance (`n_toys`, `rng_seed`,
   `n_threshold_toys`, `threshold_percentile`, both seeds, scan binary path,
   input config path, n_successful_toys, n_threshold_toys_used,
   `sigma_threshold_per_mass` graph, `q_target_per_mass` graph) and, if
   `save_per_toy_curves == true`, a `TTree` `band_per_toy` with branches
   `(toy_idx, mchi, sigma_UL_asy, sigma_UL_toy)` so we can re-quantile
   externally without re-running.

#### Implementation notes

- **Reuse, don't reimplement.** The minimization, profiling, bisection, and constraint
  term all live in the scan binary; the band binary never touches them. Each toy is
  literally one `system()` call to the same binary the user runs today.
- **Subprocess invocation:** `std::system` is acceptable for a first version
  (synchronous, blocks the worker thread). The thread pool gives us parallelism.
  If we ever need finer control (timeout, kill), switch to `popen`/`fork+exec` or
  `boost::process` later.
- **JSON manipulation:** reuse the project's existing `nlohmann::json` dependency.
  Stripping the `band` block is `cfg.erase("band")`; injecting fields is direct
  assignment.
- **ROOT reading:** standard `TFile::Open` + `Get<TGraph>(...)`.
- **Reproducibility:** with `n_workers > 1` we still get bit-identical results because
  each toy's RNG is seeded independently (`rng_seed + t`), not from a shared stream.
- **Robustness:** if a toy subprocess exits non-zero or its ROOT file is missing the
  UL graph, log the toy index and the path to its log file. If
  `abort_on_failure == false` (default), skip it and continue. Quantiles are computed
  on the surviving set; if survival fraction drops below `min_success_fraction`,
  error out before writing the band ROOT.

### Step 3 — Plotter extension to draw the band

**File:** `apps/ccdarksens_plot_dmelectron_limit.cc`

- Add a CLI flag `--band <band.root>` that, when present, autodetects which
  threshold variants are in the file (looks for `_asymptotic` and `_toy_mc`
  suffixed graphs). At least one must be present.
- Optional `--band-mode <asymptotic|toy_mc|both>` (default `both` if both
  present, else whichever one is) for explicit selection.
- Drawing order (back-to-front), single-mode: yellow ±2σ band → green ±1σ
  band → black dashed median expected → black solid observed (from the
  existing `upper_limit_sigma_e_mchi_graph` of the nominal scan).
- Drawing order (back-to-front), `both` mode: 2σ_asymptotic (faded yellow) →
  1σ_asymptotic (faded green) → 2σ_toy_mc (saturated yellow) → 1σ_toy_mc
  (saturated green) → median dashed black (toy_mc) → median dotted gray
  (asymptotic) → observed solid black. The faded vs saturated palette makes
  the difference between threshold modes visually obvious without a busy
  legend.
- Legend: "Median expected, asymptotic threshold", "Median expected, toy-MC
  threshold", "± 1σ expected", "± 2σ expected", "Observed (90 % CL)".
- Color recipe (single-mode default): 2σ = `kOrange-4`, 1σ = `kGreen+1`,
  median = black dashed, observed = black solid. Matches DAMIC-M, ATLAS, CMS
  conventions.

This works with any analysis as long as the band ROOT exists (it always will, since
the graph names are fixed by the band tool).

### Step 4 — Add a `band` block to the existing scan configs (no new config files)

We don't ship parallel `band_<analysis>.json` files. The user's existing scan config
already specifies the signal models and the background to sample from; we just append
a `band` block to it. Same JSON, same scan command works as before, and now
`ccdarksens_band <same_json>` produces the band.

Concretely for SRDM:

- Edit `configs/scan_srdm_pattern_csv_per_mass_pydme_match.json` to add the
  `band` block shown in Step 2. Suggested progression:
  1. **Smoke test** (`threshold: "asymptotic"`, `n_toys = 200`,
     `n_workers = 8`): ~1 min, validates the orchestrator and the asymptotic
     band shape.
  2. **Threshold test** (`threshold: "both"`, `n_toys = 200`,
     `n_threshold_toys = 10000`): ~2–3 min, validates Phase 0 / Phase 1
     plumbing and produces a first comparison plot of the two band variants.
  3. **Publication run** (`threshold: "both"`, `n_toys = 1000`,
     `n_threshold_toys = 10000`): ~10–15 min on 8 workers; both bands quoted
     with full quantile resolution.

For other analyses, once their scan binary is verified to satisfy the §2 contract
(Step 1), the user just adds the same kind of `band` block to whatever scan config
they already use for that analysis.

If the user prefers keeping the band config separated for cleanliness — e.g. so the
"observed limit only" run doesn't carry an unused `band` block — they can keep two
copies of the same config differing only by the presence of `band`. Both work.

### Step 5 — Documentation & validation

- Add a section to `docs/Limit_curve_1kgyear_checklist.md` referring to this app.
- Add a short README at `apps/README_band.md` describing the contract from §2 so
  authors of future scan binaries know what they need to provide.
- Validation per analysis (run once at adoption):
  - Sanity 1: at very small `K` (e.g. K = 5), keep `keep_per_toy_outputs == true`,
    inspect the per-toy curves directly, verify the median is sensible.
  - Sanity 2: at `K = 1000`, the **median expected** curve should sit close to but
    slightly tighter than the **observed** curve in the pattern-11-dominated mass
    region for SRDM (because `D_obs = 144 > B_pat[0] ≈ 141.5` is a +2.5 upward
    fluctuation).
  - Sanity 3: ±1σ band width at a representative mass should be a small fraction of
    a decade (Poisson σ on the dominant template propagates roughly linearly to σ_UL).
  - Sanity 4 (MC convergence): regenerate band with `n_toys = 200` vs `n_toys = 1000`
    and verify the median moves by ≪ 1σ band width.
  - Sanity 5 (threshold comparison, with `band.threshold = "both"`): in the
    SRDM regime where the (1,1) bin dominates, the asymptotic and toy-MC
    bands should overlap to within a few percent in σ at every mass. A larger
    discrepancy at a specific mass is a real physics signal (Wilks failing
    there) and is exactly what we want to *measure* with this comparison, not
    a bug. Document any per-mass deviation > 10 % as an analysis caveat.
  - Sanity 6 (q_target sanity): when the threshold ROOT is written, the
    `q_target_per_mass` values should sit in `[1.0, 3.0]` for a 90 %
    one-sided UL. Values < 1 indicate too few threshold toys for the upper
    tail; values > 3 indicate either a deep non-Wilks regime or a buggy fit
    at that mass. Increase `n_threshold_toys` if the per-mass MC error on
    `c_μ` exceeds ~5 %.

---

## 4. Verification: physics-input parity with the reference frameworks

Cross-check of every physics input in `DataSimulator.ipynb` (Carlos) and
`LBCLimit_from_Antoine.ipynb` (Antoine) against our SRDM scan config and
`apps/ccdarksens_scan_srdm_pattern_csv.cc`. The conclusion is that the band tool
needs **no additional physics inputs** beyond what the existing scan binary
already reads, because each toy reuses the unchanged scan binary on a copy of the
same config — only `observed_counts` is replaced with a per-toy Poisson sample.

### What the reference notebooks actually compute (and what they don't)

Important framing note: **neither reference notebook actually produces a
sensitivity band**. We use them as templates for two distinct ingredients
that the band tool combines:

| Notebook | What it computes | Toy generation source | Role of `D_obs` |
|---|---|---|---|
| `DataSimulator.ipynb` (Carlos) | Counts-level cross-check (cell 21–22). 10 000 pattern-count toys plotted as histograms vs the observed counts; *no* `σ_UL`, *no* limit setting, *no* band. | `Poisson(Bp_i + Br_i)` — **background model** (cell 21). | Drawn as a single reference line on top of the toy-counts histogram (cell 22 `axvline(144, …)`). |
| `LBCLimit_from_Antoine.ipynb` | One **observed** upper limit using a toy-MC critical value for the test statistic (the `Prob(muS)` function). No outer band loop — `data` is fixed at `Data0 = [144, 0, 0, 1, 0, 0]`. | `Poisson(μS·s_i + Bkg_i)` — **`H_μ` signal+background** at the trial signal strength being tested. | Used as the comparator `T_mu(D_obs, …)` against which the toy `T_μ` distribution is compared to find the 90 % CL UL. |
| **Our band tool (planned)** | Brazilian band: median expected `σ_UL` ± 1σ / ± 2σ envelope across `K_outer ~ 1000` outer toys. | Outer (band) toys: `Poisson(Bp_i + Br_i)` — **background model**, by analogy to Carlos cell 21. Inner (test-statistic) toys (`threshold = "toy_mc"` only): `Poisson(μ·s_i + b_i)` — **`H_μ`**, exactly Antoine. | Plotted on top as a separate observed-UL curve from a *separate* nominal scan run; never an input to the band itself. |

Both reference notebooks are consistent with the discipline "draw model-based
toys, plot the observed result as a separate overlay" — Carlos's pattern-count
toys come from the background model, Antoine's test-statistic toys come from
`H_μ`, and neither notebook ever uses `D_obs` to *generate* toys. The band
tool inherits exactly this discipline. The new ingredient that neither
notebook supplies is the **outer toy loop that produces a distribution of
`σ_UL` curves**, which is what the Brazilian band quantiles. The phrasing
"Carlos-style band" / "Antoine-style band" used loosely in earlier drafts of
this document was inaccurate; what those labels really refer to is just the
choice between asymptotic (Carlos / pydme) and toy-MC (Antoine) calibration
of `q_target`. Both styles are available via `band.threshold`.

### Where each physics input lives

| Physics input | DataSimulator (Carlos) | LBCLimit (Antoine) | Our config / scan | Where it enters per toy |
|---|---|---|---|---|
| Per-pixel read noise σ | `read_noise = 0.16` (cells 2–4) | upstream of `Bkg` | absorbed in `background_Bp` | **upstream — band tool does not see it** |
| Per-pixel dark current λ | `lamb = 3e-4` (cells 2–4) | upstream of `Bkg` | absorbed in `background_Bp` | upstream — same |
| CCD geometry | rows × cols × numbCCDs × loop_sizes | aggregated | absorbed in `Bp / Br` | upstream — same |
| `N_img` (image count for prior) | sum(loop_sizes) ≈ 20333 | LBC dataset value | `constrain_n_bins = 4450` | per-toy fit reads it |
| `N_pix` total | `Nusedpix = 1.998e9` | LBC dataset value | absorbed in `Bp / Br` | upstream |
| LBC mask/livetime γ | `gamma = 0.47445…` (cell 17) | absent (uses `Texp` only) | absorbed in `Bp / Br` | upstream |
| Pattern detection eff `epspat` | upstream `Background_efficiencies.csv` | inline `[0.908,…]` (stored, unused) | absorbed in `Bp / Br / Spat` | upstream |
| `pto` matrix (e⁻ count → pattern) | upstream `Background_efficiencies.csv` | inline arrays | absorbed in Spat CSVs | upstream |
| `Bp` per pattern | `[141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]` | `[141.4, 4.24e-2, 4.24e-2, 7.7e-6, 2.2e-5, 3.6e-6]` (different background model) | matches Carlos | **Poisson mean for toys** (band tool reads `run.background_Bp`) |
| `Br` per pattern | `[0.039, 0.039, 0.016, 0.052, 0.011, 0.035]` (radioactive Geant4) | folded into `Bkg_p` | matches Carlos | **Poisson mean for toys** (band tool reads `run.background_Br`) |
| Nuisance `θ` | implicit (per-CCD dark-current fits) | none | `[0.5, 10]` multiplicative on `Br` | per-toy fit unchanged |
| Constraint on `θ` | implicit | none | `−θ·Br + 98·N_img·ln(θ·Br/N_img)` | per-toy fit unchanged |
| QCDark signal γ (charge yield) | upstream `df` | upstream `df` | `outputs/qcdark_srdm_gamma_0p4645605` filtered Spat | per-toy fit reads same Spat |
| Exposure | computed per CCD | `Texp = 1.257 kg·day` | `mass_kg × livetime_days = 1.3 kg·day` | per-toy fit unchanged |
| Test statistic | none | `T_mu` Cowan PRC 2011, eq. (14)/(16), no constraint | `q_μ = 2(NLL − NLL_min)` PLR + asymptotic χ² + constraint | per-toy fit unchanged |

### Three subtleties to make explicit

1. **Antoine vs Carlos use different `Bkg` numbers and different likelihoods.**
   Antoine's `Bkg = [141.4, 4.24e-2, 4.24e-2, 7.7e-6, 2.2e-5, 3.6e-6]` is smaller in
   non-(1,1) bins than Carlos's `Bp + Br`, and Antoine's `T_mu` has no nuisance and
   no constraint term — it's a pure Poisson product across patterns with one
   signal-strength parameter. We have already chosen to follow Carlos / pydme for
   the observed-limit calculation (≈8 % average agreement with the pydme reference
   confirms this). The band tool inherits that choice exactly — no new decision is
   introduced. The Antoine notebook is provided primarily as the **methodology**
   template (toy-MC sampling for sensitivity), not for the numerical inputs.

2. **Toy generation is pattern-level Poisson, *not* image-level resimulation.**
   DataSimulator cells 2–4 do per-pixel image simulation with
   `read_noise + dark_current + pattern_finding` to *generate* `Bp`. The band tool
   does **not** redo any of that — it samples directly at pattern level via
   `D^{(t)}_i ~ Poisson(Bp_i + Br_i)`. This is the right thing for a frequentist
   sensitivity band: hold the model fixed and sample the Poisson statistical
   fluctuation of the count given that mean. Sampling at the image level would
   double-count the noise that is *already in* `Bp`.

3. **Toys hold `θ = 1` (nominal).** The Poisson mean per pattern is `Bp + 1·Br`.
   We do **not** sample `θ` from its prior before sampling counts. This is the
   standard frequentist B-only nominal band, and it exactly matches DataSimulator
   cell 21's `Poisson(Bp_i) + Poisson(Br_i)` (which is identically distributed to
   `Poisson(Bp_i + Br_i)` for independent Poisson processes). The constraint term
   on `θ` still operates **inside the per-toy fit** exactly as it does for the
   observed limit; only the toy-data *generation* step ignores `θ`-spread. The
   "B-only post-fit" variant from §1 (toys at `Bp + θ̂·Br`) is a future option;
   the current architecture supports it as a swap of the toy-mean expression.

### What the `band` block must specify — and what it must not

The `band` block in the input config should specify **only sampling-control
parameters**. Anything else would either be redundant (already in `run`) or
double-count physics that is already baked into `Bp / Br / Spat`.

| In `band` block | NOT in `band` block (and why) |
|---|---|
| `scan_binary` (path to existing scan binary) | `background_Bp / Br` — already in `run` |
| `n_toys` | `read_noise`, `dark_current` — upstream of pattern level |
| `rng_seed` | pattern detection efficiencies — already in `Bp / Br / Spat` |
| `n_workers` | livetime, detector mass — already in `experiment` / `detector` |
| `outdir` | constraint parameters — already in `run.constrain_*` |
| (optional) `quantiles`, `save_per_toy_curves`, `tmp_dir`, `keep_per_toy_outputs`, dotted-path overrides for `bp / br` keys | signal models / Spat — already in `srdm_signal_csv` and read by the scan binary |

Bottom line: **the band tool's only added physics is the choice of B-only
nominal Poisson sampling at the pattern level**. Every other knob — including
all noise, efficiency, exposure, and constraint parameters — is unchanged from
the existing scan run and is reused per toy via the same `scan_binary` invocation.

## 5. Files added / changed (summary)

**Added**

- `apps/ccdarksens_band.cc` — the orchestrator binary (single new executable;
  drives Phase 0 / Phase 1 / Phase 2 as in Step 2).
- `apps/README_band.md` — author-facing contract documentation including the
  base contract (1a) and the threshold extension (1b).

**Modified (existing scan configs, Step 4)**

- `configs/scan_srdm_pattern_csv_per_mass_pydme_match.json` — add an optional
  `band` block (with optional `threshold: "asymptotic" | "toy_mc" | "both"`).
  Same pattern for any other analysis that wants a band. **No new config
  files.**

**Modified (per scan binary, Step 1)**

Base contract (Step 1a, all six binaries):
- Parse `run.observed_counts` and use it as the data vector when present.
- Write `upper_limit_sigma_e_mchi_graph` as a `TGraph` of exact (m_χ, σ_UL)
  pairs.
- Uniform ~30-line patch.

Threshold extension (Step 1b, scan-by-scan; opt-in per analysis):
- In scan mode: parse optional `run.q_target_lookup_path`, read TGraph
  `q_target_per_mass`, and use that value at the bisection crossing instead
  of the asymptotic `(Φ⁻¹(cl))²` constant.
- New `run.mode = "threshold_toys"` driver that runs `n_threshold_toys`
  Poisson(s+b) sub-toys per mass at the configured `(m_χ, σ_threshold(m_χ))`,
  computes `q_μ` for each via the existing `MinimizeOverScaleMinuit` /
  `MinimizeOverSigmaAndTheta`, and writes a TGraph `q_target_per_mass` of
  `(m_χ, c_μ)` to `<outdir>/qtarget_threshold.root`.
- Uniform ~60-line patch per binary; needed only for binaries where toy-MC
  bands are desired. Binaries without it still get the asymptotic band.

**Modified (plotter)**

- `apps/ccdarksens_plot_dmelectron_limit.cc` — `--band` flag with autodetection
  of `_asymptotic` and `_toy_mc` graph variants; optional `--band-mode` for
  explicit selection.

**Build**

- `CMakeLists.txt` — new executable `ccdarksens_band`.

**Output ROOT (per analysis)**

`outputs/band/<analysis_label>/band.root`:

| Object | Always | Asymptotic mode | Toy-MC mode |
|---|---|---|---|
| `median_expected_sigma_e_mchi_<mode>` | one per mode | `_asymptotic` | `_toy_mc` |
| `band_1sigma_sigma_e_mchi_<mode>`     | one per mode | `_asymptotic` | `_toy_mc` |
| `band_2sigma_sigma_e_mchi_<mode>`     | one per mode | `_asymptotic` | `_toy_mc` |
| `q_target_per_mass`                   | only in toy_mc / both | — | TGraph |
| `sigma_threshold_per_mass`            | only in toy_mc / both | — | TGraph |
| `meta` (TTree)                        | always       | run provenance | + threshold provenance |
| `band_per_toy` (TTree)                | optional     | branch `sigma_UL_asy` | branch `sigma_UL_toy` |

For `threshold: "both"` runs, all of the above are present in a single ROOT
file, enabling direct overlay plots of the two threshold modes side-by-side.
