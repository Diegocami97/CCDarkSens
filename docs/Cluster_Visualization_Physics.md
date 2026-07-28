# Cluster Visualization — Physics and Mathematics

**Document:** `docs/Cluster_Visualization_Physics.md`  
**Script:** `utils/plot_cluster_visualization.py`  
**Author:** Diego Venegas-Vargas  
**Date:** 2026-06-09

---

## Purpose

`plot_cluster_visualization.py` simulates the 2D pixel image that a single
dark-matter–electron recoil leaves on a skipper-CCD detector.  It is a
**self-contained analytical computation** — it does not read any ROOT output
or scan result.  Given a recoil energy $E_r$, a depth $z$, and a material
described by $(E_\text{gap},\,\varepsilon_h)$, it answers:

> *What integer charge pattern does a skipper CCD record for this event?*

The pipeline has six stages:

```
(E_r, z, E_gap, ε_h)
        │
        ▼
 1. Charge yield          ⟨n_e⟩ from ionisation model
        │
        ▼
 2. Lateral diffusion     σ_xy(z, E_r)  from DiffusionPhysics.hh
        │
        ▼
 3. Pixel fractions       f_ij = ∫∫_pixel  G(x,y; σ_xy) dx dy
        │
        ▼
 4. Expected charge       q_ij = ⟨n_e⟩ · f_ij
        │
        ▼
 5. Readout noise         q̃_ij = q_ij + N(0, σ_ro²)
        │
        ▼
 6. Integer image         image_ij = round(q̃_ij)
```

Each stage is described in full detail below.

---

## Stage 1 — Charge yield ⟨n_e⟩

### Physical picture

When a DM particle scatters off an electron, it deposits a recoil energy
$E_r$.  Part of this energy goes into ionising the first electron-hole pair
(overcoming the band gap $E_\text{gap}$); the rest is shared among additional
pairs at an average cost of $\varepsilon_h$ per pair.

### Formula

Convention 1 — *gap energy consumed by the first pair*:

$$
\langle n_e \rangle = \frac{\max(E_r - E_\text{gap},\; 0)}{\varepsilon_h}
$$

This is consistent with Convention 1 in `plot_ne_vs_Er_v2.py`.  The
$\max(\cdot, 0)$ encodes the threshold: no ionisation occurs below the band
gap.

### Evaluated for the four scenarios in the script

| Material | $E_\text{gap}$ [eV] | $\varepsilon_h$ [eV] | $E_r = 4$ eV | $E_r = 2$ eV |
|---|---|---|---|---|
| Si (reference) | 1.2 | 3.8 | $(4.0-1.2)/3.8 = 0.74$ | $(2.0-1.2)/3.8 = 0.21$ |
| Low-gap | 0.1 | 0.5 | $(4.0-0.1)/0.5 = 7.8$ | $(2.0-0.1)/0.5 = 3.8$ |

The Si row at $E_r = 2$ eV gives $\langle n_e \rangle = 0.21$, which is below
the single-carrier detection threshold.  The resulting image is dominated
entirely by readout noise — this is the *intentional* physics message: a 2 eV
event is invisible in silicon but produces a 4-carrier cluster in a lower-gap
material.

---

## Stage 2 — Lateral diffusion width σ_xy

### Physical picture

After ionisation, the free electrons are swept toward the pixel layer by the
applied drift field.  During this drift they undergo transverse Brownian
diffusion — random thermal kicks from lattice scattering — spreading into a
growing lateral cloud.  The width of this cloud at the pixel surface depends
on the depth $z$ at which the interaction occurred and on the local electric
field profile.

### Motivation for the Gaussian ansatz

Each electron executes a random walk in the transverse plane as it drifts
through the bulk.  By the **central limit theorem**, the sum of many
independent random scattering displacements converges to a Gaussian
distribution regardless of the microscopic scattering details.  The transverse
position of a single electron at the pixel layer is therefore well described
by:

$$
p(x, y) = \frac{1}{2\pi\sigma_{xy}^2}\exp\!\left(-\frac{x^2+y^2}{2\sigma_{xy}^2}\right)
$$

where $\sigma_{xy}$ is the RMS lateral displacement accumulated over the full
drift path from depth $z$ to the surface.

### Formula (DiffusionPhysics.hh)

The closed-form expression used by CCDarkSens — and mirrored exactly in the
script — is:

$$
\sigma_{xy}(z,\, E_r) = \sqrt{-A\,\ln\!\left(1 - b\,z\right)}\;\cdot\;
\bigl(\alpha + \beta\, E_r^{(\text{keV})}\bigr)
\quad [\mu\text{m}]
$$

**Parameters** (validated in `CLAUDE.md`, must match scan configs):

| Symbol | Value | Units | Meaning |
|---|---|---|---|
| $A$ | 803.25 | µm² | diffusion amplitude |
| $b$ | $6.5\times10^{-4}$ | µm⁻¹ | field-depth coefficient |
| $\alpha$ | 1.0 | — | dimensionless lateral scale |
| $\beta$ | 0.0 | µm keV⁻¹ | energy dependence of diffusion |

With $\beta = 0$ (current config), the formula is purely depth-dependent:

$$
\sigma_{xy}(z) = \alpha\,\sqrt{-A\,\ln(1 - b\,z)}
$$

### Origin of the $\sqrt{-A\ln(1-bz)}$ form

For a **uniform drift field** $\mathcal{E}$, the variance of the transverse
displacement after drifting a distance $z$ grows linearly:

$$
\sigma_{xy}^2 = 2 D_T \frac{z}{v_d}
$$

where $D_T$ is the transverse diffusion coefficient and $v_d$ is the drift
velocity.  In a non-uniform field $\mathcal{E}(z)$, the local diffusion rate
varies with depth.  Integrating over the actual DAMIC-M field profile and
fitting to a closed form yields the $-A\ln(1-bz)$ argument, which approaches
$A\cdot bz$ (linear in $z$) for small $bz$, recovering the uniform-field limit.

### Evaluated at the three representative depths

| $z$ [µm] | $\sigma_{xy}$ [µm] | $\sigma_{xy}/p$ | Interpretation |
|---|---|---|---|
| 50 (shallow) | $\approx 5.2$ | 0.35 | cloud smaller than one pixel |
| 337 (mid) | $\approx 14.1$ | 0.94 | cloud comparable to pixel pitch |
| 620 (deep) | $\approx 20.4$ | 1.36 | cloud larger than one pixel |

where $p = 15\,\mu\text{m}$ is the pixel pitch.

---

## Stage 3 — Pixel fraction integral

### Setup

The 2D Gaussian charge cloud is centred at the middle of the central pixel of
a $9\times9$ grid.  Pixel $(i,j)$ occupies the area
$[x_i, x_{i+1}]\times[y_j, y_{j+1}]$ with edges at:

$$
x_k = \left(k - 4.5\right)\times p, \qquad k = 0,1,\ldots,9
$$

(and likewise for $y$), so the central pixel straddles $[-p/2,\,+p/2]$ in both
dimensions.

### The integral

The fraction of total charge landing in pixel $(i,j)$ is:

$$
f_{ij} = \int_{x_i}^{x_{i+1}}\int_{y_j}^{y_{j+1}}
\frac{1}{2\pi\sigma_{xy}^2}
\exp\!\left(-\frac{x^2+y^2}{2\sigma_{xy}^2}\right)
dx\,dy
$$

### Factorisation

Because the integrand is separable — it is a product of a function of $x$ only
and a function of $y$ only — the 2D integral splits into a product of two 1D
integrals:

$$
f_{ij} =
\underbrace{\int_{x_i}^{x_{i+1}}
\frac{1}{\sqrt{2\pi}\,\sigma_{xy}}
e^{-x^2/2\sigma_{xy}^2}\,dx}_{I_x}
\;\times\;
\underbrace{\int_{y_j}^{y_{j+1}}
\frac{1}{\sqrt{2\pi}\,\sigma_{xy}}
e^{-y^2/2\sigma_{xy}^2}\,dy}_{I_y}
$$

### Evaluating $I_x$ with the Gaussian CDF

Define the standard Gaussian cumulative distribution function:

$$
\Phi(t) = \frac{1}{\sqrt{2\pi}}\int_{-\infty}^{t} e^{-u^2/2}\,du
= \frac{1}{2}\left[1 + \text{erf}\!\left(\frac{t}{\sqrt{2}}\right)\right]
$$

where the error function is:

$$
\text{erf}(t) = \frac{2}{\sqrt{\pi}}\int_0^t e^{-u^2}\,du
$$

$\Phi(t)$ is the antiderivative of the standard Gaussian density, so:

$$
I_x = \int_{x_i}^{x_{i+1}}
\frac{1}{\sqrt{2\pi}\,\sigma_{xy}}
e^{-x^2/2\sigma_{xy}^2}\,dx
= \Phi\!\left(\frac{x_{i+1}}{\sigma_{xy}}\right)
- \Phi\!\left(\frac{x_i}{\sigma_{xy}}\right)
$$

and identically for $I_y$.  Therefore:

$$
\boxed{
f_{ij}
= \left[\Phi\!\left(\frac{x_{i+1}}{\sigma_{xy}}\right) - \Phi\!\left(\frac{x_i}{\sigma_{xy}}\right)\right]
\times
\left[\Phi\!\left(\frac{y_{j+1}}{\sigma_{xy}}\right) - \Phi\!\left(\frac{y_j}{\sigma_{xy}}\right)\right]
}
$$

This is implemented in `_gauss_pixel_fractions` using Python's `math.erf`:

```python
cdf = 0.5 * (1.0 + np.array([erf(e / (sqrt(2.0) * sigma_um)) for e in edges_um]))
frac_1d = np.diff(cdf)          # differences of CDF at pixel edges
return np.outer(frac_1d, frac_1d)   # outer product → 2D pixel fractions
```

Note: `erf(e / (sqrt(2) * σ))` evaluates $\text{erf}(t/\sqrt{2})$ at
$t = e/\sigma$, which by the identity above equals $2\Phi(e/\sigma) - 1$.

### Normalisation check

By construction, summing over all pixels recovers unity:

$$
\sum_{i,j} f_{ij}
= \left[\Phi(+\infty) - \Phi(-\infty)\right]^2 = 1^2 = 1
$$

In practice, for a $9\times9$ grid with $\sigma_{xy} \lesssim 20\,\mu\text{m}$
and $p = 15\,\mu\text{m}$, the tails outside the window contain
$< 10^{-4}$ of the charge.

### Charge containment as a function of $\sigma_{xy}/p$

| $\sigma_{xy}/p$ | Central pixel $f_{00}$ | $3\times3$ window | Regime |
|---|---|---|---|
| 0.33 | ≈ 97% | ≈ 100% | shallow Si — single-pixel event |
| 0.94 | ≈ 52% | ≈ 96% | mid-depth — charge starts spilling |
| 1.36 | ≈ 34% | ≈ 87% | deep — clear multi-pixel topology |

---

## Stage 4 — Expected charge per pixel

The expected (mean) number of electrons in pixel $(i,j)$ is simply the total
charge times the fraction reaching that pixel:

$$
q_{ij} = \langle n_e \rangle \cdot f_{ij}
$$

This is the noiseless signal image.  It is a real-valued (non-integer) array.

---

## Stage 5 — Readout noise

A skipper CCD measures charge in each pixel independently.  The measurement
adds Gaussian readout noise with RMS $\sigma_\text{ro}$:

$$
\tilde{q}_{ij} = q_{ij} + \eta_{ij}, \qquad \eta_{ij} \overset{\text{i.i.d.}}{\sim} \mathcal{N}(0,\, \sigma_\text{ro}^2)
$$

**Validated parameter** (CLAUDE.md):

$$
\sigma_\text{ro} = 0.16\,e^-
$$

The noise is independent across pixels (i.i.d.) and independent of the signal
charge.  In the script this is:

```python
noisy = expected + rng.normal(0.0, SIGMA_READOUT_E, size=expected.shape)
```

Note that dark current is **not** added here.  This is intentional: the script
visualises cluster *topology* (how signal charge is distributed spatially) not
a full detector image.  For the topology comparison the dark current —
of order $10^{-3}$–$10^{-4}\,e^-$/pixel/day — is negligible compared to
$\sigma_\text{ro}$.

---

## Stage 6 — Integer image (skipper resolution)

A skipper CCD achieves single-electron resolution by reading each pixel many
times and averaging, suppressing readout noise to $\sigma_\text{ro} \ll 0.5\,e^-$.
The resulting measurement is rounded to the nearest integer:

$$
\text{image}_{ij} = \text{round}\!\left(\tilde{q}_{ij}\right) \in \mathbb{Z}
$$

This is the final pixel image shown in the figure.

---

## Full expression for one pixel

Combining all stages, the probability distribution of the integer charge
$k_{ij}$ in pixel $(i,j)$ is:

$$
P(k_{ij} = k) = \int_{k-0.5}^{k+0.5}
\frac{1}{\sqrt{2\pi}\,\sigma_\text{ro}}
\exp\!\left(-\frac{(q - \langle n_e\rangle\, f_{ij})^2}{2\sigma_\text{ro}^2}\right)
dq
$$

$$
= \Phi\!\left(\frac{k + 0.5 - \langle n_e\rangle f_{ij}}{\sigma_\text{ro}}\right)
- \Phi\!\left(\frac{k - 0.5 - \langle n_e\rangle f_{ij}}{\sigma_\text{ro}}\right)
$$

i.e., the probability of observing integer $k$ is the Gaussian probability
mass in the interval $[k-0.5,\,k+0.5]$ centred on the expected charge
$\langle n_e\rangle f_{ij}$.

---

## The four scenarios and their physics message

```
              E_r = 4 eV                    E_r = 2 eV
         ┌────────────────────────┬────────────────────────┐
   Si    │  ⟨n_e⟩ = 0.74         │  ⟨n_e⟩ = 0.21         │
         │  barely 1 carrier      │  BELOW THRESHOLD       │
         │  no cluster topology   │  noise-only image      │
         ├────────────────────────┼────────────────────────┤
  Low-   │  ⟨n_e⟩ = 7.8          │  ⟨n_e⟩ = 3.8          │
   gap   │  rich multi-pixel      │  clear visible cluster │
         │  cluster at mid/deep   │  grows with depth      │
         └────────────────────────┴────────────────────────┘
```

The key comparison is the **right column**: at $E_r = 2$ eV, the Si detector
records nothing (the event is completely masked by readout noise), while the
low-gap detector records a 4-carrier cluster whose topology grows with depth
and is clearly distinguishable from background.  This is the primary
motivation for exploring low-band-gap semiconductor targets for sub-eV
dark matter searches.

---

## Connection to the CCDarkSens C++ framework

| Script element | C++ counterpart |
|---|---|
| `compute_sigma_xy_um(z, E)` | `ccdarksens::ComputeSigmaXYUm` in `include/ccdarksens/response/DiffusionPhysics.hh` |
| `_gauss_pixel_fractions(σ)` | `EfficiencyMC` pixel-level MC (numerical equivalent via sampling) |
| `SIGMA_READOUT_E = 0.16` | `sigma_readout_e` in `efficiency_mc` config block |
| `_scenario_ne(sc)` | `ChargeIonization::FoldToNe` (same convention, no Fano fluctuations here) |
| `PIXEL_SIZE_UM = 15.0` | `pixel_size_um` in `detector` config block |

The script uses deterministic analytic integrals where the C++ framework uses
Monte Carlo sampling, but both evaluate the same underlying physics model.
