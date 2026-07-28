# Regenerating QCDark2 dielectric functions (band gap via scissor)

> **Full workflow (recommended):** [qcdark2_dielectric_workflow.md](qcdark2_dielectric_workflow.md) — step-by-step pipeline, what **8×8×8** means, `Si_comp.h5`, interim vs production configs, **scissor vs DFT** (§13), **low-gap pheno DM study** (§14), and **next steps** (§15).  
> **Ionization coupling (plan):** [band_gap_pheno_ionization.md](band_gap_pheno_ionization.md) — `eh_pair_eV` and scaled `p100K_table.csv` per scissor gap.

CCDarkSens `qcdark2` rates use a precomputed RPA dielectric function ε(ω, **q**) in HDF5 (`epsilon`, `q`, `E`, attrs `M_cell`, `V_cell`, `dE`). The band gap is **not** a free parameter at rate-generation time; it is set when ε is built from DFT.

## Mechanism

QCDark2 applies a **scissor correction** after DFT: all conduction-band energies are shifted so the fundamental gap equals `scissor_bandgap` (eV) in the input file. The dielectric function is then computed from those MO energies.

- DFT (8×8×8 k-grid, basis, xc, …) is stored under `<save_path>/DFT_resources/` and **reused** when only `scissor_bandgap` or `name` changes.
- Each run writes `<save_path>/<name>_resources/epsilon.hdf5`.

Reference: `scissor_bandgap` in `qcdark2/dielectric_pyscf/dft_routines.py` (`convert_to_eV_and_scissor`).

## Shipped `Si_comp.h5`

The bundled composite table (`dielectric_functions/composite/Si_comp.h5`) combines:

- **LFE**, 8×8×8 k-grid, q ≲ 8 (units in file ≈ α·mₑ momentum convention), `q_shift ≈ 0.066286`
- **no LFE**, 8×8×8 k-grid, q up to ~25, `q_shift ≈ 0.01`

Templates in `configs/qcdark2/` mirror these segments. Full composite merge (splice LFE + nolfe q ranges) is not automated here yet; use one segment for band-gap studies or merge manually.

## Quick smoke test

```bash
cd /Users/diegovenegasvargas/Documents/CCDarkSens
python3 utils/qcdark2_regenerate_epsilon.py --input configs/qcdark2/Si_smoke.in
```

Output: `data/qcdark2_epsilon/Si/Si_smoke_nolfe_q1.h5`

## Production-style run (one band gap)

```bash
python3 utils/qcdark2_regenerate_epsilon.py \
  --template configs/qcdark2/Si_lfe8q.in --scissor 1.2
```

Then (optional high-q tail, same scissor):

```bash
python3 utils/qcdark2_regenerate_epsilon.py \
  --template configs/qcdark2/Si_nolfe25q.in --scissor 1.2
```

## Band-gap scan (reuses DFT)

```bash
python3 utils/qcdark2_regenerate_epsilon.py \
  --template configs/qcdark2/Si_lfe8q.in --scissor 1.1 1.2 1.3
```

Point rate grids at the packaged file:

```json
"epsilon_h5": "data/qcdark2_epsilon/Si/Si_lfe8q_gap1.2.h5"
```

## Direct QCDark2 invocation

```bash
cd /Users/diegovenegasvargas/Documents/Software/QCDark2
.venv/bin/python3 -m qcdark2.dielectric_pyscf /path/to/input.in
```

## Expectations

| Stage | Wall time (rough) |
|--------|-------------------|
| Smoke (`4³` k, `q_max=1`, no LFE) | minutes |
| `Si_lfe8q` (8³ k, LFE, q≤8) | hours–days (machine dependent) |
| `Si_nolfe25q` | similar or longer |

Use `mpi = True` in the `.in` file on a cluster (see QCDark2 docs).
