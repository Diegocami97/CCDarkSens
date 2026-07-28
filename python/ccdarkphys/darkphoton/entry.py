"""
DarkELF entry point for hidden-photon (dark photon) absorption.

Wraps ``darkelf.R_absorption()``. Absorption deposits the full hidden-photon
rest energy as a single ionization event — the signal is monochromatic at
E = m_A', not a continuous spectrum like DM-electron scattering. darkelf only
returns the *total* rate (events/kg/yr) at a given mass; there is no native
dR/dE to sample.

To slot into RateTable::MakeTH1D (which interpolates linearly between CSV
points and returns zero outside the tabulated range — see
src/io/RateTable.cc), the line is represented as a boxcar: zero, then a flat
plateau of height (rate / binsize_eV) spanning [mA_eV - binsize_eV/2,
mA_eV + binsize_eV/2], then zero again. Integrating this plateau over a
bin of width binsize_eV reproduces the correct total rate. ``binsize_eV``
must therefore match (or be no coarser than) the bin width the consuming
model.Emin_eV/Emax_eV/nbins ultimately uses downstream, or the line's power
will be mis-binned — same consistency requirement QCDark2's entry.py
enforces for its own dE.

``band_gap_eV`` overrides darkelf's tabulated E_gap for the loaded material
(a plain attribute override — E_gap is only ever used as a scalar threshold
inside darkelf, never to reshape the loaded ELF table).
Pass None to use the material's default gap from its darkelf YAML. There is
no eh_pair_eV parameter here — absorption deposits the full m_A' as a single
line; the DM-electron charge-yield model (eps_h) does not apply.

Custom dielectric path (``eps_csv``)
--------------------------------------
When ``eps_csv`` is provided, darkelf is driven with a user-supplied dielectric
tensor rather than a built-in material file. The CSV must have a header line
starting with ``#`` and columns::

    omega_eV, re_eps_xx, im_eps_xx[, re_eps_yy, im_eps_yy, re_eps_zz, im_eps_zz]

yy and zz columns are optional; when absent they default to the xx values
(isotropic approximation). The polarization-averaged loss function is::

    W_avg(ω) = (1/3) Σ_i  Im[ε_i(ω)] / (Re[ε_i(ω)]² + Im[ε_i(ω)]²)

darkelf is called with an in-darkelf-tree temp file that encodes W_avg as
Re[ε_eff]=0, Im[ε_eff]=1/W_avg — darkelf then computes Im[-1/ε_eff] = W_avg
exactly. A density correction (rho_Si / density_g_cm3) is applied to the
returned rate to account for the custom material density. After the call the
temp files are removed.

``density_g_cm3`` and ``band_gap_eV`` are required when ``eps_csv`` is set.

Public:
    compute_dRdE(material, mA_eV, epsilon, band_gap_eV=None,
                 binsize_eV=0.1, eps_csv=None, density_g_cm3=None,
                 darkelf_dir=None, darkelf_kwargs=None)
      -> dict(E_eV, dRdE_kg_year_eV, meta)

CLI example (built-in material):
    PYTHONPATH=python python3 -m ccdarkphys.darkphoton.entry \
      --material Si --mA_eV 10 --epsilon 1e-13 \
      --darkelf_dir /path/to/DarkELF \
      --out_csv data/darkphoton_rates/Si/dRdE_Si_absorption_m10.000000_e1e-13.csv

CLI example (custom dielectric):
    PYTHONPATH=python python3 -m ccdarkphys.darkphoton.entry \
      --material custom --mA_eV 0.1 --epsilon 1e-13 \
      --eps_csv data/eu5in2sb6_eps.csv --density_g_cm3 6.77 --band_gap_eV 0.06 \
      --darkelf_dir /path/to/DarkELF \
      --out_csv data/darkphoton_rates/Eu5In2Sb6/dRdE_Eu5In2Sb6_absorption_m0.100000_e1e-13.csv
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import uuid

import numpy as np

from ccdarkphys.common import constants as QEC
from ccdarkphys.common import io as CIO

# Density of Si used in darkelf's built-in Si YAML (g/cm³).
# Applied as a correction factor when driving darkelf with a custom ELF:
#   R_custom = R_darkelf_Si × (RHO_SI / density_g_cm3)
_RHO_SI_DARKELF = 2.329

_ENV_DARKELF_DIR = "CCDARK_SENS_DARKELF_DIR"

# darkelf source-file mapping per material. Same convention as
# data_path/rates_dir lookups elsewhere: extend this dict to add materials.
_DARKELF_FILES = {
    # Si using mermin k-grid only (no optical-limit file — featureless Drude-like ELF).
    # eps_electron_optical_filename set to empty string so darkelf never loads the
    # optical file even if Si_eps_electron_opticallimit.dat exists on disk.
    # NOTE: darkelf extrapolates k outside the mermin grid range (k_min~37 eV/c) for
    # m_A' < 37 eV, so rates in the 1-30 eV window are from an extrapolated model.
    "si": dict(target="Si", filename="Si_mermin.dat", phonon_filename="Si_epsphonon_data6K.dat",
               eps_electron_optical_filename="Si_eps_electron_opticallimit_DISABLED.dat"),
    # Si using measured optical constants at k=0, temperature-corrected to 130 K
    # (CCD operating temperature). Source: Edwards/Palik Handbook (1997) [Ref. 86 of
    # DAMIC-M PRL 2025]; Rajkanan et al. (1979) temperature correction [Ref. 87].
    # E_g(130K)=1.146 eV; E1 peak blueshifted to ~3.465 eV, E2 to ~4.335 eV.
    # This is the method used in the DAMIC-M PRL (2025) for the hidden photon limit.
    # Generate file first: python3 utils/build_si_optical_limit.py
    "si_optical": dict(
        target="Si",
        filename="Si_mermin.dat",
        phonon_filename="Si_epsphonon_data6K.dat",
        eps_electron_optical_filename="Si_eps_electron_opticallimit.dat",
    ),
    # Hypothetical narrow-gap material: E_gap=0.34 eV, rho=8 g/cm³, sigma_DC=300 Ω⁻¹cm⁻¹.
    # Two ε₁ brackets — see docs/DarkPhoton_Absorption_HypotheticalMaterial.md.
    # Generate files first: python3 utils/build_darkphoton_hypothetical_material.py --darkelf_dir <path>
    "hypmat_unscreened": dict(
        target="HypMat_0p34",
        filename="",
        phonon_filename="",
        eps_electron_optical_filename="HypMat_0p34_eps_electron_opticallimit_eps1_1.dat",
    ),
    "hypmat_screened": dict(
        target="HypMat_0p34",
        filename="",
        phonon_filename="",
        eps_electron_optical_filename="HypMat_0p34_eps_electron_opticallimit_eps1_12.dat",
    ),
    # QCDark2-informed ELF: Si measured optical dielectric with band_gap_eV set to the
    # scissors-shifted direct gap (2.1 eV, derived from Si_fast_gap0p34.h5 onset).
    # Uses Si's real interband structure + plasmon (~17 eV) as a proxy for SrCd₂Sb₂.
    # This variant is intermediate between the featureless Drude model and the true
    # SrCd₂Sb₂ ELF (which requires AiiDA wavefunction files, bands_workchain_pk=9491).
    # Does NOT replace hypmat_unscreened / hypmat_screened — all three coexist.
    # Generate file first: python3 utils/build_darkphoton_hypmat_qcdark2_elf.py --darkelf_dir <path>
    # band_gap_eV=2.1 is passed at init time (not encoded in the .dat file).
    "hypmat_qcdark2": dict(
        target="HypMat_0p34",
        filename="",
        phonon_filename="",
        eps_electron_optical_filename="HypMat_0p34_qcdark2_eps_electron_opticallimit.dat",
    ),
}


def _load_eps_csv(path: str) -> dict[str, np.ndarray]:
    """
    Read a custom dielectric tensor CSV.

    Required columns : omega_eV, re_eps_xx, im_eps_xx
    Optional columns : re_eps_yy, im_eps_yy, re_eps_zz, im_eps_zz
                       (default to xx values when absent → isotropic)

    Returns dict with keys: omega, re_xx, im_xx, re_yy, im_yy, re_zz, im_zz
    """
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                candidate = stripped.lstrip("# ")
                # Accept only lines whose first token looks like a column name
                # (contains "omega" or "eps" or "eV") — skip free-text comment lines.
                if any(kw in candidate.lower() for kw in ("omega", "eps", "_ev")):
                    parts = candidate.split(",")
                    # Strip any leading label like "Columns: "
                    parts[0] = parts[0].split(":")[-1]
                    header = [h.strip() for h in parts]
                continue
            rows.append([float(x) for x in stripped.split(",")])

    if not rows:
        raise ValueError(f"eps_csv {path!r} contains no data rows.")

    data = np.array(rows)
    if header is None or len(header) < data.shape[1]:
        # fallback: positional
        header = ["omega_eV", "re_eps_xx", "im_eps_xx",
                  "re_eps_yy", "im_eps_yy", "re_eps_zz", "im_eps_zz"][:data.shape[1]]

    col = {name: data[:, i] for i, name in enumerate(header)}

    omega = col["omega_eV"]
    re_xx = col["re_eps_xx"]
    im_xx = col["im_eps_xx"]
    re_yy = col.get("re_eps_yy", re_xx)
    im_yy = col.get("im_eps_yy", im_xx)
    re_zz = col.get("re_eps_zz", re_xx)
    im_zz = col.get("im_eps_zz", im_xx)

    return dict(omega=omega, re_xx=re_xx, im_xx=im_xx,
                re_yy=re_yy, im_yy=im_yy, re_zz=re_zz, im_zz=im_zz)


def _loss_function(re_eps: np.ndarray, im_eps: np.ndarray) -> np.ndarray:
    """Im[-1/ε] = Im[ε] / (Re[ε]² + Im[ε]²).  Returns 0 where |ε|² ≈ 0."""
    denom = re_eps**2 + im_eps**2
    return np.where(denom > 0, im_eps / denom, 0.0)


def _pol_avg_loss(eps_data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Return (omega, W_avg) where W_avg = (1/3)(W_xx + W_yy + W_zz)."""
    w_xx = _loss_function(eps_data["re_xx"], eps_data["im_xx"])
    w_yy = _loss_function(eps_data["re_yy"], eps_data["im_yy"])
    w_zz = _loss_function(eps_data["re_zz"], eps_data["im_zz"])
    return eps_data["omega"], (w_xx + w_yy + w_zz) / 3.0


def _write_custom_darkelf_files(darkelf_dir: str, target: str,
                                omega: np.ndarray, w_avg: np.ndarray,
                                density_g_cm3: float, band_gap_eV: float) -> str:
    """
    Write a minimal darkelf material directory for the custom ELF.

    Encoding trick: Re[ε_eff]=0, Im[ε_eff]=1/W_avg so that darkelf computes
        Im[-1/ε_eff] = Im[ε_eff] / (0² + Im[ε_eff]²) = (1/W_avg)/(1/W_avg²) = W_avg.

    Returns the optical-limit dat filename (basename only, as expected by darkelf).
    """
    data_dir = os.path.join(darkelf_dir, "data", target)
    os.makedirs(data_dir, exist_ok=True)

    # YAML — only rhoT and E_gap matter for R_absorption
    yaml_path = os.path.join(data_dir, f"{target}.yaml")
    yaml_lines = [
        f"rhoT : {density_g_cm3}  # g/cm³ — custom material",
        f"E_gap : {band_gap_eV}  # eV — custom material band gap",
        "e0 : 3.6",
        "A : 100.0",
        "omegap : 10.0",
        "lattice_spacing : 5.0",
        "cLAkms : 5.0",
        "cTAkms : 3.0",
        "ombar : 0.06",
        "Enl_list : []",
        "LOvec : [0.06]",
        "atoms : ['Custom']",
        "unitcell : {'Custom': {'A': 100.0, 'mult': 1}}",
        f"# CCDarkSens custom dielectric: rho={density_g_cm3} g/cm3, E_gap={band_gap_eV} eV",
    ]
    with open(yaml_path, "w") as f:
        f.write("\n".join(yaml_lines) + "\n")

    # Optical-limit dat file: [omega, Re[ε_eff]=0, Im[ε_eff]=1/W_avg]
    dat_name = f"{target}_eps_electron_opticallimit.dat"
    dat_path = os.path.join(data_dir, dat_name)
    with open(dat_path, "w") as f:
        f.write(f"# CCDarkSens custom ELF: W_avg encoded as Re=0, Im=1/W_avg\n")
        # Set below-gap entries to a tiny non-zero value to avoid divide-by-zero
        for w, im_val_raw in zip(omega, 1.0 / np.where(w_avg > 1e-30, w_avg, 1e-30)):
            f.write(f"{w:.6e}  {0.0:.6e}  {im_val_raw:.6e}\n")

    return dat_name


def _resolve_darkelf_dir(explicit: str | None) -> str:
    if explicit:
        return os.path.abspath(os.path.expanduser(explicit))
    env = os.environ.get(_ENV_DARKELF_DIR)
    if env:
        return os.path.abspath(os.path.expanduser(env.strip()))
    raise FileNotFoundError(
        f"darkelf source directory is required. Pass darkelf_dir=... or set {_ENV_DARKELF_DIR}."
    )


def _load_darkelf(darkelf_dir: str):
    if darkelf_dir not in sys.path:
        sys.path.insert(0, darkelf_dir)
    try:
        from darkelf import darkelf as darkelf_cls
    except ImportError as exc:
        raise ImportError(
            f"darkelf package not importable from {darkelf_dir!r}. "
            f"Pass darkelf_dir=... pointing at the darkelf repo root, or set {_ENV_DARKELF_DIR}."
        ) from exc
    return darkelf_cls


def compute_dRdE(
    material: str,
    mA_eV: float,
    epsilon: float,
    band_gap_eV: float | None = None,
    binsize_eV: float = 0.1,
    *,
    eps_csv: str | None = None,
    density_g_cm3: float | None = None,
    darkelf_dir: str | None = None,
    darkelf_kwargs: dict | None = None,
) -> dict:
    """
    Same outward contract as other CCDarkSens backends.

    Returns:
      dict with ``E_eV``, ``dRdE_kg_year_eV`` (events / kg / year / eV), ``meta``.
      The spectrum is a 4-point boxcar centered at mA_eV (see module docstring).

    When ``eps_csv`` is provided (custom dielectric path):
      - ``density_g_cm3`` is required.
      - ``band_gap_eV`` is required (used as absorption threshold).
      - ``material`` is used only as a label; built-in darkelf mappings are ignored.
      - darkelf is still used for the rate formula; a temp material dir is created
        in ``darkelf_dir/data/`` and removed after the call.
    """
    darkelf_path = _resolve_darkelf_dir(darkelf_dir)
    darkelf_cls = _load_darkelf(darkelf_path)

    # ── Custom dielectric path ────────────────────────────────────────────────
    if eps_csv is not None:
        if density_g_cm3 is None:
            raise ValueError("density_g_cm3 is required when eps_csv is set.")
        if band_gap_eV is None:
            raise ValueError("band_gap_eV is required when eps_csv is set.")

        eps_data = _load_eps_csv(eps_csv)
        omega, w_avg = _pol_avg_loss(eps_data)

        # Zero loss below band gap (darkelf's E_gap threshold may not cover the
        # full digitized range; enforce it explicitly in the loss function array).
        w_avg = np.where(omega >= float(band_gap_eV), w_avg, 0.0)

        # Write temp darkelf material dir (unique name to allow parallel runs)
        temp_target = f"CCDarkSens_custom_{uuid.uuid4().hex[:8]}"
        temp_dir = os.path.join(darkelf_path, "data", temp_target)
        try:
            dat_name = _write_custom_darkelf_files(
                darkelf_path, temp_target,
                omega, w_avg,
                float(density_g_cm3), float(band_gap_eV),
            )
            file_kwargs = dict(
                target=temp_target,
                filename="",
                phonon_filename="",
                eps_electron_optical_filename=dat_name,
            )
            obj = darkelf_cls(mX=float(mA_eV),
                              v0kms=238.0, vekms=253.7, vesckms=544.0,
                              **file_kwargs)
            obj.rhoX = float(QEC.rho_X_eVcm3)
            obj.E_gap = float(band_gap_eV)
            rate_kg_yr = float(obj.R_absorption(kappa=float(epsilon)))
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)

        aniso_note = "anisotropic (xx,yy,zz)" if "re_eps_yy" in eps_data else "isotropic (xx only)"
        meta = {
            "material": material,
            "mediator": "absorption",
            "eps_csv": eps_csv,
            "anisotropy": aniso_note,
            "density_g_cm3": float(density_g_cm3),
            "mA_eV": float(mA_eV),
            "epsilon": float(epsilon),
            "binsize_eV": float(binsize_eV),
            "rate_total_kg_yr": rate_kg_yr,
            "band_gap_eV": float(band_gap_eV),
            "rhoX_eVcm3": float(QEC.rho_X_eVcm3),
            "darkelf_dir": darkelf_path,
        }

    # ── Built-in darkelf material path ───────────────────────────────────────
    else:
        key = material.lower()
        if key not in _DARKELF_FILES:
            raise NotImplementedError(
                f"No darkelf file mapping for material={material!r}. "
                f"Pass eps_csv= to use a custom dielectric CSV."
            )

        file_kwargs = dict(_DARKELF_FILES[key])
        if darkelf_kwargs:
            file_kwargs.update(darkelf_kwargs)

        # Halo parameters from Baxter et al. 2021 (EPJC 81, 907)
        obj = darkelf_cls(mX=float(mA_eV),
                          v0kms=238.0, vekms=253.7, vesckms=544.0,
                          **file_kwargs)
        obj.rhoX = float(QEC.rho_X_eVcm3)
        if band_gap_eV is not None:
            obj.E_gap = float(band_gap_eV)
        rate_kg_yr = float(obj.R_absorption(kappa=float(epsilon)))

        meta = {
            "material": material,
            "mediator": "absorption",
            "mA_eV": float(mA_eV),
            "epsilon": float(epsilon),
            "binsize_eV": float(binsize_eV),
            "rate_total_kg_yr": rate_kg_yr,
            "E_gap_darkelf_eV": float(getattr(obj, "E_gap", float("nan"))),
            "rhoX_eVcm3": float(QEC.rho_X_eVcm3),
            "darkelf_dir": darkelf_path,
        }

    # ── Shared boxcar output ──────────────────────────────────────────────────
    half = 0.5 * float(binsize_eV)
    edge = 1e-3 * float(binsize_eV)
    height = rate_kg_yr / float(binsize_eV)
    e_points = np.array(
        [mA_eV - half - edge, mA_eV - half, mA_eV + half, mA_eV + half + edge],
        dtype=float,
    )
    drde_points = np.array([0.0, height, height, 0.0], dtype=float)
    return {
        "E_eV": e_points,
        "dRdE_kg_year_eV": drde_points,
        "meta": meta,
    }


def _header_lines(meta: dict) -> list:
    lines = [
        "# Differential Rates computed with CCDarkSens (DarkELF-Absorption entry)",
        f"# material = {meta['material']}, mediator = {meta['mediator']}",
        f"# mA' (eV) = {meta['mA_eV']}",
        f"# epsilon (kinetic mixing) = {meta['epsilon']}",
        f"# total absorption rate (events/kg/yr) = {meta['rate_total_kg_yr']:.6e}",
        f"# binsize_eV used to represent the monochromatic line = {meta['binsize_eV']}",
        f"# rhoX (eV/cm^3) = {meta['rhoX_eVcm3']}",
        f"# darkelf_dir = {meta['darkelf_dir']}",
    ]
    if "eps_csv" in meta:
        lines += [
            f"# custom eps_csv = {meta['eps_csv']}",
            f"# anisotropy handling = {meta['anisotropy']}",
            f"# density_g_cm3 = {meta['density_g_cm3']}",
            f"# band_gap_eV (threshold) = {meta['band_gap_eV']}",
        ]
    else:
        lines.append(f"# darkelf E_gap (eV) = {meta.get('E_gap_darkelf_eV', 'N/A')}")
    lines += [
        "# Output units: dR/dE in events / kg / year / eV (boxcar line, see module docstring)",
        "# Columns: E (eV), dRdE (events/kg/year/eV)",
    ]
    return lines


def _cli():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--material", default="Si")
    ap.add_argument("--mA_eV", type=float, required=True)
    ap.add_argument("--epsilon", type=float, required=True)
    ap.add_argument("--band_gap_eV", type=float, default=None,
                     help="Override darkelf's tabulated E_gap (or required threshold for --eps_csv).")
    ap.add_argument("--binsize_eV", type=float, default=0.1)
    ap.add_argument("--eps_csv", default=None,
                    help="Path to custom dielectric tensor CSV (see module docstring for format).")
    ap.add_argument("--density_g_cm3", type=float, default=None,
                    help="Material density in g/cm³ — required when --eps_csv is set.")
    ap.add_argument("--darkelf_dir", default=None)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    res = compute_dRdE(
        material=args.material,
        mA_eV=args.mA_eV,
        epsilon=args.epsilon,
        band_gap_eV=args.band_gap_eV,
        binsize_eV=args.binsize_eV,
        eps_csv=args.eps_csv,
        density_g_cm3=args.density_g_cm3,
        darkelf_dir=args.darkelf_dir,
    )
    CIO.write_csv_generic(args.out_csv, res["E_eV"], res["dRdE_kg_year_eV"], _header_lines(res["meta"]))
    print(f"[darkphoton] wrote {args.out_csv}  (N={len(res['E_eV'])}, "
          f"total rate={res['meta']['rate_total_kg_yr']:.4e} events/kg/yr)")


if __name__ == "__main__":
    _cli()
