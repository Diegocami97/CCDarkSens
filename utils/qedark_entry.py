#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
QEDark entry point for DM–electron scattering in Silicon (Skipper CCDs).

This module is intentionally built ONLY from the QEdark source you provided:
  - QEdark_constants.py     (physics constants)
  - DM_halo_dist.py         (SHM / velocity integrals)
  - Si_f2.txt               (material/crystal form-factor data from the repo)
and the equations derived in QEdark_f2.ipynb (to be transcribed into the
'RATE KERNEL' section below).

Public API:
    compute_dRdE(material, mediator, mchi_eV, sigma_e_cm2, halo,
                 band_gap_eV=1.2, eh_pair_eV=3.8, binsize_eV=0.1)
      -> dict(E_eV=array, dRdE_g_day_eV=array)

CLI:
    python utils/qedark_entry.py \
        --material Si --mediator heavy \
        --mchi_MeV 10 --sigma_e_cm2 1e-37 \
        --v0_kms 220 --vE_kms 232 --vesc_kms 544 \
        --band_gap_eV 1.2 --eh_pair_eV 3.8 --binsize_eV 0.1 \
        --out_csv data/qedark_rates/Si/heavy/dRdE_Si_heavy_m10.000000_s1e-37.csv

Units policy:
  - Input mass: eV (CLI accepts MeV then multiplies by 1e6)
  - Cross-section: cm^2
  - Halo velocities: km/s in CLI/JSON → converted to cm/s internally
  - Output differential rate: events / g / day / eV
"""

from __future__ import annotations
import os
import sys
import argparse
import numpy as np

# ---- QEDark source you provided ----
try:
    import QEdark_constants as QEC
except Exception as e:
    raise ImportError("Could not import QEdark_constants.py (make sure it is on PYTHONPATH).") from e

try:
    import DM_halo_dist as HALO
except Exception as e:
    raise ImportError("Could not import DM_halo_dist.py (make sure it is on PYTHONPATH).") from e


# -------------------------------
# Utilities & data-file handling
# -------------------------------

def _qedark_data_dir() -> str:
    """
    Resolve the directory containing material data files (e.g., Si_f2.txt).
    Search order:
      1) env QEDARK_DATA_DIR
      2) repo-relative: ./data/qedark
      3) alongside this file: ../data/qedark
    """
    env = os.environ.get("QEDARK_DATA_DIR", "").strip()
    if env:
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(os.getcwd(), "data", "qedark"),
        os.path.join(here, "..", "data", "qedark"),
        os.path.join(here, "data", "qedark"),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return c
    # fallback to current dir
    return os.getcwd()


def load_material_si_f2(filename: str = "Si_f2.txt") -> dict:
    """
    Load Silicon form-factor data table (the exact file you provided).
    Returns dict with arrays and minimal metadata.

    Expected file format:
      - Plain text table with columns describing the squared crystal form factor
        and any auxiliary grids required by the notebook derivation.
      - Since formats vary across repos, we parse flexibly:
          * ignore leading '#' comments
          * split on whitespace or commas
          * collect numeric columns into a 2D array

    NOTE: The 'RATE KERNEL' below must interpret these columns according
          to the notebook’s convention. We do NOT invent physics here.
    """
    data_path = os.path.join(_qedark_data_dir(), filename)
    if not os.path.isfile(data_path):
        raise FileNotFoundError(
            f"Material file '{filename}' not found in '{_qedark_data_dir()}'. "
            "Set QEDARK_DATA_DIR or place the file under data/qedark/."
        )

    # Load numeric rows (skip comments/blank lines)
    rows = []
    with open(data_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # tolerate CSV or whitespace
            parts = [p for p in s.replace(",", " ").split() if p]
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                # tolerate headers or malformed lines by skipping
                continue

    if not rows:
        raise RuntimeError(f"No numeric data parsed from {data_path}.")

    arr = np.asarray(rows, dtype=float)
    # The exact meaning of columns depends on the notebook definitions.
    # We return a generic dict to be interpreted in the RATE KERNEL.
    return {
        "table": arr,              # raw numerical table (N x M)
        "nrows": arr.shape[0],
        "ncols": arr.shape[1],
        "source_path": data_path,
    }


# -------------------------------
# SHM halo helper(s)
# -------------------------------

def _to_cms(x_kms: float) -> float:
    return float(x_kms) * 1.0e5

def eta_SHM(vmin_cm_s: np.ndarray,
            v0_cm_s: float,
            vE_cm_s: float,
            vesc_cm_s: float) -> np.ndarray:
    """
    Standard Halo Model (Maxwellian truncated at v_esc) velocity integral:
        η(vmin) = ∫_{|v|>vmin}^{v<vesc} f(v + vE) / v  d^3v
    Returns array matching vmin_cm_s.

    For robustness we delegate to DM_halo_dist.py if it exposes a suitable API;
    otherwise we implement a minimal closed-form (erf) version as a fallback.

    The HALO module you provided will be preferred if it has 'eta_shm' or similar.
    """
    # Try to use provided DM_halo_dist.py if a direct function exists
    for name in ("eta_shm", "eta_SHM", "eta"):
        if hasattr(HALO, name):
            fn = getattr(HALO, name)
            return np.asarray(fn(vmin_cm_s, v0_cm_s, vE_cm_s, vesc_cm_s), dtype=float)

    # Fallback: minimal closed-form (dimensionless) using standard expressions.
    # Reference: e.g. Lewin & Smith, arXiv:astro-ph/9603039; appendix for η(vmin)
    from math import sqrt, pi
    from scipy.special import erf

    v0 = v0_cm_s
    vE = vE_cm_s
    vesc = vesc_cm_s

    N = (erf(vesc/v0) - 2/np.sqrt(np.pi)*(vesc/v0)*np.exp(-(vesc/v0)**2))
    x = vmin_cm_s
    a = (x - vE)/v0
    b = (x + vE)/v0
    term = (erf(b) - erf(a)) - (2/np.sqrt(np.pi))*( (np.exp(-a*a) - np.exp(-b*b)) * v0/(2*vE) )
    eta = (1.0 / (2.0*vE*N)) * term
    eta[x > (vesc + vE)] = 0.0
    eta = np.maximum(eta, 0.0)
    return eta


# -------------------------------
# Public API: compute_dRdE
# -------------------------------

_MEDIATOR_TO_INDEX = {
    "heavy": 0, "massive": 0, "0": 0, 0: 0,
    "light": 2, "massless": 2, "2": 2, 2: 2,
}

def compute_dRdE(material: str,
                 mediator: str,
                 mchi_eV: float,
                 sigma_e_cm2: float,
                 halo: dict,
                 band_gap_eV: float = 1.2,
                 eh_pair_eV: float = 3.8,
                 binsize_eV: float = 0.1) -> dict:
    """
    Compute differential rate dR/dE_e for DM–e scattering in Silicon using ONLY
    the QEDark source provided (constants, SHM halo, material form-factor file).

    Parameters
    ----------
    material : str
        Currently 'Si' supported (expects Si_f2.txt).
    mediator : str
        'heavy'/'massive' (n=0) or 'light'/'massless' (n=2).
    mchi_eV : float
        Dark-matter mass in eV.
    sigma_e_cm2 : float
        Electron reference cross section (cm^2).
    halo : dict
        {'v0_kms':..., 'vE_kms':..., 'vesc_kms':...} or same keys in cm/s as
        'v0_cm_s','vE_cm_s','vesc_cm_s'.
    band_gap_eV : float
        Silicon band gap (eV). Below this, rates should be ~0.
    eh_pair_eV : float
        Mean energy to create an e–h pair (eV). Passed for consistency.
    binsize_eV : float
        Desired energy sampling (eV). The material table may dictate a native
        grid; we will align to it or bin to ~binsize_eV.

    Returns
    -------
    dict with:
        'E_eV' : np.ndarray
        'dRdE_g_day_eV' : np.ndarray   (events / g / day / eV)
    """

    # --- Validate & normalize inputs ---
    if material.lower() not in ("si", "silicon"):
        raise NotImplementedError("Only Silicon is supported at the moment (material='Si').")
    nFDM = _MEDIATOR_TO_INDEX.get(mediator, None)
    if nFDM is None:
        raise ValueError("mediator must be one of: heavy/massive/0 or light/massless/2.")

    # Halo velocities to cm/s
    if "v0_cm_s" in halo:
        v0_cm_s = float(halo["v0_cm_s"]); vE_cm_s = float(halo["vE_cm_s"]); vesc_cm_s = float(halo["vesc_cm_s"])
    else:
        v0_cm_s = _to_cms(halo["v0_kms"]); vE_cm_s = _to_cms(halo["vE_kms"]); vesc_cm_s = _to_cms(halo["vesc_kms"])

    # --- Load material table (exact file you provided) ---
    mat = load_material_si_f2("Si_f2.txt")
    tab = mat["table"]    # (N x M) numeric
    # NOTE: The notebook defines what each column means (q-grid, energy bins, |f_crystal|^2, etc.)
    # We'll interpret these in the RATE KERNEL section below.

    # ---- RATE KERNEL (insert notebook equations here) ----
    #
    # The exact steps are dictated by QEdark_f2.ipynb. This block should:
    #   1) Build an energy grid E_eV (≥ band_gap_eV), spaced ~binsize_eV (or native from table).
    #   2) For each energy bin, compute vmin(E) and the velocity integral η(vmin) using eta_SHM().
    #   3) Combine η(vmin), the mediator form factor (nFDM), and the crystal form-factor content
    #      from Si_f2.txt to obtain the differential rate.
    #   4) Normalize to units of events / g / day / eV.
    #
    # IMPORTANT:
    #   - DO NOT invent physics shortcuts here. Transcribe from the notebook.
    #   - Keep the final shape E_eV, dRdE matching 1D arrays.
    #
    # Below is a scaffold to make this drop-in easy and keep the API stable.

    # --- 1) Energy grid ---
    # If the material table directly provides an energy column, prefer it;
    # otherwise define a simple grid starting at band_gap_eV.
    # Here we assume col0 is energy [eV] IF the table is 2+ columns and first makes sense.
    col0 = tab[:, 0]
    if np.all(col0 >= 0) and (np.nanmax(col0) > band_gap_eV):
        # assume an energy-like first column
        E_native = np.asarray(col0, dtype=float)
        E_native = E_native[np.isfinite(E_native)]
        E_native.sort()
        # restrict to physical region above band gap
        E_eV = E_native[E_native >= max(band_gap_eV, 0.0)]
    else:
        # fallback grid (we will replace with notebook-native once we parse Si_f2 format)
        E_max_guess = float(np.nanmax(col0)) if np.isfinite(np.nanmax(col0)) else 300.0
        if not np.isfinite(E_max_guess) or E_max_guess < band_gap_eV + 1.0:
            E_max_guess = 300.0
        nbin = int(np.ceil((E_max_guess - band_gap_eV)/max(binsize_eV, 1e-3)))
        E_eV = band_gap_eV + np.arange(nbin, dtype=float)*binsize_eV

    # --- 2) Velocity integral η(vmin) ---
    # The mapping vmin(E) depends on the scattering kinematics in the notebook.
    # We put a placeholder function to be replaced by the exact expression.
    def vmin_of_E_electron_placeholder(E):
        # TODO: replace with the exact vmin(E) from the notebook derivation.
        # Placeholder returns huge vmin so that η ≈ 0; avoids unphysical numbers before we finalize.
        return np.full_like(E, 1e12, dtype=float)  # cm/s
    vmin = vmin_of_E_electron_placeholder(E_eV)
    eta_vals = eta_SHM(vmin, v0_cm_s, vE_cm_s, vesc_cm_s)

    # --- 3) Crystal form factor & mediator dependence ---
    # Si_f2.txt presumably encodes |f_crystal|^2 over some grid. We need the exact mapping
    # (from the notebook) to combine it with η(vmin) and mediator power nFDM (0 or 2).
    #
    # For now, we create a strictly-positive placeholder shape to keep the plumbing testable.
    # Replace 'shape_factor' with the true F_crystal^2(E, ...).
    shape_factor = np.ones_like(E_eV, dtype=float)

    # Mediator factor: F_DM^2 ∝ q^{-2nFDM}. The q(E) mapping is part of the notebook.
    # We'll absorb it into the 'shape_factor' once we finalize the expression.
    # For now, we keep mediator as a neutral scale (set to 1.0).
    mediator_factor = np.ones_like(E_eV, dtype=float)

    # --- 4) Normalization to events / g / day / eV ---
    # dR/dE ∝ (rho_X / mchi) * sigma_e * η(vmin) * [shape & mediator factors]
    # Constants & numeric prefactors belong to the notebook kernel.
    rho_X_eVcm3 = 0.3e9  # default; could become a parameter if the notebook uses a different value
    prefactor = (rho_X_eVcm3 / float(mchi_eV)) * float(sigma_e_cm2)
    dRdE_native_per_g_per_s_per_eV = prefactor * eta_vals * shape_factor * mediator_factor

    # Convert s → day; (we already made rate per g, not kg)
    dRdE_g_day_eV = dRdE_native_per_g_per_s_per_eV * 86400.0

    # Numerical hygiene
    dRdE_g_day_eV = np.where(np.isfinite(dRdE_g_day_eV), dRdE_g_day_eV, 0.0)
    dRdE_g_day_eV = np.maximum(dRdE_g_day_eV, 0.0)

    # Zero-out below band gap explicitly
    dRdE_g_day_eV[E_eV < band_gap_eV] = 0.0

    return {"E_eV": E_eV, "dRdE_g_day_eV": dRdE_g_day_eV}


# -------------------------------
# CSV writer & CLI
# -------------------------------

CSV_HEADER = [
    "# Differential Rates computed with QEDark entry (ccdarksens)",
    "# material = {material}, mediator = {mediator}, table = {table_path}",
    "# halo (cm/s): v0={v0_cm_s}, vE={vE_cm_s}, vesc={vesc_cm_s}",
    "# mX (eV) = {mchi_eV}",
    "# sigma_e (cm^2) = {sigma_e_cm2}",
    "# Output units: dR/dE in events / g / day / eV",
    "# Columns: E (eV), dRdE (events/g/day/eV)"
]

def _write_csv(path: str, E: np.ndarray, R: np.ndarray, meta: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for line in CSV_HEADER:
            f.write(line.format(**meta) + "\n")
        f.write("E,dRdE\n")
        for e, r in zip(E, R):
            f.write(f"{e:.8g},{r:.10g}\n")

def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--material", default="Si")
    ap.add_argument("--mediator", required=True, choices=["heavy","massive","light","massless"])
    ap.add_argument("--mchi_MeV", type=float, required=True)
    ap.add_argument("--sigma_e_cm2", type=float, required=True)
    ap.add_argument("--v0_kms", type=float, default=220.0)
    ap.add_argument("--vE_kms", type=float, default=232.0)
    ap.add_argument("--vesc_kms", type=float, default=544.0)
    ap.add_argument("--band_gap_eV", type=float, default=1.2)
    ap.add_argument("--eh_pair_eV", type=float, default=3.8)
    ap.add_argument("--binsize_eV", type=float, default=0.1)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    halo = {"v0_kms": args.v0_kms, "vE_kms": args.vE_kms, "vesc_kms": args.vesc_kms}

    out = compute_dRdE(
        material=args.material,
        mediator=args.mediator,
        mchi_eV=args.mchi_MeV * 1.0e6,
        sigma_e_cm2=args.sigma_e_cm2,
        halo=halo,
        band_gap_eV=args.band_gap_eV,
        eh_pair_eV=args.eh_pair_eV,
        binsize_eV=args.binsize_eV,
    )
    E = out["E_eV"]; R = out["dRdE_g_day_eV"]

    meta = dict(
        material=args.material, mediator=args.mediator,
        table_path=os.path.join(_qedark_data_dir(), "Si_f2.txt"),
        v0_cm_s=_to_cms(args.v0_kms), vE_cm_s=_to_cms(args.vE_kms), vesc_cm_s=_to_cms(args.vesc_kms),
        mchi_eV=args.mchi_MeV*1e6, sigma_e_cm2=args.sigma_e_cm2,
    )
    _write_csv(args.out_csv, E, R, meta)
    print(f"[qedark_entry] wrote {args.out_csv}  (N={len(E)})")


if __name__ == "__main__":
    main_cli()
