"""
Convert Tom Arbaugh's optical spectra .dat files (scalar or SOC) into the
CCDarkSens custom dielectric CSV format expected by entry.py (_load_eps_csv).

Input columns:  energy[eV]  eps1  eps2  sigma1[S/m]  absorption[cm^-1]  reflectivity
Output columns: omega_eV, re_eps_xx, im_eps_xx

Band-gap zeroing is handled downstream by entry.py (band_gap_eV config key),
so we do NOT zero eps2 here — we just pass eps1/eps2 through cleanly.

Usage:
    python3 utils/convert_optical_spectra.py \
        "/path/to/scalar_optical_spectra (1).dat" \
        data/custom_eps/srcd2sb2_eps_scalar.csv

    python3 utils/convert_optical_spectra.py \
        "/path/to/soc_optical_spectra (1).dat" \
        data/custom_eps/srcd2sb2_eps_soc.csv
"""

import sys
import numpy as np

def convert(dat_path: str, out_path: str) -> None:
    energies, eps1s, eps2s = [], [], []
    with open(dat_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                e   = float(parts[0])
                re  = float(parts[1])
                im  = float(parts[2])
            except ValueError:
                continue  # skip nan rows
            if not (np.isfinite(e) and np.isfinite(re) and np.isfinite(im)):
                continue
            energies.append(e)
            eps1s.append(re)
            eps2s.append(im)

    if not energies:
        raise RuntimeError(f"No valid rows found in {dat_path!r}")

    tag = "scalar" if "scalar" in dat_path.lower() else "soc"
    with open(out_path, "w") as f:
        f.write(f"# SrCd2Sb2 dielectric tensor from DFT ({tag} calculation, Tom Arbaugh)\n")
        f.write(f"# Source: {dat_path}\n")
        f.write("# Columns: omega_eV, re_eps_xx, im_eps_xx\n")
        for e, r, i in zip(energies, eps1s, eps2s):
            f.write(f"{e:.15e},{r:.15e},{i:.15e}\n")

    print(f"Wrote {len(energies)} points → {out_path}")
    print(f"  Energy range: {energies[0]:.4f} – {energies[-1]:.4f} eV")
    print(f"  eps2 range:   {min(eps2s):.4f} – {max(eps2s):.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.dat> <output.csv>")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
