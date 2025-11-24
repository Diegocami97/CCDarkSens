#!/usr/bin/env python3
# utils/qedark_generate.py
"""
Generate a single QEDark dR/dE CSV using the bridge.

Example:
  python3 utils/qedark_generate.py \
    --module qedark_bridge --function compute_dRdE \
    --material Si --mediator heavy \
    --mchi_MeV 10 --sigma_e_cm2 1e-37 \
    --v0_kms 220 --vE_kms 232 --vesc_kms 544 \
    --emin_eV 1 --emax_eV 300 --estep_eV 0.1 \
    --band_gap_eV 1.2 --eh_pair_eV 3.8 --binsize_eV 0.1 \
    --out_csv data/qedark_rates/Si/heavy/dRdE_Si_heavy_m10.000000_s1e-37.csv
"""
import argparse, importlib, os, sys
import numpy as np

CSV_HEADER = [
    "# Differential Rates computed with QEDark",
    "# v0,vE,vesc = [{v0_cm_s},{vE_cm_s},{vesc_cm_s}] cm/s",
    "# mX = {mchi_eV}",
    "# sigma_e = {sigma_e_cm2}",
    "# rates in units of evts/g/day/eV",
    "# Ee in units of eV"
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--module", required=True)
    ap.add_argument("--function", required=True)
    ap.add_argument("--material", required=True)
    ap.add_argument("--mediator", required=True, choices=["heavy","light","massive","massless"])
    ap.add_argument("--mchi_MeV", type=float, required=True)
    ap.add_argument("--sigma_e_cm2", required=True)  # keep as str for filename fidelity
    ap.add_argument("--v0_kms", type=float, required=True)
    ap.add_argument("--vE_kms", type=float, required=True)
    ap.add_argument("--vesc_kms", type=float, required=True)
    ap.add_argument("--emin_eV", type=float, required=True)
    ap.add_argument("--emax_eV", type=float, required=True)
    ap.add_argument("--estep_eV", type=float, required=True)
    ap.add_argument("--band_gap_eV", type=float, default=1.2)
    ap.add_argument("--eh_pair_eV", type=float, default=3.8)
    ap.add_argument("--binsize_eV", type=float, default=0.1)
    ap.add_argument("--out_csv", required=True)
    args = ap.parse_args()

    # import bridge
    try:
        mod = importlib.import_module(args.module)
    except Exception as e:
        print(f"[qedark_generate] ERROR importing module '{args.module}': {e}", file=sys.stderr)
        sys.exit(2)
    if not hasattr(mod, args.function):
        print(f"[qedark_generate] ERROR: function '{args.function}' not found in '{args.module}'", file=sys.stderr)
        sys.exit(2)
    func = getattr(mod, args.function)

    # Units & grids
    mchi_eV = args.mchi_MeV * 1.0e6
    v0_cm_s   = args.v0_kms  * 1.0e5
    vE_cm_s   = args.vE_kms  * 1.0e5
    vesc_cm_s = args.vesc_kms* 1.0e5

    # Energy grid for metadata (the bridge uses its native binning; we still pass E_eV)
    E = np.arange(args.emin_eV, args.emax_eV + 0.5*args.estep_eV, args.estep_eV)

    # Call bridge
    out = func(
        material=args.material,
        mediator=args.mediator,
        mchi_eV=mchi_eV,
        sigma_e_cm2=float(args.sigma_e_cm2),
        E_eV=E,
        halo=dict(v0_cm_s=v0_cm_s, vE_cm_s=vE_cm_s, vesc_cm_s=vesc_cm_s),
        band_gap=args.band_gap_eV,
        eh_E=args.eh_pair_eV,
        binsize_eV=args.binsize_eV,
    )

    Ee   = out["E_eV"]
    dRdE = out["dRdE_g_day_eV"]
    if Ee.shape != dRdE.shape:
        print("[qedark_generate] ERROR: E and dRdE shapes do not match.", file=sys.stderr)
        sys.exit(5)

    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    with open(args.out_csv, "w") as f:
        for line in CSV_HEADER:
            f.write(line.format(
                v0_cm_s=v0_cm_s, vE_cm_s=vE_cm_s, vesc_cm_s=vesc_cm_s,
                mchi_eV=mchi_eV, sigma_e_cm2=args.sigma_e_cm2
            ) + "\n")
        f.write("E,dRdE\n")
        for e, r in zip(Ee, dRdE):
            f.write(f"{e:.16g},{r:.16g}\n")

    print(f"[qedark_generate] wrote {args.out_csv}")

if __name__ == "__main__":
    main()
