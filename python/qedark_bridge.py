# ============================================================================
#  CCDarkSens — qedark_bridge
#  Thin adapter between the DAMIC-M collaboration QEDark4DAMIC class and the
#  CCDarkSens ccdarkphys rate-generation pipeline.
#
#  Provides a single compute_dRdE() function with the same signature as the
#  ccdarkphys.qedark.entry module so that both backends are interchangeable.
#
#  Requires: QEDark4DAMIC installed in the active Python environment.
#  See python/ccdarkphys/qedark/entry.py for the self-contained fallback.
#
#  Author: Diego Venegas-Vargas
# ============================================================================
"""
Adapter from QEDark4DAMIC (DAMIC-M collaboration class) to the CCDarkSens
ccdarkphys rate-generation interface.

Usage:
    from qedark_bridge import compute_dRdE
    result = compute_dRdE("Si", "heavy", mchi_eV=1e6, sigma_e_cm2=1e-37,
                           E_eV=E_array, halo=halo_dict)
    # result["dRdE_g_day_eV"]  — events / g / day / eV

The mediator string maps to QEDark nFDM: "heavy"/"massive" → 0, "light"/"massless" → 2.
Halo dict must contain v0_cm_s, vE_cm_s, vesc_cm_s (all in cm/s).
"""
import numpy as np

# Use the canonical class used by the collaboration.
# If your code exposes QEDark4DM instead, just swap the import.
from QEDark4DAMIC import QEDark4DAMIC

_MEDIATOR_TO_nFDM = {
    "heavy": 0, "massive": 0, "0": 0, 0: 0,
    "light": 2, "massless": 2, "2": 2, 2: 2,
}

def compute_dRdE(
    material: str,
    mediator: str,
    mchi_eV: float,
    sigma_e_cm2: float,
    E_eV: np.ndarray,
    halo: dict,
    band_gap: float = 1.2,
    eh_E: float = 3.8,
    binsize_eV: float = 0.1,
    rho_X_eVcm3: float = 0.3e9,
):
    """
    Returns:
      dict with {"E_eV": array, "dRdE_g_day_eV": array}

    Conventions:
      * mediator -> nFDM mapping: heavy/massive -> 0, light/massless -> 2
      * halo expects cm/s: keys v0_cm_s, vE_cm_s, vesc_cm_s
      * output units: events/g/day/eV
      * if the class ctor does not accept xsec_e, we scale linearly by sigma_e_cm2
    """
    nFDM = _MEDIATOR_TO_nFDM.get(mediator, 0)

    v0_cm_s   = halo.get("v0_cm_s")
    vE_cm_s   = halo.get("vE_cm_s")
    vesc_cm_s = halo.get("vesc_cm_s")
    if any(v is None for v in (v0_cm_s, vE_cm_s, vesc_cm_s)):
        raise ValueError("Halo must include v0_cm_s, vE_cm_s, vesc_cm_s (cm/s).")

    # Instantiate collab model; try passing cross-section to ctor if supported.
    supports_xsec_in_ctor = False
    try:
        qed = QEDark4DAMIC(
            rho_X=rho_X_eVcm3,
            v0=v0_cm_s, vE=vE_cm_s, vesc=vesc_cm_s,
            band_gap=band_gap, eh_E=eh_E,
            xsec_e=sigma_e_cm2,
        )
        supports_xsec_in_ctor = True
    except TypeError:
        qed = QEDark4DAMIC(
            rho_X=rho_X_eVcm3,
            v0=v0_cm_s, vE=vE_cm_s, vesc=vesc_cm_s,
            band_gap=band_gap, eh_E=eh_E,
        )

    # The collab API returns native (bin-centered) Ee and differential rates per bin.
    vparams = [v0_cm_s, vE_cm_s, vesc_cm_s]
    Ee_list, dRdE_list = qed.dRdnearray(material, float(mchi_eV), float(binsize_eV), nFDM, 'shm', vparams)

    Ee = np.asarray(Ee_list, dtype=float)
    dRdE_native = np.asarray(dRdE_list, dtype=float)

    # Match collab script convention: events/g/day/eV
    dRdE_g_day_eV = dRdE_native / 365.25 / 1000.0

    if not supports_xsec_in_ctor:
        # If the class did not absorb sigma in ctor, scale linearly here.
        dRdE_g_day_eV *= float(sigma_e_cm2)

    return {"E_eV": Ee, "dRdE_g_day_eV": dRdE_g_day_eV}
