"""
Shared physical constants for CCDarkSens Python physics modules.
Values are standard; if your original QEdark_constants.py differs,
adjust here as needed.
"""

# Fine-structure constant
alpha = 1/137.035999084

# Electron mass (eV)
me_eV = 510998.95  # 0.51099895 MeV

# Useful conversions
sec_per_day = 86400.0
cm2_per_m2 = 1.0e4
eV_per_keV = 1.0e3

# --- extra QEDark notebook constants ---
amu_kg       = 1.66053906660e-27       # kg
rho_X_eVcm3  = 0.3e9                   # local DM density in eV/cm^3
ccms         = 2.99792458e10           # speed of light [cm/s]
sec_per_year = 365.25 * 86400.0        # s / year
