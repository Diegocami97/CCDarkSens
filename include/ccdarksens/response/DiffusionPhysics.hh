// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  DiffusionPhysics.hh -- I define ComputeSigmaXYUm, the single shared
//  formula for the lateral diffusion width sigma_xy(z, E) = sqrt(-A ln(1 - b
//  z)) * (alpha + beta E_keV). Diffusion and ChargeTransport both use it, so
//  the formula exists in exactly one place.
// ===========================================================================

#pragma once
#include <algorithm>
#include <cmath>

namespace ccdarksens {

/// σ_xy(z,E) = √(-A log(1 - b·z)) · (α + β·E_keV)
///
/// Returns σ_xy in micrometers. Returns 0.0 when z is out of the valid range
/// (i.e. when 1 - b·z ≤ 0), which avoids silent NaN propagation.
// Arguments: depth z [um], energy E [eV], and the diffusion constants A [um^2], b [1/um], alpha, beta [1/keV]. Returns sigma_xy [um].
inline double ComputeSigmaXYUm(double z_um, double E_eV,
                                double A_um2, double b_umInv,
                                double alpha, double beta_per_keV) {
    const double E_keV  = std::max(0.0, E_eV) * 1e-3;
    const double inside = 1.0 - b_umInv * z_um;
    if (inside <= 0.0) return 0.0;
    return std::sqrt(-A_um2 * std::log(inside)) * (alpha + beta_per_keV * E_keV);
}

} // namespace ccdarksens
