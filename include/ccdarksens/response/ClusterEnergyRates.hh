// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterEnergyRates.hh -- Header for folding a dR/dE_true spectrum through
//  a ClusterFitMC kernel into per-E_reco-bin rates.
// ===========================================================================

#pragma once

#include <vector>

#include "ccdarksens/response/ClusterFitMC.hh"

class TH1D;

namespace ccdarksens {

/**
 * Fold a dR/dE_true DENSITY histogram (events/(kg*year*eV), e.g. from
 * RateTable::MakeTH1D) through a ClusterFitMC kernel into rates per E_reco
 * bin -- the continuous-energy analogue of FoldNeToPatternRates
 * (PatternRates.hh).
 *
 * K.Etrue_grid_eV is a SPARSE (typically ~16, often log-spaced) set of
 * points at which the kernel's per-point response K.K[i][:] was actually
 * simulated (see ClusterFitMC::BuildKernel) -- each point stands in for a
 * whole neighborhood of true energies, not just itself. So for each point i:
 *   counts_i = (integral of h_Etrue_density over point i's assigned
 *              neighborhood, i.e. from the geometric midpoint with its
 *              lower neighbor to the geometric midpoint with its upper
 *              neighbor) * exposure_kg_year
 *   rate_reco[j] += sum_i counts_i * K.K[i][j]
 *
 * h_Etrue_density does NOT need the same binning as K.Etrue_grid_eV -- it
 * can be arbitrarily finer (e.g. RateTable::MakeTH1D with many bins), since
 * each point's full neighborhood is genuinely integrated, not point-sampled.
 * FIXED (previously a real bug, see docs/ClusterFitMC_Design.md Sec. 6.8):
 * an earlier version multiplied the point-sampled density by the INPUT
 * histogram's own (fine) bin width instead of the neighborhood width
 * between sparse kernel points, silently discarding the density in between
 * consecutive kernel points -- verified to capture only ~10% of the true
 * total rate on a real spectrum. The neighborhood-integral form above is
 * the correct quadrature for a sparse point-sampled kernel.
 *
 * Unlike the DM-electron channel, there is no separate ionization/charge
 * step before this fold: Phases 1-2 already produce a quenched
 * electron-equivalent spectrum directly, so ClusterFitMC's kernel operates
 * straight on E_true (= E_ee) with no intermediate n_e-equivalent stage.
 */
std::vector<double> FoldEtrueToErecoRates(
    const TH1D& h_Etrue_density,
    double exposure_kg_year,
    const KernelMatrix& K);

}  // namespace ccdarksens
