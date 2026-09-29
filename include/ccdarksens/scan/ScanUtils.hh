// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ScanUtils.hh -- I declare the small PLR-scan helpers that the scan apps
//  share: the upper-limit crossing on a monotonized q(sigma) grid, log-
//  linear interpolation of the signal in sigma, monotonization of q, and the
//  bracket-and-bisect upper-limit search. They used to be copy-pasted
//  between the scan apps.
// ===========================================================================

#pragma once
/// Shared PLR scan utilities used by scan apps.
///
/// These functions capture the repeated q(σ) monotonicization, UL crossing,
/// signal interpolation, and bracket-and-bisect logic that was copy-pasted
/// between ccdarksens_scan_dmelectron_pattern and ccdarksens_scan_srdm_pattern_csv.

#include <functional>
#include <vector>

namespace ccdarksens::scan {

/// Upper-limit crossing on a monotonized q(σ) grid.
/// Interpolates log-linearly in σ at the first crossing of target_q.
/// Falls back to sigma_list.back() if no crossing is found.
double UlFromQMonoCrossing(const std::vector<double>& sigma_list,
                           const std::vector<double>& q_mono,
                           double target_q);

/// Log-linear interpolation of a pre-computed S_grid at an arbitrary log10(sigma).
/// S_grid[i] is the signal bin vector at sigma_list[i].
/// Clamps to the grid edges when log10_sigma is out of range.
std::vector<double> InterpolateSignal(
    const std::vector<std::vector<double>>& S_grid,
    const std::vector<double>& sigma_list,
    double log10_sigma);

/// Convert raw NLL values to a monotonized q_mu vector (running max, floor at 0).
/// q_mu[k] = max(0, 2*(nll[k] - nll_min)); q_mono is running max over k.
std::vector<double> MonotonizeQ(const std::vector<double>& nll_values,
                                double nll_min);

/// Bracket-and-bisect upper limit in log10(sigma) space.
///
/// Expands hi from log10_seed upward until q_mu_at(hi) >= target_q, then
/// bisects [lo, hi] to find the crossing. Returns log10(sigma_UL).
/// Returns log10_hi if bracketing fails to reach target_q.
///
/// Parameters:
///   q_mu_at     — callable: double(double log10_sigma) → q_mu value (floored at 0)
///   log10_lo    — lower bound of the search range
///   log10_hi    — upper bound of the search range
///   log10_seed  — starting point for bracket expansion (>= log10_lo)
///   target_q    — PLR threshold (e.g. 2.706 for 90% CL asymptotic)
double BisectUpperLimit(std::function<double(double)> q_mu_at,
                        double log10_lo,
                        double log10_hi,
                        double log10_seed,
                        double target_q);

} // namespace ccdarksens::scan
