// ============================================================================
//  CCDarkSens — PCDCalculator
//  Folds n_e/B spectra through P(q|n_e) and builds n_e reconstruction kernels from PCD tables.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once

#include <map>
#include <memory>
#include <vector>

#include <TH1D.h>

namespace ccdarksens {

/**
 * PCDCalculator:
 *
 *  Takes:
 *    - A table of P(q | n_e) histograms (from PCDBasedResponse)
 *    - Signal spectrum   S(n_e)
 *    - Background spectrum B(n_e)
 *
 *  And computes:
 *
 *    λ_q^signal    = Σ_{n_e} S(n_e) * P(q | n_e)
 *    λ_q^background = Σ_{n_e} B(n_e) * P(q | n_e)
 *
 *  The output is a pair of TH1D histograms in q-space, ready for likelihoods.
 */
class PCDCalculator {
public:
  PCDCalculator() = default;

  /**
   * Fold S(n_e) and B(n_e) with P(q | n_e).
   *
   * Inputs:
   *   - pcd_table: map<n_e, TH1D*> from PCDBasedResponse (normalized)
   *   - Sn:        signal spectrum S(n_e), indexed by n_e
   *   - Bn:        background spectrum B(n_e), indexed by n_e
   *
   * Output:
   *   - A pair: { signal_q_hist, background_q_hist } where each is a
   *     TH1D with identical binning to the P(q | n_e) histograms.
   */
  std::pair<std::unique_ptr<TH1D>, std::unique_ptr<TH1D>>
  FoldSpectra(const std::map<int, std::unique_ptr<TH1D>>& pcd_table,
              const std::map<int, double>& Sn,
              const std::map<int, double>& Bn) const;

  // =====================================================================
  // NEW: n_e reconstruction kernel and folding (pattern-space support)
  // =====================================================================

  /// Kernel type: K[i_true][j_obs] = P(n_obs = j | n_true = i)
  using NeKernel = std::vector<std::vector<double>>;

  /**
   * Build an n_e reconstruction kernel P(n_obs | n_true) from the
   * PCD table P(q | n_true).
   *
   *  - ne_min, ne_max: range of n_e indices to consider (both true and obs).
   *  - We map q-bin centers to nearest integer n_obs via std::lround(q),
   *    restrict to [ne_min, ne_max], and accumulate probabilities.
   *  - Each row (fixed n_true) is renormalized to sum to 1 if >0.
   */
  NeKernel BuildNeKernelFromPCD(const std::map<int, std::unique_ptr<TH1D>>& pcd_table,
                              int ne_min, int ne_max,
                              double sigma_res,
                              double Dqmin = 4.4,
                              double Dqmax = 3.125) const;

  /**
   * Fold a true n_e histogram h_true through the given NeKernel.
   *
   *  - h_true: TH1D with bin centers at integer n_e in [ne_min, ne_max]
   *  - kernel: K[n_true][n_obs] as returned by BuildNeKernelFromPCD
   *
   * Returns a new TH1D with the same binning as h_true, containing
   * S_rec(n_obs) = Σ_{n_true} S_true(n_true) * P(n_obs | n_true).
   */
  std::unique_ptr<TH1D> FoldNeSpectrum(TH1D& h_true,
                                       const NeKernel& kernel,
                                       int ne_min, int ne_max,
                                       const std::string& name = "S_rec_ne") const;

private:
  // No member data; all logic is functional.
};

} // namespace ccdarksens
