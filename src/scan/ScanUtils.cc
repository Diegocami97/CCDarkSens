// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ScanUtils.cc -- I implement the shared PLR-scan helpers declared in
//  ScanUtils.hh (upper-limit crossing, signal interpolation, q
//  monotonization and the bracket-and-bisect upper-limit search).
// ===========================================================================

#include "ccdarksens/scan/ScanUtils.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ccdarksens::scan {

// ----------------------------------------------------------------------------
// UlFromQMonoCrossing
//   Upper limit from a monotonized q(sigma) grid: at the first grid interval where q
//   crosses target_q I interpolate log-linearly in sigma. If there is no crossing:
//   when q is above the target everywhere I return the smallest sigma; otherwise the
//   first grid sigma at or after the minimum of q that reaches the target, and
//   finally the largest sigma. Empty or mismatched input returns 0 or the last sigma.
// ----------------------------------------------------------------------------
double UlFromQMonoCrossing(const std::vector<double>& sigma_list,
                           const std::vector<double>& q_mono,
                           double target_q)
{
  if (sigma_list.empty() || q_mono.size() != sigma_list.size())
    return sigma_list.empty() ? 0.0 : sigma_list.back();
  const std::size_t Ny = sigma_list.size();
  for (std::size_t k = 0; k + 1 < Ny; ++k) {
    const double q1 = q_mono[k];
    const double q2 = q_mono[k + 1];
    if (q1 < target_q && q2 >= target_q && q2 > q1) {
      const double s1 = sigma_list[k];
      const double s2 = sigma_list[k + 1];
      if (s1 <= 0.0 || s2 <= 0.0) break;
      const double t = (target_q - q1) / (q2 - q1);
      return std::pow(10.0, std::log10(s1) + t * (std::log10(s2) - std::log10(s1)));
    }
  }
  const double q_min = *std::min_element(q_mono.begin(), q_mono.end());
  if (q_min >= target_q) return sigma_list.front();
  const int k_best = static_cast<int>(
      std::min_element(q_mono.begin(), q_mono.end()) - q_mono.begin());
  for (int k = k_best; k < static_cast<int>(Ny); ++k)
    if (q_mono[static_cast<std::size_t>(k)] >= target_q)
      return sigma_list[static_cast<std::size_t>(k)];
  return sigma_list.back();
}

// ----------------------------------------------------------------------------
// InterpolateSignal
//   Signal vector at an arbitrary log10(sigma), linearly interpolated in log10(sigma)
//   between the two neighboring grid points; outside the grid I return the edge vector.
// ----------------------------------------------------------------------------
std::vector<double> InterpolateSignal(
    const std::vector<std::vector<double>>& S_grid,
    const std::vector<double>& sigma_list,
    double log10_sigma)
{
  const std::size_t n = S_grid.empty() ? 0u : S_grid[0].size();
  std::vector<double> out(n, 0.0);
  if (S_grid.size() < 2u) return S_grid.empty() ? out : S_grid[0];
  const double log10_lo = std::log10(sigma_list.front());
  const double log10_hi = std::log10(sigma_list.back());
  const double log10_s  = std::max(log10_lo, std::min(log10_hi, log10_sigma));
  for (std::size_t j = 0; j + 1 < sigma_list.size(); ++j) {
    const double l0 = std::log10(sigma_list[j]);
    const double l1 = std::log10(sigma_list[j + 1]);
    if (log10_s >= l0 && log10_s <= l1) {
      const double t = (l1 - l0) > 1e-300 ? (log10_s - l0) / (l1 - l0) : 0.0;
      for (std::size_t b = 0; b < n; ++b)
        out[b] = S_grid[j][b] + t * (S_grid[j + 1][b] - S_grid[j][b]);
      return out;
    }
  }
  if (log10_s <= log10_lo) return S_grid.front();
  return S_grid.back();
}

// ----------------------------------------------------------------------------
// MonotonizeQ
//   q[k] = max(0, 2*(nll[k] - nll_min)), then a running maximum over k so that q(sigma)
//   never decreases.
// ----------------------------------------------------------------------------
std::vector<double> MonotonizeQ(const std::vector<double>& nll_values,
                                double nll_min)
{
  std::vector<double> q_mono(nll_values.size());
  if (nll_values.empty()) return q_mono;
  const double q0 = 2.0 * (nll_values[0] - nll_min);
  q_mono[0] = (q0 > 0.0) ? q0 : 0.0;
  for (std::size_t k = 1; k < nll_values.size(); ++k) {
    const double q = 2.0 * (nll_values[k] - nll_min);
    q_mono[k] = std::max(q_mono[k - 1], (q > 0.0) ? q : 0.0);
  }
  return q_mono;
}

// ----------------------------------------------------------------------------
// BisectUpperLimit
//   Bracket-and-bisect in log10(sigma): starting at the seed I step upward (step 0.3,
//   doubling, up to 12 times) until q_mu reaches target_q, then bisect (24 iterations at
//   most) down to a tolerance of max(1e-3, 1% of the range) in log10(sigma) or 0.01 in q.
//   Returns log10(sigma_UL), or log10_hi if the target is never reached.
// ----------------------------------------------------------------------------
double BisectUpperLimit(std::function<double(double)> q_mu_at,
                        double log10_lo,
                        double log10_hi,
                        double log10_seed,
                        double target_q)
{
  constexpr double brack_step  = 0.3;  // initial upward step of the bracket [log10 sigma]
  constexpr int    max_expand  = 12;  // maximum number of bracket expansions
  constexpr int    max_iter    = 24;  // maximum number of bisection iterations
  constexpr double q_tol       = 0.01;  // stop when q_mu is within this of target_q
  const double     tol_x       = std::max(1e-3, 0.01 * (log10_hi - log10_lo));

  double lo   = std::max(log10_lo, log10_seed);
  double hi   = lo;
  double step = brack_step;

  for (int i = 0; i < max_expand && hi < log10_hi - 1e-12; ++i) {
    hi = std::min(hi + step, log10_hi);
    if (q_mu_at(hi) >= target_q) break;
    step *= 2.0;
  }

  if (q_mu_at(hi) < target_q) return log10_hi;

  double left = lo, right = hi;
  for (int it = 0; it < max_iter; ++it) {
    const double mid   = 0.5 * (left + right);
    const double q_mid = q_mu_at(mid);
    if (q_mid >= target_q) right = mid;
    else                   left  = mid;
    if (std::abs(right - left) < tol_x || std::abs(q_mid - target_q) < q_tol) break;
  }
  return right;
}

} // namespace ccdarksens::scan
