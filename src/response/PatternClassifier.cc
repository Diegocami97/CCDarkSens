// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternClassifier.cc -- Identifies hit patterns in 1D rows or 2D images
//  using Gaussian Pm/Pmn/Pmnl statistics matching the reference notebook
//  logic.
// ===========================================================================

#include "ccdarksens/response/PatternClassifier.hh"

#include <algorithm>
#include <cmath>
#include <limits>

namespace ccdarksens {

namespace {
  // Standard normal CDF: P(X <= x) for X ~ N(mu, sigma^2)
  inline double norm_cdf(double x, double mu, double sigma) {
    if (sigma <= 0.0) return (x >= mu) ? 1.0 : 0.0;
    const double z = (x - mu) / (sigma * std::sqrt(2.0));
    return 0.5 * (1.0 + std::erf(z));
  }
  // -log(p) for p in (0,1], clamped to avoid log(0)
  inline double neg_log_clamped(double p) {
    const double eps = 1e-300;
    if (p <= eps) return -std::log(eps);
    return -std::log(p);
  }
}

// =============================================================================
//  Constructor
// =============================================================================
PatternClassifier::PatternClassifier(const PatternClassifierConfig& cfg)
  : cfg_(cfg)
{}

// =============================================================================
//  Small helpers
// =============================================================================

// pixel above threshold?
bool PatternClassifier::is_hit(double q) const {
  return q >= cfg_.Qmin_e;
}

// Clamp + round q → int electrons (used only for local patch fallback)
int PatternClassifier::approx_int_e(double q) const {
  if (q < cfg_.Qmin_e) return 0;
  int n = static_cast<int>(std::lround(q));
  if (n < 1) n = 1;
  if (n > cfg_.max_e_per_pixel) n = cfg_.max_e_per_pixel;
  return n;
}

// =============================================================================
//  Pattern-statistic metrics: Pm / Pmn / Pmnl (same as notebook)
//  Pm = -log(norm.cdf(q,m,sigma)), Pmn = min over (m,n) perm of -log(cdf*cdf),
//  Pmnl = -log(cdf*cdf*cdf). Lower = better; pass when stat < thr.
// =============================================================================

double PatternClassifier::pattern_stat_M(const std::vector<double>& pix,
                                         int m) const
{
  if (pix.size() < 1) return std::numeric_limits<double>::infinity();

  const double q = pix[0];
  const double sigma = cfg_.sigma_res_e;
  if (sigma <= 0.0) return (q >= m) ? 0.0 : std::numeric_limits<double>::infinity();

  double p = norm_cdf(q, static_cast<double>(m), sigma);
  return neg_log_clamped(p);
}

double PatternClassifier::pattern_stat_MN(const std::vector<double>& pix,
                                          int m, int n) const
{
  if (pix.size() < 2) return std::numeric_limits<double>::infinity();

  const double sigma = cfg_.sigma_res_e;
  const double q0 = pix[0], q1 = pix[1];
  if (sigma <= 0.0) {
    double d1 = std::abs(q0 - m) + std::abs(q1 - n);
    double d2 = std::abs(q0 - n) + std::abs(q1 - m);
    return std::min(d1, d2);
  }

  double p1 = norm_cdf(q0, static_cast<double>(m), sigma) * norm_cdf(q1, static_cast<double>(n), sigma);
  double p2 = norm_cdf(q0, static_cast<double>(n), sigma) * norm_cdf(q1, static_cast<double>(m), sigma);
  return std::min(neg_log_clamped(p1), neg_log_clamped(p2));
}

double PatternClassifier::pattern_stat_MNL(const std::vector<double>& pix,
                                          int m, int n, int l) const
{
  if (pix.size() < 3) return std::numeric_limits<double>::infinity();

  const double sigma = cfg_.sigma_res_e;
  const double q0 = pix[0], q1 = pix[1], q2 = pix[2];
  if (sigma <= 0.0) {
    return std::abs(pix[0] - m) + std::abs(pix[1] - n) + std::abs(pix[2] - l);
  }

  double p = norm_cdf(q0, static_cast<double>(m), sigma)
           * norm_cdf(q1, static_cast<double>(n), sigma)
           * norm_cdf(q2, static_cast<double>(l), sigma);
  return neg_log_clamped(p);
}

// =============================================================================
//  Identify M pattern from a seed
// =============================================================================
PatternResult PatternClassifier::identify_M(const std::vector<double>& row,
                                            int seed_idx) const
{
  PatternResult pr;

  if (seed_idx < 0 || seed_idx >= (int)row.size()) {
    return pr;
  }

  double q = row[seed_idx];
  if (!is_hit(q)) return pr;

  // Optional: round(q) clamped to 1..5 (matches reference data/Background_efficiencies.csv diagonal)
  if (cfg_.single_pixel_use_round) {
    int m = static_cast<int>(std::lround(q));
    if (m < 1) {
      if (cfg_.allow_pattern_zero) {
        pr.valid = true;
        pr.start_idx = seed_idx;
        pr.cluster_size = 1;
        pr.total_charge_e = q;
        pr.label.q = {0};
        bool left_ok  = (seed_idx == 0) || !is_neighbor_hit(row[seed_idx - 1]);
        bool right_ok = (seed_idx == (int)row.size() - 1) || !is_neighbor_hit(row[seed_idx + 1]);
        pr.label.isolated = (left_ok && right_ok);
      }
      return pr;
    }
    if (m > cfg_.max_e_per_pixel) m = cfg_.max_e_per_pixel;
    pr.valid = true;
    pr.start_idx = seed_idx;
    pr.cluster_size = 1;
    pr.total_charge_e = q;
    pr.label.q = {m};
    bool left_ok  = (seed_idx == 0) || !is_neighbor_hit(row[seed_idx - 1]);
    bool right_ok = (seed_idx == (int)row.size() - 1) || !is_neighbor_hit(row[seed_idx + 1]);
    pr.label.isolated = (left_ok && right_ok);
    return pr;
  }

  // Notebook: sequential threshold check. Return [m] when Pm(q,1)<thr,...,Pm(q,m)<thr and Pm(q,m+1)>=thr.
  // So return largest m (1..5) such that pattern_stat_M(q, m) < thr_M; use m=6 to decide [5] vs none.
  if (pattern_stat_M({q}, 1) >= cfg_.thr_M) {
    if (cfg_.allow_pattern_zero) {
      pr.valid = true;
      pr.start_idx = seed_idx;
      pr.cluster_size = 1;
      pr.total_charge_e = q;
      pr.label.q = {0};
      bool left_ok  = (seed_idx == 0) || !is_neighbor_hit(row[seed_idx - 1]);
      bool right_ok = (seed_idx == (int)row.size() - 1) || !is_neighbor_hit(row[seed_idx + 1]);
      pr.label.isolated = (left_ok && right_ok);
    }
    return pr;
  }

  int best_m = 1;
  for (int m = 1; m <= cfg_.max_e_per_pixel; ++m) {
    if (pattern_stat_M({q}, m) >= cfg_.thr_M)
      break;
    best_m = m;
  }
  // Notebook checks Pm(q,6) to decide [5] vs none: if all Pm(1)..Pm(5) < thr and Pm(6) < thr, return none/[0]
  const int m_next = best_m + 1;
  if (best_m == cfg_.max_e_per_pixel && m_next <= 6) {
    if (pattern_stat_M({q}, m_next) < cfg_.thr_M)
      best_m = 0;  // will return invalid or [0] below
  }

  if (best_m == 0) {
    if (cfg_.allow_pattern_zero) {
      pr.valid = true;
      pr.start_idx = seed_idx;
      pr.cluster_size = 1;
      pr.total_charge_e = q;
      pr.label.q = {0};
      bool left_ok  = (seed_idx == 0) || !is_neighbor_hit(row[seed_idx - 1]);
      bool right_ok = (seed_idx == (int)row.size() - 1) || !is_neighbor_hit(row[seed_idx + 1]);
      pr.label.isolated = (left_ok && right_ok);
    }
    return pr;
  }

  pr.valid = true;
  pr.start_idx = seed_idx;
  pr.cluster_size = 1;
  pr.total_charge_e = q;
  pr.label.q = {best_m};

  // isolation: require neighbors below threshold
  // bool left_ok = (seed_idx == 0) || !is_hit(row[seed_idx - 1]);
  // bool right_ok = (seed_idx == (int)row.size() - 1) || !is_hit(row[seed_idx + 1]);
  // pr.label.isolated = (left_ok && right_ok);

  bool left_ok  = (seed_idx == 0) || !is_neighbor_hit(row[seed_idx - 1]);
  bool right_ok = (seed_idx == (int)row.size() - 1) || !is_neighbor_hit(row[seed_idx + 1]);
  pr.label.isolated = (left_ok && right_ok);

  return pr;
}

// =============================================================================
//  Identify MN pattern from a seed
// =============================================================================
PatternResult PatternClassifier::identify_MN(const std::vector<double>& row,
                                             int seed_idx) const
{
  PatternResult pr;
  if (!cfg_.enable_MN) return pr;

  int N = row.size();
  if (seed_idx < 0 || seed_idx >= N) return pr;

  // We assume pixel1 = seed, pixel2 = seed+1
  if (seed_idx + 1 >= N) return pr;
  if (!is_hit(row[seed_idx]) || !is_hit(row[seed_idx+1])) return pr;

  double best_stat = std::numeric_limits<double>::infinity();
  int best_m = 0, best_n = 0;

  for (int m = 1; m <= cfg_.max_e_per_pixel; ++m) {
    for (int n = 1; n <= cfg_.max_e_per_pixel; ++n) {
      double S = pattern_stat_MN({row[seed_idx], row[seed_idx+1]}, m, n);
      if (S < best_stat) {
        best_stat = S;
        best_m = m;
        best_n = n;
      }
    }
  }

  if (best_stat > cfg_.thr_MN) return pr;

  pr.valid = true;
  pr.start_idx = seed_idx;
  pr.cluster_size = 2;
  pr.total_charge_e = row[seed_idx] + row[seed_idx+1];
  pr.label.q = {best_m, best_n};

  // isolation: pixel before and after should be empty
  // bool left_ok = (seed_idx == 0) || !is_hit(row[seed_idx - 1]);
  // bool right_ok = (seed_idx + 2 > N-1) || !is_hit(row[seed_idx + 2]);
  // pr.label.isolated = (left_ok && right_ok);

  bool left_ok  = (seed_idx == 0) || !is_neighbor_hit(row[seed_idx - 1]);
  bool right_ok = (seed_idx + 2 > N-1) || !is_neighbor_hit(row[seed_idx + 2]);
  pr.label.isolated = (left_ok && right_ok);

  return pr;
}

// =============================================================================
//  Identify MNL pattern from a seed
// =============================================================================
PatternResult PatternClassifier::identify_MNL(const std::vector<double>& row,
                                              int seed_idx) const
{
  PatternResult pr;
  if (!cfg_.enable_MNL) return pr;

  int N = row.size();
  if (seed_idx < 0 || seed_idx >= N) return pr;

  // assume three consecutive pixels: seed, seed+1, seed+2
  if (seed_idx + 2 >= N) return pr;
  if (!is_hit(row[seed_idx]) ||
      !is_hit(row[seed_idx+1]) ||
      !is_hit(row[seed_idx+2])) return pr;

  double best_stat = std::numeric_limits<double>::infinity();
  int best_m = 0, best_n = 0, best_l = 0;

  // Notebook: Pmnl = min over all permutations of (m,n,l); same stat formula
  const std::vector<double> pix3 = {row[seed_idx], row[seed_idx+1], row[seed_idx+2]};
  for (int m = 1; m <= cfg_.max_e_per_pixel; ++m) {
    for (int n = 1; n <= cfg_.max_e_per_pixel; ++n) {
      for (int l = 1; l <= cfg_.max_e_per_pixel; ++l) {
        double S = std::min({
          pattern_stat_MNL(pix3, m, n, l),
          pattern_stat_MNL(pix3, m, l, n),
          pattern_stat_MNL(pix3, n, m, l),
          pattern_stat_MNL(pix3, n, l, m),
          pattern_stat_MNL(pix3, l, m, n),
          pattern_stat_MNL(pix3, l, n, m)
        });
        if (S < best_stat) {
          best_stat = S;
          best_m = m;
          best_n = n;
          best_l = l;
        }
      }
    }
  }

  if (best_stat > cfg_.thr_MNL) return pr;

  pr.valid = true;
  pr.start_idx = seed_idx;
  pr.cluster_size = 3;
  pr.total_charge_e = row[seed_idx] + row[seed_idx+1] + row[seed_idx+2];
  pr.label.q = {best_m, best_n, best_l};

  // isolation
  // bool left_ok = (seed_idx == 0) || !is_hit(row[seed_idx - 1]);
  // bool right_ok = (seed_idx + 3 > N-1) || !is_hit(row[seed_idx + 3]);
  // pr.label.isolated = (left_ok && right_ok);
  bool left_ok  = (seed_idx == 0) || !is_neighbor_hit(row[seed_idx - 1]);
  bool right_ok = (seed_idx + 3 > N-1) || !is_neighbor_hit(row[seed_idx + 3]);
  pr.label.isolated = (left_ok && right_ok);    
  return pr;
}

// =============================================================================
//  Decide what pattern appears at seed_idx
//
//  This is a faithful C++ port of the Python identify_pattern() sequential
//  decision tree (efficiencies.py).  The key property of that tree is:
//
//    Each nested "if Pmn(m,1) < thr2" call tests whether the two-pixel charge
//    is STILL compatible with a HIGHER hypothesis (m+1 total electrons).
//    Pattern identity is determined by the HIGHEST m for which the test passes,
//    NOT by which (m,n) minimises the statistic.
//
//  The original argmin implementation incorrectly returned (1,1) for events
//  that should be (2,1), (3,1), etc.
// =============================================================================
PatternResult PatternClassifier::classify_from_seed(const std::vector<double>& row,
                                                    int seed_idx) const
{
  const int N = static_cast<int>(row.size());
  if (seed_idx < 0 || seed_idx >= N || !is_seed(row[seed_idx]))
    return PatternResult{};

  const double thr  = cfg_.Qmin_e;
  const double thr2 = cfg_.thr_MN;
  const double thr3 = cfg_.thr_MNL;

  // Python: cl[cl < thr] = 0  (pixels below threshold contribute as 0 to stats)
  auto zbt = [thr](double q) -> double { return (q >= thr) ? q : 0.0; };

  const double cl1 = zbt(row[seed_idx]);
  const double cl2 = (seed_idx + 1 < N) ? zbt(row[seed_idx + 1]) : 0.0;
  const double cl3 = (seed_idx + 2 < N) ? zbt(row[seed_idx + 2]) : 0.0;

  // Raw (un-zeroed) neighbour values — used only for isolation checks.
  const double raw_left = (seed_idx - 1 >= 0) ? row[seed_idx - 1] : 0.0;
  const double raw_r2   = (seed_idx + 2 < N)  ? row[seed_idx + 2] : 0.0;
  const double raw_r3   = (seed_idx + 3 < N)  ? row[seed_idx + 3] : 0.0;

  // Isolation flags (Python: cl[0]<thr, cl[3]<thr, cl[4]<thr)
  const bool left_ok   = (raw_left < thr);
  const bool right2_ok = (raw_r2   < thr);   // right neighbour of 2-pixel span
  const bool right3_ok = (raw_r3   < thr);   // right neighbour of 3-pixel span

  // PMN: min(-log prod) over all permutations of (m,n) for pixels (cl1,cl2).
  // Returns +inf if enable_MN is false.
  auto PMN = [&](int m, int n) -> double {
    if (!cfg_.enable_MN) return std::numeric_limits<double>::infinity();
    return pattern_stat_MN({cl1, cl2}, m, n);
  };

  // PMNL: min(-log prod) over all 6 permutations of (m,n,l) for (cl1,cl2,cl3).
  // Returns +inf if enable_MNL is false.
  auto PMNL = [&](int m, int n, int l) -> double {
    if (!cfg_.enable_MNL) return std::numeric_limits<double>::infinity();
    int p[3] = {m, n, l};
    std::sort(p, p + 3);
    double best = std::numeric_limits<double>::infinity();
    do {
      double s = pattern_stat_MNL({cl1, cl2, cl3}, p[0], p[1], p[2]);
      if (s < best) best = s;
    } while (std::next_permutation(p, p + 3));
    return best;
  };

  const PatternResult kInvalid{};

  // Build a valid 2-pixel result; returns invalid if not isolated.
  auto make2 = [&](int m, int n) -> PatternResult {
    if (!left_ok || !right2_ok) return kInvalid;
    PatternResult r;
    r.valid           = true;
    r.start_idx       = seed_idx;
    r.cluster_size    = 2;
    r.total_charge_e  = row[seed_idx] + ((seed_idx+1 < N) ? row[seed_idx+1] : 0.0);
    r.label.q         = {m, n};
    r.label.isolated  = true;
    return r;
  };

  // Build a valid 3-pixel result; returns invalid if not isolated.
  auto make3 = [&](int m, int n, int l) -> PatternResult {
    if (!left_ok || !right3_ok) return kInvalid;
    PatternResult r;
    r.valid           = true;
    r.start_idx       = seed_idx;
    r.cluster_size    = 3;
    r.total_charge_e  = row[seed_idx]
                      + ((seed_idx+1 < N) ? row[seed_idx+1] : 0.0)
                      + ((seed_idx+2 < N) ? row[seed_idx+2] : 0.0);
    r.label.q         = {m, n, l};
    r.label.isolated  = true;
    return r;
  };

  // =========================================================================
  // Python identify_pattern() decision tree — translated verbatim.
  //
  // The tree is entered at the top-level Pmn(1,1) check and descends by
  // testing progressively higher-charge hypotheses.  Each level that passes
  // means the charge is consistent with more electrons, narrowing the label.
  // =========================================================================

  if (PMN(1, 1) < thr2) {
    // ---- 2-pixel branch ---------------------------------------------------
    if (PMN(2, 1) < thr2) {
      if (PMN(3, 1) < thr2) {
        if (PMN(4, 1) < thr2) {
          if (PMN(5, 1) < thr2) {
            return kInvalid;            // total charge too high
          }
          return make2(4, 1);           // [4,1]
        }
        // Pmn(4,1) fails
        if (PMN(3, 2) < thr2)
          return make2(3, 2);           // [3,2]
        if (PMNL(3, 1, 1) < thr3) {
          if (PMNL(3, 1, 2) < thr3)
            return kInvalid;
          return make3(3, 1, 1);        // [3,1,1]
        }
        return make2(3, 1);             // [3,1]

      }
      // Pmn(3,1) fails
      if (PMN(2, 2) < thr2) {
        if (PMNL(2, 2, 1) < thr3) {
          if (PMNL(2, 1, 3) < thr3) return kInvalid;
          if (PMNL(2, 2, 2) < thr3) return kInvalid;
          return make3(2, 2, 1);        // [2,2,1]
        }
        if (PMN(2, 3) < thr2) {
          if (PMNL(2, 3, 1) < thr3) return kInvalid;
          if (PMN(2, 4) < thr2)     return kInvalid;
          if (PMN(3, 3) < thr2)     return kInvalid;
          return make2(2, 3);           // [2,3]  (rare in practice)
        }
        return make2(2, 2);             // [2,2]
      }
      if (PMNL(2, 1, 1) < thr3) {
        if (PMNL(2, 1, 2) < thr3) {
          if (PMNL(2, 1, 3) < thr3) return kInvalid;
          return make3(2, 2, 1);        // [2,2,1]  (via 2,1,2 route)
        }
        return make3(2, 1, 1);          // [2,1,1]
      }
      return make2(2, 1);               // [2,1]

    }
    // Pmn(2,1) fails — (1,1) or 3-pixel upgrade
    if (PMNL(1, 1, 1) < thr3) {
      if (PMNL(1, 1, 2) < thr3) {
        if (PMNL(1, 1, 3) < thr3) {
          if (PMNL(1, 1, 4) < thr3) return kInvalid;
          return make3(3, 1, 1);        // [3,1,1]
        }
        return make3(2, 1, 1);          // [2,1,1]
      }
      return make3(1, 1, 1);            // [1,1,1]
    }
    return make2(1, 1);                 // [1,1]

  }

  // ---- Single-pixel branch: cl[0]<thr AND cl[2]<thr ----------------------
  // cl[2] is the pixel right of the seed (must be empty)
  {
    const bool right1_ok = (seed_idx + 1 >= N) || (row[seed_idx + 1] < thr);
    if (!left_ok || !right1_ok) return kInvalid;
    return identify_M(row, seed_idx);
  }
}

// =============================================================================
//  Main row-scanner
// =============================================================================
std::vector<PatternResult>
PatternClassifier::ScanRow(const std::vector<double>& row) const
{
  std::vector<PatternResult> results;
  int N = row.size();
  int i = 0;

  while (i < N) {
    // Find next seed pixel
    while (i < N && !is_hit(row[i])) {
      i++;
    }
    if (i >= N) break;

    // classify from this seed
    PatternResult pr = classify_from_seed(row, i);

    if (pr.valid) {
      results.push_back(pr);
      i = pr.start_idx + pr.cluster_size; // jump past cluster
    } else {
      i++;  // failed to classify, move on
    }
  }

  return results;
}

// =============================================================================
//  Scan 2D image with isolation (notebook-style: middle row + above/below empty)
// =============================================================================
PatternResult PatternClassifier::ScanImage2DWithIsolation(
    const std::vector<double>& row_above,
    const std::vector<double>& row_middle,
    const std::vector<double>& row_below,
    double isolation_threshold) const
{
  const int n = static_cast<int>(row_middle.size());
  if (n < 5 || (int)row_above.size() != n || (int)row_below.size() != n)
    return PatternResult();

  const double thr = (isolation_threshold >= 0.0) ? isolation_threshold : cfg_.Qmin_e;

  // Match notebook scan_image / scan_image_background: classify first, then check isolation
  // only for the columns the pattern actually spans (1 col for 1-pixel, 2 for 2-pixel, 3 for 3-pixel).
  for (int i = 0; i <= n - 4; ++i) {
    if (!is_seed(row_middle[static_cast<std::size_t>(i)])) continue;

    // Build segment (left_neighbor, row[i], row[i+1], row[i+2], row[i+3]) = 5 elements
    std::vector<double> seg(5);
    seg[0] = (i == 0) ? 0.0 : row_middle[static_cast<std::size_t>(i - 1)];
    for (int j = 0; j < 4; ++j)
      seg[static_cast<std::size_t>(j + 1)] = row_middle[static_cast<std::size_t>(i + j)];
    PatternResult pr = classify_from_seed(seg, 1);
    if (!pr.valid) continue;

    // Isolation: above/below below threshold only for pattern span (notebook: 1 col for [m], 2 for [m,n], 3 for [m,n,l])
    const int span = pr.cluster_size;
    bool isolated = true;
    for (int j = 0; j < span && isolated; ++j) {
      const std::size_t k = static_cast<std::size_t>(i + j);
      if (row_above[k] >= thr || row_below[k] >= thr) isolated = false;
    }
    if (!isolated) continue;

    // Notebook: for single-pixel reject if charge >= qmax (they use image_ext[i+1] < self.qmax)
    if (pr.label.q.size() == 1u && pr.label.q[0] > 0 &&
        pr.total_charge_e >= cfg_.Qmax_e)
      continue;

    return pr;
  }
  return PatternResult();
}

// =============================================================================
//  Scan ALL patterns in a 2D image (Python scan_image() equivalent)
//
//  Differences from ScanImage2DWithIsolation:
//    - Collects every isolated, qmax-passing pattern (not just the first)
//    - Applies qmax to ALL pattern sizes (1-, 2-, 3-pixel)
//    - Advances past each accepted pattern span (like Python i += len(var))
//  This matches the reference CSV normalization: efficiency = count / Nsims
//  can exceed 1 when multiple patterns appear in one image.
// =============================================================================
std::vector<PatternResult> PatternClassifier::ScanAllImage2DWithIsolation(
    const std::vector<double>& row_above,
    const std::vector<double>& row_middle,
    const std::vector<double>& row_below,
    double isolation_threshold) const
{
  const int ncols = static_cast<int>(row_middle.size());
  if (ncols < 5 ||
      static_cast<int>(row_above.size())  != ncols ||
      static_cast<int>(row_below.size())  != ncols)
    return {};

  const double thr = (isolation_threshold >= 0.0) ? isolation_threshold : cfg_.Qmin_e;
  std::vector<PatternResult> results;

  int i = 0;
  while (i < ncols) {
    if (row_middle[static_cast<std::size_t>(i)] < thr) { ++i; continue; }

    // Build seg[0..4]: Python image_ext[i:i+5]
    //   seg[0] = left neighbour  (cl[0])
    //   seg[1] = row_middle[i]   (cl[1] = seed)
    //   seg[2] = row_middle[i+1] (cl[2])
    //   seg[3] = row_middle[i+2] (cl[3])
    //   seg[4] = row_middle[i+3] (cl[4])
    std::vector<double> seg(5, 0.0);
    seg[0] = (i > 0)          ? row_middle[static_cast<std::size_t>(i - 1)] : 0.0;
    seg[1] = row_middle[static_cast<std::size_t>(i)];
    seg[2] = (i + 1 < ncols)  ? row_middle[static_cast<std::size_t>(i + 1)] : 0.0;
    seg[3] = (i + 2 < ncols)  ? row_middle[static_cast<std::size_t>(i + 2)] : 0.0;
    seg[4] = (i + 3 < ncols)  ? row_middle[static_cast<std::size_t>(i + 3)] : 0.0;

    PatternResult pr = classify_from_seed(seg, 1);

    if (!pr.valid) { ++i; continue; }

    const int span = pr.cluster_size;

    // qmax check: Python single→q<qmax, multi→sum<qmax; unified as total_charge_e
    if (pr.total_charge_e >= cfg_.Qmax_e) { ++i; continue; }

    // Above/below isolation: every column in the pattern span must be below thr
    // in both row_above and row_below  (Python: image_ext_up/down check)
    bool abv_ok = true;
    for (int j = 0; j < span && abv_ok; ++j) {
      const int col = i + j;
      if (col < ncols) {
        const auto k = static_cast<std::size_t>(col);
        if (row_above[k] >= thr || row_below[k] >= thr)
          abv_ok = false;
      }
    }
    if (!abv_ok) { ++i; continue; }

    // Pattern accepted — fix start_idx to absolute column position
    pr.start_idx = i;
    results.push_back(pr);
    i += span;   // Python: i += len(var)
  }
  return results;
}

// =============================================================================
//  Backward-compatibility: classify a local patch (e.g. 3×1)
// =============================================================================
PatternResult
PatternClassifier::ClassifyLocalPatch(const std::vector<double>& qpix,
                                      int NX, int NY) const
{
  // if NY > 1, flatten row-major, ignore y structure
  if (NX <= 0 || NY <= 0) return PatternResult();

  std::vector<double> row;
  row.reserve(NX*NY);

  for (int iy = 0; iy < NY; ++iy) {
    for (int ix = 0; ix < NX; ++ix) {
      int idx = iy*NX + ix;
      if (idx < (int)qpix.size())
        row.push_back(qpix[idx]);
    }
  }

  auto vec = ScanRow(row);
  if (!vec.empty()) return vec[0];
  return PatternResult();
}

// =============================================================================
//  ostream operator<< for PatternLabel
// =============================================================================
std::ostream& operator<<(std::ostream& os, const PatternLabel& lab) {
  os << "(";
  for (size_t i = 0; i < lab.q.size(); ++i) {
    os << lab.q[i];
    if (i + 1 < lab.q.size()) os << ",";
  }
  os << ")";
  if (!lab.isolated) os << " [non-isolated]";
  return os;
}

} // namespace ccdarksens
