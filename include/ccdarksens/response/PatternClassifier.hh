// ============================================================================
//  CCDarkSens — PatternClassifier
//  Header for pattern-label types, classifier configuration thresholds, and row/image scanning interfaces.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once

#include <cstddef>
#include <ostream>
#include <string>
#include <vector>

namespace ccdarksens {

// -----------------------------------------------------------------------------
// Configuration for pattern classification
// -----------------------------------------------------------------------------
struct PatternClassifierConfig {
  // Minimum charge (in electrons) to consider a pixel as "hit"
  // Reference: qmin = 3.75 * read_noise = 3.75 * 0.16 = 0.60 (efficiencies.py)
  double Qmin_e    = 0.6;

  // Maximum charge (in electrons) for neighbor pixels to be considered "empty"
  // (isolation check); reference uses the same threshold as qmin
  double neighbor_Qmax_e = 0.6;

  // Maximum charge (in electrons) we care to encode in integer labels
  double Qmax_e    = 10.0;

  // Residual (readout) noise sigma in electrons for pattern tests
  double sigma_res_e = 0.16;

  // Enable / disable MN and MNL pattern hypotheses
  bool enable_MN   = true;
  bool enable_MNL  = true;

  // Maximum integer electrons per pixel we consider in hypotheses
  int  max_e_per_pixel = 5;

  // Thresholds on the pattern-identification variables
  // (Pm/Pmn/Pmnl: -log(norm.cdf), same as notebook; pass when stat < thr)
  // Reference: thrm=3.5, thrmn=4.0, thrmnl=5.5 (efficiencies.py defaults)
  double thr_M     = 3.5;   // single-pixel compatibility threshold
  double thr_MN    = 4.0;   // 2-pixel patterns
  double thr_MNL   = 5.5;   // 3-pixel patterns

  /// If true, single-pixel branch can return pattern (0) when charge is below 1-e
  /// (notebook identify_pattern_background). Used for background efficiency matrix.
  bool allow_pattern_zero = false;

  /// If true, single-pixel branch uses round(q) clamped to 1..5 instead of
  /// sequential Pm threshold (matches reference Background_efficiencies.csv diagonal).
  bool single_pixel_use_round = false;
};

// -----------------------------------------------------------------------------
// Semantic pattern label: (m), (m,n), (m,n,l), ...
// -----------------------------------------------------------------------------
struct PatternLabel {
  std::vector<int> q;   ///< integer electrons per pixel (e.g. {2}, {3,1}, {3,1,1})
  bool isolated = true; ///< true if pattern is isolated in the row

  bool operator==(const PatternLabel& other) const {
    return (isolated == other.isolated && q == other.q);
  }

  bool operator<(const PatternLabel& other) const {
    if (isolated != other.isolated) {
      return isolated < other.isolated;
    }
    if (q.size() != other.q.size()) {
      return q.size() < other.q.size();
    }
    for (std::size_t i = 0; i < q.size(); ++i) {
      if (q[i] != other.q[i]) return q[i] < other.q[i];
    }
    return false;
  }
};

std::ostream& operator<<(std::ostream& os, const PatternLabel& lab);

// -----------------------------------------------------------------------------
// Result of classifying a cluster / pattern in a row
// -----------------------------------------------------------------------------
struct PatternResult {
  PatternLabel label;      ///< semantic pattern label
  int          start_idx;  ///< index of first pixel of pattern in the row
  int          cluster_size; ///< number of above-threshold pixels in cluster
  double       total_charge_e; ///< sum of charges in involved pixels
  bool         valid;      ///< false if no acceptable pattern was found

  PatternResult()
    : start_idx(-1), cluster_size(0),
      total_charge_e(0.0), valid(false) {}
};

// -----------------------------------------------------------------------------
// PatternClassifier: DAMIC-like pattern ID on 1D rows + local patches
// -----------------------------------------------------------------------------
class PatternClassifier {
public:
  explicit PatternClassifier(const PatternClassifierConfig& cfg);
  ~PatternClassifier() = default;

  const PatternClassifierConfig& config() const { return cfg_; }

  // -------------------------------------------------------------------------
  // High-level API: classify a full 1D row (charges in electrons)
  //
  // Returns a list of patterns (clusters) found along the row, in order.
  // This is the main interface to be used by EfficiencyMC and any analysis
  // that works in "row space" like the Python pattern-efficiency codes.
  // -------------------------------------------------------------------------
  std::vector<PatternResult>
  ScanRow(const std::vector<double>& row) const;

  /// Scan 2D image (3 rows: above, middle, below). Uses middle row only for
  /// pattern finding but requires isolation: pixels in row_above and row_below
  /// (same column range) must be below isolation_threshold (notebook-style).
  /// Returns the single best pattern (max total_charge_e) found in the event, or invalid.
  PatternResult ScanImage2DWithIsolation(
      const std::vector<double>& row_above,
      const std::vector<double>& row_middle,
      const std::vector<double>& row_below,
      double isolation_threshold = -1.0) const;  // <0 => use cfg_.Qmin_e

  /// Scan ALL patterns in the 2D image (3 rows: above, middle, below),
  /// matching the Python scan_image() behaviour: every isolated, qmax-passing
  /// pattern in the middle row is collected and returned in scan order.
  /// Patterns that fail qmax (total_charge_e >= cfg_.Qmax_e) or above/below
  /// isolation are skipped; in-row isolation is handled by classify_from_seed.
  std::vector<PatternResult> ScanAllImage2DWithIsolation(
      const std::vector<double>& row_above,
      const std::vector<double>& row_middle,
      const std::vector<double>& row_below,
      double isolation_threshold = -1.0) const;

  // -------------------------------------------------------------------------
  // Backward-compatible API: classify a local patch (e.g. 3x1)
  //
  // This is a thin wrapper that internally builds a temporary 1D row
  // from the patch and calls ScanRow(), returning the first valid pattern
  // found. Useful while we gradually migrate all callers to ScanRow().
  // -------------------------------------------------------------------------
  PatternResult
  ClassifyLocalPatch(const std::vector<double>& qpix,
                     int NX, int NY) const;

private:
  PatternClassifierConfig cfg_;

  // Helpers: basic predicates
  bool is_hit(double q) const;
  bool is_seed(double q) const { return q >= cfg_.Qmin_e; }
  bool is_neighbor_hit(double q) const { return q >= cfg_.neighbor_Qmax_e; }
  int  approx_int_e(double q) const;

  // Internal: pattern-identification metrics similar to Pmn/Pmnl
  double pattern_stat_M(const std::vector<double>& pix,
                        int m) const;

  double pattern_stat_MN(const std::vector<double>& pix,
                         int m, int n) const;

  double pattern_stat_MNL(const std::vector<double>& pix,
                          int m, int n, int l) const;

  // Try to identify a single-pixel pattern starting at seed_idx
  PatternResult identify_M(const std::vector<double>& row,
                           int seed_idx) const;

  // Try to identify a 2-pixel pattern (MN) starting at seed_idx
  PatternResult identify_MN(const std::vector<double>& row,
                            int seed_idx) const;

  // Try to identify a 3-pixel pattern (MNL) starting at seed_idx
  PatternResult identify_MNL(const std::vector<double>& row,
                             int seed_idx) const;

  // Core logic used by ScanRow(): given a seed pixel index, decide
  // which pattern (M / MN / MNL) best matches, including isolation tests.
  PatternResult classify_from_seed(const std::vector<double>& row,
                                   int seed_idx) const;
};

} // namespace ccdarksens

