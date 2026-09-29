// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  EfficiencyMC.hh -- Header for EfficiencyMC configuration, pattern-table
//  construction, and ε(n_e) application to spectra.
// ===========================================================================

#pragma once

#include <map>
#include <memory>
#include <ostream>
#include <vector>

#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/PixelSimulator.hh"
#include "ccdarksens/response/PatternClassifier.hh"

class TH1D;

namespace ccdarksens {

class PatternImageGenerator;

/**
 * EfficiencyMC configuration:
 *
 *  - Number of MC trials per n_e
 *  - Row-segment length for pattern simulation (1D)
 *  - PixelSimulator settings (diffusion, noise, dark current)
 *  - Accepted pattern labels that define the efficiency ε(n_e)
 *
 * EfficiencyMC will:
 *   1) Simulate a 1D row segment with ne_true electrons at Ee_eV.
 *   2) Use PatternClassifier::ScanRow(row) to find patterns.
 *   3) Build P(PatternLabel | n_e).
 *   4) Compute ε(n_e) by summing prob. over accepted_labels.
 */
struct EfficiencyMCConfig {
  /// Trials for each n_e
  int ne_trials = 50000;          // 50k MC events per n_e

  /// Length of the 1D row segment in pixels used in the MC
  /// (must be consistent with pix_cfg.nx and ny=1).
  int row_length = 32;            // 32-pixel rows

  /// Pixel simulator configuration (should typically have
  ///  mode = PixelSimMode::RowSegment,
  ///  nx   = row_length,
  ///  ny   = 1).
  PixelSimulatorConfig pix_cfg;

  /// Accepted patterns used to define epsilon(n_e).
  /// For example: { {1}, {2}, {1,1}, {2,1} } etc.
  std::vector<PatternLabel> accepted_labels;

  /// Random seed (if needed for additional internal RNG)
  unsigned long seed = 12345;
  /// Diagnostic only — not used in production scans. When true, Poisson DC is
  /// overlaid on the 1D row segment before classification. For DC-dependent
  /// efficiency studies, use the 2D image path (PatternImageGenerator::lambda_dc)
  /// which models DC pileup in the correct 2D spatial context.
  bool include_dc_pileup = false;
};

/// Decode pattern ID to digit list (e.g. 11 → [1,1], 211 → [2,1,1]).
inline std::vector<int> DecodePatternCode(int code) {
  if (code <= 0) return {1};
  std::vector<int> digits;
  while (code) {
    digits.insert(digits.begin(), code % 10);
    code /= 10;
  }
  return digits;
}

/**
 * EfficiencyMC computes:
 *
 *   1) Full pattern probability table:
 *        pattern_table_[n_e][label] = P(label | n_e)
 *
 *   2) Epsilon(n_e) = sum_{label in accepted_labels} P(label | n_e)
 *
 * using ChargeTransport + PixelSimulator + PatternClassifier (row-based),
 * or optionally via 2D image simulation (PatternImageGenerator path).
 */
class EfficiencyMC {
public:
  // Constructor: settings, the shared charge-transport model and the pattern classifier.
  EfficiencyMC(const EfficiencyMCConfig& cfg,
            std::shared_ptr<ChargeTransport> ct,
            std::shared_ptr<PatternClassifier> classifier);

  /// Compute and return epsilon(n_e) in the range [ne_min, ne_max]
  /// at a reference energy Ee_eV.
  ///
  /// This will:
  ///   - Build the internal pattern_table_ via MC,
  ///   - Sum over cfg_.accepted_labels to fill a TH1D with ε(n_e).
  std::unique_ptr<TH1D> PrecomputeEpsilon(int ne_min, int ne_max,
                                          double Ee_eV = 0.0);

  /// Access the full pattern table:
  ///   pattern_table_[n_e][label] = normalized P(label | n_e).
  const std::map<int, std::map<PatternLabel, double>>&
  GetPatternTable() const {
    return pattern_table_;
  }

  /// Compute ε(E, n_e) on an energy grid.
  ///  - E_grid: list of reference energies (eV)
  ///  - ne_min, ne_max: electron range
  /// Fills internal tables:
  ///  - energy_grid_eV_[iE]
  ///  - epsilon_Ene_[iE][ne - ne_min]
  void PrecomputeEpsilonVsEnergy(const std::vector<double>& E_grid,
                                 int ne_min, int ne_max);

  /// Convenience: is there at least one accepted label whose total charge == ne?
  bool IsAcceptedPattern(int ne) const;

  /// Direct access to config (read-only)
  const EfficiencyMCConfig& config() const { return cfg_; }

  /// Run n_examples MC trials for the given n_e at Ee_eV and print the
  /// pixel row charges and the classified pattern to out (for debugging;
  /// notebook-style "how the image looks" for one case).
  void PrintExampleRows(int ne, double Ee_eV, int n_examples, std::ostream& out);

  /// When set, BuildPatternTable uses 2D image simulation (notebook-style:
  /// generate_image_E + scan with isolation) instead of 1D row. Call before
  /// PrecomputeEpsilon / BuildPatternTable.
  void SetPatternImageGenerator(std::shared_ptr<PatternImageGenerator> gen);

  /// E grid and ε(E,n_e) accessors (used by apps / pipeline)
  const std::vector<double>& energy_grid_eV() const { return energy_grid_eV_; }
  const std::vector<std::vector<double>>& epsilon_Ene() const { return epsilon_Ene_; }
  std::unique_ptr<TH1D>
  PrecomputeEpsilonWithPatternEff(int ne_min, int ne_max, double Ee_eV,
                                  const std::map<std::pair<int,int>, double>& pattern_eff);

private:
  EfficiencyMCConfig cfg_;  // settings
  std::shared_ptr<ChargeTransport>   ct_;  // depth / diffusion sampler
  std::shared_ptr<PatternClassifier> classifier_;  // turns pixel charges into pattern labels
  std::shared_ptr<PatternImageGenerator> img_gen_;  ///< optional: 2D image path for efficiency

  /// Probability table: P(pattern label | n_e).
  std::map<int, std::map<PatternLabel, double>> pattern_table_;

  /// Workhorse: build pattern_table_ for ne in [ne_min, ne_max] at Ee_eV.
  void BuildPatternTable(int ne_min, int ne_max, double Ee_eV);

  /// Energy grid and ε(E, n_e) cache
  std::vector<double>                energy_grid_eV_;   // size = NE
  std::vector<std::vector<double>>   epsilon_Ene_;      // [iE][ne - ne_min]
};

} // namespace ccdarksens
