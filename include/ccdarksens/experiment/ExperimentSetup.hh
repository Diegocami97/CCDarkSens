// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ExperimentSetup.hh -- I declare the experiment configuration (mode,
//  livetime, duty cycle, n_e binning, ROI selection, observable) and the
//  ExperimentSummary that the scan apps consume. ExperimentSetup turns one
//  into the other, including the exposure in kg*year.
// ===========================================================================

#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace ccdarksens {

// How the data set is produced: real data, the Asimov data set, or toy pseudo-experiments.
enum class ExperimentMode { Observed, Asimov, Toys };

// Inclusive range of n_e bins [ne_min, ne_max] used for the spectra.
struct BinningNE {
  int ne_min = 0;
  int ne_max = 0; // inclusive
};

// ----------------------------------------------------------------------------
// ExperimentConfig
//   Raw experiment settings as read from the JSON config.
// ----------------------------------------------------------------------------
struct ExperimentConfig {
  ExperimentMode mode = ExperimentMode::Asimov;  // data-set type
  double livetime_days = 0.0;  // total livetime [days]
  double duty_cycle    = 1.0;  // fraction of the livetime that is usable, (0,1]
  BinningNE binning;  // n_e range of the spectra
  std::vector<int> roi_bins;  // n_e bins entering the likelihood (empty = all bins in binning)
  /// Pattern IDs for likelihood when observable_bins == "pattern" (e.g. 11, 21, 111, 31, 22, 211).
  std::vector<int> pattern_roi;
  /// Observable for likelihood: "n_e" (default) or "pattern".
  std::string observable_bins = "n_e";
};

// ----------------------------------------------------------------------------
// ExperimentSummary
//   Derived, validated quantities the rest of the framework uses: exposure
//   in kg*year, the final ROI bin list, the observable and the RNG seed.
// ----------------------------------------------------------------------------
struct ExperimentSummary {
  double exposure_kg_year = 0.0;  // detector mass * livetime * duty cycle, in kg*year
  BinningNE binning;  // n_e range of the spectra
  std::vector<int> roi_bins;  // sorted, de-duplicated ROI n_e bins
  std::vector<int> pattern_roi;  // pattern IDs entering the likelihood in pattern space
  std::string observable_bins = "n_e";
  uint64_t rng_seed_used = 0;  // seed actually used for toys
  std::string mode_string;  // "observed" | "asimov" | "toys"
};

// ----------------------------------------------------------------------------
// ExperimentSetup
//   I validate an ExperimentConfig and produce its ExperimentSummary.
// ----------------------------------------------------------------------------
class ExperimentSetup {
public:
  ExperimentSetup(ExperimentConfig cfg, double detector_mass_kg, uint64_t rng_seed);
  ExperimentSummary prepare_summary() const;

private:
  ExperimentConfig cfg_;  // configuration I was built with
  double           detector_mass_kg_;  // active target mass [kg]
  uint64_t         rng_seed_;  // seed passed through to the summary
};

} // namespace ccdarksens
