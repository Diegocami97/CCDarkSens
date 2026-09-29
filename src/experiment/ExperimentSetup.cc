// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ExperimentSetup.cc -- I compute the experiment summary (exposure in
//  kg*year, ROI bin list, mode string) from the JSON experiment settings and
//  the detector mass.
// ===========================================================================

#include "ccdarksens/experiment/ExperimentSetup.hh"
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace ccdarksens {

// Human-readable name of an ExperimentMode (used in logs and output files).
static std::string mode_to_string(ExperimentMode m) {
  switch (m) {
    case ExperimentMode::Observed: return "observed";
    case ExperimentMode::Asimov:   return "asimov";
    case ExperimentMode::Toys:     return "toys";
  }
  return "unknown";
}

// ----------------------------------------------------------------------------
// ExperimentSetup::ExperimentSetup
//   I store the configuration and reject unphysical settings:
//   livetime <= 0, duty_cycle outside (0,1], or ne_max <= ne_min
//   (std::invalid_argument).
// ----------------------------------------------------------------------------
ExperimentSetup::ExperimentSetup(ExperimentConfig cfg,
                                 double detector_mass_kg,
                                 uint64_t rng_seed)
  : cfg_(std::move(cfg)), detector_mass_kg_(detector_mass_kg), rng_seed_(rng_seed)
{
  if (cfg_.livetime_days <= 0.0)
    throw std::invalid_argument("livetime_days must be > 0");
  if (cfg_.duty_cycle <= 0.0 || cfg_.duty_cycle > 1.0)
    throw std::invalid_argument("duty_cycle must be in (0,1]");
  if (cfg_.binning.ne_max <= cfg_.binning.ne_min)
    throw std::invalid_argument("ne_max must be > ne_min");
}

// ----------------------------------------------------------------------------
// ExperimentSetup::prepare_summary
//   I compute
//       exposure [kg*year] = livetime_days * duty_cycle * mass_kg / 365.25
//   and build the ROI list: every n_e in [ne_min, ne_max] if the config gave
//   none, otherwise the config list sorted and de-duplicated. The pattern ROI
//   and the observable name are passed through unchanged.
// ----------------------------------------------------------------------------
ExperimentSummary ExperimentSetup::prepare_summary() const {
  ExperimentSummary s;
  // Exposure in kg·year: (livetime_days * duty_cycle) * mass_kg / days_per_year
  constexpr double days_per_year = 365.25;
  s.exposure_kg_year = cfg_.livetime_days * cfg_.duty_cycle * detector_mass_kg_ / days_per_year;
  s.binning = cfg_.binning;
  s.rng_seed_used = rng_seed_;
  s.mode_string = mode_to_string(cfg_.mode);

  if (cfg_.roi_bins.empty()) {
    s.roi_bins.reserve(static_cast<size_t>(cfg_.binning.ne_max - cfg_.binning.ne_min + 1));
    for (int ne = cfg_.binning.ne_min; ne <= cfg_.binning.ne_max; ++ne)
      s.roi_bins.push_back(ne);
  } else {
    s.roi_bins = cfg_.roi_bins;
    std::sort(s.roi_bins.begin(), s.roi_bins.end());
    s.roi_bins.erase(std::unique(s.roi_bins.begin(), s.roi_bins.end()), s.roi_bins.end());
  }
  s.pattern_roi      = cfg_.pattern_roi;
  s.observable_bins  = cfg_.observable_bins;
  return s;
}

} // namespace ccdarksens
