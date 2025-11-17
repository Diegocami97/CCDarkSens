#include "ccdarksens/experiment/ExperimentSetup.hh"
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace ccdarksens {

static std::string mode_to_string(ExperimentMode m) {
  switch (m) {
    case ExperimentMode::Observed: return "observed";
    case ExperimentMode::Asimov:   return "asimov";
    case ExperimentMode::Toys:     return "toys";
  }
  return "unknown";
}

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

ExperimentSummary ExperimentSetup::prepare_summary() const {
  ExperimentSummary s;
  // s.exposure_kg_day = detector_mass_kg_ * cfg_.livetime_days * cfg_.duty_cycle;
  s.exposure_kg_year = cfg_.livetime_days * cfg_.duty_cycle * detector_mass_kg_ / 365; // NEED to fix for later
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
  return s;
}

} // namespace ccdarksens
