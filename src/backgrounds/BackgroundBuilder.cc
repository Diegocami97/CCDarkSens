#include "ccdarksens/backgrounds/BackgroundBuilder.hh"
#include "ccdarksens/backgrounds/PoissonDarkCurrent.hh"
#include "ccdarksens/response/PatternEfficiency.hh"

#include <TH1D.h>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace ccdarksens {

BackgroundBuilder::BackgroundBuilder(int rows, int cols, double active_fraction,
                                     int ne_min, int ne_max)
: rows_(rows), cols_(cols), ne_min_(ne_min), ne_max_(ne_max), active_fraction_(active_fraction)
{
  if (rows<=0 || cols<=0) throw std::invalid_argument("rows/cols must be >0");
  if (!(active_fraction>0.0 && active_fraction<=1.0))
    throw std::invalid_argument("active_fraction must be in (0,1]");
  if (ne_max_ < ne_min_) throw std::invalid_argument("ne_max<ne_min");

  const double npix = static_cast<double>(rows_) * static_cast<double>(cols_) * active_fraction_;
  n_active_pixels_ = static_cast<long long>(std::llround(npix));
}

void BackgroundBuilder::SetTiming(double livetime_days, double duty_cycle, const TimingConfig& tcfg) {
  if (livetime_days <= 0.0) throw std::invalid_argument("livetime_days must be > 0");
  if (!(duty_cycle > 0.0 && duty_cycle <= 1.0))
    throw std::invalid_argument("duty_cycle must be in (0,1]");

  if (!tcfg.n_exposures_override.has_value() && tcfg.exposure_time_s <= 0.0)
    throw std::invalid_argument("exposure_time_s must be >0 if n_exposures is not overridden");

  livetime_days_ = livetime_days;
  duty_cycle_    = duty_cycle;
  tcfg_          = tcfg;
}

void BackgroundBuilder::SetDarkCurrent(const DarkCurrentConfig& dccfg) {
  if (dccfg.lambda_e_per_pix_per_exposure < 0.0)
    throw std::invalid_argument("lambda must be >= 0");
  dccfg_ = dccfg;
}

std::unique_ptr<TH1D> BackgroundBuilder::BuildBkgAsimov() {
  // decide number of exposures
  long long n_exposures = 0;
  if (tcfg_.n_exposures_override.has_value()) {
    n_exposures = static_cast<long long>(*tcfg_.n_exposures_override);
  } else {
    const double live_s = livetime_days_ * 86400.0 * duty_cycle_;
    n_exposures = static_cast<long long>(std::floor(live_s / tcfg_.exposure_time_s));
  }
  if (n_exposures <= 0) throw std::runtime_error("Derived n_exposures <= 0");
  n_exposures_used_ = n_exposures;

  // Build per-pixel Poisson histogram (unit normalized), then scale
  PoissonDarkCurrentBackground pix_pois(dccfg_.lambda_e_per_pix_per_exposure, /*norm=*/1.0);
  auto h = pix_pois.MakeHist(ne_min_, ne_max_);

  // Scale to all pixels and all exposures, then apply global norm
  const double scale = static_cast<double>(n_active_pixels_) * static_cast<double>(n_exposures) * dccfg_.norm_scale;
  h->Scale(scale);

  // Apply pattern efficiency if provided
  if (pe_) pe_->Apply(*h);

  return h;
}

} // namespace ccdarksens
