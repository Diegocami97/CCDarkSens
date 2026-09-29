// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  BackgroundBuilder.cc -- I build the Asimov background n_e spectrum from a
//  per-pixel Poisson dark current, scaled by the number of active pixels and
//  the number of exposures, with an optional pattern-efficiency applied on
//  top.
// ===========================================================================

#include "ccdarksens/backgrounds/BackgroundBuilder.hh"
#include "ccdarksens/backgrounds/PoissonDarkCurrent.hh"
#include "ccdarksens/response/PatternEfficiency.hh"

#include <TH1D.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <map>



namespace ccdarksens {

// ----------------------------------------------------------------------------
// BackgroundBuilder::BackgroundBuilder
//   I validate the geometry (rows, cols > 0; active_fraction in (0,1];
//   ne_max >= ne_min) and compute the number of active pixels
//   = round(rows * cols * active_fraction).
//   Throws std::invalid_argument on bad input.
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// BackgroundBuilder::SetTiming
//   I store the livetime [days], the duty cycle and the exposure timing.
//   Either n_exposures_override or a positive exposure_time_s is required so
//   that I can derive the number of exposures later on.
//   Throws std::invalid_argument if livetime <= 0, duty_cycle is outside
//   (0,1], or neither timing source is available.
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// BackgroundBuilder::SetDarkCurrent
//   I store the dark-current settings. lambda_e_per_pix_per_year must be >= 0.
// ----------------------------------------------------------------------------
void BackgroundBuilder::SetDarkCurrent(const DarkCurrentConfig& dccfg) {
  if (dccfg.lambda_e_per_pix_per_year < 0.0)
    throw std::invalid_argument("lambda_e_per_pix_per_year must be >= 0");
  dccfg_ = dccfg;
}

// ----------------------------------------------------------------------------
// BackgroundBuilder::BuildBkgAsimov
//   I build the Asimov background histogram B_obs(n_e):
//     1) number of exposures = override, or floor(livetime*86400*duty/t_exp);
//     2) lambda per exposure = lambda_year * t_exp / (365.25 d in seconds);
//     3) unit-normalized single-pixel Poisson spectrum on [ne_min, ne_max];
//     4) scale by (active pixels) * (exposures) * norm_scale;
//     5) apply the pattern efficiency, if one was set.
//   Returns a new histogram owned by the caller.
//   Throws std::runtime_error if the derived number of exposures is <= 0.
// ----------------------------------------------------------------------------
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

  // Convert λ_year -> λ_exp (e-/pix/exposure)
  const double year_s = 365.25 * 86400.0;
  const double lambda_per_exp =
    (tcfg_.exposure_time_s > 0.0)
      ? dccfg_.lambda_e_per_pix_per_year * (tcfg_.exposure_time_s / year_s)
      : 0.0;

  // Build per-pixel Poisson histogram (unit normalized), then scale
  PoissonDarkCurrentBackground pix_pois(lambda_per_exp, /*norm=*/1.0);
  auto h = pix_pois.MakeHist(ne_min_, ne_max_);

  // Scale to all pixels and all exposures, then apply global norm
  const double scale =
    static_cast<double>(n_active_pixels_) *
    static_cast<double>(n_exposures) *
    dccfg_.norm_scale;
  h->Scale(scale);

  // Apply pattern efficiency if provided
  if (pe_){
    std::cout << "[BackgroundBuilder::BuildBkgAsimov] Applying pattern efficiency to background\n";
    pe_->Apply(*h);
  }

  return h;
}

std::unique_ptr<TH1D>
// ----------------------------------------------------------------------------
// BackgroundBuilder::BuildBkgAsimov_EDependent
//   Energy-dependent entry point. The dark-current background has no true
//   energy, so I ignore the energy-dependent efficiency grid here (the signal
//   already carries the full E-dependent efficiency) and just return
//   BuildBkgAsimov().
// ----------------------------------------------------------------------------
BackgroundBuilder::BuildBkgAsimov_EDependent(
    const std::vector<double>& /*E_grid_eV*/,
    const std::vector<std::vector<double>>& /*eps_Ene*/) 
{
    // For now, just ignore the E-dependent pattern info
    // and use the existing Asimov background builder.
    // Signal is already treated with full E-dependent efficiency.
    return BuildBkgAsimov();
}


} // namespace ccdarksens
