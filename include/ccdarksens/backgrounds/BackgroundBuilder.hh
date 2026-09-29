// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  BackgroundBuilder.hh -- I declare the timing / dark-current / flat-
//  background configuration structs and the BackgroundBuilder class.
//  BackgroundBuilder turns them into the Asimov background histogram B(n_e):
//  a per-pixel Poisson dark-current spectrum scaled by the number of active
//  pixels and the number of exposures, optionally multiplied by a pattern
//  efficiency.
// ===========================================================================

#pragma once
#include <memory>
#include <optional>
#include <string>

class TH1D;

namespace ccdarksens {

// NOTE: We now interpret lambda_e_per_pix_per_exposure as a *yearly*
//       rate (electrons / pixel / year). It is converted internally
//       to a per-exposure mean using exposure_time_s.
// ----------------------------------------------------------------------------
// DarkCurrentConfig
//   Dark-current settings for the n_e-space background. I keep lambda as a
//   yearly rate (e-/pixel/year); BackgroundBuilder converts it to a mean per
//   exposure using TimingConfig::exposure_time_s.
// ----------------------------------------------------------------------------
struct DarkCurrentConfig {
  // Dark current: mean electrons / pixel / YEAR.
  // This will be converted internally to electrons / pixel / exposure
  // using the exposure_time_s from TimingConfig.
  // Store λ in e-/pixel/year
  double lambda_e_per_pix_per_year = 0.0;  // mean dark-current electrons per pixel per YEAR (0 = no dark current)
  double norm_scale = 1.0;                    // global scale factor (usually 1)
};

// ----------------------------------------------------------------------------
// FlatBackgroundConfig
//   A flat (energy-independent) radiogenic background in d.r.u.
//   = events/(kg*year*keV), defined between Emin_eV and Emax_eV.
// ----------------------------------------------------------------------------
struct FlatBackgroundConfig {
  bool   enabled               = false;  // turn flat component on/off
  double rate_per_kg_year_keV  = 0.0;    // events / (kg · year · keV)
  double Emin_eV               = 0.0;    // energy range for the flat spectrum
  double Emax_eV               = 0.0;
};

// ----------------------------------------------------------------------------
// BackgroundsConfig
//   Container for every background component I currently model. New
//   components (radiogenic lines, etc.) get added here.
// ----------------------------------------------------------------------------
struct BackgroundsConfig {
  DarkCurrentConfig dark;  // per-pixel Poisson dark current
  FlatBackgroundConfig flat;  // flat Compton/radiogenic spectrum
  // future: add other background components here
};

// ----------------------------------------------------------------------------
// TimingConfig
//   How the exposure is chopped into readouts. The number of exposures is
//   either given explicitly (n_exposures_override) or derived from the
//   livetime, duty cycle and exposure_time_s.
// ----------------------------------------------------------------------------
struct TimingConfig {
  std::optional<int>    n_exposures_override; // if set, use this explicitly
  double                exposure_time_s = 0.0; // required to derive n_exposures if not overridden
};

class PatternEfficiency; // forward

// ----------------------------------------------------------------------------
// BackgroundBuilder
//   I build the Asimov background B_obs(n_e) for a CCD:
//       B(n_e) = Npix_active * Nexposures * norm_scale * Poisson(n_e; lambda_exp)
//   and optionally apply a PatternEfficiency to it. Usage: construct with the
//   geometry, call SetTiming() and SetDarkCurrent(), then BuildBkgAsimov().
// ----------------------------------------------------------------------------
class BackgroundBuilder {
public:
  // Constructor. rows x cols is the pixel array, active_fraction the live fraction
  // of it, and [ne_min, ne_max] the n_e range of the histogram I will build.
  BackgroundBuilder(int rows, int cols, double active_fraction,
                    int ne_min, int ne_max);

  // Configure timing and dark current
  void SetTiming(double livetime_days, double duty_cycle, const TimingConfig& tcfg);
  void SetDarkCurrent(const DarkCurrentConfig& dccfg);

  // Optionally set a pattern efficiency (will be applied to Bkg Asimov)
  void SetPatternEfficiency(std::shared_ptr<PatternEfficiency> pe) { pe_ = std::move(pe); }

  // Build Asimov background histogram B_obs(ne)
  // Npix_active = rows*cols*active_fraction
  // Nexp = n_exposures_override ? value : floor(livetime_days*86400*duty_cycle/exposure_time_s)
  std::unique_ptr<TH1D> BuildBkgAsimov();

  // For logging/inspection
  long long NActivePixels() const noexcept { return n_active_pixels_; }
  long long NExposuresUsed() const noexcept { return n_exposures_used_; }

  // // Effective λ per pixel per exposure used internally
  // double LambdaPerExposure() const {
  //   if (tcfg_.exposure_time_s <= 0.0) return 0.0;
  //   // convert from e-/pix/year -> e-/pix/exposure
  //   constexpr double seconds_per_year = 365.25 * 86400.0;
  //   return dccfg_.lambda_e_per_pix_per_exposure *
  //          (tcfg_.exposure_time_s / seconds_per_year);
  // }
  // Derived: mean e- per pixel per exposure (given exposure_time_s)
  double LambdaPerExposure() const {
    const double year_s = 365.25 * 86400.0;
    if (tcfg_.exposure_time_s <= 0.0) return 0.0;
    return dccfg_.lambda_e_per_pix_per_year *
           (tcfg_.exposure_time_s / year_s);
  }

  // NEW: E-dependent wrapper, currently just forwards to BuildBkgAsimov()
  std::unique_ptr<TH1D> BuildBkgAsimov_EDependent(
    const std::vector<double>& E_grid_eV,
    const std::vector<std::vector<double>>& eps_Ene) ;

private:
  int rows_, cols_, ne_min_, ne_max_;  // geometry and n_e histogram range
  double active_fraction_;  // fraction of pixels that are live
  long long n_active_pixels_ = 0;  // rows*cols*active_fraction, rounded
  long long n_exposures_used_ = 0;  // exposure count of the last BuildBkgAsimov() call

  double livetime_days_ = 0.0;  // total livetime [days]
  double duty_cycle_    = 1.0;  // fraction of the livetime spent taking data
  TimingConfig tcfg_;  // exposure timing settings
  DarkCurrentConfig dccfg_;  // dark-current settings
  std::shared_ptr<PatternEfficiency> pe_;  // optional efficiency applied to B (null = none)
};

} // namespace ccdarksens
