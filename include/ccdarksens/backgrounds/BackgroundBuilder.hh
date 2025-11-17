#pragma once
#include <memory>
#include <optional>
#include <string>

class TH1D;

namespace ccdarksens {

// NOTE: We now interpret lambda_e_per_pix_per_exposure as a *yearly*
//       rate (electrons / pixel / year). It is converted internally
//       to a per-exposure mean using exposure_time_s.
struct DarkCurrentConfig {
  // Dark current: mean electrons / pixel / YEAR.
  // This will be converted internally to electrons / pixel / exposure
  // using the exposure_time_s from TimingConfig.
  // Store λ in e-/pixel/year
  double lambda_e_per_pix_per_year = 0.0;
  double norm_scale = 1.0;                    // global scale factor (usually 1)
};

struct FlatBackgroundConfig {
  bool   enabled               = false;  // turn flat component on/off
  double rate_per_kg_year_keV  = 0.0;    // events / (kg · year · keV)
  double Emin_eV               = 0.0;    // energy range for the flat spectrum
  double Emax_eV               = 0.0;
};

struct BackgroundsConfig {
  DarkCurrentConfig dark;
  FlatBackgroundConfig flat;
  // future: add other background components here
};

struct TimingConfig {
  std::optional<int>    n_exposures_override; // if set, use this explicitly
  double                exposure_time_s = 0.0; // required to derive n_exposures if not overridden
};

class PatternEfficiency; // forward

class BackgroundBuilder {
public:
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

private:
  int rows_, cols_, ne_min_, ne_max_;
  double active_fraction_;
  long long n_active_pixels_ = 0;
  long long n_exposures_used_ = 0;

  double livetime_days_ = 0.0;
  double duty_cycle_    = 1.0;
  TimingConfig tcfg_;
  DarkCurrentConfig dccfg_;
  std::shared_ptr<PatternEfficiency> pe_;
};

} // namespace ccdarksens
