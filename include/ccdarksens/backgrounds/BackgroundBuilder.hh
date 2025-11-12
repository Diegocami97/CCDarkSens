#pragma once
#include <memory>
#include <optional>
#include <string>

class TH1D;

namespace ccdarksens {

struct DarkCurrentConfig {
  double lambda_e_per_pix_per_exposure = 0.0; // mean e- per pixel per image
  double norm_scale = 1.0;                    // global scale factor (usually 1)
};

struct BackgroundsConfig {
  DarkCurrentConfig dark;
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

  // Optionally set a pattern efficiency (will be applied to B)
  void SetPatternEfficiency(std::shared_ptr<PatternEfficiency> pe) { pe_ = std::move(pe); }

  // Build Asimov background histogram B_obs(ne)
  // Npix_active = rows*cols*active_fraction
  // Nexp = n_exposures_override ? value : floor(livetime_days*86400*duty_cycle/exposure_time_s)
  std::unique_ptr<TH1D> BuildBkgAsimov() ;

  // For logging/inspection
  long long NActivePixels() const noexcept { return n_active_pixels_; }
  long long NExposuresUsed() const noexcept { return n_exposures_used_; }
  double    LambdaPerExposure()   const { return dccfg_.lambda_e_per_pix_per_exposure; }

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
