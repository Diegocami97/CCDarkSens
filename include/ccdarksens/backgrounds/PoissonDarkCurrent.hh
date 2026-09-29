// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PoissonDarkCurrent.hh -- I declare PoissonDarkCurrentBackground: the
//  single-pixel, single-exposure Poisson dark-current spectrum P(n_e |
//  lambda) as a ROOT histogram. BackgroundBuilder scales it up to the whole
//  detector.
// ===========================================================================

#pragma once
#include <memory>
#include <string>

class TH1D;

namespace ccdarksens {

// Simple Poisson background in n_e for a *single pixel/exposure*.
// Here lambda_e is the mean number of electrons per pixel per exposure.
// ----------------------------------------------------------------------------
// PoissonDarkCurrentBackground
//   Poisson n_e distribution for ONE pixel in ONE exposure, with lambda_e the
//   mean number of electrons per pixel per exposure. MakeHist() returns
//   norm * P(n_e) on a unit-width n_e binning.
// ----------------------------------------------------------------------------
class PoissonDarkCurrentBackground {
public:
  PoissonDarkCurrentBackground(double lambda_e,
                               double norm=1.0,
                               std::string name="PoissonDarkCurrent")
  : lambda_(lambda_e), norm_(norm), name_(std::move(name)) {}

  std::unique_ptr<TH1D> MakeHist(int ne_min, int ne_max) const;

  double lambda() const { return lambda_; }
  double normalization() const { return norm_; }
  void set_lambda(double l) { lambda_ = l; }
  void set_normalization(double n) { norm_ = n; }
  const std::string& name() const { return name_; }

private:
  double lambda_ = 0.0;  // mean e-/pix/exposure
  double norm_   = 1.0;  // overall multiplicative normalization of the histogram
  std::string name_;  // ROOT histogram name
};

} // namespace ccdarksens
