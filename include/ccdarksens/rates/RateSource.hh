// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  RateSource.hh -- I declare RateSource: helpers that create and fill
//  synthetic dR/dE histograms (flat and mono-energetic lines) for tests and
//  for placeholder spectra.
// ===========================================================================

#pragma once
#include <memory>

class TH1D;

namespace ccdarksens {

// ----------------------------------------------------------------------------
// RateSource
//   Static helpers for synthetic spectra: a linear-energy histogram, a flat
//   spectrum, and a single mono-energetic line.
// ----------------------------------------------------------------------------
class RateSource {
public:
  static std::unique_ptr<TH1D> MakeLinearEnergyHist(double emin_eV, double emax_eV, int nbins,
                                                    const char* name="dRdE");
  static void FillFlat(TH1D& dRdE, double norm_per_kg_day);
  static void FillMonoLine(TH1D& dRdE, double E0_eV, double norm_per_kg_day);
};

} // namespace ccdarksens
