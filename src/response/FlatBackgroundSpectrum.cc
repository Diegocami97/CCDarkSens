// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  FlatBackgroundSpectrum.cc -- Flat (energy-independent) dR/dE background
//  spectrum, channel-agnostic.
// ===========================================================================

#include "ccdarksens/response/FlatBackgroundSpectrum.hh"

#include <TH1D.h>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// MakeFlatDrdeSpectrum
//   Histogram of nbins between Emin_eV and Emax_eV with every bin set to the
//   flat rate converted from events/(kg*year*keV) to events/(kg*year*eV).
//   The caller supplies the ROOT name (it must be unique).
// ----------------------------------------------------------------------------
std::unique_ptr<TH1D> MakeFlatDrdeSpectrum(double rate_per_kg_year_keV,
                                            double Emin_eV, double Emax_eV, int nbins,
                                            const std::string& name) {
  auto h = std::make_unique<TH1D>(name.c_str(), "flat background;E [eV];events/(kg*year*eV)",
                                   nbins, Emin_eV, Emax_eV);
  const double flat_rate_per_eV = rate_per_kg_year_keV / 1000.0;  // keV -> eV, same conversion as the reference app
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
    h->SetBinContent(ib, flat_rate_per_eV);
  }
  return h;
}

}  // namespace ccdarksens
