// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  RateSource.cc -- I create and fill synthetic dR/dE histograms (flat and
//  mono-energetic line) for tests and placeholder spectra.
// ===========================================================================

#include "ccdarksens/rates/RateSource.hh"
#include <TH1D.h>
#include <stdexcept>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// RateSource::MakeLinearEnergyHist
//   Empty histogram with nbins uniform bins in [emin, emax] eV.
//   Throws std::invalid_argument if nbins <= 0 or emax <= emin.
// ----------------------------------------------------------------------------
std::unique_ptr<TH1D> RateSource::MakeLinearEnergyHist(double emin, double emax, int nbins,
                                                       const char* name) {
  if (nbins<=0 || emax<=emin) throw std::invalid_argument("bad energy hist");
  return std::make_unique<TH1D>(name, "dR/dE;E_{e} [eV];events/(kg day eV)", nbins, emin, emax);
}

// ----------------------------------------------------------------------------
// RateSource::FillFlat
//   Fill every bin with the same value so that the histogram integrates to
//   norm_per_kg_day over its full energy range.
// ----------------------------------------------------------------------------
void RateSource::FillFlat(TH1D& dRdE, double norm_per_kg_day) {
  const int nb = dRdE.GetNbinsX();
  const double width_total = dRdE.GetXaxis()->GetXmax() - dRdE.GetXaxis()->GetXmin();
  const double c = (width_total>0) ? (norm_per_kg_day / width_total) : 0.0;
  for (int i=1;i<=nb;++i) dRdE.SetBinContent(i, c);
}

// ----------------------------------------------------------------------------
// RateSource::FillMonoLine
//   Put a mono-energetic line at E0_eV into the bin containing it, scaled so
//   that the bin integral equals norm_per_kg_day (E0 outside the range is
//   ignored).
// ----------------------------------------------------------------------------
void RateSource::FillMonoLine(TH1D& dRdE, double E0_eV, double norm_per_kg_day) {
  const int ibin = dRdE.FindBin(E0_eV);
  if (ibin<1 || ibin>dRdE.GetNbinsX()) return;
  const double dE = dRdE.GetBinWidth(ibin);
  dRdE.SetBinContent(ibin, (dE>0) ? (norm_per_kg_day / dE) : 0.0);
}

} // namespace ccdarksens
