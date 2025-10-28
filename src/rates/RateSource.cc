#include "ccdarksens/rates/RateSource.hh"
#include <TH1D.h>
#include <stdexcept>

namespace ccdarksens {

std::unique_ptr<TH1D> RateSource::MakeLinearEnergyHist(double emin, double emax, int nbins,
                                                       const char* name) {
  if (nbins<=0 || emax<=emin) throw std::invalid_argument("bad energy hist");
  return std::make_unique<TH1D>(name, "dR/dE;E_{e} [eV];events/(kg day eV)", nbins, emin, emax);
}

void RateSource::FillFlat(TH1D& dRdE, double norm_per_kg_day) {
  const int nb = dRdE.GetNbinsX();
  const double width_total = dRdE.GetXaxis()->GetXmax() - dRdE.GetXaxis()->GetXmin();
  const double c = (width_total>0) ? (norm_per_kg_day / width_total) : 0.0;
  for (int i=1;i<=nb;++i) dRdE.SetBinContent(i, c);
}

void RateSource::FillMonoLine(TH1D& dRdE, double E0_eV, double norm_per_kg_day) {
  const int ibin = dRdE.FindBin(E0_eV);
  if (ibin<1 || ibin>dRdE.GetNbinsX()) return;
  const double dE = dRdE.GetBinWidth(ibin);
  dRdE.SetBinContent(ibin, (dE>0) ? (norm_per_kg_day / dE) : 0.0);
}

} // namespace ccdarksens
