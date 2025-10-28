#pragma once
#include <memory>

class TH1D;

namespace ccdarksens {

class RateSource {
public:
  static std::unique_ptr<TH1D> MakeLinearEnergyHist(double emin_eV, double emax_eV, int nbins,
                                                    const char* name="dRdE");
  static void FillFlat(TH1D& dRdE, double norm_per_kg_day);
  static void FillMonoLine(TH1D& dRdE, double E0_eV, double norm_per_kg_day);
};

} // namespace ccdarksens
