#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/Diffusion.hh"
#include <TH1D.h>

namespace ccdarksens {

std::unique_ptr<TH1D> DetectorResponsePipeline::Apply(const TH1D& dRdE,
                                                      double exposure_kg_year,
                                                      int ne_min, int ne_max,
                                                      double Ee_ref_eV) const {
  auto h = ion_->FoldToNe(dRdE, exposure_kg_year, ne_min, ne_max);
  if (diff_) diff_->Apply(*h, Ee_ref_eV);
  if (pe_)   pe_->Apply(*h);
  return h;
}

} // namespace ccdarksens
