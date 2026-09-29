// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  NeSpaceResponseFold.cc -- Single-pixel n_e-space ResponseFold:
//  ChargeIonization::FoldToNe + epsilon(n_e).
// ===========================================================================

#include "ccdarksens/response/NeSpaceResponseFold.hh"

#include <algorithm>
#include <TH1D.h>

#include "ccdarksens/response/ChargeIonization.hh"

namespace ccdarksens {

// Constructor: keep the ionization model and the configuration.
NeSpaceResponseFold::NeSpaceResponseFold(std::shared_ptr<ChargeIonization> ion,
                                          NeSpaceResponseFoldConfig cfg)
    : ion_(std::move(ion)), cfg_(std::move(cfg)) {}

// ----------------------------------------------------------------------------
// NeSpaceResponseFold::Fold
//   Signal per ROI n_e bin: fold dR/dE to n_e with the ionization table, then
//   multiply each ROI bin by its efficiency (clamped to [0,1]; a missing n_e counts as 0).
// ----------------------------------------------------------------------------
std::vector<double> NeSpaceResponseFold::Fold(const TH1D& dRdE_density,
                                               double exposure_kg_year) const {
  auto S_true = ion_->FoldToNe(dRdE_density, exposure_kg_year, cfg_.ne_min, cfg_.ne_max);

  std::vector<double> out;
  out.reserve(cfg_.roi_bins.size());
  for (int ne : cfg_.roi_bins) {
    const int bin = S_true->FindBin(static_cast<double>(ne));
    const double s_true = S_true->GetBinContent(bin);
    auto it = cfg_.eps_ne.find(ne);
    const double eps = (it != cfg_.eps_ne.end()) ? std::clamp(it->second, 0.0, 1.0) : 0.0;
    out.push_back(s_true * eps);
  }
  return out;
}

// ----------------------------------------------------------------------------
// NeSpaceResponseFold::FoldNeBackground
//   Same per-n_e efficiency applied to a dark-current histogram that is already in n_e.
// ----------------------------------------------------------------------------
std::vector<double> NeSpaceResponseFold::FoldNeBackground(const TH1D& B_ne_asimov) const {
  std::vector<double> out;
  out.reserve(cfg_.roi_bins.size());
  for (int ne : cfg_.roi_bins) {
    const int bin = B_ne_asimov.GetXaxis()->FindBin(static_cast<double>(ne));
    const double b_true = B_ne_asimov.GetBinContent(bin);
    auto it = cfg_.eps_ne.find(ne);
    const double eps = (it != cfg_.eps_ne.end()) ? std::clamp(it->second, 0.0, 1.0) : 0.0;
    out.push_back(b_true * eps);
  }
  return out;
}

// ----------------------------------------------------------------------------
// NeSpaceResponseFold::BinLabels
//   Labels "ne=<n>", one per ROI bin.
// ----------------------------------------------------------------------------
std::vector<std::string> NeSpaceResponseFold::BinLabels() const {
  std::vector<std::string> labels;
  labels.reserve(cfg_.roi_bins.size());
  for (int ne : cfg_.roi_bins) labels.push_back("ne=" + std::to_string(ne));
  return labels;
}

}  // namespace ccdarksens
