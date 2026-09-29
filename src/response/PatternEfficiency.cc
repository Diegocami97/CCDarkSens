// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternEfficiency.cc -- Stores and multiplies a per-bin ε(n_e) histogram
//  into target n_e spectra.
// ===========================================================================

#include "ccdarksens/response/PatternEfficiency.hh"
#include <TH1D.h>
#include <algorithm>
#include <stdexcept>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// PatternEfficiency::SetEfficiencyHist
//   Store a clone of the efficiency histogram.
// ----------------------------------------------------------------------------
void PatternEfficiency::SetEfficiencyHist(const TH1D& epsilon_ne) {
  eps_ = std::unique_ptr<TH1D>(static_cast<TH1D*>(epsilon_ne.Clone("eps_ne")));
}

// ----------------------------------------------------------------------------
// PatternEfficiency::Apply
//   Multiply every bin of target_ne by epsilon(n_e), clamped to [0,1].
//   Does nothing if no efficiency was set. Throws std::invalid_argument if the
//   two histograms have different bin counts.
// ----------------------------------------------------------------------------
void PatternEfficiency::Apply(TH1D& target_ne) const {
  if (!eps_) return;
  if (target_ne.GetNbinsX() != eps_->GetNbinsX())
    throw std::invalid_argument("PatternEfficiency: bin count mismatch");
  for (int i=1;i<=target_ne.GetNbinsX();++i) {
    const double e = std::clamp(eps_->GetBinContent(i), 0.0, 1.0);
    target_ne.SetBinContent(i, target_ne.GetBinContent(i) * e);
  }
}

} // namespace ccdarksens
