#pragma once
#include <memory>

class TH1D;

namespace ccdarksens {
class ChargeIonization;
class PatternEfficiency;
class Diffusion;

class DetectorResponsePipeline {
public:
  explicit DetectorResponsePipeline(std::shared_ptr<ChargeIonization> ion)
  : ion_(std::move(ion)) {}

  void SetDiffusion(std::shared_ptr<Diffusion> diff) { diff_ = std::move(diff); }
  void SetPatternEfficiency(std::shared_ptr<PatternEfficiency> pe) { pe_ = std::move(pe); }

  // Ee_ref_eV: characteristic recoil energy for diffusion sigma (use 0 if beta=0)
  std::unique_ptr<TH1D> Apply(const TH1D& dRdE,
                              double exposure_kg_day,
                              int ne_min, int ne_max,
                              double Ee_ref_eV = 0.0) const;

private:
  std::shared_ptr<ChargeIonization>  ion_;
  std::shared_ptr<Diffusion>         diff_;
  std::shared_ptr<PatternEfficiency> pe_;
};

} // namespace ccdarksens
