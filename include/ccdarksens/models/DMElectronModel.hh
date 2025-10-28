#pragma once
#include <vector>
#include <string>

namespace ccdarksens {

struct DMParams {
  double mass_GeV = 1.0;
  double sigma_ref = 1e-38;
  std::string mediator = "heavy";
};

struct EnergyGrid {
  std::vector<double> edges;
  std::vector<double> centers;
};

class DMElectronModel {
public:
  explicit DMElectronModel(DMParams p);
  const DMParams& params() const noexcept { return params_; }

  // compute differential rate [events / (kg·day·eV)]
  std::vector<double> differential_rate(const EnergyGrid& grid) const;

private:
  DMParams params_;
};

} // namespace ccdarksens
