#pragma once
#include <memory>

class TH1D;

namespace ccdarksens {

class PatternEfficiency {
public:
  void SetEfficiencyHist(const TH1D& epsilon_ne); // clones input
  void Apply(TH1D& target_ne) const;              // multiplies bin-by-bin

private:
  std::unique_ptr<TH1D> eps_;
};

} // namespace ccdarksens
