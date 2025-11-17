#pragma once

#include "ccdarksens/stats/ITestStatistic.hh"

namespace ccdarksens::stats {

/**
 * Poisson binned likelihood-based test statistic.
 *
 * EvaluateNLL:
 *   -ln L = sum_i [ mu_i - n_i ln(mu_i) ]   (up to additive constants)
 *
 * EvaluateRatio:
 *   q = -2 ln [ L(data | model_test) / L(data | model_null) ]
 *     = 2 [ NLL(data | model_test) - NLL(data | model_null) ]
 *
 * In the Asimov background-only setup (data=B, null=B, test=B+S),
 * q reduces to the familiar
 *   q = 2 * sum_i [ S_i - B_i * ln(1 + S_i / B_i) ].
 */
class PoissonAsimovPLR : public ITestStatistic {
public:
  PoissonAsimovPLR() = default;
  ~PoissonAsimovPLR() override = default;

  double EvaluateNLL(const std::vector<double>& data,
                     const std::vector<double>& model) const override;

  double EvaluateRatio(const std::vector<double>& data,
                       const std::vector<double>& model_test,
                       const std::vector<double>& model_null) const override;
};

} // namespace ccdarksens::stats
