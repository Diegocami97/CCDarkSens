#pragma once

#include <vector>

namespace ccdarksens::stats {

/**
 * Abstract interface for binned test statistics.
 *
 * Works on generic vectors of bin counts:
 *  - data[i]       : observed counts in bin i
 *  - model[i]      : expected counts in bin i
 *  - model_test[i] : expected counts under "test" hypothesis
 *  - model_null[i] : expected counts under "null" hypothesis
 */
class ITestStatistic {
public:
  virtual ~ITestStatistic() = default;

  /// Negative log-likelihood: -ln L(data | model)
  /// Implementations may ignore additive constants.
  virtual double EvaluateNLL(const std::vector<double>& data,
                             const std::vector<double>& model) const = 0;

  /// Likelihood ratio:
  ///   q = -2 ln [ L(data | model_test) / L(data | model_null) ]
  virtual double EvaluateRatio(const std::vector<double>& data,
                               const std::vector<double>& model_test,
                               const std::vector<double>& model_null) const = 0;
};

} // namespace ccdarksens::stats
