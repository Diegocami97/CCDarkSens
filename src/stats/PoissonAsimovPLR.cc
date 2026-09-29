// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PoissonAsimovPLR.cc -- Computes Poisson -ln L and the q = -2ΔlnL ratio
//  used for Asimov and grid-scan hypothesis tests.
// ===========================================================================

#include "ccdarksens/stats/PoissonAsimovPLR.hh"
#include "ccdarksens/stats/StatsUtils.hh"

#include <cmath>
#include <stdexcept>

namespace ccdarksens::stats {

// ----------------------------------------------------------------------------
// PoissonAsimovPLR::EvaluateNLL
//   -ln L = sum_i [ mu_i - n_i ln(mu_i) ], dropping the constant ln(n_i!).
//   Bins with model <= 0 add nothing when data is also 0 and a 1e9 penalty
//   otherwise. Throws std::invalid_argument if the sizes differ.
// ----------------------------------------------------------------------------
double PoissonAsimovPLR::EvaluateNLL(const std::vector<double>& data,
                                     const std::vector<double>& model) const {
  if (data.size() != model.size()) {
    throw std::invalid_argument("PoissonAsimovPLR::EvaluateNLL: data/model size mismatch");
  }

  const std::size_t n = data.size();
  double nll = 0.0;

  for (std::size_t i = 0; i < n; ++i) {
    const double n_i = data[i];
    const double mu_i = model[i];

    if (mu_i <= 0.0) {
      if (n_i <= 0.0) {
        // both zero: likelihood term is constant; contribute nothing
        continue;
      }
      // data>0 but mu<=0: impossible; give a large penalty.
      nll += 1e9;
      continue;
    }

    // -ln L_i = mu_i - n_i ln(mu_i)   (up to constant ln(n_i!))
    nll += mu_i - n_i * safe_log(mu_i);
  }

  return nll;
}

// ----------------------------------------------------------------------------
// PoissonAsimovPLR::EvaluateRatio
//   q = 2*(NLL_test - NLL_null), floored at 0 for numerical safety.
// ----------------------------------------------------------------------------
double PoissonAsimovPLR::EvaluateRatio(const std::vector<double>& data,
                                       const std::vector<double>& model_test,
                                       const std::vector<double>& model_null) const {
  if (data.size() != model_test.size() || data.size() != model_null.size()) {
    throw std::invalid_argument("PoissonAsimovPLR::EvaluateRatio: size mismatch");
  }

  const double nll_test = EvaluateNLL(data, model_test);
  const double nll_null = EvaluateNLL(data, model_null);

  // q = -2 ln (L_test / L_null) = 2 (NLL_test - NLL_null)
  const double q = 2.0 * (nll_test - nll_null);
  return q < 0.0 ? 0.0 : q;  // numeric safety: should be >= 0
}

} // namespace ccdarksens::stats
