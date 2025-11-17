#include "ccdarksens/backgrounds/PoissonDarkCurrent.hh"
#include <TH1D.h>
#include <cmath>
#include <vector>

namespace ccdarksens {

static double pois_pmf(int n, double lambda) {
  if (lambda <= 0.0) return (n==0) ? 1.0 : 0.0;
  double logp = -lambda + n*std::log(lambda);
  double lf = 0.0;
  for (int k=2;k<=n;++k) lf += std::log(double(k));
  return std::exp(logp - lf);
}

std::unique_ptr<TH1D> PoissonDarkCurrentBackground::MakeHist(int ne_min, int ne_max) const {
  const int nbin = ne_max - ne_min + 1;
  std::vector<double> edges(nbin+1);
  for (int i=0;i<=nbin;++i) edges[i] = (ne_min - 0.5) + i;

  auto h = std::make_unique<TH1D>(name_.c_str(), (name_+";n_{e};counts").c_str(),
                                  nbin, edges.data());
  for (int n=ne_min; n<=ne_max; ++n) {
    const double p = (n>=0) ? pois_pmf(n, lambda_) : 0.0;
    h->SetBinContent(h->FindBin(n), norm_ * p);
  }
  return h;
}

} // namespace ccdarksens
