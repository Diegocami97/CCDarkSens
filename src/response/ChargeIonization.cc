#include "ccdarksens/response/ChargeIonization.hh"
#include <TH1D.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ccdarksens {

static bool is_comment_or_empty(const std::string& s) {
  for (char c : s) { if (!std::isspace((unsigned char)c)) return c=='#'; }
  return true;
}

ChargeIonization::ChargeIonization(std::string table_csv) { LoadCSV_(table_csv); }

void ChargeIonization::LoadCSV_(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("ChargeIonization: cannot open " + path);

  std::vector<double> E;
  std::vector<std::vector<double>> cols;

  std::string line;
  while (std::getline(in, line)) {
    if (is_comment_or_empty(line)) continue;
    std::istringstream ss(line);
    std::string cell;
    std::vector<double> vals;
    while (std::getline(ss, cell, ',')) {
      if (!cell.empty() && (cell.back()=='\r' || cell.back()=='\n')) cell.pop_back();
      vals.push_back(std::stod(cell));
    }
    if (vals.size()<2) continue;
    E.push_back(vals[0]);
    if (cols.empty()) cols.resize(vals.size()-1);
    for (size_t j=1;j<vals.size();++j) cols[j-1].push_back(vals[j]);
  }
  if (E.size()<2 || cols.empty()) throw std::runtime_error("ChargeIonization: malformed table");

  // Normalize rows for robustness
  for (size_t i=0;i<E.size();++i) {
    double s=0.0; for (auto& c: cols) s += c[i];
    if (s>0) for (auto& c: cols) c[i] /= s;
  }
  pn_given_E_.clear();
  for (auto& c : cols) pn_given_E_.push_back({E, c});
}

double ChargeIonization::InterpLinearClamped(const std::vector<double>& x,
                                             const std::vector<double>& y,
                                             double xq) {
  if (xq <= x.front()) return std::clamp(y.front(), 0.0, 1.0);
  if (xq >= x.back())  return std::clamp(y.back(),  0.0, 1.0);
  auto it = std::upper_bound(x.begin(), x.end(), xq);
  size_t i1 = size_t(it - x.begin());
  size_t i0 = i1 - 1;
  double t = (xq - x[i0]) / (x[i1] - x[i0]);
  double v = (1.0 - t)*y[i0] + t*y[i1];
  return std::clamp(v, 0.0, 1.0);
}

std::unique_ptr<TH1D> ChargeIonization::FoldToNe(const TH1D& dRdE,
                                                  double exposure_kg_year,
                                                  int ne_min, int ne_max) const {
  if (ne_max < ne_min) throw std::invalid_argument("ne_max < ne_min");
  const int k = MaxNeFromTable();
  const int nbins = dRdE.GetNbinsX();
  if (nbins <= 0) throw std::invalid_argument("dRdE has no bins");

  const int n_out = ne_max - ne_min + 1;
  std::vector<double> edges(n_out+1);
  for (int i=0;i<=n_out;++i) edges[i] = (ne_min - 0.5) + i;

  auto h = std::make_unique<TH1D>("S_ne","Signal in n_{e};n_{e};counts",
                                  n_out, edges.data());

  for (int i=1;i<=nbins;++i) {
    const double Ei   = dRdE.GetBinCenter(i);
    const double dEi  = dRdE.GetBinWidth(i);
    const double rate = dRdE.GetBinContent(i);  // events/(kg year eV)
    // const double counts = rate * dEi * exposure_kg_year;
    const double counts = rate * exposure_kg_year * dEi;
    if (counts <= 0) continue;

    for (int n=1; n<=k; ++n) {
      const auto& E  = pn_given_E_[size_t(n-1)].first;
      const auto& Pn = pn_given_E_[size_t(n-1)].second;
      const double p = InterpLinearClamped(E, Pn, Ei);
      const double contrib = counts * p;
      if (n >= ne_min && n <= ne_max)
        h->AddBinContent(h->FindBin(n), contrib);
    }
  }
  return h;
}

} // namespace ccdarksens
