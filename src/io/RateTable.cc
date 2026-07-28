// ============================================================================
//  CCDarkSens — RateTable
//  Loads two-column QEDark-style rate CSVs and interpolates dR/dE onto a uniform ROOT energy histogram in events/(kg·year·eV).
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include "ccdarksens/io/RateTable.hh"
#include <TH1D.h>

#include <cctype>
#include <fstream>
#include <sstream>
#include <string>

namespace ccdarksens {

namespace {
inline bool is_comment_or_empty(const std::string& s) {
  for (char c : s) { if (c == '#') return true; if (!std::isspace(static_cast<unsigned char>(c))) return false; }
  return true;
}

inline void trim(std::string& s) {
  size_t i = 0; while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  size_t j = s.size(); while (j > i && std::isspace(static_cast<unsigned char>(s[j-1]))) --j;
  s = s.substr(i, j - i);
}

inline bool split_line(const std::string& line, double& a, double& b) {
  // try comma first
  {
    std::stringstream ss(line);
    std::string s0, s1;
    if (std::getline(ss, s0, ',') && std::getline(ss, s1, ',')) {
      trim(s0); trim(s1);
      try { a = std::stod(s0); b = std::stod(s1); return true; } catch (...) {}
    }
  }
  // fallback: whitespace
  {
    std::stringstream ss(line);
    if (ss >> a >> b) return true;
  }
  return false;
}
} // namespace

bool RateTable::LoadCSV(const std::string& path) {
  E_eV_.clear(); R_kg_year_eV_.clear(); meta_ = RateMeta{};

  std::ifstream in(path);
  if (!in) return false;

  std::string line;
  bool saw_header = false;
  while (std::getline(in, line)) {
    if (line.empty() || is_comment_or_empty(line)) continue;
    if (!saw_header) { saw_header = true; continue; } // skip header row
    double x=0, y=0;
    if (!split_line(line, x, y)) continue;
    E_eV_.push_back(x);
    R_kg_year_eV_.push_back(y);
  }
  return (!E_eV_.empty() && E_eV_.size() == R_kg_year_eV_.size());
}

std::unique_ptr<TH1D> RateTable::MakeTH1D(const std::string& name,
                                          double Emin_eV, double Emax_eV,
                                          int nbins) const {
  auto h = std::make_unique<TH1D>(name.c_str(), name.c_str(), nbins, Emin_eV, Emax_eV);
  h->Sumw2();
  if (E_eV_.empty()) return h;

  // NOTE: bin-center sampling (will be replaced with proper trapezoid integration later).
  const int n = static_cast<int>(E_eV_.size());

  for (int i = 1; i <= nbins; ++i) {
    const double Ec = h->GetBinCenter(i);

    if (Ec <= E_eV_.front()) { h->SetBinContent(i, R_kg_year_eV_.front()); continue; }
    if (Ec >= E_eV_.back())  { h->SetBinContent(i, 0.0); continue; }

    // Linear interpolation at bin center
    size_t lo = 0, hi = static_cast<size_t>(n) - 1;
    while (hi - lo > 1) {
      size_t mid = (lo + hi) / 2;
      if (E_eV_[mid] <= Ec) lo = mid; else hi = mid;
    }
    const double x0 = E_eV_[lo], x1 = E_eV_[hi];
    const double y0 = R_kg_year_eV_[lo], y1 = R_kg_year_eV_[hi];
    h->SetBinContent(i, y0 + (Ec - x0) / (x1 - x0) * (y1 - y0));
  }

  return h;
}

} // namespace ccdarksens
