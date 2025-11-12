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

  const double g_to_kg = 1.0 / 1000.0;

  for (int i = 1; i <= nbins; ++i) {
    const double Ec = h->GetBinCenter(i);

    // Binary search for bracketing indices
    if (Ec <= E_eV_.front()) {
      // h->SetBinContent(i, R_g_day_eV_.front() * g_to_kg); // old method, no need for scaling now
      h->SetBinContent(i, R_kg_year_eV_.front() );
      continue;
    }
    if (Ec >= E_eV_.back()) {
      // h->SetBinContent(i, R_g_day_eV_.back() * g_to_kg);
      h->SetBinContent(i, R_kg_year_eV_.back()); // old method, no need for scaling now
      continue;
    }

    size_t lo = 0, hi = E_eV_.size() - 1;
    while (hi - lo > 1) {
      const size_t mid = (lo + hi) / 2;
      if (E_eV_[mid] <= Ec) lo = mid; else hi = mid;
    }
    const double x0 = E_eV_[lo], x1 = E_eV_[hi];
    const double y0 = R_kg_year_eV_[lo], y1 = R_kg_year_eV_[hi];
    const double t = (Ec - x0) / (x1 - x0);
    const double y = y0 + t * (y1 - y0);

    // h->SetBinContent(i, y * g_to_kg);
    h->SetBinContent(i, y);
  }

  return h;
}

} // namespace ccdarksens
