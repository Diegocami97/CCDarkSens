// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  DarkPhotonModel.cc -- I resolve and load the hidden-photon absorption
//  rate-table CSV for one (mass, coupling) grid point and turn it into a
//  dR/dE histogram.
// ===========================================================================

#include "ccdarksens/model/DarkPhotonModel.hh"

#include <TH1D.h>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <string>
#include <cctype>

namespace fs = std::filesystem;

namespace ccdarksens {

// Out-of-line destructor so RateTable only has to be a complete type in this file.
DarkPhotonModel::~DarkPhotonModel() = default;

// Format a mass with exactly six decimals so it matches the rate-file names written by the grid generators.
static std::string format_mA_6f(double x) {
  std::ostringstream os;
  os.setf(std::ios::fmtflags(0), std::ios::floatfield);
  os << std::fixed << std::setprecision(6) << x;
  return os.str();
}

// ----------------------------------------------------------------------------
// DarkPhotonModel::ResolvePath_
//   I substitute {material}, {mediator}, {mA_eV} and {epsilon} into the filename template and prepend rates_dir.
//   Returns the full path of the CSV for this grid point.
// ----------------------------------------------------------------------------
std::string DarkPhotonModel::ResolvePath_() const {
  // Simple token replacement for {material},{mediator},{mA_eV},{epsilon}
  std::string fname = cfg_.filename_template;

  auto repl = [&](const std::string& token, const std::string& with){
    size_t pos = 0;
    while ((pos = fname.find(token, pos)) != std::string::npos) {
      fname.replace(pos, token.size(), with);
      pos += with.size();
    }
  };

  repl("{material}",  cfg_.material);
  repl("{mediator}",  cfg_.mediator);
  repl("{mA_eV}",     format_mA_6f(cfg_.mA_eV));
  repl("{epsilon}",   cfg_.epsilon_ref.empty() ? cfg_.epsilon : cfg_.epsilon_ref);

  fs::path p = fs::path(cfg_.rates_dir) / fname;
  return p.string();
}

// ----------------------------------------------------------------------------
// DarkPhotonModel::Configure
//   I store the config, resolve the CSV path and load the table.
//   Returns false if the file could not be read.
// ----------------------------------------------------------------------------
bool DarkPhotonModel::Configure(const DarkPhotonConfig& c) {
  cfg_ = c;
  table_ = std::make_unique<RateTable>();
  const auto path = ResolvePath_();
  if (!table_->LoadCSV(path)) {
    return false;
  }
  return true;
}

// ----------------------------------------------------------------------------
// DarkPhotonModel::MakeSpectrum_E
//   I return dR/dE_e [events/(kg*year*eV)] as a histogram of nbins between
//   Emin_eV and Emax_eV, named dRdE__mA_<m>__eps_<e> (sanitized so it is a legal ROOT
//   name). Returns nullptr if Configure() has not succeeded.
//   If epsilon_ref is set, I read the file at epsilon_ref and multiply the rate
//   by (epsilon/epsilon_ref)^2 -- the absorption rate scales as epsilon^2.
// ----------------------------------------------------------------------------
std::unique_ptr<TH1D> DarkPhotonModel::MakeSpectrum_E() const {
  if (!table_) return nullptr;
  auto sanitize = [](std::string s) {
    for (char& c : s) {
      if (c == '.') c = 'p';
      else if (c == '-') c = 'm';
      else if (c == '+') c = 'p';
      else if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_')) c = '_';
    }
    return s;
  };
  const std::string name = sanitize("dRdE__mA_" + format_mA_6f(cfg_.mA_eV) + "__eps_" + cfg_.epsilon);
  auto h = table_->MakeTH1D(name.c_str(), cfg_.Emin_eV, cfg_.Emax_eV, cfg_.nbins);
  if (h && !cfg_.epsilon_ref.empty() && !cfg_.epsilon.empty()) {
    const double eps     = std::stod(cfg_.epsilon);
    const double eps_ref = std::stod(cfg_.epsilon_ref);
    const double scale   = std::pow(eps / eps_ref, 2.0);
    for (int b = 1; b <= h->GetNbinsX(); ++b)
      h->SetBinContent(b, h->GetBinContent(b) * scale);
  }
  return h;
}

} // namespace ccdarksens
