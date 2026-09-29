// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  DMElectronModel.cc -- I resolve and load the DM-electron scattering rate-
//  table CSV for one (mass, coupling) grid point and turn it into a dR/dE
//  histogram.
// ===========================================================================

#include "ccdarksens/model/DMElectronModel.hh"
// #include "ccdarksens/io/RateTable.hh"

#include <TH1D.h>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <string>
#include <cctype>

namespace fs = std::filesystem;

namespace ccdarksens {

// Out-of-line destructor so RateTable only has to be a complete type in this file.
DMElectronModel::~DMElectronModel() = default;

// Format a mass with exactly six decimals so it matches the rate-file names written by the grid generators.
static std::string format_mchi_6f(double x) {
  std::ostringstream os;
  os.setf(std::ios::fmtflags(0), std::ios::floatfield);
  os << std::fixed << std::setprecision(6) << x;
  return os.str();
}

// ----------------------------------------------------------------------------
// DMElectronModel::ResolvePath_
//   I substitute {material}, {mediator}, {mchi_MeV} and {sigma_e_cm2} into the filename template and prepend rates_dir.
//   Returns the full path of the CSV for this grid point.
// ----------------------------------------------------------------------------
std::string DMElectronModel::ResolvePath_() const {
  // Simple token replacement for {material},{mediator},{mchi_MeV},{sigma_e_cm2}
  std::string fname = cfg_.filename_template;

  auto repl = [&](const std::string& token, const std::string& with){
    size_t pos = 0;
    while ((pos = fname.find(token, pos)) != std::string::npos) {
      fname.replace(pos, token.size(), with);
      pos += with.size();
    }
  };

  repl("{material}",     cfg_.material);
  repl("{mediator}",     cfg_.mediator);
  repl("{mchi_MeV}",     format_mchi_6f(cfg_.mchi_MeV));
  repl("{sigma_e_cm2}",  cfg_.sigma_e_cm2);

  fs::path p = fs::path(cfg_.rates_dir) / fname;
  return p.string();
}

// ----------------------------------------------------------------------------
// DMElectronModel::Configure
//   I store the config, resolve the CSV path and load the table.
//   Returns false if the file could not be read.
// ----------------------------------------------------------------------------
bool DMElectronModel::Configure(const DMElectronConfig& c) {
  cfg_ = c;
  table_ = std::make_unique<RateTable>();
  const auto path = ResolvePath_();
  if (!table_->LoadCSV(path)) {
    return false;
  }
  return true;
}

// ----------------------------------------------------------------------------
// DMElectronModel::MakeSpectrum_E
//   I return dR/dE_e [events/(kg*year*eV)] as a histogram of nbins between
//   Emin_eV and Emax_eV, named dRdE__mchi_<m>__sigma_<s> (sanitized so it is a legal ROOT
//   name). Returns nullptr if Configure() has not succeeded.
// ----------------------------------------------------------------------------
std::unique_ptr<TH1D> DMElectronModel::MakeSpectrum_E() const {
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
  const std::string name = sanitize("dRdE__mchi_" + format_mchi_6f(cfg_.mchi_MeV) + "__sigma_" + cfg_.sigma_e_cm2);
  return table_->MakeTH1D(name.c_str(), cfg_.Emin_eV, cfg_.Emax_eV, cfg_.nbins);
}

} // namespace ccdarksens
