// ============================================================================
//  CCDarkSens — DMElectronModel
//  Resolves and loads QEDark/QCDark rate CSV paths for a given (material, mediator, mχ, σe) and builds a binned dR/dE TH1D spectrum.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

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

DMElectronModel::~DMElectronModel() = default;

static std::string format_mchi_6f(double x) {
  std::ostringstream os;
  os.setf(std::ios::fmtflags(0), std::ios::floatfield);
  os << std::fixed << std::setprecision(6) << x;
  return os.str();
}

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

bool DMElectronModel::Configure(const DMElectronConfig& c) {
  cfg_ = c;
  table_ = std::make_unique<RateTable>();
  const auto path = ResolvePath_();
  if (!table_->LoadCSV(path)) {
    return false;
  }
  return true;
}

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
