#pragma once
#include <memory>
#include <string>
#include "ccdarksens/io/RateTable.hh" 

class TH1D;

namespace ccdarksens {

struct DMElectronConfig {
  std::string material;           // "Si"
  std::string mediator;           // "heavy" | "massless"
  std::string rates_dir;          // e.g. data/qedark_rates/Si/heavy
  std::string filename_template;  // e.g. dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv
  double      mchi_MeV = 0.0;     // used in filename
  std::string sigma_e_cm2;        // keep as string to match filename exactly

  // output spectrum binning
  double Emin_eV = 0.0;
  double Emax_eV = 20.0;
  int    nbins   = 200;
};

// class RateTable;

class DMElectronModel {
public:
  DMElectronModel() = default;
  ~DMElectronModel();                 // <-- declare only (no = default here)
  bool Configure(const DMElectronConfig& cfg);

  /// dR/dE in events / kg / day / eV
  std::unique_ptr<TH1D> MakeSpectrum_E() const;

  const DMElectronConfig& cfg() const noexcept { return cfg_; }

private:
  std::string ResolvePath_() const;
  DMElectronConfig              cfg_{};
  std::unique_ptr<RateTable>    table_;
};

} // namespace ccdarksens
