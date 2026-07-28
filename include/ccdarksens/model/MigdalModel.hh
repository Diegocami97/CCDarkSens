// ============================================================================
//  CCDarkSens — MigdalModel
//  Header for the Migdal-effect ionization model: consumes a DMNucleonConfig
//  and loads the resulting dR/dE_e rate table.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once
#include <memory>
#include <string>
#include "ccdarksens/io/RateTable.hh"
#include "ccdarksens/model/DMNucleonConfig.hh"

class TH1D;

namespace ccdarksens {

class MigdalModel {
public:
  MigdalModel() = default;
  ~MigdalModel();
  bool Configure(const DMNucleonConfig& cfg);

  /// dR/dE_e in events / kg / day / eV
  std::unique_ptr<TH1D> MakeSpectrum_E() const;

  const DMNucleonConfig& cfg() const noexcept { return cfg_; }

private:
  std::string ResolvePath_() const;
  DMNucleonConfig                cfg_{};
  std::unique_ptr<RateTable>     table_;
};

} // namespace ccdarksens
