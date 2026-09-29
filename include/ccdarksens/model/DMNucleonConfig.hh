// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  DMNucleonConfig.hh -- I define DMNucleonConfig: the DM-nucleus coupling
//  settings (target nucleus, mediator, rate-file location, grid point,
//  output binning) shared by the Migdal and WIMP-nucleon models.
// ===========================================================================

#pragma once
#include <string>

namespace ccdarksens {

// Describes the DM <-> nucleus coupling itself, independent of which
// observable (Migdal ionization spectrum today, nuclear-recoil spectrum for
// a future WIMP search) is derived from it. `mediator` here names the
// DM-nucleon mediator and is a distinct namespace from the DM-electron
// `mediator` field used by DMElectronConfig/DarkPhotonConfig — the same
// string values ("heavy"/"light") describe a different vertex.
struct DMNucleonConfig {
  std::string target_nucleus;     // e.g. "Si28" (nuclear target, distinct from "material")
  int         A = 28;             // mass number
  int         Z = 14;             // atomic number
  std::string mediator;           // "heavy" | "light" (DM-nucleon mediator)
  std::string rates_dir;          // e.g. data/migdal_rates/Si/heavy
  std::string filename_template;  // e.g. dRdE_{target_nucleus}_{mediator}_m{mchi_MeV}_s{sigma_n_cm2}.csv
  double      mchi_MeV = 0.0;     // DM mass, used in filename
  std::string sigma_n_cm2;        // DM-nucleon cross section, kept as string to match filename exactly

  // output spectrum binning (whatever energy observable the consuming model produces)
  double Emin_eV = 0.0;  // lower edge of the output spectrum [eV]
  double Emax_eV = 20.0;  // upper edge of the output spectrum [eV]
  int    nbins   = 200;  // number of output bins
};

} // namespace ccdarksens
