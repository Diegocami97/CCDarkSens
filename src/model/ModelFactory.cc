// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ModelFactory.cc -- I dispatch on the model type in the config and build
//  the corresponding signal spectrum. Adding a new model type only needs a
//  new branch here.
// ===========================================================================

#include "ccdarksens/model/ModelFactory.hh"

#include <TH1D.h>

#include "ccdarksens/model/DMElectronModel.hh"
#include "ccdarksens/model/DarkPhotonModel.hh"
#include "ccdarksens/model/MigdalModel.hh"
#include "ccdarksens/model/WimpNucleonModel.hh"

namespace ccdarksens {

// ----------------------------------------------------------------------------
// MakeSignalSpectrumE
//   I copy the relevant ModelJSON fields into the config struct of the
//   selected model (dark_photon, migdal, wimp_nucleon, otherwise dm_electron),
//   configure it for the given mass and coupling string, and return its
//   dR/dE spectrum [events/(kg*year*eV)]. If ok is non-null it is set to
//   whether the rate CSV was found and loaded.
// ----------------------------------------------------------------------------
std::unique_ptr<TH1D> MakeSignalSpectrumE(const ModelJSON& mj,
                                           double mass_val,
                                           const std::string& coupling_str,
                                           bool* ok)
{
  bool good = false;
  std::unique_ptr<TH1D> h;

  if (mj.type == "dark_photon") {
    DarkPhotonConfig c;
    c.material          = mj.material;
    c.mediator          = mj.mediator;
    c.rates_dir         = mj.rates_dir;
    c.filename_template = mj.filename_template;
    c.Emin_eV           = mj.Emin_eV;
    c.Emax_eV           = mj.Emax_eV;
    c.nbins             = mj.nbins;
    c.mA_eV             = mass_val;
    c.epsilon           = coupling_str;
    c.epsilon_ref       = mj.epsilon_ref;
    DarkPhotonModel m;
    good = m.Configure(c);
    h    = m.MakeSpectrum_E();

  } else if (mj.type == "migdal") {
    DMNucleonConfig c;
    c.target_nucleus    = mj.target_nucleus;
    c.A                 = mj.nuclear_A;
    c.Z                 = mj.nuclear_Z;
    c.mediator          = mj.mediator;
    c.rates_dir         = mj.rates_dir;
    c.filename_template = mj.filename_template;
    c.Emin_eV           = mj.Emin_eV;
    c.Emax_eV           = mj.Emax_eV;
    c.nbins             = mj.nbins;
    c.mchi_MeV          = mass_val;
    c.sigma_n_cm2       = coupling_str;
    MigdalModel m;
    good = m.Configure(c);
    h    = m.MakeSpectrum_E();

  } else if (mj.type == "wimp_nucleon") {
    DMNucleonConfig c;
    c.target_nucleus    = mj.target_nucleus;
    c.A                 = mj.nuclear_A;
    c.Z                 = mj.nuclear_Z;
    c.mediator          = mj.mediator;
    c.rates_dir         = mj.rates_dir;
    c.filename_template = mj.filename_template;
    c.Emin_eV           = mj.Emin_eV;
    c.Emax_eV           = mj.Emax_eV;
    c.nbins             = mj.nbins;
    c.mchi_MeV          = mass_val;
    c.sigma_n_cm2       = coupling_str;
    WimpNucleonModel m;
    good = m.Configure(c);
    h    = m.MakeSpectrum_E();

  } else { // "dm_electron" (default)
    DMElectronConfig c;
    c.material          = mj.material;
    c.mediator          = mj.mediator;
    c.rates_dir         = mj.rates_dir;
    c.filename_template = mj.filename_template;
    c.Emin_eV           = mj.Emin_eV;
    c.Emax_eV           = mj.Emax_eV;
    c.nbins             = mj.nbins;
    c.mchi_MeV          = mass_val;
    c.sigma_e_cm2       = coupling_str;
    DMElectronModel m;
    good = m.Configure(c);
    h    = m.MakeSpectrum_E();
  }

  if (ok) *ok = good;
  return h;
}

} // namespace ccdarksens
