// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ModelFactory.hh -- I declare MakeSignalSpectrumE, the single entry point
//  that builds the dR/dE signal spectrum for any model type (dm_electron,
//  dark_photon, migdal, wimp_nucleon) from a ModelJSON block and one (mass,
//  coupling) grid point.
// ===========================================================================

#pragma once
/// Factory for building a dR/dE signal spectrum from a ModelJSON config.
///
/// Dispatches on mj.type ("dm_electron", "dark_photon", "migdal",
/// "wimp_nucleon") and returns a TH1D* (caller owns it). Adding a new model
/// type requires one new branch in src/model/ModelFactory.cc only — no scan
/// app changes needed.

#include "ccdarksens/io/ConfigManager.hh"
#include <memory>

class TH1D;

namespace ccdarksens {

/// Build dR/dE(E) [ev/kg/day/eV] for a single (mass, coupling) grid point.
///
/// mass_val     — physical mass in the model's native units
///                  dm_electron : mchi_MeV
///                  dark_photon : mA_eV
///                  migdal      : mchi_MeV
/// coupling_str — coupling constant as a formatted string matching filename template
///                  dm_electron : sigma_e_cm2
///                  dark_photon : epsilon
///                  migdal      : sigma_n_cm2
/// ok           — optional: set to true when the rate CSV was found and loaded
///
/// Returns a valid (possibly all-zero) TH1D even on load failure — callers rely
/// on this for dummy/flat-background spectra when a grid point is missing.
std::unique_ptr<TH1D> MakeSignalSpectrumE(const ModelJSON& mj,
                                           double mass_val,
                                           const std::string& coupling_str,
                                           bool* ok = nullptr);

} // namespace ccdarksens
