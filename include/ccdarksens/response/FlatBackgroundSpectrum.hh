// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  FlatBackgroundSpectrum.hh -- Header for building a flat (energy-
//  independent) dR/dE background spectrum, channel-agnostic.
// ===========================================================================

#pragma once

#include <memory>
#include <string>

class TH1D;

namespace ccdarksens {

/**
 * Builds a flat dR/dE_true density spectrum (events/(kg*year*eV)) over
 * [Emin_eV, Emax_eV] in nbins bins -- the same flat-Compton-background
 * construction the DM-electron reference app does inline
 * (ccdarksens_scan_dmelectron_pattern.cc:729-739: borrow a ModelFactory-
 * shaped TH1D, Reset("ICES"), fill every bin with rate_per_kg_year_keV/1000),
 * generalized to take its own binning directly rather than borrowing a
 * ModelJSON's -- so it works for analysis spaces (cluster-energy) that have
 * no ModelJSON-shaped spectrum to borrow the binning from.
 *
 * rate_per_kg_year_keV: events/(kg*year*keV) -- the standard d.r.u.
 * convention BackgroundJSON::flat_bkg_norm_per_kg_year already uses.
 * Internally converted to per-eV (divide by 1000) and applied uniformly to
 * every bin, exactly matching the reference app's own conversion.
 *
 * name: ROOT histogram name -- caller-supplied, matching RateTable::MakeTH1D's
 * convention. Every TH1 auto-registers in ROOT's global directory by name;
 * a repeated fixed name across multiple calls (e.g. once per analysis space,
 * or once per rate value in a comparison) triggers "Replacing existing TH1"
 * warnings and a real double-ownership risk, not just cosmetic noise --
 * give each call a distinct name.
 */
std::unique_ptr<TH1D> MakeFlatDrdeSpectrum(double rate_per_kg_year_keV,
                                            double Emin_eV, double Emax_eV, int nbins,
                                            const std::string& name = "dRdE_flat");

}  // namespace ccdarksens
