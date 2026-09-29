// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  BackgroundFactory.cc -- Builds background vectors per
//  run.background_source, reusing BackgroundBuilder (dark current) and
//  MakeFlatDrdeSpectrum (flat d.r.u.) exactly as the reference app does,
//  generalized over ResponseFold.
// ===========================================================================

#include "ccdarksens/response/BackgroundFactory.hh"

#include <TH1D.h>

#include <cmath>
#include <stdexcept>

#include "ccdarksens/backgrounds/BackgroundBuilder.hh"
#include "ccdarksens/response/FlatBackgroundSpectrum.hh"
#include "ccdarksens/response/ResponseFold.hh"

namespace ccdarksens {

namespace {

// ----------------------------------------------------------------------------
// MakeDcFlatMigrationBackground
//   Background for background_source = "dc_flat_migration". I add up to three
//   pieces, each already folded into the response bins:
//     1) Poisson dark current from BackgroundBuilder, folded n_e -> output bins
//        (only when the fold supports it, i.e. pattern / n_e space);
//     2) the flat d.r.u. spectrum, folded through the signal response, or -- if
//        response.background_efficiency_csv is set -- weighted by the background's
//        own efficiency curve;
//     3) the low-energy excess (LEE), if enabled (WIMP-nucleon channel only).
//   Bp is returned as zeros and Br as the total, so Bp_theta_Br also works.
// ----------------------------------------------------------------------------
BackgroundResult MakeDcFlatMigrationBackground(const ConfigManager& cfg,
                                                const ExperimentSummary& summary,
                                                const ResponseFold& fold,
                                                const ResponseFactoryResult& response) {
  // Dark current has no meaning for analysis spaces with no pixel/n_e model
  // (cluster_energy) -- degrade gracefully to flat-only there rather than
  // throwing, since "dc_flat_migration" is also the default background_source
  // and shouldn't require every channel to fake a DC contribution.
  std::vector<double> B_dc;  // dark-current contribution per output bin (stays empty if unsupported)
  if (fold.SupportsDarkCurrentBackground()) {
    const auto& det = cfg.detector();

    BackgroundBuilder bld(det.geometry().rows, det.geometry().cols, det.geometry().active_fraction,
                           response.ne_min_bkg, response.ne_max);

    TimingConfig tcfg;
    tcfg.exposure_time_s = cfg.timing().exposure_time_s;
    if (cfg.timing().n_exposures_override.has_value()) {
      tcfg.n_exposures_override = cfg.timing().n_exposures_override;
    }
    bld.SetTiming(cfg.experiment_cfg().livetime_days, cfg.experiment_cfg().duty_cycle, tcfg);

    DarkCurrentConfig dcc;
    dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    dcc.norm_scale                = cfg.backgrounds().norm_scale;
    bld.SetDarkCurrent(dcc);

    auto B_dc_ne = bld.BuildBkgAsimov();
    B_dc = fold.FoldNeBackground(*B_dc_ne);
  }

  // Flat background: same energy binning as the signal model spectrum
  // (cfg.model().Emin_eV/Emax_eV/nbins), matching the reference app's own
  // construction exactly (it borrows MakeSignalSpectrumE's binning rather
  // than BackgroundJSON's separate flat_bkg_Emin_eV/Emax_eV/nbins fields).
  const auto& mj = cfg.model();
  const double flat_rate_dru = cfg.backgrounds().has_flat_bkg ? cfg.backgrounds().flat_bkg_norm_per_kg_year : 0.0;
  std::vector<double> B_flat;  // flat-spectrum contribution per output bin
  if (!cfg.response().background_efficiency_csv.empty()) {
    // Background's own detection-efficiency curve, distinct from the
    // signal's kernel (see MakeClusterEnergyFlatBackgroundWithEfficiency
    // docstring / docs/ClusterFitMC_Design.md Sec. 6.9).
    const auto eff_table = LoadBackgroundEfficiencyTable(cfg.response().background_efficiency_csv);
    B_flat = MakeClusterEnergyFlatBackgroundWithEfficiency(fold, summary.exposure_kg_year, flat_rate_dru,
                                                            eff_table)
                 .B_pat;
  } else {
    auto dRdE_flat = MakeFlatDrdeSpectrum(flat_rate_dru, mj.Emin_eV, mj.Emax_eV, mj.nbins,
                                           "dRdE_flat_bkg_factory");
    B_flat = fold.Fold(*dRdE_flat, summary.exposure_kg_year);
  }

  // Low-energy excess (Tier A, fixed shape): additive on top of the flat
  // Compton term, gated by run.background.lee.enabled -- see
  // docs/LowEnergyExcess_Design.md. cluster_energy (WIMP-nucleon) only;
  // ErecoEdgesEV() throws a clear error if called for an analysis space
  // that doesn't support it, same contract MakeClusterEnergyFlatBackgroundWithEfficiency
  // already relies on.
  std::vector<double> B_lee;  // low-energy-excess contribution per output bin
  if (cfg.backgrounds().has_lee_bkg) {
    B_lee = MakeClusterEnergyLEEBackground(fold, summary.exposure_kg_year,
                                            cfg.backgrounds().lee_rate_per_kg_day,
                                            cfg.backgrounds().lee_decay_energy_eV)
                .B_pat;
  }

  const std::size_t n = fold.NumBins();
  BackgroundResult result;
  result.B_pat.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double b_dc = (i < B_dc.size()) ? B_dc[i] : 0.0;
    const double b_flat = (i < B_flat.size()) ? B_flat[i] : 0.0;
    const double b_lee = (i < B_lee.size()) ? B_lee[i] : 0.0;
    result.B_pat[i] = b_dc + b_flat + b_lee;
  }
  result.Bp.assign(n, 0.0);
  result.Br = result.B_pat;
  return result;
}

// ----------------------------------------------------------------------------
// MakeBpBrTemplateBackground
//   Background for background_source = "bp_br_template": I take Bp and Br
//   straight from run.background_Bp / run.background_Br and set B = Bp + Br.
//   Throws std::runtime_error if their lengths differ from the number of
//   response bins.
// ----------------------------------------------------------------------------
BackgroundResult MakeBpBrTemplateBackground(const ConfigManager& cfg, const ResponseFold& fold) {
  const auto& run = cfg.run();
  const std::size_t n = fold.NumBins();
  if (run.background_Bp.size() != n || run.background_Br.size() != n) {
    throw std::runtime_error(
        "BackgroundFactory: background_source is bp_br_template but background_Bp/Br size (" +
        std::to_string(run.background_Bp.size()) + "/" + std::to_string(run.background_Br.size()) +
        ") != response fold's bin count (" + std::to_string(n) + ")");
  }
  BackgroundResult result;
  result.Bp = run.background_Bp;
  result.Br = run.background_Br;
  result.B_pat.resize(n);
  for (std::size_t i = 0; i < n; ++i) result.B_pat[i] = result.Bp[i] + result.Br[i];
  return result;
}

}  // namespace

// ----------------------------------------------------------------------------
// MakeBackground
//   Public entry point: I pick the construction from run.background_source
//   ("bp_br_template" or, by default, "dc_flat_migration") and return the
//   nominal background B_pat plus its Bp / Br split.
// ----------------------------------------------------------------------------
BackgroundResult MakeBackground(const ConfigManager& cfg,
                                 const ExperimentSummary& summary,
                                 const ResponseFold& fold,
                                 const ResponseFactoryResult& response) {
  if (cfg.run().background_source == "bp_br_template") {
    return MakeBpBrTemplateBackground(cfg, fold);
  }
  return MakeDcFlatMigrationBackground(cfg, summary, fold, response);
}

// ----------------------------------------------------------------------------
// MakeClusterEnergyFlatBackground
//   Flat background for one cluster_energy channel: I build a flat dR/dE with
//   the channel's own rate (d.r.u.) between Emin_eV and Emax_eV and fold it
//   through fold.Fold(). Bp = 0, Br = B_pat.
// ----------------------------------------------------------------------------
BackgroundResult MakeClusterEnergyFlatBackground(const ResponseFold& fold,
                                                  double exposure_kg_year,
                                                  double Emin_eV, double Emax_eV, int nbins,
                                                  double flat_rate_dru) {
  auto dRdE_flat = MakeFlatDrdeSpectrum(flat_rate_dru, Emin_eV, Emax_eV, nbins,
                                         "dRdE_flat_bkg_channel");
  std::vector<double> B_flat = fold.Fold(*dRdE_flat, exposure_kg_year);

  const std::size_t n = fold.NumBins();
  BackgroundResult result;
  result.B_pat.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    result.B_pat[i] = (i < B_flat.size()) ? B_flat[i] : 0.0;
  }
  result.Bp.assign(n, 0.0);
  result.Br = result.B_pat;
  return result;
}

// ----------------------------------------------------------------------------
// MakeClusterEnergyFlatBackgroundWithEfficiency
//   Flat Compton background weighted by the background's own efficiency curve.
//   In every E_reco bin: B = rate[d.r.u.] * width[keV] * exposure * <eff>, where <eff>
//   is the bin average of the curve from Simpson's rule with 9 sub-points (the
//   curve rises steeply near threshold, so a bin-centre sample would be off).
// ----------------------------------------------------------------------------
BackgroundResult MakeClusterEnergyFlatBackgroundWithEfficiency(
    const ResponseFold& fold,
    double exposure_kg_year,
    double flat_rate_dru,
    const BackgroundEfficiencyTable& eff_table) {
  const auto edges_eV = fold.ErecoEdgesEV();
  const std::size_t n = fold.NumBins();

  BackgroundResult result;
  result.B_pat.assign(n, 0.0);
  for (std::size_t i = 0; i < n && i + 1 < edges_eV.size(); ++i) {
    const double e_lo_kev = edges_eV[i] / 1000.0;
    const double e_hi_kev = edges_eV[i + 1] / 1000.0;
    const double width_kev = e_hi_kev - e_lo_kev;
    // Bin-averaged efficiency via Simpson's rule (9 sub-points), not a
    // single point-sample at the bin center -- the digitized curve rises
    // very steeply near threshold (over just ~50 eV, one bin width here),
    // where a center point-sample can differ from the true bin average by
    // tens of percent (see docs/ClusterFitMC_Design.md Sec. 6.9 follow-up).
    constexpr int kSubPoints = 9;  // must be odd for Simpson's rule
    double sum = 0.0;
    for (int k = 0; k < kSubPoints; ++k) {
      const double e_k = e_lo_kev + (e_hi_kev - e_lo_kev) * static_cast<double>(k) / (kSubPoints - 1);
      const double w = (k == 0 || k == kSubPoints - 1) ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0);
      sum += w * EvaluateBackgroundEfficiency(eff_table, e_k);
    }
    const double eff_avg = sum / (3.0 * (kSubPoints - 1));
    result.B_pat[i] = flat_rate_dru * width_kev * exposure_kg_year * eff_avg;
  }
  result.Bp.assign(n, 0.0);
  result.Br = result.B_pat;
  return result;
}

// ----------------------------------------------------------------------------
// MakeClusterEnergyLEEBackground
//   Low-energy-excess background: dR/dE = rate * (1/eps) * exp(-E/eps). I
//   integrate it exactly over each E_reco bin,
//       B = rate[/kg/day] * exposure[kg*day] * (exp(-E_lo/eps) - exp(-E_hi/eps)),
//   so no numerical quadrature is needed.
// ----------------------------------------------------------------------------
BackgroundResult MakeClusterEnergyLEEBackground(const ResponseFold& fold,
                                                 double exposure_kg_year,
                                                 double rate_per_kg_day,
                                                 double epsilon_eV) {
  const auto edges_eV = fold.ErecoEdgesEV();
  const std::size_t n = fold.NumBins();
  const double exposure_kg_day = exposure_kg_year * 365.25;

  BackgroundResult result;
  result.B_pat.assign(n, 0.0);
  for (std::size_t i = 0; i < n && i + 1 < edges_eV.size(); ++i) {
    // Exact integral of the normalized exponential dR/dE = (1/eps)*exp(-E/eps)
    // over [e_lo, e_hi] -- closed form, no Simpson's-rule sampling needed
    // (contrast MakeClusterEnergyFlatBackgroundWithEfficiency's digitized
    // curve, which has no closed form).
    const double e_lo = edges_eV[i];
    const double e_hi = edges_eV[i + 1];
    result.B_pat[i] = rate_per_kg_day * exposure_kg_day *
                       (std::exp(-e_lo / epsilon_eV) - std::exp(-e_hi / epsilon_eV));
  }
  result.Bp.assign(n, 0.0);
  result.Br = result.B_pat;
  return result;
}

}  // namespace ccdarksens
