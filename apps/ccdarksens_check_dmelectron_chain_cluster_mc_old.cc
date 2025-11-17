// apps/ccdarksens_check_dmelectron_chain_cluster_mc.cc
#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/detector/Detector.hh"

#include "ccdarksens/model/DMElectronModel.hh"
#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/response/ClusterMC.hh"
#include "ccdarksens/backgrounds/BackgroundBuilder.hh"

#include <TH1D.h>
#include <TFile.h>
#include <filesystem>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace ccdarksens;

// helper to sum over ROI in n_e histogram
static double SumROI(const TH1D& h, const std::vector<int>& roi, int ne_min) {
  double s = 0.0;
  for (int ne : roi) {
    const int b = (ne - ne_min + 1);  // 1-based ROOT bin index
    if (b >= 1 && b <= h.GetNbinsX()) s += h.GetBinContent(b);
  }
  return s;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: ccdarksens_check_dmelectron_chain_cluster_mc <config.json>\n";
    return 1;
  }

  try {
    ccdarksens::ConfigManager cfg(argv[1]);
    cfg.parse();

    // --- Experiment setup (same as sensitivity_cluster_mc) ---
    const auto& det = cfg.detector();
    const double mass_kg = det.mass_kg();
    ccdarksens::ExperimentSetup setup(cfg.experiment_cfg(), mass_kg, cfg.run().rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[CHECK-CHAIN] Run: " << cfg.run().label << "\n"
              << "  Mass [kg]: " << mass_kg << "\n"
              << "  Exposure [kg·year]: " << summary.exposure_kg_year << "\n"
              << "  Duty cycle: " << cfg.experiment_cfg().duty_cycle << "\n"
              << "  Binning n_e: [" << ne_min << "," << ne_max << "]\n";

    // --- 1) Load DM–electron dR/dE from QE-Dark (via DMElectronModel) ---
    const auto& mj = cfg.model();
    ccdarksens::DMElectronModel dm_model;
    ccdarksens::DMElectronConfig mc;
    mc.material          = mj.material;
    mc.mediator          = mj.mediator;
    mc.rates_dir         = mj.rates_dir;
    mc.filename_template = mj.filename_template;
    mc.mchi_MeV          = mj.mchi_MeV;
    mc.sigma_e_cm2       = mj.sigma_e_cm2;
    mc.Emin_eV           = mj.Emin_eV;
    mc.Emax_eV           = mj.Emax_eV;
    mc.nbins             = mj.nbins;

    if (!dm_model.Configure(mc)) {
      std::cerr << "[error] DMElectronModel::Configure failed (could not load QE-Dark table)\n";
      return 1;
    }

    // dRdE in units events / kg / year / eV
    auto dRdE = dm_model.MakeSpectrum_E();
    if (!dRdE) {
      std::cerr << "[error] DMElectronModel::MakeSpectrum_E returned null\n";
      return 1;
    }

    std::cout << "  ∫ dE dR/dE (per kg·year): " << dRdE->Integral("width") << "\n";

    // --- 2) Detector response: ionization + diffusion + ε_MC ---
    auto ion = std::make_shared<ccdarksens::ChargeIonization>("data/p100K_table.csv");

    auto diff = std::make_shared<ccdarksens::Diffusion>(
      /*A_um2      */ 803.25,
      /*b_umInv    */ 6.5e-4,
      /*alpha      */ 1.0,
      /*beta/keV   */ 0.0,
      /*thickness  */ 670.0,   // um
      /*sigmaR_e   */ 0.08,
      /*sigmaC_e   */ 0.16
    );

    // Geometry & ClusterMC config from JSON
    ccdarksens::GeometryLite geo;
    geo.rows           = det.geometry().rows;
    geo.cols           = det.geometry().cols;
    geo.pixel_size_um  = det.geometry().pixel_size_um;
    geo.thickness_um   = det.geometry().thickness_mm * 1000.0;
    geo.active_fraction= det.geometry().active_fraction;

    // ClusterMC config from JSON
    ccdarksens::ClusterMCConfig cmc;
    cmc.n_events_per_ne = cfg.response().cmc.n_events_per_ne;
    cmc.rng_seed        = cfg.response().cmc.rng_seed;
    cmc.sigma_readout_e = cfg.response().cmc.sigma_readout_e;
    cmc.Qmin_e          = cfg.response().cmc.Qmin_e;
    cmc.Qmax_e          = cfg.response().cmc.Qmax_e;
    cmc.diffusion = {
      cfg.response().cmc.A_um2,
      cfg.response().cmc.b_umInv,
      cfg.response().cmc.alpha,
      cfg.response().cmc.beta_per_keV
    };
    cmc.rows_bin        = cfg.response().cmc.rows_bin;
    cmc.cols_bin        = cfg.response().cmc.cols_bin;
    cmc.pileup_with_dc  = cfg.response().cmc.pileup_with_dc;

    // Convert λ_year → λ_per_exposure for DC pileup
    {
      const double lambda_year = cfg.backgrounds().lambda_e_per_pix_per_year;
      const double year_s      = 365.25 * 86400.0;
      const double exp_s       = cfg.timing().exposure_time_s;
      double lambda_exp = 0.0;
      if (exp_s > 0.0) {
        lambda_exp = lambda_year * (exp_s / year_s);
      }
      cmc.lambda_dc_per_pix_per_exposure =
        cmc.pileup_with_dc ? lambda_exp : 0.0;
    }
    cmc.half_window_pix = 2;
    cmc.enable_MN  = true;
    cmc.enable_MNL = true;

    // ClusterMC & ε_MC
    ccdarksens::ClusterMC clusterMC(cmc, geo);
    auto eps_mc = clusterMC.PrecomputeEpsilon(ne_min, ne_max);

    auto pe = std::make_shared<ccdarksens::PatternEfficiency>();
    pe->SetEfficiencyHist(*eps_mc);

    ccdarksens::DetectorResponsePipeline pipe(ion);
    pipe.SetDiffusion(diff);
    pipe.SetPatternEfficiency(pe);

    // S_raw, S_obs in counts per n_e bin for this exposure
    auto S_raw = ion->FoldToNe(*dRdE, summary.exposure_kg_year, ne_min, ne_max);
    auto S_obs = pipe.Apply(*dRdE, summary.exposure_kg_year, ne_min, ne_max, 50.0);

    // --- 3) Background Asimov: pure dark current, same ε_MC ---
    ccdarksens::BackgroundBuilder bld(det.geometry().rows,
                          det.geometry().cols,
                          det.geometry().active_fraction,
                          ne_min, ne_max);

    TimingConfig tcfg;
    tcfg.exposure_time_s = cfg.timing().exposure_time_s;
    if (cfg.timing().n_exposures_override.has_value())
      tcfg.n_exposures_override = cfg.timing().n_exposures_override;
    bld.SetTiming(cfg.experiment_cfg().livetime_days,
                  cfg.experiment_cfg().duty_cycle,
                  tcfg);

    DarkCurrentConfig dcc;
    // store per-year λ here; BackgroundBuilder converts internally
    dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    dcc.norm_scale = cfg.backgrounds().norm_scale;
    bld.SetDarkCurrent(dcc);

    // pattern-efficiency for background, etc. unchanged
    auto pe_bkg = std::make_shared<PatternEfficiency>();
    pe_bkg->SetEfficiencyHist(*eps_mc);
    bld.SetPatternEfficiency(pe_bkg);

    auto B_asimov = bld.BuildBkgAsimov();

    const double B_total_pixexp =
      static_cast<double>(bld.NActivePixels()) *
      static_cast<double>(bld.NExposuresUsed());

    // n_e = 0 bin after ε
    int b0 = (0 - ne_min + 1);
    double B_zero_after_eps = 0.0;
    if (b0 >= 1 && b0 <= B_asimov->GetNbinsX()) {
      B_zero_after_eps = B_asimov->GetBinContent(b0);
    }

    // n_e >= 1 (non-zero after ε)
    int b1   = (std::max(1, ne_min) - ne_min + 1);
    int bmax = B_asimov->GetNbinsX();
    double B_nonzero_after_eps =
      (b1 <= bmax) ? B_asimov->Integral(b1, bmax) : 0.0;

    std::cout << "  B_total_pixexp (no eps): " << B_total_pixexp << "\n"
              << "  B_zero_after_eps       : " << B_zero_after_eps << "\n"
              << "  B_nonzero_after_eps    : " << B_nonzero_after_eps << "\n";

    // --- 4) ROI yields (signal + background) ---
    const auto roi_bins = cfg.experiment_cfg().roi_bins;
    const double S_ROI = SumROI(*S_obs,    roi_bins, ne_min);
    const double B_ROI = SumROI(*B_asimov, roi_bins, ne_min);

    std::cout << "\n[checks]\n"
              << "  ∫ dE dR/dE (per kg·year): " << dRdE->Integral("width") << "\n"
              << "  S_raw integral: " << S_raw->Integral() << "\n"
              << "  S_obs integral: " << S_obs->Integral() << "\n"
              << "  B_asimov integral: " << B_asimov->Integral() << "\n"
              << "  ROI bins: ";
    for (size_t i = 0; i < roi_bins.size(); ++i) {
      std::cout << roi_bins[i]
                << (i + 1 < roi_bins.size() ? "," : "");
    }
    std::cout << "\n  S_ROI: " << S_ROI
              << "   B_ROI: " << B_ROI << "\n";

    // --- 5) Output ROOT file with all pieces ---
    std::filesystem::create_directories(cfg.run().outdir);
    const std::string fname = cfg.run().outdir + "/check_dmelectron_chain_cluster_mc.root";
    TFile fout(fname.c_str(), "RECREATE");
    dRdE->Write("dRdE_E");
    eps_mc->Write("eps_mc");
    S_raw->Write("S_ne_raw");
    S_obs->Write("S_ne_obs");
    B_asimov->Write("B_asimov_ne");
    fout.Close();

    std::cout << "\n[output] " << fname << "\n";
    return 0;

  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 2;
  }
}
