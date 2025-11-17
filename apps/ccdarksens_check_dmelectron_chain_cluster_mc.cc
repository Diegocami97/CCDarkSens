#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <sstream>

#include <TH1D.h>
#include <TFile.h>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/model/DMElectronModel.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/ClusterMC.hh"
#include "ccdarksens/backgrounds/BackgroundBuilder.hh"

using namespace ccdarksens;

// Simple helper to sum over ROI bins (same pattern as sensitivity apps)
static double SumROI(TH1D& h,
                     const std::vector<int>& roi_bins,
                     int ne_min)
{
  double s = 0.0;
  for (int ne : roi_bins) {
    const int bin = h.FindBin(ne);
    s += h.GetBinContent(bin);
  }
  return s;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json\n";
    return 1;
  }

  try {
    const std::string config_path = argv[1];
    ConfigManager cfg(config_path);
    cfg.parse();

    const auto& run = cfg.run();
    std::cout << "[run] label=" << run.label
              << " outdir='" << run.outdir << "'\n";

    // --- Experiment setup (consistent with other apps) ---
    const auto& det      = cfg.detector();
    const double mass_kg = det.mass_kg();

    ExperimentSetup setup(cfg.experiment_cfg(), mass_kg, run.rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[CHECK-CHAIN] Run: " << run.label << "\n"
              << "  Mass [kg]: " << mass_kg << "\n"
              << "  Exposure [kg·year]: " << summary.exposure_kg_year << "\n"
              << "  Duty cycle: " << cfg.experiment_cfg().duty_cycle << "\n"
              << "  Binning n_e: [" << ne_min << "," << ne_max << "]\n";

    // ------------------------------------------------------------------
    // 1) Load DM–electron dR/dE from QE-Dark (via DMElectronModel)
    //    Units: events / kg / year / eV
    // ------------------------------------------------------------------
    const auto& mj = cfg.model();
    DMElectronConfig mc;
    mc.material          = mj.material;
    mc.mediator          = mj.mediator;
    mc.rates_dir         = mj.rates_dir;
    mc.filename_template = mj.filename_template;
    mc.mchi_MeV          = mj.mchi_MeV;
    mc.sigma_e_cm2       = mj.sigma_e_cm2;
    mc.Emin_eV           = mj.Emin_eV;
    mc.Emax_eV           = mj.Emax_eV;
    mc.nbins             = mj.nbins;

    DMElectronModel dm_model;
    if (!dm_model.Configure(mc)) {
      std::cerr << "[error] DMElectronModel::Configure failed (could not load QE-Dark table)\n";
      return 1;
    }

    auto dRdE = dm_model.MakeSpectrum_E();  // events / kg / year / eV
    if (!dRdE) {
      std::cerr << "[error] DMElectronModel::MakeSpectrum_E returned null\n";
      return 1;
    }

    std::cout << "  ∫ dE dR/dE (per kg·year): "
              << dRdE->Integral("width") << "\n";

    // ------------------------------------------------------------------
    // 2) Detector response: ionization + diffusion + ε_MC from ClusterMC
    // ------------------------------------------------------------------
    auto ion = std::make_shared<ChargeIonization>("data/p100K_table.csv");

    auto diff = std::make_shared<Diffusion>(
      /*A_um2      */ 803.25,
      /*b_umInv    */ 6.5e-4,
      /*alpha      */ 1.0,
      /*beta/keV   */ 0.0,
      /*thickness  */ 670.0,   // um
      /*sigmaR_e   */ 0.08,
      /*sigmaC_e   */ 0.16
    );

    // Geometry for ClusterMC (same pattern as mvp_dm_electron_min)
    GeometryLite geo;
    geo.rows            = det.geometry().rows;
    geo.cols            = det.geometry().cols;
    geo.pixel_size_um   = det.geometry().pixel_size_um;
    geo.thickness_um    = det.geometry().thickness_mm * 1000.0;
    geo.active_fraction = det.geometry().active_fraction;

    ClusterMCConfig cmc;
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

    // λ from config is in e-/pix/year; convert ONLY for pileup in ClusterMC:
    {
      const double lambda_year = cfg.backgrounds().lambda_e_per_pix_per_year; // e-/pix/year
      const double year_s      = 365.25 * 86400.0;
      const double exp_s       = cfg.timing().exposure_time_s;
      double lambda_exp = 0.0; // e-/pix/exposure

      if (exp_s > 0.0) {
        lambda_exp = lambda_year * (exp_s / year_s);
      }

      cmc.lambda_dc_per_pix_per_exposure =
        cmc.pileup_with_dc ? lambda_exp : 0.0;
    }

    cmc.half_window_pix = 2;
    cmc.enable_MN  = true;
    cmc.enable_MNL = true;

    ClusterMC clusterMC(cmc, geo);
    auto eps_mc = clusterMC.PrecomputeEpsilon(ne_min, ne_max);

    auto pe = std::make_shared<PatternEfficiency>();
    pe->SetEfficiencyHist(*eps_mc);

    DetectorResponsePipeline pipe(ion);
    pipe.SetDiffusion(diff);
    pipe.SetPatternEfficiency(pe);

    // S_raw and S_obs: counts per n_e bin for this exposure
    auto S_raw = ion->FoldToNe(*dRdE, summary.exposure_kg_year,
                               ne_min, ne_max);
    auto S_obs = pipe.Apply(*dRdE, summary.exposure_kg_year,
                            ne_min, ne_max, /*Emax_keV_for_MC*/ 50.0);
                            
    // -------------------------------------------------
    // (2b) Flat energy background: 1 event/(kg year keV)
    // -------------------------------------------------

    // Define the flat level (events / (kg · year · keV)).
    // You can later move this into JSON; for now it's hard-coded.
    const double flat_rate_per_kg_year_keV = cfg.backgrounds().has_flat_bkg ?
        cfg.backgrounds().flat_bkg_norm_per_kg_year : 0.0;

    // Convert to events / (kg · year · eV):
    const double flat_rate_per_kg_year_eV =
        flat_rate_per_kg_year_keV / 1000.0;  // 1 keV = 1000 eV

    // Build a flat dR/dE histogram with the SAME binning as the DM spectrum.
    auto dR_flat =
        std::unique_ptr<TH1D>(static_cast<TH1D*>(dRdE->Clone("dRdE_flat")));
    dR_flat->Reset("ICES");

    const int nbins_E = dR_flat->GetNbinsX();
    for (int ib = 1; ib <= nbins_E; ++ib) {
      dR_flat->SetBinContent(ib, flat_rate_per_kg_year_eV);
    }

    // Fold flat dR/dE -> n_e using the SAME chain as the signal:
    auto Bflat_raw = ion->FoldToNe(*dR_flat, summary.exposure_kg_year, ne_min, ne_max);
    auto Bflat_obs = pipe.Apply(*dR_flat, summary.exposure_kg_year,
                                ne_min, ne_max,
                                /*Emax_keV_for_MC*/ 50.0);

    // ------------------------------------------------------------------
    // 3) Dark-current Asimov background using BackgroundBuilder
    //    Same ε_MC efficiency histogram
    // ------------------------------------------------------------------
    BackgroundBuilder bld(det.geometry().rows,
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
    // Here λ is interpreted as e-/pix/year and converted internally by BackgroundBuilder
    dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    dcc.norm_scale = cfg.backgrounds().norm_scale;
    bld.SetDarkCurrent(dcc);

    auto pe_bkg = std::make_shared<PatternEfficiency>();
    pe_bkg->SetEfficiencyHist(*eps_mc);
    bld.SetPatternEfficiency(pe_bkg);

    auto B_dc_asimov = bld.BuildBkgAsimov();

    // Total Asimov background: dark current + flat energy
    // Combine DC Asimov background with the flat component
    auto B_asimov = std::unique_ptr<TH1D>(
        static_cast<TH1D*>(B_dc_asimov->Clone("B_asimov_total")));
    B_asimov->Add(Bflat_obs.get());

    const double B_total_pixexp =
      static_cast<double>(bld.NActivePixels()) *
      static_cast<double>(bld.NExposuresUsed());

    // n_e = 0 bin after ε
    int b0 = (0 - ne_min + 1);
    double B_zero_after_eps = 0.0;
    if (b0 >= 1 && b0 <= B_asimov->GetNbinsX()) {
      B_zero_after_eps = B_asimov->GetBinContent(b0);
    }

    // n_e >= 1 bins after ε
    int b1   = (std::max(1, ne_min) - ne_min + 1);
    int bmax = B_asimov->GetNbinsX();
    double B_nonzero_after_eps =
      (b1 <= bmax) ? B_asimov->Integral(b1, bmax) : 0.0;

    std::cout << "  B_total_pixexp (no eps): " << B_total_pixexp << "\n"
              << "  B_zero_after_eps       : " << B_zero_after_eps << "\n"
              << "  B_nonzero_after_eps    : " << B_nonzero_after_eps << "\n";

    // ------------------------------------------------------------------
    // 4) ROI yields and basic checks
    // ------------------------------------------------------------------
    // const auto& roi_bins = cfg.experiment_cfg().roi_bins;
    // const double S_ROI = SumROI(*S_obs, *const_cast<std::vector<int>*>(&roi_bins), ne_min);
    // const double B_ROI = SumROI(*B_asimov, *const_cast<std::vector<int>*>(&roi_bins), ne_min);
    const auto& roi_bins = cfg.experiment_cfg().roi_bins;
    const double S_ROI      = SumROI(*S_obs,       roi_bins, ne_min);
    const double Bdc_ROI    = SumROI(*B_dc_asimov, roi_bins, ne_min);
    const double Bflat_ROI  = SumROI(*Bflat_obs,   roi_bins, ne_min);
    const double Btot_ROI   = SumROI(*B_asimov,    roi_bins, ne_min);

    std::ostringstream roi_str;
    for (size_t i = 0; i < roi_bins.size(); ++i) {
      if (i) roi_str << ",";
      roi_str << roi_bins[i];
    }

    std::cout << "\n[checks]\n"
              << "  ∫ dE dR/dE (per kg·year): " << dRdE->Integral("width") << "\n"
              << "  S_raw integral: " << S_raw->Integral() << "\n"
              << "  S_obs integral: " << S_obs->Integral() << "\n"
              << "  B_dc_asimov integral: " << B_dc_asimov->Integral() << "\n"
              << "  B_flat integral (after response): " << Bflat_obs->Integral() << "\n"
              << "  B_total integral: " << B_asimov->Integral() << "\n"
              << "  ROI bins: " << roi_str.str() << "\n"
              << "  S_ROI: " << S_ROI
              << "   Bdc_ROI: " << Bdc_ROI
              << "   Bflat_ROI: " << Bflat_ROI
              << "   Btot_ROI: " << Btot_ROI
              << "\n";

    // for (size_t i = 0; i < roi_bins.size(); ++i) {
    //   std::cout << roi_bins[i] << (i + 1 < roi_bins.size() ? "," : "");
    // }
    // std::cout << "\n  S_ROI: " << S_ROI << "   B_ROI: " << B_ROI << "\n";

    // ------------------------------------------------------------------
    // 5) Write everything to a ROOT file in run.outdir
    // ------------------------------------------------------------------
    if (!run.outdir.empty()) {
      std::filesystem::create_directories(run.outdir);
      const std::string fname = run.outdir + "/check_dmelectron_chain_cluster_mc.root";

      TFile fout(fname.c_str(), "RECREATE");
      dRdE->Write("dRdE_E");
      dR_flat->Write("dRdE_flat_E");
      eps_mc->Write("eps_mc");
      S_raw->Write("S_ne_raw");
      S_obs->Write("S_ne_obs");
      B_dc_asimov->Write("B_dc_ne");
      Bflat_obs->Write("B_flat_ne");
      B_asimov->Write("B_asimov_total_ne");
      fout.Close();

      std::cout << "[output] wrote " << fname << "\n";
    } else {
      std::cout << "[output] run.outdir is empty; not writing ROOT file.\n";
    }

  } catch (const std::exception& e) {
    std::cerr << "[fatal] exception: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
