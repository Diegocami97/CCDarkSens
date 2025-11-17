#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/detector/Detector.hh"

#include "ccdarksens/rates/RateSource.hh"
#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/backgrounds/BackgroundBuilder.hh"

#include <TH1D.h>
#include <TFile.h>
#include <filesystem>
#include <iostream>
#include <vector>
#include <algorithm>

static std::unique_ptr<TH1D> MakeFlatEfficiency(int ne_min, int ne_max, double eps) {
  const int nbin = ne_max - ne_min + 1;
  std::vector<double> edges(nbin + 1);
  for (int i=0;i<=nbin;++i) edges[i] = (ne_min - 0.5) + i;
  auto h = std::make_unique<TH1D>("eps_ne","Pattern efficiency; n_{e}; #epsilon", nbin, edges.data());
  for (int b=1;b<=nbin;++b) h->SetBinContent(b, std::clamp(eps,0.0,1.0));
  return h;
}

static double SumROI(const TH1D& h, const std::vector<int>& roi, int ne_min) {
   double s = 0.0;
   for (int ne : roi) {
     const int b = (ne - ne_min + 1);  // 1-based ROOT bin index
     if (b >= 1 && b <= h.GetNbinsX()) s += h.GetBinContent(b);
   }
   return s;
}

int main(int argc, char** argv){
  if (argc<2){ std::cerr<<"usage: ccdarksens_sensitivity_fast <config.json>\n"; return 1; }

  try {
    ccdarksens::ConfigManager cfg(argv[1]); cfg.parse();
    const auto& det = cfg.detector();
    const double mass_kg = det.mass_kg();
    ccdarksens::ExperimentSetup setup(cfg.experiment_cfg(), mass_kg, cfg.run().rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[FAST] Run: " << cfg.run().label << "\n"
              << "  Mass [kg]: " << mass_kg << "\n"
              << "  Exposure [kg·day]: " << summary.exposure_kg_year << "\n"
              << "  Duty cycle: " << cfg.experiment_cfg().duty_cycle << "\n"
              << "  Binning n_e: [" << ne_min << "," << ne_max << "]\n";

    // --- dR/dE (synthetic test) ---
    auto dRdE = ccdarksens::RateSource::MakeLinearEnergyHist(10.0, 100.0, 180, "dRdE_flat");
    ccdarksens::RateSource::FillFlat(*dRdE, 1.0);

    // --- Response modules (FAST path) ---
    auto ion = std::make_shared<ccdarksens::ChargeIonization>("data/p100K_table.csv");
    auto diff = std::make_shared<ccdarksens::Diffusion>(
      803.25, 6.5e-4, 1.0, 0.0, 670.0, 0.08, 0.16
    );

    std::unique_ptr<TH1D> eps_hist;
    if (cfg.backgrounds().has_flat_eps) eps_hist = MakeFlatEfficiency(ne_min, ne_max, cfg.backgrounds().flat_eps);
    else                                eps_hist = MakeFlatEfficiency(ne_min, ne_max, 1.0);

    auto pe = std::make_shared<ccdarksens::PatternEfficiency>();
    pe->SetEfficiencyHist(*eps_hist);

    ccdarksens::DetectorResponsePipeline pipe(ion);
    pipe.SetDiffusion(diff);
    pipe.SetPatternEfficiency(pe);

    auto S_raw = ion->FoldToNe(*dRdE, summary.exposure_kg_year, ne_min, ne_max);
    auto S_obs = pipe.Apply(*dRdE, summary.exposure_kg_year, ne_min, ne_max, 50.0);

    // --- Background Asimov (use same ε) ---
    ccdarksens::BackgroundBuilder bld(det.geometry().rows, det.geometry().cols, det.geometry().active_fraction, ne_min, ne_max);
    ccdarksens::TimingConfig tcfg;
    tcfg.exposure_time_s = cfg.timing().exposure_time_s;
    if (cfg.timing().n_exposures_override.has_value()) tcfg.n_exposures_override = cfg.timing().n_exposures_override;
    bld.SetTiming(cfg.experiment_cfg().livetime_days, cfg.experiment_cfg().duty_cycle, tcfg);

    ccdarksens::DarkCurrentConfig dcc;
    // Dark current is now configured in e-/pixel/year
    dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    dcc.norm_scale = cfg.backgrounds().norm_scale;
    bld.SetDarkCurrent(dcc);

    auto pe_bkg = std::make_shared<ccdarksens::PatternEfficiency>();
    pe_bkg->SetEfficiencyHist(*eps_hist);
    bld.SetPatternEfficiency(pe_bkg);

    auto B_asimov = bld.BuildBkgAsimov();
    // --- Background rollups (consistent across paths) ---
    const double B_total_pixexp = static_cast<double>(bld.NActivePixels()) * static_cast<double>(bld.NExposuresUsed());

    // n_e = 0 bin (after epsilon)
    int b0 = (0 - ne_min + 1);
    double B_zero_after_eps = 0.0;
    if (b0 >= 1 && b0 <= B_asimov->GetNbinsX()) {
    B_zero_after_eps = B_asimov->GetBinContent(b0);
    }

    // non-zero after eps (n_e >= 1)
    int b1 = (std::max(1, ne_min) - ne_min + 1);
    int bmax = B_asimov->GetNbinsX();
    double B_nonzero_after_eps = (b1 <= bmax) ? B_asimov->Integral(b1, bmax) : 0.0;

    std::cout << "  B_total_pixexp (no eps): " << B_total_pixexp << "\n"
            << "  B_zero_after_eps       : " << B_zero_after_eps << "\n"
            << "  B_nonzero_after_eps    : " << B_nonzero_after_eps << "\n";


    // --- ROI yields ---
    const auto roi_bins = cfg.experiment_cfg().roi_bins; // vector<int>
    const double S_ROI = SumROI(*S_obs, roi_bins, ne_min);
    const double B_ROI = SumROI(*B_asimov, roi_bins, ne_min);

    std::cout << "\n[checks]\n"
              << "  ∫ dE dR/dE (per kg·day): " << dRdE->Integral("width") << "\n"
              << "  S_raw integral: " << S_raw->Integral() << "\n"
              << "  S_obs integral: " << S_obs->Integral() << "\n"
              << "  B_asimov integral: " << B_asimov->Integral() << "\n"
              << "  ROI bins: ";
    for (size_t i=0;i<roi_bins.size();++i){ std::cout << roi_bins[i] << (i+1<roi_bins.size()? ",":""); }
    std::cout << "\n  S_ROI: " << S_ROI << "   B_ROI: " << B_ROI << "\n";

    std::filesystem::create_directories(cfg.run().outdir);
    TFile fout((cfg.run().outdir + "/fast_check.root").c_str(),"RECREATE");
    dRdE->Write(); S_raw->Write("S_ne_raw"); S_obs->Write("S_ne_obs"); B_asimov->Write("B_asimov_ne");
    fout.Close();
    std::cout << "\n[output] " << (cfg.run().outdir + "/fast_check.root") << "\n";
    return 0;
  } catch (const std::exception& e){
    std::cerr << "ERROR: " << e.what() << "\n"; return 2;
  }
}
