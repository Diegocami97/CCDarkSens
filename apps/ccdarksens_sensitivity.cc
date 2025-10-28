// apps/ccdarksens_sensitivity.cc
#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"

#include "ccdarksens/rates/RateSource.hh"
#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/backgrounds/PoissonDarkCurrent.hh"

#include <TH1D.h>
#include <TFile.h>
#include "TSystem.h"
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

// helper: build a flat epsilon(ne) TH1D aligned to [ne_min, ne_max] integer bins
static std::unique_ptr<TH1D> MakeFlatEfficiency(int ne_min, int ne_max, double epsilon) {
  const int nbin = ne_max - ne_min + 1;
  std::vector<double> edges(nbin + 1);
  for (int i = 0; i <= nbin; ++i) edges[i] = (ne_min - 0.5) + i;

  auto h = std::make_unique<TH1D>("eps_ne", "Pattern efficiency; n_{e}; #epsilon", nbin, edges.data());
  for (int b = 1; b <= nbin; ++b) h->SetBinContent(b, std::clamp(epsilon, 0.0, 1.0));
  return h;
}

// helper: print first N bins of a TH1D with integer centers
static void PrintFirstBins(const TH1D& h, int n_to_print = 10) {
  const int nb = h.GetNbinsX();
  const int n = std::min(nb, n_to_print);
  std::cout << "    n_e : counts\n";
  for (int i = 1; i <= n; ++i) {
    const int ne = static_cast<int>(std::lround(h.GetBinCenter(i)));
    std::cout << "    " << ne << " : " << h.GetBinContent(i) << "\n";
  }
  if (nb > n) std::cout << "    ... (" << (nb - n) << " more bins)\n";
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "usage: ccdarksens_sensitivity <config.json>\n";
    return 1;
  }

  try {
    // -------------------- Load config & prepare exposure/binning --------------------
    ccdarksens::ConfigManager cfg(argv[1]);
    cfg.parse();

    const auto& det = cfg.detector();
    const double mass_kg = det.mass_kg();

    ccdarksens::ExperimentSetup setup(cfg.experiment_cfg(), mass_kg, cfg.run().rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[ccdarksens] Run: " << cfg.run().label << "\n"
              << "  Mass [kg]: " << mass_kg << "\n"
              << "  Exposure [kg·day]: " << summary.exposure_kg_day << "\n"
              << "  Mode: " << summary.mode_string << "\n"
              << "  Binning n_e: [" << ne_min << ", " << ne_max << "]\n";

    // -------------------- Build a synthetic dR/dE --------------------
    // Example 1: flat spectrum with integral = 1 event/(kg·day) in [10,100] eV
    auto dRdE = ccdarksens::RateSource::MakeLinearEnergyHist(10.0, 100.0, 180, "dRdE_flat");
    ccdarksens::RateSource::FillFlat(*dRdE, /*integral*/ 1.0);

    // (Optionally, try a mono line)
    // auto dRdE = ccdarksens::RateSource::MakeLinearEnergyHist(10.0, 100.0, 180, "dRdE_line");
    // ccdarksens::RateSource::FillMonoLine(*dRdE, /*E0_eV*/ 30.0, /*integral*/ 1.0);

    // -------------------- Response modules --------------------
    // Charge-ionization table: CSV with columns Er_eV, P1, P2, ..., Pk
    // Place your table at data/p100K_table.csv (or adjust the path)
    auto ion = std::make_shared<ccdarksens::ChargeIonization>("data/p100K_table.csv");

    // Diffusion: parameters reflecting your Python model (tune kappa/readout as needed)
    auto diff = std::make_shared<ccdarksens::Diffusion>(
      /*A_um2*/          803.25,
      /*b_umInv*/        6.5e-4,
      /*alpha*/          1.0,
      /*beta_per_keV*/   0.0,      // set >0 only if you need E-dependence in sigma_xy
      /*thickness_um*/   670.0,    // 0.67 mm
      /*kappa_e_per_um*/ 0.08,     // mapping µm lateral spread -> e- smearing (tune!)
      /*sigma_readout_e*/0.16      // readout noise in e-
    );

    // Pattern efficiency: start with a flat 95% efficiency (you can later load from file)
    auto eps = MakeFlatEfficiency(ne_min, ne_max, 0.95);
    auto pe  = std::make_shared<ccdarksens::PatternEfficiency>();
    pe->SetEfficiencyHist(*eps);

    // Pipeline
    ccdarksens::DetectorResponsePipeline pipe(ion);
    pipe.SetDiffusion(diff);
    pipe.SetPatternEfficiency(pe);

    // Apply response. Ee_ref_eV is the characteristic energy for diffusion σ; with beta=0 it's irrelevant.
    double Ee_ref_eV = 50.0;
    auto S_raw = ion->FoldToNe(*dRdE, summary.exposure_kg_day, ne_min, ne_max);
    auto S_obs = pipe.Apply(*dRdE, summary.exposure_kg_day, ne_min, ne_max, Ee_ref_eV);

    // -------------------- Background: Poisson dark current --------------------
    // Mean lambda in electrons; norm is an overall scaling (set 1.0 for unit-normalized PMF)
    ccdarksens::PoissonDarkCurrentBackground bkg(/*lambda_e*/ 1.5e-4, /*norm*/ 1.0);
    auto B = bkg.MakeHist(ne_min, ne_max);

    // -------------------- Sanity prints --------------------
    const double sum_dRdE = dRdE->Integral("width"); // ∫ dE dR/dE == 1.0
    const double S_raw_sum = S_raw->Integral();
    const double S_obs_sum = S_obs->Integral();
    const double B_sum     = B->Integral();

    std::cout << "\n[checks]\n";
    std::cout << "  ∫ dE dR/dE (per kg·day):  " << sum_dRdE << "\n";
    std::cout << "  Signal before diffusion/eff (counts): " << S_raw_sum << "\n";
    std::cout << "  Signal after  diffusion+eff (counts): " << S_obs_sum << "\n";
    std::cout << "  Background Poisson sum (norm):        " << B_sum << "\n";

    std::cout << "\n[first bins of S_raw]\n";
    PrintFirstBins(*S_raw, 12);
    std::cout << "\n[first bins of S_obs]\n";
    PrintFirstBins(*S_obs, 12);
    std::cout << "\n[first bins of B]\n";
    PrintFirstBins(*B, 12);

    // -------------------- Save to ROOT file --------------------
    const std::string outroot = cfg.run().outdir + "/mvp_check.root";
    gSystem->mkdir(cfg.run().outdir.c_str(), kTRUE); // ensure outdir exists
    TFile fout(outroot.c_str(), "RECREATE");
    dRdE->Write();
    S_raw->Write("S_ne_raw");
    S_obs->Write("S_ne_obs");
    B->Write("B_ne");
    fout.Close();
    std::cout << "\n[output] wrote histograms to: " << outroot << "\n";

    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 2;
  }
}
