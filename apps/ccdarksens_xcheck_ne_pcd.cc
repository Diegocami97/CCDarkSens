// ============================================================================
//  CCDarkSens — ccdarksens_xcheck_ne_pcd
//  Cross-check executable that validates PCD→n_e→pattern folding by writing S_true, kernels, S_rec, and S_obs histograms to a debug ROOT file.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <map>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <nlohmann/json.hpp>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/model/DMElectronModel.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"

#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/PixelSimulator.hh"
#include "ccdarksens/response/PatternClassifier.hh"
#include "ccdarksens/response/EfficiencyMC.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/response/PCDBasedResponse.hh"
#include "ccdarksens/response/PCDCalculator.hh"

using nlohmann::json;
using namespace ccdarksens;

// Simple helper to get first value from grid spec
static double first_grid_value(const json& spec)
{
  if (spec.is_object()) {
    if (spec.contains("values") && spec["values"].is_array()
        && !spec["values"].empty()) {
      return spec["values"][0].get<double>();
    }
    if (spec.contains("linspace")) {
      return spec["linspace"].at("start").get<double>();
    }
    if (spec.contains("logspace")) {
      double a = spec["logspace"].at("start_exp").get<double>();
      return std::pow(10.0, a);
    }
  } else if (spec.is_array() && !spec.empty()) {
    return spec[0].get<double>();
  }
  throw std::runtime_error("first_grid_value: empty or invalid spec");
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json\n";
    return 1;
  }

  const std::string config_path = argv[1];

  try {
    // -------------------------------------------------------------------------
    // Parse configuration (same pattern as main scan app)
    // -------------------------------------------------------------------------
    ConfigManager cfg(config_path);
    cfg.parse();

    json jroot;
    {
      std::ifstream jf(config_path);
      if (!jf) {
        throw std::runtime_error("Cannot open config file " + config_path);
      }
      jf >> jroot;
    }

    const auto& run = cfg.run();
    const auto& det = cfg.detector();

    std::string outdir = run.outdir.empty()
                           ? std::string("outputs/xcheck_ne_pcd")
                           : (run.outdir + "/xcheck_ne_pcd");
    std::filesystem::create_directories(outdir);

    // -------------------------------------------------------------------------
    // Experiment setup
    // -------------------------------------------------------------------------
    ExperimentSetup setup(cfg.experiment_cfg(), det.mass_kg(), run.rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[xcheck] Exposure = " << summary.exposure_kg_year
              << " kg·year, n_e range [" << ne_min << "," << ne_max << "]\n";

    // -------------------------------------------------------------------------
    // Detector / charge-transport configs
    // -------------------------------------------------------------------------
    auto ion = std::make_shared<ChargeIonization>("data/p100K_table.csv");

    // Dark current per exposure (as in your main app)
    const double lambda_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    const double year_s = 365.25 * 86400.0;
    const double exp_time_s = cfg.timing().exposure_time_s;
    const double lambda_per_exp =
        (exp_time_s > 0.0 && year_s > 0.0)
            ? lambda_per_year * (exp_time_s / year_s)
            : 0.0;

    const auto& emj = cfg.response().emc;

    ChargeTransportConfig ct_cfg;
    ct_cfg.thickness_um   = det.geometry().thickness_mm * 1000.0;
    ct_cfg.A_um2          = emj.A_um2;
    ct_cfg.b_umInv        = emj.b_umInv;
    ct_cfg.alpha          = emj.alpha;
    ct_cfg.beta_per_keV   = emj.beta_per_keV;
    ct_cfg.rng_seed       = emj.rng_seed;

    auto ct = std::make_shared<ChargeTransport>(ct_cfg);

    // -------------------------------------------------------------------------
    // EfficiencyMC configuration (row segment)
    // -------------------------------------------------------------------------
    EfficiencyMCConfig emc_cfg;
    emc_cfg.ne_trials = emj.n_events_per_ne;

    emc_cfg.row_length = emj.row_length;
    emc_cfg.pix_cfg.mode           = PixelSimMode::RowSegment;
    emc_cfg.pix_cfg.nx             = emc_cfg.row_length;
    emc_cfg.pix_cfg.ny             = 1;
    emc_cfg.pix_cfg.pixel_size_um  = det.geometry().pixel_size_um;
    emc_cfg.pix_cfg.lambda_dc      = lambda_per_exp;
    emc_cfg.pix_cfg.sigma_readout_e= emj.sigma_readout_e;
    emc_cfg.pix_cfg.rng_seed       = emj.rng_seed;
    emc_cfg.seed                   = emj.rng_seed;

    // Accepted pattern labels derived from experiment.pattern_roi
    emc_cfg.accepted_labels.clear();
    for (int code : summary.pattern_roi) {
      PatternLabel lab;
      lab.isolated = true;
      lab.q = ccdarksens::DecodePatternCode(code);
      if (!lab.q.empty()) emc_cfg.accepted_labels.push_back(lab);
    }
    if (emc_cfg.accepted_labels.empty()) {
      PatternLabel lab;
      lab.isolated = true;
      lab.q = {1};
      emc_cfg.accepted_labels.push_back(lab);
    }

    std::cout << "[xcheck] Accepted pattern labels (from pattern_roi):\n";
    for (std::size_t i = 0; i < emc_cfg.accepted_labels.size(); ++i) {
      std::cout << "  " << i << ": q = { ";
      for (int qv : emc_cfg.accepted_labels[i].q) {
        std::cout << qv << " ";
      }
      std::cout << "}  isolated = "
                << (emc_cfg.accepted_labels[i].isolated ? "true" : "false")
                << "\n";
    }

    // PatternClassifier config
    PatternClassifierConfig pcc;
    const auto& pcc_temp = cfg.response().pattern_classifier;
    pcc.Qmin_e    = pcc_temp.Qmin_e;
    pcc.Qmax_e    = pcc_temp.Qmax_e;
    pcc.enable_MN = pcc_temp.enable_MN;
    pcc.enable_MNL= pcc_temp.enable_MNL;
    pcc.sigma_res_e   = pcc_temp.sigma_res_e;
    pcc.max_e_per_pixel = pcc_temp.max_e_per_pixel;
    pcc.thr_M     = pcc_temp.thr_M;
    pcc.thr_MN    = pcc_temp.thr_MN;
    pcc.thr_MNL   = pcc_temp.thr_MNL;

    auto classifier = std::make_shared<PatternClassifier>(pcc);
    EfficiencyMC emc(emc_cfg, ct, classifier);

    const double Ee_ref_eV = 50.0;
    auto eps_mc = emc.PrecomputeEpsilon(ne_min, ne_max, Ee_ref_eV);

    auto pe = std::make_shared<PatternEfficiency>();
    pe->SetEfficiencyHist(*eps_mc);

    std::cout << "[xcheck] EfficiencyMC epsilon(n_e) at Ee_ref = "
              << Ee_ref_eV << " eV\n";
    for (int ne = ne_min; ne <= ne_max; ++ne) {
      int bin = eps_mc->FindBin(ne);
      double eps = eps_mc->GetBinContent(bin);
      std::cout << "  n_e = " << std::setw(2) << ne
                << "  eps = " << std::fixed << std::setprecision(3) << eps
                << "\n";
    }

    // -------------------------------------------------------------------------
    // PCD response configuration
    // -------------------------------------------------------------------------
    PixelSimulatorConfig pix_cfg_pcd;
    pix_cfg_pcd.mode           = PixelSimMode::LocalPatch;
    pix_cfg_pcd.nx             = 5;
    pix_cfg_pcd.ny             = 5;
    pix_cfg_pcd.pixel_size_um  = det.geometry().pixel_size_um;
    pix_cfg_pcd.lambda_dc      = lambda_per_exp;
    pix_cfg_pcd.sigma_readout_e= emj.sigma_readout_e;
    pix_cfg_pcd.rng_seed       = emj.rng_seed;

    PCDResponseConfig pcd_cfg;
    pcd_cfg.q_min     = cfg.response().pcd.q_min;
    pcd_cfg.q_max     = cfg.response().pcd.q_max;
    pcd_cfg.nbins     = cfg.response().pcd.nbins;
    pcd_cfg.mc_trials = cfg.response().pcd.mc_trials;
    pcd_cfg.pix_cfg   = pix_cfg_pcd;

    auto pcdresp = std::make_shared<PCDBasedResponse>(pcd_cfg, ct);
    auto pcdcalc = std::make_shared<PCDCalculator>();

    // -------------------------------------------------------------------------
    // DetectorResponsePipeline wiring
    // -------------------------------------------------------------------------
    DetectorResponsePipeline pipe(ion);

    auto diff = std::make_shared<Diffusion>(
        ct_cfg.A_um2, ct_cfg.b_umInv, ct_cfg.alpha, ct_cfg.beta_per_keV,
        ct_cfg.thickness_um, 0.08, emj.sigma_readout_e);
    pipe.SetDiffusion(diff);
    pipe.SetPatternEfficiency(pe);
    pipe.SetPCDBasedResponse(pcdresp);
    pipe.SetPCDCalculator(pcdcalc);
    pipe.SetPCDKernelParams(cfg.response().pcd.sigma_res_e, cfg.response().pcd.Dqmin, cfg.response().pcd.Dqmax);

    // -------------------------------------------------------------------------
    // Choose a single DM point from the grid for cross-check
    // -------------------------------------------------------------------------
    const auto& mj = cfg.model();

    DMElectronConfig mc_base;
    mc_base.material          = mj.material;
    mc_base.mediator          = mj.mediator;
    mc_base.rates_dir         = mj.rates_dir;
    mc_base.filename_template = mj.filename_template;
    mc_base.Emin_eV           = mj.Emin_eV;
    mc_base.Emax_eV           = mj.Emax_eV;
    mc_base.nbins             = mj.nbins;

    const auto& jgrid = jroot["model"]["grid"];
    double mchi = first_grid_value(jgrid.at("mchi_MeV"));
    double sigma_val = 0.0;
    {
      const auto& s_spec = jgrid.at("sigma_e_cm2");
      if (s_spec.contains("logspace")) {
        double start_exp = s_spec["logspace"]["start_exp"].get<double>();
        sigma_val = std::pow(10.0, start_exp);
      } else {
        sigma_val = first_grid_value(s_spec);
      }
    }

    std::string fmt_sigma = ".1e";
    if (jgrid.contains("format")) {
      fmt_sigma = jgrid["format"].value("sigma", std::string(".1e"));
    }

    auto format_sigma = [&](double sigma, const std::string& fmt) {
      int prec = 6;
      if (!fmt.empty() && fmt.front() == '.' &&
          (fmt.back() == 'e' || fmt.back() == 'E')) {
        try {
          prec = std::stoi(fmt.substr(1, fmt.size() - 2));
        } catch (...) {
          prec = 6;
        }
      }
      std::ostringstream ss;
      ss.setf(std::ios::scientific);
      ss << std::setprecision(prec) << sigma;
      return ss.str();
    };

    DMElectronConfig mc = mc_base;
    mc.mchi_MeV    = mchi;
    mc.sigma_e_cm2 = format_sigma(sigma_val, fmt_sigma);

    DMElectronModel dm_sig;
    if (!dm_sig.Configure(mc)) {
      throw std::runtime_error("Failed to configure DM model for xcheck.");
    }

    auto dRdE_sig = dm_sig.MakeSpectrum_E();

    double dRdE_int = dRdE_sig->Integral("width");
    double N_raw = dRdE_int * summary.exposure_kg_year;

    std::cout << "[xcheck] Using DM point mchi=" << mchi
              << " MeV, sigma=" << sigma_val
              << " cm^2\n";
    std::cout << "  ∫ dR/dE dE = " << dRdE_int
              << " events/(kg·year)\n";
    std::cout << "  N_raw = " << N_raw << " events\n";

    // -------------------------------------------------------------------------
    // Build S_true(n_e) directly from ChargeIonization
    // -------------------------------------------------------------------------
    auto S_true_ne = ion->FoldToNe(*dRdE_sig, summary.exposure_kg_year,
                                   ne_min, ne_max);

    // -------------------------------------------------------------------------
    // Build P(q | n_true) via PCDBasedResponse and P(n_obs | n_true)
    // -------------------------------------------------------------------------
    const auto& pcd_table =
      pcdresp->BuildPCDTable(ne_min, ne_max, Ee_ref_eV);

    double sigma_res = pcc.sigma_res_e;   // or pix_cfg_pcd.sigma_readout_e
    double Dqmin = 4.4;
    double Dqmax = 3.125;

    auto kernel_ne =
    pcdcalc->BuildNeKernelFromPCD(pcd_table, ne_min, ne_max,
                                    sigma_res, Dqmin, Dqmax);

    std::cout << "[xcheck] P(n_obs | n_true) rows for n_true=1..5:\n";
    for (int ne_true = std::max(ne_min, 1); ne_true <= std::min(ne_max, 5); ++ne_true) {
      int i_true = ne_true - ne_min;
      if (i_true < 0 || i_true >= (int)kernel_ne.size()) continue;
      const auto& row = kernel_ne[(std::size_t)i_true];

      double row_sum = 0.0;
      for (double v : row) row_sum += v;

      std::cout << "  n_true = " << ne_true << "  row_sum = "
                << std::fixed << std::setprecision(3) << row_sum << "  [";
      for (int ne_obs = ne_min; ne_obs <= ne_max; ++ne_obs) {
        int j_obs = ne_obs - ne_min;
        if (j_obs < 0 || j_obs >= (int)row.size()) continue;
        double P = row[(std::size_t)j_obs];
        if (P > 1e-3) {
          std::cout << " (" << ne_obs << " -> " << std::setprecision(3) << P << ")";
        }
      }
      std::cout << " ]\n";
    }

    // -------------------------------------------------------------------------
    // Fold S_true → S_rec(n_e) via the kernel
    // -------------------------------------------------------------------------
    auto S_rec_ne = pcdcalc->FoldNeSpectrum(*S_true_ne, kernel_ne,
                                            ne_min, ne_max,
                                            "S_rec_ne");

    // -------------------------------------------------------------------------
    // Full pattern-space signal via pipeline (uses same machinery)
    // -------------------------------------------------------------------------
    pipe.SetAnalysisSpace(AnalysisSpace::Pattern);
    auto S_obs_ne = pipe.Apply(*dRdE_sig, summary.exposure_kg_year,
                               ne_min, ne_max, Ee_ref_eV);

    // -------------------------------------------------------------------------
    // Pure PCD-space signal
    // -------------------------------------------------------------------------
    pipe.SetAnalysisSpace(AnalysisSpace::PCD);
    auto S_obs_q = pipe.Apply(*dRdE_sig, summary.exposure_kg_year,
                              ne_min, ne_max, Ee_ref_eV);

    // -------------------------------------------------------------------------
    // Write everything to ROOT file
    // -------------------------------------------------------------------------
    std::string out_path = outdir + "/xcheck_ne_pcd.root";
    TFile fout(out_path.c_str(), "RECREATE");
    if (!fout.IsOpen()) {
      std::cerr << "[xcheck] ERROR: cannot create output file " << out_path << "\n";
      return 1;
    }

    S_true_ne->SetName("S_true_ne");
    S_true_ne->Write();
    S_rec_ne->Write();          // "S_rec_ne"
    S_obs_ne->SetName("S_obs_ne");
    S_obs_ne->Write();
    S_obs_q->SetName("S_obs_q");
    S_obs_q->Write();

    eps_mc->SetName("pattern_efficiency");
    eps_mc->Write();

    // Also store P(q|n_true) as a TH2 if you like (optional)
    TH2D h_pq("Pq_ne_true", ";n_{e}^{true};q [e^{-}];P(q|n_{e}^{true})",
              ne_max - ne_min + 1, ne_min - 0.5, ne_max + 0.5,
              cfg.response().pcd.nbins,
              cfg.response().pcd.q_min, cfg.response().pcd.q_max);

    for (const auto& kv : pcd_table) {
      int ne_true = kv.first;
      const TH1D* hq = kv.second.get();
      if (!hq) continue;
      int ix = h_pq.GetXaxis()->FindBin(ne_true);
      for (int ib = 1; ib <= hq->GetNbinsX(); ++ib) {
        double qcenter = hq->GetXaxis()->GetBinCenter(ib);
        double Pq      = hq->GetBinContent(ib);
        int iy = h_pq.GetYaxis()->FindBin(qcenter);
        h_pq.SetBinContent(ix, iy, Pq);
      }
    }
    h_pq.Write();

    fout.Close();
    std::cout << "[xcheck] Done. Output written to " << out_path << "\n";
  }
  catch (const std::exception& ex) {
    std::cerr << "[xcheck] ERROR: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
