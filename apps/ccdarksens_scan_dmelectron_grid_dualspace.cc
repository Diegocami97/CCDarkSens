// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ccdarksens_scan_dmelectron_grid_dualspace.cc -- Dual-space DM-electron
//  grid scan that computes test statistics in both reconstructed n_e and PCD
//  q observables on the same (mχ, σe) grid.
// ===========================================================================

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

#include "ccdarksens/backgrounds/BackgroundBuilder.hh"

using nlohmann::json;
using namespace ccdarksens;

#include "ccdarksens/utils/AppUtils.hh"
static auto expand_axis = [](const json& spec, const std::string& kind) {
  return ccdarksens::utils::ExpandAxis(spec, kind);
};

// ----------------------------------------------------------------------------
// format_sigma
//   Coupling value -> the string used in the rate-file names; a format like ".3e" gives 3 digits in scientific notation, anything else falls back to 6.
// ----------------------------------------------------------------------------
static std::string format_sigma(double sigma, const std::string& fmt)
{
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
}

// ----------------------------------------------------------------------------
// make_edges_from_centers
//   Histogram bin edges for bin centres c: midpoints between neighbors, outer edges extended by half a step (a single centre gets +/-50%).
// ----------------------------------------------------------------------------
static std::vector<double> make_edges_from_centers(const std::vector<double>& c)
{
  const std::size_t N = c.size();
  std::vector<double> edges(N + 1);
  if (N == 0) return edges;

  if (N == 1) {
    const double w = std::abs(c[0]) > 0 ? std::abs(c[0]) * 0.5 : 0.5;
    edges[0] = c[0] - w;
    edges[1] = c[0] + w;
    return edges;
  }

  edges[0] = c[0] - 0.5 * (c[1] - c[0]);
  for (std::size_t i = 1; i < N; ++i) {
    edges[i] = 0.5 * (c[i - 1] + c[i]);
  }
  edges[N] = c[N - 1] + 0.5 * (c[N - 1] - c[N - 2]);
  return edges;
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// main
//   Dual-space DM-electron scan: on the same (m_chi, sigma_e) grid I compute the Asimov q in
//   pattern space or in PCD space (reconstructed charge). Usage: <program> config.json. Steps:
//   parse the config; set up the experiment and the detector-response components
//   (ChargeTransport, EfficiencyMC, PCD response, DetectorResponsePipeline); build the dark-current
//   and flat backgrounds in the chosen space; scan the grid; write the results.
// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json\n";
    return 1;
  }

  const std::string config_path = argv[1];

  try {
    // -------------------------------------------------------------------------
    // Parse configuration
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

    if (!run.outdir.empty()) {
      std::filesystem::create_directories(run.outdir);
    }

    // Decide analysis space
    std::string analysis_space_str = "pattern";
    if (jroot.contains("response") &&
        jroot["response"].contains("analysis_space")) {
      analysis_space_str = jroot["response"]["analysis_space"].get<std::string>();
    } else {
      analysis_space_str = cfg.response().analysis_space;
    }
    const bool use_pcd_space = (analysis_space_str == "pcd");

    std::cout << "[scan] analysis_space = " << analysis_space_str
              << " (use_pcd_space=" << (use_pcd_space ? "true" : "false") << ")\n";

    // -------------------------------------------------------------------------
    // Experiment setup
    // -------------------------------------------------------------------------
    ExperimentSetup setup(cfg.experiment_cfg(), det.mass_kg(), run.rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[scan] Exposure = " << summary.exposure_kg_year
              << " kg·year, n_e range [" << ne_min << "," << ne_max << "]\n";

    // -------------------------------------------------------------------------
    // Detector response components
    // -------------------------------------------------------------------------
    auto ion = std::make_shared<ChargeIonization>("data/p100K_table.csv");

    // Background dark current rate (for both EfficiencyMC & PCD)
    const double lambda_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    const double year_s = 365.25 * 86400.0;
    const double exp_time_s = cfg.timing().exposure_time_s;
    const double lambda_per_exp =
        (exp_time_s > 0.0 && year_s > 0.0)
            ? lambda_per_year * (exp_time_s / year_s)
            : 0.0;

    // ---------------- ChargeTransport (physics of diffusion, depth, etc.) ----
    const auto& emj = cfg.response().emc;

    ChargeTransportConfig ct_cfg;
    ct_cfg.thickness_um   = det.geometry().thickness_mm * 1000.0;
    ct_cfg.A_um2          = emj.A_um2;
    ct_cfg.b_umInv        = emj.b_umInv;
    ct_cfg.alpha          = emj.alpha;
    ct_cfg.beta_per_keV   = emj.beta_per_keV;
    ct_cfg.rng_seed       = emj.rng_seed;

    auto ct = std::make_shared<ChargeTransport>(ct_cfg);

   // ---------------- EfficiencyMC (uses ct + PixelSimulator + PatternClassifier)
    EfficiencyMCConfig emc_cfg;
    emc_cfg.ne_trials = emj.n_events_per_ne;

    emc_cfg.row_length = emj.row_length;

    emc_cfg.pix_cfg.mode           = PixelSimMode::RowSegment;  // <--- row-based mode
    emc_cfg.pix_cfg.nx             = emc_cfg.row_length;
    emc_cfg.pix_cfg.ny             = 1;
    emc_cfg.pix_cfg.pixel_size_um  = det.geometry().pixel_size_um;
    emc_cfg.pix_cfg.lambda_dc      = lambda_per_exp;
    emc_cfg.pix_cfg.sigma_readout_e= emj.sigma_readout_e;
    emc_cfg.pix_cfg.rng_seed       = emj.rng_seed;

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

    for (std::size_t i = 0; i < emc_cfg.accepted_labels.size(); ++i) {
      std::cout << "[dualspace] Accepted pattern " << i << ": q = { ";
      for (int qv : emc_cfg.accepted_labels[i].q) std::cout << qv << " ";
      std::cout << "}  isolated = " << (emc_cfg.accepted_labels[i].isolated ? "true" : "false") << "\n";
    }

    emc_cfg.seed = emj.rng_seed;

    // PatternClassifier config from JSON
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
    // plus any new thresholds / sigma_res_e I added

    auto classifier = std::make_shared<PatternClassifier>(pcc);
    EfficiencyMC emc(emc_cfg, ct, classifier);


    const double Ee_ref_eV = 50.0;
    auto eps_mc = emc.PrecomputeEpsilon(ne_min, ne_max, Ee_ref_eV);

    auto pe = std::make_shared<PatternEfficiency>();
    pe->SetEfficiencyHist(*eps_mc);

    std::cout << "[XCheck] EfficiencyMC epsilon(n_e) at Ee_ref=" << Ee_ref_eV << " eV\n";
    for (int ne = ne_min; ne <= ne_max; ++ne) {
      int bin = eps_mc->FindBin(ne);
      double eps = eps_mc->GetBinContent(bin);
      std::cout << "  n_e = " << std::setw(2) << ne
                << "  eps = " << std::fixed << std::setprecision(3) << eps << "\n";
    }


    // ---------------- PCD response (ChargeTransport + PixelSimulator table) ---
    PixelSimulatorConfig pix_cfg_pcd;
    pix_cfg_pcd.mode           = PixelSimMode::LocalPatch;
    pix_cfg_pcd.nx             = 5;  // slightly larger patch for PCD
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

    // ---------------- DetectorResponsePipeline wiring ------------------------
    DetectorResponsePipeline pipe(ion);
    // Diffusion object is still used for any fast/simple branch if needed
    auto diff = std::make_shared<Diffusion>(
        ct_cfg.A_um2, ct_cfg.b_umInv, ct_cfg.alpha, ct_cfg.beta_per_keV,
        ct_cfg.thickness_um, 0.08, emj.sigma_readout_e);
    pipe.SetDiffusion(diff);

    pipe.SetPatternEfficiency(pe);
    pipe.SetPCDBasedResponse(pcdresp);
    pipe.SetPCDCalculator(pcdcalc);
    pipe.SetPCDKernelParams(cfg.response().pcd.sigma_res_e, cfg.response().pcd.Dqmin, cfg.response().pcd.Dqmax);

    pipe.SetAnalysisSpace(use_pcd_space
                          ? AnalysisSpace::PCD
                          : AnalysisSpace::Pattern);

    // -------------------------------------------------------------------------
    // Build backgrounds: Dark current always in n_e, then optionally folded
    // -------------------------------------------------------------------------
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
    dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    dcc.norm_scale                = cfg.backgrounds().norm_scale;
    bld.SetDarkCurrent(dcc);

    // Pattern mode: apply pattern ε in BackgroundBuilder via PatternEfficiency
    if (!use_pcd_space) {
      auto pe_bkg = std::make_shared<PatternEfficiency>();
      pe_bkg->SetEfficiencyHist(*eps_mc);
      bld.SetPatternEfficiency(pe_bkg);
    }

    // DC background in n_e
    auto B_dc_ne = bld.BuildBkgAsimov();  // TH1D in n_e

    // Convert to map<int,double> for PCD folding
    std::map<int,double> B_dc_map;
    for (int ne = ne_min; ne <= ne_max; ++ne) {
      int bin = B_dc_ne->FindBin(ne);
      B_dc_map[ne] = B_dc_ne->GetBinContent(bin);
    }

    // -------------------------------------------------------------------------
    // Flat background in energy → n_e or q via pipeline
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

    // Use dummy DM config to get dE binning, then set flat rate per bin
    DMElectronConfig mc_dummy = mc_base;
    mc_dummy.mchi_MeV         = 1.0;
    mc_dummy.sigma_e_cm2      = "1e-40";

    DMElectronModel dm_dummy;
    dm_dummy.Configure(mc_dummy);
    auto dRdE_flat = dm_dummy.MakeSpectrum_E();
    dRdE_flat->Reset("ICES");

    double flat_rate_per_eV = 0.0;
    if (cfg.backgrounds().has_flat_bkg) {
      flat_rate_per_eV =
          cfg.backgrounds().flat_bkg_norm_per_kg_year / 1000.0;
    }
    for (int ib = 1; ib <= dRdE_flat->GetNbinsX(); ++ib) {
      dRdE_flat->SetBinContent(ib, flat_rate_per_eV);
    }

    std::unique_ptr<TH1D> B_flat_ne;
    std::unique_ptr<TH1D> B_flat_q;

    if (!use_pcd_space) {
      pipe.SetAnalysisSpace(AnalysisSpace::Pattern);
      B_flat_ne = pipe.Apply(*dRdE_flat, summary.exposure_kg_year,
                             ne_min, ne_max, Ee_ref_eV);
    } else {
      pipe.SetAnalysisSpace(AnalysisSpace::PCD);
      B_flat_q = pipe.Apply(*dRdE_flat, summary.exposure_kg_year,
                            ne_min, ne_max, Ee_ref_eV);
    }

    // -------------------------------------------------------------------------
    // DC background: pattern vs PCD
    // -------------------------------------------------------------------------
    std::unique_ptr<TH1D> B_dc_q;

    if (!use_pcd_space) {
      // In pattern mode, B_dc_ne already includes ε_pattern if PatternEfficiency set
      B_dc_q = std::make_unique<TH1D>(*B_dc_ne); // name doesn't matter here
    } else {
      // In PCD mode, fold B_dc_true(n_e) with P(q | n_e)
      const auto& pcd_table =
          pcdresp->BuildPCDTable(ne_min, ne_max, Ee_ref_eV);
      std::cout << "[XCheck] PCD table normalization P(q|n_e)\n";
        for (int ne = ne_min; ne <= ne_max; ++ne) {
        auto it = pcd_table.find(ne);
        if (it == pcd_table.end()) continue;
        const TH1D& h = *(it->second);
        double sum = h.Integral();
        std::cout << "  n_e = " << std::setw(2) << ne
                    << "  ∑_q P(q|n_e) = " << std::fixed << std::setprecision(3) << sum
                    << "\n";
        }    

      std::map<int,double> S_empty;
      auto folded_dc = pcdcalc->FoldSpectra(pcd_table, S_empty, B_dc_map);
      B_dc_q = std::move(folded_dc.second);
    }

    // -------------------------------------------------------------------------
    // Total background histogram used in likelihood
    // -------------------------------------------------------------------------
    std::unique_ptr<TH1D> B_tot;

    if (!use_pcd_space) {
      B_tot = std::make_unique<TH1D>(*B_dc_ne);
      if (B_flat_ne) B_tot->Add(B_flat_ne.get());
    } else {
      B_tot = std::make_unique<TH1D>(*B_dc_q);
      if (B_flat_q) B_tot->Add(B_flat_q.get());
    }

    std::cout << "\n[XCheck] Background summary: "
          << (use_pcd_space ? "PCD space (q)" : "Pattern space (n_e)") << "\n";

    double B_int = B_tot->Integral();
    std::cout << "  Total B_tot integral = " << B_int << " events\n";

    if (!use_pcd_space) {
    // pattern space: report per-ROI-bin n_e
    double B_roi = 0.0;
    for (int ne : summary.roi_bins) {
        int bin = B_tot->FindBin(ne);
        double b = B_tot->GetBinContent(bin);
        B_roi += b;
        std::cout << "    n_e = " << std::setw(2) << ne
                << "  B(n_e) = " << b << "\n";
    }
    std::cout << "  ROI sum B = " << B_roi << "\n";
    } else {
    // PCD: just show a few q bins around the peak
    const int nbq = B_tot->GetNbinsX();
    std::cout << "  Sample B(q) bins:\n";
    for (int ib = 1; ib <= nbq; ib += std::max(1, nbq/10)) {
        double qcenter = B_tot->GetXaxis()->GetBinCenter(ib);
        double b = B_tot->GetBinContent(ib);
        std::cout << "    q = " << std::setw(6) << std::setprecision(2) << qcenter
                << "  B(q) = " << std::scientific << b << "\n";
    }
    }
    std::cout << "\n";

    // -------------------------------------------------------------------------
    // Build (mchi, sigma) grid
    // -------------------------------------------------------------------------
    const auto& jgrid = jroot["model"]["grid"];
    auto mchi_list   = expand_axis(jgrid.at("mchi_MeV"), "mchi_MeV");
    auto sigma_list  = expand_axis(jgrid.at("sigma_e_cm2"), "sigma_e_cm2");

    if (mchi_list.empty() || sigma_list.empty()) {
      throw std::runtime_error("Empty mchi_MeV or sigma_e_cm2 grid.");
    }

    std::string fmt_sigma = ".1e";
    if (jgrid.contains("format")) {
      fmt_sigma = jgrid["format"].value("sigma", std::string(".1e"));
    }

    auto edges_m = make_edges_from_centers(mchi_list);
    auto edges_s = make_edges_from_centers(sigma_list);

    TH2D h_q("q_mchi_sigma",
             ";m_{#chi} [MeV];#sigma_{e} [cm^{2}];q",
             static_cast<int>(mchi_list.size()), edges_m.data(),
             static_cast<int>(sigma_list.size()), edges_s.data());

    // -------------------------------------------------------------------------
    // Scan grid and compute Asimov q in pattern or PCD space
    // -------------------------------------------------------------------------
    int idx_mchi = 0;
    for (double mchi : mchi_list) {
      int idx_sigma = 0;
      for (double sigma_val : sigma_list) {

        DMElectronConfig mc = mc_base;
        mc.mchi_MeV    = mchi;
        mc.sigma_e_cm2 = format_sigma(sigma_val, fmt_sigma);

        // Configure DM model
        //
        DMElectronModel dm_sig;
        if (!dm_sig.Configure(mc)) { // failed to configure means no rates available
          std::cerr << "[scan] WARNING: failed to configure DM model for "
                    << "mchi=" << mchi
                    << ", sigma=" << sigma_val << "\n";
          h_q.SetBinContent(idx_mchi + 1, idx_sigma + 1, 0.0);
          ++idx_sigma;
          continue;
        }

        std::cout << "[scan] Computing DM signal for mchi=" << mchi
                  << " MeV, sigma=" << sigma_val << " cm^2\n";
        auto dRdE_sig = dm_sig.MakeSpectrum_E();

        pipe.SetAnalysisSpace(use_pcd_space
                              ? AnalysisSpace::PCD
                              : AnalysisSpace::Pattern);

        auto S_obs = pipe.Apply(*dRdE_sig, summary.exposure_kg_year,
                                ne_min, ne_max, Ee_ref_eV);

        if (idx_mchi == 0 && idx_sigma == 0) {
            double dRdE_int = dRdE_sig->Integral("width"); // events/(kg·year)
            double N_raw = dRdE_int * summary.exposure_kg_year;

            double S_int = S_obs->Integral();
            std::cout << "[XCheck] First grid point mchi=" << mchi
                        << " MeV, sigma=" << sigma_val << " cm^2\n";
            std::cout << "  ∫ dR/dE dE = " << dRdE_int
                        << " events/(kg·year)\n";
            std::cout << "  N_raw = " << N_raw << " events\n";
            std::cout << "  ∫ S_obs(" << (use_pcd_space ? "q" : "n_e")
                        << ") = " << S_int << " events\n";
        }

        double q_ts = 0.0;

        if (!use_pcd_space) {
          // Pattern-space likelihood: sum over ROI n_e bins
          for (int ne : summary.roi_bins) {
            int bin = S_obs->FindBin(ne);
            double s = S_obs->GetBinContent(bin);
            double b = B_tot->GetBinContent(bin);

            if (b <= 0.0)        q_ts += 2.0 * s;
            else if (s > 0.0)    q_ts += 2.0 * (s - b * std::log(1.0 + s / b));
          }
        } else {
          // PCD-space likelihood: sum over all q-bins (could restrict to q-ROI)
          const int nbq = S_obs->GetNbinsX();
          for (int ib = 1; ib <= nbq; ++ib) {
            double s = S_obs->GetBinContent(ib);
            double b = B_tot->GetBinContent(ib);

            if (b <= 0.0)        q_ts += 2.0 * s;
            else if (s > 0.0)    q_ts += 2.0 * (s - b * std::log(1.0 + s / b));
          }
        }

        h_q.SetBinContent(idx_mchi + 1, idx_sigma + 1, q_ts);

        if (run.verbosity >= 2) {
          std::cout << "[scan] mchi=" << mchi
                    << " MeV, sigma=" << sigma_val
                    << " cm^2 → q=" << q_ts << "\n";
        }

        ++idx_sigma;
      }
      ++idx_mchi;
    }

    // -------------------------------------------------------------------------
    // Write output
    // -------------------------------------------------------------------------
    std::string out_path = run.outdir + "/scan_dmelectron_grid_dualspace.root";
    TFile fout(out_path.c_str(), "RECREATE");
    if (!fout.IsOpen()) {
      std::cerr << "[scan] ERROR: cannot create output file " << out_path << "\n";
      return 1;
    }

    B_tot->Write(use_pcd_space ? "B_tot_q" : "B_tot_ne");
    eps_mc->Write("pattern_efficiency");
    h_q.Write("q_mchi_sigma");
    fout.Close();

    std::cout << "[scan] Done. Output written to " << out_path << "\n";
  }
  catch (const std::exception& ex) {
    std::cerr << "[scan] ERROR: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
