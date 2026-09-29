// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ccdarksens_scan_dmelectron_pattern_Edep.cc -- Pattern-mode DM-electron
//  grid scan variant that uses energy-dependent detector folding via
//  DetectorResponsePipeline::ApplyEDependent.
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
//   Pattern-space DM-electron scan variant that folds the signal with the energy-dependent efficiency
//   (DetectorResponsePipeline::ApplyEDependent). Usage: <program> config.json. Steps: parse the
//   config; set up the experiment and the detector-response components; build the dark-current and
//   flat backgrounds in n_e; scan the (m_chi, sigma_e) grid computing the Asimov q; write the results.
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

    // -------------------------------------------------------------------------
    // Experiment setup
    // -------------------------------------------------------------------------
    ExperimentSetup setup(cfg.experiment_cfg(), det.mass_kg(), run.rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[scan-pattern] Exposure = " << summary.exposure_kg_year
              << " kg·year, n_e range [" << ne_min << "," << ne_max << "]\n";

    // -------------------------------------------------------------------------
    // Detector response components
    // -------------------------------------------------------------------------
    auto ion = std::make_shared<ChargeIonization>("data/p100K_table.csv");

    // Dark current per exposure
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

    // ---------------- EfficiencyMC (row-based) ----------------
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

    std::cout << "[scan-pattern] Accepted pattern labels (from pattern_roi):\n";
    for (std::size_t i = 0; i < emc_cfg.accepted_labels.size(); ++i) {
      std::cout << "  " << i << ": q = { ";
      for (int qv : emc_cfg.accepted_labels[i].q) std::cout << qv << " ";
      std::cout << "}  isolated = "
                << (emc_cfg.accepted_labels[i].isolated ? "true" : "false")
                << "\n";
    }

    // PatternClassifier
    PatternClassifierConfig pcc;
    const auto& pcc_temp = cfg.response().pattern_classifier;
    pcc.Qmin_e        = pcc_temp.Qmin_e;
    pcc.Qmax_e        = pcc_temp.Qmax_e;
    pcc.enable_MN     = pcc_temp.enable_MN;
    pcc.enable_MNL    = pcc_temp.enable_MNL;
    pcc.sigma_res_e   = pcc_temp.sigma_res_e;
    pcc.max_e_per_pixel = pcc_temp.max_e_per_pixel;
    pcc.thr_M         = pcc_temp.thr_M;
    pcc.thr_MN        = pcc_temp.thr_MN;
    pcc.thr_MNL       = pcc_temp.thr_MNL;

    auto classifier = std::make_shared<PatternClassifier>(pcc);

    // Build EfficiencyMC as a shared_ptr so we can give it to the pipeline
    auto emc_ptr = std::make_shared<EfficiencyMC>(emc_cfg, ct, classifier);

    //--------------------------------------------------------------
    // Energy grid for energy-dependent pattern efficiency
    //--------------------------------------------------------------
    std::vector<double> E_grid_eV = {2,4,5,7,10,15,20,30,40,50,80,100};
    emc_ptr->PrecomputeEpsilonVsEnergy(E_grid_eV, ne_min, ne_max);

    // Optional debug print
    std::cout << "[scan-pattern] EfficiencyMC epsilon(n_e, E) table:\n";
    for (size_t iE = 0; iE < E_grid_eV.size(); ++iE) {
        std::cout << "  E = " << E_grid_eV[iE] << " eV: ";
        const auto& row = emc_ptr->epsilon_Ene()[iE];
        for (int ne = ne_min; ne <= ne_max; ++ne) {
            double eps = row[ne - ne_min];
            std::cout << " ne=" << ne << "->" << std::fixed << std::setprecision(3)
                      << eps << " ";
        }
        std::cout << "\n";
    }

    // ---------------- PCD response (for detector smearing) -------------------
    PixelSimulatorConfig pix_cfg_pcd;
    pix_cfg_pcd.mode            = PixelSimMode::LocalPatch;
    pix_cfg_pcd.nx              = 5;
    pix_cfg_pcd.ny              = 5;
    pix_cfg_pcd.pixel_size_um   = det.geometry().pixel_size_um;
    pix_cfg_pcd.lambda_dc       = lambda_per_exp;
    pix_cfg_pcd.sigma_readout_e = emj.sigma_readout_e;
    pix_cfg_pcd.rng_seed        = emj.rng_seed;

    PCDResponseConfig pcd_cfg;
    pcd_cfg.q_min     = cfg.response().pcd.q_min;
    pcd_cfg.q_max     = cfg.response().pcd.q_max;
    pcd_cfg.nbins     = cfg.response().pcd.nbins;
    pcd_cfg.mc_trials = cfg.response().pcd.mc_trials;
    pcd_cfg.pix_cfg   = pix_cfg_pcd;

    auto pcdresp = std::make_shared<PCDBasedResponse>(pcd_cfg, ct);
    auto pcdcalc = std::make_shared<PCDCalculator>();

    // ---------------- DetectorResponsePipeline wiring (Pattern mode) ---------
    DetectorResponsePipeline pipe(ion);

    auto diff = std::make_shared<Diffusion>(
        ct_cfg.A_um2, ct_cfg.b_umInv, ct_cfg.alpha, ct_cfg.beta_per_keV,
        ct_cfg.thickness_um, 0.08, emj.sigma_readout_e);
    pipe.SetDiffusion(diff);
    pipe.SetPCDBasedResponse(pcdresp);
    pipe.SetPCDCalculator(pcdcalc);
    pipe.SetPCDKernelParams(cfg.response().pcd.sigma_res_e, cfg.response().pcd.Dqmin, cfg.response().pcd.Dqmax);
    pipe.SetAnalysisSpace(AnalysisSpace::Pattern);
    pipe.EnableEDependentPattern(true);

    // give the pipeline access to EfficiencyMC (epsilon(E, n_e))
    pipe.SetEfficiencyMC(emc_ptr);

    // -------------------------------------------------------------------------
    // Backgrounds: build DC + flat in n_e (pattern space)
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

    // // DC background in pattern space (builder can optionally use PatternEfficiency)
    // auto pe_bkg = std::make_shared<PatternEfficiency>();
    // pe_bkg->SetEfficiencyHist(*eps_mc);
    // bld.SetPatternEfficiency(pe_bkg);

    // auto B_dc_ne = bld.BuildBkgAsimov();  // TH1D in n_e
    auto B_dc_ne = bld.BuildBkgAsimov_EDependent(
    emc_ptr->energy_grid_eV(),
    emc_ptr->epsilon_Ene()
    );  // TH1D in n_e

    // Flat background: build via pipeline in pattern space
    const auto& mj = cfg.model();
    DMElectronConfig mc_base;
    mc_base.material          = mj.material;
    mc_base.mediator          = mj.mediator;
    mc_base.rates_dir         = mj.rates_dir;
    mc_base.filename_template = mj.filename_template;
    mc_base.Emin_eV           = mj.Emin_eV;
    mc_base.Emax_eV           = mj.Emax_eV;
    mc_base.nbins             = mj.nbins;

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

    pipe.SetAnalysisSpace(AnalysisSpace::Pattern);
    // auto B_flat_ne = pipe.Apply(*dRdE_flat, summary.exposure_kg_year,
    //                             ne_min, ne_max, Ee_ref_eV);

    auto B_flat_ne = pipe.ApplyEDependent(
      *dRdE_flat,
      summary.exposure_kg_year,
      ne_min,
      ne_max,
      emc_ptr->energy_grid_eV(),
      emc_ptr->epsilon_Ene()
    );

    // Total background in n_e
    auto B_tot = std::make_unique<TH1D>(*B_dc_ne);
    if (B_flat_ne) B_tot->Add(B_flat_ne.get());

    std::cout << "\n[scan-pattern] Background summary (pattern space n_e)\n";
    double B_int = B_tot->Integral();
    std::cout << "  Total B_tot integral = " << B_int << " events\n";
    double B_roi = 0.0;
    for (int ne : summary.roi_bins) {
      int bin = B_tot->FindBin(ne);
      double b = B_tot->GetBinContent(bin);
      B_roi += b;
      std::cout << "    n_e = " << std::setw(2) << ne
                << "  B(n_e) = " << b << "\n";
    }
    std::cout << "  ROI sum B = " << B_roi << "\n\n";

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

    TH2D h_q("q_mchi_sigma_pattern",
             ";m_{#chi} [MeV];#sigma_{e} [cm^{2}];q_{Asimov} (pattern)",
             static_cast<int>(mchi_list.size()), edges_m.data(),
             static_cast<int>(sigma_list.size()), edges_s.data());

    // -------------------------------------------------------------------------
    // Scan grid and compute Asimov q in pattern space
    // -------------------------------------------------------------------------
    DMElectronConfig mc_sig_base = mc_base;
    std::cout << std::scientific << std::setprecision(8);

    int idx_mchi = 0;
    for (double mchi : mchi_list) {
      int idx_sigma = 0;
      for (double sigma_val : sigma_list) {

        DMElectronConfig mc = mc_sig_base;
        mc.mchi_MeV    = mchi;
        mc.sigma_e_cm2 = format_sigma(sigma_val, fmt_sigma);

        DMElectronModel dm_sig;
        if (!dm_sig.Configure(mc)) {
          std::cerr << "[scan-pattern] WARNING: failed to configure DM model for "
                    << "mchi=" << mchi
                    << ", sigma=" << sigma_val << "\n";
          h_q.SetBinContent(idx_mchi + 1, idx_sigma + 1, 0.0);
          ++idx_sigma;
          continue;
        }
        // std::cout << "[scan-pattern] Computing DM signal for mchi=" << mchi
        //           << " MeV, sigma=" << sigma_val << " cm^2\n";

        auto dRdE_sig = dm_sig.MakeSpectrum_E();

        pipe.SetAnalysisSpace(AnalysisSpace::Pattern);
        auto S_obs = pipe.ApplyEDependent(
            *dRdE_sig,
            summary.exposure_kg_year,
            ne_min,
            ne_max,
            emc_ptr->energy_grid_eV(),
            emc_ptr->epsilon_Ene()
        );

        if (idx_mchi == 0 && idx_sigma == 0) {
          double dRdE_int = dRdE_sig->Integral("width"); // events/(kg·year)
          double N_raw = dRdE_int * summary.exposure_kg_year;
          double S_int = S_obs->Integral();
          std::cout << "[scan-pattern] First grid point mchi=" << mchi
                    << " MeV, sigma=" << sigma_val << " cm^2\n";
          std::cout << "  ∫ dR/dE dE = " << dRdE_int
                    << " events/(kg·year)\n";
          std::cout << "  N_raw = " << N_raw << " events\n";
          std::cout << "  ∫ S_obs(n_e) = " << S_int << " events\n";
        }

        // After computing S_obs and before q_ts calculation:
        if (mchi == 0.53 /*or use the actual variable */) {
            std::cout << "[debug] mchi=" << mchi << " sigma=" << sigma_val << "\n";
            for (int ne = ne_min; ne <= ne_max; ++ne) {
                double sval = S_obs->GetBinContent(S_obs->FindBin(ne));
                if (ne >= 2 && ne <= 5) {
                    std::cout << "  S_obs(ne=" << ne << ") = " << sval << "\n";
                }
            }
        }

        // Asimov q-statistic in n_e pattern space (ROI only)
        double q_ts = 0.0;
        for (int ne : summary.roi_bins) {
          int bin = S_obs->FindBin(ne);
          double s = S_obs->GetBinContent(bin);
          double b = B_tot->GetBinContent(bin);

          if (b <= 0.0)        q_ts += 2.0 * s;
          else if (s > 0.0)    q_ts += 2.0 * (s - b * std::log(1.0 + s / b));
        }

        h_q.SetBinContent(idx_mchi + 1, idx_sigma + 1, q_ts);

        if (run.verbosity >= 2) {
          std::cout << "[scan-pattern] mchi=" << mchi
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
    std::string out_path = run.outdir + "/scan_dmelectron_pattern.root";
    TFile fout(out_path.c_str(), "RECREATE");
    if (!fout.IsOpen()) {
      std::cerr << "[scan-pattern] ERROR: cannot create output file " << out_path << "\n";
      return 1;
    }

    B_tot->Write("B_tot_ne");
    // eps_mc->Write("pattern_efficiency");
    h_q.Write("q_mchi_sigma_pattern");
    fout.Close();

    std::cout << "[scan-pattern] Done. Output written to " << out_path << "\n";
  }
  catch (const std::exception& ex) {
    std::cerr << "[scan-pattern] ERROR: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
