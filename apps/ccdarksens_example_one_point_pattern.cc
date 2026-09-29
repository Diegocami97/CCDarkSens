// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ccdarksens_example_one_point_pattern.cc -- Single (mχ, σe) diagnostic run
//  of the full pattern pipeline (S_obs→S_pat, backgrounds, PLR) without a
//  full grid scan.
// ===========================================================================

#include <algorithm>
#include <chrono>
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
#include <set>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TParameter.h>
#include <TCanvas.h>
#include <TLegend.h>

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
#include "ccdarksens/response/PatternImageGenerator.hh"
#include "ccdarksens/response/PatternRates.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/response/PCDBasedResponse.hh"
#include "ccdarksens/response/PCDCalculator.hh"
#include "ccdarksens/stats/TestStatisticFactory.hh"
#include "ccdarksens/stats/ProfileLikelihood.hh"
#include "ccdarksens/backgrounds/BackgroundBuilder.hh"

using nlohmann::json;
using namespace ccdarksens;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------
#include "ccdarksens/utils/AppUtils.hh"
static auto expand_axis = [](const json& spec, const std::string& kind) {
  return ccdarksens::utils::ExpandAxis(spec, kind);
};

// ----------------------------------------------------------------------------
// EfficiencyMcJsonBlock
//   Pointer to response.efficiency_mc, or to the old key response.pattern_mc, or nullptr if neither exists.
// ----------------------------------------------------------------------------
static const json* EfficiencyMcJsonBlock(const json& jroot) {
  if (!jroot.contains("response")) return nullptr;
  const auto& jresp = jroot["response"];
  if (jresp.contains("efficiency_mc")) return &jresp["efficiency_mc"];
  if (jresp.contains("pattern_mc")) return &jresp["pattern_mc"];
  return nullptr;
}

// ----------------------------------------------------------------------------
// format_sigma
//   Coupling value -> the string used in the rate-file names; a format like ".3e" gives 3 digits in scientific notation, anything else falls back to 6.
// ----------------------------------------------------------------------------
static std::string format_sigma(double sigma, const std::string& fmt)
{
  int prec = 6;
  if (!fmt.empty() && fmt.size() >= 3 && fmt.front() == '.' &&
      (fmt.back() == 'e' || fmt.back() == 'E')) {
    try { prec = std::stoi(fmt.substr(1, fmt.size() - 2)); } catch (...) {}
  }
  std::ostringstream ss;
  ss.setf(std::ios::scientific);
  ss << std::setprecision(prec) << sigma;
  return ss.str();
}

/// Load observed counts from CSV: one row of comma-separated numbers (or one number per line). Returns empty on error.
static std::vector<double> load_data_csv(const std::string& path, std::size_t expected_size)
{
  std::ifstream in(path);
  if (!in) return {};
  std::vector<double> out;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, ',')) {
      try {
        out.push_back(std::stod(cell));
      } catch (...) { break; }
    }
    if (out.size() >= expected_size) break;
  }
  return out;
}

/// Load observed counts from ROOT (histogram D_pat). Units: counts per pattern, same order as pattern_roi.
static std::vector<double> load_data_root(const std::string& path, std::size_t expected_size)
{
  TFile f(path.c_str(), "READ");
  if (!f.IsOpen()) return {};
  TH1D* h = nullptr;
  f.GetObject("D_pat", h);
  if (!h) return {};
  const int nb = h->GetNbinsX();
  std::vector<double> out;
  out.reserve(static_cast<std::size_t>(nb));
  for (int i = 1; i <= nb; ++i) out.push_back(h->GetBinContent(i));
  return out;
}

// ----------------------------------------------------------------------------
// get_exposure_from_data_file
//   Read the exposure_kg_year parameter stored in a ROOT data file. Returns false (and leaves *out alone) if the file is not ROOT or has no such parameter.
// ----------------------------------------------------------------------------
static bool get_exposure_from_data_file(const std::string& path, double* out)
{
  if (!out || path.size() < 6 || path.compare(path.size() - 5, 5, ".root") != 0) return false;
  TFile f(path.c_str(), "READ");
  if (!f.IsOpen()) return false;
  TParameter<double>* p = nullptr;
  f.GetObject("exposure_kg_year", p);
  if (!p) return false;
  *out = p->GetVal();
  return true;
}

// Load observed counts: a .root path goes to load_data_root, anything else is read as CSV.
static std::vector<double> load_data(const std::string& path, std::size_t expected_size)
{
  if (path.size() >= 5 && path.compare(path.size() - 5, 5, ".root") == 0)
    return load_data_root(path, expected_size);
  return load_data_csv(path, expected_size);
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// main
//   Single-point diagnostic of the full pattern pipeline, without a grid scan.
//   Usage: <program> config.json [mchi_MeV] [sigma_e_cm2] (the first grid point is used if
//   the mass and coupling are omitted). Steps:
//     1) build the detector response (ionization, charge transport, EfficiencyMC, classifier)
//        and the pattern efficiency table;
//     2) build the background model (dark current, flat spectrum, Bp/Br or migration);
//     3) compute the signal spectrum at the requested point and fold it into pattern space,
//        then evaluate the profile likelihood;
//     4) write the spectra and results to a ROOT file.
// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json [mchi_MeV] [sigma_e_cm2]\n"
              << "  If mchi/sigma omitted, first grid point from config is used.\n";
    return 1;
  }

  const std::string config_path = argv[1];
  double mchi_MeV = 0.0;
  double sigma_e_cm2 = 0.0;
  bool use_cli_point = (argc >= 4);
  if (use_cli_point) {
    try {
      mchi_MeV = std::stod(argv[2]);
      sigma_e_cm2 = std::stod(argv[3]);
    } catch (...) {
      std::cerr << "Invalid mchi or sigma_e argument.\n";
      return 1;
    }
  }

  try {
    ConfigManager cfg(config_path);
    cfg.parse();

    json jroot;
    {
      std::ifstream jf(config_path);
      if (!jf) throw std::runtime_error("Cannot open config file " + config_path);
      jf >> jroot;
    }

    const auto& run = cfg.run();
    const auto& det = cfg.detector();
    if (!run.outdir.empty()) std::filesystem::create_directories(run.outdir);

    ExperimentSetup setup(cfg.experiment_cfg(), det.mass_kg(), run.rng_seed);
    auto summary = setup.prepare_summary();
    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;
    const bool use_pattern_bins = (summary.observable_bins == "pattern");
    // For background we include n_e=0 when in pattern space (dark current has n_e=0)
    const int ne_min_bkg = use_pattern_bins ? std::min(ne_min, 0) : ne_min;

    if (!use_cli_point) {
      if (jroot.contains("model") && jroot["model"].contains("example_point")) {
        const auto& ep = jroot["model"]["example_point"];
        mchi_MeV = ep.value("mchi_MeV", 0.5);
        sigma_e_cm2 = ep.value("sigma_e_cm2", 1e-27);
        std::cout << "[example-one] Using example_point from config: mchi=" << mchi_MeV
                  << " MeV, sigma=" << sigma_e_cm2 << " cm^2\n";
      } else {
        const auto& jgrid = jroot["model"]["grid"];
        auto mchi_list = expand_axis(jgrid.at("mchi_MeV"), "mchi_MeV");
        auto sigma_list = expand_axis(jgrid.at("sigma_e_cm2"), "sigma_e_cm2");
        if (mchi_list.empty() || sigma_list.empty())
          throw std::runtime_error("Empty grid and no mchi/sigma given on command line.");
        mchi_MeV = mchi_list[0];
        sigma_e_cm2 = sigma_list[0];
        std::cout << "[example-one] Using first grid point: mchi=" << mchi_MeV
                  << " MeV, sigma=" << sigma_e_cm2 << " cm^2\n";
      }
    }

    std::string fmt_sigma = ".1e";
    if (jroot.contains("model") && jroot["model"].contains("grid") &&
        jroot["model"]["grid"].contains("format")) {
      fmt_sigma = jroot["model"]["grid"]["format"].value("sigma", std::string(".1e"));
    }

    // =========================================================================
    // 1. Build detector response pipeline
    // =========================================================================
    const auto& ion_cfg = cfg.response().charge_ionization;
    auto ion = std::make_shared<ChargeIonization>(ion_cfg.table_csv);
    std::cout << "[example] ChargeIonization: " << ion_cfg.table_csv << "\n";
    const auto& emj = cfg.response().emc;

    ChargeTransportConfig ct_cfg;
    ct_cfg.thickness_um = det.geometry().thickness_mm * 1000.0;
    ct_cfg.A_um2 = emj.A_um2;
    ct_cfg.b_umInv = emj.b_umInv;
    ct_cfg.alpha = emj.alpha;
    ct_cfg.beta_per_keV = emj.beta_per_keV;
    ct_cfg.rng_seed = emj.rng_seed;
    auto ct = std::make_shared<ChargeTransport>(ct_cfg);

    EfficiencyMCConfig emc_cfg;
    emc_cfg.ne_trials = emj.n_events_per_ne;
    if (const json* j_emc = EfficiencyMcJsonBlock(jroot);
        j_emc && j_emc->contains("row_length")) {
      emc_cfg.row_length = (*j_emc)["row_length"].get<int>();
    } else {
      emc_cfg.row_length = emj.row_length;
    }
    emc_cfg.pix_cfg.mode = PixelSimMode::RowSegment;
    emc_cfg.pix_cfg.nx = emc_cfg.row_length;
    emc_cfg.pix_cfg.ny = 1;
    emc_cfg.pix_cfg.pixel_size_um = det.geometry().pixel_size_um;
    emc_cfg.pix_cfg.lambda_dc = 0.0;
    emc_cfg.pix_cfg.sigma_readout_e = emj.sigma_readout_e;
    emc_cfg.pix_cfg.rng_seed = emj.rng_seed;
    emc_cfg.seed = emj.rng_seed;
    // Accepted labels follow the analysis space (see scan app for rationale):
    // pattern mode -> multi-pixel pattern_roi; n_e mode -> single-pixel patterns
    // implied by roi_bins (one pixel with n_e electrons).
    emc_cfg.accepted_labels.clear();
    if (use_pattern_bins) {
      for (int code : summary.pattern_roi) {
        PatternLabel lab;
        lab.isolated = true;
        lab.q = ccdarksens::DecodePatternCode(code);
        if (!lab.q.empty()) emc_cfg.accepted_labels.push_back(lab);
      }
    } else {
      for (int ne : summary.roi_bins) {
        if (ne <= 0) continue;
        PatternLabel lab;
        lab.isolated = true;
        lab.q = {ne};  // single pixel holding n_e electrons (valid for n_e <= 9)
        emc_cfg.accepted_labels.push_back(lab);
      }
    }
    if (emc_cfg.accepted_labels.empty()) {
      PatternLabel lab;
      lab.isolated = true;
      lab.q = {1};
      emc_cfg.accepted_labels.push_back(lab);
    }

    std::map<std::pair<int, int>, double> pattern_eff_map;
    std::string efficiency_source_desc;  // "CSV: path" or "EfficiencyMC (N trials per n_e)"
    if (const json* j_emc = EfficiencyMcJsonBlock(jroot);
        j_emc && j_emc->contains("efficiency_csv")) {
      const std::string eff_csv_path = (*j_emc)["efficiency_csv"].get<std::string>();
      efficiency_source_desc = "CSV: " + eff_csv_path;
      std::ifstream eff_in(eff_csv_path);
      if (!eff_in.is_open())
        throw std::runtime_error("Failed to open pattern efficiency CSV: " + eff_csv_path);
      std::string line;
      while (std::getline(eff_in, line)) {
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front())))
          line.erase(line.begin());
        if (line.empty() || line[0] == '#' ||
            line.rfind("pattern", 0) == 0) continue;
        std::stringstream ss(line);
        std::string col;
        int pattern_code = 0, ne_val = 0;
        double eff_val = 0.0;
        if (!std::getline(ss, col, ',')) continue;
        try { pattern_code = std::stoi(col); } catch (...) { continue; }
        if (!std::getline(ss, col, ',')) continue;
        try { ne_val = std::stoi(col); } catch (...) { continue; }
        if (!std::getline(ss, col, ',')) continue;
        try { eff_val = std::stod(col); } catch (...) { continue; }
        pattern_eff_map[{pattern_code, ne_val}] = eff_val;
      }
      std::cout << "[example-one] Loaded " << pattern_eff_map.size()
                << " (pattern, ne) entries from CSV.\n";
    } else {
      std::cout << "[example-one] No efficiency_csv; will use EfficiencyMC table if pattern bins.\n";
    }

    PatternClassifierConfig pcc;
    const auto& pcc_temp = cfg.response().pattern_classifier;
    pcc.Qmin_e = pcc_temp.Qmin_e;
    pcc.neighbor_Qmax_e = pcc_temp.neighbor_Qmax_e;
    pcc.Qmax_e = pcc_temp.Qmax_e;
    pcc.enable_MN = pcc_temp.enable_MN;
    pcc.enable_MNL = pcc_temp.enable_MNL;
    pcc.sigma_res_e = pcc_temp.sigma_res_e;
    pcc.max_e_per_pixel = pcc_temp.max_e_per_pixel;
    pcc.thr_M = pcc_temp.thr_M;
    pcc.thr_MN = pcc_temp.thr_MN;
    pcc.thr_MNL = pcc_temp.thr_MNL;
    pcc.allow_pattern_zero = pcc_temp.allow_pattern_zero;
    pcc.single_pixel_use_round = pcc_temp.single_pixel_use_round;
    auto classifier = std::make_shared<PatternClassifier>(pcc);
    auto emc_ptr = std::make_shared<EfficiencyMC>(emc_cfg, ct, classifier);

    // Optional: build pattern efficiency from 2D image + isolation (notebook-style)
    if (emj.use_2d_image_efficiency) {
      const auto& pimg = cfg.response().pattern_image;
      const bool from_json = (pimg.nrows_binned > 0 && pimg.ncols > 0 && pimg.row_binning > 0);
      const bool from_detector = !from_json && (det.geometry().rows > 0 && det.geometry().cols > 0 && pimg.row_binning > 0 && pimg.col_binning > 0);
      if (from_json || from_detector) {
        PatternImageConfig img_cfg;
        if (from_json) {
          img_cfg.nrows_binned = pimg.nrows_binned;
          img_cfg.ncols        = pimg.ncols;
        } else {
          img_cfg.raw_rows = det.geometry().rows;
          img_cfg.raw_cols = det.geometry().cols;
        }
        img_cfg.row_binning    = pimg.row_binning;
        img_cfg.col_binning    = pimg.col_binning;
        img_cfg.pixel_size_um  = pimg.pixel_size_um;
        img_cfg.sigma_readout_e= pimg.sigma_readout_e;
        img_cfg.lambda_dc      = pimg.lambda_dc;
        img_cfg.rng_seed       = pimg.rng_seed;
        img_cfg.include_dark_current = (pimg.lambda_dc > 0.0);
        auto img_gen_eff = std::make_shared<PatternImageGenerator>(img_cfg, ct);
        emc_ptr->SetPatternImageGenerator(img_gen_eff);
        std::cout << "[example-one] Pattern efficiency will be computed from 2D image (notebook-style, with isolation).\n";
        // Output example 2D pattern image (same generator used for efficiency)
        const double Ee_img_eV = 50.0;
        const int example_ne_img = 5;
        auto image_2d_eff = img_gen_eff->GenerateImage(example_ne_img, Ee_img_eV);
        TH2D h2img_eff("h2_pattern_image_eff",
                       "2D binned pattern image (n_e=5);col;row",
                       img_gen_eff->Ncols(), 0.5, img_gen_eff->Ncols() + 0.5,
                       img_gen_eff->NrowsBinned(), 0.5, img_gen_eff->NrowsBinned() + 0.5);
        for (int r = 0; r < img_gen_eff->NrowsBinned(); ++r)
          for (int c = 0; c < img_gen_eff->Ncols(); ++c)
            h2img_eff.SetBinContent(c + 1, img_gen_eff->NrowsBinned() - r,
                                    image_2d_eff[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]);
        std::string img_png_eff = run.outdir + "/example_2d_pattern_image.png";
        TCanvas cimg_eff("c_2d_eff", "2D pattern image", 800, 300);
        h2img_eff.Draw("COLZ");
        cimg_eff.SaveAs(img_png_eff.c_str());
        std::cout << "[example-one] 2D pattern image (n_e=" << example_ne_img << ") written to " << img_png_eff << "\n";
      }
    }

    const double Ee_ref_eV = 50.0;
    std::unique_ptr<TH1D> h_eps_ne;
    auto t_eff_start = std::chrono::steady_clock::now();
    if (!pattern_eff_map.empty()) {
      h_eps_ne = emc_ptr->PrecomputeEpsilonWithPatternEff(
          ne_min_bkg, ne_max, Ee_ref_eV, pattern_eff_map);
    } else {
      h_eps_ne = emc_ptr->PrecomputeEpsilon(ne_min_bkg, ne_max, Ee_ref_eV);
    }
    auto t_eff_end = std::chrono::steady_clock::now();
    double t_eff_s = std::chrono::duration<double>(t_eff_end - t_eff_start).count();
    std::cout << "[example-one] Pattern efficiency computation took " << std::fixed << std::setprecision(2) << t_eff_s << " s\n";

    PixelSimulatorConfig pix_cfg_pcd;
    pix_cfg_pcd.mode = PixelSimMode::LocalPatch;
    pix_cfg_pcd.nx = 5;
    pix_cfg_pcd.ny = 5;
    pix_cfg_pcd.pixel_size_um = det.geometry().pixel_size_um;
    pix_cfg_pcd.lambda_dc = 0.0;
    pix_cfg_pcd.sigma_readout_e = emj.sigma_readout_e;
    pix_cfg_pcd.rng_seed = emj.rng_seed;
    PCDResponseConfig pcd_cfg;
    pcd_cfg.q_min = cfg.response().pcd.q_min;
    pcd_cfg.q_max = cfg.response().pcd.q_max;
    pcd_cfg.nbins = cfg.response().pcd.nbins;
    pcd_cfg.mc_trials = cfg.response().pcd.mc_trials;
    pcd_cfg.pix_cfg = pix_cfg_pcd;
    auto pcdresp = std::make_shared<PCDBasedResponse>(pcd_cfg, ct);
    auto pcdcalc = std::make_shared<PCDCalculator>();

    DetectorResponsePipeline pipe(ion);
    // Diffusion only used when PCD is not available; with PCD, kernel has diffusion.
    auto diff = std::make_shared<Diffusion>(
        ct_cfg.A_um2, ct_cfg.b_umInv, ct_cfg.alpha, ct_cfg.beta_per_keV,
        ct_cfg.thickness_um, 0.08, emj.sigma_readout_e);
    pipe.SetDiffusion(diff);
    pipe.SetPCDBasedResponse(pcdresp);
    pipe.SetPCDCalculator(pcdcalc);
    pipe.SetPCDKernelParams(cfg.response().pcd.sigma_res_e, cfg.response().pcd.Dqmin, cfg.response().pcd.Dqmax);
    pipe.SetAnalysisSpace(AnalysisSpace::Pattern);
    auto pe = std::make_shared<PatternEfficiency>();
    pe->SetEfficiencyHist(*h_eps_ne);
    pipe.SetPatternEfficiency(pe);
    pipe.SetSkipPatternEfficiency(use_pattern_bins);
    pipe.EnableEDependentPattern(false);

    // =========================================================================
    // 2. Build background model
    // =========================================================================
    BackgroundBuilder bld(det.geometry().rows, det.geometry().cols,
                          det.geometry().active_fraction, ne_min_bkg, ne_max);
    TimingConfig tcfg;
    tcfg.exposure_time_s = cfg.timing().exposure_time_s;
    if (cfg.timing().n_exposures_override.has_value())
      tcfg.n_exposures_override = cfg.timing().n_exposures_override;
    bld.SetTiming(cfg.experiment_cfg().livetime_days,
                  cfg.experiment_cfg().duty_cycle, tcfg);
    DarkCurrentConfig dcc;
    dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    dcc.norm_scale = cfg.backgrounds().norm_scale;
    bld.SetDarkCurrent(dcc);
    if (!use_pattern_bins) bld.SetPatternEfficiency(pe);

    auto B_dc_ne = bld.BuildBkgAsimov();
    const auto& mj = cfg.model();
    DMElectronConfig mc_base;
    mc_base.material = mj.material;
    mc_base.mediator = mj.mediator;
    mc_base.rates_dir = mj.rates_dir;
    mc_base.filename_template = mj.filename_template;
    mc_base.Emin_eV = mj.Emin_eV;
    mc_base.Emax_eV = mj.Emax_eV;
    mc_base.nbins = mj.nbins;
    DMElectronConfig mc_dummy = mc_base;
    mc_dummy.mchi_MeV = 1.0;
    mc_dummy.sigma_e_cm2 = "1e-40";
    DMElectronModel dm_dummy;
    dm_dummy.Configure(mc_dummy);
    auto dRdE_flat = dm_dummy.MakeSpectrum_E();
    dRdE_flat->Reset("ICES");
    double flat_rate_per_eV = 0.0;
    if (cfg.backgrounds().has_flat_bkg)
      flat_rate_per_eV = cfg.backgrounds().flat_bkg_norm_per_kg_year / 1000.0;
    for (int ib = 1; ib <= dRdE_flat->GetNbinsX(); ++ib)
      dRdE_flat->SetBinContent(ib, flat_rate_per_eV);
    // In pattern mode use true n_e for background: ε(pattern|n_true) already has diffusion.
    std::unique_ptr<TH1D> B_flat_ne;
    if (use_pattern_bins) {
      B_flat_ne = ion->FoldToNe(*dRdE_flat, summary.exposure_kg_year, ne_min_bkg, ne_max);
    } else {
      B_flat_ne = pipe.Apply(*dRdE_flat, summary.exposure_kg_year, ne_min_bkg, ne_max, Ee_ref_eV);
    }
    auto B_tot = std::make_unique<TH1D>(*B_dc_ne);
    if (B_flat_ne) B_tot->Add(B_flat_ne.get());

    std::vector<double> B_pat;
    std::vector<std::unique_ptr<TH1D>> charge_histos;  // optional: charge distribution per n_e
    std::vector<std::string> sanity_csv_lines;          // optional: event,ne,pattern_id,q1,q2,...
    if (use_pattern_bins) {
      // Total background in n_e before pattern folding (pattern-space option only)
      std::cout << "[example-one] Total background B_tot(n_e) before pattern fold (includes n_e=0):\n";
      std::cout << std::scientific;
      for (int ne = ne_min_bkg; ne <= ne_max; ++ne) {
        int bin = B_tot->FindBin(static_cast<double>(ne));
        double b = B_tot->GetBinContent(bin);
        std::cout << "  n_e = " << std::setw(2) << ne << "  B_tot = " << b << "\n";
      }
      std::cout << "  integral = " << B_tot->Integral() << "\n";
      std::cout << std::defaultfloat;

      if (summary.pattern_roi.empty())
        throw std::runtime_error("observable_bins is 'pattern' but pattern_roi is empty.");
      if (pattern_eff_map.empty()) {
        const auto& table = emc_ptr->GetPatternTable();
        for (const auto& ne_entry : table) {
          int ne = ne_entry.first;
          for (const auto& label_prob : ne_entry.second) {
            int code = 0;
            for (int d : label_prob.first.q) code = code * 10 + d;
            pattern_eff_map[{code, ne}] = label_prob.second;
          }
        }
        efficiency_source_desc = "EfficiencyMC (BuildPatternTable, "
            + std::to_string(emc_ptr->config().ne_trials) + " trials per n_e)";
        std::cout << "[example-one] Filled pattern efficiencies from EfficiencyMC ("
                  << pattern_eff_map.size() << " entries).\n";
        // Report which pattern_roi are requested vs which have non-zero efficiency from MC
        std::cout << "[example-one] pattern_roi (all requested):";
        for (int pid : summary.pattern_roi) std::cout << " " << pid;
        std::cout << "\n";
        std::set<int> produced;
        for (const auto& kv : pattern_eff_map) {
          if (kv.second > 0.0) produced.insert(kv.first.first);
        }
        std::cout << "[example-one] Patterns with non-zero P(pattern|n_e) from MC:";
        for (int p : produced) std::cout << " " << p;
        if (produced.empty()) std::cout << " (none)";
        std::cout << "\n";
        std::set<int> requested(summary.pattern_roi.begin(), summary.pattern_roi.end());
        bool missing = false;
        for (int p : requested) {
          if (produced.find(p) == produced.end()) {
            if (!missing) std::cout << "[example-one] Patterns in pattern_roi with 0 efficiency (not produced by MC):";
            std::cout << " " << p;
            missing = true;
          }
        }
        if (missing) std::cout << "\n";
      }
      if (pattern_eff_map.empty())
        throw std::runtime_error("No pattern efficiencies (no CSV and EfficiencyMC table empty).");

      // Where efficiencies are computed (notebook-style printout)
      std::cout << "[example-one] Efficiencies computed from: " << efficiency_source_desc << "\n";

      // Write efficiency-per-pattern table to CSV (pattern_roi × roi_bins to match reference)
      std::string eff_csv_out = run.outdir + "/efficiency_per_pattern.csv";
      std::ofstream eff_out(eff_csv_out);
      if (eff_out.is_open()) {
        eff_out << "pattern,ne,Efficiency\n";
        for (int pid : summary.pattern_roi) {
          for (int ne : summary.roi_bins) {
            auto it = pattern_eff_map.find({pid, ne});
            double eff = (it != pattern_eff_map.end()) ? it->second : 0.0;
            eff_out << pid << "," << ne << "," << std::scientific << eff << "\n";
          }
        }
        eff_out.close();
        std::cout << "[example-one] Efficiency table written to " << eff_csv_out << " (pattern_roi × roi_bins)\n";
      }

      // Verbose: print efficiency table P(pattern|n_e) like the notebook
      if (run.verbosity >= 2) {
        std::cout << "[example-one] Efficiency per pattern P(pattern|n_e) (includes n_e=0):\n  n_e \\ pattern";
        for (int p : summary.pattern_roi) std::cout << "\t" << p;
        std::cout << "\n";
        for (int ne = ne_min_bkg; ne <= ne_max; ++ne) {
          std::cout << "  " << ne;
          for (int p : summary.pattern_roi) {
            auto it = pattern_eff_map.find({p, ne});
            double v = (it != pattern_eff_map.end()) ? it->second : 0.0;
            std::cout << "\t" << std::fixed << std::setprecision(4) << v;
          }
          std::cout << "\n";
        }
      }

      // Optional: E-dependent epsilon CSV (Generate_efficiencies with E grid)
      bool output_epsilon_E_csv = (jroot.contains("debug") && jroot["debug"].contains("output_epsilon_E_csv") && jroot["debug"]["output_epsilon_E_csv"].get<bool>());
      if (output_epsilon_E_csv && emc_ptr) {
        std::vector<double> E_grid = {0.0, 50.0, 100.0};
        if (jroot["debug"].contains("output_epsilon_E_grid") && jroot["debug"]["output_epsilon_E_grid"].is_array()) {
          E_grid.clear();
          for (const auto& v : jroot["debug"]["output_epsilon_E_grid"])
            E_grid.push_back(v.get<double>());
        }
        emc_ptr->PrecomputeEpsilonVsEnergy(E_grid, ne_min, ne_max);
        const auto& eg = emc_ptr->energy_grid_eV();
        const auto& eene = emc_ptr->epsilon_Ene();
        std::string eps_E_path = run.outdir + "/epsilon_vs_E.csv";
        std::ofstream ef(eps_E_path);
        if (ef.is_open()) {
          ef << "E_eV,n_e,epsilon_total\n";
          for (size_t iE = 0; iE < eg.size() && iE < eene.size(); ++iE) {
            for (int ne = ne_min; ne <= ne_max; ++ne) {
              int j = ne - ne_min;
              if (j >= 0 && j < (int)eene[iE].size())
                ef << std::scientific << eg[iE] << "," << ne << "," << eene[iE][static_cast<size_t>(j)] << "\n";
            }
          }
          ef.close();
          std::cout << "[example-one] E-dependent epsilon CSV written to " << eps_E_path << "\n";
        }
      }

      // 2D binned image (notebook-style): generate, plot PNG, scan middle row (image size from detector, binning from pattern_image)
      const auto& pimg = cfg.response().pattern_image;
      const bool from_json = (pimg.nrows_binned > 0 && pimg.ncols > 0 && pimg.row_binning > 0);
      const bool from_detector = !from_json && (det.geometry().rows > 0 && det.geometry().cols > 0 && pimg.row_binning > 0 && pimg.col_binning > 0);
      bool do_2d_image = from_json || from_detector;
      if (do_2d_image) {
        PatternImageConfig img_cfg;
        if (from_json) {
          img_cfg.nrows_binned = pimg.nrows_binned;
          img_cfg.ncols        = pimg.ncols;
        } else {
          img_cfg.raw_rows = det.geometry().rows;
          img_cfg.raw_cols = det.geometry().cols;
        }
        img_cfg.row_binning    = pimg.row_binning;
        img_cfg.col_binning    = pimg.col_binning;
        img_cfg.pixel_size_um  = pimg.pixel_size_um;
        img_cfg.sigma_readout_e= pimg.sigma_readout_e;
        img_cfg.lambda_dc      = pimg.lambda_dc;
        img_cfg.rng_seed       = pimg.rng_seed;
        img_cfg.include_dark_current = (pimg.lambda_dc > 0.0);
        auto img_gen = std::make_shared<PatternImageGenerator>(img_cfg, ct);
        const int example_ne = 5;
        auto image_2d = img_gen->GenerateImage(example_ne, Ee_ref_eV);
        TH2D h2img("h2_pattern_image",
                   "2D binned pattern image;col;row",
                   img_gen->Ncols(), 0.5, img_gen->Ncols() + 0.5,
                   img_gen->NrowsBinned(), 0.5, img_gen->NrowsBinned() + 0.5);
        for (int r = 0; r < img_gen->NrowsBinned(); ++r)
          for (int c = 0; c < img_gen->Ncols(); ++c)
            h2img.SetBinContent(c + 1, img_gen->NrowsBinned() - r, image_2d[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)]);
        std::string img_png = run.outdir + "/example_2d_pattern_image.png";
        TCanvas cimg("c_2d", "2D pattern image", 800, 300);
        h2img.Draw("COLZ");
        cimg.SaveAs(img_png.c_str());
        std::cout << "[example-one] 2D pattern image (n_e=" << example_ne << ") written to " << img_png << "\n";
        // Scan middle row for pattern (notebook scan_image style)
        const int mid_row = img_gen->NrowsBinned() / 2;
        const std::vector<double>& middle_row = image_2d[static_cast<std::size_t>(mid_row)];
        auto scan_results = classifier->ScanRow(middle_row);
        const PatternResult* best_scan = nullptr;
        for (const auto& pr : scan_results) {
          if (!pr.valid) continue;
          if (!best_scan || pr.total_charge_e > best_scan->total_charge_e) best_scan = &pr;
        }
        if (best_scan) {
          int code = 0;
          for (int d : best_scan->label.q) code = code * 10 + d;
          std::cout << "[example-one] 2D image middle-row scan: pattern " << code << ", charge " << best_scan->total_charge_e << "\n";
        } else {
          std::cout << "[example-one] 2D image middle-row scan: no pattern found.\n";
        }
        // Optional: simulate_cluster(b,c,d) and classify (notebook simulate_cluster)
        bool do_sim_cluster = (jroot.contains("debug") && jroot["debug"].contains("simulate_cluster") && jroot["debug"]["simulate_cluster"].get<bool>());
        if (do_sim_cluster && img_gen) {
          double sb = 1.0, sc = 1.0, sd = 0.0;
          if (jroot["debug"].contains("simulate_cluster_bcd")) {
            const auto& bcd = jroot["debug"]["simulate_cluster_bcd"];
            if (bcd.is_array() && bcd.size() >= 3) { sb = bcd[0]; sc = bcd[1]; sd = bcd[2]; }
          }
          auto cl = img_gen->SimulateCluster(sb, sc, sd);
          std::vector<double> mid_row_5(cl[1].begin(), cl[1].end());  // middle row, 5 cols
          auto cl_results = classifier->ScanRow(mid_row_5);
          const PatternResult* best_cl = nullptr;
          for (const auto& pr : cl_results) {
            if (!pr.valid) continue;
            if (!best_cl || pr.total_charge_e > best_cl->total_charge_e) best_cl = &pr;
          }
          std::cout << "[example-one] simulate_cluster(" << sb << "," << sc << "," << sd << ") middle row =";
          for (double v : mid_row_5) std::cout << " " << std::fixed << std::setprecision(3) << v;
          std::cout << "\n";
          if (best_cl) {
            int code = 0;
            for (int d : best_cl->label.q) code = code * 10 + d;
            std::cout << "[example-one]   classified pattern " << code << "\n";
          }
        }
        // Optional: charge distribution histos per n_e (Ploting_simulation_specific_E style)
        bool output_charge_histos = (jroot.contains("debug") && jroot["debug"].contains("output_charge_histos") && jroot["debug"]["output_charge_histos"].get<bool>());
        if (output_charge_histos && img_gen) {
          const int ncharge_trials = jroot["debug"].value("output_charge_histos_n", 2000);
          for (int ne = ne_min; ne <= ne_max; ++ne) {
            auto h = std::make_unique<TH1D>(("charge_ne" + std::to_string(ne)).c_str(),
                                           ("Charge distribution n_e=" + std::to_string(ne)).c_str(),
                                           80, 0.0, 8.0);
            for (int t = 0; t < ncharge_trials; ++t) {
              auto img = img_gen->GenerateImage(ne, Ee_ref_eV);
              const int mid_row = img_gen->NrowsBinned() / 2;
              const std::vector<double>& row = img[static_cast<std::size_t>(mid_row)];
              auto res = classifier->ScanRow(row);
              const PatternResult* best = nullptr;
              for (const auto& pr : res) {
                if (!pr.valid) continue;
                if (!best || pr.total_charge_e > best->total_charge_e) best = &pr;
              }
              if (best) h->Fill(best->total_charge_e);
            }
            charge_histos.push_back(std::move(h));
          }
          std::cout << "[example-one] Charge histos per n_e (" << (ne_max - ne_min + 1) << " histos, " << ncharge_trials << " trials each) will be written to ROOT.\n";
        }
        // Optional: sanity-check CSV (event, n_e, pattern_id, total_charge_e)
        bool output_sanity_csv = (jroot.contains("debug") && jroot["debug"].contains("output_sanity_csv") && jroot["debug"]["output_sanity_csv"].get<bool>());
        if (output_sanity_csv && img_gen) {
          const int nsanity = jroot["debug"].value("output_sanity_csv_n", 500);
          sanity_csv_lines.push_back("event,n_e,pattern_id,total_charge_e");
          int ev_id = 0;
          for (int ne = ne_min; ne <= ne_max; ++ne) {
            const int per_ne = std::max(1, nsanity / (ne_max - ne_min + 1));
            for (int t = 0; t < per_ne; ++t) {
              auto img = img_gen->GenerateImage(ne, Ee_ref_eV);
              const int mid_row = img_gen->NrowsBinned() / 2;
              const std::vector<double>& row = img[static_cast<std::size_t>(mid_row)];
              auto res = classifier->ScanRow(row);
              const PatternResult* best = nullptr;
              for (const auto& pr : res) {
                if (!pr.valid) continue;
                if (!best || pr.total_charge_e > best->total_charge_e) best = &pr;
              }
              std::ostringstream line;
              line << ev_id << "," << ne << ",";
              if (best) {
                int code = 0;
                for (int d : best->label.q) code = code * 10 + d;
                line << code << "," << std::fixed << std::setprecision(4) << best->total_charge_e;
              } else {
                line << "0,0";
              }
              sanity_csv_lines.push_back(line.str());
              ++ev_id;
            }
          }
          std::string sanity_path = run.outdir + "/sanity_check.csv";
          std::ofstream sf(sanity_path);
          if (sf.is_open()) {
            for (const auto& l : sanity_csv_lines) sf << l << "\n";
            sf.close();
            std::cout << "[example-one] Sanity CSV written to " << sanity_path << "\n";
          }
        }
      }

      if (run.background_source == "bp_br_template") {
        if (run.background_Bp.size() != summary.pattern_roi.size() ||
            run.background_Br.size() != summary.pattern_roi.size())
          throw std::runtime_error("background_source is bp_br_template but background_Bp/Br size != pattern_roi size");
        B_pat.resize(summary.pattern_roi.size());
        for (size_t i = 0; i < B_pat.size(); ++i)
          B_pat[i] = run.background_Bp[i] + run.background_Br[i];
        std::cout << "[example-one] B_pat from pydme-style Bp+Br template (background_source=bp_br_template).\n";
      } else {
      const std::string& bkg_eff_csv = cfg.backgrounds().background_efficiency_csv;
      if (!bkg_eff_csv.empty()) {
        // Background fold via migration matrix P(identified | true pattern)
        std::map<std::pair<int, int>, double> migration_map;
        std::ifstream mig_in(bkg_eff_csv);
        if (!mig_in.is_open())
          throw std::runtime_error("Failed to open background efficiency CSV: " + bkg_eff_csv);
        std::string header_line;
        if (!std::getline(mig_in, header_line))
          throw std::runtime_error("Empty background efficiency CSV: " + bkg_eff_csv);
        std::vector<int> iden_codes;  // column index -> identified pattern code (from eff_0, eff_1, ...)
        {
          std::stringstream hs(header_line);
          std::string tok;
          bool first = true;
          while (std::getline(hs, tok, ',')) {
            while (!tok.empty() && std::isspace(static_cast<unsigned char>(tok.back()))) tok.pop_back();
            while (!tok.empty() && std::isspace(static_cast<unsigned char>(tok.front()))) tok.erase(tok.begin());
            if (first) { first = false; continue; }
            if (tok.empty() || tok.rfind("eff_", 0) != 0) continue;
            try {
              int code = std::stoi(tok.substr(4));
              iden_codes.push_back(code);
            } catch (...) { continue; }
          }
        }
        std::string row_line;
        while (std::getline(mig_in, row_line)) {
          while (!row_line.empty() && std::isspace(static_cast<unsigned char>(row_line.front())))
            row_line.erase(row_line.begin());
          if (row_line.empty()) continue;
          std::stringstream rs(row_line);
          std::string first_col;
          if (!std::getline(rs, first_col, ',')) continue;
          int true_pat = 0;
          try { true_pat = std::stoi(first_col); } catch (...) { continue; }
          for (size_t j = 0; j < iden_codes.size(); ++j) {
            std::string cell;
            if (!std::getline(rs, cell, ',')) break;
            try {
              double p = std::stod(cell);
              migration_map[{true_pat, iden_codes[j]}] = p;
            } catch (...) {}
          }
        }
        std::cout << "[example-one] Loaded background migration matrix: " << migration_map.size()
                  << " (true_pat, iden_pat) entries from " << bkg_eff_csv << "\n";

        // B_true_pat: rate per true (ideal) pattern. Map B_tot(n_e) to single-pixel patterns (0)-(5).
        std::map<int, double> B_true_pat;
        for (int ne = 0; ne <= 5; ++ne) {
          int bin = B_tot->FindBin(static_cast<double>(ne));
          double rate = B_tot->GetBinContent(bin);
          if (rate != 0.0) B_true_pat[ne] = rate;  // pattern code ne = single-pixel n_e electrons
        }
        // Multi-pixel true patterns could be added here from a separate model; for DC we use single-pixel only.

        B_pat = ccdarksens::FoldBackgroundWithMigration(
            B_true_pat, summary.pattern_roi, migration_map);
        std::cout << "[example-one] B_pat from migration (background_efficiency_csv).\n";
      } else {
        B_pat = ccdarksens::FoldNeToPatternRates(
            *B_tot, ne_min_bkg, ne_max, summary.pattern_roi, pattern_eff_map);
      }
      }
      if (run.background_source != "bp_br_template")
        std::cout << "[example-one] B_tot(n_e) integral = " << std::scientific << B_tot->Integral() << "\n";
      std::cout << "[example-one] B_pat (pattern rates):";
      for (size_t i = 0; i < B_pat.size(); ++i)
        std::cout << " " << summary.pattern_roi[i] << "=" << std::scientific << B_pat[i];
      std::cout << "\n";
    }

    // =========================================================================
    // 3. Compute signal spectrum for the requested (mχ, σe) point
    // =========================================================================
    DMElectronConfig mc = mc_base;
    mc.mchi_MeV = mchi_MeV;
    mc.sigma_e_cm2 = format_sigma(sigma_e_cm2, fmt_sigma);
    DMElectronModel dm_sig;
    if (!dm_sig.Configure(mc)) {
      std::cerr << "[example-one] Failed to configure DM model for mchi=" << mchi_MeV
                << ", sigma=" << sigma_e_cm2 << "\n";
      return 1;
    }
    auto dRdE_sig = dm_sig.MakeSpectrum_E();
    pipe.SetAnalysisSpace(AnalysisSpace::Pattern);
    // ε(pattern|n_true) encodes full detector; exposure always from config (JSON).
    auto S_true = ion->FoldToNe(*dRdE_sig, summary.exposure_kg_year, ne_min, ne_max);
    auto S_obs = pipe.Apply(*dRdE_sig, summary.exposure_kg_year, ne_min, ne_max, Ee_ref_eV);

    std::cout << "\n[example-one] === Single point: mchi = " << mchi_MeV
              << " MeV, sigma = " << format_sigma(sigma_e_cm2, fmt_sigma) << " cm^2 ===\n";
    std::cout << "  S_obs(n_e):\n";
    for (int ne = ne_min; ne <= ne_max; ++ne) {
      int bin = S_obs->FindBin(static_cast<double>(ne));
      double v = S_obs->GetBinContent(bin);
      if (v != 0.0) std::cout << "    n_e = " << ne << "  S_obs = " << v << "\n";
    }
    std::cout << "  Integral S_obs = " << S_obs->Integral() << "\n";

    double q_ts = 0.0;
    if (use_pattern_bins) {
      std::vector<double> S_pat = ccdarksens::FoldNeToPatternRates(
          *S_true, ne_min, ne_max, summary.pattern_roi, pattern_eff_map);
      std::vector<double> model_test(B_pat.size());
      for (size_t i = 0; i < B_pat.size(); ++i) model_test[i] = S_pat[i] + B_pat[i];

      std::vector<double> data_pat;
      if (!run.data_path.empty()) {
        data_pat = load_data(run.data_path, B_pat.size());
        if (data_pat.size() < B_pat.size())
          throw std::runtime_error("data_path \"" + run.data_path + "\" has " + std::to_string(data_pat.size()) + " values, need " + std::to_string(B_pat.size()));
        data_pat.resize(B_pat.size());
      }
      if (!data_pat.empty()) {
        std::cout << "  Pattern bins (pattern_roi | D_pat | B_pat):\n";
        std::cout << std::scientific;
        for (size_t i = 0; i < summary.pattern_roi.size(); ++i) {
          std::cout << "    pattern " << summary.pattern_roi[i]
                    << "  D_pat = " << data_pat[i]
                    << "  B_pat = " << B_pat[i] << "\n";
        }
        std::cout << std::defaultfloat;
      } else {
        std::cout << "  Pattern bins (pattern_roi | S_pat | B_pat | model_test = S+B):\n";
        std::cout << std::scientific;
        for (size_t i = 0; i < summary.pattern_roi.size(); ++i) {
          std::cout << "    pattern " << summary.pattern_roi[i]
                    << "  S_pat = " << S_pat[i]
                    << "  B_pat = " << B_pat[i]
                    << "  model_test = " << model_test[i] << "\n";
        }
        std::cout << std::defaultfloat;
      }

      ccdarksens::stats::StatisticsConfig stats_cfg;
      stats_cfg.test_stat = "PLR";
      stats_cfg.verbosity = 0;
      auto ts = ccdarksens::stats::MakeTestStatistic(stats_cfg);
      if (!ts) throw std::runtime_error("Failed to construct PLR.");
      q_ts = ts->EvaluateRatio(B_pat, model_test, B_pat);
      std::cout << "  q (Asimov PLR) = " << q_ts << "\n";

      if (run.use_profile_likelihood) {
        ccdarksens::stats::ProfileLikelihood pl;
        std::vector<double> S_null, S_test;
        if (run.single_bin_likelihood) {
          double d_sum = 0, b_sum = 0, s_sum = 0;
          for (double v : B_pat) b_sum += v;
          if (!data_pat.empty()) {
            for (std::size_t i = 0; i < data_pat.size(); ++i) d_sum += data_pat[i];
          } else if (run.data_path.empty()) {
            d_sum = b_sum;
          } else {
            auto data = load_data(run.data_path, B_pat.size());
            if (data.size() < B_pat.size())
              throw std::runtime_error("data_path \"" + run.data_path + "\" has " + std::to_string(data.size()) + " values, need " + std::to_string(B_pat.size()));
            for (std::size_t i = 0; i < B_pat.size(); ++i) d_sum += data[i];
          }
            for (double v : S_pat) s_sum += v;
          pl.SetData({d_sum});
          pl.SetBTemplate({b_sum});
          S_null = {0.0};
          S_test = {s_sum};
          if (run.verbosity >= 1)
            std::cout << "  [profile] Single-bin, data_sum=" << d_sum << " B_sum=" << b_sum << "\n";
        } else {
          if (!data_pat.empty()) {
            pl.SetData(data_pat);
            double data_exposure = 0.0;
            if (get_exposure_from_data_file(run.data_path, &data_exposure) && data_exposure > 0.0) {
              const double cfg_exposure = summary.exposure_kg_year;
              const double rel = std::abs(data_exposure - cfg_exposure) / cfg_exposure;
              if (rel > 0.01)
                std::cout << "  [profile] WARNING: data file exposure_kg_year = " << data_exposure
                          << " but config gives " << cfg_exposure << " (relative diff " << (rel * 100) << "%). "
                          << "Match experiment.livetime_days and detector mass for correct S_pat/B_pat.\n";
              else if (run.verbosity >= 1)
                std::cout << "  [profile] Data exposure_kg_year = " << data_exposure << " (matches config).\n";
            }
            if (run.verbosity >= 1)
              std::cout << "  [profile] Using real data from " << run.data_path << "\n";
          } else {
            pl.SetData(B_pat);
            if (run.verbosity >= 1)
              std::cout << "  [profile] Asimov data (data = B)\n";
          }
          S_null.assign(B_pat.size(), 0.0);
          S_test = S_pat;
          std::string bg_model = run.background_model;
          if (bg_model == "Bp_theta_Br" || bg_model == "Bp_theta_br") {
            if (run.background_Bp.size() == B_pat.size() && run.background_Br.size() == B_pat.size()) {
              pl.SetBpBr(run.background_Bp, run.background_Br);
              if (run.verbosity >= 1)
                std::cout << "  [profile] Background: B = Bp + theta*Br (from config)\n";
            } else {
              std::vector<double> z(B_pat.size(), 0.0);
              pl.SetBpBr(z, B_pat);
              if (run.verbosity >= 1)
                std::cout << "  [profile] Background: B = theta*Br, Br = B_pat\n";
            }
            pl.SetConstrainPriorStrength(run.constrain_prior_strength);
            if (run.verbosity >= 1 && run.constrain_prior_strength > 0)
              std::cout << "  [profile] Constrain: pydme L(theta), prior_strength=" << run.constrain_prior_strength << "\n";
          } else {
            pl.SetBTemplate(B_pat);
          }
        }
        if (run.constrain_scale_prior_mean && run.constrain_scale_prior_sigma && *run.constrain_scale_prior_sigma > 0) {
          const double mean = *run.constrain_scale_prior_mean;
          const double sig = *run.constrain_scale_prior_sigma;
          pl.SetConstrain([mean, sig](double scale) {
            return 0.5 * std::pow((scale - mean) / sig, 2);
          });
          if (run.verbosity >= 1)
            std::cout << "  [profile] Constrain: Gaussian prior on scale (mean=" << mean << ", sigma=" << sig << ")\n";
        }
        double pl_lo = 0.01, pl_hi = 10.0;
        if (pl.UseBpBr()) { pl_lo = run.theta_lo; pl_hi = run.theta_hi; }
        const double q_profile = pl.EvaluateRatio(S_null, S_test, pl_lo, pl_hi);
        std::cout << "  q (profile) = " << q_profile << "\n";
      }

      // =====================================================================
      // 4. Write ROOT output
      // =====================================================================
      std::string out_path = run.outdir + "/example_one_point_pattern.root";
      TFile fout(out_path.c_str(), "RECREATE");
      if (!fout.IsOpen()) {
        std::cerr << "[example-one] Cannot create " << out_path << "\n";
      } else {
        S_obs->Write("S_obs_ne");
        if (run.background_source != "bp_br_template")
          B_tot->Write("B_tot_ne");
        const int np = static_cast<int>(summary.pattern_roi.size());
        TH1D h_B_pat("B_pat_validation",
                     "Background rate per pattern;pattern bin;rate",
                     np, 0.5, np + 0.5);
        for (int i = 0; i < np; ++i) {
          h_B_pat.SetBinContent(i + 1, B_pat[i]);
          h_B_pat.GetXaxis()->SetBinLabel(i + 1, std::to_string(summary.pattern_roi[i]).c_str());
        }
        h_B_pat.Write();

        if (!run.data_path.empty()) {
          auto data = load_data(run.data_path, B_pat.size());
          if (data.size() >= static_cast<std::size_t>(np)) {
            data.resize(np);
            TH1D h_D_pat("D_pat", "Observed counts per pattern (data);pattern bin;counts",
                         np, 0.5, np + 0.5);
            for (int i = 0; i < np; ++i) {
              h_D_pat.SetBinContent(i + 1, data[i]);
              h_D_pat.GetXaxis()->SetBinLabel(i + 1, std::to_string(summary.pattern_roi[i]).c_str());
            }
            h_D_pat.Write();
          }
        }
        // Always write S_pat and B_pat (signal from model at example point, for comparison with data)
        TH1D h_S_pat("S_pat_validation",
                     "Signal rate per pattern (example point);pattern bin;rate",
                     np, 0.5, np + 0.5);
        for (int i = 0; i < np; ++i) {
          h_S_pat.SetBinContent(i + 1, S_pat[i]);
          h_S_pat.GetXaxis()->SetBinLabel(i + 1, std::to_string(summary.pattern_roi[i]).c_str());
        }
        h_S_pat.Write();
        for (auto& h : charge_histos)
          if (h) h->Write();

        // Rate-per-pattern plot (PNG): data + B when data_path set, else S + B
        std::string png_path = run.outdir + "/example_one_point_pattern_rates_per_pattern.png";
        TCanvas c1("c_rates", "Rate per pattern", 800, 600);
        if (!run.data_path.empty()) {
          auto data = load_data(run.data_path, B_pat.size());
          if (data.size() >= static_cast<std::size_t>(np)) {
            data.resize(np);
            TH1D h_D_pat("D_pat_plot", "Observed counts per pattern (data);pattern bin;counts",
                         np, 0.5, np + 0.5);
            for (int i = 0; i < np; ++i) {
              h_D_pat.SetBinContent(i + 1, data[i]);
              h_D_pat.GetXaxis()->SetBinLabel(i + 1, std::to_string(summary.pattern_roi[i]).c_str());
            }
            h_D_pat.SetFillColorAlpha(38, 0.7);
            h_D_pat.SetLineColor(38);
            h_B_pat.SetFillColorAlpha(46, 0.7);
            h_B_pat.SetLineColor(46);
            double ymax = std::max(h_D_pat.GetMaximum(), h_B_pat.GetMaximum()) * 1.15;
            if (ymax <= 0.0) ymax = 1.0;
            h_D_pat.GetYaxis()->SetRangeUser(0.0, ymax);
            h_D_pat.Draw("BAR");
            h_B_pat.Draw("BAR SAME");
            TLegend leg(0.70, 0.75, 0.88, 0.88);
            leg.AddEntry(&h_D_pat, "Data (D_{pat})", "f");
            leg.AddEntry(&h_B_pat, "Background (B_{pat})", "f");
            leg.Draw();
          }
        } else {
          TH1D h_S_pat("S_pat_plot", "Signal rate per pattern;pattern bin;rate", np, 0.5, np + 0.5);
          for (int i = 0; i < np; ++i) {
            h_S_pat.SetBinContent(i + 1, S_pat[i]);
            h_S_pat.GetXaxis()->SetBinLabel(i + 1, std::to_string(summary.pattern_roi[i]).c_str());
          }
          h_S_pat.SetFillColorAlpha(38, 0.7);
          h_S_pat.SetLineColor(38);
          h_B_pat.SetFillColorAlpha(46, 0.7);
          h_B_pat.SetLineColor(46);
          double ymax = std::max(h_S_pat.GetMaximum(), h_B_pat.GetMaximum()) * 1.15;
          if (ymax <= 0.0) ymax = 1.0;
          h_S_pat.GetYaxis()->SetRangeUser(0.0, ymax);
          h_S_pat.Draw("BAR");
          h_B_pat.Draw("BAR SAME");
          TLegend leg(0.70, 0.75, 0.88, 0.88);
          leg.AddEntry(&h_S_pat, "Signal (S_{pat})", "f");
          leg.AddEntry(&h_B_pat, "Background (B_{pat})", "f");
          leg.Draw();
        }
        c1.SaveAs(png_path.c_str());
        std::cout << "  Rate-per-pattern plot written to " << png_path << "\n";

        fout.Close();
        std::cout << "  Output written to " << out_path << "\n";
      }
    } else {
      std::vector<double> data, model_null, model_test;
      for (int ne : summary.roi_bins) {
        int bin = S_obs->FindBin(ne);
        double s = S_obs->GetBinContent(bin);
        double b = B_tot->GetBinContent(bin);
        data.push_back(b);
        model_null.push_back(b);
        model_test.push_back(b + s);
      }
      ccdarksens::stats::StatisticsConfig stats_cfg;
      stats_cfg.test_stat = "PLR";
      stats_cfg.verbosity = 0;
      auto ts = ccdarksens::stats::MakeTestStatistic(stats_cfg);
      if (ts) q_ts = ts->EvaluateRatio(data, model_test, model_null);
      std::cout << "  q (n_e ROI PLR) = " << q_ts << "\n";
    }

    std::cout << "[example-one] Done.\n";
  } catch (const std::exception& ex) {
    std::cerr << "[example-one] ERROR: " << ex.what() << "\n";
    return 1;
  }
  return 0;
}
