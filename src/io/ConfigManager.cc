// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ConfigManager.cc -- I parse the framework JSON config into the typed run,
//  detector, experiment, background, timing, response and model structs
//  declared in ConfigManager.hh. Missing optional keys fall back to the
//  defaults declared in those structs.
// ===========================================================================

#include "ccdarksens/io/ConfigManager.hh"
#include <fstream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <cmath>

#include <nlohmann/json.hpp>

namespace ccdarksens {

// using nlohmann::json;

namespace {

// Fills a ClusterFitMCJSON from a "cluster_fit_mc" JSON block. Shared by the
// single-channel response.cluster_fit_mc path and each entry of
// response.channels[] (joint-likelihood channels), so both stay in sync.
void ParseClusterFitMCJSON_(const nlohmann::json& jp, ClusterFitMCJSON& p) {
  p.window_nx       = jp.value("window_nx",       p.window_nx);
  p.window_ny       = jp.value("window_ny",       p.window_ny);
  p.pixel_size_um   = jp.value("pixel_size_um",   p.pixel_size_um);
  p.sigma_readout_e = jp.value("sigma_readout_e", p.sigma_readout_e);
  p.lambda_dc_per_pixel = jp.value("lambda_dc_per_pixel", p.lambda_dc_per_pixel);  // diagnostic dark-current knob (0 = off)
  p.lambda_dc_calib_per_pixel = jp.value("lambda_dc_calib_per_pixel", p.lambda_dc_calib_per_pixel);  // diagnostic: DC used for the calibration toys only (<0 = same as above)
  p.one_dimensional = jp.value("one_dimensional", p.one_dimensional);
  if (jp.contains("diffusion")) {
    const auto& jd = jp.at("diffusion");
    p.A_um2        = jd.value("A_um2",        p.A_um2);
    p.b_umInv      = jd.value("b_umInv",      p.b_umInv);
    p.alpha        = jd.value("alpha",        p.alpha);
    p.beta_per_keV = jd.value("beta_per_keV", p.beta_per_keV);
    p.thickness_um = jd.value("thickness_um", p.thickness_um);
  } else {
    p.A_um2        = jp.value("A_um2",        p.A_um2);
    p.b_umInv      = jp.value("b_umInv",      p.b_umInv);
    p.alpha        = jp.value("alpha",        p.alpha);
    p.beta_per_keV = jp.value("beta_per_keV", p.beta_per_keV);
    p.thickness_um = jp.value("thickness_um", p.thickness_um);
  }
  p.fit_method       = jp.value("fit_method",       p.fit_method);
  p.sigma_xy_lo_px   = jp.value("sigma_xy_lo_px",   p.sigma_xy_lo_px);
  p.sigma_xy_hi_px   = jp.value("sigma_xy_hi_px",   p.sigma_xy_hi_px);
  p.n_toys           = jp.value("n_toys",           p.n_toys);
  p.target_tail_prob = jp.value("target_tail_prob", p.target_tail_prob);
  p.calib_rng_seed   = jp.value("calib_rng_seed",   p.calib_rng_seed);
  p.ne_trials_per_point = jp.value("ne_trials_per_point", p.ne_trials_per_point);
  p.sigma_xy_fid_min_px = jp.value("sigma_xy_fid_min_px", p.sigma_xy_fid_min_px);
  p.sigma_xy_fid_max_px = jp.value("sigma_xy_fid_max_px", p.sigma_xy_fid_max_px);
  p.eh_pair_eV       = jp.value("eh_pair_eV",       p.eh_pair_eV);
  p.fano_factor      = jp.value("fano_factor",      p.fano_factor);
  p.rng_seed         = jp.value("rng_seed",         p.rng_seed);
  p.Etrue_min_eV     = jp.value("Etrue_min_eV",     p.Etrue_min_eV);
  p.Etrue_max_eV     = jp.value("Etrue_max_eV",     p.Etrue_max_eV);
  p.Etrue_npoints    = jp.value("Etrue_npoints",    p.Etrue_npoints);
  p.Ereco_min_eV     = jp.value("Ereco_min_eV",     p.Ereco_min_eV);
  p.Ereco_max_eV     = jp.value("Ereco_max_eV",     p.Ereco_max_eV);
  p.Ereco_nbins      = jp.value("Ereco_nbins",      p.Ereco_nbins);
}

}  // namespace

// ----------------------------------------------------------------------------
// ConfigManager::ConfigManager
//   I only remember the path; the file is read in parse().
// ----------------------------------------------------------------------------
ConfigManager::ConfigManager(std::string path) : path_(std::move(path)) {}

// ----------------------------------------------------------------------------
// ConfigManager::parse
//   I read the JSON file and hand each top-level block to its parser.
//   Only "run" is mandatory; detector, experiment, backgrounds, response and
//   model are parsed if present (parse_backgrounds_ also fills the timing
//   settings). Throws std::runtime_error if the file cannot be opened, and
//   nlohmann::json exceptions on malformed JSON or a missing required key.
// ----------------------------------------------------------------------------
void ConfigManager::parse() {
  std::ifstream in(path_);
  if (!in) throw std::runtime_error("Cannot open config: " + path_);
  nlohmann::json j;
  in >> j;

  parse_run_(j.at("run"));
  if (j.contains("detector"))   parse_detector_(j.at("detector"));
  if (j.contains("experiment")) parse_experiment_(j.at("experiment"));
  if (j.contains("backgrounds")) parse_backgrounds_(j.at("backgrounds"));
  if (j.contains("response"))   parse_response_(j.at("response"));
  if (j.contains("model"))      parse_model_(j.at("model"));
}

// ----------------------------------------------------------------------------
// ConfigManager::parse_run_
//   I fill the RunHeader from the "run" block: label/output, confidence level,
//   the likelihood options (profile likelihood, minimizer, constraints,
//   background_source / background_model with their Bp/Br vectors), the
//   pydme-matching switches, and the band-tool "mode" with its
//   threshold_toys sub-block. The Gaussian scale prior is only enabled when
//   both its mean and its sigma are given.
// ----------------------------------------------------------------------------
void ConfigManager::parse_run_(const nlohmann::json& j) {
  run_.label     = j.value("label", std::string{});
  run_.outdir    = j.value("outdir", std::string{});
  run_.cl        = j.value("cl", 0.90);
  run_.test_stat = j.value("test_stat", "PLR");
  run_.n_toys    = j.value("n_toys", 0);
  run_.rng_seed  = j.value("rng_seed", 12345ULL);
  run_.verbosity = j.value("verbosity", 1);
  run_.exposure_kg_year = j.value("exposure_kg_year", 0.0);
  run_.dump_point_spectra_root = j.value("dump_point_spectra_root", false);
  run_.use_profile_likelihood = j.value("use_profile_likelihood", false);
  run_.data_path = j.value("data_path", std::string{});
  run_.single_bin_likelihood = j.value("single_bin_likelihood", false);
  if (j.contains("constrain_scale_prior_mean") && j.contains("constrain_scale_prior_sigma")) {
    run_.constrain_scale_prior_mean = j.at("constrain_scale_prior_mean").get<double>();
    run_.constrain_scale_prior_sigma = j.at("constrain_scale_prior_sigma").get<double>();
  }
  run_.background_source = j.value("background_source", std::string("dc_flat_migration"));
  run_.background_model = j.value("background_model", std::string("scale"));
  if (j.contains("background_Bp") && j["background_Bp"].is_array()) {
    run_.background_Bp.clear();
    for (const auto& v : j["background_Bp"]) run_.background_Bp.push_back(v.get<double>());
  }
  if (j.contains("background_Br") && j["background_Br"].is_array()) {
    run_.background_Br.clear();
    for (const auto& v : j["background_Br"]) run_.background_Br.push_back(v.get<double>());
  }
  run_.constrain_prior_strength = j.value("constrain_prior_strength", 0.0);
  run_.constrain_use_gamma_sign = j.value("constrain_use_gamma_sign", false);
  run_.constrain_use_tau_weighted = j.value("constrain_use_tau_weighted", false);
  run_.constrain_n_bins = j.value("constrain_n_bins", 1);
  run_.theta_lo = j.value("theta_lo", 0.5);
  run_.theta_hi = j.value("theta_hi", 10.0);
  run_.profile_likelihood_plot_mchi = j.value("profile_likelihood_plot_mchi", 0.0);
  run_.profile_minimizer = j.value("profile_minimizer", std::string("brent"));
  run_.pydme_style_ul = j.value("pydme_style_ul", false);
  run_.smooth_ul_envelope = j.value("smooth_ul_envelope", false);
  run_.mode = j.value("mode", std::string("scan"));
  run_.q_target_lookup_path = j.value("q_target_lookup_path", std::string{});
  if (j.contains("threshold_toys")) {
    const auto& jt = j.at("threshold_toys");
    run_.threshold_toys.sigma_threshold_graph_path =
        jt.value("sigma_threshold_graph_path", std::string{});
    run_.threshold_toys.n_threshold_toys = jt.value("n_threshold_toys", static_cast<long long>(10000));
    run_.threshold_toys.percentile = jt.value("percentile", 0.90);
    run_.threshold_toys.rng_seed = jt.value("rng_seed", static_cast<uint64_t>(23456ULL));
  }
}

// ----------------------------------------------------------------------------
// ConfigManager::parse_detector_
//   I build the Detector from the "detector" block. rows, cols, pixel_size_um,
//   thickness_mm, target_element and density_g_cm3 are required; active_fraction
//   (default 1), Z, A and a mass_kg override (null/absent = compute the mass
//   from the geometry) are optional.
// ----------------------------------------------------------------------------
void ConfigManager::parse_detector_(const nlohmann::json& j) {
  DetectorGeometry g;
  g.rows = j.at("rows").get<int>();
  g.cols = j.at("cols").get<int>();
  g.pixel_size_um = j.at("pixel_size_um").get<double>();
  g.thickness_mm  = j.at("thickness_mm").get<double>();
  g.active_fraction = j.value("active_fraction", 1.0);

  TargetMaterial m;
  m.element = j.at("target_element").get<std::string>();
  m.density_g_cm3 = j.at("density_g_cm3").get<double>();
  m.Z = j.value("Z", 0);
  m.A = j.value("A", 0.0);

  std::optional<double> mass_override;
  if (j.contains("mass_kg") && !j.at("mass_kg").is_null())
    mass_override = j.at("mass_kg").get<double>();

  detector_ = std::make_unique<Detector>(g, m, mass_override);
}

// Convert experiment.mode ("observed" | "asimov" | "toys") into the enum; anything else throws std::invalid_argument.
static ExperimentMode parse_mode(const std::string& s) {
  if (s == "observed") return ExperimentMode::Observed;
  if (s == "asimov")   return ExperimentMode::Asimov;
  if (s == "toys")     return ExperimentMode::Toys;
  throw std::invalid_argument("experiment.mode must be observed/asimov/toys");
}

// ----------------------------------------------------------------------------
// ConfigManager::parse_experiment_
//   I fill the ExperimentConfig from the "experiment" block. mode, livetime_days
//   and binning{ne_min, ne_max} are required; duty_cycle, roi_bins, pattern_roi
//   and observable_bins are optional.
// ----------------------------------------------------------------------------
void ConfigManager::parse_experiment_(const nlohmann::json& j) {
  exp_cfg_.mode = parse_mode(j.at("mode").get<std::string>());
  exp_cfg_.livetime_days = j.at("livetime_days").get<double>();
  exp_cfg_.duty_cycle    = j.value("duty_cycle", 1.0);

  BinningNE b;
  const auto& jb = j.at("binning");
  b.ne_min = jb.at("ne_min").get<int>();
  b.ne_max = jb.at("ne_max").get<int>();
  exp_cfg_.binning = b;

  if (j.contains("roi_bins"))
    exp_cfg_.roi_bins = j.at("roi_bins").get<std::vector<int>>();
  if (j.contains("pattern_roi"))
    exp_cfg_.pattern_roi = j.at("pattern_roi").get<std::vector<int>>();
  if (j.contains("observable_bins"))
    exp_cfg_.observable_bins = j.at("observable_bins").get<std::string>();
}

// ----------------------------------------------------------------------------
// ConfigManager::parse_response_
//   I fill the ResponseJSON from the "response" block:
//     - mode and analysis_space (the space defaults to "pcd" for mode "pcd",
//       otherwise "pattern"; "cluster_energy" is set explicitly for the WIMP channel);
//     - the pcd block, and the legacy cluster_mc block (kept only for old configs);
//     - efficiency_mc (the old key pattern_mc is still accepted);
//     - cluster_fit_mc, the background-efficiency CSV and the optional joint
//       "channels" list of the WIMP-nucleon channel;
//     - pattern_classifier, whose shared values default to those of efficiency_mc;
//     - charge_ionization and pattern_image.
//   If there is no efficiency_mc/pattern_mc block, I copy the legacy cluster_mc
//   values into emc so that old configs keep working.
// ----------------------------------------------------------------------------
void ConfigManager::parse_response_(const nlohmann::json& jr) {
  // Read detector-response backend mode
  response_.mode = jr.value("mode", std::string("fast"));
  // Decide default analysis_space from mode:
  //   mode == "pcd"  → default to PCD space
  //   otherwise      → default to pattern/n_e space
  const std::string default_space =
      (response_.mode == "pcd") ? std::string("pcd") : std::string("pattern");

  // Allow optional override in JSON: "analysis_space": "pattern" | "pcd"
  response_.analysis_space = jr.value("analysis_space", default_space);

  // --- PCD block: configure q-binning and MC trials for P(q|n_e) ---
  if (jr.contains("pcd")) {
    const auto& jp = jr.at("pcd");
    response_.pcd.q_min       = jp.value("q_min",       response_.pcd.q_min);
    response_.pcd.q_max       = jp.value("q_max",       response_.pcd.q_max);
    response_.pcd.nbins       = jp.value("nbins",       response_.pcd.nbins);
    response_.pcd.mc_trials   = jp.value("mc_trials",   response_.pcd.mc_trials);
    response_.pcd.sigma_res_e = jp.value("sigma_res_e", response_.pcd.sigma_res_e);
    response_.pcd.Dqmin       = jp.value("Dqmin",       response_.pcd.Dqmin);
    response_.pcd.Dqmax       = jp.value("Dqmax",       response_.pcd.Dqmax);
  }

  if (jr.contains("cluster_mc")) {
    const auto& jc = jr.at("cluster_mc");
    response_.cmc.n_events_per_ne = jc.value("n_events_per_ne", 20000);
    response_.cmc.sigma_readout_e = jc.value("sigma_readout_e", 0.16);
    response_.cmc.Qmin_e          = jc.value("Qmin_e", 0.3);
    response_.cmc.Qmax_e          = jc.value("Qmax_e", 4.0);
    if (jc.contains("diffusion")) {
      const auto& jd = jc.at("diffusion");
      response_.cmc.A_um2       = jd.value("A_um2", 803.25);
      response_.cmc.b_umInv     = jd.value("b_umInv", 6.5e-4);
      response_.cmc.alpha       = jd.value("alpha", 1.0);
      response_.cmc.beta_per_keV= jd.value("beta_per_keV", 0.0);
    }
    if (jc.contains("binning")) {
      const auto& jb = jc.at("binning");
      response_.cmc.rows_bin = jb.value("rows_bin", 1);
      response_.cmc.cols_bin = jb.value("cols_bin", 1);
    }
    response_.cmc.pileup_with_dc = jc.value("pileup_with_dc", false);
    response_.cmc.rng_seed       = jc.value("rng_seed", static_cast<uint64_t>(987654321ULL));
  }

  // --- EfficiencyMC block (canonical key: efficiency_mc; legacy alias: pattern_mc) ---
  const nlohmann::json* j_emc = nullptr;
  if (jr.contains("efficiency_mc"))
    j_emc = &jr.at("efficiency_mc");
  else if (jr.contains("pattern_mc"))
    j_emc = &jr.at("pattern_mc");

  if (j_emc) {
    const auto& jp = *j_emc;
    auto& p = response_.emc;

    p.n_events_per_ne = jp.value("n_events_per_ne", p.n_events_per_ne);
    p.sigma_readout_e = jp.value("sigma_readout_e", p.sigma_readout_e);
    p.Qmin_e          = jp.value("Qmin_e",          p.Qmin_e);
    p.Qmax_e          = jp.value("Qmax_e",          p.Qmax_e);
    if (jp.contains("diffusion")) {
      const auto& jd = jp.at("diffusion");
      p.A_um2        = jd.value("A_um2",        p.A_um2);
      p.b_umInv      = jd.value("b_umInv",      p.b_umInv);
      p.alpha        = jd.value("alpha",        p.alpha);
      p.beta_per_keV = jd.value("beta_per_keV", p.beta_per_keV);
    } else {
      p.A_um2        = jp.value("A_um2",        p.A_um2);
      p.b_umInv      = jp.value("b_umInv",      p.b_umInv);
      p.alpha        = jp.value("alpha",        p.alpha);
      p.beta_per_keV = jp.value("beta_per_keV", p.beta_per_keV);
    }
    p.enable_MN       = jp.value("enable_MN",       p.enable_MN);
    p.enable_MNL      = jp.value("enable_MNL",      p.enable_MNL);
    p.rng_seed        = jp.value("rng_seed",        p.rng_seed);
    p.row_length = jp.value("row_length", p.row_length);
    p.use_2d_image_efficiency = jp.value("use_2d_image_efficiency", p.use_2d_image_efficiency);
    p.efficiency_csv           = jp.value("efficiency_csv",           p.efficiency_csv);
    p.efficiency_csv_reference = jp.value("efficiency_csv_reference", p.efficiency_csv_reference);
    // accepted_labels are derived from experiment.pattern_roi in the apps
  }
  // --- ClusterFitMC block (WIMP-nucleus SI channel Phase 3/4) ---
  if (jr.contains("cluster_fit_mc")) {
    ParseClusterFitMCJSON_(jr.at("cluster_fit_mc"), response_.cluster_fit_mc);
  }
  response_.background_efficiency_csv =
      jr.value("background_efficiency_csv", response_.background_efficiency_csv);
  // --- Joint-likelihood channels (e.g. 1x1 + 1x100); empty by default ---
  if (jr.contains("channels")) {
    for (const auto& jc : jr.at("channels")) {
      ResponseChannelJSON ch;
      ch.label = jc.value("label", std::string());
      ch.flat_background_norm_per_kg_year_keV =
          jc.value("flat_background_norm_per_kg_year_keV", 0.0);
      ch.background_efficiency_csv = jc.value("background_efficiency_csv", std::string());
      if (jc.contains("cluster_fit_mc")) {
        ParseClusterFitMCJSON_(jc.at("cluster_fit_mc"), ch.cluster_fit_mc);
      }
      response_.channels.push_back(std::move(ch));
    }
  }

  // --- Pattern Classifier: inherit shared params from efficiency_mc when not set ---
  if (jr.contains("pattern_classifier")) {
    const auto& jc = jr.at("pattern_classifier");
    auto& pc = response_.pattern_classifier;
    // Shared with efficiency_mc: default from emc so config can specify once under efficiency_mc
    pc.Qmin_e          = jc.value("Qmin_e",          response_.emc.Qmin_e);
    pc.neighbor_Qmax_e  = jc.value("neighbor_Qmax_e", pc.neighbor_Qmax_e);
    pc.Qmax_e          = jc.value("Qmax_e",          response_.emc.Qmax_e);
    pc.sigma_res_e     = jc.value("sigma_res_e",    response_.emc.sigma_readout_e);
    pc.enable_MN       = jc.value("enable_MN",       response_.emc.enable_MN);
    pc.enable_MNL      = jc.value("enable_MNL",      response_.emc.enable_MNL);
    pc.max_e_per_pixel = jc.value("max_e_per_pixel", pc.max_e_per_pixel);
    pc.thr_M           = jc.value("thr_M",           pc.thr_M);
    pc.thr_MN          = jc.value("thr_MN",           pc.thr_MN);
    pc.thr_MNL         = jc.value("thr_MNL",         pc.thr_MNL);
    pc.allow_pattern_zero = jc.value("allow_pattern_zero", pc.allow_pattern_zero);
    pc.single_pixel_use_round = jc.value("single_pixel_use_round", pc.single_pixel_use_round);
  }

  // 2D image size and binning for pattern efficiency (notebook-style)
  if (jr.contains("charge_ionization")) {
    const auto& ji = jr.at("charge_ionization");
    auto& ci = response_.charge_ionization;
    ci.table_csv   = ji.value("table_csv", ci.table_csv);
    ci.band_gap_eV = ji.value("band_gap_eV", ci.band_gap_eV);
    ci.eh_pair_eV  = ji.value("eh_pair_eV", ci.eh_pair_eV);
    if (ji.contains("scenario"))
      ci.scenario = ji.at("scenario").get<std::string>();
  }

  if (jr.contains("pattern_image")) {
    const auto& pi = jr.at("pattern_image");
    auto& pimg = response_.pattern_image;
    pimg.nrows_binned   = pi.value("nrows_binned",   pimg.nrows_binned);
    pimg.ncols          = pi.value("ncols",         pimg.ncols);
    pimg.row_binning    = pi.value("row_binning",   pimg.row_binning);
    pimg.col_binning    = pi.value("col_binning",   pimg.col_binning);
    pimg.pixel_size_um  = pi.value("pixel_size_um", pimg.pixel_size_um);
    pimg.sigma_readout_e= pi.value("sigma_readout_e", pimg.sigma_readout_e);
    pimg.lambda_dc      = pi.value("lambda_dc",     pimg.lambda_dc);
    pimg.rng_seed       = pi.value("rng_seed",      static_cast<uint64_t>(pimg.rng_seed));
  }

  if (!j_emc) {
    // Copy from ClusterMC for backward compatibility when no efficiency_mc / pattern_mc block
    response_.emc.n_events_per_ne = response_.cmc.n_events_per_ne;
    response_.emc.sigma_readout_e = response_.cmc.sigma_readout_e;
    response_.emc.Qmin_e          = response_.cmc.Qmin_e;
    response_.emc.Qmax_e          = response_.cmc.Qmax_e;
    response_.emc.A_um2           = response_.cmc.A_um2;
    response_.emc.b_umInv         = response_.cmc.b_umInv;
    response_.emc.alpha           = response_.cmc.alpha;
    response_.emc.beta_per_keV    = response_.cmc.beta_per_keV;
    response_.emc.rng_seed        = response_.cmc.rng_seed;
  }
  




}


// Earlier version of parse_backgrounds_ (per-exposure dark-current key). I keep it commented out only for reference; the active version is below.
// void ConfigManager::parse_backgrounds_(const nlohmann::json& jb) {
//   if (jb.contains("dark_current")) {
//     const auto& jd = jb.at("dark_current");
//     bkg_.lambda_e_per_pix_per_exposure = jd.value("lambda_e_per_pix_per_exposure", 0.0);
//     bkg_.norm_scale = jd.value("norm_scale", 1.0);
//   }
//   if (jb.contains("pattern_efficiency")) {
//     const auto& je = jb.at("pattern_efficiency");
//     const std::string type = je.value("type","flat");
//     if (type == "flat") {
//       bkg_.has_flat_eps = true;
//       bkg_.flat_eps = je.value("epsilon", 1.0);
//     }
//   }

//   // timing knobs live under experiment (but expose via timing_)
//   // Find them from the already-parsed experiment block in the original nlohmann::json:
//   // (we still have j in this scope? If not, read from jb's parent; simplest: require backgrounds to also include timing)
//   // For clarity now, we’ll pull from a sibling "experiment" via the stored exp_cfg_ — but we need exposure_time_s & n_exposures.
//   // Easiest: require they sit under "experiment" and re-open here:
//   // (I could move this logic to parse_experiment_ later.)

//   // Try parent: we can't access parent here; instead, read from a convenience duplication in backgrounds if present:
//   if (jb.contains("timing")) {
//     const auto& jt = jb.at("timing");
//     if (jt.contains("n_exposures") && !jt.at("n_exposures").is_null())
//       timing_.n_exposures_override = jt.at("n_exposures").get<int>();
//     timing_.exposure_time_s = jt.value("exposure_time_s", 0.0);
//   }
// }

// ----------------------------------------------------------------------------
// ConfigManager::parse_backgrounds_
//   I fill the background and timing settings from the "backgrounds" block:
//     - dark_current: lambda_e_per_pix_per_year (preferred) or the older
//       lambda_e_per_pix_per_exposure, which I convert to a yearly rate with
//       the exposure time; plus norm_scale;
//     - pattern_efficiency (flat type only);
//     - timing (exposure_time_s and an optional n_exposures override);
//     - flat_background (dR/dE in events/(kg*year*keV) between Emin and Emax);
//     - background_efficiency_csv;
//     - lee: the low-energy excess used by the WIMP-nucleon channel, with two
//       published presets (damic_snolab_2020, damic_snolab_2023_skipper) that
//       explicit rate/decay-energy keys override.
// ----------------------------------------------------------------------------
void ConfigManager::parse_backgrounds_(const nlohmann::json& jb) {
  if (jb.contains("dark_current")) {
    const auto& jd = jb.at("dark_current");

    // New preferred key: lambda_e_per_pix_per_year
    if (jd.contains("lambda_e_per_pix_per_year")) {
      bkg_.lambda_e_per_pix_per_year =
        jd.value("lambda_e_per_pix_per_year", 0.0);
    } else if (jd.contains("lambda_e_per_pix_per_exposure")) {
      // Optional backward compatibility
      const double lam_exp =
        jd.value("lambda_e_per_pix_per_exposure", 0.0);
      const double year_s = 365.25 * 86400.0;
      const double exp_s  = jb.contains("timing")
                            ? jb.at("timing").value("exposure_time_s", 0.0)
                            : 0.0;
      if (exp_s > 0.0) {
        bkg_.lambda_e_per_pix_per_year = lam_exp * (year_s / exp_s);
      } else {
        bkg_.lambda_e_per_pix_per_year = 0.0;
      }
    } else {
      bkg_.lambda_e_per_pix_per_year = 0.0;
    }

    bkg_.norm_scale = jd.value("norm_scale", 1.0);
  }

  if (jb.contains("pattern_efficiency")) {
    const auto& je = jb.at("pattern_efficiency");
    const std::string type = je.value("type","flat");
    if (type == "flat") {
      bkg_.has_flat_eps = true;
      bkg_.flat_eps = je.value("epsilon", 1.0);
    }
  }

  if (jb.contains("timing")) {
    const auto& jt = jb.at("timing");
    if (jt.contains("n_exposures") && !jt.at("n_exposures").is_null())
      timing_.n_exposures_override = jt.at("n_exposures").get<int>();
    timing_.exposure_time_s = jt.value("exposure_time_s", 0.0);
  }

  // ---- flat background (energy spectrum) ----
  //
  // JSON block:
  //   "flat_background": {
  //     "norm_per_kg_year_keV": 1.0,
  //     "Emin_eV": 0.0,
  //     "Emax_eV": 20000.0,
  //     "nbins": 200
  //   }
  //
  if (jb.contains("flat_background")) {
    const auto& jf = jb.at("flat_background");
    bkg_.has_flat_bkg = true;

    bkg_.flat_bkg_norm_per_kg_year =
        jf.value("norm_per_kg_year_keV", 0.0);  // events/(kg·year·keV)
    bkg_.flat_bkg_Emin_eV = jf.value("Emin_eV", 0.0);
    bkg_.flat_bkg_Emax_eV = jf.value("Emax_eV", 0.0);
    bkg_.flat_bkg_nbins   = jf.value("nbins", 0);
  } else {
    bkg_.has_flat_bkg = false;
    bkg_.flat_bkg_norm_per_kg_year = 0.0;
    bkg_.flat_bkg_Emin_eV = 0.0;
    bkg_.flat_bkg_Emax_eV = 0.0;
    bkg_.flat_bkg_nbins   = 0;
  }

  if (jb.contains("background_efficiency_csv")) {
    bkg_.background_efficiency_csv = jb.at("background_efficiency_csv").get<std::string>();
  } else {
    bkg_.background_efficiency_csv.clear();
  }

  // ---- low-energy excess (LEE) background ----
  //
  // JSON block:
  //   "lee": {
  //     "enabled": true,
  //     "preset": "damic_snolab_2023_skipper",   // or "damic_snolab_2020"
  //     "rate_per_kg_day": 10.0,                  // overrides preset
  //     "decay_energy_eV": 89.0                   // overrides preset
  //   }
  //
  // Presets carry the published fits (see docs/LowEnergyExcess_Design.md):
  // damic_snolab_2020 = arXiv:2007.15622 (5.1 events/kg/day, eps=67 eV);
  // damic_snolab_2023_skipper = arXiv:2306.01717 (10.0, eps=89 eV), default.
  // Rate is per kg-day of DAMIC's own silicon detector with unconfirmed
  // physical origin -- using this for a non-DAMIC target material borrows
  // the shape, it is not a derived scaling.
  if (jb.contains("lee")) {
    const auto& jl = jb.at("lee");
    bkg_.has_lee_bkg = jl.value("enabled", false);

    const std::string preset = jl.value("preset", std::string("damic_snolab_2023_skipper"));
    double preset_rate = 10.0, preset_eps = 89.0;  // damic_snolab_2023_skipper
    if (preset == "damic_snolab_2020") {
      preset_rate = 5.1;
      preset_eps  = 67.0;
    }
    bkg_.lee_rate_per_kg_day = jl.value("rate_per_kg_day", preset_rate);
    bkg_.lee_decay_energy_eV = jl.value("decay_energy_eV", preset_eps);
  } else {
    bkg_.has_lee_bkg = false;
    bkg_.lee_rate_per_kg_day = 0.0;
    bkg_.lee_decay_energy_eV = 0.0;
  }
}


// ---- NEW: helper to expand QE-Dark-style axis spec into a list of doubles ----
//
// Accepts any of:
//   - plain array: [0.5, 1.0, 2.0]
//   - {"values": [...]}
//   - {"linspace": { "start": ..., "stop": ..., "num": N, "endpoint": true/false }}
//   - {"logspace": { "start_exp": ..., "stop_exp": ..., "num": N, "endpoint": true/false }}
//
// ----------------------------------------------------------------------------
// expand_axis_
//   I expand one grid-axis spec into an explicit list of numbers. Accepted
//   forms are a plain array, {"values": [...]}, {"linspace": {...}} and
//   {"logspace": {...}} (log10 exponents start_exp/stop_exp, or literal
//   start/stop > 0). Values keep the order in which they were appended (no sorting
//   or de-duplication). Throws std::runtime_error for a non-positive literal
//   logspace bound.
// ----------------------------------------------------------------------------
static void expand_axis_(const nlohmann::json& jaxis,
                         std::vector<double>& out_values)
{
  out_values.clear();

  auto append_values = [&](const nlohmann::json& arr) {
    for (const auto& x : arr) {
      if (!x.is_number()) continue;
      out_values.push_back(x.get<double>());
    }
  };

  // Plain array
  if (jaxis.is_array()) {
    append_values(jaxis);
  }

  // "values": [...]
  if (jaxis.contains("values")) {
    append_values(jaxis.at("values"));
  }

  // "linspace": {start, stop, num, endpoint}
  if (jaxis.contains("linspace")) {
    const auto& jl    = jaxis.at("linspace");
    const double a    = jl.at("start").get<double>();
    const double b    = jl.at("stop").get<double>();
    const int    num  = jl.at("num").get<int>();
    const bool   endp = jl.value("endpoint", true);
    if (num > 0) {
      if (num == 1) {
        out_values.push_back(a);
      } else {
        const double span = (b - a);
        const double step = endp ? (span / (num - 1)) : (span / num);
        for (int i = 0; i < num; ++i) {
          out_values.push_back(a + i * step);
        }
      }
    }
  }

  // "logspace": {start_exp, stop_exp, num, endpoint}
  if (jaxis.contains("logspace")) {
    const auto& jl    = jaxis.at("logspace");
    double aexp = 0.0;
    double bexp = 0.0;
    // Support both:
    //  - logspace in log10 space: {start_exp, stop_exp}
    //  - logspace in literal space (e.g. MeV): {start, stop}
    if (jl.contains("start_exp") && jl.contains("stop_exp")) {
      aexp = jl.at("start_exp").get<double>();
      bexp = jl.at("stop_exp").get<double>();
    } else {
      const double start = jl.at("start").get<double>();
      const double stop  = jl.at("stop").get<double>();
      if (start <= 0.0 || stop <= 0.0) {
        throw std::runtime_error("ConfigManager logspace(start/stop) requires start,stop > 0");
      }
      aexp = std::log10(start);
      bexp = std::log10(stop);
    }
    const int    num  = jl.at("num").get<int>();
    const bool   endp = jl.value("endpoint", true);
    if (num > 0) {
      if (num == 1) {
        out_values.push_back(std::pow(10.0, aexp));
      } else {
        const double span_exp = (bexp - aexp);
        const double step_exp = endp ? (span_exp / (num - 1)) : (span_exp / num);
        for (int i = 0; i < num; ++i) {
          const double e = aexp + i * step_exp;
          out_values.push_back(std::pow(10.0, e));
        }
      }
    }
  }

  // We keep the values in the order they were appended; no sorting or dedup.
}



// ----------------------------------------------------------------------------
// ConfigManager::parse_model_
//   I fill the ModelJSON from the "model" block; every key has a default so
//   I never throw here. I read the DM-electron / dark-photon fields, the DM-nucleon
//   fields (target nucleus, A, Z, sigma_n_cm2) and epsilon_ref, plus the optional
//   grid: the mchi_MeV and sigma_e_cm2 axes and the "options.format" hints.
//   NOTE: the sigma_n_cm2 and epsilon axes are not expanded here; the scan apps
//   read those directly from the raw JSON (ExpandAxisFromGrid).
// ----------------------------------------------------------------------------
void ConfigManager::parse_model_(const nlohmann::json& jm) {
  // Use .value(...) with defaults so we never throw if keys are missing.
  model_.type              = jm.value("type", "dm_electron");
  model_.material          = jm.value("material", "Si");
  model_.mediator          = jm.value("mediator", "heavy");
  model_.rates_dir         = jm.value("rates_dir", "data/qedark_rates/Si/heavy");
  model_.filename_template = jm.value("filename_template",
                           "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv");
  model_.mchi_MeV          = jm.value("mchi_MeV", 3.0);
  model_.sigma_e_cm2       = jm.value("sigma_e_cm2", std::string("1e-37"));
  model_.Emin_eV           = jm.value("Emin_eV", 0.0);
  model_.Emax_eV           = jm.value("Emax_eV", 20.0);
  model_.nbins             = jm.value("nbins", 200);

  // DM-nucleon coupling fields (used when type == "migdal")
  model_.target_nucleus    = jm.value("target_nucleus", std::string("Si28"));
  model_.nuclear_A         = jm.value("nuclear_A", 28);
  model_.nuclear_Z         = jm.value("nuclear_Z", 14);
  model_.sigma_n_cm2       = jm.value("sigma_n_cm2", std::string("1e-40"));
  model_.epsilon_ref       = jm.value("epsilon_ref", std::string(""));

  // ---- NEW: optional QE-Dark style grid ----
  model_.has_grid = false;
  model_.grid_mchi.values.clear();
  model_.grid_mchi.format.clear();
  model_.grid_sigma.values.clear();
  model_.grid_sigma.format.clear();

  if (jm.contains("grid")) {
    const auto& jg = jm.at("grid");
    model_.has_grid = true;

    if (jg.contains("mchi_MeV")) {
      expand_axis_(jg.at("mchi_MeV"), model_.grid_mchi.values);
    }
    if (jg.contains("sigma_e_cm2")) {
      expand_axis_(jg.at("sigma_e_cm2"), model_.grid_sigma.values);
    }

    // Optional: reuse QE-Dark formatting hints if present:
    //   "options": { "format": { "mchi": ".6f", "sigma": ".1e" } }
    if (jm.contains("options")) {
      const auto& jo = jm.at("options");
      if (jo.contains("format")) {
        const auto& jf = jo.at("format");
        model_.grid_mchi.format  = jf.value("mchi", std::string{});
        model_.grid_sigma.format = jf.value("sigma", std::string{});
      }
    }
  }
}


} // namespace ccdarksens
