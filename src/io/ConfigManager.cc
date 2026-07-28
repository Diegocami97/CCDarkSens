// ============================================================================
//  CCDarkSens — ConfigManager
//  Parses the framework JSON config into typed run, detector, experiment, background, response, and model structures.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include "ccdarksens/io/ConfigManager.hh"
#include <fstream>
#include <memory>
#include <stdexcept>
#include <utility>
#include <cmath>

#include <nlohmann/json.hpp>

namespace ccdarksens {

// using nlohmann::json;

ConfigManager::ConfigManager(std::string path) : path_(std::move(path)) {}

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
}

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

static ExperimentMode parse_mode(const std::string& s) {
  if (s == "observed") return ExperimentMode::Observed;
  if (s == "asimov")   return ExperimentMode::Asimov;
  if (s == "toys")     return ExperimentMode::Toys;
  throw std::invalid_argument("experiment.mode must be observed/asimov/toys");
}

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
    if (jp.contains("binning")) {
      const auto& jb = jp.at("binning");
      p.rows_bin = jb.value("rows_bin", p.rows_bin);
      p.cols_bin = jb.value("cols_bin", p.cols_bin);
    } else {
      p.rows_bin = jp.value("rows_bin", p.rows_bin);
      p.cols_bin = jp.value("cols_bin", p.cols_bin);
    }
    p.half_window_pix = jp.value("half_window_pix", p.half_window_pix);
    p.pileup_with_dc  = jp.value("pileup_with_dc",  p.pileup_with_dc);
    p.enable_MN       = jp.value("enable_MN",       p.enable_MN);
    p.enable_MNL      = jp.value("enable_MNL",      p.enable_MNL);
    p.rng_seed        = jp.value("rng_seed",        p.rng_seed);
    p.row_length = jp.value("row_length", p.row_length);
    p.use_2d_image_efficiency = jp.value("use_2d_image_efficiency", p.use_2d_image_efficiency);
    // accepted_labels are derived from experiment.pattern_roi in the apps
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
    response_.emc.rows_bin        = response_.cmc.rows_bin;
    response_.emc.cols_bin        = response_.cmc.cols_bin;
    response_.emc.pileup_with_dc  = response_.cmc.pileup_with_dc;
    response_.emc.rng_seed        = response_.cmc.rng_seed;
  }
  




}


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
//   // (If you want, we can move this logic to parse_experiment_ later.)

//   // Try parent: we can't access parent here; instead, read from a convenience duplication in backgrounds if present:
//   if (jb.contains("timing")) {
//     const auto& jt = jb.at("timing");
//     if (jt.contains("n_exposures") && !jt.at("n_exposures").is_null())
//       timing_.n_exposures_override = jt.at("n_exposures").get<int>();
//     timing_.exposure_time_s = jt.value("exposure_time_s", 0.0);
//   }
// }

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
}


// ---- NEW: helper to expand QE-Dark-style axis spec into a list of doubles ----
//
// Accepts any of:
//   - plain array: [0.5, 1.0, 2.0]
//   - {"values": [...]}
//   - {"linspace": { "start": ..., "stop": ..., "num": N, "endpoint": true/false }}
//   - {"logspace": { "start_exp": ..., "stop_exp": ..., "num": N, "endpoint": true/false }}
//
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
