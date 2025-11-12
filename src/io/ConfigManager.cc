#include "ccdarksens/io/ConfigManager.hh"
#include <fstream>
#include <memory>
#include <stdexcept>
#include <utility>

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
  if (j.contains("detector"))     parse_detector_(j.at("detector"));
  if (j.contains("experiment"))   parse_experiment_(j.at("experiment"));
  if (j.contains("backgrounds")) parse_backgrounds_(j.at("backgrounds"));
  if (j.contains("response")) parse_response_(j.at("response"));
  if (j.contains("model")) parse_model_(j.at("model"));
}

void ConfigManager::parse_run_(const nlohmann::json& j) {
  // run_.label     = j.at("label").get<std::string>();
  // run_.outdir    = j.at("outdir").get<std::string>();
  run_.label     = j.value("label", std::string{});
  run_.outdir    = j.value("outdir", std::string{});
  run_.cl        = j.value("cl", 0.90);
  run_.test_stat = j.value("test_stat", "PLR");
  run_.n_toys    = j.value("n_toys", 0);
  run_.rng_seed  = j.value("rng_seed", 12345ULL);
  run_.verbosity = j.value("verbosity", 1);
  run_.exposure_kg_year = j.value("exposure_kg_year", 0.0);
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
}

void ConfigManager::parse_response_(const nlohmann::json& jr) {
  response_.mode = jr.value("mode", std::string("fast"));

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
}



void ConfigManager::parse_backgrounds_(const nlohmann::json& jb) {
  if (jb.contains("dark_current")) {
    const auto& jd = jb.at("dark_current");
    bkg_.lambda_e_per_pix_per_exposure = jd.value("lambda_e_per_pix_per_exposure", 0.0);
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

  // timing knobs live under experiment (but expose via timing_)
  // Find them from the already-parsed experiment block in the original nlohmann::json:
  // (we still have j in this scope? If not, read from jb's parent; simplest: require backgrounds to also include timing)
  // For clarity now, we’ll pull from a sibling "experiment" via the stored exp_cfg_ — but we need exposure_time_s & n_exposures.
  // Easiest: require they sit under "experiment" and re-open here:
  // (If you want, we can move this logic to parse_experiment_ later.)

  // Try parent: we can't access parent here; instead, read from a convenience duplication in backgrounds if present:
  if (jb.contains("timing")) {
    const auto& jt = jb.at("timing");
    if (jt.contains("n_exposures") && !jt.at("n_exposures").is_null())
      timing_.n_exposures_override = jt.at("n_exposures").get<int>();
    timing_.exposure_time_s = jt.value("exposure_time_s", 0.0);
  }
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
}


} // namespace ccdarksens
