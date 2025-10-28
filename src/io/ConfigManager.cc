#include "ccdarksens/io/ConfigManager.hh"
#include <fstream>
#include <memory>
#include <stdexcept>
#include <utility>

#include <nlohmann/json.hpp>

namespace ccdarksens {

using nlohmann::json;

ConfigManager::ConfigManager(std::string path) : path_(std::move(path)) {}

void ConfigManager::parse() {
  std::ifstream in(path_);
  if (!in) throw std::runtime_error("Cannot open config: " + path_);
  json j;  // or: nlohmann::json j;
  in >> j;

  parse_run_(j.at("run"));
  parse_detector_(j.at("detector"));
  parse_experiment_(j.at("experiment"));
}

void ConfigManager::parse_run_(const json& j) {
  run_.label     = j.at("label").get<std::string>();
  run_.outdir    = j.at("outdir").get<std::string>();
  run_.cl        = j.value("cl", 0.90);
  run_.test_stat = j.value("test_stat", "PLR");
  run_.n_toys    = j.value("n_toys", 0);
  run_.rng_seed  = j.value("rng_seed", 12345ULL);
  run_.verbosity = j.value("verbosity", 1);
}

void ConfigManager::parse_detector_(const json& j) {
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

void ConfigManager::parse_experiment_(const json& j) {
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

} // namespace ccdarksens
