#pragma once
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json_fwd.hpp>

#include "ccdarksens/detector/Detector.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"

// namespace nlohmann { class json; }

namespace ccdarksens {

struct RunHeader {
  std::string label;
  std::string outdir;
  double      cl = 0.90;
  std::string test_stat;   // "PLR" | "CLs"
  int         n_toys = 0;
  uint64_t    rng_seed = 12345;
  int         verbosity = 1;
};

class ConfigManager {
public:
  explicit ConfigManager(std::string path);
  void parse();

  const RunHeader&        run()           const noexcept { return run_; }
  const Detector&         detector()      const noexcept { return *detector_; }
  const ExperimentConfig& experiment_cfg() const noexcept { return exp_cfg_; }

private:
  std::string  path_;
  RunHeader    run_;
  std::unique_ptr<Detector> detector_;
  ExperimentConfig exp_cfg_;

  void parse_run_(const nlohmann::json& j);
  void parse_detector_(const nlohmann::json& j);
  void parse_experiment_(const nlohmann::json& j);
};

} // namespace ccdarksens
