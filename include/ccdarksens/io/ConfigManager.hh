#pragma once
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include <optional>

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
  double exposure_kg_year = 0.0;
};

struct BackgroundJSON {
  // Dark current in e-/pixel/year
  double lambda_e_per_pix_per_year = 0.0;
  double norm_scale = 1.0;

  // pattern efficiency (flat for MVP)
  bool   has_flat_eps = false;
  double flat_eps = 1.0;

  // flat energy background (dR/dE) in events/(kg·year·keV)
  bool   has_flat_bkg = false;
  double flat_bkg_norm_per_kg_year = 0.0;  // events / (kg·year·keV)
  double flat_bkg_Emin_eV = 0.0;
  double flat_bkg_Emax_eV = 0.0;
  int    flat_bkg_nbins   = 0;
};

struct TimingJSON {
  std::optional<int> n_exposures_override;
  double exposure_time_s = 0.0;
};

struct ClusterMCJSON {
  int    n_events_per_ne = 20000;
  double sigma_readout_e = 0.16;
  double Qmin_e = 0.3;
  double Qmax_e = 4.0;
  double A_um2 = 803.25;
  double b_umInv = 6.5e-4;
  double alpha = 1.0;
  double beta_per_keV = 0.0;
  int rows_bin = 1;
  int cols_bin = 1;
  bool pileup_with_dc = false;
  uint64_t rng_seed = 987654321;
};

struct ResponseJSON {
  std::string mode = "fast"; // "fast" or "cluster_mc"
  ClusterMCJSON cmc;
};

struct ModelJSON {
  std::string type;              // e.g. "dm_electron"
  std::string material;          // "Si"
  std::string mediator;          // "heavy" | "massless"
  std::string rates_dir;         // e.g. data/qedark_rates/Si/heavy
  std::string filename_template; // e.g. dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv
  double      mchi_MeV = 0.0;    // used to resolve filename
  std::string sigma_e_cm2;       // keep string to match filename format exactly
  double      Emin_eV = 0.0;     // spectrum binning (internal representation)
  double      Emax_eV = 20.0;
  int         nbins   = 200;
};





class ConfigManager {
public:
  explicit ConfigManager(std::string path);
  void parse();

  const RunHeader&        run()           const noexcept { return run_; }
  const Detector&         detector()      const noexcept { return *detector_; }
  const ExperimentConfig& experiment_cfg() const noexcept { return exp_cfg_; }
  const BackgroundJSON& backgrounds() const noexcept { return bkg_; }
  const TimingJSON&     timing() const noexcept { return timing_; }
  const ResponseJSON& response() const noexcept { return response_; }
  const ModelJSON& model() const noexcept { return model_; }


private:
  std::string  path_;
  RunHeader    run_;
  std::unique_ptr<Detector> detector_;
  ExperimentConfig exp_cfg_;
  BackgroundJSON bkg_;
  TimingJSON     timing_;
  ResponseJSON   response_;
  ModelJSON      model_;
  

  void parse_run_(const nlohmann::json& j);
  void parse_detector_(const nlohmann::json& j);
  void parse_experiment_(const nlohmann::json& j);
  void parse_backgrounds_(const nlohmann::json& j);
  void parse_timing_(const nlohmann::json& j);
  void parse_response_(const nlohmann::json& j);
  void parse_model_(const nlohmann::json& j);

};

} // namespace ccdarksens
