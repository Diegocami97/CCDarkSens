// ============================================================================
//  CCDarkSens — ConfigManager
//  Header for JSON configuration parsing and typed accessors for all framework run/detector/experiment settings.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "ccdarksens/detector/Detector.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"

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
  /// If true, dump per-grid-point spectra (dRdE, S_true, and S_obs/S_pat) into the output ROOT file.
  bool dump_point_spectra_root = false;
  /// If true, use profile likelihood (B(theta)=scale*B_template, minimize over scale) for q; matches pydme.
  bool use_profile_likelihood = false;
  /// Optional path to CSV with observed counts per pattern bin (same order as pattern_roi). If empty and use_profile_likelihood, use Asimov (data = B).
  std::string data_path;
  /// If true, collapse pattern bins to one Poisson bin for NLL (pydme single_bin option).
  bool single_bin_likelihood = false;
  /// Optional Gaussian prior on background scale: -log prior = 0.5*((scale-mean)/sigma)^2. If both set, constrain is applied.
  std::optional<double> constrain_scale_prior_mean;
  std::optional<double> constrain_scale_prior_sigma;
  /// How to obtain B_pat: "dc_flat_migration" (default, DC+flat then migration or signal epsilon) or "bp_br_template" (pydme-style B = Bp + Br from config, no DC/flat).
  std::string background_source = "dc_flat_migration";
  /// Background model for profile: "scale" (B = scale*B_template) or "Bp_theta_Br" (B = Bp + theta*Br, pydme).
  std::string background_model = "scale";
  /// For Bp_theta_Br: Bp and Br per pattern (same order as pattern_roi). Required when background_source is "bp_br_template".
  std::vector<double> background_Bp;
  std::vector<double> background_Br;
  /// Pydme-style constrain: L = sum_i (-theta*Br_i + prior_strength*ln(theta*Br_i)). 0 = off, 98 = pydme default.
  double constrain_prior_strength = 0.0;
  /// When true, use Gamma-prior sign (+theta*Br - N*ln(theta*Br)) so prior pulls theta toward N/Br; when false, match pydme sign.
  bool constrain_use_gamma_sign = false;
  /// When true, use pydme tau-weighted constraint (tau=N/sum(Br), Nrc_weighted per bin); prior mode at theta=1.
  bool constrain_use_tau_weighted = false;
  /// When > 1, use pydme multi-bin constraint: prior_strength*n_bins*ln(theta*Br_i/n_bins). Set to len(gamma) in data (e.g. 4450 for Final_Combined_Image_Data.csv) to match pydme SRDM. 1 = single-bin form.
  int constrain_n_bins = 1;
  /// Bounds for theta in Bp_theta_Br mode (log10(xsec) or scale in scale mode use scale_lo/hi in code). Defaults match pydme.
  double theta_lo = 0.5;
  double theta_hi = 10.0;
  /// If > 0 and profile likelihood is used, write a profile-likelihood plot (q_mu vs log10(sigma_e)) at this m_chi [MeV] (closest grid point). 0 = no plot.
  double profile_likelihood_plot_mchi = 0.0;
  /// Minimizer for profiling over theta: "brent" (default, 1D Brent) or "minuit" (ROOT Minuit2). Same NLL; minuit matches pydme's minimizer.
  std::string profile_minimizer = "brent";
  /// When true, match pydme UL: accept 2D fit at boundary and use 2D nll_min for q_mu/UL (same denominator and bracket start as pydme).
  bool pydme_style_ul = false;
  /// When true, cap local maxima in the stored TH1D upper limit (legacy post-scan smoother). Default false: raw bisection UL; downstream should prefer upper_limit_sigma_e_mchi_graph.
  bool smooth_ul_envelope = false;
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

  /// Path to Background_efficiencies.csv (P(identified | true pattern)); empty = fold B_tot(n_e) with signal ε(pattern|n_e)
  std::string background_efficiency_csv;
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

struct EfficiencyMCJSON {
  int    n_events_per_ne = 20000;
  double sigma_readout_e = 0.16;
  double Qmin_e = 0.3;
  double Qmax_e = 4.0;

  // diffusion
  double A_um2 = 803.25;
  double b_umInv = 6.5e-4;
  double alpha = 1.0;
  double beta_per_keV = 0.0;

  // --- LEGACY (still load but unused by EfficiencyMC) ---
  int    half_window_pix = 1;
  int    rows_bin = 1;
  int    cols_bin = 1;
  bool   pileup_with_dc = false;

  bool   enable_MN  = true;
  bool   enable_MNL = true;
  uint64_t rng_seed = 987654321ULL;

  // --- NEW ---
  int row_length = 50;   ///< length of 1D row segment for EfficiencyMC
  bool use_2d_image_efficiency = false;  ///< if true, build P(pattern|n_e) from 2D image + isolation (notebook-style)

  // Pattern list: list of { q: [1,2,1], isolated: true }
  std::vector<std::vector<int>> accepted_pattern_q;  // flat q-lists
  std::vector<bool> accepted_pattern_isolated;       // parallel vector
};

struct PatternClassifierJSON {
    double Qmin_e           = 0.3;
    double neighbor_Qmax_e  = 0.3;  ///< max charge for "empty" neighbor (notebook: use qmin for left/right isolation)
    double Qmax_e           = 4.0;
    double sigma_res_e      = 0.16;
    bool   enable_MN        = true;
    bool   enable_MNL       = true;
    int    max_e_per_pixel  = 5;
    double thr_M            = 4.0;
    double thr_MN           = 4.0;
    double thr_MNL          = 5.5;
    bool   allow_pattern_zero = false;  ///< allow pattern (0) in single-pixel branch (background efficiency)
    bool   single_pixel_use_round = false;  ///< single-pixel = round(q) clamped 1..5 (match reference CSV)
};

struct PCDJSON {
  double q_min       = 0.0;
  double q_max       = 20.0;
  int    nbins       = 200;
  int    mc_trials   = 50000;
  double sigma_res_e = 0.21;  ///< readout sigma for PCD ne-kernel (e-)
  double Dqmin       = 0.5;   ///< lower charge window for PCD ne-kernel (e-)
  double Dqmax       = 0.5;   ///< upper charge window for PCD ne-kernel (e-)
};

/// 2D binned image geometry for pattern efficiency (notebook-style).
/// Used when generating 2D images or when binning is required for efficiency computation.
struct ChargeIonizationJSON {
  std::string table_csv = "data/p100K_table.csv";
  double band_gap_eV = 1.2;   ///< pheno metadata (table already encodes physics)
  double eh_pair_eV = 3.8;
  std::string scenario;       ///< e.g. "D-equal", "B-thresh", "ref"
};

struct PatternImageJSON {
  int nrows_binned = 3;    ///< number of rows after row binning (e.g. 3)
  int ncols        = 50;   ///< number of columns after column binning (e.g. 50)
  int row_binning  = 100;  ///< row binning factor → ny_raw = nrows_binned * row_binning
  int col_binning  = 1;    ///< column binning factor → nx_raw = ncols * col_binning (1 = no col bin)
  double pixel_size_um   = 15.0;
  double sigma_readout_e = 0.21;
  double lambda_dc      = 0.0;   ///< dark current (0 = off)
  uint64_t rng_seed     = 987654321ULL;
};

struct ResponseJSON {
  // "fast"     = analytic / table ε(n_e)
  // "cluster_mc" = legacy ClusterMC backend
  // "pattern"  = EfficiencyMC backend
  // "pcd"      = reserved for future PCD backend
  std::string mode = "fast"; // "fast" or "cluster_mc"
  // analysis_space selects which space the likelihood is defined in:
  //   "pattern" (default) → n_e / pattern-efficiency space
  //   "pcd"               → pixel-charge distribution space
  //
  // By default this is derived from `mode`:
  //   mode == "pcd"  → analysis_space = "pcd"
  //   otherwise      → analysis_space = "pattern"
  std::string analysis_space = "pattern"; // "pattern" or "pixel"
  ClusterMCJSON cmc;
  EfficiencyMCJSON emc;  
  // NEW
  PatternClassifierJSON pattern_classifier;
  PCDJSON       pcd;
  /// 2D image size and binning for pattern efficiency (optional)
  PatternImageJSON pattern_image;
  ChargeIonizationJSON charge_ionization;
};

// ---- NEW: grid axis helper for model.grid ----
struct GridAxisJSON {
  std::vector<double> values;   // expanded numeric grid
  std::string format;           // optional: printf-like, e.g. ".6f", ".1e"
};

struct ModelJSON {
  std::string type;              // e.g. "dm_electron"
  std::string material;          // "Si"
  std::string mediator;          // "heavy" | "massless"
  std::string rates_dir;         // e.g. data/qedark_rates/Si/heavy
  std::string filename_template; // e.g. dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv
  double      mchi_MeV = 0.0;    // used to resolve filename (single-point apps)
  std::string sigma_e_cm2;       // keep string to match filename format exactly
  double      Emin_eV = 0.0;     // spectrum binning (internal representation)
  double      Emax_eV = 20.0;
  int         nbins   = 200;

  // ---- NEW: QE-Dark style grid specification for scans ----
  bool        has_grid = false;
  GridAxisJSON grid_mchi;        // values in MeV
  GridAxisJSON grid_sigma;       // values in cm^2 (numeric)
};


class ConfigManager {
public:
  explicit ConfigManager(std::string path);
  void parse();

  const RunHeader&        run()           const noexcept { return run_; }
  const Detector&         detector()      const noexcept { return *detector_; }
  const ExperimentConfig& experiment_cfg() const noexcept { return exp_cfg_; }
  const BackgroundJSON&   backgrounds()   const noexcept { return bkg_; }
  const TimingJSON&       timing()        const noexcept { return timing_; }
  const ResponseJSON&     response()      const noexcept { return response_; }
  const ModelJSON&        model()         const noexcept { return model_; }


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
