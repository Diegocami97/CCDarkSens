// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ConfigManager.hh -- I declare the typed structs that hold every setting
//  read from the JSON config (run header, backgrounds, timing,
//  response/efficiency MC, cluster fit, model grid) and the ConfigManager
//  class that parses a config file into them and hands out read-only
//  accessors.
// ===========================================================================

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

// ----------------------------------------------------------------------------
// RunHeader
//   The "run" block: labels, output directory, confidence level, statistical
//   method options and the background-model choice. Most fields are documented
//   where they are declared.
// ----------------------------------------------------------------------------
struct RunHeader {
  std::string label;  // run label, used in output file names
  std::string outdir;  // output directory
  double      cl = 0.90;  // confidence level of the upper limit
  std::string test_stat;   // "PLR" | "CLs"
  int         n_toys = 0;  // number of toys (toy-based statistics only)
  uint64_t    rng_seed = 12345;  // seed for every random-number generator in the run
  int         verbosity = 1;  // 0 = quiet, higher = more logging
  double exposure_kg_year = 0.0;  // optional "exposure_kg_year" from the run block (default 0); the exposure the scans use comes from ExperimentSetup
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

  /// Band-tool extension (ccdarksens_band.cc's toy-MC threshold calibration),
  /// ported from ccdarksens_scan_srdm_pattern_csv.cc. "scan" (default) or
  /// "threshold_toys" -- the latter skips UL extraction and instead generates
  /// per-mass Poisson sub-toys at sigma_threshold(m_chi), writing a
  /// q_target_per_mass TGraph to <outdir>/qtarget_threshold.root.
  std::string mode = "scan";
  /// When non-empty, open this ROOT file's "q_target_per_mass" TGraph and use
  /// it as a per-mass PLR rejection threshold instead of the asymptotic
  /// constant (used by toy-MC band runs, produced by a prior threshold_toys run).
  std::string q_target_lookup_path;
  struct ThresholdToysJSON {
    std::string sigma_threshold_graph_path;
    long long n_threshold_toys = 10000;
    double percentile = 0.90;
    uint64_t rng_seed = 23456ULL;
  };
  ThresholdToysJSON threshold_toys;
};

// ----------------------------------------------------------------------------
// BackgroundJSON
//   The "backgrounds" block: dark current, the pattern efficiency, the flat
//   radiogenic spectrum, the background-efficiency CSV and the low-energy
//   excess (LEE) used by the WIMP-nucleon channel.
// ----------------------------------------------------------------------------
struct BackgroundJSON {
  // Dark current in e-/pixel/year
  double lambda_e_per_pix_per_year = 0.0;
  double norm_scale = 1.0;  // global multiplicative scale on the dark-current spectrum

  // pattern efficiency (flat for MVP)
  bool   has_flat_eps = false;  // true if a flat efficiency was given
  double flat_eps = 1.0;  // flat pattern efficiency

  // flat energy background (dR/dE) in events/(kg·year·keV)
  bool   has_flat_bkg = false;  // true if a flat energy background was given
  double flat_bkg_norm_per_kg_year = 0.0;  // events / (kg·year·keV)
  double flat_bkg_Emin_eV = 0.0;  // lower edge of the flat spectrum [eV]
  double flat_bkg_Emax_eV = 0.0;  // upper edge of the flat spectrum [eV]
  int    flat_bkg_nbins   = 0;  // number of bins of the flat spectrum

  /// Path to Background_efficiencies.csv (P(identified | true pattern)); empty = fold B_tot(n_e) with signal ε(pattern|n_e)
  std::string background_efficiency_csv;

  /// Low-energy excess (LEE): DAMIC SNOLAB's unexplained bulk ionization
  /// excess, dR/dE = rate*(1/eps)*exp(-E/eps). WIMP-nucleon (cluster_energy)
  /// channel only -- see docs/LowEnergyExcess_Design.md. Tier A: fixed
  /// shape, no free nuisance -- run with has_lee_bkg on/off for a
  /// conservative-vs-baseline bracket, not a profiled component.
  bool   has_lee_bkg = false;
  double lee_rate_per_kg_day = 0.0;   // events / (kg*day)
  double lee_decay_energy_eV = 0.0;   // exponential decay energy epsilon
};

// ----------------------------------------------------------------------------
// TimingJSON
//   The "backgrounds.timing" block: length of one exposure and an optional
//   explicit exposure count (otherwise derived from the livetime).
// ----------------------------------------------------------------------------
struct TimingJSON {
  std::optional<int> n_exposures_override;  // if set, use this many exposures instead of deriving them
  double exposure_time_s = 0.0;  // duration of one readout/exposure [s]
};

// ----------------------------------------------------------------------------
// ClusterMCJSON
//   LEGACY. The old ClusterMC backend no longer exists; I keep this struct
//   only so that old configs that still have a "cluster_mc" block are read and
//   copied into EfficiencyMCJSON instead of silently breaking.
// ----------------------------------------------------------------------------
struct ClusterMCJSON {
  int    n_events_per_ne = 20000;  // MC events per n_e value
  double sigma_readout_e = 0.16;  // readout noise per pixel [e-]
  double Qmin_e = 0.3;  // lower charge threshold [e-]
  double Qmax_e = 4.0;  // upper charge threshold [e-]
  double A_um2 = 803.25;  // diffusion constant A in sigma_xy^2 = -A ln(1 - b z) [um^2]
  double b_umInv = 6.5e-4;  // diffusion constant b [1/um]
  double alpha = 1.0;  // energy-dependence offset of the diffusion width
  double beta_per_keV = 0.0;  // energy-dependence slope of the diffusion width [1/keV]
  int rows_bin = 1;  // row binning factor
  int cols_bin = 1;  // column binning factor
  bool pileup_with_dc = false;  // include dark-current pile-up in the simulated events
  uint64_t rng_seed = 987654321;  // RNG seed
};

// ----------------------------------------------------------------------------
// EfficiencyMCJSON
//   The "response.efficiency_mc" block: settings of EfficiencyMC, which
//   computes P(pattern | n_e) (the detection efficiency per n_e) by Monte Carlo.
//   It is dark-current free by design; DC enters only as a background.
// ----------------------------------------------------------------------------
struct EfficiencyMCJSON {
  int    n_events_per_ne = 20000;  // MC events per n_e value
  double sigma_readout_e = 0.16;  // readout noise per pixel [e-]
  double Qmin_e = 0.3;  // lower charge threshold [e-]
  double Qmax_e = 4.0;  // upper charge threshold [e-]

  // diffusion
  double A_um2 = 803.25;  // diffusion constant A [um^2]
  double b_umInv = 6.5e-4;  // diffusion constant b [1/um]
  double alpha = 1.0;  // energy-dependence offset of the diffusion width
  double beta_per_keV = 0.0;  // energy-dependence slope of the diffusion width [1/keV]

  bool   enable_MN  = true;  // allow the two-pixel (M+N) pattern class
  bool   enable_MNL = true;  // allow the three-pixel (M+N+L) pattern class
  uint64_t rng_seed = 987654321ULL;  // RNG seed

  // --- NEW ---
  int row_length = 50;   ///< length of 1D row segment for EfficiencyMC
  bool use_2d_image_efficiency = false;  ///< if true, build P(pattern|n_e) from 2D image + isolation (notebook-style)

  // Pattern list: list of { q: [1,2,1], isolated: true }
  std::vector<std::vector<int>> accepted_pattern_q;  // flat q-lists
  std::vector<bool> accepted_pattern_isolated;       // parallel vector

  /// Optional (pattern, ne) -> efficiency CSV, overriding the MC-generated
  /// pattern table. Path resolved relative to the config file's directory.
  /// When set, this is the production path used to reproduce pydme/DAMIC-M
  /// results (see e.g. configs/scan_dmelectron_pattern_pydme_exact.json).
  std::string efficiency_csv;
  /// Optional overlay CSV applied on top of efficiency_csv (or the MC table
  /// when efficiency_csv is empty), same (pattern, ne, efficiency) format.
  std::string efficiency_csv_reference;
};

// ----------------------------------------------------------------------------
// PatternClassifierJSON
//   The "response.pattern_classifier" block: charge thresholds that decide
//   which pixel-charge configuration counts as which pattern label.
// ----------------------------------------------------------------------------
struct PatternClassifierJSON {
    double Qmin_e           = 0.3;  // minimum pixel charge to count as occupied [e-]
    double neighbor_Qmax_e  = 0.3;  ///< max charge for "empty" neighbor (notebook: use qmin for left/right isolation)
    double Qmax_e           = 4.0;  // maximum pixel charge of a single-pixel pattern [e-]
    double sigma_res_e      = 0.16;  // readout noise used to scale the thresholds [e-]
    bool   enable_MN        = true;  // enable the two-pixel class
    bool   enable_MNL       = true;  // enable the three-pixel class
    int    max_e_per_pixel  = 5;  // largest charge per pixel handled by the single-pixel branch [e-]
    double thr_M            = 4.0;  // threshold for the M pixel, in units of sigma_res_e
    double thr_MN           = 4.0;  // threshold for the M+N pair, in units of sigma_res_e
    double thr_MNL          = 5.5;  // threshold for the M+N+L triple, in units of sigma_res_e
    bool   allow_pattern_zero = false;  ///< allow pattern (0) in single-pixel branch (background efficiency)
    bool   single_pixel_use_round = false;  ///< single-pixel = round(q) clamped 1..5 (match reference CSV)
};

// ----------------------------------------------------------------------------
// PCDJSON
//   The "response.pcd" block: pixel-charge-distribution response (a
//   reserved alternative analysis space).
// ----------------------------------------------------------------------------
struct PCDJSON {
  double q_min       = 0.0;  // lowest charge bin edge [e-]
  double q_max       = 20.0;  // highest charge bin edge [e-]
  int    nbins       = 200;  // number of charge bins
  int    mc_trials   = 50000;  // Monte Carlo trials per n_e
  double sigma_res_e = 0.21;  ///< readout sigma for PCD ne-kernel (e-)
  double Dqmin       = 0.5;   ///< lower charge window for PCD ne-kernel (e-)
  double Dqmax       = 0.5;   ///< upper charge window for PCD ne-kernel (e-)
};

/// 2D binned image geometry for pattern efficiency (notebook-style).
/// Used when generating 2D images or when binning is required for efficiency computation.
// ----------------------------------------------------------------------------
// ChargeIonizationJSON
//   The "response.charge_ionization" block: which P(n_e | E) table to use and
//   the band gap / electron-hole pair energy it was generated with.
// ----------------------------------------------------------------------------
struct ChargeIonizationJSON {
  std::string table_csv = "data/p100K_table.csv";  // P(n_e | E) table
  double band_gap_eV = 1.2;   ///< pheno metadata (table already encodes physics)
  double eh_pair_eV = 3.8;  // mean energy per electron-hole pair [eV]
  std::string scenario;       ///< e.g. "D-equal", "B-thresh", "ref"
};

// ----------------------------------------------------------------------------
// PatternImageJSON
//   The "response.pattern_image" block: geometry of the 2D binned image used by
//   the diagnostic 2D pattern-efficiency path (dark current can be switched on
//   here to study the efficiency degradation).
// ----------------------------------------------------------------------------
struct PatternImageJSON {
  int nrows_binned = 3;    ///< number of rows after row binning (e.g. 3)
  int ncols        = 50;   ///< number of columns after column binning (e.g. 50)
  int row_binning  = 100;  ///< row binning factor → ny_raw = nrows_binned * row_binning
  int col_binning  = 1;    ///< column binning factor → nx_raw = ncols * col_binning (1 = no col bin)
  double pixel_size_um   = 15.0;  // pixel pitch [um]
  double sigma_readout_e = 0.21;  // readout noise per pixel [e-]
  double lambda_dc      = 0.0;   ///< dark current (0 = off)
  uint64_t rng_seed     = 987654321ULL;  // RNG seed
};

/// Config for the WIMP-nucleus SI channel's Phase 3/4 reconstruction
/// (NoiseTailCalibrator + ClusterFitMC). See docs/ClusterFitMC_Design.md.
// ----------------------------------------------------------------------------
// ClusterFitMCJSON
//   The "response.cluster_fit_mc" block for the WIMP-nucleon channel: pixel
//   window, noise, diffusion, fit-engine, noise-tail calibration and kernel
//   settings for ClusterFitMC (see docs/ClusterFitMC_Design.md).
// ----------------------------------------------------------------------------
struct ClusterFitMCJSON {
  // Pixel window geometry
  int    window_nx = 15;  // fit-window width [pixels]
  int    window_ny = 15;  // fit-window height [pixels]
  double pixel_size_um = 15.0;  // pixel pitch [um]
  double sigma_readout_e = 0.16;  // readout noise per pixel [e-]

  /// Diagnostic only: mean Poisson dark-current e-/pixel/image added to every
  /// pixel in both the noise-tail calibration toys and the kernel-building
  /// trials. Default 0 = off (production; DC is not part of this channel's
  /// model).
  double lambda_dc_per_pixel = 0.0;

  /// Diagnostic only: lambda used for the noise-tail calibration toys alone.
  /// Negative (default) = same as lambda_dc_per_pixel (self-consistent cut).
  /// Set to 0 with lambda_dc_per_pixel > 0 to mimic calibrating the cut from
  /// DC-free blank images and applying it to data that contains DC.
  double lambda_dc_calib_per_pixel = -1.0;

  /// When true, this channel is the DAMIC 1x100 readout mode: a 1D Gaussian
  /// fit over a single row (window_ny should be 1) whose truth is the
  /// y-marginal of the 2D charge cloud, with every electron's charge summed
  /// into that row regardless of y (PixelSimulatorConfig::collapse_y) --
  /// not a literal 2D fit on a ny=1 window (see ClusterFitEngine.hh /
  /// PixelSimulator.hh for why that would be wrong). Default false = the
  /// normal 1x1, 2D fit.
  bool one_dimensional = false;

  // Diffusion (ChargeTransport) -- mirrors EfficiencyMCJSON's diffusion sub-block
  double A_um2 = 803.25;  // diffusion constant A [um^2]
  double b_umInv = 6.5e-4;  // diffusion constant b [1/um]
  double alpha = 1.0;  // energy-dependence offset of the diffusion width
  double beta_per_keV = 0.0;  // energy-dependence slope of the diffusion width [1/keV]
  double thickness_um = 670.0;  // CCD thickness [um]

  // Fit engine (ClusterFitEngine)
  std::string fit_method = "nelder_mead";  ///< "nelder_mead" | "minuit2"
  double sigma_xy_lo_px = 0.1;  // lower bound of the fitted Gaussian width [pixels]
  double sigma_xy_hi_px = 2.0;  // upper bound of the fitted Gaussian width [pixels]

  // Phase 3: noise-tail calibration (NoiseTailCalibrator)
  int      n_toys = 100000;  // number of pure-noise toys for the noise-tail calibration
  double   target_tail_prob = 1e-3;  // fraction of pure-noise toys allowed to pass the Delta-LL cut
  uint64_t calib_rng_seed = 987654321ULL;  // seed of the calibration toys

  // Phase 4: kernel building (ClusterFitMC). eh_pair_eV/fano_factor are
  // PhysRevD.94.082006's own stated values -- see docs/ClusterFitMC_Design.md
  // §4.2 for the citation; fano_factor is explicitly a config knob (not a
  // hardcoded constant) because the paper itself treats it as uncertain for
  // nuclear recoils and sweeps it as a systematic.
  int    ne_trials_per_point = 5000;  // forward-simulated events per true-energy point
  double sigma_xy_fid_min_px = 0.35;  // lower edge of the fiducial (surface-rejection) cut on the fitted width [pixels]
  double sigma_xy_fid_max_px = 1.22;  // upper edge of the fiducial cut on the fitted width [pixels]
  double eh_pair_eV = 3.77;  // energy per electron-hole pair [eV]
  double fano_factor = 0.133;  // Fano factor of the charge generation
  uint64_t rng_seed = 987654321ULL;  // seed of the kernel-building events

  // E_true / E_reco grid (log-spaced true-energy points; linear E_reco bin edges)
  double Etrue_min_eV = 10.0;  // lowest true-energy grid point [eV]
  double Etrue_max_eV = 300.0;  // highest true-energy grid point [eV]
  int    Etrue_npoints = 12;  // number of log-spaced true-energy points
  double Ereco_min_eV = 0.0;  // lowest reconstructed-energy bin edge [eV]
  double Ereco_max_eV = 400.0;  // highest reconstructed-energy bin edge [eV]
  int    Ereco_nbins = 40;  // number of reconstructed-energy bins
};

/// One channel of a joint-likelihood response.channels[] array (see
/// ResponseJSON::channels below) -- e.g. DAMIC's 1x1 and 1x100 readout modes.
/// Each channel gets its own ResponseFold (built from cluster_fit_mc) and its
/// own flat-background rate; ccdarksens_scan_generic.cc folds the same
/// physical signal spectrum through each channel independently and sums
/// their profile-likelihood NLLs (see docs/ClusterFitMC_Design.md).
struct ResponseChannelJSON {
  std::string label;  // channel name, e.g. "1x1" or "1x100"
  ClusterFitMCJSON cluster_fit_mc;  // this channel's own reconstruction settings
  /// This channel's own flat Compton background rate, events/(kg*year*keV)
  /// -- replaces backgrounds.flat_background.norm_per_kg_year_keV for this
  /// channel (each channel has its own measured background rate).
  double flat_background_norm_per_kg_year_keV = 0.0;
  /// Optional path to a (E_keV, efficiency) CSV giving this channel's own
  /// background detection-efficiency curve, distinct from the signal's
  /// kernel-based efficiency (e.g. PhysRevD.94.082006 Fig. 9's dashed
  /// lines -- see docs/ClusterFitMC_Design.md Sec. 6.9). When set, the
  /// background is computed via MakeClusterEnergyFlatBackgroundWithEfficiency
  /// instead of folding through the signal kernel. Empty (default) keeps
  /// the previous behavior.
  std::string background_efficiency_csv;
};

// ----------------------------------------------------------------------------
// ResponseJSON
//   The "response" block: which analysis space and detector-response backend to
//   use, plus the settings of each backend (efficiency MC, pattern classifier,
//   pixel-charge distribution, cluster fit) and the optional joint-likelihood
//   channels.
// ----------------------------------------------------------------------------
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
  ClusterMCJSON cmc;  // legacy block; copied into emc only when no efficiency_mc block exists
  EfficiencyMCJSON emc;  // EfficiencyMC settings (the .emc accessor)
  // NEW
  PatternClassifierJSON pattern_classifier;  // pattern label thresholds
  PCDJSON       pcd;  // pixel-charge-distribution settings
  /// 2D image size and binning for pattern efficiency (optional)
  PatternImageJSON pattern_image;
  ChargeIonizationJSON charge_ionization;  // P(n_e | E) table settings
  ClusterFitMCJSON cluster_fit_mc;  // WIMP-nucleon cluster reconstruction settings
  /// Optional path to a (E_keV, efficiency) CSV giving the background's own
  /// detection-efficiency curve for the single-channel cluster_energy path
  /// (mirrors ResponseChannelJSON::background_efficiency_csv for the
  /// multi-channel case below). Empty (default) keeps the previous
  /// behavior (background folded through the signal kernel).
  std::string background_efficiency_csv;
  /// Joint-likelihood channels (e.g. DAMIC's 1x1 + 1x100 readout modes),
  /// cluster_energy analysis_space only. Empty (the default) means the
  /// single-channel path above is used exactly as before -- adding this
  /// field changes no existing config's behavior. See
  /// docs/ClusterFitMC_Design.md and ccdarksens_scan_generic.cc's joint
  /// branch (gated on !channels.empty()).
  std::vector<ResponseChannelJSON> channels;
};

// ---- NEW: grid axis helper for model.grid ----
// ----------------------------------------------------------------------------
// GridAxisJSON
//   One axis of the (mass, coupling) scan grid, expanded to explicit numbers
//   plus an optional printf-like format used to build the rate-file names.
// ----------------------------------------------------------------------------
struct GridAxisJSON {
  std::vector<double> values;   // expanded numeric grid
  std::string format;           // optional: printf-like, e.g. ".6f", ".1e"
};

// ----------------------------------------------------------------------------
// ModelJSON
//   The "model" block: which signal model to use, where its rate tables live,
//   the output spectrum binning and the (mass, coupling) scan grid.
// ----------------------------------------------------------------------------
struct ModelJSON {
  std::string type;              // "dm_electron" (default) | "dark_photon" | "migdal"
  std::string material;          // "Si" — electronic target (dm_electron, dark_photon)
  std::string mediator;          // "heavy" | "massless" — meaning depends on model.type:
                                  //   dm_electron/dark_photon: DM-electron mediator
                                  //   migdal: DM-nucleon mediator (separate namespace, same strings)
  std::string rates_dir;         // e.g. data/qedark_rates/Si/heavy
  std::string filename_template; // e.g. dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv
  double      mchi_MeV = 0.0;    // used to resolve filename (single-point apps)
  std::string sigma_e_cm2;       // keep string to match filename format exactly
  double      Emin_eV = 0.0;     // spectrum binning (internal representation)
  double      Emax_eV = 20.0;
  int         nbins   = 200;

  // ---- DM-nucleon coupling fields (model.type == "migdal"; reused by a future
  //      elastic nuclear-recoil/WIMP model — see DMNucleonConfig.hh) ----
  std::string target_nucleus = "Si28"; // nuclear target, distinct from electronic `material`
  int         nuclear_A = 28;          // mass number
  int         nuclear_Z = 14;          // atomic number
  std::string sigma_n_cm2;             // DM-nucleon cross section (kept as string for filename match)
  std::string epsilon_ref = "";        // dark_photon: if set, load rate at epsilon_ref, scale by (ε/ε_ref)²

  // ---- NEW: QE-Dark style grid specification for scans ----
  bool        has_grid = false;
  GridAxisJSON grid_mchi;        // values in MeV
  GridAxisJSON grid_sigma;       // values in cm^2 (numeric) — sigma_e_cm2, epsilon, or sigma_n_cm2 depending on type
};


// ----------------------------------------------------------------------------
// ConfigManager
//   I read one JSON config file and expose it as typed, read-only structs.
//   Usage:  ConfigManager cfg(path);  cfg.parse();  cfg.run(); cfg.model(); ...
//   parse() throws if the file cannot be opened or a required block is missing.
// ----------------------------------------------------------------------------
class ConfigManager {
public:
  // Remember the config path; nothing is read until parse() is called.
  explicit ConfigManager(std::string path);
  // Read the JSON file and fill every block (run, detector, experiment, backgrounds, timing, response, model).
  void parse();

  const RunHeader&        run()           const noexcept { return run_; }
  const Detector&         detector()      const noexcept { return *detector_; }
  const ExperimentConfig& experiment_cfg() const noexcept { return exp_cfg_; }
  const BackgroundJSON&   backgrounds()   const noexcept { return bkg_; }
  const TimingJSON&       timing()        const noexcept { return timing_; }
  const ResponseJSON&     response()      const noexcept { return response_; }
  const ModelJSON&        model()         const noexcept { return model_; }


private:
  std::string  path_;  // config file path
  RunHeader    run_;  // run block
  std::unique_ptr<Detector> detector_;  // detector built from the geometry/material blocks
  ExperimentConfig exp_cfg_;  // experiment block
  BackgroundJSON bkg_;  // backgrounds block
  TimingJSON     timing_;  // timing block
  ResponseJSON   response_;  // response block
  ModelJSON      model_;  // model block
  

  // One parser per config block; each fills the matching member from the JSON object.
  void parse_run_(const nlohmann::json& j);
  void parse_detector_(const nlohmann::json& j);
  void parse_experiment_(const nlohmann::json& j);
  void parse_backgrounds_(const nlohmann::json& j);
  void parse_timing_(const nlohmann::json& j);
  void parse_response_(const nlohmann::json& j);
  void parse_model_(const nlohmann::json& j);

};

} // namespace ccdarksens
