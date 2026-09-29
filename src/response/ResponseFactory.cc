// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ResponseFactory.cc -- Builds a ResponseFold from config, replicating the
//  reference app's own
//  ChargeIonization/ChargeTransport/PatternClassifier/EfficiencyMC and
//  ClusterFitMC construction sequences.
// ===========================================================================

#include "ccdarksens/response/ResponseFactory.hh"

#include <TH1D.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/ClusterEnergyResponseFold.hh"
#include "ccdarksens/response/ClusterFitEngine.hh"
#include "ccdarksens/response/ClusterFitMC.hh"
#include "ccdarksens/response/EfficiencyMC.hh"
#include "ccdarksens/response/NeSpaceResponseFold.hh"
#include "ccdarksens/response/NoiseTailCalibrator.hh"
#include "ccdarksens/response/PatternClassifier.hh"
#include "ccdarksens/response/PatternResponseFold.hh"

namespace ccdarksens {

namespace {

// Resolve a config-relative path the same way the reference app does:
// non-absolute paths are taken relative to (config file's directory) / "..",
// i.e. paths like "data/..." resolve from the repo root when the config
// lives in configs/.
std::string ResolveConfigRelativePath(const std::string& path, const std::string& config_path) {
  if (path.empty()) return path;
  std::filesystem::path p(path);
  if (p.is_absolute()) return p.string();
  std::filesystem::path config_dir = std::filesystem::absolute(std::filesystem::path(config_path)).parent_path();
  return (config_dir / ".." / path).lexically_normal().string();
}

// Parse a (pattern_code, ne, efficiency) CSV, skipping comments/header rows,
// same tolerant row format as the reference app's inline parser.
void LoadPatternEffCsv(const std::string& resolved_path,
                        std::map<std::pair<int, int>, double>* out,
                        std::size_t* n_rows_loaded) {
  std::ifstream in(resolved_path);
  if (!in.is_open()) {
    throw std::runtime_error("ResponseFactory: failed to open pattern efficiency CSV: " + resolved_path);
  }
  std::string line;
  std::size_t n = 0;
  while (std::getline(in, line)) {
    while (!line.empty() && std::isspace(static_cast<unsigned char>(line.front()))) line.erase(line.begin());
    if (line.empty() || line[0] == '#') continue;
    if (line.rfind("pattern", 0) == 0) continue;  // header row

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

    (*out)[{pattern_code, ne_val}] = eff_val;
    ++n;
  }
  if (n_rows_loaded) *n_rows_loaded = n;
}

// n log-spaced values from lo to hi inclusive (empty if n <= 0 or a bound is not positive).
std::vector<double> LogSpace(double lo, double hi, int n) {
  std::vector<double> out;
  if (n <= 0 || lo <= 0.0 || hi <= 0.0) return out;
  out.reserve(static_cast<std::size_t>(n));
  const double log_lo = std::log10(lo);
  const double log_hi = std::log10(hi);
  for (int i = 0; i < n; ++i) {
    const double t = (n == 1) ? 0.0 : static_cast<double>(i) / static_cast<double>(n - 1);
    out.push_back(std::pow(10.0, log_lo + t * (log_hi - log_lo)));
  }
  return out;
}

}  // namespace

// ----------------------------------------------------------------------------
// MakeClusterEnergyFoldFromJSON
//   Build the cluster_energy fold for one ClusterFitMCJSON block:
//     1) calibrate the Delta LL noise-tail cut with pure-noise toys (using
//        lambda_dc_calib_per_pixel if it is set, else lambda_dc_per_pixel);
//     2) build the charge-transport model;
//     3) build the kernel K[E_true][E_reco] with ClusterFitMC (always with the
//        physical lambda_dc_per_pixel), on log-spaced true energies and linear
//        E_reco bins;
//     4) wrap it in ClusterEnergyResponseFold.
//   This is the expensive step (minutes), done once before the grid loop.
// ----------------------------------------------------------------------------
ResponseFactoryResult MakeClusterEnergyFoldFromJSON(const ClusterFitMCJSON& jcfg) {
  NoiseTailCalibratorConfig calib_cfg;
  calib_cfg.pix_cfg.nx = jcfg.window_nx;
  calib_cfg.pix_cfg.ny = jcfg.window_ny;
  calib_cfg.pix_cfg.pixel_size_um = jcfg.pixel_size_um;
  calib_cfg.pix_cfg.sigma_readout_e = jcfg.sigma_readout_e;
  calib_cfg.pix_cfg.lambda_dc = (jcfg.lambda_dc_calib_per_pixel >= 0.0)
                                    ? jcfg.lambda_dc_calib_per_pixel
                                    : jcfg.lambda_dc_per_pixel;
  calib_cfg.pix_cfg.rng_seed = jcfg.calib_rng_seed;
  calib_cfg.pix_cfg.collapse_y = jcfg.one_dimensional;
  calib_cfg.fit_cfg.sigma_pix_e = jcfg.sigma_readout_e;
  calib_cfg.fit_cfg.one_dimensional = jcfg.one_dimensional;
  calib_cfg.fit_cfg.method = (jcfg.fit_method == "minuit2")
                                 ? ClusterFitConfig::Method::kMinuit2
                                 : ClusterFitConfig::Method::kNelderMead;
  calib_cfg.fit_cfg.sigma_xy_lo_px = jcfg.sigma_xy_lo_px;
  calib_cfg.fit_cfg.sigma_xy_hi_px = jcfg.sigma_xy_hi_px;
  calib_cfg.n_toys = jcfg.n_toys;
  calib_cfg.target_tail_prob = jcfg.target_tail_prob;

  NoiseTailCalibrator calibrator(calib_cfg);
  const auto calib_result = calibrator.Run();

  ChargeTransportConfig ct_cfg;
  ct_cfg.thickness_um = jcfg.thickness_um;
  ct_cfg.A_um2 = jcfg.A_um2;
  ct_cfg.b_umInv = jcfg.b_umInv;
  ct_cfg.alpha = jcfg.alpha;
  ct_cfg.beta_per_keV = jcfg.beta_per_keV;
  ct_cfg.rng_seed = jcfg.rng_seed;
  auto ct = std::make_shared<ChargeTransport>(ct_cfg);

  ClusterFitMCConfig mc_cfg;
  mc_cfg.pix_cfg = calib_cfg.pix_cfg;
  mc_cfg.pix_cfg.lambda_dc = jcfg.lambda_dc_per_pixel;  // kernel trials always use the physical DC, whatever the calibration used
  mc_cfg.fit_cfg = calib_cfg.fit_cfg;
  mc_cfg.delta_ll_cut = calib_result.delta_ll_cut;
  mc_cfg.sigma_xy_fid_min_px = jcfg.sigma_xy_fid_min_px;
  mc_cfg.sigma_xy_fid_max_px = jcfg.sigma_xy_fid_max_px;
  mc_cfg.eh_pair_eV = jcfg.eh_pair_eV;
  mc_cfg.fano_factor = jcfg.fano_factor;
  mc_cfg.ne_trials_per_point = jcfg.ne_trials_per_point;
  mc_cfg.rng_seed = jcfg.rng_seed;

  ClusterFitMC cluster_mc(mc_cfg, ct);

  const auto Etrue_grid = LogSpace(jcfg.Etrue_min_eV, jcfg.Etrue_max_eV, jcfg.Etrue_npoints);
  std::vector<double> Ereco_edges;
  Ereco_edges.reserve(static_cast<std::size_t>(jcfg.Ereco_nbins) + 1);
  for (int i = 0; i <= jcfg.Ereco_nbins; ++i) {
    Ereco_edges.push_back(jcfg.Ereco_min_eV +
        i * (jcfg.Ereco_max_eV - jcfg.Ereco_min_eV) / jcfg.Ereco_nbins);
  }

  auto kernel = cluster_mc.BuildKernel(Etrue_grid, Ereco_edges);

  ResponseFactoryResult result;
  result.fold = std::make_unique<ClusterEnergyResponseFold>(std::move(kernel));
  result.ion = nullptr;
  result.ne_min_bkg = 0;
  result.ne_max = 0;
  return result;
}

namespace {

// ----------------------------------------------------------------------------
// MakePatternOrNeSpaceFold
//   Build the pattern-space or n_e-space fold:
//     1) ionization table and ChargeTransport from the config;
//     2) EfficiencyMC with the accepted labels for the chosen observable
//        (multi-pixel patterns from pattern_roi, or single-pixel ones from roi_bins);
//     3) the pattern classifier from response.pattern_classifier;
//     4) pattern_eff_map: from efficiency_csv if given, otherwise from the MC
//        pattern table, with efficiency_csv_reference overlaid on top;
//     5) n_e >= 10 is treated as fully efficient;
//     6) wrap everything in PatternResponseFold or NeSpaceResponseFold.
//   Throws std::runtime_error if no efficiencies are available or the pattern ROI is empty.
// ----------------------------------------------------------------------------
ResponseFactoryResult MakePatternOrNeSpaceFold(const ConfigManager& cfg,
                                                const ExperimentSummary& summary,
                                                const std::string& config_path) {
  const auto& det = cfg.detector();
  const auto& emj = cfg.response().emc;
  const bool use_pattern_bins = (summary.observable_bins == "pattern");

  const int ne_min = summary.binning.ne_min;
  const int ne_max = summary.binning.ne_max;
  const int ne_min_bkg = use_pattern_bins ? std::min(ne_min, 0) : ne_min;

  auto ion = std::make_shared<ChargeIonization>(cfg.response().charge_ionization.table_csv);

  ChargeTransportConfig ct_cfg;
  ct_cfg.thickness_um = det.geometry().thickness_mm * 1000.0;
  ct_cfg.A_um2        = emj.A_um2;
  ct_cfg.b_umInv      = emj.b_umInv;
  ct_cfg.alpha        = emj.alpha;
  ct_cfg.beta_per_keV = emj.beta_per_keV;
  ct_cfg.rng_seed     = emj.rng_seed;
  auto ct = std::make_shared<ChargeTransport>(ct_cfg);

  EfficiencyMCConfig emc_cfg;
  emc_cfg.ne_trials = emj.n_events_per_ne;
  emc_cfg.row_length = emj.row_length;
  emc_cfg.pix_cfg.mode            = PixelSimMode::RowSegment;
  emc_cfg.pix_cfg.nx              = emc_cfg.row_length;
  emc_cfg.pix_cfg.ny              = 1;
  emc_cfg.pix_cfg.pixel_size_um   = det.geometry().pixel_size_um;
  emc_cfg.pix_cfg.lambda_dc       = 0.0;
  emc_cfg.pix_cfg.sigma_readout_e = emj.sigma_readout_e;
  emc_cfg.pix_cfg.rng_seed        = emj.rng_seed;
  emc_cfg.seed                    = emj.rng_seed;

  // Accepted labels follow the analysis space, matching the reference app:
  // pattern mode uses the multi-pixel SRDM patterns in pattern_roi; n_e mode
  // uses the single-pixel patterns implied by roi_bins.
  if (use_pattern_bins) {
    for (int code : summary.pattern_roi) {
      PatternLabel lab;
      lab.isolated = true;
      lab.q = DecodePatternCode(code);
      if (!lab.q.empty()) emc_cfg.accepted_labels.push_back(lab);
    }
  } else {
    for (int ne : summary.roi_bins) {
      if (ne <= 0) continue;
      PatternLabel lab;
      lab.isolated = true;
      lab.q = {ne};
      emc_cfg.accepted_labels.push_back(lab);
    }
  }
  if (emc_cfg.accepted_labels.empty()) {
    PatternLabel lab;
    lab.isolated = true;
    lab.q = {1};
    emc_cfg.accepted_labels.push_back(lab);
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

  // pattern_eff_map: (pattern_code, ne) -> efficiency. Two sources, same as
  // the reference app: efficiency_csv (production path -- see
  // configs/scan_dmelectron_pattern_pydme_exact.json) replaces the table
  // entirely; otherwise it's filled from EfficiencyMC's own MC-generated
  // pattern table. efficiency_csv_reference always overlays on top, whichever
  // base was used.
  std::map<std::pair<int, int>, double> pattern_eff_map;
  if (!emj.efficiency_csv.empty()) {
    const std::string resolved = ResolveConfigRelativePath(emj.efficiency_csv, config_path);
    LoadPatternEffCsv(resolved, &pattern_eff_map, nullptr);
    if (pattern_eff_map.empty()) {
      throw std::runtime_error("ResponseFactory: efficiency_csv has no valid data: " + resolved);
    }
  }

  const double Ee_ref_eV = 50.0;  // reference energy for the efficiency MC [eV]
  std::unique_ptr<TH1D> h_eps_ne;
  if (!pattern_eff_map.empty()) {
    h_eps_ne = emc_ptr->PrecomputeEpsilonWithPatternEff(ne_min_bkg, ne_max, Ee_ref_eV, pattern_eff_map);
  } else {
    h_eps_ne = emc_ptr->PrecomputeEpsilon(ne_min_bkg, ne_max, Ee_ref_eV);

    const auto& table = emc_ptr->GetPatternTable();
    for (const auto& ne_entry : table) {
      const int ne = ne_entry.first;
      for (const auto& label_prob : ne_entry.second) {
        int code = 0;
        for (int d : label_prob.first.q) code = code * 10 + d;
        pattern_eff_map[{code, ne}] = label_prob.second;
      }
    }
    if (pattern_eff_map.empty()) {
      throw std::runtime_error(
          "ResponseFactory: no pattern efficiencies (no efficiency_csv and EfficiencyMC table empty).");
    }
  }

  if (!emj.efficiency_csv_reference.empty()) {
    const std::string resolved = ResolveConfigRelativePath(emj.efficiency_csv_reference, config_path);
    std::size_t n_ref = 0;
    LoadPatternEffCsv(resolved, &pattern_eff_map, &n_ref);
  }

  ResponseFactoryResult result;
  result.ion = ion;
  result.ne_min_bkg = ne_min_bkg;
  result.ne_max = ne_max;

  if (use_pattern_bins) {
    if (summary.pattern_roi.empty()) {
      throw std::runtime_error("ResponseFactory: observable_bins is 'pattern' but experiment.pattern_roi is empty.");
    }
    PatternResponseFoldConfig fold_cfg;
    fold_cfg.ne_min = ne_min;
    fold_cfg.ne_max = ne_max;
    fold_cfg.pattern_roi = summary.pattern_roi;
    fold_cfg.pattern_eff_map = std::move(pattern_eff_map);
    result.fold = std::make_unique<PatternResponseFold>(ion, std::move(fold_cfg));
  } else {
    NeSpaceResponseFoldConfig fold_cfg;
    fold_cfg.ne_min = ne_min;
    fold_cfg.ne_max = ne_max;
    fold_cfg.roi_bins = summary.roi_bins;
    for (int ne : summary.roi_bins) {
      const int bin = h_eps_ne->FindBin(static_cast<double>(ne));
      fold_cfg.eps_ne[ne] = h_eps_ne->GetBinContent(bin);
    }
    result.fold = std::make_unique<NeSpaceResponseFold>(ion, std::move(fold_cfg));
  }
  return result;
}

}  // namespace

// ----------------------------------------------------------------------------
// MakeResponseFold
//   Public entry point: pick the fold from response.analysis_space
//   ("cluster_energy", otherwise the pattern / n_e space set by
//   experiment.observable_bins). "pcd" is not supported and throws.
// ----------------------------------------------------------------------------
ResponseFactoryResult MakeResponseFold(const ConfigManager& cfg,
                                        const ExperimentSummary& summary,
                                        const std::string& config_path) {
  const std::string& analysis_space = cfg.response().analysis_space;
  if (analysis_space == "cluster_energy") {
    return MakeClusterEnergyFoldFromJSON(cfg.response().cluster_fit_mc);
  }
  if (analysis_space == "pcd") {
    throw std::runtime_error(
        "ResponseFactory: analysis_space 'pcd' has no working fold to mirror "
        "(the reference app's signal path bypasses its own PCD wiring too).");
  }
  return MakePatternOrNeSpaceFold(cfg, summary, config_path);
}

}  // namespace ccdarksens
