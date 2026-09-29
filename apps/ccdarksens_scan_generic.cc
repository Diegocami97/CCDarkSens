// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_scan_generic.cc
//  Generic, config-driven multi-channel sensitivity scan over a (mchi,
//  coupling) grid, producing a sigma_UL(mchi) exclusion curve for ANY
//  model/analysis-space combination the framework supports (dm_electron,
//  dark_photon, migdal, wimp_nucleon) x (pattern, n_e, cluster_energy).
//
//  Built on top of the same already-validated primitives the reference app
//  (ccdarksens_scan_dmelectron_pattern.cc, NOT modified, read-only reference)
//  uses -- ModelFactory, ResponseFactory/ResponseFold, BackgroundFactory,
//  ProfileLikelihood, ScanUtils -- generalized over ResponseFold so the same
//  Phase A / Phase B / Phase C loop works for every analysis space.
//
//  Output object names match the reference app's own so ccdarksens_band and
//  ccdarksens_plot_dmelectron_limit (both unmodified) work unchanged:
//    q_mchi_sigma (TH2D), upper_limit_sigma_e_mchi (TH1D),
//    upper_limit_sigma_e_mchi_graph (TGraph, the one both consumers require),
//    exposure_kg_year (TParameter<double>).
//
//  Usage: ccdarksens_scan_generic <config.json>
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TParameter.h>
#include <cctype>
#include <map>
#include <random>

#include <nlohmann/json.hpp>

#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/model/ModelFactory.hh"
#include "ccdarksens/response/BackgroundFactory.hh"
#include "ccdarksens/response/ResponseFactory.hh"
#include "ccdarksens/scan/ScanUtils.hh"
#include "ccdarksens/stats/ProfileLikelihood.hh"
#include "ccdarksens/utils/AppUtils.hh"

using nlohmann::json;
using namespace ccdarksens;

namespace {

// ----------------------------------------------------------------------------
// FormatSigma
//   Coupling value -> the string used in the rate-file names. A format like ".3e"
//   gives 3 digits in scientific notation; anything else falls back to 6.
// ----------------------------------------------------------------------------
std::string FormatSigma(double sigma, const std::string& fmt) {
  int prec = 6;
  if (!fmt.empty() && fmt.front() == '.' && (fmt.back() == 'e' || fmt.back() == 'E')) {
    try {
      prec = std::stoi(fmt.substr(1, fmt.size() - 2));
    } catch (...) {
      prec = 6;
    }
  }
  std::ostringstream ss;
  ss.setf(std::ios::scientific);
  ss << std::setprecision(prec) << sigma;
  return ss.str();
}

// ----------------------------------------------------------------------------
// ExpandAxisFromGrid
//   Expand the first of the given keys that exists in model.grid into a list of values
//   (through utils::ExpandAxis). Throws if none of the keys is present.
// ----------------------------------------------------------------------------
std::vector<double> ExpandAxisFromGrid(const json& jgrid, std::initializer_list<std::string> keys) {
  for (const auto& k : keys) {
    if (jgrid.contains(k)) return ccdarksens::utils::ExpandAxis(jgrid.at(k), k);
  }
  throw std::runtime_error("model.grid: none of the expected keys found");
}

// ----------------------------------------------------------------------------
// MakeEdgesFromCenters
//   Bin edges for a histogram whose bin centres are c: midpoints between neighbors, with
//   the outer edges extended by half a step. (A single centre gets +/-50% around it.)
// ----------------------------------------------------------------------------
std::vector<double> MakeEdgesFromCenters(const std::vector<double>& c) {
  const std::size_t N = c.size();
  std::vector<double> edges(N + 1);
  if (N == 0) return edges;
  if (N == 1) {
    const double w = std::abs(c[0]) > 0 ? std::abs(c[0]) * 0.5 : 0.5;
    edges[0] = c[0] - w;
    edges[1] = c[0] + w;
    return edges;
  }
  edges[0] = c[0] - 0.5 * (c[1] - c[0]);
  for (std::size_t i = 1; i < N; ++i) edges[i] = 0.5 * (c[i - 1] + c[i]);
  edges[N] = c[N - 1] + 0.5 * (c[N - 1] - c[N - 2]);
  return edges;
}

// ----------------------------------------------------------------------------
// ObservedCountsFromConfig
//   Observed counts given inline as run.observed_counts, or an empty vector if the key is
//   missing. Throws if the length differs from the number of response bins.
// ----------------------------------------------------------------------------
std::vector<double> ObservedCountsFromConfig(const json& jroot, std::size_t expected_size) {
  if (!jroot.contains("run")) return {};
  const auto& jr = jroot["run"];
  if (!jr.contains("observed_counts") || jr["observed_counts"].is_null()) return {};
  auto counts = jr["observed_counts"].get<std::vector<double>>();
  if (counts.empty()) return {};
  if (counts.size() != expected_size) {
    throw std::runtime_error("run.observed_counts size (" + std::to_string(counts.size()) +
                              ") must match the response fold's bin count (" + std::to_string(expected_size) + ")");
  }
  return counts;
}

// Sum of all elements of v.
double SumVec(const std::vector<double>& v) {
  double s = 0.0;
  for (double x : v) s += x;
  return s;
}

// Ported near-verbatim from the reference app's own data-loading helpers
// (ccdarksens_scan_dmelectron_pattern.cc) so run.data_path works the same
// way here: real observed counts, from CSV (flat row) or ROOT (D_pat
// histogram, mapped by bin label to roi_ids when labels are present).
std::vector<double> LoadDataCsv(const std::string& path, std::size_t expected_size) {
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
      } catch (...) {
        break;
      }
    }
    if (out.size() >= expected_size) break;
  }
  return out;
}

// Integer pattern id from a histogram bin label, or -1 if the label is empty or not a number.
int ParseDataBinLabel(const char* label) {
  if (!label || !label[0]) return -1;
  try {
    return std::stoi(label);
  } catch (...) {
    return -1;
  }
}

// ----------------------------------------------------------------------------
// LoadDataRoot
//   Read the D_pat histogram from a ROOT file. If its bins are labeled with pattern ids, I
//   return the counts in the order of roi_ids (a missing id counts as 0); otherwise the bins
//   in file order. Returns an empty vector if the file or the histogram is missing.
// ----------------------------------------------------------------------------
std::vector<double> LoadDataRoot(const std::string& path, const std::vector<int>& roi_ids) {
  TFile f(path.c_str(), "READ");
  if (!f.IsOpen()) return {};
  TH1D* h = nullptr;
  f.GetObject("D_pat", h);
  if (!h) return {};
  const int nb = h->GetNbinsX();

  std::map<int, double> by_id;
  bool has_labels = false;
  for (int i = 1; i <= nb; ++i) {
    const char* lbl = h->GetXaxis()->GetBinLabel(i);
    const int id = ParseDataBinLabel(lbl);
    if (id >= 0) {
      has_labels = true;
      by_id[id] = h->GetBinContent(i);
    }
  }
  if (!roi_ids.empty() && has_labels) {
    std::vector<double> out;
    out.reserve(roi_ids.size());
    for (int pid : roi_ids) {
      auto it = by_id.find(pid);
      out.push_back(it != by_id.end() ? it->second : 0.0);
    }
    return out;
  }
  std::vector<double> out;
  out.reserve(static_cast<std::size_t>(nb));
  for (int i = 1; i <= nb; ++i) out.push_back(h->GetBinContent(i));
  return out;
}

// Load observed counts: a ROOT file (D_pat) goes to LoadDataRoot, anything else is read as CSV.
std::vector<double> LoadData(const std::string& path, std::size_t expected_size, const std::vector<int>& roi_ids) {
  if (path.size() >= 5 && path.compare(path.size() - 5, 5, ".root") != 0) return LoadDataCsv(path, expected_size);
  return LoadDataRoot(path, roi_ids);
}

// ---- Joint-likelihood scan (response.channels[], e.g. DAMIC's 1x1 + 1x100
// readout modes) -- a separate, self-contained path from the single-channel
// Phase A/B/C below, so the existing (parity-gated) single-channel path is
// never touched when response.channels is non-empty. Each channel is folded
// and profiled completely independently (its own ResponseFold, its own
// ProfileLikelihood/background); at each grid point the two channels' NLLs
// are simply summed before MonotonizeQ/UlFromQMonoCrossing -- this is exact,
// not an approximation, because each channel's background nuisance
// parameter only ever appears in that channel's own NLL term (see
// docs/ClusterFitMC_Design.md for the derivation). Scope for this pass:
// profile_minimizer brent/minuit only, background_model=scale only (each
// channel's own flat rate via response.channels[].flat_background_norm_per_kg_year_keV),
// Asimov data only (data_k = B_k) -- matches what the DAMIC 2016 joint
// config actually needs; see the plan's explicit scope cuts for why.
int RunJointChannelScan(ConfigManager& cfg, const json& jroot) {
  const auto& run = cfg.run();
  if (run.mode == "threshold_toys") {
    throw std::runtime_error("scan-generic: run.mode=threshold_toys is not supported for "
                              "joint-channel (response.channels) configs.");
  }
  if (run.background_model == "Bp_theta_Br" || run.background_model == "Bp_theta_br") {
    throw std::runtime_error("scan-generic: background_model=Bp_theta_Br is not supported for "
                              "joint-channel configs yet (scale only).");
  }
  if (run.profile_minimizer != "brent" && run.profile_minimizer != "minuit") {
    throw std::runtime_error("scan-generic: profile_minimizer '" + run.profile_minimizer +
                              "' is not supported for joint-channel configs yet (brent/minuit only).");
  }
  if (!run.data_path.empty() && run.verbosity >= 1) {
    std::cout << "[scan-generic][joint] NOTE: run.data_path is ignored for joint-channel configs "
                 "(Asimov data = background per channel only, this pass).\n";
  }

  ExperimentSetup setup(cfg.experiment_cfg(), cfg.detector().mass_kg(), run.rng_seed);
  auto summary = setup.prepare_summary();
  const bool single_bin = run.single_bin_likelihood;

  struct Channel {
    std::string label;
    ResponseFactoryResult response;
    ccdarksens::stats::ProfileLikelihood pl;
  };

  std::vector<Channel> channels;
  for (const auto& ch_json : cfg.response().channels) {
    Channel ch;
    ch.label = ch_json.label;
    ch.response = MakeClusterEnergyFoldFromJSON(ch_json.cluster_fit_mc);
    ResponseFold& fold = *ch.response.fold;
    BackgroundResult bg;
    if (!ch_json.background_efficiency_csv.empty()) {
      // Background's own detection-efficiency curve, distinct from the
      // signal's kernel (e.g. PhysRevD.94.082006 Fig. 9's dashed lines) --
      // see docs/ClusterFitMC_Design.md Sec. 6.9.
      const auto eff_table = LoadBackgroundEfficiencyTable(ch_json.background_efficiency_csv);
      bg = MakeClusterEnergyFlatBackgroundWithEfficiency(
          fold, summary.exposure_kg_year, ch_json.flat_background_norm_per_kg_year_keV, eff_table);
    } else {
      bg = MakeClusterEnergyFlatBackground(fold, summary.exposure_kg_year, cfg.model().Emin_eV,
                                            cfg.model().Emax_eV, cfg.model().nbins,
                                            ch_json.flat_background_norm_per_kg_year_keV);
    }
    std::vector<double> B_pl = single_bin ? std::vector<double>{SumVec(bg.B_pat)} : bg.B_pat;
    ch.pl.SetData(B_pl);  // Asimov: data = background
    ch.pl.SetBTemplate(B_pl);
    std::cout << "[scan-generic][joint] channel=" << ch.label << " bins=" << fold.NumBins()
              << " flat_bkg_rate=" << ch_json.flat_background_norm_per_kg_year_keV
              << " background_efficiency_csv=" << (ch_json.background_efficiency_csv.empty() ? "(none)" : ch_json.background_efficiency_csv)
              << " total_B=" << SumVec(bg.B_pat) << "\n";
    channels.push_back(std::move(ch));
  }
  if (channels.empty()) throw std::runtime_error("scan-generic: response.channels is empty.");

  const double profile_param_lo = 0.01, profile_param_hi = 10.0;
  const double target_q = std::pow(TMath::NormQuantile(run.cl), 2);  // q threshold for the one-sided limit at confidence level cl

  const auto& mj = cfg.model();
  const auto& jgrid = jroot.at("model").at("grid");
  std::vector<double> mchi_list, sigma_list;  // mass axis and coupling axis of the scan
  if (mj.type == "migdal" || mj.type == "wimp_nucleon") {
    mchi_list = ExpandAxisFromGrid(jgrid, {"mchi_MeV"});
    sigma_list = ExpandAxisFromGrid(jgrid, {"sigma_n_cm2", "sigma_e_cm2"});
  } else {
    mchi_list = ExpandAxisFromGrid(jgrid, {"mchi_MeV"});
    sigma_list = ExpandAxisFromGrid(jgrid, {"sigma_e_cm2"});
  }
  if (mchi_list.empty() || sigma_list.empty()) throw std::runtime_error("Empty mass or coupling grid.");
  const std::string fmt_sigma = jgrid.contains("format") ? jgrid["format"].value("sigma", std::string(".1e")) : ".1e";

  auto edges_m = MakeEdgesFromCenters(mchi_list);
  auto edges_s = MakeEdgesFromCenters(sigma_list);
  TH2D h_q("q_mchi_sigma", ";m_{#chi} [MeV];coupling;q_{Asimov}", static_cast<int>(mchi_list.size()), edges_m.data(),
           static_cast<int>(sigma_list.size()), edges_s.data());
  TH1D h_upper_limit("upper_limit_sigma_e_mchi", ";m_{#chi} [MeV];coupling (joint CL)",
                      static_cast<int>(mchi_list.size()), edges_m.data());

  int idx_mchi = 0;
  for (double mchi : mchi_list) {
    std::vector<double> nll_values;
    nll_values.reserve(sigma_list.size());
    for (double sigma_val : sigma_list) {
      bool sig_ok = false;
      auto dRdE_sig = MakeSignalSpectrumE(mj, mchi, FormatSigma(sigma_val, fmt_sigma), &sig_ok);
      if (!sig_ok && run.verbosity >= 1) {
        std::cout << "[scan-generic][joint] WARNING: failed to load rate table for mchi=" << mchi
                   << ", coupling=" << sigma_val << "; using all-zero signal.\n";
      }
      double nll_sum = 0.0;
      for (auto& ch : channels) {
        std::vector<double> S = ch.response.fold->Fold(*dRdE_sig, summary.exposure_kg_year);
        std::vector<double> S_pl = single_bin ? std::vector<double>{SumVec(S)} : S;
        const double nll = (run.profile_minimizer == "minuit")
            ? ch.pl.MinimizeOverScaleMinuit(S_pl, profile_param_lo, profile_param_hi).second
            : ch.pl.MinimizeOverScale(S_pl, profile_param_lo, profile_param_hi).second;
        nll_sum += nll;
      }
      nll_values.push_back(nll_sum);
    }
    const double nll_min = *std::min_element(nll_values.begin(), nll_values.end());
    const auto q_mono = ccdarksens::scan::MonotonizeQ(nll_values, nll_min);
    const double ul_sigma = ccdarksens::scan::UlFromQMonoCrossing(sigma_list, q_mono, target_q);
    for (int k = 0; k < static_cast<int>(nll_values.size()); ++k) {
      if (q_mono[static_cast<std::size_t>(k)] <= 0.0) continue;
      h_q.SetBinContent(idx_mchi + 1, k + 1, q_mono[static_cast<std::size_t>(k)]);
    }
    h_upper_limit.SetBinContent(idx_mchi + 1, ul_sigma);
    if (run.verbosity >= 1) {
      std::cout << "[scan-generic][joint] mchi=" << mchi << " -> upper_limit=" << ul_sigma << "\n";
    }
    ++idx_mchi;
  }

  if (run.smooth_ul_envelope) {
    constexpr int kEnvW = 5;
    const int n_mchi = h_upper_limit.GetNbinsX();
    std::vector<double> raw(static_cast<std::size_t>(n_mchi));
    for (int i = 0; i < n_mchi; ++i) raw[static_cast<std::size_t>(i)] = h_upper_limit.GetBinContent(i + 1);
    std::vector<double> env(static_cast<std::size_t>(n_mchi), 0.0);
    for (int i = 0; i < n_mchi; ++i) {
      const int lo = std::max(0, i - kEnvW);
      const int hi = std::min(n_mchi - 1, i + kEnvW);
      double win_min = raw[static_cast<std::size_t>(i)];
      for (int j = lo; j <= hi; ++j) {
        const double v = raw[static_cast<std::size_t>(j)];
        if (v > 0.0 && (win_min <= 0.0 || v < win_min)) win_min = v;
      }
      env[static_cast<std::size_t>(i)] = win_min;
    }
    for (int i = 0; i < n_mchi; ++i) h_upper_limit.SetBinContent(i + 1, env[static_cast<std::size_t>(i)]);
  }

  const int n_mchi = h_upper_limit.GetNbinsX();
  TGraph g_upper_limit(n_mchi);
  for (int i = 0; i < n_mchi; ++i) {
    g_upper_limit.SetPoint(i, mchi_list[static_cast<std::size_t>(i)], h_upper_limit.GetBinContent(i + 1));
  }
  g_upper_limit.SetName("upper_limit_sigma_e_mchi_graph");

  const std::string out_path = run.outdir + "/scan_generic.root";
  TFile fout(out_path.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    throw std::runtime_error("[scan-generic][joint] cannot create output file " + out_path);
  }
  h_q.Write("q_mchi_sigma");
  h_upper_limit.Write("upper_limit_sigma_e_mchi");
  g_upper_limit.Write("upper_limit_sigma_e_mchi_graph");
  TParameter<double> p_exp("exposure_kg_year", summary.exposure_kg_year);
  p_exp.Write();
  fout.Close();

  std::cout << "[scan-generic][joint] Done. Output written to " << out_path << "\n";
  return 0;
}

}  // namespace

// ----------------------------------------------------------------------------
// main
//   Generic scan over the (mass, coupling) grid of the config. Usage: <program> config.json.
//     Setup      - parse the config, then dispatch to RunJointChannelScan when
//                  response.channels is set;
//     Phase A    - build once: experiment summary, response fold and background, the
//                  observed data (run.observed_counts > run.data_path > Asimov), and the
//                  ProfileLikelihood in scale or Bp + theta*Br mode;
//     Grid       - mass and coupling axes for the model type, plus an optional per-mass
//                  q threshold (run.q_target_lookup_path);
//     Threshold  - run.mode = "threshold_toys": toy-MC calibration of that threshold, then exit;
//     Phase B    - for every mass, profile the NLL along the coupling axis (optionally
//                  with the 2D pydme fit and continuous bisection), monotonize q and
//                  read off the upper limit;
//     Phase C    - optional envelope smoothing and writing <outdir>/scan_generic.root
//                  (q_mchi_sigma, upper_limit_sigma_e_mchi[_graph], exposure_kg_year).
//   Returns 0 on success, 1 on any error.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json\n";
    return 1;
  }
  const std::string config_path = argv[1];

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
    if (!run.outdir.empty()) std::filesystem::create_directories(run.outdir);

    if (!cfg.response().channels.empty()) {
      return RunJointChannelScan(cfg, jroot);
    }

    // ---- Phase A: build once, before the grid loop ----
    ExperimentSetup setup(cfg.experiment_cfg(), cfg.detector().mass_kg(), run.rng_seed);
    auto summary = setup.prepare_summary();
    const bool use_pattern_bins = (summary.observable_bins == "pattern");

    std::cout << "[scan-generic] analysis_space=" << cfg.response().analysis_space
              << " observable_bins=" << summary.observable_bins
              << " exposure=" << summary.exposure_kg_year << " kg*year\n";

    auto response = MakeResponseFold(cfg, summary, config_path);
    ResponseFold& fold = *response.fold;
    const std::size_t n_bins_full = fold.NumBins();  // output bins of the fold, before any single-bin collapse

    auto background = MakeBackground(cfg, summary, fold, response);

    // Observed data (raw, before single_bin_likelihood collapse). Precedence
    // matches the reference app exactly: run.observed_counts (inline config
    // array) overrides run.data_path (CSV/ROOT file), which overrides Asimov
    // (data = background).
    const auto& roi_ids = use_pattern_bins ? summary.pattern_roi : summary.roi_bins;
    std::vector<double> data_raw = ObservedCountsFromConfig(jroot, n_bins_full);  // observed counts per bin (empty = not given)
    if (data_raw.empty() && !run.data_path.empty()) {
      data_raw = LoadData(run.data_path, n_bins_full, roi_ids);
      if (data_raw.size() < n_bins_full) {
        throw std::runtime_error("data_path \"" + run.data_path + "\" has " + std::to_string(data_raw.size()) +
                                  " values, need " + std::to_string(n_bins_full));
      }
      data_raw.resize(n_bins_full);
      if (run.verbosity >= 1) std::cout << "[scan-generic] Real data loaded from " << run.data_path << "\n";
    }
    if (data_raw.empty()) data_raw = background.B_pat;

    const bool single_bin = run.single_bin_likelihood;
    const std::size_t n_bins_pl = single_bin ? 1u : n_bins_full;  // bins actually entering the likelihood

    std::vector<double> data_pl = single_bin ? std::vector<double>{SumVec(data_raw)} : data_raw;
    std::vector<double> B_pat_pl = single_bin ? std::vector<double>{SumVec(background.B_pat)} : background.B_pat;
    std::vector<double> Bp_pl = single_bin ? std::vector<double>{SumVec(background.Bp)} : background.Bp;
    std::vector<double> Br_pl = single_bin ? std::vector<double>{SumVec(background.Br)} : background.Br;

    ccdarksens::stats::ProfileLikelihood profile_pl;  // likelihood with the background nuisance parameter
    profile_pl.SetData(data_pl);
    const bool use_bp_br = (run.background_model == "Bp_theta_Br" || run.background_model == "Bp_theta_br");
    if (use_bp_br) {
      profile_pl.SetBpBr(Bp_pl, Br_pl);
      profile_pl.SetConstrainPriorStrength(run.constrain_prior_strength);
      profile_pl.SetConstrainGammaSign(run.constrain_use_gamma_sign);
      profile_pl.SetConstrainUseTauWeighted(run.constrain_use_tau_weighted);
      profile_pl.SetConstrainNBins(run.constrain_n_bins);
      profile_pl.SetAcceptBoundary2DMinimum(run.pydme_style_ul);
    } else {
      profile_pl.SetBTemplate(B_pat_pl);
    }
    double profile_param_lo = 0.01, profile_param_hi = 10.0;  // bounds of the nuisance parameter (scale or theta)
    if (profile_pl.UseBpBr()) {
      profile_param_lo = run.theta_lo;
      profile_param_hi = run.theta_hi;
    }
    const std::vector<double> profile_S_null(n_bins_pl, 0.0);  // background-only signal vector
    const double target_q = std::pow(TMath::NormQuantile(run.cl), 2);  // q threshold for the one-sided limit at confidence level cl

    // ---- Grid ----
    const auto& mj = cfg.model();
    const auto& jgrid = jroot.at("model").at("grid");
    std::vector<double> mchi_list, sigma_list;  // mass axis and coupling axis of the scan
    if (mj.type == "dark_photon") {
      mchi_list = ExpandAxisFromGrid(jgrid, {"mA_eV", "mchi_MeV"});
      sigma_list = ExpandAxisFromGrid(jgrid, {"epsilon", "sigma_e_cm2"});
    } else if (mj.type == "migdal" || mj.type == "wimp_nucleon") {
      mchi_list = ExpandAxisFromGrid(jgrid, {"mchi_MeV"});
      sigma_list = ExpandAxisFromGrid(jgrid, {"sigma_n_cm2", "sigma_e_cm2"});
    } else {
      mchi_list = ExpandAxisFromGrid(jgrid, {"mchi_MeV"});
      sigma_list = ExpandAxisFromGrid(jgrid, {"sigma_e_cm2"});
    }
    if (mchi_list.empty() || sigma_list.empty()) throw std::runtime_error("Empty mass or coupling grid.");
    std::string fmt_sigma = jgrid.contains("format") ? jgrid["format"].value("sigma", std::string(".1e")) : ".1e";

    // ---- run.q_target_lookup_path: optional per-mass PLR threshold override,
    // ported from ccdarksens_scan_srdm_pattern_csv.cc so band.cc's toy-MC
    // calibration works against this app too. ----
    std::map<double, double> q_target_by_mass;  // per-mass PLR threshold from a previous threshold_toys run (empty = use target_q)
    if (!run.q_target_lookup_path.empty()) {
      TFile fql(run.q_target_lookup_path.c_str(), "READ");
      if (!fql.IsOpen() || fql.IsZombie()) {
        throw std::runtime_error("Cannot open q_target_lookup_path: " + run.q_target_lookup_path);
      }
      TGraph* g_qt = dynamic_cast<TGraph*>(fql.Get("q_target_per_mass"));
      if (!g_qt) {
        throw std::runtime_error("q_target_lookup_path missing TGraph 'q_target_per_mass': " + run.q_target_lookup_path);
      }
      for (int i = 0; i < g_qt->GetN(); ++i) {
        double x = 0.0, y = 0.0;
        g_qt->GetPoint(i, x, y);
        q_target_by_mass[x] = y;
      }
      fql.Close();
      if (run.verbosity >= 1) {
        std::cout << "[scan-generic] Loaded q_target lookup with " << q_target_by_mass.size() << " entries from "
                  << run.q_target_lookup_path << "\n";
      }
    }
    auto target_q_for_mass = [&](double mchi) -> double {
      if (q_target_by_mass.empty()) return target_q;
      const double rel_tol_mass = 1e-9;
      for (const auto& kv : q_target_by_mass) {
        if (std::abs(kv.first - mchi) <= rel_tol_mass * std::max(1.0, std::abs(mchi))) return kv.second;
      }
      double best = target_q, best_diff = std::numeric_limits<double>::infinity();
      for (const auto& kv : q_target_by_mass) {
        const double d = std::abs(kv.first - mchi);
        if (d < best_diff) {
          best_diff = d;
          best = kv.second;
        }
      }
      return best;
    };

    // ---- run.mode == "threshold_toys": band.cc's toy-MC threshold
    // calibration phase (early exit before the normal grid loop), ported
    // near-verbatim from ccdarksens_scan_srdm_pattern_csv.cc. ----
    if (run.mode == "threshold_toys") {
      const auto& tt = run.threshold_toys;
      if (tt.sigma_threshold_graph_path.empty()) {
        throw std::runtime_error("run.threshold_toys.sigma_threshold_graph_path is required");
      }
      if (tt.n_threshold_toys < 100) {
        throw std::runtime_error("threshold_toys.n_threshold_toys must be >= 100 (got " +
                                  std::to_string(tt.n_threshold_toys) + ")");
      }
      if (tt.percentile <= 0.0 || tt.percentile >= 1.0) {
        throw std::runtime_error("threshold_toys.percentile must be in (0,1), got " + std::to_string(tt.percentile));
      }

      std::map<double, double> sigma_thr_by_mass;
      {
        TFile fst(tt.sigma_threshold_graph_path.c_str(), "READ");
        if (!fst.IsOpen() || fst.IsZombie()) {
          throw std::runtime_error("Cannot open sigma_threshold_graph_path: " + tt.sigma_threshold_graph_path);
        }
        TGraph* g_st = dynamic_cast<TGraph*>(fst.Get("sigma_threshold_per_mass"));
        if (!g_st) {
          throw std::runtime_error("sigma_threshold_graph_path missing TGraph 'sigma_threshold_per_mass': " +
                                    tt.sigma_threshold_graph_path);
        }
        for (int i = 0; i < g_st->GetN(); ++i) {
          double x = 0.0, y = 0.0;
          g_st->GetPoint(i, x, y);
          sigma_thr_by_mass[x] = y;
        }
        fst.Close();
      }
      auto sigma_thr_for_mass = [&](double mchi) -> double {
        const double rel_tol_mass = 1e-9;
        for (const auto& kv : sigma_thr_by_mass) {
          if (std::abs(kv.first - mchi) <= rel_tol_mass * std::max(1.0, std::abs(mchi))) return kv.second;
        }
        double best = -1.0, best_diff = std::numeric_limits<double>::infinity();
        for (const auto& kv : sigma_thr_by_mass) {
          const double d = std::abs(kv.first - mchi);
          if (d < best_diff) {
            best_diff = d;
            best = kv.second;
          }
        }
        return best;
      };

      if (run.verbosity >= 1) {
        std::cout << "[scan-generic] mode=threshold_toys n_threshold_toys=" << tt.n_threshold_toys
                  << " percentile=" << tt.percentile << " seed=" << tt.rng_seed << "\n";
      }

      std::vector<double> qt_mchi, qt_value;
      qt_mchi.reserve(mchi_list.size());
      qt_value.reserve(mchi_list.size());

      int ix = 0;
      for (double mchi : mchi_list) {
        const double sigma_thr = sigma_thr_for_mass(mchi);
        if (!(sigma_thr > 0.0)) {
          throw std::runtime_error("sigma_threshold for mchi=" + std::to_string(mchi) + " is non-positive in " +
                                    tt.sigma_threshold_graph_path);
        }

        std::vector<std::vector<double>> S_grid;
        S_grid.reserve(sigma_list.size());
        for (double sigma_val : sigma_list) {
          bool sig_ok = false;
          auto dRdE_sig = MakeSignalSpectrumE(mj, mchi, FormatSigma(sigma_val, fmt_sigma), &sig_ok);
          std::vector<double> S = fold.Fold(*dRdE_sig, summary.exposure_kg_year);
          S_grid.push_back(single_bin ? std::vector<double>{SumVec(S)} : S);
        }
        const double log10_lo = std::log10(sigma_list.front());
        const double log10_hi = std::log10(sigma_list.back());
        const double log10_sigma_thr_seed = std::log10(sigma_thr);
        std::vector<double> S_thr = ccdarksens::scan::InterpolateSignal(S_grid, sigma_list, log10_sigma_thr_seed);

        std::vector<double> lambda(n_bins_pl, 0.0);
        for (std::size_t i = 0; i < lambda.size(); ++i) {
          lambda[i] = S_thr[i] + Bp_pl[i] + Br_pl[i];
          if (lambda[i] < 0.0) throw std::runtime_error("Negative Poisson mean encountered in threshold_toys.");
        }

        auto S_from_log10 = [&S_grid, &sigma_list, log10_lo, log10_hi](double log10_s) -> std::vector<double> {
          const std::size_t n = S_grid.empty() ? 0u : S_grid[0].size();
          std::vector<double> out(n, 0.0);
          if (S_grid.size() < 2u) return S_grid.empty() ? out : S_grid[0];
          const double log10_s_clamp = std::max(log10_lo, std::min(log10_hi, log10_s));
          for (std::size_t j = 0; j + 1 < sigma_list.size(); ++j) {
            const double l0 = std::log10(sigma_list[j]);
            const double l1 = std::log10(sigma_list[j + 1]);
            if (log10_s_clamp >= l0 && log10_s_clamp <= l1) {
              const double t = (l1 - l0) > 1e-300 ? (log10_s_clamp - l0) / (l1 - l0) : 0.0;
              for (std::size_t b = 0; b < n; ++b) out[b] = S_grid[j][b] + t * (S_grid[j + 1][b] - S_grid[j][b]);
              return out;
            }
          }
          if (log10_s_clamp <= std::log10(sigma_list[0])) return S_grid[0];
          return S_grid.back();
        };

        const bool do_2d_fit_here =
            use_pattern_bins && (run.profile_minimizer == "minuit2d" || run.profile_minimizer == "pydme") &&
            sigma_list.size() >= 2u;

        std::seed_seq seq{static_cast<std::uint32_t>(tt.rng_seed & 0xFFFFFFFFu),
                          static_cast<std::uint32_t>((tt.rng_seed >> 32) & 0xFFFFFFFFu),
                          static_cast<std::uint32_t>(ix), 0xC0FFEEu};
        std::mt19937_64 rng(seq);
        std::vector<std::poisson_distribution<long long>> poisson_per_bin;
        poisson_per_bin.reserve(lambda.size());
        for (double mu : lambda) poisson_per_bin.emplace_back(mu);

        std::vector<double> q_values(static_cast<std::size_t>(tt.n_threshold_toys), 0.0);
        for (long long k = 0; k < tt.n_threshold_toys; ++k) {
          std::vector<double> nk(lambda.size(), 0.0);
          for (std::size_t i = 0; i < nk.size(); ++i) nk[i] = static_cast<double>(poisson_per_bin[i](rng));

          profile_pl.SetData(nk);
          const auto pr_top = profile_pl.MinimizeOverScaleMinuit(S_thr, profile_param_lo, profile_param_hi);
          const double theta_hat_top = pr_top.first;
          const double nll_top = pr_top.second;

          double nll_glob = nll_top;
          if (do_2d_fit_here && S_grid.size() >= 2u) {
            auto res = profile_pl.MinimizeOverSigmaAndTheta(S_from_log10, log10_lo, log10_hi, profile_param_lo,
                                                              profile_param_hi, log10_sigma_thr_seed, theta_hat_top);
            if (res.ok && res.nll_min < nll_glob) nll_glob = res.nll_min;
          }
          const double q = 2.0 * (nll_top - nll_glob);
          q_values[static_cast<std::size_t>(k)] = (q > 0.0) ? q : 0.0;
        }

        const std::size_t q_idx =
            static_cast<std::size_t>(tt.percentile * static_cast<double>(tt.n_threshold_toys - 1));
        std::nth_element(q_values.begin(), q_values.begin() + static_cast<std::ptrdiff_t>(q_idx), q_values.end());
        const double c_mu = q_values[q_idx];
        qt_mchi.push_back(mchi);
        qt_value.push_back(c_mu);
        if (run.verbosity >= 1) {
          std::cout << "[scan-generic] threshold_toys mchi=" << mchi << " sigma_thr=" << sigma_thr
                    << " q_target(p=" << tt.percentile << ")=" << c_mu << "\n";
        }
        // restore SetData to the nominal data for consistency, in case this
        // profile_pl instance is reused (it currently is not, but be explicit).
        profile_pl.SetData(data_pl);
        ++ix;
      }

      const std::string thr_out_path = run.outdir + "/qtarget_threshold.root";
      TFile fout_thr(thr_out_path.c_str(), "RECREATE");
      if (!fout_thr.IsOpen()) throw std::runtime_error("Failed to create output file: " + thr_out_path);
      TGraph g_qt(static_cast<int>(qt_mchi.size()));
      g_qt.SetName("q_target_per_mass");
      g_qt.SetTitle(";m_{#chi} [MeV];q_{#mu, target} (toy-MC threshold)");
      for (int i = 0; i < static_cast<int>(qt_mchi.size()); ++i) {
        g_qt.SetPoint(i, qt_mchi[static_cast<std::size_t>(i)], qt_value[static_cast<std::size_t>(i)]);
      }
      g_qt.Write("q_target_per_mass");
      TParameter<double>("threshold_percentile", tt.percentile).Write();
      TParameter<double>("n_threshold_toys", static_cast<double>(tt.n_threshold_toys)).Write();
      fout_thr.Close();
      std::cout << "[scan-generic] Wrote " << thr_out_path << "\n";
      return 0;
    }

    auto edges_m = MakeEdgesFromCenters(mchi_list);
    auto edges_s = MakeEdgesFromCenters(sigma_list);
    TH2D h_q("q_mchi_sigma", ";m_{#chi} [MeV];coupling;q_{Asimov}", static_cast<int>(mchi_list.size()), edges_m.data(),
             static_cast<int>(sigma_list.size()), edges_s.data());
    TH1D h_upper_limit("upper_limit_sigma_e_mchi", ";m_{#chi} [MeV];coupling (90% CL)",
                        static_cast<int>(mchi_list.size()), edges_m.data());

    // ---- Phase B: grid loop (1D brent/minuit, and pattern-space-only 2D
    // minuit2d/pydme pre-fit + pydme continuous bisection -- matching the
    // reference app, which only runs the 2D fit for pattern-space) ----
    const bool pydme_mode = (run.profile_minimizer == "pydme");  // true = also run the continuous pydme bisection
    std::vector<double> ul_pydme_bisection_vals(mchi_list.size(), -1.0);  // pydme-style upper limit per mass (-1 = not available)

    int idx_mchi = 0;
    for (double mchi : mchi_list) {
      double nll_null;
      if (run.profile_minimizer == "minuit" || run.profile_minimizer == "minuit2d" ||
          run.profile_minimizer == "pydme") {
        nll_null = profile_pl.MinimizeOverScaleMinuit(profile_S_null, profile_param_lo, profile_param_hi).second;
      } else {
        nll_null = profile_pl.MinimizeOverScale(profile_S_null, profile_param_lo, profile_param_hi).second;
      }

      // Pydme-style: one 2D Minuit fit over (log10(sigma), theta) per mchi,
      // giving the true global NLL minimum; the per-sigma S_grid computed
      // here is then reused (not re-folded) inside the sigma loop below.
      const bool do_2d_fit = use_pattern_bins &&
          (run.profile_minimizer == "minuit2d" || run.profile_minimizer == "pydme") && sigma_list.size() >= 2u;
      double nll_min_2d = std::numeric_limits<double>::max();
      bool use_2d_nll_min = false;
      double log10_sigma_hat_pydme = std::log10(sigma_list.front());
      std::vector<std::vector<double>> S_grid;
      if (do_2d_fit) {
        S_grid.reserve(sigma_list.size());
        for (double sigma_val : sigma_list) {
          bool sig_ok = false;
          auto dRdE_sig = MakeSignalSpectrumE(mj, mchi, FormatSigma(sigma_val, fmt_sigma), &sig_ok);
          std::vector<double> S = fold.Fold(*dRdE_sig, summary.exposure_kg_year);
          S_grid.push_back(single_bin ? std::vector<double>{SumVec(S)} : S);
        }
        const double log10_lo = std::log10(sigma_list.front());
        const double log10_hi = std::log10(sigma_list.back());
        auto S_from_log10 = [&S_grid, &sigma_list, log10_lo, log10_hi](double log10_s) -> std::vector<double> {
          const std::size_t n = S_grid.empty() ? 0u : S_grid[0].size();
          std::vector<double> out(n, 0.0);
          if (S_grid.size() < 2u) return S_grid.empty() ? out : S_grid[0];
          const double log10_s_clamp = std::max(log10_lo, std::min(log10_hi, log10_s));
          for (std::size_t j = 0; j + 1 < sigma_list.size(); ++j) {
            const double l0 = std::log10(sigma_list[j]);
            const double l1 = std::log10(sigma_list[j + 1]);
            if (log10_s_clamp >= l0 && log10_s_clamp <= l1) {
              const double t = (l1 - l0) > 1e-300 ? (log10_s_clamp - l0) / (l1 - l0) : 0.0;
              for (std::size_t b = 0; b < n; ++b) out[b] = S_grid[j][b] + t * (S_grid[j + 1][b] - S_grid[j][b]);
              return out;
            }
          }
          if (log10_s_clamp <= std::log10(sigma_list[0])) return S_grid[0];
          return S_grid.back();
        };
        auto result = profile_pl.MinimizeOverSigmaAndTheta(S_from_log10, log10_lo, log10_hi, profile_param_lo,
                                                             profile_param_hi);
        if (result.ok) {
          nll_min_2d = result.nll_min;
          use_2d_nll_min = true;
          log10_sigma_hat_pydme = result.log10_sigma_hat;
          if (run.verbosity >= 1) {
            std::cout << "[scan-generic] mchi=" << mchi << " " << run.profile_minimizer
                      << ": log10_sigma_hat=" << result.log10_sigma_hat << " theta_hat=" << result.theta_hat
                      << " nll_min=" << result.nll_min << "\n";
          }
        }
      }

      std::vector<double> nll_values;
      nll_values.reserve(sigma_list.size());
      int idx_sigma = 0;
      for (double sigma_val : sigma_list) {
        if (do_2d_fit && idx_sigma < static_cast<int>(S_grid.size())) {
          const std::vector<double>& S_for_pl = S_grid[static_cast<std::size_t>(idx_sigma)];
          const double nll_m = profile_pl.MinimizeOverScaleMinuit(S_for_pl, profile_param_lo, profile_param_hi).second;
          const double nll_b = profile_pl.MinimizeOverScale(S_for_pl, profile_param_lo, profile_param_hi).second;
          nll_values.push_back(std::min(nll_m, nll_b));
          ++idx_sigma;
          continue;
        }

        bool sig_ok = false;
        auto dRdE_sig = MakeSignalSpectrumE(mj, mchi, FormatSigma(sigma_val, fmt_sigma), &sig_ok);
        if (!sig_ok && run.verbosity >= 1) {
          std::cout << "[scan-generic] WARNING: failed to load rate table for mchi=" << mchi
                     << ", coupling=" << sigma_val << "; using all-zero signal.\n";
        }
        // MakeSignalSpectrumE always returns a valid (possibly all-zero) TH1D
        // even on load failure -- fold it through normally rather than
        // special-casing; an all-zero signal correctly yields q=0 there.
        std::vector<double> S = fold.Fold(*dRdE_sig, summary.exposure_kg_year);
        std::vector<double> S_for_pl = single_bin ? std::vector<double>{SumVec(S)} : S;

        double nll;
        if (use_pattern_bins) {
          // Matches the reference app's pattern-space branch: "minuit" also
          // cross-checks against Brent and takes the better (never-worse) fit.
          if (run.profile_minimizer == "minuit") {
            const double nll_m = profile_pl.MinimizeOverScaleMinuit(S_for_pl, profile_param_lo, profile_param_hi).second;
            const double nll_b = profile_pl.MinimizeOverScale(S_for_pl, profile_param_lo, profile_param_hi).second;
            nll = std::min(nll_m, nll_b);
          } else {
            nll = profile_pl.MinimizeOverScale(S_for_pl, profile_param_lo, profile_param_hi).second;
          }
        } else {
          // n_e-space branch: matches the reference app exactly -- "minuit"/
          // "pydme" use Minuit only here (no Brent cross-check).
          if (run.profile_minimizer == "minuit" || run.profile_minimizer == "pydme") {
            nll = profile_pl.MinimizeOverScaleMinuit(S_for_pl, profile_param_lo, profile_param_hi).second;
          } else {
            nll = profile_pl.MinimizeOverScale(S_for_pl, profile_param_lo, profile_param_hi).second;
          }
        }
        nll_values.push_back(nll);
        ++idx_sigma;
      }

      const double nll_min_grid = *std::min_element(nll_values.begin(), nll_values.end());
      // Use the 2D global minimum when pydme_style_ul requests it and the 2D
      // fit succeeded; else the grid minimum, to avoid a spurious flat curve
      // (matches the reference app exactly).
      const double nll_min = (run.pydme_style_ul && use_2d_nll_min) ? nll_min_2d : nll_min_grid;
      const double target_q_mchi = target_q_for_mass(mchi);
      const auto q_mono = ccdarksens::scan::MonotonizeQ(nll_values, nll_min);
      const double ul_sigma = ccdarksens::scan::UlFromQMonoCrossing(sigma_list, q_mono, target_q_mchi);

      double ul_pydme_bisection = -1.0;
      if (pydme_mode && S_grid.size() >= 2u) {
        const double log10_lo = std::log10(sigma_list.front());
        const double log10_hi = std::log10(sigma_list.back());
        auto q_mu_at = [&](double log10_sig) {
          std::vector<double> S = ccdarksens::scan::InterpolateSignal(S_grid, sigma_list, log10_sig);
          const double nll = profile_pl.MinimizeOverScaleMinuit(S, profile_param_lo, profile_param_hi).second;
          return 2.0 * (nll - nll_min);
        };
        const int k_best = static_cast<int>(std::min_element(nll_values.begin(), nll_values.end()) - nll_values.begin());
        const double log10_sigma_best = std::log10(sigma_list[static_cast<std::size_t>(k_best)]);
        const double log10_seed = use_2d_nll_min ? std::max({log10_sigma_hat_pydme, log10_lo, log10_sigma_best})
                                                  : std::max(log10_lo, log10_sigma_best);
        const double log10_ul =
            ccdarksens::scan::BisectUpperLimit(q_mu_at, log10_lo, log10_hi, log10_seed, target_q_mchi);
        if (log10_ul < log10_hi) ul_pydme_bisection = std::pow(10.0, log10_ul);
      }
      if (ul_pydme_bisection > 0.0) ul_pydme_bisection_vals[static_cast<std::size_t>(idx_mchi)] = ul_pydme_bisection;

      for (int k = 0; k < static_cast<int>(nll_values.size()); ++k) {
        if (q_mono[static_cast<std::size_t>(k)] <= 0.0) continue;
        h_q.SetBinContent(idx_mchi + 1, k + 1, q_mono[static_cast<std::size_t>(k)]);
      }
      h_upper_limit.SetBinContent(idx_mchi + 1, ul_sigma);
      if (run.verbosity >= 1) {
        std::cout << "[scan-generic] mchi=" << mchi << " -> upper_limit=" << ul_sigma << "\n";
      }
      ++idx_mchi;
    }

    // ---- Phase C: post-processing + output ----
    if (run.smooth_ul_envelope) {
      constexpr int kEnvW = 5;
      const int n_mchi = h_upper_limit.GetNbinsX();
      std::vector<double> raw(static_cast<std::size_t>(n_mchi));
      for (int i = 0; i < n_mchi; ++i) raw[static_cast<std::size_t>(i)] = h_upper_limit.GetBinContent(i + 1);
      std::vector<double> env(static_cast<std::size_t>(n_mchi), 0.0);
      for (int i = 0; i < n_mchi; ++i) {
        const int lo = std::max(0, i - kEnvW);
        const int hi = std::min(n_mchi - 1, i + kEnvW);
        double win_min = raw[static_cast<std::size_t>(i)];
        for (int j = lo; j <= hi; ++j) {
          const double v = raw[static_cast<std::size_t>(j)];
          if (v > 0.0 && (win_min <= 0.0 || v < win_min)) win_min = v;
        }
        env[static_cast<std::size_t>(i)] = win_min;
      }
      for (int i = 0; i < n_mchi; ++i) h_upper_limit.SetBinContent(i + 1, env[static_cast<std::size_t>(i)]);
    }

    const int n_mchi = h_upper_limit.GetNbinsX();
    TGraph g_upper_limit(n_mchi);
    for (int i = 0; i < n_mchi; ++i) {
      g_upper_limit.SetPoint(i, mchi_list[static_cast<std::size_t>(i)], h_upper_limit.GetBinContent(i + 1));
    }
    g_upper_limit.SetName("upper_limit_sigma_e_mchi_graph");

    std::vector<double> pydme_m, pydme_s;
    for (int i = 0; i < n_mchi; ++i) {
      const double s = ul_pydme_bisection_vals[static_cast<std::size_t>(i)];
      if (s > 0.0) {
        pydme_m.push_back(mchi_list[static_cast<std::size_t>(i)]);
        pydme_s.push_back(s);
      }
    }
    TGraph g_pydme_bisection(static_cast<int>(pydme_m.size()));
    if (!pydme_m.empty()) {
      g_pydme_bisection.SetName("upper_limit_sigma_e_mchi_pydme_bisection");
      for (int i = 0; i < static_cast<int>(pydme_m.size()); ++i) {
        g_pydme_bisection.SetPoint(i, pydme_m[static_cast<std::size_t>(i)], pydme_s[static_cast<std::size_t>(i)]);
      }
    }

    const std::string out_path = run.outdir + "/scan_generic.root";
    TFile fout(out_path.c_str(), "RECREATE");
    if (!fout.IsOpen()) {
      std::cerr << "[scan-generic] ERROR: cannot create output file " << out_path << "\n";
      return 1;
    }
    h_q.Write("q_mchi_sigma");
    h_upper_limit.Write("upper_limit_sigma_e_mchi");
    g_upper_limit.Write("upper_limit_sigma_e_mchi_graph");
    if (!pydme_m.empty()) g_pydme_bisection.Write("upper_limit_sigma_e_mchi_pydme_bisection");
    TParameter<double> p_exp("exposure_kg_year", summary.exposure_kg_year);
    p_exp.Write();
    fout.Close();

    std::cout << "[scan-generic] Done. Output written to " << out_path << "\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "[scan-generic] ERROR: " << e.what() << "\n";
    return 1;
  }
}
