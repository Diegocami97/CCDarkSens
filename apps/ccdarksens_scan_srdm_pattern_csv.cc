// ============================================================================
//  CCDarkSens — ccdarksens_scan_srdm_pattern_csv
//  Grid scan over SRDM pattern-space signal CSVs that evaluates profile-likelihood q(mχ,σ) and extracts upper limits, including a threshold-toys mode for band workflows.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TParameter.h>
#include <TMath.h>

#include <nlohmann/json.hpp>

#include "ccdarksens/stats/ProfileLikelihood.hh"

using nlohmann::json;

namespace {

static std::vector<double> make_edges_from_centers(const std::vector<double>& c) {
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
  for (std::size_t i = 1; i < N; ++i) {
    edges[i] = 0.5 * (c[i - 1] + c[i]);
  }
  edges[N] = c[N - 1] + 0.5 * (c[N - 1] - c[N - 2]);
  return edges;
}

static std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::string token;
  std::stringstream ss(line);
  while (std::getline(ss, token, ',')) out.push_back(token);
  return out;
}

static bool approx_equal(double a, double b, double rel_tol = 1e-12) {
  // Works well for scientific-notation sigmas around 1e-38.
  const double denom = std::max(1e-300, std::abs(a));
  return std::abs(a - b) <= rel_tol * denom;
}

struct SrdmMassData {
  double mchi = 0.0;
  std::vector<double> sigmas;                // xsec values (ascending)
  std::vector<std::vector<double>> Spat;   // Spat[i] corresponds to sigmas[i]
};

static int find_column(const std::vector<std::string>& header, const std::string& col_name) {
  for (std::size_t i = 0; i < header.size(); ++i) {
    if (header[i] == col_name) return static_cast<int>(i);
  }
  return -1;
}

static std::optional<std::pair<double, std::string>> parse_mchi_and_tag(
    const std::string& filename,
    const std::regex& mchi_regex,
    std::string* full_match_out = nullptr) {
  std::smatch m;
  if (!std::regex_search(filename, m, mchi_regex)) return std::nullopt;
  if (m.size() < 2) return std::nullopt;

  try {
    const double mchi = std::stod(m[1].str());
    if (full_match_out) {
      // m[0] is the full match; tag is used only for debugging/logging.
      *full_match_out = m[0].str();
    }
    return std::make_pair(mchi, filename);
  } catch (...) {
    return std::nullopt;
  }
}

static SrdmMassData load_srdm_csv(
    const std::string& csv_path,
    double mchi,
    const std::vector<int>& pattern_roi,
    double signal_rate_scale) {
  std::ifstream in(csv_path);
  if (!in.is_open()) throw std::runtime_error("Failed to open SRDM CSV: " + csv_path);

  std::string header_line;
  if (!std::getline(in, header_line)) throw std::runtime_error("Empty CSV: " + csv_path);

  const auto header = split_csv_line(header_line);
  if (header.empty()) throw std::runtime_error("Missing CSV header: " + csv_path);

  const int xsec_col = find_column(header, "xsec");
  if (xsec_col < 0) throw std::runtime_error("CSV missing 'xsec' column: " + csv_path);

  // Build column indices for each requested pattern bin, in pattern_roi order.
  std::vector<int> spat_cols;
  spat_cols.reserve(pattern_roi.size());
  for (int pid : pattern_roi) {
    const std::string col = "S" + std::to_string(pid);
    const int idx = find_column(header, col);
    if (idx < 0) {
      throw std::runtime_error("CSV " + csv_path + " missing column '" + col +
                               "' for pattern_roi entry " + std::to_string(pid));
    }
    spat_cols.push_back(idx);
  }

  SrdmMassData d;
  d.mchi = mchi;

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    const auto toks = split_csv_line(line);
    if (toks.size() <= static_cast<std::size_t>(xsec_col)) continue;

    double sigma = 0.0;
    try {
      sigma = std::stod(toks[xsec_col]);
    } catch (...) {
      continue;
    }

    std::vector<double> svec(pattern_roi.size(), 0.0);
    for (std::size_t j = 0; j < pattern_roi.size(); ++j) {
      const int col = spat_cols[j];
      if (toks.size() <= static_cast<std::size_t>(col)) {
        throw std::runtime_error("CSV " + csv_path + " malformed row (too few columns).");
      }
      try {
        svec[j] = std::stod(toks[col]) * signal_rate_scale;
      } catch (...) {
        svec[j] = 0.0;
      }
    }

    d.sigmas.push_back(sigma);
    d.Spat.push_back(std::move(svec));
  }

  // Ensure sigma order (CSV should already be ascending, but make it deterministic).
  if (d.sigmas.size() != d.Spat.size()) throw std::runtime_error("Internal SRDM CSV parse error.");

  std::vector<std::size_t> order(d.sigmas.size());
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return d.sigmas[a] < d.sigmas[b]; });

  std::vector<double> sig_sorted;
  std::vector<std::vector<double>> spat_sorted;
  sig_sorted.reserve(d.sigmas.size());
  spat_sorted.reserve(d.Spat.size());
  for (std::size_t idx : order) {
    sig_sorted.push_back(d.sigmas[idx]);
    spat_sorted.push_back(std::move(d.Spat[idx]));
  }
  d.sigmas = std::move(sig_sorted);
  d.Spat = std::move(spat_sorted);

  return d;
}

static int find_sigma_index(const std::vector<double>& sigmas, double target, double rel_tol = 1e-12) {
  for (std::size_t i = 0; i < sigmas.size(); ++i) {
    if (approx_equal(sigmas[i], target, rel_tol)) return static_cast<int>(i);
  }
  return -1;
}

// Interpolate Spat between neighboring sigma points in log10(sigma) space.
// - Exact matches (within rel_tol) return the stored template.
// - Out-of-range sigmas are clamped to the nearest available template.
static std::vector<double> interpolate_spat_on_sigma(
    const SrdmMassData& mass,
    double sigma,
    double rel_tol = 1e-12) {
  if (sigma <= 0.0) throw std::runtime_error("Sigma must be > 0 for interpolation.");
  if (mass.sigmas.empty() || mass.Spat.empty()) return {};
  if (mass.sigmas.size() != mass.Spat.size()) {
    throw std::runtime_error("Internal error: mass.sigmas/mass.Spat size mismatch.");
  }

  const auto& sigmas = mass.sigmas;

  const int idx_exact = find_sigma_index(sigmas, sigma, rel_tol);
  if (idx_exact >= 0) return mass.Spat[static_cast<std::size_t>(idx_exact)];

  // If only one point exists, nothing to interpolate.
  if (sigmas.size() == 1) return mass.Spat.front();

  if (sigma <= sigmas.front()) return mass.Spat.front();
  if (sigma >= sigmas.back()) return mass.Spat.back();

  // Find bracketing points: sigmas[i_low] <= sigma <= sigmas[i_high]
  const auto it_upper = std::upper_bound(sigmas.begin(), sigmas.end(), sigma);
  const std::size_t i_high = static_cast<std::size_t>(it_upper - sigmas.begin());
  if (i_high == 0u || i_high >= sigmas.size()) {
    // Should be impossible due to the clamping checks above, but keep it safe.
    return (sigma < sigmas.front()) ? mass.Spat.front() : mass.Spat.back();
  }
  const std::size_t i_low = i_high - 1u;

  const double s_low = sigmas[i_low];
  const double s_high = sigmas[i_high];
  if (s_low <= 0.0 || s_high <= 0.0) throw std::runtime_error("Non-positive sigma encountered in templates.");

  const double log10_s = std::log10(sigma);
  const double log10_low = std::log10(s_low);
  const double log10_high = std::log10(s_high);

  const double denom = (log10_high - log10_low);
  const double t = (std::abs(denom) > 1e-300) ? ((log10_s - log10_low) / denom) : 0.0;

  const auto& Spat_low = mass.Spat[i_low];
  const auto& Spat_high = mass.Spat[i_high];
  if (Spat_low.size() != Spat_high.size()) {
    throw std::runtime_error("Internal error: Spat template size mismatch across sigmas.");
  }

  std::vector<double> out(Spat_low.size(), 0.0);
  for (std::size_t b = 0; b < out.size(); ++b) {
    out[b] = Spat_low[b] + t * (Spat_high[b] - Spat_low[b]);
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json\n";
    return 1;
  }

  const std::string config_path = argv[1];
  try {

    // =========================================================================
    // 1. Parse JSON configuration
    // =========================================================================
    std::ifstream jf(config_path);
    if (!jf.is_open()) throw std::runtime_error("Cannot open config file: " + config_path);
    json jroot;
    jf >> jroot;

    const auto run = jroot.at("run");
    const auto exp = jroot.at("experiment");
    const auto det = jroot.at("detector");
    const auto sig = jroot.at("srdm_signal_csv");

    const bool use_profile_likelihood = run.value("use_profile_likelihood", false);
    if (!use_profile_likelihood)
      throw std::runtime_error("This SRDM CSV app currently requires run.use_profile_likelihood=true");

    const std::string outdir = run.at("outdir").get<std::string>();
    const double cl = run.value("cl", 0.9);
    const int verbosity = run.value("verbosity", 1);

    // Optional: dump Spat (pattern-space expected signal) for one (mchi, sigma)
    // point to validate units / templates.
    const double dump_spat_point_mchi_MeV = run.value("dump_spat_point_mchi_MeV", -1.0);
    const double dump_spat_point_sigma_e_cm2 = run.value("dump_spat_point_sigma_e_cm2", -1.0);
    const bool dump_spat_point =
        (dump_spat_point_mchi_MeV > 0.0 && dump_spat_point_sigma_e_cm2 > 0.0);

    const std::string background_source = run.value("background_source", std::string("bp_br_template"));
    const std::string background_model = run.value("background_model", std::string("Bp_theta_Br"));
    if (background_source != "bp_br_template" || background_model != "Bp_theta_Br") {
      throw std::runtime_error(
          "This SRDM CSV app only supports background_source='bp_br_template' and background_model='Bp_theta_Br'.");
    }

    const std::vector<int> pattern_roi = exp.at("pattern_roi").get<std::vector<int>>();
    if (pattern_roi.empty()) throw std::runtime_error("experiment.pattern_roi is empty.");

    // Background templates: Bp and Br per pattern bin (same order as pattern_roi).
    const std::vector<double> background_Bp = run.at("background_Bp").get<std::vector<double>>();
    const std::vector<double> background_Br = run.at("background_Br").get<std::vector<double>>();
    if (background_Bp.size() != pattern_roi.size() || background_Br.size() != pattern_roi.size()) {
      throw std::invalid_argument("background_Bp/background_Br size must match experiment.pattern_roi size.");
    }

    // Optional: observed pattern counts (one per pattern_roi bin) to run an observed
    // limit instead of an Asimov sensitivity. When omitted, data is set to Bp+Br
    // (theta=1 nominal) so q_mu corresponds to the median expected sensitivity.
    std::vector<double> observed_counts;
    bool use_observed_data = false;
    if (run.contains("observed_counts") && !run.at("observed_counts").is_null()) {
      observed_counts = run.at("observed_counts").get<std::vector<double>>();
      if (!observed_counts.empty()) {
        if (observed_counts.size() != pattern_roi.size()) {
          throw std::invalid_argument(
              "run.observed_counts size must match experiment.pattern_roi size.");
        }
        use_observed_data = true;
      }
    }

    const double theta_lo = run.value("theta_lo", 0.5);
    const double theta_hi = run.value("theta_hi", 10.0);

    const double constrain_prior_strength = run.value("constrain_prior_strength", 0.0);
    const bool constrain_use_gamma_sign = run.value("constrain_use_gamma_sign", false);
    const bool constrain_use_tau_weighted = run.value("constrain_use_tau_weighted", false);
    const int constrain_n_bins = run.value("constrain_n_bins", 1);

    const std::string profile_minimizer = run.value("profile_minimizer", std::string("pydme"));
    const bool pydme_mode = (profile_minimizer == "pydme");
    const bool pydme_style_ul = run.value("pydme_style_ul", pydme_mode);

    // Band-tool extension (optional). Two new keys, both with zero-impact defaults:
    //   run.mode (default "scan"): when "threshold_toys", the binary skips UL extraction
    //     and instead generates Poisson(s+b) sub-toys per mass at sigma_threshold(m_chi)
    //     and writes a TGraph q_target_per_mass to <outdir>/qtarget_threshold.root.
    //   run.q_target_lookup_path (default ""): when non-empty, open that ROOT file and
    //     use the TGraph "q_target_per_mass" as the per-mass PLR rejection threshold
    //     instead of the asymptotic constant (Phi^-1(cl))^2. Used by toy-MC band runs.
    const std::string run_mode = run.value("mode", std::string("scan"));
    if (run_mode != "scan" && run_mode != "threshold_toys") {
      throw std::runtime_error(
          "run.mode must be 'scan' or 'threshold_toys', got: " + run_mode);
    }
    const std::string q_target_lookup_path =
        run.value("q_target_lookup_path", std::string(""));

    const std::string sigma_grid_mode = run.value("sigma_grid_mode", std::string("intersection"));
    const bool use_union_sigma = (sigma_grid_mode == "union");
    const bool use_per_mass_sigma = (sigma_grid_mode == "per_mass");
    if (sigma_grid_mode != "intersection" && sigma_grid_mode != "union" && sigma_grid_mode != "per_mass") {
      throw std::runtime_error(
          "run.sigma_grid_mode must be 'intersection', 'union', or 'per_mass', got: " + sigma_grid_mode);
    }

    // Exposure for normalization bookkeeping (data is Asimov so not used directly here).
    const double livetime_days = exp.value("livetime_days", 0.0);
    const double duty_cycle = exp.value("duty_cycle", 1.0);
    const double mass_kg = det.value("mass_kg", 0.0);
    const double exposure_kg_year = (livetime_days * duty_cycle * mass_kg) / 365.25;

    // SRDM CSV signal values are rates in "gram/day" convention. Convert to expected
    // counts by multiplying by (exposure_kg_day * 1000 g/kg).
    const double exposure_kg_day_default = livetime_days * duty_cycle * mass_kg;
    const double signal_exposure_kg_day = sig.value("signal_exposure_kg_day", exposure_kg_day_default);
    const bool rates_are_per_gram_per_day = sig.value("rates_are_per_gram_per_day", true);

    // SRDM CSV inputs.
    const std::string input_dir = sig.at("input_dir").get<std::string>();
    const std::string filename_regex_str = sig.value(
        "filename_regex",
        std::string(R"(pattern_signal_summed_mX([0-9]+(?:\.[0-9]+)?)_full_QCD)"));
    const double signal_rate_scale = sig.value("signal_rate_scale", 1.0);
    const double signal_rate_scale_total =
        rates_are_per_gram_per_day ? signal_rate_scale * signal_exposure_kg_day * 1000.0
                                    : signal_rate_scale * signal_exposure_kg_day;

    const std::regex mchi_regex(filename_regex_str);
    if (!std::filesystem::exists(input_dir))
      throw std::runtime_error("SRDM input_dir does not exist: " + input_dir);

    // =========================================================================
    // 2. Discover and load per-mass signal CSV files
    // =========================================================================
    std::vector<std::pair<double, std::string>> csvs_by_mass;
    for (const auto& ent : std::filesystem::directory_iterator(input_dir)) {
      if (!ent.is_regular_file()) continue;
      const auto path = ent.path();
      const std::string fname = path.filename().string();

      auto parsed = parse_mchi_and_tag(fname, mchi_regex);
      if (!parsed) continue;

      const double mchi = parsed->first;
      csvs_by_mass.push_back({mchi, path.string()});
    }
    if (csvs_by_mass.empty()) {
      throw std::runtime_error("No SRDM CSVs found in " + input_dir +
                               " matching filename_regex=" + filename_regex_str);
    }

    // Sort by mchi and load each CSV.
    std::sort(csvs_by_mass.begin(), csvs_by_mass.end(),
              [](const auto& a, const auto& b) { return a.first < b.first; });

    std::vector<SrdmMassData> masses;
    masses.reserve(csvs_by_mass.size());
    for (const auto& [mchi, path] : csvs_by_mass) {
      if (verbosity >= 1) {
        std::cout << "[scan-srdm-csv] Loading mchi=" << std::scientific << mchi
                  << " from " << path << "\n";
      }
      masses.push_back(load_srdm_csv(path, mchi, pattern_roi, signal_rate_scale_total));
    }

    if (masses.empty()) throw std::runtime_error("No masses loaded.");

    const double rel_tol = 1e-12;
    // sigma_list_common: σ values used for the scan when mode is intersection or union (same grid for every mass).
    // sigma_y_axis: σ bin centers for TH2D q(mchi, σ); for per_mass it is the union of all masses' grids.
    std::vector<double> sigma_list_common;
    std::vector<double> sigma_y_axis;

    if (use_per_mass_sigma) {
      for (const auto& md : masses) {
        for (double s : md.sigmas) sigma_y_axis.push_back(s);
      }
      if (sigma_y_axis.empty()) throw std::runtime_error("Sigma union (for TH2D axis) is empty.");

      std::sort(sigma_y_axis.begin(), sigma_y_axis.end());
      std::vector<double> sigma_y_dedup;
      sigma_y_dedup.reserve(sigma_y_axis.size());
      for (double s : sigma_y_axis) {
        if (sigma_y_dedup.empty() || !approx_equal(s, sigma_y_dedup.back(), rel_tol)) {
          sigma_y_dedup.push_back(s);
        }
      }
      sigma_y_axis = std::move(sigma_y_dedup);
      if (sigma_y_axis.empty()) throw std::runtime_error("Sigma union axis is empty after deduplication.");
    } else if (use_union_sigma) {
      // Union of sigma grids across masses (interpolate missing per-mass templates).
      for (const auto& md : masses) {
        sigma_list_common.insert(sigma_list_common.end(), md.sigmas.begin(), md.sigmas.end());
      }
      if (sigma_list_common.empty()) throw std::runtime_error("Sigma union is empty.");

      std::sort(sigma_list_common.begin(), sigma_list_common.end());

      std::vector<double> sigma_list_dedup;
      sigma_list_dedup.reserve(sigma_list_common.size());
      for (double s : sigma_list_common) {
        if (sigma_list_dedup.empty() || !approx_equal(s, sigma_list_dedup.back(), rel_tol)) {
          sigma_list_dedup.push_back(s);
        }
      }
      sigma_list_common = std::move(sigma_list_dedup);

      if (sigma_list_common.empty()) throw std::runtime_error("Sigma union is empty after deduplication.");
      sigma_y_axis = sigma_list_common;
    } else {
      // Intersection: common σ grid for sensitivity / reproduction (exact CSV rows only).
      const auto& sig0 = masses.front().sigmas;
      sigma_list_common.reserve(sig0.size());
      for (double s : sig0) {
        bool ok = true;
        for (std::size_t im = 1; im < masses.size(); ++im) {
          if (find_sigma_index(masses[im].sigmas, s, rel_tol) < 0) {
            ok = false;
            break;
          }
        }
        if (ok) sigma_list_common.push_back(s);
      }
      if (sigma_list_common.empty()) throw std::runtime_error("Sigma intersection is empty.");
      std::sort(sigma_list_common.begin(), sigma_list_common.end());
      sigma_y_axis = sigma_list_common;
    }

    // mchi axis centers.
    std::vector<double> mchi_list;
    mchi_list.reserve(masses.size());
    for (const auto& md : masses) mchi_list.push_back(md.mchi);
    std::sort(mchi_list.begin(), mchi_list.end());

    const auto edges_m = make_edges_from_centers(mchi_list);
    const auto edges_s = make_edges_from_centers(sigma_y_axis);

    const int Nx = static_cast<int>(mchi_list.size());
    const int Ny = static_cast<int>(sigma_y_axis.size());
    const std::string q_title_data_tag = use_observed_data ? "observed" : "Asimov";
    TH2D hq("q_mchi_sigma_pattern",
            (";m_{#chi} [MeV];#sigma_{e} [cm^{2}] (q_{#mu} " + q_title_data_tag + ")").c_str(),
            Nx, edges_m.data(), Ny, edges_s.data());

    TH1D h_upper_limit("upper_limit_sigma_e_mchi",
                        ";m_{#chi} [MeV];#sigma_{e} [cm^{2}] (90% CL)",
                        Nx, edges_m.data());

    // =========================================================================
    // 3. Build background template and choose data vector
    // =========================================================================
    // Nominal total background template (Bp + Br at theta=1) — always saved as a
    // diagnostic reference, even when running an observed-data fit.
    std::vector<double> B_pat(pattern_roi.size(), 0.0);
    for (std::size_t i = 0; i < pattern_roi.size(); ++i) {
      B_pat[i] = background_Bp[i] + background_Br[i];
    }

    TH1D h_Btot("Btot",
                 ";pattern_roi (bin);B_{tot} = B_p + B_r (theta=1 nominal)",
                 static_cast<int>(pattern_roi.size()), 0.5,
                 static_cast<double>(pattern_roi.size()) + 0.5);
    for (int ib = 0; ib < static_cast<int>(pattern_roi.size()); ++ib) {
      h_Btot.SetBinContent(ib + 1, B_pat[static_cast<std::size_t>(ib)]);
      h_Btot.GetXaxis()->SetBinLabel(ib + 1, std::to_string(pattern_roi[ib]).c_str());
    }

    // Choose the data vector for the profile likelihood: observed counts when
    // provided (observed limit), else B_pat (Asimov sensitivity).
    const std::vector<double>& data_pat = use_observed_data ? observed_counts : B_pat;
    if (verbosity >= 1) {
      std::cout << "[scan-srdm-csv] data mode: " << (use_observed_data ? "observed" : "Asimov (Bp+Br)") << " | D = [";
      for (std::size_t i = 0; i < data_pat.size(); ++i) {
        std::cout << data_pat[i] << (i + 1 < data_pat.size() ? ", " : "");
      }
      std::cout << "]\n";
    }

    auto profile_pl = std::make_unique<ccdarksens::stats::ProfileLikelihood>();
    profile_pl->SetData(data_pat);
    profile_pl->SetBpBr(background_Bp, background_Br);
    profile_pl->SetConstrainPriorStrength(constrain_prior_strength);
    profile_pl->SetConstrainGammaSign(constrain_use_gamma_sign);
    profile_pl->SetConstrainUseTauWeighted(constrain_use_tau_weighted);
    profile_pl->SetConstrainNBins(constrain_n_bins);
    profile_pl->SetAcceptBoundary2DMinimum(pydme_style_ul);

    const bool do_2d_fit = (profile_minimizer == "pydme" || profile_minimizer == "minuit2d");
    if (verbosity >= 1) {
      std::cout << "[scan-srdm-csv] sigma_grid_mode=" << sigma_grid_mode;
      if (use_per_mass_sigma) {
        std::cout << "  TH2D_Ybins(union)=" << sigma_y_axis.size() << " (scan uses each mass's own σ grid)\n";
      } else {
        std::cout << "  N_sigma=" << sigma_y_axis.size() << ", masses: " << mchi_list.size()
                  << ", do_2d_fit=" << (do_2d_fit ? "true" : "false") << "\n";
      }
    }

    const double target_q_asymptotic = std::pow(TMath::NormQuantile(cl), 2);
    // target_q_asymptotic is the asymptotic (chi^2 / Wilks) PLR rejection threshold for
    // the configured CL. When run.q_target_lookup_path is non-empty, target_q_for_mass()
    // overrides it on a per-mass basis using the TGraph from the band tool's threshold
    // phase; otherwise the constant is returned for every mass.

    std::map<double, double> q_target_by_mass;
    if (!q_target_lookup_path.empty()) {
      TFile fql(q_target_lookup_path.c_str(), "READ");
      if (!fql.IsOpen() || fql.IsZombie()) {
        throw std::runtime_error(
            "Cannot open q_target_lookup_path: " + q_target_lookup_path);
      }
      TGraph* g_qt = dynamic_cast<TGraph*>(fql.Get("q_target_per_mass"));
      if (!g_qt) {
        throw std::runtime_error(
            "q_target_lookup_path missing TGraph 'q_target_per_mass': " + q_target_lookup_path);
      }
      const int Nq = g_qt->GetN();
      for (int i = 0; i < Nq; ++i) {
        double x = 0.0, y = 0.0;
        g_qt->GetPoint(i, x, y);
        q_target_by_mass[x] = y;
      }
      fql.Close();
      if (verbosity >= 1) {
        std::cout << "[scan-srdm-csv] Loaded q_target lookup with " << Nq
                  << " entries from " << q_target_lookup_path << "\n";
      }
    }

    auto target_q_for_mass = [&](double mchi) -> double {
      if (q_target_by_mass.empty()) return target_q_asymptotic;
      // Exact (within mass-grid tolerance) match preferred; fall back to nearest neighbor.
      const double rel_tol_mass = 1e-9;
      for (const auto& kv : q_target_by_mass) {
        if (std::abs(kv.first - mchi) <=
            rel_tol_mass * std::max(1.0, std::abs(mchi))) {
          return kv.second;
        }
      }
      double best_q = target_q_asymptotic;
      double best_diff = std::numeric_limits<double>::infinity();
      for (const auto& kv : q_target_by_mass) {
        const double d = std::abs(kv.first - mchi);
        if (d < best_diff) {
          best_diff = d;
          best_q = kv.second;
        }
      }
      return best_q;
    };
    // =========================================================================
    // 4. Initialise output histograms and per-mass grid structures
    // =========================================================================
    // Fill q histogram.
    const bool have_multiple_masses = masses.size() > 1;
    (void)have_multiple_masses;

    std::unique_ptr<TH1D> h_spat_dump;
    bool spat_dumped = false;
    if (dump_spat_point) {
      const int np = static_cast<int>(pattern_roi.size());
      h_spat_dump = std::make_unique<TH1D>(
          "Spat_dump",
          ";pattern_roi (bin);S_{pat} from SRDM CSV (expected counts, Asimov template)",
          np, 0.5, np + 0.5);
      for (int ib = 0; ib < np; ++ib) {
        h_spat_dump->GetXaxis()->SetBinLabel(ib + 1,
                                             std::to_string(pattern_roi[ib]).c_str());
      }
    }

    // Map mchi -> mass data index for stable access even if csv order differs from mchi_list order.
    std::map<double, std::size_t> mass_index_by_value;
    for (std::size_t i = 0; i < masses.size(); ++i) mass_index_by_value[masses[i].mchi] = i;

    // =========================================================================
    // 5. Threshold-toys mode (band tool Phase 1 — early exit)
    //    run.mode == "threshold_toys": generate per-mass q_target lookup.
    //
    // For each mass:
    //   sigma_thr = sigma_threshold_per_mass(m_chi)   // from Phase 0 graph
    //   S_thr     = Spat at (m_chi, sigma_thr)        // CSV / log10 interpolation
    //   lambda    = S_thr + Bp + theta_nominal * Br   // theta_nominal = 1
    //   for k = 1..n_threshold_toys:
    //       n^{(k)} ~ Poisson(lambda)                 // independent RNG stream
    //       SetData(n^{(k)})
    //       nll_top  = MinimizeOverScaleMinuit(S_thr, theta_lo, theta_hi).second
    //       nll_glob = MinimizeOverSigmaAndTheta(...,
    //                      log10(sigma_thr), theta_hat_top).nll_min     // seeded 2D fit
    //       q^{(k)}  = max(0, 2 * (nll_top - nll_glob))
    //   c_mu(m_chi) = quantile(q[], threshold_percentile)
    //
    // Important: the 2D Simplex is *seeded* at (log10(σ_thr), θ̂_top) rather
    // than the box midpoint. Because σ_thr is by construction the asymptotic
    // 90% UL — which sits high in the σ scan range — a midpoint-seeded Simplex
    // routinely fails to find the global minimum (especially with the boundary
    // rejection guard). The seeded fit is fast (one Minuit run per toy) and
    // robust. SetAcceptBoundary2DMinimum(pydme_style_ul) ensures we don't
    // discard valid boundary minima.
    //
    // Reuses the *exact* ProfileLikelihood instance and minimizers used in scan
    // mode — no new likelihood code, only a different driver loop.
    // -------------------------------------------------------------------------
    if (run_mode == "threshold_toys") {
      const auto thr = run.at("threshold_toys");
      const std::string sigma_thr_path =
          thr.at("sigma_threshold_graph_path").get<std::string>();
      const long long n_threshold_toys =
          thr.value("n_threshold_toys", static_cast<long long>(10000));
      const double threshold_percentile = thr.value("percentile", 0.90);
      const std::uint64_t threshold_rng_seed =
          thr.value("rng_seed", static_cast<std::uint64_t>(23456ULL));

      if (n_threshold_toys < 100) {
        throw std::runtime_error(
            "threshold_toys.n_threshold_toys must be >= 100 (got " +
            std::to_string(n_threshold_toys) + ").");
      }
      if (threshold_percentile <= 0.0 || threshold_percentile >= 1.0) {
        throw std::runtime_error(
            "threshold_toys.percentile must be in (0,1), got " +
            std::to_string(threshold_percentile));
      }

      std::map<double, double> sigma_thr_by_mass;
      {
        TFile fst(sigma_thr_path.c_str(), "READ");
        if (!fst.IsOpen() || fst.IsZombie()) {
          throw std::runtime_error(
              "Cannot open sigma_threshold_graph_path: " + sigma_thr_path);
        }
        TGraph* g_st = dynamic_cast<TGraph*>(fst.Get("sigma_threshold_per_mass"));
        if (!g_st) {
          throw std::runtime_error(
              "sigma_threshold_graph_path missing TGraph 'sigma_threshold_per_mass': " +
              sigma_thr_path);
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
          if (std::abs(kv.first - mchi) <=
              rel_tol_mass * std::max(1.0, std::abs(mchi))) {
            return kv.second;
          }
        }
        double best = -1.0;
        double best_diff = std::numeric_limits<double>::infinity();
        for (const auto& kv : sigma_thr_by_mass) {
          const double d = std::abs(kv.first - mchi);
          if (d < best_diff) {
            best_diff = d;
            best = kv.second;
          }
        }
        return best;
      };

      if (verbosity >= 1) {
        std::cout << "[scan-srdm-csv] mode=threshold_toys n_threshold_toys="
                  << n_threshold_toys << " percentile=" << threshold_percentile
                  << " seed=" << threshold_rng_seed << "\n";
      }

      std::vector<double> qt_mchi;
      std::vector<double> qt_value;
      qt_mchi.reserve(mchi_list.size());
      qt_value.reserve(mchi_list.size());

      for (int ix = 0; ix < Nx; ++ix) {
        const double mchi = mchi_list[ix];
        const auto it_im = mass_index_by_value.find(mchi);
        if (it_im == mass_index_by_value.end()) continue;
        const std::size_t im = it_im->second;

        const double sigma_thr = sigma_thr_for_mass(mchi);
        if (!(sigma_thr > 0.0)) {
          throw std::runtime_error(
              "sigma_threshold for mchi=" + std::to_string(mchi) +
              " is non-positive in " + sigma_thr_path);
        }

        // Spat at sigma_thr: log10-interpolated from this mass's CSV templates.
        const std::vector<double> S_thr =
            interpolate_spat_on_sigma(masses[im], sigma_thr, rel_tol);
        if (S_thr.size() != pattern_roi.size()) {
          throw std::runtime_error(
              "Internal error: S_thr size mismatch for mchi=" + std::to_string(mchi));
        }

        // lambda_i = S_thr_i + Bp_i + theta_nominal * Br_i (theta_nominal = 1)
        std::vector<double> lambda(pattern_roi.size(), 0.0);
        for (std::size_t i = 0; i < lambda.size(); ++i) {
          lambda[i] = S_thr[i] + background_Bp[i] + background_Br[i];
          if (lambda[i] < 0.0) {
            throw std::runtime_error(
                "Negative Poisson mean encountered in threshold_toys.");
          }
        }

        // S_from_log10_sigma callable used by the global 2D minimizer.
        const double log10_lo = std::log10(masses[im].sigmas.front());
        const double log10_hi = std::log10(masses[im].sigmas.back());
        const double log10_sigma_thr_seed = std::log10(sigma_thr);
        auto S_from_log10 = [&masses, im, rel_tol](double log10_s) -> std::vector<double> {
          return interpolate_spat_on_sigma(masses[im], std::pow(10.0, log10_s), rel_tol);
        };

        // Independent RNG stream per mass for reproducibility under different
        // mass orderings: derive a stream key from (seed, ix).
        std::seed_seq seq{static_cast<std::uint32_t>(threshold_rng_seed & 0xFFFFFFFFu),
                          static_cast<std::uint32_t>((threshold_rng_seed >> 32) & 0xFFFFFFFFu),
                          static_cast<std::uint32_t>(ix),
                          0xC0FFEEu};
        std::mt19937_64 rng(seq);

        std::vector<double> q_values(static_cast<std::size_t>(n_threshold_toys), 0.0);
        std::vector<std::poisson_distribution<long long>> poisson_per_bin;
        poisson_per_bin.reserve(lambda.size());
        for (double mu : lambda) poisson_per_bin.emplace_back(mu);

        for (long long k = 0; k < n_threshold_toys; ++k) {
          std::vector<double> nk(pattern_roi.size(), 0.0);
          for (std::size_t i = 0; i < nk.size(); ++i) {
            nk[i] = static_cast<double>(poisson_per_bin[i](rng));
          }

          profile_pl->SetData(nk);

          // Constrained NLL: profile θ at fixed σ = σ_thr.
          const auto pr_top =
              profile_pl->MinimizeOverScaleMinuit(S_thr, theta_lo, theta_hi);
          const double theta_hat_top = pr_top.first;
          const double nll_top = pr_top.second;

          // Unconstrained NLL: 2D Simplex seeded at (log10(σ_thr), θ̂_top),
          // which lies very near the global min for toys generated at σ_thr.
          double nll_glob = nll_top;
          double dbg_log10_sig_hat = log10_sigma_thr_seed;
          double dbg_theta_hat = theta_hat_top;
          bool dbg_ok = false;
          if (do_2d_fit && masses[im].sigmas.size() >= 2u) {
            auto res = profile_pl->MinimizeOverSigmaAndTheta(
                S_from_log10, log10_lo, log10_hi, theta_lo, theta_hi,
                log10_sigma_thr_seed, theta_hat_top);
            if (res.ok && res.nll_min < nll_glob) {
              nll_glob = res.nll_min;
              dbg_log10_sig_hat = res.log10_sigma_hat;
              dbg_theta_hat = res.theta_hat;
              dbg_ok = true;
            }
          }

          const double q = 2.0 * (nll_top - nll_glob);
          q_values[static_cast<std::size_t>(k)] = (q > 0.0) ? q : 0.0;

          if (verbosity >= 2 && k < 5) {
            std::cout << "[scan-srdm-csv]   thr_toy mchi=" << std::scientific << mchi
                      << " k=" << k
                      << " nll_top=" << std::setprecision(8) << nll_top
                      << " theta_top=" << theta_hat_top
                      << " nll_glob=" << nll_glob
                      << " q=" << q
                      << " log10_sig_hat=" << dbg_log10_sig_hat
                      << " theta_hat=" << dbg_theta_hat
                      << " 2d_ok=" << dbg_ok
                      << "\n";
          }
        }

        const std::size_t q_idx =
            static_cast<std::size_t>(threshold_percentile *
                                     static_cast<double>(n_threshold_toys - 1));
        std::nth_element(q_values.begin(), q_values.begin() + q_idx, q_values.end());
        const double c_mu = q_values[q_idx];

        qt_mchi.push_back(mchi);
        qt_value.push_back(c_mu);

        if (verbosity >= 1) {
          std::cout << "[scan-srdm-csv] threshold_toys mchi=" << std::scientific << mchi
                    << " sigma_thr=" << sigma_thr
                    << " q_target(p=" << threshold_percentile << ")="
                    << std::setprecision(6) << c_mu << "\n";
        }
      }

      std::filesystem::create_directories(outdir);
      const std::string thr_out_path = outdir + "/qtarget_threshold.root";
      TFile fout_thr(thr_out_path.c_str(), "RECREATE");
      if (!fout_thr.IsOpen()) {
        throw std::runtime_error("Failed to create output file: " + thr_out_path);
      }
      TGraph g_qt(static_cast<int>(qt_mchi.size()));
      g_qt.SetName("q_target_per_mass");
      g_qt.SetTitle(";m_{#chi} [MeV];q_{#mu, target} (toy-MC threshold)");
      for (int i = 0; i < static_cast<int>(qt_mchi.size()); ++i) {
        g_qt.SetPoint(i, qt_mchi[static_cast<std::size_t>(i)],
                      qt_value[static_cast<std::size_t>(i)]);
      }
      g_qt.Write("q_target_per_mass");
      TParameter<double>("threshold_percentile", threshold_percentile).Write();
      TParameter<double>("n_threshold_toys",
                         static_cast<double>(n_threshold_toys)).Write();
      fout_thr.Close();

      if (verbosity >= 1) {
        std::cout << "[scan-srdm-csv] Wrote " << thr_out_path << "\n";
      }
      return 0;
    }

    double t_br_mchi = 0.0;
    double t_br_sigma = 0.0;
    double t_br_q_mono = 0.0;
    double t_br_nll = 0.0;
    double t_br_theta_hat = 0.0;
    std::unique_ptr<TTree> tree_q_scan_per_mass;
    if (use_per_mass_sigma) {
      tree_q_scan_per_mass = std::make_unique<TTree>(
          "q_scan_per_mass",
          "Per-mass xsec grid: mchi, sigma, q_mono, nll, theta_hat (Asimov PLR at CSV points)");
      tree_q_scan_per_mass->Branch("mchi_MeV", &t_br_mchi);
      tree_q_scan_per_mass->Branch("sigma_e_cm2", &t_br_sigma);
      tree_q_scan_per_mass->Branch("q_mono", &t_br_q_mono);
      tree_q_scan_per_mass->Branch("nll", &t_br_nll);
      tree_q_scan_per_mass->Branch("theta_hat", &t_br_theta_hat);
    }

    // =========================================================================
    // 6. Main grid scan — compute q(mχ, σ) and extract 90% CL upper limits
    // =========================================================================
    for (int ix = 0; ix < Nx; ++ix) {
      const double mchi = mchi_list[ix];
      const auto it = mass_index_by_value.find(mchi);
      if (it == mass_index_by_value.end()) continue;
      const std::size_t im = it->second;

      const std::vector<double>& sig_scan =
          use_per_mass_sigma ? masses[im].sigmas : sigma_list_common;
      if (sig_scan.empty()) throw std::runtime_error("Empty sigma scan grid for a mass.");

      // Build S_grid on sig_scan (union uses interpolation onto common sig_scan points).
      std::vector<std::vector<double>> S_grid;
      S_grid.reserve(sig_scan.size());
      for (double s : sig_scan) {
        if (use_union_sigma && !use_per_mass_sigma) {
          S_grid.push_back(interpolate_spat_on_sigma(masses[im], s, rel_tol));
        } else {
          const int idx_s = find_sigma_index(masses[im].sigmas, s, rel_tol);
          if (idx_s < 0) {
            throw std::runtime_error("Sigma missing on mass template grid (unexpected).");
          }
          S_grid.push_back(masses[im].Spat[static_cast<std::size_t>(idx_s)]);
        }
      }

      // nll(sigma) and nll_min.
      std::vector<double> nll_values;
      nll_values.reserve(sig_scan.size());
      std::vector<double> theta_hat_values;
      theta_hat_values.reserve(sig_scan.size());

      double nll_min_grid = std::numeric_limits<double>::infinity();
      double nll_min_2d = nll_min_grid;
      bool use_nll_min_2d = false;
      double log10_sigma_hat_pydme = std::log10(sig_scan.front());  // for pydme UL bracketing

      if (do_2d_fit && sig_scan.size() >= 2u) {
        const double log10_lo = std::log10(sig_scan.front());
        const double log10_hi = std::log10(sig_scan.back());
        auto S_from_log10 = [&S_grid, &sig_scan, log10_lo, log10_hi](double log10_s) -> std::vector<double> {
          const std::size_t n = S_grid.empty() ? 0u : S_grid[0].size();
          if (S_grid.size() < 2u) return S_grid.empty() ? std::vector<double>(n, 0.0) : S_grid[0];

          const double log10_s_clamp = std::max(log10_lo, std::min(log10_hi, log10_s));
          for (std::size_t j = 0; j + 1 < sig_scan.size(); ++j) {
            const double l0 = std::log10(sig_scan[j]);
            const double l1 = std::log10(sig_scan[j + 1]);
            if (log10_s_clamp >= l0 && log10_s_clamp <= l1) {
              const double t = (l1 - l0) > 1e-300 ? (log10_s_clamp - l0) / (l1 - l0) : 0.0;
              std::vector<double> out(n, 0.0);
              for (std::size_t b = 0; b < n; ++b) out[b] = S_grid[j][b] + t * (S_grid[j + 1][b] - S_grid[j][b]);
              return out;
            }
          }
          if (log10_s_clamp <= std::log10(sig_scan.front())) return S_grid.front();
          return S_grid.back();
        };

        auto res = profile_pl->MinimizeOverSigmaAndTheta(
            S_from_log10, log10_lo, log10_hi, theta_lo, theta_hi);
        if (res.ok) {
          nll_min_2d = res.nll_min;
          use_nll_min_2d = true;
          log10_sigma_hat_pydme = res.log10_sigma_hat;
          if (verbosity >= 2) {
            std::cout << "[scan-srdm-csv] mchi=" << mchi << " 2D-fit nll_min=" << res.nll_min
                      << " log10_sigma_hat=" << res.log10_sigma_hat << " theta_hat=" << res.theta_hat << "\n";
          }
        }
      }

      // Compute nll(sigma) at discrete sigma points (used for q entries).
      for (std::size_t k = 0; k < sig_scan.size(); ++k) {
        const auto pr = profile_pl->MinimizeOverScaleMinuit(S_grid[k], theta_lo, theta_hi);
        const double theta_hat = pr.first;
        const double nll = pr.second;
        nll_values.push_back(nll);
        theta_hat_values.push_back(theta_hat);

        if (dump_spat_point && !spat_dumped &&
            approx_equal(mchi, dump_spat_point_mchi_MeV) &&
            approx_equal(sig_scan[k], dump_spat_point_sigma_e_cm2)) {
          if (!h_spat_dump) throw std::runtime_error("Internal error: dump_spat_point enabled but h_spat_dump=nullptr");
          const auto& spat = S_grid[k];
          if (spat.size() != pattern_roi.size()) {
            throw std::runtime_error("Internal error: Spat dump size mismatch");
          }
          for (std::size_t b = 0; b < spat.size(); ++b) {
            h_spat_dump->SetBinContent(static_cast<int>(b) + 1, spat[b]);
          }
          spat_dumped = true;
        }

        nll_min_grid = std::min(nll_min_grid, nll);
      }

      const bool use_2d_nll_min = (pydme_style_ul && use_nll_min_2d);
      const double nll_min_final = use_2d_nll_min ? nll_min_2d : nll_min_grid;
      if (!std::isfinite(nll_min_final)) throw std::runtime_error("Non-finite nll_min encountered.");

      // Monotonicize q(σ): running max so UL crossing is well-defined and reduces sawtooth
      std::vector<double> q_mu(nll_values.size()), q_mono(nll_values.size());
      for (std::size_t k = 0; k < nll_values.size(); ++k) {
        const double q = 2.0 * (nll_values[k] - nll_min_final);
        q_mu[k] = (q > 0.0) ? q : 0.0;
      }
      q_mono[0] = q_mu[0];
      for (std::size_t k = 1; k < nll_values.size(); ++k)
        q_mono[k] = std::max(q_mono[k - 1], q_mu[k]);

      // Compute upper limit σ_UL (smallest sigma with q_mu >= target_q).
      // target_q is per-mass when q_target_lookup_path was provided (toy-MC threshold);
      // otherwise it falls back to the asymptotic constant for every mass.
      const double target_q = target_q_for_mass(mchi);
      double ul_sigma = sig_scan.back();
      const double log10_lo = std::log10(sig_scan.front());
      const double log10_hi = std::log10(sig_scan.back());
      if (pydme_mode && S_grid.size() >= 2u) {
        auto S_from_log10 = [&S_grid, &sig_scan, log10_lo, log10_hi](double log10_s) -> std::vector<double> {
          const std::size_t n = S_grid.empty() ? 0u : S_grid[0].size();
          std::vector<double> out(n, 0.0);
          if (S_grid.size() < 2u) return S_grid.empty() ? out : S_grid[0];
          const double log10_s_clamp = std::max(log10_lo, std::min(log10_hi, log10_s));
          for (std::size_t j = 0; j + 1 < sig_scan.size(); ++j) {
            const double l0 = std::log10(sig_scan[j]);
            const double l1 = std::log10(sig_scan[j + 1]);
            if (log10_s_clamp >= l0 && log10_s_clamp <= l1) {
              const double t = (l1 - l0) > 1e-300 ? (log10_s_clamp - l0) / (l1 - l0) : 0.0;
              for (std::size_t b = 0; b < n; ++b)
                out[b] = S_grid[j][b] + t * (S_grid[j + 1][b] - S_grid[j][b]);
              return out;
            }
          }
          if (log10_s_clamp <= std::log10(sig_scan.front())) return S_grid.front();
          return S_grid.back();
        };

        auto q_mu_at = [&](double log10_sig) {
          std::vector<double> S = S_from_log10(log10_sig);
          const double nll = profile_pl->MinimizeOverScaleMinuit(S, theta_lo, theta_hi).second;
          const double q = 2.0 * (nll - nll_min_final);
          return q > 0.0 ? q : 0.0;
        };

        const int k_best_pydme =
            static_cast<int>(std::min_element(nll_values.begin(), nll_values.end()) - nll_values.begin());
        const double log10_sigma_best = std::log10(sig_scan[static_cast<std::size_t>(k_best_pydme)]);
        const double lo = use_2d_nll_min ? std::max({log10_sigma_hat_pydme, log10_lo, log10_sigma_best})
                                         : std::max(log10_lo, log10_sigma_best);
        double hi = lo;

        const double ul_brack_step = 0.3;
        const int ul_max_expand = 12;
        const int ul_max_iter = 24;
        const double ul_q_tol = 0.01;
        const double tol_x = std::max(1e-3, 0.01 * (log10_hi - log10_lo));
        double step = ul_brack_step;
        int brack_count = 0;
        for (int _ = 0; _ < ul_max_expand && hi < log10_hi - 1e-12; ++_) {
          brack_count++;
          hi = std::min(hi + step, log10_hi);
          if (q_mu_at(hi) >= target_q) break;
          step *= 2.0;
        }

        if (q_mu_at(hi) >= target_q) {
          double left = lo, right = hi;
          for (int it = 0; it < ul_max_iter; ++it) {
            const double mid = 0.5 * (left + right);
            const double q_mid = q_mu_at(mid);
            if (q_mid >= target_q) right = mid;
            else left = mid;
            if (std::abs(right - left) < tol_x || std::abs(q_mid - target_q) < ul_q_tol) break;
          }
          ul_sigma = std::pow(10.0, right);
          if (verbosity >= 1) {
            std::cout << "[scan-srdm-csv] mchi=" << std::scientific << mchi << " MeV pydme UL: log10(sigma)="
                      << std::log10(ul_sigma) << " (bracketing steps=" << brack_count << ", bisection)\n";
          }
        } else {
          if (verbosity >= 1) {
            const double q_at_hi = q_mu_at(log10_hi);
            std::cout << "[scan-srdm-csv] mchi=" << std::scientific << mchi
                      << " MeV pydme UL: bracketing did not reach target_q (q at log10_sigma_hi="
                      << log10_hi << " is " << q_at_hi << " < " << target_q << "); using grid.\n";
          }
          const int k_best =
              static_cast<int>(std::min_element(nll_values.begin(), nll_values.end()) - nll_values.begin());
          for (int k = k_best; k < static_cast<int>(nll_values.size()); ++k) {
            if (q_mono[k] >= target_q) {
              ul_sigma = sig_scan[static_cast<std::size_t>(k)];
              break;
            }
          }
        }
      } else {
        const int k_best =
            static_cast<int>(std::min_element(nll_values.begin(), nll_values.end()) - nll_values.begin());
        for (int k = k_best; k < static_cast<int>(nll_values.size()); ++k) {
          if (q_mono[k] >= target_q) {
            ul_sigma = sig_scan[static_cast<std::size_t>(k)];
            break;
          }
        }
      }

      for (std::size_t k = 0; k < sig_scan.size(); ++k) {
        int iy_bin = static_cast<int>(k);
        if (use_per_mass_sigma) {
          iy_bin = find_sigma_index(sigma_y_axis, sig_scan[k], rel_tol);
          if (iy_bin < 0) {
            throw std::runtime_error("per_mass: sigma not found on union TH2D Y axis (unexpected).");
          }
        }
        hq.SetBinContent(ix + 1, iy_bin + 1, q_mono[k]);

        if (tree_q_scan_per_mass) {
          t_br_mchi = mchi;
          t_br_sigma = sig_scan[k];
          t_br_q_mono = q_mono[k];
          t_br_nll = nll_values[k];
          t_br_theta_hat = theta_hat_values[k];
          tree_q_scan_per_mass->Fill();
        }

        if (verbosity >= 1) {
          std::cout << "[scan-srdm-csv] mchi=" << std::scientific << mchi << " MeV"
                    << " sigma=" << std::scientific << sig_scan[k] << " cm^2"
                    << " theta_hat=" << std::setprecision(6) << theta_hat_values[k]
                    << " nll=" << std::setprecision(8) << nll_values[k]
                    << " q_mono=" << std::setprecision(8) << q_mono[k] << "\n";
        }
      }
      h_upper_limit.SetBinContent(ix + 1, ul_sigma);
    }

    // =========================================================================
    // 7. Post-processing and ROOT output
    // =========================================================================
    // Optional: remove upward spikes from the limit curve (envelope smoothing).
    {
      const int n_mchi_lim = h_upper_limit.GetNbinsX();
      std::vector<double> sigma_lim(n_mchi_lim);
      for (int i = 0; i < n_mchi_lim; ++i) sigma_lim[i] = h_upper_limit.GetBinContent(i + 1);
      for (int pass = 0; pass < 10; ++pass) {
        bool changed = false;
        for (int i = 1; i < n_mchi_lim - 1; ++i) {
          const double s_prev = sigma_lim[i - 1];
          double& s_curr = sigma_lim[i];
          const double s_next = sigma_lim[i + 1];
          if (s_curr <= 0.0 || s_prev <= 0.0 || s_next <= 0.0) continue;
          if (s_curr > s_prev && s_curr > s_next) {
            const double cap = std::max(s_prev, s_next);
            if (s_curr > cap) {
              s_curr = cap;
              changed = true;
            }
          }
        }
        if (!changed) break;
      }
      for (int i = 0; i < n_mchi_lim; ++i) h_upper_limit.SetBinContent(i + 1, sigma_lim[i]);
    }

    // Build a TGraph mirror of the upper-limit curve at the *exact* mchi values.
    // The TH1D bin centers don't equal the original masses for log-spaced grids
    // (because make_edges_from_centers only places bin centers at the inputs
    //  when the inputs are evenly spaced), so the TH1D alone is unsuitable for
    //  consumers that read X positions back. Downstream tools should prefer
    //  this TGraph when present.
    TGraph g_upper_limit(static_cast<int>(mchi_list.size()));
    g_upper_limit.SetName("upper_limit_sigma_e_mchi_graph");
    g_upper_limit.SetTitle(";m_{#chi} [MeV];#sigma_{e} [cm^{2}] (90% CL)");
    for (int i = 0; i < static_cast<int>(mchi_list.size()); ++i) {
      g_upper_limit.SetPoint(i, mchi_list[static_cast<std::size_t>(i)],
                             h_upper_limit.GetBinContent(i + 1));
    }

    // Write output.
    std::filesystem::create_directories(outdir);
    const std::string out_path = outdir + "/scan_srdm_pattern_csv.root";
    TFile fout(out_path.c_str(), "RECREATE");
    if (!fout.IsOpen()) {
      throw std::runtime_error("Failed to create output file: " + out_path);
    }
    hq.Write("q_mchi_sigma_pattern");
    h_upper_limit.Write("upper_limit_sigma_e_mchi");
    g_upper_limit.Write("upper_limit_sigma_e_mchi_graph");
    h_Btot.Write("Btot");

    {
      const std::string data_title = use_observed_data
                                         ? ";pattern_roi (bin);D_{obs} (observed counts)"
                                         : ";pattern_roi (bin);D = B_{p} + B_{r} (Asimov, theta=1)";
      TH1D h_data("Data",
                  data_title.c_str(),
                  static_cast<int>(pattern_roi.size()), 0.5,
                  static_cast<double>(pattern_roi.size()) + 0.5);
      for (int ib = 0; ib < static_cast<int>(pattern_roi.size()); ++ib) {
        h_data.SetBinContent(ib + 1, data_pat[static_cast<std::size_t>(ib)]);
        h_data.GetXaxis()->SetBinLabel(ib + 1, std::to_string(pattern_roi[ib]).c_str());
      }
      h_data.Write("Data");
    }

    if (tree_q_scan_per_mass) tree_q_scan_per_mass->Write();
    if (h_spat_dump) h_spat_dump->Write();

    fout.cd();
    TParameter<double>("exposure_kg_year", exposure_kg_year).Write();
    TParameter<double>("signal_exposure_kg_day", signal_exposure_kg_day).Write();
    TParameter<double>("rates_are_per_gram_per_day", rates_are_per_gram_per_day ? 1.0 : 0.0).Write();
    TParameter<double>("signal_rate_scale_total", signal_rate_scale_total).Write();
    TParameter<double>("signal_rate_scale", signal_rate_scale).Write();

    fout.Close();
    if (verbosity >= 1) {
      std::cout << "[scan-srdm-csv] Wrote " << out_path << "\n";
    }
  } catch (const std::exception& e) {
    std::cerr << "[scan-srdm-csv] ERROR: " << e.what() << "\n";
    return 1;
  }

  return 0;
}

