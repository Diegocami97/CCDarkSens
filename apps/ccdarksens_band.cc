// ============================================================================
//  CCDarkSens — ccdarksens_band
//  Orchestrates repeated subprocess calls to a configured scan binary to build median and ±1σ/±2σ expected sensitivity bands from Poisson background-only toys.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#endif

#include <TDirectory.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TParameter.h>
#include <TError.h>

#include <nlohmann/json.hpp>

using nlohmann::json;
namespace fs = std::filesystem;

namespace {

// -----------------------------------------------------------------------------
// JSON helpers
// -----------------------------------------------------------------------------

// Resolve a dotted-path JSON key like "run.background_Bp" into a const ref.
// Throws if the path doesn't exist or the leaf isn't an array of numbers.
static const json& resolve_dotted(const json& root, const std::string& dotted) {
  std::vector<std::string> parts;
  {
    std::stringstream ss(dotted);
    std::string token;
    while (std::getline(ss, token, '.')) parts.push_back(token);
  }
  const json* cur = &root;
  for (const auto& p : parts) {
    if (!cur->is_object() || !cur->contains(p)) {
      throw std::runtime_error("Dotted key not found: " + dotted);
    }
    cur = &(cur->at(p));
  }
  return *cur;
}

static std::vector<double> resolve_dotted_double_vec(const json& root,
                                                     const std::string& dotted) {
  const json& leaf = resolve_dotted(root, dotted);
  if (!leaf.is_array()) {
    throw std::runtime_error("Dotted key is not an array: " + dotted);
  }
  std::vector<double> out;
  out.reserve(leaf.size());
  for (const auto& v : leaf) out.push_back(v.get<double>());
  return out;
}

// Write a JSON object to a file (pretty-printed for human inspection on failure).
static void write_json(const fs::path& p, const json& j) {
  std::ofstream out(p);
  if (!out.is_open()) {
    throw std::runtime_error("Cannot open for write: " + p.string());
  }
  out << std::setw(2) << j << "\n";
}

// -----------------------------------------------------------------------------
// ROOT helpers (graph harvesting + writing)
// -----------------------------------------------------------------------------

struct UlCurve {
  std::vector<double> mchi;
  std::vector<double> sigma_ul;
};

// Open a ROOT file and read a TGraph by name. Returns its (x,y) points.
//
// Restores gDirectory after closing the read file. Otherwise destroying the
// TFile resets gDirectory to gROOT and a caller that still has an output TFile
// open will see Write() fail with "current directory is not associated with a
// file" (e.g. embedding q_target_per_mass into band.root).
static UlCurve read_tgraph(const fs::path& root_path,
                           const std::string& graph_name) {
  // Suppress ROOT's noisy "no such file" / "no such object" Info/Warning
  // messages for failed reads; we want clean status logging from this binary.
  const Int_t prev_level = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kError;

  TDirectory* const dir_before = gDirectory;
  UlCurve out;
  try {
    std::unique_ptr<TFile> f(TFile::Open(root_path.c_str(), "READ"));
    gErrorIgnoreLevel = prev_level;

    if (!f || f->IsZombie()) {
      throw std::runtime_error("Cannot open ROOT file: " + root_path.string());
    }
    TGraph* g = dynamic_cast<TGraph*>(f->Get(graph_name.c_str()));
    if (!g) {
      throw std::runtime_error("ROOT file " + root_path.string() +
                               " missing TGraph '" + graph_name + "'");
    }
    const int N = g->GetN();
    out.mchi.reserve(N);
    out.sigma_ul.reserve(N);
    for (int i = 0; i < N; ++i) {
      double x = 0.0, y = 0.0;
      g->GetPoint(i, x, y);
      out.mchi.push_back(x);
      out.sigma_ul.push_back(y);
    }
  } catch (...) {
    gErrorIgnoreLevel = prev_level;
    if (dir_before) {
      dir_before->cd();
    }
    throw;
  }
  if (dir_before) {
    dir_before->cd();
  }
  return out;
}

static TGraph make_tgraph(const std::string& name, const std::string& title,
                          const std::vector<double>& xs,
                          const std::vector<double>& ys) {
  TGraph g(static_cast<int>(xs.size()));
  g.SetName(name.c_str());
  g.SetTitle(title.c_str());
  for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
    g.SetPoint(i, xs[static_cast<std::size_t>(i)], ys[static_cast<std::size_t>(i)]);
  }
  return g;
}

// -----------------------------------------------------------------------------
// Subprocess invocation
// -----------------------------------------------------------------------------

struct ScanInvocationResult {
  int exit_code = 0;
  std::string log_path;
};

// Synchronously invoke `<scan_binary> <cfg_path> > <log_path> 2>&1` via
// std::system. Returns the exit code and the log path so we can show it on
// failure. The user's PATH is honored.
static ScanInvocationResult invoke_scan(const std::string& scan_binary,
                                        const fs::path& cfg_path,
                                        const fs::path& log_path) {
  // Quote both args: simple shell-safe path quoting (paths under tmp_dir or the
  // user's working dir; no shell metachars expected in normal use).
  auto sh_quote = [](const std::string& s) -> std::string {
    std::string out;
    out.reserve(s.size() + 2);
    out += "'";
    for (char c : s) {
      if (c == '\'') out += "'\\''";
      else out += c;
    }
    out += "'";
    return out;
  };

  const std::string cmd = sh_quote(scan_binary) + " " + sh_quote(cfg_path.string()) +
                          " > " + sh_quote(log_path.string()) + " 2>&1";
  ScanInvocationResult r;
  r.exit_code = std::system(cmd.c_str());
#if defined(WIFEXITED) && defined(WEXITSTATUS)
  // std::system on POSIX returns the wait status; pull out the actual exit
  // code if normal exit, otherwise propagate as nonzero.
  if (WIFEXITED(r.exit_code)) {
    r.exit_code = WEXITSTATUS(r.exit_code);
  }
#endif
  r.log_path = log_path.string();
  return r;
}

// -----------------------------------------------------------------------------
// Quantile extraction
// -----------------------------------------------------------------------------

static double quantile_inplace(std::vector<double>& v, double q) {
  if (v.empty()) return std::nan("");
  if (q <= 0.0) return *std::min_element(v.begin(), v.end());
  if (q >= 1.0) return *std::max_element(v.begin(), v.end());
  const std::size_t n = v.size();
  // Linear interpolation between adjacent ranks (numpy "linear" interpolation).
  const double pos = q * static_cast<double>(n - 1);
  const std::size_t lo = static_cast<std::size_t>(std::floor(pos));
  const std::size_t hi = static_cast<std::size_t>(std::ceil(pos));
  std::nth_element(v.begin(), v.begin() + lo, v.end());
  const double y_lo = v[lo];
  if (lo == hi) return y_lo;
  // After nth_element, elements after `lo` are >= v[lo]; need partial sort up
  // to `hi` to safely take v[hi].
  std::nth_element(v.begin() + lo + 1, v.begin() + hi, v.end());
  const double y_hi = v[hi];
  const double frac = pos - static_cast<double>(lo);
  return y_lo + frac * (y_hi - y_lo);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json [--stop-after phase0|phase1|phase2]\n"
              << "  Config must contain a 'band' block with at least:\n"
              << "    scan_binary (string), n_toys (int), outdir (string)\n"
              << "  --stop-after phase0  run Phase 0 only (Asimov σ_threshold; needs toy_mc|both).\n"
              << "  --stop-after phase1  run Phases 0+1 only (writes qtarget_threshold.root).\n"
              << "  --stop-after phase2  full band (default; same as omitting the flag).\n";
    return 1;
  }

  const fs::path config_path = argv[1];
  std::string stop_after;  // empty or "phase2" = run outer toy loop; "phase0"/"phase1" = early exit
  for (int ai = 2; ai < argc; ++ai) {
    const std::string arg(argv[ai]);
    if (arg == "--stop-after" && ai + 1 < argc) {
      stop_after = argv[++ai];
      for (auto& c : stop_after) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
      }
      continue;
    }
    std::cerr << "[band] ERROR: Unknown CLI argument: " << arg
              << " (expected optional --stop-after <phase0|phase1|phase2>)\n";
    return 1;
  }
  if (stop_after.empty()) {
    stop_after = "phase2";
  }
  if (stop_after != "phase0" && stop_after != "phase1" && stop_after != "phase2") {
    std::cerr << "[band] ERROR: --stop-after must be phase0, phase1, or phase2 (got: "
              << stop_after << ")\n";
    return 1;
  }

  const auto t_start = std::chrono::steady_clock::now();

  try {
    json cfg_user;
    {
      std::ifstream cf(config_path);
      if (!cf.is_open()) {
        throw std::runtime_error("Cannot open config: " + config_path.string());
      }
      cf >> cfg_user;
    }

    if (!cfg_user.contains("band") || !cfg_user.at("band").is_object()) {
      throw std::runtime_error(
          "Config has no 'band' block. Add it (see docs/Plan_sensitivity_band_app.md §2).");
    }
    const json band = cfg_user.at("band");

    // ----------------------- band block parsing ------------------------------
    const std::string scan_binary =
        band.at("scan_binary").get<std::string>();
    const long long n_toys = band.at("n_toys").get<long long>();
    if (n_toys < 1) throw std::runtime_error("band.n_toys must be >= 1.");
    const fs::path band_outdir = fs::path(band.at("outdir").get<std::string>());

    const std::uint64_t rng_seed =
        band.value("rng_seed", static_cast<std::uint64_t>(12345ULL));
    int n_workers = band.value("n_workers", 0);
    if (n_workers <= 0) {
      n_workers = static_cast<int>(std::thread::hardware_concurrency());
      if (n_workers <= 0) n_workers = 1;
    }
    n_workers = std::min<int>(n_workers, static_cast<int>(n_toys));

    const std::string ul_graph_name =
        band.value("ul_graph_name", std::string("upper_limit_sigma_e_mchi_graph"));
    const std::string bp_key =
        band.value("background_bp_key", std::string("run.background_Bp"));
    const std::string br_key =
        band.value("background_br_key", std::string("run.background_Br"));

    const bool save_per_toy_curves = band.value("save_per_toy_curves", false);
    const bool keep_per_toy_outputs = band.value("keep_per_toy_outputs", false);
    const bool abort_on_failure = band.value("abort_on_failure", false);
    const double min_success_fraction = band.value("min_success_fraction", 0.95);

    const std::string threshold_mode = band.value("threshold", std::string("asymptotic"));
    if (threshold_mode != "asymptotic" && threshold_mode != "toy_mc" &&
        threshold_mode != "both") {
      throw std::runtime_error(
          "band.threshold must be 'asymptotic' | 'toy_mc' | 'both', got: " + threshold_mode);
    }
    const bool do_asy = (threshold_mode == "asymptotic" || threshold_mode == "both");
    const bool do_toy = (threshold_mode == "toy_mc" || threshold_mode == "both");

    if ((stop_after == "phase0" || stop_after == "phase1") && !do_toy) {
      throw std::runtime_error(
          "--stop-after " + stop_after +
          " requires band.threshold 'toy_mc' or 'both' (got '" + threshold_mode + "').");
    }

    const long long n_threshold_toys =
        band.value("n_threshold_toys", static_cast<long long>(10000));
    const double threshold_percentile = band.value("threshold_percentile", 0.90);
    const std::uint64_t threshold_seed =
        band.value("threshold_seed", static_cast<std::uint64_t>(23456ULL));

    // Default tmp_dir = <outdir>/_toys
    const fs::path tmp_dir = band.contains("tmp_dir")
                                 ? fs::path(band.at("tmp_dir").get<std::string>())
                                 : (band_outdir / "_toys");

    // Quantile config (defaults: 0.025 / 0.16 / 0.50 / 0.84 / 0.975).
    double q_med = 0.50, q_l1 = 0.16, q_h1 = 0.84, q_l2 = 0.025, q_h2 = 0.975;
    if (band.contains("quantiles") && band.at("quantiles").is_object()) {
      const auto& qb = band.at("quantiles");
      q_med = qb.value("median", q_med);
      q_l1  = qb.value("low_1sigma", q_l1);
      q_h1  = qb.value("high_1sigma", q_h1);
      q_l2  = qb.value("low_2sigma", q_l2);
      q_h2  = qb.value("high_2sigma", q_h2);
    }

    fs::create_directories(band_outdir);
    fs::create_directories(tmp_dir);

    // ----------------------- resolve sampling means --------------------------
    const std::vector<double> Bp = resolve_dotted_double_vec(cfg_user, bp_key);
    const std::vector<double> Br = resolve_dotted_double_vec(cfg_user, br_key);
    if (Bp.size() != Br.size() || Bp.empty()) {
      throw std::runtime_error("Bp / Br vectors empty or size mismatch.");
    }
    std::vector<double> lambda(Bp.size(), 0.0);
    for (std::size_t i = 0; i < Bp.size(); ++i) {
      if (Bp[i] < 0.0 || Br[i] < 0.0) {
        throw std::runtime_error("Negative Bp / Br entry encountered.");
      }
      lambda[i] = Bp[i] + Br[i];
    }

    std::cout << "[band] config         : " << config_path << "\n"
              << "[band] scan_binary    : " << scan_binary << "\n"
              << "[band] threshold      : " << threshold_mode << "\n"
              << "[band] n_toys         : " << n_toys << "\n"
              << "[band] n_workers      : " << n_workers << "\n"
              << "[band] outdir         : " << band_outdir << "\n"
              << "[band] tmp_dir        : " << tmp_dir << "\n"
              << "[band] N_pattern (Bp) : " << Bp.size() << "\n"
              << "[band] rng_seed       : " << rng_seed << "\n";
    if (do_toy) {
      std::cout << "[band] threshold_seed : " << threshold_seed << "\n"
                << "[band] n_threshold    : " << n_threshold_toys << "\n"
                << "[band] threshold_pct  : " << threshold_percentile << "\n";
    }

    // -------------------------------------------------------------------------
    // Phase 0 — Asimov asymptotic UL (only for toy_mc / both)
    // -------------------------------------------------------------------------
    fs::path sigma_threshold_root;
    UlCurve sigma_threshold_curve;
    if (do_toy) {
      const fs::path phase0_dir = tmp_dir / "phase0";
      fs::create_directories(phase0_dir);

      json cfg_phase0 = cfg_user;
      cfg_phase0.erase("band");
      cfg_phase0["run"].erase("observed_counts");
      cfg_phase0["run"].erase("mode");
      cfg_phase0["run"].erase("q_target_lookup_path");
      cfg_phase0["run"].erase("threshold_toys");
      cfg_phase0["run"]["outdir"] = phase0_dir.string();
      // Increase verbosity floor of 1 so we get the limit lines in the log.
      if (!cfg_phase0["run"].contains("verbosity")) {
        cfg_phase0["run"]["verbosity"] = 1;
      }

      const fs::path phase0_cfg = phase0_dir / "config.json";
      const fs::path phase0_log = phase0_dir / "scan.log";
      write_json(phase0_cfg, cfg_phase0);

      std::cout << "[band] Phase 0 — Asimov asymptotic UL → " << phase0_dir << "\n";
      const ScanInvocationResult r0 = invoke_scan(scan_binary, phase0_cfg, phase0_log);
      if (r0.exit_code != 0) {
        throw std::runtime_error("Phase 0 scan failed (exit=" +
                                 std::to_string(r0.exit_code) +
                                 "). See log: " + phase0_log.string());
      }

      // The scan's output ROOT is named after its analysis (e.g.
      // scan_srdm_pattern_csv.root). Find the unique *.root in phase0_dir
      // containing ul_graph_name.
      fs::path scan_out_root;
      for (const auto& entry : fs::directory_iterator(phase0_dir)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".root") continue;
        try {
          (void)read_tgraph(entry.path(), ul_graph_name);
          scan_out_root = entry.path();
          break;
        } catch (...) {
          // not the right ROOT; keep looking
        }
      }
      if (scan_out_root.empty()) {
        throw std::runtime_error(
            "Phase 0: no ROOT in " + phase0_dir.string() +
            " contains TGraph '" + ul_graph_name +
            "'. Check scan log: " + phase0_log.string());
      }
      sigma_threshold_curve = read_tgraph(scan_out_root, ul_graph_name);
      std::cout << "[band]    σ_threshold loaded from " << scan_out_root
                << " (N_mass=" << sigma_threshold_curve.mchi.size() << ")\n";

      // Save as sigma_threshold_per_mass.
      sigma_threshold_root = tmp_dir / "sigma_threshold.root";
      {
        TFile fst(sigma_threshold_root.c_str(), "RECREATE");
        TGraph g_st = make_tgraph(
            "sigma_threshold_per_mass",
            ";m_{#chi} [MeV];#sigma_{threshold} (Asimov asymptotic UL) [cm^{2}]",
            sigma_threshold_curve.mchi, sigma_threshold_curve.sigma_ul);
        g_st.Write("sigma_threshold_per_mass");
        fst.Close();
      }
      std::cout << "[band]    wrote " << sigma_threshold_root << "\n";

      if (stop_after == "phase0") {
        const auto t_end = std::chrono::steady_clock::now();
        const double secs =
            std::chrono::duration<double>(t_end - t_start).count();
        std::cout << "[band] --stop-after phase0: finished in " << std::fixed
                  << std::setprecision(1) << secs << " s\n";
        return 0;
      }
    }

    // -------------------------------------------------------------------------
    // Phase 1 — q_target threshold toys (only for toy_mc / both)
    // -------------------------------------------------------------------------
    fs::path qtarget_root;
    if (do_toy) {
      const fs::path phase1_dir = tmp_dir / "phase1";
      fs::create_directories(phase1_dir);

      json cfg_phase1 = cfg_user;
      cfg_phase1.erase("band");
      cfg_phase1["run"].erase("observed_counts");
      cfg_phase1["run"]["mode"] = "threshold_toys";
      cfg_phase1["run"]["outdir"] = phase1_dir.string();
      cfg_phase1["run"]["threshold_toys"] = {
          {"sigma_threshold_graph_path", sigma_threshold_root.string()},
          {"n_threshold_toys", n_threshold_toys},
          {"percentile", threshold_percentile},
          {"rng_seed", threshold_seed}};
      if (!cfg_phase1["run"].contains("verbosity")) {
        cfg_phase1["run"]["verbosity"] = 1;
      }

      const fs::path phase1_cfg = phase1_dir / "config.json";
      const fs::path phase1_log = phase1_dir / "scan.log";
      write_json(phase1_cfg, cfg_phase1);

      std::cout << "[band] Phase 1 — q_target threshold toys → " << phase1_dir
                << "  (N=" << n_threshold_toys << " sub-toys per mass)\n";
      const ScanInvocationResult r1 = invoke_scan(scan_binary, phase1_cfg, phase1_log);
      if (r1.exit_code != 0) {
        throw std::runtime_error("Phase 1 threshold_toys failed (exit=" +
                                 std::to_string(r1.exit_code) +
                                 "). See log: " + phase1_log.string());
      }
      qtarget_root = phase1_dir / "qtarget_threshold.root";
      if (!fs::exists(qtarget_root)) {
        throw std::runtime_error(
            "Phase 1 finished but " + qtarget_root.string() + " was not produced.");
      }
      const UlCurve qt = read_tgraph(qtarget_root, "q_target_per_mass");
      std::cout << "[band]    q_target_per_mass loaded (N=" << qt.mchi.size()
                << "). Range: [" << *std::min_element(qt.sigma_ul.begin(), qt.sigma_ul.end())
                << ", "
                << *std::max_element(qt.sigma_ul.begin(), qt.sigma_ul.end()) << "]\n";

      if (stop_after == "phase1") {
        const auto t_end = std::chrono::steady_clock::now();
        const double secs =
            std::chrono::duration<double>(t_end - t_start).count();
        std::cout << "[band] --stop-after phase1: finished in " << std::fixed
                  << std::setprecision(1) << secs << " s\n";
        std::cout << "[band]    q_target lookup: " << qtarget_root << "\n";
        return 0;
      }
    }

    // -------------------------------------------------------------------------
    // Phase 2 — outer toy loop (always)
    //
    // For each toy t = 0..n_toys-1:
    //   1. Seed mt19937_64 with (rng_seed XOR t-mixing) for reproducible
    //      independent streams.
    //   2. Sample D^(t)_i ~ Poisson(lambda_i).
    //   3. Build per-toy config: stripped of band block, with observed_counts
    //      = D^(t), outdir = tmp_dir/toy_<t>/{asy,toy}.
    //   4. Invoke the scan binary; harvest σ_UL curve from
    //      upper_limit_sigma_e_mchi_graph in its output ROOT.
    //
    // For threshold="both", each toy runs the scan binary twice on the same
    // D^(t): once without q_target_lookup_path (asymptotic), once with
    // q_target_lookup_path attached (toy_mc).
    // -------------------------------------------------------------------------
    std::cout << "[band] Phase 2 — outer toy loop  (K=" << n_toys
              << ", workers=" << n_workers << ")\n";

    // Per-toy results, indexed by (mode_idx, toy_idx) where mode_idx 0=asy, 1=toy.
    // Each entry is a UlCurve. Empty = failed toy.
    struct ToyResult {
      bool ok_asy = false;
      bool ok_toy = false;
      UlCurve curve_asy;
      UlCurve curve_toy;
      std::string fail_reason;
    };
    std::vector<ToyResult> results(static_cast<std::size_t>(n_toys));

    auto run_one_toy = [&](long long t) {
      const fs::path toy_dir = tmp_dir / ("toy_" + std::to_string(t));
      fs::create_directories(toy_dir);

      // Independent stream per toy: combine seed with toy index via seed_seq.
      std::seed_seq seq{static_cast<std::uint32_t>(rng_seed & 0xFFFFFFFFu),
                        static_cast<std::uint32_t>((rng_seed >> 32) & 0xFFFFFFFFu),
                        static_cast<std::uint32_t>(t & 0xFFFFFFFFu),
                        static_cast<std::uint32_t>((t >> 32) & 0xFFFFFFFFu),
                        0xBADCAFEu};
      std::mt19937_64 rng(seq);

      std::vector<double> D_toy(lambda.size(), 0.0);
      for (std::size_t i = 0; i < lambda.size(); ++i) {
        std::poisson_distribution<long long> pd(lambda[i]);
        D_toy[i] = static_cast<double>(pd(rng));
      }

      // Common per-toy config: strip band, mode/threshold-toys; inject toy data.
      json cfg_base = cfg_user;
      cfg_base.erase("band");
      cfg_base["run"].erase("mode");
      cfg_base["run"].erase("threshold_toys");
      cfg_base["run"].erase("q_target_lookup_path");
      // observed_counts becomes the toy realization for this iteration.
      cfg_base["run"]["observed_counts"] = D_toy;
      // Quiet down per-toy logs (still go to scan.log).
      cfg_base["run"]["verbosity"] = 0;

      auto run_pass = [&](const std::string& tag,
                          bool attach_q_target_lookup,
                          UlCurve& out_curve,
                          bool& out_ok,
                          std::string& out_fail) -> void {
        const fs::path pass_dir = toy_dir / tag;
        fs::create_directories(pass_dir);

        json cfg_pass = cfg_base;
        cfg_pass["run"]["outdir"] = pass_dir.string();
        if (attach_q_target_lookup) {
          cfg_pass["run"]["q_target_lookup_path"] = qtarget_root.string();
        }

        const fs::path pass_cfg = pass_dir / "config.json";
        const fs::path pass_log = pass_dir / "scan.log";
        try {
          write_json(pass_cfg, cfg_pass);
        } catch (const std::exception& ex) {
          out_ok = false;
          out_fail = std::string("write_json: ") + ex.what();
          return;
        }
        const ScanInvocationResult rp = invoke_scan(scan_binary, pass_cfg, pass_log);
        if (rp.exit_code != 0) {
          out_ok = false;
          out_fail = "scan exit=" + std::to_string(rp.exit_code) +
                     " log=" + pass_log.string();
          return;
        }

        // Find the scan's output ROOT in pass_dir (one *.root containing UL graph).
        fs::path scan_out;
        for (const auto& entry : fs::directory_iterator(pass_dir)) {
          if (!entry.is_regular_file()) continue;
          if (entry.path().extension() != ".root") continue;
          try {
            out_curve = read_tgraph(entry.path(), ul_graph_name);
            scan_out = entry.path();
            break;
          } catch (...) {
          }
        }
        if (scan_out.empty()) {
          out_ok = false;
          out_fail = "no ROOT containing '" + ul_graph_name + "' in " +
                     pass_dir.string();
          return;
        }
        out_ok = true;
      };

      ToyResult r;
      if (do_asy) {
        run_pass("asymptotic", /*attach_q_target_lookup=*/false,
                 r.curve_asy, r.ok_asy, r.fail_reason);
      } else {
        r.ok_asy = false;
      }
      if (do_toy) {
        std::string fr;
        run_pass("toy_mc", /*attach_q_target_lookup=*/true,
                 r.curve_toy, r.ok_toy, fr);
        if (!fr.empty()) r.fail_reason = (r.fail_reason.empty() ? fr : r.fail_reason + "; " + fr);
      } else {
        r.ok_toy = false;
      }
      results[static_cast<std::size_t>(t)] = std::move(r);

      if (!keep_per_toy_outputs) {
        std::error_code ec;
        fs::remove_all(toy_dir, ec);
      }
    };

    // Thread pool: a simple atomic-counter work-stealing scheme.
    std::atomic<long long> next_toy{0};
    std::mutex log_mu;
    std::atomic<long long> n_done{0};
    std::atomic<long long> n_failed{0};
    auto worker = [&]() {
      while (true) {
        const long long t = next_toy.fetch_add(1);
        if (t >= n_toys) return;
        try {
          run_one_toy(t);
        } catch (const std::exception& ex) {
          ToyResult r;
          r.ok_asy = false;
          r.ok_toy = false;
          r.fail_reason = std::string("exception: ") + ex.what();
          results[static_cast<std::size_t>(t)] = std::move(r);
        }
        const long long done = n_done.fetch_add(1) + 1;
        const auto& res = results[static_cast<std::size_t>(t)];
        const bool any_pass_ok = (do_asy ? res.ok_asy : true) &&
                                 (do_toy ? res.ok_toy : true);
        if (!any_pass_ok) n_failed.fetch_add(1);
        if (done % std::max<long long>(1, n_toys / 20) == 0 || done == n_toys) {
          std::lock_guard<std::mutex> lk(log_mu);
          std::cout << "[band]    progress " << done << " / " << n_toys
                    << " (failed=" << n_failed.load() << ")\n"
                    << std::flush;
        }
      }
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(n_workers));
    for (int w = 0; w < n_workers; ++w) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    long long n_ok_asy = 0, n_ok_toy = 0;
    for (const auto& r : results) {
      if (r.ok_asy) ++n_ok_asy;
      if (r.ok_toy) ++n_ok_toy;
    }
    std::cout << "[band]    completed: ok_asy=" << n_ok_asy
              << " / " << (do_asy ? n_toys : 0)
              << "  ok_toy=" << n_ok_toy
              << " / " << (do_toy ? n_toys : 0) << "\n";

    if (do_asy) {
      const double frac = static_cast<double>(n_ok_asy) / static_cast<double>(n_toys);
      if (frac < min_success_fraction) {
        if (abort_on_failure) {
          throw std::runtime_error(
              "Asymptotic survival fraction " + std::to_string(frac) +
              " < min_success_fraction=" + std::to_string(min_success_fraction));
        } else {
          std::cerr << "[band] WARNING: asymptotic survival " << frac
                    << " < " << min_success_fraction << "\n";
        }
      }
    }
    if (do_toy) {
      const double frac = static_cast<double>(n_ok_toy) / static_cast<double>(n_toys);
      if (frac < min_success_fraction) {
        if (abort_on_failure) {
          throw std::runtime_error(
              "Toy-MC survival fraction " + std::to_string(frac) +
              " < min_success_fraction=" + std::to_string(min_success_fraction));
        } else {
          std::cerr << "[band] WARNING: toy_mc survival " << frac
                    << " < " << min_success_fraction << "\n";
        }
      }
    }

    // -------------------------------------------------------------------------
    // Quantiles per mass
    //
    // Build a unified mass grid as the union of all surviving toys' mchi
    // points. Each per-toy curve is sampled at the master grid via nearest-
    // mass match; the contract requires every per-toy curve uses the same
    // mchi grid (they all come from the same scan binary on the same config),
    // so this match is always exact in the supported case.
    // -------------------------------------------------------------------------
    auto pivot_and_quantile = [&](bool use_asy) -> std::vector<TGraph> {
      std::vector<TGraph> out;
      std::vector<double> master_mchi;
      bool master_set = false;

      // Find master grid from the first surviving toy.
      for (const auto& r : results) {
        const bool ok = use_asy ? r.ok_asy : r.ok_toy;
        if (!ok) continue;
        const auto& c = use_asy ? r.curve_asy : r.curve_toy;
        master_mchi = c.mchi;
        master_set = true;
        break;
      }
      if (!master_set) return out;

      const std::size_t M = master_mchi.size();
      // per-mass vector of σ_UL across surviving toys
      std::vector<std::vector<double>> per_mass(M);
      for (auto& v : per_mass) v.reserve(static_cast<std::size_t>(n_toys));

      // Verify and accumulate.
      auto match_idx = [&](const std::vector<double>& mchi_curve,
                           std::size_t i_master) -> int {
        const double m = master_mchi[i_master];
        const double tol = 1e-9 * std::max(1.0, std::abs(m));
        for (std::size_t j = 0; j < mchi_curve.size(); ++j) {
          if (std::abs(mchi_curve[j] - m) <= tol) return static_cast<int>(j);
        }
        return -1;
      };

      for (const auto& r : results) {
        const bool ok = use_asy ? r.ok_asy : r.ok_toy;
        if (!ok) continue;
        const auto& c = use_asy ? r.curve_asy : r.curve_toy;
        for (std::size_t i = 0; i < M; ++i) {
          const int j = match_idx(c.mchi, i);
          if (j < 0) continue;
          const double s = c.sigma_ul[static_cast<std::size_t>(j)];
          if (s > 0.0 && std::isfinite(s)) per_mass[i].push_back(s);
        }
      }

      std::vector<double> y_med(M), y_l1(M), y_h1(M), y_l2(M), y_h2(M);
      for (std::size_t i = 0; i < M; ++i) {
        std::vector<double>& v = per_mass[i];
        if (v.empty()) {
          y_med[i] = y_l1[i] = y_h1[i] = y_l2[i] = y_h2[i] = 0.0;
          continue;
        }
        // We need each quantile, but quantile_inplace does partial sorting
        // that mutates the vector — copy once per quantile.
        auto qv = [&](double q) -> double {
          std::vector<double> tmp(v);
          return quantile_inplace(tmp, q);
        };
        y_med[i] = qv(q_med);
        y_l1[i]  = qv(q_l1);
        y_h1[i]  = qv(q_h1);
        y_l2[i]  = qv(q_l2);
        y_h2[i]  = qv(q_h2);
      }

      const std::string suf = use_asy ? "_asymptotic" : "_toy_mc";
      out.push_back(make_tgraph(
          "median_expected_sigma_e_mchi" + suf,
          ";m_{#chi} [MeV];median expected #sigma_{UL} [cm^{2}]",
          master_mchi, y_med));
      out.push_back(make_tgraph(
          "band_1sigma_low_sigma_e_mchi" + suf,
          ";m_{#chi} [MeV];#sigma_{UL} 16 % quantile [cm^{2}]",
          master_mchi, y_l1));
      out.push_back(make_tgraph(
          "band_1sigma_high_sigma_e_mchi" + suf,
          ";m_{#chi} [MeV];#sigma_{UL} 84 % quantile [cm^{2}]",
          master_mchi, y_h1));
      out.push_back(make_tgraph(
          "band_2sigma_low_sigma_e_mchi" + suf,
          ";m_{#chi} [MeV];#sigma_{UL} 2.5 % quantile [cm^{2}]",
          master_mchi, y_l2));
      out.push_back(make_tgraph(
          "band_2sigma_high_sigma_e_mchi" + suf,
          ";m_{#chi} [MeV];#sigma_{UL} 97.5 % quantile [cm^{2}]",
          master_mchi, y_h2));
      return out;
    };

    std::vector<TGraph> graphs_asy = do_asy ? pivot_and_quantile(true) : std::vector<TGraph>{};
    std::vector<TGraph> graphs_toy = do_toy ? pivot_and_quantile(false) : std::vector<TGraph>{};

    // -------------------------------------------------------------------------
    // Write band.root
    // -------------------------------------------------------------------------
    const fs::path band_root = band_outdir / "band.root";
    {
      TFile fout(band_root.c_str(), "RECREATE");
      if (!fout.IsOpen()) {
        throw std::runtime_error("Cannot create band ROOT: " + band_root.string());
      }
      for (auto& g : graphs_asy) g.Write();
      for (auto& g : graphs_toy) g.Write();

      // Phase artifacts as side-graphs (helpful for plotting / debugging).
      if (!sigma_threshold_curve.mchi.empty()) {
        TGraph g_st = make_tgraph(
            "sigma_threshold_per_mass",
            ";m_{#chi} [MeV];#sigma_{threshold} (Asimov asymptotic UL) [cm^{2}]",
            sigma_threshold_curve.mchi, sigma_threshold_curve.sigma_ul);
        g_st.Write("sigma_threshold_per_mass");
      }
      if (!qtarget_root.empty() && fs::exists(qtarget_root)) {
        const UlCurve qt = read_tgraph(qtarget_root, "q_target_per_mass");
        TGraph g_qt = make_tgraph(
            "q_target_per_mass",
            ";m_{#chi} [MeV];q_{#mu, target} (toy-MC threshold)",
            qt.mchi, qt.sigma_ul);
        g_qt.Write("q_target_per_mass");
      }

      // Provenance metadata.
      TParameter<double>("n_toys", static_cast<double>(n_toys)).Write();
      TParameter<double>("rng_seed", static_cast<double>(rng_seed)).Write();
      TParameter<double>("n_workers", static_cast<double>(n_workers)).Write();
      TParameter<double>("n_ok_asy", static_cast<double>(n_ok_asy)).Write();
      TParameter<double>("n_ok_toy", static_cast<double>(n_ok_toy)).Write();
      TParameter<double>("threshold_asymptotic_run",
                         do_asy ? 1.0 : 0.0).Write();
      TParameter<double>("threshold_toy_mc_run",
                         do_toy ? 1.0 : 0.0).Write();
      if (do_toy) {
        TParameter<double>("n_threshold_toys",
                           static_cast<double>(n_threshold_toys)).Write();
        TParameter<double>("threshold_percentile", threshold_percentile).Write();
        TParameter<double>("threshold_seed",
                           static_cast<double>(threshold_seed)).Write();
      }

      // Save per-toy curves if requested (TTree with branches indexed by toy).
      if (save_per_toy_curves) {
        long long t_idx = 0;
        double t_mchi = 0.0;
        double t_asy = 0.0;
        double t_toy = 0.0;
        TTree band_per_toy("band_per_toy", "Per-toy σ_UL points");
        band_per_toy.Branch("toy_idx", &t_idx);
        band_per_toy.Branch("mchi_MeV", &t_mchi);
        band_per_toy.Branch("sigma_UL_asy", &t_asy);
        band_per_toy.Branch("sigma_UL_toy", &t_toy);
        for (std::size_t t = 0; t < results.size(); ++t) {
          t_idx = static_cast<long long>(t);
          const auto& r = results[t];
          // Use whichever curve has the master grid for mass listing.
          const auto& m_grid = r.ok_asy ? r.curve_asy.mchi : r.curve_toy.mchi;
          for (std::size_t i = 0; i < m_grid.size(); ++i) {
            t_mchi = m_grid[i];
            t_asy = (r.ok_asy && i < r.curve_asy.sigma_ul.size())
                        ? r.curve_asy.sigma_ul[i] : 0.0;
            t_toy = (r.ok_toy && i < r.curve_toy.sigma_ul.size())
                        ? r.curve_toy.sigma_ul[i] : 0.0;
            band_per_toy.Fill();
          }
        }
        band_per_toy.Write();
      }

      fout.Close();
    }

    // Optional cleanup of the toy scratch dir (Phase 0 / Phase 1 artifacts kept
    // so that downstream tools can re-read sigma_threshold / q_target if they
    // want without redoing the scan).
    if (!keep_per_toy_outputs) {
      // Per-toy dirs already removed inside run_one_toy; remove only empty
      // children of tmp_dir to keep Phase 0 / Phase 1 around.
      std::error_code ec;
      for (const auto& entry : fs::directory_iterator(tmp_dir, ec)) {
        if (ec) break;
        if (entry.is_directory() && entry.path().filename().string().rfind("toy_", 0) == 0) {
          std::error_code rc;
          fs::remove_all(entry.path(), rc);
        }
      }
    }

    const auto t_end = std::chrono::steady_clock::now();
    const double secs =
        std::chrono::duration<double>(t_end - t_start).count();
    std::cout << "[band] DONE in " << std::fixed << std::setprecision(1) << secs
              << " s   →   " << band_root << "\n";
  } catch (const std::exception& ex) {
    std::cerr << "[band] ERROR: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
