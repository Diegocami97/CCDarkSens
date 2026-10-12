// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternTables.cc -- Generation, reading, writing and comparison of the
//  signal efficiency table P(pattern | n_e) and the background identification
//  matrix B[identified | injected]; see PatternTables.hh.
// ===========================================================================

#include "ccdarksens/response/PatternTables.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "ccdarksens/response/ChargeTransport.hh"

namespace ccdarksens {

namespace {

// ----------------------------------------------------------------------------
// PatternString
//   Pattern digits to a string, e.g. {2,1,1} -> "211"; used to sort the rows of
//   the background matrix the way the reference table is sorted.
// ----------------------------------------------------------------------------
std::string PatternString(const std::vector<int>& p) {
  std::string s;
  for (int d : p) s += std::to_string(d);
  return s;
}

// ----------------------------------------------------------------------------
// WarnIfDifferent
//   Print a warning when a config value differs from the one the shipped
//   reference background table was validated with.
// ----------------------------------------------------------------------------
void WarnIfDifferent(const char* name, double got, double want) {
  if (std::abs(got - want) > 1e-9)
    std::cerr << "[bg-eff] WARNING: " << name << " = " << got
              << " differs from the validated value " << want << "\n";
}

} // namespace

// ----------------------------------------------------------------------------
// PatternCode
//   Pattern digits to integer code, e.g. {2,1,1} -> 211.
// ----------------------------------------------------------------------------
int PatternCode(const std::vector<int>& p) {
  int code = 0;
  for (int d : p) code = code * 10 + d;
  return code;
}

// ----------------------------------------------------------------------------
// EnumerateIdealPatterns
//   The empty pattern (0) first, then every one-, two- and three-pixel pattern
//   with 1..5 e- per pixel and total charge <= max_sum, ordered by length and
//   then lexicographically (notebook pattern_simulation_background).
// ----------------------------------------------------------------------------
std::vector<std::vector<int>> EnumerateIdealPatterns(int max_sum) {
  std::vector<std::vector<int>> out;
  out.push_back({0});
  std::set<std::vector<int>> uniq;
  for (int a = 1; a <= 5 && a <= max_sum; ++a) uniq.insert({a});
  for (int a = 1; a <= 5; ++a)
    for (int b = 1; b <= 5; ++b)
      if (a + b <= max_sum) uniq.insert({a, b});
  for (int a = 1; a <= 5; ++a)
    for (int b = 1; b <= 5; ++b)
      for (int c = 1; c <= 5; ++c)
        if (a + b + c <= max_sum) uniq.insert({a, b, c});
  std::vector<std::vector<int>> rest(uniq.begin(), uniq.end());
  std::sort(rest.begin(), rest.end(), [](const std::vector<int>& x, const std::vector<int>& y) {
    if (x.size() != y.size()) return x.size() < y.size();
    return x < y;
  });
  for (const auto& p : rest) out.push_back(p);
  return out;
}

// ----------------------------------------------------------------------------
// MakePatternTableSetup
//   Build the charge transport, the image generator and the classifier from the
//   response blocks of the config. The signal table uses the 3 x 50 binned image;
//   the background matrix uses the 3 x 5 cluster geometry and the classifier
//   flags of the reference background table.
// ----------------------------------------------------------------------------
PatternTableSetup MakePatternTableSetup(const ConfigManager& cfg, bool background_mode) {
  const auto& det = cfg.detector();
  const auto& emj = cfg.response().emc;
  const auto& pimg = cfg.response().pattern_image;
  const auto& pcc_j = cfg.response().pattern_classifier;

  ChargeTransportConfig ct_cfg;
  ct_cfg.thickness_um = det.geometry().thickness_mm * 1000.0;
  ct_cfg.A_um2 = emj.A_um2;
  ct_cfg.b_umInv = emj.b_umInv;
  ct_cfg.alpha = emj.alpha;
  ct_cfg.beta_per_keV = emj.beta_per_keV;
  ct_cfg.rng_seed = emj.rng_seed;
  auto ct = std::make_shared<ChargeTransport>(ct_cfg);

  PatternImageConfig img_cfg;
  if (background_mode) {
    const bool from_json = (pimg.nrows_binned > 0 && pimg.ncols > 0 && pimg.row_binning > 0);
    const bool from_detector = !from_json && (det.geometry().rows > 0 && det.geometry().cols > 0 &&
                                              pimg.row_binning > 0 && pimg.col_binning > 0);
    if (from_json) {
      img_cfg.nrows_binned = pimg.nrows_binned;
      img_cfg.ncols = pimg.ncols;
    } else if (from_detector) {
      img_cfg.raw_rows = det.geometry().rows;
      img_cfg.raw_cols = det.geometry().cols;
    } else {
      img_cfg.nrows_binned = pimg.nrows_binned > 0 ? pimg.nrows_binned : 3;
      img_cfg.ncols = pimg.ncols > 0 ? pimg.ncols : 50;
    }
    img_cfg.row_binning = pimg.row_binning > 0 ? pimg.row_binning : 100;
    img_cfg.col_binning = pimg.col_binning > 0 ? pimg.col_binning : 1;
    img_cfg.pixel_size_um = pimg.pixel_size_um;
    img_cfg.sigma_readout_e = pimg.sigma_readout_e;
    img_cfg.lambda_dc = pimg.lambda_dc;
    img_cfg.rng_seed = pimg.rng_seed;
  } else {
    img_cfg.nrows_binned = pimg.nrows_binned > 0 ? pimg.nrows_binned : 3;
    img_cfg.ncols = pimg.ncols > 0 ? pimg.ncols : 50;
    img_cfg.row_binning = pimg.row_binning > 0 ? pimg.row_binning : 100;
    img_cfg.col_binning = pimg.col_binning > 0 ? pimg.col_binning : 1;
    img_cfg.pixel_size_um = pimg.pixel_size_um;
    img_cfg.sigma_readout_e = pimg.sigma_readout_e;
    img_cfg.lambda_dc = pimg.lambda_dc;
    img_cfg.include_dark_current = (pimg.lambda_dc > 0.0);
    img_cfg.rng_seed = pimg.rng_seed;
  }

  PatternClassifierConfig pcc;
  pcc.Qmin_e = pcc_j.Qmin_e;
  pcc.neighbor_Qmax_e = pcc_j.neighbor_Qmax_e;
  pcc.Qmax_e = pcc_j.Qmax_e;
  pcc.enable_MN = pcc_j.enable_MN;
  pcc.enable_MNL = pcc_j.enable_MNL;
  pcc.sigma_res_e = pcc_j.sigma_res_e;
  pcc.max_e_per_pixel = pcc_j.max_e_per_pixel;
  pcc.thr_M = pcc_j.thr_M;
  pcc.thr_MN = pcc_j.thr_MN;
  pcc.thr_MNL = pcc_j.thr_MNL;
  pcc.selection = pcc_j.selection;
  pcc.chain_max_pair = pcc_j.chain_max_pair;
  pcc.chain_allow_22 = pcc_j.chain_allow_22;
  pcc.chain_horizontal_isolation = pcc_j.chain_horizontal_isolation;
  pcc.chain_first_only = pcc_j.chain_first_only;
  if (background_mode) {
    pcc.allow_pattern_zero = true;         // the identification matrix has an empty-pattern row
    pcc.single_pixel_use_round = false;    // the reference uses the plain threshold test
    WarnIfDifferent("pattern_image.sigma_readout_e", pimg.sigma_readout_e, 0.16);
    WarnIfDifferent("pattern_classifier.Qmin_e", pcc.Qmin_e, 0.60);
    WarnIfDifferent("pattern_classifier.neighbor_Qmax_e", pcc.neighbor_Qmax_e, 0.60);
    WarnIfDifferent("pattern_classifier.Qmax_e", pcc.Qmax_e, 5.6);
    WarnIfDifferent("pattern_classifier.thr_M", pcc.thr_M, 3.5);
    WarnIfDifferent("pattern_classifier.thr_MN", pcc.thr_MN, 4.0);
    WarnIfDifferent("pattern_classifier.thr_MNL", pcc.thr_MNL, 5.5);
  } else {
    pcc.allow_pattern_zero = pcc_j.allow_pattern_zero;
    pcc.single_pixel_use_round = pcc_j.single_pixel_use_round;
  }

  PatternTableSetup s;
  s.img_gen = std::make_shared<PatternImageGenerator>(img_cfg, ct);
  s.classifier = std::make_shared<PatternClassifier>(pcc);
  s.selection = pcc.selection;
  s.alpha = emj.alpha;
  s.beta_per_keV = emj.beta_per_keV;
  s.A_um2 = emj.A_um2;
  s.b_umInv = emj.b_umInv;
  s.lambda_dc = pimg.lambda_dc;
  s.sigma_readout_e = pimg.sigma_readout_e;
  return s;
}

// ----------------------------------------------------------------------------
// GenerateSignalTable
//   For n_e = 1..ne_max, generate Nsims images with n_e diffused electrons
//   (Ee = 0, so only alpha matters), scan the middle row and count how many times
//   each pattern is found. Several patterns can be found in one image, so a row of
//   probabilities can add to more than 1, as in the reference table.
// ----------------------------------------------------------------------------
PatternTable GenerateSignalTable(const PatternTableSetup& setup, int Nsims, int ne_max) {
  std::map<int, std::map<int, long>> counts;  // ne -> pattern code -> identifications
  std::set<int> codes;
  for (int ne = 1; ne <= ne_max; ++ne) {
    for (int it = 0; it < Nsims; ++it) {
      const auto image = setup.img_gen->GenerateImage(ne, 0.0);
      if (image.size() < 3) continue;
      const auto prs = setup.classifier->ScanAllImage2DWithIsolation(image[0], image[1], image[2], -1.0);
      for (const auto& pr : prs) {
        const int code = PatternCode(pr.label.q);
        counts[ne][code] += 1;
        codes.insert(code);
      }
    }
  }
  PatternTable table;
  for (int code : codes)
    for (int ne = 1; ne <= ne_max; ++ne) {
      const auto it = counts[ne].find(code);
      const long c = (it != counts[ne].end()) ? it->second : 0;
      table[code][ne] = static_cast<double>(c) / static_cast<double>(Nsims);
    }
  return table;
}

// ----------------------------------------------------------------------------
// GenerateBackgroundMatrix
//   For every ideal pattern I simulate Nsims noisy 3x5 clusters, classify each
//   with the 2D isolation scan and count which patterns are identified; a cluster
//   with no pattern counts as the empty pattern {0}.
// ----------------------------------------------------------------------------
PatternTable GenerateBackgroundMatrix(const PatternTableSetup& setup, int Nsims) {
  const std::vector<std::vector<int>> ideals = EnumerateIdealPatterns(5);
  PatternTable matrix;
  for (const auto& ideal : ideals) {
    double b = 0, c = 0, d = 0;
    if (ideal.size() >= 1) b = static_cast<double>(ideal[0]);
    if (ideal.size() >= 2) c = static_cast<double>(ideal[1]);
    if (ideal.size() >= 3) d = static_cast<double>(ideal[2]);
    std::map<int, int> counts;  // identified code -> number of identifications
    for (int t = 0; t < Nsims; ++t) {
      const auto cl = setup.img_gen->SimulateCluster(b, c, d);
      const auto prs = setup.classifier->ScanAllImage2DWithIsolation(cl[0], cl[1], cl[2], -1.0);
      if (prs.empty()) {
        counts[0]++;
      } else {
        for (const auto& pr : prs) counts[PatternCode(pr.label.q)]++;
      }
    }
    const int inj_code = PatternCode(ideal);
    for (const auto& id_pat : ideals) {
      const int id_code = PatternCode(id_pat);
      matrix[id_code][inj_code] =
          (counts.count(id_code) ? counts[id_code] : 0) / static_cast<double>(Nsims);
    }
  }
  return matrix;
}

// ----------------------------------------------------------------------------
// WriteSignalCsv
//   Header comment with the generation parameters, then pattern,ne,Efficiency rows
//   ordered by n_e and then by ascending pattern code (5 decimals).
// ----------------------------------------------------------------------------
bool WriteSignalCsv(const std::string& path, const PatternTable& table,
                    const PatternTableSetup& setup, int Nsims, int ne_max) {
  std::ofstream of(path);
  if (!of.is_open()) return false;
  of << "# Simulation Parameters: Nsims=" << Nsims << ", alpha=" << setup.alpha
     << ", beta=" << setup.beta_per_keV << ", A=" << setup.A_um2 << ",b=" << setup.b_umInv
     << ",lambda=" << setup.lambda_dc << ", selection=" << setup.selection << "\n";
  of << "pattern,ne,Efficiency\n";
  of << std::fixed << std::setprecision(5);
  for (int ne = 1; ne <= ne_max; ++ne) {
    for (const auto& [code, row] : table) {
      const auto it = row.find(ne);
      of << code << "," << ne << "," << (it != row.end() ? it->second : 0.0) << "\n";
    }
  }
  return true;
}

// ----------------------------------------------------------------------------
// WriteBackgroundCsv
//   iden_pat,eff_<injected>,... with one row per identified pattern in string-sorted
//   order (the order of the reference table) and one column per ideal pattern.
// ----------------------------------------------------------------------------
bool WriteBackgroundCsv(const std::string& path, const PatternTable& matrix) {
  std::ofstream of(path);
  if (!of.is_open()) return false;
  const std::vector<std::vector<int>> ideals = EnumerateIdealPatterns(5);
  of << "iden_pat";
  for (const auto& ideal : ideals) of << ",eff_" << PatternString(ideal);
  of << "\n";
  std::vector<std::vector<int>> rows_sorted = ideals;
  std::sort(rows_sorted.begin(), rows_sorted.end(), [](const std::vector<int>& x, const std::vector<int>& y) {
    return PatternString(x) < PatternString(y);
  });
  for (const auto& id_pat : rows_sorted) {
    of << PatternString(id_pat);
    const auto row_it = matrix.find(PatternCode(id_pat));
    for (const auto& inj : ideals) {
      double v = 0.0;
      if (row_it != matrix.end()) {
        const auto it = row_it->second.find(PatternCode(inj));
        if (it != row_it->second.end()) v = it->second;
      }
      of << "," << std::scientific << v;
    }
    of << "\n";
  }
  return true;
}

// ----------------------------------------------------------------------------
// ReadSignalCsv
//   pattern,ne,Efficiency rows; lines starting with '#' are comments, from which I
//   take the "Nsims=" value when asked; duplicate (pattern, ne) rows are summed.
// ----------------------------------------------------------------------------
PatternTable ReadSignalCsv(const std::string& path, long* nsims_out) {
  PatternTable table;
  if (nsims_out) *nsims_out = 0;
  std::ifstream f(path);
  if (!f) {
    std::cerr << "[tables] Cannot open " << path << "\n";
    return table;
  }
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') {
      const auto pos = line.find("Nsims=");
      if (pos != std::string::npos && nsims_out && *nsims_out == 0) {
        try { *nsims_out = std::stol(line.substr(pos + 6)); } catch (...) {}
      }
      continue;
    }
    if (line.rfind("pattern", 0) == 0) continue;
    std::istringstream ss(line);
    std::string tok;
    int pat = 0, ne = 0;
    double eff = 0.0;
    if (!std::getline(ss, tok, ',')) continue;
    try { pat = std::stoi(tok); } catch (...) { continue; }
    if (!std::getline(ss, tok, ',')) continue;
    try { ne = std::stoi(tok); } catch (...) { continue; }
    if (!std::getline(ss, tok, ',')) continue;
    try { eff = std::stod(tok); } catch (...) { continue; }
    table[pat][ne] += eff;
  }
  return table;
}

// ----------------------------------------------------------------------------
// ReadBackgroundCsv
//   Header iden_pat,eff_<injected>,...; one row per identified pattern.
// ----------------------------------------------------------------------------
PatternTable ReadBackgroundCsv(const std::string& path) {
  PatternTable matrix;
  std::ifstream f(path);
  if (!f) {
    std::cerr << "[tables] Cannot open " << path << "\n";
    return matrix;
  }
  std::string line;
  if (!std::getline(f, line)) return matrix;
  std::vector<int> injected;
  {
    std::istringstream ss(line);
    std::string tok;
    std::getline(ss, tok, ',');
    while (std::getline(ss, tok, ',')) {
      const auto pos = tok.rfind('_');
      try { injected.push_back(std::stoi(tok.substr(pos + 1))); } catch (...) { injected.push_back(-1); }
    }
  }
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    std::istringstream ss(line);
    std::string tok;
    if (!std::getline(ss, tok, ',')) continue;
    int iden = 0;
    try { iden = std::stoi(tok); } catch (...) { continue; }
    for (int inj : injected) {
      if (!std::getline(ss, tok, ',')) break;
      double v = 0.0;
      try { v = std::stod(tok); } catch (...) {}
      matrix[iden][inj] = v;
    }
  }
  return matrix;
}

// ----------------------------------------------------------------------------
// PooledPull
//   (p_new - p_ref) / sigma, with sigma from the pooled estimate of the two
//   binomial proportions; zero when the pooled variance vanishes.
// ----------------------------------------------------------------------------
double PooledPull(double p_new, double n_new, double p_ref, double n_ref) {
  const double k = p_new * n_new + p_ref * n_ref;
  const double ph = k / (n_new + n_ref);
  const double var = ph * (1.0 - ph) * (1.0 / n_new + 1.0 / n_ref);
  return var <= 0.0 ? 0.0 : (p_new - p_ref) / std::sqrt(var);
}

// ----------------------------------------------------------------------------
// AverageTables
//   Mean over independent runs of the same table; an entry missing from a run is 0.
// ----------------------------------------------------------------------------
PatternTable AverageTables(const std::vector<PatternTable>& tables) {
  PatternTable out;
  if (tables.empty()) return out;
  for (const auto& t : tables)
    for (const auto& [a, row] : t)
      for (const auto& [b, v] : row) out[a][b] += v;
  const double n = static_cast<double>(tables.size());
  for (auto& [a, row] : out)
    for (auto& [b, v] : row) v /= n;
  return out;
}

// ----------------------------------------------------------------------------
// CompareTables
//   Union of the keys of both tables; pull of each entry; summary over the entries
//   that are non-zero in either table.
// ----------------------------------------------------------------------------
TableComparison CompareTables(const PatternTable& ref, const PatternTable& neu,
                              double n_ref, double n_new, const std::set<int>& patterns) {
  std::set<std::pair<int, int>> keys;
  for (const auto& [a, row] : ref) for (const auto& kv : row) keys.insert({a, kv.first});
  for (const auto& [a, row] : neu) for (const auto& kv : row) keys.insert({a, kv.first});
  auto value = [](const PatternTable& t, int a, int b) {
    const auto ia = t.find(a);
    if (ia == t.end()) return 0.0;
    const auto ib = ia->second.find(b);
    return ib == ia->second.end() ? 0.0 : ib->second;
  };
  TableComparison cmp;
  double sum_sq = 0.0;
  for (const auto& [a, b] : keys) {
    if (!patterns.empty() && !patterns.count(a)) continue;
    TableComparisonEntry e;
    e.key_a = a;
    e.key_b = b;
    e.ref = value(ref, a, b);
    e.neu = value(neu, a, b);
    e.pull = PooledPull(e.neu, n_new, e.ref, n_ref);
    cmp.entries.push_back(e);
    if (e.ref > 0.0 || e.neu > 0.0) {
      ++cmp.n_nonzero;
      sum_sq += e.pull * e.pull;
      if (std::abs(e.pull) > 3.0) ++cmp.n_beyond_3;
      if (std::abs(e.pull) > 2.0) ++cmp.n_beyond_2;
      cmp.max_abs_pull = std::max(cmp.max_abs_pull, std::abs(e.pull));
    }
  }
  if (cmp.n_nonzero > 0) cmp.rms_pull = std::sqrt(sum_sq / cmp.n_nonzero);
  return cmp;
}

// ----------------------------------------------------------------------------
// PrintComparison
//   Summary of a comparison: entry counts, RMS pull, entries beyond 2 and 3 sigma
//   against the number expected from a faithful reproduction, the per-n_e sums for
//   a signal table, and the eight largest pulls.
// ----------------------------------------------------------------------------
void PrintComparison(std::ostream& os, const TableComparison& cmp, bool signal_layout) {
  os << (signal_layout ? "signal" : "background") << ": " << cmp.entries.size() << " entries, "
     << cmp.n_nonzero << " non-zero in either table\n";
  os << std::fixed << std::setprecision(2);
  os << "  RMS pull (non-zero entries): " << cmp.rms_pull << "   chi2/n: " << cmp.rms_pull * cmp.rms_pull << "\n";
  os << std::setprecision(1);
  os << "  entries beyond 3 sigma: " << cmp.n_beyond_3 << " (expected " << 0.0027 * cmp.n_nonzero
     << ")   beyond 2 sigma: " << cmp.n_beyond_2 << " (expected " << 0.0455 * cmp.n_nonzero
     << ")   max |pull|: " << std::setprecision(2) << cmp.max_abs_pull << "\n";
  if (signal_layout) {
    struct Acc { double ref = 0.0, neu = 0.0, sq = 0.0; int n = 0; };
    std::map<int, Acc> by_ne;
    for (const auto& e : cmp.entries) {
      if (!(e.ref > 0.0 || e.neu > 0.0)) continue;
      auto& a = by_ne[e.key_b];
      a.ref += e.ref; a.neu += e.neu; a.sq += e.pull * e.pull; ++a.n;
    }
    os << "  per n_e:  sum(ref)  sum(new)  RMS pull\n";
    for (const auto& [ne, a] : by_ne)
      os << "    " << std::setw(2) << ne << "   " << std::setw(8) << std::setprecision(4) << a.ref << " "
         << std::setw(8) << a.neu << "  " << std::setw(8) << std::setprecision(2)
         << std::sqrt(a.sq / a.n) << "\n";
  }
  std::vector<TableComparisonEntry> worst;
  for (const auto& e : cmp.entries)
    if (e.ref > 0.0 || e.neu > 0.0) worst.push_back(e);
  std::sort(worst.begin(), worst.end(), [](const TableComparisonEntry& x, const TableComparisonEntry& y) {
    return std::abs(x.pull) > std::abs(y.pull);
  });
  os << "  largest pulls:\n";
  os << std::setprecision(5) << std::defaultfloat;
  for (std::size_t i = 0; i < worst.size() && i < 8; ++i) {
    const auto& e = worst[i];
    os << "    (" << e.key_a << ", " << e.key_b << "): ref " << e.ref << "  new " << e.neu
       << "  pull " << std::showpos << std::fixed << std::setprecision(2) << e.pull
       << std::noshowpos << std::defaultfloat << std::setprecision(5) << "\n";
  }
}

// ----------------------------------------------------------------------------
// WriteComparisonCsv
//   One row per compared entry with the reference value, the new value, the
//   difference and the pull.
// ----------------------------------------------------------------------------
bool WriteComparisonCsv(const std::string& path, const TableComparison& cmp, bool signal_layout) {
  std::ofstream of(path);
  if (!of.is_open()) return false;
  of << (signal_layout ? "pattern,ne" : "iden_pat,inj_pat") << ",ref,new,diff,pull\n";
  of << std::setprecision(6);
  for (const auto& e : cmp.entries)
    of << e.key_a << "," << e.key_b << "," << e.ref << "," << e.neu << ","
       << (e.neu - e.ref) << "," << e.pull << "\n";
  return true;
}

// ----------------------------------------------------------------------------
// ParseTableCheckArgs
//   Positional values first (n_trials, reference csv), then the --flags described
//   in PatternTables.hh. A flag with a list (--new, --patterns) takes every
//   following token up to the next flag.
// ----------------------------------------------------------------------------
TableCheckArgs ParseTableCheckArgs(int argc, char** argv) {
  TableCheckArgs a;
  auto is_flag = [](const std::string& t) { return t.size() > 2 && t[0] == '-' && t[1] == '-'; };
  for (int i = 1; i < argc; ++i) {
    const std::string t = argv[i];
    if (!is_flag(t)) { a.positional.push_back(t); continue; }
    if (t == "--config" && i + 1 < argc) { a.config = argv[++i]; }
    else if (t == "--out" && i + 1 < argc) { a.out_csv = argv[++i]; }
    else if (t == "--n-ref" && i + 1 < argc) { a.n_ref = std::stod(argv[++i]); }
    else if (t == "--n-new" && i + 1 < argc) { a.n_new = std::stod(argv[++i]); }
    else if (t == "--new") {
      while (i + 1 < argc && !is_flag(argv[i + 1])) a.new_files.push_back(argv[++i]);
    } else if (t == "--patterns") {
      while (i + 1 < argc && !is_flag(argv[i + 1])) a.patterns.insert(std::stoi(argv[++i]));
    } else {
      std::cerr << "[tables] unknown or incomplete option " << t << "\n";
    }
  }
  return a;
}

namespace {

// ----------------------------------------------------------------------------
// FinishComparison
//   Print the comparison, write the optional CSV and decide the exit code.
// ----------------------------------------------------------------------------
int FinishComparison(const TableComparison& cmp, bool signal_layout, const TableCheckArgs& args) {
  PrintComparison(std::cout, cmp, signal_layout);
  if (!args.out_csv.empty()) {
    if (WriteComparisonCsv(args.out_csv, cmp, signal_layout)) std::cout << "  wrote " << args.out_csv << "\n";
    else std::cerr << "  cannot write " << args.out_csv << "\n";
  }
  const bool pass = (cmp.rms_pull <= 1.5 && cmp.max_abs_pull <= 5.0);
  std::cout << "  RESULT: " << (pass ? "PASS" : "FAIL")
            << " (RMS pull <= 1.5 and no entry beyond 5 sigma)\n";
  return pass ? 0 : 1;
}

} // namespace

// ----------------------------------------------------------------------------
// RunCsvTableComparison
//   Read the reference and the new table(s), average the new ones, take the trial
//   counts from the flags or the CSV headers, and compare.
// ----------------------------------------------------------------------------
int RunCsvTableComparison(bool signal_layout, const std::string& ref_csv, const TableCheckArgs& args) {
  long nsims_ref = 0;
  PatternTable ref = signal_layout ? ReadSignalCsv(ref_csv, &nsims_ref) : ReadBackgroundCsv(ref_csv);
  if (ref.empty()) { std::cerr << "[tables] reference table is empty: " << ref_csv << "\n"; return 2; }
  std::vector<PatternTable> runs;
  long nsims_new = 0;
  for (const auto& f : args.new_files) {
    long n = 0;
    PatternTable t = signal_layout ? ReadSignalCsv(f, &n) : ReadBackgroundCsv(f);
    if (t.empty()) { std::cerr << "[tables] new table is empty: " << f << "\n"; return 2; }
    if (nsims_new == 0) nsims_new = n;
    runs.push_back(std::move(t));
  }
  const double n_ref = args.n_ref > 0 ? args.n_ref : (nsims_ref > 0 ? nsims_ref : (signal_layout ? 0.0 : 1e5));
  const double n_one = args.n_new > 0 ? args.n_new : static_cast<double>(nsims_new);
  if (n_ref <= 0 || n_one <= 0) {
    std::cerr << "[tables] give --n-ref and --n-new (the CSV header has no Nsims= value)\n";
    return 2;
  }
  const PatternTable neu = AverageTables(runs);
  std::cout << "reference " << ref_csv << "  (" << n_ref << " trials/entry)\n";
  std::cout << "new       " << args.new_files.size() << " file(s), " << n_one << " trials each\n";
  const TableComparison cmp =
      CompareTables(ref, neu, n_ref, n_one * static_cast<double>(runs.size()), args.patterns);
  return FinishComparison(cmp, signal_layout, args);
}

// ----------------------------------------------------------------------------
// RunGeneratedTableComparison
//   Compare a table generated in this process with the reference CSV.
// ----------------------------------------------------------------------------
int RunGeneratedTableComparison(bool signal_layout, const std::string& ref_csv,
                                const PatternTable& generated, double n_new,
                                const TableCheckArgs& args) {
  long nsims_ref = 0;
  PatternTable ref = signal_layout ? ReadSignalCsv(ref_csv, &nsims_ref) : ReadBackgroundCsv(ref_csv);
  if (ref.empty()) { std::cerr << "[tables] reference table is empty: " << ref_csv << "\n"; return 2; }
  const double n_ref = args.n_ref > 0 ? args.n_ref : (nsims_ref > 0 ? nsims_ref : (signal_layout ? 0.0 : 1e5));
  if (n_ref <= 0) {
    std::cerr << "[tables] give --n-ref (the reference header has no Nsims= value)\n";
    return 2;
  }
  std::cout << "reference " << ref_csv << "  (" << n_ref << " trials/entry)\n";
  std::cout << "generated " << n_new << " trials/entry\n";
  const TableComparison cmp = CompareTables(ref, generated, n_ref, n_new, args.patterns);
  return FinishComparison(cmp, signal_layout, args);
}

} // namespace ccdarksens
