// ============================================================================
//  CCDarkSens — ccdarksens_pattern_background_eff
//  Monte Carlo app that simulates ideal patterns through ChargeTransport and PatternClassifier and writes background identification-efficiency CSVs.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/PatternClassifier.hh"
#include "ccdarksens/response/PatternImageGenerator.hh"

using namespace ccdarksens;

// Enumerate ideal patterns: (0) first (background), then (1)..(5), (1,1),..., (1,1,1),... with sum <= max_sum (notebook pattern_simulation_background).
static std::vector<std::vector<int>> enumerate_ideal_patterns(int max_sum) {
  std::vector<std::vector<int>> out;
  out.push_back({0});  // background / empty
  std::set<std::vector<int>> uniq;
  for (int a = 1; a <= 5 && a <= max_sum; ++a)
    uniq.insert({a});
  for (int a = 1; a <= 5; ++a) {
    for (int b = 1; b <= 5; ++b) {
      if (a + b > max_sum) continue;
      uniq.insert({a, b});
    }
  }
  for (int a = 1; a <= 5; ++a) {
    for (int b = 1; b <= 5; ++b) {
      for (int c = 1; c <= 5; ++c) {
        if (a + b + c > max_sum) continue;
        uniq.insert({a, b, c});
      }
    }
  }
  std::vector<std::vector<int>> rest(uniq.begin(), uniq.end());
  std::sort(rest.begin(), rest.end(), [](const std::vector<int>& x, const std::vector<int>& y) {
    if (x.size() != y.size()) return x.size() < y.size();
    return x < y;
  });
  for (const auto& p : rest) out.push_back(p);
  return out;
}

static int pattern_to_code(const std::vector<int>& p) {
  int code = 0;
  for (int d : p) code = code * 10 + d;
  return code;
}

static std::string pattern_to_str(const std::vector<int>& p) {
  std::string s;
  for (int d : p) s += std::to_string(d);
  return s;
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json [Nsims] [outdir]\n";
    return 1;
  }
  const std::string config_path = argv[1];
  int Nsims = 10000;
  if (argc >= 3) Nsims = std::max(1, std::stoi(argv[2]));
  std::string outdir = ".";
  if (argc >= 4) outdir = argv[3];

  try {
    ConfigManager cfg(config_path);
    cfg.parse();
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
    const bool from_json = (pimg.nrows_binned > 0 && pimg.ncols > 0 && pimg.row_binning > 0);
    const bool from_detector = !from_json && (det.geometry().rows > 0 && det.geometry().cols > 0 && pimg.row_binning > 0 && pimg.col_binning > 0);
    if (from_json) {
      img_cfg.nrows_binned = pimg.nrows_binned;
      img_cfg.ncols        = pimg.ncols;
    } else if (from_detector) {
      img_cfg.raw_rows = det.geometry().rows;
      img_cfg.raw_cols = det.geometry().cols;
    } else {
      img_cfg.nrows_binned = pimg.nrows_binned > 0 ? pimg.nrows_binned : 3;
      img_cfg.ncols        = pimg.ncols > 0 ? pimg.ncols : 50;
    }
    img_cfg.row_binning    = pimg.row_binning > 0 ? pimg.row_binning : 100;
    img_cfg.col_binning    = pimg.col_binning > 0 ? pimg.col_binning : 1;
    img_cfg.pixel_size_um  = pimg.pixel_size_um;
    img_cfg.sigma_readout_e= pimg.sigma_readout_e;
    img_cfg.lambda_dc      = pimg.lambda_dc;
    img_cfg.rng_seed       = pimg.rng_seed;
    auto img_gen = std::make_shared<PatternImageGenerator>(img_cfg, ct);

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
    pcc.allow_pattern_zero = true;  // notebook identify_pattern_background: allow pattern (0)
    pcc.single_pixel_use_round = true;  // round(q) for single-pixel to match reference diagonal (2)->(2), etc.
    auto classifier = std::make_shared<PatternClassifier>(pcc);

    std::vector<std::vector<int>> ideals = enumerate_ideal_patterns(5);
    std::map<int, std::map<int, double>> eff;  // ideal_code -> (identified_code -> efficiency)

    for (const auto& ideal : ideals) {
      double b = 0, c = 0, d = 0;
      if (ideal.size() >= 1) b = static_cast<double>(ideal[0]);
      if (ideal.size() >= 2) c = static_cast<double>(ideal[1]);
      if (ideal.size() >= 3) d = static_cast<double>(ideal[2]);
      std::map<int, int> counts;
      for (int t = 0; t < Nsims; ++t) {
        auto cl = img_gen->SimulateCluster(b, c, d);
        // 2D scan with isolation (notebook scan_image_background uses 3-row image)
        auto res = classifier->ScanImage2DWithIsolation(cl[0], cl[1], cl[2], -1.0);
        int code = 0;
        if (res.valid)
          code = pattern_to_code(res.label.q);
        counts[code]++;
      }
      int ideal_code = pattern_to_code(ideal);
      for (const auto& id_ideal : ideals) {
        int id_code = pattern_to_code(id_ideal);
        double e = (counts.count(id_code) ? counts[id_code] : 0) / static_cast<double>(Nsims);
        eff[ideal_code][id_code] = e;
      }
    }

    std::string csv_path = outdir + "/Background_efficiencies.csv";
    std::ofstream of(csv_path);
    if (!of.is_open()) {
      std::cerr << "Cannot write " << csv_path << "\n";
      return 1;
    }
    of << "iden_pat";
    for (const auto& ideal : ideals)
      of << ",eff_" << pattern_to_str(ideal);
    of << "\n";
    for (const auto& ideal : ideals) {
      of << pattern_to_str(ideal);
      int ideal_code = pattern_to_code(ideal);
      const auto& row = eff[ideal_code];
      for (const auto& id_ideal : ideals) {
        int id_code = pattern_to_code(id_ideal);
        auto it = row.find(id_code);
        of << "," << std::scientific << (it != row.end() ? it->second : 0.0);
      }
      of << "\n";
    }
    of.close();
    std::cout << "Wrote " << csv_path << " (Nsims=" << Nsims << ", " << ideals.size() << " ideal patterns).\n";
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
