// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_validate_response_factory.cc
//  Slice 4 validation: diff ResponseFactory + BackgroundFactory's B_pat
//  against the reference app's (ccdarksens_scan_dmelectron_pattern.cc) own
//  printed B_pat values, for both background_source modes
//  ("dc_flat_migration"-style default, and "bp_br_template").
//
//  Usage: ccdarksens_validate_response_factory <config.json> <pattern_id>:<B_pat_ref> [...]
//             [--dump-sb=<path.csv> --dump-sb-mchi=<mass_val> --dump-sb-sigma=<coupling_str>]
//
//  --dump-sb (with --dump-sb-mchi/--dump-sb-sigma) folds the config's
//  model.* signal spectrum at one representative grid point through the
//  same ResponseFold + BackgroundFactory this app already builds, and
//  writes S(bin)/B(bin) per response bin -- for plotting the folded
//  spectra stage of the pipeline. Only applies to cluster_energy /
//  pattern / n_e analysis spaces that expose per-bin S and B (the same
//  ones the rest of this app already handles).
// ===========================================================================

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TH1D.h"

#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/model/ModelFactory.hh"
#include "ccdarksens/response/BackgroundFactory.hh"
#include "ccdarksens/response/ResponseFactory.hh"

using namespace ccdarksens;

namespace {
// ----------------------------------------------------------------------------
// FlagValue
//   Value of the first command-line argument that starts with prefix (e.g. "--dump-sb="), or an empty string.
// ----------------------------------------------------------------------------
std::string FlagValue(int argc, char** argv, const std::string& prefix) {
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg.rfind(prefix, 0) == 0) return arg.substr(prefix.size());
  }
  return "";
}
}  // namespace

// ----------------------------------------------------------------------------
// main
//   Check ResponseFactory/BackgroundFactory against reference numbers: build the fold
//   and background from the config, then compare the computed B_pat of each ROI bin with the
//   <pattern_id>:<B_pat_ref> values given on the command line (relative difference below
//   1e-6 counts as a match). With --dump-sb I also write the folded S and B per bin for one
//   signal point. Returns 0 if every given reference matches.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "Usage: %s <config.json> <pattern_id>:<B_pat_ref> [...]\n", argv[0]);
    return 1;
  }
  const std::string config_path = argv[1];
  const std::string dump_sb_path = FlagValue(argc, argv, "--dump-sb=");
  const std::string dump_sb_mchi = FlagValue(argc, argv, "--dump-sb-mchi=");
  const std::string dump_sb_sigma = FlagValue(argc, argv, "--dump-sb-sigma=");

  std::map<int, double> ref_by_pattern;  // reference B_pat per pattern (or n_e) bin from the command line
  for (int i = 2; i < argc; ++i) {
    std::string tok = argv[i];
    if (tok.rfind("--", 0) == 0) continue;
    auto pos = tok.find(':');
    if (pos == std::string::npos) continue;
    int pid = std::stoi(tok.substr(0, pos));
    double val = std::stod(tok.substr(pos + 1));
    ref_by_pattern[pid] = val;
  }

  ConfigManager cfg(config_path);
  cfg.parse();
  ExperimentSetup setup(cfg.experiment_cfg(), cfg.detector().mass_kg(), cfg.run().rng_seed);
  auto summary = setup.prepare_summary();

  auto rf = MakeResponseFold(cfg, summary, config_path);
  auto bg = MakeBackground(cfg, summary, *rf.fold, rf);

  std::printf("[validate-response-factory] background_source=%s background_model=%s\n",
              cfg.run().background_source.c_str(), cfg.run().background_model.c_str());

  if (!dump_sb_path.empty()) {
    if (dump_sb_mchi.empty() || dump_sb_sigma.empty()) {
      std::fprintf(stderr, "--dump-sb requires --dump-sb-mchi=<mass_val> and --dump-sb-sigma=<coupling_str>\n");
      return 2;
    }
    const double mchi = std::stod(dump_sb_mchi);
    bool sig_ok = false;
    auto dRdE_sig = MakeSignalSpectrumE(cfg.model(), mchi, dump_sb_sigma, &sig_ok);
    std::vector<double> S = rf.fold->Fold(*dRdE_sig, summary.exposure_kg_year);
    const auto& B = bg.B_pat;

    std::ofstream out(dump_sb_path);
    out << "# folded S(bin)/B(bin), model.type=" << cfg.model().type
        << " mchi=" << mchi << " coupling=" << dump_sb_sigma
        << " exposure_kg_year=" << summary.exposure_kg_year
        << " signal_csv_found=" << (sig_ok ? "1" : "0") << "\n";
    auto labels = rf.fold->BinLabels();
    out << "bin_index,bin_label,S,B\n";
    for (std::size_t i = 0; i < S.size(); ++i) {
      const double b = (i < B.size()) ? B[i] : 0.0;
      const std::string lbl = (i < labels.size()) ? labels[i] : "";
      out << i << "," << lbl << "," << S[i] << "," << b << "\n";
    }
    std::printf("[validate-response-factory] wrote %s (%zu bins, signal_csv_found=%s)\n",
                dump_sb_path.c_str(), S.size(), sig_ok ? "yes" : "no");
  }

  bool all_ok = true;  // stays true while every reference matches
  const auto& roi = summary.pattern_roi.empty() ? summary.roi_bins : summary.pattern_roi;
  for (std::size_t i = 0; i < roi.size() && i < bg.B_pat.size(); ++i) {
    const int pid = roi[i];
    const double computed = bg.B_pat[i];
    auto it = ref_by_pattern.find(pid);
    if (it == ref_by_pattern.end()) {
      std::printf("  bin %d: computed B_pat = %.8e (no reference given)\n", pid, computed);
      continue;
    }
    const double ref = it->second;
    const double rel = (ref != 0.0) ? std::fabs(computed - ref) / std::fabs(ref) : std::fabs(computed);
    const bool ok = rel < 1e-6;
    all_ok = all_ok && ok;
    std::printf("  bin %d: computed = %.8e  ref = %.8e  rel_diff = %.3e  %s\n", pid, computed, ref, rel,
                ok ? "OK" : "MISMATCH");
  }

  std::printf(all_ok ? "[validate-response-factory] PASS\n" : "[validate-response-factory] FAIL\n");
  return all_ok ? 0 : 1;
}
