// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ccdarksens_pattern_background_eff.cc -- Monte Carlo app that simulates
//  ideal patterns through ChargeTransport and PatternClassifier and writes
//  background identification-efficiency CSVs.
// ===========================================================================

#include <algorithm>
#include <iostream>
#include <string>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/response/PatternTables.hh"

using namespace ccdarksens;

// ----------------------------------------------------------------------------
// main
//   Monte Carlo of the background pattern-identification matrix. For every ideal
//   pattern I simulate Nsims noisy 3x5 clusters, classify each with the 2D isolation
//   scan and count which patterns are identified (a cluster with none counts as
//   the empty pattern {0}), then write B[identified | injected] to
//   <outdir>/Background_efficiencies.csv in the layout of data/Background_efficiencies.csv:
//   one row per identified pattern (string-sorted), one column eff_<injected> per ideal pattern.
//   I warn when the config's noise or thresholds differ from the values the shipped table was
//   validated with (sigma_ro 0.16, Qmin 0.6, Qmax 5.6, thresholds 3.5/4.0/5.5). The generation
//   itself is GenerateBackgroundMatrix in PatternTables.cc, shared with
//   ccdarksens_validate_background_efficiency.
//   Usage: <program> config.json [Nsims=100000] [outdir=.]. Returns 0 on success, 1 on error.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json [Nsims] [outdir]\n";
    return 1;
  }
  const std::string config_path = argv[1];
  int Nsims = 100000;  // simulated clusters per ideal pattern
  if (argc >= 3) Nsims = std::max(1, std::stoi(argv[2]));
  std::string outdir = ".";
  if (argc >= 4) outdir = argv[3];

  try {
    ConfigManager cfg(config_path);
    cfg.parse();
    const PatternTableSetup setup = MakePatternTableSetup(cfg, /*background_mode=*/true);
    const PatternTable matrix = GenerateBackgroundMatrix(setup, Nsims);

    const std::string csv_path = outdir + "/Background_efficiencies.csv";
    if (!WriteBackgroundCsv(csv_path, matrix)) {
      std::cerr << "Cannot write " << csv_path << "\n";
      return 1;
    }
    std::cout << "Wrote " << csv_path << " (Nsims=" << Nsims << ", "
              << EnumerateIdealPatterns(5).size() << " ideal patterns).\n";
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
