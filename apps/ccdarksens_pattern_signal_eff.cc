// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ccdarksens_pattern_signal_eff.cc -- Monte Carlo app that deposits n_e
//  electrons in a simulated 3x50 binned image, runs the configured pattern
//  selection and writes P(pattern | n_e) in the layout of
//  data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv.
// ===========================================================================

#include <algorithm>
#include <iostream>
#include <string>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/response/PatternTables.hh"

using namespace ccdarksens;

// ----------------------------------------------------------------------------
// main
//   For each n_e = 1..ne_max I generate Nsims images with n_e diffused electrons
//   (Ee = 0, so only alpha matters), scan the middle row with the selection set in
//   response.pattern_classifier.selection ("tree" or "chain"), and count how many
//   times each pattern is found (the generation itself is GenerateSignalTable in
//   PatternTables.cc, shared with ccdarksens_validate_pattern_efficiency).
//   Writes <outdir>/Efficiencies_patterns.csv.
//   Usage: <program> config.json [Nsims=100000] [outdir=.] [ne_max=5]. Returns 0 on success, 1 on error.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json [Nsims] [outdir] [ne_max]\n";
    return 1;
  }
  const std::string config_path = argv[1];
  int Nsims = 100000;
  if (argc >= 3) Nsims = std::max(1, std::stoi(argv[2]));
  std::string outdir = ".";
  if (argc >= 4) outdir = argv[3];
  int ne_max = 5;
  if (argc >= 5) ne_max = std::max(1, std::stoi(argv[4]));

  try {
    ConfigManager cfg(config_path);
    cfg.parse();
    const PatternTableSetup setup = MakePatternTableSetup(cfg, /*background_mode=*/false);

    std::cout << "[signal-eff] selection=" << setup.selection << "  alpha=" << setup.alpha
              << "  lambda_dc=" << setup.lambda_dc << "  sigma_ro=" << setup.sigma_readout_e
              << "  Nsims=" << Nsims << "\n";

    const PatternTable table = GenerateSignalTable(setup, Nsims, ne_max);

    const std::string csv_path = outdir + "/Efficiencies_patterns.csv";
    if (!WriteSignalCsv(csv_path, table, setup, Nsims, ne_max)) {
      std::cerr << "Cannot write " << csv_path << "\n";
      return 1;
    }
    std::cout << "Wrote " << csv_path << " (" << table.size() << " patterns, n_e = 1.." << ne_max << ").\n";
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << "\n";
    return 1;
  }
  return 0;
}
