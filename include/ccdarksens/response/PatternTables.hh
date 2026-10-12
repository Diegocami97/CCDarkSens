// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  PatternTables.hh -- Generation, reading, writing and comparison of the two
//  pattern efficiency tables: the signal table P(pattern | n_e) and the
//  background identification matrix B[identified | injected]. The generator
//  apps and the validators share this code so that one config produces and
//  checks the same table.
// ===========================================================================

#pragma once

#include <iosfwd>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/response/PatternClassifier.hh"
#include "ccdarksens/response/PatternImageGenerator.hh"

namespace ccdarksens {

/// Table keyed by two integers: [pattern code][n_e] for the signal table and
/// [identified code][injected code] for the background matrix.
using PatternTable = std::map<int, std::map<int, double>>;

// ----------------------------------------------------------------------------
// PatternTableSetup
//   Everything the two generators need, built once from a config: the image
//   generator, the classifier and the parameters I print in the CSV header.
// ----------------------------------------------------------------------------
struct PatternTableSetup {
  std::shared_ptr<PatternImageGenerator> img_gen;  ///< 2D image / cluster generator
  std::shared_ptr<PatternClassifier>     classifier; ///< pattern selection (tree or chain)
  std::string selection;      ///< "tree" or "chain", as set in the config
  double alpha = 0.0;         ///< diffusion alpha of the config
  double beta_per_keV = 0.0;  ///< diffusion beta of the config
  double A_um2 = 0.0;         ///< diffusion A of the config
  double b_umInv = 0.0;       ///< diffusion b of the config
  double lambda_dc = 0.0;     ///< dark-current pile-up per binned pixel
  double sigma_readout_e = 0.0; ///< readout noise of the generated images
};

/// Build the setup from a parsed config. In background mode I use the cluster
/// geometry and classifier flags of the reference background table (empty pattern
/// allowed, no single-pixel rounding) and warn when the noise or thresholds differ
/// from the values that table was validated with.
PatternTableSetup MakePatternTableSetup(const ConfigManager& cfg, bool background_mode);

/// Simulate Nsims images with n_e = 1..ne_max electrons and return P(pattern | n_e)
/// as table[pattern][n_e]; every pattern seen at least once has an entry for every n_e.
PatternTable GenerateSignalTable(const PatternTableSetup& setup, int Nsims, int ne_max);

/// Simulate Nsims noisy 3x5 clusters for each ideal pattern and return
/// B[identified | injected] as table[identified][injected]. A cluster with no
/// pattern counts as the empty pattern {0}.
PatternTable GenerateBackgroundMatrix(const PatternTableSetup& setup, int Nsims);

/// All ideal patterns with total charge <= max_sum: the empty pattern first, then
/// the one-, two- and three-pixel patterns with 1..5 e- per pixel.
std::vector<std::vector<int>> EnumerateIdealPatterns(int max_sum);

/// Pattern digits to integer code, e.g. {2,1,1} -> 211.
int PatternCode(const std::vector<int>& p);

/// Write a signal table in the layout of data/Efficiencies_patterns_*.csv.
bool WriteSignalCsv(const std::string& path, const PatternTable& table,
                    const PatternTableSetup& setup, int Nsims, int ne_max);

/// Write a background matrix in the layout of data/Background_efficiencies.csv
/// (rows = identified pattern, string-sorted; columns = eff_<injected>).
bool WriteBackgroundCsv(const std::string& path, const PatternTable& matrix);

/// Read a signal table (pattern,ne,Efficiency). Duplicate rows are summed. If
/// nsims_out is given it receives the "Nsims=" value of the header comment, or 0.
PatternTable ReadSignalCsv(const std::string& path, long* nsims_out = nullptr);

/// Read a background matrix (iden_pat,eff_<injected>,...).
PatternTable ReadBackgroundCsv(const std::string& path);

// ----------------------------------------------------------------------------
// TableComparison
//   Entry-by-entry comparison with a pooled two-sample binomial test.
// ----------------------------------------------------------------------------
struct TableComparisonEntry {
  int key_a = 0;      ///< pattern code (signal) or identified code (background)
  int key_b = 0;      ///< n_e (signal) or injected code (background)
  double ref = 0.0;   ///< reference value
  double neu = 0.0;   ///< new value
  double pull = 0.0;  ///< (new - ref) / pooled sigma
};

struct TableComparison {
  std::vector<TableComparisonEntry> entries; ///< union of both tables
  int n_nonzero = 0;       ///< entries non-zero in either table
  double rms_pull = 0.0;   ///< RMS pull over the non-zero entries
  int n_beyond_3 = 0;      ///< non-zero entries with |pull| > 3
  int n_beyond_2 = 0;      ///< non-zero entries with |pull| > 2
  double max_abs_pull = 0.0;
};

/// Pooled two-sample binomial pull of two probabilities from n_new and n_ref trials.
double PooledPull(double p_new, double n_new, double p_ref, double n_ref);

/// Compare two tables. If patterns is not empty only entries whose first key is in
/// it are used (signal tables).
TableComparison CompareTables(const PatternTable& ref, const PatternTable& neu,
                              double n_ref, double n_new,
                              const std::set<int>& patterns = {});

/// Average several tables of equal meaning (independent runs); missing entries are 0.
PatternTable AverageTables(const std::vector<PatternTable>& tables);

/// Print the summary of a comparison; signal_layout adds the per-n_e sums.
void PrintComparison(std::ostream& os, const TableComparison& cmp, bool signal_layout);

/// Write the per-entry comparison as CSV.
bool WriteComparisonCsv(const std::string& path, const TableComparison& cmp, bool signal_layout);

// ----------------------------------------------------------------------------
// TableCheckArgs
//   Command line shared by the two validators: up to two positional values
//   (n_trials, reference csv) kept from the original usage, plus
//     --config cfg.json        generate the table from this config (the same one the
//                              generator app uses) instead of the built-in constants
//     --new f1.csv [f2.csv..]  compare existing tables with the reference, no simulation
//     --n-ref N --n-new N      trials behind each reference / new entry (default: the
//                              Nsims= value in the CSV header)
//     --patterns c1 c2 ..      signal only: restrict to these pattern codes
//     --out file.csv           write the per-entry comparison
// ----------------------------------------------------------------------------
struct TableCheckArgs {
  std::vector<std::string> positional;  ///< n_trials, reference csv path
  std::string config;                   ///< --config
  std::vector<std::string> new_files;   ///< --new
  double n_ref = 0.0;                   ///< --n-ref (0 = from the header or the default)
  double n_new = 0.0;                   ///< --n-new (0 = from the header)
  std::set<int> patterns;               ///< --patterns
  std::string out_csv;                  ///< --out
};

/// Parse the validator command line described above.
TableCheckArgs ParseTableCheckArgs(int argc, char** argv);

/// Compare the CSV tables of args.new_files (averaged) with ref_csv and print the
/// summary. Returns 0 if the RMS pull is <= 1.5 and no entry is beyond 5 sigma,
/// 1 otherwise, 2 on a read error.
int RunCsvTableComparison(bool signal_layout, const std::string& ref_csv, const TableCheckArgs& args);

/// Same comparison for a table generated in memory with n_new trials per entry.
int RunGeneratedTableComparison(bool signal_layout, const std::string& ref_csv,
                                const PatternTable& generated, double n_new,
                                const TableCheckArgs& args);

} // namespace ccdarksens
