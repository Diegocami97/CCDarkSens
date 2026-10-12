// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ccdarksens_classify_pair.cc -- Hand check of the pattern identification:
//  takes the charges of consecutive pixels of one row, prints the p_mn
//  statistics of the data-selection note and the label the framework's own
//  PatternClassifier returns for the "chain" and "tree" selections.
// ===========================================================================

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ccdarksens/response/PatternClassifier.hh"

using namespace ccdarksens;

// ----------------------------------------------------------------------------
// PrintUsage
//   One-screen description of the command line.
// ----------------------------------------------------------------------------
static void PrintUsage(const char* prog) {
  std::cerr
      << "Usage: " << prog << " q1 q2 [q3 ...] (--sigma S | --bound B) [options]\n"
      << "  q1 q2 ...        calibrated charges (e-) of consecutive pixels of one row\n"
      << "  --sigma S        readout noise of the CCD (e-); default 0.16\n"
      << "  --bound B        1 e- bound of the CCD (= 3.75 sigma), alternative to --sigma\n"
      << "  --qmax Q         maximum total charge of a pattern (default 5.6 e-)\n"
      << "  --above a1 a2..  charges of the row above the pixels (default: empty)\n"
      << "  --below b1 b2..  charges of the row below the pixels (default: empty)\n"
      << "  --left L --right R  charge of the pixel to the left / right of the list\n"
      << "Prints p_mn of the first pair, the chain steps and the label returned by the\n"
      << "classifier for selection = chain (published) and tree (notebook).\n";
}

// ----------------------------------------------------------------------------
// LabelString
//   Pattern digits as {3,1}; "none" when the scan found no pattern.
// ----------------------------------------------------------------------------
static std::string LabelString(const std::vector<int>& q) {
  std::string s = "{";
  for (std::size_t i = 0; i < q.size(); ++i) s += (i ? "," : "") + std::to_string(q[i]);
  return s + "}";
}

// ----------------------------------------------------------------------------
// main
//   Build a 3-row image whose middle row holds the given pixels (padded with
//   empty pixels), run the classifier with each selection and print the result
//   next to the statistics that decide it. Returns 0 on success, 1 on bad input.
// ----------------------------------------------------------------------------
int main(int argc, char** argv) {
  std::vector<double> q;       // pixels of the middle row
  std::vector<double> above, below;
  double sigma = 0.16;         // readout noise (e-)
  double qmax = 5.6;           // maximum pattern charge (e-)
  double left = 0.0, right = 0.0;
  int collecting = 0;          // 0 = pixels, 1 = above, 2 = below
  for (int i = 1; i < argc; ++i) {
    const std::string t = argv[i];
    if (t == "-h" || t == "--help") { PrintUsage(argv[0]); return 0; }
    if (t == "--sigma" && i + 1 < argc) { sigma = std::atof(argv[++i]); collecting = 0; continue; }
    if (t == "--bound" && i + 1 < argc) { sigma = std::atof(argv[++i]) / 3.75; collecting = 0; continue; }
    if (t == "--qmax" && i + 1 < argc) { qmax = std::atof(argv[++i]); collecting = 0; continue; }
    if (t == "--left" && i + 1 < argc) { left = std::atof(argv[++i]); collecting = 0; continue; }
    if (t == "--right" && i + 1 < argc) { right = std::atof(argv[++i]); collecting = 0; continue; }
    if (t == "--above") { collecting = 1; continue; }
    if (t == "--below") { collecting = 2; continue; }
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (end == t.c_str() || *end != '\0') { std::cerr << "Unrecognised argument: " << t << "\n"; return 1; }
    (collecting == 1 ? above : collecting == 2 ? below : q).push_back(v);
  }
  if (q.size() < 2 || sigma <= 0.0) { PrintUsage(argv[0]); return 1; }

  // Layout of the middle row: 3 empty pixels, optional left pixel, the list, optional right pixel, 3 empty.
  const int pad = 3;
  std::vector<double> mid(pad, 0.0);
  mid.push_back(left);
  const int first = static_cast<int>(mid.size());
  mid.insert(mid.end(), q.begin(), q.end());
  mid.push_back(right);
  mid.insert(mid.end(), pad, 0.0);
  std::vector<double> up(mid.size(), 0.0), dn(mid.size(), 0.0);
  for (std::size_t k = 0; k < above.size() && first + k < up.size(); ++k) up[first + k] = above[k];
  for (std::size_t k = 0; k < below.size() && first + k < dn.size(); ++k) dn[first + k] = below[k];

  PatternClassifierConfig cfg;
  cfg.sigma_res_e = sigma;
  cfg.Qmin_e = 3.75 * sigma;
  cfg.neighbor_Qmax_e = cfg.Qmin_e;
  cfg.Qmax_e = qmax;
  cfg.thr_M = 3.5;
  cfg.thr_MN = 4.0;
  cfg.thr_MNL = 5.5;
  cfg.max_e_per_pixel = 5;
  cfg.enable_MN = true;
  cfg.enable_MNL = true;
  cfg.chain_first_only = false;

  PatternClassifier stat(cfg);
  std::cout << std::fixed << std::setprecision(4);
  std::cout << "sigma = " << sigma << " e-   Qmin = " << cfg.Qmin_e << " e-   Qmax = " << qmax
            << " e-   pixels =";
  double total = 0.0;
  for (double v : q) { std::cout << " " << std::setprecision(2) << v; total += v; }
  std::cout << std::setprecision(4) << "   total = " << total << " e-\n\n";

  const double a = q[0], b = q[1];
  std::cout << "p_mn = -log[cdf(q1;m) cdf(q2;n)] for the first pair (" << std::setprecision(2) << a << ", " << b
            << ")\n" << std::setprecision(2);
  const int hyp[][2] = {{1, 1}, {2, 1}, {1, 2}, {3, 1}, {1, 3}, {2, 2}, {3, 2}, {2, 3}};
  for (const auto& h : hyp)
    std::cout << "   (" << h[0] << "," << h[1] << ")  " << std::setw(8) << stat.StatPairOrdered(a, b, h[0], h[1]) << "\n";
  const double p11 = stat.StatPairOrdered(a, b, 1, 1);
  const double p21 = stat.StatPair(a, b, 2, 1);
  const double p31 = stat.StatPair(a, b, 3, 1);
  const double p22 = stat.StatPairOrdered(a, b, 2, 2);
  std::cout << "chain: p_11 = " << p11 << "  p_21 = min(p_21,p_12) = " << p21
            << "  p_31 = min(p_31,p_13) = " << p31 << "  p_22 = " << p22
            << "   (pairs pass below " << cfg.thr_MN << ", triples below " << cfg.thr_MNL << ")\n";
  if (q.size() >= 3) {
    std::cout << "triple (first three pixels): p_111 = " << stat.StatTriple(q[0], q[1], q[2], 1, 1, 1)
              << "  p_211 = " << std::min(stat.StatTriple(q[0], q[1], q[2], 2, 1, 1),
                                          stat.StatTriple(q[0], q[1], q[2], 1, 2, 1))
              << "\n";
  }
  std::cout << "\n";

  for (const std::string sel : {"chain", "tree"}) {
    cfg.selection = sel;
    PatternClassifier clf(cfg);
    const auto res = clf.ScanAllImage2DWithIsolation(up, mid, dn, -1.0);
    std::cout << "selection " << std::left << std::setw(6) << sel << ": ";
    if (res.empty()) std::cout << "no pattern\n";
    else {
      for (std::size_t k = 0; k < res.size(); ++k)
        std::cout << (k ? "  " : "") << LabelString(res[k].label.q) << " (start pixel "
                  << (res[k].start_idx - first) << ", total " << std::setprecision(2) << res[k].total_charge_e
                  << " e-)";
      std::cout << "\n";
    }
    std::cout << std::right;
  }
  return 0;
}
