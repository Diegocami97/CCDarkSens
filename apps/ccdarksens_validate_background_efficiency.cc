// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_validate_background_efficiency.cc
//  Validates the background confusion matrix against the collaboration
//  reference:  data/Background_efficiencies.csv
//
//  For each injected pattern (a,b,c) the Python reference generates:
//    simulate_cluster(a,b,c) → 3×5 noisy cluster →
//    scan_image_background()  → identified pattern (including {0} = noise)
//  and records B[identified][injected] = count / Nsims.
//
//  This app reproduces that logic using:
//    PatternImageGenerator::SimulateCluster(a,b,c)  (identical to Python)
//    PatternClassifier with allow_pattern_zero=true
//    ScanAllImage2DWithIsolation on the 3×5 cluster
//    → if empty, count as pattern {0}
//
//  Usage:
//    ccdarksens_validate_background_efficiency [n_trials] [ref_csv_path]
//
//  Defaults: n_trials = 100000, ref = data/Background_efficiencies.csv
// ===========================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/PatternClassifier.hh"
#include "ccdarksens/response/PatternImageGenerator.hh"

using namespace ccdarksens;

// ---------------------------------------------------------------------------
// Encode/decode pattern codes  (e.g. {2,1,1} → 211)
// ---------------------------------------------------------------------------
static int encode_code(const std::vector<int>& q) {
    int code = 0;
    for (int d : q) code = code * 10 + d;
    return code;
}

// All injected patterns: (a,b,c) for simulate_cluster(a,b,c).
// Matches Python: all_pats = [(0,)] + combos of 1..5 up to sum≤5, size 1..3.
static std::vector<std::tuple<int,int,int,int>> build_injected_patterns() {
    // Returns {code, a, b, c} where code = integer label for the column header
    std::vector<std::tuple<int,int,int,int>> result;
    // (0,) → code=0, a=b=c=0
    result.emplace_back(0, 0, 0, 0);
    // size=1: (m,) for m=1..5
    for (int m = 1; m <= 5; ++m)
        result.emplace_back(m, m, 0, 0);
    // size=2: (m,n) for m,n in [1,5], m+n<=5
    for (int m = 1; m <= 5; ++m)
        for (int n = 1; n <= 5; ++n)
            if (m + n <= 5)
                result.emplace_back(m*10 + n, m, n, 0);
    // size=3: (m,n,l) for m,n,l in [1,5], m+n+l<=5
    for (int m = 1; m <= 5; ++m)
        for (int n = 1; n <= 5; ++n)
            for (int l = 1; l <= 5; ++l)
                if (m + n + l <= 5)
                    result.emplace_back(m*100 + n*10 + l, m, n, l);
    return result;
}

// ---------------------------------------------------------------------------
// Load reference CSV → ref[iden_code][inj_code] = efficiency
// ---------------------------------------------------------------------------
static std::map<int, std::map<int,double>>
load_background_csv(const std::string& path,
                    std::vector<int>& inj_codes_out,
                    std::vector<int>& iden_codes_out)
{
    std::map<int, std::map<int,double>> ref;
    std::ifstream f(path);
    if (!f) {
        std::cerr << "[validate_bg] Cannot open: " << path << "\n";
        return ref;
    }

    std::string line;
    // Parse header row: iden_pat,eff_0,eff_1,...
    if (!std::getline(f, line)) return ref;
    {
        std::istringstream ss(line);
        std::string tok;
        std::getline(ss, tok, ',');  // "iden_pat"
        inj_codes_out.clear();
        while (std::getline(ss, tok, ',')) {
            // tok = "eff_0", "eff_11", etc.
            if (tok.size() > 4)
                try { inj_codes_out.push_back(std::stoi(tok.substr(4))); }
                catch (...) {}
        }
    }

    iden_codes_out.clear();
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        if (!std::getline(ss, tok, ',')) continue;
        int iden_code = 0;
        try { iden_code = std::stoi(tok); } catch (...) { continue; }
        iden_codes_out.push_back(iden_code);
        for (int inj : inj_codes_out) {
            if (!std::getline(ss, tok, ',')) break;
            double val = 0.0;
            try { val = std::stod(tok); } catch (...) {}
            ref[iden_code][inj] = val;
        }
    }
    return ref;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // -----------------------------------------------------------------------
    // Parameters — must match the Python reference (same thresholds/noise)
    // -----------------------------------------------------------------------
    int n_trials = 100000;
    if (argc >= 2) { try { n_trials = std::stoi(argv[1]); } catch (...) {} }

    std::string ref_csv = "data/Background_efficiencies.csv";
    if (argc >= 3) ref_csv = argv[2];

    const double sigma_ro  = 0.16;
    const double Qmin_e    = 0.60;    // 3.75 * 0.16
    const double Qmax_e    = 5.6;
    const double thr_M     = 3.5;
    const double thr_MN    = 4.0;
    const double thr_MNL   = 5.5;
    const double pixel_um  = 15.0;
    const uint64_t rng_seed = 123456789ULL;

    // -----------------------------------------------------------------------
    // Build components
    // -----------------------------------------------------------------------
    // ChargeTransport is needed for PatternImageGenerator constructor but not
    // used during SimulateCluster (no diffusion sampling).  Provide minimal cfg.
    ChargeTransportConfig ct_cfg;
    ct_cfg.thickness_um = 669.0;
    ct_cfg.A_um2        = 803.25;
    ct_cfg.b_umInv      = 0.00065;
    ct_cfg.alpha        = 1.0;
    ct_cfg.rng_seed     = rng_seed;
    auto ct = std::make_shared<ChargeTransport>(ct_cfg);

    // PatternImageGenerator — 3×5 config, just for SimulateCluster
    PatternImageConfig img_cfg;
    img_cfg.nrows_binned    = 3;
    img_cfg.ncols           = 5;
    img_cfg.row_binning     = 1;   // not used by SimulateCluster
    img_cfg.col_binning     = 1;
    img_cfg.pixel_size_um   = pixel_um;
    img_cfg.sigma_readout_e = sigma_ro;
    img_cfg.rng_seed        = rng_seed;
    auto img_gen = std::make_shared<PatternImageGenerator>(img_cfg, ct);

    // PatternClassifier with allow_pattern_zero=true (background mode)
    PatternClassifierConfig pcc;
    pcc.Qmin_e          = Qmin_e;
    pcc.neighbor_Qmax_e = Qmin_e;
    pcc.Qmax_e          = Qmax_e;
    pcc.sigma_res_e     = sigma_ro;
    pcc.thr_M           = thr_M;
    pcc.thr_MN          = thr_MN;
    pcc.thr_MNL         = thr_MNL;
    pcc.max_e_per_pixel = 5;
    pcc.enable_MN       = true;
    pcc.enable_MNL      = true;
    pcc.allow_pattern_zero = true;   // ← key: identifies sub-threshold as pattern {0}
    auto classifier = std::make_shared<PatternClassifier>(pcc);

    // -----------------------------------------------------------------------
    // Load reference
    // -----------------------------------------------------------------------
    std::vector<int> ref_inj_codes, ref_iden_codes;
    auto ref = load_background_csv(ref_csv, ref_inj_codes, ref_iden_codes);
    if (ref.empty()) {
        std::cerr << "[validate_bg] Failed to load reference CSV.\n";
        return 1;
    }
    std::cout << "=================================================================\n";
    std::cout << "  Background efficiency validation vs collaboration reference\n";
    std::cout << "=================================================================\n";
    std::cout << "  Reference CSV : " << ref_csv << "\n";
    std::cout << "  MC trials/pat : " << n_trials << "\n";
    std::cout << "  Noise         : σ_ro=" << sigma_ro << " e⁻\n";
    std::cout << "  Thresholds    : Qmin=" << Qmin_e << "  thr_M=" << thr_M
              << "  thr_MN=" << thr_MN << "  thr_MNL=" << thr_MNL << "\n";
    std::cout << "  N_ref_trials  : 100000 (deduced from min non-zero value 1e-5)\n";
    std::cout << "-----------------------------------------------------------------\n";

    const int n_ref_trials = 100000;

    // -----------------------------------------------------------------------
    // Run background MC for each injected pattern
    // -----------------------------------------------------------------------
    auto injected_pats = build_injected_patterns();

    // cpp_counts[inj_code][iden_code] = count
    std::map<int, std::map<int, std::size_t>> cpp_counts;

    for (const auto& [inj_code, a, b, c] : injected_pats) {
        auto& counts = cpp_counts[inj_code];
        counts[0] = 0;  // ensure {0} entry exists

        for (int trial = 0; trial < n_trials; ++trial) {
            // simulate_cluster(a, b, c) → 3×5 noisy cluster
            auto cluster = img_gen->SimulateCluster(
                static_cast<double>(a),
                static_cast<double>(b),
                static_cast<double>(c));

            const std::vector<double>& row_above  = cluster[0];
            const std::vector<double>& row_middle = cluster[1];
            const std::vector<double>& row_below  = cluster[2];

            // scan_image_background: collect all patterns, else {0}
            auto prs = classifier->ScanAllImage2DWithIsolation(
                row_above, row_middle, row_below, -1.0);

            if (prs.empty()) {
                counts[0] += 1;
            } else {
                for (const auto& pr : prs) {
                    int code = encode_code(pr.label.q);
                    counts[code] += 1;
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Comparison
    // -----------------------------------------------------------------------
    // Pull: diff / sqrt(σ²_cpp + σ²_ref)
    auto combined_sigma = [&](double p_ref) -> double {
        const double p = std::max(p_ref, 1e-6);
        const double v_cpp = p * (1.0 - p) / static_cast<double>(n_trials);
        const double v_ref = p * (1.0 - p) / static_cast<double>(n_ref_trials);
        return std::sqrt(v_cpp + v_ref);
    };

    std::cout << "\n";
    std::cout << std::left
              << std::setw(10) << "iden_pat"
              << std::setw(10) << "inj_pat"
              << std::setw(12) << "C++ eff"
              << std::setw(12) << "ref eff"
              << std::setw(12) << "diff"
              << std::setw(8)  << "pull"
              << "\n";
    std::cout << std::string(64, '-') << "\n";

    int n_checks = 0, n_fail = 0;
    double max_pull = 0.0, sum_sq = 0.0;

    // Only compare entries where ref_val >= 1e-4 to avoid comparing negligible entries
    const double min_ref_val = 1e-4;

    for (int inj_code : ref_inj_codes) {
        const auto& counts = cpp_counts[inj_code];
        bool first_in_col = true;

        for (int iden_code : ref_iden_codes) {
            double ref_val = 0.0;
            {
                auto it1 = ref.find(iden_code);
                if (it1 != ref.end()) {
                    auto it2 = it1->second.find(inj_code);
                    if (it2 != it1->second.end()) ref_val = it2->second;
                }
            }

            std::size_t cnt = 0;
            {
                auto it = counts.find(iden_code);
                if (it != counts.end()) cnt = it->second;
            }
            double cpp_val = static_cast<double>(cnt) / static_cast<double>(n_trials);

            // Skip very small entries to keep output manageable
            if (ref_val < min_ref_val && cpp_val < min_ref_val) continue;

            double diff  = cpp_val - ref_val;
            double sigma = combined_sigma(ref_val);
            double pull  = std::abs(diff) / sigma;
            bool   fail  = (pull > 3.0);

            ++n_checks;
            if (fail) ++n_fail;
            max_pull  = std::max(max_pull, pull);
            sum_sq   += pull * pull;

            if (first_in_col) { std::cout << "\n"; first_in_col = false; }

            std::cout << std::left
                      << std::setw(10) << iden_code
                      << std::setw(10) << inj_code
                      << std::fixed << std::setprecision(5)
                      << std::setw(12) << cpp_val
                      << std::setw(12) << ref_val
                      << std::setw(12) << diff
                      << std::setprecision(1)
                      << std::setw(8)  << pull
                      << (fail ? "  *** FAIL" : "")
                      << "\n";
        }
    }

    // -----------------------------------------------------------------------
    // Summary
    // -----------------------------------------------------------------------
    std::cout << "\n" << std::string(64, '=') << "\n";
    std::cout << "  Entries compared (ref >= " << min_ref_val << ") : " << n_checks << "\n";
    std::cout << "  Failures >3σ              : " << n_fail   << "\n";
    std::cout << "  Max pull                  : " << std::fixed << std::setprecision(1)
              << max_pull << " σ\n";
    if (n_checks > 0)
        std::cout << "  RMS pull                  : " << std::setprecision(2)
                  << std::sqrt(sum_sq / n_checks) << " σ\n";
    std::cout << "  (σ = combined C++ + ref,  N_cpp=" << n_trials
              << ",  N_ref=" << n_ref_trials << ")\n";
    std::cout << std::string(64, '=') << "\n";

    if (n_fail == 0)
        std::cout << "  RESULT: PASS — all entries within 3σ of reference\n";
    else
        std::cout << "  RESULT: FAIL — " << n_fail << " entries deviate > 3σ\n";
    std::cout << std::string(64, '=') << "\n";

    return (n_fail == 0) ? 0 : 1;
}
