// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_validate_pattern_efficiency.cc
//  Validates EfficiencyMC against the collaboration reference:
//    data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv
//
//  Runs EfficiencyMC (2D image path, matching Python generate_image_E) and
//  computes P(pattern | n_e) for n_e = 1..5.  Prints a side-by-side comparison
//  with the reference CSV and reports pull = |diff| / σ_MC.
//
//  Usage:
//    ccdarksens_validate_pattern_efficiency [n_trials_per_ne] [ref_csv_path]
//
//  Defaults: n_trials = 200000, ref = data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv
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
#include <tuple>
#include <vector>

#include <TH1D.h>

#include "ccdarksens/response/ChargeTransport.hh"
#include "ccdarksens/response/EfficiencyMC.hh"
#include "ccdarksens/response/PatternClassifier.hh"
#include "ccdarksens/response/PatternImageGenerator.hh"

using namespace ccdarksens;

// ---------------------------------------------------------------------------
// Load reference CSV → map[pattern_code][n_e] = efficiency
// ---------------------------------------------------------------------------
static std::map<int, std::map<int, double>>
load_reference_csv(const std::string& path)
{
    std::map<int, std::map<int, double>> ref;
    std::ifstream f(path);
    if (!f) { std::cerr << "[validate] Cannot open reference: " << path << "\n"; return ref; }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("pattern", 0) == 0) continue;  // header
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
        ref[pat][ne] = eff;
    }
    return ref;
}

// ---------------------------------------------------------------------------
// Map PatternLabel → integer code (e.g. {2,1} → 21)
// ---------------------------------------------------------------------------
static int label_to_code(const PatternLabel& lab) {
    int code = 0;
    for (int d : lab.q) code = code * 10 + d;
    return code;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // -----------------------------------------------------------------------
    // Parameters — match the reference CSV header exactly:
    //   Nsims=1000000, alpha=1, beta=0, A=803.25, b=0.00065, lambda=0.00041
    //   read_noise=0.16, qmin=3.75*0.16=0.60, qmax=5.6
    //   thrm=3.5, thrmn=4.0, thrmnl=5.5
    //   poison=True (DC pileup ON)
    // -----------------------------------------------------------------------
    const int ne_min = 1;
    const int ne_max = 5;

    int ne_trials = 200000;
    if (argc >= 2) { try { ne_trials = std::stoi(argv[1]); } catch (...) {} }

    std::string ref_csv = "data/Efficiencies_patterns_Nsims1000000_DCTrue_alpha1.csv";
    if (argc >= 3) ref_csv = argv[2];

    // Physics
    const double A_um2        = 803.25;
    const double b_umInv      = 0.00065;
    const double alpha_diff   = 1.0;
    const double beta_per_keV = 0.0;
    const double Ee_eV        = 0.0;       // reference uses E=0
    const double pixel_um     = 15.0;
    const double thickness_um = 669.0;     // z_max in Python

    // Noise / thresholds
    const double sigma_ro  = 0.16;
    const double lambda_dc = 0.00041;      // e-/pixel/exposure, DC pileup ON
    const double Qmin_e    = 0.60;         // 3.75 * 0.16
    const double Qmax_e    = 5.6;          // Python default qmax
    const double thr_M     = 3.5;
    const double thr_MN    = 4.0;
    const double thr_MNL   = 5.5;

    const uint64_t rng_seed = 987654321ULL;

    // Pattern codes present in the reference CSV (in display order)
    const std::vector<int> ref_pattern_codes = {
        1, 2, 3, 4, 5,
        11, 41, 32, 311, 31,
        221, 22, 23, 211, 111, 21
    };

    // -----------------------------------------------------------------------
    // Build components
    // -----------------------------------------------------------------------
    ChargeTransportConfig ct_cfg;
    ct_cfg.thickness_um  = thickness_um;
    ct_cfg.A_um2         = A_um2;
    ct_cfg.b_umInv       = b_umInv;
    ct_cfg.alpha         = alpha_diff;
    ct_cfg.beta_per_keV  = beta_per_keV;
    ct_cfg.rng_seed      = rng_seed;
    auto ct = std::make_shared<ChargeTransport>(ct_cfg);

    PatternClassifierConfig pcc;
    pcc.Qmin_e          = Qmin_e;
    pcc.neighbor_Qmax_e = Qmin_e;   // isolation threshold = qmin (matches Python)
    pcc.Qmax_e          = Qmax_e;
    pcc.sigma_res_e     = sigma_ro;
    pcc.thr_M           = thr_M;
    pcc.thr_MN          = thr_MN;
    pcc.thr_MNL         = thr_MNL;
    pcc.max_e_per_pixel = 5;
    pcc.enable_MN       = true;
    pcc.enable_MNL      = true;
    auto classifier = std::make_shared<PatternClassifier>(pcc);

    // 2D image generator — matches Python generate_image_E:
    //   image shape (3*bining, 50) → binned to (3, 50), bining=100
    PatternImageConfig img_cfg;
    img_cfg.nrows_binned        = 3;
    img_cfg.ncols               = 50;
    img_cfg.row_binning         = 100;
    img_cfg.col_binning         = 1;
    img_cfg.pixel_size_um       = pixel_um;
    img_cfg.sigma_readout_e     = sigma_ro;
    img_cfg.lambda_dc           = lambda_dc;
    img_cfg.include_dark_current= true;
    img_cfg.randomize_center    = true;  // match Python x0=uniform(16,35), y0=uniform(16,85)
    img_cfg.center_margin_pix   = 15.0;
    img_cfg.rng_seed            = rng_seed;
    auto img_gen = std::make_shared<PatternImageGenerator>(img_cfg, ct);

    // EfficiencyMC config — accepted_labels = all patterns in reference
    EfficiencyMCConfig emc_cfg;
    emc_cfg.ne_trials         = ne_trials;
    emc_cfg.row_length        = 50;
    emc_cfg.pix_cfg.mode      = PixelSimMode::RowSegment;
    emc_cfg.pix_cfg.nx        = 50;
    emc_cfg.pix_cfg.ny        = 1;
    emc_cfg.pix_cfg.pixel_size_um    = pixel_um;
    emc_cfg.pix_cfg.sigma_readout_e  = sigma_ro;
    emc_cfg.pix_cfg.lambda_dc        = lambda_dc;
    emc_cfg.include_dc_pileup = true;
    emc_cfg.seed              = rng_seed;

    // Accept every pattern that appears in the reference
    const std::vector<std::vector<int>> all_q = {
        {1},{2},{3},{4},{5},
        {1,1},{4,1},{3,2},{3,1,1},{3,1},
        {2,2,1},{2,2},{2,3},{2,1,1},{1,1,1},{2,1}
    };
    for (const auto& q : all_q) {
        PatternLabel lab;
        lab.q = q;
        lab.isolated = true;
        emc_cfg.accepted_labels.push_back(lab);
    }

    auto emc = std::make_shared<EfficiencyMC>(emc_cfg, ct, classifier);
    emc->SetPatternImageGenerator(img_gen);  // use 2D image path

    // -----------------------------------------------------------------------
    // Run MC
    // -----------------------------------------------------------------------
    std::cout << "=================================================================\n";
    std::cout << "  EfficiencyMC validation vs collaboration reference\n";
    std::cout << "=================================================================\n";
    std::cout << "  Reference CSV : " << ref_csv << "\n";
    std::cout << "  MC trials/n_e : " << ne_trials << "\n";
    std::cout << "  Diffusion     : A=" << A_um2 << " µm², b=" << b_umInv << " µm⁻¹"
              << "  alpha=" << alpha_diff << "  beta=" << beta_per_keV << "\n";
    std::cout << "  Noise/DC      : σ_ro=" << sigma_ro << " e⁻"
              << "  λ_dc=" << lambda_dc << " e⁻/pix/exp (pileup ON)\n";
    std::cout << "  Thresholds    : Qmin=" << Qmin_e << "  thr_M=" << thr_M
              << "  thr_MN=" << thr_MN << "  thr_MNL=" << thr_MNL << "\n";
    std::cout << "  Image         : 3×50 (row_bin=100), pixel=" << pixel_um << " µm\n";
    std::cout << "-----------------------------------------------------------------\n";
    std::cout << "Running 2D-image EfficiencyMC...\n\n";
    std::flush(std::cout);

    emc->PrecomputeEpsilon(ne_min, ne_max, Ee_eV);
    const auto& table = emc->GetPatternTable();

    // -----------------------------------------------------------------------
    // Load reference
    // -----------------------------------------------------------------------
    auto ref = load_reference_csv(ref_csv);
    if (ref.empty()) {
        std::cerr << "[validate] Failed to load reference CSV. Aborting.\n";
        return 1;
    }
    std::cout << "Reference loaded: " << ref.size() << " patterns.\n\n";

    // -----------------------------------------------------------------------
    // Comparison table
    // -----------------------------------------------------------------------
    // Combined uncertainty: σ_comb = sqrt(σ²_cpp + σ²_ref)
    // Both C++ MC and reference use ne_trials and N_ref=1000000 trials respectively.
    const int n_ref_trials = 1000000;
    auto mc_sigma = [&](double p_ref) -> double {
        const double p = std::max(p_ref, 1e-6);
        const double var_cpp = p * (1.0 - p) / static_cast<double>(ne_trials);
        const double var_ref = p * (1.0 - p) / static_cast<double>(n_ref_trials);
        return std::sqrt(var_cpp + var_ref);
    };

    std::cout << std::left
              << std::setw(5)  << "n_e"
              << std::setw(8)  << "pat"
              << std::setw(12) << "C++ eff"
              << std::setw(12) << "ref eff"
              << std::setw(12) << "diff"
              << std::setw(8)  << "pull σ"
              << "\n";
    std::cout << std::string(57, '-') << "\n";

    int n_checks = 0, n_fail = 0;
    double max_pull = 0.0;
    double sum_sq_pull = 0.0;

    for (int ne = ne_min; ne <= ne_max; ++ne) {
        auto it_ne = table.find(ne);
        bool printed_ne = false;

        for (int code : ref_pattern_codes) {
            auto it_ref = ref.find(code);
            double ref_val = 0.0;
            if (it_ref != ref.end()) {
                auto it_ne_ref = it_ref->second.find(ne);
                if (it_ne_ref != it_ref->second.end())
                    ref_val = it_ne_ref->second;
            }

            // Sum isolated + non-isolated variants: reference Python doesn't split on this
            double cpp_val = 0.0;
            if (it_ne != table.end()) {
                for (const auto& [lab, prob] : it_ne->second) {
                    if (label_to_code(lab) == code) cpp_val += prob;
                }
            }

            // Skip truly negligible entries in both
            if (ref_val < 1e-5 && cpp_val < 1e-5) continue;

            double diff  = cpp_val - ref_val;
            double sigma = mc_sigma(ref_val);
            double pull  = std::abs(diff) / sigma;
            bool   fail  = (pull > 3.0);

            if (fail) ++n_fail;
            ++n_checks;
            max_pull    = std::max(max_pull, pull);
            sum_sq_pull += pull * pull;

            if (!printed_ne) {
                std::cout << "\n";
                printed_ne = true;
            }

            std::cout << std::left
                      << std::setw(5)  << ne
                      << std::setw(8)  << code
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
    std::cout << "\n" << std::string(57, '=') << "\n";
    std::cout << "  Entries compared : " << n_checks << "\n";
    std::cout << "  Failures >3σ     : " << n_fail   << "\n";
    std::cout << "  Max pull         : " << std::fixed << std::setprecision(1)
              << max_pull << " σ\n";
    if (n_checks > 0)
        std::cout << "  RMS pull         : " << std::setprecision(2)
                  << std::sqrt(sum_sq_pull / n_checks) << " σ\n";
    std::cout << "  (σ = sqrt(p*(1-p)/" << ne_trials << " + p*(1-p)/" << n_ref_trials << ") — combined C++ + ref uncertainty)\n";
    std::cout << std::string(57, '=') << "\n";

    if (n_fail == 0)
        std::cout << "  RESULT: PASS — all entries within 3σ of reference\n";
    else
        std::cout << "  RESULT: FAIL — " << n_fail << " entries deviate > 3σ\n";
    std::cout << std::string(57, '=') << "\n";

    return (n_fail == 0) ? 0 : 1;
}
