// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: ccdarksens_compute_brc.cc
//  Computes the random-coincidence background B^rc_p for each identified
//  pattern p using the formula from Table I of the DAMIC-M SRDM paper:
//
//    B^rc_p = N_total × Σ_q [ Π_i Poisson(q_i | λ_{q_i}) ] × B[p|q]
//
//  where:
//    q      = injected pattern (charge tuple, one entry per pixel)
//    λ_k    = per-pixel Poisson mean for k-electron dark current events
//               (computed from θ_k: λ_k = θ_k × 1e-4)
//    B[p|q] = background confusion matrix (from Background_efficiencies.csv)
//    N_total = total number of pixels = sum(Nusedpix) over all images
//
//  The six reported patterns are: {11}, {21}, {111}, {31}, {22}, {211}
//  (codes 11, 21, 111, 31, 22, 211 in the CSV's iden_pat column).
//
//  Usage:
//    ccdarksens_compute_brc [theta1] [theta2] [theta3] [N_total]
//                           [beff_csv] [image_csv]
//
//  Defaults:
//    theta1 = 3.0   (→ lambda1 = 3e-4 e-/pixel)
//    theta2 = 2.0   (→ lambda2 = 2e-4 e-/pixel)
//    theta3 = 2.0   (→ lambda3 = 2e-4 e-/pixel)
//    N_total computed from image_csv (sum of Nusedpix), or 1853807441 if not found
//    beff_csv = data/Background_efficiencies.csv
//    image_csv = data/Final_Combined_Image_Data.csv
// ===========================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Poisson PMF: P(k | lambda) = e^{-lambda} * lambda^k / k!
// Returns 0 for k<0 or lambda<0.
// ---------------------------------------------------------------------------
static double poisson_pmf(int k, double lambda) {
    if (k < 0 || lambda < 0.0) return 0.0;
    if (k == 0) return std::exp(-lambda);
    // log P = -lambda + k*log(lambda) - log(k!)
    double logp = -lambda + k * std::log(lambda);
    for (int i = 1; i <= k; ++i) logp -= std::log(static_cast<double>(i));
    return std::exp(logp);
}

// ---------------------------------------------------------------------------
// Load Background_efficiencies.csv → B[iden_code][inj_code]
// Header row: iden_pat,eff_0,eff_11,...
// Data rows:  iden_code,val0,val1,...
// ---------------------------------------------------------------------------
static std::map<int, std::map<int,double>>
load_beff(const std::string& path, std::vector<int>& inj_codes)
{
    std::map<int, std::map<int,double>> B;
    std::ifstream f(path);
    if (!f) { std::cerr << "[brc] Cannot open: " << path << "\n"; return B; }

    std::string line;
    // Header
    if (!std::getline(f, line)) return B;
    {
        std::istringstream ss(line);
        std::string tok;
        std::getline(ss, tok, ',');  // "iden_pat"
        inj_codes.clear();
        while (std::getline(ss, tok, ','))
            try { inj_codes.push_back(std::stoi(tok.substr(4))); } catch (...) {}
    }

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        if (!std::getline(ss, tok, ',')) continue;
        int iden = 0;
        try { iden = std::stoi(tok); } catch (...) { continue; }
        for (int inj : inj_codes) {
            if (!std::getline(ss, tok, ',')) break;
            try { B[iden][inj] = std::stod(tok); } catch (...) {}
        }
    }
    return B;
}

// ---------------------------------------------------------------------------
// Load Final_Combined_Image_Data.csv → sum of Nusedpix column
// ---------------------------------------------------------------------------
static long long load_ntotal(const std::string& path)
{
    std::ifstream f(path);
    if (!f) return -1;
    std::string line;
    if (!std::getline(f, line)) return -1;

    // Find column index of "Nusedpix"
    std::vector<std::string> headers;
    {
        std::istringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) headers.push_back(tok);
    }
    int col = -1;
    for (int i = 0; i < (int)headers.size(); ++i)
        if (headers[i] == "Nusedpix") { col = i; break; }
    if (col < 0) return -1;

    long long total = 0;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string tok;
        for (int i = 0; i <= col; ++i)
            if (!std::getline(ss, tok, ',')) { tok = ""; break; }
        try { total += std::stoll(tok); } catch (...) {}
    }
    return total;
}

// ---------------------------------------------------------------------------
// Decode integer pattern code → digit list (e.g. 211 → {2,1,1})
// Special: code 0 → {0}
// ---------------------------------------------------------------------------
static std::vector<int> decode(int code) {
    if (code <= 0) return {0};
    std::vector<int> digits;
    while (code) { digits.insert(digits.begin(), code % 10); code /= 10; }
    return digits;
}

// ---------------------------------------------------------------------------
// All injected patterns that appear as columns in the CSV
// These are the same 26 entries from build_injected_patterns() in the
// validate_background_efficiency app.
// ---------------------------------------------------------------------------
static std::vector<int> all_injected_codes()
{
    std::vector<int> v;
    v.push_back(0);
    for (int m = 1; m <= 5; ++m) v.push_back(m);
    for (int m = 1; m <= 5; ++m)
        for (int n = 1; n <= 5; ++n)
            if (m + n <= 5) v.push_back(m*10+n);
    for (int m = 1; m <= 5; ++m)
        for (int n = 1; n <= 5; ++n)
            for (int l = 1; l <= 5; ++l)
                if (m+n+l <= 5) v.push_back(m*100+n*10+l);
    return v;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // -----------------------------------------------------------------------
    // Parameters
    // -----------------------------------------------------------------------
    double theta1 = 3.0, theta2 = 2.0, theta3 = 2.0;
    long long N_total_arg = -1;
    std::string beff_csv  = "data/Background_efficiencies.csv";
    std::string image_csv = "data/Final_Combined_Image_Data.csv";

    if (argc >= 2) { try { theta1 = std::stod(argv[1]); } catch (...) {} }
    if (argc >= 3) { try { theta2 = std::stod(argv[2]); } catch (...) {} }
    if (argc >= 4) { try { theta3 = std::stod(argv[3]); } catch (...) {} }
    if (argc >= 5) { try { N_total_arg = std::stoll(argv[4]); } catch (...) {} }
    if (argc >= 6) beff_csv  = argv[5];
    if (argc >= 7) image_csv = argv[6];

    // λ_k = θ_k × 1e-4
    const double lam1 = theta1 * 1e-4;
    const double lam2 = theta2 * 1e-4;
    const double lam3 = theta3 * 1e-4;

    // λ for k-electron pixel: k=0→0, k=1→lam1, k=2→lam2, k≥3→lam3
    auto lambda_for = [&](int k) -> double {
        if (k <= 0) return 0.0;
        if (k == 1) return lam1;
        if (k == 2) return lam2;
        return lam3;  // k >= 3
    };

    // -----------------------------------------------------------------------
    // N_total
    // -----------------------------------------------------------------------
    long long N_total = (N_total_arg > 0) ? N_total_arg : load_ntotal(image_csv);
    if (N_total <= 0) {
        N_total = 1853807441LL;  // fallback from paper
        std::cerr << "[brc] Could not load N_total from CSV; using hardcoded "
                  << N_total << "\n";
    }

    // -----------------------------------------------------------------------
    // Load B[iden|inj]
    // -----------------------------------------------------------------------
    std::vector<int> inj_codes;
    auto B = load_beff(beff_csv, inj_codes);
    if (B.empty()) { std::cerr << "[brc] Failed to load B matrix. Aborting.\n"; return 1; }

    // -----------------------------------------------------------------------
    // Print header
    // -----------------------------------------------------------------------
    std::cout << "============================================================\n";
    std::cout << "  B^rc_p computation\n";
    std::cout << "  Formula: B^rc_p = N_total × Σ_q [ Π Poisson(q_i|λ_i) ] × B[p|q]\n";
    std::cout << "============================================================\n";
    std::cout << "  θ = (" << theta1 << ", " << theta2 << ", " << theta3 << ")\n";
    std::cout << "  λ = (" << lam1   << ", " << lam2   << ", " << lam3   << ") e-/pixel\n";
    std::cout << "  N_total = " << N_total << " pixels\n";
    std::cout << "  B matrix: " << beff_csv << "\n";
    std::cout << "------------------------------------------------------------\n\n";

    // -----------------------------------------------------------------------
    // Compute B^rc_p for each identified pattern
    // The 6 patterns reported in the paper (Table I):
    // -----------------------------------------------------------------------
    const std::vector<int> report_pats = {11, 21, 111, 31, 22, 211};

    // All injected pattern codes (columns of B matrix)
    auto all_inj = all_injected_codes();

    std::cout << std::left << std::setw(8) << "pattern"
              << std::right << std::setw(16) << "B^rc_p"
              << "\n";
    std::cout << std::string(26, '-') << "\n";

    for (int p : report_pats) {
        double brc = 0.0;

        for (int q_code : all_inj) {
            // B[p|q] — probability that injected pattern q is identified as p
            auto it1 = B.find(p);
            if (it1 == B.end()) continue;
            auto it2 = it1->second.find(q_code);
            if (it2 == it1->second.end()) continue;
            double Bpq = it2->second;
            if (Bpq == 0.0) continue;

            // P(q) = Π_i Poisson(q_i | λ_{q_i})
            std::vector<int> digits = decode(q_code);

            // For the (0,) pattern: a single sub-threshold pixel; its
            // Poisson weight is the probability of seeing ≥1 pixel with 0
            // charge above threshold — effectively weight 1 (it's the null
            // background level). The notebook uses P(0|0)=1 for this entry.
            double Pq = 1.0;
            for (int qi : digits) {
                double lam = lambda_for(qi);
                Pq *= poisson_pmf(qi, lam);
            }

            brc += Pq * Bpq;
        }

        brc *= static_cast<double>(N_total);

        std::cout << std::left  << std::setw(8) << p
                  << std::right << std::fixed << std::setprecision(4)
                  << std::setw(16) << brc
                  << "\n";
    }

    // -----------------------------------------------------------------------
    // Also print a detailed breakdown for the dominant pattern {11}
    // -----------------------------------------------------------------------
    std::cout << "\n------------------------------------------------------------\n";
    std::cout << "  Detailed breakdown for pattern p={11}:\n";
    std::cout << "  (showing injected patterns q with contribution > 0.001)\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << std::left  << std::setw(8)  << "q"
              << std::right << std::setw(12) << "P(q)"
              << std::setw(12) << "B[11|q]"
              << std::setw(16) << "contribution"
              << "\n";
    std::cout << std::string(48, '-') << "\n";

    double brc_11_total = 0.0;
    for (int q_code : all_inj) {
        auto it1 = B.find(11);
        if (it1 == B.end()) continue;
        auto it2 = it1->second.find(q_code);
        if (it2 == it1->second.end()) continue;
        double Bpq = it2->second;
        if (Bpq == 0.0) continue;

        std::vector<int> digits = decode(q_code);
        double Pq = 1.0;
        for (int qi : digits) Pq *= poisson_pmf(qi, lambda_for(qi));

        double contrib = Pq * Bpq * static_cast<double>(N_total);
        brc_11_total += contrib;

        if (contrib > 0.001) {
            std::cout << std::left  << std::setw(8) << q_code
                      << std::right << std::scientific << std::setprecision(3)
                      << std::setw(12) << Pq
                      << std::setw(12) << Bpq
                      << std::fixed << std::setprecision(4)
                      << std::setw(16) << contrib
                      << "\n";
        }
    }
    std::cout << std::string(48, '-') << "\n";
    std::cout << std::left << std::setw(32) << "  Total B^rc_{11}:"
              << std::right << std::fixed << std::setprecision(4)
              << std::setw(16) << brc_11_total << "\n";

    std::cout << "\n============================================================\n";
    std::cout << "  Reference (Table I): B^rc_{11} = 141.4\n";
    std::cout << "============================================================\n";

    return 0;
}
