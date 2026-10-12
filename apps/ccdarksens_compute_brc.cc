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
//    ccdarksens_compute_brc --config rates.json [--lambda1 L] [--lambda2 L] [--lambda3 L]
//                           [--n-pixels N] [--segments table.csv] [--matrix B.csv] [--out f.json]
//
//  Defaults:
//    theta1 = 3.0   (→ lambda1 = 3e-4 e-/pixel)
//    theta2 = 2.0   (→ lambda2 = 2e-4 e-/pixel)
//    theta3 = 2.0   (→ lambda3 = 2e-4 e-/pixel)
//    N_total computed from image_csv (sum of Nusedpix), or 1853807441 if not found
//    beff_csv = data/Background_efficiencies.csv
//    image_csv = data/Final_Combined_Image_Data.csv
//
//  The config is the file of configs/examples/background_pattern_rates.json (blocks
//  "rates", "pixels", "migration_matrix_csv", "patterns", "breakdown_pattern",
//  "output_json"); a flag overrides the matching entry. A segments table (columns
//  n_pixels, lambda_1 and optionally lambda_2, lambda_3plus; one row per image, CCD or
//  column group) makes the sum run over segments with their own rates, as in Eq. (4) of
//  the SRDM paper. The six B^rc_p are also printed as a background_Bp line for a scan
//  config and written to --out (or "output_json") as JSON.
// ===========================================================================

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

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
// Segment: pixel readouts that share the same per-pixel rates (an image, a CCD,
// a group of columns). The sum over segments is Eq. (4) of the SRDM paper; one
// segment is the uniform-rate formula of pydme's Background_pattern_.
// ---------------------------------------------------------------------------
struct Segment {
    double n_pixels = 0.0;   // pixel readouts in the segment
    double lam[4] = {0.0, 0.0, 0.0, 0.0};  // mean number of k-electron events per readout, k = 1, 2, >=3
};

// ---------------------------------------------------------------------------
// load_segments
//   Read a CSV with columns n_pixels and lambda_1 (and optionally lambda_2,
//   lambda_3plus); a missing rate column falls back to the default rates.
// ---------------------------------------------------------------------------
static std::vector<Segment>
load_segments(const std::string& path, const double def[4])
{
    std::vector<Segment> segs;
    std::ifstream f(path);
    if (!f) { std::cerr << "[brc] Cannot open: " << path << "\n"; return segs; }
    std::string line;
    if (!std::getline(f, line)) return segs;
    std::vector<std::string> hdr;
    {
        std::istringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) hdr.push_back(tok);
    }
    auto col = [&](const char* name) {
        for (std::size_t i = 0; i < hdr.size(); ++i) if (hdr[i] == name) return static_cast<int>(i);
        return -1;
    };
    const int c_n = col("n_pixels"), c1 = col("lambda_1"), c2 = col("lambda_2"), c3 = col("lambda_3plus");
    if (c_n < 0) { std::cerr << "[brc] segments table needs an n_pixels column\n"; return segs; }
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> v;
        std::istringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) v.push_back(tok);
        auto get = [&](int c, double fallback) {
            if (c < 0 || c >= static_cast<int>(v.size()) || v[c].empty()) return fallback;
            try { return std::stod(v[c]); } catch (...) { return fallback; }
        };
        Segment s;
        s.n_pixels = get(c_n, 0.0);
        s.lam[1] = get(c1, def[1]);
        s.lam[2] = get(c2, def[2]);
        s.lam[3] = get(c3, def[3]);
        segs.push_back(s);
    }
    return segs;
}

// ---------------------------------------------------------------------------
// pattern_weight
//   P(q) = Π_i Poisson(q_i | λ_{q_i}) for the ideal pattern code q in a segment;
//   the empty pattern has weight 1, and k >= 3 electrons use the 3+ rate.
// ---------------------------------------------------------------------------
static double pattern_weight(int q_code, const Segment& seg)
{
    double w = 1.0;
    for (int qi : decode(q_code)) {
        const double lam = (qi <= 0) ? 0.0 : seg.lam[std::min(qi, 3)];
        w *= poisson_pmf(qi, lam);
    }
    return w;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // -----------------------------------------------------------------------
    // Parameters: positional values keep the original meaning, flags override
    // -----------------------------------------------------------------------
    double theta1 = 3.0, theta2 = 2.0, theta3 = 2.0;
    long long N_total_arg = -1;
    std::string beff_csv  = "data/Background_efficiencies.csv";
    std::string image_csv = "data/Final_Combined_Image_Data.csv";
    std::string config_path, segments_path, out_json;
    std::vector<int> report_pats = {11, 21, 111, 31, 22, 211};
    int breakdown_pat = 11;
    double flag_l1 = -1, flag_l2 = -1, flag_l3 = -1;
    std::string flag_matrix;          // --matrix, applied after the config so a flag wins
    long long flag_npix = -1;         // --n-pixels, same

    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        const std::string t = argv[i];
        auto next = [&](double& dst) { if (i + 1 < argc) dst = std::atof(argv[++i]); };
        if (t == "--config" && i + 1 < argc) config_path = argv[++i];
        else if (t == "--segments" && i + 1 < argc) segments_path = argv[++i];
        else if (t == "--matrix" && i + 1 < argc) flag_matrix = argv[++i];
        else if (t == "--out" && i + 1 < argc) out_json = argv[++i];
        else if (t == "--lambda1") next(flag_l1);
        else if (t == "--lambda2") next(flag_l2);
        else if (t == "--lambda3") next(flag_l3);
        else if (t == "--n-pixels" && i + 1 < argc) flag_npix = static_cast<long long>(std::atof(argv[++i]));
        else if (t == "-h" || t == "--help") {
            std::cout << "Usage: ccdarksens_compute_brc [theta1] [theta2] [theta3] [N_total] [beff_csv] [image_csv]\n"
                      << "       ccdarksens_compute_brc --config rates.json [--lambda1 L] [--lambda2 L] [--lambda3 L]\n"
                      << "                              [--n-pixels N] [--segments table.csv] [--matrix B.csv] [--out f.json]\n"
                      << "  A flag overrides the matching entry of the config.\n";
            return 0;
        }
        else if (t.rfind("--", 0) == 0) { std::cerr << "[brc] unknown option " << t << "\n"; return 1; }
        else positional.push_back(t);
    }

    if (!config_path.empty()) {
        std::ifstream cf(config_path);
        if (!cf) { std::cerr << "[brc] Cannot open config: " << config_path << "\n"; return 1; }
        nlohmann::json j;
        cf >> j;
        if (j.contains("rates")) {
            theta1 = j["rates"].value("lambda_1", theta1 * 1e-4) * 1e4;
            theta2 = j["rates"].value("lambda_2", theta2 * 1e-4) * 1e4;
            theta3 = j["rates"].value("lambda_3plus", theta3 * 1e-4) * 1e4;
        }
        if (j.contains("pixels")) {
            const auto& px = j["pixels"];
            if (px.contains("n_pixel_readouts") && !px["n_pixel_readouts"].is_null())
                N_total_arg = static_cast<long long>(px["n_pixel_readouts"].get<double>());
            image_csv = px.value("image_table_csv", image_csv);
            if (segments_path.empty() && px.contains("segments_csv") && !px["segments_csv"].is_null())
                segments_path = px["segments_csv"].get<std::string>();
        }
        beff_csv = j.value("migration_matrix_csv", beff_csv);
        if (j.contains("patterns")) report_pats = j["patterns"].get<std::vector<int>>();
        breakdown_pat = j.value("breakdown_pattern", report_pats.empty() ? 11 : report_pats[0]);
        if (out_json.empty() && j.contains("output_json") && !j["output_json"].is_null())
            out_json = j["output_json"].get<std::string>();
    }
    if (config_path.empty()) {
        if (positional.size() >= 1) { try { theta1 = std::stod(positional[0]); } catch (...) {} }
        if (positional.size() >= 2) { try { theta2 = std::stod(positional[1]); } catch (...) {} }
        if (positional.size() >= 3) { try { theta3 = std::stod(positional[2]); } catch (...) {} }
        if (positional.size() >= 4) { try { N_total_arg = std::stoll(positional[3]); } catch (...) {} }
        if (positional.size() >= 5) beff_csv  = positional[4];
        if (positional.size() >= 6) image_csv = positional[5];
    }
    if (!flag_matrix.empty()) beff_csv = flag_matrix;
    if (flag_npix > 0) N_total_arg = flag_npix;
    if (flag_l1 >= 0) theta1 = flag_l1 * 1e4;
    if (flag_l2 >= 0) theta2 = flag_l2 * 1e4;
    if (flag_l3 >= 0) theta3 = flag_l3 * 1e4;

    // λ_k = θ_k × 1e-4
    const double lam1 = theta1 * 1e-4;
    const double lam2 = theta2 * 1e-4;
    const double lam3 = theta3 * 1e-4;
    const double lam_default[4] = {0.0, lam1, lam2, lam3};

    // -----------------------------------------------------------------------
    // Segments: one uniform segment with N_total pixels, or the rows of a table
    // -----------------------------------------------------------------------
    std::vector<Segment> segments;
    long long N_total = 0;
    if (!segments_path.empty()) {
        segments = load_segments(segments_path, lam_default);
        if (segments.empty()) { std::cerr << "[brc] no segments loaded. Aborting.\n"; return 1; }
        double n = 0.0;
        for (const auto& sg : segments) n += sg.n_pixels;
        N_total = static_cast<long long>(std::llround(n));
    } else {
        N_total = (N_total_arg > 0) ? N_total_arg : load_ntotal(image_csv);
        if (N_total <= 0) {
            N_total = 1853807441LL;  // fallback from paper
            std::cerr << "[brc] Could not load N_total from CSV; using hardcoded "
                      << N_total << "\n";
        }
        Segment one;
        one.n_pixels = static_cast<double>(N_total);
        one.lam[1] = lam1; one.lam[2] = lam2; one.lam[3] = lam3;
        segments.push_back(one);
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
    if (!segments_path.empty()) {
        std::cout << "  rates per segment from " << segments_path << " (" << segments.size() << " segments)\n";
        for (int k = 1; k <= 3; ++k) {
            double mean = 0.0;
            for (const auto& sg : segments) mean += sg.n_pixels * sg.lam[k];
            std::cout << "  pixel-weighted mean λ_" << k << (k == 3 ? "+" : " ") << " = "
                      << mean / static_cast<double>(N_total) << " e-/pixel\n";
        }
    } else {
        std::cout << "  θ = (" << theta1 << ", " << theta2 << ", " << theta3 << ")\n";
        std::cout << "  λ = (" << lam1   << ", " << lam2   << ", " << lam3   << ") e-/pixel\n";
    }
    std::cout << "  N_total = " << N_total << " pixels\n";
    std::cout << "  B matrix: " << beff_csv << "\n";
    std::cout << "------------------------------------------------------------\n\n";

    // -----------------------------------------------------------------------
    // Compute B^rc_p for each identified pattern (default: the six of Table I)
    // -----------------------------------------------------------------------
    auto all_inj = all_injected_codes();

    std::cout << std::left << std::setw(8) << "pattern"
              << std::right << std::setw(16) << "B^rc_p"
              << std::setw(16) << "diagonal q=p"
              << std::setw(16) << "misidentified"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    std::vector<double> brc_values;
    for (int p : report_pats) {
        double brc = 0.0;
        double diag = 0.0;  // q = p term: pattern p injected and identified correctly

        for (int q_code : all_inj) {
            auto it1 = B.find(p);
            if (it1 == B.end()) continue;
            auto it2 = it1->second.find(q_code);
            if (it2 == it1->second.end()) continue;
            const double Bpq = it2->second;
            if (Bpq == 0.0) continue;

            // N_s-weighted P(q) summed over the segments
            double weighted = 0.0;
            for (const auto& sg : segments) weighted += sg.n_pixels * pattern_weight(q_code, sg);

            brc += weighted * Bpq;
            if (q_code == p) diag = weighted * Bpq;
        }
        brc_values.push_back(brc);

        std::cout << std::left  << std::setw(8) << p
                  << std::right << std::scientific << std::setprecision(4)
                  << std::setw(16) << brc
                  << std::setw(16) << diag
                  << std::setw(16) << (brc - diag)
                  << "\n";
    }

    // -----------------------------------------------------------------------
    // Detailed breakdown for one pattern (default {11})
    // -----------------------------------------------------------------------
    std::cout << "\n------------------------------------------------------------\n";
    std::cout << "  Detailed breakdown for pattern p={" << breakdown_pat << "}:\n";
    std::cout << "  (showing injected patterns q with contribution > 0.001)\n";
    std::cout << "------------------------------------------------------------\n";
    std::cout << std::left  << std::setw(8)  << "q"
              << std::right << std::setw(12) << "P(q)"
              << std::setw(12) << ("B[" + std::to_string(breakdown_pat) + "|q]")
              << std::setw(16) << "contribution"
              << "\n";
    std::cout << std::string(48, '-') << "\n";

    double brc_bd_total = 0.0;
    for (int q_code : all_inj) {
        auto it1 = B.find(breakdown_pat);
        if (it1 == B.end()) continue;
        auto it2 = it1->second.find(q_code);
        if (it2 == it1->second.end()) continue;
        const double Bpq = it2->second;
        if (Bpq == 0.0) continue;

        double weighted = 0.0;
        for (const auto& sg : segments) weighted += sg.n_pixels * pattern_weight(q_code, sg);
        const double Pq = weighted / static_cast<double>(N_total);  // pixel-averaged P(q)

        const double contrib = weighted * Bpq;
        brc_bd_total += contrib;

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
    std::cout << std::left << std::setw(32) << ("  Total B^rc_{" + std::to_string(breakdown_pat) + "}:")
              << std::right << std::fixed << std::setprecision(4)
              << std::setw(16) << brc_bd_total << "\n";

    // -----------------------------------------------------------------------
    // background_Bp line for a scan config, and optional JSON file
    // -----------------------------------------------------------------------
    nlohmann::json bp = nlohmann::json::array();
    for (double v : brc_values) {
        std::ostringstream os;
        os << std::setprecision(6) << v;
        bp.push_back(std::stod(os.str()));
    }
    std::cout << "\n------------------------------------------------------------\n";
    std::cout << "  Paste into a scan config under \"run\" (pattern order";
    for (int p : report_pats) std::cout << " " << p;
    std::cout << "):\n  " << nlohmann::json{{"background_Bp", bp}}.dump() << "\n";
    if (!out_json.empty()) {
        std::ofstream of(out_json);
        if (!of) { std::cerr << "[brc] Cannot write " << out_json << "\n"; return 1; }
        nlohmann::json out;
        out["patterns"] = report_pats;
        out["n_pixel_readouts"] = N_total;
        out["background_Bp"] = bp;
        of << out.dump(2) << "\n";
        std::cout << "  Wrote " << out_json << "\n";
    }

    std::cout << "\n============================================================\n";
    std::cout << "  Published pydme Bp = [141.4, 0.111, 0.042, 0.019, 2.5e-5, 5.8e-5]\n";
    std::cout << "  for {11,21,111,31,22,211}. Its {11} entry equals the diagonal term\n";
    std::cout << "  (q = p) at theta1 = 2.90; the 'misidentified' column above is extra.\n";
    std::cout << "  The other five Bp entries are not one consistent theta-fit of this\n";
    std::cout << "  formula with the shipped B[p|q], so do not expect this app to return them.\n";
    std::cout << "============================================================\n";

    return 0;
}
