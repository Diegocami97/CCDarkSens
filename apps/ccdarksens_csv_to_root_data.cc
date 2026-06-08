// ============================================================================
//  CCDarkSens — ccdarksens_csv_to_root_data
//  Converts observed pattern-count CSVs into ROOT TH1D data histograms (D_pat) with optional exposure/mass metadata for likelihood inputs.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>

#include <TH1D.h>
#include <TFile.h>
#include <TParameter.h>

#include "ccdarksens/io/ConfigManager.hh"

// Trim whitespace
static std::string trim(const std::string& s) {
  auto start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  auto end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end == std::string::npos ? std::string::npos : end - start + 1);
}

// Parse string to double; return 0.0 on failure (e.g. date columns)
static double parse_double(const std::string& s) {
  if (s.empty()) return 0.0;
  char* end = nullptr;
  double v = std::strtod(s.c_str(), &end);
  return (end != nullptr && end != s.c_str()) ? v : 0.0;
}

// Default mass per pixel (kg): Si 15 µm × 15 µm × 0.67 mm, rho = 2.33 g/cm³ (LBC/pydme style)
static const double k_default_mass_per_pixel_kg = 3.51e-10;

// Parse pattern_roi from config; return empty if no config or no pattern_roi
static std::vector<int> get_pattern_roi_from_config(const std::string& config_path) {
  try {
    ccdarksens::ConfigManager cfg(config_path);
    cfg.parse();
    const auto& roi = cfg.experiment_cfg().pattern_roi;
    return std::vector<int>(roi.begin(), roi.end());
  } catch (...) {
    return {};
  }
}

// Optional exposure/mass from CSV (when texp and Nusedpix/Npix columns present).
struct ExposureMass {
  double exposure_kg_year = 0.0;
  double mass_kg = 0.0;           // effective mass (time-weighted), kg
  double livetime_days = 0.0;     // total livetime, days (sum of texp)
  bool has_value = false;
};

// Detect CSV format and parse into (pattern_ids, counts).
// pattern_roi_hint: if non-empty, use this order for output; otherwise use order from CSV.
// out_exposure_mass: if non-null, filled when CSV has texp and Nusedpix/Npix columns.
static bool parse_csv(const std::string& path,
                      std::vector<int>& pattern_ids,
                      std::vector<double>& counts,
                      const std::vector<int>& pattern_roi_hint,
                      ExposureMass* out_exposure_mass = nullptr) {
  std::ifstream in(path);
  if (!in) {
    std::cerr << "[csv2root] Cannot open " << path << "\n";
    return false;
  }
  pattern_ids.clear();
  counts.clear();

  std::string line;
  if (!std::getline(in, line)) return false;
  line = trim(line);
  while (!line.empty() && static_cast<unsigned char>(line[0]) < 32) line.erase(0, 1);

  // Check for header with Count_Candidate_* or "pattern" + count column
  if (line.find("Count_Candidate_") != std::string::npos ||
      (line.find("pattern") != std::string::npos && (line.find("count") != std::string::npos || line.find("Count") != std::string::npos))) {
    // Parse header
    std::vector<std::string> headers;
    std::istringstream hs(line);
    std::string tok;
    while (std::getline(hs, tok, ',')) {
      std::string t = trim(tok);
      // Strip BOM or other leading non-printable from first column
      while (!t.empty() && static_cast<unsigned char>(t[0]) < 32) t.erase(0, 1);
      headers.push_back(t);
    }
    bool has_count_candidate = false;
    for (const auto& h : headers) {
      if (h.find("Count_Candidate_") != std::string::npos) { has_count_candidate = true; break; }
    }
    if (has_count_candidate) {
      // Build list of (column_index, pattern_id) for Count_Candidate_* columns
      std::vector<std::pair<size_t, int>> pat_columns;
      const std::string prefix("Count_Candidate_");
      for (size_t i = 0; i < headers.size(); ++i) {
        std::string h = headers[i];
        auto pos = h.find(prefix);
        if (pos == std::string::npos) {
          // Also accept common typo Count_Candiate_
          pos = h.find("Count_Candiate_");
          if (pos != std::string::npos) h = "Count_Candidate_" + h.substr(14);
        }
        if (pos != std::string::npos) {
          std::string num = h.substr(pos + prefix.size());
          try {
            int pat = std::stoi(num);
            pat_columns.push_back({i, pat});
          } catch (...) {}
        }
      }
      if (pat_columns.empty()) {
        std::cerr << "[csv2root] No Count_Candidate_* columns found.\n";
        return false;
      }
      pattern_ids.clear();
      for (const auto& p : pat_columns) pattern_ids.push_back(p.second);
      counts.assign(pat_columns.size(), 0.0);

      // Find optional columns for exposure/mass (texp in days, Nusedpix or Npix)
      int texp_col = -1;
      int npix_col = -1;
      for (size_t i = 0; i < headers.size(); ++i) {
        if (headers[i] == "texp") texp_col = static_cast<int>(i);
        else if (headers[i] == "Nusedpix" || headers[i] == "Npix") npix_col = static_cast<int>(i);
      }
      double exposure_kg_days = 0.0;
      double sum_texp_days = 0.0;
      const double mass_per_pixel_kg = k_default_mass_per_pixel_kg;

      // Integrate over time: sum counts from every data row; optionally accumulate exposure
      int n_rows = 0;
      while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        std::istringstream ds(line);
        std::vector<std::string> cells;
        while (std::getline(ds, tok, ',')) cells.push_back(trim(tok));
        for (size_t k = 0; k < pat_columns.size(); ++k) {
          size_t col = pat_columns[k].first;
          if (col < cells.size()) counts[k] += parse_double(cells[col]);
        }
        if (out_exposure_mass && texp_col >= 0 && npix_col >= 0 &&
            static_cast<size_t>(texp_col) < cells.size() &&
            static_cast<size_t>(npix_col) < cells.size()) {
          double t = parse_double(cells[texp_col]);
          double n = parse_double(cells[npix_col]);
          if (t > 0 && n > 0) {
            exposure_kg_days += t * n * mass_per_pixel_kg;
            sum_texp_days += t;
          }
        }
        ++n_rows;
      }
      if (n_rows == 0) {
        std::cerr << "[csv2root] No data rows after Count_Candidate_* header.\n";
        return false;
      }
      if (n_rows > 1) {
        std::cout << "[csv2root] Summed " << n_rows << " time rows → total candidates per pattern.\n";
      }
      if (out_exposure_mass && sum_texp_days > 0 && exposure_kg_days > 0) {
        out_exposure_mass->exposure_kg_year = exposure_kg_days / 365.25;  // match ExperimentSetup (kg·year)
        out_exposure_mass->mass_kg = exposure_kg_days / sum_texp_days;
        out_exposure_mass->livetime_days = sum_texp_days;
        out_exposure_mass->has_value = true;
        std::cout << "[csv2root] From texp and Nusedpix: livetime = " << out_exposure_mass->livetime_days
                  << " days, effective mass = " << out_exposure_mass->mass_kg << " kg, exposure = "
                  << out_exposure_mass->exposure_kg_year << " kg·year.\n";
      }
    } else {
      // pattern,count style: need first data row
      if (!std::getline(in, line)) {
        std::cerr << "[csv2root] Missing data row after header.\n";
        return false;
      }
      line = trim(line);
      std::istringstream ds(line);
      std::vector<double> row;
      while (std::getline(ds, tok, ',')) {
        try {
          row.push_back(std::stod(trim(tok)));
        } catch (...) {
          row.push_back(0.0);
        }
      }
      // pattern,count style
      int pat_col = -1, count_col = -1;
      for (size_t i = 0; i < headers.size(); ++i) {
        std::string h = headers[i];
        if (h == "pattern") pat_col = static_cast<int>(i);
        else if (h == "count" || h == "Count") count_col = static_cast<int>(i);
      }
      if (pat_col < 0 || count_col < 0 || static_cast<size_t>(pat_col) >= row.size() || static_cast<size_t>(count_col) >= row.size()) {
        std::cerr << "[csv2root] Need 'pattern' and 'count' columns.\n";
        return false;
      }
      pattern_ids.push_back(static_cast<int>(row[pat_col]));
      counts.push_back(row[count_col]);
      // Read rest of rows
      while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) break;
        std::istringstream rs(line);
        std::vector<double> r;
        while (std::getline(rs, tok, ',')) {
          try { r.push_back(std::stod(trim(tok))); } catch (...) { r.push_back(0); }
        }
        if (r.size() >= static_cast<size_t>(std::max(pat_col, count_col) + 1)) {
          pattern_ids.push_back(static_cast<int>(r[pat_col]));
          counts.push_back(r[count_col]);
        }
      }
    }
  } else {
    // No header: one or two columns
    std::vector<double> first_row;
    {
      std::istringstream ss(line);
      std::string cell;
      while (std::getline(ss, cell, ',')) {
        try {
          first_row.push_back(std::stod(trim(cell)));
        } catch (...) { break; }
      }
    }
    if (first_row.size() >= 2) {
      // Two columns: pattern, count (maybe more rows)
      pattern_ids.push_back(static_cast<int>(first_row[0]));
      counts.push_back(first_row[1]);
      while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) break;
        std::istringstream ss(line);
        std::string a, b;
        if (std::getline(ss, a, ',') && std::getline(ss, b, ',')) {
          try {
            pattern_ids.push_back(std::stoi(trim(a)));
            counts.push_back(std::stod(trim(b)));
          } catch (...) {}
        }
      }
    } else {
      // Single row of counts (pattern order from config or 1..N)
      for (double v : first_row) counts.push_back(v);
      if (!pattern_roi_hint.empty() && pattern_roi_hint.size() == counts.size()) {
        pattern_ids = pattern_roi_hint;
      } else {
        for (size_t i = 0; i < counts.size(); ++i)
          pattern_ids.push_back(static_cast<int>(i) + 1);
      }
    }
  }

  if (counts.empty()) {
    std::cerr << "[csv2root] No data parsed.\n";
    return false;
  }

  // If we have pattern_roi_hint and parsed (pattern, count) pairs, reorder to pattern_roi
  if (!pattern_roi_hint.empty() && pattern_ids.size() == counts.size()) {
    std::map<int, double> pat_to_count;
    for (size_t i = 0; i < pattern_ids.size(); ++i)
      pat_to_count[pattern_ids[i]] = counts[i];
    pattern_ids = pattern_roi_hint;
    counts.clear();
    for (int pid : pattern_ids) {
      auto it = pat_to_count.find(pid);
      counts.push_back(it != pat_to_count.end() ? it->second : 0.0);
    }
  }

  return true;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " <data.csv> <output.root> [config.json]\n";
    std::cerr << "  config.json optional: used to set pattern_roi order and labels (experiment.pattern_roi).\n";
    return 1;
  }
  const std::string csv_path(argv[1]);
  const std::string root_path(argv[2]);
  const std::string config_path(argc >= 4 ? argv[3] : "");

  std::vector<int> pattern_roi = get_pattern_roi_from_config(config_path);
  if (!config_path.empty() && pattern_roi.empty())
    std::cout << "[csv2root] Warning: no experiment.pattern_roi in config; using CSV order.\n";

  std::vector<int> pattern_ids;
  std::vector<double> counts;
  ExposureMass em;
  if (!parse_csv(csv_path, pattern_ids, counts, pattern_roi, &em)) {
    return 1;
  }

  const int np = static_cast<int>(counts.size());
  TH1D h_data("D_pat",
              "Observed counts per pattern;pattern bin;counts",
              np, 0.5, np + 0.5);
  for (int i = 0; i < np; ++i) {
    h_data.SetBinContent(i + 1, counts[i]);
    h_data.GetXaxis()->SetBinLabel(i + 1, std::to_string(pattern_ids[i]).c_str());
  }

  TFile fout(root_path.c_str(), "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[csv2root] Cannot create " << root_path << "\n";
    return 1;
  }
  h_data.Write();
  if (em.has_value) {
    TParameter<double> p_exp("exposure_kg_year", em.exposure_kg_year);
    TParameter<double> p_mass("mass_kg", em.mass_kg);
    TParameter<double> p_livetime("livetime_days", em.livetime_days);
    p_exp.Write();
    p_mass.Write();
    p_livetime.Write();
  }
  fout.Close();

  std::cout << "[csv2root] Wrote " << np << " pattern bins to " << root_path << " (histogram D_pat, same style as S_pat_validation).\n";
  if (em.has_value) {
    std::cout << "[csv2root] Also wrote livetime_days, mass_kg (effective mass), and exposure_kg_year (from texp and Nusedpix).\n";
  }
  return 0;
}
