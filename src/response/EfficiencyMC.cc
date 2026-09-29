// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  EfficiencyMC.cc -- Monte Carlo builder of P(pattern|n_e) and ε(n_e) by
//  simulating rows or images and scanning them with PatternClassifier.
// ===========================================================================

#include "ccdarksens/response/EfficiencyMC.hh"
#include "ccdarksens/response/PatternImageGenerator.hh"

#include <TH1D.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace ccdarksens {

// Constructor: keep the settings and the shared helpers (the charge-transport model owns its own RNG).
EfficiencyMC::EfficiencyMC(const EfficiencyMCConfig& cfg,
                     std::shared_ptr<ChargeTransport> ct,
                     std::shared_ptr<PatternClassifier> classifier)
  : cfg_(cfg),
    ct_(std::move(ct)),
    classifier_(std::move(classifier))
{
  // ct_ owns its RNG; cfg_.seed is there in case I later want extra RNG here.
}

// Switch the pattern table to the 2D image path by handing over an image generator.
void EfficiencyMC::SetPatternImageGenerator(std::shared_ptr<PatternImageGenerator> gen) {
  img_gen_ = std::move(gen);
}

// -----------------------------------------------------------------------------
// Helper: encode a PatternLabel into a unique integer code for easy printing
// -----------------------------------------------------------------------------
int EncodePatternCode(const std::vector<int>& q) {
  int code = 0;
  for (int d : q) {
    code = code * 10 + d;
  }
  return code;
}

// -----------------------------------------------------------------------------
// Build pattern_table_[n_e][label] = P(label | n_e)
// by simulating cfg_.ne_trials events for each n_e in [ne_min, ne_max].
// -----------------------------------------------------------------------------
void EfficiencyMC::BuildPatternTable(int ne_min, int ne_max, double Ee_eV)
{
  pattern_table_.clear();

  // --- 2D image path (notebook-style: generate_image_E + scan with isolation) ---
  if (img_gen_) {
    for (int ne_true = ne_min; ne_true <= ne_max; ++ne_true) {
      std::map<PatternLabel, std::size_t> counts;  // how many images gave each pattern label
      std::size_t n_trials = 0;  // images actually simulated
      for (int it = 0; it < cfg_.ne_trials; ++it) {
        auto image_2d = img_gen_->GenerateImage(ne_true, Ee_eV);
        const int nrows = img_gen_->NrowsBinned();
        const int ncols = img_gen_->Ncols();
        if (nrows < 3 || ncols < 5) continue;
        const std::vector<double>& row_above  = image_2d[0];
        const std::vector<double>& row_middle = image_2d[1];
        const std::vector<double>& row_below  = image_2d[2];
        // Python scan_image() collects ALL patterns per image; count each.
        auto all_prs = classifier_->ScanAllImage2DWithIsolation(
            row_above, row_middle, row_below, -1.0);
        for (const auto& pr : all_prs)
          counts[pr.label] += 1;
        n_trials++;
      }
      std::map<PatternLabel, double> plabels;
      if (n_trials > 0) {
        for (auto& kv : counts) {
          plabels[kv.first] =
              static_cast<double>(kv.second) / static_cast<double>(n_trials);
        }
      }
      pattern_table_[ne_true] = std::move(plabels);
    }
    return;
  }

  // --- 1D row path (original) ---
  PixelSimulator pixSim(cfg_.pix_cfg);
  const int row_len = cfg_.row_length;

  if (row_len <= 0) {
    throw std::runtime_error("EfficiencyMC::BuildPatternTable: row_length <= 0");
  }

  if (cfg_.pix_cfg.nx < row_len || cfg_.pix_cfg.ny != 1) {
    std::cerr << "[EfficiencyMC] WARNING: pix_cfg {nx,ny} = {"
              << cfg_.pix_cfg.nx << "," << cfg_.pix_cfg.ny
              << "} not fully consistent with row_length=" << row_len
              << " (expected nx>=row_length, ny=1)\n";
  }

  for (int ne_true = ne_min; ne_true <= ne_max; ++ne_true) {

    std::map<PatternLabel, std::size_t> counts;  // how many trials gave each pattern label
    std::size_t n_trials = 0;  // trials that produced a usable row

    for (int it = 0; it < cfg_.ne_trials; ++it) {

      pixSim.Reset();

      std::vector<double> xs(ne_true), ys(ne_true);
      const double z_um = ct_->SampleDepthUm();
      const double sigma_xy_um = ct_->SigmaXYUm(z_um, Ee_eV);
      ct_->SampleCloudXY(0.0, 0.0, sigma_xy_um, ne_true, xs, ys);

      for (int i = 0; i < ne_true; ++i) {
        pixSim.DepositElectron(xs[i], ys[i]);
      }

      if (cfg_.include_dc_pileup) {
        pixSim.AddDarkCurrent();
      }
      pixSim.AddReadoutNoise();

      const std::vector<double>& qpix = pixSim.PixelCharges();

      if ((int)qpix.size() < row_len) {
        continue;
      }

      std::vector<double> row(qpix.begin(), qpix.begin() + row_len);

      auto patterns = classifier_->ScanRow(row);

      const PatternResult* best = nullptr;
      for (const auto& pr : patterns) {
        if (!pr.valid) continue;
        if (!best || pr.total_charge_e > best->total_charge_e) {
          best = &pr;
        }
      }

      if (best) {
        counts[best->label] += 1;
      }

      n_trials++;
    } // end trials

    // Normalize to get P(label | n_e)
    std::map<PatternLabel, double> plabels;
    if (n_trials > 0) {
      for (auto& kv : counts) {
        plabels[kv.first] =
            static_cast<double>(kv.second) / static_cast<double>(n_trials);
      }
    }

    pattern_table_[ne_true] = std::move(plabels);
  }
}

// -----------------------------------------------------------------------------
// Print example row(s) for one n_e (notebook-style: "how the image looks")
// -----------------------------------------------------------------------------
void EfficiencyMC::PrintExampleRows(int ne, double Ee_eV, int n_examples,
                                 std::ostream& out)
{
  if (!ct_ || !classifier_) return;
  PixelSimulator pixSim(cfg_.pix_cfg);
  const int row_len = cfg_.row_length;
  if (row_len <= 0 || ne <= 0) return;

  out << "[EfficiencyMC] Example row(s) for n_e = " << ne
      << ", E_ref = " << Ee_eV << " eV (" << n_examples << " trial(s)), row_length = " << row_len << " px\n";
  if (row_len <= 5) {
    out << "  (short row: pattern (1) is common when electrons pile in one pixel; for multi-pixel patterns 11, 21, etc., set efficiency_mc.row_length in config, e.g. 32)\n";
  }

  for (int ex = 0; ex < n_examples; ++ex) {
    pixSim.Reset();
    std::vector<double> xs(ne), ys(ne);
    const double z_um = ct_->SampleDepthUm();
    const double sigma_xy_um = ct_->SigmaXYUm(z_um, Ee_eV);
    ct_->SampleCloudXY(0.0, 0.0, sigma_xy_um, ne, xs, ys);
    for (int i = 0; i < ne; ++i) {
      pixSim.DepositElectron(xs[i], ys[i]);
    }
    if (cfg_.include_dc_pileup) pixSim.AddDarkCurrent();
    pixSim.AddReadoutNoise();
    const std::vector<double>& qpix = pixSim.PixelCharges();
    if ((int)qpix.size() < row_len) continue;

    std::vector<double> row(qpix.begin(), qpix.begin() + row_len);
    auto patterns = classifier_->ScanRow(row);
    const PatternResult* best = nullptr;
    for (const auto& pr : patterns) {
      if (!pr.valid) continue;
      if (!best || pr.total_charge_e > best->total_charge_e) best = &pr;
    }

    out << "  trial " << (ex + 1) << "  row_charges(e) =";
    for (int i = 0; i < row_len; ++i) {
      out << " " << std::fixed << std::setprecision(3) << row[i];
    }
    out << "\n";
    if (best) {
      out << "           classified pattern = " << best->label
          << "  (charge = " << std::fixed << std::setprecision(3)
          << best->total_charge_e << " e)\n";
    } else {
      out << "           classified pattern = (none)\n";
    }
  }
}

// -----------------------------------------------------------------------------
// Precompute epsilon(n_e) with pattern-efficiency weights from CSV:
//   epsilon(n_e) = sum_{label in accepted_labels} P(label | n_e) * Eff_csv(code, n_e)
// where code = EncodePatternCode(label.q).
//
// Fast path (efficiency_csv provided): skip MC entirely and use CSV values
// directly as eps(ne) = sum_{accepted label L} Eff_csv(L, ne).
// -----------------------------------------------------------------------------
std::unique_ptr<TH1D>
EfficiencyMC::PrecomputeEpsilonWithPatternEff(
    int ne_min, int ne_max, double Ee_eV,
    const std::map<std::pair<int,int>, double>& pattern_eff)
{
  if (!ct_ || !classifier_) {
    throw std::runtime_error(
      "EfficiencyMC::PrecomputeEpsilonWithPatternEff: missing ChargeTransport or PatternClassifier");
  }

  const int nbins = ne_max - ne_min + 1;
  auto h = std::make_unique<TH1D>("eps_mc_csv",
                                  "EfficiencyMC Efficiency with CSV weights",
                                  nbins, ne_min - 0.5, ne_max + 0.5);
  h->Sumw2();

  // CSV-only fast path: efficiency_csv provided → use CSV values directly as eps(ne), skip MC
  if (!pattern_eff.empty()) {
    for (int ne = ne_min; ne <= ne_max; ++ne) {
      double eps_val = 0.0;
      for (const auto& lab : cfg_.accepted_labels) {
        int code = EncodePatternCode(lab.q);
        auto it_w = pattern_eff.find({code, ne});
        if (it_w != pattern_eff.end()) eps_val += it_w->second;
      }
      h->SetBinContent(ne - ne_min + 1, eps_val);
    }
    return h;
  }

  // Full MC path: build pattern table, then weight by CSV
  BuildPatternTable(ne_min, ne_max, Ee_eV);

  // Loop over true n_e
  for (int ne = ne_min; ne <= ne_max; ++ne) {

    auto it = pattern_table_.find(ne);
    if (it == pattern_table_.end()) {
      h->SetBinContent(ne - ne_min + 1, 0.0);
      h->SetBinError  (ne - ne_min + 1, 0.0);
      continue;
    }

    const auto& plabels = it->second; // map<PatternLabel,double> = P(L | ne)

    double eps_val = 0.0;

    // Sum over accepted pattern labels, but now each label is weighted by Eff_csv(code, ne)
    for (const auto& lab : cfg_.accepted_labels) {

      auto ip = plabels.find(lab);
      if (ip == plabels.end()) continue;

      // Pattern label → integer code (1,2,11,21,...)
      int code = EncodePatternCode(lab.q);

      double eff_csv = 1.0;
      auto key  = std::make_pair(code, ne);
      auto it_w = pattern_eff.find(key);
      if (it_w != pattern_eff.end()) {
        eff_csv = it_w->second;
      } else {
        // If there is no entry in the CSV, eff_csv could be 1.0 or 0.0.
        // Start with 1.0 to avoid artificially killing patterns that weren't tabulated.
        eff_csv = 1.0;
      }

      eps_val += ip->second * eff_csv;
    }

    h->SetBinContent(ne - ne_min + 1, eps_val);
    h->SetBinError  (ne - ne_min + 1, 0.0);
  }

  return h;
}


// -----------------------------------------------------------------------------
// Precompute epsilon(n_e) in a TH1D.
//   epsilon(n_e) = sum_{label in cfg_.accepted_labels} P(label | n_e).
// -----------------------------------------------------------------------------
std::unique_ptr<TH1D> EfficiencyMC::PrecomputeEpsilon(int ne_min, int ne_max,
                                                   double Ee_eV)
{
  if (!ct_ || !classifier_) {
    throw std::runtime_error(
      "EfficiencyMC::PrecomputeEpsilon: missing ChargeTransport or PatternClassifier");
  }

  // Build the pattern table via MC
  BuildPatternTable(ne_min, ne_max, Ee_eV);

  // Optional debug printout
  for (int ne = ne_min; ne <= ne_max; ++ne) {
    auto it = pattern_table_.find(ne);
    if (it == pattern_table_.end()) continue;

    std::cerr << "[EfficiencyMC] n_e = " << ne << "\n";
    for (const auto& kv : it->second) {
      std::cerr << "  label=" << kv.first
                << "  P=" << kv.second << "\n";
    }
  }

  // Prepare output histogram
  const int nbins = ne_max - ne_min + 1;
  auto h = std::make_unique<TH1D>("eps_mc", "EfficiencyMC Efficiency",
                                  nbins, ne_min - 0.5, ne_max + 0.5);
  h->Sumw2();

  // Fill bins with ε(n_e)
  for (int ne = ne_min; ne <= ne_max; ++ne) {

    auto it = pattern_table_.find(ne);
    if (it == pattern_table_.end()) {
      h->SetBinContent(ne - ne_min + 1, 0.0);
      h->SetBinError  (ne - ne_min + 1, 0.0);
      continue;
    }

    const auto& plabels = it->second;

    double eps_val = 0.0;
    // Sum over accepted pattern labels
    for (const auto& lab : cfg_.accepted_labels) {
      auto ip = plabels.find(lab);
      if (ip != plabels.end()) {
        eps_val += ip->second;
      }
    }

    h->SetBinContent(ne - ne_min + 1, eps_val);
    h->SetBinError  (ne - ne_min + 1, 0.0);
  }

  return h;
}

// -----------------------------------------------------------------------------
// Energy-dependent ε(E, n_e).
//   - E_grid: list of E_ref in eV.
//   - For each E_ref, we rebuild the pattern_table_ and sum over accepted labels.
// -----------------------------------------------------------------------------
void EfficiencyMC::PrecomputeEpsilonVsEnergy(const std::vector<double>& E_grid,
                                          int ne_min, int ne_max)
{
    if (!ct_ || !classifier_) {
      throw std::runtime_error(
        "EfficiencyMC::PrecomputeEpsilonVsEnergy: missing ChargeTransport or PatternClassifier");
    }

    energy_grid_eV_ = E_grid;
    int NE = static_cast<int>(E_grid.size());
    int Nn = ne_max - ne_min + 1;

    epsilon_Ene_.assign(NE, std::vector<double>(Nn, 0.0));

    for (int iE = 0; iE < NE; ++iE) {
        double Eref = E_grid[iE];

        // Rebuild pattern table at this energy
        BuildPatternTable(ne_min, ne_max, Eref);

        // For each n_e, sum probabilities over accepted labels
        for (int ne = ne_min; ne <= ne_max; ++ne) {
            auto it = pattern_table_.find(ne);
            if (it == pattern_table_.end()) continue;

            const auto& plabels = it->second;
            double eps_val = 0.0;
            for (const auto& acc_lab : cfg_.accepted_labels) {
                auto ip = plabels.find(acc_lab);
                if (ip != plabels.end()) {
                    eps_val += ip->second;
                }
            }
            epsilon_Ene_[iE][ne - ne_min] = eps_val;
        }
    }
}

// ----------------------------------------------------------------------------
// EfficiencyMC::IsAcceptedPattern
//   True if at least one accepted, isolated pattern label has a total charge
//   (the sum of its digits) equal to ne.
// ----------------------------------------------------------------------------
bool EfficiencyMC::IsAcceptedPattern(int ne) const {
    for (const auto& lab : cfg_.accepted_labels) {
        int qsum = 0;
        for (int q : lab.q) {
            qsum += q;
        }
        if (qsum == ne && lab.isolated) {
            return true;
        }
    }
    return false;
}

} // namespace ccdarksens
