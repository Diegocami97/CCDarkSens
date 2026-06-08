// ============================================================================
//  CCDarkSens — ccdarksens_scan_dmelectron_toy
//  Minimal background-free DM-e toy scan in charge space using flat efficiency and a fixed “3 events” discovery criterion over a hardcoded mass grid.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include <TH1D.h>
#include <TFile.h>
#include <TGraph.h>

#include "ccdarksens/model/DMElectronModel.hh"
#include "ccdarksens/response/ChargeIonization.hh"

using namespace ccdarksens;

// ----------------------------------------------------------------------------
// Helper: build a DMElectronConfig for given mchi, sigma
// ----------------------------------------------------------------------------
static DMElectronConfig MakeDMConfig(double mchi_MeV, const std::string& sigma_str)
{
  DMElectronConfig mc;
  mc.material          = "Si";
  mc.mediator          = "heavy";
  mc.rates_dir         = "data/qedark_rates/Si/heavy/after_eta_fix";
  mc.filename_template = "dRdE_{material}_{mediator}_m{mchi_MeV}_s{sigma_e_cm2}.csv";

  // Energy grid: match your pattern config (0–20 eV, 200 bins) :contentReference[oaicite:1]{index=1}
  mc.Emin_eV           = 0.0;
  mc.Emax_eV           = 20.0;
  mc.nbins             = 200;

  mc.mchi_MeV          = mchi_MeV;
  mc.sigma_e_cm2       = sigma_str;

  return mc;
}

// ----------------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // Hardcode some basic "experiment" settings, independent of any JSON.
  const double exposure_kg_year = 1.0;   // exposure in kg·year (can change)
  const double N_lim            = 3.0;   // "3 events" criterion (background-free)
  const double sigma0           = 1.4e-40; // reference cross-section [cm^2]
  const std::string sigma0_str  = "1.4e-40";

  // Electron-ROI: n_e = 1–5 (like typical DM-e analyses)
  const int ne_min = 0;                  // ChargeIonization can start at 0
  const int ne_max = 20;                 // must cover at least up to ROI max
  const std::vector<int> roi_bins = {1, 2, 3, 4, 5};

  // Simple, flat detection efficiency (Mike-like toy)
  const double eps_flat = 0.95;

  // DM mass grid [MeV] — you can adjust as needed
  std::vector<double> mchi_list = {0.53, 1.0, 3.0, 10.0, 100.0, 1000.0};

  std::cout << "=== CCDarkSens DM-e Toy Scan (Charge Space, Bkg-Free) ===\n";
  std::cout << "Exposure       : " << exposure_kg_year << " kg·year\n";
  std::cout << "N_lim          : " << N_lim << " events (background-free)\n";
  std::cout << "sigma0         : " << sigma0 << " cm^2 (reference)\n";
  std::cout << "ROI (n_e)      : { ";
  for (int ne : roi_bins) std::cout << ne << " ";
  std::cout << "}\n";
  std::cout << "Flat epsilon   : " << eps_flat << "\n\n";

  // --------------------------------------------------------------------------
  // Charge Ionization model: P(n_e | E)
  // Uses your existing table path (same as pattern app) :contentReference[oaicite:2]{index=2}
  // --------------------------------------------------------------------------
  auto ion = std::make_shared<ChargeIonization>("data/p100K_table.csv");

  // Output containers for the limit curve
  std::vector<double> mchi_vals;
  std::vector<double> sigma_lim_vals;

  std::cout << std::scientific << std::setprecision(6);

  for (double mchi : mchi_list) {
    // ------------------------------------------------------------------------
    // 1) Configure DM model at (mchi, sigma0), build dR/dE
    // ------------------------------------------------------------------------
    DMElectronConfig mc = MakeDMConfig(mchi, sigma0_str);
    DMElectronModel dm;
    std::cout << "  [toy-scan] Using rates_dir=" << mc.rates_dir
          << ", mchi=" << mc.mchi_MeV
          << ", sigma=" << mc.sigma_e_cm2 << "\n";

    if (!dm.Configure(mc)) {
      std::cerr << "[toy-scan] WARNING: failed to configure DM model for mchi="
                << mchi << " MeV, sigma=" << sigma0_str << " cm^2\n";
      mchi_vals.push_back(mchi);
      sigma_lim_vals.push_back(std::numeric_limits<double>::quiet_NaN());
      continue;
    }

    auto dRdE = dm.MakeSpectrum_E();  // TH1D: dR/dE [events/(kg·year·eV)]

    // Optional debug: integrated rate
    double dRdE_int = dRdE->Integral("width"); // events / (kg·year)
    double N_raw = dRdE_int * exposure_kg_year;

    std::cout << "[toy-scan] mchi=" << mchi << " MeV\n";
    std::cout << "  ∫ dR/dE dE = " << dRdE_int
              << " events/(kg·year), N_raw=" << N_raw << " events\n";

    // ------------------------------------------------------------------------
    // 2) Fold dR/dE → S_true(n_e) via ChargeIonization
    //    S_true(n_e) is the *true* number of DM events in each n_e bin.
    // ------------------------------------------------------------------------
    auto h_ne_true = ion->FoldToNe(*dRdE, exposure_kg_year, ne_min, ne_max);
    // Apply flat, n_e-independent detection efficiency
    h_ne_true->Scale(eps_flat);

    // Debug: print S_true(n_e) in ROI at sigma0
    std::cout << "  S_true(n_e) in ROI [sigma0=" << sigma0 << " cm^2]:\n";
    for (int ne : roi_bins) {
      int bin = h_ne_true->FindBin(static_cast<double>(ne));
      double s = h_ne_true->GetBinContent(bin);
      std::cout << "    n_e=" << ne << " -> S = " << s << " events\n";
    }

    // ------------------------------------------------------------------------
    // 3) Compute total signal in ROI at sigma0
    // ------------------------------------------------------------------------
    double S_ROI_sigma0 = 0.0;
    for (int ne : roi_bins) {
      int bin = h_ne_true->FindBin(static_cast<double>(ne));
      S_ROI_sigma0 += h_ne_true->GetBinContent(bin);
    }

    std::cout << "  S_ROI(sigma0) = " << S_ROI_sigma0 << " events in ROI\n";

    double sigma_lim = std::numeric_limits<double>::quiet_NaN();
    if (S_ROI_sigma0 > 0.0) {
      // Background-free limit: sigma_lim = N_lim / S_ROI(sigma0) * sigma0
      sigma_lim = (N_lim / S_ROI_sigma0) * sigma0;
      std::cout << "  => sigma_lim (bkg-free, N_lim=" << N_lim
                << ") = " << sigma_lim << " cm^2\n\n";
    } else {
      std::cout << "  => S_ROI(sigma0)=0 → sigma_lim undefined (set to NaN)\n\n";
    }

    mchi_vals.push_back(mchi);
    sigma_lim_vals.push_back(sigma_lim);
  }

  // --------------------------------------------------------------------------
  // Save limit curve to ROOT file for plotting
  // --------------------------------------------------------------------------
  TGraph gr(static_cast<int>(mchi_vals.size()),
            mchi_vals.data(), sigma_lim_vals.data());
  gr.SetName("g_sigma_lim_toy");
  gr.SetTitle("Toy DM-e limit (charge space, bkg-free);m_{#chi} [MeV];#sigma_{e} [cm^{2}]");

  TFile fout("scan_dmelectron_toy.root", "RECREATE");
  if (!fout.IsOpen()) {
    std::cerr << "[toy-scan] ERROR: cannot create output file scan_dmelectron_toy.root\n";
    return 1;
  }
  gr.Write();
  fout.Close();

  std::cout << "=== Toy scan finished. Output written to scan_dmelectron_toy.root ===\n";

  return 0;
}
