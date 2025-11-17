#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/detector/Detector.hh"

#include "ccdarksens/rates/RateSource.hh"
#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/backgrounds/PoissonDarkCurrent.hh"
#include "ccdarksens/backgrounds/BackgroundBuilder.hh"

#include <TH1D.h>
#include <TFile.h>
#include <filesystem>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>

// helper: build a flat epsilon(ne) TH1D aligned to [ne_min, ne_max] integer bins
static std::unique_ptr<TH1D> MakeFlatEfficiency(int ne_min, int ne_max, double epsilon) {
    const int nbin = ne_max - ne_min + 1;
    std::vector<double> edges(nbin + 1);
    for (int i = 0; i <= nbin; ++i) edges[i] = (ne_min - 0.5) + i;

    auto h = std::make_unique<TH1D>("eps_ne", "Pattern efficiency; n_{e}; #epsilon", nbin, edges.data());
    for (int b = 1; b <= nbin; ++b) h->SetBinContent(b, std::clamp(epsilon, 0.0, 1.0));
    return h;
}

// helper: print first N bins of a TH1D with integer centers
static void PrintFirstBins(const TH1D& h, int n_to_print = 10) {
    const int nb = h.GetNbinsX();
    const int n = std::min(nb, n_to_print);
    std::cout << "    n_e : counts\n";
    for (int i = 1; i <= n; ++i) {
        const int ne = static_cast<int>(std::lround(h.GetBinCenter(i)));
        std::cout << "    " << ne << " : " << h.GetBinContent(i) << "\n";
    }
    if (nb > n) std::cout << "    ... (" << (nb - n) << " more bins)\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: ccdarksens_sensitivity <config.json>\n";
        return 1;
    }

    try {
        // -------------------- Load config & experiment --------------------
        ccdarksens::ConfigManager cfg(argv[1]);
        cfg.parse();

        const auto& det = cfg.detector();
        const double mass_kg = det.mass_kg();

        ccdarksens::ExperimentSetup setup(cfg.experiment_cfg(), mass_kg, cfg.run().rng_seed);
        auto summary = setup.prepare_summary();

        const int ne_min = summary.binning.ne_min;
        const int ne_max = summary.binning.ne_max;

        std::cout << "[ccdarksens] Run: " << cfg.run().label << "\n"
                  << "  Mass [kg]: " << mass_kg << "\n"
                  << "  Exposure [kg·day]: " << summary.exposure_kg_year << "\n"
                  << "  Mode: " << summary.mode_string << "\n"
                  << "  Binning n_e: [" << ne_min << ", " << ne_max << "]\n";

        // -------------------- Synthetic dR/dE for a quick signal check (optional) --------------------
        // This is just to verify the response pipeline; sensitivity defaults to background-only.
        auto dRdE = ccdarksens::RateSource::MakeLinearEnergyHist(10.0, 100.0, 180, "dRdE_flat");
        ccdarksens::RateSource::FillFlat(*dRdE, /*integral per kg·day*/ 1.0);

        // -------------------- Response modules (for signal path) --------------------
        auto ion = std::make_shared<ccdarksens::ChargeIonization>("data/p100K_table.csv");
        auto diff = std::make_shared<ccdarksens::Diffusion>(
            /*A_um2*/ 803.25,
            /*b_umInv*/ 6.5e-4,
            /*alpha*/ 1.0,
            /*beta_per_keV*/ 0.0,
            /*thickness_um*/ 670.0,
            /*kappa_e_per_um*/ 0.08,
            /*sigma_readout_e*/ 0.16
        );

        // pattern efficiency for signal (use same ε as background unless specified otherwise)
        std::unique_ptr<TH1D> eps_sig;
        if (cfg.backgrounds().has_flat_eps) {
            eps_sig = MakeFlatEfficiency(ne_min, ne_max, cfg.backgrounds().flat_eps);
        } else {
            eps_sig = MakeFlatEfficiency(ne_min, ne_max, 1.0);
        }
        auto pe = std::make_shared<ccdarksens::PatternEfficiency>();
        pe->SetEfficiencyHist(*eps_sig);

        ccdarksens::DetectorResponsePipeline pipe(ion);
        pipe.SetDiffusion(diff);
        pipe.SetPatternEfficiency(pe);

        // Build signal histograms (for inspection; NOT used in background-only sensitivity)
        auto S_raw = ion->FoldToNe(*dRdE, summary.exposure_kg_year, ne_min, ne_max);
        double Ee_ref_eV = 50.0; // only matters if beta_per_keV > 0
        auto S_obs = pipe.Apply(*dRdE, summary.exposure_kg_year, ne_min, ne_max, Ee_ref_eV);

        // -------------------- Background-only Asimov builder --------------------
        ccdarksens::BackgroundBuilder bld(det.geometry().rows,
                                          det.geometry().cols,
                                          det.geometry().active_fraction,
                                          ne_min, ne_max);

        // Timing from experiment & timing shim in backgrounds
        ccdarksens::TimingConfig tcfg;
        const auto& timing_json = cfg.timing();
        tcfg.exposure_time_s = timing_json.exposure_time_s;
        if (timing_json.n_exposures_override.has_value())
            tcfg.n_exposures_override = timing_json.n_exposures_override;

        bld.SetTiming(cfg.experiment_cfg().livetime_days,
                      cfg.experiment_cfg().duty_cycle,
                      tcfg);

        // Dark current settings
        ccdarksens::DarkCurrentConfig dcc;
        // Dark current is now configured in e-/pixel/year
        dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
        dcc.norm_scale = cfg.backgrounds().norm_scale;
        bld.SetDarkCurrent(dcc);

        // Pattern efficiency on the background (flat ε for MVP)
        if (cfg.backgrounds().has_flat_eps) {
            auto eps_bkg = MakeFlatEfficiency(ne_min, ne_max, cfg.backgrounds().flat_eps);
            auto pe_bkg  = std::make_shared<ccdarksens::PatternEfficiency>();
            pe_bkg->SetEfficiencyHist(*eps_bkg);
            bld.SetPatternEfficiency(pe_bkg);
        }

        // Build Asimov background
        auto B_asimov = bld.BuildBkgAsimov();

        // -------------------- Sanity checks --------------------
        std::cout << "\n[checks]\n";
        std::cout << "  ∫ dE dR/dE (per kg·day): " << dRdE->Integral("width") << "\n";
        std::cout << "  S_raw integral (counts): " << S_raw->Integral() << "\n";
        std::cout << "  S_obs integral (counts): " << S_obs->Integral() << "\n";
        std::cout << "  B_asimov integral (counts): " << B_asimov->Integral() << "\n";

        std::cout << "\n[background-asimov]\n";
        std::cout << "  N_active_pixels : " << bld.NActivePixels() << "\n";
        std::cout << "  N_exposures     : " << bld.NExposuresUsed() << "\n";
        std::cout << "  lambda (e-/pix/exposure): " << bld.LambdaPerExposure() << "\n";

        // -------------------- Write ROOT outputs --------------------
        std::filesystem::create_directories(cfg.run().outdir);
        std::string outpath = cfg.run().outdir + "/mvp_check.root";
        TFile fout(outpath.c_str(), "RECREATE");
        dRdE->Write();
        S_raw->Write("S_ne_raw");
        S_obs->Write("S_ne_obs");
        B_asimov->Write("B_asimov_ne");
        fout.Close();

        std::cout << "\n[output] histograms written to: " << outpath << "\n";
        std::cout << "        (For background-only Asimov sensitivity, use B_asimov_ne as the pseudo-data.)\n";

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}
