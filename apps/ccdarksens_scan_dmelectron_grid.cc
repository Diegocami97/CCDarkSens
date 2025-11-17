// ccdarksens_scan_dmelectron_grid.cc
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <sstream>


#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>

#include <nlohmann/json.hpp>

#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/model/DMElectronModel.hh"
#include "ccdarksens/experiment/ExperimentSetup.hh"
#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/DetectorResponsePipeline.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/ClusterMC.hh"
#include "ccdarksens/backgrounds/BackgroundBuilder.hh"

using nlohmann::json;
using namespace ccdarksens;

// ----------------------------------------------------------------------
// Small helpers to reuse QEdark-style grid description
// ----------------------------------------------------------------------

// C++ version of the Python _expand_axis() from qedark_generate_grid.py
static std::vector<double> expand_axis(const json& spec, const std::string& kind)
{
  std::vector<double> vals;

  if (spec.is_object()) {
    if (spec.contains("values")) {
      for (const auto& v : spec.at("values")) {
        vals.push_back(static_cast<double>(v.get<double>()));
      }
    }
    if (spec.contains("linspace")) {
      const auto& p   = spec.at("linspace");
      const double a  = p.at("start").get<double>();
      const double b  = p.at("stop").get<double>();
      const int    n  = p.at("num").get<int>();
      const bool endpoint = p.value("endpoint", true);

      if (n <= 0) {
        throw std::runtime_error("linspace.num must be >0 for axis " + kind);
      }
      if (n == 1) {
        vals.push_back(a);
      } else {
        const int nsteps = endpoint ? (n - 1) : n;
        const double step = (b - a) / static_cast<double>(nsteps);
        for (int i = 0; i < n; ++i) {
          vals.push_back(a + i * step);
        }
      }
    }
    if (spec.contains("logspace")) {
      const auto& p  = spec.at("logspace");
      const double a = p.at("start_exp").get<double>();
      const double b = p.at("stop_exp").get<double>();
      const int    n = p.at("num").get<int>();
      const bool endpoint = p.value("endpoint", true);

      if (n <= 0) {
        throw std::runtime_error("logspace.num must be >0 for axis " + kind);
      }
      if (n == 1) {
        vals.push_back(std::pow(10.0, a));
      } else {
        const int nsteps = endpoint ? (n - 1) : n;
        const double step = (b - a) / static_cast<double>(nsteps);
        for (int i = 0; i < n; ++i) {
          const double e = a + i * step;
          vals.push_back(std::pow(10.0, e));
        }
      }
    }
  } else if (spec.is_array()) {
    for (const auto& v : spec) vals.push_back(v.get<double>());
  } else {
    throw std::runtime_error("Grid spec for '" + kind +
                             "' must be dict or list in JSON.");
  }

  // uniq-preserve (same behavior as _uniq_preserve)
  std::vector<double> out;
  out.reserve(vals.size());
  for (double v : vals) {
    if (std::find(out.begin(), out.end(), v) == out.end()) {
      out.push_back(v);
    }
  }
  return out;
}

// Turn a sigma value into a filename string using a format like ".1e"
static std::string format_sigma(double sigma, const std::string& fmt)
{
  int prec = 6; // default
  if (!fmt.empty() && fmt.front() == '.' && (fmt.back() == 'e' || fmt.back() == 'E')) {
    try {
      prec = std::stoi(fmt.substr(1, fmt.size() - 2));
    } catch (...) {
      prec = 6;
    }
  }
  std::ostringstream ss;
  ss.setf(std::ios::scientific);
  ss << std::setprecision(prec) << sigma;
  return ss.str();
}

// Build edges from a sorted list of centers so each point gets its own bin
static std::vector<double> make_edges_from_centers(const std::vector<double>& c)
{
  const std::size_t N = c.size();
  std::vector<double> edges(N + 1);
  if (N == 0) return edges;

  if (N == 1) {
    const double w = std::abs(c[0]) > 0 ? std::abs(c[0]) * 0.5 : 0.5;
    edges[0] = c[0] - w;
    edges[1] = c[0] + w;
    return edges;
  }

  edges[0] = c[0] - 0.5 * (c[1] - c[0]);
  for (std::size_t i = 1; i < N; ++i) {
    edges[i] = 0.5 * (c[i - 1] + c[i]);
  }
  edges[N] = c[N - 1] + 0.5 * (c[N - 1] - c[N - 2]);
  return edges;
}

// Sum over ROI bins
static double SumROI(TH1D& h,
                     const std::vector<int>& roi_bins,
                     int /*ne_min*/)
{
  double s = 0.0;
  for (int ne : roi_bins) {
    const int bin = h.FindBin(ne);
    s += h.GetBinContent(bin);
  }
  return s;
}

// ----------------------------------------------------------------------
int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " config.json\n";
    return 1;
  }

  const std::string config_path = argv[1];

  try {
    // ------------------------------------------------------------------
    // 0) Parse config (framework-wide knobs via ConfigManager)
    // ------------------------------------------------------------------
    ConfigManager cfg(config_path);
    cfg.parse();

    const auto& run = cfg.run();
    const auto& det = cfg.detector();

    std::cout << "[scan] label=" << run.label
              << " outdir='" << run.outdir << "'\n";

    // Outdir
    if (!run.outdir.empty()) {
      std::filesystem::create_directories(run.outdir);
    }

    // Experiment setup
    const double mass_kg = det.mass_kg();
    ExperimentSetup setup(cfg.experiment_cfg(), mass_kg, run.rng_seed);
    auto summary = setup.prepare_summary();

    const int ne_min = summary.binning.ne_min;
    const int ne_max = summary.binning.ne_max;

    std::cout << "[scan] Mass [kg]: " << mass_kg << "\n"
              << "       Exposure [kg·year]: " << summary.exposure_kg_year << "\n"
              << "       n_e binning: [" << ne_min << "," << ne_max << "]\n";

    // ------------------------------------------------------------------
    // 1) Build detector response pipeline (same as check_dmelectron_chain_cluster_mc)
    // ------------------------------------------------------------------
    const auto& mj = cfg.model();

    // DMElectronModel config will be filled per grid point, but we keep
    // everything that does NOT change here.
    DMElectronConfig mc_base;
    mc_base.material          = mj.material;
    mc_base.mediator          = mj.mediator;
    mc_base.rates_dir         = mj.rates_dir;
    mc_base.filename_template = mj.filename_template;
    mc_base.Emin_eV           = mj.Emin_eV;
    mc_base.Emax_eV           = mj.Emax_eV;
    mc_base.nbins             = mj.nbins;

    // ChargeIonization
    auto ion = std::make_shared<ChargeIonization>("data/p100K_table.csv");

    // Diffusion (same MVP numbers as in your check app)
    auto diff = std::make_shared<Diffusion>(
      803.25,   // A_um2
      6.5e-4,   // b_umInv
      1.0,      // alpha
      0.0,      // beta_per_keV
      det.geometry().thickness_mm * 1000.0, // thickness_um
      0.08,     // sigma_readoutR_e (MVP)
      0.16      // sigma_readoutC_e (MVP)
    );

    // Geometry & ClusterMC
    GeometryLite geo;
    geo.rows         = det.geometry().rows;
    geo.cols         = det.geometry().cols;
    geo.pixel_size_um= det.geometry().pixel_size_um;
    geo.thickness_um = det.geometry().thickness_mm * 1000.0;

    ClusterMCConfig cmc;
    cmc.n_events_per_ne = cfg.response().cmc.n_events_per_ne;
    cmc.rng_seed        = cfg.response().cmc.rng_seed;
    cmc.sigma_readout_e = cfg.response().cmc.sigma_readout_e;
    cmc.Qmin_e          = cfg.response().cmc.Qmin_e;
    cmc.Qmax_e          = cfg.response().cmc.Qmax_e;
    cmc.diffusion = {
      cfg.response().cmc.A_um2,
      cfg.response().cmc.b_umInv,
      cfg.response().cmc.alpha,
      cfg.response().cmc.beta_per_keV
    };
    cmc.rows_bin        = cfg.response().cmc.rows_bin;
    cmc.cols_bin        = cfg.response().cmc.cols_bin;
    cmc.pileup_with_dc  = cfg.response().cmc.pileup_with_dc;

    // Convert λ_year -> λ_exp (for pileup inside ClusterMC)
    {
      const double lambda_year = cfg.backgrounds().lambda_e_per_pix_per_year;
      const double year_s      = 365.25 * 86400.0;
      const double exp_s       = cfg.timing().exposure_time_s;
      double lambda_exp = 0.0;
      if (exp_s > 0.0) lambda_exp = lambda_year * (exp_s / year_s);
      cmc.lambda_dc_per_pix_per_exposure =
        cmc.pileup_with_dc ? lambda_exp : 0.0;
    }

    cmc.half_window_pix = 2;
    cmc.enable_MN  = true;
    cmc.enable_MNL = true;

    ClusterMC clusterMC(cmc, geo);
    auto eps_mc = clusterMC.PrecomputeEpsilon(ne_min, ne_max);

    auto pe = std::make_shared<PatternEfficiency>();
    pe->SetEfficiencyHist(*eps_mc);

    DetectorResponsePipeline pipe(ion);
    pipe.SetDiffusion(diff);
    pipe.SetPatternEfficiency(pe);

    // ------------------------------------------------------------------
    // 2) Asimov background: dark current + flat energy bkg
    // ------------------------------------------------------------------

    // 2a) Flat energy background using same dR/dE binning as QE-Dark spectrum
    //     We don't have dR/dE yet, so temporarily build a dummy, then clone
    //     its binning later when we configure DMElectronModel for the FIRST grid point.
    //
    // For cleanliness, we'll delay building the flat spectrum until we have
    // a real dR/dE. Here we just prepare the Asimov dark current first.

    BackgroundBuilder bld(det.geometry().rows,
                          det.geometry().cols,
                          det.geometry().active_fraction,
                          ne_min, ne_max);

    TimingConfig tcfg;
    tcfg.exposure_time_s = cfg.timing().exposure_time_s;
    if (cfg.timing().n_exposures_override.has_value())
      tcfg.n_exposures_override = cfg.timing().n_exposures_override;

    bld.SetTiming(cfg.experiment_cfg().livetime_days,
                  cfg.experiment_cfg().duty_cycle,
                  tcfg);

    DarkCurrentConfig dcc;
    dcc.lambda_e_per_pix_per_year = cfg.backgrounds().lambda_e_per_pix_per_year;
    dcc.norm_scale                = cfg.backgrounds().norm_scale;
    bld.SetDarkCurrent(dcc);

    auto pe_bkg = std::make_shared<PatternEfficiency>();
    pe_bkg->SetEfficiencyHist(*eps_mc);
    bld.SetPatternEfficiency(pe_bkg);

    auto B_dc_asimov = bld.BuildBkgAsimov(); // n_e spectrum (counts) after ε

    // We'll add the flat component once we know the QE-Dark binning (below).

    // ------------------------------------------------------------------
    // 3) Read grid description directly from JSON (QEdark-style)
    // ------------------------------------------------------------------
    std::ifstream jfin(config_path);
    if (!jfin) {
      throw std::runtime_error("Cannot reopen config for grid parsing: " + config_path);
    }
    json jroot;
    jfin >> jroot;

    if (!jroot.contains("model") || !jroot["model"].contains("grid")) {
      throw std::runtime_error("Config must contain model.grid for the scan.");
    }

    const auto& jgrid = jroot["model"]["grid"];

    auto mchi_list = expand_axis(jgrid.at("mchi_MeV"), "mchi_MeV");
    auto sigma_list = expand_axis(jgrid.at("sigma_e_cm2"), "sigma_e_cm2");

    if (mchi_list.empty() || sigma_list.empty()) {
      throw std::runtime_error("Empty mchi_MeV or sigma_e_cm2 grid.");
    }

    // Optional format block: { "mchi": ".6f", "sigma": ".1e" }
    std::string fmt_sigma = ".1e";
    if (jgrid.contains("format") && jgrid["format"].is_object()) {
      const auto& jf = jgrid["format"];
      fmt_sigma = jf.value("sigma", std::string(".1e"));
      // (mchi format is handled by DMElectronModel itself; we assume .6f.)
    }

    std::cout << "[scan] Grid: "
              << "N_mchi=" << mchi_list.size()
              << "  N_sigma=" << sigma_list.size() << "\n";

    // Prepare 2D histogram for q(mchi, sigma)
    std::vector<double> edges_m = make_edges_from_centers(mchi_list);
    std::vector<double> edges_s = make_edges_from_centers(sigma_list);

    TH2D h_q("q_mchi_sigma",
             ";m_{#chi} [MeV];#sigma_{e} [cm^{2}];q = -2 ln(L/L_{0})",
             static_cast<int>(mchi_list.size()), edges_m.data(),
             static_cast<int>(sigma_list.size()), edges_s.data());

    // ------------------------------------------------------------------
    // 4) Need one reference dR/dE to define flat background binning
    //    Use the FIRST grid point to get the QE-Dark energy binning.
    // ------------------------------------------------------------------
    {
      DMElectronConfig mc_first = mc_base;
      mc_first.mchi_MeV    = mchi_list.front();
      mc_first.sigma_e_cm2 = format_sigma(sigma_list.front(), fmt_sigma);

      DMElectronModel dm_first;
      if (!dm_first.Configure(mc_first)) {
        throw std::runtime_error("DMElectronModel::Configure failed for first grid point.");
      }
      auto dRdE_ref = dm_first.MakeSpectrum_E();
      dRdE_ref->SetName("dRdE_reference");

      // Flat background level in events/(kg·year·keV)
      const double flat_rate_per_kg_year_keV = cfg.backgrounds().has_flat_bkg
        ? cfg.backgrounds().flat_bkg_norm_per_kg_year
        : 0.0;

      const double flat_rate_per_kg_year_eV = flat_rate_per_kg_year_keV / 1000.0;

      auto dR_flat =
        std::unique_ptr<TH1D>(static_cast<TH1D*>(dRdE_ref->Clone("dRdE_flat")));
      dR_flat->Reset("ICES");

      const int nbins_E = dR_flat->GetNbinsX();
      for (int ib = 1; ib <= nbins_E; ++ib) {
        dR_flat->SetBinContent(ib, flat_rate_per_kg_year_eV);
      }

      auto Bflat_obs = pipe.Apply(*dR_flat, summary.exposure_kg_year,
                                  ne_min, ne_max,
                                  /*Emax_keV_for_MC*/ 50.0);

      // Total Asimov background = dark current + flat component
      auto B_asimov = std::unique_ptr<TH1D>(
        static_cast<TH1D*>(B_dc_asimov->Clone("B_asimov_total")));
      B_asimov->Add(Bflat_obs.get());

      // Keep this histogram for the whole scan (same for all grid points)
      B_dc_asimov.swap(B_asimov);
    }

    TH1D& B_asimov = *B_dc_asimov; // alias

    // Make vector of background counts in ROI
    std::vector<int> roi_bins = summary.roi_bins;
    std::cout << "[scan] ROI bins: ";
    for (std::size_t i = 0; i < roi_bins.size(); ++i) {
      if (i) std::cout << ",";
      std::cout << roi_bins[i];
    }
    std::cout << "\n";

    // ------------------------------------------------------------------
    // 5) Scan over grid, compute Asimov q(mchi, sigma)
    // ------------------------------------------------------------------
    int idx_mchi = 0;
    for (double mchi : mchi_list) {
      int idx_sigma = 0;
      for (double sigma_val : sigma_list) {
        // Configure DMElectronModel for this grid point
        DMElectronConfig mc = mc_base;
        mc.mchi_MeV    = mchi;
        mc.sigma_e_cm2 = format_sigma(sigma_val, fmt_sigma);

        DMElectronModel dm_model;
        if (!dm_model.Configure(mc)) {
          std::cerr << "[scan] WARNING: failed to load QE-Dark table for "
                    << "mchi=" << mchi
                    << " sigma=" << sigma_val
                    << " (file name uses sigma_str=" << mc.sigma_e_cm2 << ")\n";
          // Fill with a large q? Or just skip (NaN)? For now set q=0 and continue.
          h_q.SetBinContent(1 + idx_mchi, 1 + idx_sigma, 0.0);
          ++idx_sigma;
          continue;
        }

        auto dRdE = dm_model.MakeSpectrum_E();
        dRdE->SetName("dRdE_signal");

        // Signal in n_e after full response, for this exposure
        auto S_obs = pipe.Apply(*dRdE, summary.exposure_kg_year,
                                ne_min, ne_max,
                                /*Emax_keV_for_MC*/ 50.0);

        // Compute Asimov PLR q = 2 sum_i [ s_i - b_i log(1 + s_i/b_i) ]
        double q = 0.0;
        for (int ne : roi_bins) {
          const int bin = B_asimov.FindBin(ne);
          const double b = B_asimov.GetBinContent(bin);
          const double s = S_obs->GetBinContent(bin);

          if (b <= 0.0) {
            // If background is zero, the Asimov formula tends to 2*s
            // (since ln(1 + s/b) ~ ln(s/b) -> large). MVP: use 2*s.
            q += 2.0 * s;
          } else if (s > 0.0) {
            q += 2.0 * (s - b * std::log(1.0 + s / b));
          }
          // if s=0, contribution is 0.
        }

        h_q.SetBinContent(1 + idx_mchi, 1 + idx_sigma, q);

        if (run.verbosity >= 2) {
          std::cout << "  [scan] mchi=" << mchi
                    << " MeV, sigma=" << sigma_val
                    << " cm^2  -> q=" << q << "\n";
        }

        ++idx_sigma;
      }
      ++idx_mchi;
    }

    // ------------------------------------------------------------------
    // 6) Save output
    // ------------------------------------------------------------------
    const std::string out_path =
      run.outdir.empty()
      ? "scan_dmelectron_grid.root"
      : (run.outdir + "/scan_dmelectron_grid.root");

    TFile fout(out_path.c_str(), "RECREATE");
    if (!fout.IsOpen()) {
      std::cerr << "[scan] ERROR: cannot create output file " << out_path << "\n";
      return 1;
    }

    // Save background and a few summaries
    B_asimov.Write("B_asimov");
    h_q.Write("q_mchi_sigma");

    fout.Write();
    fout.Close();

    std::cout << "[scan] Wrote " << out_path << "\n";
    std::cout << "[scan] Done.\n";

  } catch (const std::exception& ex) {
    std::cerr << "[scan] ERROR: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
