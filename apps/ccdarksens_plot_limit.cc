// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ccdarksens_plot_limit.cc -- ROOT plotting executable that draws DM-
//  electron 90% CL upper-limit curves from scan ROOT outputs and can overlay
//  Brazilian-band envelopes from ccdarksens_band.
//
//  This is my fork of ccdarksens_plot_dmelectron_limit.cc. It adds the --wimp
//  mode (WIMP-nucleon spin-independent literature overlay); the original file
//  stays as the reference for the DM-electron plots (see
//  docs/GenericScanApp_Design.md).
// ===========================================================================

#include <cmath>
#include <cctype>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TLatex.h>
#include <TH1D.h>
#include <TH1.h>
#include <TH2D.h>
#include <TList.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TGraphSmooth.h>
#include <TTree.h>
#include <TError.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TApplication.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TSystem.h>
#include <fstream>
#include <sstream>
#include <algorithm>  // if not already included (for min_element/max_element)

#include <memory>
#include <set>
#include <map>
#include <limits>
#include <cmath>



#include <vector>
#include <limits>
#include <cmath>

// ----------------------------------------------------------------------------
// EnsureOutputDirForFile
//   Create the directory part of file_path (recursively) if the path contains one.
// ----------------------------------------------------------------------------
static void EnsureOutputDirForFile(const std::string& file_path) {
  const std::size_t slash = file_path.find_last_of('/');
  if (slash == std::string::npos || slash == 0) return;
  gSystem->mkdir(file_path.substr(0, slash).c_str(), true);
}

// ----------------------------------------------------------------------------
// ResolveOutputPath
//   Absolute version of a relative output path (prefixed with the current directory), used for log messages.
// ----------------------------------------------------------------------------
static std::string ResolveOutputPath(const std::string& rel) {
  if (rel.empty() || rel[0] == '/') return rel;
  return std::string(gSystem->pwd()) + "/" + rel;
}

// ----------------------------------------------------------------------------
// SanitizeName
//   Turn a label into a legal ROOT object name: every character that is not alphanumeric becomes '_', repeated underscores collapse and trailing ones are removed.
// ----------------------------------------------------------------------------
static std::string SanitizeName(const std::string& s) {
  std::string out;
  for (char c : s) out += (std::isalnum(static_cast<unsigned char>(c)) || c == '_') ? c : '_';
  // Collapse runs of underscores
  std::string result;
  bool prev_under = false;
  for (char c : out) {
    if (c == '_' && prev_under) continue;
    result += c;
    prev_under = (c == '_');
  }
  // Strip trailing underscore
  while (!result.empty() && result.back() == '_') result.pop_back();
  return result;
}

// ----------------------------------------------------------------------------
// SaveMainLimitOutputs
//   Save the finished plot: the canvas as a PDF, and a ROOT file holding every limit curve
//   (renamed from its label), the q histogram, the literature curves and the store-only curves.
//   Returns true if the PDF exists afterwards.
// ----------------------------------------------------------------------------
static bool SaveMainLimitOutputs(TCanvas* c,
                                 const std::vector<TGraph*>& glimits,
                                 TH2D* hq,
                                 const std::string& out_pdf,
                                 const std::string& out_root,
                                 const char* phase,
                                 const std::vector<std::pair<TGraph*, std::string>>& literature = {},
                                 const std::vector<std::pair<TGraph*, std::string>>& store_only = {},
                                 const std::vector<std::string>& limit_labels = {}) {
  if (!c) {
    std::cerr << "[limit] (" << phase << ") Canvas unavailable; skip save.\n";
    return false;
  }
  EnsureOutputDirForFile(out_pdf);
  EnsureOutputDirForFile(out_root);
  c->cd();
  c->Update();
  c->SaveAs(out_pdf.c_str());
  {
    TFile fout(out_root.c_str(), "RECREATE");
    for (std::size_t ig = 0; ig < glimits.size(); ++ig) {
      TGraph* g = glimits[ig];
      if (!g) continue;
      if (ig < limit_labels.size() && !limit_labels[ig].empty()) {
        const std::string safe = SanitizeName(limit_labels[ig]);
        if (!safe.empty()) g->SetName(safe.c_str());
      }
      g->Write();
    }
    if (hq) hq->Write();
    for (const auto& [g, name] : literature) {
      if (g) { g->SetName(name.c_str()); g->Write(); }
    }
    for (const auto& [g, name] : store_only) {
      if (g) { g->SetName(name.c_str()); g->Write(); }
    }
    fout.Close();
  }
  const bool pdf_ok =
      !gSystem->AccessPathName(out_pdf.c_str(), kReadPermission);
  std::cout << "[limit] (" << phase << ") "
            << (pdf_ok ? "Saved" : "WARNING: PDF may be missing") << ": "
            << ResolveOutputPath(out_pdf) << "\n"
            << "[limit]           "
            << ResolveOutputPath(out_root) << "\n";
  return pdf_ok;
}

// ----------------------------------------------------------------------------
// CsvEscape
//   Quote a CSV field (doubling any quotes) if it contains a comma, quote or line break.
// ----------------------------------------------------------------------------
static std::string CsvEscape(const std::string& s) {
  bool need_quotes = false;
  for (char c : s) {
    if (c == '"' || c == ',' || c == '\n' || c == '\r') {
      need_quotes = true;
      break;
    }
  }
  if (!need_quotes) return s;
  std::string out = "\"";
  for (char c : s) {
    if (c == '"') out += "\"\"";
    else out += c;
  }
  out += '"';
  return out;
}

// ----------------------------------------------------------------------------
// SaveLimitContoursCsv
//   Write the limit curves as CSV (mchi_MeV, sigma_e_cm2; a label column when there is more than
//   one curve), skipping non-finite or non-positive points. Returns false if nothing was written.
// ----------------------------------------------------------------------------
static bool SaveLimitContoursCsv(const std::string& path,
                                 const std::vector<TGraph*>& glimits,
                                 const std::vector<std::string>& labels,
                                 double q_threshold) {
  if (path.empty()) return false;
  if (glimits.empty()) {
    std::cerr << "[limit] ERROR: no limit curves to write to CSV.\n";
    return false;
  }
  EnsureOutputDirForFile(path);
  std::ofstream out(path);
  if (!out) {
    std::cerr << "[limit] ERROR: cannot open CSV for write: " << path << "\n";
    return false;
  }
  out << "# Upper-limit contour: mchi_MeV vs sigma_e_cm2 (dark_photon: mA_eV vs epsilon)\n";
  out << "# q_threshold=" << q_threshold << "\n";
  const bool multi = (glimits.size() > 1);
  if (multi) {
    out << "label,mchi_MeV,sigma_e_cm2\n";
  } else {
    out << "mchi_MeV,sigma_e_cm2\n";
  }
  std::size_t n_written = 0;
  for (std::size_t ig = 0; ig < glimits.size(); ++ig) {
    TGraph* g = glimits[ig];
    if (!g) continue;
    const std::string lbl =
        (ig < labels.size() && !labels[ig].empty())
            ? labels[ig]
            : ("curve_" + std::to_string(ig));
    const int n = g->GetN();
    for (int i = 0; i < n; ++i) {
      double mchi = 0.0, sigma = 0.0;
      g->GetPoint(i, mchi, sigma);
      if (!(sigma > 0.0 && std::isfinite(mchi) && std::isfinite(sigma))) continue;
      if (multi) {
        out << CsvEscape(lbl) << ',';
      }
      out << std::scientific << std::setprecision(8) << mchi << ','
          << std::scientific << std::setprecision(8) << sigma << '\n';
      ++n_written;
    }
  }
  if (n_written == 0) {
    std::cerr << "[limit] ERROR: CSV export produced no points.\n";
    return false;
  }
  std::cout << "[limit] Saved contour CSV (" << n_written << " points): "
            << ResolveOutputPath(path) << "\n";
  return true;
}

// ----------------------------------------------------------------------------
// MakeFilledBetween
//   Closed polygon that fills the region between two curves over [x_min, x_max]: n_samples points
//   log-spaced in x along each curve, with the given fill colour and transparency.
// ----------------------------------------------------------------------------
TGraph* MakeFilledBetween(const TGraph* g_bottom,
                          const TGraph* g_top,
                          double x_min,
                          double x_max,
                          int n_samples,
                          int fill_color,
                          double alpha = 0.25)
{
    if (!g_bottom || !g_top || n_samples < 2) return nullptr;

    TGraph* g_band = new TGraph(2 * n_samples);
    g_band->SetName("g_filled_between");
    g_band->SetTitle("");

    // sample in log-x, since my axes are log
    double logxmin = std::log10(x_min);
    double logxmax = std::log10(x_max);

    // upper edge (g_top) from left -> right
    for (int i = 0; i < n_samples; ++i) {
        double logx = logxmin + (logxmax - logxmin) * i / (n_samples - 1);
        double x    = std::pow(10.0, logx);

        double y_top    = g_top->Eval(x, 0, "S");
        double y_bottom = g_bottom->Eval(x, 0, "S");

        // safety: if something goes crazy, skip that point
        if (y_top <= 0 || y_bottom <= 0) {
            y_top    = 0;
            y_bottom = 0;
        }

        g_band->SetPoint(i, x, y_top);
    }

    // lower edge (g_bottom) from right -> left
    for (int i = 0; i < n_samples; ++i) {
        int j = n_samples - 1 - i;
        double logx = logxmin + (logxmax - logxmin) * j / (n_samples - 1);
        double x    = std::pow(10.0, logx);

        double y_bottom = g_bottom->Eval(x, 0, "S");
        if (y_bottom <= 0) y_bottom = 0;

        g_band->SetPoint(n_samples + i, x, y_bottom);
    }

    g_band->SetFillStyle(1001);
    g_band->SetFillColorAlpha(fill_color, alpha);
    g_band->SetLineWidth(0);

    return g_band;
}

// Closed ribbons for Brazilian bands: use only native TGraph knots where both
// edges are finite and strictly positive. This avoids TGraph::Eval extrapolation
// beyond the scan mass range and linear interpolation through y=0 placeholders
// (empty toy bins in ccdarksens_band), which produced vertical spikes and false
// high-mass tails on sparse per-mass grids.
static std::vector<TGraph*> MakeBandRibbonPieces(const TGraph* g_bottom,
                                                 const TGraph* g_top)
{
    std::vector<TGraph*> out;
    if (!g_bottom || !g_top) return out;
    const int n = g_bottom->GetN();
    if (n < 2 || g_top->GetN() != n) return out;

    auto point_ok = [&](int i) -> bool {
        double xb, yb, xt, yt;
        g_bottom->GetPoint(i, xb, yb);
        g_top->GetPoint(i, xt, yt);
        if (!(std::isfinite(xb) && std::isfinite(yb) && std::isfinite(xt) &&
              std::isfinite(yt)))
            return false;
        const double xtol =
            1e-9 * (std::fabs(xb) + std::fabs(xt) + 1.0);
        if (std::fabs(xb - xt) > xtol) return false;
        return yb > 0.0 && yt > 0.0;
    };

    int run_lo = -1;
    auto close_run = [&](int run_hi) {
        if (run_lo < 0) return;
        const int lo = run_lo;
        const int hi = run_hi;
        run_lo = -1;
        if (hi < lo + 1) return;

        const int npt = hi - lo + 1;
        TGraph* g = new TGraph(2 * npt);
        g->SetName("g_band_ribbon");
        g->SetTitle("");
        for (int k = 0; k < npt; ++k) {
            double x, yb, yt;
            g_bottom->GetPoint(lo + k, x, yb);
            g_top->GetPoint(lo + k, x, yt);
            g->SetPoint(k, x, yt);
        }
        for (int k = 0; k < npt; ++k) {
            const int i = hi - k;
            double x, yb, yt;
            g_bottom->GetPoint(i, x, yb);
            g_top->GetPoint(i, x, yt);
            (void)yt;
            g->SetPoint(npt + k, x, yb);
        }
        out.push_back(g);
    };

    for (int i = 0; i < n; ++i) {
        if (!point_ok(i)) {
            if (run_lo >= 0) close_run(i - 1);
            continue;
        }
        if (run_lo < 0) run_lo = i;
    }
    if (run_lo >= 0) close_run(n - 1);
    return out;
}

// Log-log plots: keep only finite knots with x>0, y>0 so Draw("L") does not
// traverse y=0 placeholders from ccdarksens_band (which breaks or vanishes on
// log y) and so the polyline matches the ribbon support.
static TGraph* FilterGraphPositiveLogSafe(const TGraph* g, const char* name)
{
    if (!g) return nullptr;
    std::vector<double> xv, yv;
    for (int i = 0; i < g->GetN(); ++i) {
        double x = 0.0, y = 0.0;
        g->GetPoint(i, x, y);
        if (std::isfinite(x) && std::isfinite(y) && x > 0.0 && y > 0.0) {
            xv.push_back(x);
            yv.push_back(y);
        }
    }
    if (xv.size() < 2) return nullptr;
    TGraph* out = new TGraph(static_cast<int>(xv.size()));
    for (std::size_t i = 0; i < xv.size(); ++i)
        out->SetPoint(static_cast<int>(i), xv[i], yv[i]);
    if (name && name[0] != '\0') out->SetName(name);
    return out;
}

// ----------------------------------------------------------------------------
// median_sorted_copy
//   Median of v (NaN if empty); v is passed by value and sorted internally.
// ----------------------------------------------------------------------------
static double median_sorted_copy(std::vector<double> v)
{
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    if (n % 2u) return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Read band.root TTree band_per_toy (requires save_per_toy_curves=true in band
// config) and write one panel per mass: empirical #sigma_{UL} from outer toys.
// which: "toy" | "asy" | "both" (asy = sigma_UL_asy per toy).
// scan_root_for_obs_ul: if non-empty, load upper_limit_sigma_e_mchi_graph and
// draw one vertical line per panel at the observed UL (#sigma) for that m_chi.
static int draw_band_per_toy_sigma_hists(const std::string& band_root_path,
                                         const std::string& which_lc,
                                         const std::string& scan_root_for_obs_ul)
{
    TFile f(band_root_path.c_str(), "READ");
    if (!f.IsOpen() || f.IsZombie()) {
        std::cerr << "[limit][band-per-toy] cannot open " << band_root_path << "\n";
        return 1;
    }
    TTree* tr = dynamic_cast<TTree*>(f.Get("band_per_toy"));
    if (!tr) {
        std::cerr << "[limit][band-per-toy] TTree 'band_per_toy' missing in "
                  << band_root_path
                  << " — set band.save_per_toy_curves=true and re-run ccdarksens_band.\n";
        return 1;
    }
    long long br_toy_idx = 0;
    double br_mchi = 0.0;
    double br_asy = 0.0;
    double br_toy = 0.0;
    tr->SetBranchAddress("toy_idx", &br_toy_idx);
    tr->SetBranchAddress("mchi_MeV", &br_mchi);
    tr->SetBranchAddress("sigma_UL_asy", &br_asy);
    tr->SetBranchAddress("sigma_UL_toy", &br_toy);

    std::map<double, std::vector<double>> toy_by_m;
    std::map<double, std::vector<double>> asy_by_m;
    const Long64_t nent = tr->GetEntries();
    for (Long64_t i = 0; i < nent; ++i) {
        tr->GetEntry(i);
        if (br_toy > 0.0 && std::isfinite(br_toy)) toy_by_m[br_mchi].push_back(br_toy);
        if (br_asy > 0.0 && std::isfinite(br_asy)) asy_by_m[br_mchi].push_back(br_asy);
    }
    f.Close();

    const bool draw_toy =
        (which_lc == "toy" || which_lc == "both");
    const bool draw_asy =
        (which_lc == "asy" || which_lc == "asymptotic" || which_lc == "both");
    if (!draw_toy && !draw_asy) {
        std::cerr << "[limit][band-per-toy] unknown --band-per-toy-which '"
                  << which_lc << "' (use toy, asy, or both).\n";
        return 1;
    }

    std::vector<double> masses;
    if (draw_toy && draw_asy) {
        std::set<double> ms;
        for (const auto& kv : toy_by_m) ms.insert(kv.first);
        for (const auto& kv : asy_by_m) ms.insert(kv.first);
        for (double m : ms) masses.push_back(m);
    } else if (draw_toy) {
        for (const auto& kv : toy_by_m) masses.push_back(kv.first);
    } else {
        for (const auto& kv : asy_by_m) masses.push_back(kv.first);
    }

    if (masses.empty()) {
        std::cerr << "[limit][band-per-toy] no positive #sigma_{UL} entries found.\n";
        return 1;
    }

    std::unique_ptr<TGraph> g_obs_ul;
    if (!scan_root_for_obs_ul.empty()) {
        TFile fs(scan_root_for_obs_ul.c_str(), "READ");
        if (!fs.IsOpen() || fs.IsZombie()) {
            std::cerr << "[limit][band-per-toy] cannot open scan file for observed UL: "
                      << scan_root_for_obs_ul << "\n";
        } else {
            auto* g_in =
                dynamic_cast<TGraph*>(fs.Get("upper_limit_sigma_e_mchi_graph"));
            if (!g_in) {
                std::cerr << "[limit][band-per-toy] missing "
                             "upper_limit_sigma_e_mchi_graph in "
                          << scan_root_for_obs_ul << " — no observed UL lines.\n";
            } else {
                TObject* cl = g_in->Clone("g_obs_ul_band_per_toy");
                g_obs_ul.reset(dynamic_cast<TGraph*>(cl));
                if (!g_obs_ul) delete cl;
            }
        }
    }
    if (g_obs_ul) {
        const int gn = g_obs_ul->GetN();
        double gx0 = 0.0, gy0 = 0.0, gxm = 0.0, gym = 0.0;
        if (gn > 0) {
            g_obs_ul->GetPoint(0, gx0, gy0);
            g_obs_ul->GetPoint(gn - 1, gxm, gym);
        }
        std::cout << "[limit][band-per-toy] Loaded Obs. #sigma_{UL} (data) graph: N="
                  << gn << " points, m_{#chi} in [" << gx0 << ", " << gxm
                  << "] MeV from\n  " << scan_root_for_obs_ul << "\n";
    } else if (!scan_root_for_obs_ul.empty()) {
        std::cout << "[limit][band-per-toy] No upper_limit_sigma_e_mchi_graph (obs UL lines disabled).\n";
    }

    const int np = static_cast<int>(masses.size());
    const int ncol = std::max(1, static_cast<int>(std::ceil(std::sqrt(static_cast<double>(np)))));
    const int nrow = (np + ncol - 1) / ncol;

    TCanvas* cpt =
        new TCanvas("c_band_per_toy", "Per-toy #sigma_{UL} (band outer toys)",
                    520 * ncol, 380 * nrow);
    cpt->Divide(ncol, nrow);

    // TGraph vertical segment: TLine/DrawClone often fails on log-x subpads of a
    // divided canvas; heap TGraph + "L SAME" stays on the pad until SaveAs.
    auto draw_world_vline = [](double xv, double y0, double y1, Color_t col,
                               Style_t sty, int wid) {
        if (!gPad || !std::isfinite(xv) || !std::isfinite(y0) || !std::isfinite(y1))
            return;
        auto* g = new TGraph(2);
        g->SetPoint(0, xv, y0);
        g->SetPoint(1, xv, y1);
        g->SetLineColor(col);
        g->SetLineStyle(sty);
        g->SetLineWidth(wid);
        g->Draw("L SAME");
    };

    // Use per-pad x ranges so low/high masses remain readable and do not collapse
    // against global scan edges in a single shared axis.

    std::vector<TH1D*> written_for_root;

    for (int ip = 0; ip < np; ++ip) {
        const double mchi = masses[static_cast<std::size_t>(ip)];
        cpt->cd(ip + 1);
        gPad->SetLogx(1);
        gPad->SetLeftMargin(0.11);
        gPad->SetBottomMargin(0.12);

        std::vector<double> vt = toy_by_m[mchi];
        std::vector<double> va = asy_by_m[mchi];

        double smin = std::numeric_limits<double>::infinity();
        double smax = 0.0;
        auto widen = [&](const std::vector<double>& v) {
            for (double s : v) {
                if (!(s > 0.0) || !std::isfinite(s)) continue;
                smin = std::min(smin, s);
                smax = std::max(smax, s);
            }
        };
        double sig_obs = 0.0;
        bool have_sig_obs = false;
        if (g_obs_ul) {
            sig_obs = g_obs_ul->Eval(mchi);
            if (std::isfinite(sig_obs) && sig_obs > 0.0) have_sig_obs = true;
            std::cout << std::scientific << std::setprecision(6)
                      << "[limit][band-per-toy] m_{#chi}=" << mchi
                      << " MeV  Obs.#sigma_{UL}(data)=" << sig_obs
                      << " cm^{2}  have_obs=" << (have_sig_obs ? "yes" : "no")
                      << "\n";
        }
        if (draw_toy) widen(vt);
        if (draw_asy) widen(va);
        if (have_sig_obs) {
            smin = std::min(smin, sig_obs);
            smax = std::max(smax, sig_obs);
        }
        if (!(smax > 0.0) || !std::isfinite(smin)) {
            TLatex tx;
            tx.SetNDC();
            tx.SetTextSize(0.06);
            tx.DrawLatex(0.15, 0.5, "no entries");
            continue;
        }

        const int ntoy_eff = draw_toy ? static_cast<int>(vt.size()) : 0;
        const int nasy_eff = draw_asy ? static_cast<int>(va.size()) : 0;
        const int n_for_bins = std::max(ntoy_eff, nasy_eff);
        const int nb_local = std::min(
            180, std::max(24, static_cast<int>(std::lround(2.0 * std::sqrt(
                                    static_cast<double>(std::max(n_for_bins, 1)))))));
        double x_lo = smin * 0.88;
        double x_hi = smax * 1.12;
        if (!(x_lo > 0.0) || x_lo >= x_hi) {
            x_lo = smax * 1e-4;
            x_hi = smax * 1.12;
        }
        // Extra margin on log-x so vertical lines do not sit on the frame.
        x_lo /= 1.03;
        x_hi *= 1.03;

        TH1D* h_toy = nullptr;
        TH1D* h_asy = nullptr;

        if (draw_toy && !vt.empty()) {
            h_toy = new TH1D(Form("h_band_toy_%d", ip),
                              ";#bar{#sigma}_{e} [cm^{2}];toys / bin",
                              nb_local, x_lo, x_hi);
            for (double s : vt) h_toy->Fill(s);
            h_toy->SetLineColor(kBlue + 2);
            h_toy->SetLineWidth(2);
            h_toy->SetFillColorAlpha(kBlue + 1, 0.35);
            h_toy->SetMinimum(0);
            h_toy->Draw("HIST");
        }
        if (draw_asy && !va.empty()) {
            h_asy = new TH1D(Form("h_band_asy_%d", ip),
                             ";#bar{#sigma}_{e} [cm^{2}];toys / bin",
                             nb_local, x_lo, x_hi);
            for (double s : va) h_asy->Fill(s);
            h_asy->SetLineColor(kRed + 1);
            h_asy->SetLineWidth(2);
            h_asy->SetFillStyle(0);
            h_asy->SetMinimum(0);
            const char* opt = (h_toy && h_toy->GetEntries() > 0) ? "HIST SAME" : "HIST";
            h_asy->Draw(opt);
        }

        // Keep each pad on its local range to avoid compressed low-mass panels.
        if (h_toy) h_toy->GetXaxis()->SetRangeUser(x_lo, x_hi);
        if (h_asy) h_asy->GetXaxis()->SetRangeUser(x_lo, x_hi);

        TLatex cap;
        cap.SetNDC();
        cap.SetTextSize(0.045);
        cap.DrawLatex(0.12, 0.88, Form("m_{#chi} = %.4g MeV", mchi));
        if (draw_toy && ntoy_eff > 0) {
            TLatex nlab;
            nlab.SetNDC();
            nlab.SetTextSize(0.035);
            nlab.DrawLatex(0.12, 0.82, Form("N_{toy}=%d (toy-MC #sigma_{UL})", ntoy_eff));
        }
        if (draw_asy && nasy_eff > 0) {
            TLatex nlab2;
            nlab2.SetNDC();
            nlab2.SetTextSize(0.035);
            nlab2.DrawLatex(0.12, draw_toy && ntoy_eff > 0 ? 0.76 : 0.82,
                           Form("N=%d (asympt. pass #sigma_{UL})", nasy_eff));
        }
        if (have_sig_obs) {
            TLatex tobs;
            tobs.SetNDC();
            tobs.SetTextColor(kGreen + 2);
            tobs.SetTextSize(0.03);
            tobs.DrawLatex(0.12, 0.68, "Obs. #sigma_{UL} (data)");
        }

        gPad->Modified();
        gPad->Update();
        const double y_u0 = gPad->GetUymin();
        const double y_u1 = gPad->GetUymax();

        bool drew_median = false;
        if (draw_toy && !vt.empty()) {
            const double med = median_sorted_copy(vt);
            if (std::isfinite(med) && med > 0.0) {
                draw_world_vline(med, y_u0, y_u1, kBlack, 2, 2);
                drew_median = true;
            }
        }
        if (!drew_median && draw_asy && !va.empty()) {
            const double med = median_sorted_copy(va);
            if (std::isfinite(med) && med > 0.0)
                draw_world_vline(med, y_u0, y_u1, kBlack, 2, 2);
        }
        if (have_sig_obs) {
            std::cout << std::scientific << std::setprecision(6)
                      << "[limit][band-per-toy]   -> vertical line at #bar{#sigma}_{e}="
                      << sig_obs << " cm^{2}, y=" << y_u0 << ".." << y_u1
                      << "  pad=" << (gPad ? gPad->GetName() : "?") << "\n";
            draw_world_vline(sig_obs, y_u0, y_u1, kGreen + 2, 1, 3);
        }

        if (h_toy) {
            h_toy->SetDirectory(nullptr);
            written_for_root.push_back(h_toy);
        }
        if (h_asy) {
            h_asy->SetDirectory(nullptr);
            written_for_root.push_back(h_asy);
        }
    }

    cpt->cd(0);
    cpt->Update();
    if (!gROOT->IsBatch()) {
        cpt->Draw(); // keep heap canvas + open window until app.Run() returns
    }
    const std::string out_pdf = "dme_band_per_toy_sigma.pdf";
    cpt->SaveAs(out_pdf.c_str());
    if (gROOT->IsBatch()) {
        delete cpt;
    }

    TFile fout("dme_band_per_toy_sigma.root", "RECREATE");
    for (TH1D* h : written_for_root) {
        if (h) h->Write();
    }
    fout.Close();

    std::cout << std::defaultfloat << std::setprecision(6);
    std::cout << "[limit][band-per-toy] Saved: " << out_pdf << "\n"
              << "[limit][band-per-toy]         dme_band_per_toy_sigma.root\n";
    return 0;
}

// TGraph* MakeExactEnvelope(const std::vector<const TGraph*>& graphs,
//                           double xmin, double xmax,
//                           int n_samples = 2000)
// {
//     if (graphs.empty()) return nullptr;

//     TGraph* g_env = new TGraph();
//     g_env->SetName("g_env");
//     g_env->SetTitle("");

//     double logxmin = std::log10(xmin);
//     double logxmax = std::log10(xmax);

//     int ip = 0;
//     for (int i = 0; i < n_samples; ++i) {
//         double logx = logxmin + (logxmax - logxmin) * i / (n_samples - 1);
//         double x    = std::pow(10.0, logx);

//         double y_min = std::numeric_limits<double>::infinity();

//         for (const TGraph* g : graphs) {
//             if (!g || g->GetN() < 2) continue;

//             // only use graphs that cover this x
//             double gx0, gy0, gx1, gy1;
//             g->GetPoint(0, gx0, gy0);
//             g->GetPoint(g->GetN()-1, gx1, gy1);
//             if (x < std::min(gx0,gx1) || x > std::max(gx0,gx1)) continue;

//             double y = g->Eval(x, 0, "S");  // linear interpolation
//             if (y > 0.0 && y < y_min) y_min = y;
//         }

//         if (y_min < std::numeric_limits<double>::infinity()) {
//             g_env->SetPoint(ip++, x, y_min);
//         }
//     }

//     return g_env;
// }

TGraph* MakeLowerEnvelope(const std::vector<const TGraph*>& graphs,
                          double xmin, double xmax,
                          int n_samples = 2000)
{
    if (graphs.empty()) return nullptr;

    TGraph* g_env = new TGraph();
    g_env->SetName("g_env");
    g_env->SetTitle("");

    double logxmin = std::log10(xmin);
    double logxmax = std::log10(xmax);

    int ip = 0;

    for (int i = 0; i < n_samples; ++i) {
        double logx = logxmin + (logxmax - logxmin) * i / (n_samples - 1);
        double x    = std::pow(10.0, logx);

        double y_min = std::numeric_limits<double>::infinity();

        for (const TGraph* g : graphs) {
            if (!g || g->GetN() < 2) continue;

            // --- true x-range of this graph ---
            double gx_min =  std::numeric_limits<double>::infinity();
            double gx_max = -std::numeric_limits<double>::infinity();
            for (int j = 0; j < g->GetN(); ++j) {
                double gx, gy;
                g->GetPoint(j, gx, gy);
                if (gx < gx_min) gx_min = gx;
                if (gx > gx_max) gx_max = gx;
            }

            if (x < gx_min || x > gx_max) continue;  // x outside this curve

            double y = g->Eval(x, 0, "S");  // linear interpolation
            if (y > 0.0 && y < y_min) {
                y_min = y;
            }
        }

        if (y_min < std::numeric_limits<double>::infinity()) {
            g_env->SetPoint(ip++, x, y_min);
        }
    }

    return g_env;
}

// ----------------------------------------------------------------------------
// MakeExactEnvelope
//   Lower envelope of several exclusion curves: at every x of any curve I take the smallest positive
//   y among the curves that cover that x (interpolated). Returns a new graph owned by the caller.
// ----------------------------------------------------------------------------
TGraph* MakeExactEnvelope(const std::vector<const TGraph*>& graphs)
{
    std::set<double> all_x;

    // Collect all x-values from all curves
    for (const TGraph* g : graphs) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);
            all_x.insert(x);
        }
    }

    // Build envelope graph
    TGraph* env = new TGraph();
    env->SetName("g_env");
    env->SetTitle("");

    int ip = 0;

    for (double x : all_x) {
        double y_min = std::numeric_limits<double>::infinity();

        for (const TGraph* g : graphs) {
            if (!g) continue;

            // Check if x is inside the graph range
            double gx0, gy0, gx1, gy1;
            g->GetPoint(0, gx0, gy0);
            g->GetPoint(g->GetN()-1, gx1, gy1);

            if (x < std::min(gx0,gx1) || x > std::max(gx0,gx1)) 
                continue;

            double y = g->Eval(x, 0, "S");
            if (y > 0.0 && y < y_min)
                y_min = y;
        }

        if (y_min < std::numeric_limits<double>::infinity()) {
            env->SetPoint(ip++, x, y_min);
        }
    }

    return env;
}



// ----------------------------------------------------------------------------
// EnforceMonotonicEnvelope
//   Remove upward wiggles: going along the points, y is never allowed to rise above the previous value (in place).
// ----------------------------------------------------------------------------
void EnforceMonotonicEnvelope(TGraph* g)
{
    if (!g) return;
    int n = g->GetN();
    if (n < 3) return;

    double prev_x, prev_y;
    g->GetPoint(0, prev_x, prev_y);

    for (int i = 1; i < n; ++i) {
        double x, y;
        g->GetPoint(i, x, y);

        // enforce monotonic behaviour (no oscillatory dips)
        if (y > prev_y) {
            y = prev_y; // clamp downward only
            g->SetPoint(i, x, y);
        }

        prev_y = y;
    }
}

// ----------------------------------------------------------------------------
// LoadExclusionCSV
//   Read an exclusion curve from a two-column text file ("x, y", "x y", ";" or tab separated; comment
//   and header lines are skipped), scale the columns and style the line. Returns nullptr if the file
//   cannot be opened or has no valid points.
// ----------------------------------------------------------------------------
TGraph* LoadExclusionCSV(const std::string& csv_path,
                         int line_color = kGray+2,
                         int line_style = 1,
                         int line_width = 2,
                        double scale_x = 1.0,
                        double scale_y = 1.0)
{
  std::ifstream in(csv_path);
  if (!in.is_open()) {
    std::cerr << "[exclusion] ERROR: cannot open CSV file " << csv_path << "\n";
    return nullptr;
  }

  std::vector<double> xs;
  std::vector<double> ys;
  std::string line;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    // Skip comments or header lines starting with '#' or letters
    if (line[0] == '#' || std::isalpha(static_cast<unsigned char>(line[0]))) {
      continue;
    }

    std::istringstream ss(line);
    double x = 0.0, y = 0.0;
    char sep = 0;

    // Very simple "x, y" or "x y" parser
    if (!(ss >> x)) continue;
    if (ss.peek() == ',' || ss.peek() == ';' || ss.peek() == '\t') {
      ss >> sep;
    }
    if (!(ss >> y)) continue;

    xs.push_back(x);
    ys.push_back(y);
  }

  if (xs.empty()) {
    std::cerr << "[exclusion] WARNING: no valid points found in " << csv_path << "\n";
    return nullptr;
  }

  auto* g = new TGraph(static_cast<int>(xs.size()));
  g->SetName(("g_excl_" + csv_path).c_str());
  g->SetTitle("");

  for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
    g->SetPoint(i, xs[i]*scale_x, ys[i]*scale_y);
  }

  g->SetLineColor(line_color);
  g->SetLineStyle(line_style);
  g->SetLineWidth(line_width);

  return g;
}

// Dark photon literature CSVs use several column layouts; mass is in keV, ε in
// column k/limit.  x_col/y_col are 0-based field indices after comma-splitting.
// True if the first significant character of the line is a letter (a header row).
static bool LineLooksLikeCSVHeader(const std::string& line)
{
  for (char c : line) {
    if (c == ',' || c == ' ' || c == '\t' || c == '\r') continue;
    return std::isalpha(static_cast<unsigned char>(c));
  }
  return false;
}

// ----------------------------------------------------------------------------
// SplitCSVDoubles
//   Split a line at commas into numbers, silently skipping empty or non-numeric fields.
// ----------------------------------------------------------------------------
static std::vector<double> SplitCSVDoubles(const std::string& line)
{
  std::vector<double> vals;
  std::string field;
  std::stringstream ss(line);
  while (std::getline(ss, field, ',')) {
    if (field.empty()) continue;
    try {
      vals.push_back(std::stod(field));
    } catch (...) {
      // skip non-numeric tokens
    }
  }
  return vals;
}

TGraph* LoadDarkPhotonLimitCSV(const std::string& csv_path,
                             int x_col,
                             int y_col,
                             double x_scale,
                             Color_t line_color = kGray + 2,
                             Style_t line_style = 1,
                             Width_t line_width = 2,
                             double y_scale = 1.0)
{
  std::ifstream in(csv_path);
  if (!in.is_open()) {
    std::cerr << "[exclusion] ERROR: cannot open CSV file " << csv_path << "\n";
    return nullptr;
  }

  std::vector<double> xs;
  std::vector<double> ys;
  std::string line;

  while (std::getline(in, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;
    if (LineLooksLikeCSVHeader(line)) continue;

    const std::vector<double> fields = SplitCSVDoubles(line);
    if (x_col < 0 || y_col < 0 ||
        static_cast<std::size_t>(x_col) >= fields.size() ||
        static_cast<std::size_t>(y_col) >= fields.size()) {
      continue;
    }

    const double x = fields[static_cast<std::size_t>(x_col)] * x_scale;
    const double y = fields[static_cast<std::size_t>(y_col)] * y_scale;
    if (!(x > 0.0 && y > 0.0 && std::isfinite(x) && std::isfinite(y))) continue;

    xs.push_back(x);
    ys.push_back(y);
  }

  if (xs.empty()) {
    std::cerr << "[exclusion] WARNING: no valid points found in " << csv_path << "\n";
    return nullptr;
  }

  auto* g = new TGraph(static_cast<int>(xs.size()));
  g->SetName(("g_dp_" + csv_path).c_str());
  g->SetTitle("");
  for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
    g->SetPoint(i, xs[i], ys[i]);
  }
  g->SetLineColor(line_color);
  g->SetLineStyle(line_style);
  g->SetLineWidth(line_width);
  return g;
}

// The XENONnT HP bracket as two curves (lo, hi) and a band (TGraphAsymmErrors) between them.
struct DarkPhotonXenonNtBracket {
  TGraph* lo = nullptr;
  TGraph* hi = nullptr;
  TGraphAsymmErrors* band = nullptr;
};

// ----------------------------------------------------------------------------
// LoadDarkPhotonXenonNtHPBracket
//   Read a CSV with mass, lower and upper epsilon columns (mass multiplied by mass_scale; rows with
//   non-positive values skipped) and build the lower/upper curves plus the band between them.
// ----------------------------------------------------------------------------
static DarkPhotonXenonNtBracket LoadDarkPhotonXenonNtHPBracket(
    const std::string& csv_path,
    double mass_scale,
    Color_t color,
    Style_t line_style = 7,
    Width_t line_width = 3)
{
  DarkPhotonXenonNtBracket out;
  std::ifstream in(csv_path);
  if (!in.is_open()) {
    std::cerr << "[exclusion] ERROR: cannot open CSV file " << csv_path << "\n";
    return out;
  }

  std::vector<double> xs;
  std::vector<double> ylo;
  std::vector<double> yhi;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    if (LineLooksLikeCSVHeader(line)) continue;
    const std::vector<double> fields = SplitCSVDoubles(line);
    if (fields.size() < 3) continue;
    const double x = fields[0] * mass_scale;
    const double lo = fields[1];
    const double hi = fields[2];
    if (!(x > 0.0 && lo > 0.0 && hi > 0.0)) continue;
    xs.push_back(x);
    ylo.push_back(lo);
    yhi.push_back(hi);
  }
  if (xs.empty()) {
    std::cerr << "[exclusion] WARNING: no valid XENONnT HP bracket points in "
              << csv_path << "\n";
    return out;
  }

  out.lo = new TGraph(static_cast<int>(xs.size()));
  out.hi = new TGraph(static_cast<int>(xs.size()));
  out.band = new TGraphAsymmErrors(static_cast<int>(xs.size()));
  for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
    out.lo->SetPoint(i, xs[i], ylo[i]);
    out.hi->SetPoint(i, xs[i], yhi[i]);
    const double ymid = 0.5 * (ylo[i] + yhi[i]);
    const double dy = 0.5 * std::abs(yhi[i] - ylo[i]);
    out.band->SetPoint(i, xs[i], ymid);
    out.band->SetPointError(i, 0.0, 0.0, dy, dy);
  }

  for (TGraph* g : {out.lo, out.hi}) {
    g->SetLineColor(color);
    g->SetLineStyle(line_style);
    g->SetLineWidth(line_width);
  }
  out.band->SetFillColor(color);
  out.band->SetFillStyle(3005);
  out.band->SetLineColor(color);
  out.band->SetLineWidth(0);
  return out;
}



// ----------------------------------------------------------------------------
// LoadMassLimitDAT
//   Read a "mass,limit" file (first line is a header), multiply the limit by yscale and style the line.
// ----------------------------------------------------------------------------
TGraph* LoadMassLimitDAT(const std::string& filename,
                         double yscale = 1,   // multiply Y column by this factor
                         Color_t  color  = kBlack,
                         Style_t  lstyle = 1,
                         Width_t  lwidth = 2)
{
  std::ifstream in(filename);
  if (!in.is_open()) {
    Error("LoadMassLimitDAT", "Cannot open file '%s'", filename.c_str());
    return nullptr;
  }

  std::vector<double> vx, vy;
  std::string line;

  bool first = true;

  while (std::getline(in, line)) {
    if (line.empty()) continue;

    // Skip header line: "mass_mev,limit"
    if (first) {
      first = false;
      continue;
    }

    std::stringstream ss(line);

    double mass = 0.0;
    double limit = 0.0;
    char comma;

    // Format: mass,limit
    if (!(ss >> mass)) continue;
    if (ss.peek() == ',' || ss.peek() == ';') ss >> comma;
    if (!(ss >> limit)) continue;

    vx.push_back(mass);
    vy.push_back(limit * yscale);   // apply scaling here
  }

  auto gr = new TGraph(vx.size());
  for (int i = 0; i < (int)vx.size(); i++)
    gr->SetPoint(i, vx[i], vy[i]);

  gr->SetLineColor(color);
  gr->SetLineStyle(lstyle);
  gr->SetLineWidth(lwidth);

  return gr;
}

// need a function that reads in from txt file and makes a graph

// ----------------------------------------------------------------------------
// LoadTxtToTGraph
//   Read an "x,y" text file (comment and blank lines skipped). x values above 1e4 are treated
//   as masses in eV and converted to MeV. Returns nullptr if the file cannot be opened.
// ----------------------------------------------------------------------------
TGraph* LoadTxtToTGraph(const std::string& filename,
                        Color_t  color  = kBlack,
                        Style_t  lstyle = 1,
                        Width_t  lwidth = 2)
{
  std::ifstream in(filename);
  if (!in.is_open()) {
    Error("LoadTxtToTGraph", "Cannot open file '%s'", filename.c_str());
    return nullptr;
  }

  std::vector<double> vx, vy;
  std::string line;

  while (std::getline(in, line)) {
    if (line.empty())            continue;
    if (line[0] == '#')          continue;   // skip comments
    if (line.find_first_not_of(" \t\r\n") == std::string::npos)
      continue;                              // skip whitespace-only

    std::stringstream ss(line);

    double x = 0.0;
    double y = 0.0;
    char comma;

    // Expect "x,y" (e.g. "1.0e+06,2.0e-28")
    if (!(ss >> x)) continue;               // failed to read x → skip line

    // Optional comma (or other separator)
    if (ss.peek() == ',' || ss.peek() == ';')
      ss >> comma;

    if (!(ss >> y)) continue;               // failed to read y → skip line

    // mass_ev column (paper export): convert to MeV for axis
    if (x > 1.0e4) x *= 1.0e-6;

    vx.push_back(x);
    vy.push_back(y);
  }

  auto gr = new TGraph(static_cast<int>(vx.size()));
  for (int i = 0; i < static_cast<int>(vx.size()); ++i) {
    gr->SetPoint(i, vx[i], vy[i]);
  }

  gr->SetLineColor(color);
  gr->SetLineStyle(lstyle);
  gr->SetLineWidth(lwidth);

  return gr;
}

// Set the visible axis window without TGraph::SetLimits (which blocks ROOT GUI zoom).
// Set the visible axis window of the current pad through the frame histogram.
static void ApplyPadAxisRangeUser(double x_lo, double x_hi, double y_lo, double y_hi)
{
  if (!gPad) return;
  TH1* frame = nullptr;
  TList* prims = gPad->GetListOfPrimitives();
  if (prims) {
    TIter next(prims);
    TObject* obj = nullptr;
    while ((obj = next())) {
      if (obj->InheritsFrom(TH1::Class())) {
        frame = static_cast<TH1*>(obj);
        break;
      }
    }
  }
  if (!frame) return;
  if (x_lo > 0.0 && x_hi > x_lo) {
    frame->GetXaxis()->SetRangeUser(x_lo, x_hi);
  }
  if (y_lo > 0.0 && y_hi > y_lo) {
    frame->GetYaxis()->SetRangeUser(y_lo, y_hi);
  }
  gPad->Modified();
  gPad->Update();
}

// ----------------------------------------------------------------------------
// UpdateRangesFromGraph
//   Widen the running (xmin, xmax, ymin, ymax) so that it contains every point of g.
// ----------------------------------------------------------------------------
void UpdateRangesFromGraph(double& xmin, double& xmax,
                           double& ymin, double& ymax,
                           const TGraph* g)
{
  if (!g) return;
  int n = g->GetN();
  for (int i = 0; i < n; ++i) {
    double x, y;
    g->GetPoint(i, x, y);
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  }
}

// ----------------------------------------------------------------------------
// UpdateRangesFromGraphWindow
//   Same as UpdateRangesFromGraph but only for points with x inside [xwin_lo, xwin_hi].
// ----------------------------------------------------------------------------
static void UpdateRangesFromGraphWindow(double& xmin, double& xmax,
                                      double& ymin, double& ymax,
                                      const TGraph* g,
                                      double xwin_lo, double xwin_hi)
{
  if (!g) return;
  for (int i = 0; i < g->GetN(); ++i) {
    double x, y;
    g->GetPoint(i, x, y);
    if (x < xwin_lo || x > xwin_hi) continue;
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  }
}

// ----------------------------------------------------------------------------
// MakeFilledBand
//   Polygon that fills the area between curve g and the level y_bottom (points sorted in x), with transparency alpha and no outline.
// ----------------------------------------------------------------------------
TGraph* MakeFilledBand(const TGraph* g,
                       double y_bottom,
                       int fill_color,
                       double alpha = 0.25)
{
  if (!g) return nullptr;

  const int n = g->GetN();
  if (n < 2) return nullptr;

  // Copy points and sort in x (just to be safe)
  std::vector<std::pair<double,double>> pts(n);
  double x, y;
  for (int i = 0; i < n; ++i) {
    g->GetPoint(i, x, y);
    pts[i] = std::make_pair(x, y);
  }
  std::sort(pts.begin(), pts.end(),
            [](const auto& a, const auto& b){ return a.first < b.first; });

  // Build closed polygon: curve + two points at y_bottom
  TGraph* gf = new TGraph(n + 2);
  for (int i = 0; i < n; ++i) {
    gf->SetPoint(i, pts[i].first, pts[i].second);
  }
  gf->SetPoint(n,   pts.back().first, y_bottom);  // bottom-right
  gf->SetPoint(n+1, pts.front().first, y_bottom); // bottom-left

  gf->SetName((std::string(g->GetName()) + "_fill").c_str());
  gf->SetTitle("");

  gf->SetFillStyle(1001);
  gf->SetFillColorAlpha(fill_color, alpha);
  gf->SetLineWidth(0); // no outline for the fill

  return gf;
}

// ----------------------------------------------------------------------------
// MakeFilledBandAbove
//   Polygon that fills the area between curve g and the level y_top over [x_left, x_right] (default:
//   the range of g); where the range goes beyond the curve the gap is filled up to y_top.
// ----------------------------------------------------------------------------
TGraph* MakeFilledBandAbove(const TGraph* g,
                            double y_top,
                            int fill_color,
                            double alpha = 0.25,
                            double x_left  = -1.0,
                            double x_right = -1.0)
{
  if (!g) return nullptr;

  const int n = g->GetN();
  if (n < 2) return nullptr;

  std::vector<std::pair<double,double>> pts(n);
  double x, y;
  for (int i = 0; i < n; ++i) {
    g->GetPoint(i, x, y);
    pts[i] = std::make_pair(x, y);
  }
  std::sort(pts.begin(), pts.end(),
            [](const auto& a, const auto& b){ return a.first < b.first; });

  const double xl = (x_left  > 0.0) ? x_left  : pts.front().first;
  const double xr = (x_right > 0.0) ? x_right : pts.back().first;

  // Build polygon that fills between the curve and y_top over [xl, xr].
  // When xl/xr extend beyond the curve, the gap region is filled to y_top.
  //
  // Vertex order:
  //   curve drop-in at left (only if xl < pts.front().x)
  //   curve points left → right
  //   curve drop-out at right (only if xr > pts.back().x)
  //   close along top from xr → xl
  bool extend_left  = (xl < pts.front().first);
  bool extend_right = (xr > pts.back().first);
  int extra = (extend_left ? 1 : 0) + (extend_right ? 1 : 0) + 2;  // +2 for top corners
  TGraph* gf = new TGraph(n + extra);

  int ip = 0;
  if (extend_left)  gf->SetPoint(ip++, pts.front().first, y_top);  // drop from top to curve start
  for (int i = 0; i < n; ++i) gf->SetPoint(ip++, pts[i].first, pts[i].second);
  if (extend_right) gf->SetPoint(ip++, pts.back().first, y_top);   // rise from curve end to top
  gf->SetPoint(ip++, xr, y_top);  // top-right
  gf->SetPoint(ip++, xl, y_top);  // top-left (ROOT closes back to first point)

  gf->SetName((std::string(g->GetName()) + "_fill_above").c_str());
  gf->SetTitle("");

  gf->SetFillStyle(1001);
  gf->SetFillColorAlpha(fill_color, alpha);
  gf->SetLineWidth(0);

  return gf;
}

// ----------------------------------------------------------------------------
// LabelGraphAtX
//   Text label placed next to the point of g nearest to x_label (compared in log x), 25% above the curve.
// ----------------------------------------------------------------------------
TLatex* LabelGraphAtX(const TGraph* g,
                      double x_label,
                      const char* text,
                      int color,
                      double size = 0.035)
{
    if (!g) return nullptr;

    // Find nearest point in graph
    int n = g->GetN();
    double best_dx = 1e99;
    double best_x = 0, best_y = 0;

    for (int i = 0; i < n; ++i) {
        double x, y;
        g->GetPoint(i, x, y);
        double dx = fabs(log(x) - log(x_label));  // compare in log-space
        if (dx < best_dx) {
            best_dx = dx;
            best_x = x;
            best_y = y;
        }
    }

    // Slight vertical offset (above curve)
    best_y *= 1.25;

    TLatex* t = new TLatex(best_x, best_y, text);
    t->SetTextColor(color);
    t->SetTextSize(size);
    t->SetTextAlign(12); // left-aligned
    t->SetTextFont(42);
    t->SetNDC(false);

    return t;
}

// graphs: all exclusion curves that define the union exclusion
// xmin,xmax: mass range for the band (same as my plot axes)
// y_top: top of fill (same as my plot max)
// fill_color/alpha: style
TGraph* MakeUnionBandAbove(const std::vector<const TGraph*>& graphs,
                           double xmin, double xmax,
                           double y_top,
                           int fill_color,
                           double alpha = 0.25)
{
    if (graphs.empty()) return nullptr;

    // 1) Build the lower envelope sampled densely in log x
    const int N = 4000;  // can push to 6000 if needed
    TGraph* g_env = new TGraph();
    g_env->SetName("g_env_union");
    g_env->SetTitle("");

    double logxmin = std::log10(xmin);
    double logxmax = std::log10(xmax);

    int ip = 0;
    for (int i = 0; i < N; ++i) {
        double logx = logxmin + (logxmax - logxmin) * i / (N - 1);
        double x    = std::pow(10.0, logx);

        double y_min = std::numeric_limits<double>::infinity();

        for (const TGraph* g : graphs) {
            if (!g || g->GetN() < 2) continue;

            double gx0, gy0, gx1, gy1;
            g->GetPoint(0, gx0, gy0);
            g->GetPoint(g->GetN() - 1, gx1, gy1);

            if (x < std::min(gx0, gx1) || x > std::max(gx0, gx1))
                continue;

            double y = g->Eval(x, 0, "S");  // linear interpolation
            if (y > 0.0 && y < y_min)
                y_min = y;
        }

        if (y_min < std::numeric_limits<double>::infinity()) {
            // 2) tiny downward fudge so we never sit above the true best curve
            y_min *= 0.98;  // try 0.99 for a tighter fill
            g_env->SetPoint(ip++, x, y_min);
        }
    }

    // 3) Turn that envelope into a filled band up to y_top
    TGraph* g_fill = MakeFilledBandAbove(g_env, y_top, fill_color, alpha);
    g_fill->SetName("g_union_fill");

    return g_fill;
}



// Extract a 90% CL upper-limit TGraph from a scan ROOT file by scanning the q TH2D
// column-by-column (same as --from-qhist), so band curves are consistent with the central line.
TGraph* ExtractLimitCurveFromRoot(const std::string& path, double q_thr = 2.71) {
    TFile* f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "[dc-band] ERROR: cannot open " << path << "\n";
        return nullptr;
    }
    TH2D* h = dynamic_cast<TH2D*>(f->Get("q_mchi_sigma_pattern"));
    if (!h) h = dynamic_cast<TH2D*>(f->Get("q_mchi_sigma"));
    if (!h) {
        std::cerr << "[dc-band] ERROR: no q TH2D in " << path << "\n";
        f->Close(); return nullptr;
    }
    h->SetDirectory(nullptr);
    f->Close();
    std::vector<double> xs, ys;
    for (int ix = 1; ix <= h->GetNbinsX(); ++ix) {
        double mchi = h->GetXaxis()->GetBinCenter(ix);
        double sigma_ul = -1;
        for (int iy = 1; iy <= h->GetNbinsY(); ++iy) {
            if (h->GetBinContent(ix, iy) > q_thr) {
                sigma_ul = h->GetYaxis()->GetBinCenter(iy); break;
            }
        }
        if (sigma_ul > 0) { xs.push_back(mchi); ys.push_back(sigma_ul); }
    }
    if (xs.empty()) return nullptr;
    return new TGraph(xs.size(), xs.data(), ys.data());
}

// Log-linear interpolation of TGraph at x.
static double InterpGraphLogLog(TGraph* g, double x) {
    int n = g->GetN();
    double* gx = g->GetX(); double* gy = g->GetY();
    if (x <= gx[0])   return gy[0];
    if (x >= gx[n-1]) return gy[n-1];
    for (int i = 0; i < n-1; ++i) {
        if (x >= gx[i] && x <= gx[i+1]) {
            double t = (std::log(x)-std::log(gx[i])) / (std::log(gx[i+1])-std::log(gx[i]));
            return std::exp(std::log(gy[i]) + t*(std::log(gy[i+1])-std::log(gy[i])));
        }
    }
    return gy[n-1];
}

// Build a closed filled-polygon band whose bottom = point-wise min and top = point-wise max
// over g_lo, g_hi, and g_cen. Uses the INTERSECTION x-range of all three curves and clips
// y values to [ymin_clip, ymax_clip] to avoid ROOT hatching artefacts from off-axis polygons.
TGraph* MakeDcBandPolygon(TGraph* g_lo, TGraph* g_hi, TGraph* g_cen,
                           double ymin_clip = 1e-40, double ymax_clip = 1e40) {
    // Intersection of x-ranges: only build where all three curves have coverage
    auto xrange = [](TGraph* g) -> std::pair<double,double> {
        return {g->GetX()[0], g->GetX()[g->GetN()-1]};
    };
    double xlo = std::max({xrange(g_lo).first,  xrange(g_hi).first,  xrange(g_cen).first});
    double xhi = std::min({xrange(g_lo).second, xrange(g_hi).second, xrange(g_cen).second});
    if (xlo >= xhi) return nullptr;

    std::vector<double> allx;
    for (TGraph* g : {g_lo, g_hi, g_cen}) {
        for (int i = 0; i < g->GetN(); ++i) {
            double x, y; g->GetPoint(i, x, y);
            if (x >= xlo && x <= xhi) allx.push_back(x);
        }
    }
    std::sort(allx.begin(), allx.end());
    allx.erase(std::unique(allx.begin(), allx.end()), allx.end());

    int n = allx.size();
    std::vector<double> lo(n), hi(n);
    for (int i = 0; i < n; ++i) {
        double x = allx[i];
        double a = InterpGraphLogLog(g_lo,  x);
        double b = InterpGraphLogLog(g_hi,  x);
        double c = InterpGraphLogLog(g_cen, x);
        lo[i] = std::max(std::min({a,b,c}), ymin_clip);
        hi[i] = std::min(std::max({a,b,c}), ymax_clip);
    }
    std::vector<double> bx, by;
    for (int i = 0; i < n; ++i)    { bx.push_back(allx[i]); by.push_back(lo[i]); }
    for (int i = n-1; i >= 0; --i) { bx.push_back(allx[i]); by.push_back(hi[i]); }
    bx.push_back(allx[0]); by.push_back(lo[0]);
    return new TGraph(bx.size(), bx.data(), by.data());
}

// ----------------------------------------------------------------------------
// main
//   Plot the exclusion limit from one or more scan ROOT files (each optionally followed by a label),
//   with options for dark-photon and Migdal plots, the Brazilian band from ccdarksens_band, DC bands,
//   literature curves, axis ranges and output paths (see the usage text for the flags). Flow: open the
//   files and get q(m_chi, sigma), extract the limit curves, overlay the band and external exclusions,
//   set the axis ranges from all curves, draw the plot and write the PDF, ROOT and CSV outputs.
// ----------------------------------------------------------------------------
int main(int argc, char** argv)
{
  // ---------------------------------------------------------------------------
  // CLI notes (latest update):
  //   - Supports multiple input ROOT files, each optionally followed by a label:
  //       scan1.root "label1" scan2.root "label2" ...
  //   - Supports --from-qhist to force UL extraction from q(mchi,sigma)
  //     histogram even when upper_limit_sigma_e_mchi exists.
  //   - Falls back to legacy single-file behavior when only one input file is
  //     provided (or inferred from argv[1]).
  // ---------------------------------------------------------------------------



  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " scan_dmelectron_grid.root [\"label\"] [q_threshold] [mediator]\n"
              << "  Optional flags:\n"
              << "    --from-qhist    Build limit from q(mchi,sigma) histogram even if upper_limit_sigma_e_mchi exists.\n"
              << "    --draw-both     Draw stored UL and pydme bisection diagnostic (or q-map if bisection graph absent).\n"
              << "    --dark-photon   Dark photon absorption plot (m_{A'} vs #varepsilon); "
              << "loads data/previous_limits/dark_photon/ (.csv and .txt).\n"
              << "    --dc-band-lo <root>  (dark-photon) Lower band edge scan (best DC). Combined with --dc-band-hi,\n"
              << "                         draws a hatched band and only the central (positional) line.\n"
              << "    --dc-band-hi <root>  Upper band edge scan (worst DC).\n"
              << "    --migdal        Migdal effect plot (m_{chi} vs sigma_n); "
              << "loads data/previous_limits/Migdal_heavy/ (combine with --heavy or --light).\n"
              << "    --wimp          WIMP-nucleon SI elastic plot (m_{chi} vs sigma_n); "
              << "loads data/previous_limits/WIMP/DAMIC_2016_SNOLAB_*.csv; always implies --from-qhist.\n"
              << "    --wimp-minimal  Same as --wimp, but draws only the observed curve (black) and our scan curve(s) (red).\n"
              << "    --band <band.root>   Overlay Brazilian band (median + ribbons) from ccdarksens_band.\n"
              << "    --band-mode asymptotic|toy_mc|both   Band threshold variant (default: auto).\n"
              << "    --band-per-toy-hists   Also write dme_band_per_toy_sigma.pdf (+ .root) from TTree band_per_toy\n"
              << "                           (requires band.save_per_toy_curves=true when running ccdarksens_band).\n"
              << "    --band-per-toy-which toy|asy|both   Which per-toy #sigma_{UL} branch(es) to histogram (default: toy).\n"
              << "    --title <text>  Window / plot title (also drawn at top of main pad).\n"
              << "    --out-pdf <path>  Output PDF for main DM-e limit canvas (default: dme_limit_curve.pdf).\n"
              << "    --out-root <path> Output ROOT for main limit graphs (default: dme_limit_curve.root).\n"
              << "    --out-csv <path>  Write limit contour(s) as CSV (mchi_MeV, sigma_e_cm2; label column if multiple curves).\n"
              << "    --reference-scan <root> <legend>  Solid reference curve (e.g. Si 1.2 eV / 3.8 eV);\n"
              << "                      reuses an existing input if the path matches.\n"
              << "    --show-damic    Overlay canonical paper-export QEdark limit (ScienceRun2024_results-1).\n"
              << "    --plain-legend  Use CLI labels in the legend without (q-map) / (UL, q-grid) suffixes.\n"
              << "    --legend-left   Place legend on the left (compact, two-line labels).\n"
              << "    --legend-right  Place legend on the right (compact, two-line labels).\n"
              << "    --batch         Non-interactive: save immediately, no GUI (no axis editing).\n"
              << "  Without --batch: PDF/ROOT are written before the GUI opens; close the window to\n"
              << "  overwrite with final axis ranges (or use --batch for a one-shot save).\n"
              << "  Multiple files:\n"
              << "    scan1.root \"label1\" scan2.root \"label2\" ... [q_threshold] [mediator]\n"
              << "  mediator: 'heavy', 'light', or 'intermediate' (default: heavy). Chooses literature curves.\n"
              << "    --intermediate  Intermediate mediator plot; suppresses all literature curves.\n";
    return 1;
  }

  bool force_from_qhist = false;
  bool draw_both = false;
  bool dark_photon = false;
  bool migdal_mode = false;
  bool wimp_mode = false;
  bool wimp_minimal = false;  // --wimp-minimal: only observed (black) + our scan curve(s) (red)
  double q_thr = 2.71; // default ~90% CL, 1 dof
  std::string mediator = "heavy";

  // --band <band.root>: optional sensitivity-band file produced by
  // build/ccdarksens_band. When provided, the plotter overlays the median
  // expected and ± 1σ / ± 2σ Brazilian band on the standard exclusion plot
  // BEFORE drawing the observed limit curve(s).
  // --band-mode <asymptotic|toy_mc|both>: which threshold variant to draw.
  // Default = autodetect (if both present, draw both; else whichever is).
  std::string band_path;
  std::string band_mode_arg;
  bool band_per_toy_hists = false;
  std::string band_per_toy_which = "toy";
  std::string plot_title;
  std::string out_pdf_cli;
  std::string out_root_cli;
  std::string out_csv_cli;
  bool show_damic = false;
  bool plain_legend = false;
  bool legend_left = false;
  bool legend_right = false;
  bool smooth_curves = false;
  std::string dc_band_lo_path;  // --dc-band-lo: best DC scan (band lower edge)
  std::string dc_band_hi_path;  // --dc-band-hi: worst DC scan (band upper edge)

  std::vector<std::string> in_paths;
  std::vector<std::string> legend_labels;
  std::string reference_scan_path;
  std::string reference_scan_label;
  int reference_input_index = -1;
  size_t reference_graph_index = static_cast<size_t>(-1);
  // --store-only <rootfile> <label>: extract limit curve and write to ROOT output, but do not draw
  std::vector<std::pair<std::string, std::string>> store_only_paths;

  // Axis range overrides (optional; auto-determined from data if not set)
  double axis_x_min = -1, axis_x_max = -1, axis_y_min = -1, axis_y_max = -1;

  auto looks_like_root = [](const std::string& s) -> bool {
    if (s.size() < 5) return false;
    const std::string tail = s.substr(s.size() - 5);
    return (tail == ".root" || tail == ".ROOT");
  };
  auto lower_copy = [](std::string s) -> std::string {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
  };
  auto is_number_token = [](const std::string& tok, double& out) -> bool {
    try {
      size_t pos = 0;
      out = std::stod(tok, &pos);
      return pos == tok.size();
    } catch (...) {
      return false;
    }
  };
  auto stem_label = [](const std::string& p) -> std::string {
    std::string base = p;
    const auto slash = base.find_last_of('/');
    if (slash != std::string::npos) base = base.substr(slash + 1);
    const auto dot = base.find_last_of('.');
    if (dot != std::string::npos) base = base.substr(0, dot);
    return base.empty() ? std::string("scan") : base;
  };
  auto paths_same = [](const std::string& a, const std::string& b) -> bool {
    if (a == b) return true;
    if (a.size() >= b.size() &&
        a.compare(a.size() - b.size(), b.size(), b) == 0) {
      return true;
    }
    if (b.size() >= a.size() &&
        b.compare(b.size() - a.size(), a.size(), a) == 0) {
      return true;
    }
    return false;
  };

  // Robust parsing so calls like:
  //   <rootfile> "Diego" 2.71 heavy --from-qhist
  // don't crash by attempting std::stod("Diego").
  for (int i = 1; i < argc; ++i) {
    const std::string tok(argv[i]);
    const std::string tok_lc = lower_copy(tok);

    if (tok_lc == "--from-qhist") {
      force_from_qhist = true;
      continue;
    }
    if (tok_lc == "--draw-both") {
      draw_both = true;
      continue;
    }
    if (tok_lc == "--dark-photon") {
      dark_photon = true;
      continue;
    }
    if (tok_lc == "--migdal") {
      migdal_mode = true;
      continue;
    }
    if (tok_lc == "--wimp") {
      wimp_mode = true;
      force_from_qhist = true;  // --wimp always reads the UL from q(mchi,sigma)
      continue;
    }
    if (tok_lc == "--wimp-minimal") {
      wimp_mode = true;
      wimp_minimal = true;
      force_from_qhist = true;
      continue;
    }
    if (tok_lc == "--x-min" && i + 1 < argc) { axis_x_min = std::stod(argv[++i]); continue; }
    if (tok_lc == "--x-max" && i + 1 < argc) { axis_x_max = std::stod(argv[++i]); continue; }
    if (tok_lc == "--y-min" && i + 1 < argc) { axis_y_min = std::stod(argv[++i]); continue; }
    if (tok_lc == "--y-max" && i + 1 < argc) { axis_y_max = std::stod(argv[++i]); continue; }
    if (tok_lc == "--band") {
      if (i + 1 < argc) {
        band_path = argv[++i];
      }
      continue;
    }
    if (tok_lc == "--band-mode") {
      if (i + 1 < argc) {
        band_mode_arg = lower_copy(argv[++i]);
      }
      continue;
    }
    if (tok_lc == "--band-per-toy-hists") {
      band_per_toy_hists = true;
      continue;
    }
    if (tok_lc == "--band-per-toy-which") {
      if (i + 1 < argc) {
        band_per_toy_which = lower_copy(argv[++i]);
      }
      continue;
    }
    if (tok_lc == "--title") {
      if (i + 1 < argc) {
        plot_title = argv[++i];
      }
      continue;
    }
    if (tok_lc == "--out-pdf") {
      if (i + 1 < argc) {
        out_pdf_cli = argv[++i];
      }
      continue;
    }
    if (tok_lc == "--out-root") {
      if (i + 1 < argc) {
        out_root_cli = argv[++i];
      }
      continue;
    }
    if (tok_lc == "--out-csv") {
      if (i + 1 < argc) {
        out_csv_cli = argv[++i];
      }
      continue;
    }
    if (tok_lc == "--show-damic") {
      show_damic = true;
      continue;
    }
    if (tok_lc == "--plain-legend") {
      plain_legend = true;
      continue;
    }
    if (tok_lc == "--smooth-curve") {
      smooth_curves = true;
      continue;
    }
    if (tok_lc == "--legend-left") {
      legend_left = true;
      continue;
    }
    if (tok_lc == "--legend-right") {
      legend_right = true;
      continue;
    }
    if (tok_lc == "--dc-band-lo") {
      if (i + 1 < argc) dc_band_lo_path = argv[++i];
      continue;
    }
    if (tok_lc == "--dc-band-hi") {
      if (i + 1 < argc) dc_band_hi_path = argv[++i];
      continue;
    }
    if (tok_lc == "--reference-scan") {
      if (i + 2 < argc) {
        reference_scan_path = argv[++i];
        reference_scan_label = argv[++i];
      } else {
        std::cerr << "[limit] ERROR: --reference-scan requires <rootfile> <legend>\n";
        return 1;
      }
      continue;
    }
    if (tok_lc == "--store-only") {
      if (i + 2 < argc) {
        std::string so_path  = argv[++i];
        std::string so_label = argv[++i];
        store_only_paths.emplace_back(so_path, so_label);
      } else {
        std::cerr << "[limit] ERROR: --store-only requires <rootfile> <label>\n";
        return 1;
      }
      continue;
    }
    if (tok_lc == "heavy" || tok_lc == "--heavy") {
      mediator = "heavy";
      continue;
    }
    if (tok_lc == "intermediate" || tok_lc == "--intermediate") {
      mediator = "intermediate";
      continue;
    }
    if (tok_lc == "light" || tok_lc == "--light") {
      mediator = "light";
      continue;
    }

    if (looks_like_root(tok)) {
      in_paths.push_back(tok);

      std::string lbl = stem_label(tok);
      if (i + 1 < argc) {
        const std::string nxt(argv[i + 1]);
        const std::string nxt_lc = lower_copy(nxt);
        double tmp = 0.0;
        const bool nxt_is_flag = !nxt.empty() && nxt[0] == '-';
        const bool nxt_is_mediator = (nxt_lc == "heavy" || nxt_lc == "light");
        const bool nxt_is_root = looks_like_root(nxt);
        const bool nxt_is_number = is_number_token(nxt, tmp);

        if (!nxt_is_flag && !nxt_is_mediator && !nxt_is_root && !nxt_is_number) {
          lbl = nxt;
          ++i;  // consume label token
        }
      }
      legend_labels.push_back(lbl);
      continue;
    }

    double v = 0.0;
    if (is_number_token(tok, v)) {
      q_thr = v;
      continue;
    }
  }

  // Fallback for the legacy calling convention: first arg is the file.
  if (in_paths.empty()) {
    in_paths.push_back(argv[1]);
    legend_labels.push_back((argc >= 3) ? std::string(argv[2]) : std::string("scan"));
  }

  if (!reference_scan_path.empty()) {
    bool found = false;
    for (size_t i = 0; i < in_paths.size(); ++i) {
      if (paths_same(in_paths[i], reference_scan_path)) {
        reference_input_index = static_cast<int>(i);
        if (!reference_scan_label.empty()) legend_labels[i] = reference_scan_label;
        found = true;
        break;
      }
    }
    if (!found) {
      in_paths.push_back(reference_scan_path);
      legend_labels.push_back(reference_scan_label.empty() ? std::string("Si reference")
                                                           : reference_scan_label);
      reference_input_index = static_cast<int>(in_paths.size()) - 1;
    }
  }

  if (band_per_toy_hists && band_path.empty()) {
    std::cerr << "[limit] ERROR: --band-per-toy-hists requires --band <band.root>\n";
    return 1;
  }

  if (mediator != "heavy" && mediator != "light" && mediator != "intermediate") {
    std::cerr << "[limit] ERROR: mediator must be 'heavy', 'light', or 'intermediate', got '" << mediator << "'\n";
    return 1;
  }

  std::cout << "[limit] Inputs (" << in_paths.size() << "):\n";
  for (size_t i = 0; i < in_paths.size(); ++i) {
    const bool is_ref = (static_cast<int>(i) == reference_input_index);
    std::cout << "  - " << in_paths[i] << " (label=\"" << legend_labels[i] << "\""
              << (is_ref ? ", reference" : "") << ")\n";
  }
  std::cout << "[limit] q_threshold = " << q_thr << "\n";
  std::cout << "[limit] mediator = " << mediator << " (literature curves)\n";
  if (draw_both) {
    std::cout << "[limit] Drawing stored UL (q-grid) and pydme bisection diagnostic (or q-map fallback).\n";
  } else if (force_from_qhist) {
    std::cout << "[limit] Forcing curve build from q(mchi,sigma) histogram.\n";
  }
  if (band_per_toy_hists) {
    std::cout << "[limit] Will write per-toy #sigma_{UL} histograms (which=" << band_per_toy_which
              << ") after the main figure.\n";
  }

  // Optional: --batch saves PDF/ROOT immediately without GUI. Otherwise
  // canvases stay open for axis editing; outputs are written when windows close.
  // Strip flags (and their values) before TApplication parses argv.
  bool batch = false;
  {
    auto skip_value_flag = [&](const std::string& a_lc) -> bool {
      return a_lc == "--band" || a_lc == "--band-mode" ||
             a_lc == "--band-per-toy-which" || a_lc == "--title" ||
             a_lc == "--out-pdf" || a_lc == "--out-root" || a_lc == "--out-csv" ||
             a_lc == "--reference-scan";
    };
    int dst = 1;
    for (int src = 1; src < argc; ++src) {
      const std::string a(argv[src]);
      const std::string a_lc = lower_copy(a);
      if (a_lc == "--batch") {
        batch = true;
        continue;
      }
      if (a_lc == "--from-qhist" || a_lc == "--draw-both" ||
          a_lc == "--band-per-toy-hists" || a_lc == "--show-damic" ||
          a_lc == "--plain-legend" || a_lc == "--legend-left" ||
          a_lc == "--legend-right" || a_lc == "--dark-photon" ||
          a_lc == "--migdal" || a_lc == "--wimp" || a_lc == "--wimp-minimal" ||
          a_lc == "--intermediate") {
        continue;
      }
      if (skip_value_flag(a_lc)) {
        if (src + 1 < argc) ++src;
        continue;
      }
      argv[dst++] = argv[src];
    }
    argc = dst;
  }

  // 1) Create the ROOT application (must be before creating canvases)
  TApplication app("ccdarksens_plot_dmelectron_limit", &argc, argv);
  if (batch) {
    gROOT->SetBatch(kTRUE);
  }

  // ------------------------------------------------------------------
  // 1) Open file and get q(mchi, sigma) histogram
  // ------------------------------------------------------------------
  std::vector<double> mchi_vals;
  std::vector<double> sigma_lim_vals;
  std::vector<TGraph*> glimits;
  TGraph* glimit = nullptr; // first curve (used by legacy plotting code below)
  TH2D* hq = nullptr;       // first histogram (used for saving below)
  std::vector<std::string> legend_entries_per_graph;

  for (size_t idx_file = 0; idx_file < in_paths.size(); ++idx_file) {
    const std::string in_path = in_paths[idx_file];

    TFile* fin = TFile::Open(in_path.c_str(), "READ");
    if (!fin || fin->IsZombie()) {
      std::cerr << "[limit] ERROR: cannot open file " << in_path << "\n";
      return 1;
    }

    TH2D* hq_local = dynamic_cast<TH2D*>(fin->Get("q_mchi_sigma_pattern"));
    if (!hq_local) hq_local = dynamic_cast<TH2D*>(fin->Get("q_mchi_sigma"));
    if (!hq_local) {
      std::cerr << "[limit] ERROR: histogram 'q_mchi_sigma_pattern' or 'q_mchi_sigma' not found in "
                << in_path << "\n";
      fin->Close();
      return 1;
    }
    hq_local->SetDirectory(nullptr);

    TH1D* hul = dynamic_cast<TH1D*>(fin->Get("upper_limit_sigma_e_mchi"));
    if (hul) hul->SetDirectory(nullptr);

    // Preferred over hul when present: a TGraph stores the exact mchi grid,
    // bypassing TH1D bin-center artifacts on log-spaced (uneven) X grids.
    TGraph* g_ul_in = dynamic_cast<TGraph*>(fin->Get("upper_limit_sigma_e_mchi_graph"));
    TGraph* g_ul_owned = nullptr;
    if (g_ul_in) {
      g_ul_owned = dynamic_cast<TGraph*>(g_ul_in->Clone(
          ("g_upper_limit_in_" + std::to_string(idx_file)).c_str()));
    }

    TGraph* g_pydme_in = dynamic_cast<TGraph*>(fin->Get("upper_limit_sigma_e_mchi_pydme_bisection"));
    TGraph* g_pydme_owned = nullptr;
    if (g_pydme_in) {
      g_pydme_owned = dynamic_cast<TGraph*>(g_pydme_in->Clone(
          ("g_pydme_bisection_in_" + std::to_string(idx_file)).c_str()));
    }

    fin->Close();

    const bool have_ul_source = (g_ul_owned != nullptr) || (hul != nullptr);
    const bool use_hul = (have_ul_source && !force_from_qhist);
    const bool build_ul = draw_both ? have_ul_source : use_hul;
    const bool build_pydme_diag = draw_both && (g_pydme_owned != nullptr);
    const bool build_q = draw_both ? !build_pydme_diag : !use_hul;

    bool any_points = false;

    auto legend_with_suffix = [&](const char* suffix) -> std::string {
      return plain_legend ? legend_labels[idx_file]
                          : legend_labels[idx_file] + suffix;
    };

    auto build_and_push_ul = [&]() {
      std::vector<double> mchi_vals_local;
      std::vector<double> sigma_lim_vals_local;
      if (g_ul_owned) {
        const int N = g_ul_owned->GetN();
        for (int i = 0; i < N; ++i) {
          double x = 0.0, y = 0.0;
          g_ul_owned->GetPoint(i, x, y);
          if (y > 0.0) {
            mchi_vals_local.push_back(x);
            sigma_lim_vals_local.push_back(y);
          }
        }
      } else if (hul) {
        const int Nx = hul->GetNbinsX();
        for (int ix = 1; ix <= Nx; ++ix) {
          const double mchi = hul->GetXaxis()->GetBinCenter(ix);
          const double sigma_lim = hul->GetBinContent(ix);
          if (sigma_lim > 0.0) {
            mchi_vals_local.push_back(mchi);
            sigma_lim_vals_local.push_back(sigma_lim);
          }
        }
      }
      if (mchi_vals_local.empty()) return;

      const int Npoints = static_cast<int>(mchi_vals_local.size());
      TGraph* g_ul = new TGraph(Npoints);
      g_ul->SetName(("g_dme_limit_" + std::to_string(idx_file) + "_ul").c_str());
      g_ul->SetTitle("");
      for (int i = 0; i < Npoints; ++i) {
        g_ul->SetPoint(i, mchi_vals_local[i], sigma_lim_vals_local[i]);
      }
      glimits.push_back(g_ul);
      legend_entries_per_graph.push_back(legend_with_suffix(" (UL, q-grid)"));
      mchi_vals.insert(mchi_vals.end(), mchi_vals_local.begin(), mchi_vals_local.end());
      sigma_lim_vals.insert(sigma_lim_vals.end(), sigma_lim_vals_local.begin(), sigma_lim_vals_local.end());
      if (idx_file == 0 && !glimit) glimit = g_ul;
      any_points = true;
    };

    auto build_and_push_pydme = [&]() {
      if (!g_pydme_owned) return;
      std::vector<double> mchi_vals_local;
      std::vector<double> sigma_lim_vals_local;
      const int N = g_pydme_owned->GetN();
      for (int i = 0; i < N; ++i) {
        double x = 0.0, y = 0.0;
        g_pydme_owned->GetPoint(i, x, y);
        if (y > 0.0) {
          mchi_vals_local.push_back(x);
          sigma_lim_vals_local.push_back(y);
        }
      }
      if (mchi_vals_local.empty()) return;

      const int Npoints = static_cast<int>(mchi_vals_local.size());
      TGraph* g_pb = new TGraph(Npoints);
      g_pb->SetName(("g_dme_limit_" + std::to_string(idx_file) + "_pydme").c_str());
      g_pb->SetTitle("");
      for (int i = 0; i < Npoints; ++i)
        g_pb->SetPoint(i, mchi_vals_local[i], sigma_lim_vals_local[i]);
      glimits.push_back(g_pb);
      legend_entries_per_graph.push_back(legend_with_suffix(" (pydme bisection)"));
      mchi_vals.insert(mchi_vals.end(), mchi_vals_local.begin(), mchi_vals_local.end());
      sigma_lim_vals.insert(sigma_lim_vals.end(), sigma_lim_vals_local.begin(), sigma_lim_vals_local.end());
      if (idx_file == 0 && !glimit) glimit = g_pb;
      any_points = true;
    };

    auto build_and_push_q = [&]() {
      std::vector<double> mchi_vals_local;
      std::vector<double> sigma_lim_vals_local;

      const int Nx = hq_local->GetNbinsX();
      const int Ny = hq_local->GetNbinsY();

      for (int ix = 1; ix <= Nx; ++ix) {
        const double mchi = hq_local->GetXaxis()->GetBinCenter(ix);

        // Extract q vs sigma for this column
        std::vector<double> sig(Ny), q(Ny);
        for (int iy = 1; iy <= Ny; ++iy) {
          sig[iy - 1] = hq_local->GetYaxis()->GetBinCenter(iy);
          q[iy - 1] = hq_local->GetBinContent(ix, iy);
        }

        // Find first crossing q(i) < q_thr <= q(i+1)
        double sigma_lim = -1.0;
        for (int iy = 0; iy < Ny - 1; ++iy) {
          const double q1 = q[iy];
          const double q2 = q[iy + 1];

          if (q1 < q_thr && q2 >= q_thr && q2 > q1) {
            const double s1 = sig[iy];
            const double s2 = sig[iy + 1];
            if (s1 <= 0.0 || s2 <= 0.0) break;

            const double logS1 = std::log10(s1);
            const double logS2 = std::log10(s2);
            const double t = (q_thr - q1) / (q2 - q1);
            const double logSlim = logS1 + t * (logS2 - logS1);
            sigma_lim = std::pow(10.0, logSlim);
            break;
          }
        }

        if (sigma_lim > 0.0) {
          mchi_vals_local.push_back(mchi);
          sigma_lim_vals_local.push_back(sigma_lim);
        } else {
          const double q_min = *std::min_element(q.begin(), q.end());
          if (q_min >= q_thr) {
            // excludes all scanned sigma -> use lowest scanned sigma
            mchi_vals_local.push_back(mchi);
            sigma_lim_vals_local.push_back(sig.front());
          }
        }
      }

      if (mchi_vals_local.empty()) return;

      const int Npoints = static_cast<int>(mchi_vals_local.size());
      TGraph* g_q = new TGraph(Npoints);
      g_q->SetName(("g_dme_limit_" + std::to_string(idx_file) + "_q").c_str());
      g_q->SetTitle("");
      for (int i = 0; i < Npoints; ++i) {
        g_q->SetPoint(i, mchi_vals_local[i], sigma_lim_vals_local[i]);
      }
      glimits.push_back(g_q);
      legend_entries_per_graph.push_back(legend_with_suffix(" (q-map)"));
      mchi_vals.insert(mchi_vals.end(), mchi_vals_local.begin(), mchi_vals_local.end());
      sigma_lim_vals.insert(sigma_lim_vals.end(), sigma_lim_vals_local.begin(), sigma_lim_vals_local.end());
      if (idx_file == 0 && !glimit) glimit = g_q;
      any_points = true;
    };

    if (build_ul) {
      build_and_push_ul();
    }

    if (build_pydme_diag) {
      build_and_push_pydme();
      if (g_pydme_owned) {
        delete g_pydme_owned;
        g_pydme_owned = nullptr;
      }
    } else if (build_q) {
      build_and_push_q();
      if (!any_points && have_ul_source && force_from_qhist) {
        std::cerr << "[limit] WARNING: q-map empty or no q=" << q_thr
                  << " crossings in " << in_path
                  << "; falling back to upper_limit_sigma_e_mchi(_graph).\n";
        build_and_push_ul();
      }
      if (static_cast<int>(idx_file) == reference_input_index) {
        reference_graph_index = glimits.size() - 1;
      }
    }

    if (hul) {
      delete hul;
      hul = nullptr;
    }
    if (g_ul_owned) {
      delete g_ul_owned;
      g_ul_owned = nullptr;
    }

    if (!any_points) {
      std::cerr << "[limit] ERROR: no valid limit points were found for " << in_path << "\n";
      continue;
    }

    if (idx_file == 0) hq = hq_local;
  }

  if (!glimit) {
    std::cerr << "[limit] ERROR: no valid curves were built.\n";
    return 1;
  }

  if (!out_csv_cli.empty()) {
    if (!SaveLimitContoursCsv(out_csv_cli, glimits, legend_entries_per_graph,
                              q_thr)) {
      return 1;
    }
  }

  // ------------------------------------------------------------------
  // Optional sensitivity band (--band <band.root>)
  //
  // Schema produced by build/ccdarksens_band:
  //   median_expected_sigma_e_mchi{_asymptotic|_toy_mc}
  //   band_1sigma_low_sigma_e_mchi{_asymptotic|_toy_mc}
  //   band_1sigma_high_sigma_e_mchi{_asymptotic|_toy_mc}
  //   band_2sigma_low_sigma_e_mchi{_asymptotic|_toy_mc}
  //   band_2sigma_high_sigma_e_mchi{_asymptotic|_toy_mc}
  //
  // We load both variants if available; --band-mode selects which to draw
  // (default = "both" if both sets are present; otherwise the present one).
  // ------------------------------------------------------------------
  struct BandSet {
    TGraph* median = nullptr;
    TGraph* l1 = nullptr;
    TGraph* h1 = nullptr;
    TGraph* l2 = nullptr;
    TGraph* h2 = nullptr;
    bool complete() const { return median && l1 && h1 && l2 && h2; }
  };
  BandSet band_asy, band_toy;
  bool draw_band_asy = false, draw_band_toy = false;

  if (!band_path.empty()) {
    TFile* fb = TFile::Open(band_path.c_str(), "READ");
    if (!fb || fb->IsZombie()) {
      std::cerr << "[limit][band] WARNING: cannot open band file: " << band_path
                << "\n";
      if (fb) { fb->Close(); delete fb; }
    } else {
      auto load = [&](const std::string& base, const std::string& suf) -> TGraph* {
        TGraph* g = dynamic_cast<TGraph*>(fb->Get((base + suf).c_str()));
        if (!g) return nullptr;
        return static_cast<TGraph*>(g->Clone((base + suf + "_plot").c_str()));
      };
      auto load_set = [&](const std::string& suf) -> BandSet {
        BandSet b;
        b.median = load("median_expected_sigma_e_mchi", suf);
        b.l1     = load("band_1sigma_low_sigma_e_mchi", suf);
        b.h1     = load("band_1sigma_high_sigma_e_mchi", suf);
        b.l2     = load("band_2sigma_low_sigma_e_mchi", suf);
        b.h2     = load("band_2sigma_high_sigma_e_mchi", suf);
        return b;
      };
      band_asy = load_set("_asymptotic");
      band_toy = load_set("_toy_mc");
      fb->Close();
      delete fb;

      const std::string mode = band_mode_arg.empty() ? std::string("auto")
                                                     : band_mode_arg;
      if (mode == "asymptotic") {
        draw_band_asy = band_asy.complete();
      } else if (mode == "toy_mc") {
        draw_band_toy = band_toy.complete();
      } else if (mode == "both") {
        draw_band_asy = band_asy.complete();
        draw_band_toy = band_toy.complete();
      } else {
        // auto: prefer both if available; otherwise whichever is.
        draw_band_asy = band_asy.complete();
        draw_band_toy = band_toy.complete();
      }
      std::cout << "[limit][band] " << band_path
                << "  asymptotic=" << (draw_band_asy ? "ON" : "off")
                << "  toy_mc=" << (draw_band_toy ? "ON" : "off") << "\n";
      if (!draw_band_asy && !draw_band_toy) {
        std::cerr << "[limit][band] WARNING: neither band variant complete in "
                  << band_path << "; nothing drawn.\n";
      }
    }
  }

  // int c_damic2025  = TColor::GetColor(210,  90,  80);   // soft red
  // int c_sensei  = TColor::GetColor(210,  95,  85);   // soft red
  // int c_supercdms  = TColor::GetColor(210,  95,  85);   // soft red
  // int c_darkside  = TColor::GetColor(210,  95,  85);   // soft red
  // int c_panda4t  = TColor::GetColor(210,  95,  85);   // soft red
  // int c_xenonnt  = TColor::GetColor(210,  95,  85);   // soft red

  int c_damic2025  = TColor::GetColor(150, 175, 140);;   // soft red
// 
  int c_srdm       = TColor::GetColor(100, 100, 150);   // bluish gray
  int c_sensei     = TColor::GetColor( 40,  80, 200);   // deep blue
  int c_supercdms  = TColor::GetColor(160,  70, 180);   // purple
  int c_darkside   = TColor::GetColor( 70, 160,  70);   // green
  int c_panda4t    = TColor::GetColor(210, 140,  60);   // orange
  int c_xenonnt    = TColor::GetColor(100, 100, 150);   // bluish gray

  // int c_freezein   = TColor::GetColor(130, 160, 210);   // light blue

  // int c_freezein = TColor::GetColor(110,80,130);
  // int c_freezein = TColor::GetColor(70,90,140);
  // int c_freezein = TColor::GetColor(30,110,105);
  // int c_freezein = TColor::GetColor(90, 60, 130);
  // int c_freezein = TColor::GetColor(110, 70, 160);
  // int c_freezein = TColor::GetColor(120, 90, 150);
  // int c_freezein = TColor::GetColor(150, 40, 50);
  // int c_freezein = TColor::GetColor(170, 55, 70);
  int c_freezein = TColor::GetColor(180, 40, 55);


  // ------------------------------------------------------------------
  // 4) Load external exclusions / theory curves
  // ------------------------------------------------------------------
  TGraph* g_damic_2025      = nullptr;
  TGraph* g_sensei          = nullptr;
  TGraph* g_supercdms       = nullptr;
  TGraph* g_darkside50      = nullptr;
  TGraph* g_panda4T         = nullptr;
  TGraph* g_xenonnT         = nullptr;
  TGraph* g_xenonnT_hi        = nullptr;
  TGraphAsymmErrors* g_xenonnT_bracket = nullptr;
  TGraph* g_damicm_mike     = nullptr;
  TGraph* g_model           = nullptr;
  TGraph* g_solar_reflected = nullptr;
  TGraph* g_srdm            = nullptr;
  TGraph* g_wimp_nominal    = nullptr;
  TGraph* g_wimp_lower      = nullptr;
  TGraph* g_wimp_upper      = nullptr;
  TGraph* g_wimp_damicm_proj = nullptr;
  TGraph* g_wimp_snolab_uncertain = nullptr;
  TGraph* g_wimp_radomir = nullptr;

  if (migdal_mode) {
    // --- Migdal effect: σ_n (DM-nucleon) literature ---
    // DAMIC-M (2025 LBC) is a heavy-mediator-only digitization -- do not show
    // it on light-mediator plots, since we have no light-mediator LBC
    // reference curve to display in its place (data/previous_limits/Migdal_light/
    // is empty). The four external experiment curves below are kept for both
    // mediators as a rough visual reference.
    const std::string migdal_base = "data/previous_limits/Migdal_heavy/";
    if (mediator == "heavy") {
      g_damic_2025 = LoadTxtToTGraph(
          migdal_base + "DAMIC-M_2025_DMn-Migdal_heavymediator.txt", kBlack, 1, 3);
    }
    g_xenonnT    = LoadMassLimitDAT(
        migdal_base + "XENON1T-Migdal.dat",    1e-39, c_xenonnt,  1, 2);
    g_panda4T    = LoadExclusionCSV(
        migdal_base + "Pandax_Migdal_2023.csv", c_panda4t, 1, 2, 1000.0, 1.0);
    g_darkside50 = LoadMassLimitDAT(
        migdal_base + "Darkside-Migdal.dat",    1e-39, c_darkside, 1, 2);
    g_sensei     = LoadMassLimitDAT(
        migdal_base + "SENSEI-Migdal.dat",      1e-35, c_sensei,   1, 2);

    int n_loaded = 0;
    if (g_damic_2025) ++n_loaded;
    if (g_xenonnT)    ++n_loaded;
    if (g_panda4T)    ++n_loaded;
    if (g_darkside50) ++n_loaded;
    if (g_sensei)     ++n_loaded;
    std::cout << "[limit][migdal] Loaded " << n_loaded
              << " literature curve(s) from " << migdal_base << "\n";
  } else if (wimp_mode) {
    // --- WIMP-nucleon SI elastic literature: DAMIC at SNOLAB, 0.6 kg-day
    // (PhysRevD.94.082006 / arXiv:1607.07410). See docs/ClusterFitMC_Design.md
    // Sec.6 for the full provenance/labeling history. Primary trusted curve is
    // observed_fig11_diego -- a direct digitization of the paper's own Fig. 11
    // (the single "This work, 0.6 kg d" curve, unambiguous -- no separate
    // observed/expected lines to confuse in that figure). radomir's curve and
    // the nominal/lower/upper triplet do NOT match it well at low mass and are
    // kept only as secondary, less-certain references. CSV columns are
    // (mchi_GeV, sigma_n_cm2); scale_x=1000 converts to this app's MeV axis.
    const std::string wimp_base = "data/previous_limits/WIMP/";
    g_damic_2025 = LoadExclusionCSV(
        wimp_base + "DAMIC_2016_SNOLAB_observed_fig11_diego.csv", kRed + 1, 1, 3, 1000.0, 1.0);
    // nominal: still WIMPyCCD-sourced (uncertain provenance) -- kept for
    // context only, not independently confirmed. lower/upper: directly
    // digitized by the user from the actual Fig. 11 PDF (the paper's own
    // "expected sensitivity +-1 sigma" band around the red observed line,
    // per its caption) -- same trusted-source standard as g_damic_2025,
    // superseding the earlier WIMPyCCD-sourced lower/upper files.
    g_wimp_nominal = LoadExclusionCSV(
        wimp_base + "DAMIC_2016_SNOLAB_expected_nominal.csv", kGray + 2, 3, 2, 1000.0, 1.0);
    g_wimp_lower = LoadExclusionCSV(
        wimp_base + "DAMIC_2016_SNOLAB_expected_lower_diego.csv", kRed + 1, 2, 3, 1000.0, 1.0);
    g_wimp_upper = LoadExclusionCSV(
        wimp_base + "DAMIC_2016_SNOLAB_expected_upper_diego.csv", kRed + 1, 2, 3, 1000.0, 1.0);
    // Context only: DAMIC-M projected sensitivity (1 kg-year) -- a DIFFERENT,
    // newer/larger experiment, not part of the 2016-reproduction comparison itself.
    g_wimp_damicm_proj = LoadExclusionCSV(
        wimp_base + "DAMIC-M_projected_1kgyear.csv", kAzure + 1, 2, 2, 1000.0, 1.0);
    // Same WIMPyCCD digitized/ folder, but UNCERTAIN provenance (sits ~10x below even
    // the expected-band lower edge, implausible for a real result) -- included on request,
    // drawn distinctly (magenta dash-dot) so it is never mistaken for the trusted observed curve.
    g_wimp_snolab_uncertain = LoadExclusionCSV(
        wimp_base + "DAMIC_2016_SNOLAB_damic_at_snolab_uncertain.csv", kMagenta + 1, 9, 2, 1000.0, 1.0);
    // Demoted to secondary reference: does not match observed_fig11_diego well at low
    // mass (~20x off at 1.3-2 GeV) despite tracking it more closely at high mass --
    // drawn thin/muted so it is never mistaken for the primary trusted curve.
    g_wimp_radomir = LoadExclusionCSV(
        wimp_base + "DAMIC_2016_SNOLAB_observed_radomir.csv", kTeal + 2, 2, 1, 1000.0, 1.0);

    int n_loaded = 0;
    if (g_damic_2025)            ++n_loaded;
    if (g_wimp_nominal)          ++n_loaded;
    if (g_wimp_lower)            ++n_loaded;
    if (g_wimp_upper)            ++n_loaded;
    if (g_wimp_damicm_proj)      ++n_loaded;
    if (g_wimp_snolab_uncertain) ++n_loaded;
    if (g_wimp_radomir)          ++n_loaded;
    std::cout << "[limit][wimp] Loaded " << n_loaded
              << " literature curve(s) from " << wimp_base << "\n";
  } else if (!dark_photon) {
    // --- DM-electron literature curves ---
    if (mediator == "intermediate") {
      // No published literature curves for intermediate mediator — skip all loading.
      std::cout << "[limit] intermediate mediator: no literature curves loaded.\n";
    } else {
    const std::string limits_base = "data/previous_limits/" + mediator + "_mediator/";
    const std::string damic_paper_export_heavy =
        "collab_frameworks/pydme/analysis/DailyModulation/LBC-Sep2024/paper_figures/data/"
        "LBC_results/ScienceRun2024_results-Pattern/DAMIC-M_2025_QEDark_DMe_heavymediator.txt";
    const std::string damic_paper_export_light =
        "collab_frameworks/pydme/analysis/DailyModulation/LBC-Sep2024/paper_figures/data/"
        "LBC_results/ScienceRun2024_results-Pattern/DAMIC-M_2025_QEDark_DMe_ulightmediator.txt";
    const std::string damic_file = (mediator == "heavy")
      ? "DAMIC-M_2025_QEDark_DMe_heavymediator.txt"
      : "DAMIC-M_2025_QEDark_DMe_ulightmediator.txt";
    const std::string freeze_file = (mediator == "heavy")
      ? "Freeze-Out_F1_CMS-community.csv"
      : "Freeze_in_limit.csv";
    const std::string srdm_file = (mediator == "heavy")
      ? "SRDM_limit_heavy.dat"
      : "SRDM_limit_light.dat";

    {
      const std::string damic_path =
          (mediator == "heavy") ? damic_paper_export_heavy
          : (mediator == "light") ? damic_paper_export_light
                                  : (limits_base + damic_file);
      g_damic_2025 = LoadTxtToTGraph(damic_path, kBlack, 1, 3);
      if (!g_damic_2025) {
        std::cerr << "[limit] WARNING: cannot load paper-export reference from "
                  << damic_path << "\n";
      }
    }
    g_sensei       = LoadExclusionCSV(limits_base + "SENSEI.csv", c_xenonnt, 1, 2);
    g_supercdms    = LoadExclusionCSV(limits_base + "SuperCDMS.csv", c_xenonnt, 1, 2);
    g_darkside50   = LoadExclusionCSV(limits_base + "DarkSide50.csv", c_xenonnt, 1, 2);
    g_panda4T      = LoadExclusionCSV(limits_base + "Panda4T.csv", c_xenonnt, 1, 2);
    g_xenonnT      = LoadExclusionCSV(limits_base + "Xenon.csv", c_xenonnt, 1, 2);
    g_damicm_mike  = LoadExclusionCSV(limits_base + "damic-m_1kgyear_mike.csv",
                                      kCyan, (mediator == "heavy") ? 2 : 1, 3);
    g_model        = LoadExclusionCSV(limits_base + freeze_file, c_freezein, 1, 3);
    g_solar_reflected = LoadMassLimitDAT(limits_base + srdm_file, 1e-38, kOrange+2, 1, 2);
    g_srdm         = LoadExclusionCSV(
        "data/previous_limits/srdm/carlos_srdm_ulm_limit.csv", c_freezein, 1, 3);
    } // end heavy/light literature block
  } else {
    // --- Dark photon absorption literature (ScienceRun2024 notebook HP figure) ---
    const std::string dp_base = "data/previous_limits/dark_photon/";
    constexpr double kKeVtoEv = 1000.0;
    // Paolo color scheme; direct-detection literature fills share one gray
    const Color_t kHpXenon      = kViolet - 5;
    const Color_t kHpSuperCDMS  = kGray + 1;

    g_damic_2025 = LoadDarkPhotonLimitCSV(
        dp_base + "DAMIC-M_2025_DAMICmodel_HP.txt", 0, 1, 1.0, kBlack, 1, 3);
    g_panda4T = LoadDarkPhotonLimitCSV(
        dp_base + "SCDMS_HP2024.txt", 0, 1, 1.0, kHpSuperCDMS, 7, 2);
    const DarkPhotonXenonNtBracket xnt = LoadDarkPhotonXenonNtHPBracket(
        dp_base + "dark_photon_xenonnt_hp_bracket.csv", kKeVtoEv, kHpXenon, 7, 3);
    g_xenonnT = xnt.lo;
    g_xenonnT_hi = xnt.hi;
    g_xenonnT_bracket = xnt.band;
    g_darkside50 = LoadDarkPhotonLimitCSV(
        dp_base + "dark_photon_stellar_limits_combined.csv", 1, 2, kKeVtoEv,
        kGreen + 2, 3, 2);

    int n_loaded = 0;
    if (g_damic_2025) ++n_loaded;
    if (g_panda4T) ++n_loaded;
    if (g_xenonnT) ++n_loaded;
    if (g_darkside50) ++n_loaded;
    if (g_model) ++n_loaded;
    std::cout << "[limit][dark-photon] Loaded " << n_loaded
              << " literature curve(s) from " << dp_base
              << " (notebook HP set + stellar limits)\n";
  }



  // Expand ranges using all non-null graphs
  double xmin = *std::min_element(mchi_vals.begin(), mchi_vals.end());
  double xmax = *std::max_element(mchi_vals.begin(), mchi_vals.end());
  double ymin = *std::min_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
  double ymax = *std::max_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
  if (migdal_mode) {
    if (g_damic_2025) UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_xenonnT);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_panda4T);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_darkside50);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_sensei);
  } else if (wimp_mode) {
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025);
    // Expected +-1sigma band (now directly digitized, trusted) stays visible
    // even in --wimp-minimal; only the other secondary/context curves drop.
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_wimp_lower);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_wimp_upper);
    if (!wimp_minimal) {
      UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_wimp_nominal);
      UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_wimp_snolab_uncertain);
      UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_wimp_radomir);
    }
  } else if (!dark_photon) {
    if (g_damic_2025) UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_sensei);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_supercdms);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_darkside50);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_panda4T);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_xenonnT);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_model);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_solar_reflected);
  }

  // For dark photon: freeze xmax from scan data before literature curves can shrink it.
  // Literature curves (LBC, stellar) may not extend as far as our scan grid.
  const double xmax_data = xmax;

  // y-bounds of the plot (must match glimit->SetMinimum/Maximum)
  // For dark photon, only update y from g_srdm, not x (literature curve may end at ~20 eV).
  if (dark_photon) {
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_srdm);
    xmax = xmax_data;  // restore data-driven xmax
  } else if (!migdal_mode) {
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_srdm);
  }

  // Expand x/y ranges to include any drawn band envelopes (use the 2σ outer
  // edges so the full Brazilian-band ribbon is visible).
  if (draw_band_asy) {
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, band_asy.l2);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, band_asy.h2);
  }
  if (draw_band_toy) {
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, band_toy.l2);
    UpdateRangesFromGraph(xmin, xmax, ymin, ymax, band_toy.h2);
  }

  double y_bottom = ymin * 0.3;
  double y_top    = ymax * 3.0;
  if (dark_photon) {
    // Stellar rows reach ~1e-7 eV; ensure literature fills cover the full axis top.
    constexpr double kDpLiteratureFillTop = 1e-6;
    y_top = std::max(y_top, kDpLiteratureFillTop);
  } else if (migdal_mode) {
    // σ_n limits: SENSEI/XENON1T can be as high as ~1e-25 at low mass — cover top.
    y_top = std::max(y_top, 1e-24);
  }

  // CLI axis overrides: apply after auto-range computation
  if (axis_x_min > 0) xmin     = axis_x_min;
  if (axis_x_max > 0) xmax     = axis_x_max;
  if (axis_y_min > 0) y_bottom = axis_y_min;
  if (axis_y_max > 0) y_top    = axis_y_max;

  // ---- Experimental contours: fill ABOVE the curve ----
  // TGraph* g_damic_2025_fill = MakeFilledBandAbove(
  //     g_damic_2025, y_top,
  //     g_damic_2025 ? g_damic_2025->GetLineColor() : c_damic2025, 0.25);

  TGraph* g_damic_2025_fill = nullptr;
  TGraph* g_sensei_fill = nullptr;
  TGraph* g_supercdms_fill = nullptr;
  TGraph* g_darkside50_fill = nullptr;
  TGraph* g_panda4T_fill = nullptr;
  TGraph* g_xenonnT_fill = nullptr;
  TGraph* g_direct_det_env = nullptr;
  TGraph* g_direct_det_fill = nullptr;
  TGraph* g_model_fill = nullptr;
  TGraph* g_migdal_env = nullptr;
  TGraph* g_migdal_env_fill = nullptr;
  if (migdal_mode) {
    // Build lower envelope of all available Migdal limits for the shaded exclusion band.
    std::vector<const TGraph*> migdal_curves;
    if (g_damic_2025)  migdal_curves.push_back(g_damic_2025);
    if (g_xenonnT)     migdal_curves.push_back(g_xenonnT);
    if (g_panda4T)     migdal_curves.push_back(g_panda4T);
    if (g_darkside50)  migdal_curves.push_back(g_darkside50);
    if (g_sensei)      migdal_curves.push_back(g_sensei);
    if (!migdal_curves.empty()) {
      g_migdal_env = MakeLowerEnvelope(migdal_curves, xmin, xmax_data);
      g_migdal_env_fill = MakeFilledBandAbove(g_migdal_env, y_top, kGray + 1, 0.25);
    }
  } else if (dark_photon) {
    std::vector<const TGraph*> direct_det_curves;
    if (g_damic_2025) direct_det_curves.push_back(g_damic_2025);
    if (g_panda4T) direct_det_curves.push_back(g_panda4T);
    if (g_xenonnT) direct_det_curves.push_back(g_xenonnT);
    if (!direct_det_curves.empty()) {
      const double kDpEnvMassMin = 0.01;
      const double kDpEnvMassMax = xmax_data;  // driven by scan data range
      g_direct_det_env = MakeLowerEnvelope(direct_det_curves, kDpEnvMassMin, kDpEnvMassMax);
      g_direct_det_fill = MakeFilledBandAbove(g_direct_det_env, y_top, kGray + 1, 0.25);
    }
    g_darkside50_fill = MakeFilledBandAbove(
        g_darkside50, y_top, kGreen + 2, 0.30);
  } else {
    if (g_damic_2025) {
      g_damic_2025_fill = MakeFilledBandAbove(g_damic_2025, y_top, c_damic2025, 0.25);
    }
    g_sensei_fill = MakeFilledBandAbove(
        g_sensei, y_top,
        g_sensei ? g_sensei->GetLineColor() : c_sensei, 0.25);
    g_supercdms_fill = MakeFilledBandAbove(
        g_supercdms, y_top,
        g_supercdms ? g_supercdms->GetLineColor() : c_supercdms, 0.25);
    g_darkside50_fill = MakeFilledBandAbove(
        g_darkside50, y_top,
        g_darkside50 ? g_darkside50->GetLineColor() : c_darkside, 0.25);
    g_panda4T_fill = MakeFilledBandAbove(
        g_panda4T, y_top,
        g_panda4T ? g_panda4T->GetLineColor() : c_panda4t, 0.25);
    g_xenonnT_fill = MakeFilledBandAbove(
        g_xenonnT, y_top,
        g_xenonnT ? g_xenonnT->GetLineColor() : c_xenonnt, 0.25);
  }

  TGraph* g_solar_reflected_fill = MakeFilledBandAbove(
      g_solar_reflected, y_top,
      g_solar_reflected ? g_solar_reflected->GetLineColor() : kOrange+2, 0.25);

  if (!dark_photon && !migdal_mode && !wimp_mode) {
    // Heavy: red fill below freeze-out target, no solid line
    // Light: solid red line only
    if (mediator == "heavy")
      g_model_fill = MakeFilledBand(g_model, y_bottom, kRed - 4, 0.35);
  }

  // (Optional) if I ever draw g_damicm_mike and want it shaded too:
  // TGraph* g_damicm_mike_fill = MakeFilledBand(
  //     g_damicm_mike, y_bottom, g_damicm_mike->GetLineColor(), 0.25);



  // ------------------------------------------------------------------
  // 5) Make a quick publication-style plot
  // ------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.025);
  gStyle->SetLabelSize(0.035,"xyz");
  gStyle->SetTitleSize(0.035,"xyz");
  gStyle->SetTitleOffset(1.2,"y");
  gStyle->SetTitleOffset(1.1,"x");
  gStyle->SetTickLength(0.02,"x");
  gStyle->SetTickLength(0.02,"y");


  const std::string canvas_title =
      plot_title.empty()
          ? (dark_photon  ? std::string("Dark photon limit")
             : migdal_mode ? std::string("Migdal limit")
             : wimp_mode   ? std::string("WIMP-nucleon SI limit")
                           : std::string("DM-e limit"))
          : plot_title;
  TCanvas* c = new TCanvas("c_limit", canvas_title.c_str(), 900, 700);
  c->SetLogx();
  c->SetLogy();
  c->SetTicks(1,1);
  c->SetEditable(kTRUE);
  c->SetLeftMargin(0.14);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(plot_title.empty() ? 0.06 : 0.10);

  const std::vector<Color_t> curve_colors = {
      kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 1, kCyan + 1};
  // Si reference (1.2 eV / 3.8 eV): gold — distinct from freeze-in/out (maroon) and dashed sweeps.
  const Color_t kSiReferenceColor = static_cast<Color_t>(TColor::GetColor(218, 165, 32));
  const std::vector<Style_t> curve_styles = {2, 1, 3, 4, 5, 6, 7};
  for (size_t i = 0; i < glimits.size(); ++i) {
    if (!glimits[i]) continue;
    if (smooth_curves) {
      TGraphSmooth gs("gs");
      TGraph* tmp = gs.SmoothLowess(glimits[i], "", 0.15);
      TGraph* smoothed = new TGraph(*tmp);
      smoothed->SetName(glimits[i]->GetName());
      glimits[i] = smoothed;
    }
    if (i == reference_graph_index || (wimp_mode && glimits.size() == 1)) {
      // wimp_mode with a single result curve: black dashed against the
      // paper's now-red observed/band curves, per the user's requested
      // color scheme (paper = red, our result = black).
      glimits[i]->SetLineColor(kBlack);
      glimits[i]->SetLineStyle(kDashed);
      glimits[i]->SetLineWidth(3);
    } else if (wimp_mode && glimits.size() == 2) {
      // wimp_mode with exactly two result curves (e.g. two DAMIC-M exposure
      // scenarios): black + blue, per explicit request.
      glimits[i]->SetLineColor(i == 0 ? kBlack : kBlue + 1);
      glimits[i]->SetLineStyle(kDashed);
      glimits[i]->SetLineWidth(3);
    } else if (wimp_mode && glimits.size() == 3) {
      // wimp_mode with exactly three result curves (e.g. two DAMIC-M
      // exposure scenarios plus the DAMIC 2016 reproduction): black, blue,
      // green, per explicit request.
      static const std::vector<Color_t> three_colors = {kBlack, kBlue + 1, kGreen + 2};
      glimits[i]->SetLineColor(three_colors[i]);
      glimits[i]->SetLineStyle(kDashed);
      glimits[i]->SetLineWidth(3);
    } else {
      glimits[i]->SetLineWidth(3);
      glimits[i]->SetLineColor(curve_colors[i % curve_colors.size()]);
      glimits[i]->SetLineStyle(kDashed);
    }
  }

  // ------------------------------------------------------------------
  // 4) Axis ranges from BOTH our limit curve AND external exclusions
  // ------------------------------------------------------------------
//   double xmin = *std::min_element(mchi_vals.begin(), mchi_vals.end());
//   double xmax = *std::max_element(mchi_vals.begin(), mchi_vals.end());
//   double ymin = *std::min_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
//   double ymax = *std::max_element(sigma_lim_vals.begin(), sigma_lim_vals.end());

  // Load external exclusion(s) from CSV
  // IMPORTANT: use the correct path where the CSV actually lives
  // e.g. "./SRDM_XENON1T-s20_ulightmediator.csv"

//   glimit->GetXaxis()->SetLimits(xmin * 0.8, xmax );
//   glimit->SetMinimum(ymin * 0.3);
//   glimit->SetMaximum(ymax * 3.0);

//   // Draw our limit first (defines axes)
//   glimit->Draw("AL");

//   // Overlay external exclusion(s)
//   // Overlay external exclusions
//   if (g_damic_2025)  g_damic_2025->Draw("L SAME");
//   if (g_sensei)      g_sensei->Draw("L SAME");
//   if (g_supercdms)   g_supercdms->Draw("L SAME");
//   if (g_darkside50)  g_darkside50->Draw("L SAME");
//   if (g_panda4T)     g_panda4T->Draw("L SAME");
//   // if (g_damicm_mike) g_damicm_mike->Draw("L SAME");
//   // Overlay theory curves
// //   if (g_freeze_in)   g_freeze_in->Draw("L SAME");
//   if (g_model)  g_model->Draw("L SAME");

// Pre-computed lower envelope of all direct-detection limits (incl. DAMIC-M 2025).
// Generated by utils/build_dme_direct_detection_envelope.py.
// The DAMIC-M 2025 solid line is drawn separately on top.
TGraph* g_env_electron      = nullptr;
TGraph* g_env_fill_electron = nullptr;
if (!dark_photon && !migdal_mode && !wimp_mode && mediator != "intermediate") {
  const std::string env_csv = "data/previous_limits/" + mediator + "_mediator/direct_detection_envelope.csv";
  g_env_electron = LoadExclusionCSV(env_csv, kGray + 1, 1, 0);
  if (g_env_electron) {
    g_env_fill_electron = MakeFilledBandAbove(g_env_electron, y_top, kGray + 1, 0.25);
  }
}


  const char* kXtitle = dark_photon ? "m_{A'} [eV]"     : "m_{#chi} [MeV/c^{2}]";
  const char* kYtitle = dark_photon ? "#varepsilon"
                      : (migdal_mode || wimp_mode) ? "#bar{#sigma}_{n} [cm^{2}]"
                                    : "#bar{#sigma}_{e} [cm^{2}]";
  // x-axis: for dark photon/wimp use data-driven range; for DM-electron/migdal use fixed
  // publication window. CLI overrides (--x-min/--x-max/--y-min/--y-max) take precedence.
  // wimp mode needs data-driven range -- the WIMP-nucleon mass grid (0.5-10 GeV) would be
  // mostly clipped by the DM-electron/migdal fixed window (0.1-1000 MeV).
  const double kMassPlotMin = (axis_x_min > 0) ? axis_x_min : ((dark_photon || wimp_mode) ? xmin  : 1e-1);
  const double kMassPlotMax = (axis_x_max > 0) ? axis_x_max : ((dark_photon || wimp_mode) ? xmax  : 1e3);
  const double y_plot_max =
      (axis_y_max > 0) ? axis_y_max : (dark_photon ? std::max(ymax * 3.0, 1e-6) : ymax * 3.0);

  // View range for interactive GUI: scan limits only (not literature to 10^11 MeV).
  const double x_view_lo = kMassPlotMin;
  const double x_view_hi = kMassPlotMax;
  const double y_view_lo = (axis_y_min > 0) ? axis_y_min : y_bottom;
  const double y_view_hi = y_plot_max;

  glimit->GetXaxis()->SetTitle(kXtitle);
  glimit->GetYaxis()->SetTitle(kYtitle);
  glimit->GetXaxis()->SetTitleSize(0.05);
  glimit->GetYaxis()->SetTitleSize(0.05);
  glimit->GetXaxis()->SetLabelSize(0.04);
  glimit->GetYaxis()->SetLabelSize(0.04);
  if (batch) {
    glimit->GetXaxis()->SetLimits(kMassPlotMin, kMassPlotMax);
    glimit->SetMinimum(y_bottom);
    glimit->SetMaximum(y_plot_max);
  }
  glimit->Draw("AL");

  // ---- Brazilian sensitivity band(s) (drawn first so curves overlay them) ----
  // 2σ uses HEP-canonical "orange"; 1σ uses "green". When both variants are
  // drawn, asymptotic ribbons are drawn first at low alpha; toy_mc fills on top.
  TGraph* band_fill_asy_2s = nullptr;
  TGraph* band_fill_asy_1s = nullptr;
  TGraph* band_fill_toy_2s = nullptr;
  TGraph* band_fill_toy_1s = nullptr;
  TGraph* band_toy_median_vis = nullptr;
  TGraph* band_asy_median_vis = nullptr;

  auto draw_band_ribbons = [&](const TGraph* lo2, const TGraph* hi2,
                               const TGraph* lo1, const TGraph* hi1,
                               int color2, double alpha2, int color1,
                               double alpha1) -> std::pair<TGraph*, TGraph*> {
      TGraph* leg2 = nullptr;
      TGraph* leg1 = nullptr;
      for (TGraph* g : MakeBandRibbonPieces(lo2, hi2)) {
          g->SetFillStyle(1001);
          g->SetFillColorAlpha(color2, alpha2);
          g->SetLineWidth(0);
          g->Draw("F SAME");
          if (!leg2) leg2 = g;
      }
      for (TGraph* g : MakeBandRibbonPieces(lo1, hi1)) {
          g->SetFillStyle(1001);
          g->SetFillColorAlpha(color1, alpha1);
          g->SetLineWidth(0);
          g->Draw("F SAME");
          if (!leg1) leg1 = g;
      }
      return {leg2, leg1};
  };

  // Draw asymptotic ribbons under toy MC when both are enabled so the
  // low-alpha asymptotic wash does not sit on top of the toy fills.
  if (draw_band_asy && draw_band_toy) {
    const double a2 = 0.20;
    const double a1 = 0.20;
    const auto pr_asy =
        draw_band_ribbons(band_asy.l2, band_asy.h2, band_asy.l1, band_asy.h1,
                          kOrange, a2, kGreen + 1, a1);
    band_fill_asy_2s = pr_asy.first;
    band_fill_asy_1s = pr_asy.second;
    const auto pr_toy =
        draw_band_ribbons(band_toy.l2, band_toy.h2, band_toy.l1, band_toy.h1,
                          kOrange, 0.55, kGreen + 1, 0.55);
    band_fill_toy_2s = pr_toy.first;
    band_fill_toy_1s = pr_toy.second;
  } else {
    if (draw_band_toy) {
      const auto pr_toy =
          draw_band_ribbons(band_toy.l2, band_toy.h2, band_toy.l1, band_toy.h1,
                            kOrange, 0.55, kGreen + 1, 0.55);
      band_fill_toy_2s = pr_toy.first;
      band_fill_toy_1s = pr_toy.second;
    }
    if (draw_band_asy) {
      const auto pr_asy =
          draw_band_ribbons(band_asy.l2, band_asy.h2, band_asy.l1, band_asy.h1,
                            kOrange, 0.55, kGreen + 1, 0.55);
      band_fill_asy_2s = pr_asy.first;
      band_fill_asy_1s = pr_asy.second;
    }
  }

  // Median lines are drawn after limit curves (see below) so they are not
  // hidden by filled regions, SRDM (also black), or scan lines.

  // ---- DC-sweep band (dark photon: hatched fill spanning lo..hi DC, central line only) ----
  TGraph* g_dc_band = nullptr;
  const bool draw_dc_band = dark_photon && !dc_band_lo_path.empty() && !dc_band_hi_path.empty();
  if (draw_dc_band) {
    TGraph* g_dc_lo  = ExtractLimitCurveFromRoot(dc_band_lo_path,  q_thr);
    TGraph* g_dc_hi  = ExtractLimitCurveFromRoot(dc_band_hi_path,  q_thr);
    TGraph* g_dc_cen = glimits.empty() ? nullptr : glimits[0];
    if (g_dc_lo && g_dc_hi && g_dc_cen) {
        g_dc_band = MakeDcBandPolygon(g_dc_lo, g_dc_hi, g_dc_cen, y_bottom, y_plot_max);
        g_dc_band->SetFillColor(kRed - 4);
        g_dc_band->SetFillStyle(3001);
        g_dc_band->SetLineColor(0);
    } else {
        std::cerr << "[dc-band] WARNING: failed to build band polygon (check --dc-band-lo/hi paths)\n";
    }
  }

  // ---- Filled regions (behind lines) ----
  if (migdal_mode) {
    if (g_migdal_env_fill) g_migdal_env_fill->Draw("F SAME");
  } else if (wimp_mode) {
    // No filled envelope for WIMP mode -- just the four literature lines below.
  } else if (dark_photon) {
    if (g_darkside50_fill)  g_darkside50_fill->Draw("F SAME");
    if (g_direct_det_fill)  g_direct_det_fill->Draw("F SAME");
    if (g_dc_band)          g_dc_band->Draw("F SAME");
  } else {
    if (g_env_fill_electron)    g_env_fill_electron->Draw("F SAME");
    if (g_solar_reflected_fill) g_solar_reflected_fill->Draw("F SAME");
    if (g_model_fill)           g_model_fill->Draw("F SAME");
  }

  // ---- Line contours on top ----
  if (migdal_mode) {
    if (g_xenonnT)    g_xenonnT->Draw("L SAME");
    if (g_panda4T)    g_panda4T->Draw("L SAME");
    if (g_darkside50) g_darkside50->Draw("L SAME");
    if (g_sensei)     g_sensei->Draw("L SAME");
    if (g_damic_2025) g_damic_2025->Draw("L SAME");
  } else if (wimp_mode) {
    if (!wimp_minimal) {
      if (g_wimp_damicm_proj)      g_wimp_damicm_proj->Draw("L SAME");      // context only, different experiment
      if (g_wimp_snolab_uncertain) g_wimp_snolab_uncertain->Draw("L SAME"); // uncertain provenance, see comment above
      if (g_wimp_radomir)          g_wimp_radomir->Draw("L SAME");          // demoted secondary reference
      if (g_wimp_nominal) g_wimp_nominal->Draw("L SAME");
    }
    // Expected +-1sigma band: directly digitized, trusted -- stays visible
    // even in --wimp-minimal.
    if (g_wimp_lower)   g_wimp_lower->Draw("L SAME");
    if (g_wimp_upper)   g_wimp_upper->Draw("L SAME");
    if (g_damic_2025)   g_damic_2025->Draw("L SAME");  // observed (Fig. 11, directly verified), drawn last/on top
  } else if (dark_photon) {
    if (g_damic_2025)  g_damic_2025->Draw("L SAME");
  } else {
    if (show_damic && g_damic_2025) g_damic_2025->Draw("L SAME");
    if (g_damic_2025)  g_damic_2025->Draw("L SAME");
    if (mediator != "heavy" && g_model) g_model->Draw("L SAME");
  }

  // Limit curves: scan sweeps first, Si reference on top (solid red).
  // In DC-band mode (dark photon): only draw the central curve (glimits[0]) as dashed red.
  if (draw_dc_band) {
    if (!glimits.empty() && glimits[0]) {
      glimits[0]->SetLineColor(kRed + 1);
      glimits[0]->SetLineStyle(2);
      glimits[0]->SetLineWidth(3);
      glimits[0]->Draw("L SAME");
    }
    // Still draw the reference scan (e.g. DAMIC-M 1 kg-yr Si projection)
    if (reference_graph_index < glimits.size() && glimits[reference_graph_index]) {
      glimits[reference_graph_index]->Draw("L SAME");
    }
  } else {
    for (size_t i = 0; i < glimits.size(); ++i) {
      if (!glimits[i] || i == reference_graph_index) continue;
      glimits[i]->Draw("L SAME");
    }
    if (reference_graph_index < glimits.size() && glimits[reference_graph_index]) {
      glimits[reference_graph_index]->Draw("L SAME");
    }
  }

  // Band medians on top: log-safe copy skips y=0 knots. Toy MC median uses a
  // color distinct from scan UL curves (often black dashed Asimov) which sit
  // just under this draw and can match the median in (x,y) — same black reads
  // as one line even when z-ordered on top.
  if (draw_band_toy && band_toy.median) {
    band_toy_median_vis = FilterGraphPositiveLogSafe(band_toy.median,
                                                     "g_median_toy_mc_logsafe");
    TGraph* gmed = band_toy_median_vis ? band_toy_median_vis : band_toy.median;
    gmed->SetLineColor(kBlue + 2);
    gmed->SetLineStyle(kDashed);
    gmed->SetLineWidth(2);
    gmed->Draw("L SAME");
  }
  if (draw_band_asy && band_asy.median) {
    band_asy_median_vis = FilterGraphPositiveLogSafe(band_asy.median,
                                                     "g_median_asymptotic_logsafe");
    TGraph* gmed = band_asy_median_vis ? band_asy_median_vis : band_asy.median;
    gmed->SetLineColor(draw_band_toy ? kGray + 2 : kBlack);
    gmed->SetLineStyle(draw_band_toy ? kDotted : kDashed);
    gmed->SetLineWidth(2);
    gmed->Draw("L SAME");
  }

  const bool legend_compact = legend_left || legend_right;
  TLegend* leg = legend_left
                     ? new TLegend(0.11, 0.50, 0.40, 0.90)
                     : legend_right
                           ? new TLegend(0.52, 0.43, 0.89, 0.88)
                           : new TLegend(0.55, 0.45, 0.82, 0.82);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(legend_compact ? 0.024 : 0.03);
  if (legend_compact) {
    leg->SetEntrySeparation(0.35);
  }

  // Sensitivity-band legend entries (only when bands are actually drawn).
  if (draw_band_toy) {
    if (band_toy.median) {
      TGraph* le_med = band_toy_median_vis ? band_toy_median_vis : band_toy.median;
      leg->AddEntry(le_med,
                    draw_band_asy ? "Median expected (toy MC)"
                                  : "Median expected",
                    "l");
    }
    if (band_fill_toy_1s)
      leg->AddEntry(band_fill_toy_1s,
                    draw_band_asy ? "#pm 1#sigma (toy MC)" : "#pm 1#sigma",
                    "f");
    if (band_fill_toy_2s)
      leg->AddEntry(band_fill_toy_2s,
                    draw_band_asy ? "#pm 2#sigma (toy MC)" : "#pm 2#sigma",
                    "f");
  }
  if (draw_band_asy) {
    if (band_asy.median) {
      TGraph* le_med = band_asy_median_vis ? band_asy_median_vis : band_asy.median;
      leg->AddEntry(le_med,
                    draw_band_toy ? "Median expected (asymptotic)"
                                  : "Median expected",
                    "l");
    }
    if (!draw_band_toy && band_fill_asy_1s)
      leg->AddEntry(band_fill_asy_1s, "#pm 1#sigma", "f");
    if (!draw_band_toy && band_fill_asy_2s)
      leg->AddEntry(band_fill_asy_2s, "#pm 2#sigma", "f");
  }

  if (migdal_mode) {
    // Migdal legend order:
    //   1. Excluded region envelope (fill)
    //   2. Individual experiments
    //   3. Our scan curve(s) — last input first (heaviest DC → most sensitive)
    if (g_migdal_env_fill) leg->AddEntry(g_migdal_env_fill, "Direct Detection limits", "f");
    if (g_xenonnT)    leg->AddEntry(g_xenonnT,    "XENON1T (2019, Migdal)", "l");
    if (g_panda4T)    leg->AddEntry(g_panda4T,    "PandaX-4T (2023, Migdal)", "l");
    if (g_darkside50) leg->AddEntry(g_darkside50, "DarkSide-50 (Migdal)", "l");
    if (g_sensei)     leg->AddEntry(g_sensei,     "SENSEI (Migdal)", "l");
    if (g_damic_2025) leg->AddEntry(g_damic_2025, "DAMIC-M (2025 LBC, Migdal)", "l");
    for (int i = static_cast<int>(glimits.size()) - 1; i >= 0; --i) {
      if (!glimits[i]) continue;
      const std::string lbl =
          (i < static_cast<int>(legend_entries_per_graph.size())
               ? legend_entries_per_graph[i] : std::string("scan"));
      leg->AddEntry(glimits[i], lbl.c_str(), "l");
    }
    const char* fdm_label = (mediator == "heavy") ? "#bf{F_{DM} = 1}"
                                                   : "#bf{F_{DM} = ( #alpha m_{N} / q )^{2}}";
    leg->AddEntry((TObject*)0, fdm_label, "");
  } else if (wimp_mode) {
    // WIMP-nucleon SI legend order:
    //   1. Observed limit (DAMIC 2016 SNOLAB, Fig. 11, digitized directly -- trusted)
    //   2. Expected median + +-1sigma band edges
    //   3. Context: DAMIC-M projected (different experiment)
    //   4. Secondary/uncertain references
    //   5. Our scan curve(s)
    if (g_damic_2025)   leg->AddEntry(g_damic_2025,   "DAMIC 2016 SNOLAB, Fig. 11 (observed)", "l");
    if (g_wimp_lower)   leg->AddEntry(g_wimp_lower,   "DAMIC 2016 expected #pm1#sigma, Fig. 11 (digitized)", "l");
    if (!wimp_minimal) {
      if (g_wimp_nominal) leg->AddEntry(g_wimp_nominal, "DAMIC 2016 expected median (uncertain prov.)", "l");
      if (g_wimp_damicm_proj) leg->AddEntry(g_wimp_damicm_proj, "DAMIC-M projected, 1 kg-yr (different exp.)", "l");
      if (g_wimp_radomir) leg->AddEntry(g_wimp_radomir, "radomir's curve (secondary, uncertain)", "l");
      if (g_wimp_snolab_uncertain) leg->AddEntry(g_wimp_snolab_uncertain, "\"damic at snolab.csv\" (uncertain prov.)", "l");
    }
    for (int i = static_cast<int>(glimits.size()) - 1; i >= 0; --i) {
      if (!glimits[i]) continue;
      const std::string lbl =
          (i < static_cast<int>(legend_entries_per_graph.size())
               ? legend_entries_per_graph[i] : std::string("scan"));
      leg->AddEntry(glimits[i], lbl.c_str(), "l");
    }
  } else if (dark_photon) {
    // Dark photon legend order:
    //   1. Direct Detection limits (fill)
    //   2. Stellar limits (fill)
    //   3. DAMIC-M 2025 LBC (solid line)
    //   4. DC band (fill) + central line, or CLI curves in reverse order
    if (g_direct_det_fill) leg->AddEntry(g_direct_det_fill, "Direct Detection limits", "f");
    if (g_darkside50_fill) leg->AddEntry(g_darkside50_fill, "Stellar limits", "f");
    else if (g_darkside50) leg->AddEntry(g_darkside50,      "Stellar limits", "l");
    if (g_damic_2025)      leg->AddEntry(g_damic_2025,      "DAMIC-M (2025 LBC)", "l");
    if (draw_dc_band && g_dc_band && !glimits.empty() && glimits[0]) {
      const std::string cen_lbl =
          legend_entries_per_graph.empty() ? std::string("SrCd_{2}Sb_{2}, 1 kg-yr, DC = 10^{-3} e/pix/d")
                                           : legend_entries_per_graph[0];
      leg->AddEntry(glimits[0], cen_lbl.c_str(), "l");
      leg->AddEntry(g_dc_band,  "Band: DC #in [10^{-6}, 10^{-1}] e/pix/d", "f");
      if (reference_graph_index < glimits.size() && glimits[reference_graph_index]) {
        const std::string ref_lbl =
            (reference_graph_index < legend_entries_per_graph.size())
                ? legend_entries_per_graph[reference_graph_index] : std::string("DAMIC-M (1 kg-yr)");
        leg->AddEntry(glimits[reference_graph_index], ref_lbl.c_str(), "l");
      }
    } else {
      for (int i = static_cast<int>(glimits.size()) - 1; i >= 0; --i) {
        if (!glimits[i]) continue;
        const std::string lbl =
            (i < static_cast<int>(legend_entries_per_graph.size())
                 ? legend_entries_per_graph[i] : std::string("scan"));
        leg->AddEntry(glimits[i], lbl.c_str(), "l");
      }
    }
    leg->AddEntry((TObject*)0, "#bf{Hidden photon}", "");
  } else {
    if (g_env_fill_electron)    leg->AddEntry(g_env_fill_electron, "Direct Detection limits", "f");
    if (g_solar_reflected_fill) leg->AddEntry(g_solar_reflected_fill, "Stellar Limits", "f");
    if (g_damic_2025)           leg->AddEntry(g_damic_2025, "DAMIC-M (2025)", "l");
    if (mediator == "heavy") {
      if (g_model_fill) leg->AddEntry(g_model_fill, "Freeze-out target", "f");
    } else {
      if (g_model) leg->AddEntry(g_model, "Freeze-in target", "l");
    }
    if (mediator == "heavy") {
      leg->AddEntry((TObject*)0, "#bf{F_{DM} = 1}", "");
    } else {
      leg->AddEntry((TObject*)0, "#bf{F_{DM} = ( #alpha  m_{e} / q )^{2}}", "");
    }
    for (size_t i = 0; i < glimits.size(); ++i) {
      if (!glimits[i]) continue;
      const std::string lbl =
          (i < legend_entries_per_graph.size()
               ? legend_entries_per_graph[i] : std::string("scan"));
      leg->AddEntry(glimits[i], lbl.c_str(), "l");
    }
  }

  leg->Draw();

  if (!plot_title.empty()) {
    TLatex title_tex;
    title_tex.SetNDC();
    title_tex.SetTextAlign(22);
    title_tex.SetTextFont(42);
    title_tex.SetTextSize(0.045);
    title_tex.DrawLatex(0.5, 0.97, plot_title.c_str());
  }

  TLatex latex;
  latex.SetTextFont(42);
  latex.SetTextSize(0.035);
  latex.SetNDC(false);  // use axis coordinates, not 0–1

  // DAMIC-M 1 kg-year (solid black)
  latex.SetTextColor(kBlack);
  // latex.DrawLatex(3.0, 2e-40, "DAMIC-M, Pattern Analysis-CCDarkSens");

  if (show_damic && g_damic_2025) {
    latex.SetTextColor(g_damic_2025->GetLineColor());
    latex.DrawLatex(7.0, 2e-37, "paper export (ScienceRun2024)");
  }

  // // SENSEI
  // latex.SetTextColor(g_sensei->GetLineColor());
  // latex.DrawLatex(4.0, 7e-36, "SENSEI (2025)");

  // // SuperCDMS
  // latex.SetTextColor(g_supercdms->GetLineColor());
  // latex.DrawLatex(8.0, 5e-33, "SuperCDMS (2025)");

  // // DarkSide-50
  // latex.SetTextColor(g_darkside50->GetLineColor());
  // latex.DrawLatex(40.0, 3e-35, "DarkSide-50 (2023)");

  // // PandaX-4T
  // latex.SetTextColor(g_panda4T->GetLineColor());
  // latex.DrawLatex(120.0, 8e-36, "PandaX-4T (2023)");

  // // Xenon-nT
  // latex.SetTextColor(g_xenonnT->GetLineColor());
  // latex.DrawLatex(100.0, 3e-34, "XENON-1T/nT (2019,2025)");

  c->Update();
  // Enforce final visible window for both interactive and batch exports.
  ApplyPadAxisRangeUser(x_view_lo, x_view_hi, y_view_lo, y_view_hi);

  const bool multi_inputs = (glimits.size() > 1);
  const std::string out_pdf = out_pdf_cli.empty()
      ? (multi_inputs ? "dme_limit_curve_multi.pdf" : "dme_limit_curve.pdf")
      : out_pdf_cli;
  const std::string out_root = out_root_cli.empty()
      ? (multi_inputs ? "dme_limit_curve_multi.root" : "dme_limit_curve.root")
      : out_root_cli;

  if (band_per_toy_hists) {
    const std::string obs_scan =
        in_paths.empty() ? std::string() : in_paths.front();
    (void)draw_band_per_toy_sigma_hists(band_path, band_per_toy_which, obs_scan);
  }

  // return 0;





  // Legacy Migdal comparison canvas (batch DM-e only — skip for dark photon, migdal, and wimp mode).
  TCanvas* c_migdal_limit = nullptr;
  if (batch && !dark_photon && !migdal_mode && !wimp_mode && mediator != "intermediate") {
  // Now we make a new plot for the migdal case
  // ------------------------------------------------------------------

  TGraph* g_damic_2025_migdal   = LoadTxtToTGraph("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/Migdal/DAMIC-M_2025_DMn-Migdal_heavymediator.txt",
                                            kBlack, 1, 3);

  TGraph* g_damic_wimp = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/WIMP/damicm_fit.csv",
                                            kBlack, 2, 3,1000,1);

  TGraph* g_xenonnt_migdal    = LoadMassLimitDAT("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/Migdal/XENON1T-Migdal.dat",
                                            1e-39, c_xenonnt, 1, 2);

  TGraph* g_xenonnt_wimp    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/WIMP/XENONnT_2025.csv",
                                            c_xenonnt, 1, 2,1000,1);

  // TGraph* g_panda4T_migdal    = LoadMassLimitDAT("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/Migdal/PANDA4T-Migdal.dat",
  //                                           1e-39, c_xenonnt, 1, 2);

  TGraph* g_panda4T_migdal    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/Migdal/Pandax_Migdal_2023.csv",
                                            c_xenonnt, 1, 2,1000,1);

  TGraph* g_panda4T_wimp    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/WIMP/PandaX4T_2025.csv",
                                            c_xenonnt, 1, 2, 1000,1);

  TGraph* g_LZ_wimp    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/WIMP/LZ_2025.csv",
                                            c_xenonnt, 1, 2, 1000,1);

  TGraph* g_darkside50_wimp   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/WIMP/DarkSide50_2025.csv",
                                            c_xenonnt, 1, 2, 1000,1);

  TGraph* g_darkside50_migdal   = LoadMassLimitDAT("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/Migdal/Darkside-Migdal.dat",
                                            1e-39, c_xenonnt, 1, 2);

   TGraph* g_sensei_migdal   = LoadMassLimitDAT("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/Migdal/SENSEI-Migdal.dat",
                                            1e-35, c_xenonnt, 1, 2);

  TGraph* g_pandax_2022_wimp   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/WIMP/PandaX_4T_Wimp_2022.csv",
                                            c_xenonnt, 1, 2, 1000,1);

   TGraph* g_xenon1t_2020_wimp   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/WIMP/Xenon_1T_Wimp_2020.csv",
                                            c_xenonnt, 1, 2, 1000,1);

                                

  

  // TGraph* g_sensei_migdal  = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/SENSEI.csv",
  //                                           c_xenonnt, 2, 3);

  // TGraph* g_supercdms_migdal    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/SuperCDMS.csv",
  //                                           c_xenonnt, 2, 3);

  // TGraph* g_darkside50_migdal   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/DarkSide50.csv",
  //                                           c_xenonnt, 2, 3);

  // TGraph* g_panda4T_migdal      = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/Panda4T.csv",
  //                                           c_xenonnt, 2, 3);

  // TGraph* g_xenonnT_migdal    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/Xenon.csv",
  //                                           c_xenonnt, 2, 3);

  
  TGraph* g_damic_2025_migdal_projection = g_damic_2025_migdal
      ? (TGraph*) g_damic_2025_migdal->Clone("g_damic_2025_migdal_projection") : nullptr;
  double F_improvement = 3000;

  if (g_damic_2025_migdal_projection) {
    for (int i = 0; i < g_damic_2025_migdal_projection->GetN(); ++i) {
        double x, y;
        g_damic_2025_migdal_projection->GetPoint(i, x, y);
        g_damic_2025_migdal_projection->SetPoint(i, x, y / F_improvement);
    }
  }

  // Include projections for WIMP searches with DAMIC-M
  TGraph* g_damic_2025_wimp_det_effects = (TGraph*) g_damic_wimp->Clone("g_damic_2025_wimp_det_effects");
  TGraph* g_damic_2025_wimp_diego = (TGraph*) g_damic_wimp->Clone("g_damic_2025_wimp_diego");
  double F_change_wimp = 4;
  for (int i = 0; i < g_damic_2025_wimp_det_effects->GetN(); ++i) {
      double x, y;
      g_damic_2025_wimp_det_effects->GetPoint(i, x, y);
      g_damic_2025_wimp_det_effects->SetPoint(i, x, (y * F_change_wimp));  // factor > 1 => better sensitivity

      g_damic_2025_wimp_diego->GetPoint(i, x, y);
      g_damic_2025_wimp_diego->SetPoint(i, x, y);  // factor > 1 => better sensitivity

      // std::cout<<"[limit] DAMIC-M WIMP proj: mchi = " << x << " MeV/c², sigma_lim = " << (y * F_change_wimp) << " cm² (det effects), "
      //          << y << " cm² (Diego's estimate)" << std::endl;
  }


   xmin = *std::min_element(mchi_vals.begin(), mchi_vals.end());
   xmax = *std::max_element(mchi_vals.begin(), mchi_vals.end());
   ymin = *std::min_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
   ymax = *std::max_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025_migdal);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_wimp);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025_wimp_det_effects);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025_wimp_diego);
  
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_xenonnt_migdal);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_xenonnt_wimp);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_panda4T_migdal);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_panda4T_wimp);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_LZ_wimp);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_darkside50_wimp);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_darkside50_migdal);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_sensei_migdal);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025_migdal_projection);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_pandax_2022_wimp);
  // UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_xenon1t_2020_wimp);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_sensei);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_supercdms);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_darkside50);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_panda4T);
// //   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_sensei_2025);
// //   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_freeze_in);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_model);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damicm_mike);
  // xmax = 4e3; // limit x-axis to 1 GeV for better visualization

  // y-bounds of the plot (must match glimit->SetMinimum/Maximum)
 y_bottom = ymin * 0.3;
 y_top    = ymax * 3.0;

//  y_bottom = 1e-43 * 0.3;


  TGraph* g_damic_2025_migdal_fill = MakeFilledBandAbove(
      g_damic_2025_migdal, y_top,
      c_damic2025, 0.20);

  TGraph* g_xenonnt_migdal_fill = MakeFilledBandAbove(
      g_xenonnt_migdal, y_top,
      c_damic2025, 0.25);

  TGraph* g_panda4T_migdal_fill = MakeFilledBandAbove(
      g_panda4T_migdal, y_top,
      c_damic2025, 0.25);

  TGraph* g_damic_wimp_fill = MakeFilledBandAbove(
      g_damic_wimp, y_top,
      c_damic2025, 0.25);

  TGraph* g_xenonnt_wimp_fill = MakeFilledBandAbove(
      g_xenonnt_wimp, y_top,
      c_damic2025, 0.25);

  TGraph* g_panda4T_wimp_fill = MakeFilledBandAbove(
      g_panda4T_wimp, y_top,
      c_damic2025, 0.25);

  TGraph* g_LZ_wimp_fill = MakeFilledBandAbove(
      g_LZ_wimp, y_top,
      c_damic2025, 0.25);

  TGraph* g_darkside50_wimp_fill = MakeFilledBandAbove(
      g_darkside50_wimp, y_top,
      c_damic2025, 0.25);

  TGraph* g_darkside50_migdal_fill = MakeFilledBandAbove(
      g_darkside50_migdal, y_top,
      c_damic2025, 0.25);

  TGraph* g_sensei_migdal_fill = MakeFilledBandAbove(
      g_sensei_migdal, y_top,
      c_damic2025, 0.25);

  TGraph* g_damic_2025_migdal_projection_fill = MakeFilledBandAbove(
      g_damic_2025_migdal_projection, y_top,
      c_damic2025, 0.25);

  TGraph* g_pandax_2022_wimp_fill = MakeFilledBandAbove(
      g_pandax_2022_wimp, y_top,
      c_damic2025, 0.25);

  TGraph* g_xenon1t_2020_wimp_fill = MakeFilledBandAbove(
      g_xenon1t_2020_wimp, y_top,
      c_damic2025, 0.25);

  // y-bounds
// y_bottom = ymin * 0.3;
// y_top    = ymax * 3.0;

// collect all curves that should define the excluded region
std::vector<const TGraph*> curves_for_envelope = {
    
    // g_damic_wimp,
    g_xenonnt_migdal,
    // g_xenonnt_wimp,
    // g_damic_2025_migdal,
    g_panda4T_migdal,
    g_panda4T_wimp,
    g_LZ_wimp,
    g_darkside50_wimp,
    g_darkside50_migdal,
    g_sensei_migdal,
    g_pandax_2022_wimp,
    g_xenon1t_2020_wimp
};

double x_patch_min = 2.0;    // MeV/c^2, tune these by eye
double x_patch_max = 500.0;  // MeV/c^2

TGraph* g_migdal_patch = MakeFilledBetween(
    g_damic_2025_migdal,   // bottom curve
    g_sensei_migdal,       // top curve
    x_patch_min,
    x_patch_max,
    800,                   // dense sampling
    c_damic2025,           // same green
    0.25                   // same alpha
);

TGraph* g_env = MakeLowerEnvelope(curves_for_envelope, xmin, xmax);
TGraph* g_env_fill = MakeFilledBandAbove(g_env, y_top, c_damic2025, 0.25);
// // EnforceMonotonicEnvelope(g_env);

// // single filled band from the envelope
// TGraph* g_env_fill = nullptr;
// if (g_env) {
  g_env_fill = MakeFilledBandAbove(g_env, y_top, c_damic2025, 0.25);
  
// }

  

  

  


  c_migdal_limit = new TCanvas("c_migdal_limit", "Migdal limit", 900, 700);
  c_migdal_limit->SetLogx();
  c_migdal_limit->SetLogy();
  c_migdal_limit->SetTicks(1,1);
  c_migdal_limit->SetEditable(kTRUE);
  c_migdal_limit->SetLeftMargin(0.14);
  c_migdal_limit->SetRightMargin(0.04);
  c_migdal_limit->SetBottomMargin(0.12);
  c_migdal_limit->SetTopMargin(0.06);

  g_damic_wimp->SetLineWidth(3);
  g_damic_wimp->SetLineColor(kBlack);
  g_damic_wimp->SetLineStyle(1);

  g_damic_2025_wimp_det_effects->SetLineWidth(3);
  g_damic_2025_wimp_det_effects->SetLineColor(kBlue+2);
  g_damic_2025_wimp_det_effects->SetLineStyle(1);

  g_damic_2025_wimp_diego->SetLineWidth(3);
  g_damic_2025_wimp_diego->SetLineColor(kRed+2);
  g_damic_2025_wimp_diego->SetLineStyle(1);



  const char* kXtitleM = "m_{#chi} [MeV/c^{2}]";
  const char* kYtitleM = "#bar{#sigma}_{n} [cm^{2}]";
  g_damic_wimp->GetXaxis()->SetTitle(kXtitleM);
  g_damic_wimp->GetYaxis()->SetTitle(kYtitleM);
  g_damic_wimp->GetXaxis()->SetTitleSize(0.05);
  g_damic_wimp->GetYaxis()->SetTitleSize(0.05);
  g_damic_wimp->GetXaxis()->SetLabelSize(0.04);
  g_damic_wimp->GetYaxis()->SetLabelSize(0.04);
  g_damic_wimp->GetXaxis()->SetLimits(kMassPlotMin, kMassPlotMax);
  g_damic_wimp->SetMinimum(y_bottom);
  g_damic_wimp->SetMaximum(y_plot_max);
  // Draw axes (and a first version of the line)
  // g_damic_2025_migdal_projection->Draw("AL");
  // ---- Filled region (ONE band) ----
  // if (g_env_fill) g_env_fill->Draw("F SAME");
  // if (g_migdal_patch)  g_migdal_patch->Draw("F SAME");

  // ---- Filled regions (behind lines) ----
  // if (g_damic_2025_migdal_projection_fill)  g_damic_2025_migdal_projection_fill->Draw("F SAME");
  // if (g_damic_2025_migdal_fill)  g_damic_2025_migdal_fill->Draw("F SAME");
  // if (g_xenonnt_migdal_fill)  g_xenonnt_migdal_fill->Draw("F SAME");
  // if (g_panda4T_migdal_fill)  g_panda4T_migdal_fill->Draw("F SAME");
  // if (g_damic_2025_migdal_fill)  g_damic_2025_migdal_fill->Draw("F SAME");

  // // if (g_xenonnt_wimp_fill)  g_xenonnt_wimp_fill->Draw("F SAME");
  // // if (g_panda4T_wimp_fill)  g_panda4T_wimp_fill->Draw("F SAME");
  // if (g_LZ_wimp_fill)  g_LZ_wimp_fill->Draw("F SAME");
  // if (g_darkside50_wimp_fill)  g_darkside50_wimp_fill->Draw("F SAME");
  // if (g_xenonnt_migdal_fill)  g_xenonnt_migdal_fill->Draw("F SAME");
  // if (g_sensei_fill)      g_sensei_fill->Draw("F SAME");
  // if (g_supercdms_fill)   g_supercdms_fill->Draw("F SAME");
  // if (g_darkside50_fill)  g_darkside50_fill->Draw("F SAME");
  // if (g_panda4T_fill)     g_panda4T_fill->Draw("F SAME");
  // if (g_xenonnT_fill)     g_xenonnT_fill->Draw("F SAME");
  // if (g_model_fill)       g_model_fill->Draw("F SAME");
  // if (g_damicm_mike_fill) g_damicm_mike_fill->Draw("F SAME");

  // ---- Line contours on top ----
  // if (g_damic_2025_migdal)  g_damic_2025->Draw("L SAME");
  // if (g_sensei)      g_sensei->Draw("L SAME");
  // if (g_supercdms)   g_supercdms->Draw("L SAME");
  // if (g_darkside50)  g_darkside50->Draw("L SAME");
  // if (g_panda4T)     g_panda4T->Draw("L SAME");
  // if (g_xenonnT)     g_xenonnT->Draw("L SAME");
  // if (g_model)       g_model->Draw("L SAME");
  // g_damic_2025_migdal->Draw("L SAME");
  g_damic_wimp->Draw("AL");
  g_damic_2025_wimp_det_effects->Draw("L SAME");
  g_damic_2025_wimp_diego->Draw("L SAME");
  
  // g_xenonnt_migdal->Draw("L SAME");
  // g_xenonnt_wimp->Draw("L SAME");
  // g_panda4T_migdal->Draw("L SAME");
  // g_darkside50_wimp->Draw("L SAME");
  // g_panda4T_wimp->Draw("L SAME");
  // g_LZ_wimp->Draw("L SAME");
  // g_sensei_migdal->Draw("L SAME");
  // g_pandax_2022_wimp->Draw("L SAME");
  // g_xenon1t_2020_wimp->Draw("L SAME");
  // g_darkside50_migdal->Draw("L SAME");
  
  // if (g_damicm_mike) g_damicm_mike->Draw("L SAME");

  // Finally, bring our solid DAMIC-M curve to the very front
  // g_damic_2025_migdal_projection->Draw("L SAME");

  TLegend* leg_wimp = new TLegend(0.65, 0.55, 0.92, 0.92); // move it to avoid overlapping curves
  leg_wimp->SetBorderSize(0);
  leg_wimp->SetFillStyle(0);
  leg_wimp->SetTextSize(0.03);
  leg_wimp->AddEntry(g_damic_wimp,      "DAMIC-M, 1 kg-year (Cinyu)", "l");
  leg_wimp->AddEntry(g_damic_2025_wimp_diego,      "DAMIC-M, 1 kg-year (Diego)", "l");
  leg_wimp->AddEntry(g_damic_2025_wimp_det_effects,      "DAMIC-M, 1 kg-year (Diego w/ det. effects)", "l");
  leg_wimp->Draw();

  gPad->RedrawAxis();

  
  TLatex latex_migdal;
  latex_migdal.SetTextFont(42);
  latex_migdal.SetTextSize(0.035);
  latex_migdal.SetNDC(false);  // use axis coordinates, not 0–1
  
  
  
  c_migdal_limit->Update();

  const std::string migdal_out_pdf = "dme_migdal_limit_curve.pdf";
  const std::string migdal_out_root = "dme_migdal_limit_curve.root";
  c_migdal_limit->SaveAs(migdal_out_pdf.c_str());
  TFile fout_migdal(migdal_out_root.c_str(), "RECREATE");
  if (g_damic_2025_migdal) g_damic_2025_migdal->Write();
  if (hq) hq->Write();
  fout_migdal.Close();
  std::cout << "[limit] Saved: " << migdal_out_pdf << "\n"
            << "[limit]         " << migdal_out_root << "\n";
  }  // batch: migdal canvas

  // Store-only curves: extracted from ROOT files, written to output but not drawn
  std::vector<std::pair<TGraph*, std::string>> store_only_graphs;
  for (const auto& [so_path, so_label] : store_only_paths) {
    TGraph* g = ExtractLimitCurveFromRoot(so_path, q_thr);
    if (g) {
      // sanitize label into a ROOT-safe name
      std::string root_name = so_label;
      for (char& ch : root_name) if (!std::isalnum(static_cast<unsigned char>(ch))) ch = '_';
      store_only_graphs.emplace_back(g, root_name);
    } else {
      std::cerr << "[limit] WARNING: --store-only could not extract limit curve from " << so_path << "\n";
    }
  }

  // Literature curves to persist alongside scan limits
  std::vector<std::pair<TGraph*, std::string>> literature_curves;
  if (migdal_mode) {
    if (g_damic_2025)  literature_curves.push_back({g_damic_2025,  "lit_damic_2025_migdal"});
    if (g_xenonnT)     literature_curves.push_back({g_xenonnT,     "lit_xenon1t_migdal"});
    if (g_panda4T)     literature_curves.push_back({g_panda4T,     "lit_pandax_migdal"});
    if (g_darkside50)  literature_curves.push_back({g_darkside50,  "lit_darkside_migdal"});
    if (g_sensei)      literature_curves.push_back({g_sensei,      "lit_sensei_migdal"});
    if (g_migdal_env)  literature_curves.push_back({g_migdal_env,  "lit_migdal_envelope"});
  } else if (wimp_mode) {
    if (g_damic_2025)   literature_curves.push_back({g_damic_2025,   "lit_damic2016_observed"});
    if (g_wimp_nominal) literature_curves.push_back({g_wimp_nominal, "lit_damic2016_expected_nominal"});
    if (g_wimp_lower)   literature_curves.push_back({g_wimp_lower,   "lit_damic2016_expected_lower"});
    if (g_wimp_upper)   literature_curves.push_back({g_wimp_upper,   "lit_damic2016_expected_upper"});
    if (g_wimp_damicm_proj) literature_curves.push_back({g_wimp_damicm_proj, "lit_damicm_projected_1kgyear"});
    if (g_wimp_snolab_uncertain) literature_curves.push_back({g_wimp_snolab_uncertain, "lit_damic2016_snolab_uncertain"});
    if (g_wimp_radomir) literature_curves.push_back({g_wimp_radomir, "lit_damic2016_radomir_secondary"});
  } else if (!dark_photon) {
    if (g_damic_2025)        literature_curves.push_back({g_damic_2025,        "lit_damic_2025"});
    if (g_env_electron)      literature_curves.push_back({g_env_electron,      "lit_dd_envelope"});
    if (g_solar_reflected)   literature_curves.push_back({g_solar_reflected,   "lit_stellar"});
    if (g_model)             literature_curves.push_back({g_model,             mediator == "heavy" ? "lit_freeze_out" : "lit_freeze_in"});
  }

  if (!batch) {
    c->SetBit(kCanDelete, false);
    SaveMainLimitOutputs(c, glimits, hq, out_pdf, out_root, "initial", literature_curves, store_only_graphs, legend_entries_per_graph);
    std::cout << "[limit] Interactive mode:\n"
              << "  - Double-click an axis to set log-scale min/max (best for zoom).\n"
              << "  - Or use View menu -> Zoom / zoom tool, then drag a box on the pad.\n"
              << "  - Close the window when done to refresh PDF/ROOT with final axes.\n";
    app.Run();
    SaveMainLimitOutputs(c, glimits, hq, out_pdf, out_root, "final", literature_curves, store_only_graphs, legend_entries_per_graph);
  } else {
    SaveMainLimitOutputs(c, glimits, hq, out_pdf, out_root, "batch", literature_curves, store_only_graphs, legend_entries_per_graph);
  }

  return 0;





  







}
