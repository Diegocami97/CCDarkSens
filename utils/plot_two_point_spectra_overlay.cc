// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  plot_two_point_spectra_overlay.cc -- ROOT macro that overlays dRdE and
//  n_e/pattern spectra from a two-point dump ROOT file into comparison PDFs
//  under outplots/.
// ===========================================================================

#include <TFile.h>
#include <TKey.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>

#include <regex>
#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <limits>

// ROOT macro usage:
//   root -l -b -q 'utils/plot_two_point_spectra_overlay.cc("outputs/two_point_spectra_dump_ne/scan_dmelectron_pattern.root")'
void plot_two_point_spectra_overlay(const char* infile_c = "")
{
  const std::string infile = (infile_c && infile_c[0]) ? infile_c : std::string{};
  if (infile.empty()) {
    std::cerr << "plot_two_point_spectra_overlay: missing input ROOT file\n";
    return;
  }

  gStyle->SetOptStat(0);
  gSystem->mkdir("outplots", /*recursive=*/true);

  TFile* f = TFile::Open(infile.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "plot_two_point_spectra_overlay: failed to open " << infile << "\n";
    return;
  }

  auto starts_with = [](const std::string& s, const std::string& p) -> bool {
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
  };

  // Collect suffix tags from dRdE__<tag> histograms.
  const std::string prefix_dRdE = "dRdE__";
  std::set<std::string> pointTags;

  TIter itKeys(f->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(itKeys())) {
    const std::string name = key->GetName();
    if (starts_with(name, prefix_dRdE)) {
      pointTags.insert(name.substr(prefix_dRdE.size())); // everything after dRdE__
    }
  }

  if (pointTags.empty()) {
    std::cerr << "plot_two_point_spectra_overlay: no dRdE__* histograms found in file\n";
    return;
  }

  std::cout << "plot_two_point_spectra_overlay: found " << pointTags.size()
            << " point tags in " << infile << "\n";
  for (const auto& tag : pointTags) {
    std::cout << "  tag=" << tag << "\n";
  }

  auto format_trim = [](double v, int precision) -> std::string {
    if (!std::isfinite(v)) return "nan";
    if (v == 0.0) return "0";
    std::ostringstream ss;
    // General formatting (no fixed) so 10/100 don't show decimals.
    ss << std::setprecision(precision) << v;
    std::string s = ss.str();
    // Trim trailing zeros after decimal if present.
    auto dot = s.find('.');
    if (dot != std::string::npos) {
      while (!s.empty() && s.back() == '0') s.pop_back();
      if (!s.empty() && s.back() == '.') s.pop_back();
    }
    return s.empty() ? "0" : s;
  };

  // Convert sanitized tags back to readable numbers.
  auto mchiTagToLabel = [](const std::string& mchi_tag) -> std::string {
    // mchi_tag example: "0p500000"
    std::string x = mchi_tag;
    for (char& c : x) if (c == 'p') c = '.';
    return x;
  };

  auto sigmaTagToLabel = [](std::string sigma_tag) -> std::string {
    // sigma_tag example: "5p0em37" (from sanitize_root_name and format_sigma)
    // em -> e- and ep -> e+
    sigma_tag = std::regex_replace(sigma_tag, std::regex("em"), "e-");
    sigma_tag = std::regex_replace(sigma_tag, std::regex("ep"), "e+");
    for (char& c : sigma_tag) if (c == 'p') c = '.';
    return sigma_tag;
  };

  auto tagToLegend = [&](const std::string& tag) -> std::string {
    // Expected tag: mchi_<mchi>__sigma_<sigma>
    std::smatch m;
    const std::regex re(R"(mchi_([^_]+)__sigma_([^_]+))");
    if (std::regex_search(tag, m, re) && m.size() >= 3) {
      const std::string mchi_raw = mchiTagToLabel(m[1].str()); // e.g. "10.000000"
      const std::string sigma_raw = sigmaTagToLabel(m[2].str()); // e.g. "5.0e-37"

      double mchi_val = 0.0;
      double sigma_val = 0.0;
      try { mchi_val = std::stod(mchi_raw); } catch (...) { mchi_val = 0.0; }
      try { sigma_val = std::stod(sigma_raw); } catch (...) { sigma_val = 0.0; }

      const std::string mchi_s = format_trim(mchi_val, /*precision=*/10);

      std::ostringstream ss;
      ss << std::setprecision(6) << std::scientific << sigma_val; // e.g. 5.000000e-37
      std::string sigma_s = ss.str();
      // Trim trailing zeros in mantissa: "5.00000e-37" -> "5e-37"
      auto epos = sigma_s.find('e');
      if (epos != std::string::npos) {
        std::string mant = sigma_s.substr(0, epos);
        std::string expn = sigma_s.substr(epos); // includes 'e'
        auto dot = mant.find('.');
        if (dot != std::string::npos) {
          while (!mant.empty() && mant.back() == '0') mant.pop_back();
          if (!mant.empty() && mant.back() == '.') mant.pop_back();
        }
        sigma_s = mant + expn;
      }

      return "m_{#chi}=" + mchi_s + " MeV, #sigma_{e}=" + sigma_s + " cm^{2}";
    }
    return tag;
  };

  auto basename = [](const std::string& path) -> std::string {
    auto slash = path.find_last_of('/');
    const std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
    auto dot = name.find_last_of('.');
    return (dot == std::string::npos) ? name : name.substr(0, dot);
  };

  const std::string base = basename(infile);

  // Families to plot: only those actually present in the file.
  const std::vector<std::string> families = {
    "dRdE__",
    "S_true_ne__",
    "S_obs_ne__",
    "S_pat__"
  };

  auto family_present = [&](const std::string& famPrefix) -> bool {
    // just test first tag
    const std::string anyTag = *pointTags.begin();
    const std::string hname = famPrefix + anyTag;
    return f->Get(hname.c_str()) != nullptr;
  };

  int color_i = 1;
  for (const auto& fam : families) {
    if (!family_present(fam)) continue;

    const std::string famLabel =
      (fam == "dRdE__") ? "dRdE" :
      (fam == "S_true_ne__") ? "S_true_n_e" :
      (fam == "S_obs_ne__") ? "S_obs_n_e" :
      (fam == "S_pat__") ? "S_pat_pattern" : fam;

    TCanvas* c = new TCanvas(Form("c_%s_%s", famLabel.c_str(), base.c_str()), famLabel.c_str(), 900, 650);
    c->SetLogy(true);

    TLegend* leg = new TLegend(0.52, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    struct Item {
      TH1* h = nullptr;
      std::string tag;
      double max = 0.0;
    };

    // Clone histograms and compute y-range using the maximum of the two histograms.
    std::vector<Item> items;
    double max_y = 0.0;
    double min_pos_y = std::numeric_limits<double>::infinity();

    int clone_idx = 0;
    for (const auto& tag : pointTags) {
      const std::string hname = fam + tag;
      TH1* h_in = dynamic_cast<TH1*>(f->Get(hname.c_str()));
      if (!h_in) continue;

      TH1* h = static_cast<TH1*>(h_in->Clone(Form("%s_clone_%d", hname.c_str(), clone_idx)));
      h->SetDirectory(nullptr);

      const double hmax = h->GetMaximum();
      if (hmax > max_y) max_y = hmax;

      for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        const double v = h->GetBinContent(ib);
        if (v > 0.0 && v < min_pos_y) min_pos_y = v;
      }

      items.push_back({h, tag, hmax});
      ++clone_idx;
    }

    if (items.empty()) {
      delete c;
      continue;
    }

    if (items.size() < pointTags.size()) {
      std::cout << "plot_two_point_spectra_overlay: WARNING found only "
                << items.size() << "/" << pointTags.size()
                << " histograms for family '" << fam << "' in " << infile << "\n";
    }

    // Sort by maximum descending so the largest histogram is drawn first.
    std::sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
      return a.max > b.max;
    });

    // Extract into hs vector in sorted order.
    std::vector<TH1*> hs;
    hs.reserve(items.size());

    // Assign colors/styles in the sorted order (max first => red, second => blue).
    int color_i_local = 1;
    for (auto& it : items) {
      TH1* h = it.h;

      const int color = (color_i_local == 1) ? kRed+1 : kBlue+1;
      h->SetLineColor(color);
      h->SetMarkerColor(color);
      h->SetMarkerStyle(20 + color_i_local);
      h->SetLineWidth(2);

      // Produce a lighter fill color by blending the color with white.
      const TColor* col = gROOT->GetColor(color);
      Float_t r = 0.0f, gcol = 0.0f, b = 0.0f;
      col->GetRGB(r, gcol, b);
      const Float_t blendFrac = 0.7f;
      const Float_t rr = r + blendFrac * (1.0f - r);
      const Float_t gg = gcol + blendFrac * (1.0f - gcol);
      const Float_t bb = b + blendFrac * (1.0f - b);
      const Int_t lighterFillColor = TColor::GetColor(rr, gg, bb);
      h->SetFillStyle(3001); // semi-transparent
      h->SetFillColor(lighterFillColor);

      leg->AddEntry(h, tagToLegend(it.tag).c_str(), "l");
      hs.push_back(h);

      ++color_i_local;
      if (color_i_local > 6) color_i_local = 1; // safety
    }

    // Set explicit y-range for consistent log-scale plots.
    if (max_y > 0.0) {
      TH1* h0 = hs.front();
      h0->SetMaximum(max_y);
      if (std::isfinite(min_pos_y) && min_pos_y > 0.0) {
        h0->SetMinimum(min_pos_y * 0.5);
      } else {
        h0->SetMinimum(1e-300);
      }
    }

    // Explicit axis units:
    // - dRdE is a differential rate table: events/(kg·year·eV)
    // - S_true/S_obs are expected event counts after applying exposure and folding/response.
    if (!hs.empty()) {
      TH1* h0 = hs.front();
      if (fam == "dRdE__") {
        h0->GetYaxis()->SetTitle("dR/dE [events/(kg #times year #times eV)]");
      } else if (fam == "S_true_ne__") {
        h0->GetYaxis()->SetTitle("S_{true}(n_{e}) [Counts]");
      } else if (fam == "S_obs_ne__") {
        h0->GetYaxis()->SetTitle("S_{obs}(n_{e}) [Counts]");
      } else if (fam == "S_pat__") {
        h0->GetYaxis()->SetTitle("S_{pat} [Counts]");
      }
      gPad->Modified();
      gPad->Update();
    }

    // Draw all points (use 'hist' for robust visibility on log plots).
    bool firstDraw = true;
    for (auto* h : hs) {
      if (firstDraw) {
        h->Draw("hist");
        firstDraw = false;
      } else {
        h->Draw("hist same");
      }
    }
    leg->Draw();
    c->Update();

    const std::string outpng = "outplots/overlay_" + famLabel + "_" + base + ".png";
    c->SaveAs(outpng.c_str());
    delete c;
  }
}
