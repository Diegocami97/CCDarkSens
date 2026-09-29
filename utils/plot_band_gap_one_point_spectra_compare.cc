// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: plot_band_gap_one_point_spectra_compare.cc
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  File: plot_band_gap_one_point_spectra_compare.cc
//  Overlay dR/dE and S_true(n_e) from one-point band-gap spectra dumps.
//
//  Usage:
//    root -l -b -q 'utils/plot_band_gap_one_point_spectra_compare.cc("outputs/band_gap_one_point_spectra")'
// ===========================================================================

#include <TFile.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TLine.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

struct CaseSpectra {
  std::string id;
  TH1* dRdE = nullptr;
  TH1* S_true = nullptr;
};

// ----------------------------------------------------------------------------
// Basename
//   File name of a path without the directory and without the extension.
// ----------------------------------------------------------------------------
std::string Basename(const std::string& path)
{
  const auto slash = path.find_last_of('/');
  const std::string name = (slash == std::string::npos) ? path : path.substr(slash + 1);
  const auto dot = name.find_last_of('.');
  return (dot == std::string::npos) ? name : name.substr(0, dot);
}

// ----------------------------------------------------------------------------
// FindFirstHist
//   First histogram of the file whose name starts with prefix, or nullptr.
// ----------------------------------------------------------------------------
TH1* FindFirstHist(TFile* f, const char* prefix)
{
  if (!f) return nullptr;
  TIter it(f->GetListOfKeys());
  while (auto* key = static_cast<TKey*>(it())) {
    const std::string name = key->GetName();
    if (name.rfind(prefix, 0) == 0) {
      return dynamic_cast<TH1*>(f->Get(name.c_str()));
    }
  }
  return nullptr;
}

// ----------------------------------------------------------------------------
// LoadCase
//   Load the dR/dE and S_true(n_e) histograms of one case from its ROOT file into out; returns false (with a message) if the file or a histogram is missing.
// ----------------------------------------------------------------------------
bool LoadCase(const std::string& case_id, const std::string& root_path, CaseSpectra& out)
{
  TFile* f = TFile::Open(root_path.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "plot_band_gap_one_point_spectra_compare: cannot open " << root_path << "\n";
    return false;
  }
  out.id = case_id;
  TH1* hd = FindFirstHist(f, "dRdE__");
  TH1* hs = FindFirstHist(f, "S_true_ne__");
  if (!hd || !hs) {
    std::cerr << "plot_band_gap_one_point_spectra_compare: missing dRdE or S_true in "
              << root_path << "\n";
    f->Close();
    return false;
  }
  out.dRdE = static_cast<TH1*>(hd->Clone((case_id + "__dRdE").c_str()));
  out.S_true = static_cast<TH1*>(hs->Clone((case_id + "__S_true").c_str()));
  out.dRdE->SetDirectory(nullptr);
  out.S_true->SetDirectory(nullptr);
  f->Close();
  return true;
}

// ----------------------------------------------------------------------------
// StyleHist
//   Line colour, style and width 2, no markers, and the x-axis title of a histogram.
// ----------------------------------------------------------------------------
void StyleHist(TH1* h, int color, int style, const char* axis_title)
{
  h->SetLineColor(color);
  h->SetLineWidth(2);
  h->SetLineStyle(style);
  h->SetMarkerStyle(0);
  h->GetXaxis()->SetTitle(axis_title);
}

// ----------------------------------------------------------------------------
// IntegralWidth
//   Sum of bin content times bin width (0 for a null histogram).
// ----------------------------------------------------------------------------
double IntegralWidth(TH1* h)
{
  if (!h) return 0.0;
  double sum = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    sum += h->GetBinContent(i) * h->GetBinWidth(i);
  }
  return sum;
}

// ----------------------------------------------------------------------------
// GapEvFromTag
//   Band gap in eV from a tag such as "gap0p7" ('p' is the decimal point); -1 if the tag is malformed.
// ----------------------------------------------------------------------------
double GapEvFromTag(const std::string& gap_tag)
{
  if (gap_tag == "gap1p2") return 1.2;
  if (gap_tag.rfind("gap", 0) != 0 || gap_tag.size() < 4) return -1.0;
  std::string s = gap_tag.substr(3);
  for (char& c : s) {
    if (c == 'p') c = '.';
  }
  return std::atof(s.c_str());
}

// ----------------------------------------------------------------------------
// GapTagFromId
//   Gap tag of a case id: the part before the first underscore.
// ----------------------------------------------------------------------------
std::string GapTagFromId(const std::string& case_id)
{
  const auto us = case_id.find('_');
  return (us == std::string::npos) ? case_id : case_id.substr(0, us);
}

// ----------------------------------------------------------------------------
// ScenarioFromId
//   Scenario of a case id: the part after the first underscore (empty if none).
// ----------------------------------------------------------------------------
std::string ScenarioFromId(const std::string& case_id)
{
  const auto us = case_id.find('_');
  return (us == std::string::npos) ? "" : case_id.substr(us + 1);
}

// ----------------------------------------------------------------------------
// GapLabelEv
//   Band gap as a label such as "0.7 eV" (one decimal, 1.2 eV kept as 1.2).
// ----------------------------------------------------------------------------
std::string GapLabelEv(double gap_eV)
{
  char buf[32];
  if (std::fabs(gap_eV - 1.2) < 1e-9)
    std::snprintf(buf, sizeof(buf), "1.2");
  else
    std::snprintf(buf, sizeof(buf), "%.1f", gap_eV);
  return std::string(buf) + " eV";
}

// ----------------------------------------------------------------------------
// EhLabelEv
//   Electron-hole pair energy as a label with one decimal.
// ----------------------------------------------------------------------------
std::string EhLabelEv(double eh_eV)
{
  char buf[32];
  if (std::fabs(eh_eV - 3.8) < 1e-9)
    std::snprintf(buf, sizeof(buf), "3.8");
  else if (std::fabs(eh_eV - 1.2) < 1e-9)
    std::snprintf(buf, sizeof(buf), "1.2");
  else
    std::snprintf(buf, sizeof(buf), "%.1f", eh_eV);
  return std::string(buf);
}

// ROOT histogram titles: ASCII + TLatex (E_gap, epsilon_h). No UTF-8 punctuation.
std::string EgapEpsilonLegendLabel(double gap_eV, double eh_eV)
{
  return "E_{gap} = " + GapLabelEv(gap_eV) + " eV, #varepsilon_{h} = " + EhLabelEv(eh_eV) +
         " eV";
}

constexpr double kEhBThreshEv = 3.8;
constexpr double kNeZoomLo = 0.5;   // bin edges for n_e = 1
constexpr double kNeZoomHi = 5.5;   // bin edges for n_e = 5

// ----------------------------------------------------------------------------
// EhPairEvForCase
//   eps_h of a case: equal to the gap for D-equal, the fixed B-thresh value for B-thresh, -1 otherwise.
// ----------------------------------------------------------------------------
double EhPairEvForCase(const CaseSpectra& cs)
{
  const std::string scen = ScenarioFromId(cs.id);
  const double gap = GapEvFromTag(GapTagFromId(cs.id));
  if (scen == "D-equal") return gap;
  if (scen == "B-thresh") return kEhBThreshEv;
  return -1.0;
}

// ----------------------------------------------------------------------------
// IonizationLegendLabel
//   Legend label with the (E_gap, eps_h) pair of a case.
// ----------------------------------------------------------------------------
std::string IonizationLegendLabel(const CaseSpectra& cs)
{
  const double gap = GapEvFromTag(GapTagFromId(cs.id));
  const double eh = EhPairEvForCase(cs);
  return EgapEpsilonLegendLabel(gap, eh);
}

// ----------------------------------------------------------------------------
// SetLogYRangeFromHists
//   Set a logarithmic y range for a group of histograms from their largest and smallest positive values (optionally within a bin range); the top is ymax_scale times the maximum.
// ----------------------------------------------------------------------------
void SetLogYRangeFromHists(const std::vector<TH1*>& hists, double ymax_scale = 3.0,
                           int xbin_first = 0, int xbin_last = 0)
{
  double ymax = 0.0;
  double ymin_pos = std::numeric_limits<double>::infinity();
  for (TH1* h : hists) {
    if (!h) continue;
    const int ib0 = (xbin_first > 0) ? xbin_first : 1;
    const int ib1 = (xbin_last > 0) ? xbin_last : h->GetNbinsX();
    for (int ib = ib0; ib <= ib1; ++ib) {
      const double v = h->GetBinContent(ib);
      if (v > ymax) ymax = v;
      if (v > 0.0 && v < ymin_pos) ymin_pos = v;
    }
  }
  if (ymax <= 0.0) return;
  const double ymin = (ymin_pos < std::numeric_limits<double>::infinity()) ? ymin_pos * 0.3
                                                                           : ymax * 1e-6;
  for (TH1* h : hists) {
    if (h) h->GetYaxis()->SetRangeUser(ymin, ymax * ymax_scale);
  }
}

}  // namespace

// ----------------------------------------------------------------------------
// plot_band_gap_one_point_spectra_compare
//   ROOT macro: read every case directory of the one-point spectra dump (default outputs/band_gap_one_point_spectra) and overlay dR/dE and S_true(n_e) of all cases in PDFs under outplots/band_gap_one_point_spectra/.
// ----------------------------------------------------------------------------
void plot_band_gap_one_point_spectra_compare(const char* results_dir_c = "")
{
  gStyle->SetOptStat(0);
  const std::string results_dir =
      (results_dir_c && results_dir_c[0]) ? results_dir_c : "outputs/band_gap_one_point_spectra";

  gSystem->mkdir("outplots", true);
  gSystem->mkdir("outplots/band_gap_one_point_spectra", true);

  std::vector<CaseSpectra> cases;
  TSystemDirectory dir("cases", results_dir.c_str());
  TList* entries = dir.GetListOfFiles();
  if (!entries) {
    std::cerr << "plot_band_gap_one_point_spectra_compare: no directory " << results_dir << "\n";
    return;
  }

  entries->Sort();
  TIter it(entries);
  while (auto* ent = static_cast<TSystemFile*>(it())) {
    const std::string name = ent->GetName();
    if (name == "." || name == "..") continue;
    if (!ent->IsDirectory()) continue;
    // Only standard case ids: gap0p1_D-equal, gap1p2_B-thresh, etc.
    if (name.find("gap") != 0) continue;

    const std::string root_path = results_dir + "/" + name + "/scan_dmelectron_pattern.root";
    if (gSystem->AccessPathName(root_path.c_str())) continue;

    CaseSpectra cs;
    if (LoadCase(name, root_path, cs)) {
      cases.push_back(cs);
      std::cout << "loaded case " << name << "\n";
    }
  }

  if (cases.empty()) {
    std::cerr << "plot_band_gap_one_point_spectra_compare: no cases with spectra found\n";
    return;
  }

  std::sort(cases.begin(), cases.end(),
            [](const CaseSpectra& a, const CaseSpectra& b) { return a.id < b.id; });

  // Summary table
  {
    const std::string tbl = "outplots/band_gap_one_point_spectra/integrals_summary.txt";
    FILE* fp = fopen(tbl.c_str(), "w");
    if (fp) {
      fprintf(fp, "# case  integral_dRdE_dE [events/kg·yr]  integral_S_true [counts]\n");
      for (const auto& c : cases) {
        fprintf(fp, "%s  %.6e  %.6e\n", c.id.c_str(), IntegralWidth(c.dRdE),
                IntegralWidth(c.S_true));
      }
      fclose(fp);
      std::cout << "wrote " << tbl << "\n";
    }
  }

  struct OverlayOpts {
    double x_lo = 0.0;
    double x_hi = 0.0;
    const char* hist_title = nullptr;
    double leg_x1 = 0.50;
    double leg_y1 = 0.50;
    double leg_x2 = 0.88;
    double leg_y2 = 0.88;
  };

  auto draw_overlay = [&](const char* which, const char* ytitle, const char* outfile,
                          const std::vector<CaseSpectra>& subset,
                          const std::function<std::string(const CaseSpectra&)>& legend_label,
                          const OverlayOpts& opts = OverlayOpts{}) {
    if (subset.empty()) return;
    TCanvas c("c", which, 1000, 700);
    c.SetLogy();

    const char* xtitle = (std::string(which) == "dRdE") ? "E_{e} [eV]" : "n_{e}";
    const bool zoom_x = (opts.x_hi > opts.x_lo);
    int xbin_first = 0;
    int xbin_last = 0;

    bool first = true;
    int color = 1;
    std::vector<TH1*> drawn;
    auto* leg = new TLegend(opts.leg_x1, opts.leg_y1, opts.leg_x2, opts.leg_y2);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(zoom_x ? 0.030 : 0.032);

    for (const auto& cs : subset) {
      TH1* h = (std::string(which) == "dRdE") ? cs.dRdE : cs.S_true;
      if (!h) continue;
      StyleHist(h, color, 1, xtitle);
      h->GetYaxis()->SetTitle(ytitle);
      if (opts.hist_title) h->SetTitle(opts.hist_title);
      if (zoom_x) {
        h->GetXaxis()->SetRangeUser(opts.x_lo, opts.x_hi);
        if (xbin_first == 0) {
          xbin_first = h->GetXaxis()->FindBin(opts.x_lo + 1e-6);
          xbin_last = h->GetXaxis()->FindBin(opts.x_hi - 1e-6);
        }
      }
      h->Draw(first ? "HIST" : "HIST SAME");
      leg->AddEntry(h, legend_label(cs).c_str(), "l");
      drawn.push_back(h);
      first = false;
      ++color;
      if (color == 5) color = 6;
      if (color > 9) color = 1;
    }
    if (first) {
      delete leg;
      return;
    }

    SetLogYRangeFromHists(drawn, 3.0, xbin_first, xbin_last);
    leg->Draw();
    c.SaveAs(outfile);
    delete leg;
    std::cout << "wrote " << outfile << "\n";
  };

  auto by_scenario = [&](const std::string& scen) {
    std::vector<CaseSpectra> out;
    for (const auto& cs : cases) {
      if (ScenarioFromId(cs.id) == scen) out.push_back(cs);
    }
    std::sort(out.begin(), out.end(), [](const CaseSpectra& a, const CaseSpectra& b) {
      return GapEvFromTag(GapTagFromId(a.id)) < GapEvFromTag(GapTagFromId(b.id));
    });
    return out;
  };

  const auto label_gap = [](const CaseSpectra& cs) {
    const double gap = GapEvFromTag(GapTagFromId(cs.id));
    return EgapEpsilonLegendLabel(gap, EhPairEvForCase(cs));
  };
  const auto label_ionization = [](const CaseSpectra& cs) { return IonizationLegendLabel(cs); };

  const std::vector<CaseSpectra> d_equal = by_scenario("D-equal");
  const std::vector<CaseSpectra> b_thresh = by_scenario("B-thresh");

  OverlayOpts ne_zoom;
  ne_zoom.x_lo = kNeZoomLo;
  ne_zoom.x_hi = kNeZoomHi;
  ne_zoom.hist_title = "S_{true}(n_{e}); 1 #leq n_{e} #leq 5";
  ne_zoom.leg_x1 = 0.48;
  ne_zoom.leg_y1 = 0.42;
  ne_zoom.leg_x2 = 0.89;
  ne_zoom.leg_y2 = 0.88;

  draw_overlay("dRdE", "dR/dE [events/(kg yr eV)]",
               "outplots/band_gap_one_point_spectra/compare_dRdE_D-equal_all_gaps.pdf",
               d_equal, label_gap);
  draw_overlay("S_true", "expected counts",
               "outplots/band_gap_one_point_spectra/compare_S_true_D-equal_all_gaps.pdf",
               d_equal, label_gap);
  draw_overlay("dRdE", "dR/dE [events/(kg yr eV)]",
               "outplots/band_gap_one_point_spectra/compare_dRdE_B-thresh_all_gaps.pdf",
               b_thresh, label_gap);
  draw_overlay("S_true", "expected counts",
               "outplots/band_gap_one_point_spectra/compare_S_true_B-thresh_all_gaps.pdf",
               b_thresh, label_gap);

  // Ionization-table comparison in n_e = 1..5 (low-ne structure)
  draw_overlay("S_true", "expected counts",
               "outplots/band_gap_one_point_spectra/compare_ionization_D-equal_all_gaps_ne1to5.pdf",
               d_equal, label_ionization, ne_zoom);
  draw_overlay("S_true", "expected counts",
               "outplots/band_gap_one_point_spectra/compare_ionization_B-thresh_all_gaps_ne1to5.pdf",
               b_thresh, label_ionization, ne_zoom);
  draw_overlay("S_true", "expected counts",
               "outplots/band_gap_one_point_spectra/compare_ionization_all_tables_ne1to5.pdf",
               cases, label_ionization, ne_zoom);

  // Per gap: D-equal vs B-thresh on same axes (before / after panels)
  std::vector<std::string> gap_tags;
  for (const auto& cs : cases) {
    const std::string gt = GapTagFromId(cs.id);
    if (std::find(gap_tags.begin(), gap_tags.end(), gt) == gap_tags.end()) {
      gap_tags.push_back(gt);
    }
  }
  std::sort(gap_tags.begin(), gap_tags.end(),
            [](const std::string& a, const std::string& b) {
              return GapEvFromTag(a) < GapEvFromTag(b);
            });

  for (const std::string& gt : gap_tags) {
    const CaseSpectra* d_case = nullptr;
    const CaseSpectra* b_case = nullptr;
    for (const auto& cs : cases) {
      if (GapTagFromId(cs.id) != gt) continue;
      if (ScenarioFromId(cs.id) == "D-equal") d_case = &cs;
      if (ScenarioFromId(cs.id) == "B-thresh") b_case = &cs;
    }
    if (!d_case || !b_case) {
      std::cerr << "plot_band_gap_one_point_spectra_compare: skip " << gt
                << " (missing D-equal or B-thresh)\n";
      continue;
    }

    const double gap_eV = GapEvFromTag(gt);
    const std::string gap_lbl = GapLabelEv(gap_eV);
    const std::string leg_d = EgapEpsilonLegendLabel(gap_eV, gap_eV);
    const std::string leg_b = EgapEpsilonLegendLabel(gap_eV, kEhBThreshEv);

    TCanvas c(("c_cmp_" + gt).c_str(), gt.c_str(), 1100, 500);
    c.Divide(2, 1);
    std::vector<TLegend*> pad_legs;

    auto draw_panel = [&](int pad, TH1* hd, TH1* hb, const std::string& panel_title,
                          const char* xtitle, const char* ytitle) {
      c.cd(pad);
      gPad->SetLogy();
      StyleHist(hd, kBlue + 1, 1, xtitle);
      StyleHist(hb, kRed + 1, 1, xtitle);
      hd->GetYaxis()->SetTitle(ytitle);
      hb->GetYaxis()->SetTitle(ytitle);
      hd->SetTitle(panel_title.c_str());
      hd->Draw("HIST");
      hb->Draw("HIST SAME");
      SetLogYRangeFromHists({hd, hb});
      auto* leg = new TLegend(0.48, 0.62, 0.89, 0.88);
      leg->SetBorderSize(1);
      leg->SetFillStyle(1001);
      leg->SetFillColor(kWhite);
      leg->SetTextSize(0.038);
      leg->AddEntry(hd, leg_d.c_str(), "l");
      leg->AddEntry(hb, leg_b.c_str(), "l");
      leg->Draw();
      pad_legs.push_back(leg);
      gPad->Modified();
      gPad->Update();
    };

    const std::string title_gap = "E_{gap} = " + gap_lbl + " eV";
    draw_panel(1, d_case->dRdE, b_case->dRdE,
               "Before ionization: dR/dE; " + title_gap, "E_{e} [eV]",
               "dR/dE [events/(kg yr eV)]");
    draw_panel(2, d_case->S_true, b_case->S_true,
               "After ionization: S_{true}(n_{e}); " + title_gap, "n_{e}", "expected counts");

    const std::string out =
        "outplots/band_gap_one_point_spectra/compare_D-equal_vs_B-thresh__" + gt + ".pdf";
    c.SaveAs(out.c_str());
    for (TLegend* leg : pad_legs) delete leg;
    std::cout << "wrote " << out << "\n";

    // Single-panel ionization compare at low n_e (same scissor gap, two p100K tables)
    std::vector<CaseSpectra> pair = {*d_case, *b_case};
    const std::string gap_title = "S_{true}(n_{e}); E_{gap} = " + gap_lbl + " eV";
    OverlayOpts gap_zoom = ne_zoom;
    gap_zoom.hist_title = gap_title.c_str();
    const std::string out_ne =
        "outplots/band_gap_one_point_spectra/compare_ionization_ne1to5__" + gt + ".pdf";
    draw_overlay("S_true", "expected counts", out_ne.c_str(), pair, label_ionization, gap_zoom);
  }

  // Per-case before / after (two pads)
  for (const auto& cs : cases) {
    TCanvas c(("c_" + cs.id).c_str(), cs.id.c_str(), 1100, 500);
    c.Divide(2, 1);

    c.cd(1);
    gPad->SetLogy();
    StyleHist(cs.dRdE, kBlack, 1, "E_{e} [eV]");
    cs.dRdE->SetTitle(("Before ionization: dR/dE — " + cs.id).c_str());
    cs.dRdE->GetYaxis()->SetTitle("dR/dE [events/(kg yr eV)]");
    cs.dRdE->Draw("HIST");

    c.cd(2);
    gPad->SetLogy();
    StyleHist(cs.S_true, kBlue + 1, 1, "n_{e}");
    cs.S_true->SetTitle(("After ionization: S_{true}(n_{e}) — " + cs.id).c_str());
    cs.S_true->GetYaxis()->SetTitle("expected counts");
    cs.S_true->Draw("HIST");

    const std::string out =
        "outplots/band_gap_one_point_spectra/before_after__" + cs.id + ".pdf";
    c.SaveAs(out.c_str());
    std::cout << "wrote " << out << "\n";
  }

  for (auto& cs : cases) {
    delete cs.dRdE;
    delete cs.S_true;
  }
}
