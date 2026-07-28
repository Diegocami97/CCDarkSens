// plot_darkphoton_band.C
// Dark photon projection: central line (DC=1e-3) + hatched band (DC=1e-5 to DC=1e-1)
// Run: root -b -q utils/plot_darkphoton_band.C

#include "TFile.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

// Extract 90% CL upper limit curve from q(mchi, sigma) histogram.
// Returns TGraph of (mchi, sigma_UL) scanning each mchi column for q > q_threshold.
TGraph* ExtractLimitCurve(const std::string& root_path, double q_threshold = 2.71) {
    TFile* f = TFile::Open(root_path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        printf("ERROR: cannot open %s\n", root_path.c_str());
        return nullptr;
    }

    // Try common histogram names
    TH2D* h = nullptr;
    for (const char* name : {"q_hist_mchi_sigma", "q_mchi_sigma", "h_q", "qhist"}) {
        h = (TH2D*) f->Get(name);
        if (h) break;
    }
    if (!h) {
        // List keys and try first TH2
        TIter next(f->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*) next())) {
            TObject* obj = key->ReadObj();
            if (obj->InheritsFrom("TH2")) { h = (TH2D*) obj; break; }
        }
    }
    if (!h) { printf("ERROR: no TH2 found in %s\n", root_path.c_str()); f->Close(); return nullptr; }

    std::vector<double> xs, ys;
    int nx = h->GetNbinsX(), ny = h->GetNbinsY();
    for (int ix = 1; ix <= nx; ++ix) {
        double mchi = h->GetXaxis()->GetBinCenter(ix);
        double sigma_ul = -1;
        // Scan from low sigma upward; first bin where q > threshold is the UL
        for (int iy = 1; iy <= ny; ++iy) {
            double q = h->GetBinContent(ix, iy);
            if (q > q_threshold) { sigma_ul = h->GetYaxis()->GetBinCenter(iy); break; }
        }
        if (sigma_ul > 0) { xs.push_back(mchi); ys.push_back(sigma_ul); }
    }
    f->Close();
    if (xs.empty()) return nullptr;
    return new TGraph(xs.size(), xs.data(), ys.data());
}

// Interpolate TGraph at x (log-linear interpolation).
double InterpGraph(TGraph* g, double x) {
    int n = g->GetN();
    double* xs = g->GetX(); double* ys = g->GetY();
    if (x <= xs[0])   return ys[0];
    if (x >= xs[n-1]) return ys[n-1];
    for (int i = 0; i < n-1; ++i) {
        if (x >= xs[i] && x <= xs[i+1]) {
            double t = (std::log(x) - std::log(xs[i])) / (std::log(xs[i+1]) - std::log(xs[i]));
            return std::exp(std::log(ys[i]) + t * (std::log(ys[i+1]) - std::log(ys[i])));
        }
    }
    return ys[n-1];
}

// Build a filled band TGraph using a common mass grid, taking point-wise min/max.
// This ensures the central line is always inside the band.
TGraph* MakeBand(TGraph* g_lo, TGraph* g_hi, TGraph* g_cen) {
    // Union of all x-points
    std::vector<double> allx;
    for (int i = 0; i < g_lo->GetN();  ++i) { double x,y; g_lo->GetPoint(i,x,y);  allx.push_back(x); }
    for (int i = 0; i < g_hi->GetN();  ++i) { double x,y; g_hi->GetPoint(i,x,y);  allx.push_back(x); }
    for (int i = 0; i < g_cen->GetN(); ++i) { double x,y; g_cen->GetPoint(i,x,y); allx.push_back(x); }
    std::sort(allx.begin(), allx.end());
    allx.erase(std::unique(allx.begin(), allx.end()), allx.end());

    std::vector<double> lo_y, hi_y;
    for (double x : allx) {
        double ylo  = InterpGraph(g_lo,  x);
        double yhi  = InterpGraph(g_hi,  x);
        double ycen = InterpGraph(g_cen, x);
        lo_y.push_back(std::min({ylo, yhi, ycen}));
        hi_y.push_back(std::max({ylo, yhi, ycen}));
    }

    int n = allx.size();
    std::vector<double> bx, by;
    for (int i = 0; i < n; ++i)   { bx.push_back(allx[i]); by.push_back(lo_y[i]); }
    for (int i = n-1; i >= 0; --i) { bx.push_back(allx[i]); by.push_back(hi_y[i]); }
    bx.push_back(allx[0]); by.push_back(lo_y[0]);
    return new TGraph(bx.size(), bx.data(), by.data());
}

void plot_darkphoton_band() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const double Q_THR = 2.71;

    // Files: band edges (DC=1e-5 best, DC=1e-1 worst) and central line (DC=1e-3)
    const char* f_lo  = "outputs/darkphoton/hypmat_unscreened_ne/scan_dmelectron_pattern.root";       // DC=1e-5
    const char* f_c1  = "outputs/darkphoton/hypmat_unscreened_ne_dc1e4/scan_dmelectron_pattern.root"; // DC=1e-4
    const char* f_cen = "outputs/darkphoton/hypmat_unscreened_ne_dc1e3/scan_dmelectron_pattern.root"; // DC=1e-3 (central)
    const char* f_c2  = "outputs/darkphoton/hypmat_unscreened_ne_dc1e2/scan_dmelectron_pattern.root"; // DC=1e-2
    const char* f_hi  = "outputs/darkphoton/hypmat_unscreened_ne_dc1e1/scan_dmelectron_pattern.root"; // DC=1e-1

    TGraph* g_lo  = ExtractLimitCurve(f_lo,  Q_THR);
    TGraph* g_cen = ExtractLimitCurve(f_cen, Q_THR);
    TGraph* g_hi  = ExtractLimitCurve(f_hi,  Q_THR);

    if (!g_lo || !g_cen || !g_hi) { printf("ERROR: failed to extract limit curves\n"); return; }

    // --- Build band polygon ---
    TGraph* g_band = MakeBand(g_lo, g_hi, g_cen);

    // --- Styling ---
    const int band_color = kRed - 4;

    g_band->SetFillColor(band_color);
    g_band->SetFillStyle(3001);
    g_band->SetLineColor(0);  // no border line

    g_cen->SetLineColor(kRed + 1);
    g_cen->SetLineWidth(3);
    g_cen->SetLineStyle(2);  // dashed

    // --- Canvas ---
    TCanvas* c = new TCanvas("c_dp_band", "Dark photon band", 800, 650);
    c->SetLogx(); c->SetLogy();
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
    c->SetRightMargin(0.04); c->SetTopMargin(0.05);

    // Frame
    double xmin = 0.34, xmax = 25.0, ymin = 1e-17, ymax = 1e-12;
    TH2D* frame = new TH2D("frame","", 100, xmin, xmax, 100, ymin, ymax);
    frame->GetXaxis()->SetTitle("m_{A'} [eV]");
    frame->GetYaxis()->SetTitle("#varepsilon");
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetLabelSize(0.04);
    frame->GetYaxis()->SetLabelSize(0.04);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->Draw("AXIS");

    g_band->Draw("F SAME");
    g_cen->Draw("L SAME");

    // --- Legend ---
    TLegend* leg = new TLegend(0.38, 0.70, 0.96, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
    leg->AddEntry(g_cen,  "SrCd_{2}Sb_{2}, 1 kg-yr, DC = 10^{-3} e/pix/d", "L");
    leg->AddEntry(g_band, "Band: DC #in [10^{-5}, 10^{-1}] e/pix/d", "F");
    leg->Draw();

    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.04); lat.SetTextColor(kGray+2);
    lat.DrawLatex(0.55, 0.62, "HypMat (SrCd_{2}Sb_{2}), 1 kg-yr");

    c->RedrawAxis();
    c->SaveAs("outputs/darkphoton/hypmat_band_dc1e3_center.pdf");
    printf("Saved: outputs/darkphoton/hypmat_band_dc1e3_center.pdf\n");
}
