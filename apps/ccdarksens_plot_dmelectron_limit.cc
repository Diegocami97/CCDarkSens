#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TAxis.h>
#include <fstream>
#include <sstream>
#include <algorithm>  // if not already included (for min_element/max_element)

TGraph* LoadExclusionCSV(const std::string& csv_path,
                         int line_color = kGray+2,
                         int line_style = 1,
                         int line_width = 2)
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
    g->SetPoint(i, xs[i], ys[i]);
  }

  g->SetLineColor(line_color);
  g->SetLineStyle(line_style);
  g->SetLineWidth(line_width);

  return g;
}


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


int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0]
              << " scan_dmelectron_grid.root [q_threshold]\n";
    return 1;
  }

  const std::string in_path = argv[1];
  double q_thr = 2.71; // default ~90% CL, 1 dof

  if (argc >= 3) {
    q_thr = std::stod(argv[2]);
  }

  std::cout << "[limit] Input file: " << in_path << "\n";
  std::cout << "[limit] q_threshold = " << q_thr << "\n";

  // ------------------------------------------------------------------
  // 1) Open file and get q(mchi, sigma) histogram
  // ------------------------------------------------------------------
  TFile* fin = TFile::Open(in_path.c_str(), "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "[limit] ERROR: cannot open file " << in_path << "\n";
    return 1;
  }

  TH2D* hq = dynamic_cast<TH2D*>(fin->Get("q_mchi_sigma"));
  if (!hq) {
    std::cerr << "[limit] ERROR: histogram 'q_mchi_sigma' not found.\n";
    fin->Close();
    return 1;
  }

  hq->SetDirectory(nullptr); // decouple from file
  fin->Close();

  const int Nx = hq->GetNbinsX();
  const int Ny = hq->GetNbinsY();

  std::cout << "[limit] q_mchi_sigma: Nx=" << Nx << "  Ny=" << Ny << "\n";

  // ------------------------------------------------------------------
  // 2) For each mchi-bin (x), find sigma where q crosses q_thr
  // ------------------------------------------------------------------
  std::vector<double> mchi_vals;
  std::vector<double> sigma_lim_vals;

  for (int ix = 1; ix <= Nx; ++ix) {
    const double mchi = hq->GetXaxis()->GetBinCenter(ix);

    // Extract q vs sigma for this column
    std::vector<double> sig(Ny), q(Ny);
    for (int iy = 1; iy <= Ny; ++iy) {
      sig[iy - 1] = hq->GetYaxis()->GetBinCenter(iy);
      q[iy - 1]   = hq->GetBinContent(ix, iy);
    }

    // Sanity check: (roughly) monotonic q vs sigma? Print warning if not.
    bool monotonic_warned = false;
    for (int iy = 1; iy < Ny; ++iy) {
      if (q[iy] + 0.05 * std::fabs(q[iy]) < q[iy - 1]) {
        if (!monotonic_warned) {
          std::cerr << "[limit] WARNING: q not monotonic in sigma for mchi="
                    << mchi << " (iy=" << iy << ")\n";
          monotonic_warned = true;
        }
      }
    }

    // Find first crossing q(i) < q_thr <= q(i+1)
    double sigma_lim = -1.0;
    for (int iy = 0; iy < Ny - 1; ++iy) {
      const double q1 = q[iy];
      const double q2 = q[iy + 1];

      if (q1 < q_thr && q2 >= q_thr && q2 > q1) {
        // Interpolate in log10(sigma) between iy and iy+1
        const double s1 = sig[iy];
        const double s2 = sig[iy + 1];
        const double logS1 = std::log10(s1);
        const double logS2 = std::log10(s2);

        const double t = (q_thr - q1) / (q2 - q1);
        const double logSlim = logS1 + t * (logS2 - logS1);
        sigma_lim = std::pow(10.0, logSlim);
        std::cout << "[XCheck] mchi=" << mchi
            << " crossing between sigma=" << s1 << " (q=" << q1 << ")"
            << " and sigma=" << s2 << " (q=" << q2 << ")\n"
            << "         -> sigma_lim ~ " << sigma_lim << "\n";
        break;
      }
    }

    if (sigma_lim > 0.0) {
      mchi_vals.push_back(mchi);
      sigma_lim_vals.push_back(sigma_lim);
    } else {
      // No crossing:
      // - if all q < q_thr -> not sensitive at this mchi in scanned range
      // - if all q >= q_thr -> excludes down to lowest sigma scanned
      const double q_min = *std::min_element(q.begin(), q.end());
      const double q_max = *std::max_element(q.begin(), q.end());

      if (q_max < q_thr) {
        std::cout << "[limit] mchi=" << mchi
                  << " MeV: no exclusion within scanned sigma range "
                     "(max q=" << q_max << ")\n";
        // XCheck: also we could print S_ROI and B_ROI if we decide to have them saved elsewhere,
        // or at least the sigma range:
        std::cout << "[XCheck]  sigma range = [" << sig.front()
                    << ", " << sig.back() << "]\n";
      } else if (q_min >= q_thr) {
            std::cout << "[limit] mchi=" << mchi
                        << " MeV: excludes all scanned sigma (min q=" << q_min
                        << "). Using lowest sigma bin as limit.\n";
            std::cout << "[XCheck]  sigma_min=" << sig.front()
                        << "  q(sigma_min)=" << q.front()
                        << "  q(sigma_max)=" << q.back() << "\n";

            mchi_vals.push_back(mchi);
            sigma_lim_vals.push_back(sig.front());
      } else {
        std::cout << "[limit] mchi=" << mchi
                  << " MeV: crossing ambiguous (q-range "
                  << q_min << "–" << q_max << "). Skipping.\n";
      }
    }
  }

  if (mchi_vals.empty()) {
    std::cerr << "[limit] ERROR: no valid limit points were found.\n";
    return 1;
  }

  // ------------------------------------------------------------------
  // 3) Build TGraph of the exclusion curve
  // ------------------------------------------------------------------
  const int Npoints = static_cast<int>(mchi_vals.size());
  TGraph* glimit = new TGraph(Npoints);
  glimit->SetName("g_dme_limit");
  glimit->SetTitle("");

  for (int i = 0; i < Npoints; ++i) {
    glimit->SetPoint(i, mchi_vals[i], sigma_lim_vals[i]);


  }

    // ------------------------------------------------------------------
  // 4) Load external exclusions / theory curves from CSV
  // ------------------------------------------------------------------
  // Adjust these paths to your local directory structure:
  TGraph* g_damic_2025   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/heavy_mediator/DAMIC-M_this_work_bold_black_.csv",
                                            kRed, 2, 3);

  TGraph* g_sensei       = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/heavy_mediator/SENSEI.csv",
                                            kBlue+1, 2, 3);

  TGraph* g_supercdms    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/heavy_mediator/SuperCDMS.csv",
                                            kMagenta+2, 2, 3);

  TGraph* g_darkside50   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/heavy_mediator/DarkSide50.csv",
                                            kGreen+2, 2, 3);

  TGraph* g_panda4T      = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/heavy_mediator/Panda4T.csv",
                                            kOrange+7, 2, 3);

  TGraph* g_damicm_mike      = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/heavy_mediator/damic-m_1kgyear_mike.csv",
                                            kCyan, 2, 3);


  // Theory targets: freeze-in / freeze-out
  TGraph* g_model    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/heavy_mediator/freezeout.csv",
                                            kGray+2, 7, 3);  // dotted-ish



//   TGraph* g_damic_2025   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/DAMIC-M_this_work_bold_black_.csv",
//                                             kRed, 2, 3);

//   TGraph* g_sensei       = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/SENSEI.csv",
//                                             kBlue+1, 2, 3);

//   TGraph* g_supercdms    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/SuperCDMS.csv",
//                                             kMagenta+2, 2, 3);

//   TGraph* g_darkside50   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/DarkSide50.csv",
//                                             kGreen+2, 2, 3);

//   TGraph* g_panda4T      = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/Panda4T.csv",
//                                             kOrange+7, 2, 3);

//   TGraph* g_damicm_mike      = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/damic-m_1kgyear_mike.csv",
//                                             kCyan, 2, 3);

// //   // Theory targets: freeze-in / freeze-out
// //   TGraph* g_freeze_in    = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/dme_ulm_freeze_in.csv",
// //                                             kGray+2, 3, 3);  // dotted-ish

//   TGraph* g_model   = LoadExclusionCSV("/Users/diegovenegasvargas/Documents/CCDarkSens/data/previous_limits/light_mediator/freezein.csv",
//                                             kGray+2, 7, 3);  // another style




    

  // Expand ranges using all non-null graphs
  double xmin = *std::min_element(mchi_vals.begin(), mchi_vals.end());
  double xmax = *std::max_element(mchi_vals.begin(), mchi_vals.end());
  double ymin = *std::min_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
  double ymax = *std::max_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damic_2025);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_sensei);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_supercdms);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_darkside50);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_panda4T);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_sensei_2025);
//   UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_freeze_in);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_model);
  UpdateRangesFromGraph(xmin, xmax, ymin, ymax, g_damicm_mike);
//   xmax = 10e3; // limit x-axis to 1 GeV for better visualization


  // ------------------------------------------------------------------
  // 5) Make a quick publication-style plot
  // ------------------------------------------------------------------
  gStyle->SetOptStat(0);

//   TCanvas* c = new TCanvas("c_limit", "DM-e limit", 900, 700);
//   c->SetLogx();
//   c->SetLogy();
//   c->SetTicks(1,1);
//   c->SetLeftMargin(0.14);
//   c->SetRightMargin(0.04);
//   c->SetBottomMargin(0.12);
//   c->SetTopMargin(0.06);

//   glimit->SetLineWidth(3);
//   glimit->SetLineColor(kBlack);
//   glimit->SetLineStyle(1);

//   glimit->GetXaxis()->SetTitle("m_{#chi} (MeV/c^{2})");
//   glimit->GetYaxis()->SetTitle("#sigma_{e} (cm^{2})");

//   glimit->GetXaxis()->SetTitleSize(0.05);
//   glimit->GetYaxis()->SetTitleSize(0.05);
//   glimit->GetXaxis()->SetLabelSize(0.045);
//   glimit->GetYaxis()->SetLabelSize(0.045);

//   // Adjust plot ranges nice-ish
//   double xmin = mchi_vals.front();
//   double xmax = mchi_vals.back();
//   double ymin = *std::min_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
//   double ymax = *std::max_element(sigma_lim_vals.begin(), sigma_lim_vals.end());

//   glimit->GetXaxis()->SetLimits(xmin * 0.8, xmax * 1.2);
//   glimit->SetMinimum(ymin * 0.3);
//   glimit->SetMaximum(ymax * 3.0);

//   glimit->Draw("AL");

//   // Simple legend, you can customise text later
//   TLegend* leg = new TLegend(0.18,0.78,0.55,0.90);
//   leg->SetBorderSize(0);
//   leg->SetFillStyle(0);
//   leg->SetTextSize(0.04);
//   leg->AddEntry(glimit, "DAMIC-M, this work (QE-Dark chain), Ultra-Light Mediator", "l");
//   leg->Draw();

//   c->Update();

  TCanvas* c = new TCanvas("c_limit", "DM-e limit", 900, 700);
  c->SetLogx();
  c->SetLogy();
  c->SetTicks(1,1);
  c->SetLeftMargin(0.14);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.06);

  glimit->SetLineWidth(3);
  glimit->SetLineColor(kBlack);
  glimit->SetLineStyle(1);

  glimit->GetXaxis()->SetTitle("m_{#chi} [MeV/c^{2}]");
  glimit->GetYaxis()->SetTitle("#sigma_{e} [cm^{2}]");
  glimit->GetXaxis()->SetTitleSize(0.05);
  glimit->GetYaxis()->SetTitleSize(0.05);
  glimit->GetXaxis()->SetLabelSize(0.04);
  glimit->GetYaxis()->SetLabelSize(0.04);

  // ------------------------------------------------------------------
  // 4) Axis ranges from BOTH our limit curve AND external exclusions
  // ------------------------------------------------------------------
//   double xmin = *std::min_element(mchi_vals.begin(), mchi_vals.end());
//   double xmax = *std::max_element(mchi_vals.begin(), mchi_vals.end());
//   double ymin = *std::min_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
//   double ymax = *std::max_element(sigma_lim_vals.begin(), sigma_lim_vals.end());

  // Load external exclusion(s) from CSV
  // IMPORTANT: use the correct path where your CSV actually lives
  // e.g. "./SRDM_XENON1T-s20_ulightmediator.csv"

  glimit->GetXaxis()->SetLimits(xmin * 0.8, xmax );
  glimit->SetMinimum(ymin * 0.3);
  glimit->SetMaximum(ymax * 3.0);

  // Draw our limit first (defines axes)
  glimit->Draw("AL");

  // Overlay external exclusion(s)
  // Overlay external exclusions
  if (g_damic_2025)  g_damic_2025->Draw("L SAME");
//   if (g_sensei)      g_sensei->Draw("L SAME");
//   if (g_supercdms)   g_supercdms->Draw("L SAME");
//   if (g_darkside50)  g_darkside50->Draw("L SAME");
//   if (g_panda4T)     g_panda4T->Draw("L SAME");
  if (g_damicm_mike) g_damicm_mike->Draw("L SAME");
  // Overlay theory curves
//   if (g_freeze_in)   g_freeze_in->Draw("L SAME");
  if (g_model)  g_model->Draw("L SAME");

  TLegend* leg = new TLegend(0.65, 0.55, 0.92, 0.92); // move it to avoid overlapping curves
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);

  leg->AddEntry(glimit,      "DAMIC-M, 1 kg-year", "l");

  if (g_damic_2025)  leg->AddEntry(g_damic_2025,  "DAMIC-M (2025)",          "l");
//   if (g_sensei)      leg->AddEntry(g_sensei,      "SENSEI",                  "l");
//   if (g_supercdms)   leg->AddEntry(g_supercdms,   "SuperCDMS",               "l");
//   if (g_darkside50)  leg->AddEntry(g_darkside50,  "DarkSide-50",             "l");
//   if (g_panda4T)     leg->AddEntry(g_panda4T,     "PandaX-4T",               "l");
  if (g_damicm_mike) leg->AddEntry(g_damicm_mike, "DAMIC-M (1 kg-year, Mike)", "l");

//   if (g_freeze_in)   leg->AddEntry(g_freeze_in,   "Freeze-in target",        "l");
//   if (g_model)  leg->AddEntry(g_model,  "Freeze-in target",       "l");
  if (g_model)  leg->AddEntry(g_model,  "Freeze-out target",       "l");

  leg->AddEntry((TObject*)0, "#bf{F_{DM} = 1}", "");
//   leg->AddEntry((TObject*)0, "#bf{F_{DM} = ( #alpha  m_{e} / q )^{2}}", "");

  leg->Draw();
  c->Update();


  // Save a couple of outputs
  std::string out_pdf = "dme_limit_curve.pdf";
  std::string out_root= "dme_limit_curve.root";
  c->SaveAs(out_pdf.c_str());

  TFile fout(out_root.c_str(), "RECREATE");
  glimit->Write();
  hq->Write(); // keep the map too
  fout.Close();

  std::cout << "[limit] Saved: " << out_pdf  << "\n"
            << "[limit]         " << out_root << "\n";

  return 0;
}
