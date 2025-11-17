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
      } else if (q_min >= q_thr) {
        std::cout << "[limit] mchi=" << mchi
                  << " MeV: excludes all scanned sigma (min q=" << q_min
                  << "). Using lowest sigma bin as limit.\n";
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
  // 4) Make a quick publication-style plot
  // ------------------------------------------------------------------
  gStyle->SetOptStat(0);

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

  glimit->GetXaxis()->SetTitle("m_{#chi} (MeV/c^{2})");
  glimit->GetYaxis()->SetTitle("#sigma_{e} (cm^{2})");

  glimit->GetXaxis()->SetTitleSize(0.05);
  glimit->GetYaxis()->SetTitleSize(0.05);
  glimit->GetXaxis()->SetLabelSize(0.045);
  glimit->GetYaxis()->SetLabelSize(0.045);

  // Adjust plot ranges nice-ish
  double xmin = mchi_vals.front();
  double xmax = mchi_vals.back();
  double ymin = *std::min_element(sigma_lim_vals.begin(), sigma_lim_vals.end());
  double ymax = *std::max_element(sigma_lim_vals.begin(), sigma_lim_vals.end());

  glimit->GetXaxis()->SetLimits(xmin * 0.8, xmax * 1.2);
  glimit->SetMinimum(ymin * 0.3);
  glimit->SetMaximum(ymax * 3.0);

  glimit->Draw("AL");

  // Simple legend, you can customise text later
  TLegend* leg = new TLegend(0.18,0.78,0.55,0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(glimit, "DAMIC-M, this work (QE-Dark chain)", "l");
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
