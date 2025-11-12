// Run with:
//   root -l -q 'plot_p100K_table.C("data/p100K_table.csv")'
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TColor.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cctype>

static bool is_comment_or_empty_line(const std::string& s) {
  for (char c : s) { if (!std::isspace((unsigned char)c)) return c=='#'; }
  return true;
}

// --- custom 20-color palette (based on matplotlib "tab20") ---
static std::vector<int> make_custom_colors() {
  std::vector<std::array<double,3>> rgb = {
    {31,119,180}, {255,127,14}, {44,160,44}, {214,39,40}, {148,103,189},
    {140,86,75},  {227,119,194}, {127,127,127}, {188,189,34}, {23,190,207},
    {174,199,232}, {255,187,120}, {152,223,138}, {255,152,150}, {197,176,213},
    {196,156,148}, {247,182,210}, {199,199,199}, {219,219,141}, {158,218,229}
  };
  std::vector<int> cols;
  cols.reserve(rgb.size());
  for (size_t i=0;i<rgb.size();++i) {
    // cols.push_back(TColor::GetColor(rgb[i][0]/255., rgb[i][1]/255., rgb[i][2]/255.));
    cols.push_back(
    TColor::GetColor(
        (Float_t)(rgb[i][0]/255.0),
        (Float_t)(rgb[i][1]/255.0),
        (Float_t)(rgb[i][2]/255.0)
    )
    );

  }
  return cols;
}

void plot_p100K_table(const char* path = "data/p100K_table.csv") {
  std::ifstream in(path);
  if (!in) { Error("plot_p100K_table","Cannot open %s", path); return; }

  std::vector<double> E;
  std::vector<std::vector<double>> P;
  std::string line; int ncols=-1;

  while (std::getline(in,line)) {
    if (is_comment_or_empty_line(line)) continue;
    std::istringstream ss(line);
    std::string cell;
    std::vector<double> vals;
    while (std::getline(ss, cell, ',')) {
      if (!cell.empty() && (cell.back()=='\r'||cell.back()=='\n')) cell.pop_back();
      vals.push_back(cell.empty()?0.0:std::stod(cell));
    }
    if (vals.size()<2) continue;
    if (ncols<0) ncols=(int)vals.size();
    if ((int)vals.size()!=ncols) continue;
    E.push_back(vals[0]);
    P.emplace_back(vals.begin()+1, vals.end());
  }
  if (E.size()<2) { Error("plot_p100K_table","No valid rows in %s", path); return; }

  const int N=(int)E.size();
  const int K=(int)P.front().size();
  const int Kshow=std::min(20,K);

  // --- mean n_e(E) ---
  std::vector<double> mean_ne(N,0.0);
  for (int i=0;i<N;++i){
    double s=0.0; for (int j=0;j<K;++j) s+=(j+1)*P[i][j];
    mean_ne[i]=s;
  }

  // --- style ---
  gStyle->SetOptStat(0);
  auto colors = make_custom_colors();

  // --- probability graphs ---
  std::vector<TGraph*> prob_graphs; prob_graphs.reserve(K);
  for (int j=0;j<K;++j){
    auto* g=new TGraph(N);
    for (int i=0;i<N;++i) g->SetPoint(i,E[i],P[i][j]);
    g->SetLineWidth(2);
    if (j<Kshow){
      g->SetLineColor(colors[j % colors.size()]);
      g->SetMarkerColor(colors[j % colors.size()]);
    } else {
      g->SetLineColor(kGray+2);
      g->SetLineStyle(7);
    }
    prob_graphs.push_back(g);
  }

  auto* mg=new TMultiGraph();
  for(auto* g:prob_graphs) mg->Add(g,"L");

  auto* c1=new TCanvas("c1","P(ne|E)",1100,750);
  c1->cd();
//   gPad->SetLogx(true);
//   gPad->SetLogy(true);
  mg->Draw("A");
  mg->GetXaxis()->SetTitle("Recoil energy E_{r} [eV]");
  mg->GetYaxis()->SetTitle("P(n_{e} | E_{r})");
  mg->SetTitle("Ionization probability P(n_{e}|E_{r})");

  auto* leg=new TLegend(0.63,0.43,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (int j=0;j<Kshow;++j)
    leg->AddEntry(prob_graphs[j],Form("n_{e} = %d",j+1),"l");
  if (K>Kshow) leg->AddEntry(prob_graphs[Kshow],"others","l");
  leg->Draw();

  // --- mean ne ---
  auto* gMean=new TGraph(N);
  for(int i=0;i<N;++i) gMean->SetPoint(i,E[i],mean_ne[i]);
  gMean->SetLineColor(kMagenta+1);
  gMean->SetLineWidth(3);

  auto* c2=new TCanvas("c2","<ne>(E)",900,600);
  c2->cd();
  gPad->SetLogx(true);
  gPad->SetLogy(true);
  gMean->Draw("AL");
  gMean->GetXaxis()->SetTitle("Recoil energy E_{r} [eV]");
  gMean->GetYaxis()->SetTitle("#LT n_{e} #GT");
  gMean->SetTitle("Mean electron yield #LT n_{e} #GT vs E_{r}");

  c1->SaveAs("outplots/p100K_Pne_vs_E.pdf");
  c2->SaveAs("outplots/p100K_mean_ne_vs_E.pdf");


  // --- write to ROOT only ---
  TFile fout("outputs/p100K_table_plots.root","RECREATE");
  c1->Write();
  c2->Write();
  mg->Write("Pne_multigraph");
  for(int j=0;j<K;++j) prob_graphs[j]->Write(Form("Pne_ne%d",j+1));
  gMean->Write("mean_ne");
  fout.Close();

  printf("[ok] Wrote all ROOT objects to outputs/p100K_table_plots.root\n");
}
