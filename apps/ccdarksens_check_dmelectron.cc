#include "ccdarksens/io/ConfigManager.hh"
#include "ccdarksens/model/DMElectronModel.hh"
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <iostream>
#include <filesystem>

using namespace ccdarksens;
namespace fs = std::filesystem;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <config.json>\n";
    return 1;
  }

  ConfigManager cfg(argv[1]);
  cfg.parse();

  std::string outdir = cfg.run().outdir;
  if (outdir.empty())
    outdir = "outputs/" + (cfg.run().label.empty() ? std::string("unnamed")
                                                   : cfg.run().label);

  std::error_code ec;
  std::filesystem::create_directories(outdir, ec);
  if (ec) {
    std::cerr << "[error] failed to create outdir '" << outdir
              << "': " << ec.message() << "\n";
    return 1;
  }

  std::cout << "[run] label=" << cfg.run().label
            << " outdir='" << outdir << "'\n";


  const auto& m = cfg.model();                // new accessor we added
  DMElectronModel model;
  DMElectronConfig mc;
  mc.material          = m.material;
  mc.mediator          = m.mediator;
  mc.rates_dir         = m.rates_dir;
  mc.filename_template = m.filename_template;
  mc.mchi_MeV          = m.mchi_MeV;
  mc.sigma_e_cm2       = m.sigma_e_cm2;
  mc.Emin_eV           = m.Emin_eV;
  mc.Emax_eV           = m.Emax_eV;
  mc.nbins             = m.nbins;

  if (!model.Configure(mc)) {
    std::cerr << "[error] failed to load dRdE file\n";
    return 1;
  }

  std::unique_ptr<TH1D> h = model.MakeSpectrum_E();
  if (!h) { std::cerr << "[error] no histogram\n"; return 1; }

  h->Scale(cfg.run().exposure_kg_year);        // counts/eV for this exposure
  std::cout << "[check] Integral = " << h->Integral("width") << " counts\n";

  gStyle->SetOptStat(0);
  TCanvas c("c","dR/dE check",800,600);
  h->SetLineColor(kAzure+2);
  h->SetLineWidth(2);
  h->GetXaxis()->SetTitle("E_{e} [eV]");
  h->GetYaxis()->SetTitle("dR/dE [counts / eV]");
  h->Draw("hist");
  c.SaveAs((fs::path(cfg.run().outdir)/"dRdE_check.pdf").c_str());
//   c.SaveAs("(fs::path(cfg.run().outdir)"/"dRdE_check.pdf").c_str());

  return 0;
}
