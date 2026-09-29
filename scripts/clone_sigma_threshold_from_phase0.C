// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  clone_sigma_threshold_from_phase0.C -- One-shot helper: I copy the
//  Phase-0 upper-limit TGraph into sigma_threshold.root so that the
//  threshold_toys mode can read it.
// ===========================================================================

// One-shot: copy Phase-0 UL TGraph into sigma_threshold.root for threshold_toys.
void clone_sigma_threshold_from_phase0(
    const char* phase0_root =
        "outputs/band_srdm_pattern_csv_per_mass_pydme_match/_toys/phase0/"
        "scan_srdm_pattern_csv.root",
    const char* out_root =
        "outputs/band_srdm_pattern_csv_per_mass_pydme_match/_toys/"
        "sigma_threshold.root") {
  TFile f0(phase0_root, "READ");
  if (!f0.IsOpen()) {
    Error("clone_sigma_threshold_from_phase0", "cannot open %s", phase0_root);
    return;
  }
  TGraph* g = dynamic_cast<TGraph*>(f0.Get("upper_limit_sigma_e_mchi_graph"));
  if (!g) {
    Error("clone_sigma_threshold_from_phase0",
          "missing upper_limit_sigma_e_mchi_graph in %s", phase0_root);
    return;
  }
  TFile fo(out_root, "RECREATE");
  TGraph* gc = dynamic_cast<TGraph*>(g->Clone("sigma_threshold_per_mass"));
  gc->SetTitle(";m_{#chi} [MeV];#sigma_{threshold} [cm^{2}]");
  gc->Write("sigma_threshold_per_mass");
  fo.Close();
  f0.Close();
  Info("clone_sigma_threshold_from_phase0", "wrote %s", out_root);
}
