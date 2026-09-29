// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  DetectorResponsePipeline.cc -- Master response orchestrator that maps
//  dR/dE to pattern-space S_obs(n_e) or PCD-space S_obs(q), including
//  optional PCD reconstruction kernels.
// ===========================================================================

#include "ccdarksens/response/DetectorResponsePipeline.hh"

#include <algorithm>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <TH1D.h>

namespace ccdarksens {

/**
 * MASTER detector-response application.
 *
 * Pattern mode:
 *   - Ionization:       S(E) → S(n_e)
 *   - If PCD machinery present:
 *        build P(q | n_e_true) → P(n_obs | n_true) → S_rec(n_obs)
 *     else:
 *        optional Diffusion in n_e
 *   - PatternEfficiency: ε(n_e_obs)
 *   - Result:           S_obs(n_e_obs)
 *
 * PCD mode:
 *   - Ionization:        S(E) → S(n_e)
 *   - PCDBasedResponse:  build P(q | n_e)
 *   - PCDCalculator:     fold S(n_e) with P(q | n_e)
 *   - Result:            S_obs(q)
 */
std::unique_ptr<TH1D>
DetectorResponsePipeline::Apply(const TH1D& dRdE,
                                double exposure_kg_year,
                                int ne_min, int ne_max,
                                double Ee_ref_eV) const
{
  if (!ion_) {
    throw std::runtime_error("DetectorResponsePipeline::Apply: ChargeIonization not set.");
  }

  //
  // ────────────────────────────────────────────────────────────────
  // PATTERN MODE
  // ────────────────────────────────────────────────────────────────
  //
  if (analysis_space_ == AnalysisSpace::Pattern) {

    //
    // 1) dR/dE → S_true(n_e) expected true electron counts
    //
    auto h_ne_true = ion_->FoldToNe(dRdE, exposure_kg_year, ne_min, ne_max);

    //
    // 2) If PCD machinery is available, use it to build a
    //    reconstruction kernel P(n_obs | n_true) from P(q | n_true),
    //    then fold S_true → S_rec. Otherwise, fall back to the old
    //    behavior (just S_true, with optional histogram-level diffusion).
    //
    std::unique_ptr<TH1D> h_ne_rec;

    const bool have_pcd_kernel = (pcd_response_ && pcd_calc_);

    // When PCD is used, diffusion is already inside P(q|n_e) from the full
    // detector MC (ChargeTransport, PixelSimulator). We never call diff_->Apply()
    // here, so diffusion is applied only once. When PCD is not used, we apply
    // histogram-level diff_ in the else branch only.
    if (have_pcd_kernel) {
      //
      // 2a) Build P(q | n_e_true) with full detector MC (includes diffusion)
      //
      const auto& pcd_table =
        pcd_response_->BuildPCDTable(ne_min, ne_max, Ee_ref_eV);

      //
      // 2b) Build P(n_obs | n_true) from the PCD table
      //
      auto kernel_ne =
        pcd_calc_->BuildNeKernelFromPCD(pcd_table, ne_min, ne_max,
                                        pcd_sigma_res_e_, pcd_Dqmin_, pcd_Dqmax_);

      //
      // 2c) Fold S_true(n_true) → S_rec(n_obs)
      //
      h_ne_rec = pcd_calc_->FoldNeSpectrum(*h_ne_true, kernel_ne,
                                           ne_min, ne_max, "S_rec_ne");

      // NOTE: we do NOT apply additional n_e-level diffusion here,
      // because the PCD-based kernel already encodes detector smearing
      // (diffusion, DC, readout noise, etc.).
    } else {
      // Optional: apply histogram-level diffusion in n_e space
      h_ne_rec = std::make_unique<TH1D>(*h_ne_true); // copy S_true

      if (diff_) {
        diff_->Apply(*h_ne_rec, Ee_ref_eV);
      }
    }

    //
    // 3) Apply pattern efficiency ε(n_e_obs) unless skipped (pattern-space fold applies it once)
    //
    if (!pe_ && !skip_pattern_efficiency_) {
      throw std::runtime_error(
        "DetectorResponsePipeline::Apply (Pattern): PatternEfficiency not set.");
    }
    if (pe_ && !skip_pattern_efficiency_) {
      pe_->Apply(*h_ne_rec);
    }

    //
    // Done — h_ne_rec is now the observable S_obs(n_e_obs)
    //
    return h_ne_rec;
  }

  //
  // ────────────────────────────────────────────────────────────────
  // PCD MODE
  // ────────────────────────────────────────────────────────────────
  //
  else if (analysis_space_ == AnalysisSpace::PCD) {

    //
    // 1) dR/dE → S(n_e)
    //
    if (ne_min < 0 || ne_max < ne_min) {
      throw std::runtime_error(
        "DetectorResponsePipeline::Apply (PCD): invalid n_e range.");
    }

    auto h_ne = ion_->FoldToNe(dRdE, exposure_kg_year,
                               ne_min, ne_max);

    //
    // Convert S(n_e) histogram → map<int,double>
    //
    std::map<int, double> S_map;
    for (int ne = ne_min; ne <= ne_max; ++ne) {
      int bin = h_ne->FindBin(ne);
      S_map[ne] = h_ne->GetBinContent(bin);
    }

    //
    // 2) Build P(q | n_e)
    //
    if (!pcd_response_) {
      throw std::runtime_error(
        "DetectorResponsePipeline::Apply (PCD): PCDBasedResponse not set.");
    }

    const auto& pcd_table =
      pcd_response_->BuildPCDTable(ne_min, ne_max, Ee_ref_eV);

    //
    // 3) Fold S(n_e) with P(q | n_e)
    //
    if (!pcd_calc_) {
      throw std::runtime_error(
        "DetectorResponsePipeline::Apply (PCD): PCDCalculator not set.");
    }

    // Background folding is not handled here. This function only returns
    // the SIGNAL observable. Background folding will be done by BackgroundBuilder.
    std::map<int, double> B_empty;

    auto folded_pair = pcd_calc_->FoldSpectra(pcd_table, S_map, B_empty);

    //
    // Return the SIGNAL q-spectrum
    //
    return std::move(folded_pair.first);
  }

  //
  // ────────────────────────────────────────────────────────────────
  // SHOULD NOT REACH HERE
  // ────────────────────────────────────────────────────────────────
  //
  throw std::runtime_error(
      "DetectorResponsePipeline::Apply: Unknown analysis space.");
}

// -----------------------------------------------------------------------------
// ApplyEDependent: full triple convolution in energy, n_true, and reconstruction,
// then pattern efficiency ε(E, n_obs).
// -----------------------------------------------------------------------------
std::unique_ptr<TH1D>
DetectorResponsePipeline::ApplyEDependent(TH1D& dRdE,
                                          double exposure_kg_year,
                                          int ne_min, int ne_max,
                                          const std::vector<double>& E_grid,
                                          const std::vector<std::vector<double>>& eps_Ene) const
{
    if (!ion_)
        throw std::runtime_error("ApplyEDependent: missing ChargeIonization.");

    if (!pcd_response_ || !pcd_calc_)
        throw std::runtime_error("ApplyEDependent: missing PCD machinery.");

    if (!emc_)
        throw std::runtime_error("ApplyEDependent: missing EfficiencyMC (epsilon).");

    const int Nn = ne_max - ne_min + 1;  // number of n_e bins

    auto h = std::make_unique<TH1D>("S_obs_ne",
                                    "Observed S(n_e);n_e;counts",
                                    Nn, ne_min - 0.5, ne_max + 0.5);
    h->Sumw2();

    // Build P(q | n_true) once (assumed E–independent)
    const auto& pcd_table =
        pcd_response_->BuildPCDTable(ne_min, ne_max, /*dummy E*/ 0.0);

    // Build P(n_obs | n_true) once
    auto kernel =
        pcd_calc_->BuildNeKernelFromPCD(pcd_table, ne_min, ne_max,
                                        pcd_sigma_res_e_, pcd_Dqmin_, pcd_Dqmax_);

    // Helper lambda: interpolate eps_Ene(E) in energy
    auto interpolate_eps_row =
        [&](double E, std::vector<double>& eps_row_out)
    {
        eps_row_out.assign(Nn, 0.0);

        // If E is below the first grid point, use the first row
        if (E <= E_grid.front()) {
            eps_row_out = eps_Ene.front();
            return;
        }
        // If E is above the last grid point, use the last row
        if (E >= E_grid.back()) {
            eps_row_out = eps_Ene.back();
            return;
        }

        // Find bracketing indices i_lo, i_hi such that
        // E_grid[i_lo] <= E < E_grid[i_hi]
        auto it_hi = std::upper_bound(E_grid.begin(), E_grid.end(), E);
        const int i_hi = static_cast<int>(std::distance(E_grid.begin(), it_hi));
        const int i_lo = i_hi - 1;

        const double E_lo = E_grid[i_lo];
        const double E_hi = E_grid[i_hi];
        const double t = (E - E_lo) / (E_hi - E_lo);

        const auto& row_lo = eps_Ene[i_lo];
        const auto& row_hi = eps_Ene[i_hi];

        for (int io = 0; io < Nn; ++io) {
            eps_row_out[io] = (1.0 - t) * row_lo[io] + t * row_hi[io];
        }
    };

    // MAIN LOOP OVER dR/dE BINS
    const int nbins = dRdE.GetNbinsX();
    std::vector<double> eps_row(Nn, 0.0);  // epsilon(n_obs) interpolated at the current energy

    for (int ib = 1; ib <= nbins; ++ib) {
        const double E         = dRdE.GetBinCenter(ib);
        const double dRdE_val  = dRdE.GetBinContent(ib);
        const double dE        = dRdE.GetBinWidth(ib);

        double weight = dRdE_val * exposure_kg_year * dE;
        if (weight <= 0.0) continue;

        // 1) P(n_true | E)
        auto Pntrue = ion_->ProbNeGivenE(E, ne_min, ne_max);

        // 2) epsilon(n_obs, E) via interpolation on precomputed table
        interpolate_eps_row(E, eps_row);

        // 3) Triple convolution
        for (int n_true = ne_min; n_true <= ne_max; ++n_true) {
            const int it = n_true - ne_min;
            const double ptrue = Pntrue[it];
            if (ptrue <= 0.0) continue;

            for (int n_obs = ne_min; n_obs <= ne_max; ++n_obs) {
                const int io = n_obs - ne_min;
                const double P_rec = kernel[it][io];
                if (P_rec <= 0.0) continue;

                const double eps = eps_row[io];
                if (eps <= 0.0) continue;

                const double contrib = weight * ptrue * P_rec * eps;
                h->AddBinContent(h->FindBin(n_obs), contrib);
            }
        }
    }

    return h;
}

} // namespace ccdarksens
