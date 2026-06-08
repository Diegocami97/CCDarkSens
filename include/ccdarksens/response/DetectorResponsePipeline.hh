// ============================================================================
//  CCDarkSens — DetectorResponsePipeline
//  Header for the unified pattern/PCD detector-response pipeline and E-dependent folding entry points.
//
//  Author: Diego Venegas-Vargas
// ============================================================================

#pragma once

#include <map>
#include <memory>
#include <vector>
#include <utility>

#include <TH1D.h>

#include "ccdarksens/response/ChargeIonization.hh"
#include "ccdarksens/response/Diffusion.hh"
#include "ccdarksens/response/PatternEfficiency.hh"
#include "ccdarksens/response/EfficiencyMC.hh"
#include "ccdarksens/response/PCDBasedResponse.hh"
#include "ccdarksens/response/PCDCalculator.hh"

namespace ccdarksens {

/// Analysis space selector: pattern-based or PCD-based.
enum class AnalysisSpace {
  Pattern,
  PCD
};

/**
 * DetectorResponsePipeline:
 *
 *   MASTER wrapper controlling the detector response applied to S(E)
 *   to obtain an observable-space distribution for the likelihood.
 *
 *   Two primary modes:
 *
 *   1) Pattern space (n_e histogram)
 *      - Ionization → (optionally PCD-based reconstruction) → PatternEfficiency
 *      - Produces TH1D in n_e (1..N)
 *
 *   2) PCD space (q histogram)
 *      - Ionization → PCDBasedResponse (MC P(q|n_e)) → PCDCalculator
 *      - Produces TH1D in q (continuous charge)
 *
 *   In addition, ApplyEDependent implements an E-dependent triple convolution
 *   used by your dedicated pattern-scan app.
 */
class DetectorResponsePipeline {
public:
  DetectorResponsePipeline(std::shared_ptr<ChargeIonization> ion)
  : ion_(std::move(ion)) {}

  //
  // Setter for analysis mode
  //
  void SetAnalysisSpace(AnalysisSpace space) {
    analysis_space_ = space;
  }

  //
  // Setters for pattern-mode components
  //
  void SetDiffusion(std::shared_ptr<Diffusion> diff) {
    diff_ = std::move(diff);
  }

  void SetPatternEfficiency(std::shared_ptr<PatternEfficiency> pe) {
    pe_ = std::move(pe);
  }

  /// When true, Apply (pattern mode) returns S_rec without applying ε(n_e).
  /// Use when folding to pattern bins so acceptance is applied only in the fold.
  void SetSkipPatternEfficiency(bool skip) { skip_pattern_efficiency_ = skip; }

  void SetEfficiencyMC(const std::shared_ptr<EfficiencyMC>& emc) { emc_ = emc; }

  //
  // Setters for PCD-mode components
  //
  void SetPCDBasedResponse(std::shared_ptr<PCDBasedResponse> pcdresp) {
    pcd_response_ = std::move(pcdresp);
  }

  void SetPCDCalculator(std::shared_ptr<PCDCalculator> pcdcalc) {
    pcd_calc_ = std::move(pcdcalc);
  }

  /// Configure BuildNeKernelFromPCD parameters (defaults match historical hard-coded values).
  void SetPCDKernelParams(double sigma_res_e, double Dqmin, double Dqmax) {
    pcd_sigma_res_e_ = sigma_res_e;
    pcd_Dqmin_       = Dqmin;
    pcd_Dqmax_       = Dqmax;
  }

  /**
   * Apply the chosen detector response to the input S(E).
   *
   * For Pattern mode:
   *   Returns a TH1D in n_e bins.
   *
   * For PCD mode:
   *   Returns a TH1D in q bins.
   *
   * ne_min/ne_max apply only to pattern mode.
   */
  std::unique_ptr<TH1D> Apply(const TH1D& dRdE,
                              double exposure_kg_year,
                              int ne_min, int ne_max,
                              double Ee_ref_eV = 0.0) const;

  /// Getter for analysis mode
  AnalysisSpace GetAnalysisSpace() const { return analysis_space_; }

  // Flag for apps that want to use the dedicated E-dependent logic.
  // (Your scan app is already calling ApplyEDependent directly.)
  bool use_Edependent_pattern_ = false;
  void EnableEDependentPattern(bool v) { use_Edependent_pattern_ = v; }

  /// E-dependent triple convolution used by your pattern-scan app.
  std::unique_ptr<TH1D> ApplyEDependent( TH1D& dRdE,
                                         double exposure_kg_year,
                                         int ne_min, int ne_max,
                                         const std::vector<double>& E_grid,
                                         const std::vector<std::vector<double>>& eps_Ene) const;

private:
  // Core required component for both modes:
  std::shared_ptr<ChargeIonization> ion_;

  // Pattern-mode components:
  std::shared_ptr<Diffusion>         diff_;
  std::shared_ptr<PatternEfficiency> pe_;
  std::shared_ptr<EfficiencyMC>      emc_;

  // PCD-mode components:
  std::shared_ptr<PCDBasedResponse>  pcd_response_;
  std::shared_ptr<PCDCalculator>     pcd_calc_;

  // Selected analysis space
  AnalysisSpace analysis_space_ = AnalysisSpace::Pattern;

  bool skip_pattern_efficiency_ = false;

  // PCD kernel parameters (used in BuildNeKernelFromPCD calls)
  double pcd_sigma_res_e_ = 0.21;
  double pcd_Dqmin_       = 0.5;
  double pcd_Dqmax_       = 0.5;
};

} // namespace ccdarksens
