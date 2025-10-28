#pragma once
#include <cstddef>

class TH1D;

namespace ccdarksens {

// Parametric diffusion model + fast Gaussian rebinning in n_e space.
// sigma_xy(z,E) = sqrt(-A log(1 - b z)) * (alpha + beta * E_keV)
//
// Parameters:
//  - A_um2, b_umInv, alpha, beta_per_keV: match the Python model
//  - thickness_um: CCD thickness (for z-averaging sigma_xy)
//  - kappa_e_per_um: maps lateral sigma_xy [um] -> sigma_e [electrons]
//  - sigma_readout_e: extra Gaussian readout noise in electrons (quadrature)
class Diffusion {
public:
  Diffusion(double A_um2       = 803.25,
            double b_umInv     = 6.5e-4,
            double alpha       = 1.0,
            double beta_per_keV= 0.0,
            double thickness_um= 670.0,     // example: 0.67 mm
            double kappa_e_per_um = 0.0,    // set >0 to enable geometry->e- mapping
            double sigma_readout_e = 0.0)   // additional noise
  : A_um2_(A_um2), b_umInv_(b_umInv), alpha_(alpha), beta_per_keV_(beta_per_keV),
    thickness_um_(thickness_um), kappa_e_per_um_(kappa_e_per_um),
    sigma_readout_e_(sigma_readout_e) {}

  // Compute <sigma_e> by z-averaging sigma_xy over [0, thickness_um]
  // Ee_eV is the characteristic electron-recoil energy (use 0 if beta=0).
  double ComputeSigmaElectrons(double Ee_eV, std::size_t nz_steps = 512) const;

  // Convolve integer-binned TH1D (n_e bins centered at integers) with a
  // Gaussian whose mean is each source bin's n_e and stddev = ComputeSigmaElectrons(Ee_eV).
  // Uses CDF between [q-0.5, q+0.5] so probability mass is conserved.
  void Apply(TH1D& h_ne, double Ee_eV) const;

  // Accessors
  void set_kappa_e_per_um(double k)     { kappa_e_per_um_ = k; }
  void set_sigma_readout_e(double s)    { sigma_readout_e_ = s; }
  void set_thickness_um(double t)       { thickness_um_ = t; }

  double A_um2()        const { return A_um2_; }
  double b_umInv()      const { return b_umInv_; }
  double alpha()        const { return alpha_; }
  double beta_per_keV() const { return beta_per_keV_; }

private:
  // sigma_xy(z, Ee) in micrometers (Ee in eV; beta is per keV)
  double sigma_xy_um_(double z_um, double Ee_eV) const;

  static double normal_cdf_(double x); // Phi(x)
  static double interval_prob_(double mu, double sigma, double a, double b);

  double A_um2_;
  double b_umInv_;
  double alpha_;
  double beta_per_keV_;
  double thickness_um_;
  double kappa_e_per_um_;
  double sigma_readout_e_;
};

} // namespace ccdarksens
