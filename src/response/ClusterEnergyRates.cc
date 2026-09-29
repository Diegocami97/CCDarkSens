// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  ClusterEnergyRates.cc -- Folds a dR/dE_true density spectrum through a
//  ClusterFitMC kernel into per-E_reco-bin rates.
// ===========================================================================

#include "ccdarksens/response/ClusterEnergyRates.hh"
#include <TH1D.h>

#include <algorithm>
#include <cmath>

namespace ccdarksens {

namespace {

// Geometric-midpoint edges for a (typically log-spaced) point grid: edge[i]
// is the boundary between point i-1's and point i's assigned neighborhood.
// The two outer edges extend by the same log-ratio as the adjacent step, so
// every point (including the first/last) gets a neighborhood of nonzero
// width. Falls back to the arithmetic midpoint for any pair of points that
// aren't both positive (geometric mean undefined there).
// Returns n+1 edges around n grid points (details in the comment just above).
std::vector<double> NeighborhoodEdges(const std::vector<double>& points) {
  const std::size_t n = points.size();
  std::vector<double> edges(n + 1, 0.0);
  if (n == 0) return edges;
  if (n == 1) {
    const double w = std::abs(points[0]) > 0.0 ? std::abs(points[0]) * 0.5 : 0.5;
    edges[0] = points[0] - w;
    edges[1] = points[0] + w;
    return edges;
  }

  auto midpoint = [](double a, double b) {
    return (a > 0.0 && b > 0.0) ? std::sqrt(a * b) : 0.5 * (a + b);
  };

  for (std::size_t i = 1; i < n; ++i) edges[i] = midpoint(points[i - 1], points[i]);
  // Extend the outer edges by the same ratio/offset as the adjacent interior gap.
  edges[0] = (points[0] > 0.0 && edges[1] > 0.0 && points[0] != edges[1])
                 ? points[0] * (points[0] / edges[1])
                 : points[0] - (edges[1] - points[0]);
  edges[n] = (points[n - 1] > 0.0 && edges[n - 1] > 0.0 && points[n - 1] != edges[n - 1])
                 ? points[n - 1] * (points[n - 1] / edges[n - 1])
                 : points[n - 1] + (points[n - 1] - edges[n - 1]);
  edges[0] = std::max(edges[0], 0.0);
  return edges;
}

}  // namespace

// ----------------------------------------------------------------------------
// FoldEtrueToErecoRates
//   I fold dR/dE_true through the kernel K[E_true][E_reco]. For every kernel
//   point I integrate the input density over that point's own neighborhood
//   (bins between the geometric midpoints to its neighbors, the shared boundary
//   bin going to the upper point), multiply by the exposure, and spread the
//   counts over E_reco with the kernel row. See ClusterEnergyRates.hh.
// ----------------------------------------------------------------------------
std::vector<double> FoldEtrueToErecoRates(
    const TH1D& h_Etrue_density,
    double exposure_kg_year,
    const KernelMatrix& K) {
  const std::size_t n_ereco = K.Ereco_edges_eV.empty() ? 0 : K.Ereco_edges_eV.size() - 1;
  std::vector<double> rate_reco(n_ereco, 0.0);
  if (K.Etrue_grid_eV.empty()) return rate_reco;

  TH1D& h = const_cast<TH1D&>(h_Etrue_density);
  const auto edges = NeighborhoodEdges(K.Etrue_grid_eV);
  const int nbins_h = h.GetNbinsX();

  for (std::size_t iE = 0; iE < K.Etrue_grid_eV.size(); ++iE) {
    // Integrate the density over this point's full assigned neighborhood
    // [edges[iE], edges[iE+1]) -- not just the single fine bin nearest the
    // point -- so every part of the input spectrum is counted exactly once
    // across the sparse kernel grid (see header docstring).
    int bin_lo = h.FindBin(edges[iE]);
    int bin_hi = h.FindBin(edges[iE + 1]);
    // Exclude the boundary bin from this point's neighborhood when it isn't
    // the last point, so the next point's [bin_lo] picks it up instead --
    // avoids double-counting the shared boundary bin between neighbors.
    if (iE + 1 < K.Etrue_grid_eV.size() && bin_hi > bin_lo) bin_hi -= 1;
    bin_lo = std::max(bin_lo, 1);
    bin_hi = std::min(bin_hi, nbins_h);

    const double counts_i =
        (bin_hi >= bin_lo) ? h.Integral(bin_lo, bin_hi, "width") * exposure_kg_year : 0.0;

    const auto& row = K.K[iE];
    for (std::size_t j = 0; j < n_ereco; ++j) {
      rate_reco[j] += counts_i * row[j];
    }
  }
  return rate_reco;
}

}  // namespace ccdarksens
