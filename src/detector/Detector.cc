// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  Detector.cc -- I validate the detector geometry/material and compute the
//  active target mass from the geometry (or return the mass override from
//  the config).
// ===========================================================================

#include "ccdarksens/detector/Detector.hh"
#include <cmath>
#include <stdexcept>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// Detector::Detector
//   I store the geometry, material and optional mass override, and reject
//   anything unphysical (non-positive rows/cols/pitch/thickness/density,
//   active_fraction outside (0,1], empty element name) with
//   std::invalid_argument.
// ----------------------------------------------------------------------------
Detector::Detector(DetectorGeometry geom, TargetMaterial mat, std::optional<double> mass_override)
  : geom_(std::move(geom)), mat_(std::move(mat)), mass_override_kg_(mass_override)
{
  if (geom_.rows <= 0 || geom_.cols <= 0)
    throw std::invalid_argument("rows/cols must be > 0");
  if (geom_.pixel_size_um <= 0.0)
    throw std::invalid_argument("pixel_size_um must be > 0");
  if (geom_.thickness_mm <= 0.0)
    throw std::invalid_argument("thickness_mm must be > 0");
  if (geom_.active_fraction <= 0.0 || geom_.active_fraction > 1.0)
    throw std::invalid_argument("active_fraction must be (0,1]");
  if (mat_.element.empty())
    throw std::invalid_argument("target_element required");
  if (mat_.density_g_cm3 <= 0.0)
    throw std::invalid_argument("density_g_cm3 > 0");
}

// ----------------------------------------------------------------------------
// Detector::compute_mass_from_geometry_kg
//   mass = density * (cols*pitch) * (rows*pitch) * thickness * active_fraction,
//   with the pitch converted from um to cm and the thickness from mm to cm.
//   Returns kg.
// ----------------------------------------------------------------------------
double Detector::compute_mass_from_geometry_kg() const {
  const double pix_cm      = geom_.pixel_size_um * 1e-4; // microns to cm
  const double thickness_cm= geom_.thickness_mm  * 0.1; // mm to cm
  const double width_cm    = geom_.cols * pix_cm; // cm
  const double height_cm   = geom_.rows * pix_cm;
  const double area_cm2    = width_cm * height_cm;
  const double volume_cm3  = area_cm2 * thickness_cm * geom_.active_fraction;
  const double mass_g      = mat_.density_g_cm3 * volume_cm3;
  return mass_g * 1e-3;
}

// ----------------------------------------------------------------------------
// Detector::mass_kg
//   Active target mass [kg]: the config override if I was given one,
//   otherwise the geometric mass.
// ----------------------------------------------------------------------------
double Detector::mass_kg() const {
  if (mass_override_kg_) return *mass_override_kg_;
  return compute_mass_from_geometry_kg();
}

} // namespace ccdarksens
