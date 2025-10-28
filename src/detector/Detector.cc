#include "ccdarksens/detector/Detector.hh"
#include <cmath>
#include <stdexcept>

namespace ccdarksens {

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

double Detector::mass_kg() const {
  if (mass_override_kg_) return *mass_override_kg_;
  return compute_mass_from_geometry_kg();
}

} // namespace ccdarksens
