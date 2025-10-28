#pragma once
#include <optional>
#include <string>

namespace ccdarksens {

struct DetectorGeometry {
  int    rows = 0;
  int    cols = 0;
  double pixel_size_um = 0.0;   // micrometers
  double thickness_mm  = 0.0;   // millimeters
  double active_fraction = 1.0; // 0..1
};

struct TargetMaterial {
  std::string element;          // "Si", "Ge", ...
  int         Z = 0;            // optional (0 = unknown)
  double      A = 0.0;          // optional (0 = unknown)
  double      density_g_cm3 = 0.0;
};

class Detector {
public:
  Detector(DetectorGeometry geom, TargetMaterial mat, std::optional<double> mass_kg_override);

  const DetectorGeometry& geometry() const noexcept { return geom_; }
  const TargetMaterial&   material() const noexcept { return mat_; }

  double mass_kg() const;  // Computes or returns override

private:
  double compute_mass_from_geometry_kg() const;

  DetectorGeometry        geom_;
  TargetMaterial          mat_;
  std::optional<double>   mass_override_kg_;
};

} // namespace ccdarksens
