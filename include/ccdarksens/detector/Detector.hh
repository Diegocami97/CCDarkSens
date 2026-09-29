// ===========================================================================
//  Diego Venegas-Vargas
//  DAMIC-M collaboration
//  CCDarkSens Framework
//
//  Detector.hh -- I declare the detector description: the CCD geometry
//  (DetectorGeometry), the target material (TargetMaterial) and the Detector
//  class that combines them and returns the active target mass in kg.
// ===========================================================================

#pragma once
#include <optional>
#include <string>

namespace ccdarksens {

// ----------------------------------------------------------------------------
// DetectorGeometry
//   Physical layout of the CCD: pixel array, pixel pitch, thickness and the
//   live (active) fraction of the array.
// ----------------------------------------------------------------------------
struct DetectorGeometry {
  int    rows = 0;  // number of pixel rows
  int    cols = 0;  // number of pixel columns
  double pixel_size_um = 0.0;   // micrometers
  double thickness_mm  = 0.0;   // millimeters
  double active_fraction = 1.0; // 0..1
};

// ----------------------------------------------------------------------------
// TargetMaterial
//   Target element and bulk density. Z and A are optional (0 = unknown) and
//   are only needed by the nuclear-recoil channels.
// ----------------------------------------------------------------------------
struct TargetMaterial {
  std::string element;          // "Si", "Ge", ...
  int         Z = 0;            // optional (0 = unknown)
  double      A = 0.0;          // optional (0 = unknown)
  double      density_g_cm3 = 0.0;  // bulk density [g/cm^3]
};

// ----------------------------------------------------------------------------
// Detector
//   I validate a geometry + material pair and provide the active target mass.
//   The mass is either computed from the geometry (pixels x pitch^2 x
//   thickness x active_fraction x density) or taken from an explicit override.
// ----------------------------------------------------------------------------
class Detector {
public:
  // Constructor. Throws std::invalid_argument if the geometry or material is unphysical.
  Detector(DetectorGeometry geom, TargetMaterial mat, std::optional<double> mass_kg_override);

  const DetectorGeometry& geometry() const noexcept { return geom_; }
  const TargetMaterial&   material() const noexcept { return mat_; }

  double mass_kg() const;  // Computes or returns override

private:
  // Mass [kg] from geometry and density (used when there is no override).
  double compute_mass_from_geometry_kg() const;

  DetectorGeometry        geom_;  // pixel array / thickness / active fraction
  TargetMaterial          mat_;  // target element and density
  std::optional<double>   mass_override_kg_;  // if set, mass_kg() returns this instead of the geometric mass
};

} // namespace ccdarksens
