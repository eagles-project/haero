#ifndef HAERO_MAM4_WATER_UPTAKE_HPP
#define HAERO_MAM4_WATER_UPTAKE_HPP

#include "haero/process.hpp"

namespace haero {

/// @class Mam4WaterUptake
/// This process implements MAM4's treatment of water uptake, a diagnostic
/// aerosol process that updates the following diagnostic variables:
/// * `aero_water`, the per-mode aerosol water mixing ratio [kg/kg]
/// * `total_aero_water`, the aerosol water mixing ratio summed over modes [kg/kg]
/// * `mean_wet_diameter`, the per-mode geometric mean wet diameter [m]
/// * `wet_density`, the density of wet air [kg / m**3] <-- ???
class Mam4WaterUptake: public FDiagnosticProcess {
  public:

  /// Creates a process that simulates water uptake.
  /// @param [in] use_bisection If true, solves the Kohler equation using
  ///                           the bisection method. Otherwise, solves the
  ///                           Kohler equation using an algebraic method.
  explicit Mam4WaterUptake(bool use_bisection = false);
};

}

#endif
