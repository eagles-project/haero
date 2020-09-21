#ifndef HAERO_ATMOSPHERE_HPP
#define HAERO_ATMOSPHERE_HPP

#include "haero/haero_config.hpp"
#include "haero/physical_constants.hpp"

namespace haero {

/// This type defines ambient atmospheric conditions for supported models.
struct AtmosphericConditions {
  /// The model used to generate the atmospheric conditions.
  /// * uniform: all vertical levels have identical atmospheric conditions
  /// * hydrostatic: the vertical profile of the atmospheric conditions is
  ///                determined using the hydrostatic equilibrium approximation.
  enum { uniform, hydrostatic } model;
  /// Parameters specific to each model.
  union {
    /// Uniform model.
    struct {
      /// Mean molecular weight of air [kg/mol].
      Real mu;
      /// Scaled atmospheric height [m].
      Real H;
      /// Pressure [Pa].
      Real p0;
      /// Temperature [K].
      Real T0;
      /// Clear-sky relative humidity [-]
      Real phi0;
      /// Cloud fraction [-]
      Real N0;
    } uniform;
    /// Hydrostatic model (we assume a constant lapse rate).
    struct {
      /// Reference pressure at z = 0 [Pa].
      Real p0;
      /// Reference temperature at z = 0 [K].
      Real T0;
      /// Temperature lapse rate @f$ \Gamma = -\partial T / \partial z @f$ [K/m]
      Real lapse_rate;
    } hydrostatic;
  } params;
};

AtmosphericConditions stable_hydrostatic(const Real& refp=p_std_pa,
                                         const Real& reft = 295.0,
                                         const Real& gamma=0.5*dry_adiabatic_lapse_rate_K_per_m);

AtmosphericConditions neutrally_stable_hydrostatic(const Real& refp=p_std_pa,
                                         const Real& reft = 295.0);

AtmosphericConditions unstable_hydrostatic(const Real& refp=p_std_pa,
                                         const Real& reft = 295.0,
                                         const Real& gamma=1.5*dry_adiabatic_lapse_rate_K_per_m);


} // namespace haero
#endif
