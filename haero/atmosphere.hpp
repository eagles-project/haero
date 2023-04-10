// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HAERO_ATMOSPHERE_HPP
#define HAERO_ATMOSPHERE_HPP

#include <haero/haero.hpp>

namespace haero {

/// @class Atmosphere
/// This type stores atmospheric state variables inherited from a host model.
class Atmosphere final {
  // number of vertical levels
  int num_levels_;

public:
  /// default constructor (individual views must be set manually)
  Atmosphere() = default;

  // Copy construction and assignment are supported for moving data between
  // host and device, and for populating multi-column views.
  Atmosphere(const Atmosphere &) = default;
  Atmosphere &operator=(const Atmosphere &) = default;

  /// destructor, valid on both host and device
  KOKKOS_FUNCTION
  ~Atmosphere() {}

  // views storing atmospheric state data for a single vertical column

  /// temperature [K]
  ColumnView temperature;

  /// pressure [Pa]
  ColumnView pressure;

  /// water vapor mass mixing ratio [kg vapor/kg dry air]
  ColumnView vapor_mixing_ratio;

  /// height at the midpoint of each vertical level [m]
  ColumnView height;

  /// hydroѕtatic "pressure thickness" defined as the difference in hydrostatic
  /// pressure levels between the interfaces bounding a vertical level [Pa]
  ColumnView hydrostatic_dp;

  /// cloud fraction [-]
  ColumnView cloud_fraction;

  /// vertical updraft velocity used for ice nucleation [m/s]
  ColumnView updraft_vel_ice_nucleation;

  // column-specific planetary boundary height [m]
  Real planetary_boundary_height;

  /// returns the number of vertical levels per column in the system
  KOKKOS_INLINE_FUNCTION
  int num_levels() const { return num_levels_; }

  /// Sets the planetary boundary height [m].
  KOKKOS_INLINE_FUNCTION
  void set_planetary_boundary_height(Real pblh) {
    planetary_boundary_height = pblh;
  }

  /// Returns true iff all atmospheric quantities are nonnegative, using the
  /// given thread team to parallelize the check.
  KOKKOS_INLINE_FUNCTION
  bool quantities_nonnegative(const ThreadTeam &team) const {
    const int nk = num_levels();
    int violations = 0;
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, nk),
        KOKKOS_CLASS_LAMBDA(int k, int &violation) {
          if ((temperature(k) < 0) || (pressure(k) < 0) ||
              (vapor_mixing_ratio(k) < 0)) {
            violation = 1;
          }
        },
        violations);
    return (violations == 0);
  }
};

} // namespace haero

#endif
