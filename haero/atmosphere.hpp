#ifndef HAERO_ATMOSPHERE_HPP
#define HAERO_ATMOSPHERE_HPP

#include "haero/haero.hpp"

namespace haero {

/// @class Atmosphere
/// This type stores atmospheric state variables inherited from a host model.
class Atmosphere final {
  public:

  /// Creates an Atmosphere that stores unmanaged views of atmospheric column
  /// data owned and managed by the atmosphere host model.
  /// @param [in] num_levels the number of vertical levels per column stored by
  ///                        the state
  /// @param [in] temp A view of temperature column data [K] managed by the host
  ///                  model
  /// @param [in] press A view of pressure column data [Pa] managed by the host
  ///                   model
  /// @param [in] rel_hum A view of relative_humidity column data [-] managed
  ///                     by the host model
  /// @param [in] ht A view of height column data [m] on level interfaces,
  ///                managed by the host model
  /// @param [in] pblh The column-specific planetary boundary height [m],
  ///                  computed by the host model
  Atmosphere(int num_levels,
             const ColumnView temp,
             const ColumnView press,
             const ColumnView rel_hum,
             const ColumnView ht,
             Real pblh);

  /// Destructor.
  KOKKOS_FUNCTION
  ~Atmosphere();

  /// Returns the number of vertical levels per column in the system.
  KOKKOS_INLINE_FUNCTION
  int num_levels() const { return num_levels_; }

  /// Returns atmospheric temperature data [K], defined at each vertical level
  KOKKOS_INLINE_FUNCTION
  const ColumnView temperature() const { return temperature_; }

  /// Returns atmospheric pressure data [Pa], defined at each vertical level.
  KOKKOS_INLINE_FUNCTION
  const ColumnView pressure() const { return pressure_; }

  /// Returns relative humidity data [-], defined at each vertical level.
  KOKKOS_INLINE_FUNCTION
  const ColumnView relative_humidity() const { return relative_humidity_; }

  /// Returns height data [m], defined at interfaces between vertical levels.
  KOKKOS_INLINE_FUNCTION
  const ColumnView height() const { return height_; }

  /// Returns the planetary boundary height [m].
  KOKKOS_INLINE_FUNCTION
  const Real planetary_boundary_height() const { return pblh_; }

  private:

  // Number of vertical levels.
  const int num_levels_;

  // Views.
  const ColumnView temperature_;
  const ColumnView pressure_;
  const ColumnView relative_humidity_;
  const ColumnView height_;

  // Planetary boundary height.
  Real pblh_;
};

}

#endif
