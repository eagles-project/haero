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
  /// @param [in] press A view of total (moist) pressure column data [Pa]
  ///                   managed by the host model
  /// @param [in] qv A view of water vapor mixing ratio column data
  ///                [kg vapor/kg dry air] managed by the host model
  /// @param [in] ht A view of height column data [m] on level interfaces,
  ///                managed by the host model
  /// @param [in] pblh The column-specific planetary boundary height [m],
  ///                  computed by the host model
  /// @param [in] pdel The hydrostatic "pressure thickness" defined as the
  ///                  difference in hydrostatic pressure levels at interfaces
  ///                  bounding each vertical level [Pa]
  Atmosphere(int num_levels, const ColumnView temp, const ColumnView press,
             const ColumnView qv, const ColumnView ht, const ColumnView pdel,
             Real pblh);

  Atmosphere() = default;
  /// Destructor.
  KOKKOS_FUNCTION
  ~Atmosphere() {}

  // Views.
  ColumnView temperature;
  ColumnView pressure;
  ColumnView vapor_mixing_ratio;
  ColumnView height;
  ColumnView hydrostatic_dp;

  // Planetary boundary height.
  Real planetary_boundary_height;

  /// Returns the number of vertical levels per column in the system.
  KOKKOS_INLINE_FUNCTION
  int num_levels() const { return num_levels_; }

  /// Sets the planetary boundary height [m].
  KOKKOS_INLINE_FUNCTION
  void set_planetary_boundary_height(Real pblh) {
    planetary_boundary_height = pblh;
  }

 private:
  // Number of vertical levels.
  int num_levels_;
};

}  // namespace haero

#endif
