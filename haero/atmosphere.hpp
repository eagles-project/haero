#ifndef HAERO_ATMOSPHERE_HPP
#define HAERO_ATMOSPHERE_HPP

#include "haero/haero.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

/// @class Atmosphere
/// This type stores atmospheric state variables inherited from a host model.
class Atmosphere final {
  public:

  /// This type represents an array mapping a vertical level index to a pack.
  /// The vertical level(s) are identified by the index.
  /// Our views are unmanaged in general, to allow a host model to assume
  /// responsibility for managing resources.
  using ColumnView = ekat::Unmanaged<Kokkos::View<PackType*> >;

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
  Atmosphere(int num_levels,
             const Kokkos::View<PackType*>& temp,
             const Kokkos::View<PackType*>& press,
             const Kokkos::View<PackType*>& rel_hum,
             const Kokkos::View<PackType*>& ht);

  /// Destructor.
  ~Atmosphere();

  /// Returns the number of vertical levels per column in the system.
  int num_levels() const { return num_levels_; }

  /// Returns atmospheric temperature data [K], defined at each vertical level
  const ColumnView& temperature() const { return temperature_; }

  /// Returns atmospheric pressure data [Pa], defined at each vertical level.
  const ColumnView& pressure() const { return pressure_; }

  /// Returns relative humidity data [-], defined at each vertical level.
  const ColumnView& relative_humidity() const { return relative_humidity_; }

  /// Returns height data [m], defined at interfaces between vertical levels.
  const ColumnView& height() const { return height_; }

  private:

  // Number of vertical levels.
  const int num_levels_;

  // Unmanaged views.
  const ColumnView temperature_;
  const ColumnView pressure_;
  const ColumnView relative_humidity_;
  const ColumnView height_;
};

}

#endif
