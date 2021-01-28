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

  /// This is the device on which the Atmosphere stores its data.
  using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;

  /// This type represents vectorizable packs of Reals of length HAERO_PACK_SIZE.
  using PackType = ekat::Pack<Real, HAERO_PACK_SIZE>;

  /// This type represents a multidimensional array mapping a column and
  /// vertical level index to a pack.
  /// * The column is identified by the index i.
  /// * The vertical level(s) are identified by the index k.
  /// So view[i][k] yields the desired pack.
  /// Our views are unmanaged in general, to allow a host model to assume
  /// responsibility for managing resources.
  using ColumnView = ekat::Unmanaged<Kokkos::View<PackType**> >;

  /// Creates an Atmosphere that stores unmanaged views of atmospheric column
  /// data owned and managed by the atmosphere host model.
  /// @param [in] num_columns the number of vertical columns stored by the state
  /// @param [in] num_levels the number of vertical levels per column stored by
  ///                        the state
  /// @param [in] temp A view of temperature data [K] managed by the host
  ///                  model
  /// @param [in] press A view of pressure data [Pa] managed by the host
  ///                   model
  /// @param [in] rel_hum A view of relative_humidity data [-] managed
  ///                     by the host model
  /// @param [in] ht A view of height data [m] on level interfaces, managed
  ///                by the host model
  Atmosphere(int num_columns,
             int num_levels,
             const Kokkos::View<PackType**>& temp,
             const Kokkos::View<PackType**>& press,
             const Kokkos::View<PackType**>& rel_hum,
             const Kokkos::View<PackType**>& ht):
    num_columns_(num_columns),
    num_levels_(num_levels),
    temperature_(temp),
    pressure_(press),
    relative_humidity_(rel_hum),
    height_(ht) {}

  /// Destructor.
  ~Atmosphere() {}

  /// Returns the number of independent atmospheric columns in the system.
  int num_columns() const { return num_columns_; }

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

  // Number of columns / vertical levels.
  int num_columns_;
  int num_levels_;

  // Unmanaged views
  ColumnView temperature_;
  ColumnView pressure_;
  ColumnView relative_humidity_;
  ColumnView height_;
};

}

#endif
