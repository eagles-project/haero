#ifndef HAERO_PHYS_COLUMN_HPP
#define HAERO_PHYS_COLUMN_HPP

#include "haero_config.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"

namespace haero {

struct phys_column {
  using kokkos_device_types = ekat::KokkosTypes<ekat::DefaultDevice>;
  using kokkos_host_types = ekat::KokkosTypes<ekat::HostDevice>;
  using real_pack_type = ekat::Pack<Real,HAERO_PACK_SIZE>;

  using view_1d = kokkos_device_types::view_1d<real_pack_type>;
};


} // namespace haero
#endif
