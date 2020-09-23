#ifndef HAERO_COLUMN_BASE_HPP
#define HAERO_COLUMN_BASE_HPP

#include "haero/haero_config.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

class ColumnBase {
  public :
  using kokkos_device_types = ekat::KokkosTypes<ekat::DefaultDevice>;
  using kokkos_host_types = ekat::KokkosTypes<ekat::HostDevice>;
  using real_pack_type = ekat::Pack<Real,HAERO_PACK_SIZE>;
  using mask_type = ekat::Mask<HAERO_PACK_SIZE>;

  using view_1d = kokkos_device_types::view_1d<real_pack_type>;
  using modal_var_view = kokkos_device_types::view_2d<Real>;
  using mask_view = kokkos_device_types::view_1d<mask_type>;

  using pack_info = ekat::PackInfo<HAERO_PACK_SIZE>;

  int nlev;
  int npack_mid;
  int npack_interface;

  ColumnBase() = delete;

  protected:
    ColumnBase(const int& nl) : nlev(nl), npack_mid(pack_info::num_packs(nl)),
      npack_interface(pack_info::num_packs(nl+1)) {};
};

} // namespace haero
#endif
