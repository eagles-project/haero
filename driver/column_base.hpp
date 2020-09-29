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
  using view_2d = kokkos_device_types::view_2d<real_pack_type>;
  using modal_var_view = kokkos_device_types::view_2d<Real>;
  using mask_view = kokkos_device_types::view_1d<mask_type>;

  using pack_info = ekat::PackInfo<HAERO_PACK_SIZE>;

  inline int num_levels() const {return m_num_levels;}
  inline int num_packs_lev() const {return m_num_packs_lev;}
  inline int num_packs_int() const {return m_num_packs_int;}

  ColumnBase() = delete;

  protected:
    ColumnBase(const int& nl) : m_num_levels(nl), m_num_packs_lev(pack_info::num_packs(nl)),
      m_num_packs_int(pack_info::num_packs(nl+1)) {}

  int m_num_levels;
  int m_num_packs_lev;
  int m_num_packs_int;
};

} // namespace haero
#endif
