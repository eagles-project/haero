#ifndef HAERO_COLUMN_BASE_HPP
#define HAERO_COLUMN_BASE_HPP

#include "haero/haero_config.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

using index_pair = std::pair<int,int>;

KOKKOS_INLINE_FUNCTION
index_pair level_to_pack_ind(const int& level_index, const int& npacks) {
  index_pair result;
  result.first = level_index/npacks;
  result.second = level_index%npacks;
  return result;
}

KOKKOS_INLINE_FUNCTION
int pack_to_level_ind(const int& pack_num, const int& pack_ind) {
  return HAERO_PACK_SIZE*pack_num + pack_ind;
}

KOKKOS_INLINE_FUNCTION
int pack_to_level_ind(const index_pair& pind) {
  return pack_to_level_ind(pind.first, pind.second);
}

class column_base {
  public :
  using kokkos_device_types = ekat::KokkosTypes<ekat::DefaultDevice>;
  using kokkos_host_types = ekat::KokkosTypes<ekat::HostDevice>;
  using real_pack_type = ekat::Pack<Real,HAERO_PACK_SIZE>;
  using mask_type = ekat::Mask<HAERO_PACK_SIZE>;

  using view_1d = kokkos_device_types::view_1d<real_pack_type>;
  using mask_view = kokkos_device_types::view_1d<mask_type>;

  protected:
    column_base() {}

};

} // namespace haero
#endif
