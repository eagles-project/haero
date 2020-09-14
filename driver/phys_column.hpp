#ifndef HAERO_PHYS_COLUMN_HPP
#define HAERO_PHYS_COLUMN_HPP

#include "haero/haero_config.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include <map>
#include <string>
#include <vector>

namespace haero {

struct phys_column {
  using kokkos_device_types = ekat::KokkosTypes<ekat::DefaultDevice>;
  using kokkos_host_types = ekat::KokkosTypes<ekat::HostDevice>;
  using real_pack_type = ekat::Pack<Real,HAERO_PACK_SIZE>;

  using view_1d = kokkos_device_types::view_1d<real_pack_type>;

  phys_column() = delete;

  phys_column(const int& nlevs);

  phys_column(const int& nlevs, const std::vector<std::string>& var_names);

  std::map<std::string, view_1d> level_vars;
  std::map<std::string, view_1d> interface_vars;

  int npack_mid;
  int npack_interface;
  int nlev;
};


} // namespace haero
#endif
