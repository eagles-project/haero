#include "haero/haero_config.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include <iostream>
#include <string>
#include "catch2/catch.hpp"

using namespace haero;
using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
using real_pack_type = ekat::Pack<Real, HAERO_PACK_SIZE>;
using view_1d = kokkos_device_type::view_1d<real_pack_type>;

TEST_CASE("ncwriter", "") {
  const std::string fname = "test_file.nc";
  NcWriter ncf(fname);

  REQUIRE (ncf.get_ndims() == 1);
  REQUIRE (ncf.get_nvars() == 0);

  const int nlev = 10;
  ncf.add_level_dims(nlev);

  REQUIRE (ncf.get_ndims() == 3);

  const int nmodes = 4;
  ncf.add_mode_dim(nmodes);

  REQUIRE (ncf.get_ndims() == 4);

  ncf.define_time_var<Real>();

  view_1d test_level_var("test_level_var", nlev);
  view_1d test_interface_var("interface_var", nlev+1);

  ncf.define_level_var("level_test_var", ekat::units::Units::nondimensional(), test_level_var);
  REQUIRE (ncf.get_nvars() == 2);
  REQUIRE (ncf.get_varid(test_level_var) != NC_EBADID);

  ncf.define_interface_var("interface_test_var", ekat::units::Units::nondimensional(), test_interface_var);
  REQUIRE (ncf.get_nvars() == 3);
  REQUIRE (ncf.get_varid(test_interface_var) != NC_EBADID);

  ncf.close();
  std::cout << ncf.info_string();
}
