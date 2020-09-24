#include "haero/haero_config.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include "catch2/catch.hpp"

using namespace haero;
using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
using real_pack_type = ekat::Pack<Real, HAERO_PACK_SIZE>;
using view_1d = kokkos_device_type::view_1d<real_pack_type>;

TEST_CASE("ncwriter", "") {

  /// Create a new netcdf file
  const std::string fname = "test_file.nc";
  NcWriter ncf(fname);

  /** The constructor automatically defines a dimension for time.

    no other dimensions and no variables are defined so far.
  */
  REQUIRE (ncf.get_ndims() == 1);
  REQUIRE (ncf.get_nvars() == 0);

  /**
    Add a dimension for level midpoints and a dimension for level interfaces.
  */
  const int nlev = 10;
  ncf.add_level_dims(nlev);
  REQUIRE (ncf.get_ndims() == 3);

  /**
    Add a dimension for modes

    @todo We may not keep this
  */
  const int nmodes = 4;
  ncf.add_mode_dim(nmodes);
  REQUIRE (ncf.get_ndims() == 4);

  /**
    Define the coordinate variable for the time dimension.

    This automatically initializes the first time value to zero.
  */
  ncf.define_time_var();

  /**
    Create views and define variables from their metadata.

    Note: No data is copied at this stage, so we don't have to worry about
    which memory space it's in.
  */
  view_1d test_level_var("test_level_var", nlev);
  view_1d test_interface_var("interface_var", nlev+1);
  ncf.define_level_var("level_test_var", ekat::units::Units::nondimensional(), test_level_var);
  REQUIRE (ncf.get_nvars() == 2);
  REQUIRE (ncf.get_varid("level_test_var") != NC_EBADID);

  ncf.define_interface_var("interface_test_var", ekat::units::Units::nondimensional(), test_interface_var);
  REQUIRE (ncf.get_nvars() == 3);
  REQUIRE (ncf.get_varid("interface_test_var") != NC_EBADID);


  /**
    The default data type for data arrays is std::vector; it's always on the host.

    Hence, we can write data from a std::vector to the nc file.
  */
  ncf.define_level_var("plus_minus_one", ekat::units::Units::nondimensional());
  std::vector<Real> pm1(nlev);
  for (int i=0; i<nlev; ++i) {
    pm1[i] = std::pow(-1, i);
  }

  ncf.add_level_variable_data("plus_minus_one", 0, pm1);

  /**
    Here we try to add data to a time step index that is out of bounds.
  */
  REQUIRE_THROWS (ncf.add_level_variable_data("plus_minus_one", 1, pm1));
  /**
    Add the time index and try again.
  */
  const Real t=1.0;
  ncf.add_time_value(t);
  REQUIRE_NOTHROW (ncf.add_level_variable_data("plus_minus_one", 1, pm1));

  ncf.close();
  std::cout << ncf.info_string();
}
