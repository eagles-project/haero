#include "driver/column_netcdf.hpp"
#include <iostream>
#include <string>
#include "catch2/catch.hpp"

using namespace haero;

TEST_CASE("column_netcdf", "") {
  const std::string fname = "test_file.nc";
  NetCDFFileHandler ncf(fname);

  REQUIRE (ncf.get_ndims() == 0);
  REQUIRE (ncf.get_nvars() == 0);

  const int nlev = 10;
  ncf.add_level_dims(nlev);

  REQUIRE (ncf.get_ndims() == 2);

  const int nmodes = 4;
  ncf.add_mode_dim(nmodes);

  REQUIRE (ncf.get_ndims() == 3);

  ncf.add_time_dim();

  REQUIRE (ncf.get_ndims() == 4);

  ncf.close_file();
  std::cout << ncf.info_string();
}
