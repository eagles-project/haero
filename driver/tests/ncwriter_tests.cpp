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
using view_2d = kokkos_device_type::view_2d<real_pack_type>;
using view_3d = kokkos_device_type::view_3d<real_pack_type>;

TEST_CASE("ncwriter", "") {

  /// Create a new netcdf file
  const std::string fname = "test_file.nc";
  NcWriter ncf(fname);
  const NcWriter::text_att_type file_desc = std::make_pair("NcWriter_unit_test",
    "Testing purposes only; some variables will be undefined.");
  ncf.add_file_attribute(file_desc);

  /** The constructor automatically defines a dimension for time.

    no other dimensions and no variables are defined so far.
  */
  REQUIRE (ncf.get_ndims() == 1);
  REQUIRE (ncf.get_nvars() == 0);

  const int ncol = 2;
  ncf.add_column_dim(ncol);
  REQUIRE (ncf.get_ndims() == 2);

  /**
    Add a dimension for level midpoints and a dimension for level interfaces.
  */
  const int nlev = 10;
  ncf.add_level_dims(nlev);
  REQUIRE (ncf.get_ndims() == 4);

  /**
    Add a dimension for modes

    @todo We may not keep this
  */
  const int nmodes = 4;
  const std::vector<std::string> mode_names = {"aitken", "accumulation", "black_carbon", "coarse"};
  std::vector<std::pair<Real,Real>> diameter_minmax;
  diameter_minmax.push_back(std::make_pair<Real>(10E-6, 30E-6));
  diameter_minmax.push_back(std::make_pair<Real>(30E-6, 300E-6));
  diameter_minmax.push_back(std::make_pair<Real>(20E-6, 200E-6));
  diameter_minmax.push_back(std::make_pair<Real>(200E-6, 4000E-6));
  ncf.add_mode_dim(nmodes, mode_names, diameter_minmax);
  REQUIRE (ncf.get_ndims() == 5);

  /**
    Define the coordinate variable for the time dimension.

    This automatically initializes the first time value to zero.
  */
  ncf.define_time_var();

  /**
    Create views for basic level and interface variables,
      and define netCDF variables from their metadata.

    Note: No data is copied at this stage, so we don't have to worry about
    which memory space it's in.
  */
  view_2d test_level_var("test_level_var", ncol, nlev);
  view_2d test_interface_var("interface_var", ncol, nlev+1);

  ncf.define_level_var("level_test_var", ekat::units::Units::nondimensional(), test_level_var);
  REQUIRE (ncf.get_varid("level_test_var") != NC_EBADID);

  ncf.define_interface_var("plus_minus_two", ekat::units::Units::nondimensional(), test_interface_var);
  REQUIRE (ncf.get_varid("plus_minus_two") != NC_EBADID);

  /**
    Create views for modal aerosols, and define corresponding netCDF variables.
  */
  view_3d test_modal_var("plus_minus_modenum", ncol, nmodes, nlev);
  ncf.define_modal_var("test_aerosol_from_view", ekat::units::pow(ekat::units::kg, -1), test_modal_var);
  /**
    Initialize a single column's worth of data on the host, copy to device
  */
  auto modal_col0 = Kokkos::subview(test_modal_var, 0, Kokkos::ALL, Kokkos::ALL);
  auto hmodalc0 = Kokkos::create_mirror_view(modal_col0);
  for (int i = 0; i<nmodes; ++i ) {
    for (int j = 0; j<nlev; ++j) {
      hmodalc0(i,j) = i * std::pow(-1,j);
    }
  }
  Kokkos::deep_copy(modal_col0, hmodalc0);
  /**
    Add the data to each column in the netCDF file.
  */
  ncf.add_variable_data("test_aerosol_from_view", 0, 0, test_modal_var);

  /// Test the stripped-down version of modal variable definition
  ncf.define_modal_var("test_aerosol", ekat::units::pow(ekat::units::kg, -1));

  /**
    The default data type for data arrays is a 1D std::vector; it's always on the host.

    Hence, we can write data from a std::vector to the nc file.
  */
  ncf.define_level_var("plus_minus_one", ekat::units::Units::nondimensional());
  std::vector<Real> pm1(nlev);
  for (int i=0; i<nlev; ++i) {
    pm1[i] = std::pow(-1, i);
  }
  for (int i=0; i<ncol; ++i) {
    ncf.add_level_variable_data("plus_minus_one", 0, i, pm1);
  }

  /**
    Test writing to netCDF file from a view
  */
  auto hpm2 = Kokkos::create_mirror_view(test_interface_var);
  for (int i=0; i<nlev+1; ++i) {
    hpm2(0,i) = 2*std::pow(-1,i);
    hpm2(1,i) = 2*std::pow(-1,i);
  }
  Kokkos::deep_copy(test_interface_var, hpm2);
  for (int i=0; i<ncol; ++i) {
    ncf.add_variable_data("plus_minus_two", 0, i, test_interface_var);
  }


  /**
    Here we try to add data to a time step index that is out of bounds.
  */
  REQUIRE_THROWS (ncf.add_level_variable_data("plus_minus_one", 1, 0, pm1));
  /**
    Add the time index and try again.
  */
  const Real t=1.0;
  ncf.add_time_value(t);
  REQUIRE_NOTHROW (ncf.add_level_variable_data("plus_minus_one", 1, 0, pm1));

  ncf.close();
  std::cout << ncf.info_string();
}
