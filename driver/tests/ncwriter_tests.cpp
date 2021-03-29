#include "haero/haero_config.hpp"
#include "haero/mode.hpp"
#include "haero/aerosol_species.hpp"
#include "driver/ncwriter.hpp"
#include "driver/ncwriter_impl.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include <iostream>
#include <string>
#include <cmath>
#include "catch2/catch.hpp"

using namespace haero;
using namespace haero::driver;
using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
using real_pack_type = ekat::Pack<Real, HAERO_PACK_SIZE>;
using view_2d = kokkos_device_type::view_2d<real_pack_type>;
using view_1d = kokkos_device_type::view_1d<real_pack_type>;

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

  /**
    Add a dimension for level midpoints and a dimension for level interfaces.
  */
  const int nlev = 10;
  ncf.add_level_dims(nlev);
  REQUIRE (ncf.get_ndims() == 3);

  /**
    Add a dimension for modes
  */
  std::vector<Mode> modes = create_mam4_modes();
  int nmodes = modes.size();
  ncf.add_mode_dim(modes);
  REQUIRE (ncf.get_ndims() == 4);

  /**
    Add a dimension for species
  */
  std::vector<AerosolSpecies> species = create_mam4_aerosol_species();
  int nspec = species.size();
  ncf.add_aerosol_dim(species);
  REQUIRE(ncf.get_ndims() == 5);
  REQUIRE(ncf.num_aerosols() == nspec);

  std::vector<GasSpecies> gases = create_mam4_gas_species();
  int ngas = gases.size();
  ncf.add_gas_dim(gases);
  REQUIRE(ncf.get_ndims() == 6);
  REQUIRE(ncf.num_gases() == ngas);

  /**
    Define the coordinate variable for the time dimension.

    This automatically initializes the first time value to zero.
  */
  ncf.define_time_var();

  /**
    Define a variable that is just one number (0-dimensional)
  */
  ncf.define_scalar_var("scalar_var_zero", ekat::units::Units::nondimensional());

  /**
    Create views for basic level and interface variables,
      and define netCDF variables from their metadata.

    Note: No data is copied at this stage, so we don't have to worry about
    which memory space it's in.
  */
  view_1d test_level_var("test_level_var", pack_info::num_packs(nlev));
  view_1d test_interface_var("interface_var", pack_info::num_packs(nlev+1));

  ncf.define_level_var("level_test_var", ekat::units::Units::nondimensional(), test_level_var);
  REQUIRE (ncf.get_varid("level_test_var") != NC_EBADID);

  ncf.define_interface_var("plus_minus_two", ekat::units::Units::nondimensional(), test_interface_var);
  REQUIRE (ncf.get_varid("plus_minus_two") != NC_EBADID);

  /**
    Create views for modal aerosols, and define corresponding netCDF variables.
  */
  view_2d test_modal_var("plus_minus_modenum", nmodes, pack_info::num_packs(nlev));
  /**
    Initialize a single column's worth of data on the host, copy to device
  */
  auto hmodal = Kokkos::create_mirror_view(test_modal_var);
  for (int i = 0; i<nmodes; ++i ) {
    for (int j = 0; j<nlev; ++j) {
      hmodal(i, pack_info::pack_idx(j))[pack_info::vec_idx(j)] = i * std::pow(-1,j);
    }
  }
  Kokkos::deep_copy(test_modal_var, hmodal);

  /// Test the stripped-down version of modal variable definition
  ncf.define_level_var("level_var_subview0", ekat::units::Units::nondimensional(),
    Kokkos::subview(test_modal_var,0,Kokkos::ALL));

  /**
    The default data type for data arrays is a 1D std::vector; it's always on the host.

    Hence, we can write data from a std::vector to the nc file.
  */
  ncf.define_level_var("plus_minus_one", ekat::units::Units::nondimensional());
  std::vector<Real> pm1(nlev);
  for (int i=0; i<nlev; ++i) {
    pm1[i] = std::pow(-1, i);
  }
  ncf.add_level_variable_data("plus_minus_one", 0, pm1);

  /**
    Test writing to netCDF file from a view
  */
  int time_idx = 0;
  int mode_idx = 0;
  int spec_idx = 0;
  auto hpm2 = Kokkos::create_mirror_view(test_interface_var);
  for (int i=0; i<nlev+1; ++i) {
    hpm2(pack_info::pack_idx(i))[pack_info::vec_idx(i)] = 2*std::pow(-1,i);
  }
  Kokkos::deep_copy(test_interface_var, hpm2);
  ncf.add_variable_data("plus_minus_two", time_idx, mode_idx, spec_idx, test_interface_var);


  /**
    Here we try to add data to a time step index that is out of bounds.
  */
  ++time_idx;
  REQUIRE_THROWS (ncf.add_level_variable_data("plus_minus_one", time_idx, pm1));
  /**
    Add the time index and try again.
  */
  const Real t=1.0;
  ncf.add_time_value(t);
  REQUIRE_NOTHROW (ncf.add_level_variable_data("plus_minus_one", time_idx,pm1));
  ncf.add_variable_data("plus_minus_two", time_idx, mode_idx, spec_idx, test_interface_var);

  ncf.close();
  std::cout << ncf.info_string();
}
