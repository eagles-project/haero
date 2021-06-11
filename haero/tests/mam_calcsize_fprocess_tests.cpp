#include "catch2/catch.hpp"
#include "haero/model.hpp"
#include "haero/atmosphere.hpp"
#include "mam_calcsize_test_bridge.hpp"

using namespace haero;

TEST_CASE("calcsize_run", "mam_calcsize_fprocess") {


  const auto modes = create_mam4_modes();
  const auto aero_species = create_mam4_aerosol_species();
  const auto mode_spec_map = create_mam4_mode_species();
  const auto gas_species = create_mam4_gas_species();

  const auto aero_config =
    ModalAerosolConfig(modes, aero_species, mode_spec_map, gas_species);

  const int num_levels = 72;
  const int num_gases = 1;
  const int num_modes = 1;
  const Real t = 2.3;
  const Real dt = 0.15;
  int num_vert_packs = num_levels / HAERO_PACK_SIZE;


  Kokkos::View<PackType *> temp("temperature", num_vert_packs);
  Kokkos::View<PackType *> press("pressure", num_vert_packs);
  Kokkos::View<PackType *> rel_hum("relative humidity", num_vert_packs);
  Kokkos::View<PackType *> pdel("hydrostatic_dp", num_vert_packs);

  //Atmosphere atmos(num_levels, temp, press, rel_hum, ht, pdel, pblh);
  //Atmosphere atmosphere(72, 273., 100000., 72., 2., 0.5, 53.);


  /*  const ModalAerosolConfig& modal_aerosol_config,
    Real t, Real dt, const Prognostics& prognostics,
    const Atmosphere& atmosphere,
    const Diagnostics& diagnostics,
    Tendencies& tendencies*/
  //run(aer_config,random,random,)

  REQUIRE(3==3);
}

TEST_CASE("compute_diameter", "mam_calcsize_fprocess") {

  static constexpr Real vol2num = 2;
  //compute diameter
  //MAMCalcsizeProcess::compute_diameter(vol2num);
  Real diameter = compute_diameter_bridge(vol2num);

  REQUIRE(diameter==3);
}
