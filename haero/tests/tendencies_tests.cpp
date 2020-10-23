#include "haero/tendencies.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>

using namespace haero;

TEST_CASE("tendencies_ctor", "") {

  // Create a set of prognostics and add some modes and species to it.
  Prognostics progs(10, 72);
  Mode mode1("mode1", 1e-3, 1e-1, 1e-4, 0.0, 0.0);
  std::vector<Species> mode1_species;
  mode1_species.push_back(Species("Sulfate", "SO4"));
  progs.add_aerosol_mode(mode1, mode1_species);

  // Try to create tendencies for these prognostics before assembly.
  REQUIRE_THROWS(Tendencies(progs));

  // Assemble the prognostic variables, create tendencies for it, and make sure
  // the vitals match up.
  progs.assemble();
  Tendencies tends(progs);
  REQUIRE(tends.num_aerosol_modes() == progs.num_aerosol_modes());
  for (int m = 0; m < tends.num_aerosol_modes(); ++m) {
    REQUIRE(tends.num_aerosol_species(m) == progs.num_aerosol_species(m));
  }
  REQUIRE(tends.num_gas_species() == progs.num_gas_species());
  REQUIRE(tends.num_columns() == progs.num_columns());
  REQUIRE(tends.num_levels() == progs.num_levels());
  for (int m = 0; m < tends.num_aerosol_modes(); ++m) {
    const auto& tends_int_aeros = tends.interstitial_aerosols(m);
    const auto& progs_int_aeros = progs.interstitial_aerosols(m);
    REQUIRE(tends_int_aeros.extent(0) == progs_int_aeros.extent(0));
    REQUIRE(tends_int_aeros.extent(1) == progs_int_aeros.extent(1));
    REQUIRE(tends_int_aeros.extent(2) == progs_int_aeros.extent(2));

    const auto& tends_cld_aeros = tends.cloudborne_aerosols(m);
    const auto& progs_cld_aeros = progs.cloudborne_aerosols(m);
    REQUIRE(tends_cld_aeros.extent(0) == progs_cld_aeros.extent(0));
    REQUIRE(tends_cld_aeros.extent(1) == progs_cld_aeros.extent(1));
    REQUIRE(tends_cld_aeros.extent(2) == progs_cld_aeros.extent(2));
  }

  const auto& tends_modal_num_densities = tends.modal_num_densities();
  const auto& progs_modal_num_densities = progs.modal_num_densities();
  REQUIRE(tends_modal_num_densities.extent(0) == progs_modal_num_densities.extent(0));
  REQUIRE(tends_modal_num_densities.extent(1) == progs_modal_num_densities.extent(1));
  REQUIRE(tends_modal_num_densities.extent(2) == progs_modal_num_densities.extent(2));

  REQUIRE(tends.num_levels() == progs.num_levels());
}

