#include "available_diagnostics.hpp"
#include "haero/modal_aerosol_config.hpp"
#include "haero/math_helpers.hpp"
#include "haero/mode.hpp"
#include "haero/aerosol_species.hpp"
#include "catch2/catch.hpp"
#include <iostream>

using namespace haero;

TEST_CASE ("wet_radius_diagnostic", "") {

  using KohlerPoly = KohlerPolynomial<PackType>;

  // setup aerosol config
  const auto config = create_simple_test_config();
  std::cout << config.info_string();

  REQUIRE(config.num_modes() == 2);
  const int nmodes = config.num_modes();
  REQUIRE(config.num_gases() == 1);
  REQUIRE(config.max_species_per_mode() == 2);
  REQUIRE(config.num_aerosol_populations == 3);
  REQUIRE(config.aerosol_mode_index("test_mode0") == 0);
  REQUIRE(config.aerosol_mode_index("test_mode1") == 1);
  REQUIRE(config.aerosol_species_index(0, "TS0") == 0);
  REQUIRE(config.aerosol_species_index(0, "TS1") == -1);
  REQUIRE(config.aerosol_species_index(1, "TS0") == 0);
  REQUIRE(config.aerosol_species_index(1, "TS1") == 1);

  for (int m=0; m<config.num_modes(); ++m) {
    for (int s=0; s<config.max_species_per_mode(); ++s) {
      std::cout << "config.population_index(" << m << ", " << s << ") = " << config.population_index(m,s) << "\n";
    }
  }
  const auto pop_inds = config.create_population_indices_view();

  // setup host model stuff
  const int nlev = 8;
  const int npacks = PackInfo::num_packs(nlev);

  ColumnView relative_humidity("rel_humidity", npacks);
  auto h_relative_humidity = Kokkos::create_mirror_view(relative_humidity);
  const Real drelh = (KohlerPoly::rel_humidity_max - KohlerPoly::rel_humidity_min)/(nlev-1);
  for (int k=0; k<nlev; ++k) {
    const int pack_idx = PackInfo::pack_idx(k);
    const int vec_idx = PackInfo::vec_idx(k);
    const Real rh = KohlerPoly::rel_humidity_min + k * drelh;

    h_relative_humidity(pack_idx)[vec_idx] = rh;
  }
  Kokkos::deep_copy(relative_humidity, h_relative_humidity);

  SpeciesColumnView q_aero("aerosol_mass_mixing_ratios", config.num_aerosol_populations, npacks);
  auto h_q_aero = Kokkos::create_mirror_view(q_aero);
  for (int p=0; p<config.num_aerosol_populations; ++p) {
    for (int k=0; k<nlev; ++k) {
      const int pack_idx = PackInfo::pack_idx(k);
      const int vec_idx = PackInfo::vec_idx(k);
      h_q_aero(p,pack_idx)[vec_idx] = (p + k + 1)*1e-6;
    }
  }
  Kokkos::deep_copy(q_aero, h_q_aero);
  ModalColumnView num_ratios("aerosol_mode_number_mixing_ratios", config.num_modes(), npacks);
  auto h_num_ratios = Kokkos::create_mirror_view(num_ratios);
  for (int m=0; m<nmodes; ++m) {
    for (int k=0; k<nlev; ++k) {
      const int pack_idx = PackInfo::pack_idx(k);
      const int vec_idx = PackInfo::vec_idx(k);
      h_num_ratios(m, pack_idx)[vec_idx] = (k+1)*(m+1)*1e4;
    }
  }
  Kokkos::deep_copy(num_ratios, h_num_ratios);
  DeviceType::view_1d<AerosolSpecies> aerosols_in_mode("aerosols_in_mode", config.max_species_per_mode());

  // modal averages
  ModalColumnView mode_wet_radius("mode_wet_radius", nmodes, npacks);
  ModalColumnView mode_mean_particle_dry_volume("mode_mean_particle_dry_volume", nmodes, npacks);
  ModalColumnView mode_hygroscopicity("mode_hygroscopicity", nmodes, npacks);
  ModalColumnView mode_dry_particle_radius("mode_dry_particle_radius", nmodes, npacks);

  for (int m=0; m<nmodes; ++m) {

    config.aerosol_species_for_mode(m, aerosols_in_mode);

    Kokkos::parallel_for(npacks, ModalMeanParticleVolume(
      Kokkos::subview(mode_mean_particle_dry_volume, m, Kokkos::ALL),
      q_aero, num_ratios, aerosols_in_mode,
      Kokkos::subview(pop_inds, m, Kokkos::ALL), config.h_n_species_per_mode(m), m));

    Kokkos::parallel_for(npacks, ModalHygroscopicity(
      Kokkos::subview(mode_hygroscopicity, m, Kokkos::ALL),
      q_aero, aerosols_in_mode,
      Kokkos::subview(pop_inds, m, Kokkos::ALL), config.h_n_species_per_mode(m)));

    Kokkos::parallel_for(npacks, KOKKOS_LAMBDA (const int pack_idx) {
      mode_dry_particle_radius(m,pack_idx) =
        0.5*modal_mean_particle_diameter(mode_mean_particle_dry_volume(m, pack_idx),
          config.d_aerosol_modes(m).log_sigma);
    });

    Kokkos::parallel_for(npacks, ModalWetRadius(
      Kokkos::subview(mode_wet_radius, m, Kokkos::ALL),
      Kokkos::subview(mode_hygroscopicity, m, Kokkos::ALL),
      Kokkos::subview(mode_dry_particle_radius, m, Kokkos::ALL),
      relative_humidity));

  }

  // verification
  const auto mode0spec = config.aerosol_species_for_mode(0);
  auto h_mode_hyg = Kokkos::create_mirror_view(mode_hygroscopicity);
  auto h_mode_dry_vol = Kokkos::create_mirror_view(mode_mean_particle_dry_volume);
  Kokkos::deep_copy(h_mode_hyg, mode_hygroscopicity);
  Kokkos::deep_copy(h_mode_dry_vol, mode_mean_particle_dry_volume);
  for (int pack_idx=0; pack_idx<npacks; ++pack_idx) {
    REQUIRE( FloatingPoint<PackType>::equiv(h_mode_hyg(0,pack_idx), PackType(mode0spec[0].hygroscopicity)) );
    REQUIRE( FloatingPoint<PackType>::equiv(h_mode_dry_vol(0,pack_idx),  q_aero(0,pack_idx) / mode0spec[0].density / h_num_ratios(0,pack_idx)));
  }

  const auto mode1spec = config.aerosol_species_for_mode(1);
  for (int pack_idx=0; pack_idx<npacks; ++pack_idx) {
    const PackType mode_vol_mr = q_aero(1,pack_idx)/mode1spec[0].density + q_aero(2,pack_idx)/mode1spec[1].density;
    const PackType mhyg = (q_aero(1,pack_idx)*mode1spec[0].hygroscopicity/mode1spec[0].density +
        q_aero(2,pack_idx)*mode1spec[1].hygroscopicity/mode1spec[1].density) / mode_vol_mr;
    REQUIRE( FloatingPoint<PackType>::equiv(h_mode_dry_vol(1,pack_idx), mode_vol_mr / h_num_ratios(1,pack_idx)));
    REQUIRE( FloatingPoint<PackType>::equiv(h_mode_hyg(1,pack_idx), mhyg) );

  }

}
