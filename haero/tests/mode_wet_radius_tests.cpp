#include <iomanip>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/aerosol_species.hpp"
#include "haero/check.hpp"
#include "haero/diagnostics/kohler_solve.hpp"
#include "haero/diagnostics/mode_dry_particle_volume.hpp"
#include "haero/diagnostics/mode_hygroscopicity.hpp"
#include "haero/diagnostics/mode_wet_radius.hpp"
#include "haero/math.hpp"
#include "haero/modal_aerosol_config.hpp"
#include "haero/mode.hpp"

using namespace haero;

TEST_CASE("wet_radius_diagnostic", "") {
  using KohlerPoly = KohlerPolynomial<ekat::Pack<double, HAERO_PACK_SIZE>>;

  /** setup aerosol config

    This ModalAerosolConfig is simple enough to compute modal averages by hand,
    which we use for verification.
  */
  const auto config = create_simple_test_config();
  std::cout << config.info_string();

  /// Tests for create_simple_test_config();
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

  /// Write population indices to console
  for (int m = 0; m < config.num_modes(); ++m) {
    for (int s = 0; s < config.max_species_per_mode(); ++s) {
      std::cout << "config.population_index(" << m << ", " << s
                << ") = " << config.population_index(m, s) << "\n";
    }
  }
  /// Create a device view of population indices
  const auto pop_inds = config.create_population_indices_view();

  // setup host model stuff
  const int nlev = 40;
  const int npacks = PackInfo::num_packs(nlev);

  // define relative humidity for input data
  ColumnView relative_humidity("rel_humidity", npacks);
  auto h_relative_humidity = Kokkos::create_mirror_view(relative_humidity);
  const Real drelh =
      (KohlerPoly::rel_humidity_max - KohlerPoly::rel_humidity_min) /
      (nlev - 1);
  for (int k = 0; k < nlev; ++k) {
    const int pack_idx = PackInfo::pack_idx(k);
    const int vec_idx = PackInfo::vec_idx(k);
    const Real rh = KohlerPoly::rel_humidity_min + k * drelh;

    h_relative_humidity(pack_idx)[vec_idx] = rh;
  }
  Kokkos::deep_copy(relative_humidity, h_relative_humidity);

  // define mass mixing ratios for input data
  SpeciesColumnView q_aero("aerosol_mass_mixing_ratios",
                           config.num_aerosol_populations, npacks);
  auto h_q_aero = Kokkos::create_mirror_view(q_aero);
  for (int p = 0; p < config.num_aerosol_populations; ++p) {
    for (int k = 0; k < nlev; ++k) {
      const int pack_idx = PackInfo::pack_idx(k);
      const int vec_idx = PackInfo::vec_idx(k);
      h_q_aero(p, pack_idx)[vec_idx] = (p + 2 * k + 1) * 1e-6;
    }
  }
  Kokkos::deep_copy(q_aero, h_q_aero);

  // define mode number mixing ratios for input data
  ModalColumnView num_ratios("aerosol_mode_number_mixing_ratios",
                             config.num_modes(), npacks);
  auto h_num_ratios = Kokkos::create_mirror_view(num_ratios);
  for (int m = 0; m < nmodes; ++m) {
    for (int k = 0; k < nlev; ++k) {
      const int pack_idx = PackInfo::pack_idx(k);
      const int vec_idx = PackInfo::vec_idx(k);
      h_num_ratios(m, pack_idx)[vec_idx] = (2 * k + 1) * (m + 1) * 1e6;
    }
  }
  Kokkos::deep_copy(num_ratios, h_num_ratios);

  // Allocate a view to hold AerosolSpecies for each mode, individually
  DeviceType::view_2d<AerosolSpecies> aerosols_in_mode(
      "aerosols_in_mode", nmodes, config.max_species_per_mode());
  Kokkos::parallel_for(
      nmodes, KOKKOS_LAMBDA(const int m) {
        config.aerosol_species_for_mode(
            m, Kokkos::subview(aerosols_in_mode, m, Kokkos::ALL));
      });

  // modal averages
  ModalColumnView mode_wet_radius("mode_wet_radius", nmodes, npacks);
  ModalColumnView mode_mean_particle_dry_volume("mode_mean_particle_dry_volume",
                                                nmodes, npacks);
  ModalColumnView mode_hygroscopicity("mode_hygroscopicity", nmodes, npacks);
  ModalColumnView mode_dry_particle_radius("mode_dry_particle_radius", nmodes,
                                           npacks);

  std::cout << "starting kernels\n";
  // on host, loop over each mode
  for (int m = 0; m < nmodes; ++m) {
    // compute the modal mean particle volume
    Kokkos::parallel_for(
        npacks,
        ModalMeanParticleVolume(
            Kokkos::subview(mode_mean_particle_dry_volume, m, Kokkos::ALL),
            q_aero, num_ratios,
            Kokkos::subview(aerosols_in_mode, m, Kokkos::ALL),
            Kokkos::subview(pop_inds, m, Kokkos::ALL),
            config.h_n_species_per_mode(m), m));
    std::cout << "\tdry particle volume ready\n";

    // compute the modal mean hygroscopicity
    Kokkos::parallel_for(
        npacks, ModalHygroscopicity(
                    Kokkos::subview(mode_hygroscopicity, m, Kokkos::ALL),
                    q_aero, Kokkos::subview(aerosols_in_mode, m, Kokkos::ALL),
                    Kokkos::subview(pop_inds, m, Kokkos::ALL),
                    config.h_n_species_per_mode(m)));
    std::cout << "\thygroscopicity ready\n";

    // compute the modal mean dry particle *radius* (not diameter)
    Kokkos::parallel_for(
        npacks, KOKKOS_LAMBDA(const int pack_idx) {
          mode_dry_particle_radius(m, pack_idx) =
              0.5 *
              config.d_aerosol_modes[m].mean_particle_diameter_from_volume(
                  mode_mean_particle_dry_volume(m, pack_idx));
        });
    std::cout << "\tdry particle radius ready\n";

    // compute the modal avg wet radius
    Kokkos::parallel_for(
        npacks, ModalWetRadius(
                    Kokkos::subview(mode_wet_radius, m, Kokkos::ALL),
                    Kokkos::subview(mode_hygroscopicity, m, Kokkos::ALL),
                    Kokkos::subview(mode_dry_particle_radius, m, Kokkos::ALL),
                    relative_humidity, config.h_aerosol_modes[m]));
    std::cout << "\twet radius ready\n";
  }

  std::cout << "kernel evaluations complete.\n" << std::endl;

  // for verification: copy kernel results to host
  auto h_mode_hyg = Kokkos::create_mirror_view(mode_hygroscopicity);
  auto h_mode_dry_vol =
      Kokkos::create_mirror_view(mode_mean_particle_dry_volume);
  auto h_mode_dry_radius = Kokkos::create_mirror_view(mode_dry_particle_radius);
  Kokkos::deep_copy(h_mode_hyg, mode_hygroscopicity);
  Kokkos::deep_copy(h_mode_dry_vol, mode_mean_particle_dry_volume);
  Kokkos::deep_copy(h_mode_dry_radius, mode_dry_particle_radius);
  auto h_mode_wet_radius = Kokkos::create_mirror_view(mode_wet_radius);
  Kokkos::deep_copy(h_mode_wet_radius, mode_wet_radius);

  // Mode 0 has only 1 species, so its modal averages should match that single
  // species
  const auto mode0spec = config.aerosol_species_for_mode(0);
  for (int pack_idx = 0; pack_idx < npacks; ++pack_idx) {
    // check hygroscopicity matches exact value
    REQUIRE(FloatingPoint<PackType>::equiv(
        h_mode_hyg(0, pack_idx), PackType(mode0spec[0].hygroscopicity)));
    // check particle dry volume matches exact value
    REQUIRE(FloatingPoint<PackType>::equiv(h_mode_dry_vol(0, pack_idx),
                                           h_q_aero(0, pack_idx) /
                                               mode0spec[0].density /
                                               h_num_ratios(0, pack_idx)));

    // check that wet radius >= dry radius always
    REQUIRE(Check<PackType>::is_greater_or_equal(
        h_mode_wet_radius(0, pack_idx), h_mode_dry_radius(0, pack_idx)));
  }
  std::cout << "Tested mode 0 hygroscopicity == exact value\n";
  std::cout << "Tested mode 0 dry_volume == exact value \n";
  std::cout << "Tested mode 0 wet radius >= dry radius\n";

  // Mode 1 has 2 species, so we can easily compute the exact mass-weighted
  // modal averages
  const auto mode1spec = config.aerosol_species_for_mode(1);
  for (int pack_idx = 0; pack_idx < npacks; ++pack_idx) {
    const PackType mode_vol_mr = h_q_aero(1, pack_idx) / mode1spec[0].density +
                                 h_q_aero(2, pack_idx) / mode1spec[1].density;
    const PackType mhyg = (h_q_aero(1, pack_idx) * mode1spec[0].hygroscopicity /
                               mode1spec[0].density +
                           h_q_aero(2, pack_idx) * mode1spec[1].hygroscopicity /
                               mode1spec[1].density) /
                          mode_vol_mr;
    // check modal avg particle dry volume
    REQUIRE(FloatingPoint<PackType>::equiv(
        h_mode_dry_vol(1, pack_idx), mode_vol_mr / h_num_ratios(1, pack_idx)));
    // check modal avg hygroscopicity
    REQUIRE(FloatingPoint<PackType>::equiv(h_mode_hyg(1, pack_idx), mhyg));
    // check that wet radius >= dry radius always
    REQUIRE(Check<PackType>::is_greater_or_equal(
        h_mode_wet_radius(1, pack_idx), h_mode_dry_radius(1, pack_idx)));
  }

  std::cout << "Tested mode 1 hygroscopicity == exact value\n";
  std::cout << "Tested mode 1 dry_volume == exact value \n";
  std::cout << "Tested mode 1 wet radius >= dry radius\n";

  for (int m = 0; m < nmodes; ++m) {
    std::cout << "Mode " << m << ":\n";
    std::cout << "\t" << std::setw(HAERO_PACK_SIZE * 18) << "dry radius"
              << std::setw(HAERO_PACK_SIZE * 18) << "wet radius"
              << std::setw(HAERO_PACK_SIZE * 18) << " wet - dry\n";
    for (int pack_idx = 0; pack_idx < npacks; ++pack_idx) {
      std::cout << std::setw(HAERO_PACK_SIZE * 18)
                << h_mode_dry_radius(m, pack_idx)
                << std::setw(HAERO_PACK_SIZE * 18)
                << h_mode_wet_radius(m, pack_idx)
                << std::setw(HAERO_PACK_SIZE * 18)
                << h_mode_wet_radius(m, pack_idx) -
                       h_mode_dry_radius(m, pack_idx)
                << "\n";
    }
  }
}
