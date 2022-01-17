#include "haero/processes/mam_calcsize_process.hpp"

#include <Kokkos_Core.hpp>
#include <cmath>

namespace haero {

MAMCalcsizeProcess::MAMCalcsizeProcess()
    : DeviceAerosolProcess<MAMCalcsizeProcess>("MAMCalcsizeProcess") {}
//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMCalcsizeProcess::init_(const ModalAerosolConfig &modal_aerosol_config) {
  nmodes = modal_aerosol_config.num_aerosol_modes();

  num_populations = modal_aerosol_config.num_aerosol_populations;

  Kokkos::resize(num_mode_species, nmodes);
  {
    auto nm_spec = Kokkos::create_mirror_view(num_mode_species);
    for (int i = 0; i < nmodes; i++)
      nm_spec[i] = modal_aerosol_config.aerosol_species_for_mode(i).size();

    max_nspec = nm_spec(0);
    Kokkos::resize(density, max_nspec);

    for (int i = 0; i < nm_spec.extent(0); i++)
      max_nspec = std::max(max_nspec, static_cast<std::size_t>(nm_spec(i)));

    Kokkos::deep_copy(num_mode_species, nm_spec);
  }

  Kokkos::resize(spec_density, max_nspec * nmodes);
  {
    const auto h_spec_density = Kokkos::create_mirror_view(spec_density);
    const auto all_species = modal_aerosol_config.aerosol_species;
    for (int i = 0; i < all_species.size(); i++)
      h_spec_density(i) = all_species[i].density;
    Kokkos::deep_copy(spec_density, h_spec_density);
  }

  Kokkos::resize(population_offsets, nmodes);
  {
    const auto pop_offsets = Kokkos::create_mirror_view(population_offsets);
    for (int i = 0; i < nmodes; i++)
      pop_offsets[i] = modal_aerosol_config.population_index(i, 0);
    Kokkos::deep_copy(population_offsets, pop_offsets);
  }

  Kokkos::resize(v2nmin_nmodes, nmodes);
  Kokkos::resize(v2nmax_nmodes, nmodes);
  Kokkos::resize(dgnmin_nmodes, nmodes);
  Kokkos::resize(dgnmax_nmodes, nmodes);
  Kokkos::resize(common_factor_nmodes, nmodes);

  {
    const auto h_v2nmin_nmodes = Kokkos::create_mirror_view(v2nmin_nmodes);
    const auto h_v2nmax_nmodes = Kokkos::create_mirror_view(v2nmax_nmodes);
    const auto h_dgnmin_nmodes = Kokkos::create_mirror_view(dgnmin_nmodes);
    const auto h_dgnmax_nmodes = Kokkos::create_mirror_view(dgnmax_nmodes);
    const auto h_common_factor_nmodes =
        Kokkos::create_mirror_view(common_factor_nmodes);

    for (int i = 0; i < nmodes; i++) {
      const auto &mode = modal_aerosol_config.aerosol_modes[i];
      using T = decltype(v2nmin_nmodes)::value_type;
      h_v2nmin_nmodes[i] = mode.min_vol_to_num_ratio<T>();
      h_v2nmax_nmodes[i] = mode.max_vol_to_num_ratio<T>();
      h_dgnmin_nmodes[i] = mode.min_diameter;
      h_dgnmax_nmodes[i] = mode.max_diameter;
      h_common_factor_nmodes[i] =
          std::exp(4.5 * std::log(std::pow(mode.mean_std_dev, 2.0))) *
          Constants::pi_sixth;
    }

    Kokkos::deep_copy(v2nmin_nmodes, h_v2nmin_nmodes);
    Kokkos::deep_copy(v2nmax_nmodes, h_v2nmax_nmodes);
    Kokkos::deep_copy(dgnmin_nmodes, h_dgnmin_nmodes);
    Kokkos::deep_copy(dgnmax_nmodes, h_dgnmax_nmodes);
    Kokkos::deep_copy(common_factor_nmodes, h_common_factor_nmodes);
  }

  aitken_idx = modal_aerosol_config.aerosol_mode_index("aitken");
  accum_idx = modal_aerosol_config.aerosol_mode_index("accum");
}

}  // namespace haero
