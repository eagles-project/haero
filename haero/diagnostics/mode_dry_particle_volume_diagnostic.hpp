#ifndef HAERO_MODE_DRY_PARTICLE_VOLUME_DIAGNOSTIC_HPP
#define HAERO_MODE_DRY_PARTICLE_VOLUME_DIAGNOSTIC_HPP

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "haero/aerosol_species.hpp"
#include "haero/haero.hpp"

namespace haero {

/** @brief This struct computes a reduction across the species in a mode to find
  the modal mean particle volume.

  The mean particle volume is proportional to the 3rd moment of the mode's
  log-normal PDF.
*/
struct ModalMeanParticleVolume {
  ColumnView modal_mean_volume;                                // output
  SpeciesColumnView mass_mixing_ratios;                        // input
  ModalColumnView number_mixing_ratios;                        // input
  DeviceType::view_1d<const AerosolSpecies> aerosols_in_mode;  // input
  DeviceType::view_1d<const int> population_indices;           // input
  int n_aerosols_in_mode;                                      // input
  int mode_idx;                                                // input

  KOKKOS_INLINE_FUNCTION
  ModalMeanParticleVolume(ColumnView vol, const SpeciesColumnView mass_mrs,
                          const ModalColumnView number_mrs,
                          DeviceType::view_1d<AerosolSpecies> aeros,
                          DeviceType::view_1d<int> pops, const int n,
                          const int m)
      : modal_mean_volume(vol),
        mass_mixing_ratios(mass_mrs),
        number_mixing_ratios(number_mrs),
        aerosols_in_mode(aeros),
        population_indices(pops),
        n_aerosols_in_mode(n),
        mode_idx(m) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int pack_idx) const {
    PackType vol_mixing_ratio(0);  // m3 aerosol per kg air
    for (int i = 0; i < n_aerosols_in_mode; ++i) {
      const int s = population_indices(i);
      vol_mixing_ratio +=
          mass_mixing_ratios(s, pack_idx) / aerosols_in_mode(i).density;
    }
    modal_mean_volume(pack_idx) =
        vol_mixing_ratio / number_mixing_ratios(mode_idx, pack_idx);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename DeviceType::MemberType& column_team) const {
    // const int col_idx = column_team.league_rank();
    const int pack_idx = column_team.team_rank();
    PackType vol_mixing_ratio(0);
    Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(column_team, n_aerosols_in_mode),
        [=](const int& i, PackType& update) {
          const int s = population_indices(i);
          update +=
              mass_mixing_ratios(s, pack_idx) / aerosols_in_mode(i).density;
        },
        vol_mixing_ratio);
    modal_mean_volume(pack_idx) =
        vol_mixing_ratio / number_mixing_ratios(mode_idx, pack_idx);
  }
};

}  // namespace haero
#endif
