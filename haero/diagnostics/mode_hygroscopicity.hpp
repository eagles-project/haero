#ifndef HAERO_MODAL_MEAN_HYGROSCOPICITY_HPP
#define HAERO_MODAL_MEAN_HYGROSCOPICITY_HPP

#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "haero/aerosol_species.hpp"
#include "haero/haero.hpp"

namespace haero {

struct ModalHygroscopicity {
  ColumnView hygroscopicity;
  SpeciesColumnView mass_mixing_ratios;
  DeviceType::view_1d<const AerosolSpecies> aerosols_in_mode;
  DeviceType::view_1d<const int> population_indices;
  int n_aerosols_in_mode;

  KOKKOS_INLINE_FUNCTION
  ModalHygroscopicity(ColumnView hyg, const SpeciesColumnView mass_mrs,
                      const DeviceType::view_1d<AerosolSpecies> aeros,
                      const DeviceType::view_1d<int> pops, const int n)
      : hygroscopicity(hyg),
        mass_mixing_ratios(mass_mrs),
        aerosols_in_mode(aeros),
        population_indices(pops),
        n_aerosols_in_mode(n) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int pack_idx) const {
    PackType hyg(0);
    PackType vol_mix_ratio(0);
    for (int i = 0; i < n_aerosols_in_mode; ++i) {
      const int s = population_indices(i);
      hyg += mass_mixing_ratios(s, pack_idx) *
             aerosols_in_mode(i).hygroscopicity / aerosols_in_mode(i).density;
      vol_mix_ratio +=
          mass_mixing_ratios(s, pack_idx) / aerosols_in_mode(i).density;
    }
    hygroscopicity(pack_idx) = hyg / vol_mix_ratio;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename DeviceType::MemberType& column_team) const {
    const int pack_idx = column_team.team_rank();
    PackType hyg(0);
    PackType vol_mix_ratio(0);
    const auto thread_policy =
        Kokkos::TeamThreadRange(column_team, n_aerosols_in_mode);
    Kokkos::parallel_reduce(
        thread_policy,
        [=](const int i, PackType& update) {
          const int s = population_indices(i);
          update += mass_mixing_ratios(s, pack_idx) *
                    aerosols_in_mode(i).hygroscopicity /
                    aerosols_in_mode(i).density;
        },
        hyg);
    Kokkos::parallel_reduce(
        thread_policy,
        [=](const int i, PackType& update) {
          const int s = population_indices(i);
          update +=
              mass_mixing_ratios(s, pack_idx) / aerosols_in_mode(i).density;
        },
        vol_mix_ratio);
    column_team.team_barrier();
    hygroscopicity(pack_idx) = hyg / vol_mix_ratio;
  }
};

}  // namespace haero
#endif
