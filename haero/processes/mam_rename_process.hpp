#ifndef HAERO_MAM_RENAME_PROCESS_HPP
#define HAERO_MAM_RENAME_PROCESS_HPP

#include "haero/aerosol_process.hpp"
#include "haero_config.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include "kokkos/Kokkos_Vector.hpp"

namespace haero {

namespace {

using haero::ColumnView;
using haero::SpeciesColumnView;

KOKKOS_INLINE_FUNCTION
void dryvolume_change(const TeamType& team, const int imode,
                      SpeciesColumnView q_vmr, SpeciesColumnView q_del_growth,
                      const int start_species_index,
                      const int end_species_index, ColumnView dryvol,
                      ColumnView deldryvol, ColumnView mass_2_vol) {
  printf("(%d,%d)\n", start_species_index, end_species_index);
  Kokkos::parallel_scan(
      Kokkos::TeamThreadRange(team, start_species_index, end_species_index),
      KOKKOS_LAMBDA(const std::size_t ispec, PackType& update,
                    const bool final_pass) {
        update += q_del_growth(imode, ispec) * mass_2_vol(ispec);
        if (final_pass) {
          deldryvol(ispec) = update;
        }
      });
  Kokkos::parallel_scan(
      Kokkos::TeamThreadRange(team, start_species_index, end_species_index),
      KOKKOS_LAMBDA(const std::size_t ispec, PackType& update,
                    const bool final_pass) {
        update += q_vmr(imode, ispec) * mass_2_vol(ispec);
        if (final_pass) {
          //segfault
          dryvol(ispec) = update - deldryvol(ispec);
        }
      });
}

KOKKOS_INLINE_FUNCTION
void compute_dryvol_change_in_src_mode(
    int num_modes, const view_1d_int_type dest_mode_of_mode_mapping,
    bool is_cloudy, const SpeciesColumnView qi_vmr,
    SpeciesColumnView qi_del_growth, const SpeciesColumnView qcld_vmr,
    SpeciesColumnView qcld_del_growth, ColumnView dryvol_a,
    ColumnView deldryvol_a, ColumnView dryvol_c, ColumnView deldryvol_c,
    view_1d_int_type population_offsets, int num_populations,
    const ColumnView mass_2_vol, view_1d_int_type start_species_for_mode,
    view_1d_int_type end_species_for_mode, view_1d_int_type stencil,
    const TeamType& team) {
  Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, num_modes),
      KOKKOS_LAMBDA(const int src_mode) {
        // Create mask such that functions below only run for source modes
        // that have a destination mode
        stencil(src_mode) =
            static_cast<int>(dest_mode_of_mode_mapping(src_mode) > 0);

        start_species_for_mode(src_mode) = population_offsets(src_mode);

        // find start and end index of species in this mode in the
        // "population" array The indices are same for interstitial and
        // cloudborne species
        end_species_for_mode(src_mode) = population_offsets(src_mode + 1);
      });

  for (int src_mode = 0; src_mode < num_modes; src_mode++) {
    if (!stencil(src_mode)) continue;
    dryvolume_change(
        team, src_mode, qi_vmr, qi_del_growth, start_species_for_mode(src_mode),
        end_species_for_mode(src_mode), dryvol_a, deldryvol_a, mass_2_vol);
  }

  if (is_cloudy) {
    for (int src_mode = 0; src_mode < num_modes; src_mode++) {
      if (!stencil(src_mode)) continue;
      dryvolume_change(team, src_mode, qcld_vmr, qcld_del_growth,
                       start_species_for_mode(src_mode),
                       end_species_for_mode(src_mode), dryvol_c, deldryvol_c,
                       mass_2_vol);
    }
  }
}
}  // namespace

/// \brief Bindings for the rename subroutine
class MAMRenameProcess final : public DeviceAerosolProcess<MAMRenameProcess> {
 public:
  MAMRenameProcess()
      : DeviceAerosolProcess<MAMRenameProcess>(RenameProcess,
                                               "MAMRenameProcess"),
        is_cloudy{true} {}

  KOKKOS_INLINE_FUNCTION
  ~MAMRenameProcess() {}

 protected:
  //------------------------------------------------------------------------
  //                                Overrides
  //------------------------------------------------------------------------

  void init_(const ModalAerosolConfig& modal_aerosol_config) override;

  KOKKOS_INLINE_FUNCTION
  void run_(const TeamType& team, Real t, Real dt,
            const Prognostics& prognostics, const Atmosphere& atmosphere,
            const Diagnostics& diagnostics,
            const Tendencies& tendencies) const override {
    const auto& qi_vmr = prognostics.interstitial_aerosols;
    const auto& qcld_vmr = prognostics.cloud_aerosols;

    Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, max_aer), KOKKOS_LAMBDA(const int i) {
          for (int j = 0; j < num_modes; j++) {
            qi_del_growth(j, i) = 0;
            qcld_del_growth(j, i) = 0;
          }
        });

    compute_dryvol_change_in_src_mode(
        num_modes, dest_mode_of_mode_mapping, is_cloudy, qi_vmr, qi_del_growth,
        qcld_vmr, qcld_del_growth, dryvol_a, deldryvol_a, dryvol_c, deldryvol_c,
        population_offsets, num_populations, mass_2_vol, start_species_for_mode,
        end_species_for_mode, stencil, team);
  }

 private:
  void find_renaming_pairs_(
      const ModalAerosolConfig& config,
      view_1d_int_type dest_mode_of_mode_mapping, std::size_t& num_pairs,
      view_1d_pack_type size_factor, view_1d_pack_type fmode_dist_tail_fac,
      view_1d_pack_type volume2num_lo_relaxed,
      view_1d_pack_type volume2num_hi_relaxed,
      view_1d_pack_type ln_diameter_tail_fac, view_1d_pack_type diameter_cutoff,
      view_1d_pack_type ln_dia_cutoff, view_1d_pack_type diameter_belowcutoff,
      view_1d_pack_type dryvol_smallest) const;

 private:
  // TODO: These variable names are ambiguous and ought to be updated in another
  // PR. The Calcsize references many of these variables.
  ColumnView dgnumlo;
  ColumnView dgnumhi;
  ColumnView dgnum;
  ColumnView alnsg;
  ColumnView mass_2_vol;

  // FIXME: naming
  SpeciesColumnView qi_del_growth;
  SpeciesColumnView qcld_del_growth;
  ColumnView dryvol_a;
  ColumnView dryvol_c;
  ColumnView deldryvol_a;
  ColumnView deldryvol_c;

  view_1d_int_type population_offsets;
  view_1d_int_type start_species_for_mode;
  view_1d_int_type end_species_for_mode;
  view_1d_int_type stencil;
  view_1d_int_type dest_mode_of_mode_mapping;

  int max_aer;
  int num_populations;
  int naer;
  int num_modes;

  bool is_cloudy;
};

}  // namespace haero

#endif
