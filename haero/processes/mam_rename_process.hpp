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

template <std::size_t Rank, Kokkos::Iterate Layout = Kokkos::Iterate::Left>
using RangeHelper = Kokkos::MDRangePolicy<Kokkos::Rank<Rank, Layout>>;

KOKKOS_INLINE_FUNCTION
void dryvolume_change(const int imode, SpeciesColumnView q_vmr,
                      SpeciesColumnView q_del_growth,
                      const int start_species_index,
                      const int end_species_index, PackType& dryvol,
                      PackType& deldryvol, const ColumnView& mass_2_vol) {
  PackType tmp_dryvol{0.}, tmp_del_dryvol{0.};
  for (int ispec = start_species_index; ispec < end_species_index; ispec++) {
    tmp_dryvol = tmp_dryvol + q_vmr(ispec, imode) * mass_2_vol(ispec);
    tmp_del_dryvol =
        tmp_del_dryvol + q_del_growth(ispec, imode) * mass_2_vol(ispec);
  }

  dryvol = tmp_dryvol - tmp_del_dryvol;
  deldryvol = tmp_del_dryvol;
}

KOKKOS_INLINE_FUNCTION
void compute_dryvol_change_in_src_mode(
    int num_modes, const view_1d_int_type& dest_mode_of_mode_mapping,
    bool is_cloudy, SpeciesColumnView qi_vmr, SpeciesColumnView qi_del_growth,
    SpeciesColumnView qcld_vmr, SpeciesColumnView qcld_del_growth,
    ColumnView& dryvol_a, ColumnView& deldryvol_a, ColumnView& dryvol_c,
    ColumnView& deldryvol_c, view_1d_int_type population_offsets,
    int num_populations, const ColumnView& mass_2_vol) {
  // FIXME: This array seems to only be assigned to once, and it's taken
  // directly from population_offsets. Will this ever be calculated differently?
  view_1d_int_type start_species_for_mode("start_species_for_mode", num_modes);
  Kokkos::deep_copy(start_species_for_mode, population_offsets);

  view_1d_int_type end_species_for_mode("end_species_for_mode", num_modes);
  Kokkos::parallel_for(
      num_modes, KOKKOS_LAMBDA(const int src_mode) {
        // find start and end index of species in this mode in the "population"
        // array The indices are same for interstitial and cloudborne species
        end_species_for_mode(src_mode) =
            (src_mode == num_modes) ? num_populations
                                    : population_offsets(src_mode + 1) - 1;
      });

  Kokkos::parallel_for(
      num_modes, KOKKOS_LAMBDA(const int src_mode) {
        const auto dest_mode = dest_mode_of_mode_mapping(src_mode);

        // FIXME: How sparse is dest_mode_of_mode_mapping? If it's sparse,
        // should we generate tuples of source/destination modes as a
        // first-pass? This sort of early return could kill performance if it
        // happens often.
        if (dest_mode <= 0) return;

        dryvolume_change(src_mode, qi_vmr, qi_del_growth,
                         start_species_for_mode(src_mode),
                         end_species_for_mode(src_mode), dryvol_a(src_mode),
                         deldryvol_a(src_mode), mass_2_vol);
      });

  if (is_cloudy) {
    Kokkos::parallel_for(
        num_modes, KOKKOS_LAMBDA(const int src_mode) {
          const auto dest_mode = dest_mode_of_mode_mapping(src_mode);

          // FIXME: see above
          if (dest_mode <= 0) return;

          dryvolume_change(src_mode, qi_vmr, qi_del_growth,
                           start_species_for_mode(src_mode),
                           end_species_for_mode(src_mode), dryvol_c(src_mode),
                           deldryvol_c(src_mode), mass_2_vol);
        });
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

    // FIXME: naming
    SpeciesColumnView qi_del_growth("intersitial growth", max_aer, num_modes);
    SpeciesColumnView qcld_del_growth("cloudborne growth", max_aer, num_modes);

    Kokkos::parallel_for(
        RangeHelper<2>({{0, 0}}, {{max_aer, num_modes}}),
        KOKKOS_LAMBDA(const int i, const int j) {
          qi_del_growth(i, j) = 0;
          qcld_del_growth(i, j) = 0;
        });

    ColumnView dryvol_a("dryvol_a", num_modes);
    ColumnView dryvol_c("dryvol_c", num_modes);
    ColumnView deldryvol_a("deldryvol_a", num_modes);
    ColumnView deldryvol_c("deldryvol_c", num_modes);

    compute_dryvol_change_in_src_mode(
        num_modes, dest_mode_of_mode_mapping, is_cloudy, qi_vmr, qi_del_growth,
        qcld_vmr, qcld_del_growth, dryvol_a, deldryvol_a, dryvol_c, deldryvol_c,
        population_offsets, num_populations, mass_2_vol);
  }

 private:
  void find_renaming_pairs_(const ModalAerosolConfig& config,
                            view_1d_int_type& dest_mode_of_mode_mapping,
                            std::size_t& num_pairs,
                            view_1d_pack_type& size_factor,
                            view_1d_pack_type& fmode_dist_tail_fac,
                            view_1d_pack_type& volume2num_lo_relaxed,
                            view_1d_pack_type& volume2num_hi_relaxed,
                            view_1d_pack_type& ln_diameter_tail_fac,
                            view_1d_pack_type& diameter_cutoff,
                            view_1d_pack_type& ln_dia_cutoff,
                            view_1d_pack_type& diameter_belowcutoff,
                            view_1d_pack_type& dryvol_smallest) const;

 private:
  // TODO: These variable names are ambiguous and ought to be updated in another
  // PR. The Calcsize references many of these variables.
  ColumnView dgnumlo;
  ColumnView dgnumhi;
  ColumnView dgnum;
  ColumnView alnsg;
  ColumnView mass_2_vol;
  view_1d_int_type population_offsets;

  int max_aer;
  int num_populations;
  int naer;
  int num_modes;

  view_1d_int_type dest_mode_of_mode_mapping;

  bool is_cloudy;
};

}  // namespace haero

#endif
