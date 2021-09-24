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

template<std::size_t Rank, Kokkos::Iterate Layout=Kokkos::Iterate::Left>
using RangeHelper = Kokkos::MDRangePolicy<Kokkos::Rank<Rank, Layout>>;

  // integer,  intent(in):: num_modes ! total number of modes
  // integer,  intent(in):: dest_mode_of_mode(:) ! destination mode for a mode
  // logical,  intent(in) :: iscldy ! true if it is a cloudy cell
  // real(wp), intent(in) :: qi_vmr(:,:)           ! mass mixing ratios (mmr) [kmol/kmol]
  // real(wp), intent(in) :: qi_del_growth(:,:) !growth in mmr [kmol/kmol]
  // real(wp), intent(in), optional :: qcld_vmr(:,:)
  // real(wp), intent(in), optional :: qcld_del_growth(:,:)
  // !intent-outs
  // real(wp), intent(out) :: dryvol_a(:), dryvol_c(:)       !dry volumes (before growth) [m3/kmol-air]
  // real(wp), intent(out) :: deldryvol_a(:), deldryvol_c(:) !change in dry volumes [m3/kmol-air]
KOKKOS_INLINE_FUNCTION
void compute_dryvol_change_in_src_mode(int num_modes, int dest_mode_of_mode,
    bool is_cloudy, SpeciesColumnView qi_vmr, SpeciesColumnView qi_del_growth,
    SpeciesColumnView qcld_vmr, SpeciesColumnView qcld_del_growth,
    ColumnView dryvol_a, ColumnView deldryvol_a,
    ColumnView dryvol_c, ColumnView deldryvol_c) {
  // 
}

}

/// \brief Bindings for the rename subroutine
class MAMRenameProcess final : public DeviceAerosolProcess<MAMRenameProcess> {
 public:

  using integral_type = int;
  using size_type = std::size_t;

  MAMRenameProcess();

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

    auto qi_vmr = prognostics.interstitial_aerosols();
    auto qcld_vmr = prognostics.cloud_aerosols();

    // FIXME: naming
    SpeciesColumnView qi_del_growth("intersitial growth", max_aer, num_modes);
    SpeciesColumnView qcld_del_growth("cloudborne growth", max_aer, num_modes);

    Kokkos::parallel_for(RangeHelper<2>({{0, 0}}, {{max_aer, num_modes}}),
        KOKKOS_LAMBDA(const int i, const int j) {
          qi_del_growth(i, j) = 0;
          qcld_del_growth(i, j) = 0;
        });

    // have:
    // - dest_mode_of_mode_mapping
    // - num_modes
    // - dgnumlo;
    // - dgnumhi;
    // - dgnum;
    // - alnsg;
    // - population_offsets;
    // - max_aer;
    // - num_populations;
    // - naer;
    // - num_modes;

    // compute_dryvol_change_in_src_mode(num_modes, dest_mode_of_mode_mapping, is_cloudy, qi_vmr, qi_del_growth)

    // Callsite -- dont use
    // call compute_dryvol_change_in_src_mode(num_modes, dest_mode_of_mode, &              !input
    //     iscldy_subarea, q_interstitial, qaer_del_grow4rnam, q_cloudborne, qaercw_del_grow4rnam, & !input
    //     dryvol_a, deldryvol_a, dryvol_c, deldryvol_c)                                     !output
    //
    // Def site -- use these names
    // compute_dryvol_change_in_src_mode(
    //     num_modes, dest_mode_of_mode, iscldy, qi_vmr, qi_del_growth, qcld_vmr,
    //     qcld_del_growth, dryvol_a, deldryvol_a, dryvol_c, deldryvol_c);
  }

 private:
  void find_renaming_pairs_(const ModalAerosolConfig& config,
                            view_1d_int_type& dest_mode_of_mode_mapping,
                            size_type& num_pairs, view_1d_pack_type& size_factor,
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
  ColumnView population_offsets;

  integral_type max_aer;
  integral_type num_populations;
  integral_type naer;
  integral_type num_modes;

  view_1d_int_type dest_mode_of_mode_mapping;

  bool is_cloudy;
};

}  // namespace haero

#endif
