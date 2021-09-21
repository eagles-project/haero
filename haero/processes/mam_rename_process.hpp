#ifndef HAERO_MAM_RENAME_PROCESS_HPP
#define HAERO_MAM_RENAME_PROCESS_HPP

#include "haero/aerosol_process.hpp"
#include "haero_config.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include "kokkos/Kokkos_Vector.hpp"

namespace {

// KOKKOS_INLINE_FUNCTION
// void compute_dryvol_change_in_src_mode(int num_modes, int dest_mode_of_mode,
//     bool is_cloudy, qi_vmr, qi_del_growth, qcld_vmr, qcld_del_growth,
//     dryvol_a, deldryvol_a, dryvol_c, deldryvol_c) {
// 
// }

}

namespace haero {

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
  view_1d_pack_type dgnumlo;
  view_1d_pack_type dgnumhi;
  view_1d_pack_type dgnum;
  view_1d_pack_type alnsg;
  view_1d_pack_type population_offsets;

  integral_type max_aer;
  integral_type num_populations;
  integral_type naer;
  integral_type num_modes;
};

}  // namespace haero

#endif
