#ifndef HAERO_MAM_RENAME_PROCESS_HPP
#define HAERO_MAM_RENAME_PROCESS_HPP

#include "haero/aerosol_process.hpp"
#include "haero_config.hpp"
#include "kokkos/Kokkos_Core.hpp"
#include "kokkos/Kokkos_Vector.hpp"

namespace haero {

/// \brief Bindings for the rename subroutine
class MAMRenameProcess final : public DeviceAerosolProcess<MAMRenameProcess> {
 public:
  template <typename T>
  using Container = Kokkos::vector<T>;

  using Integral = int;
  using Size = std::size_t;

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
            const Tendencies& tendencies) const override {}

 private:
  void find_renaming_pairs_(
      const ModalAerosolConfig& config,
      const Container<Integral>& dest_mode_of_mode_mapping, Size& num_pairs,
      Container<Real>& sz_factor, Container<Real>& fmode_dist_tail_fac,
      Container<Real>& v2n_lo_rlx, Container<Real>& v2n_hi_rlx,
      Container<Real>& ln_diameter_tail_fac, Container<Real>& diameter_cutoff,
      Container<Real>& ln_dia_cutoff, Container<Real>& diameter_belowcutoff,
      Container<Real>& dryvol_smallest) const;

 private:
  // TODO: These variable names are ambiguous and ought to be updated in another
  // PR. The Calcsize references many of these variables.
  Container<Real> dgnumlo;
  Container<Real> dgnumhi;
  Container<Real> dgnum;
  Container<Real> alnsg;
};

}  // namespace haero

#endif
