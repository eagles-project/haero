#ifndef HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP
#define HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP

#include "haero/aerosol_process.hpp"

namespace haero
{

/// \brief Bindings for the rename subroutine
class MAMRenameProcess final : public DeviceAerosolProcess<MAMRenameProcess> {

public:
  template <typename T>
  using Container = std::vector<T>;

  using Integral = int;
  using Size = std::size_t;

  MAMRenameProcess();

  KOKKOS_INLINE_FUNCTION
  ~MAMRenameProcess() {}

 protected:
  //------------------------------------------------------------------------
  //                                Overrides
  //------------------------------------------------------------------------

  virtual void init_(const ModalAerosolConfig &modal_aerosol_config) override;

  KOKKOS_FUNCTION
  virtual void run_(Real t, Real dt, const Prognostics &prognostics,
            const Atmosphere &atmosphere, const Diagnostics &diagnostics,
            const Tendencies &tendencies) const override;

private:
  void find_renaming_pairs_(
      Size nmodes,
      const Container<Integral>& dest_mode_of_mode_mapping,
      Size& num_pairs,
      Container<Real>& sz_factor,
      Container<Real>& fmode_dist_tail_fac,
      Container<Real>& v2n_lo_rlx,
      Container<Real>& v2n_hi_rlx,
      Container<Real>& ln_diameter_tail_fac,
      Container<Real>& diameter_cutoff,
      Container<Real>& ln_dia_cutoff,
      Container<Real>& diameter_belowcutoff,
      Container<Real>& dryvol_smallest) const;

private:
  std::vector<Real> dgnumlo;
  std::vector<Real> dgnumhi;
  std::vector<Real> dgnum;
  std::vector<Real> alnsg;
};

}  // namespace haero

#endif
