#ifndef HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP
#define HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP

#include "haero/aerosol_process.hpp"

namespace haero
{

/// \brief Bindings for the rename subroutine
class MAMRenameProcess final : public DeviceAerosolProcess<MAMRenameProcess> {
 public:
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

  void find_renaming_pairs_(const std::size_t nmodes,
                            const std::vector<Real>& dest_mode_of_mode,
                            const std::size_t num_pairs,
                            const std::vector<Real>& sz_factor,
                            const std::vector<Real>& fmode_dist_tail_fac,
                            const std::vector<Real>& v2n_lo_rlx,
                            const std::vector<Real>& v2n_hi_rlx,
                            const std::vector<Real>& ln_diameter_tail_fac,
                            const std::vector<Real>& diameter_cutoff,
                            const std::vector<Real>& ln_dia_cutoff,
                            const std::vector<Real>& diameter_belowcutoff,
                            const std::vector<Real>& dryvol_smallest) const;

private:

  std::vector<Real> dgnumlo;
  std::vector<Real> dgnumhi;
  std::vector<Real> dgnum;
  std::vector<Real> alnsg;
};

}  // namespace haero

#endif
