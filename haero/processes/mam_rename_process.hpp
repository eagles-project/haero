#ifndef HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP
#define HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP

#include "haero/aerosol_process.hpp"

namespace haero {

/// \brief Bindings for the rename subroutine
class MAMRenameProcess final : public AerosolProcess {
 public:
  MAMRenameProcess();

  KOKKOS_INLINE_FUNCTION
  ~MAMRenameProcess() {}

 protected:
  //------------------------------------------------------------------------
  //                                Overrides
  //------------------------------------------------------------------------

  void init_(const ModalAerosolConfig &modal_aerosol_config) override;

  KOKKOS_FUNCTION
  void run_(Real t, Real dt, const Prognostics &prognostics,
            const Atmosphere &atmosphere, const Diagnostics &diagnostics,
            Tendencies &tendencies) const override {}
};

}  // namespace haero

#endif
