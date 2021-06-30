#ifndef HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP
#define HAERO_MAM_RENAME_SUBAREA_PROCESS_HPP

#include "haero/aerosol_process.hpp"

namespace haero {

  /// \brief Bindings for the rename_subarea subroutine
  class MAMRenameSubareaProcess : public AerosolProcess {
  public:
    MAMRenameSubareaProcess();

    KOKKOS_INLINE_FUNCTION
    virtual ~MAMRenameSubareaProcess() {}

    virtual void init(const ModalAerosolConfig &modal_aerosol_config) override;

    KOKKOS_FUNCTION
    virtual void run(const ModalAerosolConfig &modal_aerosol_config, Real t,
                     Real dt, const Prognostics &prognostics,
                     const Atmosphere &atmosphere, const Diagnostics &diagnostics,
                     Tendencies &tendencies) const override;
  };

} // namespace haero

#endif
