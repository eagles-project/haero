#ifndef HAERO_MAM_CALCSIZE_PROCESS_HPP
#define HAERO_MAM_CALCSIZE_PROCESS_HPP

#include "haero/aerosol_process.hpp"


namespace haero {

  class MAMCalcsizeProcess : public AerosolProcess {

  public:
    MAMCalcsizeProcess();

    /// Destructor.
    KOKKOS_INLINE_FUNCTION
    virtual ~MAMCalcsizeProcess() {}

    virtual void init(const ModalAerosolConfig &modal_aerosol_config) override;

    KOKKOS_FUNCTION
    virtual void run(const ModalAerosolConfig &modal_aerosol_config, Real t,
                     Real dt, const Prognostics &prognostics,
                     const Atmosphere &atmosphere, const Diagnostics &diagnostics,
                     Tendencies &tendencies) const override;

    };

}  // namespace haero

#endif
