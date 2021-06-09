#include "haero/processes/mam_calcsize_process.hpp"

#include <cmath>

namespace haero {

  MAMCalcsizeProcess::MAMCalcsizeProcess()
    : AerosolProcess(CalcsizeProcess, "MAMCalcsizeProcess") {}

  //------------------------------------------------------------------------
  //                                Accessors
  //------------------------------------------------------------------------

  void MAMCalcsizeProcess::init(
                                  const ModalAerosolConfig& modal_aerosol_config) {}

  KOKKOS_FUNCTION
  void MAMCalcsizeProcess::run(const ModalAerosolConfig& modal_aerosol_config,
                               Real t, Real dt, const Prognostics& prognostics,
                               const Atmosphere& atmosphere,
                               const Diagnostics& diagnostics,
                               Tendencies& tendencies) const {}

}  // namespace haero
