#include "haero/processes/mam_rename_subarea_process.hpp"

namespace haero {

  MAMRenameSubareaProcess::MAMRenameSubareaProcess()
    : AerosolProcess(RenameSubareaProcess, "MAMRenameSubareaProcess") {}

  void MAMCalcsizeProcess::init(
                                  const ModalAerosolConfig& modal_aerosol_config) {}

  KOKKOS_FUNCTION
  void MAMCalcsizeProcess::run(const ModalAerosolConfig& modal_aerosol_config,
                               Real t, Real dt, const Prognostics& prognostics,
                               const Atmosphere& atmosphere,
                               const Diagnostics& diagnostics,
                               Tendencies& tendencies) const {}

} // namespace haero
