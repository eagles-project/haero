#include "haero/processes/mam_rename_subarea_process.hpp"
#include <iostream>

namespace haero {

  MAMRenameSubareaProcess::MAMRenameSubareaProcess()
    : AerosolProcess(RenameSubareaProcess, "MAMRenameSubareaProcess") {}

  void MAMRenameSubareaProcess::init(
                                  const ModalAerosolConfig& modal_aerosol_config) {}

  KOKKOS_FUNCTION
  void MAMRenameSubareaProcess::run(const ModalAerosolConfig& modal_aerosol_config,
                               Real t, Real dt, const Prognostics& prognostics,
                               const Atmosphere& atmosphere,
                               const Diagnostics& diagnostics,
                               Tendencies& tendencies) const {
    std::cout << "Running in here\n";
  }

} // namespace haero
