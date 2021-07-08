#include "haero/processes/mam_rename_process.hpp"
#include <iostream>

namespace haero {

  MAMRenameProcess::MAMRenameProcess()
    : AerosolProcess(RenameProcess, "MAMRenameProcess") {}

  void MAMRenameProcess::init(
                                  const ModalAerosolConfig& modal_aerosol_config) {}

  KOKKOS_FUNCTION
  void MAMRenameProcess::run(const ModalAerosolConfig& modal_aerosol_config,
                               Real t, Real dt, const Prognostics& prognostics,
                               const Atmosphere& atmosphere,
                               const Diagnostics& diagnostics,
                               Tendencies& tendencies) const {
    std::cout << "Running in here\n";
  }

} // namespace haero
