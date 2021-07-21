#include "haero/processes/mam_rename_process.hpp"

#include <iostream>

namespace haero {

MAMRenameProcess::MAMRenameProcess()
    : AerosolProcess(RenameProcess, "MAMRenameProcess") {}

void MAMRenameProcess::init_(const ModalAerosolConfig& modal_aerosol_config) {}

}  // namespace haero
