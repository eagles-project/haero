#include "haero/processes/mam_rename_process.hpp"

#include <iostream>

namespace haero {

MAMRenameProcess::MAMRenameProcess()
    : DeviceAerosolProcess<MAMRenameProcess>(RenameProcess,
                                             "MAMRenameProcess") {}

void MAMRenameProcess::init_(const ModalAerosolConfig& config) {}

}  // namespace haero
