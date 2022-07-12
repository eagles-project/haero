
#include "haero/processes/mam_nucleation_process.hpp"

#include <cmath>

namespace haero {

MAMNucleationProcess::MAMNucleationProcess()
    : DeviceAerosolProcess<MAMNucleationProcess>("MAMNucleationProcess") {}

//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMNucleationProcess::init_(
    const ModalAerosolConfig& modal_aerosol_config) {}

}  // namespace haero
