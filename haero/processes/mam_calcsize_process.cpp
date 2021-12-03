#include "haero/processes/mam_calcsize_process.hpp"

#include <cmath>

namespace haero {

MAMCalcsizeProcess::MAMCalcsizeProcess()
    : DeviceAerosolProcess<MAMCalcsizeProcess>("MAMCalcsizeProcess") {}
//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMCalcsizeProcess::init_(const ModalAerosolConfig& modal_aerosol_config) {
}

}  // namespace haero
