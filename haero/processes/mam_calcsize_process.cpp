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

}  // namespace haero
