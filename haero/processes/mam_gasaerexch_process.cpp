
#include "haero/processes/mam_gasaerexch_process.hpp"

#include <cmath>

namespace haero {

MAMGasAerosolExchangeProcess::MAMGasAerosolExchangeProcess()
    : DeviceAerosolProcess(NucleationProcess, "MAMGasAerosolExchangeProcess") {}

//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMGasAerosolExchangeProcess::init_(
    const ModalAerosolConfig& modal_aerosol_config) {}

}  // namespace haero
