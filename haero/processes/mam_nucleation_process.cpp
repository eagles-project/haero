
#include "haero/processes/mam_nucleation_process.hpp"

#include <cmath>

namespace haero {

MAMNucleationProcess::MAMNucleationProcess()
    : AerosolProcess(NucleationProcess, "MAMNucleationProcess") {}

//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMNucleationProcess::init(
    const ModalAerosolConfig& modal_aerosol_config) {}

}  // namespace haero
