
#include <cmath>

#include "haero/processes/mam_nucleation_process.hpp"


namespace haero {


MAMNucleationProcess::MAMNucleationProcess():
  AerosolProcess(NucleationProcess, "MAMNucleationProcess") {}


//------------------------------------------------------------------------
//                                Accessors
//------------------------------------------------------------------------

void MAMNucleationProcess::init(const ModalAerosolConfig& modal_aerosol_config) {}

KOKKOS_FUNCTION
void MAMNucleationProcess::run(const ModalAerosolConfig& modal_aerosol_config,
                   Real t, Real dt,
                   const Prognostics& prognostics,
                   const Atmosphere& atmosphere,
                   const Diagnostics& diagnostics,
                   Tendencies& tendencies) const {}

}
