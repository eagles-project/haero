#include "haero/model.hpp"
#include "mam4_nucleation_fprocess.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines that implement this process.
extern void mam4_nucleation_init();
extern void mam4_nucleation_run(Real t, Real dt, void* progs, void* atm, void* diags, void* tends);
extern void mam4_nucleation_finalize();

} // extern "C"

namespace haero {

Mam4NucleationFProcess::Mam4NucleationFProcess():
  FPrognosticProcess(haero::NucleationProcess, "MAM4 Nucleation (Fortran)",
                     mam4_nucleation_init, mam4_nucleation_run,
                     mam4_nucleation_finalize) {
}

} // end haero namespace

