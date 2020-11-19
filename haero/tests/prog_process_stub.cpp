#include "haero/model.hpp"
#include "prog_process_stub.hpp"

using haero::Real;

extern "C" {

// Fortran creation function.
extern void prog_stub_init();
extern void prog_stub_run(Real t, Real dt, void* progs, void* diags, void* tends);
extern void prog_stub_finalize();

} // extern "C"

namespace haero {

ProgProcessStub::ProgProcessStub():
  FPrognosticProcess(haero::ActivationProcess, "Prognostic process stub (Fortran)",
                     prog_stub_init, prog_stub_run, prog_stub_finalize) {
}

} // end haero namespace

