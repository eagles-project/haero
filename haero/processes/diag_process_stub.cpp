#include "haero/model.hpp"
#include "haero/processes/diag_process_stub.hpp"

using haero::Real;

extern "C" {

// Fortran creation function.
extern void diag_stub_init();
extern void diag_stub_update(Real t, void* progs, void* diags);
extern void diag_stub_finalize();

} // extern "C"

namespace haero {

DiagProcessStub::DiagProcessStub():
  FDiagnosticProcess(haero::WaterUptakeProcess, "Diagnostic process stub (Fortran)",
                     {}, {}, diag_stub_init, diag_stub_update, diag_stub_finalize) {
}

} // end haero namespace

