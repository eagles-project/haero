#include "diag_fprocess_stub.hpp"

#include "haero/model.hpp"

using haero::Real;

extern "C" {

// Fortran creation function.
extern void diag_stub_init();
extern void diag_stub_update(Real t, void* progs, void* atm, void* diags);
extern void diag_stub_finalize();

}  // extern "C"

namespace haero {

DiagFProcessStub::DiagFProcessStub()
    : FDiagnosticProcess(haero::WaterUptakeProcess,
                         "Diagnostic process stub (Fortran)",
                         {"temperature"},  // atm variables
                         {},  // modal aerosol-species-specific variables
                         {"pressure"},  // gas-species-specific variables
                         {"pressure"},  // mode-specific variables
                         diag_stub_init, diag_stub_update, diag_stub_finalize) {
}

}  // namespace haero
