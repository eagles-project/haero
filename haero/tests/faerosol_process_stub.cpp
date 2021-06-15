#include "faerosol_process_stub.hpp"

#include "haero/model.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines that implement this process.
extern void process_stub_init();
extern void process_stub_run(Real t, Real dt, void* progs, void* atm,
                             void* diags, void* tends);
extern void process_stub_finalize();
extern void process_stub_set_integer_param(const char* name, int value);
extern void process_stub_set_logical_param(const char* name, bool value);
extern void process_stub_set_real_param(const char* name, Real value);

}  // extern "C"

namespace haero {

FAerosolProcessStub::FAerosolProcessStub()
    : FAerosolProcess(haero::ActivationProcess,
                      "Aerosol process stub (Fortran)", process_stub_init,
                      process_stub_run, process_stub_finalize,
                      process_stub_set_integer_param,
                      process_stub_set_logical_param,
                      process_stub_set_real_param) {
}

}  // namespace haero
