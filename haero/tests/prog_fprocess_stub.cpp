#include "haero/model.hpp"
#include "prog_fprocess_stub.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines that implement this process.
extern void prog_stub_init();
extern void prog_stub_run(Real t, Real dt, void* progs, void* atm, void* diags, void* tends);
extern void prog_stub_finalize();

// Decay rate. Since Fortran-backed processes cannot be run in more than one
// instance of a model, we can safely use this as a global variable.
Real decay_rate_;

// A C function that the Fortran implementation can call to obtain the given
// exponential decay rate.
Real prog_stub_decay_rate() {
  return decay_rate_;
}

} // extern "C"

namespace haero {

ProgFProcessStub::ProgFProcessStub(Real decay_rate):
  FPrognosticProcess(haero::ActivationProcess, "Prognostic process stub (Fortran)",
                     prog_stub_init, prog_stub_run, prog_stub_finalize) {
  // No positive decay rates!
  EKAT_ASSERT_MSG(decay_rate <= 0.0, "decay_rate must be non-positive!");

  // Set the decay rate global variable to make it available to Fortran.
  decay_rate_ = decay_rate;
}

} // end haero namespace

