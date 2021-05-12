#include "faerosol_process_stub.hpp"

#include "haero/model.hpp"

using haero::Real;

extern "C" {

// Fortran subroutines that implement this process.
extern void process_stub_init();
extern void process_stub_run(Real t, Real dt, void* progs, void* atm,
                             void* diags, void* tends);
extern void process_stub_finalize();

// Decay rate. Since Fortran-backed processes cannot be run in more than one
// instance of a model, we can safely use this as a global variable.
Real decay_rate_;

// A C function that the Fortran implementation can call to obtain the given
// exponential decay rate.
Real process_stub_decay_rate() { return decay_rate_; }

}  // extern "C"

namespace haero {

FAerosolProcessStub::FAerosolProcessStub(Real decay_rate)
    : FAerosolProcess(haero::ActivationProcess,
                      "Aerosol process stub (Fortran)", process_stub_init,
                      process_stub_run, process_stub_finalize) {
  // No positive decay rates!
  EKAT_ASSERT_MSG(decay_rate <= 0.0, "decay_rate must be non-positive!");

  // Set the decay rate global variable to make it available to Fortran.
  decay_rate_ = decay_rate;
}

}  // namespace haero
