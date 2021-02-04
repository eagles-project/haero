#ifndef HAERO_PROG_FPROCESS_STUB_HPP
#define HAERO_PROG_FPROCESS_STUB_HPP

#include "haero/process.hpp"

namespace haero {

/// This prognostic process implements a simple exponential decay model for
/// the transfer of cloudborne aerosols to interstitial aerosols in variable-
/// mode systems. See prog_process_stub.f90 for a complete description.
class ProgFProcessStub: public FPrognosticProcess {
  public:

  /// Create a prognostic process with the specified exponential decay rate.
  explicit ProgFProcessStub(Real decay_rate);
};

} // end haero namespace

#endif
