#ifndef HAERO_FAEROSOL_PROCESS_STUB_HPP
#define HAERO_FAEROSOL_PROCESS_STUB_HPP

#include "haero/aerosol_process.hpp"

namespace haero {

/// This prognostic process implements a simple exponential decay model for
/// the transfer of cloudborne aerosols to interstitial aerosols in variable-
/// mode systems. See prog_process_stub.f90 for a complete description.
class FAerosolProcessStub: public FAerosolProcess {
  public:

  /// Create a prognostic process with the specified exponential decay rate.
  explicit FAerosolProcessStub(Real decay_rate);
};

} // end haero namespace

#endif
