#ifndef HAERO_MAM4_NUCLEATION_FPROCESS_HPP
#define HAERO_MAM4_NUCLEATION_FPROCESS_HPP

#include "haero/process.hpp"

namespace haero {

/// This prognostic process implements a simple exponential decay model for
/// the transfer of cloudborne aerosols to interstitial aerosols in variable-
/// mode systems. See prog_process_stub.f90 for a complete description.
class Mam4NucleationFProcess: public FPrognosticProcess {
  public:

  /// Create a MAM4 nucleation process implemented in Fortran. This
  /// implementation is based on that found in modal_aero_newnuc.F90 in
  /// MAM4.
  Mam4NucleationFProcess();
};

} // end haero namespace

#endif
