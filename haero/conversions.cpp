#include "haero/conversions.hpp"

using haero::Real;

// Interoperable C functions for converting data in Fortran.
// See haero_conversions.F90 for details on how these functions are used.
extern "C" {

Real relative_humidity_from_vapor_mixing_ratio_f(Real qv, Real p, Real T) {
  return haero::conversions::relative_humidity_from_vapor_mixing_ratio(qv, p, T);
}

}  // extern "C"

