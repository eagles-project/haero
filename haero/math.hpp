#ifndef HAERO_MATH_HPP
#define HAERO_MATH_HPP

#include <haero/haero.hpp>

#include <ekat/ekat_assert.hpp>
#include <ekat/ekat_scalar_traits.hpp>

#ifndef __CUDACC__
#include <cmath>
#endif

namespace haero {

// bring in std:: or CUDA C++ functions
#ifdef __CUDACC__
using ::cbrt;
using ::erf;
using ::exp;
using ::expm1;
using ::isinf;
using ::isnan;
using ::log;
using ::log10;
using ::pow;
using ::sqrt;
using ::tanh;
using ::tgamma;
#else
using std::cbrt;
using std::erf;
using std::exp;
using std::expm1;
using std::isinf;
using std::isnan;
using std::log;
using std::log10;
using std::pow;
using std::sqrt;
using std::tanh;
using std::tgamma;
#endif

//// bring in kokkos math functions
// using namespace Kokkos::Experimental;

/// Returns the minimum of a and b.
KOKKOS_INLINE_FUNCTION
Real min(Real a, Real b) { return (a < b) ? a : b; }

// Returns the maximum of a and b.
KOKKOS_INLINE_FUNCTION
Real max(Real a, Real b) { return (a > b) ? a : b; }

KOKKOS_INLINE_FUNCTION
Real square(Real x) { return x * x; }

KOKKOS_INLINE_FUNCTION
Real cube(Real x) { return x * x * x; }

} // namespace haero
#endif
