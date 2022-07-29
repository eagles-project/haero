#ifndef HAERO_MATH_HPP
#define HAERO_MATH_HPP

#include <haero/haero.hpp>

#include <ekat/ekat_pack_math.hpp>

#ifndef __CUDACC__
#include <cmath>
#endif

namespace haero {

// bring in ekat math functions
using ekat::abs;
using ekat::exp;
using ekat::expm1;
using ekat::log;
using ekat::log10;
using ekat::tgamma;
using ekat::sqrt;
using ekat::cbrt;
using ekat::tanh;
using ekat::erf;

// bring in std:: or CUDA C++ functions
#ifndef __CUDACC__
using ::exp;
using ::log;
using ::sqrt;
#else
using std::exp;
using std::log;
using std::sqrt;
#endif

//// bring in kokkos math functions
//using namespace Kokkos::Experimental;

// min functions for Haero's Real type and for mixed Real/Pack arguments
KOKKOS_INLINE_FUNCTION
Real min(Real a, Real b) {
  return (a < b) ? a : b;
}

KOKKOS_INLINE_FUNCTION
PackType min(Real a, const PackType& b) {
  PackType result = b;
  result.set(a < b, a);
  return result;
}

KOKKOS_INLINE_FUNCTION
PackType min(const PackType& a, Real b) {
  PackType result = a;
  result.set(b < a, b);
  return result;
}

// max functions for Haero's Real type and for mixed Real/Pack arguments
KOKKOS_INLINE_FUNCTION
Real max(Real a, Real b) {
  return (a > b) ? a : b;
}

KOKKOS_INLINE_FUNCTION
PackType max(Real a, const PackType& b) {
  PackType result = b;
  result.set(a > b, a);
  return result;
}

KOKKOS_INLINE_FUNCTION
PackType max(const PackType& a, Real b) {
  PackType result = a;
  result.set(b > a, b);
  return result;
}

}  // end namespace haero

#endif
