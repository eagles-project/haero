// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HAERO_MATH_HPP
#define HAERO_MATH_HPP

#include <haero/haero.hpp>

#include <ekat/ekat_assert.hpp>
#include <ekat/ekat_scalar_traits.hpp>

#include <Kokkos_NumericTraits.hpp>

#ifndef __CUDACC__
#include <cmath>
#include <limits>
#endif

namespace haero {

// bring in std:: or CUDA C++ functions
#ifdef __CUDACC__
using ::atan;
using ::cbrt;
using ::cos;
using ::erf;
using ::erfc;
using ::exp;
using ::expm1;
using ::isinf;
using ::isnan;
using ::log;
using ::log10;
using ::pow;
using ::round;
using ::sin;
using ::sqrt;
using ::tanh;
using ::tgamma;
KOKKOS_INLINE_FUNCTION constexpr Real min() {
  return Kokkos::Experimental::norm_min_v<Real>;
}
KOKKOS_INLINE_FUNCTION constexpr Real max() {
  return Kokkos::Experimental::finite_max_v<Real>;
}
KOKKOS_INLINE_FUNCTION constexpr Real epsilon() {
  return Kokkos::Experimental::epsilon_v<Real>;
}
#else
using std::atan;
using std::cbrt;
using std::cos;
using std::erf;
using std::erfc;
using std::exp;
using std::expm1;
using std::isinf;
using std::isnan;
using std::log;
using std::log10;
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::tanh;
using std::tgamma;
KOKKOS_INLINE_FUNCTION constexpr Real min() {
  return std::numeric_limits<Real>::min();
}
KOKKOS_INLINE_FUNCTION constexpr Real max() {
  return std::numeric_limits<Real>::max();
}
KOKKOS_INLINE_FUNCTION constexpr Real epsilon() {
  return std::numeric_limits<Real>::epsilon();
}
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
