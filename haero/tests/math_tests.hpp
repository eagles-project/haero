// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HAERO_MATH_TESTS_HPP
#define HAERO_MATH_TESTS_HPP

#include "haero/floating_point.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"

#include "haero/math.hpp"
#include "haero/root_finders.hpp"

namespace haero {
namespace math {

/** @brief Cubic Legendre polynomial, @f$ P_3(x) @f$

  This struct is used for unit tests and to demonstrate the required
  interface for scalar functions to be used with rootfinding algorithms.

  Each scalar function must implement the T operator() (const T x) method that
  returns the function's value at x.

  If it is to be used by the Newton solver, it must also implement the
  derivative(const T x) method, which returns f'(x).

  For an application example, see KohlerPolynomial.
*/
template <typename T> struct LegendreCubic {
  static_assert(std::is_floating_point<
                    typename ekat::ScalarTraits<T>::scalar_type>::value,
                "floating point type");

  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  KOKKOS_INLINE_FUNCTION
  LegendreCubic() {}

  KOKKOS_INLINE_FUNCTION
  T operator()(const T &xin) const { return 0.5 * (-3 * xin + 5 * cube(xin)); }

  /// @brief Derivative of the cubic Legendre polynomial, @f$ \frac{d
  /// P_3}{dx}(x) @f$
  KOKKOS_INLINE_FUNCTION
  T derivative(const T &xin) const { return 0.5 * (-3 + 15 * square(xin)); }
};

/** @brief Quartic Legendre polynomial, @f$ P_4(x) @f$

  This struct is used for unit tests and to demonstrate the required
  interface for scalar functions to be used with rootfinding algorithms.

  Each scalar function must implement the T operator() (const T x) method.

  If it is to be used by the Newton solver, it must also implement the
  derivative(const T x) method, which returns f'(x).

  For an application example, see KohlerPolynomial.
*/
template <typename T> struct LegendreQuartic {
  static_assert(std::is_floating_point<
                    typename ekat::ScalarTraits<T>::scalar_type>::value,
                "floating point type.");

  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  KOKKOS_INLINE_FUNCTION
  LegendreQuartic() {}

  /// @brief Quartic Legendre polynomial, @f$ P_4(x) @f$
  KOKKOS_INLINE_FUNCTION
  T operator()(const T xin) const {
    return 0.125 * (3 - 30 * square(xin) + 35 * pow(xin, 4));
  }

  /// @brief Quartic Legendre polynomial, @f$ \frac{d P_4}{dx}(x) @f$
  KOKKOS_INLINE_FUNCTION
  T derivative(const T xin) const {
    return 0.125 * (-60 * xin + 140 * cube(xin));
  }
};

/** @brief Quadratic polynomial.

  This struct is used for unit tests and to demonstrate the required
  interface for scalar functions to be used with rootfinding algorithms.

  Each scalar function must implement the T operator() (const T x) method.

  If it is to be used by the Newton solver, it must also implement the
  derivative(const T x) method, which returns f'(x).
*/
template <typename T> struct MonicParabola {
  static_assert(std::is_floating_point<
                    typename ekat::ScalarTraits<T>::scalar_type>::value,
                "floating point type.");
  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  Real a;
  Real b;

  KOKKOS_INLINE_FUNCTION
  MonicParabola() : a(0), b(3) {}

  KOKKOS_INLINE_FUNCTION
  T operator()(const T xin) const {
    const T aa(a);
    const T bb(b);
    return square(xin) + aa * xin + bb;
  }

  KOKKOS_INLINE_FUNCTION
  T derivative(const T xin) const {
    const T aa(a);
    return 2 * xin + aa;
  }
};

} // namespace math
} // namespace haero
#endif
