#ifndef HAERO_FLOATING_POINT_UTILS_HPP
#define HAERO_FLOATING_POINT_UTILS_HPP

#include "haero/haero_config.hpp"
#include "ekat/ekat_assert.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

/**
  struct for help with common floating point operations

*/
template <typename T=Real>
struct FloatingPoint {
  static_assert(std::is_floating_point<T>::value, "floating point type required.");

  /// Default tolerance for floating point comparisons
  static constexpr T zero_tol = (std::is_same<T,float>::value) ? 1.0E-7 : 1.0E-13;

  /// Define floating point zero by @f$\lvert x \rvert < \epsilon_{tol}@f$
  KOKKOS_INLINE_FUNCTION
  static bool zero(const T x, const T tol=zero_tol) {
    EKAT_ASSERT(tol>0);
    return std::abs(x) < tol;
  }

  /// Define floating point equivalence by @f$\lvert x_0 - x_1 \rvert < \epsilon_{tol}@f$
  KOKKOS_INLINE_FUNCTION
  static bool equiv(const T x0, const T x1, const T tol=zero_tol) {
    EKAT_ASSERT(tol>0);
    return std::abs(x0-x1) < tol;
  }

  /** Define floating point in bounds as @f$ l - \epsilon_{tol} < x < u + \epsilon_{tol}@f$

    @param [in] x
    @param [in] lower lower bound @f$l@f$
    @param [in] upper upper bound @f$u@f$
    @param [in] tol tolerance @f$\epsilon_{tol}@f$.
  */
  KOKKOS_INLINE_FUNCTION
  static bool in_bounds(const T x, const T lower, const T upper, const T tol = zero_tol) {
    EKAT_ASSERT(tol>0);
    return (x >= (lower - tol) && x <= (upper + tol));
  }

  /** multiplier for safe division by x, @f$ \frac{1}{x} \approx \frac{x}{x^2 + \epsilon_{tol}^2}@f$

    For use with removable singularities.

    @warning the return value is a multiplication factor, not a divisor

    @param [in] x
    @param [in] tol
    @return @f$\frac{1}{x} \approx \frac{x}{x^2 + \epsilon_{tol}^2@f$
  */
  KOKKOS_INLINE_FUNCTION
  static T safe_denominator(const T x, const T tol=zero_tol) {
    EKAT_ASSERT(tol>0);
    return x / (x*x + tol*tol);
  }
};


} // namespace haero
#endif
