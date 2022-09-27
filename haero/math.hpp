#ifndef HAERO_MATH_HPP
#define HAERO_MATH_HPP

#include <haero/check.hpp>
#include <haero/constants.hpp>
#include <haero/floating_point.hpp>
#include <haero/haero.hpp>

#include <ekat/ekat_assert.hpp>
#include <ekat/ekat_pack.hpp>
#include <ekat/ekat_pack_math.hpp>
#include <ekat/ekat_scalar_traits.hpp>

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
#ifdef __CUDACC__
using ::exp;
using ::log;
using ::sqrt;
using ::isnan;
#else
using std::exp;
using std::log;
using std::sqrt;
using std::isnan;
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

template <typename T>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!ekat::ScalarTraits<T>::is_simd, T>::type
    square(T x) {
  static_assert(std::is_arithmetic<T>::value, "square: bad type");
  return x * x;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<!ekat::ScalarTraits<T>::is_simd, T>::type
    cube(T x) {
  static_assert(std::is_arithmetic<T>::value, "arithmetic type");
  return x * x * x;
}

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
template <typename T>
struct LegendreCubic {
  static_assert(std::is_floating_point<
                    typename ekat::ScalarTraits<T>::scalar_type>::value,
                "floating point type");

  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  KOKKOS_INLINE_FUNCTION
  LegendreCubic() {}

  KOKKOS_INLINE_FUNCTION
  T operator()(const T& xin) const { return 0.5 * (-3 * xin + 5 * cube(xin)); }

  /// @brief Derivative of the cubic Legendre polynomial, @f$ \frac{d
  /// P_3}{dx}(x) @f$
  KOKKOS_INLINE_FUNCTION
  T derivative(const T& xin) const { return 0.5 * (-3 + 15 * square(xin)); }
};

/** @brief Quartic Legendre polynomial, @f$ P_4(x) @f$

  This struct is used for unit tests and to demonstrate the required
  interface for scalar functions to be used with rootfinding algorithms.

  Each scalar function must implement the T operator() (const T x) method.

  If it is to be used by the Newton solver, it must also implement the
  derivative(const T x) method, which returns f'(x).

  For an application example, see KohlerPolynomial.
*/
template <typename T>
struct LegendreQuartic {
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
template <typename T>
struct MonicParabola {
  static_assert(std::is_floating_point<typename ekat::ScalarTraits<T>::scalar_type>::value,
    "floating point type.");
  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  static constexpr Real a = 0;
  static constexpr Real b = 3;

  KOKKOS_INLINE_FUNCTION
  MonicParabola() {}

  KOKKOS_INLINE_FUNCTION
  T operator() (const T xin) const {
    return square(xin) + a * xin + b;
  }

  KOKKOS_INLINE_FUNCTION
  T derivative(const T xin) const {
    return 2*xin + a;
  }
};

/** @brief Scalar rootfinding algorithm that employs Newton's method.
  This algorithm has a quadratic convergence rate but it is not guaranteed to
  converge. It may converge to an incorrect root if a poor initial guess is
  chosen. It requires both function values and derivative values.
*/
// namespace {
// template <typename T>
// KOKKOS_INLINE_FUNCTION Real scalarize(const T pack) {
//   return pack[0];
// }
// template <>
// KOKKOS_INLINE_FUNCTION Real scalarize(const Real r) {
//   return r;
// }
// }  // namespace
template <typename ScalarFunction>
struct NewtonSolver {
  using value_type = typename ScalarFunction::value_type;

  /// maximum number of iterations allowed per solve
  static constexpr int max_iter = 200;
  /// solution
  value_type xroot;
  /// tolerance
  Real conv_tol;
  /// iteration counter
  int counter;
  /// absolute difference between two consecutive iterations
  value_type iter_diff;
  /// Scalar function whose root we need
  const ScalarFunction& f;
  /// true if a failure condition is met
  bool fail;

  /** @brief Constructor.

    @param [in] x0 Initial guess for the rootfinding algorithm
    @param [in] a0 Unused. Required to have unified interface for all solvers.
    @param [in] b0 Unused. Required to have unified interface for all solvers.
    @param [in] tol Convergence tolerance that determines when a root is found.
    @param [in] fn ScalarFunction instance whose root needs to be found
  */
  KOKKOS_INLINE_FUNCTION
  NewtonSolver(const value_type x0, const value_type a0, const value_type b0, const Real& tol,
                     const ScalarFunction& fn)
      : xroot(x0),
        conv_tol(tol),
        counter(0),
        iter_diff(std::numeric_limits<Real>::max()),
        f(fn),
        fail(false) {}

  /// Solves for the root.  Prints a warning message if the convergence
  /// tolerance is not met before the maximum number of iterations is achieved.
  KOKKOS_INLINE_FUNCTION
  value_type solve() {return solve_impl<value_type>();}

  template <typename VT>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<std::is_floating_point<VT>::value,
                              value_type>::type
      solve_impl() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      const value_type xnp1 = xroot - f(xroot) / f.derivative(xroot);
      iter_diff = abs(xnp1 - xroot);
      keep_going = !(FloatingPoint<value_type>::zero(iter_diff, conv_tol));
      EKAT_KERNEL_ASSERT_MSG(counter <= max_iter, "NewtonSolver: max iterations");
      if (counter > max_iter) {
        keep_going = false;
        fail = true;
      }
      if (isnan(xnp1)) {
        keep_going = false;
        fail = true;
      }
      xroot = xnp1;
    }
    return xroot;
  }

  template <typename VT>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<VT::packtag, value_type>::type
  solve_impl() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      const value_type xnp1 = xroot - f(xroot) / f.derivative(xroot);
      iter_diff = abs(xnp1 - xroot);
      keep_going = !(FloatingPoint<value_type>::zero(iter_diff, conv_tol));
      EKAT_KERNEL_ASSERT_MSG(counter <= max_iter, "NewtonSolver: max iterations");
      if (counter > max_iter) {
        keep_going = false;
        fail = true;
      }
      if (isnan(xnp1).any()) {
        keep_going = false;
        fail = true;
      }
      xroot = xnp1;
    }
    return xroot;
  }
};

/** padding factor stops BracketedNewtonSolver from getting stuck at an endpoint

  Larger values tend to require more iterations for the solver to converge.
*/
static constexpr Real bracket_pad_factor = 1e-7;

/** @brief Newton iterations protected by a bracket that prevents Newton from
going outside its bounds.

  This solver is a compromise between the speed of Newton's method and the
robustness of the bisection method. In addition to the requirement that the
initial interval contains a root, this solver requires that function value at
the initial interval endpoints have opposite sign.
*/
template <typename ScalarFunction>
struct BracketedNewtonSolver {
  using value_type = typename ScalarFunction::value_type;
  static constexpr int max_iter = 200;
  value_type xroot;
  value_type a;
  value_type b;
  Real conv_tol;
  value_type fa;
  value_type fx;
  value_type fb;
  const ScalarFunction& f;
  int counter;
  value_type iter_diff;
  /// true if a failure condition is met
  bool fail;

  /** @brief Constructor.

    @param [in] x0 x0 Initial guess for the rootfinding algorithm
    @param [in] a0 Left endpoint of the initial interval
    @param [in] b0 Right endpoint of the initial interval
    @param [in] tol Convergence tolerance that determines when a root is found.
    @param [in] fn ScalarFunction instance whose root needs to be found
  */
  KOKKOS_INLINE_FUNCTION
  BracketedNewtonSolver(const value_type x0, const value_type a0, const value_type b0,
                        const Real tol, const ScalarFunction& fn)
      : xroot(x0),
        a(a0),
        b(b0),
        conv_tol(tol),
        fa(fn(value_type(a0))),
        fx(fn(x0)),
        fb(fn(value_type(b0))),
        f(fn),
        counter(0),
        iter_diff(std::numeric_limits<Real>::max()),
        fail(false) {
    EKAT_KERNEL_ASSERT(Check<value_type>::is_positive(b - a));
    EKAT_KERNEL_ASSERT(Check<value_type>::is_negative(fa * fb));
  }

  KOKKOS_INLINE_FUNCTION
  value_type solve() { return solve_impl<value_type>(); }

  template <typename VT>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<std::is_floating_point<VT>::value,
                              value_type>::type
      solve_impl() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      // newton step
      value_type x = xroot - fx / f.derivative(xroot);
      // safeguard: require x to be inside current bracket
      // assure progress: guard against tiny steps
      const Real pad_fac = bracket_pad_factor;
      const Real pad = pad_fac * (b - a);
      if (!FloatingPoint<Real>::in_bounds(x, a + pad, b - pad)) {
        x = 0.5 * (a + b);
      }
      fx = f(x);
      // update bracket
      if (fx * fa > 0) {
        a = x;
        fa = fx;
      } else {
        b = x;
        fb = fx;
      }
      // check convergence
      iter_diff = abs(x - xroot);
      keep_going = !FloatingPoint<value_type>::zero(iter_diff, conv_tol);
      // prevent infinite loops
      EKAT_KERNEL_ASSERT_MSG(counter <= max_iter, "BracketedNewtonSolver: max iterations");
      if (counter > max_iter) {
        keep_going = false;
        fail = true;
      }
      if (isnan(x)) {
        keep_going = false;
        fail = true;
      }
      xroot = x;
    }
    return xroot;
  }

  template <typename VT>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<VT::packtag, value_type>::type
  solve_impl() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      // newton step
      VT x = xroot - fx / f.derivative(xroot);
      // safeguard: require x to be inside current bracket
      // assure progress: guard against tiny steps
      const Real pad_fac = bracket_pad_factor;
      const VT pad = pad_fac * (b - a);
      const auto mbnd = (x >= (a + pad)) and (x <= (b - pad));
      const auto mid = 0.5*(a + b);
      x.set(mbnd, x, mid);
      fx = f(x);
      // update bracket
      const auto msign = (fx * fa > 0);
      a.set(msign, x, a);
      fa.set(msign, fx, fa);
      b.set(msign, b, x);
      fb.set(msign, fb, fx);
      // check convergence
      iter_diff = abs(x - xroot);
      keep_going = !FloatingPoint<VT>::zero(iter_diff, conv_tol);
      // prevent infinite loops
      EKAT_KERNEL_ASSERT_MSG(counter <= max_iter, "NewtonSolver: max iterations");
      if (counter > max_iter) {
        keep_going = false;
        fail = true;
      }
      if (isnan(x).any()) {
        keep_going = false;
        fail = true;
      }
      xroot = x;
    }
    return xroot;
  }
};

/** @brief Scalar rootfinding algorithm that employs the recursive bisection
  method. This method has only a linear convergence rate, but it is guaranteed
  to converge if the initial interval contains a root.

  The Legendre polynomials above demonstrate the required ScalarFunction
  interface.

  This solver requires the initial interval to contain a root.

  For an application example, see KohlerPolynomial.

*/
template <typename ScalarFunction>
struct BisectionSolver {
  using value_type = typename ScalarFunction::value_type;

  /// maximum number of iterations allowed
  static constexpr int max_iter = 200;
  /// solution
  value_type xroot;
  /// left endpoint of root search interval
  value_type a;
  /// right endpoint of root search interval
  value_type b;
  /// tolerance
  Real conv_tol;
  /// function value at left endpoint of search interval
  value_type fa;
  /// next iteration solution
  value_type xnp1;
  /// iteration counter
  int counter;
  /// width of the current search interval
  value_type iter_diff;
  /// scalar function whose root we need
  const ScalarFunction& f;
  /// true if a failure condition is met
  bool fail;

  /** @brief Constructor.

    @param [in] x0 Unused. Required to have unified interface for all solvers.
    @param [in] a0 Left endpoint of the initial interval
    @param [in] b0 Right endpoint of the initial interval
    @param [in] tol Convergence tolerance that determines when a root is found.
    @param [in] fn ScalarFunction instance whose root needs to be found
  */
  KOKKOS_INLINE_FUNCTION
  BisectionSolver(const value_type x0, const value_type a0, const value_type b0, const Real tol,
                  const ScalarFunction& fn)
      : xroot(0.5 * (a0 + b0)),
        a(a0),
        b(b0),
        conv_tol(tol),
        fa(fn(a0)),
        xnp1(0.5 * (a0 + b0)),
        counter(0),
        iter_diff(b0 - a0),
        f(fn),
        fail(false) {}

  /// Solves for the root.  Prints a warning message if the convergence
  /// tolerance is not met before the maximum number of iterations is achieved.
  KOKKOS_INLINE_FUNCTION
  value_type solve() { return solve_impl<value_type>(); }

  template <typename VT>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<std::is_floating_point<VT>::value,
                              value_type>::type
      solve_impl() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      const value_type fx = f(xroot);
      xnp1 = 0.5 * (a + b);
      if (fx * fa < 0) {
        b = xroot;
      } else {
        a = xroot;
        fa = fx;
      }
      iter_diff = b - a;
      xroot = xnp1;
      keep_going = !(FloatingPoint<value_type>::zero(iter_diff, conv_tol));
      EKAT_KERNEL_ASSERT_MSG(counter <= max_iter, "BisectionSolver: max iterations");
      if (counter > max_iter) {
        keep_going = false;
        fail = true;
      }
      if (isnan(xnp1)) {
        keep_going = false;
        fail = true;
      }
    }
    return xroot;
  }

  template <typename VT>
  KOKKOS_INLINE_FUNCTION typename std::enable_if<VT::packtag, value_type>::type
  solve_impl() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      xnp1 = 0.5 * (a + b);
      const value_type fx = f(xroot);
      const auto m = (fx * fa < 0);
      b.set(m, xroot);
      a.set(!m, xroot);
      fa.set(!m, fx);
      iter_diff = b - a;
      xroot = xnp1;
      keep_going = !(FloatingPoint<value_type>::zero(iter_diff, conv_tol));
      EKAT_KERNEL_ASSERT_MSG(counter <= max_iter, "BisectionSolver: max iterations");
      if (counter > max_iter) {
        keep_going = false;
        fail = true;
      }
      if (isnan(xnp1).any()) {
        keep_going = false;
        fail = true;
      }
    }
    return xroot;
  }
};

}  // namespace math

}  // namespace haero
#endif
