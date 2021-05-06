#ifndef HAERO_MATH_HELPERS_HPP
#define HAERO_MATH_HELPERS_HPP

#include "haero.hpp"
#include "floating_point.hpp"
#include "physical_constants.hpp"
#include <cmath>
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_scalar_traits.hpp"

namespace haero {

template <typename T> KOKKOS_INLINE_FUNCTION
typename std::enable_if<!ekat::ScalarTraits<T>::is_simd, T>::type square(T x) {
  static_assert(std::is_arithmetic<T>::value, "square: bad type");
  return x*x;
}

template <typename T> KOKKOS_INLINE_FUNCTION
typename std::enable_if<!ekat::ScalarTraits<T>::is_simd,T>::type cube(T x) {
  static_assert(std::is_arithmetic<T>::value, "arithmetic type");
  return x*x*x;
}

template <typename T> KOKKOS_INLINE_FUNCTION
T sphere_radius_from_volume(const T vol) {
  static_assert(std::is_floating_point<typename ekat::ScalarTraits<T>::scalar_type>::value, "arithmetic type");
  static constexpr Real four_thirds_pi = 4*constants::pi/3.0;
  return cbrt(vol / four_thirds_pi);
}

namespace math {

/** @brief Performs a single iteration of Newton's method for finding the root of a scalar
  equation, e.g.,  @f$ f(x) = 0 @f$.

  @param [in] xn Current solution, @f$ x_n @f$
  @param [in] fx Function value at current solution, @f$ f(x_n) @f$
  @param [in] fpx Function derivative at current solution, @f$ f'(x_n) @f$
  @return Next iterate,  @f$ x_{n+1} @f$
*/
template <typename T> KOKKOS_INLINE_FUNCTION
T next_newton_scalar_iteration(const T& xn, const T& fx, const T& fpx) {
  static_assert(std::is_floating_point<typename ekat::ScalarTraits<T>::scalar_type>::value, "floating point type.");
  return xn - fx / fpx;
}

/** @brief Performs a single iteration of the recursive bisection algorithm for finding
  roots of scalar equations @f$ f(x) = 0 @f$.

  @param [out] xnp1 Next iterate, @f$ x_{n+1} @f$
  @param [inout] an Left endpoint of current search interval, @f$ a_n @f$
  @param [inout] bn Right endpoint of current search interval, @f$ b_n @f$
  @param [inout] fan Function value at left endpoint of current search interval, @f$ f(a_n) @f$
  @param [in] xn Current solution, @f$ x_n @f$
  @param [in] fx Function value at current solution, @f$ f(x_n) @f$

*/
template <typename T> KOKKOS_INLINE_FUNCTION
void next_bisection_scalar_iteration(T& xnp1, T& an, T& bn, T& fan, const T& xn, const T& fx)
{
  static_assert(std::is_floating_point<typename ekat::ScalarTraits<T>::scalar_type>::value, "floating point type.");
  xnp1 = 0.5*(an + bn);
  if (fx * fan < 0) {
    bn = xn;
  }
  else {
    an = xn;
    fan = fx;
  }
}

template <> KOKKOS_INLINE_FUNCTION
void next_bisection_scalar_iteration<PackType>(PackType& xnp1, PackType& an, PackType& bn,
  PackType& fan, const PackType& xn, const PackType& fx) {

  xnp1 = 0.5*(an + bn);
  const auto m = (fx * fan < 0);
  const auto not_m = !m;
  ekat_masked_loop(m, s) {bn[s] = xn[s];};
  ekat_masked_loop(not_m, s) {an[s] = xn[s]; fan[s] = fx[s];};
}



/** @brief Cubic Legendre polynomial, @f$ P_3(x) @f$

  This struct is used for unit tests and to demonstrate the required
  interface for scalar functions to be used with rootfinding algorithms.

  Each scalar function must implement the operator() (const Real x) method that returns the
  function's value at x.

  If it is to be used by the Newton solver, it must also implement the
  derivative(const Real x) method, which returns f'(x).

  For an application example, see KohlerPolynomial.
*/
template <typename T>
struct LegendreCubic {
  static_assert(std::is_floating_point<typename ekat::ScalarTraits<T>::scalar_type>::value, "floating point type");

  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  KOKKOS_INLINE_FUNCTION
  LegendreCubic() {}

  KOKKOS_INLINE_FUNCTION
  T operator() (const T& xin) const {
    return 0.5*(-3*xin + 5*cube(xin));
  }

  /// @brief Derivative of the cubic Legendre polynomial, @f$ \frac{d P_3}{dx}(x) @f$
  KOKKOS_INLINE_FUNCTION
  T derivative(const T& xin) const {
    return 0.5*(-3 + 15*square(xin));
  }
};


/** @brief Quartic Legendre polynomial, @f$ P_4(x) @f$

  This struct is used for unit tests and to demonstrate the required
  interface for scalar functions to be used with rootfinding algorithms.

  Each scalar function must implement the operator() (const Real x) method.

  If it is to be used by the Newton solver, it must also implement the
  derivative(const Real x) method, which returns f'(x).

  For an application example, see KohlerPolynomial.
*/
template <typename T>
struct LegendreQuartic {
  static_assert(std::is_floating_point<typename ekat::ScalarTraits<T>::scalar_type>::value, "floating point type.");

  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  KOKKOS_INLINE_FUNCTION
  LegendreQuartic() {}

  /// @brief Quartic Legendre polynomial, @f$ P_4(x) @f$
  KOKKOS_INLINE_FUNCTION
  T operator() (const T xin) const {
    return 0.125 * (3 - 30*square(xin) + 35*pow(xin,4));
  }

  /// @brief Quartic Legendre polynomial, @f$ \frac{d P_4}{dx}(x) @f$
  KOKKOS_INLINE_FUNCTION
  T derivative(const T xin) const {
    return 0.125 * (-60*xin + 140*cube(xin));
  }
};


/** @brief Scalar rootfinding algorithm that employs Newton's method.
  This algorithm has a quadratic convergence rate but it is not guaranteed to converge.
  It may converge to an incorrect root if a poor initial guess is chosen.
  It requires both function values and derivative values.

  Template parameter ScalarFunction must implement the function_eval(const Real x)
  and derivative_eval(const Real x) methods.

  The Legendre polynomials above demonstrate the required interface.

  For an application example, see KohlerPolynomial.

*/
template <typename ScalarFunction>
struct ScalarNewtonSolver {

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

  /** @brief Constructor.

    @param [in] xinit Initial guess for the rootfinding algorithm
    @param [in] tol Convergence tolerance that determines when a root is found.
    @param [in] fn ScalarFunction instance whose root needs to be found
  */
  KOKKOS_INLINE_FUNCTION
  ScalarNewtonSolver(const value_type& xinit, const Real& tol, const ScalarFunction& fn) :
    xroot(xinit),
    conv_tol(tol),
    counter(0),
    iter_diff(std::numeric_limits<Real>::max()),
    f(fn)
    {}


  /// Solves for the root.  Prints a warning message if the convergence tolerance is not met before the maximum number of iterations is achieved.
  KOKKOS_INLINE_FUNCTION
  value_type solve() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      const value_type fxn = f(xroot);
      const value_type fprimexn = f.derivative(xroot);
      const value_type xnp1 = next_newton_scalar_iteration(xroot, fxn, fprimexn);
      iter_diff = abs(xnp1 - xroot);
      xroot = xnp1;
      keep_going = !(FloatingPoint<value_type>::zero(iter_diff, conv_tol));
      if (counter >= max_iter) {
        printf("newton solve warning: max iterations reached\n");
        keep_going = false;
      }
    }
    return xroot;
  }
};


/** @brief Scalar rootfinding algorithm that employs the recursive bisection method.
  This method has only a linear convergence rate, but it is guaranteed to converge if
  the initial interval contains a root.

  Template parameter ScalarFunction must implement the function_eval(const Real x) method.
  The derivative_eval(const Real x) method is not necessary to use this solver.

  The Legendre polynomials above demonstrate the required interface.

  For an application example, see KohlerPolynomial.

*/
template <typename ScalarFunction> struct BisectionSolver {
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

  /** @brief Constructor.

    @param [in] a0 Left endpoint of the initial interval
    @param [in] b0 Right endpoint of the initial interval
    @param [in] tol Convergence tolerance that determines when a root is found.
    @param [in] fn ScalarFunction instance whose root needs to be found
  */
  KOKKOS_INLINE_FUNCTION
  BisectionSolver(const value_type& a0, const value_type& b0, const Real tol, const ScalarFunction& fn) :
    xroot(0.5*(a0 + b0)),
    a(a0),
    b(b0),
    conv_tol(tol),
    fa(fn(a0)),
    xnp1(0.5*(a0+b0)),
    counter(0),
    iter_diff(b0-a0),
    f(fn) {}

  /// Solves for the root.  Prints a warning message if the convergence tolerance is not met before the maximum number of iterations is achieved.
  KOKKOS_INLINE_FUNCTION
  value_type solve() {
    bool keep_going = true;
    while (keep_going) {
      ++counter;
      const value_type fx = f(xroot);
      next_bisection_scalar_iteration(xnp1, a, b, fa, xroot, fx);
      iter_diff = b-a;
      xroot = xnp1;
      keep_going = !(FloatingPoint<value_type>::zero(iter_diff, conv_tol));
      if (counter >= max_iter) {
        printf("bisection solve warning: max iterations reached");
        keep_going = false;
      }
    }
    return xroot;
  }
};


} // namespace math


} // namespace haero
#endif
