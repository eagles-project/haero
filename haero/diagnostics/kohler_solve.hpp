#ifndef HAERO_KOHLER_SOLVE_HPP
#define HAERO_KOHLER_SOLVE_HPP

#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "haero/floating_point.hpp"
#include "haero/haero.hpp"
#include "haero/math.hpp"

namespace haero {
/** Coefficient accounting for Kelvin effect on droplets

   In the documentation, this constant is denoted by A.

  @todo verify that the assumptions that lead to this constant still hold at
 SCREAM resolutions (e.g., constant temperature, constant surface tension) */
static constexpr Real kelvin_droplet_effect_coeff = 0.00120746723156361711;

/** @brief Struct that represents the Kohler polynomial.

  @f$ K(r_w) = \log(s) r_w^4 - A r_w^3 + (B - \log(s))r_d^3 r_w + A r_d^3 @f$

  where r_w is the wet radius, s is relative humidity, A is the Kelvin effect
  coefficient, B is hygroscopicity, and r_d is the dry radius.

  Each instance corresponds to a separate set of coefficients, which are
  functions of the inputs.

  This struct conforms to the interface prescribed in math_helpers.hpp for
  scalar functions that are to be used with numerical rootfinding algorithms.

  The Kohler polynomial is a quartic polynomial whose variable is particle wet
  radius. Equilibrium solutions are roots of this polynomial. Algebraically,
  there are 2 complex roots and 2 real roots.  Of the real roots, one is
  positive, the other is negative.
  Physically, only the real, positive root makes sense.

  This struct is templated on scalar type so that it can be used with PackType.
  If it is used with PackType, each element of the PackType corresponds to a
  separate KohlerPolynomial, with distinct coefficients.

  @warning This polynomial is severely ill-conditioned, to the point that it is
  sensitive to order-of-operations changes caused by compiler optimization
  flags.  We therefore require double precision.

  Properties of the Kohler Polynomial that are useful to finding its roots:

  1. K(0) = kelvin_droplet_effect_coeff * cube(r_dry) > 0
  2. K(r_dry) = r_dry**4 * hygroscopicity > 0

  Properties of the Kohler Polynomial that are useful to finding its roots given
  inputs that are within the bounds defined below:

  1. K(25*r_dry) < 0

*/
template <typename T = ekat::Pack<double, HAERO_PACK_SIZE>>
struct KohlerPolynomial {
  static_assert(
      std::is_same<typename ekat::ScalarTraits<T>::scalar_type, double>::value,
      "double precision required.");

  /// Minimum value of relative humidity
  static constexpr double rel_humidity_min = 0.05;
  /// Above this relative humidity is considered saturated air, and cloud
  /// aerosol processes would apply
  static constexpr double rel_humidity_max = 0.98;
  /// Minimum hygroscopicity for E3SM v1 aerosol species
  static constexpr double hygro_min = 1e-6;
  /// Maximum hygroscopicity for E3SM v1 aerosol species
  static constexpr double hygro_max = 1.3;
  /// Minimum particle size for E3SM v1
  static constexpr double dry_radius_min_microns = 1e-3;
  /// Maximum particle size for Kohler theory
  static constexpr double dry_radius_max_microns = 30;

  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

  /// Coefficient in the Kohler polynomial
  T log_rel_humidity;
  /// Coefficient in the Kohler polynomial
  T hygroscopicity;
  /// Safe return value
  T dry_radius;
  /// Coefficient in the Kohler polynomial
  T dry_radius_cubed;

  /** Constructor. Creates 1 instance of a KohlerPolynomial.

    @param rel_h relative humidity
    @param hygro hygroscopicity
    @param dry_rad_microns particle dry radius [ 1e-6 m ]
  */
  template <typename U>
  KOKKOS_INLINE_FUNCTION KohlerPolynomial(const U& rel_h, const U& hygro,
                                          const U& dry_rad_microns)
      : log_rel_humidity(log(rel_h)),
        hygroscopicity(hygro),
        dry_radius(dry_rad_microns),
        dry_radius_cubed(cube(dry_rad_microns)) {
    EKAT_KERNEL_ASSERT(valid_inputs(T(rel_h), T(hygro), T(dry_rad_microns)));
  }

  template <typename U>
  KOKKOS_INLINE_FUNCTION KohlerPolynomial(const MaskType& m, const U& rel_h,
                                          const U& hygro,
                                          const U& dry_rad_microns)
      : log_rel_humidity(log(rel_h)),
        hygroscopicity(hygro),
        dry_radius(dry_rad_microns),
        dry_radius_cubed(cube(dry_rad_microns)) {
    EKAT_KERNEL_ASSERT(valid_inputs(m, T(rel_h), T(hygro), T(dry_rad_microns)));
  }

  /** Evaluates the Kohler polynomial.

    @f$ K(r_w) = \log(s) r_w^4 - A r_w^3 + (B - \log(s))r_d^3 r_w + A r_d^3 @f$

    where r_w is the wet radius, s is relative humidity, A is the Kelvin effect
    coefficient, B is hygroscopicity, and r_d is the dry radius.

    @param [in] Polynomial input value, wet_radius wet radius in microns [ 1e-6
    m]
    @return Polynomial value, wet_radius in microns [ 1e-6 m]
  */
  template <typename U>
  KOKKOS_INLINE_FUNCTION T operator()(const U& wet_radius) const {
    const T rwet = T(wet_radius);
    const Real kelvinA = kelvin_droplet_effect_coeff;
    const T result = (log_rel_humidity * rwet - kelvinA) * cube(rwet) +
                     ((hygroscopicity - log_rel_humidity) * rwet + kelvinA) *
                         dry_radius_cubed;
    return result;
  }

  /** @brief Evaluates the derivative of the Kohler polynomial with respect to
    wet radius

    @f$ K'(r_w) = \frac{\partial K}(\partial r_w)(r_w) @f$

    @param [in] Polynomial input value, wet radius in microns [ 1e-6 m]
    @return Polynomial slope at input value
  */
  template <typename U>
  KOKKOS_INLINE_FUNCTION T derivative(const U& wet_radius) const {
    const T rwet = T(wet_radius);
    const T wet_radius_squared = square(rwet);
    const Real kelvinA = kelvin_droplet_effect_coeff;
    const T result =
        (4 * log_rel_humidity * rwet - 3 * kelvinA) * wet_radius_squared +
        (hygroscopicity - log_rel_humidity) * dry_radius_cubed;
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  bool valid_inputs(const T& relh, const T& hyg, const T& dry_rad) const {
    return (FloatingPoint<T>::in_bounds(relh, rel_humidity_min,
                                        rel_humidity_max) and
            FloatingPoint<T>::in_bounds(hyg, hygro_min, hygro_max) and
            FloatingPoint<T>::in_bounds(dry_rad, dry_radius_min_microns,
                                        dry_radius_max_microns));
  }

  KOKKOS_INLINE_FUNCTION
  bool valid_inputs() const {
    return valid_inputs(exp(this->log_rel_humidity), this->hygroscopicity,
                        this->dry_radius);
  }

  template <typename ST = T>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<ekat::ScalarTraits<ST>::is_simd, bool>::type
      valid_inputs(const MaskType& m) const {
    const T relative_humidity = exp(this->log_rel_humidity);
    const double rhmin = rel_humidity_min;
    const double rhmax = rel_humidity_max;
    const double hmin = hygro_min;
    const double hmax = hygro_max;
    const double rmin = dry_radius_min_microns;
    const double rmax = dry_radius_max_microns;
    const double tol = FloatingPoint<double>::zero_tol;
    const auto in_bounds_mask =
        (relative_humidity >= (rhmin-tol)) && (relative_humidity <= (rhmax+tol)) &&
        (this->hygroscopicity >= (hmin-tol)) && (this->hygroscopicity <= (hmax+tol)) &&
        (this->dry_radius >= (rmin-tol)) && (this->dry_radius <= (rmax+tol));
    return (in_bounds_mask || (!m)).all();
  }

  template <typename ST = T>
  KOKKOS_INLINE_FUNCTION
      typename std::enable_if<ekat::ScalarTraits<ST>::is_simd, bool>::type
      valid_inputs(const MaskType& m, const T& relh, const T& hyg,
                   const T& dry_rad) const {
    const double rhmin = rel_humidity_min;
    const double rhmax = rel_humidity_max;
    const double hmin = hygro_min;
    const double hmax = hygro_max;
    const double rmin = dry_radius_min_microns;
    const double rmax = dry_radius_max_microns;
    const double tol = FloatingPoint<double>::zero_tol;
    const auto in_bounds_mask = (relh >= (rhmin-tol)) && (relh <= (rhmax+tol)) &&
                                (hyg >= (hmin-tol)) && (hyg <= (hmax+tol)) &&
                                (dry_rad >= (rmin-tol)) && (dry_rad <= (rmax+tol));
    return (in_bounds_mask || (!m)).all();
  }

  /** @brief Writes a string containing a Mathematica script that may be used to
    generate the verification data.

    @param [in] n the number of points in each parameter range
  */
  static std::string mathematica_verification_program(const int n);

  /** @brief Writes a string containing a Matlab script that may be used to
    generate the verification data.

    @param [in] n the number of points in each parameter range
  */
  static std::string matlab_verification_program(const int n);
};

/** @brief This functor solves for the positive root of a packed instance of the
 * KohlerPolynomial.
 */
struct KohlerNewtonSolve {
  typedef ekat::Pack<double, HAERO_PACK_SIZE> double_pack;

  double_pack relative_humidity;
  double_pack hygroscopicity;
  double_pack dry_radius_microns;
  Real conv_tol;
  int n_iter;

  KOKKOS_INLINE_FUNCTION
  KohlerNewtonSolve(const PackType& rel_h, const PackType& hygro,
                    const PackType& dry_rad, const Real tol)
      : relative_humidity(rel_h),
        hygroscopicity(hygro),
        dry_radius_microns(dry_rad),
        conv_tol(tol),
        n_iter(0) {}

  KOKKOS_INLINE_FUNCTION
  PackType operator()() {
    double_pack wet_radius_init(25 * dry_radius_microns);
    double_pack result(-1);
    const auto kpoly = KohlerPolynomial<double_pack>(
        relative_humidity, hygroscopicity, dry_radius_microns);
    auto solver = math::ScalarNewtonSolver<KohlerPolynomial<double_pack>>(
        wet_radius_init, conv_tol, kpoly);
    result = solver.solve();
    n_iter = solver.counter;
    return PackType(result);
  }
};

/** @brief This functor solves for the positive root of a packed instance of the
 * KohlerPolynomial.
 */
struct KohlerBisectionSolve {
  typedef ekat::Pack<double, HAERO_PACK_SIZE> double_pack;
  double_pack relative_humidity;
  double_pack hygroscopicity;
  double_pack dry_radius_microns;
  Real conv_tol;
  int n_iter;

  KOKKOS_INLINE_FUNCTION
  KohlerBisectionSolve(const PackType& rel_h, const PackType& hygro,
                       const PackType& dry_rad, const Real tol)
      : relative_humidity(rel_h),
        hygroscopicity(hygro),
        dry_radius_microns(dry_rad),
        conv_tol(tol),
        n_iter(0) {}

  KOKKOS_INLINE_FUNCTION
  PackType operator()() {
    double_pack left_endpt(0.9 * dry_radius_microns);
    double_pack right_endpt(200.0);
    double_pack result(-1);
    const auto kpoly = KohlerPolynomial<double_pack>(
        relative_humidity, hygroscopicity, dry_radius_microns);
    auto solver = math::BisectionSolver<KohlerPolynomial<double_pack>>(
        left_endpt, right_endpt, conv_tol, kpoly);
    result = solver.solve();
    n_iter = solver.counter;
    return PackType(result);
  }
};

/** @brief This functor solves for the positive root of a packed instance of the
 * KohlerPolynomial.
 */
struct KohlerBracketedNewtonSolve {
  typedef ekat::Pack<double, HAERO_PACK_SIZE> double_pack;
  double_pack relative_humidity;
  double_pack hygroscopicity;
  double_pack dry_radius_microns;
  Real conv_tol;
  int n_iter;

  KOKKOS_INLINE_FUNCTION
  KohlerBracketedNewtonSolve(const PackType& rel_h, const PackType& hygro,
                             const PackType& dry_rad, const Real tol)
      : relative_humidity(rel_h),
        hygroscopicity(hygro),
        dry_radius_microns(dry_rad),
        conv_tol(tol),
        n_iter(0) {}

  KOKKOS_INLINE_FUNCTION
  PackType operator()() {
    // initial bracket
    const double a0 = ekat::min(dry_radius_microns);
    const double b0 = 25 * ekat::max(dry_radius_microns);
    double_pack r0(10 * dry_radius_microns);  // initial guess
    double_pack result(-1);
    const auto kpoly = KohlerPolynomial<double_pack>(
        relative_humidity, hygroscopicity, dry_radius_microns);
    auto solver = math::BracketedNewtonSolver<KohlerPolynomial<double_pack>>(
        r0, a0, b0, conv_tol, kpoly);
    result = solver.solve();
    EKAT_KERNEL_ASSERT((result > 0).all());
    n_iter = solver.counter;
    return PackType(result);
  }
};

}  // namespace haero
#endif
