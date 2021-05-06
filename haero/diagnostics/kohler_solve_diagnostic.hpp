#ifndef HAERO_KOHLER_SOLVE_HPP
#define HAERO_KOHLER_SOLVE_HPP

#include "haero/haero.hpp"
#include "haero/math_helpers.hpp"
#include "haero/floating_point.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_math.hpp"

namespace haero {

/** @brief Struct that represents the Kohler polynomial.

  @f$ K(r_w) = \log(s) r_w^4 - A r_w^3 + (B - \log(s))r_d^3 r_w + A r_d^3 @f$

  where r_w is the wet radius, s is relative humidity, A is the Kelvin effect coefficient,
    B is hygroscopicity, and r_d is the dry radius.

  Each instance corresponds to a separate set of coefficients, which are functions of the inputs.

  This struct conforms to the interface prescribed in math_helpers.hpp for scalar functions
  that are to be used with numerical rootfinding algorithms.

  The Kohler polynomial is a quartic polynomial whose variable is particle wet radius.
  Equilibrium solutions are roots of this polynomial.
  Algebraically, there are 2 complex roots and 2 real roots.  Of the real roots, one is
  positive, the other is negative.
  Physically, only the real, positive root makes sense.

  This struct is templated on scalar type so that it can be used with PackType.
  If it is used with PackType, each element of the PackType corresponds to a separate
  KohlerPolynomial, with distinct coefficients.
*/
template <typename T>
struct KohlerPolynomial {
  static_assert(std::is_floating_point<typename ekat::ScalarTraits<T>::scalar_type>::value, "floating point type");

  /// Minimum value of relative humidity
  static constexpr Real rel_humidity_min = 0.05;
  /// Above this relative humidity is considered saturated air, and cloud aerosol processes would apply
  static constexpr Real rel_humidity_max = 0.98;
  /// Minimum hygroscopicity for E3SM v1 aerosol species
  static constexpr Real hygro_min =  1e-10;
  /// Maximum hygroscopicity for E3SM v1 aerosol species
  static constexpr Real hygro_max = 1.3;
  /// Minimum particle size for E3SM v1
  static constexpr Real dry_radius_min_microns = 1e-3;
  /// Maximum particle size for Kohler theory
  static constexpr Real dry_radius_max_microns = 30;

  using value_type = T;
  using scalar_type = typename ekat::ScalarTraits<T>::scalar_type;

   /** Coefficient accounting for Kelvin effect on droplets

    In the documentation, this constant is denoted by A.

   @todo verify that the assumptions that lead to this constant still hold at SCREAM resolutions
  (e.g., constant temperature, constant surface tension) */
  static constexpr Real kelvin_droplet_effect_coeff = 0.00120746723156361711;

  /// Coefficient in the Kohler polynomial
  T log_rel_humidity;
  /// Coefficient in the Kohler polynomial
  T hygroscopicity;
  /// Coefficient in the Kohler polynomial
  T dry_radius_cubed;

  /** Constructor. Creates 1 instance of a KohlerPolynomial.

    @param rel_h relative humidity
    @param hygro hygroscopicity
    @param dry_rad_microns particle dry radius [ 1e-6 m ]
  */
  KOKKOS_INLINE_FUNCTION
  KohlerPolynomial(const T& rel_h, const T& hygro, const T& dry_rad_microns) :
    log_rel_humidity(log(rel_h)),
    hygroscopicity(hygro),
    dry_radius_cubed(cube(dry_rad_microns)) {
    EKAT_KERNEL_ASSERT(valid_inputs(rel_h, hygro, dry_rad_microns));
    }


  /** Evaluates the Kohler polynomial.

    @f$ K(r_w) = \log(s) r_w^4 - A r_w^3 + (B - \log(s))r_d^3 r_w + A r_d^3 @f$

    where r_w is the wet radius, s is relative humidity, A is the Kelvin effect coefficient,
    B is hygroscopicity, and r_d is the dry radius.

    @param [in] Polynomial input value, wet_radius wet radius in microns [ 1e-6 m]
    @return Polynomial value, wet_radius in microns [ 1e-6 m]
  */
  KOKKOS_INLINE_FUNCTION
  T operator() (const T& wet_radius) const {
    const T wet_radius_cubed = cube(wet_radius);
    const Real kelvinA = 0.00120746723156361711;
    T result = log_rel_humidity * wet_radius * wet_radius_cubed -
      kelvinA * wet_radius_cubed +
      (hygroscopicity * dry_radius_cubed - log_rel_humidity * dry_radius_cubed)*wet_radius +
      kelvinA * dry_radius_cubed;
    const auto nans = isnan(wet_radius);
    vector_simd for (int i=0; i<HAERO_PACK_SIZE; ++i) {
      if (nans[i]) {
        result[i] = 0;
      }
    }
    return result;
  }

  /** @brief Evaluates the derivative of the Kohler polynomial with respect to wet radius

    @f$ K'(r_w) = \frac{\partial K}(\partial r_w)(r_w) @f$

    @param [in] Polynomial input value, wet radius in microns [ 1e-6 m]
    @return Polynomial slope at input value
  */
  KOKKOS_INLINE_FUNCTION
  T derivative(const T& wet_radius) const {
    const T wet_radius_squared = square(wet_radius);
    const Real kelvinA = 0.00120746723156361711;
    T result = 4*log_rel_humidity*wet_radius*wet_radius_squared -
      3*kelvinA*wet_radius_squared +
      hygroscopicity*dry_radius_cubed - log_rel_humidity*dry_radius_cubed;
    const auto nans = isnan(wet_radius);
    vector_simd for (int i=0; i<HAERO_PACK_SIZE; ++i) {
      if (nans[i]) {
        result[i] = 0;
      }
    }
    return result;
  }

  KOKKOS_INLINE_FUNCTION
  bool valid_inputs(const T& relh, const T& hyg, const T& dry_rad) {
    return (FloatingPoint<T>::in_bounds(relh, rel_humidity_min, rel_humidity_max) and
            FloatingPoint<T>::in_bounds(hyg, hygro_min, hygro_max) and
            FloatingPoint<T>::in_bounds(dry_rad, dry_radius_min_microns, dry_radius_max_microns));
  }

  static std::string mathematica_verification_program(const int n, const std::string& output_dir = "~");
};


/** @brief This functor solves for the positive root of a packed instance of the KohlerPolynomial.
*/
struct KohlerNewtonSolve {
  const PackType& relative_humidity;
  const PackType& hygroscopicity;
  const PackType& dry_radius_microns;
  Real conv_tol;
  int n_iter;

  KOKKOS_INLINE_FUNCTION
  KohlerNewtonSolve(const PackType& rel_h, const PackType& hygro, const PackType& dry_rad, const Real tol) :
    relative_humidity(rel_h),
    hygroscopicity(hygro),
    dry_radius_microns(dry_rad),
    conv_tol(tol),
    n_iter(0)
    {}

  KOKKOS_INLINE_FUNCTION
  PackType operator() () {
    /// first guess at root needs to be sufficiently positive to avoid convergence to the negative root.
    /// avg. case suggests an initial guess of ~ 10 * dry_radius is good and convergence is achieved within < 10 iterations
    /// worst case (large dry radius, near-saturated air) may require initial guesses ~ 25 * dry_radius
    /// using the worst case value for all guesses causes the number of iterations to go past 30, reducing performance.
    /// hence, we use the avg case guess and add a guard to increase the initial root if the algorithm does not converge
    /// to the Kohler polynomial's positive root.
    PackType wet_radius_init(25*dry_radius_microns);
    PackType result(-1);
    bool keep_going = true;
    while ( keep_going ) {
      auto solver = math::ScalarNewtonSolver<KohlerPolynomial<PackType>>(wet_radius_init, conv_tol,
        KohlerPolynomial<PackType>(relative_humidity, hygroscopicity, dry_radius_microns));
      result = solver.solve();
      n_iter += solver.counter;
      keep_going = (result < 0).any();
      /// iterations must converge
      EKAT_KERNEL_ASSERT(solver.counter < math::ScalarNewtonSolver<KohlerPolynomial<PackType>>::max_iter);
      /// if converged to the wrong root, try a bigger initial guess
      if (keep_going) {
        wet_radius_init *= 2;
      }
    }
    return result;
  }
};


struct KohlerBisectionSolve {
  const PackType& relative_humidity;
  const PackType& hygroscopicity;
  const PackType& dry_radius_microns;
  Real conv_tol;
  int n_iter;

  KOKKOS_INLINE_FUNCTION
  KohlerBisectionSolve(const PackType& rel_h, const PackType& hygro, const PackType& dry_rad, const Real tol) :
    relative_humidity(rel_h),
    hygroscopicity(hygro),
    dry_radius_microns(dry_rad),
    conv_tol(tol),
    n_iter(0)
    {}

  KOKKOS_INLINE_FUNCTION
  PackType operator() () {
    PackType left_endpt(0.9*dry_radius_microns);
    PackType right_endpt(200.0);
    PackType result(-1);
    auto solver = math::BisectionSolver<KohlerPolynomial<PackType>>(left_endpt, right_endpt, conv_tol,
      KohlerPolynomial<PackType>(relative_humidity, hygroscopicity, dry_radius_microns));
    result = solver.solve();
    n_iter = solver.counter;
    return result;
  }
};

} // namespace haero
#endif
