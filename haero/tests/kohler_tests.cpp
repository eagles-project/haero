#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

#include "catch2/catch.hpp"
#include "ekat/util/ekat_math_utils.hpp"
#include "haero/diagnostics/kohler_solve_diagnostic.hpp"
#include "haero/haero.hpp"
#include "haero/math_helpers.hpp"
#include "haero/utils.hpp"

using namespace haero;

struct KohlerTestInput {
  DeviceType::view_1d<PackType> relative_humidity;
  DeviceType::view_1d<PackType> hygroscopicity;
  DeviceType::view_1d<PackType> dry_radius;

  typename DeviceType::view_1d<PackType>::HostMirror h_relative_humidity;
  typename DeviceType::view_1d<PackType>::HostMirror h_hygroscopicity;
  typename DeviceType::view_1d<PackType>::HostMirror h_dry_radius;

  explicit KohlerTestInput(const int n)
      : relative_humidity("relative_humidity", PackInfo::num_packs(cube(n))),
        hygroscopicity("hygroscopicity", PackInfo::num_packs(cube(n))),
        dry_radius("dry_radius", PackInfo::num_packs(cube(n))) {
    EKAT_REQUIRE(n > 1);
    const Real drelh = (KohlerPolynomial<PackType>::rel_humidity_max -
                        KohlerPolynomial<PackType>::rel_humidity_min) /
                       (n - 1);
    const Real dhyg = (KohlerPolynomial<PackType>::hygro_max -
                       KohlerPolynomial<PackType>::hygro_min) /
                      (n - 1);
    const Real ddry = (KohlerPolynomial<PackType>::dry_radius_max_microns -
                       KohlerPolynomial<PackType>::dry_radius_min_microns) /
                      (n - 1);

    h_relative_humidity = Kokkos::create_mirror_view(relative_humidity);
    h_hygroscopicity = Kokkos::create_mirror_view(hygroscopicity);
    h_dry_radius = Kokkos::create_mirror_view(dry_radius);

    int ind = 0;
    for (int i = 0; i < n; ++i) {
      const Real rel_h =
          KohlerPolynomial<PackType>::rel_humidity_min + i * drelh;
      for (int j = 0; j < n; ++j) {
        const Real hyg = KohlerPolynomial<PackType>::hygro_min + j * dhyg;
        for (int k = 0; k < n; ++k) {
          const Real drad =
              KohlerPolynomial<PackType>::dry_radius_min_microns + k * ddry;
          const int pack_idx = PackInfo::pack_idx(ind);
          const int vec_idx = PackInfo::vec_idx(ind++);
          h_relative_humidity(pack_idx)[vec_idx] = rel_h;
          h_hygroscopicity(pack_idx)[vec_idx] = hyg;
          h_dry_radius(pack_idx)[vec_idx] = drad;
        }
      }
    }

    Kokkos::deep_copy(relative_humidity, h_relative_humidity);
    Kokkos::deep_copy(hygroscopicity, h_hygroscopicity);
    Kokkos::deep_copy(dry_radius, h_dry_radius);
  }
};

struct KohlerTestFunctor {
  DeviceType::view_1d<PackType> newton_sol;
  DeviceType::view_1d<PackType> newton_err;
  DeviceType::view_1d<int> newton_iterations;
  DeviceType::view_1d<PackType> bisection_sol;
  DeviceType::view_1d<PackType> bisection_err;
  DeviceType::view_1d<int> bisection_iterations;
  DeviceType::view_1d<PackType> bracket_sol;
  DeviceType::view_1d<PackType> bracket_err;
  DeviceType::view_1d<int> bracket_iterations;
  DeviceType::view_1d<PackType> rh_in;
  DeviceType::view_1d<PackType> hyg_in;
  DeviceType::view_1d<PackType> dry_rad;
  DeviceType::view_1d<PackType> true_sol;
  DeviceType::view_1d<MaskType> pack_masks;
  Real tol;

  KohlerTestFunctor(DeviceType::view_1d<PackType> nsol,
                    DeviceType::view_1d<PackType> nerr,
                    DeviceType::view_1d<int> niter,
                    DeviceType::view_1d<PackType> bsol,
                    DeviceType::view_1d<PackType> berr,
                    DeviceType::view_1d<int> biter,
                    DeviceType::view_1d<PackType> brksol,
                    DeviceType::view_1d<PackType> brkerr,
                    DeviceType::view_1d<int> brkiter,
                    const DeviceType::view_1d<PackType> rh,
                    const DeviceType::view_1d<PackType> hyg,
                    const DeviceType::view_1d<PackType> drad,
                    const DeviceType::view_1d<PackType> tsol,
                    const DeviceType::view_1d<MaskType> masks, const Real ctol)
      : newton_sol(nsol),
        newton_err(nerr),
        newton_iterations(niter),
        bisection_sol(bsol),
        bisection_err(berr),
        bisection_iterations(biter),
        bracket_sol(brksol),
        bracket_err(brkerr),
        bracket_iterations(brkiter),
        rh_in(rh),
        hyg_in(hyg),
        dry_rad(drad),
        true_sol(tsol),
        pack_masks(masks),
        tol(ctol) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int pack_idx) const {
    KohlerNewtonSolve newton(rh_in(pack_idx), hyg_in(pack_idx),
                             dry_rad(pack_idx), tol);
    newton_sol(pack_idx) = newton();
    newton_err(pack_idx) = abs(newton_sol(pack_idx) - true_sol(pack_idx))/true_sol(pack_idx);
    newton_iterations(pack_idx) = newton.n_iter;
    KohlerBisectionSolve bisection(rh_in(pack_idx), hyg_in(pack_idx),
                                   dry_rad(pack_idx), tol);
    bisection_sol(pack_idx) = bisection();
    bisection_err(pack_idx) = abs(bisection_sol(pack_idx) - true_sol(pack_idx))/true_sol(pack_idx);
    bisection_iterations(pack_idx) = bisection.n_iter;

    KohlerBracketedNewtonSolve bracket(rh_in(pack_idx), hyg_in(pack_idx), dry_rad(pack_idx), tol);
    bracket_sol(pack_idx) = bracket();
    bracket_err(pack_idx) = abs(bracket_sol(pack_idx) - true_sol(pack_idx))/true_sol(pack_idx);
    bracket_iterations(pack_idx) = bracket.n_iter;

    vector_simd for (int s=0; s<HAERO_PACK_SIZE; ++s) {
      if (!pack_masks(pack_idx)[s]) {
        newton_err(pack_idx)[s] = 0;
        bisection_err(pack_idx)[s] = 0;
        bracket_err(pack_idx)[s] = 0;
      }
    }

#ifndef HAERO_USE_CUDA
    if (false /*(newton_err(pack_idx) > tol).any()*/) {
      std::cout << "error exceeds tolerance at pack " << pack_idx << "\n";
      std::cout << "\tarray indices: ";
      for (int i = 0; i < HAERO_PACK_SIZE; ++i) {
        std::cout << PackInfo::array_idx(pack_idx, i) << " ";
      }
      std::cout << "\n\trh " << rh_in(pack_idx) << " hyg " << hyg_in(pack_idx)
                << " dry_rad " << dry_rad(pack_idx) << "\n";
      std::cout << "\tnewton sol = " << newton_sol(pack_idx)
                << " n_iter = " << newton_iterations(pack_idx) << "\n";
      std::cout << "\ttrue_sol = " << true_sol(pack_idx) << "\n";
    }
#endif

    if (pack_idx == rh_in.extent(0) - 1) {
      ekat_masked_loop(!pack_masks(pack_idx), s) {
        newton_err(pack_idx)[s] = 0;
        bisection_err(pack_idx)[s] = 0;
      }
    }
  }
};

template <typename ScalarType>
struct PackMaxReduce {
  DeviceType::view_1d<ekat::Pack<ScalarType, HAERO_PACK_SIZE>> view;

  explicit PackMaxReduce(
      DeviceType::view_1d<ekat::Pack<ScalarType, HAERO_PACK_SIZE>> v)
      : view(v) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int pack_idx, ScalarType& s) const {
    vector_disabled for (int i = 0; i < HAERO_PACK_SIZE; ++i) {
      s = (s > view(pack_idx)[i] ? s : view(pack_idx)[i]);
    }
  }
};

TEST_CASE("KohlerPolynomial properties", "") {
  // Generate input data
  static constexpr int N = 20;
  const KohlerTestInput test_inputs(N);

  const int num_packs = PackInfo::num_packs(cube(N));

  DeviceType::view_1d<PackType> kohler_at_zero("kohler_at_zero", num_packs);
  DeviceType::view_1d<PackType> kohler_at_rdry("kohler_at_rdry", num_packs);
  DeviceType::view_1d<PackType> kohler_at_25rdry("kohler_at_25rdry", num_packs);

  Kokkos::parallel_for(
      num_packs, KOKKOS_LAMBDA(const int pack_idx) {
        const auto poly =
            KohlerPolynomial<PackType>(test_inputs.relative_humidity(pack_idx),
                                       test_inputs.hygroscopicity(pack_idx),
                                       test_inputs.dry_radius(pack_idx));

        kohler_at_zero(pack_idx) = poly(PackType(0));
        kohler_at_rdry(pack_idx) = poly(test_inputs.dry_radius(pack_idx));
        kohler_at_25rdry(pack_idx) =
            poly(25 * test_inputs.dry_radius(pack_idx));
      });

  auto h_kohler_at_zero = Kokkos::create_mirror_view(kohler_at_zero);
  auto h_kohler_at_rdry = Kokkos::create_mirror_view(kohler_at_rdry);
  auto h_kohler_at_25rdry = Kokkos::create_mirror_view(kohler_at_25rdry);
  auto h_rdry = Kokkos::create_mirror_view(test_inputs.dry_radius);
  auto h_hygro = Kokkos::create_mirror_view(test_inputs.hygroscopicity);
  auto h_relh = Kokkos::create_mirror_view(test_inputs.relative_humidity);
  Kokkos::deep_copy(h_kohler_at_zero, kohler_at_zero);
  Kokkos::deep_copy(h_kohler_at_rdry, kohler_at_rdry);
  Kokkos::deep_copy(h_kohler_at_25rdry, kohler_at_25rdry);
  Kokkos::deep_copy(h_rdry, test_inputs.dry_radius);
  Kokkos::deep_copy(h_hygro, test_inputs.hygroscopicity);
  Kokkos::deep_copy(h_relh, test_inputs.relative_humidity);

  for (int i = 0; i < num_packs; ++i) {
    REQUIRE(FloatingPoint<PackType>::equiv(
        h_kohler_at_zero(i), kelvin_droplet_effect_coeff * cube(h_rdry(i))));

    REQUIRE( FloatingPoint<PackType>::rel(h_kohler_at_rdry(i), h_rdry(i)*cube(h_rdry(i))* h_hygro(i) , 4.75e-6) );
    REQUIRE( (h_kohler_at_rdry(i) > 0).all());
    REQUIRE((h_kohler_at_25rdry(i) < 0).all());
  }
}

TEST_CASE("KohlerSolve-verification", "") {
  /// Generate input data
  static constexpr int N = 20;
  const KohlerTestInput test_inputs(N);

  const int num_packs = PackInfo::num_packs(cube(N));

  const Real conv_tol = (std::is_same<float, Real>::value ? 1.0e-3 : 1.0e-10);

  std::cout << "generating 3-parameter sweep for the Kohler solve with "
            << cube(N) << " trials.\n";

  DeviceType::view_1d<PackType> newton_sol("newton_sol", num_packs);
  DeviceType::view_1d<PackType> newton_error("newton_error", num_packs);
  DeviceType::view_1d<int> newton_iterations("newton_iterations", num_packs);
  DeviceType::view_1d<PackType> bisection_sol("bisection_sol", num_packs);
  DeviceType::view_1d<PackType> bisection_error("bisection_error", num_packs);
  DeviceType::view_1d<int> bisection_iterations("bisection_iterations",
                                                num_packs);
  DeviceType::view_1d<PackType> true_sol("true_sol", num_packs);
  DeviceType::view_1d<MaskType> pack_masks("pack_masks", num_packs);
  DeviceType::view_1d<PackType> bracket_sol("bracket_sol", num_packs);
  DeviceType::view_1d<PackType> bracket_error("bracket_error", num_packs);
  DeviceType::view_1d<int> bracket_iterations("bracket_iterations", num_packs);

  Kokkos::parallel_for(
      num_packs, KOKKOS_LAMBDA(const int pack_idx) {
        if (pack_idx < num_packs - 1) {
          pack_masks(pack_idx) = MaskType(true);
        } else {
          vector_disabled for (int i = 0; i < PackInfo::last_vec_end(cube(N));
                               ++i) {
            pack_masks(pack_idx).set(i, true);
          }
          vector_disabled for (int i = PackInfo::last_vec_end(cube(N));
                               i < HAERO_PACK_SIZE; ++i) {
            pack_masks(pack_idx).set(i, false);
          }
        }
      });

  std::cout << "returned from padding mask init.\n";

  std::string dfile = HAERO_TEST_DATA_DIR;
  dfile += "/mm_kohler_roots.txt";
  std::cout << "reading true solutions from data file: " << dfile << "\n";
  std::ifstream mm_sols(dfile);
  auto h_true_sol = Kokkos::create_mirror_view(true_sol);
  REQUIRE(mm_sols.is_open());
  Real mmroot;
  int idx = 0;
  while (mm_sols >> mmroot) {
    const int pack_idx = PackInfo::pack_idx(idx);
    const int vec_idx = PackInfo::vec_idx(idx++);
    h_true_sol(pack_idx)[vec_idx] = mmroot;
  }
  Kokkos::deep_copy(true_sol, h_true_sol);
  mm_sols.close();

  std::cout << "launching test kernel\n";
  Kokkos::parallel_for(
      "KohlerVerificatioTest", num_packs,
      KohlerTestFunctor(newton_sol, newton_error, newton_iterations,
                        bisection_sol, bisection_error, bisection_iterations,
                        bracket_sol, bracket_error, bracket_iterations,
                        test_inputs.relative_humidity,
                        test_inputs.hygroscopicity, test_inputs.dry_radius,
                        true_sol, pack_masks, conv_tol));

  Real max_err_newton;
  Real max_err_bisection;
  Real max_err_bracket;
  int max_iter_newton;
  int max_iter_bisection;
  int max_iter_bracket;

  std::cout << "computing reductions to generate error statistics\n";
  Kokkos::parallel_reduce(num_packs, PackMaxReduce<Real>(newton_error),
                          Kokkos::Max<Real>(max_err_newton));
  Kokkos::parallel_reduce(num_packs, PackMaxReduce<Real>(bisection_error),
                          Kokkos::Max<Real>(max_err_bisection));
  Kokkos::parallel_reduce(num_packs, PackMaxReduce<Real>(bracket_error),
    Kokkos::Max<Real>(max_err_bracket));
  Kokkos::parallel_reduce(
      num_packs,
      KOKKOS_LAMBDA(const int pack_idx, int& ctr) {
        ctr = (ctr > newton_iterations(pack_idx) ? ctr
                                                 : newton_iterations(pack_idx));
      },
      Kokkos::Max<int>(max_iter_newton));
  Kokkos::parallel_reduce(
      num_packs,
      KOKKOS_LAMBDA(const int pack_idx, int& ctr) {
        ctr = (ctr > bisection_iterations(pack_idx)
                   ? ctr
                   : bisection_iterations(pack_idx));
      },
      Kokkos::Max<int>(max_iter_bisection));
  Kokkos::parallel_reduce(num_packs, KOKKOS_LAMBDA (const int pack_idx, int& ctr) {
      ctr = (ctr > bracket_iterations(pack_idx) ? ctr : bracket_iterations(pack_idx));
    }, Kokkos::Max<int>(max_iter_bracket));

  std::cout << "Newton solve:\n";
  std::cout << "\t max err = " << max_err_newton << "\n";
  std::cout << "\tmax iter = " << max_iter_newton << "\n";
  std::cout << "Bisection solve:\n";
  std::cout << "\t max err = " << max_err_bisection << "\n";
  std::cout << "\tmax iter = " << max_iter_bisection << "\n";
  std::cout << "Bracket solve:\n";
  std::cout << "\t max err = " << max_err_bracket << "\n";
  std::cout << "\tmax iter = " << max_iter_bracket << "\n";

#if HAERO_DOUBLE_PRECISION || !defined(NDEBUG)
  REQUIRE(max_err_newton < 1.5 * conv_tol);
#else
  std::cout
      << "DISABLED Newton Solver test due to single precision/release build\n";
#endif
  REQUIRE(max_err_bisection < 2.3 * conv_tol);
  REQUIRE(max_err_bracket < 2.3 * conv_tol);

  std::cout << "To generate the verification data with Mathematica, run this "
               "program:\n\n";

  std::cout << KohlerPolynomial<PackType>::mathematica_verification_program(N);

  std::cout << "\n\nTo generate the verification data with Matlab, run this "
               "program:\n\n";

  std::cout << KohlerPolynomial<PackType>::matlab_verification_program(N);
}
