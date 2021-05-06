#include "haero/haero.hpp"
#include "haero/diagnostics/kohler_solve_diagnostic.hpp"
#include "haero/math_helpers.hpp"
#include "haero/utils.hpp"
#include "ekat/util/ekat_math_utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace haero;

std::string mathematica_verification_program();

struct KohlerTestInput {
  DeviceType::view_1d<PackType> relative_humidity;
  DeviceType::view_1d<PackType> hygroscopicity;
  DeviceType::view_1d<PackType> dry_radius;

  typename DeviceType::view_1d<PackType>::HostMirror h_relative_humidity;
  typename DeviceType::view_1d<PackType>::HostMirror h_hygroscopicity;
  typename DeviceType::view_1d<PackType>::HostMirror h_dry_radius;

  explicit KohlerTestInput(const int n) :
    relative_humidity("relative_humidity", PackInfo::num_packs(cube(n))),
    hygroscopicity("hygroscopicity", PackInfo::num_packs(cube(n))),
    dry_radius("dry_radius", PackInfo::num_packs(cube(n))) {
      EKAT_REQUIRE(n>1);
      const Real drelh = (KohlerPolynomial<PackType>::rel_humidity_max -
        KohlerPolynomial<PackType>::rel_humidity_min)/(n-1);
      const Real dhyg = (KohlerPolynomial<PackType>::hygro_max -
        KohlerPolynomial<PackType>::hygro_min)/(n-1);
      const Real ddry = (KohlerPolynomial<PackType>::dry_radius_max_microns -
        KohlerPolynomial<PackType>::dry_radius_min_microns)/(n-1);

      h_relative_humidity = Kokkos::create_mirror_view(relative_humidity);
      h_hygroscopicity = Kokkos::create_mirror_view(hygroscopicity);
      h_dry_radius = Kokkos::create_mirror_view(dry_radius);

      int ind=0;
      for (int i=0; i<n; ++i) {
        const Real rel_h = KohlerPolynomial<PackType>::rel_humidity_min + i * drelh;
        for (int j=0; j<n; ++j) {
          const Real hyg = KohlerPolynomial<PackType>::hygro_min + j * dhyg;
          for (int k=0; k<n; ++k) {
            const Real drad = KohlerPolynomial<PackType>::dry_radius_min_microns + k * ddry;
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
  DeviceType::view_1d<PackType> rh_in;
  DeviceType::view_1d<PackType> hyg_in;
  DeviceType::view_1d<PackType> dry_rad;
  DeviceType::view_1d<PackType> true_sol;
  Real tol;

  KohlerTestFunctor(DeviceType::view_1d<PackType> nsol, DeviceType::view_1d<PackType> nerr,
    DeviceType::view_1d<int> niter, DeviceType::view_1d<PackType> bsol, DeviceType::view_1d<PackType> berr,
    DeviceType::view_1d<int> biter, const DeviceType::view_1d<PackType> rh,
    const DeviceType::view_1d<PackType> hyg, const DeviceType::view_1d<PackType> drad,
    const DeviceType::view_1d<PackType> tsol, const Real ctol) :
    newton_sol(nsol),
    newton_err(nerr),
    newton_iterations(niter),
    bisection_sol(bsol),
    bisection_err(berr),
    bisection_iterations(biter),
    rh_in(rh),
    hyg_in(hyg),
    dry_rad(drad),
    true_sol(tsol),
    tol(ctol) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int pack_idx) const {
    KohlerNewtonSolve newton(rh_in(pack_idx), hyg_in(pack_idx), dry_rad(pack_idx), tol);
    newton_sol(pack_idx) = newton();
    newton_err(pack_idx) = abs(newton_sol(pack_idx) - true_sol(pack_idx));
    newton_iterations(pack_idx) = newton.n_iter;
    KohlerBisectionSolve bisection(rh_in(pack_idx), hyg_in(pack_idx), dry_rad(pack_idx), tol);
    bisection_sol(pack_idx) = bisection();
    bisection_err(pack_idx) = abs(bisection_sol(pack_idx) - true_sol(pack_idx));
    bisection_iterations(pack_idx) = bisection.n_iter;
  }
};

template <typename ScalarType>
struct PackMaxReduce {
  DeviceType::view_1d<ekat::Pack<ScalarType,HAERO_PACK_SIZE>> view;

  explicit PackMaxReduce(DeviceType::view_1d<ekat::Pack<ScalarType,HAERO_PACK_SIZE>> v) :
    view(v) {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int pack_idx, ScalarType& s) const {
    vector_simd for (int i=0; i<HAERO_PACK_SIZE; ++i) {
      s = (s > view(pack_idx)[i] ? s : view(pack_idx)[i]);
    }
  }
};

TEST_CASE("KohlerSolve-verification", "") {
  /// Generate input data
  static constexpr int N = 20;
  const KohlerTestInput test_inputs(N);

  const int num_packs = PackInfo::num_packs(cube(N));

  const Real conv_tol = (std::is_same<float,Real>::value ? 1.0e-3 : 1.0e-10);

  DeviceType::view_1d<PackType> newton_sol("newton_sol", num_packs);
  DeviceType::view_1d<PackType> newton_error("newton_error", num_packs);
  DeviceType::view_1d<int> newton_iterations("newton_iterations", num_packs);
  DeviceType::view_1d<PackType> bisection_sol("bisection_sol", num_packs);
  DeviceType::view_1d<PackType> bisection_error("bisection_error", num_packs);
  DeviceType::view_1d<int> bisection_iterations("bisection_iterations", num_packs);
  DeviceType::view_1d<PackType> true_sol("true_sol", num_packs);

  std::string dfile = HAERO_TEST_DATA_DIR;
  dfile += "/mm_kohler_roots.txt";
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

  Kokkos::parallel_for("KohlerVerificatioTest", num_packs,
    KohlerTestFunctor(newton_sol, newton_error, newton_iterations,
        bisection_sol, bisection_error, bisection_iterations,
        test_inputs.relative_humidity, test_inputs.hygroscopicity, test_inputs.dry_radius,
        true_sol, conv_tol));

  Real max_err_newton;
  Real max_err_bisection;
  int max_iter_newton;
  int max_iter_bisection;

  Kokkos::parallel_reduce(num_packs, PackMaxReduce<Real>(newton_error), Kokkos::Max<Real>(max_err_newton));
  Kokkos::parallel_reduce(num_packs, PackMaxReduce<Real>(bisection_error), Kokkos::Max<Real>(max_err_bisection));
  Kokkos::parallel_reduce(num_packs, KOKKOS_LAMBDA (const int pack_idx, int& ctr) {
    ctr = (ctr > newton_iterations(pack_idx) ? ctr : newton_iterations(pack_idx));
  }, Kokkos::Max<int>(max_iter_newton));
  Kokkos::parallel_reduce(num_packs, KOKKOS_LAMBDA (const int pack_idx, int& ctr) {
    ctr = (ctr > bisection_iterations(pack_idx) ? ctr : bisection_iterations(pack_idx));
  }, Kokkos::Max<int>(max_iter_bisection));

  std::cout << "Newton solve:\n";
  std::cout << "\t max err = " << max_err_newton << "\n";
  std::cout << "\tmax iter = " << max_iter_newton << "\n";
  std::cout << "Bisection solve:\n";
  std::cout << "\t max err = " << max_err_bisection << "\n";
  std::cout << "\tmax iter = " << max_iter_bisection << "\n";

  REQUIRE(max_err_newton < conv_tol);
  REQUIRE(max_err_bisection < conv_tol);

}



