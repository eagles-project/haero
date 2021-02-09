#include "haero/model.hpp"
#include "chemUtil.hpp"
#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"

using namespace haero;
using namespace chemUtil;

TEST_CASE("toy problem chemistry", "haero_unit_tests")
{
  Kokkos::initialize();
  {

  // arguments to the constructor are: 1) directory containing chem files
                                    // 2) detail (boolean)
                                    // 3) nBatch (int)
                                    // 4) verbose (boolean)
                                    // 5) theta/latitude (real)
                                    // 6) lambda/longitude (real)
                                    // 7) initial mass X (real)
                                    // 8) initial mass X2 (real)

  // provide arbitrarily chosen inputs
  chemSolver chem_solver("toy-problem/", false, 1, false, 20.0, 37.612, 2.5e-7,
                         1.8e-6);

  // run the problem on device
  real_type_2d_view results = chem_solver.get_results();

  // create mirror view and deep copy to host
  auto results_host = Kokkos::create_mirror_view(results);
  Kokkos::deep_copy(results_host, results);

  // verify that known inputs result in known outputs
  REQUIRE(FloatingPoint<Real>::equiv(results_host(0, 0), 3.27585e-11,
                                     FloatingPoint<Real>::zero_tol));
  REQUIRE(FloatingPoint<Real>::equiv(results_host(0, 1), -1.63793e-11,
                                     FloatingPoint<Real>::zero_tol));

  } // kokkos scope
  Kokkos::finalize();
}

