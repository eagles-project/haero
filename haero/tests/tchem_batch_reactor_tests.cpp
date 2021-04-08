#include "haero/model.hpp"
#include "chemUtil.hpp"
#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"

#include <sys/stat.h>
#include <sys/types.h>

using namespace haero;
using namespace chemUtil;

TEST_CASE("toy problem chemistry", "haero_unit_tests")
{

  real_type lat = 20.0;
  real_type lon = 37.612;

  const real_type latc(20.0); //! degrees
  const real_type lonc(300.0);// ! degrees
  real_type k1 = ats<real_type>::sin(lat*PI()/180) * ats<real_type>::sin(latc*PI()/180) +
                       ats<real_type>::cos(lat*PI()/180) * ats<real_type>::cos(latc*PI()/180) *
                       ats<real_type>::cos(lon*PI()/180 - lonc*PI()/180);
  k1 = k1 > 0 ? k1 : 0;
  real_type k2 = 1;

  // Write out some test data to our current working directory.
  const char* chem_inp = R"INPUT(ELEMENTS
X /1/
END
SPECIES
X X2
END
THERM ALL
    300.000  1000.000  5000.000
X                        X  1               G   200.000  6000.000 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967000E+00                   4
X2                       X  2               G   200.000  6000.000 1000.00      1
 2.50000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00 0.00000000E+00    2
-7.45375000E+02 4.37967000E+00 2.50000000E+00 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00-7.45375000E+02 4.37967000E+00                   4
END
REACTIONS
X2=>2X      17E+0    1   1
X+X=>X2      21E+0    1   1
END
)INPUT";
  mkdir("data", 0777);
  mkdir("data/toy-problem", 0777);
  FILE* f = fopen("data/toy-problem/chem.inp", "w");
  fprintf(f, "%s", chem_inp);
  fclose(f);
  f = fopen("data/toy-problem/therm.dat", "w");
  fclose(f);

  // Run the test.
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
  const Real val00 = results_host(0, 0);
  const Real val01 = results_host(0, 1);
  REQUIRE(FloatingPoint<Real>::equiv(val00, 3.27585e-11,
                                     FloatingPoint<Real>::zero_tol));
  REQUIRE(FloatingPoint<Real>::equiv(val01, -1.63793e-11,
                                     FloatingPoint<Real>::zero_tol));

  } // kokkos scope
  Kokkos::finalize();
}

