#include "haero/model.hpp"
#include "chemUtil.hpp"
#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"

#include <sys/stat.h>
#include <sys/types.h>

using namespace haero;
using namespace chemUtil;

TEST_CASE("TChem tendency computation tests", "haero_unit_tests"){

  // lat/lon for sun zenith
  const Real latz(20.0); // ! degrees
  const Real lonz(300.0);// ! degrees

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
X2=>2X      1E+0    1   1
X+X=>X2      1E+0    1   1
END
  )INPUT";
  mkdir("data", 0777);
  mkdir("data/toy-problem", 0777);
  FILE* f = fopen("data/toy-problem/chem.inp", "w");
  fprintf(f, "%s", chem_inp);
  fclose(f);
  f = fopen("data/toy-problem/therm.dat", "w");
  fclose(f);

  SECTION("light side of terminator"){

    // lat/lon for column position
    Real lat = 20.0;
    Real lon = 37.5;
    // initial concentrations for the two species
    Real initX = 1.0e-6;
    Real initX2 = 1.0e-6;

    // calculate k1 reaction rate, based on position of column and sun's zenith
    Real k1 = sin(lat*PI()/180) * sin(latz*PI()/180) +
                       cos(lat*PI()/180) * cos(latz*PI()/180) *
                       cos(lon*PI()/180 - lonz*PI()/180);
    k1 = k1 > 0 ? k1 : 0;
    Real k2 = 1;

    // arguments to the constructor are: 1) directory containing chem files
                                      // 2) detail (boolean)
                                      // 3) nBatch (int)
                                      // 4) verbose (boolean)
                                      // 5) initial mass X (real)
                                      // 6) initial mass X2 (real)

    // provide arbitrarily chosen inputs
    chemSolver chem_solver("toy-problem/", false, 1, false, k1, k2,
                           initX, initX2);

    // run the problem on device
    real_type_2d_view results = chem_solver.get_results();

    // create mirror view and deep copy to host
    auto results_host = Kokkos::create_mirror_view(results);
    Kokkos::deep_copy(results_host, results);

    // eq (4) tendency should be positive, eq (5) negative
    const Real val00 = results_host(0, 0);
    const Real val01 = results_host(0, 1);
    REQUIRE(FloatingPoint<Real>::in_bounds(val00, 0.0, 1.0e14,
                                           FloatingPoint<Real>::zero_tol));
    REQUIRE(FloatingPoint<Real>::in_bounds(val01, -1.0e14, 0.0,
                                           FloatingPoint<Real>::zero_tol));

  }

  SECTION("dark side of terminator"){

    // lat/lon for column position
    Real lat = 20.0;
    Real lon = 50.0;
    // initial concentrations for the two species
    Real initX = 1.0e-6;
    Real initX2 = 1.0e-6;

    // calculate k1 reaction rate, based on position of column and sun's zenith
    Real k1 = sin(lat*PI()/180) * sin(latz*PI()/180) +
                       cos(lat*PI()/180) * cos(latz*PI()/180) *
                       cos(lon*PI()/180 - lonz*PI()/180);
    k1 = k1 > 0 ? k1 : 0;
    Real k2 = 1;

    // arguments to the constructor are: 1) directory containing chem files
                                      // 2) detail (boolean)
                                      // 3) nBatch (int)
                                      // 4) verbose (boolean)
                                      // 5) initial mass X (real)
                                      // 6) initial mass X2 (real)

    // provide arbitrarily chosen inputs
    chemSolver chem_solver("toy-problem/", false, 1, false, k1, k2,
                           initX, initX2);

    // run the problem on device
    real_type_2d_view results = chem_solver.get_results();

    // create mirror view and deep copy to host
    auto results_host = Kokkos::create_mirror_view(results);
    Kokkos::deep_copy(results_host, results);

    // eq (4) tendency should be negative, eq (5) positive
    const Real val00 = results_host(0, 0);
    const Real val01 = results_host(0, 1);
    REQUIRE(FloatingPoint<Real>::in_bounds(val00, -1.0e14, 0.0,
                                           FloatingPoint<Real>::zero_tol));
    REQUIRE(FloatingPoint<Real>::in_bounds(val01, 0.0, 1.0e14,
                                           FloatingPoint<Real>::zero_tol));

  }

  SECTION("zero tendencies"){

    // lat/lon for column position
    Real lat = 20.0;
    Real lon = 37.5;

    // calculate k1 reaction rate, based on position of column and sun's zenith
    Real k1 = sin(lat*PI()/180) * sin(latz*PI()/180) +
                       cos(lat*PI()/180) * cos(latz*PI()/180) *
                       cos(lon*PI()/180 - lonz*PI()/180);
    k1 = k1 > 0 ? k1 : 0;
    Real k2 = 1;

    // initial concentrations for the two species
    // calculate initial X2 such that tendencies will be zero
    Real initX = 1.0e-6;
    Real initX2 = (k2 / k1) * pow(initX, 2);

    // arguments to the constructor are: 1) directory containing chem files
                                      // 2) detail (boolean)
                                      // 3) nBatch (int)
                                      // 4) verbose (boolean)
                                      // 5) initial mass X (real)
                                      // 6) initial mass X2 (real)

    // provide arbitrarily chosen inputs
    chemSolver chem_solver("toy-problem/", false, 1, false, k1, k2,
                           initX, initX2);

    // run the problem on device
    real_type_2d_view results = chem_solver.get_results();

    // create mirror view and deep copy to host
    auto results_host = Kokkos::create_mirror_view(results);
    Kokkos::deep_copy(results_host, results);

    // verify that both tendencies are zero
    const Real val00 = results_host(0, 0);
    const Real val01 = results_host(0, 1);
    REQUIRE(FloatingPoint<Real>::zero(val00, FloatingPoint<Real>::zero_tol));
    REQUIRE(FloatingPoint<Real>::zero(val01, FloatingPoint<Real>::zero_tol));

  }
}


