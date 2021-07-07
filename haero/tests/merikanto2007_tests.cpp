#include <cstdio>

#include "catch2/catch.hpp"
#include "haero/processes/merikanto2007.hpp"

using namespace haero;

// These tests aren't really unit tests in the usual sense--they are "smoke
// tests" in that they exercise the parameterizations of Merikanto et al
// (2007) and generate plot files that can be used to compare the results of
// the parameterizations with the calculations from classical nucleation
// theory (e.g. given in Chapter 11 of the Third Edition of Seinfeld and
// Pandis, 2016).

// Compare the output of this test with Merikanto et al (2007), figure 2.
// Because NH3 affects H2SO4 nucleation, the figure has a few different sets of
// curves
TEST_CASE("merikanto2007_figure_2") {
  PackType temps[2] = {235.15, 273.15};
  PackType xi_nh3s[2][2] = {
      // temperature dependent
      {0.1, 10},   // T = 235.15 K
      {10, 1000},  // T = 273.15 K
  };
  PackType c_h2so4s[2][2][10] = {
      // temperature and NH3 dependent
      {
          // T = 235.15 K
          {2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9},  // xi = 0.1
          {1e5, 2e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9}   // xi = 10
      },
      {
          // T = 273.15 K
          {1e8, 2e8, 3e8, 4e8, 5e8, 6e8, 7e8, 8e8, 9e8, 1e9},  // xi = 10
          {5e7, 8e7, 1e8, 2e8, 5e8, 6e8, 7e8, 8e8, 9e8, 1e9}   // xi = 1000
      }};
  PackType rel_hum(0.5);
  FILE* fp = fopen("merikanto2007_figure_2.dat", "w");
  for (int i = 0; i < 2; ++i) {
    PackType temp = temps[i];
    for (int j = 0; j < 2; ++j) {
      PackType xi_nh3 = xi_nh3s[i][j];
      for (int k = 0; k < 10; ++k) {
        PackType c_h2so4 = c_h2so4s[i][j][k];
        auto log_J =
            merikanto2007::log_nucleation_rate(temp, rel_hum, c_h2so4, xi_nh3);
        fprintf(fp, "%g %g %g %g\n", temp[0], xi_nh3[0], c_h2so4[0],
                exp(log_J)[0]);
      }
    }
  }
  fclose(fp);
}
