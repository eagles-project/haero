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

        // If we're above the onset temperature, J = 0.
        auto onset_temp =
            merikanto2007::onset_temperature(rel_hum, c_h2so4, xi_nh3);
        PackType J(0);
        if (temp[0] < onset_temp[0]) {
          auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                          c_h2so4, xi_nh3);
          J = exp(log_J);
        }

        fprintf(fp, "%g %g %g %g\n", temp[0], xi_nh3[0], c_h2so4[0], J[0]);
      }
    }
  }
  fclose(fp);
}

// Compare the output of this test with Merikanto et al (2007), figure 3.
TEST_CASE("merikanto2007_figure_3") {
  PackType temps[2] = {235.15, 273.15};
  PackType c_h2so4s[2] = {1e6, 1e9};
  PackType xi_nh3s[2][2] = {// T = 235.15
                            {0.1, 100},
                            // T = 273.15
                            {10, 1000}};
  FILE* fp = fopen("merikanto2007_figure_3.dat", "w");
  for (int i = 0; i < 2; ++i) {
    PackType temp = temps[i];
    PackType c_h2so4 = c_h2so4s[i];
    for (int j = 0; j < 2; ++j) {
      PackType xi_nh3 = xi_nh3s[i][j];
      for (int k = 0; k < 10; ++k) {
        PackType rel_hum(0.05 + k * 0.1);

        // If we're above the onset temperature, J = 0.
        auto onset_temp =
            merikanto2007::onset_temperature(rel_hum, c_h2so4, xi_nh3);
        PackType J(0);
        if (temp[0] < onset_temp[0]) {
          auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                          c_h2so4, xi_nh3);
          J = exp(log_J);
        }

        fprintf(fp, "%g %g %g %g\n", temp[0], xi_nh3[0], rel_hum[0], J[0]);
      }
    }
  }
  fclose(fp);
}

TEST_CASE("merikanto2007_log_h2so4_nucleation_threshold") {
  // Compute the nucleation threshold for a given temperature and NH3 mix ratio.
  PackType temp(240), xi(10);
  auto thresh_c_h2so4 =
      exp(merikanto2007::log_h2so4_nucleation_threshold(temp, xi));

  // Compute the onset temperature to make sure we're below it.
  PackType rel_hum(0.5);
  auto onset_temp =
      merikanto2007::onset_temperature(rel_hum, thresh_c_h2so4, xi);
  REQUIRE(temp[0] < onset_temp[0]);

  // Compute the log of the nucleation rate and make sure it's zero or less.
  auto log_J =
      merikanto2007::log_nucleation_rate(temp, rel_hum, thresh_c_h2so4, xi);
  REQUIRE(log_J[0] <= 0);
}
