#include <cstdio>

#include "catch2/catch.hpp"
#include "haero/processes/vehkamaki2002.hpp"

using namespace haero;

// These tests aren't really unit tests in the usual sense--they are "smoke
// tests" in that they exercise the parameterizations of Vehkamaki et al
// (2002) and generate plot files that can be used to compare the results of
// the parameterizations with the calculations from classical nucleation
// theory (e.g. given in Chapter 11 of the Third Edition of Seinfeld and
// Pandis, 2016).

// Compare the output of this test with Vehkamaki et al (2002), figure 7.
TEST_CASE("vehkamaki2002_figure_7") {
  FILE* fp = fopen("vehkamaki2002_figure_7.dat", "w");
  for (int i = 0; i < 12; ++i) {
    PackType temp(190.15 + 10 * i);
    for (int j = 0; j < 10; ++j) {
      PackType rel_hum(0.05 + j * 0.1);
      auto c_h2so4 = vehkamaki2002::h2so4_nucleation_threshold(temp, rel_hum);
      fprintf(fp, "%g %g %g\n", temp[0], rel_hum[0], c_h2so4[0]);
    }
  }
  fclose(fp);
}

// Compare the output of this test with Seinfeld and Pandis, figure 11.11.
// Alas, we can only test one temperature value within the region of validity!
TEST_CASE("vehkamaki2002_sp_figure_11.11") {
  FILE* fp = fopen("vehkamaki2002_sp_figure_11.11.dat", "w");
  PackType temp(273.0);
  for (int i = 0; i < 9; ++i) {
    PackType rel_hum(0.2 + i * 0.1);
    auto c_h2so4 = vehkamaki2002::h2so4_nucleation_threshold(temp, rel_hum);
    auto x_crit =
        vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, temp, rel_hum);
    fprintf(fp, "%g %g\n", rel_hum[0], x_crit[0]);
  }
  fclose(fp);
}

// Compare the output of this test with Vehkamaki et al (2002), figure 8
// (dotted line).
TEST_CASE("vehkamaki2002_figure_8_dotted") {
  FILE* fp = fopen("vehkamaki2002_figure_8_dotted.dat", "w");
  PackType temp(236.0), rel_hum(0.55);
  for (int i = 0; i < 100; ++i) {
    PackType c_h2so4(1e6 + 1e7 * i);
    auto x_crit =
        vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, temp, rel_hum);
    auto J = vehkamaki2002::nucleation_rate(c_h2so4, temp, rel_hum, x_crit);
    fprintf(fp, "%g %g\n", c_h2so4[0], J[0]);
  }
  fclose(fp);
}

// Compare the output of this test with Vehkamaki et al (2002), figure 8
// (dashed line).
TEST_CASE("vehkamaki2002_figure_8_dashed") {
  FILE* fp = fopen("vehkamaki2002_figure_8_dashed.dat", "w");
  PackType temp(236.0), rel_hum(0.55);
  for (int i = 0; i < 100; ++i) {
    PackType c_h2so4(1e6 + 1e7 * i);
    auto x_crit =
        vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, temp, rel_hum);
    auto n_tot =
        vehkamaki2002::num_critical_molecules(c_h2so4, temp, rel_hum, x_crit);
    fprintf(fp, "%g %g\n", c_h2so4[0], n_tot[0]);
  }
  fclose(fp);
}
