#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_nucleation_process.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace haero;

TEST_CASE("ternary_nuc_merik2007_test", "mam_nucleation_packed") {
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0  = 987659;
  const unsigned p1  =  12373;
  long unsigned seed =  54319;
  auto random = [&]() {
    seed =  (seed * p1)%p0;
    return double(seed)/p0;
  };
  for (int i=0; i<1000; i+=PackType::n) {
    PackType  t_p;
    PackType rh_p;
    PackType c2_p;
    PackType c3_p;
    for (int p=0; p<PackType::n; ++p) {
      const double  t=  235 +   60*random();  // range 235-295
      const double rh= 0.05 +   .9*random();  // range .05-.95
      const double c2= 5.e4 + 1.e8*random();  // range 5x10^4 - 10^9 
      const double c3=  0.1 +  999*random();  // range 0.1 - 1000
      t_p [p] = t;
      rh_p[p] = rh;
      c2_p[p] = c2;
      c3_p[p] = c3;
    }

    PackType j_log_p(0), ntot_p(0), nacid_p(0), namm_p(0), r_p(0);
    MAMNucleationProcess::ternary_nuc_merik2007(t_p, rh_p, c2_p, c3_p, j_log_p, ntot_p, nacid_p, namm_p, r_p);

    for (int p=0; p<PackType::n; ++p) {
      const double j_log_pck = j_log_p[p];
      const double ntot_pck  = ntot_p [p];
      const double nacid_pck = nacid_p[p];
      const double namm_pck  = namm_p [p];
      const double r_pck     = r_p    [p];

      const double  t = t_p [p];
      const double rh = rh_p[p];
      const double c2 = c2_p[p];
      const double c3 = c3_p[p];
      double j_log_dbl = 0;
      double ntot_dbl  = 0; 
      double nacid_dbl = 0;
      double namm_dbl  = 0; 
      double r_dbl     = 0;
      MAMNucleationProcess::ternary_nuc_merik2007(t, rh, c2, c3, j_log_dbl, ntot_dbl, nacid_dbl, namm_dbl, r_dbl);

      REQUIRE(j_log_dbl == j_log_pck);
      REQUIRE(ntot_dbl  == ntot_pck);
      REQUIRE(nacid_dbl == nacid_pck);
      REQUIRE(namm_dbl  == namm_pck);
      REQUIRE(r_dbl     == r_pck);
    }
  }
}

