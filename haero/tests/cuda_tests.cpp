#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "haero/floating_point.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace haero;

TEST_CASE("cuda_log_tolerance", "cuda_tests") {
  using fp_helper = FloatingPoint<double>;
  using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;
  using SolutionView = DeviceType::view_1d<double>; 
  static const int NUM_T=10000;
  const double tolerance = 2.5e-16;
  SolutionView x_val("input_values",NUM_T);
  SolutionView y_val("log_values",NUM_T);
  auto h_x_val = Kokkos::create_mirror_view(x_val);
  auto h_y_val = Kokkos::create_mirror_view(y_val);
  for (int i=0; i<NUM_T; ++i) {
    // uniform distribution [0,1]
    h_x_val(i) = double(i+1)/NUM_T;
  }
  Kokkos::deep_copy(x_val, h_x_val);
  Kokkos::parallel_for("Log_on_GPU", NUM_T,
    KOKKOS_LAMBDA(const int i) {
      y_val(i) = std::log(x_val(i));
    }
  );
  Kokkos::deep_copy(h_y_val, y_val);
  for (int i=0; i<NUM_T; ++i) {
    REQUIRE(fp_helper::equiv(std::log(h_x_val(i)), h_y_val(i), tolerance));
  }
}

KOKKOS_INLINE_FUNCTION
double ternary_nuc_merik2007(const double t,
                             const double rh,
                             const double c2,
                             const double c3)
{
  const double
      j_log = -12.861848898625231 + 4.905527742256349*c3 - 358.2337705052991*rh - \
        0.05463019231872484*c3*t + 4.8630382337426985*rh*t + \
        0.00020258394697064567*c3*(t*t) - 0.02175548069741675*rh*(t*t) - \
        2.502406532869512e-7*c3*(t*t*t) + 0.00003212869941055865*rh*(t*t*t) - \
        4.39129415725234e6/(log(c2)*log(c2)) + (56383.93843154586*t)/(log(c2)*log(c2)) - \
        (239.835990963361*(t*t))/(log(c2)*log(c2)) + \
        (0.33765136625580167*(t*t*t))/(log(c2)*log(c2)) - \
        (629.7882041830943*rh)/((c3*c3*c3)*log(c2)) + \
        (7.772806552631709*rh*t)/((c3*c3*c3)*log(c2)) - \
        (0.031974053936299256*rh*(t*t))/((c3*c3*c3)*log(c2)) + \
        (0.00004383764128775082*rh*(t*t*t))/((c3*c3*c3)*log(c2)) + \
        1200.472096232311*log(c2) - 17.37107890065621*t*log(c2) + \
        0.08170681335921742*(t*t)*log(c2) - 0.00012534476159729881*(t*t*t)*log(c2) - \
        14.833042158178936*(log(c2)*log(c2)) + 0.2932631303555295*t*(log(c2)*log(c2)) - \
        0.0016497524241142845*(t*t)*(log(c2)*log(c2)) + \
        2.844074805239367e-6*(t*t*t)*(log(c2)*log(c2)) - 231375.56676032578*log(c3) - \
        100.21645273730675*rh*log(c3) + 2919.2852552424706*t*log(c3) + \
        0.977886555834732*rh*t*log(c3) - 12.286497122264588*(t*t)*log(c3) - \
        0.0030511783284506377*rh*(t*t)*log(c3) + \
        0.017249301826661612*(t*t*t)*log(c3) + 2.967320346100855e-6*rh*(t*t*t)*log(c3) + \
        (2.360931724951942e6*log(c3))/log(c2) - \
        (29752.130254319443*t*log(c3))/log(c2) + \
        (125.04965118142027*(t*t)*log(c3))/log(c2) - \
        (0.1752996881934318*(t*t*t)*log(c3))/log(c2) + \
        5599.912337254629*log(c2)*log(c3) - 70.70896612937771*t*log(c2)*log(c3) + \
        0.2978801613269466*(t*t)*log(c2)*log(c3) - \
        0.00041866525019504*(t*t*t)*log(c2)*log(c3) + 75061.15281456841*(log(c3)*log(c3)) - \
        931.8802278173565*t*(log(c3)*log(c3)) + 3.863266220840964*(t*t)*(log(c3)*log(c3)) - \
        0.005349472062284983*(t*t*t)*(log(c3)*log(c3)) - \
        (732006.8180571689*(log(c3)*log(c3)))/log(c2) + \
        (9100.06398573816*t*(log(c3)*log(c3)))/log(c2) - \
        (37.771091915932004*(t*t)*(log(c3)*log(c3)))/log(c2) + \
        (0.05235455395566905*(t*t*t)*(log(c3)*log(c3)))/log(c2) - \
        1911.0303773001353*log(c2)*(log(c3)*log(c3)) + \
        23.6903969622286*t*log(c2)*(log(c3)*log(c3)) - \
        0.09807872005428583*(t*t)*log(c2)*(log(c3)*log(c3)) + \
        0.00013564560238552576*(t*t*t)*log(c2)*(log(c3)*log(c3)) - \
        3180.5610833308*(log(c3)*log(c3)*log(c3)) + 39.08268568672095*t*(log(c3)*log(c3)*log(c3)) - \
        0.16048521066690752*(t*t)*(log(c3)*log(c3)*log(c3)) + \
        0.00022031380023793877*(t*t*t)*(log(c3)*log(c3)*log(c3)) + \
        (40751.075322248245*(log(c3)*log(c3)*log(c3)))/log(c2) - \
        (501.66977622013934*t*(log(c3)*log(c3)*log(c3)))/log(c2) + \
        (2.063469732254135*(t*t)*(log(c3)*log(c3)*log(c3)))/log(c2) - \
        (0.002836873785758324*(t*t*t)*(log(c3)*log(c3)*log(c3)))/log(c2) + \
        2.792313345723013*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) - \
        0.03422552111802899*t*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) + \
        0.00014019195277521142*(t*t)*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) - \
        1.9201227328396297e-7*(t*t*t)*(log(c2)*log(c2))*(log(c3)*log(c3)*log(c3)) - \
        980.923146020468*log(rh) + 10.054155220444462*t*log(rh) - \
        0.03306644502023841*(t*t)*log(rh) + 0.000034274041225891804*(t*t*t)*log(rh) + \
        (16597.75554295064*log(rh))/log(c2) - \
        (175.2365504237746*t*log(rh))/log(c2) + \
        (0.6033215603167458*(t*t)*log(rh))/log(c2) - \
        (0.0006731787599587544*(t*t*t)*log(rh))/log(c2) - \
        89.38961120336789*log(c3)*log(rh) + 1.153344219304926*t*log(c3)*log(rh) - \
        0.004954549700267233*(t*t)*log(c3)*log(rh) + \
        7.096309866238719e-6*(t*t*t)*log(c3)*log(rh) + \
        3.1712136610383244*(log(c3)*log(c3)*log(c3))*log(rh) - \
        0.037822330602328806*t*(log(c3)*log(c3)*log(c3))*log(rh) + \
        0.0001500555743561457*(t*t)*(log(c3)*log(c3)*log(c3))*log(rh) - \
        1.9828365865570703e-7*(t*t*t)*(log(c3)*log(c3)*log(c3))*log(rh);
  return j_log;
}

TEST_CASE("cuda_ternary_nuc_merik2007_debgu_verses_opt", "cuda_tests") {
  /// This test shows the difference between debug and optimized CUDA code.
  /// The above function will produce a return value of  -10.3334643256330061
  /// when compiled optimized on the GPU. Unoptimized is -10.3334643243283608
  /// But the CPU version using gcc 7.2.0 produces       -10.3334643243283608
  /// in both cases.  So it appears to get consistent GPU vs CPU results,
  /// compile both unoptimized.

  using fp_helper = FloatingPoint<double>;
  using DeviceType = ekat::KokkosTypes<ekat::DefaultDevice>;
  using SolutionView = DeviceType::view_1d<double>; 
  const double tolerance = 1.5e-9;

  const double t = 277.2392748914352;
  const double rh= 0.4482234759162829;
  const double c2= 68835279.12974012;
  const double c3= 801.8841360226556;

  SolutionView sol("return_val",1);
  auto h_sol = Kokkos::create_mirror_view(sol);
  Kokkos::parallel_for("call_on_GPU", 1,
    KOKKOS_LAMBDA(const int) {
      const double y = ternary_nuc_merik2007(t, rh, c2, c3);
      sol(0) = y;
    }
  );
  Kokkos::deep_copy(h_sol, sol);
  const double y = ternary_nuc_merik2007(t, rh, c2, c3);
  REQUIRE(fp_helper::equiv(y, h_sol(0), tolerance));
}

