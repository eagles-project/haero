#include <cmath>
#include <iomanip>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/conversions.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_gasaerexch_fprocess.hpp"
#include "haero/processes/mam_gasaerexch_process.hpp"
#include "mam_gasaerexch_test_bridge.hpp"

using namespace haero;

TEST_CASE("mam_gasaerexch_1subarea_1gas_nonvolatile",
          "mam_gasaerexch_fprocess") {
  using fp_helper = FloatingPoint<Real>;
#ifdef NDEBUG
  const Real tolerance = 1.0e-08;
#else
  const Real tolerance = 1.0e-10;
#endif

  auto aero_config = ModalAerosolConfig::create_mam4_config();
  const int num_levels = 72;

  HostDiagnostics diagnostics(aero_config, num_levels);
  MAMGasAerosolExchangeFProcess mam_gasaerexch_fprocess;
  mam_gasaerexch_fprocess.init(aero_config);  // this sets up Fortran stuff
  // MAMGasAerosolExchangeProcess mam_gasaerexch_process;("gasaerexch Test",
  //                                                        aero_config,
  //                                                        diagnostics);

  const int max_gas = 2;
  const int max_aer = 7;
  const int iaer_pom = 3;
  const int nsoa = 1;
  const int npoa = 1;
  const int npca = 4;
  const Real pstd = 101325.0;
  const Real r_universal = 314.467591;
  const int lund = 93;
  const Real dt = 1.0000000000000000;
  const Real temp = 273.00000000000000;
  const Real pmid = 100000.00000000000;
  const Real aircon = 4.4055781358372036E-002;
  const int n_mode = 4;
  const int ntot_amode = 4;
  const int max_mode = 5;
  Real qgas_cur[] = {9.6553333333333350E-011, 2.9533516044307414E-014};
  Real qgas_avg[] = {0.0000000000000000, 0.0000000000000000};
  Real qaer_cur[] = {2253176148.8415728, 22531761488.415726, 2253176.1488415725,
                     4506352297.6831455, 0.0000000000000000};
  Real qnum_cur[] = {2253176148.8415728, 22531761488.415726, 2253176.1488415725,
                     4506352297.6831455, 0.0000000000000000};
  Real uptkaer[] = {6.6080559925483628E-004, 8.1580938179609411E-004,
                    3.6489706456131803E-004, 4.5049020316212097E-004,
                    2.5940960604944770E-005, 3.2025877290055270E-005,
                    7.1043539209444609E-005, 8.7708073098079762E-005,
                    0.0000000000000000,      0.0000000000000000};
  int mode_aging_optaa[] = {0, 0, 0, 1, 0};
  int lptr2_soa_a_amode[] = {13, 20, 29, -999888777};

  mam_soaexch_1subarea_bridge(
      max_gas, max_aer, iaer_pom, nsoa, npoa, npca, pstd, r_universal, lund, dt,
      temp, pmid, aircon, n_mode, ntot_amode, max_mode, qgas_cur, qgas_avg,
      qaer_cur, qnum_cur, uptkaer, mode_aging_optaa, lptr2_soa_a_amode);
  REQUIRE(fp_helper::equiv(qgas_cur[0], 0.0, tolerance));
  REQUIRE(fp_helper::equiv(qgas_cur[1], 2.95335e-14, tolerance));
  REQUIRE(fp_helper::equiv(qgas_avg[0], 4.82767e-11, tolerance));
  REQUIRE(fp_helper::equiv(qgas_avg[1], 0.0, tolerance));
}
