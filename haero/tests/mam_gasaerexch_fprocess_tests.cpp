#include <cmath>
#include <iomanip>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/conversions.hpp"
#include "haero/floating_point.hpp"
#include "haero/model.hpp"
#include "haero/processes/mam_gasaerexch_fprocess.hpp"
#include "haero/processes/mam_gasaerexch_process.hpp"
#include "mam_gasaerexch_test_bridge.hpp"

using namespace haero;

Model* get_model_for_unit_tests(ModalAerosolConfig& aero_config) {
  const int num_levels = 72;
  static Model* model(Model::ForUnitTests(aero_config, num_levels));
  return model;
}

TEST_CASE("mam_gasaerexch_1subarea_1gas_nonvolatile",
          "mam_gasaerexch_fprocess") {
  using fp_helper = FloatingPoint<Real>;
#ifdef NDEBUG
  const Real tolerance = 1.0e-08;
#else
  const Real tolerance = 1.0e-10;
#endif
  /// Test the compute_tendencies function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
//const unsigned p0 = 987659;
//const unsigned p1 = 12373;
//long unsigned seed = 54319;
//auto random = [&]() {
//  seed = (seed * p1) % p0;
//  return Real(seed) / p0;
//};


  const auto aero_species = create_mam4_aerosol_species();
  auto aero_config = create_mam4_modal_aerosol_config();
  const int num_levels = 72;
  const int num_modes = aero_config.aerosol_modes.size();
  const int num_gases = aero_config.gas_species.size();

  std::vector<int> num_aero_species(num_modes);
  std::vector<Mode> modes = create_mam4_modes();
  std::map<std::string, std::vector<std::string>> mode_species =
      create_mam4_mode_species();
  for (int m = 0; m < num_modes; ++m) {
    num_aero_species[m] = mode_species[modes[m].name()].size();
  }

  HostDiagnostics diagnostics(num_modes, num_aero_species, num_gases,
                              num_levels);

  get_model_for_unit_tests(aero_config);
  AerosolProcessType type = CloudBorneWetRemovalProcess;
  MAMGasAerosolExchangeProcess mam_gasaerexch_process(type, "gasaerexch Test",
                                              aero_config, diagnostics);
  init_bridge();

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
  Real uptkaer[] =  {6.6080559925483628E-004, 8.1580938179609411E-004, 3.6489706456131803E-004, 
                     4.5049020316212097E-004, 2.5940960604944770E-005, 3.2025877290055270E-005, 
                     7.1043539209444609E-005, 8.7708073098079762E-005, 0.0000000000000000, 0.0000000000000000};
  int mode_aging_optaa[] = {     0,          0,          0,          1,          0};
  int lptr2_soa_a_amode[]  = {  13,         20,         29, -999888777};
  
  mam_soaexch_1subarea_bridge(
    lund, 
    dt, 
    temp, 
    pmid, 
    aircon, 
    n_mode, 
    ntot_amode, 
    max_mode, 
    qgas_cur, 
    qgas_avg, 
    qaer_cur, 
    qnum_cur, 
    uptkaer, 
    mode_aging_optaa, 
    lptr2_soa_a_amode);

    REQUIRE(fp_helper::equiv(qgas_cur[0], 9.6445411630348715E-011, tolerance));
    REQUIRE(fp_helper::equiv(qgas_cur[1], 2.9533516044307414E-014, tolerance));
    REQUIRE(fp_helper::equiv(qgas_avg[0], 9.6499372481841032E-011, tolerance));
    REQUIRE(fp_helper::equiv(qgas_avg[1], 0.0000000000000000     , tolerance));

#if 0
  for (int i = 0; i < 100; ++i) {
    const Real deltat(random());
    const PackType temp(235 + 60 * random());  // range 235-295
    const PackType pmid(
        96325 + 10000 * random());    // pressure in Pascal, sea level=101,325
    const PackType aircon(random());  // air molar concentration (kmol/m3)
    const PackType zmid(500 + 10000 * random());  // layer midpoint height (m)
    const Real pblh(1000 + 1000 * random());      // boundary layer height (m)
    const PackType relhum(0.05 + .9 * random());  // range .05-.95
    const PackType uptkrate_h2so4(
        100 * random());  // h2so4 uptake rate to aerosol (1/s)
    const PackType del_h2so4_gasprod(random());
    const PackType del_h2so4_aeruptk(random());

    view_1d_pack_type qgas_cur("qgas_cur", 4);
    view_1d_pack_type qgas_avg("qgas_avg", 4);
    view_1d_pack_type qnum_cur("qnum_cur", 4);
    view_1d_pack_type qwtr_cur("qwrt_cur", 4);
    view_2d_pack_type qaer_cur("qaer_cur", 4, 4);
    for (int i = 0; i < 4; ++i) {
      qgas_cur(i) = random();
      qgas_avg(i) = random();
      qnum_cur(i) = random();
      qwtr_cur(i) = random();
      for (int j = 0; j < 4; ++j) {
        qaer_cur(i, j) = random();
      }
    }

    const Real adjust_factor_bin_tern_ratenucl = random();
    const Real adjust_factor_pbl_ratenucl = random();
    mam_gasaerexch_process.set_param("adjust_factor_bin_tern_ratenucl",
                                     adjust_factor_bin_tern_ratenucl);
    mam_gasaerexch_process.set_param("adjust_factor_pbl_ratenucl",
                                     adjust_factor_pbl_ratenucl);
    mam_gasaerexch_process.set_param("newnuc_adjust_factor_dnaitdt", 1.0);

    PackType dndt_ait(0);
    PackType dmdt_ait(0);
    PackType dso4dt_ait(0);
    PackType dnh4dt_ait(0);
    PackType nclusterdt(0);
    mam_gasaerexch_process.compute_tendencies(
        deltat, temp, pmid, aircon, zmid, pblh, relhum, uptkrate_h2so4,
        del_h2so4_gasprod, del_h2so4_aeruptk, qgas_cur, qgas_avg, qnum_cur,
        qaer_cur, qwtr_cur, dndt_ait, dmdt_ait, dso4dt_ait, dnh4dt_ait,
        nclusterdt);
    Real dndt_ait_f = 0;
    Real dmdt_ait_f = 0;
    Real dso4dt_ait_f = 0;
    Real dnh4dt_ait_f = 0;
    Real nclusterdt_f = 0;
    compute_tendencies_bridge(
        adjust_factor_bin_tern_ratenucl, adjust_factor_pbl_ratenucl, deltat,
        temp[0], pmid[0], aircon[0], zmid[0], pblh, relhum[0],
        uptkrate_h2so4[0], del_h2so4_gasprod[0], del_h2so4_aeruptk[0],
        &qgas_cur(0)[0], &qgas_avg(0)[0], &qnum_cur(0)[0], &qaer_cur(0, 0)[0],
        &qwtr_cur(0)[0], dndt_ait_f, dmdt_ait_f, dso4dt_ait_f, dnh4dt_ait_f,
        nclusterdt_f);

    REQUIRE((fp_helper::equiv(dndt_ait[0], dndt_ait_f, tolerance) ||
             fp_helper::rel(dndt_ait[0], dndt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(dmdt_ait[0], dmdt_ait_f, tolerance) ||
             fp_helper::rel(dmdt_ait[0], dmdt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(dso4dt_ait[0], dso4dt_ait_f, tolerance) ||
             fp_helper::rel(dso4dt_ait[0], dso4dt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(dnh4dt_ait[0], dnh4dt_ait_f, tolerance) ||
             fp_helper::rel(dnh4dt_ait[0], dnh4dt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(nclusterdt[0], nclusterdt_f, tolerance) ||
             fp_helper::rel(nclusterdt[0], nclusterdt_f, tolerance)));
  }
#endif
}
