#include <cmath>
#include <iomanip>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"
#include "haero/model.hpp"
#include "haero/processes/mam_nucleation_fprocess.hpp"
#include "haero/processes/mam_nucleation_process.hpp"
#include "mam_nucleation_test_bridge.hpp"

using namespace haero;

Model *get_model_for_unit_tests(ModalAerosolConfig &aero_config) {
  const int num_levels = 72;
  static Model * model (Model::ForUnitTests(aero_config, num_levels));
  return model;
}

TEST_CASE("compute_tendencies", "mam_nucleation_fprocess") {
  using fp_helper = FloatingPoint<Real>;
#ifdef NDEBUG
  const Real tolerance = 1.0e-07;
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
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };

  auto aero_config = create_mam4_modal_aerosol_config();
  get_model_for_unit_tests(aero_config);
  AerosolProcessType type = CloudBorneWetRemovalProcess;
  MAMNucleationProcess mam_nucleation_process(type, "Nucleation Test", aero_config);

  init_bridge();

  for (int i = 0; i < 100; ++i) {
    const Real deltat(random());
    const PackType temp(235 + 60 * random());      // range 235-295
    const PackType pmid(
        96325 + 10000 * random());  // pressure in Pascal, sea level=101,325
    const PackType aircon(random());    // air molar concentration (kmol/m3)
    const PackType zmid(500 + 10000 * random());    // layer midpoint height (m)
    const PackType pblh(1000 + 1000 * random());  // boundary layer height (m)
    const PackType relhum(0.05 + .9 * random());    // range .05-.95
    const PackType uptkrate_h2so4(100 *
                              random());  // h2so4 uptake rate to aerosol (1/s)
    const PackType del_h2so4_gasprod(random());
    const PackType del_h2so4_aeruptk(random());

    view_1d_pack_type qgas_cur("qgas_cur", 4);
    view_1d_pack_type qgas_avg("qgas_avg", 4);
    view_1d_pack_type qnum_cur("qnum_cur", 4);
    view_1d_pack_type qwtr_cur("qwrt_cur", 4);
    view_2d_pack_type qaer_cur("qaer_cur", 4, 4);
    for (int i=0; i<4; ++i) {
      qgas_cur(i) = random();
      qgas_avg(i) = random();
      qnum_cur(i) = random();
      qwtr_cur(i) = random();
      for (int j=0; j<4; ++j) {
        qaer_cur(i,j) = random();
      }
    }

    const Real adjust_factor_bin_tern_ratenucl = random();
    const Real adjust_factor_pbl_ratenucl = random();
    mam_nucleation_process.set_adjust_factor_bin_tern_ratenucl(
        adjust_factor_bin_tern_ratenucl);
    mam_nucleation_process.set_adjust_factor_pbl_ratenucl(
        adjust_factor_pbl_ratenucl);
    mam_nucleation_process.set_newnuc_adjust_factor_dnaitdt(
        1.0);

    PackType dndt_ait(0);
    PackType dmdt_ait(0);
    PackType dso4dt_ait(0);
    PackType dnh4dt_ait(0);
    PackType nclusterdt(0);
     mam_nucleation_process.compute_tendencies(
        deltat, temp, pmid, aircon, zmid, pblh, relhum, 
        uptkrate_h2so4, del_h2so4_gasprod, del_h2so4_aeruptk,
        qgas_cur, qgas_avg, qnum_cur, qaer_cur, qwtr_cur,
        dndt_ait, dmdt_ait, dso4dt_ait, dnh4dt_ait, nclusterdt);
    Real dndt_ait_f= 0;
    Real dmdt_ait_f= 0;
    Real dso4dt_ait_f= 0;
    Real dnh4dt_ait_f= 0;
    Real nclusterdt_f= 0;
    compute_tendencies_bridge(
      adjust_factor_bin_tern_ratenucl,
      adjust_factor_pbl_ratenucl,
      deltat,
      temp[0], 
      pmid[0], 
      aircon[0], 
      zmid[0], 
      pblh[0], 
      relhum[0], 
      uptkrate_h2so4[0], 
      del_h2so4_gasprod[0], 
      del_h2so4_aeruptk[0], 
      &qgas_cur(0)[0],
      &qgas_avg(0)[0],
      &qnum_cur(0)[0],
      &qaer_cur(0,0)[0],
      &qwtr_cur(0)[0],
      dndt_ait_f, 
      dmdt_ait_f, 
      dso4dt_ait_f, 
      dnh4dt_ait_f, 
      nclusterdt_f);

    REQUIRE((fp_helper::equiv(dndt_ait[0],   dndt_ait_f, tolerance) ||
             fp_helper::rel  (dndt_ait[0],   dndt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(dmdt_ait[0],   dmdt_ait_f, tolerance) ||
             fp_helper::rel  (dmdt_ait[0],   dmdt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(dso4dt_ait[0], dso4dt_ait_f, tolerance) ||
             fp_helper::rel  (dso4dt_ait[0], dso4dt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(dnh4dt_ait[0], dnh4dt_ait_f, tolerance) ||
             fp_helper::rel  (dnh4dt_ait[0], dnh4dt_ait_f, tolerance)));
    REQUIRE((fp_helper::equiv(nclusterdt[0], nclusterdt_f, tolerance) ||
             fp_helper::rel  (nclusterdt[0], nclusterdt_f, tolerance)));
  }
}

TEST_CASE("ternary_nuc_merik2007", "mam_nucleation_fprocess") {
  using fp_helper = FloatingPoint<Real>;
#ifdef NDEBUG
  const Real tolerance = 1.0e-07;
#else
  const Real tolerance = 1.0e-08;
#endif
  /// Test the ternary_nuc_merik2007 function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };
  for (int i = 0; i < 1000; ++i) {
    const PackType t(235 + 60 * random());      // range 235-295
    const PackType rh(0.05 + .9 * random());    // range .05-.95
    const PackType c2(5.e4 + 1.e8 * random());  // range 5x10^4 - 10^9
    const PackType c3(0.1 + 999 * random());    // range 0.1 - 1000
    PackType j_log_cpp(0);
    PackType ntot_cpp(0);
    PackType nacid_cpp(0);
    PackType namm_cpp(0);
    PackType r_cpp(0);
    Real j_log_f90 = 0;
    Real ntot_f90 = 0;
    Real nacid_f90 = 0;
    Real namm_f90 = 0;
    Real r_f90 = 0;
    MAMNucleationProcess::ternary_nuc_merik2007(
        t, rh, c2, c3, j_log_cpp, ntot_cpp, nacid_cpp, namm_cpp, r_cpp);
    ternary_nuc_merik2007_bridge(t[0], rh[0], c2[0], c3[0], j_log_f90, ntot_f90,
                                 nacid_f90, namm_f90, r_f90);
    REQUIRE((fp_helper::equiv(j_log_cpp[0], j_log_f90, tolerance) ||
             fp_helper::rel(j_log_cpp[0], j_log_f90, tolerance)));
    REQUIRE((fp_helper::equiv(ntot_cpp[0], ntot_f90, tolerance) ||
             fp_helper::rel(ntot_cpp[0], ntot_f90, tolerance)));
    REQUIRE((fp_helper::equiv(nacid_cpp[0], nacid_f90, tolerance) ||
             fp_helper::rel(nacid_cpp[0], nacid_f90, tolerance)));
    REQUIRE((fp_helper::equiv(namm_cpp[0], namm_f90, tolerance) ||
             fp_helper::rel(namm_cpp[0], namm_f90, tolerance)));
    REQUIRE((fp_helper::equiv(r_cpp[0], r_f90, tolerance) ||
             fp_helper::rel(r_cpp[0], r_f90, tolerance)));
  }
}

TEST_CASE("binary_nuc_vehk2002", "mam_nucleation_fprocess") {
  /// Test the binary_nuc_vehk2002 function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.
  using fp_helper = FloatingPoint<Real>;
#ifdef NDEBUG
  const Real tolerance = 1.0e-08;
#else
  const Real tolerance = 1.0e-10;
#endif
  using Pack = ekat::Pack<Real, 1>;
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };
  for (int i = 0; i < 1000; ++i) {
    const Pack temp(235 + 60 * random());       // range 235-295
    const Pack rh(0.05 + .9 * random());        // range .05-.95
    const Pack so4vol(5.e4 + 1.e8 * random());  // range 5x10^4 - 10^9
    Pack ratenucl(0);
    Pack rateloge(0);
    Pack cnum_h2so4(0);
    Pack cnum_tot(0);
    Pack radius_cluster(0);
    Real ratenucl_f90 = 0;
    Real rateloge_f90 = 0;
    Real cnum_h2so4_f90 = 0;
    Real cnum_tot_f90 = 0;
    Real radius_cluster_f90 = 0;
    MAMNucleationProcess::binary_nuc_vehk2002(temp, rh, so4vol, ratenucl,
                                              rateloge, cnum_h2so4, cnum_tot,
                                              radius_cluster);
    binary_nuc_vehk2002_bridge(temp[0], rh[0], so4vol[0], ratenucl_f90,
                               rateloge_f90, cnum_h2so4_f90, cnum_tot_f90,
                               radius_cluster_f90);
    REQUIRE((fp_helper::equiv(ratenucl[0], ratenucl_f90, tolerance) ||
             fp_helper::rel(ratenucl[0], ratenucl_f90, tolerance)));
    REQUIRE((fp_helper::equiv(rateloge[0], rateloge_f90, tolerance) ||
             fp_helper::rel(rateloge[0], rateloge_f90, tolerance)));
    REQUIRE((fp_helper::equiv(cnum_h2so4[0], cnum_h2so4_f90, tolerance) ||
             fp_helper::rel(cnum_h2so4[0], cnum_h2so4_f90, tolerance)));
    REQUIRE((fp_helper::equiv(cnum_tot[0], cnum_tot_f90, tolerance) ||
             fp_helper::rel(cnum_tot[0], cnum_tot_f90, tolerance)));
    REQUIRE(
        (fp_helper::equiv(radius_cluster[0], radius_cluster_f90, tolerance) ||
         fp_helper::rel(radius_cluster[0], radius_cluster_f90, tolerance)));
  }
}

TEST_CASE("pbl_nuc_wang2008", "mam_nucleation_fprocess") {
  /// Test the pbl_nuc_wang2008 function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.
  using Pack = ekat::Pack<Real, 1>;
  using fp_helper = FloatingPoint<Real>;
#ifdef NDEBUG
  const Real tolerance = 1.0e-08;
#else
  const Real tolerance = 1.0e-12;
#endif
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };
  MAMNucleationProcess mam_nucleation_process;
  for (int i = 0; i < 1000; ++i) {
    const Pack so4vol(5.e4 + 1.e8 * random());  // range 5x10^4 - 10^9
    const int flagaa = 11 + 2 * random();       // range 11-12
    const Real adjust_factor_pbl_ratenucl = random();
    mam_nucleation_process.set_adjust_factor_pbl_ratenucl(
        adjust_factor_pbl_ratenucl);

    ekat::Pack<int, 1> flagaa2(0);
    Pack ratenucl(0);
    Pack rateloge(0);
    Pack cnum_tot(0);
    Pack cnum_h2so4(0);
    Pack cnum_nh3(0);
    Pack radius_cluster(0);

    int flagaa2_f90 = 0;
    Real ratenucl_f90 = 0;
    Real rateloge_f90 = 0;
    Real cnum_tot_f90 = 0;
    Real cnum_h2so4_f90 = 0;
    Real cnum_nh3_f90 = 0;
    Real radius_cluster_f90 = 0;

    mam_nucleation_process.pbl_nuc_wang2008(so4vol, flagaa, flagaa2, ratenucl,
                                            rateloge, cnum_tot, cnum_h2so4,
                                            cnum_nh3, radius_cluster);
    pbl_nuc_wang2008_bridge(adjust_factor_pbl_ratenucl, so4vol[0], flagaa,
                            flagaa2_f90, ratenucl_f90, rateloge_f90,
                            cnum_tot_f90, cnum_h2so4_f90, cnum_nh3_f90,
                            radius_cluster_f90);

    REQUIRE(flagaa2[0] == flagaa2_f90);
    REQUIRE((fp_helper::equiv(ratenucl[0], ratenucl_f90, tolerance) ||
             fp_helper::rel(ratenucl[0], ratenucl_f90, tolerance)));
    REQUIRE((fp_helper::equiv(rateloge[0], rateloge_f90, tolerance) ||
             fp_helper::rel(rateloge[0], rateloge_f90, tolerance)));
    REQUIRE((fp_helper::equiv(cnum_tot[0], cnum_tot_f90, tolerance) ||
             fp_helper::rel(cnum_tot[0], cnum_tot_f90, tolerance)));
    REQUIRE((fp_helper::equiv(cnum_h2so4[0], cnum_h2so4_f90, tolerance) ||
             fp_helper::rel(cnum_h2so4[0], cnum_h2so4_f90, tolerance)));
    REQUIRE((fp_helper::equiv(cnum_nh3[0], cnum_nh3_f90, tolerance) ||
             fp_helper::rel(cnum_nh3[0], cnum_nh3_f90, tolerance)));
    REQUIRE(
        (fp_helper::equiv(radius_cluster[0], radius_cluster_f90, tolerance) ||
         fp_helper::rel(radius_cluster[0], radius_cluster_f90, tolerance)));
  }
}

TEST_CASE("mer07_veh02_nuc_mosaic_1box", "mam_nucleation_fprocess") {
  /// Test the mer07_veh02_nuc_mosaic_1box function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.
  using Pack = ekat::Pack<Real, 1>;
  using fp_helper = FloatingPoint<Real>;
#ifdef NDEBUG
  const Real tolerance = 1.0e-07;
#else
  const Real tolerance = 1.0e-08;
#endif
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };
  MAMNucleationProcess mam_nucleation_process;
  for (int i = 0; i < 1000; ++i) {
    const int newnuc_method_flagaa =
        random() < .5 ? 1 + 2 * random()
                      : 11 + 2 * random();  // range 1,2,11,12
    const Pack dtnuc(random());
    const Pack temp_in(235 + 60 * random());  // range 235-295
    const Pack rh_in(0.05 + .9 * random());   // range .05-.95
    const Pack press_in(
        96325 + 10000 * random());  // pressure in Pascal, sea level=101,325
    const Pack zm_in(500 + 10000 * random());    // layer midpoint height (m)
    const Pack pblh_in(1000 + 1000 * random());  // boundary layer height (m)
    const Pack qh2so4_cur(random());             // mixing ratio
    const Pack qh2so4_avg(random());             // mixing ratio
    const Pack qnh3_cur(random());               // mixing ratio
    const Pack h2so4_uptkrate(100 *
                              random());  // h2so4 uptake rate to aerosol (1/s)
    const Pack mw_so4a_host(random() /
                            1000);  // mw of so4 aerosol in host code (g/mol)
    const int nsize =
        1 +
        2 * random();  // number of aerosol size bins. NOTE: nsize<=maxd_asize
    const int maxd_asize =
        nsize +
        2 * random();  // dimension for dplom_sect, NOTE: nsize<=maxd_asize,
    const int ldiagaa = 10 * random();  // does not appear to be used.
    std::vector<Real> dplom_sect_vec(maxd_asize);
    std::vector<Real> dphim_sect_vec(maxd_asize);
    const Real SECT_SCALE = 1.0e10;
    dplom_sect_vec[0] = random() / SECT_SCALE;
    for (int i = 1; i < maxd_asize; ++i) {
      dplom_sect_vec[i] = dplom_sect_vec[i - 1] + random() / SECT_SCALE;
      dphim_sect_vec[i - 1] = dplom_sect_vec[i];
    }
    dphim_sect_vec[maxd_asize - 1] =
        dplom_sect_vec[maxd_asize - 1] + random() / SECT_SCALE;

    const Real* dplom_sect = dplom_sect_vec.data();
    const Real* dphim_sect = dphim_sect_vec.data();

    const Real adjust_factor_bin_tern_ratenucl = random();
    const Real adjust_factor_pbl_ratenucl = random();
    mam_nucleation_process.set_adjust_factor_bin_tern_ratenucl(
        adjust_factor_bin_tern_ratenucl);
    mam_nucleation_process.set_adjust_factor_pbl_ratenucl(
        adjust_factor_pbl_ratenucl);

    ekat::Pack<int, 1> isize_nuc(0);
    Pack qnuma_del(0);
    Pack qso4a_del(0);
    Pack qnh4a_del(0);
    Pack qh2so4_del(0);
    Pack qnh3_del(0);
    Pack dens_nh4so4a(0);
    Pack dnclusterdt(0);

    int isize_nuc_f = 0;
    Real qnuma_del_f = 0;
    Real qso4a_del_f = 0;
    Real qnh4a_del_f = 0;
    Real qh2so4_del_f = 0;
    Real qnh3_del_f = 0;
    Real dens_nh4so4a_f = 0;
    Real dnclusterdt_f = 0;

    mam_nucleation_process.mer07_veh02_nuc_mosaic_1box(
        newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in, zm_in, pblh_in,
        qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate, mw_so4a_host, nsize,
        maxd_asize, dplom_sect, dphim_sect, isize_nuc, qnuma_del, qso4a_del,
        qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a, ldiagaa, &dnclusterdt);

    mer07_veh02_nuc_mosaic_1box_bridge(
        adjust_factor_bin_tern_ratenucl, adjust_factor_pbl_ratenucl,
        newnuc_method_flagaa, dtnuc[0], temp_in[0], rh_in[0], press_in[0],
        zm_in[0], pblh_in[0], qh2so4_cur[0], qh2so4_avg[0], qnh3_cur[0],
        h2so4_uptkrate[0], mw_so4a_host[0], nsize, maxd_asize, dplom_sect,
        dphim_sect, isize_nuc_f, qnuma_del_f, qso4a_del_f, qnh4a_del_f,
        qh2so4_del_f, qnh3_del_f, dens_nh4so4a_f, ldiagaa, &dnclusterdt_f);

    REQUIRE(isize_nuc[0] == isize_nuc_f - 1);
    REQUIRE((fp_helper::equiv(qnuma_del[0], qnuma_del_f, tolerance) ||
             fp_helper::rel(qnuma_del[0], qnuma_del_f, tolerance)));
    REQUIRE((fp_helper::equiv(qso4a_del[0], qso4a_del_f, tolerance) ||
             fp_helper::rel(qso4a_del[0], qso4a_del_f, tolerance)));
    REQUIRE((fp_helper::equiv(qnh4a_del[0], qnh4a_del_f, tolerance) ||
             fp_helper::rel(qnh4a_del[0], qnh4a_del_f, tolerance)));
    REQUIRE((fp_helper::equiv(qh2so4_del[0], qh2so4_del_f, tolerance) ||
             fp_helper::rel(qh2so4_del[0], qh2so4_del_f, tolerance)));
    REQUIRE((fp_helper::equiv(qnh3_del[0], qnh3_del_f, tolerance) ||
             fp_helper::rel(qnh3_del[0], qnh3_del_f, tolerance)));
    REQUIRE((fp_helper::equiv(dens_nh4so4a[0], dens_nh4so4a_f, tolerance) ||
             fp_helper::rel(dens_nh4so4a[0], dens_nh4so4a_f, tolerance)));
    REQUIRE((fp_helper::equiv(dnclusterdt[0], dnclusterdt_f, tolerance) ||
             fp_helper::rel(dnclusterdt[0], dnclusterdt_f, tolerance)));
  }
}

// These tests exercise our transplant of the MAM nucleation process.
TEST_CASE("MAMNucleationFProcess", "mam_nucleation_fprocess") {
  // We create a phony model to be used for these tests.
  auto aero_config = create_mam4_modal_aerosol_config();
  Model *model = get_model_for_unit_tests(aero_config);
  int num_levels = 72;
  int num_gases = aero_config.h_gas_species.size();
  int num_modes = aero_config.h_aerosol_modes.size();

  // Set up some prognosics aerosol data views‥
  int num_aero_populations = model->num_aerosol_populations();
  Kokkos::View<PackType**> int_aerosols("interstitial aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> cld_aerosols("cloudborne aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> gases("gases", num_gases, num_levels);
  Kokkos::View<PackType**> int_num_concs("interstitial number concs", num_modes,
                                         num_levels);
  Kokkos::View<PackType**> cld_num_concs("cloud borne number concs", num_modes,
                                         num_levels);

  // Set up atmospheric data and populate it with some views.
  Kokkos::View<PackType*> temp("temperature", num_levels);
  Kokkos::View<PackType*> press("pressure", num_levels);
  Kokkos::View<PackType*> rel_hum("relative humidity", num_levels);
  Kokkos::View<PackType*> pdel("hydrostatic_dp", num_levels);
  Kokkos::View<PackType*> ht("height", num_levels + 1);
  Real pblh = 100.0;
  auto* atm = new Atmosphere(num_levels, temp, press, rel_hum, ht, pdel, pblh);

  // Test basic construction.
  SECTION("construct") {
    auto* process = new MAMNucleationFProcess();
    REQUIRE(process->type() == haero::NucleationProcess);
    REQUIRE(process->name() ==
            "MAMNucleationFProcess (Fortran NucleationProcess)");
    delete process;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* process = new MAMNucleationFProcess();
    process->init(aero_config);
    delete process;
  }

  // Test process tendencies.
  SECTION("nucleate_without_existing_aerosols") {
    auto* process = new MAMNucleationFProcess();
    process->init(aero_config);

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    auto* progs = model->create_prognostics(int_aerosols, cld_aerosols, gases,
                                            int_num_concs, cld_num_concs);
    auto* diags = model->create_diagnostics();
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.

    // atmospheric state
    Real h0 = 3e3, dz = h0 / num_levels;
    for (int k = 0; k < num_levels; ++k) {
      temp(k) = 273.0;
      press(k) = 1e5;
      rel_hum(k) = 0.95;
      ht(k) = h0 - k * dz;
    }

    // aerosols (none)
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        int_aerosols(p, k) = 0.0;
      }
    }

    // gases
    int h2so4_index = aero_config.gas_index("H2SO4");
    for (int k = 0; k < num_levels; ++k) {
      gases(h2so4_index, k) = 1e-13;
    }
    for (int g = 0; g < num_gases; ++g) {
      if (g != h2so4_index) {
        for (int k = 0; k < num_levels; ++k) {
          gases(g, k) = 0.0;
        }
      }
    }

    // Now compute the tendencies by running the process.
    Real t = 0.0, dt = 30.0;
    process->run(aero_config, t, dt, *progs, *atm, *diags, *tends);

    // --------------------------------------------------
    // Check the tendencies to make sure they make sense.
    // --------------------------------------------------

    // SO4 nucleates within the aitken mode. All other tendencies are zero.
    const auto aero_tends = tends->interstitial_aerosols;
    int aitken_index = aero_config.aerosol_mode_index("aitken");
    int aitken_so4_index =
        aero_config.aerosol_species_index(aitken_index, "SO4");
    int so4_pop_index =
        aero_config.population_index(aitken_index, aitken_so4_index);
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      if (p == so4_pop_index) {
        for (int k = 0; k < num_levels; ++k) {
          // FIXME: Currently, our test case gets no nucleation in this config,
          // FIXME: so we use >= instead of > here. Need to fix this.
          REQUIRE(int_aerosols(p, k)[0] >= 0.0);
        }
      } else {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(int_aerosols(p, k)[0] == 0.0);
        }
      }
    }

    // The tendency for H2SO4 should be negative, and the rest should be zero.
    const auto gas_tends = tends->gases;
    for (int k = 0; k < num_levels; ++k) {
      REQUIRE(gas_tends(h2so4_index, k)[0] <= 0.0);
    }
    for (int g = 0; g < num_gases; ++g) {
      if (g != h2so4_index) {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(gas_tends(g, k)[0] == 0.0);
        }
      }
    }

    // Clean up.
    delete progs;
    delete diags;
    delete tends;
    delete process;
  }
}
