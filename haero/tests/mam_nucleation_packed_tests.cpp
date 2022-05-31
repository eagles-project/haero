#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_nucleation_process.hpp"

using namespace haero;

TEST_CASE("ternary_nuc_merik2007_test", "mam_nucleation_packed") {
  using fp_helper = FloatingPoint<double>;
#ifdef NDEBUG
  const double tolerance = 1.0e-08;
#else
  const double tolerance = 1.0e-20;
#endif
  /// Test the ternary_nuc_merik2007 function directly by calling
  /// the function with both a pack size of 1 and the default pack
  /// size and comparing the result.  The pack size of 1 has been
  /// tested against the Fortran version elsewhere so it assumed
  /// to be correct.

  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  using Pack = ekat::Pack<Real, 1>;
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return double(seed) / p0;
  };
  for (int i = 0; i < 1000; i += PackType::n) {
    PackType t_p;
    PackType rh_p;
    PackType c2_p;
    PackType c3_p;
    for (int p = 0; p < PackType::n; ++p) {
      t_p[p] = 235 + 60 * random();      // range 235-295
      rh_p[p] = 0.05 + .9 * random();    // range .05-.95
      c2_p[p] = 5.e4 + 1.e8 * random();  // range 5x10^4 - 10^9
      c3_p[p] = 0.1 + 999 * random();    // range 0.1 - 1000
    }

    PackType j_log_p(0), ntot_p(0), nacid_p(0), namm_p(0), r_p(0);
    MAMNucleationProcess::ternary_nuc_merik2007(t_p, rh_p, c2_p, c3_p, j_log_p,
                                                ntot_p, nacid_p, namm_p, r_p);

    for (int p = 0; p < PackType::n; ++p) {
      const Pack t(t_p[p]);
      const Pack rh(rh_p[p]);
      const Pack c2(c2_p[p]);
      const Pack c3(c3_p[p]);
      Pack j_log_dbl(0);
      Pack ntot_dbl(0);
      Pack nacid_dbl(0);
      Pack namm_dbl(0);
      Pack r_dbl(0);
      MAMNucleationProcess::ternary_nuc_merik2007(
          t, rh, c2, c3, j_log_dbl, ntot_dbl, nacid_dbl, namm_dbl, r_dbl);

      REQUIRE((fp_helper::equiv(j_log_dbl[0], j_log_p[p], tolerance) ||
               fp_helper::rel(j_log_dbl[0], j_log_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(ntot_dbl[0], ntot_p[p], tolerance) ||
               fp_helper::rel(ntot_dbl[0], ntot_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(nacid_dbl[0], nacid_p[p], tolerance) ||
               fp_helper::rel(nacid_dbl[0], nacid_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(namm_dbl[0], namm_p[p], tolerance) ||
               fp_helper::rel(namm_dbl[0], namm_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(r_dbl[0], r_p[p], tolerance) ||
               fp_helper::rel(r_dbl[0], r_p[p], tolerance)));
    }
  }
}

TEST_CASE("binary_nuc_vehk2002", "mam_nucleation_packed") {
  using fp_helper = FloatingPoint<double>;
#ifdef NDEBUG
  const double tolerance = 1.0e-12;
#else
  const double tolerance = 1.0e-20;
#endif
  /// Test the binary_nuc_vehk2002 function directly by calling
  /// the function with both a pack size of 1 and the default pack
  /// size and comparing the result.  The pack size of 1 has been
  /// tested against the Fortran version elsewhere so it assumed
  /// to be correct.
  using Pack = ekat::Pack<Real, 1>;
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return double(seed) / p0;
  };
  for (int i = 0; i < 1000; i += PackType::n) {
    PackType temp_p;
    PackType rh_p;
    PackType so4vol_p;
    for (int p = 0; p < PackType::n; ++p) {
      temp_p[p] = 235 + 60 * random();       // range 235-295
      rh_p[p] = 0.05 + .9 * random();        // range .05-.95
      so4vol_p[p] = 5.e4 + 1.e8 * random();  // range 5x10^4 - 10^9
    }
    PackType ratenucl_p(0), rateloge_p(0), cnum_h2so4_p(0), cnum_tot_p(0),
        radius_cluster_p(0);
    MAMNucleationProcess::binary_nuc_vehk2002(
        temp_p, rh_p, so4vol_p, ratenucl_p, rateloge_p, cnum_h2so4_p,
        cnum_tot_p, radius_cluster_p);
    for (int p = 0; p < PackType::n; ++p) {
      const Pack temp(temp_p[p]);
      const Pack rh(rh_p[p]);
      const Pack so4vol(so4vol_p[p]);
      Pack ratenucl(0);
      Pack rateloge(0);
      Pack cnum_h2so4(0);
      Pack cnum_tot(0);
      Pack radius_cluster(0);
      MAMNucleationProcess::binary_nuc_vehk2002(temp, rh, so4vol, ratenucl,
                                                rateloge, cnum_h2so4, cnum_tot,
                                                radius_cluster);
      REQUIRE((fp_helper::equiv(ratenucl[0], ratenucl_p[p], tolerance) ||
               fp_helper::rel(ratenucl[0], ratenucl_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(rateloge[0], rateloge_p[p], tolerance) ||
               fp_helper::rel(rateloge[0], rateloge_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(cnum_h2so4[0], cnum_h2so4_p[p], tolerance) ||
               fp_helper::rel(cnum_h2so4[0], cnum_h2so4_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(cnum_tot[0], cnum_tot_p[p], tolerance) ||
               fp_helper::rel(cnum_tot[0], cnum_tot_p[p], tolerance)));
      REQUIRE((
          fp_helper::equiv(radius_cluster[0], radius_cluster_p[p], tolerance) ||
          fp_helper::rel(radius_cluster[0], radius_cluster_p[p], tolerance)));
    }
  }
}

TEST_CASE("pbl_nuc_wang2008", "mam_nucleation_packed") {
  using fp_helper = FloatingPoint<double>;
#ifdef NDEBUG
  const double tolerance = 1.0e-12;
#else
  const double tolerance = 1.0e-20;
#endif
  /// Test the pbl_nuc_wang2008 function directly by calling
  /// the function with both a pack size of 1 and the default pack
  /// size and comparing the result.  The pack size of 1 has been
  /// tested against the Fortran version elsewhere so it assumed
  /// to be correct.
  using Pack = ekat::Pack<Real, 1>;
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return double(seed) / p0;
  };
  MAMNucleationProcess mam_nucleation_process;
  for (int i = 0; i < 1000; i += PackType::n) {
    PackType so4vol_p;
    const int flagaa = 11 + 2 * random();  // range 11-12
    const Real adjust_factor_pbl_ratenucl = random();
    for (int p = 0; p < PackType::n; ++p)
      so4vol_p[p] = 5.e4 + 1.e8 * random();  // range 5x10^4 - 10^9
    mam_nucleation_process.set_param("adjust_factor_pbl_ratenucl",
                                     adjust_factor_pbl_ratenucl);

    ekat::Pack<int, PackType::n> flagaa2_p(0);
    PackType ratenucl_p(0), rateloge_p(0), cnum_tot_p(0), cnum_h2so4_p(0),
        cnum_nh3_p(0), radius_cluster_p(0);
    mam_nucleation_process.pbl_nuc_wang2008(
        so4vol_p, flagaa, flagaa2_p, ratenucl_p, rateloge_p, cnum_tot_p,
        cnum_h2so4_p, cnum_nh3_p, radius_cluster_p);

    for (int p = 0; p < PackType::n; ++p) {
      const Pack so4vol(so4vol_p[p]);
      ekat::Pack<int, 1> flagaa2(0);
      Pack ratenucl(0), rateloge(0), cnum_tot(0), cnum_h2so4(0), cnum_nh3(0),
          radius_cluster(0);
      mam_nucleation_process.pbl_nuc_wang2008(so4vol, flagaa, flagaa2, ratenucl,
                                              rateloge, cnum_tot, cnum_h2so4,
                                              cnum_nh3, radius_cluster);
      REQUIRE((fp_helper::equiv(ratenucl[0], ratenucl_p[p], tolerance) ||
               fp_helper::rel(ratenucl[0], ratenucl_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(rateloge[0], rateloge_p[p], tolerance) ||
               fp_helper::rel(rateloge[0], rateloge_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(cnum_tot[0], cnum_tot_p[p], tolerance) ||
               fp_helper::rel(cnum_tot[0], cnum_tot_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(cnum_h2so4[0], cnum_h2so4_p[p], tolerance) ||
               fp_helper::rel(cnum_h2so4[0], cnum_h2so4_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(cnum_nh3[0], cnum_nh3_p[p], tolerance) ||
               fp_helper::rel(cnum_nh3[0], cnum_nh3_p[p], tolerance)));
      REQUIRE((
          fp_helper::equiv(radius_cluster[0], radius_cluster_p[p], tolerance) ||
          fp_helper::rel(radius_cluster[0], radius_cluster_p[p], tolerance)));
      REQUIRE(flagaa2[0] == flagaa2_p[p]);
    }
  }
}

TEST_CASE("mer07_veh02_nuc_mosaic_1box", "mam_nucleation_process") {
  using fp_helper = FloatingPoint<double>;
#ifdef NDEBUG
  const double tolerance = 1.0e-08;
#else
  const double tolerance = 1.0e-20;
#endif
  /// Test the mer07_veh02_nuc_mosaic_1bo function directly by calling
  /// the function with both a pack size of 1 and the default pack
  /// size and comparing the result.  The pack size of 1 has been
  /// tested against the Fortran version elsewhere so it assumed
  /// to be correct.
  using Pack = ekat::Pack<Real, 1>;
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0 = 987659;
  const unsigned p1 = 12373;
  long unsigned seed = 54319;
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return double(seed) / p0;
  };
  MAMNucleationProcess mam_nucleation_process;
  for (int i = 0; i < 1000; ++i) {
    const int newnuc_method_flagaa =
        random() < .5 ? 1 + 2 * random()
                      : 11 + 2 * random();  // range 1,2,11,12
    PackType dtnuc_p;
    PackType temp_in_p;
    PackType rh_in_p;
    PackType press_in_p;
    PackType zm_in_p;
    PackType pblh_in_p;
    PackType qh2so4_cur_p;
    PackType qh2so4_avg_p;
    PackType qnh3_cur_p;
    PackType h2so4_uptkrate_p;
    PackType mw_so4a_host_p;
    for (int p = 0; p < PackType::n; ++p) {
      dtnuc_p[p] = random();
      temp_in_p[p] = 235 + 60 * random();  // range 235-295
      rh_in_p[p] = 0.05 + .9 * random();   // range .05-.95
      press_in_p[p] =
          96325 + 10000 * random();  // pressure in Pascal, sea level=101,325
      zm_in_p[p] = 500 + 10000 * random();    // layer midpoint height (m)
      pblh_in_p[p] = 1000 + 1000 * random();  // boundary layer height (m)
      qh2so4_cur_p[p] = random();             // mixing ratio
      qh2so4_avg_p[p] = random();             // mixing ratio
      qnh3_cur_p[p] = random();               // mixing ratio
      h2so4_uptkrate_p[p] =
          100 * random();  // h2so4 uptake rate to aerosol (1/s)
      mw_so4a_host_p[p] =
          random() / 1000;  // mw of so4 aerosol in host code (g/mol)
    }
    const int nsize =
        1 +
        2 * random();  // number of aerosol size bins. NOTE: nsize<=maxd_asize
    const int maxd_asize =
        nsize +
        2 * random();  // dimension for dplom_sect, NOTE: nsize<=maxd_asize,
    const int ldiagaa = 10 * random();  // does not appear to be used.
    std::vector<Real> dplom_sect(maxd_asize);
    std::vector<Real> dphim_sect(maxd_asize);
    const Real SECT_SCALE = 1.0e10;
    dplom_sect[0] = random() / SECT_SCALE;
    for (int i = 1; i < maxd_asize; ++i) {
      dplom_sect[i] = dplom_sect[i - 1] + random() / SECT_SCALE;
      dphim_sect[i - 1] = dplom_sect[i];
    }
    dphim_sect[maxd_asize - 1] =
        dplom_sect[maxd_asize - 1] + random() / SECT_SCALE;

    const Real adjust_factor_bin_tern_ratenucl = random();
    const Real adjust_factor_pbl_ratenucl = random();
    mam_nucleation_process.set_param("adjust_factor_bin_tern_ratenucl",
                                     adjust_factor_bin_tern_ratenucl);
    mam_nucleation_process.set_param("adjust_factor_pbl_ratenucl",
                                     adjust_factor_pbl_ratenucl);

    ekat::Pack<int, PackType::n> isize_nuc_p(0);
    PackType qnuma_del_p(0);
    PackType qso4a_del_p(0);
    PackType qnh4a_del_p(0);
    PackType qh2so4_del_p(0);
    PackType qnh3_del_p(0);
    PackType dens_nh4so4a_p(0);
    PackType dnclusterdt_p(0);
    mam_nucleation_process.mer07_veh02_nuc_mosaic_1box(
        newnuc_method_flagaa, dtnuc_p, temp_in_p, rh_in_p, press_in_p, zm_in_p,
        pblh_in_p, qh2so4_cur_p, qh2so4_avg_p, qnh3_cur_p, h2so4_uptkrate_p,
        mw_so4a_host_p, nsize, maxd_asize, dplom_sect.data(), dphim_sect.data(),
        isize_nuc_p, qnuma_del_p, qso4a_del_p, qnh4a_del_p, qh2so4_del_p,
        qnh3_del_p, dens_nh4so4a_p, ldiagaa, &dnclusterdt_p);
    for (int p = 0; p < PackType::n; ++p) {
      const Pack dtnuc(dtnuc_p[p]);
      const Pack temp_in(temp_in_p[p]);
      const Pack rh_in(rh_in_p[p]);
      const Pack press_in(press_in_p[p]);
      const Pack zm_in(zm_in_p[p]);
      const Pack pblh_in(pblh_in_p[p]);
      const Pack qh2so4_cur(qh2so4_cur_p[p]);
      const Pack qh2so4_avg(qh2so4_avg_p[p]);
      const Pack qnh3_cur(qnh3_cur_p[p]);
      const Pack h2so4_uptkrate(h2so4_uptkrate_p[p]);
      const Pack mw_so4a_host(mw_so4a_host_p[p]);
      Pack qnuma_del(0);
      Pack qso4a_del(0);
      Pack qnh4a_del(0);
      Pack qh2so4_del(0);
      Pack qnh3_del(0);
      Pack dens_nh4so4a(0);
      Pack dnclusterdt(0);
      ekat::Pack<int, 1> isize_nuc(0);

      mam_nucleation_process.mer07_veh02_nuc_mosaic_1box(
          newnuc_method_flagaa, dtnuc, temp_in, rh_in, press_in, zm_in, pblh_in,
          qh2so4_cur, qh2so4_avg, qnh3_cur, h2so4_uptkrate, mw_so4a_host, nsize,
          maxd_asize, dplom_sect.data(), dphim_sect.data(), isize_nuc,
          qnuma_del, qso4a_del, qnh4a_del, qh2so4_del, qnh3_del, dens_nh4so4a,
          ldiagaa, &dnclusterdt);
      REQUIRE((fp_helper::equiv(qnuma_del[0], qnuma_del_p[p], tolerance) ||
               fp_helper::rel(qnuma_del[0], qnuma_del_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(qso4a_del[0], qso4a_del_p[p], tolerance) ||
               fp_helper::rel(qso4a_del[0], qso4a_del_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(qnh4a_del[0], qnh4a_del_p[p], tolerance) ||
               fp_helper::rel(qnh4a_del[0], qnh4a_del_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(qh2so4_del[0], qh2so4_del_p[p], tolerance) ||
               fp_helper::rel(qh2so4_del[0], qh2so4_del_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(qnh3_del[0], qnh3_del_p[p], tolerance) ||
               fp_helper::rel(qnh3_del[0], qnh3_del_p[p], tolerance)));
      REQUIRE(
          (fp_helper::equiv(dens_nh4so4a[0], dens_nh4so4a_p[p], tolerance) ||
           fp_helper::rel(dens_nh4so4a[0], dens_nh4so4a_p[p], tolerance)));
      REQUIRE((fp_helper::equiv(dnclusterdt[0], dnclusterdt_p[p], tolerance) ||
               fp_helper::rel(dnclusterdt[0], dnclusterdt_p[p], tolerance)));
      REQUIRE(isize_nuc[0] == isize_nuc_p[p]);
    }
  }
}
