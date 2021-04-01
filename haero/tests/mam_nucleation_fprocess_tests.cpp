#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_nucleation_fprocess.hpp"
#include "haero/processes/mam_nucleation_process.hpp"
#include "mam_nucleation_test_bridge.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace haero;

using namespace haero;

TEST_CASE("ternary_nuc_merik2007", "mam_nucleation_fprocess") {
  /// Test the ternary_nuc_merik2007 function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0  = 987659;
  const unsigned p1  =  12373;
  long unsigned seed =  54319;
  auto random = [&]() {
    seed =  (seed * p1)%p0;
    return double(seed)/p0;
  };
  for (int i=0; i<1000; ++i) {
    const double  t=  235 +   60*random();  // range 235-295
    const double rh= 0.05 +   .9*random();  // range .05-.95
    const double c2= 5.e4 + 1.e8*random();  // range 5x10^4 - 10^9
    const double c3=  0.1 +  999*random();  // range 0.1 - 1000
    double j_log_cpp = 0;
    double ntot_cpp = 0;
    double nacid_cpp= 0;
    double namm_cpp = 0;
    double r_cpp    = 0;
    double j_log_f90 = 0;
    double ntot_f90 = 0;
    double nacid_f90= 0;
    double namm_f90 = 0;
    double r_f90    = 0;
    MAMNucleationProcess::ternary_nuc_merik2007(t, rh, c2, c3, j_log_cpp, ntot_cpp, nacid_cpp, namm_cpp, r_cpp);
    ternary_nuc_merik2007_bridge(t, rh, c2, c3, j_log_f90, ntot_f90, nacid_f90, namm_f90, r_f90);
    REQUIRE(j_log_cpp == j_log_f90);
    REQUIRE(ntot_cpp  == ntot_f90);
    REQUIRE(nacid_cpp == nacid_f90);
    REQUIRE(namm_cpp  == namm_f90);
    REQUIRE(r_cpp     == r_f90);
  }
}

TEST_CASE("binary_nuc_vehk2002", "mam_nucleation_fprocess") {
  /// Test the binary_nuc_vehk2002 function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  using fp_helper = FloatingPoint<double>;
  const double abs_tol = 2.0e-12;
  const double rel_tol = 2.0e-12;
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0  = 987659;
  const unsigned p1  =  12373;
  long unsigned seed =  54319;
  auto random = [&]() {
    seed =  (seed * p1)%p0;
    return double(seed)/p0;
  };
  for (int i=0; i<1000; ++i) {
    const double temp   =  235 +   60*random();  // range 235-295
    const double rh     = 0.05 +   .9*random();  // range .05-.95
    const double so4vol = 5.e4 + 1.e8*random();  // range 5x10^4 - 10^9
    double ratenucl           = 0;
    double rateloge           = 0;
    double cnum_h2so4         = 0;
    double cnum_tot           = 0;
    double radius_cluster     = 0;
    double ratenucl_f90       = 0;
    double rateloge_f90       = 0;
    double cnum_h2so4_f90     = 0;
    double cnum_tot_f90       = 0;
    double radius_cluster_f90 = 0;
    MAMNucleationProcess::binary_nuc_vehk2002(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster);
    binary_nuc_vehk2002_bridge(temp, rh, so4vol, ratenucl_f90, rateloge_f90, cnum_h2so4_f90, cnum_tot_f90, radius_cluster_f90);
    REQUIRE( (fp_helper::equiv(ratenucl      , ratenucl_f90, abs_tol)       || fp_helper::rel(ratenucl   , ratenucl_f90, rel_tol)) );
    REQUIRE( (fp_helper::equiv(rateloge      ,  rateloge_f90, abs_tol)      || fp_helper::rel(rateloge   ,  rateloge_f90, rel_tol)) );
    REQUIRE( (fp_helper::equiv(cnum_h2so4    , cnum_h2so4_f90, abs_tol)     || fp_helper::rel(cnum_h2so4 , cnum_h2so4_f90, rel_tol)) );
    REQUIRE( (fp_helper::equiv(cnum_tot      , cnum_tot_f90, abs_tol)       || fp_helper::rel(cnum_tot   , cnum_tot_f90, rel_tol)) );
    REQUIRE( (fp_helper::equiv(radius_cluster, radius_cluster_f90, abs_tol) || fp_helper::rel(radius_cluster   , radius_cluster_f90, rel_tol)) );
  }
}

TEST_CASE("pbl_nuc_wang2008", "mam_nucleation_fprocess") {
  /// Test the pbl_nuc_wang2008 function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  using fp_helper = FloatingPoint<double>;
  const double abs_tol = 2.0e-12;
  const double rel_tol = 2.0e-12;
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0  = 987659;
  const unsigned p1  =  12373;
  long unsigned seed =  54319;
  auto random = [&]() {
    seed =  (seed * p1)%p0;
    return double(seed)/p0;
  };
  MAMNucleationProcess mam_nucleation_process;
  for (int i=0; i<1000; ++i) {
    const double so4vol = 5.e4 + 1.e8*random();  // range 5x10^4 - 10^9
    const int flagaa = 11 + 2*random();  // range 11-12   
    const double adjust_factor_pbl_ratenucl = random();
    mam_nucleation_process.set_adjust_factor_pbl_ratenucl(adjust_factor_pbl_ratenucl);

    int    flagaa2 = 0;
    double ratenucl              = 0;
    double rateloge              = 0;
    double cnum_tot              = 0;
    double cnum_h2so4            = 0;
    double cnum_nh3              = 0;
    double radius_cluster        = 0;

    int    flagaa2_f90 = 0;
    double ratenucl_f90              = 0;
    double rateloge_f90              = 0;
    double cnum_tot_f90              = 0;
    double cnum_h2so4_f90            = 0;
    double cnum_nh3_f90              = 0;
    double radius_cluster_f90        = 0;

    mam_nucleation_process.pbl_nuc_wang2008(so4vol, flagaa, flagaa2, ratenucl, rateloge, 
      cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster);
    pbl_nuc_wang2008_bridge(adjust_factor_pbl_ratenucl, so4vol, flagaa, flagaa2_f90, ratenucl_f90, rateloge_f90, 
      cnum_tot_f90, cnum_h2so4_f90, cnum_nh3_f90, radius_cluster_f90);

    REQUIRE( flagaa2 == flagaa2_f90);
    REQUIRE( (fp_helper::equiv(ratenucl       , ratenucl_f90        , abs_tol) || fp_helper::rel(ratenucl          , ratenucl_f90       , rel_tol)) );
    REQUIRE( (fp_helper::equiv(rateloge       , rateloge_f90        , abs_tol) || fp_helper::rel(rateloge          , rateloge_f90       , rel_tol)) );
    REQUIRE( (fp_helper::equiv(cnum_tot       , cnum_tot_f90        , abs_tol) || fp_helper::rel(cnum_tot          , cnum_tot_f90       , rel_tol)) );
    REQUIRE( (fp_helper::equiv(cnum_h2so4     , cnum_h2so4_f90      , abs_tol) || fp_helper::rel(cnum_h2so4        , cnum_h2so4_f90     , rel_tol)) );
    REQUIRE( (fp_helper::equiv(cnum_nh3       , cnum_nh3_f90        , abs_tol) || fp_helper::rel(cnum_nh3          , cnum_nh3_f90       , rel_tol)) );
    REQUIRE( (fp_helper::equiv(radius_cluster , radius_cluster_f90  , abs_tol) || fp_helper::rel(radius_cluster    , radius_cluster_f90 , rel_tol)) );
  }
}

TEST_CASE("mer07_veh02_nuc_mosaic_1box", "mam_nucleation_fprocess") {
  /// Test the mer07_veh02_nuc_mosaic_1box function directly by calling both
  /// the original Fortran version and the new C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  using fp_helper = FloatingPoint<double>;
  const double abs_tol = 5.0e-12;
  const double rel_tol = 5.0e-12;
  // Define a pseudo-random generator [0-1) that is consistent across platforms.
  // Manually checked the first 100,000 values to be unique.
  const unsigned p0  = 987659;
  const unsigned p1  =  12373;
  long unsigned seed =  54319;
  auto random = [&]() {
    seed =  (seed * p1)%p0;
    return double(seed)/p0;
  };
  MAMNucleationProcess mam_nucleation_process;
  for (int i=0; i<1000; ++i) {
    const int newnuc_method_flagaa = random() < .5 ? 1+2*random() : 11+2*random();  // range 1,2,11,12
    const double dtnuc             = random();  
    const double temp_in           = 235   +   60*random();  // range 235-295
    const double rh_in             = 0.05  +   .9*random();  // range .05-.95
    const double press_in          = 96325 + 10000*random(); // pressure in Pascal, sea level=101,325
    const double zm_in             =   500 + 10000*random(); // layer midpoint height (m)
    const double pblh_in           =  1000 +  1000*random(); // boundary layer height (m) 
    const double qh2so4_cur        = random();               // mixing ratio
    const double qh2so4_avg        = random();               // mixing ratio
    const double qnh3_cur          = random();               // mixing ratio
    const double h2so4_uptkrate    = 100*random();           // h2so4 uptake rate to aerosol (1/s)
    const double mw_so4a_host      = random()/1000;          // mw of so4 aerosol in host code (g/mol)
    const int nsize                = 1+2*random();           // number of aerosol size bins. NOTE: nsize<=maxd_asize
    const int maxd_asize           = nsize + 2*random();     // dimension for dplom_sect, NOTE: nsize<=maxd_asize,
    const int ldiagaa              = 10*random();            // does not appear to be used.
    std::vector<double> dplom_sect_vec(maxd_asize);
    std::vector<double> dphim_sect_vec(maxd_asize);
    const double SECT_SCALE = 1.0e10;
    dplom_sect_vec[0] = random()/SECT_SCALE;
    for (int i=1; i<maxd_asize; ++i) {
      dplom_sect_vec[i]   = dplom_sect_vec[i-1] + random()/SECT_SCALE;
      dphim_sect_vec[i-1] = dplom_sect_vec[i];
    }
    dphim_sect_vec[maxd_asize-1] = dplom_sect_vec[maxd_asize-1]+random()/SECT_SCALE;

    const double *dplom_sect = dplom_sect_vec.data();   
    const double *dphim_sect = dphim_sect_vec.data();   

    const double adjust_factor_bin_tern_ratenucl = random();
    const double adjust_factor_pbl_ratenucl = random();
    mam_nucleation_process.set_adjust_factor_bin_tern_ratenucl(adjust_factor_bin_tern_ratenucl);
    mam_nucleation_process.set_adjust_factor_pbl_ratenucl(adjust_factor_pbl_ratenucl);


    int    isize_nuc    = 0;
    double qnuma_del    = 0;
    double qso4a_del    = 0;
    double qnh4a_del    = 0;
    double qh2so4_del   = 0;
    double qnh3_del     = 0;
    double dens_nh4so4a = 0;
    double dnclusterdt  = 0;

    int    isize_nuc_f    = 0;
    double qnuma_del_f    = 0;
    double qso4a_del_f    = 0;
    double qnh4a_del_f    = 0;
    double qh2so4_del_f   = 0;
    double qnh3_del_f     = 0;
    double dens_nh4so4a_f = 0;
    double dnclusterdt_f  = 0;

    mam_nucleation_process.mer07_veh02_nuc_mosaic_1box(
      newnuc_method_flagaa,
      dtnuc,
      temp_in,
      rh_in,
      press_in,
      zm_in,
      pblh_in,
      qh2so4_cur,
      qh2so4_avg,
      qnh3_cur,
      h2so4_uptkrate,
      mw_so4a_host,
      nsize,
      maxd_asize,
      dplom_sect,  
      dphim_sect,  
      isize_nuc,
      qnuma_del,
      qso4a_del,
      qnh4a_del,
      qh2so4_del,
      qnh3_del,
      dens_nh4so4a,
      ldiagaa,
      &dnclusterdt);

    mer07_veh02_nuc_mosaic_1box_bridge(
      adjust_factor_bin_tern_ratenucl, 
      adjust_factor_pbl_ratenucl,
      newnuc_method_flagaa,
      dtnuc,
      temp_in,
      rh_in,
      press_in,
      zm_in,
      pblh_in,
      qh2so4_cur,
      qh2so4_avg,
      qnh3_cur,
      h2so4_uptkrate,
      mw_so4a_host,
      nsize,
      maxd_asize,
      dplom_sect,  
      dphim_sect,  
      isize_nuc_f,
      qnuma_del_f,
      qso4a_del_f,
      qnh4a_del_f,
      qh2so4_del_f,
      qnh3_del_f,
      dens_nh4so4a_f,
      ldiagaa,
      &dnclusterdt_f);

    REQUIRE(                    isize_nuc    ==  isize_nuc_f-1  );
    REQUIRE( (fp_helper::equiv( qnuma_del     ,  qnuma_del_f    , abs_tol) || fp_helper::rel( qnuma_del      ,  qnuma_del_f    , rel_tol)) );
    REQUIRE( (fp_helper::equiv( qso4a_del     ,  qso4a_del_f    , abs_tol) || fp_helper::rel( qso4a_del      ,  qso4a_del_f    , rel_tol)) );
    REQUIRE( (fp_helper::equiv( qnh4a_del     ,  qnh4a_del_f    , abs_tol) || fp_helper::rel( qnh4a_del      ,  qnh4a_del_f    , rel_tol)) );
    REQUIRE( (fp_helper::equiv( qh2so4_del    ,  qh2so4_del_f   , abs_tol) || fp_helper::rel( qh2so4_del     ,  qh2so4_del_f   , rel_tol)) );
    REQUIRE( (fp_helper::equiv( qnh3_del      ,  qnh3_del_f     , abs_tol) || fp_helper::rel( qnh3_del       ,  qnh3_del_f     , rel_tol)) );
    REQUIRE( (fp_helper::equiv( dens_nh4so4a  ,  dens_nh4so4a_f , abs_tol) || fp_helper::rel( dens_nh4so4a   ,  dens_nh4so4a_f , rel_tol)) );
    REQUIRE( (fp_helper::equiv( dnclusterdt   ,  dnclusterdt_f  , abs_tol) || fp_helper::rel( dnclusterdt    ,  dnclusterdt_f  , rel_tol)) );
  }
}

// These tests exercise our transplant of the MAM nucleation process.
TEST_CASE("MAMNucleationFProcess", "mam_nucleation_fprocess") {

  // We create a phony model to be used for these tests.
  auto aero_config = create_mam4_modal_aerosol_config();
  int num_levels = 72;
  auto* model = Model::ForUnitTests(aero_config, num_levels);
  int num_gases = aero_config.h_gas_species.size();
  int num_modes = aero_config.h_aerosol_modes.size();

  // Set up some prognosics aerosol data views‥
  int num_aero_populations = model->num_aerosol_populations();
  Kokkos::View<PackType**> int_aerosols("interstitial aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> cld_aerosols("cloudborne aerosols",
                                        num_aero_populations, num_levels);
  Kokkos::View<PackType**> gases("gases", num_gases, num_levels);
  Kokkos::View<PackType**> modal_concs("modal number concs", num_modes,
                                       num_levels);

  // Set up atmospheric data and populate it with some views.
  Kokkos::View<PackType*> temp("temperature", num_levels);
  Kokkos::View<PackType*> press("pressure", num_levels);
  Kokkos::View<PackType*> rel_hum("relative humidity", num_levels);
  Kokkos::View<PackType*> ht("height", num_levels+1);
  Real pblh = 100.0;
  auto* atm = new Atmosphere(num_levels, temp, press, rel_hum, ht, pblh);

  // Test basic construction.
  SECTION("construct") {
    auto* process = new MAMNucleationFProcess();
    REQUIRE(process->type() == haero::NucleationProcess);
    REQUIRE(process->name() == "MAMNucleationFProcess (Fortran prognostic NucleationProcess)");
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
                                            modal_concs);
    auto* diags = model->create_diagnostics();
    auto* tends = new Tendencies(*progs);

    // Set initial conditions.

    // atmospheric state
    Real h0 = 3e3, dz = h0/num_levels;
    for (int k = 0; k < num_levels; ++k) {
      temp(k) = 273.0;
      press(k) = 1e5;
      rel_hum(k) = 0.95;
      ht(k) = h0 - k*dz;
    }

    // aerosols (none)
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        int_aerosols(p, k) = 0.0;
      }
    }

    // gases
    int h2so4_index = aero_config.gas_index("H2SO4");
    printf("h2so4 index: %d\n", h2so4_index);
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
    const auto aero_tends = tends->interstitial_aerosols();
    int aitken_index = aero_config.aerosol_mode_index("aitken");
    int aitken_so4_index = aero_config.aerosol_species_index(aitken_index, "SO4");
    int so4_pop_index = aero_config.population_index(aitken_index, aitken_so4_index);
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      if (p == so4_pop_index) {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(int_aerosols(p, k)[0] > 0.0);
        }
      } else {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(int_aerosols(p, k)[0] == 0.0);
        }
      }
    }

    // The tendency for H2SO4 should be negative, and the rest should be zero.
    const auto gas_tends = tends->gases();
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

  delete model;
}

