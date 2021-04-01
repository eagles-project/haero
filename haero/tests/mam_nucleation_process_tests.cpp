#include "haero/model.hpp"
#include "haero/floating_point.hpp"
#include "haero/processes/mam_nucleation_process.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace haero;

TEST_CASE("ternary_nuc_merik2007_test", "mam_nucleation_process") {
  /// Test the ternary_nuc_merik2007 function directly by calling both
  /// the Nvidia Cuda GPU version and the CPU C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  using fp_helper = FloatingPoint<double>;
  using SolutionView = DeviceType::view_1d<double>; 
  const double tolerance = 5.0e-9;
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

    SolutionView solution("ternary_nuc_merik2007",5);
    Kokkos::parallel_for("ternary_nuc_merik2007.mam_nucleation_process", 1,
      KOKKOS_LAMBDA(const int) {
        double j_log=0, ntot=0, nacid=0, namm=0, r=0;
        MAMNucleationProcess::ternary_nuc_merik2007(t, rh, c2, c3, j_log, ntot, nacid, namm, r);
        solution(0) = j_log;
        solution(1) = ntot;
        solution(2) = nacid;
        solution(3) = namm;
        solution(4) = r;
      }
    );
    auto h_solution = Kokkos::create_mirror_view(solution);
    Kokkos::deep_copy(h_solution, solution);
    const double j_log_kok = h_solution(0);
    const double ntot_kok  = h_solution(1);
    const double nacid_kok = h_solution(2);
    const double namm_kok  = h_solution(3);
    const double r_kok     = h_solution(4);

    double j_log_cpp = 0;
    double ntot_cpp  = 0; 
    double nacid_cpp = 0;
    double namm_cpp  = 0; 
    double r_cpp     = 0;
    MAMNucleationProcess::ternary_nuc_merik2007(t, rh, c2, c3, j_log_cpp, ntot_cpp, nacid_cpp, namm_cpp, r_cpp);
    REQUIRE(fp_helper::equiv(j_log_cpp , j_log_kok, tolerance));
    REQUIRE(fp_helper::equiv(ntot_cpp  , ntot_kok,  tolerance));
    REQUIRE(fp_helper::equiv(nacid_cpp , nacid_kok, tolerance));
    REQUIRE(fp_helper::equiv(namm_cpp  , namm_kok,  tolerance));
    REQUIRE(fp_helper::equiv(r_cpp     , r_kok,     tolerance));
  }
}

TEST_CASE("binary_nuc_vehk2002", "mam_nucleation_process") {
  /// Test the binary_nuc_vehk2002 function directly by calling both
  /// the Nvidia Cuda GPU version and the CPU C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  using fp_helper = FloatingPoint<double>;
  using SolutionView = DeviceType::view_1d<double>; 
  const double tolerance = 5.0e-12;
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

    SolutionView solution("binary_nuc_vehk2002",5);
    Kokkos::parallel_for("binary_nuc_vehk2002.mam_nucleation_process", 1,
      KOKKOS_LAMBDA(const int) {
        double ratenucl=0, rateloge=0, cnum_h2so4=0, cnum_tot=0, radius_cluster=0;
        MAMNucleationProcess::binary_nuc_vehk2002(temp, rh, so4vol, ratenucl, rateloge, cnum_h2so4, cnum_tot, radius_cluster);
        solution(0) = ratenucl;
        solution(1) = rateloge;
        solution(2) = cnum_h2so4;
        solution(3) = cnum_tot;
        solution(4) = radius_cluster;
      }
    );
    auto h_solution = Kokkos::create_mirror_view(solution);
    Kokkos::deep_copy(h_solution, solution);
    const double ratenucl_kok       = h_solution(0);
    const double rateloge_kok       = h_solution(1);
    const double cnum_h2so4_kok     = h_solution(2);
    const double cnum_tot_kok       = h_solution(3);
    const double radius_cluster_kok = h_solution(4);

    double ratenucl_cpu       = 0;
    double rateloge_cpu       = 0; 
    double cnum_h2so4_cpu     = 0;
    double cnum_tot_cpu       = 0; 
    double radius_cluster_cpu = 0;
    MAMNucleationProcess::binary_nuc_vehk2002(temp, rh, so4vol, ratenucl_cpu, rateloge_cpu, cnum_h2so4_cpu, cnum_tot_cpu, radius_cluster_cpu);
    REQUIRE( (fp_helper::equiv(ratenucl_cpu       , ratenucl_kok       , tolerance) || fp_helper::rel(ratenucl_cpu       , ratenucl_kok       , tolerance)) );
    REQUIRE( (fp_helper::equiv(rateloge_cpu       , rateloge_kok       , tolerance) || fp_helper::rel(rateloge_cpu       , rateloge_kok       , tolerance)) );
    REQUIRE( (fp_helper::equiv(cnum_h2so4_cpu     , cnum_h2so4_kok     , tolerance) || fp_helper::rel(cnum_h2so4_cpu     , cnum_h2so4_kok     , tolerance)) );
    REQUIRE( (fp_helper::equiv(cnum_tot_cpu       , cnum_tot_kok       , tolerance) || fp_helper::rel(cnum_tot_cpu       , cnum_tot_kok       , tolerance)) );
    REQUIRE( (fp_helper::equiv(radius_cluster_cpu , radius_cluster_kok , tolerance) || fp_helper::rel(radius_cluster_cpu , radius_cluster_kok , tolerance)) );
  }
}


TEST_CASE("pbl_nuc_wang2008", "mam_nucleation_process") {
  /// Test the pbl_nuc_wang2008 function directly by calling both
  /// the Nvidia Cuda GPU version and the CPU C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  using fp_helper = FloatingPoint<double>;
  using SolutionView = DeviceType::view_1d<double>; 
  using FlagaaView   = DeviceType::view_1d<int>; 
  const double tolerance = 5.0e-12;
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

    SolutionView solution("pbl_nuc_wang2008",6);
    FlagaaView   flags("newnuc_method_flagaa",1);
    
    Kokkos::parallel_for("pbl_nuc_wang2008.mam_nucleation_process", 1,
      KOKKOS_LAMBDA(const int) {
        int flagaa2 = 0;
        double ratenucl=0, rateloge=0, cnum_tot=0, cnum_h2so4=0, cnum_nh3=0, radius_cluster=0;
        mam_nucleation_process.pbl_nuc_wang2008(so4vol, flagaa, flagaa2, ratenucl, rateloge,
          cnum_tot, cnum_h2so4, cnum_nh3, radius_cluster);
        solution(0) = ratenucl;        
        solution(1) = rateloge;        
        solution(2) = cnum_tot;        
        solution(3) = cnum_h2so4;      
        solution(4) = cnum_nh3;        
        solution(5) = radius_cluster; 
        flags(0)    = flagaa2;
      }
    );
    auto h_solution = Kokkos::create_mirror_view(solution);
    auto h_flags = Kokkos::create_mirror_view(flags);
    Kokkos::deep_copy(h_solution, solution);
    Kokkos::deep_copy(h_flags, flags);
    const double ratenucl_kok       = h_solution(0);
    const double rateloge_kok       = h_solution(1);
    const double cnum_tot_kok       = h_solution(2);
    const double cnum_h2so4_kok     = h_solution(3);
    const double cnum_nh3_kok       = h_solution(4);
    const double radius_cluster_kok = h_solution(5);
    const int    flagaa2_kok        = h_flags(0);

    double ratenucl_cpu       = 0;
    double rateloge_cpu       = 0; 
    double cnum_tot_cpu       = 0;
    double cnum_h2so4_cpu     = 0; 
    double cnum_nh3_cpu       = 0;
    double radius_cluster_cpu = 0;
    int    flagaa2_cpu        = 0; 

    mam_nucleation_process.pbl_nuc_wang2008(so4vol, flagaa, flagaa2_cpu, ratenucl_cpu, rateloge_cpu, cnum_tot_cpu, cnum_h2so4_cpu, cnum_nh3_cpu, radius_cluster_cpu);
    REQUIRE( (fp_helper::equiv(ratenucl_cpu      , ratenucl_kok        , tolerance) || fp_helper::rel(ratenucl_cpu      , ratenucl_kok        , tolerance)) );
    REQUIRE( (fp_helper::equiv(rateloge_cpu      , rateloge_kok        , tolerance) || fp_helper::rel(rateloge_cpu      , rateloge_kok        , tolerance)) );
    REQUIRE( (fp_helper::equiv(cnum_tot_cpu      , cnum_tot_kok        , tolerance) || fp_helper::rel(cnum_tot_cpu      , cnum_tot_kok        , tolerance)) );
    REQUIRE( (fp_helper::equiv(cnum_h2so4_cpu    , cnum_h2so4_kok      , tolerance) || fp_helper::rel(cnum_h2so4_cpu    , cnum_h2so4_kok      , tolerance)) );
    REQUIRE( (fp_helper::equiv(cnum_nh3_cpu      , cnum_nh3_kok        , tolerance) || fp_helper::rel(cnum_nh3_cpu      , cnum_nh3_kok        , tolerance)) );
    REQUIRE( (fp_helper::equiv(radius_cluster_cpu, radius_cluster_kok  , tolerance) || fp_helper::rel(radius_cluster_cpu, radius_cluster_kok  , tolerance)) );
    REQUIRE( flagaa2_cpu ==flagaa2_kok );
  }
}

TEST_CASE("mer07_veh02_nuc_mosaic_1box", "mam_nucleation_process") {
  /// Test the mer07_veh02_nuc_mosaic_1box function directly by calling both
  /// the Nvidia Cuda GPU version and the CPU C++ version and compare
  /// the result. The testing process is to generate a bunch of random
  /// input values and check the output values are close.  Differences
  /// in Fortran and C++ means the result is not identical but we hope
  /// it is within numerical round off.

  using fp_helper = FloatingPoint<double>;
  using SolutionView = DeviceType::view_1d<double>; 
  using FlagaaView   = DeviceType::view_1d<int>; 
  const double tolerance = 1.0e-08;
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
    const int ldiagaa              = 10*random();             // does not appear to be used.
    std::vector<double> dplom_sect(maxd_asize);
    std::vector<double> dphim_sect(maxd_asize);
    const double SECT_SCALE = 1.0e10;
    dplom_sect[0] = random()/SECT_SCALE;
    for (int i=1; i<maxd_asize; ++i) {
      dplom_sect[i]   = dplom_sect[i-1] + random()/SECT_SCALE;
      dphim_sect[i-1] = dplom_sect[i];
    }
    dphim_sect[maxd_asize-1] = dplom_sect[maxd_asize-1]+random()/SECT_SCALE;

    const double adjust_factor_bin_tern_ratenucl = random();
    const double adjust_factor_pbl_ratenucl = random();
    mam_nucleation_process.set_adjust_factor_bin_tern_ratenucl(adjust_factor_bin_tern_ratenucl);
    mam_nucleation_process.set_adjust_factor_pbl_ratenucl(adjust_factor_pbl_ratenucl);


    SolutionView solution("mer07_veh02_nuc_mosaic_1box",7);
    FlagaaView   flags("newnuc_method_flagaa",1);

    SolutionView dplom_sect_view("dplom_sect", dplom_sect.size());
    SolutionView dphim_sect_view("dphim_sect", dphim_sect.size());
    auto h_dplom_sect_view = Kokkos::create_mirror_view(dplom_sect_view);
    auto h_dphim_sect_view = Kokkos::create_mirror_view(dphim_sect_view);
    for (int i=0; i<maxd_asize; ++i) {
      h_dplom_sect_view(i) = dplom_sect[i];
      h_dphim_sect_view(i) = dphim_sect[i];
    }  
    Kokkos::deep_copy(dplom_sect_view, h_dplom_sect_view);
    Kokkos::deep_copy(dphim_sect_view, h_dphim_sect_view);
    
    Kokkos::parallel_for("mer07_veh02_nuc_mosaic_1box.mam_nucleation_process", 1,
      KOKKOS_LAMBDA(const int) {
        int    isize_nuc    = 0;
        double qnuma_del    = 0;
        double qso4a_del    = 0;
        double qnh4a_del    = 0;
        double qh2so4_del   = 0;
        double qnh3_del     = 0;
        double dens_nh4so4a = 0;
        double dnclusterdt  = 0;
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
          dplom_sect_view.data(),
          dphim_sect_view.data(),
          isize_nuc,
          qnuma_del,
          qso4a_del,
          qnh4a_del,
          qh2so4_del,
          qnh3_del,
          dens_nh4so4a,
          ldiagaa,
          &dnclusterdt);

        solution(0) = qnuma_del    ;        
        solution(1) = qso4a_del    ;        
        solution(2) = qnh4a_del    ;        
        solution(3) = qh2so4_del   ;      
        solution(4) = qnh3_del     ;        
        solution(5) = dens_nh4so4a ; 
        solution(6) = dnclusterdt  ;
        flags(0)    = isize_nuc    ;
      }
    );
    auto h_solution = Kokkos::create_mirror_view(solution);
    auto h_flags = Kokkos::create_mirror_view(flags);
    Kokkos::deep_copy(h_solution, solution);
    Kokkos::deep_copy(h_flags, flags);
    const double qnuma_del_kok     = h_solution(0);
    const double qso4a_del_kok     = h_solution(1);
    const double qnh4a_del_kok     = h_solution(2);
    const double qh2so4_del_kok    = h_solution(3);
    const double qnh3_del_kok      = h_solution(4);
    const double dens_nh4so4a_kok  = h_solution(5);
    const double dnclusterdt_kok   = h_solution(6);
    const int    isize_nuc_kok     = h_flags(0);

    double qnuma_del_cpu    = 0;
    double qso4a_del_cpu    = 0; 
    double qnh4a_del_cpu    = 0;
    double qh2so4_del_cpu   = 0; 
    double qnh3_del_cpu     = 0;
    double dens_nh4so4a_cpu = 0;
    double dnclusterdt_cpu  = 0;
    int    isize_nuc_cpu    = 0; 

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
          dplom_sect.data(),
          dphim_sect.data(),
          isize_nuc_cpu,
          qnuma_del_cpu,
          qso4a_del_cpu,
          qnh4a_del_cpu,
          qh2so4_del_cpu,
          qnh3_del_cpu,
          dens_nh4so4a_cpu,
          ldiagaa,
          &dnclusterdt_cpu);

    REQUIRE( (fp_helper::equiv( qnuma_del_cpu    ,  qnuma_del_kok     , tolerance) || fp_helper::rel( qnuma_del_cpu    ,  qnuma_del_kok    , tolerance)) );
    REQUIRE( (fp_helper::equiv( qso4a_del_cpu    ,  qso4a_del_kok     , tolerance) || fp_helper::rel( qso4a_del_cpu    ,  qso4a_del_kok    , tolerance)) );
    REQUIRE( (fp_helper::equiv( qnh4a_del_cpu    ,  qnh4a_del_kok     , tolerance) || fp_helper::rel( qnh4a_del_cpu    ,  qnh4a_del_kok    , tolerance)) );
    REQUIRE( (fp_helper::equiv( qh2so4_del_cpu   ,  qh2so4_del_kok    , tolerance) || fp_helper::rel( qh2so4_del_cpu   ,  qh2so4_del_kok   , tolerance)) );
    REQUIRE( (fp_helper::equiv( qnh3_del_cpu     ,  qnh3_del_kok      , tolerance) || fp_helper::rel( qnh3_del_cpu     ,  qnh3_del_kok     , tolerance)) );
    REQUIRE( (fp_helper::equiv( dens_nh4so4a_cpu ,  dens_nh4so4a_kok  , tolerance) || fp_helper::rel( dens_nh4so4a_cpu ,  dens_nh4so4a_kok , tolerance)) );
    REQUIRE( (fp_helper::equiv( dnclusterdt_cpu  ,  dnclusterdt_kok   , tolerance) || fp_helper::rel( dnclusterdt_cpu  ,  dnclusterdt_kok  , tolerance)) );
    REQUIRE( isize_nuc_cpu == isize_nuc_kok );
  }
}

// These tests exercise our transplant of the MAM nucleation process.
TEST_CASE("virtual_process_test", "mam_nucleation_process") {

  // We create a phony model to be used for these tests.
  auto aero_config = create_mam4_modal_aerosol_config();
  int num_levels = 72;
  int num_vert_packs = num_levels/HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  int num_iface_packs = (num_levels+1)/HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels+1)) {
    num_iface_packs++;
  }
  auto* model = Model::ForUnitTests(aero_config, num_levels);
  int num_gases = aero_config.h_gas_species.size();
  int num_modes = aero_config.h_aerosol_modes.size();

  // Set up some prognosics aerosol data views‥
  int num_aero_populations = model->num_aerosol_populations();
  SpeciesColumnView        int_aerosols("interstitial aerosols",
                                        num_aero_populations, num_vert_packs);
  SpeciesColumnView        cld_aerosols("cloudborne aerosols",
                                        num_aero_populations, num_vert_packs);
  SpeciesColumnView        gases("gases", num_gases, num_vert_packs);
  ModalColumnView          modal_concs("modal number concs", num_modes,
                                       num_vert_packs);

  // Set up atmospheric data and populate it with some views.
  ColumnView temp("temperature", num_vert_packs);
  ColumnView press("pressure", num_vert_packs);
  ColumnView rel_hum("relative humidity", num_vert_packs);
  ColumnView ht("height", num_iface_packs);
  Real pblh = 100.0;

  auto* atm = new Atmosphere(num_levels, temp, press, rel_hum, ht, pblh);

  // Test basic construction.
  SECTION("construct") {
    auto* process = new MAMNucleationProcess();
    REQUIRE(process->type() == haero::NucleationProcess);
    REQUIRE(process->name() == "MAMNucleationProcess");
    delete process;
  }

  // Test process initialization.
  SECTION("init_process") {
    auto* process = new MAMNucleationProcess();
    process->init(aero_config);
    delete process;
  }

  // Test process tendencies.
  SECTION("nucleate_without_existing_aerosols") {
    auto* process = new MAMNucleationProcess();
    process->init(aero_config);

    // Initialize prognostic and diagnostic variables, and construct a
    // tendencies container.
    Prognostics* progs = model->create_prognostics(int_aerosols, cld_aerosols, gases, modal_concs);
    Diagnostics* diags = model->create_diagnostics();
    Tendencies*  tends = new Tendencies(*progs);

    // Set initial conditions.

    // atmospheric state
    Real h0 = 3e3, dz = h0/num_levels;
    auto h_temp    = Kokkos::create_mirror_view(temp);
    auto h_press   = Kokkos::create_mirror_view(press);
    auto h_rel_hum = Kokkos::create_mirror_view(rel_hum);
    auto h_ht      = Kokkos::create_mirror_view(ht);
    for (int k = 0; k < num_levels; ++k) {
      h_temp(pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 273.0;
      h_press(pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 1e5;
      h_rel_hum(pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 0.95;
    }
    for (int k = 0; k < num_levels+1; ++k) {
      h_ht(pack_info::pack_idx(k))[pack_info::vec_idx(k)] = h0 - k*dz;
    }
    Kokkos::deep_copy(temp,    h_temp);
    Kokkos::deep_copy(press,   h_press);
    Kokkos::deep_copy(rel_hum, h_rel_hum);
    Kokkos::deep_copy(ht,      h_ht);

    // aerosols (none)
    auto h_int_aerosols = Kokkos::create_mirror_view(int_aerosols);
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      for (int k = 0; k < num_levels; ++k) {
        h_int_aerosols(p, pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 0.0;
      }
    }
    Kokkos::deep_copy(int_aerosols, h_int_aerosols);

    // gases
    int h2so4_index = aero_config.gas_index("H2SO4");
    printf("h2so4 index: %d\n", h2so4_index);
    auto h_gases = Kokkos::create_mirror_view(gases);
    for (int k = 0; k < num_levels; ++k) {
      h_gases(h2so4_index, pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 1e-13;
    }
    for (int g = 0; g < num_gases; ++g) {
      if (g != h2so4_index) {
        for (int k = 0; k < num_levels; ++k) {
          h_gases(g, pack_info::pack_idx(k))[pack_info::vec_idx(k)] = 0.0;
        }
      }
    }
    Kokkos::deep_copy(gases, h_gases);

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
    Kokkos::deep_copy(h_int_aerosols, int_aerosols);
    for (int p = 0; p < aero_config.num_aerosol_populations; ++p) {
      if (p == so4_pop_index) {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(h_int_aerosols(p, pack_info::pack_idx(k))[pack_info::vec_idx(k)] > 0.0);
        }
      } else {
        for (int k = 0; k < num_levels; ++k) {
          REQUIRE(h_int_aerosols(p, pack_info::pack_idx(k))[pack_info::vec_idx(k)] == 0.0);
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
