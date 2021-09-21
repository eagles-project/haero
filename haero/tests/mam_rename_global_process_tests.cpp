#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "haero/model.hpp"
#include "haero/processes/mam_rename_process.hpp"

using namespace haero;

Model *get_model_for_unit_tests(const ModalAerosolConfig &aero_config,
                                const std::size_t num_levels) {
  static Model *model(Model::ForUnitTests(aero_config, num_levels));
  return model;
}

TEST_CASE("mam_rename_global_run", "") {
  // Define an arbitrary number of columns for this rank to
  // process.  This would be out of the total number of
  // columns in the atmosphere from the whole-earth model.
  const int num_atm_columns = 1245;
  using View1D = kokkos_device_type::view_1d<PackType>;
  using View2D = kokkos_device_type::view_2d<PackType>;
  using GlobalSpeciesColumnView = kokkos_device_type::view_3d<PackType>;
  using GlobalModeColumnView = DeviceType::view_3d<PackType>;

  // Create a default aerosol configuration for testing
  // purposes only.
  ModalAerosolConfig aero_config = create_mam4_modal_aerosol_config();
  static constexpr std::size_t num_levels{72};  // number of levels
  std::size_t num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  std::size_t num_iface_packs = (num_levels + 1) / HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels + 1)) {
    num_iface_packs++;
  }

  Model *model = get_model_for_unit_tests(aero_config, num_levels);
  const std::size_t num_gases{
      aero_config.gas_species.size()};  // number of gases
  const std::size_t num_modes{
      aero_config.aerosol_modes.size()};  // number of modes

  // Set up some prognostics aerosol data views
  const int num_aero_populations = model->num_aerosol_populations();

  // Define a pseudo-random number generator [0,1) based on a few primes.
  static constexpr unsigned p0{987659};
  static constexpr unsigned p1{12373};
  long unsigned seed{54319};
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };

  // interstitial aerosols mmr [kg/kg(of air)]
  GlobalSpeciesColumnView global_int_aerosols("global interstitial aerosols",
                                              num_atm_columns,
                                              num_aero_populations,
                                              num_vert_packs);
  // cloud borne aerosols mmr [kg/kg(of air)]
  GlobalSpeciesColumnView global_cld_aerosols("global cloudborne aerosols",
                                              num_atm_columns,
                                              num_aero_populations,
                                              num_vert_packs);
  {
    // Set initial conditions
    // aerosols mass mixing ratios
    auto h_int_aerosols = Kokkos::create_mirror_view(global_int_aerosols);
    auto h_cld_aerosols = Kokkos::create_mirror_view(global_cld_aerosols);
    for (int column = 0; column < num_atm_columns; ++column) {
      for (std::size_t p = 0; p < num_aero_populations; ++p) {
        for (std::size_t k = 0; k < num_vert_packs; ++k) {
          h_int_aerosols(column, p, k) = random() * 10e-10;
          h_cld_aerosols(column, p, k) = random() * 10e-10;
        }
      }
    }
    Kokkos::deep_copy(global_int_aerosols, h_int_aerosols);
    Kokkos::deep_copy(global_cld_aerosols, h_cld_aerosols);
  }
  GlobalSpeciesColumnView global_gases("global gases", num_atm_columns, num_gases,
                                       num_vert_packs);

  GlobalModeColumnView global_int_num_concs(
      "global interstitial number concs",
      num_atm_columns,
      num_modes,       // interstitial aerosols number
      num_vert_packs); // mixing ratios [#/kg(of air)]
  GlobalModeColumnView global_cld_num_concs(
      "global cloud borne number concs",
      num_atm_columns,
      num_modes,       // cloud borne aerosols number
      num_vert_packs); // mixing ratios [#/kg(of air)]
  {
    // aerosols number mixing ratios
    auto h_int_nmrs = Kokkos::create_mirror_view(global_int_num_concs);
    auto h_cld_nmrs = Kokkos::create_mirror_view(global_cld_num_concs);
    for (int column = 0; column < num_atm_columns; ++column) {
      for (std::size_t imode = 0; imode < num_modes; ++imode) {
        for (std::size_t k = 0; k < num_vert_packs; ++k) {
          h_int_nmrs(column, imode, k) = 1e8 + random();
          h_cld_nmrs(column, imode, k) = 1e8 + random();
        }
      }
    }
    Kokkos::deep_copy(global_int_num_concs, h_int_nmrs);
    Kokkos::deep_copy(global_cld_num_concs, h_cld_nmrs);
  }
  // Set up atmospheric data and populate it.
  // Each view is vertical_levels x number_of_columns
  // These 2D Views are then sub-viewed for each
  // column to set up the classes.
  View2D global_temp("temperature", num_atm_columns, num_vert_packs);  //[K]
  {
    auto host_temp = Kokkos::create_mirror_view(global_temp);
    for (int i=0; i<num_vert_packs; ++i)
      for (int j=0; j<num_atm_columns; ++j)
        for (int s=0; s<HAERO_PACK_SIZE; ++s) 
          host_temp(j,i)[s] = 100*j + i;
    Kokkos::deep_copy(global_temp, host_temp);
  }
  View2D global_press("pressure", num_atm_columns, num_vert_packs);  //[Pa]
  View2D global_rel_hum("relative humidity", num_atm_columns, num_vert_packs);
  View2D global_qv("vapor mixing ratio", num_atm_columns, num_vert_packs);
  View2D global_pdel("hydrostatic_dp", num_atm_columns, num_vert_packs);  //[Pa]
  View2D global_ht("height", num_atm_columns, num_iface_packs);

  // Views for the Hearo classes.  These are views of class instances
  // not pointers to instances.  This means that the classes must
  // have default constructors and copy operators to define them.
  // There is one class instance for every column of the atmosphere.
  // These are device views so will have to be deep copied in order
  // to define the class instances.
  kokkos_device_type::view_1d<Diagnostics> global_diagnostics("global diags",
                                                              num_atm_columns);
  kokkos_device_type::view_1d<Atmosphere> global_atmosphere("global Atmosphere",
                                                            num_atm_columns);
  kokkos_device_type::view_1d<Prognostics> global_prognostics(
      "global Prognostics", num_atm_columns);
  kokkos_device_type::view_1d<Tendencies> global_tendencies("global Tendencies",
                                                            num_atm_columns);

  {  // Create the device copies of diagnostics,, atmosphere, prognostics and
     // tendencies
    // These are created on host and copied to device through a
    // Kokkos::deep_copy.  This works for POD and Kokkos::Views and should
    // even allow std member data to be copied to the device without
    // error but calls to std classes are not available on device code so
    // the data would not be usable.
    auto host_diags = Kokkos::create_mirror_view(global_diagnostics);
    auto host_atmos = Kokkos::create_mirror_view(global_atmosphere);
    auto host_progs = Kokkos::create_mirror_view(global_prognostics);
    auto host_tends = Kokkos::create_mirror_view(global_tendencies);

    // Loop over all of the columns to define the class data.
    for (int column = 0; column < num_atm_columns; ++column) {
      // Kokkos views can be sub-viewed in the constructor using
      // the special ALL flag and selecting a single column.
      //
      // All of these Kokkos views are Device Views!  Do not
      // try to access them on Host without a create_mirrow_view.
      // So, even though host_atmos is an instance of an
      // Atmosphere that is defined on the host, the Views
      // used to define it are Device views.
      SpeciesColumnView int_aerosols(global_int_aerosols, column, Kokkos::ALL, Kokkos::ALL);
      SpeciesColumnView cld_aerosols(global_cld_aerosols, column, Kokkos::ALL, Kokkos::ALL);
      SpeciesColumnView gases(global_gases, column, Kokkos::ALL, Kokkos::ALL);
      ModeColumnView int_num_concs(global_int_num_concs, column, Kokkos::ALL, Kokkos::ALL);
      ModeColumnView cld_num_concs(global_cld_num_concs, column, Kokkos::ALL, Kokkos::ALL);

      View1D temp(global_temp, column, Kokkos::ALL);    //[K]


      View1D press(global_press, column, Kokkos::ALL);  //[Pa]
      View1D rel_hum(global_rel_hum, column, Kokkos::ALL);
      View1D qv(global_qv, column, Kokkos::ALL);
      View1D pdel(global_pdel, column, Kokkos::ALL);  //[Pa]
      View1D ht(global_ht, column, Kokkos::ALL);
      Real pblh{100.0};  // planetary BL height [m]

      // Create the Hearo class instances from the Kokkos
      // Views.  The diagnostics and tendencies classes 
      // define and hold their own Kokkos Views.  In an acutal
      // application these two classes would be checked on
      // the host to make sure the correct Views are defined.
      Prognostics *progs = model->create_prognostics(
          int_aerosols, cld_aerosols, int_num_concs, cld_num_concs, gases);
      HostDiagnostics *diags = model->create_diagnostics();
      const Tendencies tends(*progs);
      const Atmosphere atm(num_levels, temp, press, qv, ht, pdel, pblh);

      host_atmos(column) = atm;
      host_tends(column) = tends;
      host_progs(column) = *progs;
      host_diags(column) = *diags;
      delete diags;
      delete progs;
    }
    // Now that the classes are set up on host, can
    // be copied back to device.
    Kokkos::deep_copy(global_atmosphere, host_atmos);
    Kokkos::deep_copy(global_prognostics, host_progs);
    Kokkos::deep_copy(global_tendencies, host_tends);
    Kokkos::deep_copy(global_diagnostics, host_diags);
  }

  SECTION("rename_global_run") {
    //  Now run the analysis over the global atmosphere model.
    auto process = new MAMRenameProcess();
    process->init(aero_config);
    auto d_process = process->copy_to_device();
    typedef ekat::ExeSpaceUtils<>::TeamPolicy::member_type TeamHandleType;
    const auto &teamPolicy = ekat::ExeSpaceUtils<>::get_default_team_policy(
        num_atm_columns, num_vert_packs);
    Real t = 0.0, dt = 30.0;
    // Parallel launch on GPU or over threads.
    Kokkos::parallel_for(
        teamPolicy, KOKKOS_LAMBDA(const TeamHandleType &team) {
          const int column = team.league_rank();

          // Since the Kokkos 
          Diagnostics &diags = global_diagnostics(column);
          Atmosphere &atmos = global_atmosphere(column);
          Prognostics &progs = global_prognostics(column);
          Tendencies &tends = global_tendencies(column);
          d_process->run(team, t, dt, progs, atmos, diags, tends);
        });
    AerosolProcess::delete_on_device(d_process);
  }
  auto check_temp = [&](const int col, const int lev, const PackType &p) { 
    bool check = true;
    for (int s=0; s<HAERO_PACK_SIZE; ++s) 
      check = check && p[s] == 100*col + lev;
    return check; 
  };
  {
    // Copy back to host to check results.
    auto host_diags = Kokkos::create_mirror_view(global_diagnostics);
    auto host_atmos = Kokkos::create_mirror_view(global_atmosphere);
    auto host_progs = Kokkos::create_mirror_view(global_prognostics);
    auto host_tends = Kokkos::create_mirror_view(global_tendencies);
    Kokkos::deep_copy(host_diags, global_diagnostics);
    Kokkos::deep_copy(host_atmos, global_atmosphere);
    Kokkos::deep_copy(host_progs, global_prognostics);
    Kokkos::deep_copy(host_tends, global_tendencies);
    
    for (int column = 0; column < num_atm_columns; ++column) {
      const Atmosphere &atm = host_atmos(column);
      // To view the data in a column, deep_copy the View back to host
      auto host_temp = Kokkos::create_mirror_view(atm.temperature);
      Kokkos::deep_copy(host_temp, atm.temperature);
      for (int i=0; i<num_vert_packs; ++i) {
        REQUIRE(check_temp(column, i, host_temp(i)));
      }
    }
  }
  // All the class instances have been stored in Kokkos::View and should be
  // deleted when the views go out of scope.
}
