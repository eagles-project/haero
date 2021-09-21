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
  const int num_atm_columns = 1245;
  using View1D = Kokkos::View<PackType *>;

  auto aero_config = create_mam4_modal_aerosol_config();
  static constexpr std::size_t num_levels{72};  // number of levels
  std::size_t num_vert_packs = num_levels / HAERO_PACK_SIZE;
  if (num_vert_packs * HAERO_PACK_SIZE < num_levels) {
    num_vert_packs++;
  }
  std::size_t num_iface_packs = (num_levels + 1) / HAERO_PACK_SIZE;
  if (num_iface_packs * HAERO_PACK_SIZE < (num_levels + 1)) {
    num_iface_packs++;
  }

  auto *model = get_model_for_unit_tests(aero_config, num_levels);
  const std::size_t num_gases{
      aero_config.gas_species.size()};  // number of gases
  const std::size_t num_modes{
      aero_config.aerosol_modes.size()};  // number of modes

  // Set up some prognostics aerosol data views
  const int num_aero_populations = model->num_aerosol_populations();
  using GlobalSpeciesColumnView = kokkos_device_type::view_3d<PackType>;
  using GlobalModeColumnView = DeviceType::view_3d<PackType>;

  static constexpr unsigned p0{987659};
  static constexpr unsigned p1{12373};
  long unsigned seed{54319};
  auto random = [&]() {
    seed = (seed * p1) % p0;
    return Real(seed) / p0;
  };

  GlobalSpeciesColumnView global_int_aerosols(
      "global interstitial aerosols", num_atm_columns,
      num_aero_populations,
      num_vert_packs);  // interstitial aerosols mmr [kg/kg(of air)]
  GlobalSpeciesColumnView global_cld_aerosols(
      "global cloudborne aerosols", num_atm_columns,
      num_aero_populations,
      num_vert_packs);  // cloud borne aerosols mmr [kg/kg(of air)]
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
  GlobalSpeciesColumnView global_gases("global gases",
                                       num_atm_columns,
                                       num_gases, num_vert_packs);

  GlobalModeColumnView global_int_num_concs(
      "global interstitial number concs", num_atm_columns,
      num_modes,
      num_vert_packs);  // interstitial aerosols number
                        // mixing ratios [#/kg(of air)]
  GlobalModeColumnView global_cld_num_concs(
      "global cloud borne number concs", num_atm_columns,
      num_modes,
      num_vert_packs);  // cloud borne aerosols number
                        // mixing ratios [#/kg(of air)]
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
  // Set up atmospheric data and populate it with some views.
  using View2D = Kokkos::View<PackType **>;
  View2D global_temp("temperature", num_atm_columns,
                     num_vert_packs);  //[K]
  View2D global_press("pressure", num_atm_columns,
                      num_vert_packs);  //[Pa]
  View2D global_rel_hum("relative humidity", num_atm_columns,
                        num_vert_packs);
  View2D global_qv("vapor mixing ratio", num_atm_columns,
                   num_vert_packs);
  View2D global_pdel("hydrostatic_dp", num_atm_columns,
                     num_vert_packs);  //[Pa]
  View2D global_ht("height", num_atm_columns, num_iface_packs);

  kokkos_device_type::view_1d<Diagnostics> global_diagnostics(
      "global diags", num_atm_columns);
  ;
  kokkos_device_type::view_1d<Atmosphere> global_atmosphere(
      "global Atmosphere", num_atm_columns);
  ;
  kokkos_device_type::view_1d<Prognostics> global_prognostics(
      "global Prognostics", num_atm_columns);
  ;
  kokkos_device_type::view_1d<Tendencies> global_tendencies(
      "global Tendencies", num_atm_columns);
  ;

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

    view_1d_int_type num_aero_species = vector_to_basic_1dview(
        model->calc_num_aero_species(), "Number aero species");
    for (int column = 0; column < num_atm_columns; ++column) {
      SpeciesColumnView int_aerosols = Kokkos::subview(
          global_int_aerosols, column, Kokkos::ALL, Kokkos::ALL);
      SpeciesColumnView cld_aerosols = Kokkos::subview(
          global_cld_aerosols, column, Kokkos::ALL, Kokkos::ALL);
      SpeciesColumnView gases =
          Kokkos::subview(global_gases, column, Kokkos::ALL, Kokkos::ALL);
      ;

      ModeColumnView int_num_concs = Kokkos::subview(
          global_int_num_concs, column, Kokkos::ALL, Kokkos::ALL);
      ;
      ModeColumnView cld_num_concs = Kokkos::subview(
          global_cld_num_concs, column, Kokkos::ALL, Kokkos::ALL);
      ;

      View1D temp = Kokkos::subview(global_temp, column, Kokkos::ALL);    //[K]
      View1D press = Kokkos::subview(global_press, column, Kokkos::ALL);  //[Pa]
      View1D rel_hum = Kokkos::subview(global_rel_hum, column, Kokkos::ALL);
      View1D qv = Kokkos::subview(global_qv, column, Kokkos::ALL);
      View1D pdel = Kokkos::subview(global_pdel, column, Kokkos::ALL);  //[Pa]
      View1D ht = Kokkos::subview(global_ht, column, Kokkos::ALL);
      Real pblh{100.0};  // planetary BL height [m]

      const Atmosphere atm(num_levels, temp, press, qv, ht, pdel, pblh);
      const Prognostics progs(num_aero_species.extent(0), num_aero_species,
                              num_gases, num_levels, int_aerosols, cld_aerosols,
                              int_num_concs, cld_num_concs, gases);
      const Tendencies tends(progs);
      HostDiagnostics *diags = model->create_diagnostics();

      host_atmos(column) = atm;
      host_progs(column) = progs;
      host_tends(column) = tends;
      host_diags(column) = *diags;
      delete diags;
    }
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
    Kokkos::parallel_for(
        teamPolicy, KOKKOS_LAMBDA(const TeamHandleType &team) {
          const int column = team.league_rank();

          Diagnostics &diags = global_diagnostics(column);
          Atmosphere &atmos = global_atmosphere(column);
          Prognostics &progs = global_prognostics(column);
          Tendencies &tends = global_tendencies(column);
          d_process->run(team, t, dt, progs, atmos, diags, tends);
        });
    AerosolProcess::delete_on_device(d_process);
  }
  // All the class instances have been stored in Kokkos::View and should be
  // deleted when the views go out of scope.
}
