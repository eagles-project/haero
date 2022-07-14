#include <catch2/catch.hpp>
#include <cmath>
#include <haero/mam4/mam4.hpp>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <skywalker.hpp>

using namespace haero;
using namespace skywalker;

TEST_CASE("test_constructor", "mam4_gasaerexch_process") {
  mam4::AeroConfig mam4_config;
  mam4::GasAerExchProcess process(mam4_config);
  REQUIRE(process.name() == "MAM4 gas/aersol exchange");
  REQUIRE(process.config() == mam4_config);
}

TEST_CASE("test_compute_tendencies", "mam4_gasaerexch_process") {
  int nlev = 72;
  Real pblh = 1000;
  Atmosphere atm(nlev, pblh);
  mam4::Prognostics progs(nlev);
  mam4::Diagnostics diags(nlev);
  mam4::Tendencies tends(nlev);

  mam4::AeroConfig mam4_config;
  mam4::GasAerExchProcess process(mam4_config);

  // Single-column dispatch.
  auto team_policy = TeamPolicy(1u, Kokkos::AUTO);
  Real t = 0.0, dt = 30.0;
  Kokkos::parallel_for(
      team_policy, KOKKOS_LAMBDA(const TeamType& team) {
        process.compute_tendencies(team, t, dt, atm, progs, diags, tends);
      });
}

TEST_CASE("test_multicol_compute_tendencies", "mam4_gasaerexch_process") {
  // Now we process multiple columns within a single dіspatch (mc means
  // "multi-column").
  int ncol = 8;
  DeviceType::view_1d<Atmosphere> mc_atm("mc_progs", ncol);
  DeviceType::view_1d<mam4::Prognostics> mc_progs("mc_atm", ncol);
  DeviceType::view_1d<mam4::Diagnostics> mc_diags("mc_diags", ncol);
  DeviceType::view_1d<mam4::Tendencies> mc_tends("mc_tends", ncol);
  for (int icol = 0; icol < ncol; ++icol) {
    Kokkos::parallel_for(
        "Load multi-column views", 1, KOKKOS_LAMBDA(const int) {
          int nlev = 72;
          Real pblh = 1000;
          mc_atm(icol) = Atmosphere(nlev, pblh);
          mc_progs(icol) = mam4::Prognostics(nlev);
          mc_diags(icol) = mam4::Diagnostics(nlev);
          mc_tends(icol) = mam4::Tendencies(nlev);
        });
  }

  mam4::AeroConfig mam4_config;
  mam4::GasAerExchProcess process(mam4_config);

  // Dispatch over all the above columns.
  auto team_policy = TeamPolicy(ncol, Kokkos::AUTO);
  Real t = 0.0, dt = 30.0;
  Kokkos::parallel_for(
      team_policy, KOKKOS_LAMBDA(const TeamType& team) {
        const int icol = team.league_rank();
        process.compute_tendencies(team, t, dt, mc_atm(icol), mc_progs(icol),
                                   mc_diags(icol), mc_tends(icol));
      });
}

TEST_CASE("test_gas_aer_uptkrates_1box1gas", "mam4_gasaerexch_process") {
  mam4::AeroConfig mam4_config;
  mam4::GasAerExchImpl gasaerexch;
  gasaerexch.init(mam4_config);
  REQUIRE(std::string(gasaerexch.name()) == "MAM4 gas/aersol exchange");

  const std::string input_file = "skywalker_gasaerexch_uptkrates_1box1gas.yaml";
  const std::string output_file = "skywalker_gasaerexch_uptkrates_1box1gas.py";

  // Load the ensemble. Any error encountered is fatal.
  std::unique_ptr<Ensemble> ensemble(load_ensemble(input_file, "settings"));

  Settings settings = ensemble->settings();
  if (settings.has("name")) {
    const std::string name = settings.get("name");
    if ((name != "gas_aer_uptkrates_1box1gas")) {
      std::cerr << "Invalid name: " << name << std::endl;
      exit(0);
    }
  }
  // Run the ensemble.
  try {
    ensemble->process([&](const Input& input, Output& output) {
      // Ensemble parameters
      if (!input.has("temp")) {
        std::cerr << "Required name: "
                  << "temp" << std::endl;
        exit(0);
      }
      if (!input.has_array("dgncur_awet")) {
        std::cerr << "Required name: "
                  << "dgncur_awet" << std::endl;
        exit(0);
      }
      if (!input.has_array("lnsg")) {
        std::cerr << "Required name: "
                  << "lnsg" << std::endl;
        exit(0);
      }
      if (!input.has_array("aernum")) {
        std::cerr << "Required name: "
                  << "aernum" << std::endl;
        exit(0);
      }
      const bool has_solution = input.has_array("uptkaer");

      //-------------------------------------------------------
      // Process input, do calculations, and prepare output
      //-------------------------------------------------------
      const int n_mode = 4;
      const int nghq = 2;
      const Real accom = 0.65000000000000002;
      const Real beta_inp = 0.0000000000000000;
      const Real pi = 3.1415926535897931;
      const Real r_universal = 8314.4675910000005;
      const Real mw_gas = 98.078400000000002;
      const Real mw_air = 28.966000000000001;
      const Real pmid = 100000.00000000000;
      const Real pstd = 101325.00000000000;
      const Real vol_molar_gas = 42.880000000000003;
      const Real vol_molar_air = 20.100000000000001;

      const Kokkos::Array<bool, n_mode> l_condense_to_mode = {true, true, true,
                                                              true};

      // Parse input
      Kokkos::Array<PackType, n_mode> dgncur_awet;
      Kokkos::Array<Real, n_mode> lnsg;
      {
        const std::vector<Real> array = input.get_array("dgncur_awet");
        for (size_t i = 0; i < array.size() && i < n_mode; ++i)
          dgncur_awet[i] = array[i];
      }
      {
        const std::vector<Real> array = input.get_array("lnsg");
        for (size_t i = 0; i < array.size() && i < n_mode; ++i)
          lnsg[i] = array[i];
      }
      const std::vector<Real> aernum = input.get_array("aernum");
      const Real temp = input.get("temp");
      std::vector<Real> test_uptkaer;
      if (has_solution) {
        test_uptkaer = input.get_array("uptkaer");
      }

      kokkos_device_type::view_1d<PackType> uptkaer_dev("uptkaer on device",
                                                        n_mode);
      Kokkos::parallel_for(
          "gasaerexch.gas_aer_uptkrates_1box1gas", 1, KOKKOS_LAMBDA(const int) {
            Kokkos::Array<PackType, n_mode> uptkaer;
            gasaerexch.gas_aer_uptkrates_1box1gas(
                l_condense_to_mode, temp, pmid, pstd, mw_gas, mw_air,
                vol_molar_gas, vol_molar_air, accom, r_universal, pi, beta_inp,
                nghq, dgncur_awet, lnsg, uptkaer);
            for (size_t i = 0; i < n_mode; ++i) uptkaer_dev(i) = uptkaer[i];
          });
      Kokkos::Array<PackType, n_mode> uptkaer;
      {
        auto host_view = Kokkos::create_mirror_view(uptkaer_dev);
        Kokkos::deep_copy(host_view, uptkaer_dev);
        for (size_t i = 0; i < n_mode; ++i) uptkaer[i] = host_view[i];
      }
      // Write the computed nucleation rate.
      {
        std::vector<Real> values(n_mode);
        for (size_t i = 0; i < values.size() && i < n_mode; ++i)
          values[i] = uptkaer[i][0];
        output.set("uptkaer", values);
      }
    });
    // Write out a Python module.
    ensemble->write(output_file);
  } catch (Exception& e) {
    std::cerr << ": Error: " << e.what() << std::endl;
  }
}
