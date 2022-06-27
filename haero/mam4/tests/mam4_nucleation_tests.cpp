#include <haero/mam4/mam4.hpp>

#include <catch2/catch.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace haero;

TEST_CASE("test_constructor", "mam4_nucleation_process") {
  mam4::AeroConfig mam4_config;
  mam4::NucleationProcess process(mam4_config);
  REQUIRE(process.name() == "MAM4 nucleation");
  REQUIRE(process.config() == mam4_config);
}

TEST_CASE("test_compute_tendencies", "mam4_nucleation_process") {
  int nlev = 72;
  Real pblh = 1000;
  Atmosphere atm(nlev, pblh);
  mam4::Prognostics progs(nlev);
  mam4::Diagnostics diags(nlev);
  mam4::Tendencies tends(nlev);

  mam4::AeroConfig mam4_config;
  mam4::NucleationProcess process(mam4_config);

  // Single-column dispatch.
  auto team_policy = haero::TeamPolicy(1u, Kokkos::AUTO);
  Real t = 0.0, dt = 30.0;
  Kokkos::parallel_for(team_policy, KOKKOS_LAMBDA(const TeamType& team) {
    process.compute_tendencies(team, t, dt, atm, progs, diags, tends);
  });
}

