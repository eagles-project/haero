#include "catch2/catch.hpp"
#include <Kokkos_Core.hpp>
#include <cstdio>


TEST_CASE("kokkos_team_execute", "") {
typedef Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>::member_type TeamHandleType;
const auto& teamPolicy = Kokkos::TeamPolicy<Kokkos::DefaultExecutionSpace>(9u, Kokkos::AUTO);
Kokkos::parallel_for(teamPolicy,
                    KOKKOS_LAMBDA(const TeamHandleType& team)
                    {
                      printf("League Rank:%d\n",team.league_rank());
                      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, 0u, 3u), [&] (const int& i)
                      {
                        printf("Thread in Team:%d\n",i);
                      });
                    });
}

