// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include "atmosphere.hpp"

#include <ekat/ekat_assert.hpp>

namespace haero {

Atmosphere::Atmosphere(int num_levels, Real pblh)
    : temperature("temperature", num_levels), pressure("pressure", num_levels),
      vapor_mixing_ratio("vapor mixing ratio", num_levels),
      height("height", num_levels),
      hydrostatic_dp("hydrostatic_dp", num_levels),
      cloud_fraction("cloud_fraction", num_levels),
      updraft_vel_ice_nucleation("updraft_vel_ice_nucleation", num_levels),
      planetary_boundary_height(pblh), num_levels_(num_levels) {
  EKAT_REQUIRE_MSG(num_levels > 0,
                   "Number of vertical levels must be positive");
  EKAT_REQUIRE_MSG(pblh >= 0.0,
                   "Planetary boundary height must be non-negative");
}

Atmosphere::Atmosphere(int num_levels, const ColumnView temp,
                       const ColumnView press, const ColumnView qv,
                       const ColumnView ht, const ColumnView pdel,
                       const ColumnView cloud_f, const ColumnView uv_ice_nuc,
                       Real pblh)
    : temperature(temp), pressure(press), vapor_mixing_ratio(qv), height(ht),
      hydrostatic_dp(pdel), cloud_fraction(cloud_f),
      updraft_vel_ice_nucleation(uv_ice_nuc), planetary_boundary_height(pblh),
      num_levels_(num_levels) {
  EKAT_REQUIRE_MSG(num_levels > 0,
                   "Number of vertical levels must be positive");
  EKAT_REQUIRE_MSG(pblh >= 0.0,
                   "Planetary boundary height must be non-negative");

  // Make sure the views we're given are large enough to store their data.
  EKAT_REQUIRE_MSG(temp.extent(0) >= num_levels,
                   "Temperature view must have extent >= " << num_levels);
  EKAT_REQUIRE_MSG(press.extent(0) >= num_levels,
                   "Pressure view must have extent >= " << num_levels);
  EKAT_REQUIRE_MSG(
      qv.extent(0) >= num_levels,
      "Vapor mixing ratio view must have extent >= " << num_levels);
  EKAT_REQUIRE_MSG(
      pdel.extent(0) >= num_levels,
      "Hydrostatic pressure thickness must have extent >= " << num_levels);
  EKAT_REQUIRE_MSG(ht.extent(0) >= num_levels + 1,
                   "Height view must have extent >= " << num_levels);
  EKAT_REQUIRE_MSG(cloud_f.extent(0) >= num_levels + 1,
                   "cloud fraction must have extent >= " << num_levels);
  EKAT_REQUIRE_MSG(
      uv_ice_nuc.extent(0) >= num_levels + 1,
      "updraft velocity for ice nucleation must have extent >= " << num_levels);
}

} // namespace haero
