// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include "atmosphere.hpp"

#include <ekat/ekat_assert.hpp>

namespace haero {

// this parameter should match the number of Views in the Atmosphere class
static constexpr size_t num_views = 7;

// this constructor allocates memory on-device for all of its views
Atmosphere::Atmosphere(int num_levels, Real pblh)
    : num_levels_(num_levels),
      view_storage_(reinterpret_cast<Real *>(
          Kokkos::kokkos_malloc(sizeof(Real) * num_views * num_levels))),
      temperature(&view_storage_[0], num_levels),
      pressure(&view_storage_[num_levels], num_levels),
      vapor_mixing_ratio(&view_storage_[2 * num_levels], num_levels),
      height(&view_storage_[3 * num_levels], num_levels),
      hydrostatic_dp(&view_storage_[4 * num_levels], num_levels),
      cloud_fraction(&view_storage_[5 * num_levels], num_levels),
      updraft_vel_ice_nucleation(&view_storage_[6 * num_levels], num_levels),
      planetary_boundary_height(pblh) {
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
    : num_levels_(num_levels), view_storage_(nullptr), temperature(temp),
      pressure(press), vapor_mixing_ratio(qv), height(ht), hydrostatic_dp(pdel),
      cloud_fraction(cloud_f), updraft_vel_ice_nucleation(uv_ice_nuc),
      planetary_boundary_height(pblh) {
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

Atmosphere::~Atmosphere() {
  if (view_storage_) {
    Kokkos::kokkos_free(view_storage_);
  }
}

} // namespace haero
