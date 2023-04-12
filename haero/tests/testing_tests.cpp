// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include <haero/testing.hpp>

#include <catch2/catch.hpp>

using namespace haero;

TEST_CASE("create_column_view", "") {
  ColumnView view = testing::create_column_view(72);
  REQUIRE(view.extent(0) == 72);
}

TEST_CASE("create_atmosphere", "") {
  Atmosphere atm = testing::create_atmosphere(72, 100.0);
  REQUIRE(atm.temperature.extent(0) == 72);
  REQUIRE(atm.planetary_boundary_layer_height == 100.0);
}
