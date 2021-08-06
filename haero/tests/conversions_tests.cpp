#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/conversions.hpp"
#include "haero/floating_point.hpp"

using namespace haero;

// static constexpr Real g_to_kg = 0.001;

TEST_CASE("conversions:mmr_from_number_conc") {
  Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
  PackType n(1e5);                        // number density
  PackType rho_air(1.1805);               // dry air mass density (20 C)
  auto mmr = conversions::mmr_from_number_conc(n, mu, rho_air);
  REQUIRE(FloatingPoint<PackType>::equiv(
      mmr, n * mu / (rho_air * Constants::avogadro)));
  REQUIRE(FloatingPoint<PackType>::equiv(
      n, conversions::number_conc_from_mmr(mmr, mu, rho_air)));
}

TEST_CASE("conversions:number_conc_from_mmr") {
  Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
  PackType mmr(0.2);                      // mass mixing ratio
  PackType rho_air(1.1805);               // dry air mass density (20 C)
  auto n = conversions::number_conc_from_mmr(mmr, mu, rho_air);
  REQUIRE(FloatingPoint<PackType>::equiv(
      n, mmr * (rho_air * Constants::avogadro) / mu));
  REQUIRE(FloatingPoint<PackType>::equiv(
      mmr, conversions::mmr_from_number_conc(n, mu, rho_air)));
}

TEST_CASE("conversions:mmr_from_vmr") {
  Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
  PackType vmr(0.2);                      // molar mixing ratio
  auto mmr = conversions::mmr_from_vmr(vmr, mu);
  REQUIRE(FloatingPoint<PackType>::equiv(
      mmr, vmr * mu / Constants::molec_weight_dry_air));
  REQUIRE(
      FloatingPoint<PackType>::equiv(vmr, conversions::vmr_from_mmr(mmr, mu)));
}

TEST_CASE("conversions:vmr_from_mmr") {
  Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
  PackType mmr(0.2);                      // mass mixing ratio
  auto vmr = conversions::vmr_from_mmr(mmr, mu);
  REQUIRE(FloatingPoint<PackType>::equiv(
      vmr, mmr * Constants::molec_weight_dry_air / mu));
  REQUIRE(
      FloatingPoint<PackType>::equiv(mmr, conversions::mmr_from_vmr(vmr, mu)));
}
