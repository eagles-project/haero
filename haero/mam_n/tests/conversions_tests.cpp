#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "ekat/logging/ekat_logger.hpp"
#include "haero/conversions.hpp"
#include "haero/floating_point.hpp"

using namespace haero;

TEST_CASE("conversions::mass_and_volume_mixing_ratios") {
  SECTION("mmr_from_number_conc") {
    Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
    PackType n(1e5);                        // number density
    PackType rho_air(1.1805);               // dry air mass density (20 C)
    auto mmr = conversions::mmr_from_number_conc(n, mu, rho_air);
    REQUIRE(FloatingPoint<PackType>::equiv(
        mmr, n * mu / (rho_air * Constants::avogadro)));
    REQUIRE(FloatingPoint<PackType>::equiv(
        n, conversions::number_conc_from_mmr(mmr, mu, rho_air)));
  }

  SECTION("number_conc_from_mmr") {
    Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
    PackType mmr(0.2);                      // mass mixing ratio
    PackType rho_air(1.1805);               // dry air mass density (20 C)
    auto n = conversions::number_conc_from_mmr(mmr, mu, rho_air);
    REQUIRE(FloatingPoint<PackType>::equiv(
        n, mmr * (rho_air * Constants::avogadro) / mu));
    REQUIRE(FloatingPoint<PackType>::equiv(
        mmr, conversions::mmr_from_number_conc(n, mu, rho_air)));
  }

  SECTION("mmr_from_vmr") {
    Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
    PackType vmr(0.2);                      // molar mixing ratio
    auto mmr = conversions::mmr_from_vmr(vmr, mu);
    REQUIRE(FloatingPoint<PackType>::equiv(
        mmr, vmr * mu / Constants::molec_weight_dry_air));
    REQUIRE(FloatingPoint<PackType>::equiv(vmr,
                                           conversions::vmr_from_mmr(mmr, mu)));
  }

  SECTION("vmr_from_mmr") {
    Real mu = Constants::molec_weight_so4;  // molecular weight of sulfate
    PackType mmr(0.2);                      // mass mixing ratio
    auto vmr = conversions::vmr_from_mmr(mmr, mu);
    REQUIRE(FloatingPoint<PackType>::equiv(
        vmr, mmr * Constants::molec_weight_dry_air / mu));
    REQUIRE(FloatingPoint<PackType>::equiv(mmr,
                                           conversions::mmr_from_vmr(vmr, mu)));
  }
}

TEST_CASE("conversions::water_vapor_functions") {
  ekat::Comm comm;

  ekat::logger::Logger<> logger("water_vapor_tests",
                                ekat::logger::LogLevel::debug, comm);

  const Real tol = (std::is_same<Real, double>::value ? 1e-9 : 1e-5);

  SECTION("mixing_ratio_vs_specific_humidity") {
    /*
  These tests are adapted from the Penn State METEO300 section 3.1 quiz,

  https://www.e-education.psu.edu/meteo300/node/519

  accessed 11/22/2021.
*/
    const PackType rho_vapor = 0.01;   // [kg / m3]
    const PackType rho_dry_air = 1.1;  // [kg / m3]

    // mixing ratio = kg water vapor / kg dry air
    const PackType qv = rho_vapor / rho_dry_air;
    // spec humidity = kg water vapor / kg moist air
    const PackType s = rho_vapor / (rho_dry_air + rho_vapor);
    logger.info(
        "water vapor density = {}, dry air density = {}; vapor mixing ratio qv "
        "= {}, specific humidity = {}",
        rho_vapor, rho_dry_air, qv, s);

    logger.debug("|qv - 0.00909091| = {}", abs(qv - 0.00909091));
    logger.debug("|s  - 0.00900901| = {}", abs(s - 0.00900901));

    REQUIRE(FloatingPoint<PackType>::equiv(qv, PackType(0.00909091), tol));
    REQUIRE(FloatingPoint<PackType>::equiv(s, PackType(0.00900901), tol));
  }
  SECTION("relative humidity") {
    const PackType qv = 0.021;  // kg vapor / kg dry air
    //     const PackType qvsat = 0.025; // kg vapor / kg dry air;
    const PackType P0 = 100000;   // Pa
    const PackType T0 = 293.486;  // K

    const auto rh =
        conversions::relative_humidity_from_vapor_mixing_ratio(qv, P0, T0);
    const auto es = conversions::vapor_saturation_pressure_magnus(T0);
    const auto qvsat = es / P0;

    logger.info("relative humidity with qv = 0.021 is {}", rh);
    logger.info("vapor saturation mixing ratio = {}", qvsat);
    logger.debug("|rh - 0.840019| = {}", abs(rh - 0.840019));
    logger.debug("|qvsat - 0.0249994| = {}", abs(qvsat - 0.0249994));

    REQUIRE(FloatingPoint<PackType>::equiv(rh, 0.840019, 2e-7));
    REQUIRE(FloatingPoint<PackType>::equiv(qvsat, 0.0249994, 3e-8));
  }
}
