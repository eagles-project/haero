#include "haero/haero.hpp"
#include "haero/diagnostics/kohler_solve.hpp"
#include "haero/math_helpers.hpp"
#include "haero/utils.hpp"
#include "catch2/catch.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace haero;

std::string mathematica_verification_program();

struct KohlerTestInput {
  DeviceType::view_1d<PackType> relative_humidity;
  DeviceType::view_1d<PackType> hygroscopicity;
  DeviceType::view_1d<PackType> dry_radius;

  typename DeviceType::view_1d<PackType>::HostMirror h_relative_humidity;
  typename DeviceType::view_1d<PackType>::HostMirror h_hygroscopicity;
  typename DeviceType::view_1d<PackType>::HostMirror h_dry_radius;

  KohlerTestInput(const int n) :
    relative_humidity("relative_humidity", cube(n)),
    hygroscopicity("hygroscopicity", cube(n)),
    dry_radius("dry_radius", cube(n)) {
      EKAT_REQUIRE(n>1);
      const Real drelh = (KohlerPolynomial<PackType>::rel_humidity_max -
        KohlerPolynomial<PackType>::rel_humidity_min)/(n-1);
      const Real dhyg = (KohlerPolynomial<PackType>::hygro_max -
        KohlerPolynomial<PackType>::hygro_min)/(n-1);
      const Real drad = (KohlerPolynomial<PackType>::dry_radius_max_microns -
        KohlerPolynomial<PackType>::dry_radius_min_microns)/(n-1);

      h_relative_humidity = Kokkos::create_mirror_view(relative_humidity);
      h_hygroscopicity = Kokkos::create_mirror_view(hygroscopicity);
      h_dry_radius = Kokkos::create_mirror_view(dry_radius);

      int ind=0;
      for (int i=0; i<n; ++i) {
        const Real rel_h = KohlerPolynomial<PackType>::rel_humidity_min + i * drelh;
        for (int j=0; j<n; ++j) {
          const Real hyg = KohlerPolynomial<PackType>::hygro_min + j * dhyg;
          for (int k=0; k<n; ++k) {
            const Real drad = KohlerPolynomial<PackType>::dry_radius_min_microns + k * drad;
            h_relative_humidity(ind) = rel_h;
            h_hygroscopicity(ind) = hyg;
            h_dry_radius(ind++) = drad;
          }
        }
      }

      Kokkos::deep_copy(relative_humidity, h_relative_humidity);
      Kokkos::deep_copy(hygroscopicity, h_hygroscopicity);
      Kokkos::deep_copy(dry_radius, h_dry_radius);
    }
};

TEST_CASE("KohlerSolve-Newton", "") {

}

TEST_CASE("KohlerSolve-Bisection", "") {
}

TEST_CASE("KohlerSolve-verification", "") {
  /// Generate input data
  static constexpr int N = 20;


  std::fstream mm_sols(HAERO_TEST_DATA_DIR + "/mm_kohler_roots.txt");
  REQUIRE(mm_sols.is_open());

  mm_sols.close();

}



