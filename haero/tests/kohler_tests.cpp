#include <cmath>
#include <iostream>

#include "catch2/catch.hpp"
#include "haero/haero.hpp"
#include "haero/kohler.hpp"
#include "haero/math.hpp"
#include "haero/constants.hpp"
#include "haero/floating_point.hpp"
#include "ekat/ekat_pack_math.hpp"

using namespace haero;

/**
  This struct builds an array of n input values for each
  parameter in the Kohler polynomial for use on device or host,
  for a total of n**3 trials.
*/
struct KohlerTestInput {
  DeviceType::view_1d<PackType> relative_humidity;
  DeviceType::view_1d<PackType> hygroscopicity;
  DeviceType::view_1d<PackType> dry_radius;
  static constexpr Real temperature = Constants::triple_pt_h2o;

  typename DeviceType::view_1d<PackType>::HostMirror h_relative_humidity;
  typename DeviceType::view_1d<PackType>::HostMirror h_hygroscopicity;
  typename DeviceType::view_1d<PackType>::HostMirror h_dry_radius;

  explicit KohlerTestInput(const int n) :
    relative_humidity("relative_humidity", PackInfo::num_packs(cube(n))),
    hygroscopicity("hygroscopicity", PackInfo::num_packs(cube(n))),
    dry_radius("dry_radius", PackInfo::num_packs(cube(n))) {
      EKAT_REQUIRE(n > 1);
      const Real drelh = (KohlerPolynomial<Real>::rel_humidity_max -
                          KohlerPolynomial<Real>::rel_humidity_min) /
                          (n - 1);
      const Real dhyg =  (KohlerPolynomial<Real>::hygro_max -
                          KohlerPolynomial<Real>::hygro_min) /
                          (n - 1);
      const Real ddry =  (KohlerPolynomial<Real>::dry_radius_max_microns -
                          KohlerPolynomial<Real>::dry_radius_min_microns) /
                          (n - 1);

      h_relative_humidity = Kokkos::create_mirror_view(relative_humidity);
      h_hygroscopicity = Kokkos::create_mirror_view(hygroscopicity);
      h_dry_radius = Kokkos::create_mirror_view(dry_radius);

      int ind = 0;
      for (int i=0; i<n; ++i) {
        const Real rel_h = KohlerPolynomial<Real>::rel_humidity_min + i*drelh;
        for (int j=0; j<n; ++j) {
          const Real hyg = KohlerPolynomial<Real>::hygro_min + j*dhyg;
          for (int k=0; k<n; ++k) {
            const Real dry_rad =
              KohlerPolynomial<Real>::dry_radius_min_microns + k*ddry;
            const int pack_idx = PackInfo::pack_idx(ind);
            const int vec_idx = PackInfo::vec_idx(ind++);
            h_relative_humidity(pack_idx)[vec_idx] = rel_h;
            h_hygroscopicity(pack_idx)[vec_idx] = hyg;
            h_dry_radius(pack_idx)[vec_idx] = dry_rad;
          }
        }
      }
      Kokkos::deep_copy(relative_humidity, h_relative_humidity);
      Kokkos::deep_copy(hygroscopicity, h_hygroscopicity);
      Kokkos::deep_copy(dry_radius, h_dry_radius);
    }
};
