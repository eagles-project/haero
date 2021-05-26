#ifndef HAERO_MODAL_WET_RADIUS_DIAGNOSTIC_HPP
#define HAERO_MODAL_WET_RADIUS_DIAGNOSTIC_HPP

#include "haero/floating_point.hpp"
#include "haero/haero.hpp"
#include "haero/mode.hpp"
#include "kohler_solve.hpp"

namespace haero {

struct ModalWetRadius {
  static constexpr Real meters2microns = 1e3;
  static constexpr Real microns2meters = 1e-3;
  static constexpr Real dry_radius_min_microns =
      KohlerPolynomial<double>::dry_radius_min_microns;
  ColumnView wet_radius_meters;
  ConstColumnView modal_hygroscopicity;
  ConstColumnView modal_dry_radius_meters;
  ConstColumnView relative_humidity;
  Mode mode;
  Real hysteresis_fac;
  Real tol;

  KOKKOS_INLINE_FUNCTION
  ModalWetRadius(ColumnView rwet, const ColumnView hyg, const ColumnView rdry,
                 const ColumnView rh, const Mode& m)
      : wet_radius_meters(rwet),
        modal_hygroscopicity(hyg),
        modal_dry_radius_meters(rdry),
        relative_humidity(rh),
        mode(m),
        tol(100 * FloatingPoint<Real>::zero_tol) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int pack_idx) const {
    const Real dry_rmin = dry_radius_min_microns;
    const Real to_microns = meters2microns;
    const Real to_meters = microns2meters;

    // mask = true for all rh values below crystallization_pt
    const auto rh_low =
        (relative_humidity(pack_idx) <= mode.crystallization_pt);
    // mask = true for all values between crystallization_pt and
    // deliquescence_pt
    const auto rh_mid =
        (relative_humidity(pack_idx) > mode.crystallization_pt) &&
        (relative_humidity(pack_idx) <= mode.deliquescence_pt);
    // mask = true for all particles whose dry radius lies below the
    // KohlerPolynomial's minimum
    const auto too_small =
        (to_microns * modal_dry_radius_meters(pack_idx) < dry_rmin);
    const auto use_dry_radius = (rh_low || too_small);
    const auto needs_kohler = !use_dry_radius;
    const Real hysteresis_fac =
        1 / (mode.deliquescence_pt - mode.crystallization_pt);

    // below crystallization point, or for particles that are smaller than 1 nm,
    // there is no water uptake --- the particle radius is the dry radius
    ekat_masked_loop(use_dry_radius, s) {
      wet_radius_meters(pack_idx)[s] = modal_dry_radius_meters(pack_idx)[s];
    };

    // for all other rel. humidities, we need the Kohler equation to find the
    // wet radius
    ekat_masked_loop(needs_kohler, s) {
      const auto kpoly = KohlerPolynomial<double>(
          relative_humidity(pack_idx)[s], modal_hygroscopicity(pack_idx)[s],
          to_microns * modal_dry_radius_meters(pack_idx)[s]);
      auto solver = math::ScalarNewtonSolver<KohlerPolynomial<double>>(
          25 * to_microns * modal_dry_radius_meters(pack_idx)[s], tol, kpoly);
      wet_radius_meters(pack_idx)[s] = to_meters * solver.solve();
    };

    // for relative humidities between the crystallization and deliquescence
    // points, adjust wet radius due to hysteresis.
    ekat_masked_loop(rh_mid, s) {
      const Real dry_vol = mode.mean_particle_volume_from_diameter(
          2 * modal_dry_radius_meters(pack_idx)[s]);
      Real wet_vol = mode.mean_particle_volume_from_diameter(
          2 * wet_radius_meters(pack_idx)[s]);
      EKAT_KERNEL_ASSERT(wet_vol >= dry_vol);
      const Real water_vol =
          (wet_vol - dry_vol) *
          (relative_humidity(pack_idx)[s] - mode.crystallization_pt) *
          hysteresis_fac;
      EKAT_KERNEL_ASSERT(water_vol >= 0);
      wet_vol = dry_vol + water_vol;
      wet_radius_meters(pack_idx)[s] =
          0.5 * mode.mean_particle_diameter_from_volume(wet_vol);
    };
  }
};

}  // namespace haero
#endif
