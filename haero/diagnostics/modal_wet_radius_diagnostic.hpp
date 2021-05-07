#ifndef HAERO_MODAL_WET_RADIUS_DIAGNOSTIC_HPP
#define HAERO_MODAL_WET_RADIUS_DIAGNOSTIC_HPP

#include "haero/haero.hpp"
#include "haero/mode.hpp"
#include "kohler_solve_diagnostic.hpp"
#include "haero/floating_point.hpp"

namespace haero {

struct ModalWetRadius {
  ColumnView wet_radius_meters;
  ConstColumnView modal_hygroscopicity;
  ConstColumnView modal_dry_radius_meters;
  ConstColumnView relative_humidity;
  Real tol;

  KOKKOS_INLINE_FUNCTION
  ModalWetRadius(ColumnView rwet, const ColumnView hyg, const ColumnView rdry,
    const ColumnView rh) :
    wet_radius_meters(rwet),
    modal_hygroscopicity(hyg),
    modal_dry_radius_meters(rdry),
    relative_humidity(rh),
    tol(100*FloatingPoint<Real>::zero_tol)
    {}

  KOKKOS_INLINE_FUNCTION
  void operator() (const int pack_idx) const {
    const Real meters2microns = 1e3;
    const Real microns2meters = 1e-3;
    KohlerNewtonSolve kohler(relative_humidity(pack_idx), modal_hygroscopicity(pack_idx),
      meters2microns*modal_dry_radius_meters(pack_idx), tol);
      wet_radius_meters(pack_idx) = microns2meters * kohler();
  }
};

} // namespace haero
#endif
