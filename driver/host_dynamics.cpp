#include "host_dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "ekat/ekat_assert.hpp"
#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/util/ekat_units.hpp"
#include "haero/floating_point.hpp"
#include "haero/physical_constants.hpp"
#include "haero/utils.hpp"
#include "ncwriter_impl.hpp"

namespace haero {
namespace driver {

void HostDynamics::init_from_interface_heights(std::vector<Real> z0,
                                               AtmosphericConditions& ac) {
  using namespace constants;

  EKAT_REQUIRE_MSG(z0.size() == nlev_ + 1,
                   "number of initial heights must match number of levels + 1");
  EKAT_REQUIRE_MSG(vector_is_monotone(z0),
                   "initial heights must be a monotone array.");

  const bool increasing = (z0[1] > z0[0]);
  if (increasing) std::reverse(z0.begin(), z0.end());

  ac.ztop = z0[0];
  ac.ptop = hydrostatic_pressure_at_height(ac.ztop, ac);

  auto hphi0 = Kokkos::create_mirror_view(phi0);
  auto hrho0 = Kokkos::create_mirror_view(rho0);
  auto hp = Kokkos::create_mirror_view(p);
  auto hthetav = Kokkos::create_mirror_view(thetav);
  auto hqv = Kokkos::create_mirror_view(qv);
  auto hw = Kokkos::create_mirror_view(w);
  auto hpint = Kokkos::create_mirror_view(phydro_int);

  /// set interface geopotential
  for (int k = 0; k < nlev_ + 1; ++k) {
    // Taylor et al. 2020 fig. 1 interface idx = k+1/2
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);
    hphi0(pack_idx)[vec_idx] = gravity * z0[k];
    hw(pack_idx)[vec_idx] = 0;
    hpint(pack_idx)[vec_idx] = hydrostatic_pressure_at_height(z0[k], ac);
  }

  /// set midpoint pressure, density, virtual potential temperature, water vapor
  /// mixing ratio
  for (int k = 0; k < nlev_; ++k) {
    // Taylor et al. 2020 fig. 1 level idx = k+1
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);

    const int kphalf_idx = k + 1;  // array idx of interface k + 1/2
    const int kmhalf_idx = k;      // array idx of interface k - 1/2
    const auto kphalf_pack_idx = PackInfo::pack_idx(kphalf_idx);
    const auto kphalf_vec_idx = PackInfo::vec_idx(kphalf_idx);
    const auto kmhalf_pack_idx = PackInfo::pack_idx(kmhalf_idx);
    const auto kmhalf_vec_idx = PackInfo::vec_idx(kmhalf_idx);

    const Real phimid = 0.5 * (hphi0(kmhalf_pack_idx)[kmhalf_vec_idx] +
                               hphi0(kphalf_pack_idx)[kphalf_vec_idx]);
    const Real zmid = phimid / gravity;
    const Real pres = hydrostatic_pressure_at_height(zmid, ac);
    const Real Tv = ac.Tv0 - ac.Gammav * zmid;
    hp(pack_idx)[vec_idx] = pres;
    hrho0(pack_idx)[vec_idx] = pres / (r_gas_dry_air * Tv);
    hthetav(pack_idx)[vec_idx] = Tv / exner_function(pres);
    hqv(pack_idx)[vec_idx] = water_vapor_mixing_ratio(zmid, ac);
  }

  rho0surf = AtmosphericConditions::pref / (r_gas_dry_air * ac.Tv0);
  ps = hydrostatic_pressure_at_height(0, ac);
  EKAT_ASSERT_MSG(
      FloatingPoint<Real>::rel(ps, AtmosphericConditions::pref),
      "surface pressure must equal the reference pressure, 1000 hPa.");

  Kokkos::deep_copy(w, hw);
  Kokkos::deep_copy(phi0, hphi0);
  Kokkos::deep_copy(phi, phi0);
  Kokkos::deep_copy(rho0, hrho0);
  Kokkos::deep_copy(rho, rho0);
  Kokkos::deep_copy(thetav, hthetav);
  Kokkos::deep_copy(qv, hqv);
  Kokkos::deep_copy(p, hp);

  Kokkos::deep_copy(phydro_int, hpint);

  update_thickness(ac);
}

void HostDynamics::init_from_uniform_heights(const AtmosphericConditions& ac) {
  using namespace constants;

  const Real dz = ac.ztop / nlev_;

  /// set interface geopotential and velocity
  // host mirrors preceded by "h"
  auto hphi0 = Kokkos::create_mirror_view(phi0);
  auto hw = Kokkos::create_mirror_view(w);
  auto hpint = Kokkos::create_mirror_view(phydro_int);
  for (int k = 0; k < nlev_ + 1; ++k) {
    const int pack_idx = PackInfo::pack_idx(k);
    const int vec_idx = PackInfo::vec_idx(k);
    // Taylor et al. 2020 fig. 1 interface idx = k+1/2
    const Real z = ac.ztop - k * dz;
    hphi0(pack_idx)[vec_idx] = gravity * z;
    hw(pack_idx)[vec_idx] = 0;
    hpint(pack_idx)[vec_idx] = hydrostatic_pressure_at_height(z, ac);
  }
  Kokkos::deep_copy(w, hw);
  Kokkos::deep_copy(phi0, hphi0);
  Kokkos::deep_copy(phi, phi0);
  Kokkos::deep_copy(phydro_int, hpint);

  /// set midpoint pressure, density, virtual potential temperature, water vapor
  /// mixing ratio
  // host mirrors preceded by "h"
  auto hp = Kokkos::create_mirror_view(p);
  auto hthetav = Kokkos::create_mirror_view(thetav);
  auto hqv = Kokkos::create_mirror_view(qv);
  auto hrho0 = Kokkos::create_mirror_view(rho0);
  for (int k = 0; k < nlev_; ++k) {
    // Taylor et al. 2020 fig. 1 level idx = k+1
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);

    const int kphalf_idx = k + 1;  // array idx of interface k + 1/2
    const int kmhalf_idx = k;      // array idx of interface k - 1/2
    const auto kphalf_pack_idx = PackInfo::pack_idx(kphalf_idx);
    const auto kphalf_vec_idx = PackInfo::vec_idx(kphalf_idx);
    const auto kmhalf_pack_idx = PackInfo::pack_idx(kmhalf_idx);
    const auto kmhalf_vec_idx = PackInfo::vec_idx(kmhalf_idx);

    const Real phimid = 0.5 * (hphi0(kmhalf_pack_idx)[kmhalf_vec_idx] +
                               hphi0(kphalf_pack_idx)[kphalf_vec_idx]);
    const Real zmid = phimid / gravity;
    const Real pres = hydrostatic_pressure_at_height(zmid, ac);
    const Real Tv = ac.Tv0 - ac.Gammav * zmid;
    hp(pack_idx)[vec_idx] = pres;
    hrho0(pack_idx)[vec_idx] = pres / (r_gas_dry_air * Tv);
    hthetav(pack_idx)[vec_idx] = Tv / exner_function(pres);
    hqv(pack_idx)[vec_idx] = water_vapor_mixing_ratio(zmid, ac);
  }

  rho0surf = AtmosphericConditions::pref / (r_gas_dry_air * ac.Tv0);
  ps = hydrostatic_pressure_at_height(0, ac);
  EKAT_ASSERT_MSG(
      FloatingPoint<Real>::rel(ps, AtmosphericConditions::pref),
      "surface pressure must equal the reference pressure, 1000 hPa.");

  Kokkos::deep_copy(rho0, hrho0);
  Kokkos::deep_copy(rho, rho0);
  Kokkos::deep_copy(thetav, hthetav);
  Kokkos::deep_copy(qv, hqv);
  Kokkos::deep_copy(p, hp);

  update_thickness(ac);
}

void HostDynamics::init_from_interface_pressures(std::vector<Real> p0,
                                                 AtmosphericConditions& ac) {
  using namespace constants;

  EKAT_REQUIRE_MSG(
      p0.size() == nlev_ + 1,
      "number of initial pressures must match number of levels + 1");
  EKAT_REQUIRE_MSG(vector_is_monotone(p0),
                   "initial pressures must be a monotone array.");

  const bool increasing = (p0[1] > p0[0]);
  if (!increasing) std::reverse(p0.begin(), p0.end());
  EKAT_REQUIRE_MSG(p0.back() == AtmosphericConditions::pref,
                   "surface pressure must initialize to 1000 hPa.");

  ac.ptop = p0[0];
  ac.ztop = height_at_pressure(p0[0], ac);
  // host mirrors preceded by "h"
  auto hphi0 = Kokkos::create_mirror_view(phi0);
  auto hrho0 = Kokkos::create_mirror_view(rho0);
  auto hp = Kokkos::create_mirror_view(p);
  auto hthetav = Kokkos::create_mirror_view(thetav);
  auto hqv = Kokkos::create_mirror_view(qv);
  auto hw = Kokkos::create_mirror_view(w);
  auto hpint = Kokkos::create_mirror_view(phydro_int);

  /// set interface geopotential
  for (int k = 0; k < nlev_ + 1; ++k) {
    // Taylor et al. 2020 fig. 1 interface idx = k+1/2
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);

    hw(pack_idx)[vec_idx] = 0;
    const Real z = height_at_pressure(p0[k], ac);
    hphi0(pack_idx)[vec_idx] = gravity * z;
    const Real phydro = hydrostatic_pressure_at_height(z, ac);
    hpint(pack_idx)[vec_idx] = phydro;

    // assert that p(z) = p(z(p))
    //     if (!FloatingPoint<Real>::equiv(z, height_at_pressure(phydro, ac),
    //     35*FloatingPoint<Real>::zero_tol)) {
    //       std::cout << " z = " << z << " z(p) = " <<
    //       height_at_pressure(phydro, ac) << " |diff| = "
    //         << std::abs(z - height_at_pressure(phydro, ac)) << "\n";
    //     }
    EKAT_ASSERT(
        FloatingPoint<Real>::equiv(z, height_at_pressure(phydro, ac),
                                   (std::is_same<float, Real>::value
                                        ? 2e-3
                                        : 35 * FloatingPoint<Real>::zero_tol)));
  }

  /// set midpoint pressure, density, virtual potential temperature, water vapor
  /// mixing ratio
  for (int k = 0; k < nlev_; ++k) {
    // Taylor et al. 2020 fig. 1 level idx = k+1
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);

    const int kphalf_idx = k + 1;  // array idx of interface k + 1/2
    const int kmhalf_idx = k;      // array idx of interface k - 1/2
    const auto kphalf_pack_idx = PackInfo::pack_idx(kphalf_idx);
    const auto kphalf_vec_idx = PackInfo::vec_idx(kphalf_idx);
    const auto kmhalf_pack_idx = PackInfo::pack_idx(kmhalf_idx);
    const auto kmhalf_vec_idx = PackInfo::vec_idx(kmhalf_idx);

    const Real phimid = 0.5 * (hphi0(kmhalf_pack_idx)[kmhalf_vec_idx] +
                               hphi0(kphalf_pack_idx)[kphalf_vec_idx]);
    const Real zmid = phimid / gravity;
    const Real pres = hydrostatic_pressure_at_height(zmid, ac);
    const Real Tv = ac.Tv0 - ac.Gammav * zmid;
    hp(pack_idx)[vec_idx] = pres;
    hrho0(pack_idx)[vec_idx] = pres / (r_gas_dry_air * Tv);
    hthetav(pack_idx)[vec_idx] = Tv / exner_function(pres);
    hqv(pack_idx)[vec_idx] = water_vapor_mixing_ratio(zmid, ac);
  }

  ps = p0.back();
  EKAT_ASSERT(FloatingPoint<Real>::rel(
      ps, hpint(PackInfo::last_pack_idx(
              nlev_ + 1))[PackInfo::last_vec_end(nlev_ + 1) - 1]));
  rho0surf = AtmosphericConditions::pref / (r_gas_dry_air * ac.Tv0);

  Kokkos::deep_copy(w, hw);
  Kokkos::deep_copy(phi0, hphi0);
  Kokkos::deep_copy(phi, hphi0);
  Kokkos::deep_copy(rho0, hrho0);
  Kokkos::deep_copy(rho, hrho0);
  Kokkos::deep_copy(thetav, hthetav);
  Kokkos::deep_copy(qv, hqv);
  Kokkos::deep_copy(p, hp);
  rho0surf = AtmosphericConditions::pref / (r_gas_dry_air * ac.Tv0);

  update_thickness(ac);
}

void HostDynamics::init_from_uniform_pressures(
    const AtmosphericConditions& ac) {
  using namespace constants;

  const Real delp = (AtmosphericConditions::pref - ac.ptop) / nlev_;
  std::cout << "unif. pressure dp = " << delp << "\n";

  // host mirrors preceded by "h"
  auto hphi0 = Kokkos::create_mirror_view(phi0);
  auto hw = Kokkos::create_mirror_view(w);
  auto hp = Kokkos::create_mirror_view(p);
  auto hrho0 = Kokkos::create_mirror_view(rho0);
  auto hthetav = Kokkos::create_mirror_view(thetav);
  auto hqv = Kokkos::create_mirror_view(qv);
  auto hdz = Kokkos::create_mirror_view(dz);
  auto hpint = Kokkos::create_mirror_view(phydro_int);

  // scalarized views preceded by "s"
  auto shphi0 = ekat::scalarize(hphi0);
  auto shw = ekat::scalarize(hw);
  auto shp = ekat::scalarize(hp);
  auto shrho0 = ekat::scalarize(hrho0);
  auto shthetav = ekat::scalarize(hthetav);
  auto shqv = ekat::scalarize(hqv);
  auto shdz = ekat::scalarize(hdz);
  auto shpint = ekat::scalarize(hpint);

  //   / set interface geopotential and velocity
  for (int k = 0; k < nlev_ + 1; ++k) {
    //     Taylor et al. 2020 fig. 1 interface idx = k+1/2
    const Real punif = ac.ptop + k * delp;
    const Real z = height_at_pressure(punif, ac);
    shphi0(k) = gravity * z;
    shw(k) = 0;
    shpint(k) = punif;
  }

  EKAT_ASSERT(FloatingPoint<Real>::equiv(AtmosphericConditions::pref,
                                         ac.ptop + nlev_ * delp));

  for (int k = 0; k < nlev_; ++k) {
    // Taylor et al. 2020 fig. 1 level idx = k
    const int kmhalf_idx = k;      // array idx of interface k - 1/2
    const int kphalf_idx = k + 1;  // array idx of interface k + 1/2

    const Real phimid = 0.5 * (shphi0(kmhalf_idx) + shphi0(kphalf_idx));
    const Real zmid = phimid / gravity;
    const Real pres = hydrostatic_pressure_at_height(zmid, ac);
    const Real Tv = ac.Tv0 - ac.Gammav * zmid;

    shp(k) = pres;
    shrho0(k) = pres / (r_gas_dry_air * Tv);
    shthetav(k) = Tv / exner_function(pres);
    shqv(k) = water_vapor_mixing_ratio(zmid, ac);
    shdz(k) = (shphi0(kmhalf_idx) - shphi0(kphalf_idx)) / gravity;
  }

  Kokkos::deep_copy(w, hw);
  Kokkos::deep_copy(phi0, hphi0);
  Kokkos::deep_copy(phi, phi0);
  Kokkos::deep_copy(qv, hqv);
  Kokkos::deep_copy(thetav, hthetav);
  Kokkos::deep_copy(rho0, hrho0);
  Kokkos::deep_copy(rho, rho0);
  Kokkos::deep_copy(p, hp);
  Kokkos::deep_copy(dz, hdz);
  Kokkos::deep_copy(phydro_int, hpint);

  ps = AtmosphericConditions::pref;
  rho0surf = AtmosphericConditions::pref / (r_gas_dry_air * ac.Tv0);

  update_thickness(ac);
}

void HostDynamics::update_thickness(const AtmosphericConditions& conds) {
  using namespace constants;

  Kokkos::parallel_for("HostDynamics::hydrostatic_pint", nlev_ + 1,
                       HydrostaticPressureUpdate(phi, phydro_int, conds));

  Kokkos::fence();

  // differences at level midpoints

  auto dz_local = ekat::scalarize(dz);
  auto dph = ekat::scalarize(hydrostatic_dp);
  auto phi_local = ekat::scalarize(phi);
  auto pint_local = ekat::scalarize(phydro_int);

  Kokkos::parallel_for(
      "HostDynamics::dz", nlev_, KOKKOS_LAMBDA(const int k) {
        const int kmhalf_idx = k;
        const int kphalf_idx = k + 1;
        // negative sign because levels & interfaces are indexed from model top
        // to surface, but z increases from surface to top
        dz_local(k) =
            -(phi_local(kphalf_idx) - phi_local(kmhalf_idx)) / gravity;
        dph(k) = -(pint_local(kphalf_idx) - pint_local(kmhalf_idx));
      });

  Kokkos::fence();
}

std::string HostDynamics::info_string(int tab_level) const {
  std::string tabstr = indent_string(tab_level);
  std::ostringstream ss;
  ss << tabstr << "HostDynamics info:\n";
  tabstr += "\t";
  ss << "nlev = " << nlev_ << "\n";
  ss << "ps = " << ps << "\n";
  ss << "rho0surf = " << rho0surf << "\n";
  ss << "t = " << t << "\n";
  return ss.str();
}

void HostDynamics::update(const Real newt, const AtmosphericConditions& ac) {
  using namespace constants;
  // interface update
  Kokkos::parallel_for("HostDynamics::InterfaceUpdate", nlev_ + 1,
                       DynamicsInterfaceUpdate(newt, phi, phi0, w, ac, nlev_));
  Kokkos::fence();

  const auto phi_local = ekat::scalarize(phi);
  const auto phi0_local = ekat::scalarize(phi0);
  auto rho_local = ekat::scalarize(rho);
  auto rho0_local = ekat::scalarize(rho0);
  auto thetav_local = ekat::scalarize(thetav);
  auto p_local = ekat::scalarize(p);
  Kokkos::parallel_for(
      "HostDynamics::MidpointUpdate", nlev_, KOKKOS_LAMBDA(const int k) {
        const Real phimid = 0.5 * (phi_local(k + 1) + phi_local(k));
        const Real phi0mid = 0.5 * (phi0_local(k + 1) + phi0_local(k));
        rho_local(k) = density(newt, phimid, phi0mid, rho0_local(k), ac);
        p_local(k) = pressure(rho_local(k), thetav_local(k));
      });

  Kokkos::fence();

  const Real rhosurf = density(newt, 0, 0, rho0surf, ac);
  ps = pressure(rhosurf, ac.Tv0);
  update_thickness(ac);
}

void HostDynamics::nc_init_dynamics_variables(
    NcWriter& writer, const AtmosphericConditions& conds) const {
  using att_type = NcWriter::text_att_type;
  using var_atts = std::vector<att_type>;

  // Add level and interface dimensions to file, if necessary
  if (!writer.levels_defined()) {
    writer.add_level_dims(nlev_);
  }

  // Interface variables
  const var_atts w_atts = {
      std::make_pair("cf_long_name", "upward_air_velocity"),
      std::make_pair("short_name", "w")};
  const auto w_units = ekat::units::m / ekat::units::s;

  const var_atts phi_atts = {std::make_pair("cf_long_name", "geopotential"),
                             std::make_pair("short_name", "phi")};
  const auto phi_units = ekat::units::pow(ekat::units::m, 2) *
                         ekat::units::pow(ekat::units::s, -2);

  const var_atts ph_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name",
                     "hydrostatic_air_pressure_at_interfaces"),
      std::make_pair("short_name", "phydro_int")};
  const auto ph_units = ekat::units::Pa;

  writer.define_interface_var("vertical_velocity", w_units, w, w_atts);
  writer.define_interface_var("geopotential", phi_units, phi, phi_atts);
  writer.define_interface_var("phydro_int", ph_units, phydro_int, ph_atts);

  // level variables
  const var_atts thetav_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "virtual_potential_temperature"),
      std::make_pair("short_name", "theta_v")};
  const auto thetav_units = ekat::units::K;

  const var_atts p_atts = {std::make_pair("cf_long_name", "air_pressure"),
                           std::make_pair("amip_short_name", "plev"),
                           std::make_pair("short_name", "p")};
  const auto p_units = ekat::units::Pa;

  const var_atts qv_atts = {
      std::make_pair("cf_long_name", "humidity_mixing_ratio"),
      std::make_pair("haero_long_name", "water_vapor_mass_mixing_ratio"),
      std::make_pair("short_name", "q_v")};
  const auto qv_units = ekat::units::kg / ekat::units::kg;

  const var_atts rho_atts = {std::make_pair("cf_long_name", "air_density"),
                             std::make_pair("short_name", "rho")};
  const auto rho_units = ekat::units::kg * ekat::units::pow(ekat::units::m, -3);

  const var_atts dz_atts = {
      std::make_pair("cf_long_name",
                     "atmosphere_layer_thickness_expressed_as_geopotential_"
                     "height_difference"),
      std::make_pair("short_name", "dz")};
  const auto dz_units = ekat::units::m;

  const var_atts dp_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("short_name", "pdel"),
      std::make_pair("haero_long_name",
                     "atmosphere_layer_thickness_expressed_as_hydrostatic_"
                     "pressure_difference"),
      std::make_pair("scream_name", "pseudo_density")};
  const auto dp_units = ekat::units::Pa / ekat::units::m;

  writer.define_level_var("density", rho_units, rho, rho_atts);
  writer.define_level_var("thetav", thetav_units, thetav, thetav_atts);
  writer.define_level_var("qv", qv_units, qv, qv_atts);
  writer.define_level_var("p", p_units, p, p_atts);
  writer.define_level_var("dz", dz_units, dz, dz_atts);
  writer.define_level_var("hydrostatic_dp", dp_units, hydrostatic_dp, dp_atts);

  // surface variables
  const var_atts ps_atts = {
      std::make_pair("cf_long_name", "surface_air_pressure"),
      std::make_pair("short_name", "psurf"),
      std::make_pair("amip_short_name", "ps")};
  const auto ps_units = ekat::units::Pa;
  writer.define_time_dependent_scalar_var("surface_pressure", ps_units,
                                          ps_atts);

  // scalar variables and parameters
  const var_atts p0_atts = {
      std::make_pair(
          "cf_long_name",
          "reference_air_pressure_for_atmospheric_vertical_coordinate"),
      std::make_pair("short_name", "pref")};
  writer.define_scalar_var("p0", ekat::units::Pa, p0_atts,
                           AtmosphericConditions::pref);

  const var_atts Tv0_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "reference_virtual_temperature")};
  const auto Tv0_units = ekat::units::K;
  writer.define_scalar_var("Tv0", Tv0_units, Tv0_atts, conds.Tv0);

  const var_atts Gammav_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name",
                     "initial_virtual_temperature_lapse_rate"),
      std::make_pair("short_name", "Gamma_v")};
  const auto Gammav_units = ekat::units::K / ekat::units::m;
  writer.define_scalar_var("Tv_lapse_rate", Gammav_units, Gammav_atts,
                           conds.Gammav);

  const var_atts qv0_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name",
                     "initial_surface_water_vapor_mixing_ratio"),
      std::make_pair("short_name", "qv0")};
  writer.define_scalar_var("qv0", ekat::units::Units::nondimensional(),
                           qv0_atts, conds.qv0);

  const var_atts qv1_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name",
                     "initial_water_vapor_mixing_ratio_decay_rate"),
      std::make_pair("short_name", "qv1")};
  writer.define_scalar_var("qv1", ekat::units::pow(ekat::units::m, -1),
                           qv1_atts, conds.qv1);

  const var_atts ztop_atts = {
      std::make_pair("cf_long_name", "altitude_at_top_of_atmosphere_model"),
      std::make_pair("short_name", "z_top")};
  writer.define_scalar_var("ztop", ekat::units::m, ztop_atts, conds.ztop);

  const var_atts ptop_atts = {
      std::make_pair("cf_long_name", "air_pressure_at_top_of_atmosphere_model"),
      std::make_pair("short_name", "p_top")};
  writer.define_scalar_var("ptop", ekat::units::Pa, ptop_atts, conds.ptop);

  const var_atts w0_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "host_dynamics_maximum_velocity"),
      std::make_pair("short_name", "w0")};
  writer.define_scalar_var("w0", ekat::units::m / ekat::units::s, w0_atts,
                           conds.w0);

  const var_atts tp_atts = {
      std::make_pair("cf_long_name", "null"),
      std::make_pair("haero_long_name", "host_dynamics_oscillation_period"),
      std::make_pair("short_name", "t_period")};
  writer.define_scalar_var("tperiod", ekat::units::s, tp_atts, conds.tperiod);
}

void HostDynamics::nc_write_data(NcWriter& writer,
                                 const size_t time_idx) const {
  const int null_idx = -1;

  writer.add_variable_data("vertical_velocity", time_idx, null_idx, null_idx,
                           w);
  writer.add_variable_data("geopotential", time_idx, null_idx, null_idx, phi);
  writer.add_variable_data("phydro_int", time_idx, null_idx, null_idx,
                           phydro_int);
  writer.add_variable_data("density", time_idx, null_idx, null_idx, rho);
  writer.add_variable_data("thetav", time_idx, null_idx, null_idx, thetav);
  writer.add_variable_data("qv", time_idx, null_idx, null_idx, qv);
  writer.add_variable_data("p", time_idx, null_idx, null_idx, p);
  writer.add_variable_data("hydrostatic_dp", time_idx, null_idx, null_idx,
                           hydrostatic_dp);
  writer.add_variable_data("dz", time_idx, null_idx, null_idx, dz);
  writer.add_time_dependent_scalar_value("surface_pressure", time_idx, ps);
}

Atmosphere HostDynamics::create_atmospheric_state(ColumnView temp,
                                                  ColumnView relh,
                                                  ColumnView z) const {
  using namespace constants;
  const auto p_local = p;
  const auto phi_local = phi;
  const auto thetav_local = thetav;
  const auto qv_local = qv;
  Kokkos::parallel_for(
      "HostDynamics:CreateAtmosphereLevels", nlev_, KOKKOS_LAMBDA(const int k) {
        const int pack_idx = PackInfo::pack_idx(k);
        const int vec_idx = PackInfo::vec_idx(k);

        const Real P = p_local(pack_idx)[vec_idx];
        const Real Tv = thetav_local(pack_idx)[vec_idx] * exner_function(P);
        const Real T = temperature_from_virtual_temperature(
            Tv, qv_local(pack_idx)[vec_idx]);

        const Real qvsat = qvsat_tetens(T, P);

        relh(pack_idx)[vec_idx] = qv_local(pack_idx)[vec_idx] *
                                  FloatingPoint<>::safe_denominator(qvsat);
        temp(pack_idx)[vec_idx] = T;
      });
  Kokkos::parallel_for(
      "HostDynamics::CreateAtmosphereInterfaces", nlev_ + 1,
      KOKKOS_LAMBDA(const int k) {
        const int pack_idx = PackInfo::pack_idx(k);
        const int vec_idx = PackInfo::vec_idx(k);

        z(pack_idx)[vec_idx] = phi_local(pack_idx)[vec_idx] / gravity;
      });
  Real pblh = 100.0;

  return Atmosphere(nlev_, temp, p, relh, z, hydrostatic_dp, pblh);
}

void HostDynamics::update_atmospheric_state(Atmosphere& atm) const {
  using namespace constants;
  const auto p_local = p;
  const auto phi_local = phi;
  const auto thetav_local = thetav;
  const auto qv_local = qv;

  auto temperature = atm.temperature;
  auto rel_humidity = atm.relative_humidity;
  auto level_heights = atm.height;
  Kokkos::parallel_for(
      "HostDynamics:UpdateAtmosphereLevels", nlev_, KOKKOS_LAMBDA(const int k) {
        const int pack_idx = PackInfo::pack_idx(k);
        const int vec_idx = PackInfo::vec_idx(k);

        const Real P = p_local(pack_idx)[vec_idx];
        const Real Tv = thetav_local(pack_idx)[vec_idx] * exner_function(P);
        const Real T = temperature_from_virtual_temperature(
            Tv, qv_local(pack_idx)[vec_idx]);
        const Real qvsat = qvsat_tetens(T, P);

        temperature(pack_idx)[vec_idx] = T;
        rel_humidity(pack_idx)[vec_idx] =
            qv_local(pack_idx)[vec_idx] *
            FloatingPoint<>::safe_denominator(qvsat);
      });

  Kokkos::parallel_for(
      "HostDynamics::UpdateAtmosphereInterfaces", nlev_ + 1,
      KOKKOS_LAMBDA(const int k) {
        const int pack_idx = PackInfo::pack_idx(k);
        const int vec_idx = PackInfo::vec_idx(k);

        level_heights(pack_idx)[vec_idx] =
            phi_local(pack_idx)[vec_idx] / gravity;
      });
}

}  // namespace driver
}  // namespace haero
