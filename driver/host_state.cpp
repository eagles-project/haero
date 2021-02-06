#include "host_state.hpp"
#include "ekat/ekat_assert.hpp"
#include "haero/utils.hpp"
#include "haero/floating_point.hpp"
#include "haero/physical_constants.hpp"
#include <cmath>
#include <algorithm>
#include <sstream>

namespace haero {
namespace driver {

void HostDynamics::init_from_interface_heights(std::vector<Real> z0, 
        const AtmosphericConditions& ac) {
  EKAT_REQUIRE_MSG(z0.size() == nlev_+1, "number of initial heights must match number of levels + 1");  
  EKAT_REQUIRE_MSG(vector_is_monotone(z0), "initial heights must be a monotone array.");
  
  const bool increasing = (z0[1] > z0[0]);
  if (increasing) std::reverse(z0.begin(), z0.end());
  
  auto hphi0 = Kokkos::create_mirror_view(phi0);
  auto hrho0 = Kokkos::create_mirror_view(rho0);
  auto hp = Kokkos::create_mirror_view(p);
  auto hthetav = Kokkos::create_mirror_view(thetav);
  auto hqv = Kokkos::create_mirror_view(qv);
  auto hw = Kokkos::create_mirror_view(w);
  
  /// set interface geopotential
  for (int k=0; k<nlev_+1; ++k) {
    // Taylor et al. 2020 fig. 1 interface idx = k+1/2
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);
    hphi0(pack_idx)[vec_idx] = gravity_m_per_s2 * z0[k];
    hw(pack_idx)[vec_idx] = 0;
  }

  /// set midpoint pressure, density, virtual potential temperature, water vapor mixing ratio
  for (int k=0; k<nlev_; ++k) {
    // Taylor et al. 2020 fig. 1 level idx = k+1
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);
    
    const int kphalf_idx = k+1; // array idx of interface k + 1/2
    const int kmhalf_idx = k; // array idx of interface k - 1/2
    const auto kphalf_pack_idx = PackInfo::pack_idx(kphalf_idx);
    const auto kphalf_vec_idx = PackInfo::vec_idx(kphalf_idx);
    const auto kmhalf_pack_idx = PackInfo::pack_idx(kmhalf_idx);
    const auto kmhalf_vec_idx = PackInfo::vec_idx(kmhalf_idx);
    
    const Real phimid = 0.5*(hphi0(kmhalf_pack_idx)[kmhalf_vec_idx] +
                                        hphi0(kphalf_pack_idx)[kphalf_vec_idx]);
    const Real zmid = phimid/gravity_m_per_s2;
    const Real pres = hydrostatic_pressure_at_height(zmid, ac);
    const Real Tv = ac.Tv0 - ac.Gammav * zmid;
    hp(pack_idx)[vec_idx] = pres;
    hrho0(pack_idx)[vec_idx] = pres / (r_gas_dry_air_joule_per_k_per_kg * Tv);
    hthetav(pack_idx)[vec_idx] = Tv / exner_function(pres);
    hqv(pack_idx)[vec_idx] = water_vapor_mixing_ratio(zmid, ac);
  }
  
  ps = hydrostatic_pressure_at_height(0, ac);
  EKAT_ASSERT_MSG(FloatingPoint<Real>::equiv(ps,AtmosphericConditions::pref),
    "surface pressure must equal the reference pressure, 1000 hPa.");
  
  Kokkos::deep_copy(w,hw);
  Kokkos::deep_copy(phi0,hphi0);
  Kokkos::deep_copy(phi,phi0);
  Kokkos::deep_copy(rho0,hrho0);
  Kokkos::deep_copy(rho,rho0);
  Kokkos::deep_copy(thetav,hthetav);
  Kokkos::deep_copy(qv,hqv);
  Kokkos::deep_copy(p,hp);
}

void HostDynamics::init_from_interface_pressures(std::vector<Real> p0, const AtmosphericConditions& ac) {
  EKAT_REQUIRE_MSG(p0.size() == nlev_ + 1, "number of initial pressures must match number of levels + 1");
  EKAT_REQUIRE_MSG(vector_is_monotone(p0), "initial pressures must be a monotone array.");
  
  const bool increasing = (p0[1]>p0[0]);
  if (!increasing) std::reverse(p0.begin(), p0.end());
  EKAT_REQUIRE_MSG(p0.back()== AtmosphericConditions::pref, "surface pressure must initialize to 1000 hPa.");
  
  auto hphi0 = Kokkos::create_mirror_view(phi0);
  auto hrho0 = Kokkos::create_mirror_view(rho0);
  auto hp = Kokkos::create_mirror_view(p);
  auto hthetav = Kokkos::create_mirror_view(thetav);
  auto hqv = Kokkos::create_mirror_view(qv);
  auto hw = Kokkos::create_mirror_view(w);
  
  /// set interface geopotential
  for (int k=0; k<nlev_+1; ++k) {
    // Taylor et al. 2020 fig. 1 interface idx = k+1/2
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);
  
    hw(pack_idx)[vec_idx] = 0;
    const Real z = height_at_pressure(p0[k],ac);
    hphi0(pack_idx)[vec_idx] = gravity_m_per_s2 * z;
  }
  
  /// set midpoint pressure, density, virtual potential temperature, water vapor mixing ratio
  for (int k=0; k<nlev_; ++k) {
    // Taylor et al. 2020 fig. 1 level idx = k+1
    const auto pack_idx = PackInfo::pack_idx(k);
    const auto vec_idx = PackInfo::vec_idx(k);
    
    const int kphalf_idx = k+1; // array idx of interface k + 1/2
    const int kmhalf_idx = k; // array idx of interface k - 1/2
    const auto kphalf_pack_idx = PackInfo::pack_idx(kphalf_idx);
    const auto kphalf_vec_idx = PackInfo::vec_idx(kphalf_idx);
    const auto kmhalf_pack_idx = PackInfo::pack_idx(kmhalf_idx);
    const auto kmhalf_vec_idx = PackInfo::vec_idx(kmhalf_idx);
    
    const Real phimid = 0.5*(hphi0(kmhalf_pack_idx)[kmhalf_vec_idx] +
                                        hphi0(kphalf_pack_idx)[kphalf_vec_idx]);
    const Real zmid = phimid/gravity_m_per_s2;
    const Real pres = hydrostatic_pressure_at_height(zmid, ac);
    const Real Tv = ac.Tv0 - ac.Gammav * zmid;
    hp(pack_idx)[vec_idx] = pres;
    hrho0(pack_idx)[vec_idx] = pres / (r_gas_dry_air_joule_per_k_per_kg * Tv);
    hthetav(pack_idx)[vec_idx] = Tv / exner_function(pres);
    hqv(pack_idx)[vec_idx] = water_vapor_mixing_ratio(zmid, ac);
  }
  
  ps = p0.back();
  
  Kokkos::deep_copy(w,hw);
  Kokkos::deep_copy(phi0,hphi0);
  Kokkos::deep_copy(phi,hphi0);
  Kokkos::deep_copy(rho0,hrho0);
  Kokkos::deep_copy(rho,hrho0);
  Kokkos::deep_copy(thetav,hthetav);
  Kokkos::deep_copy(qv,hqv);
  Kokkos::deep_copy(p,hp);
}

std::string HostDynamics::info_string(int tab_level) const {
  std::string tabstr = indent_string(tab_level);
  std::ostringstream ss;
  ss << tabstr << "HostDynamics info:\n";
  tabstr += "\t";
  ss << "nlev = " << nlev_ << "\n";
  ss << "ps = " << ps << "\n";
  ss << "t = " << t << "\n";
  return ss.str();
}

void HostDynamics::update(const Real t, const AtmosphericConditions& ac) {
  // interface update
  // set local variables for lambda
  auto phi_local = phi;
  auto phi0_local = phi0;
  auto w_local = w;
  Kokkos::parallel_for("HostDynamics::InterfaceUpdate", PackInfo::num_packs(nlev_+1), 
    KOKKOS_LAMBDA (const int pack_idx) {
    for (int vi=0; vi<PackInfo::vec_end(nlev_+1,pack_idx); ++vi) {
      const Real geop = geopotential(t, phi0_local(pack_idx)[vi], ac);
      phi_local(pack_idx)[vi] = geop;
      w_local(pack_idx)[vi] = velocity(t, geop, ac);
    }
  });
  
  // midpoint update
  //TODO: Phi0mid could be computed at init, then kept.
  auto rho_local = rho;
  auto rho0_local = rho0;
  auto thetav_local = thetav;
  auto p_local = p;
  Kokkos::parallel_for("HostDynamics::MidpointUpdate", PackInfo::num_packs(nlev_), 
    KOKKOS_LAMBDA (const int pack_idx) {
    for (int vi=0; vi<PackInfo::vec_end(nlev_, pack_idx); ++vi) {
      const int k = PackInfo::array_idx(pack_idx,vi);
      const int kphalf_pack = PackInfo::pack_idx(k+1);
      const int kphalf_vec = PackInfo::vec_idx(k+1);
      
      const Real phimid = 0.5*(phi_local(pack_idx)[vi] + phi_local(kphalf_pack)[kphalf_vec]);
      const Real phi0mid = 0.5*(phi0_local(pack_idx)[vi] + phi0_local(kphalf_pack)[kphalf_vec]);
      
      rho_local(pack_idx)[vi] = density(t, phimid, phi0mid, rho0_local(pack_idx)[vi], ac);
      p_local(pack_idx)[vi] = pressure(rho_local(pack_idx)[vi], thetav_local(pack_idx)[vi]);
    }
  });
}

} // namespace driver 
} // namespace haero

