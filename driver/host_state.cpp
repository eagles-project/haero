#include "host_state.hpp"
#include "ekat/ekat_assert.hpp"
#include "haero/utils.hpp"
#include "haero/physical_constants.hpp"
#include <cmath>
#include <algorithm>

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
  
  Kokkos::deep_copy(w,hw);
  Kokkos::deep_copy(phi0,hphi0);
  Kokkos::deep_copy(phi,hphi0);
  Kokkos::deep_copy(rho0,hrho0);
  Kokkos::deep_copy(rho,hrho0);
  Kokkos::deep_copy(thetav,hthetav);
  Kokkos::deep_copy(qv,hqv);
  Kokkos::deep_copy(p,hp);
}

} // namespace driver 
} // namespace haero

