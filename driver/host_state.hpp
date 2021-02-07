#ifndef HAERO_DRIVER_HOST_STATE_HPP
#define HAERO_DRIVER_HOST_STATE_HPP

#include "haero/haero_config.hpp"
#include "haero/atmosphere.hpp"
#include "haero/math_helpers.hpp"
#include "host_params.hpp"
#include "ncwriter.hpp"
#include "Kokkos_Core.hpp"
#include <string>
#include <vector>

namespace haero {
namespace driver {

class HostDynamics final {
  public: 
    using ColumnView = Kokkos::View<PackType*>;
    
    /// vertical velocity (interface variable)
    ColumnView w;
    /// geopotential (interface variable)
    ColumnView phi;
    
    /// density (midpoint variable)
    ColumnView rho;
    /// virtual potential temperature (midpoint variable)
    ColumnView thetav;
    /// water vapor mass mixing ratio (midpoint variable)
    ColumnView qv;
    /// pressure (midpoint variable)
    ColumnView p;
    
    /// elapsed time
    Real t;
    /// surface pressure
    Real ps;
    
    /** Constructor. Allocates memory, but does not initialize column data
     
      @param nl number of levels in column
    */
    HostDynamics(const int nl) :
      w("w",PackInfo::num_packs(nl+1)),
      phi("phi",PackInfo::num_packs(nl+1)),
      rho("rho",PackInfo::num_packs(nl)),
      thetav("thetav",PackInfo::num_packs(nl)),
      qv("qv",PackInfo::num_packs(nl)),
      p("p",PackInfo::num_packs(nl)),
      t(0), ps(0), nlev_(nl), 
      phi0("phi0",PackInfo::num_packs(nl+1)),
      rho0("rho0", PackInfo::num_packs(nl))
    {}
       
    HostDynamics() = delete;
    
    /// Update dynamics to new time   
    void update(const Real newt, const AtmosphericConditions& ac);
    
    /** Initialize column data at t=0 to stationary, hydrostatic balance
    
      @param z0 initial heights of interfaces
      @param ac atmospheric conditions to define parameters
    */
    void init_from_interface_heights(std::vector<Real> z0, 
        const AtmosphericConditions& ac);
    
    /** Initialize column data at t=0 to stationary, hydrostatic balance
    
      @param p0 initial pressure at interfaces
      @param ac atmospheric conditions to define parameters
    */
    void init_from_interface_pressures(std::vector<Real> p0, 
        const AtmosphericConditions& ac);
            
    /** Write basic information about *this to a string.
    */
    std::string info_string(int tab_level=0) const;
    
    void nc_init_dynamics_variables(NcWriter& writer, const AtmosphericConditions& conds) const;
    
    void nc_write_data(NcWriter& writer, const size_t time_idx) const;
    
    Atmosphere create_atmospheric_state(Kokkos::View<PackType*> temp, 
      Kokkos::View<PackType*> relh, Kokkos::View<PackType*> z) const;
    
    void update_atmospheric_state(Atmosphere& atm) const;
    
  protected:
    /// number of levels in column
    int nlev_;
    /// intial geopotential values
    ColumnView phi0;
    /// initial density values
    ColumnView rho0;
    
    Real rho0surf;
};

KOKKOS_INLINE_FUNCTION
Real geopotential(const Real t, const Real phi0, const AtmosphericConditions& ac) {
  const Real tanarg = pi * phi0 / (2*gravity_m_per_s2*ac.ztop);
  const Real exparg = ac.w0 * ac.tperiod * square(std::sin(pi*t/ac.tperiod))/(ac.ztop);
  return 2*gravity_m_per_s2*ac.ztop*std::atan(std::tan(tanarg)*std::exp(exparg)) / pi;
}

KOKKOS_INLINE_FUNCTION
Real velocity(const Real t, const Real phi, const AtmosphericConditions& ac) {
  return ac.w0 * std::sin(phi/(gravity_m_per_s2*ac.ztop))*std::sin(2*pi *t / ac.tperiod);
}

KOKKOS_INLINE_FUNCTION
Real density(const Real t, const Real phi, const Real phi0, const Real rho0, const AtmosphericConditions& ac) {
  const Real cosarg1 = pi*phi/(gravity_m_per_s2*ac.ztop);
  const Real cosarg2 = 2*pi*t/ac.tperiod;
  const Real cosarg3 = pi*phi0/(gravity_m_per_s2*ac.ztop);
  const Real exparg = std::cos(cosarg1)*std::cos(cosarg2)-std::cos(cosarg3);
  return rho0*std::exp(ac.w0*ac.tperiod*exparg/(2*ac.ztop));
}

KOKKOS_INLINE_FUNCTION
Real pressure(const Real rho, const Real thetav) {
  const Real coeff = std::pow(AtmosphericConditions::pref, -AtmosphericConditions::kappa)*
    r_gas_dry_air_joule_per_k_per_kg;
  return std::pow(coeff*rho*thetav, 1/(1-AtmosphericConditions::kappa));
}

KOKKOS_INLINE_FUNCTION
Real qvsat_tetens(const Real T, const Real p) {
  static constexpr Real half15ln10 = 17.269388197455342630;
  static constexpr Real tetens_coeff = 380.042;
  return tetens_coeff * std::exp(half15ln10*(T - 273)/(T-36)) / p;
}

} // namespace driver
} // namespace haero
#endif