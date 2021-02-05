#ifndef HAERO_DRIVER_HOST_STATE_HPP
#define HAERO_DRIVER_HOST_STATE_HPP

#include "haero/haero_config.hpp"
#include "haero/atmosphere.hpp"
#include "host_params.hpp"
#include "ncwriter.hpp"
#include "Kokkos_Core.hpp"
#include <string>
#include <vector>

namespace haero {
namespace driver {

class HostDynamics {
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
    void update(const Real newt);
    
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
        
    /** Write column data to a new time index in a netcdf file.
    
      @param writer netcdf writer
      @param time_idx index of time to add
    */
    void update_ncdata(NcWriter& writer, const size_t time_idx) const;
    
    /** Write basic information about *this to a string.
    */
    std::string info_string() const;
    
  protected:
    /// number of levels in column
    int nlev_;
    /// intial geopotential values
    ColumnView phi0;
    /// initial density values
    ColumnView rho0;
        
};

struct DynamicsUpdate {
};

}
}
#endif