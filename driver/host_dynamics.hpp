#ifndef HAERO_DRIVER_HOST_STATE_HPP
#define HAERO_DRIVER_HOST_STATE_HPP

#include "haero/haero.hpp"
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
    /// layer thickness in height (midpoint variable)
    ColumnView dz;
    /** approximate layer thickness in pressure (midpoint variable)

      "approximate" because it's calculated based on the hydrostatic assumption;
      in a non-hydrostatic atmosphere, it's an approximation.

      This is analogous to "pseudo-density" in HOMME-NH.
    */
    ColumnView hydrostatic_dp;

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
      p("plev",PackInfo::num_packs(nl)),
      dz("dzlev", PackInfo::num_packs(nl)),
      hydrostatic_dp("pseudo_density", PackInfo::num_packs(nl)),
      t(0), ps(0), nlev_(nl),
      phi0("phi0",PackInfo::num_packs(nl+1)),
      rho0("rho0", PackInfo::num_packs(nl)),
      phydro_int("hydrostatic_pressure_interface", PackInfo::num_packs(nl+1))
    {}

    HostDynamics() = delete;

    /// Update dynamics to new time
    void update(const Real newt, const AtmosphericConditions& ac);

    /** Initialize column data at t=0 to stationary, hydrostatic balance

      Calling this method will change the ztop/ptop values of AtmosphericConditions to match the
      first height value.

      @param z0 initial heights of interfaces
      @param ac atmospheric conditions to define parameters
    */
    void init_from_interface_heights(std::vector<Real> z0,
        AtmosphericConditions& ac);

    /** Initialize column data at t=0 to stationary, hydrostatic balance

      Calling this method will change the ztop/ptop values of AtmosphericConditions to match the
      first pressure value.

      @param p0 initial pressure at interfaces
      @param ac atmospheric conditions to define parameters
    */
    void init_from_interface_pressures(std::vector<Real> p0,
         AtmosphericConditions& ac);

    /** @brief initialize column data using equally-spaced levels in height coordinates.

      layer thickness dz = ac.ztop/nlev

    */
    void init_from_uniform_heights(const AtmosphericConditions& ac);


    /** @brief initialize column data using equally-spaced levels in pressure coordinates.

      layer thickness dp = ac.ptop/nlev

      Note: this layer thickness does not equate to the hydrostatic_dp layer thickness, because
      the column is only in hydrostatic balance at t = integer multiples of ac.tperiod.

    */
    void init_from_uniform_pressures(const AtmosphericConditions& ac);

    /** Write basic information about *this to a string.
    */
    std::string info_string(int tab_level=0) const;

    /** @brief Initializes dynamics variables in a netCDF file.

      @param [in/out] writer
      @param [in] conds
    */
    void nc_init_dynamics_variables(NcWriter& writer, const AtmosphericConditions& conds) const;

    /** @brief Writes dynamics variables' data to a netCDF file.

      @param [in/out] writer
      @param [in] time_idx
    */
    void nc_write_data(NcWriter& writer, const size_t time_idx) const;

    /** @brief Creates a haero::Atmosphere instance to provide dynamics input data to parameterizations.

      Note: Atmosphere member variables that are identical to existing HostDynamics ColumnView members
        are simple view copies (e.g., hydrostatic_dp and pressure). Atmosphere member variables
        that are not already ColumnViews in HostDynamics must be passed as ColumnViews to
        this function (e.g., temperature, height, and relative humidity).

      @param [in/out] temp view to store temperature (rank 1, size = nlev)
      @param [in/out] relh view to store relative humidity (rank 1, size = nlev)
      @param [in/out] z view to store level interface heights (rank 1, size = nlev + 1)
    */
    Atmosphere create_atmospheric_state(ColumnView temp,
      ColumnView relh, ColumnView z) const;

    /** @brief Updates a haero::Atmosphere instance's data (e.g., after a change in the time variable)

      @param atm
    */
    void update_atmospheric_state(Atmosphere& atm) const;

    inline int nlev() const {return nlev_;}


#ifndef HAERO_USE_CUDA // variables below are meant to be private, but must be public for gpu builds
  private:
#endif
    /// number of levels in column
    int nlev_;
    /// initial density at the surface
    Real rho0surf;
    /// intial geopotential values
    ColumnView phi0;
    /// initial density values
    ColumnView rho0;
    /// hydrostatic pressure (interface variable)
    ColumnView phydro_int;

    /** @brief compute discrete approximations of vertical derivatives using centered finite differences
      as described by Taylor et al. 2020.
    */
    void update_thickness(const AtmosphericConditions& conds);
};

/** @brief Defines the Lagrangian geopotential for HostDynamics' 1d toy model.
  This function must be called before the other dynamics functions, because they require
  the current (time t) geopotential.

  @param [in] t time [s]
  @param [in] phi0 or @f$ \phi_0@f$ initial geopotential height [m<sup>2</sup>s<sup>-2</sup]
  @param [in] ac toy model parameters
  @return geopotential at time t, @f$ \phi(t) @f$ such that $\phi(0) = \phi_0 @f$ [m<sup>2</sup>s<sup>-2</sup]
*/
KOKKOS_INLINE_FUNCTION
Real geopotential(const Real t, const Real phi0, const AtmosphericConditions& ac) {
  using namespace constants;
  const Real tanarg = pi * phi0 / (2*gravity*ac.ztop);
  const Real exparg = ac.w0 * ac.tperiod * square(std::sin(pi*t/ac.tperiod))/(ac.ztop);
  return 2*gravity*ac.ztop*std::atan(std::tan(tanarg)*std::exp(exparg)) / pi;
}

/** @brief Defines velocity for a Lagrangian geopotential surface.

  @param [in] t time [s]
  @param [in] phi geopotential at time t, @f$\phi(t)@f$ [m<sup>2</sup>s<sup>-2</sup]
  @param [in] ac toy model parameters
  @return vertical velocity w [m/s]
*/
KOKKOS_INLINE_FUNCTION
Real velocity(const Real t, const Real phi, const AtmosphericConditions& ac) {
  using namespace constants;
  return ac.w0 * std::sin(phi/(gravity*ac.ztop))*std::sin(2*pi *t / ac.tperiod);
}

/** @brief Defines the density of a Lagrangian parcel for the 1d toy model.


  @param [in] t time [s]
  @param [in] phi geopotential at time t, @f$\phi(t)@f$ [m<sup>2</sup>s<sup>-2</sup]
  @param [in] phi0 or @f$ \phi_0@f$ initial geopotential height [m<sup>2</sup>s<sup>-2</sup]
  @param [in] rho0 or @f$ \rho_0@f$ initial density [kg m<sup>-3</sup]
  @param [in] ac toy model parameters
  @return density @f$\rho(t)@f$ such that @f$\rho(0) = \rho_0@f$
*/
KOKKOS_INLINE_FUNCTION
Real density(const Real t, const Real phi, const Real phi0, const Real rho0, const AtmosphericConditions& ac) {
  using namespace constants;
  const Real cosarg1 = pi*phi/(gravity*ac.ztop);
  const Real cosarg2 = 2*pi*t/ac.tperiod;
  const Real cosarg3 = pi*phi0/(gravity*ac.ztop);
  const Real exparg = std::cos(cosarg1)*std::cos(cosarg2)-std::cos(cosarg3);
  return rho0*std::exp(ac.w0*ac.tperiod*exparg/(2*ac.ztop));
}

/** @brief Computes pressure, given density and virtual potential temperature.

  @param [in] rho density [kg m<sup>-3</sup>]
  @param [in] thetav virtual potential temperatore [K]
  @return pressure [Pa]
*/
KOKKOS_INLINE_FUNCTION
Real pressure(const Real rho, const Real thetav) {
  using namespace constants;
  const Real coeff = std::pow(AtmosphericConditions::pref, -AtmosphericConditions::kappa)*
    r_gas_dry_air;
  return std::pow(coeff*rho*thetav, 1/(1-AtmosphericConditions::kappa));
}

/** @brief Computes the water vapor saturation mixing ratio.

  Uses the Tetens equation from Soong-Ogura 1973 equation (A1) or Klemp-Wilhelmson 1978 eqn. (2.11).
  Used to define relative humidity as rh = qv / qvsat.

  @param [in] T temperature (dry air, not virtual temperature) [K]
  @param [in] Pressure [Pa]
  @return qvsat the value of the water vapor mixing ratio of a saturated parcel at the same T and P.
*/
KOKKOS_INLINE_FUNCTION
Real qvsat_tetens(const Real T, const Real p) {
  static constexpr Real half15ln10 = 17.269388197455342630;
  static constexpr Real tetens_coeff = 380.042;
  return tetens_coeff * std::exp(half15ln10*(T - 273)/(T-36)) / p;
}

struct DynamicsInterfaceUpdate {
  Real t;
  ColumnView phi;
  ConstColumnView phi0;
  ColumnView w;
  AtmosphericConditions conds;
  int nlev;

  DynamicsInterfaceUpdate(const Real newt, ColumnView phi_, ColumnView phi0_, ColumnView w_,
    const AtmosphericConditions ac, const int nl) :
    t(newt),
    phi(phi_),
    phi0(phi0_),
    w(w_),
    conds(ac),
    nlev(nl) {}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int k) const {
    using namespace constants;
    const int pack_idx = PackInfo::pack_idx(k);
    const int vec_idx = PackInfo::vec_idx(k);
    Real geop;
    if (k>0 && k<nlev) {
      geop = geopotential(t, phi0(pack_idx)[vec_idx], conds);
    }
    else {
      geop = (k==0 ? gravity*conds.ztop : 0);
    }
    phi(pack_idx)[vec_idx] = geop;
    w(pack_idx)[vec_idx] = velocity(t, geop, conds);
  }
};

struct HydrostaticPressureUpdate {
  ConstColumnView phi;
  ColumnView phydro_int;
  AtmosphericConditions conds;

  HydrostaticPressureUpdate(ColumnView phi_, ColumnView ph_, const AtmosphericConditions ac) :
    phi(phi_),
    phydro_int(ph_),
    conds(ac)
  {}

  KOKKOS_INLINE_FUNCTION
  void operator () (const int k) const {
    using namespace constants;
    const int pack_idx = PackInfo::pack_idx(k);
    const int vec_idx = PackInfo::vec_idx(k);
    phydro_int(pack_idx)[vec_idx] =
      hydrostatic_pressure_at_height(phi(pack_idx)[vec_idx]/gravity, conds);
  }

};

} // namespace driver
} // namespace haero
#endif
