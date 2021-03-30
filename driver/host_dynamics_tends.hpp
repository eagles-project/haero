#ifndef HAERO_DYNAMICS_TENDENCIES_HPP
#define HAERO_DYNAMICS_TENDENCIES_HPP

#include "haero/haero.hpp"
#include "host_dynamics.hpp"
#include "host_params.hpp"

namespace haero {
namespace driver {

/** @class DynamicsTendencies
  This type stores tendencies for the driver's 1d dynamics model, for use with time
  stepping methods.
*/
class DynamicsTendencies final {
  public:
    /// velocity tendency
    ColumnView w_tend;
    /// geopotential tendency
    ColumnView phi_tend;
    /// density tendency
    ColumnView rho_tend;
    /// virtual potential temperature tendency
    ColumnView thetav_tend;
    /// water vapor mixing ratio tendency
    ColumnView qv_tend;

    /** Constructor. Allocates memory, initializes to zero.

      @param [in] nl number of vertical levels
    */
    DynamicsTendencies(const int nl) :
      w_tend("w_tend", PackInfo::num_packs(nl+1)),
      phi_tend("phi_tend", PackInfo::num_packs(nl+1)),
      rho_tend("rho_tend", PackInfo::num_packs(nl)),
      thetav_tend("thetav_tend", PackInfo::num_packs(nl)),
      qv_tend("qv_tend", PackInfo::num_packs(nl)) {}

    DynamicsTendencies() = delete;

    /// compute all tendencies (dynamics only)
    void compute(const Real t, ConstColumnView phi, ConstColumnView rho, const AtmosphericConditions& conds);

  private:
};

/** @brief Defines the velocity tendency of the Lagrangian geopotential surface phi at time t.

  @param [in] t time
  @param [in] phi geopotential
  @param [in] conds 1d model parameters
  @return @f$dw/dt@f$
*/
KOKKOS_INLINE_FUNCTION
Real wtend(const Real t, const Real phi, const AtmosphericConditions& conds) {
  using namespace constants;
  const Real sinarg = pi * phi / (gravity * conds.ztop);
  const Real targ = 2 * pi * t / conds.tperiod;
  const Real maxamp1 = 2 * pi * conds.w0 / conds.tperiod;
  const Real maxamp2 = pi * square(conds.w0) / (2 * conds.ztop);
  return maxamp1 * std::sin(sinarg) * std::cos(targ) +
    maxamp2 * std::sin(2*sinarg) * square(std::sin(targ));
}


/** @brief Defines the geopotential tendency of the Lagrangian geopotential surface phi at time t.

  @param [in] t time
  @param [in] phi geopotential
  @param [in] conds 1d model parameters
  @return @f$d\phi/dt@f$
*/
KOKKOS_INLINE_FUNCTION
Real phitend(const Real t, const Real phi, const AtmosphericConditions& conds) {
  using namespace constants;
  const Real sinarg = pi * phi / (gravity * conds.ztop);
  const Real targ = 2 * pi * t / conds.tperiod;
  return gravity * conds.w0 * std::sin(sinarg) * std::sin(targ);
}


/** @brief Defines the density tendency of the Lagrangian geopotential surface phi at time t.

  @param [in] t time
  @param [in] phi geopotential
  @param [in] rho density
  @param [in] conds 1d model parameters
  @return @f$d\rho/dt@f$
*/
KOKKOS_INLINE_FUNCTION
Real rhotend(const Real t, const Real phi, const Real rho, const AtmosphericConditions& conds) {
  using namespace constants;
  const Real cosarg = pi * phi / (gravity * conds.ztop);
  const Real targ = 2*pi*t/conds.tperiod;
  return rho * pi * conds.w0 / conds.ztop * std::cos(cosarg) * std::sin(targ);
}

/** @brief Defines the virtual potential temperature tendency of the Lagrangian geopotential surface phi at time t.

  @return @f$d\theta_v/dt@f$
*/
KOKKOS_INLINE_FUNCTION
Real thetavtend() {return 0;}

/** @brief Defines the water vapor mixing ratio tendency of the Lagrangian geopotential surface phi at time t.

  @return @f$dq_v/dt@f$
*/
KOKKOS_INLINE_FUNCTION
Real qvtend() {return 0;}

} // namespace driver
} // namespace haero
#endif
