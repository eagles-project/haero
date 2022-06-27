#ifndef HAERO_MAM4_NUCLEATION_IMPL_HPP
#define HAERO_MAM4_NUCLEATION_IMPL_HPP

#include <haero/atmosphere.hpp>
#include <haero/mam4/aero_config.hpp>

namespace haero {
namespace mam4 {

/// @class NucleationImpl
/// This class implements MAM4's nucleation parameterization. Its structure is
/// defined by the usage of the impl_ member in the AeroProcess class in
/// ../aero_process.hpp.
class NucleationImpl {
 public:

  // name -- unique name of the process implemented by this class
  const char* name() const { return "MAM4 nucleation"; }

  // init -- initializes the implementation with MAM4's configuration
  void init(const AeroConfig& config) {
    // Nothing to do here--all config information is fixed.
  }

  // validate -- validates the given atmospheric state and prognostics against
  // assumptions made by this implementation, returning true if the states are
  // valid, false if not
  KOKKOS_INLINE_FUNCTION
  bool validate(const AeroConfig& config,
                const TeamType& team,
                const Atmosphere& atm,
                const Prognostics& progs) const {

    // Make sure relevant atmospheric quantities are physical.
    const int nk = PackInfo::num_packs(atm.num_levels());
    int violations = 0;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nk),
      KOKKOS_LAMBDA(int k, int& violation) {
        if ((atm.temperature(k) < 0).any() ||
            (atm.pressure(k) < 0).any() ||
            (atm.vapor_mixing_ratio(k) < 0).any()) {
          violation = 1;
        }
      }, violations);

    if (violations == 0) { // all clear so far
      // Check for negative mixing ratios.
      Kokkos::parallel_reduce(Kokkos::TeamThreadRange(team, nk),
        KOKKOS_LAMBDA(int k, int& violation) {
          for (int mode = 0; mode < 4; ++mode) { // check mode mmrs
            if ((progs.n_mode[mode](k) < 0).any()) {
              ++violation;
            } else {
              for (int spec = 0; spec < 7; ++spec) { // check aerosol mmrs
                if ((progs.q_aero[mode][spec](k) < 0).any()) {
                  ++violation;
                  break;
                }
              }
            }
            if (violation > 0) break;
          }
          if (violation == 0) {
            for (int gas = 0; gas < 13; ++gas) { // check gas mmrs
              if ((progs.q_gas[gas](k) < 0).any()) ++violation;
            }
          }
        }, violations);
    }
    return (violations > 0);
  }

  // compute_tendencies -- computes tendencies and updates diagnostics
  KOKKOS_INLINE_FUNCTION
  void compute_tendencies(const AeroConfig& config,
                          const TeamType& team, Real t, Real dt,
                          const Atmosphere& atm,
                          const Prognostics& progs,
                          Diagnostics& diags,
                          Tendencies& tends) const {
    // For now, just zero all the tendencies.
    const int nk = PackInfo::num_packs(atm.num_levels());
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nk),
      KOKKOS_LAMBDA(int k) {
        for (int mode = 0; mode < 4; ++mode) {
          tends.n_mode[mode](k) = 0;
          for (int spec = 0; spec < 7; ++spec) {
            tends.q_aero[mode][spec](k) = 0;
          }
        }
        for (int gas = 0; gas < 13; ++gas) {
          tends.q_gas[gas](k) = 0;
        }
      });
  }
};

} // namespace mam4
} // namespace haero

#endif
