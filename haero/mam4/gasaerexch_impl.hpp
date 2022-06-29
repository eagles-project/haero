#ifndef HAERO_MAM4_GASAEREXCH_IMPL_HPP
#define HAERO_MAM4_GASAEREXCH_IMPL_HPP

#include <haero/atmosphere.hpp>
#include <haero/mam4/aero_config.hpp>

namespace haero {
namespace mam4 {

/// @class GasAerExchImpl
/// This class implements MAM4's gas/aersol exchange  parameterization. Its structure is
/// defined by the usage of the impl_ member in the AeroProcess class in
/// ../aero_process.hpp.
class GasAerExchImpl {
  static const int num_mode   = 4;
  static const int num_gas    = 13;
  static const int num_aer    = 7;
  static const int nait       = static_cast<int>(ModeIndex::Aitken);
  static const int npca       = static_cast<int>(ModeIndex::PrimaryCarbon);
  static const int igas_h2so4 = static_cast<int>(GasId::H2SO4);
  static const int igas_nh3   = static_cast<int>(GasId::NH3);
  static const int iaer_h2so4 = static_cast<int>(GasId::H2SO4);
  static const int iaer_so4   = static_cast<int>(AeroId::SO4);
  static const int iaer_pom   = static_cast<int>(AeroId::POM);
  static const int nsoa       = static_cast<int>(AeroId::SOA);
  //static const int npoa     = static_cast<int>(AeroId::POA);

  static constexpr Real mw_h2so4 = Constants::molec_weight_h2so4;

  // ratio of gas uptake coeff w.r.t. that of h2so4
  static constexpr Real soag_h2so4_uptake_coeff_ratio = 0.81;  // for SOAG
  static constexpr Real nh3_h2so4_uptake_coeff_ratio  = 2.08;  // for NH3

  haero::kokkos_device_type::view_2d<bool> l_mode_can_contain_species;
  haero::kokkos_device_type::view_2d<bool> l_gas_condense_to_mode;
  haero::kokkos_device_type::view_1d<bool> l_mode_can_age;
  haero::kokkos_device_type::view_1d<int>  idx_gas_to_aer;
  haero::kokkos_device_type::view_1d<int>  eqn_and_numerics_category;
  haero::kokkos_device_type::view_1d<Real> uptk_rate_factor;
  haero::kokkos_device_type::view_1d<Real> qgas_cur;
  haero::kokkos_device_type::view_1d<Real> qgas_avg;
  haero::kokkos_device_type::view_1d<Real> qgas_netprod_otrproc;
  haero::kokkos_device_type::view_2d<PackType> uptkaer;
  haero::kokkos_device_type::view_1d<Real> alnsg_aer;
  haero::kokkos_device_type::view_1d<int>  mode_aging_optaa;

 public:

  // name -- unique name of the process implemented by this class
  const char* name() const { return "MAM4 gas/aersol exchange"; }

  // init -- initializes the implementation with MAM4's configuration
  void init(const AeroConfig& config) {
    Kokkos::resize(l_mode_can_contain_species, num_aer, num_mode);
    Kokkos::resize(l_gas_condense_to_mode, num_gas, num_mode);
    Kokkos::resize(l_mode_can_age, num_aer);
    Kokkos::resize(idx_gas_to_aer, num_gas);
    Kokkos::resize(eqn_and_numerics_category, num_gas);
    Kokkos::resize(uptk_rate_factor, num_gas);

    Kokkos::resize(qgas_cur, num_gas);
    Kokkos::resize(qgas_avg, num_gas);
    Kokkos::resize(qgas_netprod_otrproc, num_gas);

    Kokkos::resize(uptkaer, num_gas, num_mode); // gas to aerosol mass transfer rate (1/s)

    Kokkos::resize(alnsg_aer, num_mode);
    Kokkos::resize(mode_aging_optaa, num_mode);

    
    //------------------------------------------------------------------
    // MAM currently assumes that the uptake rate of other gases
    // are proportional to the uptake rate of sulfuric acid gas (H2SO4).
    // Here the array uptk_rate_factor contains the uptake rate ratio
    // w.r.t. that of H2SO4.
    //------------------------------------------------------------------
    // Set default to 0, meaning no uptake
    Kokkos::parallel_for(num_gas,
      KOKKOS_LAMBDA(int k) {
        uptk_rate_factor(k) = 0.0;
      });

    // H2SO4 is the ref species, so the ratio is 1
    uptk_rate_factor(igas_h2so4) = 1.0;
 
    // For NH3 
    uptk_rate_factor(igas_nh3) = nh3_h2so4_uptake_coeff_ratio;

    // For SOAG. (igas_soag_bgn and igas_soag_end are the start- and end-indices)
    // Remove use of igas_soag_bgn but keep comments in case need to be added back
    // uptk_rate_factor(igas_soag:igas_soagzz) = soag_h2so4_uptake_coeff_ratio
    //         igas_soag = ngas + 1
    //         igas_soagzz = ngas + nsoa
    //         uptk_rate_factor(igas_soag:igas_soagzz) = soag_h2so4_uptake_coeff_ratio
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
  // NOTE: that both diags and tends are const below--this means their views
  // NOTE: are fixed, but the data in those views is allowed to vary.
  KOKKOS_INLINE_FUNCTION
  void compute_tendencies(const AeroConfig& config,
                          const TeamType& team, Real t, Real dt,
                          const Atmosphere& atm,
                          const Prognostics& progs,
                          const Diagnostics& diags,
                          const Tendencies& tends) const {
    
    const int nghq = 2;  // set number of ghq points for direct ghq
    const bool l_calc_gas_uptake_coeff = config.calculate_gas_uptake_coefficient;
    //=========================================================================
    // Initialize the time-step mean gas concentration (explain why?)
    //=========================================================================
    Kokkos::parallel_for(num_gas,
      KOKKOS_LAMBDA(int k) {
        qgas_avg(k) = 0.0;
      });
    
    if (l_calc_gas_uptake_coeff) {
      const int nk = PackInfo::num_packs(atm.num_levels());
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nk),
      KOKKOS_LAMBDA(int k) {
        //=========================================================================
        // Calculate the reference uptake coefficient for all aerosol modes using
        // properties of the H2SO4 gas
        //=========================================================================
        PackType uptkaer_ref[num_mode]={0,0,0,0};  // initialize with zero (-> no uptake)
        bool l_condense_to_mode[num_mode]={true,true,true,true};   // do calcullation for ALL modes
        const int igas = igas_h2so4; // use properties of the H2SO4 gas
        const PackType &temp = atm.temperature(k);
        const PackType &pmid = atm.pressure(k);
        const PackType &dgn_awet = atm.wet_geometric_mean_diameter(k);
        const Real pstd = Constants::pressure_stp;
        const Real mw_h2so4 = Constants::molec_weight_h2so4;
        const Real mw_air = Constants::molec_weight_dry_air;
        const Real vol_molar_h2so4 = Constants::molec_weight_h2so4;
        const Real vol_molar_air = Constants::molec_weight_dry_air;
        const Real accom_coef_h2so4 = 0.65;     
        const Real r_universal = Constants::r_gas;
        const Real r_pi = Constants::pi;
        const Real beta_inp = 0; // quadrature parameter (--)
        const int ntot_amode = num_mode;
        const int nghq = config.number_gauss_points_for_integration;

        Real alnsg_aer[num_mode]={0,0,0,0};  // TODO: Figure out what this should be.
#if 0
        gas_aer_uptkrates_1box1gas(l_condense_to_mode,                           // in
                                temp, pmid, pstd, mw_h2so4, 1000*mw_air,    // in
                                vol_molar_h2so4, vol_molar_air,             // in
                                accom_coef_h2so4, r_universal, r_pi, beta_inp, // in
                                ntot_amode, nghq,                           // in
                                dgn_awet, alnsg_aer,                        // in
                                uptkaer_ref);                               // inout
#endif
        //========================================================================================
        // Assign uptake rate to each gas species and each mode using the ref. value uptkaer_ref
        // calculated above and the uptake rate factor specified as constants at the
        // beginning of the module
        //========================================================================================
        Kokkos::parallel_for(num_gas*num_mode,
        KOKKOS_LAMBDA(int k) {
          uptkaer(k/num_mode,k%num_mode) = 0.0; //default is no uptake
        });
        Kokkos::parallel_for(num_gas*num_mode,
        KOKKOS_LAMBDA(int k) {
          const int igas = k/num_mode;
          const int imode = k%num_mode;
          if (l_gas_condense_to_mode(igas,imode))
            uptkaer(igas,imode) = uptkaer_ref[imode] * uptk_rate_factor(igas);
        });

        // total uptake rate (sum of all aerosol modes) for h2so4.
        // Diagnosed for calling routine. Not used in this subroutne.
        tends.uptkrate_h2so4(k) = 0;
        for (int n=0; n<num_mode; ++n)
          tends.uptkrate_h2so4(k) += uptkaer(igas_h2so4,n);
      });
    }


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
