#ifndef HAERO_SIMPLE_NUCLEATION_PROCESS_HPP
#define HAERO_SIMPLE_NUCLEATION_PROCESS_HPP

#include <iomanip>

#include "haero/aerosol_process.hpp"
#include "haero/constants.hpp"
#include "haero/conversions.hpp"
#include "haero/processes/kerminen2002.hpp"
#include "haero/processes/merikanto2007.hpp"
#include "haero/processes/vehkamaki2002.hpp"
#include "haero/processes/wang2008.hpp"

namespace haero {

/// @class SimpleNucleationProcess
/// This aerosol process implements homogeneous binary/ternary nucleation
/// involving sulfuric acid and methane gases. It is based on classical
/// nucleation theory as parameterized by Vehkamaki et al (2002) and
/// Merikanto et al (2007). Nucleated particles are placed into the appropriate
/// mode(s) based on their computed size, or they are grown to fit the mode
/// with the minimum size.
class SimpleNucleationProcess final
    : public DeviceAerosolProcess<SimpleNucleationProcess> {
  using RealVector = kokkos_device_type::view_1d<Real>;
  using IntVector = kokkos_device_type::view_1d<int>;

  // This struct contains cross validation parameters for use with skywalker.
  struct SkywalkerParams {
    // Specified number concentration of H2SO4.
    Real c_h2so4;

    // Constructor
    SkywalkerParams():
      c_h2so4(-1.0) {}
  };

  //----------------------------------------------------------
  //                  Adjustable parameters
  //----------------------------------------------------------
  // All of the named parameters below can be set using the
  // set_param() methods.
  //----------------------------------------------------------

  /// Adjustment factor applied to nucleation rate (default: 1)
  Real nucleation_rate_factor_;

  /// Adjustment factor applied for planetary boundary layer (default: 1)
  Real pbl_factor_;

  /// Adjustment factor applied to tendency for nucleated species (default: 1)
  Real tendency_factor_;

  /// Nucleation method selection (default: 2)
  /// 2 = Vehkamaki et al (2002) binary nucleation
  /// 3 = Merikanto el al (2007) ternary nucleation
  int nucleation_method_;

  /// Planetary boundary layer (PBL) method selection (default: 0)
  /// 0 = no PBL adjustment
  /// 1 = first-order
  /// 2 = second-order
  int pbl_method_;

  //----------------------------------------------------------
  //                       Bookkeeping
  //----------------------------------------------------------
  // These metadata allow the process to work with the given
  // modal aerosol configuration without having to look stuff
  // up all the time.

  // Index of H2SO4 gas
  int igas_h2so4_;

  // Index of NH3 gas
  int igas_nh3_;

  /// Number of aerosol modes.
  int num_modes_;

  // Species and population indices of SO4 aerosol within aerosol modes
  IntVector iaer_so4_, ipop_so4_;

  // Species and population indices of NH4 aerosol within aerosol modes
  IntVector iaer_nh4_, ipop_nh4_;

  // The geometric mean particle diameters for all aerosol modes
  RealVector d_mean_aer_;

  /// The minimum particle diameters for all aerosol modes
  RealVector d_min_aer_;

  /// The maximum particle diameters for all aerosol modes
  RealVector d_max_aer_;

  // Molecular weights of H2SO4 and NH3 gases.
  Real mu_h2so4_, mu_nh3_;

  // Molecular weights of SO4 and NH4 aerosols.
  Real mu_so4_, mu_nh4_;

  // Skywalker cross-validation parameters.
  int skywalker_mode_;
  SkywalkerParams skywalker_;

 public:
  /// Constructor
  SimpleNucleationProcess();

  /// Destructor
  KOKKOS_INLINE_FUNCTION
  ~SimpleNucleationProcess() {}

  /// Copy constructor (for transferring between host and device)
  KOKKOS_INLINE_FUNCTION
  SimpleNucleationProcess(const SimpleNucleationProcess &rhs)
      : DeviceAerosolProcess<SimpleNucleationProcess>(rhs),
        nucleation_rate_factor_(rhs.nucleation_rate_factor_),
        pbl_factor_(rhs.pbl_factor_),
        tendency_factor_(rhs.tendency_factor_),
        nucleation_method_(rhs.nucleation_method_),
        pbl_method_(rhs.pbl_method_),
        igas_h2so4_(rhs.igas_h2so4_),
        igas_nh3_(rhs.igas_nh3_),
        num_modes_(rhs.num_modes_),
        iaer_so4_(rhs.iaer_so4_),
        iaer_nh4_(rhs.iaer_nh4_),
        d_mean_aer_(rhs.d_mean_aer_),
        d_min_aer_(rhs.d_min_aer_),
        d_max_aer_(rhs.d_max_aer_),
        skywalker_mode_(rhs.skywalker_mode_),
        skywalker_(rhs.skywalker_) {}

 protected:
  void init_(const ModalAerosolConfig &config) override;

  KOKKOS_INLINE_FUNCTION
  void run_(Real t, Real dt, const Prognostics &prognostics,
            const Atmosphere &atmosphere, const Diagnostics &diagnostics,
            const Tendencies &tendencies) const override {
    // Do we have any gas from which to nucleate new aerosol particles?
    if ((igas_h2so4_ == -1) or ((nucleation_method_ == 3) and (igas_nh3_ == -1))) {
      return;
    }

    // Do we track the aerosols nucleated from our gases?
    bool have_nucleated_aerosols = false;
    for (int m = 0; m < num_modes_; ++m) {
      if ((iaer_so4_(m) != -1) or (iaer_nh4_(m) != -1)) {
        have_nucleated_aerosols = true;
        break;
      }
    }
    if (not have_nucleated_aerosols) return;

    // If we're in skywalker mode, go do that thing.
    if (skywalker_mode_ == 1) {
      run_skywalker_mode_(t, dt, prognostics, atmosphere, diagnostics, tendencies);
      return;
    }

    // Calculate tendencies for aerosols/gases at each vertical level k.
    const int nk = atmosphere.temperature.extent(0);
    Kokkos::parallel_for(
        "nucleation_rate", nk, KOKKOS_LAMBDA(const int k) {
          const auto temp = atmosphere.temperature(k);
          const auto press = atmosphere.pressure(k);
          const auto qv = atmosphere.vapor_mixing_ratio(k);
          const auto h = atmosphere.height(k);
          const auto rho_d = gas_kinetics::air_mass_density(press, temp, qv);
          auto rel_hum = conversions::relative_humidity_from_vapor_mixing_ratio(
              qv, press, temp);
          const auto q_h2so4 = prognostics.gases(igas_h2so4_, k);  // mmr
          auto c_h2so4 =  // number concentration of H2SO4 [#/cc]
              1e6 * conversions::number_conc_from_mmr(q_h2so4, mu_h2so4_, rho_d);

          // Compute the base rate of nucleation using our selected method.
          PackType J;       // nucleation rate [#/cc]
          PackType r_crit;  // radius of critical cluster [nm]
          PackType n_crit;  // total # of molecules in a critical cluster [#]
          PackType n_crit_h2so4, n_crit_nh3;  // numbers of gas molecules in
                                              // the critical cluser [#]
          if (nucleation_method_ == 2) {      // binary nucleation
            auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(
                c_h2so4, temp, rel_hum);
            J = vehkamaki2002::nucleation_rate(c_h2so4, temp, rel_hum, x_crit);
            n_crit = vehkamaki2002::num_critical_molecules(c_h2so4, temp,
                                                           rel_hum, x_crit);
            n_crit_h2so4 = n_crit;
            n_crit_nh3 = 0;
            r_crit = vehkamaki2002::critical_radius(x_crit, n_crit);
          } else {
            EKAT_KERNEL_ASSERT(nucleation_method_ == 3);
            // Compute the molar/volume mixing ratio of NH3 gas [ppt].
            const auto q_nh3 = prognostics.gases(igas_nh3_, k);  // mmr
            auto xi_nh3 = 1e12 * conversions::vmr_from_mmr(q_nh3, mu_nh3_);

            auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                            c_h2so4, xi_nh3);
            J = exp(log_J);
            n_crit_h2so4 = merikanto2007::num_h2so4_molecules(log_J, temp,
                                                              c_h2so4, xi_nh3);
            n_crit_nh3 =
                merikanto2007::num_nh3_molecules(log_J, temp, c_h2so4, xi_nh3);
            n_crit = n_crit_h2so4 + n_crit_nh3;
            r_crit =
                merikanto2007::critical_radius(log_J, temp, c_h2so4, xi_nh3);
          }

          // Apply a correction for the planetary boundary layer if requested.
          if (pbl_method_ > 0) {
            PackType J_pbl(0);
            const auto within_pbl = (h <= atmosphere.planetary_boundary_height);
            if (pbl_method_ == 1) {
              J_pbl.set(within_pbl,
                        wang2008::first_order_pbl_nucleation_rate(c_h2so4));
            } else if (pbl_method_ == 2) {
              J_pbl.set(within_pbl,
                        wang2008::second_order_pbl_nucleation_rate(c_h2so4));
            }
            // If the modified nucleation rate is less than the original, J = 0.
            const auto J_pbl_lt_J = (J_pbl <= J);
            J.set(within_pbl and J_pbl_lt_J, 0);
            J.set(within_pbl and not J_pbl_lt_J, J_pbl);
            // All fresh nuclei within the planetary boundary layer are 1 nm in
            // diameter.
            r_crit.set(within_pbl and not J_pbl_lt_J, 0.5);
            // Estimate the number of particles in the critical cluѕter from its
            // mass and radius, assuming H2SO4 only.
            const auto pi_sixth = Constants::pi_sixth;
            const auto Na = Constants::avogadro;
            const auto mw_h2so4 = Constants::molec_weight_h2so4;
            static const Real rho_h2so4 = 1.8;  // density of pure H2SO4 [g/cc]
            auto d_crit = r_crit * 2e-7;  // diameter of critical cluster [cm]
            auto V_crit = cube(d_crit) * pi_sixth;  // volume [cc]
            auto m_crit = rho_h2so4 * V_crit;       // mass [g]
            n_crit_h2so4 = m_crit / mw_h2so4 * Na;
            n_crit_nh3 = 0;
            n_crit = n_crit_h2so4;
          }

          // Now determine the mode into which we place the nucleated particles.
          // Each particle large enough to fit into a mode directly is placed
          // into this mode--others are grown until they fit into the mode with
          // the smallest minimum diameter.
          // TODO

          Real initial_nuc_diam = 1;  // initial nucleus diameter [nm]
          Real final_nuc_diam = 3;    // final, grown nucleus diameter [nm]
          // TODO

          int nuc_mode = 0;  // TODO: Need to decide which mode to use for this
          Real mean_modal_diam = d_mean_aer_(nuc_mode);

          // Apply the correction of Kerminen and Kulmala (2002) to the
          // nucleation rate based on their growth.
          PackType rho_nuc;              // TODO
          Real gamma_h2so4 = 5.0 / 3.0;  // TODO: can we do better?
          PackType speed_h2so4 =
              gas_kinetics::molecular_speed(temp, mu_h2so4_, gamma_h2so4);
          PackType nuc_growth_rate = kerminen2002::nucleation_growth_rate(
              rho_nuc, c_h2so4, speed_h2so4, mu_h2so4_);
          const auto q_so4 =
              prognostics.interstitial_aerosols(ipop_so4_(nuc_mode), k);
          auto c_so4 =
              1e6 * conversions::number_conc_from_mmr(q_so4, mu_so4_, rho_d);
          static constexpr Real h2so4_accom_coeff = 0.65;
          auto apparent_J_factor = kerminen2002::apparent_nucleation_factor(
              c_so4, mean_modal_diam, h2so4_accom_coeff, temp, nuc_growth_rate,
              rho_nuc, initial_nuc_diam, final_nuc_diam);
          J *= apparent_J_factor;

          // Grow nucleated particles so they can fit into the desired mode.
          // TODO
        });
  }

  void set_param_(const std::string &name, Real value) override;
  void set_param_(const std::string &name, int value) override;

  // This method is for Jeff and Hui's nucleation cross-validation exercise.
  KOKKOS_INLINE_FUNCTION
  void run_skywalker_mode_(Real t, Real dt, const Prognostics &prognostics,
            const Atmosphere &atmosphere, const Diagnostics &diagnostics,
            const Tendencies &tendencies) const {
    // Compute the nucleation rate and apply it directly to the aitken mode.
    const int nk = atmosphere.temperature.extent(0);
    Kokkos::parallel_for(
        "nucleation rate (skywalker mode)", nk, KOKKOS_LAMBDA(const int k) {
          const auto temp = atmosphere.temperature(k);
          const auto press = atmosphere.pressure(k);
          const auto qv = atmosphere.vapor_mixing_ratio(k);
          const auto rho_d = gas_kinetics::air_mass_density(press, temp, qv);

          auto rel_hum = conversions::relative_humidity_from_vapor_mixing_ratio(
              qv, press, temp);

          // Determine the number concentration of H2SO4 gas [#/cc].
          PackType c_h2so4;
          if (skywalker_.c_h2so4 > 0.0) {
            c_h2so4 = skywalker_.c_h2so4;
          } else {
            const auto q_h2so4 = prognostics.gases(igas_h2so4_, k);  // mmr
            c_h2so4 = 1e6 * conversions::number_conc_from_mmr(q_h2so4, mu_h2so4_, rho_d);
          }

          // Compute the nucleation rate using our selected method.
          PackType J;       // nucleation rate [#/cc]
          PackType r_crit;  // radius of critical cluster [nm]
          PackType n_crit;  // total # of molecules in a critical cluster [#]
          PackType n_crit_h2so4, n_crit_nh3;  // numbers of gas molecules in
                                              // the critical cluser [#]
          // TODO: x_crit blows up for relative humidity = 0. Do we need to
          // TODO: handle this case gracefully?
          auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(
              c_h2so4, temp, rel_hum);
          if (nucleation_method_ == 2) { // binary nucleation
            J = vehkamaki2002::nucleation_rate(c_h2so4, temp, rel_hum, x_crit);
            n_crit = vehkamaki2002::num_critical_molecules(c_h2so4, temp,
                                                           rel_hum, x_crit);
            n_crit_h2so4 = n_crit;
            n_crit_nh3 = 0;
            r_crit = vehkamaki2002::critical_radius(x_crit, n_crit);
          } else { // ternary nucleation
            const auto q_nh3 = prognostics.gases(igas_nh3_, k);  // mmr
            auto xi_nh3 = 1e12 * conversions::vmr_from_mmr(q_nh3, mu_nh3_);
            auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                            c_h2so4, xi_nh3);
            J = exp(log_J);
            n_crit_h2so4 = merikanto2007::num_h2so4_molecules(log_J, temp,
                                                              c_h2so4, xi_nh3);
            n_crit_nh3 =
                merikanto2007::num_nh3_molecules(log_J, temp, c_h2so4, xi_nh3);
            n_crit = n_crit_h2so4 + n_crit_nh3;
            r_crit =
                merikanto2007::critical_radius(log_J, temp, c_h2so4, xi_nh3);
          }

          // Place the nucleation rate into the H2SO4 species of the Aitken mode.
          int nuc_mode = 1;

//          if (k == 0) {
//            printf("c_h2so4\tT\tRH\tx*\tJ\n");
//          }
//          printf("%g\t%g\t%g\t%g\t%g\n", c_h2so4[0], temp[0], rel_hum[0], x_crit[0], J[0]);

          PackType& dqdt = tendencies.interstitial_aerosols(ipop_so4_(nuc_mode), k);
          J.set(J < 0, 0);
          dqdt = J;
//          const auto mw_air = Constants::molec_weight_dry_air;
//          const auto mw_so4 = Constants::molec_weight_so4;
//          auto c_air = rho_d / mw_air;
//          PackType& dqndt = tendencies.interstitial_num_mix_ratios(nuc_mode, k);
//          PackType& dqdt = tendencies.interstitial_aerosols(ipop_so4_(nuc_mode), k);
//          PackType& dqgdt = tendencies.gases(igas_h2so4_, k);
//          dqndt = 1e6 * J / c_air; // convert to [#/m3]
//          dqdt = dqndt * mw_so4 / mw_air;
//          dqgdt = -dqdt;
        });
  }
};

}  // namespace haero

#endif
