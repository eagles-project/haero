#ifndef HAERO_SIMPLE_NUCLEATION_PROCESS_HPP
#define HAERO_SIMPLE_NUCLEATION_PROCESS_HPP

#include <iomanip>

#include "haero/aerosol_process.hpp"
#include "haero/conversions.hpp"
#include "haero/physical_constants.hpp"
#include "haero/processes/kerminen2002.hpp"
#include "haero/processes/merikanto2007.hpp"
#include "haero/processes/vehkamaki2002.hpp"
#include "haero/processes/wang2008.hpp"

namespace haero {

/// @class SimpleNucleationProcess
/// This aerosol process implements homogeneous binary/ternary nucleation
/// involving sulfuric acid and methane gases. It is based on classical
/// nucleation theory as parameterized by Vehkamaki et al (2002) and
/// Merikanto et al (2007).
class SimpleNucleationProcess final : public AerosolProcess {
  //----------------------------------------------------------
  //                  Adjustable parameters
  //----------------------------------------------------------
  // All of the named parameters below can be set using the
  // set_param() methods.
  //----------------------------------------------------------

  /// Adjustment factor applied to nucleation rate (default: 1)
  Real nucleation_rate_factor;

  /// Adjustment factor applied for planetary boundary layer (default: 1)
  Real pbl_factor;

  /// Adjustment factor applied to tendency for nucleated species (default: 1)
  Real tendency_factor;

  /// Nucleation method selection (default: 2)
  /// 2 = Vehkamaki et al (2002) binary nucleation
  /// 3 = Merikanto el al (2007) ternary nucleation
  int nucleation_method;

  /// The name of the aerosol mode into which nucleated aerosols are placed
  /// (default: "aitken")
  std::string nucleation_mode;

  /// Planetary boundary layer (PBL) method selection (default: 0)
  /// 0 = no PBL adjustment
  /// 1 = first-order
  /// 2 = second-order
  int pbl_method;

  //----------------------------------------------------------
  //                       Bookkeeping
  //----------------------------------------------------------
  // These metadata allow the process to work with the given
  // modal aerosol configuration without having to look stuff
  // up all the time.

  // Index of the nucleation mode
  int imode;

  // Index of H2SO4 gas
  int igas_h2so4;

  // Index of NH3 gas
  int igas_nh3;

  // Species and population indices of SO4 aerosol within the nucleation mode
  int iaer_so4, ipop_so4;

  // Species and population indices of NH4 aerosol within the nucleation mode
  int iaer_nh4, ipop_nh4;

  // The geometric mean particle diameters for all aerosol modes
  view_1d_scalar_type d_mean_aer;

  /// The minimum particle diameters for all aerosol modes
  view_1d_scalar_type d_min_aer;

  /// The maximum particle diameters for all aerosol modes
  view_1d_scalar_type d_max_aer;

  // Molecular weights of H2SO4 and NH3 gases.
  Real mu_h2so4, mu_nh3;

  // Molecular weights of SO4 and NH4 aerosols.
  Real mu_so4, mu_nh4;

 public:
  /// Constructor
  SimpleNucleationProcess();

  /// Destructor
  KOKKOS_INLINE_FUNCTION
  ~SimpleNucleationProcess() {}

  /// Copy constructor (for transferring between host and device)
  KOKKOS_INLINE_FUNCTION
  SimpleNucleationProcess(const SimpleNucleationProcess &rhs)
      : AerosolProcess(rhs),
        nucleation_rate_factor(rhs.nucleation_rate_factor),
        pbl_factor(rhs.pbl_factor),
        tendency_factor(rhs.tendency_factor),
        nucleation_method(rhs.nucleation_method),
        nucleation_mode("aitken"),
        pbl_method(rhs.pbl_method),
        imode(rhs.imode),
        igas_h2so4(rhs.igas_h2so4),
        igas_nh3(rhs.igas_nh3),
        iaer_so4(rhs.iaer_so4),
        iaer_nh4(rhs.iaer_nh4),
        d_mean_aer(rhs.d_mean_aer),
        d_min_aer(rhs.d_min_aer),
        d_max_aer(rhs.d_max_aer) {}

  /// not assignable
  AerosolProcess &operator=(const SimpleNucleationProcess &) = delete;

  void init(const ModalAerosolConfig &config) override;

  KOKKOS_FUNCTION
  void run(const ModalAerosolConfig &config, Real t, Real dt,
           const Prognostics &prognostics, const Atmosphere &atmosphere,
           const Diagnostics &diagnostics,
           Tendencies &tendencies) const override {
    // Do we have any gas from which to nucleate new aerosol particles?
    if ((igas_h2so4 == -1) or ((nucleation_method == 3) and (igas_nh3 == -1))) {
      return;
    }

    // Do we track the aerosols nucleated from our gases?
    if ((iaer_so4 == -1) and (iaer_nh4 == -1)) {
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
          const auto q_h2so4 = prognostics.gases(igas_h2so4, k);  // mmr
          auto c_h2so4 =  // number concentration of H2SO4 [#/cc]
              1e6 * conversions::number_conc_from_mmr(q_h2so4, mu_h2so4, rho_d);

          // Compute the base rate of nucleation using our selected method.
          PackType J;       // nucleation rate [#/cc]
          PackType r_crit;  // radius of critical cluster [nm]
          PackType n_crit;  // number of molecules in a critical cluster [#]
          if (nucleation_method == 2) {  // binary nucleation
            auto x_crit = vehkamaki2002::h2so4_critical_mole_fraction(
                c_h2so4, temp, rel_hum);
            J = vehkamaki2002::nucleation_rate(c_h2so4, temp, rel_hum, x_crit);
            n_crit = vehkamaki2002::num_critical_molecules(c_h2so4, temp,
                                                           rel_hum, x_crit);
            r_crit = vehkamaki2002::critical_radius(x_crit, n_crit);
          } else {
            EKAT_KERNEL_ASSERT(nucleation_method == 3);
            // Compute the molar/volume mixing ratio of NH3 gas [ppt].
            const auto q_nh3 = prognostics.gases(igas_nh3, k);  // mmr
            auto xi_nh3 = 1e12 * conversions::vmr_from_mmr(q_nh3, mu_nh3);

            auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                            c_h2so4, xi_nh3);
            J = exp(log_J);
            n_crit = merikanto2007::num_critical_molecules(log_J, temp, c_h2so4,
                                                           xi_nh3);
            r_crit =
                merikanto2007::critical_radius(log_J, temp, c_h2so4, xi_nh3);
          }

          // Apply a correction for the planetary boundary layer if requested.
          if (pbl_method > 0) {
            PackType J_pbl(0);
            const auto within_pbl = (h <= atmosphere.planetary_boundary_height);
            if (pbl_method == 1) {
              J_pbl.set(within_pbl,
                        wang2008::first_order_pbl_nucleation_rate(c_h2so4));
            } else if (pbl_method == 2) {
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
            // mass and radius.
            // TODO
          }

          Real initial_nuc_diam = 1;  // initial nucleus diameter [nm]
          Real final_nuc_diam = 3;    // final, grown nucleus diameter [nm]
          // TODO

          Real mean_modal_diam = d_mean_aer(imode);

          // Apply the correction of Kerminen and Kulmala (2002) to the
          // nucleation rate.
          PackType rho_nuc;              // TODO
          Real gamma_h2so4 = 5.0 / 3.0;  // TODO: can we do better?
          PackType speed_h2so4 =
              gas_kinetics::molecular_speed(temp, mu_h2so4, gamma_h2so4);
          PackType nuc_growth_rate = kerminen2002::nucleation_growth_rate(
              rho_nuc, c_h2so4, speed_h2so4, mu_h2so4);
          const auto q_so4 = prognostics.interstitial_aerosols(ipop_so4, k);
          auto c_so4 =
              1e6 * conversions::number_conc_from_mmr(q_so4, mu_so4, rho_d);
          static constexpr Real h2so4_accom_coeff = 0.65;
          auto apparent_J_factor = kerminen2002::apparent_nucleation_factor(
              c_so4, mean_modal_diam, h2so4_accom_coeff, temp, nuc_growth_rate,
              rho_nuc, initial_nuc_diam, final_nuc_diam);
          J *= apparent_J_factor;

          // Grow nucleated particles so they can fit into the desired mode.
          // TODO
        });
  }

  void set_param(const std::string &name, Real value) override;
  void set_param(const std::string &name, int value) override;
  void set_param(const std::string &name, const std::string &value) override;
};

}  // namespace haero

#endif
