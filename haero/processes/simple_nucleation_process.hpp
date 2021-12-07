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

  /// Binary nucleation as implemented by Vehkamaki et al (2002).
  static const int binary_nucleation = 2;

  /// Ternary nucleation as implemented by Merikanto et al (2007).
  static const int ternary_nucleation = 3;

  /// Planetary boundary corrections as described in Wang et al (2008).
  static const int no_pbl_correction = 0;
  static const int first_order_pbl_correction = 1;
  static const int second_order_pbl_correction = 2;

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

  /// Nucleation method selection (default: binary_nucleation)
  int nucleation_method_;

  /// Planetary boundary layer (PBL) method selection (default: none)
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

  /// Index of nucleation mode (into which nucleated particles are placed).
  int inuc_mode_;

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
  Real mw_h2so4_, mw_nh3_;

  // Molecular weights of SO4 and NH4 aerosols.
  Real mw_so4_, mw_nh4_;

 public:
  /// Constructor
  SimpleNucleationProcess();

  /// Destructor
  KOKKOS_INLINE_FUNCTION
  ~SimpleNucleationProcess() {}

  /// Copy constructor (for transferring between host and device)
  KOKKOS_INLINE_FUNCTION
  SimpleNucleationProcess(const SimpleNucleationProcess& rhs)
      : DeviceAerosolProcess<SimpleNucleationProcess>(rhs),
        nucleation_rate_factor_(rhs.nucleation_rate_factor_),
        pbl_factor_(rhs.pbl_factor_),
        tendency_factor_(rhs.tendency_factor_),
        nucleation_method_(rhs.nucleation_method_),
        pbl_method_(rhs.pbl_method_),
        igas_h2so4_(rhs.igas_h2so4_),
        igas_nh3_(rhs.igas_nh3_),
        num_modes_(rhs.num_modes_),
        inuc_mode_(rhs.inuc_mode_),
        iaer_so4_(rhs.iaer_so4_),
        ipop_so4_(rhs.ipop_so4_),
        iaer_nh4_(rhs.iaer_nh4_),
        ipop_nh4_(rhs.ipop_nh4_),
        d_mean_aer_(rhs.d_mean_aer_),
        d_min_aer_(rhs.d_min_aer_),
        d_max_aer_(rhs.d_max_aer_) {}

 private:
  void init_(const ModalAerosolConfig& config) override;

  /// Computes the base nucleation rate.
  /// @param [in] c_h2so4 the number concentration of H2SO4 gas [#/cc]
  /// @param [in] q_nh3 the mass mixing ratio of NH3 gas [kg gas/kg dry air]
  /// @param [in] temp atmospheric temperature [K]
  /// @param [in] rel_hum atmospheric relative humidity [-]
  /// @param [in] rho_d mass density of dry air [kg/m3]
  /// @param [out] J the base nucleation rate [#/cc/s]
  /// @param [out] r_crit the radius of a critical cluster [nm]
  /// @param [out] n_crit_h2so4 the number of H2SO4 gas molecules in the cluster
  /// @param [out] n_crit_nh3 the number of NH3 gas molecules in the cluster
  KOKKOS_INLINE_FUNCTION
  void compute_nucleation_rate_(const PackType& c_h2so4, const PackType& q_nh3,
                                const PackType& temp, const PackType& rel_hum,
                                const PackType& rho_d, PackType& J,
                                PackType& r_crit, PackType& n_crit_h2so4,
                                PackType& n_crit_nh3) const {
    if (nucleation_method_ == binary_nucleation) {
      auto x_crit =
          vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, temp, rel_hum);
      J = vehkamaki2002::nucleation_rate(c_h2so4, temp, rel_hum, x_crit);
      auto n_crit =
          vehkamaki2002::num_critical_molecules(c_h2so4, temp, rel_hum, x_crit);
      n_crit_h2so4 = n_crit;
      n_crit_nh3 = 0;
      r_crit = vehkamaki2002::critical_radius(x_crit, n_crit);
    } else {
      EKAT_KERNEL_ASSERT(nucleation_method_ == ternary_nucleation);
      // Compute the molar/volume mixing ratio of NH3 gas [ppt].
      auto xi_nh3 = 1e12 * conversions::vmr_from_mmr(q_nh3, mw_nh3_);

      auto log_J =
          merikanto2007::log_nucleation_rate(temp, rel_hum, c_h2so4, xi_nh3);
      J = exp(log_J);
      n_crit_h2so4 =
          merikanto2007::num_h2so4_molecules(log_J, temp, c_h2so4, xi_nh3);
      n_crit_nh3 =
          merikanto2007::num_nh3_molecules(log_J, temp, c_h2so4, xi_nh3);
      r_crit = merikanto2007::critical_radius(log_J, temp, c_h2so4, xi_nh3);
    }
  }

  /// Applies a planetary boundary layer correction to the nucleation rate J,
  /// the critical cluster radius r_crit, and the critical number concentrations
  /// of H2SO4 and NH3 gases.
  /// @param [in] c_h2so4 the number concentration of H2SO4 gas [#/cc]
  /// @param [in] h the height of the given cell(s) [m]
  /// @param [in] pblh the planetary boundary layer height [m]
  KOKKOS_INLINE_FUNCTION
  void apply_pbl_correction_(const PackType& c_h2so4, const PackType& h,
                             const PackType& pblh, PackType& J,
                             PackType& r_crit, PackType& n_crit_h2so4,
                             PackType& n_crit_nh3) const {
    if (pbl_method_ != no_pbl_correction) {
      PackType J_pbl(0);
      const auto within_pbl = (h <= pblh);
      if (pbl_method_ == first_order_pbl_correction) {
        J_pbl.set(within_pbl,
                  wang2008::first_order_pbl_nucleation_rate(c_h2so4));
      } else if (pbl_method_ == second_order_pbl_correction) {
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
      auto d_crit = r_crit * 2e-7;        // diameter of critical cluster [cm]
      auto V_crit = cube(d_crit) * pi_sixth;  // volume [cc]
      auto m_crit = rho_h2so4 * V_crit;       // mass [g]
      n_crit_h2so4 = m_crit / mw_h2so4 * Na;
      n_crit_nh3 = 0;
    }
  }

  // Given a nucleation rate and information about a critical cluster (CC) and
  // its composition, grows particles over a time change dt to produce changes
  // in mixing ratios (MR).
  KOKKOS_INLINE_FUNCTION
  void grow_particles_(const PackType& J, // nucleation rate
                       const PackType& n_crit_h2so4, // # H2SO4 particles in CC
                       const PackType& n_crit_nh3, // # NH3 particles in CC
                       const PackType& r_crit, // radius of CC
                       const Real dt, // growth period
                       const PackType& temp, // temperature
                       const PackType& rel_hum, // relative humidity
                       const PackType& c_air, // air density(?!?!?!)
                       const PackType& accom_coeff_h2so4, // accom coefficient
                       const PackType& q_h2so4, // H2SO4 gas MR
                       const PackType& q_nh3, // NH3 gas MR
                       const PackType& c_so4, // SO4 number concentration
                       PackType& dq_h2so4, // change in H2SO4 gas MR
                       PackType& dq_nh3, // change in NH3 gas MR
                       PackType& dq_so4, // change in SO4 aerosol MR
                       PackType& dq_nh4, // change in NH4 aerosol MR
                       PackType& dq_n // change in modal number MR
                       ) const {
  }

  KOKKOS_INLINE_FUNCTION
  void run_(const TeamType& team, Real t, Real dt,
            const Prognostics& prognostics, const Atmosphere& atmosphere,
            const Diagnostics& diagnostics,
            const Tendencies& tendencies) const override {
    // Do we have any gas from which to nucleate new aerosol particles?
    if ((igas_h2so4_ == -1) or
        ((nucleation_method_ == ternary_nucleation) and (igas_nh3_ == -1))) {
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

    // Below this logarithmic nucleation rate threshold, no nucleation occurs.
    const Real log_J_cutoff = -13.82;

    // Calculate tendencies for aerosols/gases at each vertical level k.
    const int nk = atmosphere.temperature.extent(0);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nk), [=](int k) {
      const auto temp = atmosphere.temperature(k);
      const auto press = atmosphere.pressure(k);
      const auto qv = atmosphere.vapor_mixing_ratio(k);
      const auto h = atmosphere.height(k);
      const Real pblh = atmosphere.planetary_boundary_height;
      const auto rho_d = gas_kinetics::air_mass_density(press, temp, qv);
      auto rel_hum = conversions::relative_humidity_from_vapor_mixing_ratio(
          qv, press, temp);
      // Mass mixing ratios (mmr, [kg gas/kg dry air]) for H2SO4 and NH3 gas.
      const auto q_h2so4 = prognostics.gases(igas_h2so4_, k);
      const PackType q_nh3 =
          (igas_nh3_ >= 0) ? prognostics.gases(igas_nh3_, k) : 0;

      // Find the number concentration for H2SO4, which is used by our
      // parameterizations.
      auto c_h2so4 =
          1e6 * conversions::number_conc_from_mmr(q_h2so4, mw_h2so4_, rho_d);

      // Compute the base rate of nucleation using our selected method.
      PackType J;                         // nucleation rate [#/cc]
      PackType r_crit;                    // radius of critical cluster [nm]
      PackType n_crit_h2so4, n_crit_nh3;  // # gas molecules in the cluser [#]
      compute_nucleation_rate_(c_h2so4, q_nh3, temp, rel_hum, rho_d, J, r_crit,
                               n_crit_h2so4, n_crit_nh3);

      // Apply a correction for the planetary boundary layer if requested.
      apply_pbl_correction_(c_h2so4, h, pblh, J, r_crit, n_crit_h2so4,
                            n_crit_nh3);
      // Apply the cutoff threshold to J.
      J.set((log(J) <= log_J_cutoff), 0);

      // Grow the nucleated particles, applying the growth factor to J.
      int p_so4 = ipop_so4_(inuc_mode_);  // population index for nucleated SO4.
      const auto q_so4 = prognostics.interstitial_aerosols(p_so4, k);
      auto c_so4 =
          1e6 * conversions::number_conc_from_mmr(q_so4, mw_so4_, rho_d);
      PackType dq_h2so4, dq_nh3, dq_so4, dq_nh4, dq_n; // mixing ratio changes
      grow_particles_(J, n_crit_h2so4, n_crit_nh3, r_crit, dt, temp, rel_hum,
          c_air, accom_coeff_h2so4, q_h2so4, q_nh3, c_h2so4, c_so4,
          dq_h2so4, dq_nh3, dq_so4, dq_nh4, dq_n);

      // Compute tendencies. All nucleated particles are placed into the
      // selected nucleation mode.

      // Number mixing ratio tendency for nucleation mode.
      dq_n *= 1e3; // convert number MR from #/mol-air to #/kmol-air
      auto dq_n_dt = dq_n / dt;
      tendencies.interstitial_num_mix_ratios(inuc_mode_, k) = dq_n_dt;

      // Compute the mass mixing ratio tendencies.
      auto mass_so4 = dq_so4 * mw_so4_; // mass xferred to SO4
      auto mass_nh4 = dq_nh4 * mw_nh4_; // mass xferred to NH4
      auto nuc_mass = mass_so4 + mass_nh4; // total nucleated mass
      // fraction of mass xferred to SO4
      auto mass_frac_so4 = max(mass_so4, 1e-35) / max(nuc_mass, 1e-35);
      auto dm_dt = max(0, nuc_mass/dt); // nucleated mass rate

      // =====================================================
      // "Various adjustments to keep the solution reasonable"
      // =====================================================

      // Apply particle size constraints to nucleated particles.
      static const Real pi = Constants::pi;
      Real M = rho_so4 * pi / 6;
      auto mass1p_lo = M * cube(d_min_aer_(inuc_mode_));
      auto mass1p_hi = M * cube(d_max_aer_(inuc_mode_));
      auto mass1p = dm_dt / dq_n_dt;
      dq_n_dt.set((mass1p < mass1p_lo), dm_dt/mass1p_lo);
      dm_dt.set((mass1p > mass1p_hi), dq_n_dt*mass1p_hi);

      // "apply adjustment factor to avoid unrealistically high aitken number
      //  concentrations in mid and upper troposphere"
      // The adjustment factors look like they're unity at the moment, so let's
      // skip this step.

      // Zero out rates too small to consider.
      auto dq_n_dt_lt_100 = (dq_n_dt < 100);
      dq_n_dt.set(dq_n_dt_lt_100, 0);
      dm_dt.set(dq_n_dt_lt_100, 0);

      // Finally, diagnose the aerosol tendencies.
      tendencies.interstitial_aerosols(p_so4, k) =
        dm_dt * mass_frac_so4 / mw_so4_;
      tendencies.interstitial_aerosols(p_so4, k) =
        dm_dt * (1.0 - mass_frac_so4) / mw_so4_;
    });
  }

  using AerosolProcess::set_param_;
  void set_param_(const std::string& name, Real value) override;
  void set_param_(const std::string& name, int value) override;
};

}  // namespace haero

#endif
