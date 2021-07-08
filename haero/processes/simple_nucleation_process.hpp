#ifndef HAERO_SIMPLE_NUCLEATION_PROCESS_HPP
#define HAERO_SIMPLE_NUCLEATION_PROCESS_HPP

#include <iomanip>

#include "haero/aerosol_process.hpp"
#include "haero/physical_constants.hpp"

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

  // Index of SO4 aerosol within the nucleation mode
  int iaer_so4;

  // The geometric mean particle diameters for all aerosol modes
  view_1d_scalar_type d_mean_aer;

  /// The minimum particle diameters for all aerosol modes
  view_1d_scalar_type d_min_aer;

  /// The maximum particle diameters for all aerosol modes
  view_1d_scalar_type d_max_aer;

 public:
  /// Constructor
  SimpleNucleationProcess();

  /// Destructor
  KOKKOS_INLINE_FUNCTION
  virtual ~SimpleNucleationProcess() {}

  /// Default copy constructor. For use in moving host instance to device.
  KOKKOS_INLINE_FUNCTION
  SimpleNucleationProcess(const SimpleNucleationProcess &rhs)
      : AerosolProcess(rhs),
        nucleation_rate_factor(rhs.nucleation_rate_factor),
        pbl_factor(rhs.pbl_factor),
        tendency_factor(rhs.tendency_factor),
        nucleation_method(rhs.nucleation_method),
        nucleation_mode("aitken"),
        pbl_method(rhs.pbl_method),
        igas_h2so4(rhs.igas_h2so4),
        igas_nh3(rhs.igas_nh3),
        iaer_so4(rhs.iaer_so4),
        d_mean_aer("mean particle diameters", 0),
        d_min_aer("minimum particle diameters", 0),
        d_max_aer("maximum particle diameters", 0) {}

  /// not assignable
  AerosolProcess &operator=(const SimpleNucleationProcess &) = delete;

  void init(const ModalAerosolConfig &config) override;

  KOKKOS_FUNCTION
  void run(const ModalAerosolConfig &config, Real t, Real dt,
           const Prognostics &prognostics, const Atmosphere &atmosphere,
           const Diagnostics &diagnostics,
           Tendencies &tendencies) const override {}

  void set_param(const std::string &name, Real value) override;
  void set_param(const std::string &name, int value) override;
  void set_param(const std::string &name, const std::string &value) override;
};

}  // namespace haero

#endif
