#ifndef HAERO_SELECTED_PROCESSES_HPP
#define HAERO_SELECTED_PROCESSES_HPP

#include "haero/process.hpp"

namespace haero {

/// @struct SelectedProcesses
/// This type specifies a set of selected aerosol processes from those available
/// to a Haero simulation. Each process is selected from a set of available
/// algorithms, or can be disabled to allow for testing. This struct must be
/// updated when a new process implementation is made available.
struct SelectedProcesses final {

  /// Default constructor (no processes selected).
  SelectedProcesses();

  /// Available process models for aerosol activation, in which interstitial
  /// aerosol particles are coated by water and converted to cloud-borne
  /// aerosol particles.
  enum Activation {
    /// No activation process
    NoActivation
  };
  Activation activation;

  /// Available process models for cloude-borne wet removal, in which
  /// cloude-borne aerosol particles are TODO
  enum CloudBorneWetRemoval {
    /// No interstitial wet removal process
    NoCloudBorneWetRemoval
  };
  CloudBorneWetRemoval cloudborne_wet_removal;

  /// Available process models for coagulation, in which interstitial aerosol
  /// particles combine to form larger interstitial particles.
  enum Coagulation {
    // No coagulation process
    NoCoagulation
  };
  /// The selected coagulation model
  Coagulation coagulation;

  /// Available process models for condensation, in which TODO
  enum Condensation {
    /// No condensation process
    NoCondensation
  };
  /// The selected condensation model
  Condensation condensation;

  /// Available process models for dry deposition, in which interstitial
  /// aerosol particles are deposited on water droplets.
  enum DryDeposition {
    /// No dry deposition process
    NoDryDeposition
  };
  /// The selected dry deposition model
  DryDeposition dry_deposition;

  /// Available process models for emissions, in which interstitial aerosol
  /// particles are added from an external source.
  enum Emissions {
    // No emissions process
    NoEmissions
  };
  /// The selected emissions model
  Emissions emissions;

  /// Available process models for intersitial wet removal, in which
  /// interstitial aerosol particles are caught by falling rain.
  enum InterstitialWetRemoval {
    /// No interstitial wet removal process
    NoInterstitialWetRemoval
  };
  /// The selected interstitial wet removal model
  InterstitialWetRemoval interstitial_wet_removal;

  /// Available process models for nucleation, in which interstitial aerosol
  /// mass is formed by gathering condensed vapor from surrounding gas.
  enum Nucleation {
    /// MAM4 legacy Fortran nucleation process
    MAM4FNucleation,
    /// No nucleation process
    NoNucleation
  };
  /// The selected nucleation model
  Nucleation nucleation;

  /// Available process models for resuspension, in which cloud-borne aerosol
  /// particles are converted to interstitial aerosols.
  enum Resuspension {
    /// No resuspension process
    NoResuspension
  };
  /// The selected resuspension model
  Resuspension resuspension;

  /// Available process models for water uptake, in which water molecules are
  /// captured by interstitial aerosol particles.
  enum WaterUptake {
    /// no water uptake model
    NoWaterUptake,
    /// Water uptake module from MAM4 in Fortran
    MAM4Fortran,
    /// Water uptake module from MAM4 in Fortran using bisection for Kohler solve
    MAM4KohlerBisectionFortran
  };
  /// The selected water uptake model.
  WaterUptake water_uptake;

};

/// Given a (prognostic) aerosol process type and a set of selected
/// parameterizations, this function creates and returns a pointer to a
/// newly-allocated PrognosticProcess instance representing a specific
/// parameterization (or implementation). The implementation of this function
/// must be updated whenever a new PrognosticProcess implementation is made
/// available.
/// @param [in] type The type of aerosol process to be selected.
/// @param [in] selections A struct containing selected aerosol processes
///                        and their implementations.
/// @returns A pointer to a newly allocated process reflecting the given
///          selections.
PrognosticProcess* select_prognostic_process(ProcessType type,
                                             const SelectedProcesses& selections);

/// Given a (diagnostic) aerosol process type and a set of selected
/// parameterizations, this function creates and returns a pointer to a
/// newly-allocated DiagnosticProcess instance representing a specific
/// parameterization (or implementation). The implementation of this function
/// must be updated whenever a new DiagnosticProcess implementation is made
/// available.
/// @param [in] type The type of aerosol process to be selected.
/// @param [in] selections A struct containing selected aerosol processes
///                        and their implementations.
/// @returns A pointer to a newly allocated process reflecting the given
///          selections.
DiagnosticProcess* select_diagnostic_process(ProcessType type,
                                             const SelectedProcesses& selections);

} // end haero namespace

#endif
