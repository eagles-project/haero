#ifndef HAERO_SELECTED_PROCESSES_HPP
#define HAERO_SELECTED_PROCESSES_HPP

#include "haero/aerosol_process.hpp"

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
    /// MAM legacy C++ nucleation process
    MAMNucleation,
#if HAERO_FORTRAN
    /// MAM legacy Fortran nucleation process
    MAMFNucleation,
#endif
    /// Simplified nucleation process
    SimpleNucleation,
    /// No nucleation process
    NoNucleation
  };
  /// The selected nucleation model
  Nucleation nucleation;

  /// Available process models for gas aerosol exchange, in which interstitial
  /// aerosol mass is formed by gathering condensed vapor from surrounding gas.
  enum GasAerosolExchange {
    /// MAM legacy C++ gas aerosol exchange process
    MAMGasAerosolExchange,
#if HAERO_FORTRAN
    /// MAM legacy Fortran gas aerosol exchange process
    MAMFGasAerosolExchange,
#endif
    /// No gas aerosol exchange process
    NoGasAerosolExchange
  };
  /// The selected gas aerosol exchange model
  GasAerosolExchange gasaerosolexchange;

  /// Available process models for calcsize, in which we compute dry diameter of
  /// the the particles and do the inter-mode mass and number transfer based on
  /// the new size (for both interstitial and cloud borne aerosol)
  enum Calcsize {
    /// MAM legacy C++ calcsize process
    MAMCalcsize,
#if HAERO_FORTRAN
    /// MAM legacy Fortran calcsize process
    MAMFCalcsize,
#endif
    /// No calcsize process
    NoCalcsize
  };
  /// The selected calcsize model
  Calcsize calcsize;

  enum Rename {
    /// MAM legacy C++ Rename process
    MAMRename,
#if HAERO_FORTRAN
    /// MAM legacy Fortran rename process
    MAMFRename,
#endif
    /// No Rename process
    NoRename,
  };

  /// The selected rename model
  Rename rename;

  /// Available process models for resuspension, in which cloud-borne aerosol
  /// particles are converted to interstitial aerosols.
  enum Resuspension {
    /// No resuspension process
    NoResuspension
  };
  /// The selected resuspension model
  Resuspension resuspension;
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
AerosolProcess* select_aerosol_process(AerosolProcessType type,
                                       const SelectedProcesses& selections);

}  // namespace haero

#endif
