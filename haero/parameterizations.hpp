#ifndef HAERO_PARAMETERIZATIONS_HPP
#define HAERO_PARAMETERIZATIONS_HPP

namespace haero {

/// @struct Parameterizations
/// This type specifies aerosol parameterizations available to a Haero
/// simulation. A set of parameterizations is a set of physical processes that
/// are approximately represented according to a specific algorithm. Each
/// process can be selected from a set of available algorithms, or can be
/// disabled to allow for testing.
struct Parameterizations final {

  /// Available process models for aerosol activation, in which interstitial
  /// aerosol particles are coated by water and converted to cloud-borne
  /// aerosol particles.
  enum Activation {
    NoActivation // no activation  model
  };
  Activation activation;

  /// Available process models for cloude-borne wet removal, in which
  /// cloude-borne aerosol particles are TODO
  enum CloudBorneWetRemoval {
    NoCloudBorneWetRemoval // no interstitial wet removal model
  };
  CloudBorneWetRemoval cloudborne_wet_removal;

  /// Available process models for coagulation, in which interstitial aerosol
  /// particles combine to form larger interstitial particles.
  enum Coagulation {
    NoCoagulation // no coagulation model
  };
  /// The selected coagulation model
  Coagulation coagulation;

  /// Available process models for condensation, in which TODO
  enum Condensation {
    NoCondensation // no condensation model
  };
  /// The selected condensation model
  Condensation condensation;

  /// Available process models for dry deposition, in which interstitial
  /// aerosol particles are deposited on water droplets.
  enum DryDeposition {
    NoDryDeposition // no dry deposition model
  };

  /// Available process models for emissions, in which interstitial aerosol
  /// particles are added from an external source.
  enum Emissions {
    NoEmissions // no emissions model
  };

  /// Available process models for intersitial wet removal, in which
  /// interstitial aerosol particles are caught by falling rain.
  enum InterstitialWetRemoval {
    NoInterstitialWetRemoval // no interstitial wet removal model
  };
  InterstitialWetRemoval interstitial_wet_removal;

  /// Available process models for nucleation, in which interstitial aerosol
  /// mass is formed by gathering condensed vapor from surrounding gas.
  enum Nucleation {
    NoNucleation // no nucleation model
  };
  /// The selected nucleation model
  Nucleation nucleation;

  /// Available process models for resuspension, in which cloud-borne aerosol
  /// particles are converted to interstitial aerosols.
  enum Resuspension {
    NoResuspension // no resuspension model
  };
  /// The selected resuspension model
  Resuspension resuspension;

  /// Available process models for water uptake, in which water molecules are
  /// captured by interstitial aerosol particles.
  enum WaterUptake {
    NoWaterUptake // no water uptake model
  };
  /// The selected water uptake model.
  WaterUptake water_uptake;

};

}

#endif
