#ifndef HAERO_PARAMETRIZATIONS_HPP
#define HAERO_PARAMETRIZATIONS_HPP

namespace haero {

/// This type specifies aerosol parametrizations available to a Haero
/// simulation. A set of parametrizations is a set of physical processes that
/// are approximately represented according to a specific algorithm. Each
/// process can be selected from a set of available algorithms, or can be
/// disabled to allow for testing.
struct Parametrizations final {

  /// Available process models for aerosol microphysics
  enum AerosolMicrophysics {
    NoAerosolMicrophysics // no aerosol microphysics model
  };
  AerosolMicrophysics aerosol_microphysics;

  /// Available process models for coagulation
  enum Coagulation {
    NoCoagulation // no coagulation model
  };
  /// The selected coagulation model
  Coagulation coagulation;

  /// Available process models for condensation
  enum Condensation {
    NoCondensation // no condensation model
  };
  /// The selected condensation model
  Condensation condensation;

  /// Available process models for gas-aerosol exchange
  enum GasAerosolExchange {
    NoGasAerosolExchange // no gas-aerosol exchange model
  };
  /// The selected gas-aerosol exchange model
  GasAerosolExchange gas_aerosol_exchange;

  /// Available process models for nucleation
  enum Nucleation {
    NoNucleation // no nucleation model
  };
  /// The selected nucleation model
  Nucleation nucleation;

};

}

#endif
