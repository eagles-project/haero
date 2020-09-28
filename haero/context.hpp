#ifndef HAERO_CONTEXT_HPP
#define HAERO_CONTEXT_HPP

#include "haero/mode.hpp"
#include "haero/species.hpp"
#include "haero/chemistry.hpp"
#include "haero/parametrizations.hpp"
#include "haero/state.hpp"

namespace haero {

/// This type stores the definition of an aerosol system to be simulated.
/// A Context is available in every parameterization, and stores the modes,
/// aerosol and gas species, and chemical reactions of interest, in addition to
/// specifications for which physical processes are to be modeled. A Context
/// can also create State objects, which store all relevant data for an aerosol
/// simulation.
class Context final {
  public:

  /// Creates a Context that supports the specified parametrizations.
  /// @param [in] parametrizations the set of parametrizations (including
  ///                              implementations) supported by the resulting
  ///                              Context
  /// @param [in] aerosol_modes a list of aerosol modes supported by the Context
  /// @param [in] aerosol_species a list of aerosol species supported by the
  ///                             Context
  /// @param [in] mode_species a map that defines the association of aerosol
  ///             species with an aerosol mode. The keys in this map are names
  ///             of aerosol modes (corresponding to those found in `modes`),
  ///             and the values are lists of symbolic names of aerosol species
  ///             (supplied in `aerosol_species`) that belong to those modes.
  /// @param [in] gas_species a list of gas species supported by the Context
  /// @param [in] gas_chemistry a ChemicalMechanism representing gas chemistry.
  Context(const Parametrizations& parametrizations,
          const std::vector<Mode>& aerosol_modes,
          const std::vector<Species>& aerosol_species,
          const std::map<std::string, std::vector<std::string> >& mode_species,
          const std::vector<Species>& gas_species,
          const ChemicalMechanism& gas_chemistry = InactiveChemicalMechanism());

  /// Context objects are not deep-copyable. They should be passed by reference.
  Context(const Context&) = delete;

  /// Destructor.
  ~Context();

  /// Context objects are not assignable either.
  Context& operator=(const Context&) = delete;

  /// Returns the parametrizations associated with this Context.
  const Parametrizations& parametrizations() const;

  /// Returns the list of aerosol modes associated with this Context.
  const std::vector<Mode>& modes() const;

  /// Returns the list of aerosol species associated with this Context.
  const std::vector<Species>& aerosol_species() const;

  /// Returns the list of gas species associated with this Context.
  const std::vector<Species>& gas_species() const;

  /// Returns the chemical mechanism for gas chemistry associated with this
  /// Context.
  const ChemicalMechanism& gas_chemistry() const;

  /// Creates a new State object that can be used with this Context.
  State create_state() const;

  private:

  // Parameterizations, modes, species, chemistry...
  Parametrizations parametrizations_;
  std::vector<Mode> modes_;
  std::vector<Species> aero_species_;
  std::vector<Species> gas_species_;
  ChemicalMechanism gas_chem_;
};

}

#endif
