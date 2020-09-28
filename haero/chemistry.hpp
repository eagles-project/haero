#ifndef HAERO_CHEMISTRY_HPP
#define HAERO_CHEMISTRY_HPP

namespace haero {

/// This type represents a chemical reaction.
class Reaction final {
  // Just a placeholder so far...
};

/// This interface represents a chemical mechanism that can equilibrate systems
/// and / model chemical kinetics where needed. Subclasses of this interface
/// can provide different solvers to implement the mechanism.
class ChemicalMechanism {
  // Just a placeholder so far...
};

/// This subclass of ChemicalMechanism provides an empty, do-nothing chemical
/// mechanism.
class InactiveChemicalMechanism: public ChemicalMechanism {
  public:

  /// Creates do-nothing chemical mechanism.
  InactiveChemicalMechanism();
};

}

#endif
