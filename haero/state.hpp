#ifndef HAERO_STATE_HPP
#define HAERO_STATE_HPP

#include <map>
#include <vector>
#include "haero/mode.hpp"
#include "haero/species.hpp"

namespace haero {

/// This type stores state information for an aerosol system. It stores
/// * mixing ratios of for each mode-specific interstitial and cloud-borne
///   aerosol species
/// * interstitial and cloud-borne number concentrations for each aerosol mode
/// * mole fractions for gas species
class State final {
  public:

  /// Creates a State using information supplied. Usually, you don't want to
  /// call this constructor explicity--create a State using a Context object
  /// instead.
  State();

  /// Destructor.
  ~State();

  private:

};

}

#endif
