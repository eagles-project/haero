#ifndef HAERO_HPP
#define HAERO_HPP

#include "haero/haero_config.hpp"

namespace haero {

// Returns haero's version string.
const char* version();

// Returns haero's git revision hash, or "unknown" if not found.
const char* revision();

// Returns true iff this build has changes that weren't committed to the git
// repo.
bool has_uncommitted_changes();

} // end namespace haero

#endif
