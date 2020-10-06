#include "phys_column.hpp"
#include "haero/physical_constants.hpp"

namespace haero {
namespace driver {

phys_column::phys_column(const int& nlev) {
  npack_mid = nlev / HAERO_PACK_SIZE;
  npack_interface = (nlev+1) / HAERO_PACK_SIZE;
}

} // namespace driver
}// namespace haero
