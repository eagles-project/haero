#include "driver.hpp"

namespace {

} // anonymous namespace

namespace haero {

void haero_driver(const std::vector<SimulationInput>& ensemble)
{
  // Our own "physics buffer" is just a set of named fields,
  // housed in a C++ map. We use vectors for storage and Kokkos
  // unmanaged views for multi-dimensional array access.
  std::map<std::string, std::vector<Real> > pbuf;
}

} // namespace haero

