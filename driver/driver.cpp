#include "standalone/driver.hpp"

// This is our all-C++ driver. It's only built if/when the driver is compiled
// with the -DMAMBOX_CXX_DRIVER flag.

#ifdef MAMBOX_CXX_DRIVER

namespace {

// Here's a home-made "physics buffer" that demonstrates how one can store
// fields with Kokkos Views.
std::map<std::string, Kokkos::View<Real> > pbuf;

// Writes output to a NetCDF file with the given path.
void write_netcdf(const std::map<std::string, std::vector<Real> >& pbuf)
{
}

} // anonymous namespace

void mambox_driver(const Sim_input_data& data)
{
  // Our own "physics buffer" is just a set of named fields,
  // housed in a C++ map. We use vectors for storage and Kokkos
  // unmanaged views for multi-dimensional array access.
  std::map<std::string, std::vector<Real> > pbuf;
}

#endif
