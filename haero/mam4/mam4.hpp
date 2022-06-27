#ifndef HAERO_MAM4_HPP
#define HAERO_MAM4_HPP

// This header makes available all MAM4 processes.
#include <haero/aero_process.hpp>
#include <haero/mam4/nucleation_impl.hpp>

namespace haero {
namespace mam4 {

using NucleationProcess = AeroProcess<AeroConfig, NucleationImpl>;

} // namespace mam4
} // namespace haero

#endif
