#ifndef HAERO_SELECT_PROCESS_HPP
#define HAERO_SELECT_PROCESS_HPP

#include "haero/aero_process.hpp"
#include "haero/parameterizations.hpp"

namespace haero {

/// Given an aerosol process type and a set of selected parameterizations,
/// this function creates and returns a pointer to a newly-allocated
/// AeroProcess instance. The implementation of this function must be updated
/// whenever a new AeroProcess implementation is made available.
/// @param [in] type The type of aerosol process to be selected.
/// @param [in] selections A struct containing selected aerosol parametrizations
///                        and their implementations.
/// @returns A pointer to a newly allocated process reflecting the given
///          selections.
AeroProcess* select_process(AeroProcessType type,
                            const Parameterizations& selections);

}

#endif
