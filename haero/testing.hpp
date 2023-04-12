// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HAERO_TESTING_HPP
#define HAERO_TESTING_HPP

#include <haero/atmosphere.hpp>

namespace haero {

// The testing namespace contains tools that are useful only in testing
// environments.
namespace testing {

/// creates an Atmosphere object that stores a column of data with the given
/// number of vertical levels and the given planetary boundary height
/// @param [in] num_levels the number of vertical levels per column stored by
///                        the state
/// @param [in] pblh The column-specific planetary boundary height [m],
///                  computed by the host model
Atmosphere create_atmosphere(int num_levels, Real pblh);

/// Creates a standalone ColumnView that uses resources allocated by a memory
/// pool.
ColumnView create_column_view(int num_levels);

/// Call this at the end of a testing session to delete all ColumnViews
/// allocated by create_column_view. This is called by Haero's implementation
/// of ekat_finalize_test_session, which is called automatically at the end of
/// each Catch2-powered unit test.
void finalize();

} // end namespace testing

} // end namespace haero

#endif
