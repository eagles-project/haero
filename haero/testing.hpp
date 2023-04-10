// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#ifndef HAERO_TESTING_HPP
#define HAERO_TESTING_HPP

#include <haero/haero.hpp>

namespace haero {

// The testing namespace contains tools that are useful only in testing
// environments.
namespace testing {

/// Creates a standalone ColumnView that uses resources allocated by a memory
/// pool.
ColumnView create_column_view(int num_levels);

} // end namespace testing

} // end namespace haero

#endif
