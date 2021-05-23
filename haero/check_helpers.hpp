#ifndef HAERO_CHECK_HELPERS_HPP
#define HAERO_CHECK_HELPERS_HPP

#include "haero/haero.hpp"

namespace haero {

/** @brief Boolean helper class for performing "Checks."

  Primarily, these checks are used inside our various *ASSERT* and *REQUIRE*
  statements, when we have to support both POD and ekat::Packs.  See
  math_helpers.hpp for examples.

  Before you add to this class, be sure to also check floating_point.hpp.  In
  that file we define a very similar set of boolean functions; however, those
  specifically account for floating point roundoff error.  These do not.
*/
template <typename ScalarType>
struct Check {
  KOKKOS_INLINE_FUNCTION
  static bool is_negative(const ScalarType& x) { return (x < 0); }

  KOKKOS_INLINE_FUNCTION
  static bool is_positive(const ScalarType& x) { return (x > 0); }
};

/** @brief Boolean helper class for performing "Checks."

  Note that the requirements for reducing an ekat::Mask to a single boolean
  are not defined --- users must specify whether the required condition holds
  for _all_ values of the pack, or _any_ value.
*/
template <typename ScalarType>
struct Check<ekat::Pack<ScalarType, HAERO_PACK_SIZE>> {
  using ValueType = ekat::Pack<ScalarType, HAERO_PACK_SIZE>;

  KOKKOS_INLINE_FUNCTION
  static bool is_negative(const ValueType& x) { return (x < 0).all(); }

  KOKKOS_INLINE_FUNCTION
  static bool is_positive(const ValueType& x) { return (x > 0).all(); }
};

}  // namespace haero
#endif
