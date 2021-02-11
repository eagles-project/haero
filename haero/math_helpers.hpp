#ifndef HAERO_MATH_HELPERS_HPP
#define HAERO_MATH_HELPERS_HPP

#include "haero/haero_config.hpp"
#include "floating_point.hpp"

namespace haero {

template <typename T> KOKKOS_INLINE_FUNCTION
T square(T x) {
  static_assert(std::is_arithmetic<T>::value, "square: bad type");
  return x*x;
}

}
#endif