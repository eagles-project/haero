#ifndef HAERO_VIEW_HELPERS_HPP
#define HAERO_VIEW_HELPERS_HPP

#include "haero/haero_config.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include <vector>

namespace haero {

#define VIEW_REAL_TYPE_IS_SP \
typename std::enable_if<std::is_same<typename ekat::ScalarTraits<typename ViewType::value_type>::scalar_type, float>::value,void>::type

#define VIEW_REAL_TYPE_IS_DP \
typename std::enable_if<std::is_same<typename ekat::ScalarTraits<typename ViewType::value_type>::scalar_type, double>::value,void>::type


} // namespace haero
#endif
