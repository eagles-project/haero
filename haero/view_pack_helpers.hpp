#ifndef HAERO_VIEW_HELPERS_HPP
#define HAERO_VIEW_HELPERS_HPP

#include "haero/haero_config.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/ekat_scalar_traits.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include <vector>

namespace haero {

//---------------------- Compile-time stuff ---------------//

/// Kokkos device definitions
using kokkos_device_type = ekat::KokkosTypes<ekat::DefaultDevice>;
using kokkos_host_type   = ekat::KokkosTypes<ekat::HostDevice>;

/// Ekat Pack definitions
using real_pack_type = ekat::Pack<Real, HAERO_PACK_SIZE>;
using pack_mask_type = ekat::Mask<HAERO_PACK_SIZE>;
using pack_info = ekat::PackInfo<HAERO_PACK_SIZE>;

/// View definitions
using view_1d_scalar_type = kokkos_device_type::view_1d<Real>;
using view_1d_pack_type = kokkos_device_type::view_1d<real_pack_type>;
using mask_view_1d_type = kokkos_device_type::view_1d<pack_mask_type>;

using view_2d_scalar_type = kokkos_device_type::view_2d<Real>;
using view_2d_pack_type = kokkos_device_type::view_2d<real_pack_type>;
using mask_view_2d_type = kokkos_device_type::view_2d<pack_mask_type>;

/// Compile-time functions
#define VIEW_REAL_TYPE_IS_SP \
typename std::enable_if<std::is_same<typename ekat::ScalarTraits<typename ViewType::value_type>::scalar_type, float>::value,void>::type

#define VIEW_REAL_TYPE_IS_DP \
typename std::enable_if<std::is_same<typename ekat::ScalarTraits<typename ViewType::value_type>::scalar_type, double>::value,void>::type

//------------------------ Run-time stuff ---------------//

/// fwd decl
template <typename VT>
typename std::enable_if<!ekat::ScalarTraits<typename VT::value_type>::is_simd, std::vector<Real>>::type
view1d_to_vector_impl(const VT& v);

/// fwd decl
template <typename VT>
typename std::enable_if<ekat::ScalarTraits<typename VT::value_type>::is_simd, std::vector<Real>>::type
view1d_to_vector_impl(const VT& v);

/** @brief Convert a std::vector<Real> to Kokkos::View<Real*>.

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_1d_scalar_type vector_to_basic_1dview(const std::vector<Real>& vector, const std::string& view_name);


/** @brief Convert a std::vector<Real> to Kokkos::View<Pack<Real,HAERO_PACK_SIZE>*>

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_1d_pack_type vector_to_packed_1dview(const std::vector<Real>& vector, const std::string& view_name);

template <typename ViewType>
std::vector<Real> view1d_to_vector(const ViewType& v) {return view1d_to_vector_impl<ViewType>(v); }

/** @brief Convert a 2d set of vectors, (`std::vector<std::vector<Real>>`) to `Kokkos::View<Real**>`

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_2d_scalar_type vector_to_basic_2dview(const std::vector<std::vector<Real>>& vectors, const std::string& view_name);

/** @brief Convert a 2d set of vectors (`std::vector<std::vector<Real>>`) to `Kokkos::View<Pack<Real,HAERO_PACK_SIZE>**>`

  view(i,j) = 1 x PACK_SIZE block

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_2d_pack_type vectors_to_row_packed_2dview(const std::vector<std::vector<Real>>& vectors, const std::string& view_name);

/** @brief Convert a 2d set of vectors (`std::vector<std::vector<Real>>`) to `Kokkos::View<Pack<Real,HAERO_PACK_SIZE>**>`

  view(i,j) = PACK_SIZE x 1 block

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_2d_pack_type vectors_to_col_packed_2dview(const std::vector<std::vector<Real>>& vectors, const std::string& view_name);


//---------------------- Impl details (don't call these directly) ---------------//

template <typename VT>
typename std::enable_if<!ekat::ScalarTraits<typename VT::value_type>::is_simd, std::vector<Real>>::type
view1d_to_vector_impl(const VT& v) {
  auto hm = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hm, v);
  std::vector<Real> result(v.extent(0));
  for (int i=0; i<v.extent(0); ++i) {
    result[i] = hm(i);
  }
  return result;
}

template <typename VT>
typename std::enable_if<ekat::ScalarTraits<typename VT::value_type>::is_simd, std::vector<Real>>::type
view1d_to_vector_impl(const VT& v) {
  auto hm = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hm, v);
  std::vector<Real> result(v.extent(0));
  for (int i=0; i<v.extent(0); ++i) {
    result[i] = hm(pack_info::pack_idx(i))[pack_info::vec_idx(i)];
  }
  return result;
}

} // namespace haero
#endif
