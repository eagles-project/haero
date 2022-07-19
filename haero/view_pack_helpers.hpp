#ifndef HAERO_VIEW_HELPERS_HPP
#define HAERO_VIEW_HELPERS_HPP

#include <Kokkos_Core.hpp>
#include <ekat/ekat_scalar_traits.hpp>
#include <type_traits>
#include <vector>

#include "haero/haero.hpp"
#include "haero/haero_config.hpp"

namespace haero {

//---------------------- Compile-time stuff ---------------//

/// Kokkos device definitions
using kokkos_device_type = DeviceType;
using kokkos_host_type = HostType;

/// Ekat Pack definitions
using real_pack_type = ekat::Pack<Real, HAERO_PACK_SIZE>;
using pack_mask_type = ekat::Mask<HAERO_PACK_SIZE>;
using pack_info = ekat::PackInfo<HAERO_PACK_SIZE>;

/// View definitions
using view_1d_scalar_type = kokkos_device_type::view_1d<Real>;
using view_1d_int_type = kokkos_device_type::view_1d<int>;
using view_1d_pack_type = kokkos_device_type::view_1d<real_pack_type>;
using mask_view_1d_type = kokkos_device_type::view_1d<pack_mask_type>;

using view_2d_scalar_type = kokkos_device_type::view_2d<Real>;
using view_2d_pack_type = kokkos_device_type::view_2d<real_pack_type>;
using mask_view_2d_type = kokkos_device_type::view_2d<pack_mask_type>;

//------------------------ Run-time stuff ---------------//

/// fwd decl
template <typename VT>
typename std::enable_if<!ekat::ScalarTraits<typename VT::value_type>::is_simd,
                        std::vector<Real>>::type
view1d_to_vector_impl(const VT& v, const int& array_length);

/// fwd decl
template <typename VT>
typename std::enable_if<ekat::ScalarTraits<typename VT::value_type>::is_simd,
                        std::vector<Real>>::type
view1d_to_vector_impl(const VT& v, const int& array_length);

/** @brief Convert a std::vector<Real> to Kokkos::View<Real*>.

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_1d_scalar_type vector_to_basic_1dview(const std::vector<Real>& vector,
                                           const std::string& view_name);

/** @brief Convert a std::vector<int> to Kokkos::View<int*>.

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_1d_int_type vector_to_basic_1dview(const std::vector<int>& vector,
                                        const std::string& view_name);

template <typename T>
kokkos_device_type::view_1d<T> vector_to_1dview(const std::vector<T>& vector,
                                                const std::string& view_name) {
  kokkos_device_type::view_1d<T> result(view_name, vector.size());
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < vector.size(); ++i) {
    hm(i) = vector[i];
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

/** @brief Convert a std::vector<Real> to
  Kokkos::View<Pack<Real,HAERO_PACK_SIZE>*>

  @param [in] vector
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_1d_pack_type vector_to_packed_1dview(const std::vector<Real>& vector,
                                          const std::string& view_name);

/** @brief Converts a 1d view (or 1d subview) to a std::vector

  @warning For packed data, the length of the output vector may not be the same
  as the extent of the view.

  @param [in] v view with ViewType::Rank = 1
  @param [in] array_length length of output vector
*/
template <typename ViewType>
std::vector<Real> view1d_to_vector(const ViewType& v,
                                   const int array_length = 0) {
  static_assert(ViewType::Rank == 1,
                "view1d_to_vector error : View must be rank 1.");
  return view1d_to_vector_impl<ViewType>(v, array_length);
}

/** @brief Convert a 2d set of vectors, (`std::vector<std::vector<Real>>`) to
  `Kokkos::View<Real**>`

  @param [in] vectors
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_2d_scalar_type vector_to_basic_2dview(
    const std::vector<std::vector<Real>>& vectors,
    const std::string& view_name);

/** @brief Convert a 2d set of vectors (`std::vector<std::vector<Real>>`) to
  `Kokkos::View<Pack<Real,HAERO_PACK_SIZE>**>`

  view(i,j) = 1 x PACK_SIZE block

  @param [in] vectors
  @param [in] view_name
  @return Kokkos view + deep copied data
*/
view_2d_pack_type vectors_to_row_packed_2dview(
    const std::vector<std::vector<Real>>& vectors,
    const std::string& view_name);

template <typename T>
kokkos_device_type::view_2d<T> vector_to_2dview(
    const std::vector<std::vector<T>>& vectors, const std::string& view_name) {
  const int mm = vectors.size();
  const int nn = vectors[0].size();
  kokkos_device_type::view_2d<T> result(view_name, mm, nn);
  auto hm = Kokkos::create_mirror_view(result);
  for (int i = 0; i < mm; ++i) {
    for (int j = 0; j < nn; ++j) {
      hm(i, j) = vectors[i][j];
    }
  }
  Kokkos::deep_copy(result, hm);
  return result;
}

//---------------------- Impl details (don't call these directly)
//---------------//

template <typename VT>
typename std::enable_if<!ekat::ScalarTraits<typename VT::value_type>::is_simd,
                        std::vector<Real>>::type
view1d_to_vector_impl(const VT& v, const int& array_length) {
  auto hm = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hm, v);
  std::vector<Real> result(v.extent(0));
  for (int i = 0; i < v.extent(0); ++i) {
    result[i] = hm(i);
  }
  return result;
}

template <typename VT>
typename std::enable_if<ekat::ScalarTraits<typename VT::value_type>::is_simd,
                        std::vector<Real>>::type
view1d_to_vector_impl(const VT& v, const int& array_length) {
  auto hm = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(hm, v);
  std::vector<Real> result(array_length);
  for (int i = 0; i < v.extent(0); ++i) {
    for (int j = 0; j < pack_info::vec_end(array_length, i); ++j) {
      result[pack_info::array_idx(i, j)] = hm(i)[j];
    }
  }
  return result;
}

}  // namespace haero
#endif
