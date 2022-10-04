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

/// View definitions
using view_1d_scalar_type = DeviceType::view_1d<Real>;
using view_1d_int_type = DeviceType::view_1d<int>;
using view_1d_pack_type = DeviceType::view_1d<PackType>;
using mask_view_1d_type = DeviceType::view_1d<MaskType>;

using view_2d_scalar_type = DeviceType::view_2d<Real>;
using view_2d_pack_type = DeviceType::view_2d<PackType>;
using mask_view_2d_type = DeviceType::view_2d<MaskType>;

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

/** @brief Initialize a view of packs to zero, preserve quiet nan entries
  in padding.

  The default Pack constructor initializes all entries to scalar_type::invalid
  (usually a quiet nan) to protect users from mishandling padded entries.

  This function initializes non-padded entries to zero, and preserves the
  invalid values for the padded entries.
*/
template <typename T>
void zero_init(DeviceType::view_1d<ekat::Pack<T,HAERO_PACK_SIZE>> view, const int n_entries) {
  const int np = PackInfo::num_packs(view.extent(0));
  if (n_entries == np) {
    // pack size = 1
    ekat::Pack<T,HAERO_PACK_SIZE> zero_pack(T(0));
    Kokkos::deep_copy(view, zero_pack);
  }
  else {
    auto h_view = Kokkos::create_mirror_view(view);
    const int last_pack_idx = PackInfo::last_pack_idx(n_entries);
    for (int i=0; i<last_pack_idx; ++i) {
      h_view(i) = ekat::Pack<T,HAERO_PACK_SIZE>(T(0));
    }
    ekat::Mask<HAERO_PACK_SIZE> last_mask(false);
    for (int i=0; i<PackInfo::last_vec_end(n_entries); ++i) {
      last_mask.set(i, true);
    }
    h_view(last_pack_idx) = ekat::Pack<T,HAERO_PACK_SIZE>(last_mask, T(0));

    Kokkos::deep_copy(view, h_view);
  }
}

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
DeviceType::view_1d<T> vector_to_1dview(const std::vector<T>& vector,
                                                const std::string& view_name) {
  DeviceType::view_1d<T> result(view_name, vector.size());
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
DeviceType::view_2d<T> vector_to_2dview(
    const std::vector<std::vector<T>>& vectors, const std::string& view_name) {
  const int mm = vectors.size();
  const int nn = vectors[0].size();
  DeviceType::view_2d<T> result(view_name, mm, nn);
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
    for (int j = 0; j < PackInfo::vec_end(array_length, i); ++j) {
      result[PackInfo::array_idx(i, j)] = hm(i)[j];
    }
  }
  return result;
}

}  // namespace haero
#endif
