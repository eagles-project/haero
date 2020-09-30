#ifndef HAERO_COLUMN_BASE_HPP
#define HAERO_COLUMN_BASE_HPP

#include "haero/haero_config.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "kokkos/Kokkos_Core.hpp"

namespace haero {

/** @brief Base class for the driver's columns.

  Provides common types and methods for using Views of Scalar Packs.

  Uses EKAT to define common kokkos types.

  Defines default view types for Haero's Driver
*/
class ColumnBase {
  public :

  /// Device Type
  using kokkos_device_types = ekat::KokkosTypes<ekat::DefaultDevice>;
  /// Host Type
  using kokkos_host_types = ekat::KokkosTypes<ekat::HostDevice>;
  /// Pack Type for simd ops on cpu
  using real_pack_type = ekat::Pack<Real,HAERO_PACK_SIZE>;
  /// Mask type for use with Packs
  using mask_type = ekat::Mask<HAERO_PACK_SIZE>;

  /// A 1d view of packs for 1d variables
  using view_1d = kokkos_device_types::view_1d<real_pack_type>;
  /// A 2d view of packs for 2d variables
  using view_2d = kokkos_device_types::view_2d<real_pack_type>;
  /** A view type that adds a dimension for number of modes.

    @todo Should we keep this?
  */
  using modal_var_view = kokkos_device_types::view_2d<Real>;
  /// A view of masks. Currently unused.
  using mask_view = kokkos_device_types::view_1d<mask_type>;

  /** @brief This is perhaps the most useful typedef in ColumnBase.

    To loop over packed views, you have to use 2 nested loops, one to index
      over packs, and one to index within packs.

      ````
      for (int pack_idx=0; pack_idx<num_packs; ++pack_idx) {
        for (int vec_idx=0; vec_idx<PACK_SIZE; ++vec_idx) {
          //********** Convert from (pack, vec) pairs to level or interface index
          const int array_index = pack_info::array_idx(pack_idx, vec_idx);
        }
      }
      ````

      But you have to be careful --- the loops as written above will also loop over
      the padded, unused indices in a Pack.  To fix this, instead of stopping the
      inner loop at `vec_idx = PACK_SIZE`, use `pack_info::vec_end(array_length, pack_idx)`
      instead, where `array_idx` = number of levels (for midpoint variables) or
      `array_idx` = number of levels + 1 (for interfaces).

    To loop over levels or interfaces and get (pack, vec) pairs, use 1 loop:

    ````
      for (int lev_ind = 0; lev_ind < num_levels; ++lev_ind) {
        //********** Convert from level or interface index to (pack, vec) pairs
        const int pack_idx = pack_info::pack_idx(lev_ind);
        cosnt int vec_idx  = pack_info::vec_idx(lev_ind);
      }
    ````

  */
  using pack_info = ekat::PackInfo<HAERO_PACK_SIZE>;

  /** @brief return the number of levels per column
  */
  inline int num_levels() const {return m_num_levels;}

  /** @brief return the number of packs used for level variables.
  */
  inline int num_packs_lev() const {return m_num_packs_lev;}

  /** @brief return the number of packs used for interface variables
  */
  inline int num_packs_int() const {return m_num_packs_int;}

  ColumnBase() = delete;

  protected:
    /** @brief Constructor.

      @param [in] nl number of vertical levels (not interfaces) in the column
    */
    ColumnBase(const int& nl) : m_num_levels(nl), m_num_packs_lev(pack_info::num_packs(nl)),
      m_num_packs_int(pack_info::num_packs(nl+1)) {}

  int m_num_levels;
  int m_num_packs_lev;
  int m_num_packs_int;
};

} // namespace haero
#endif
