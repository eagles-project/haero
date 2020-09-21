#include "column_base.hpp"
#include "dyn_column.hpp"
#include "haero/utils.hpp"
#include <sstream>

namespace haero {

dyn_column::dyn_column(const int& nl) : nlev(nl) {
  if (nlev%HAERO_PACK_SIZE == 0) {
    npack_mid = nlev/HAERO_PACK_SIZE;
  }
  else {
    npack_mid = nlev/HAERO_PACK_SIZE+1;
  }
  if ((nlev+1)%HAERO_PACK_SIZE == 0) {
    npack_interface = (nlev+1)/HAERO_PACK_SIZE;
  }
  else {
    npack_interface = (nlev+2)/HAERO_PACK_SIZE;
  }
  level_masks = mask_view("level_masks", npack_mid);
  interface_masks = mask_view("interface_masks", npack_interface);

  auto lmasks = this->level_masks;
  Kokkos::parallel_for("init_mask_vals_levels",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0}, {npack_mid,HAERO_PACK_SIZE}),
    KOKKOS_LAMBDA (const int& i, const int& j) {
      const int lev_ind = pack_to_level_ind(i,j);
      lmasks(i).set(j, !(lev_ind < nlev));
    });
  auto imasks = this->interface_masks;
  Kokkos::parallel_for("init_mask_vals_interfaces",
    Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{npack_interface, HAERO_PACK_SIZE}),
    KOKKOS_LAMBDA (const int& i, const int& j) {
      const int ifc_ind = pack_to_level_ind(i,j);
      imasks(i).set(j, !(ifc_ind < nlev+1));
    });

  level_vars.emplace("theta_v", view_1d("theta_v",npack_mid));
  level_vars.emplace("dpids", view_1d("dpids", npack_mid));
  level_vars.emplace("exner", view_1d("exner", npack_mid));
  level_vars.emplace("qv", view_1d("qv", npack_mid));

  interface_vars.emplace("w", view_1d("w", npack_interface));
  interface_vars.emplace("phi", view_1d("phi", npack_interface));
  interface_vars.emplace("pi", view_1d("pi", npack_interface));
  interface_vars.emplace("dpds", view_1d("dpds",npack_interface));
  interface_vars.emplace("mu", view_1d("mu", npack_interface));
}

std::string dyn_column::info_string(const int& tab_lev) const {
  std::ostringstream ss;
  std::string tabstr = indent_string(tab_lev);
  ss << tabstr << "haero dynamics column info:\n";
  tabstr += "\t";
  ss << tabstr << "nlev = " << nlev << '\n';
  ss << tabstr << "npack_mid = " << npack_mid << '\n';
  ss << tabstr << "npack_interface = " << npack_interface << '\n';
  ss << tabstr << "allocated level midpoint variables:\n";
  for (auto& var : level_vars) {
    ss << tabstr << '\t' << var.first << '\n';
  }
  ss << tabstr << "allocated level interface variables:\n";
  for (auto& var : interface_vars) {
    ss << tabstr << '\t' << var.first << '\n';
  }
  return ss.str();
}

} // namespace haero
