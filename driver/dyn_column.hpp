#ifndef HAERO_DYN_COLUMN_HPP
#define HAERO_DYN_COLUMN_HPP

#include "haero/haero_config.hpp"
#include "column_base.hpp"
#include <map>
#include <string>
#include <vector>

namespace haero {

class dyn_column : public column_base {
  public:

  dyn_column() = delete;

  dyn_column(const int& nl);

  std::map<std::string, view_1d> level_vars;
  mask_view level_masks;
  std::map<std::string, view_1d> interface_vars;
  mask_view interface_masks;

  int npack_mid;
  int npack_interface;
  int nlev;

  std::string info_string(const int& tab_lev=0) const ;
};


} // namespace haero
#endif
