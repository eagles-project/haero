#ifndef HAERO_PHYS_COLUMN_HPP
#define HAERO_PHYS_COLUMN_HPP

#include "haero/haero_config.hpp"
#include "column_base.hpp"
#include <map>
#include <string>
#include <vector>

namespace haero {
namespace driver {

class phys_column : public column_base {
  public:

  phys_column() = delete;

  phys_column(const int& nlevs);

  phys_column(const int& nlevs, const std::vector<std::string>& var_names);

  std::map<std::string, view_1d> level_vars;
  std::map<std::string, view_1d> interface_vars;

  int npack_mid;
  int npack_interface;
  int nlev;
};


} // namespace driver
} // namespace haero
#endif
