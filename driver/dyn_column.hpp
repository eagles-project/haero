#ifndef HAERO_DYN_COLUMN_HPP
#define HAERO_DYN_COLUMN_HPP

#include "haero/haero_config.hpp"
#include "column_base.hpp"
#include "atmosphere.hpp"
#include <map>
#include <string>
#include <vector>

namespace haero {

class DynColumn : public ColumnBase {
  public:
  DynColumn(const int& nl) : ColumnBase(nl),
    w("w", num_packs_int()),
    phi("phi", num_packs_int()),
    dpids("dpids", num_packs_lev()),
    thetav("thetav", num_packs_lev()),
    p("p", num_packs_lev()),
    qv("qv", num_packs_lev()),
    pi("pi", num_packs_lev()),
    mu("mu", num_packs_int()),
    exner("exner", num_packs_int()) {}

  void init_from_interface_heights(const AtmosphericConditions& conds, const std::vector<Real>& z_vals);

  void init_from_interface_pressures(const AtmosphericConditions& conds, const std::vector<Real>& p_vals);

  std::string info_string(const int& tab_lev=0) const ;

  view_1d w;
  view_1d phi;
  view_1d dpids;
  view_1d thetav;
  view_1d p;
  view_1d qv;
  view_1d pi;
  view_1d mu;
  view_1d exner;

  DynColumn() = delete;

  protected:
};


} // namespace haero
#endif
