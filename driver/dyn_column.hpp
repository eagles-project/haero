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
    w("w", npack_interface),
    phi("phi", npack_interface),
    dpids("dpids", npack_mid),
    thetav("thetav", npack_mid),
    p("p", npack_mid),
    qv("qv", npack_mid),
    pi("pi", npack_interface),
    mu("mu", npack_interface),
    exner("exner", npack_mid) {}

  void initialize(const AtmosphericConditions& conds);

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
