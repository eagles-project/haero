#ifndef HAERO_UTILS_HPP
#define HAERO_UTILS_HPP

#include "haero/haero_config.hpp"
#include <string>
#include <iostream>

namespace haero {

std::string indent_string(const int tab_lev);

std::string& tolower(std::string& s);

class ProgressBar {
  std::string name_;
  int niter_;
  Real freq_;
  int it_;
  Real next_;
  std::ostream& os_;

  public:
    ProgressBar(const std::string& name, const int niterations,
      const Real write_freq = 10.0, std::ostream& os = std::cout);

    void update();
};

}
#endif
