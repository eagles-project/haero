#ifndef HAERO_UTILS_HPP
#define HAERO_UTILS_HPP

#include "haero/haero_config.hpp"
#include <string>
#include <iostream>
#include <vector>

namespace haero {

/** @defgroup Utilities Miscellaneous utilities

 @{
*/

/// return a string with tab_lev copies of the indent character, '\\t'
std::string indent_string(const int tab_lev);

/// convert a string to lower case
std::string& tolower(std::string& s);

/// Get the filename extension from a string containing a filename
std::string get_filename_ext(const std::string& fname);

/** @brief Returns true if a std::vector's elements are either increasing or decreasing.

  Returns false if the vector contains duplicate values or is not in an increasing or
    decreasing order.

  @param [in] vals
*/
bool vector_is_monotone(const std::vector<Real>& vals);

/** @brief A simple progress bar to show program % completion with the console.
*/
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


/// @} defgroup utilities
} // namespace haero
#endif
