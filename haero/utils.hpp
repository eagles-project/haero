#ifndef HAERO_UTILS_HPP
#define HAERO_UTILS_HPP

#include <iostream>
#include <string>
#include <vector>

#include "haero/haero.hpp"

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

/// Write a horizontal line to a string
std::string line_delim();

/** @brief Returns true if a std::vector's elements are either increasing or
  decreasing.

  Returns false if the vector contains duplicate values or is not in an
  increasing or decreasing order.

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
}  // namespace haero
#endif
