#include "utils.hpp"
#include <sstream>

namespace haero {

std::string indent_string(const int tab_lev) {
  std::string result("");
  for (int i=0; i<tab_lev; ++i) {
    result += "\t";
  }
  return result;
}

std::string& tolower(std::string& s) {
  for (auto& c: s) {
    c = std::tolower(c);
  }
  return s;
}

}// namespace haero
