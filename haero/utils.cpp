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

std::string get_filename_ext(const std::string& fname) {
  const auto dotpos = fname.rfind('.');
  std::string result = "";
  if (dotpos != std::string::npos) {
    result = fname.substr(dotpos);
  }
  return result;
}

}// namespace haero
