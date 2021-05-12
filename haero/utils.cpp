#include "utils.hpp"

#include <sstream>

namespace haero {

std::string indent_string(const int tab_lev) {
  std::string result("");
  for (int i = 0; i < tab_lev; ++i) {
    result += "\t";
  }
  return result;
}

std::string& tolower(std::string& s) {
  for (auto& c : s) {
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

bool vector_is_monotone(const std::vector<Real>& vals) {
  bool result = true;
  if (vals.size() > 1) {
    // check increasing or decreasing based on first 2 values
    bool increasing = (vals[1] > vals[0]);
    if (increasing) {
      for (int i = 1; i < vals.size(); ++i) {
        if (vals[i] <= vals[i - 1]) {
          result = false;
          break;
        }
      }
    } else {
      for (int i = 1; i < vals.size(); ++i) {
        if (vals[i] >= vals[i - 1]) {
          result = false;
          break;
        }
      }
    }
  }
  return result;
}

}  // namespace haero
