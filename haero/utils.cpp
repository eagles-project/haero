// Copyright (c) 2021, National Technology & Engineering Solutions of Sandia,
// LLC (NTESS). Copyright (c) 2022, Battelle Memorial Institute
// SPDX-License-Identifier: BSD-3-Clause

#include "utils.hpp"

#include <ekat/ekat_assert.hpp>

#include <cstring>
#include <sstream>

namespace haero {

std::string indent_string(const int tab_lev) {
  std::string result("");
  for (int i = 0; i < tab_lev; ++i) {
    result += "\t";
  }
  return result;
}

std::string line_delim() { return "------------------------------\n"; }

std::string &tolower(std::string &s) {
  for (auto &c : s) {
    c = std::tolower(c);
  }
  return s;
}

bool is_boolean(const std::string &s) {
  return ((strcasecmp("true", s.c_str()) == 0) ||
          (strcasecmp("false", s.c_str()) == 0) ||
          (strcasecmp("yes", s.c_str()) == 0) ||
          (strcasecmp("no", s.c_str()) == 0) ||
          (strcasecmp("on", s.c_str()) == 0) ||
          (strcasecmp("off", s.c_str()) == 0));
}

bool as_boolean(const std::string &s) {
  EKAT_REQUIRE_MSG(is_boolean(s), "'" << s << "' is not a boolean!\n"
                                      << "Must be one of (case-insensitive):\n"
                                      << "* \"true\", \"false\"\n"
                                      << "* \"yes\", \"no\"\n"
                                      << "* \"on\", \"off\",\n");
  return ((strcasecmp("true", s.c_str()) == 0) ||
          (strcasecmp("yes", s.c_str()) == 0) ||
          (strcasecmp("on", s.c_str()) == 0));
}

std::string get_filename_ext(const std::string &fname) {
  const auto dotpos = fname.rfind('.');
  std::string result = "";
  if (dotpos != std::string::npos) {
    result = fname.substr(dotpos);
  }
  return result;
}

bool vector_is_monotone(const std::vector<Real> &vals) {
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

} // namespace haero
