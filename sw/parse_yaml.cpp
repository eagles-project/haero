#include "parse_yaml.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cstdarg>
#include <set>

namespace {

// Root node of the YAML file being parsed
YAML::Node root_;

// Here are C strings that have been constructed from Fortran strings.
std::set<std::string>* fortran_strings_ = nullptr;

void delete_fortran_strings() {
  // Delete all relevant Fortran string instances.
  if (fortran_strings_ != nullptr) {
    delete fortran_strings_;
    fortran_strings_ = nullptr;
  }
}

} // namespace

extern "C" {

// Returns a newly-allocated C string for the given Fortran string pointer with
// the given length. Resources for this string are managed by the running model.
const char* new_c_string(char* f_str_ptr, int f_str_len) {
  if (fortran_strings_ == nullptr) {
    fortran_strings_ = new std::set<std::string>();
    atexit(delete_fortran_strings);
  }

  // Before we allocate any more strings, check to see whether we've already
  // got this one.
  std::string f_str(f_str_ptr, f_str_len);
  auto iter = fortran_strings_->find(f_str);
  if (iter != fortran_strings_->end()) { // got it!
    return iter->c_str();
  } else { // we need a new one
    auto inserted = fortran_strings_->insert(f_str);
    return inserted.first->c_str();
  }
}

// Begin parsing the given YAML file. Halts on failure with an error message.
void parse_yaml_begin(const char* filename) {
  if (root_.IsNull()) {
    // Load the contents of the file into our root node.
    try {
      root_ = YAML::LoadFile(filename);
    } catch (std::exception& e) {
      fprintf(stderr, "Error encountered reading %s: %s\n", filename, e.what());
      exit(1);
    }
  } else {
    fprintf(stderr, "Called parse_yaml_begin() more than once!\n");
    exit(1);
  }
}

// Get the specified process name from the YAML file. Halts on failure with
// an error message.
void parse_yaml_process(const char** process_name) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse process name without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["process"] and root_["process"].IsMap()) {
    auto process = root_["process"];
    if (process["mam"]) {
      auto proc_name = process["mam"].as<std::string>();
      *process_name = new_c_string(const_cast<char*>(proc_name.c_str()), proc_name.length());
    } else {
      fprintf(stderr, "'mam' entry not found in process section!\n");
      exit(1);
    }
  } else {
    fprintf(stderr, "Did not find a valid process section!\n");
    exit(1);
  }
}

// Get the timestepping information from the YAML file. Halts on failure with
// an error message.
void parse_yaml_timestepping(double* dt, int* nsteps) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse timestepping data without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["timestepping"] and root_["timestepping"].IsMap()) {
    auto ts = root_["timestepping"];
    if (ts["dt"]) {
      *dt = ts["dt"].as<double>();
    } else {
      fprintf(stderr, "'dt' not found in timestepping section!\n");
      exit(1);
    }
    if (ts["nsteps"]) {
      *nsteps = ts["nsteps"].as<int>();
    } else {
      fprintf(stderr, "'nsteps' not found in timestepping section!\n");
      exit(1);
    }
  } else {
    fprintf(stderr, "Did not find a valid timestepping section!\n");
    exit(1);
  }
}

// Get the number of parameters to be walked. Halts on failure with an error
// message.
void parse_yaml_num_parameters(int* num_params) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse parameter data without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["parameters"] and root_["parameters"].IsMap()) {
    *num_params = static_cast<int>(root_["parameters"].size());
  } else {
    fprintf(stderr, "Did not find a valid parameters section!\n");
    exit(1);
  }
}

namespace {

int get_num_param_values(const std::string& param_name, const YAML::Node& param) {
  // If the parameter has length 3, it can either be an array of values or a
  // [start, stop, step] specification. In the former case, the values should
  // appear in ascending order.
  int len = static_cast<int>(param.size());
  if (len == 3) {
    double value0 = param[0].as<double>();
    double value1 = param[1].as<double>();
    double value2 = param[2].as<double>();
    if (not (value0 < value1)) { // this should be true in any case!
      fprintf(stderr, "Invalid values for %s: second must be greater than first.",
              param_name.c_str());
      exit(1);
    }
    if (not (value1 < value2)) { // [start, stop, step]
      return static_cast<int>((value1 - value0)/value2) + 1;
    }
  }
  return len;
}

}

// Get the name of the (1-based) nth parameter and its number of values. Halts
// on failure with an error message.
void parse_yaml_parameter(int n, const char** parameter_name, int* num_values) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse parameter data without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["parameters"] and root_["parameters"].IsMap()) {
    auto params = root_["parameters"];
    int num_params = static_cast<int>(params.size());
    if ((n < 1) || (n > num_params)) {
      fprintf(stderr, "Invalid parameter index: %d (must be between 1 and %d)\n",
        n, num_params);
      exit(1);
    } else {
      int i = 1;
      for (auto iter = params.begin(); iter != params.end(); ++iter, ++i) {
        if (i == n) {
          auto param = iter->second;
          if (param.IsSequence()) {
            auto param_name = iter->first.as<std::string>();
            *parameter_name = new_c_string(const_cast<char*>(param_name.c_str()), param_name.length());
            *num_values = get_num_param_values(param_name, param);
          } else {
            fprintf(stderr, "Parameter '%s' (%d) is not a sequence of values!\n",
              iter->first.as<std::string>().c_str(), n);
            exit(1);
          }
        }
      }
    }
  } else {
    fprintf(stderr, "Did not find a valid parameters section!\n");
    exit(1);
  }
}

// Get the values for the nth parameter.
void parse_yaml_parameter_values(int n, double* param_values) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse parameter data without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["parameters"] and root_["parameters"].IsMap()) {
    auto params = root_["parameters"];
    int num_params = static_cast<int>(params.size());
    if ((n < 1) || (n > num_params)) {
      fprintf(stderr, "Invalid parameter index: %d (must be between 1 and %d)\n",
        n, num_params);
      exit(1);
    } else {
      int i = 1;
      for (auto iter = params.begin(); iter != params.end(); ++iter, ++i) {
        if (i == n) {
          auto param = iter->second;
          if (param.IsSequence()) {
            auto param_name = iter->first.as<std::string>();
            int num_values = get_num_param_values(param_name, param);
            if (num_values != static_cast<int>(param.size())) { // [start, stop, step]
              double low = param[0].as<double>();
              double step_size = param[2].as<double>();
              for (int j = 0; j < num_values; ++j) {
                param_values[j] = low + j*step_size;
              }
            } else { // array
              for (int j = 0; j < num_values; ++j) {
                param_values[j] = param[j].as<double>();
              }
            }
          } else {
            fprintf(stderr, "Parameter '%s' (%d) is not a sequence of values!\n",
              iter->first.as<std::string>().c_str(), n);
            exit(1);
          }
        }
      }
    }
  } else {
    fprintf(stderr, "Did not find a valid parameters section!\n");
    exit(1);
  }
}

// Get the initial conditions for the atmosphere.
void parse_yaml_atmosphere(double* temp,
                           double* press,
                           double* rel_hum,
                           double* height,
                           double* cld_frac) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse atmosphere data without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["atmosphere"] and root_["atmosphere"].IsMap()) {
    auto atm = root_["atmosphere"];
    if (atm["temperature"]) {
      *temp = atm["temperature"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'temperature' in the atmosphere section!\n");
      exit(1);
    }
    if (atm["pressure"]) {
      *press = atm["pressure"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'pressure' in the atmosphere section!\n");
      exit(1);
    }
    if (atm["relative_humidity"]) {
      *rel_hum = atm["relative_humidity"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'relative_humidity' in the atmosphere section!\n");
      exit(1);
    }
    if (atm["height"]) {
      *height = atm["height"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'height' in the atmosphere section!\n");
      exit(1);
    }
    if (atm["cloud_fraction"]) {
      *cld_frac = atm["cloud_fraction"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'cloud_fraction' in the atmosphere section!\n");
      exit(1);
    }
  } else {
    fprintf(stderr, "Did not find a valid atmosphere section!\n");
    exit(1);
  }
}

// Get the initial conditions for (MAM4) aerosol prognostics.
void parse_yaml_aerosols(// modal number concentrations
                         double* nc1, double* nc2, double* nc3, double* nc4,
                         // mode 1 mixing fractions
                         double* mfso41, double* mfpom1, double* mfsoa1,
                         double* mfbc1, double* mfdst1, double* mfncl1,
                         double* mfmom1,
                         // mode 2 mixing fractions
                         double* mfso42, double* mfsoa2, double* mfncl2,
                         double* mfmom2,
                         // mode 3 mixing fractions
                         double* mfdst3, double* mfncl3, double* mfso43,
                         double* mfbc3, double* mfpom3, double* mfsoa3,
                         // mode 4 mixing fractions
                         double* mfpom4, double* mfbc4, double* mfmom4) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse aerosol data without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["aerosols"] and root_["aerosols"].IsMap()) {
    auto aerosols = root_["aerosols"];
    if (aerosols["accumulation"]) {
      auto accum = aerosols["accumulation"];
      if (accum["number_conc"]) {
        *nc1 = accum["number_conc"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'number_conc' in accumulation mode!\n");
        exit(1);
      }
      if (accum["so4"]) {
        *mfso41 = accum["so4"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'so4' aerosol in accumulation mode!\n");
        exit(1);
      }
      if (accum["pom"]) {
        *mfpom1 = accum["pom"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'pom' aerosol in accumulation mode!\n");
        exit(1);
      }
      if (accum["soa"]) {
        *mfsoa1 = accum["soa"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'soa' aerosol in accumulation mode!\n");
        exit(1);
      }
      if (accum["bc"]) {
        *mfbc1 = accum["bc"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'bc' aerosol in accumulation mode!\n");
        exit(1);
      }
      if (accum["dst"]) {
        *mfdst1 = accum["dst"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'dst' aerosol in accumulation mode!\n");
        exit(1);
      }
      if (accum["ncl"]) {
        *mfncl1 = accum["ncl"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'ncl' aerosol in accumulation mode!\n");
        exit(1);
      }
      if (accum["mom"]) {
        *mfmom1 = accum["mom"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'mom' aerosol in accumulation mode!\n");
        exit(1);
      }
    } else {
      fprintf(stderr, "Did not find 'accumulation' mode in aerosols section!\n");
      exit(1);
    }
    if (aerosols["aitken"]) {
      auto aitken = aerosols["aitken"];
      if (aitken["number_conc"]) {
        *nc2 = aitken["number_conc"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'number_conc' in aitken mode!\n");
        exit(1);
      }
      if (aitken["so4"]) {
        *mfso42 = aitken["so4"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'so4' aerosol in aitken mode!\n");
        exit(1);
      }
      if (aitken["soa"]) {
        *mfsoa2 = aitken["soa"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'soa' aerosol in aitken mode!\n");
        exit(1);
      }
      if (aitken["ncl"]) {
        *mfncl2 = aitken["ncl"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'ncl' aerosol in aitken mode!\n");
        exit(1);
      }
      if (aitken["mom"]) {
        *mfmom2 = aitken["mom"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'mom' aerosol in aitken mode!\n");
        exit(1);
      }
    } else {
      fprintf(stderr, "Did not find 'aitken' mode in aerosols section!\n");
      exit(1);
    }
    if (aerosols["coarse"]) {
      auto coarse = aerosols["coarse"];
      if (coarse["number_conc"]) {
        *nc3 = coarse["number_conc"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'number_conc' in coarse mode!\n");
        exit(1);
      }
      if (coarse["so4"]) {
        *mfso43 = coarse["so4"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'so4' aerosol in coarse mode!\n");
        exit(1);
      }
      if (coarse["pom"]) {
        *mfpom3 = coarse["pom"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'pom' aerosol in coarse mode!\n");
        exit(1);
      }
      if (coarse["soa"]) {
        *mfsoa3 = coarse["soa"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'soa' aerosol in coarse mode!\n");
        exit(1);
      }
      if (coarse["bc"]) {
        *mfbc3 = coarse["bc"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'bc' aerosol in coarse mode!\n");
        exit(1);
      }
      if (coarse["dst"]) {
        *mfdst3 = coarse["dst"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'dst' aerosol in coarse mode!\n");
        exit(1);
      }
      if (coarse["ncl"]) {
        *mfncl3 = coarse["ncl"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'ncl' aerosol in coarse mode!\n");
        exit(1);
      }
    } else {
      fprintf(stderr, "Did not find 'coarse' mode in aerosols section!\n");
      exit(1);
    }
    if (aerosols["primary_carbon"]) {
      auto pcarbon = aerosols["primary_carbon"];
      if (pcarbon["number_conc"]) {
        *nc4 = pcarbon["number_conc"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'number_conc' in primary_carbon mode!\n");
        exit(1);
      }
      if (pcarbon["pom"]) {
        *mfpom4 = pcarbon["pom"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'pom' aerosol in primary_carbon mode!\n");
        exit(1);
      }
      if (pcarbon["bc"]) {
        *mfbc4 = pcarbon["bc"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'bc' aerosol in primary_carbon mode!\n");
        exit(1);
      }
      if (pcarbon["mom"]) {
        *mfmom4 = pcarbon["mom"].as<double>();
      } else {
        fprintf(stderr, "Did not find 'mom' aerosol in primary_carbon mode!\n");
        exit(1);
      }
    } else {
      fprintf(stderr, "Did not find 'primary_carbon' mode in aerosols section!\n");
      exit(1);
    }
  } else {
    fprintf(stderr, "Did not find a valid aerosols section!\n");
    exit(1);
  }
}

// Get the initial conditions for (MAM4) gases.
void parse_yaml_gases(double* qso2, double* qh2so4, double* qsoag) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to parse gas data without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  }
  else if (root_["gases"] and root_["gases"].IsMap()) {
    auto gases = root_["gases"];
    if (gases["qso2"]) {
      *qso2 = gases["qso2"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'qso2' in the gases section!\n");
      exit(1);
    }
    if (gases["qh2so4"]) {
      *qh2so4 = gases["qh2so4"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'qh2so4' in the gases section!\n");
      exit(1);
    }
    if (gases["qsoag"]) {
      *qsoag = gases["qsoag"].as<double>();
    } else {
      fprintf(stderr, "Did not find 'qsoag' in the gases section!\n");
      exit(1);
    }
  } else {
    fprintf(stderr, "Did not find a valid gases section!\n");
    exit(1);
  }
}

// Finish parsing the YAML file.
void parse_yaml_end(void) {
  if (root_.IsNull()) {
    fprintf(stderr, "Attempted to call parse_yaml_end() without calling "
      "parse_yaml_begin()!\n");
    exit(1);
  } else {
    root_.reset();
  }
}

} // extern "C"

ParameterWalk parse_yaml(const std::string& filename) {
  ParameterWalk pw;
  return pw;
}
