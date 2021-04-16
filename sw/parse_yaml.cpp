#include "parse_yaml.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cstdarg>
#include <set>

ParameterWalk parse_yaml(const haero::ModalAerosolConfig& aerosol_config,
                         const std::string& filename) {
  using Real = haero::Real;

  ParameterWalk pw;
  try {
    auto root = YAML::LoadFile(filename);

    // process section
    if (root["process"] and root["process"].IsMap()) {
      auto process = root["process"];
      if (process["haero"]) {
        pw.process = process["mam"].as<std::string>();
      } else {
        throw YamlException("'haero' entry not found in process section!");
      }
    } else {
      throw YamlException("Did not find a valid process section!");
    }

    // timestepping section
    if (root["timestepping"] and root["timestepping"].IsMap()) {
      auto ts = root["timestepping"];
      if (ts["dt"]) {
        pw.dt = ts["dt"].as<haero::Real>();
      } else {
        throw YamlException("'dt' not found in timestepping section!\n");
      }
      if (ts["total_time"]) {
        pw.total_time = ts["total_time"].as<haero::Real>();
      } else {
        throw YamlException("'total_time' not found in timestepping section!\n");
      }
    } else {
      throw YamlException("Did not find a valid timestepping section!\n");
      exit(1);
    }

    // enѕemble section
    if (root["ensemble"] and root["ensemble"].IsMap()) {
      auto params = root["ensemble"];
      for (auto iter = params.begin(); iter != params.end(); ++iter) {
        auto param = iter->second;
        if (param.IsSequence()) {
          auto param_name = iter->first.as<std::string>();
          size_t len = param.size();
          if (len == 3) {
            Real value0 = param[0].as<Real>();
            Real value1 = param[1].as<Real>();
            Real value2 = param[2].as<Real>();
            if (not (value0 < value1)) { // this should be true in any case!
              throw YamlException(std::string("Invalid values for ") +
                  param_name + std::string(": second must be greater than first."));
            }
            if (not (value1 < value2)) { // [start, stop, step]
              len = static_cast<size_t>((value1 - value0)/value2) + 1;
              pw.ensemble[param_name].resize(len);
              for (int j = 0; j < len; ++j) {
                pw.ensemble[param_name][j] = value0 + j*value2;
              }
            } else {
              pw.ensemble[param_name] = {value0, value1, value2};
            }
          } else {
            pw.ensemble[param_name] = param.as<std::vector<Real>>();
          }
        } else if (param.IsMap()) { // could be an aerosol or gas species mmr
          for (auto iter = param.begin(); iter != param.end(); ++iter) {
            // Is this a valid aerosol mode?
            auto mode_name = iter->first.as<std::string>();
            auto mode = iter->second;
            int mode_index = aerosol_config.aerosol_mode_index(mode_name);
            if (mode_index != -1) {
              auto mode_species = aerosol_config.aerosol_species_for_mode(mode_index);
              for (auto aero_iter = mode.begin(); aero_iter != mode.end(); ++aero_iter) {
                auto aero_name = aero_iter->first.as<std::string>();
                auto aero_species = aero_iter->second;
                // Find the aerosol index within this species.
                int aero_index = aerosol_config.aerosol_species_index(mode_index, aero_name);
                if (aero_index == -1) {
                  throw YamlException(std::string("Found invalid aerosol species '") +
                    aero_name + std::string("' within mode '") + mode_name +
                    std::string("' in the ensemble section!"));
                }
                auto mmr_name = mode_name + std::string(":") + aero_name;
                pw.ensemble[mmr_name] = aero_species.as<std::vector<Real>>();
              }
            } else {
              // Maybe it's a gas species?
              auto gas_name = mode_name;
              auto gas_species = mode;
              int gas_index = aerosol_config.gas_index(gas_name);
              if (gas_index != -1) {
                pw.ensemble[gas_name] = gas_species.as<std::vector<Real>>();
              } else {
                throw YamlException(std::string("Parameter '") +
                  iter->first.as<std::string>() +
                  std::string("' is not a valid aerosol mode or gas species!\n"));
              }
            }
          }
        } else {
          throw YamlException(std::string("Parameter '") +
              iter->first.as<std::string>() +
              std::string("' is not a sequence of values or a map!\n"));
        }
      }
    } else {
      throw YamlException("Did not find a valid ensemble section!\n");
    }

    // atmosphere section
    if (root["atmosphere"] and root["atmosphere"].IsMap()) {
      auto atm = root["atmosphere"];
      if (atm["temperature"]) {
        pw.temperature = atm["temperature"].as<Real>();
      } else {
        throw YamlException("Did not find 'temperature' in the atmosphere section!\n");
      }
      if (atm["pressure"]) {
        pw.pressure = atm["pressure"].as<Real>();
      } else {
        throw YamlException("Did not find 'pressure' in the atmosphere section!\n");
      }
      if (atm["relative_humidity"]) {
        pw.relative_humidity = atm["relative_humidity"].as<Real>();
      } else {
        throw YamlException("Did not find 'relative_humidity' in the atmosphere section!\n");
      }
      if (atm["height"]) {
        pw.height = atm["height"].as<Real>();
      } else {
        throw YamlException("Did not find 'height' in the atmosphere section!\n");
      }
      if (atm["hydrostatic_dp"]) {
        pw.hydrostatic_dp = atm["hydrostatic_dp"].as<Real>();
      } else {
        throw YamlException("Did not find 'hydrostatic_dp' in the atmosphere section!\n");
      }
    } else {
      throw YamlException("Did not find a valid atmosphere section!\n");
    }

    // aerosols section
    pw.number_concs.resize(aerosol_config.h_aerosol_modes.size());
    pw.aero_mmrs.resize(aerosol_config.h_aerosol_modes.size());
    std::vector<int> found_mode(pw.number_concs.size(), 0);
    std::vector<std::vector<int>> found_aerosol(pw.aero_mmrs.size());
    if (root["aerosols"] and root["aerosols"].IsMap()) {
      auto aerosols = root["aerosols"];
      for (auto iter = aerosols.begin(); iter != aerosols.end(); ++iter) {
        auto mode_name = iter->first.as<std::string>();
        auto mode = iter->second;
        if (mode.IsMap()) {
          // Determine the mode's index.
          size_t mode_index;
          for (mode_index = 0; mode_index < aerosol_config.h_aerosol_modes.size(); ++mode_index) {
            if (mode_name == aerosol_config.h_aerosol_modes(mode_index).name()) {
              found_mode[mode_index] = 1;
              break;
            }
          }

          if (mode_index < pw.aero_mmrs.size()) {
            // Get the initial data for the aerosol species in this mode.
            auto mode_species = aerosol_config.aerosol_species_for_mode(mode_index);
            bool found_number_conc = false;
            pw.aero_mmrs[mode_index].resize(mode_species.size());
            found_aerosol[mode_index].resize(mode_species.size());
            for (auto aero_iter = mode.begin(); aero_iter != mode.end(); ++aero_iter) {
              auto aero_name = aero_iter->first.as<std::string>();
              auto aero_species = aero_iter->second;
              if (aero_name == "number_conc") { // number conc, not species name!
                pw.number_concs[mode_index] = aero_species.as<Real>();
                found_number_conc = true;
              } else {
                // Find the aerosol index within this species.
                size_t aero_index;
                for (aero_index = 0; aero_index < mode_species.size(); ++aero_index) {
                  if (aero_name == mode_species[aero_index].symbol()) {
                    found_aerosol[mode_index][aero_index] = 1;
                    break;
                  }
                }
                if (aero_index < mode_species.size()) {
                  pw.aero_mmrs[mode_index][aero_index] = aero_species.as<Real>();
                } else {
                  throw YamlException(std::string("Found invalid aerosol species '") +
                    aero_name + std::string("' within mode '") + mode_name +
                    std::string("' in the aerosols section!"));
                }
              }
              if (not found_number_conc) {
                throw YamlException(std::string("Did not find 'number_conc' in '") +
                  mode_name + std::string("' mode!"));
              }
            }
          } else {
            throw YamlException(std::string("Found invalid aerosol mode '") +
              mode_name + std::string("' in the aerosols section!"));
          }
        }
      }

      // Did we find all the modes we need?
      for (size_t m = 0; m < found_mode.size(); ++m) {
        if (found_mode[m] == 0) {
          throw YamlException(std::string("Didn't find mode '") +
            aerosol_config.h_aerosol_modes(m).name() +
            std::string("' in aerosols section."));
        } else {
          for (size_t a = 0; a < found_aerosol[m].size(); ++a) {
            if (found_aerosol[m][a] == 0) {
              auto mode_species = aerosol_config.aerosol_species_for_mode(m);
              throw YamlException(std::string("Didn't find aerosol '") +
                mode_species[a].symbol() + std::string("' in aerosols section."));
            }
          }
        }
      }
    } else {
      throw YamlException("Did not find a valid aerosols section!");
      exit(1);
    }

    // gases section
    pw.gas_mmrs.resize(aerosol_config.h_gas_species.size());
    std::vector<int> found_gas(pw.gas_mmrs.size(), 0);
    if (root["gases"] and root["gases"].IsMap()) {
      auto gases = root["gases"];
      for (auto iter = gases.begin(); iter != gases.end(); ++iter) {
        auto gas_name = iter->first.as<std::string>();
        auto gas = iter->second;
        size_t index;
        for (index = 0; index < aerosol_config.h_gas_species.size(); ++index) {
          if (gas_name == aerosol_config.h_gas_species(index).symbol()) {
            found_gas[index] = 1;
            break;
          }
        }
        if (index < pw.gas_mmrs.size()) {
          pw.gas_mmrs[index] = gas.as<Real>();
        } else {
          throw YamlException(std::string("Found invalid gas '") +
            gas_name + "' in the gases section!");
        }
      }

      // Did we find all the gases we need?
      for (size_t g = 0; g < found_gas.size(); ++g) {
        if (found_gas[g] == 0) {
          throw YamlException(std::string("Didn't find gas '") +
            aerosol_config.h_gas_species(g).symbol() +
            std::string("' in gases section."));
        }
      }
    } else {
      throw YamlException("Did not find a valid gases section!");
    }
  } catch (std::exception& e) {
    throw YamlException(e.what());
  }
  return pw;
}

