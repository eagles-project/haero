#include "parse_yaml.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cstdarg>
#include <set>

namespace {

using Real = haero::Real;
using ModalAerosolConfig = haero::ModalAerosolConfig;
using YamlException = skywalker::YamlException;

std::vector<Real> parse_value_array(const std::string& name,
                                    const YAML::Node& value) {
  size_t len = value.size();
  if (len == 3) {
    Real value0 = value[0].as<Real>();
    Real value1 = value[1].as<Real>();
    Real value2 = value[2].as<Real>();
    if (not(value0 < value1)) {  // this should be true in any case!
      throw YamlException(
          std::string("Invalid values for '") + name +
          std::string("': second must be greater than first."));
    }
    if (not(value1 < value2)) {  // [start, stop, step]
      len = static_cast<size_t>((value1 - value0) / value2) + 1;
      std::vector<Real> values(len);
      for (int j = 0; j < len; ++j) {
        values[j] = value0 + j * value2;
      }
      return values;
    } else {
      return std::vector<Real>({value0, value1, value2});
    }
  } else {
    return value.as<std::vector<Real>>();
  }
}

void parse_atm_ensemble_params(const YAML::Node& root,
                               std::map<std::string, std::vector<Real>>& params) {
  for (auto iter = root.begin(); iter != root.end(); ++iter) {
    auto param = iter->second;
    if (param.IsSequence()) {
      auto param_name = std::string("atmosphere.") + iter->first.as<std::string>();
      params[param_name] = parse_value_array(param_name, param);
    } else {
      throw YamlException(
          std::string("Parameter 'atmosphere:") + iter->first.as<std::string>() +
          std::string("' is not a sequence of values!\n"));
    }
  }
}

void parse_aero_ensemble_params(const ModalAerosolConfig& aerosol_config,
                                const YAML::Node& root,
                                std::map<std::string, std::vector<Real>>& params) {
  for (auto riter = root.begin(); riter != root.end(); ++riter) {
    auto mode_name = riter->first.as<std::string>();
    int mode_index = aerosol_config.aerosol_mode_index(mode_name, false);
    auto mode = riter->second;
    if (mode_index != -1) {  // it's an aerosol mode!
      for (auto miter = mode.begin(); miter != mode.end(); ++miter) {
        auto group_name = miter->first.as<std::string>();
        auto group = miter->second;
        if (group.IsMap()) {
          if ((group_name != "cloud") and (group_name != "interstitial")) {
            continue;
          }
          for (auto giter = group.begin(); giter != group.end(); ++giter) {
            // Is this a valid aerosol species?
            auto aero_name = giter->first.as<std::string>();
            // Find the aerosol index within this species.
            int aero_index = aerosol_config.aerosol_species_index(
                mode_index, aero_name, false);
            if (aero_index == -1) {
              throw YamlException(
                  std::string("Found invalid aerosol species '") + aero_name +
                  std::string("' within mode '") + mode_name +
                  std::string("' in the ensemble section!"));
            }
            auto mmr_name = std::string("aerosols.") + group_name +
                            std::string(".") + mode_name + std::string(".") +
                            aero_name; // also works for number_conc
            auto aero_species = giter->second;
            params[mmr_name] = parse_value_array(mmr_name, aero_species);
          }
        } else {
          throw YamlException(
              std::string("Parameter 'aerosols:") + mode_name + std::string(":") +
              group_name + std::string("' is not a valid aerosol mode!\n"));
        }
      }
    } else {
      throw YamlException(
          std::string("Parameter 'aerosols:") + mode_name +
          std::string("' is not a map!\n"));
    }
  }
}

void parse_gas_ensemble_params(const ModalAerosolConfig& aerosol_config,
                               const YAML::Node& root,
                               std::map<std::string, std::vector<Real>>& params) {
  for (auto iter = root.begin(); iter != root.end(); ++iter) {
    auto param = iter->second;
    if (param.IsSequence()) {
      auto gas_name = iter->first.as<std::string>();
      int gas_index = aerosol_config.gas_index(gas_name, false);
      if (gas_index == -1) {
        throw YamlException(
          std::string("Parameter 'gases:") + gas_name +
          std::string("' is not a valid gas species!\n"));
      }
      gas_name = std::string("gases.") + gas_name;
      params[gas_name] = parse_value_array(gas_name, param);
    } else {
      throw YamlException(
          std::string("Parameter 'gases:") + iter->first.as<std::string>() +
          std::string("' is not a sequence!\n"));
    }
  }
}

}

namespace skywalker {

ParameterWalk parse_yaml(const haero::ModalAerosolConfig& aerosol_config,
                         const std::string& filename) {

  ParameterWalk pw(aerosol_config);
  try {
    auto root = YAML::LoadFile(filename);

    // process section
    if (root["process"] and root["process"].IsMap()) {
      auto process = root["process"];
      if (process["haero"]) {
        pw.process = process["haero"].as<std::string>();
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
        pw.ref_input.dt = ts["dt"].as<Real>();
      } else {
        throw YamlException("'dt' not found in timestepping section!\n");
      }
      if (ts["total_time"]) {
        pw.ref_input.total_time = ts["total_time"].as<Real>();
      } else {
        throw YamlException(
            "'total_time' not found in timestepping section!\n");
      }
    } else {
      throw YamlException("Did not find a valid timestepping section!\n");
      exit(1);
    }

    // enѕemble section
    if (root["ensemble"] and root["ensemble"].IsMap()) {
      auto ensemble = root["ensemble"];
      for (auto eiter = ensemble.begin(); eiter != ensemble.end(); ++eiter) {
        auto group_name = eiter->first.as<std::string>();
        if ((group_name != "atmosphere") and (group_name != "gases") and
            (group_name != "aerosols")) {
          continue;
        }
        auto group = eiter->second;
        if (group_name == "atmosphere") {
          parse_atm_ensemble_params(group, pw.ensemble);
        } else if (group_name == "gases") {
          parse_gas_ensemble_params(aerosol_config, group, pw.ensemble);
        } else {
          parse_aero_ensemble_params(aerosol_config, group, pw.ensemble);
        }
      }
    } else {
      throw YamlException("Did not find a valid ensemble section!\n");
    }

    // atmosphere section
    if (root["atmosphere"] and root["atmosphere"].IsMap()) {
      auto atm = root["atmosphere"];
      if (atm["temperature"]) {
        pw.ref_input.temperature = atm["temperature"].as<Real>();
      } else {
        throw YamlException(
            "Did not find 'temperature' in the atmosphere section!\n");
      }
      if (atm["pressure"]) {
        pw.ref_input.pressure = atm["pressure"].as<Real>();
      } else {
        throw YamlException(
            "Did not find 'pressure' in the atmosphere section!\n");
      }
      if (atm["relative_humidity"]) {
        pw.ref_input.relative_humidity = atm["relative_humidity"].as<Real>();
      } else {
        throw YamlException(
            "Did not find 'relative_humidity' in the atmosphere section!\n");
      }
      if (atm["height"]) {
        pw.ref_input.height = atm["height"].as<Real>();
      } else {
        throw YamlException(
            "Did not find 'height' in the atmosphere section!\n");
      }
      if (atm["hydrostatic_dp"]) {
        pw.ref_input.hydrostatic_dp = atm["hydrostatic_dp"].as<Real>();
      } else {
        throw YamlException(
            "Did not find 'hydrostatic_dp' in the atmosphere section!\n");
      }
      if (atm["planetary_boundary_layer_height"]) {
        pw.ref_input.planetary_boundary_layer_height =
            atm["planetary_boundary_layer_height"].as<Real>();
      } else {
        throw YamlException(
            "Did not find 'planetary_boundary_layer_height' in the atmosphere "
            "section!\n");
      }
    } else {
      throw YamlException("Did not find a valid atmosphere section!\n");
    }

    // aerosols section
    int num_modes = aerosol_config.num_modes();
    pw.ref_input.interstitial_number_concs.resize(num_modes);
    pw.ref_input.cloud_number_concs.resize(num_modes);
    pw.ref_input.interstitial_aero_mmrs.resize(
        aerosol_config.num_aerosol_populations);
    pw.ref_input.cloud_aero_mmrs.resize(aerosol_config.num_aerosol_populations,
                                        0.0);
    if (root["aerosols"] and root["aerosols"].IsMap()) {
      auto aerosols = root["aerosols"];
      std::vector<int> found_mode(num_modes, 0);
      for (auto miter = aerosols.begin(); miter != aerosols.end(); ++miter) {
        auto mode_name = miter->first.as<std::string>();
        auto mode = miter->second;
        if (mode.IsMap()) {
          // Determine the mode's index.
          int mode_index = aerosol_config.aerosol_mode_index(mode_name, false);
          if (mode_index != -1) {
            found_mode[mode_index] = 1;
          }

          if (mode_index < num_modes) {
            for (auto giter = mode.begin(); giter != mode.end(); ++giter) {
              auto group_name = giter->first.as<std::string>();
              if ((group_name != "cloud") and (group_name != "interstitial")) {
                continue;
              }
              auto group = giter->second;

              // Get the initial data for the aerosol species in this mode.
              auto mode_species =
                aerosol_config.aerosol_species_for_mode(mode_index);
              bool found_number_conc = false;
              std::vector<int> found_aerosol(mode_species.size());
              for (auto aiter = group.begin(); aiter != group.end(); ++aiter) {
                auto aero_name = aiter->first.as<std::string>();
                auto aero_species = aiter->second;
                if (aero_name == "number_conc") {  // number conc, not species name!
                  if (group_name == "interstitial") {
                    pw.ref_input.interstitial_number_concs[mode_index] = aero_species.as<Real>();
                  } else {
                    pw.ref_input.cloud_number_concs[mode_index] = aero_species.as<Real>();
                  }
                  found_number_conc = true;
                } else {
                  // Find the aerosol index within this species.
                  int aero_index = aerosol_config.aerosol_species_index(
                      mode_index, aero_name, false);
                  if (aero_index == -1) {
                    throw YamlException(
                        std::string("Found invalid aerosol species '") +
                        aero_name + std::string("' within mode '") + mode_name +
                        std::string("' in the aerosols section!"));
                  } else {
                    found_aerosol[aero_index] = 1;
                    int pop_index =
                      aerosol_config.population_index(mode_index, aero_index);
                    if (group_name == "interstitial") {
                      pw.ref_input.interstitial_aero_mmrs[pop_index] =
                        aero_species.as<Real>();
                    } else {
                      pw.ref_input.cloud_aero_mmrs[pop_index] =
                        aero_species.as<Real>();
                    }
                  }
                  if (not found_number_conc) {
                    throw YamlException(
                        std::string("Did not find 'cloud_number_conc' in ") +
                        group_name + std::string("'") +
                        mode_name + std::string("' mode!"));
                  }
                }
              }
              for (size_t a = 0; a < found_aerosol.size(); ++a) {
                if (found_aerosol[a] == 0) {
                  auto mode_species = aerosol_config.aerosol_species_for_mode(mode_index);
                  throw YamlException(std::string("Didn't find aerosol '") +
                      mode_species[a].symbol() +
                      std::string("' in mode '") + mode_name +
                      std::string("' of aerosols:") + group_name + std::string(" section."));
                }
              }
            }
          } else {
            throw YamlException(std::string("Found invalid aerosol mode '") +
                mode_name +
                std::string("' in the aerosols section!"));
          }
        }
      }

      // Did we find all the modes we need?
      for (size_t m = 0; m < found_mode.size(); ++m) {
        auto mode_name = aerosol_config.h_aerosol_modes(m).name();
        if (found_mode[m] == 0) {
          throw YamlException(std::string("Didn't find mode '") + mode_name +
              std::string("' in aerosols section."));
        }
      }
    } else {
      throw YamlException("Did not find a valid aerosols section!");
      exit(1);
    }

    // gases section
    pw.ref_input.gas_mmrs.resize(aerosol_config.num_gases());
    std::vector<int> found_gas(pw.ref_input.gas_mmrs.size(), 0);
    if (root["gases"] and root["gases"].IsMap()) {
      auto gases = root["gases"];
      for (auto iter = gases.begin(); iter != gases.end(); ++iter) {
        auto gas_name = iter->first.as<std::string>();
        int gas_index = aerosol_config.gas_index(gas_name, false);
        if (gas_index != -1) {  // valid gas
          found_gas[gas_index] = 1;
          auto gas = iter->second;
          pw.ref_input.gas_mmrs[gas_index] = gas.as<Real>();
        } else {
          throw YamlException(std::string("Found invalid gas '") + gas_name +
                              "' in the gases section!");
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

}  // namespace skywalker
