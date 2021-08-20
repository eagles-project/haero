#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cstdarg>
#include <set>

#include "haero/modal_aerosol_config.hpp"
#include "skywalker.hpp"

namespace {

using Real = skywalker::Real;
using ModalAerosolConfig = haero::ModalAerosolConfig;
using ParameterWalk = skywalker::ParameterWalk;
using YamlException = skywalker::YamlException;

std::vector<Real> parse_value_array(const std::string& name,
                                    const YAML::Node& value) {
  size_t len = value.size();
  if (len == 3) {
    Real value0 = value[0].as<Real>();
    Real value1 = value[1].as<Real>();
    Real value2 = value[2].as<Real>();
    if (not(value0 < value1)) {  // this should be true in any case!
      throw YamlException(std::string("Invalid values for '") + name +
                          std::string("': second must be greater than first."));
    }
    if (not(value1 < value2)) {  // [start, stop, step]
      len = static_cast<size_t>(std::ceil((value1 - value0) / value2) + 1);
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

void parse_atm_ensemble_params(
    const YAML::Node& atm, std::map<std::string, std::vector<Real>>& params) {
  for (auto iter : atm) {
    auto param = iter.second;
    if (not param.IsSequence()) {
      throw YamlException(std::string("Parameter 'atmosphere:") +
                          iter.first.as<std::string>() +
                          std::string("' is not a sequence of values!\n"));
    }
    auto param_name = std::string("atmosphere.") + iter.first.as<std::string>();
    params[param_name] = parse_value_array(param_name, param);
  }
}

void parse_aero_ensemble_params(
    const ModalAerosolConfig& aerosol_config, const YAML::Node& aero,
    std::map<std::string, std::vector<Real>>& params) {
  for (auto aiter : aero) {
    auto mode_name = aiter.first.as<std::string>();
    int mode_index = aerosol_config.aerosol_mode_index(mode_name, false);
    auto mode = aiter.second;
    if (mode_index == -1) {  // not an aerosol mode!
      throw YamlException(std::string("Parameter 'aerosols:") + mode_name +
                          std::string("' is not a valid aerosol mode!\n"));
    }
    for (auto miter : mode) {
      auto group_name = miter.first.as<std::string>();
      auto group = miter.second;
      if (not group.IsMap()) {
        throw YamlException(std::string("Parameter 'aerosols:") + mode_name +
                            std::string(":") + group_name +
                            std::string("' is not a map!\n"));
      }
      if ((group_name != "cloud") and (group_name != "interstitial")) {
        throw YamlException(
            std::string("Іnvalid entry in parameter 'aerosols:") + mode_name +
            std::string(": ") + group_name);
      }
      for (auto giter : group) {
        // Is this a valid aerosol species?
        auto aero_name = giter.first.as<std::string>();
        // Find the aerosol index within this species.
        int aero_index =
            aerosol_config.aerosol_species_index(mode_index, aero_name, false);
        if (aero_index == -1) {
          throw YamlException(std::string("Found invalid aerosol species '") +
                              aero_name + std::string("' within mode '") +
                              mode_name +
                              std::string("' in the ensemble section!"));
        }
        auto mmr_name = std::string("aerosols.") + group_name +
                        std::string(".") + mode_name + std::string(".") +
                        aero_name;  // also works for num_mix_ratio
        auto aero_species = giter.second;
        params[mmr_name] = parse_value_array(mmr_name, aero_species);
      }
    }
  }
}

void parse_gas_ensemble_params(
    const ModalAerosolConfig& aerosol_config, const YAML::Node& gases,
    std::map<std::string, std::vector<Real>>& params) {
  for (auto iter : gases) {
    auto param = iter.second;
    if (not param.IsSequence()) {
      throw YamlException(std::string("Parameter 'gases:") +
                          iter.first.as<std::string>() +
                          std::string("' is not a sequence!\n"));
    }
    auto gas_name = iter.first.as<std::string>();
    int gas_index = aerosol_config.gas_index(gas_name, false);
    if (gas_index == -1) {
      throw YamlException(std::string("Parameter 'gases:") + gas_name +
                          std::string("' is not a valid gas species!\n"));
    }
    gas_name = std::string("gases.") + gas_name;
    params[gas_name] = parse_value_array(gas_name, param);
  }
}

void parse_process_section(const YAML::Node& process, ParameterWalk& pw) {
  // Parse the process based on the model implementation (e.g. "mam" or "haero")
  if (not process[pw.model_impl]) {
    throw YamlException(std::string("'") + pw.model_impl +
                        std::string("' entry not found in process section!"));
  }
  const auto& model_impl = process[pw.model_impl];
  if (not model_impl.IsMap()) {
    throw YamlException(std::string("'") + pw.model_impl +
                        "' in process section must be a map!");
  }
  if (not model_impl["name"]) {
    throw YamlException(std::string("'") + pw.model_impl +
                        "' in process section must have a 'name' entry!");
  }
  pw.process = model_impl["name"].as<std::string>();

  // Parse process-specific parameters.
  if (model_impl["params"]) {
    const auto& params = model_impl["params"];
    if (not params.IsMap()) {
      throw YamlException(std::string("'params' in process:") +
                          pw.model_impl + std::string(" section must be a map!"));
    }
    for (const auto& param: params) {
      const auto& name = param.first.as<std::string>();
      const auto& value = param.second.as<std::string>();
      pw.process_params[name] = value;
    }
  }
}

void parse_timestepping_section(const YAML::Node& ts, ParameterWalk& pw) {
  if (not ts["dt"]) {
    throw YamlException("'dt' not found in timestepping section!");
  }
  pw.ref_input.dt = ts["dt"].as<Real>();

  if (not ts["total_time"]) {
    throw YamlException("'total_time' not found in timestepping section!");
  }
  pw.ref_input.total_time = ts["total_time"].as<Real>();
}

void parse_ensemble_section(const YAML::Node& ensemble,
                            const ModalAerosolConfig& aerosol_config,
                            ParameterWalk& pw) {
  for (auto eiter : ensemble) {
    auto group_name = eiter.first.as<std::string>();
    if ((group_name != "atmosphere") and (group_name != "gases") and
        (group_name != "aerosols")) {
      continue;
    }
    auto group = eiter.second;
    if (group_name == "atmosphere") {
      parse_atm_ensemble_params(group, pw.ensemble);
    } else if (group_name == "gases") {
      parse_gas_ensemble_params(aerosol_config, group, pw.ensemble);
    } else {
      parse_aero_ensemble_params(aerosol_config, group, pw.ensemble);
    }
  }
}

void parse_atmosphere_section(const YAML::Node& atm, ParameterWalk& pw) {
  if (not atm["temperature"]) {
    throw YamlException(
        "Did not find 'temperature' in the atmosphere section!\n");
  }
  pw.ref_input.temperature = atm["temperature"].as<Real>();

  if (not atm["pressure"]) {
    throw YamlException("Did not find 'pressure' in the atmosphere section!\n");
  }
  pw.ref_input.pressure = atm["pressure"].as<Real>();

  if (not atm["relative_humidity"]) {
    throw YamlException(
        "Did not find 'relative_humidity' in the atmosphere section!\n");
  }
  pw.ref_input.vapor_mixing_ratio = atm["vapor_mixing_ratio"].as<Real>();

  if (not atm["height"]) {
    throw YamlException("Did not find 'height' in the atmosphere section!\n");
  }
  pw.ref_input.height = atm["height"].as<Real>();

  if (not atm["hydrostatic_dp"]) {
    throw YamlException(
        "Did not find 'hydrostatic_dp' in the atmosphere section!\n");
  }
  pw.ref_input.hydrostatic_dp = atm["hydrostatic_dp"].as<Real>();

  if (not atm["planetary_boundary_layer_height"]) {
    throw YamlException(
        "Did not find 'planetary_boundary_layer_height' in the atmosphere "
        "section!\n");
  }
  pw.ref_input.planetary_boundary_layer_height =
      atm["planetary_boundary_layer_height"].as<Real>();
}

void parse_aerosols_section(const YAML::Node& aerosols,
                            const ModalAerosolConfig& aerosol_config,
                            ParameterWalk& pw) {
  int num_modes = aerosol_config.num_modes();

  // Initialize reference input data.
  pw.ref_input.interstitial_number_mix_ratios.resize(num_modes);
  pw.ref_input.cloud_number_mix_ratios.resize(num_modes);
  pw.ref_input.interstitial_aero_mmrs.resize(
      aerosol_config.num_aerosol_populations);
  pw.ref_input.cloud_aero_mmrs.resize(aerosol_config.num_aerosol_populations,
                                      0.0);

  std::vector<int> found_mode(num_modes, 0);
  for (auto miter : aerosols) {
    auto mode_name = miter.first.as<std::string>();
    auto mode = miter.second;

    // Determine the mode's index.
    int mode_index = aerosol_config.aerosol_mode_index(mode_name, false);
    if (mode_index == -1) {
      throw YamlException(std::string("Found invalid aerosol mode '") +
                          mode_name +
                          std::string("' in the aerosols section!"));
    }
    found_mode[mode_index] = 1;

    if (not mode.IsMap()) {
      throw YamlException(std::string("aerosols:") + mode_name +
                          std::string(" is not a map!"));
    }

    EKAT_REQUIRE(mode_index < num_modes);  // aerosol_config guarantees this.
    for (auto giter : mode) {
      auto group_name = giter.first.as<std::string>();
      if ((group_name != "cloud") and (group_name != "interstitial")) {
        continue;
      }
      auto group = giter.second;

      // Get the initial data for the aerosol species in this mode.
      auto mode_species = aerosol_config.aerosol_species_for_mode(mode_index);
      bool found_num_mix_ratio = false;
      std::vector<int> found_aerosol(mode_species.size());
      for (auto aiter : group) {
        auto aero_name = aiter.first.as<std::string>();
        auto aero_species = aiter.second;
        if (aero_name == "number_mix_ratio") {  // number mix ratio, not species
          if (group_name == "interstitial") {
            pw.ref_input.interstitial_number_mix_ratios[mode_index] =
                aero_species.as<Real>();
          } else {
            pw.ref_input.cloud_number_mix_ratios[mode_index] =
                aero_species.as<Real>();
          }
          found_num_mix_ratio = true;
        } else {
          // Find the aerosol index within this species.
          int aero_index = aerosol_config.aerosol_species_index(
              mode_index, aero_name, false);
          if (aero_index == -1) {
            throw YamlException(std::string("Found invalid aerosol species '") +
                                aero_name + std::string("' within mode '") +
                                mode_name +
                                std::string("' in the aerosols section!"));
          }
          found_aerosol[aero_index] = 1;
          int pop_index =
              aerosol_config.population_index(mode_index, aero_index);
          if (group_name == "interstitial") {
            pw.ref_input.interstitial_aero_mmrs[pop_index] =
                aero_species.as<Real>();
          } else {
            pw.ref_input.cloud_aero_mmrs[pop_index] = aero_species.as<Real>();
          }
          // Did we find a num_mix_ratio value for this mode?
          if (not found_num_mix_ratio) {
            throw YamlException(std::string("Did not find 'num_mix_ratio' in ") +
                                group_name + std::string("'") + mode_name +
                                std::string("' mode!"));
          }
        }
      }

      // Did we find all the aerosols we expect in this mode?
      for (size_t a = 0; a < found_aerosol.size(); ++a) {
        if (found_aerosol[a] == 0) {
          auto mode_species =
              aerosol_config.aerosol_species_for_mode(mode_index);
          throw YamlException(std::string("Didn't find aerosol '") +
                              mode_species[a].symbol() +
                              std::string("' in mode '") + mode_name +
                              std::string("' of aerosols:") + group_name +
                              std::string(" section."));
        }
      }
    }
  }

  // Did we find all the modes we need?
  for (size_t m = 0; m < found_mode.size(); ++m) {
    auto mode_name = aerosol_config.aerosol_modes[m].name();
    if (found_mode[m] == 0) {
      throw YamlException(std::string("Didn't find mode '") + mode_name +
                          std::string("' in aerosols section."));
    }
  }
}

void parse_gases_section(const YAML::Node& gases,
                         const ModalAerosolConfig& aerosol_config,
                         ParameterWalk& pw) {
  pw.ref_input.gas_mmrs.resize(aerosol_config.num_gases());

  std::vector<int> found_gas(pw.ref_input.gas_mmrs.size(), 0);
  for (auto iter : gases) {
    auto gas_name = iter.first.as<std::string>();
    int gas_index = aerosol_config.gas_index(gas_name, false);
    if (gas_index == -1) {  // invalid gas
      throw YamlException(std::string("Found invalid gas '") + gas_name +
                          "' in the gases section!");
    }
    found_gas[gas_index] = 1;
    auto gas = iter.second;
    pw.ref_input.gas_mmrs[gas_index] = gas.as<Real>();
  }

  // Did we find all the gases we need?
  for (size_t g = 0; g < found_gas.size(); ++g) {
    if (found_gas[g] == 0) {
      throw YamlException(std::string("Didn't find gas '") +
                          aerosol_config.gas_species[g].symbol() +
                          std::string("' in gases section."));
    }
  }
}

}  // namespace

namespace skywalker {

ParameterWalk load_ensemble(const haero::ModalAerosolConfig& aerosol_config,
                            const std::string& filename,
                            const std::string& model_impl) {
  ParameterWalk pw(aerosol_config, model_impl);
  try {
    auto root = YAML::LoadFile(filename);

    if (not(root["process"] and root["process"].IsMap())) {
      throw YamlException("Did not find a valid process section!");
    }
    parse_process_section(root["process"], pw);

    if (not(root["timestepping"] and root["timestepping"].IsMap())) {
      throw YamlException("Did not find a valid timestepping section!\n");
    }
    parse_timestepping_section(root["timestepping"], pw);

    if (not(root["ensemble"] and root["ensemble"].IsMap())) {
      throw YamlException("Did not find a valid ensemble section!\n");
    }
    parse_ensemble_section(root["ensemble"], aerosol_config, pw);

    if (not(root["atmosphere"] and root["atmosphere"].IsMap())) {
      throw YamlException("Did not find a valid atmosphere section!\n");
    }
    parse_atmosphere_section(root["atmosphere"], pw);

    if (not(root["aerosols"] and root["aerosols"].IsMap())) {
      throw YamlException("Did not find a valid aerosols section!");
    }
    parse_aerosols_section(root["aerosols"], aerosol_config, pw);

    if (not(root["gases"] and root["gases"].IsMap())) {
      throw YamlException("Did not find a valid gases section!");
    }
    parse_gases_section(root["gases"], aerosol_config, pw);
  } catch (std::exception& e) {
    throw YamlException(e.what());
  }
  return pw;
}

}  // namespace skywalker
