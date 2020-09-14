#include <algorithm>
#include <yaml-cpp/yaml.h>
#include "read_yaml.hpp"

namespace haero {

namespace {

std::vector<Mode> read_modes(const YAML::Node& root)
{
  std::vector<Mode> modes;
  if (root["modes"] and root["modes"].IsMap())
  {
    auto node = root["modes"];
    for (auto iter = node.begin(); iter != node.end(); ++iter)
    {
      std::string name = iter->first.as<std::string>();
      auto mnode = iter->second;
      if (not mnode["D_min"])
        throw YamlException("mode entry has no minimum diameter (D_min).");
      else if (not mnode["D_max"])
        throw YamlException("mode entry has no maximum diameter (D_max).");
      else if (not mnode["sigma"])
        throw YamlException("mode entry has no geometric stddev (sigma).");
      else
      {
        modes.push_back(Mode(name,
                             mnode["D_min"].as<Real>(),
                             mnode["D_max"].as<Real>(),
                             mnode["sigma"].as<Real>()));
      }
    }
  }
  else
    throw YamlException("No modes section was found!");
  return modes;
}

std::vector<Species> read_aerosol_species(const YAML::Node& root)
{
  std::vector<Species> species;
  if (root["aerosols"] and root["aerosols"].IsMap())
  {
    auto node = root["aerosols"];
    for (auto iter = node.begin(); iter != node.end(); ++iter)
    {
      std::string symbol = iter->first.as<std::string>();
      auto snode = iter->second;
      if (not snode["name"])
        throw YamlException("aerosol species entry has no name.");
      else
      {
        Species s;
        s.symbol = symbol;
        s.name = snode["name"].as<std::string>();
        species.push_back(s);
      }
    }
  }
  else
    throw YamlException("No aerosols section was found!");
  return species;
}

std::vector<Species> read_gas_species(const YAML::Node& root)
{
  std::vector<Species> species;
  if (root["gases"] and root["gases"].IsMap())
  {
    auto node = root["gases"];
    for (auto iter = node.begin(); iter != node.end(); ++iter)
    {
      std::string symbol = iter->first.as<std::string>();
      auto snode = iter->second;
      if (not snode["name"])
        throw YamlException("gas species entry has no name.");
      else
      {
        Species s;
        s.symbol = symbol;
        s.name = snode["name"].as<std::string>();
        species.push_back(s);
      }
    }
  }
  else
    throw YamlException("No gases section was found!");
  return species;
}

InitialConditions read_initial_conditions(const std::vector<Mode>& modes,
                                          const std::vector<Species>& aerosols,
                                          const std::vector<Species>& gases,
                                          const YAML::Node& root)
{
  // First read all the values directly into memory.
  InitialConditions ics;
  if (root["initial_conditions"] and root["initial_conditions"].IsMap())
  {
    auto n = root["initial_conditions"];
    if (n["aerosols"] and n["aerosols"].IsMap())
    {
      auto a = n["aerosols"];
      ics.aerosols.resize(modes.size());

      // Iterate over modes, validating their names.
      for (auto a_iter = a.begin(); a_iter != a.end(); ++a_iter)
      {
        // Fetch the variable and its initial condition.
        std::string mode_name = a_iter->first.as<std::string>();
        auto mode_iter = std::find_if(modes.begin(), modes.end(),
          [&](const Mode& m) { return m.name == mode_name; });
        if (mode_iter == modes.end())
          throw YamlException("Invalid mode specified in aerosol initial conditions!");

        // The mode is valid, so compute its index.
        size_t mode_index = mode_iter - modes.begin();

        // Is the mode entry a map?
        if (a[mode_name] and a[mode_name].IsMap())
        {
          auto m = a[mode_name];

          // Iterate over the aerosol ics for the mode. These ICs are mass
          // fractions and must sum to 1.
          Real mass_frac_sum = 0.0;
          for (auto m_iter = m.begin(); m_iter != m.end(); ++m_iter)
          {
            auto aero_symbol = m_iter->first.as<std::string>();
            auto aero_iter = std::find_if(aerosols.begin(), aerosols.end(),
              [&](const Species& s) { return s.symbol == aero_symbol; });
            if (aero_iter == aerosols.end())
              throw YamlException("Invalid aerosol found in initial conditions!");
            auto aero_ic = m_iter->second.as<Real>();
            if (aero_ic < 0.0)
              throw YamlException("Negative aerosol mass fraction found!");
            else if (aero_ic > 1.0)
              throw YamlException("Invalid aerosol mass fraction (> 1) found!");
            ics.aerosols[mode_index][aero_symbol] = aero_ic;
            mass_frac_sum += aero_ic;
          }
          if (std::abs(mass_frac_sum - 1.0) > 1e-12) // TODO: Use better machinery
            throw YamlException("Aerosol mass fractions don't sum to 1!");
        }
        else
          throw YamlException("No aerosol initial conditions given for mode!");
      }
    }
    else
      throw YamlException("No aerosols subsection was found in initial_conditions!");
    if (n["gases"] and n["gases"].IsMap())
    {
      auto g = n["gases"];

      // Iterate over the gas ICs (mole fractions).
      for (auto g_iter = g.begin(); g_iter != g.end(); ++g_iter)
      {
        auto gas_symbol = g_iter->first.as<std::string>();
        auto gas_iter = std::find_if(gases.begin(), gases.end(),
              [&](const Species& s) { return s.symbol == gas_symbol; });
        if (gas_iter == gases.end())
          throw YamlException("Invalid gas found in initial conditions!");
        auto gas_ic = g_iter->second.as<Real>();
        if (gas_ic < 0.0)
          throw YamlException("Negative gas mole fraction found!");
        else if (gas_ic > 1.0)
          throw YamlException("Invalid gas mole fraction (> 1) found!");
        ics.gases[gas_symbol] = gas_ic;
      }
    }
    else
      throw YamlException("No gases subsection was found in initial_conditions!");
    if (n["modes"] and n["modes"].IsMap())
    {
      auto m = n["modes"];

      // Iterate over the mode ICs (number densities).
      for (auto m_iter = m.begin(); m_iter != m.end(); ++m_iter)
      {
        auto mode_name = m_iter->first.as<std::string>();
        auto mode_iter = std::find_if(modes.begin(), modes.end(),
              [&](const Mode& mm) { return mm.name == mode_name; });
        if (mode_iter == modes.end())
          throw YamlException("Invalid gas found in initial conditions!");
        size_t mode_index = mode_iter - modes.begin();
        auto mode_ic = m_iter->second.as<Real>();
        if (mode_ic < 0.0)
          throw YamlException("Negative mode number density found!");
        ics.modes[mode_index] = mode_ic;
      }
    }
    else
      throw YamlException("No modes subsection was found in initial_conditions!");
  }
  else
    throw YamlException("No initial_conditions section was found!");

  return ics;
}

AerosolProcesses read_physics_settings(const YAML::Node& root)
{
  AerosolProcesses settings = {true, true, true, true, true, true, true};
  if (root["physics"] and root["physics"].IsMap())
  {
    auto node = root["physics"];
    if (node["growth"] && not node["growth"].as<bool>())
      settings.growth = false;
    if (node["gas_chemistry"] && not node["gas_chemistry"].as<bool>())
      settings.gas_chemistry = false;
    if (node["cloud_chemistry"] && not node["cloud_chemistry"].as<bool>())
      settings.cloud_chemistry = false;
    if (node["gas_aerosol_exchange"] && not node["gas_aerosol_exchange"].as<bool>())
      settings.gas_aerosol_exchange = false;
    if (node["mode_merging"] && not node["mode_merging"].as<bool>())
      settings.mode_merging = false;
    if (node["nucleation"] && not node["nucleation"].as<bool>())
      settings.nucleation = false;
    if (node["coagulation"] && not node["coagulation"].as<bool>())
      settings.coagulation = false;
  }
  return settings;
}

AtmosphericConditions read_atmosphere(const YAML::Node& root)
{
  AtmosphericConditions atm;
  if (root["atmosphere"] and root["atmosphere"].IsMap())
  {
    auto a = root["atmosphere"];

    // Which model are we using?
    if (a["model"])
    {
      auto model = a["model"].as<std::string>();
      if (model == "uniform")
        atm.model = AtmosphericConditions::uniform;
      else if (model != "hydrostatic")
        atm.model = AtmosphericConditions::hydrostatic;
      else
        throw YamlException("Invalid model for atmospheric conditions!");
    }
    else
      throw YamlException("No model found in atmosphere section!");

    // Now fetch the relevant parameters.
    if (atm.model == AtmosphericConditions::uniform)
    {
      if (a["mu"])
        atm.params.uniform.mu = a["mu"].as<Real>();
      else
        throw YamlException("Missing mean molecular weight (mu) for uniform atm!");
      if (a["H"])
        atm.params.uniform.H = a["H"].as<Real>();
      else
        throw YamlException("Missing scaled height (H) for uniform atm!");
      if (a["p0"])
        atm.params.uniform.p0 = a["p0"].as<Real>();
      else
        throw YamlException("Missing pressure (p0) for uniform atm!");
      if (a["T0"])
        atm.params.uniform.T0 = a["T0"].as<Real>();
      else
        throw YamlException("Missing temperature (T0) for uniform atm!");
      if (a["phi0"])
        atm.params.uniform.phi0 = a["phi0"].as<Real>();
      else
        throw YamlException("Missing relative humidity (phi0) for uniform atm!");
      if (a["N0"])
        atm.params.uniform.N0 = a["N0"].as<Real>();
      else
        throw YamlException("Missing cloud fraction (N0) for uniform atm!");
    }
    else if (atm.model == AtmosphericConditions::hydrostatic)
    {
      // FIXME
    }
  }
  else
    throw YamlException("No atmosphere section was found!");
  return atm;
}

GridParams read_grid_params(const YAML::Node& root)
{
  GridParams grid;
  if (root["grid"] and root["grid"].IsMap())
  {
    auto g = root["grid"];
    if (g["num_columns"])
    {
      auto num_columns = g["num_columns"].as<int>();
      if (num_columns <= 0)
        throw YamlException("Non-positive num_columns found in grid section!");
      grid.num_columns = num_columns;
    }
    else
      throw YamlException("No num_columns found in grid section!");
    if (g["num_levels"])
    {
      auto num_levels = g["num_levels"].as<int>();
      if (num_levels <= 0)
        throw YamlException("Non-positive num_levels found in grid section!");
      grid.num_levels = num_levels;
    }
    else
      throw YamlException("No num_levels found in grid section!");
  }
  else
    throw YamlException("No grid section was found!");
  return grid;
}

SimpleChemistryModel read_chemistry_model(const YAML::Node& root)
{
  SimpleChemistryModel chem_model;
  if (root["chemistry"] and root["chemistry"].IsMap())
  {
    auto chem = root["chemistry"];
    if (chem["production"] and chem["production"].IsMap())
    {
      auto prod = chem["production"];
      for (auto iter = prod.begin(); iter != prod.end(); ++iter)
      {
        auto gas = iter->first.as<std::string>();
        auto rate = iter->second.as<Real>();
        chem_model.production_rates[gas] = rate;
      }
    }
  }
  return chem_model;
}

std::vector<SimulationParams> read_simulation_params(const YAML::Node& root)
{
  SimulationParams params0;
  std::vector<Real> dts;
  if (root["simulation"] and root["simulation"].IsMap())
  {
    auto sim = root["simulation"];
    if (sim["timestep"])
    {
      auto timestep = root["timestep"];
      if (timestep.IsSequence())
      {
        for (size_t i = 0; i < timestep.size(); ++i)
          dts.push_back(timestep[i].as<Real>());
      }
      else
        params0.dt = timestep.as<Real>();
    }
    if (sim["duration"])
      params0.duration = root["duration"].as<Real>();
    else
      throw YamlException("No duration found in simulation section.");

    // Output parameters.
    params0.output_dir = std::string(".");
    params0.output_prefix = std::string("haero_output");
    params0.output_freq = -1;
    if (sim["output"])
    {
      auto node = sim["output"];
      if (node["prefix"])
        params0.output_prefix = node["prefix"].as<std::string>();
      if (node["directory"])
        params0.output_dir = node["directory"].as<std::string>();
      if (node["frequency"])
        params0.output_freq = node["frequency"].as<int>();
    }
  }
  else
    throw YamlException("No simulation section was found!");

  // If we have more than one timestep, we take the cartesian product of
  // dts with params0. Otherwise we just return a vector containing params0.
  if (not dts.empty()) // multiple timesteps
  {
    std::vector<SimulationParams> params;
    for (size_t i = 0; i < dts.size(); ++i)
    {
      params[i] = params0;
      params[i].dt = dts[i];
    }
    return params;
  }
  else
    return std::vector<SimulationParams>(1, params0);
}

} // anonymous namespace

std::vector<SimulationInput> read_yaml(const std::string& filename)
{
  // Try to load the thing.
  try
  {
    auto root = YAML::LoadFile(filename);

    // Okay, it loaded. Now we fill our vector with simulation input. Because
    // this file might contain data for an ensemble and not a single simulation,
    // we first create an "input0" struct to hold all the data common to all
    // members of the ensemble. In the case that our file holds data for only
    // a single simulation, our input vector ends up containing just input0.
    SimulationInput input0;
    input0.modes = read_modes(root);
    input0.aerosols = read_aerosol_species(root);
    input0.gases = read_gas_species(root);
    input0.physics = read_physics_settings(root);
    input0.atmosphere = read_atmosphere(root);
    input0.grid = read_grid_params(root);
    input0.initial_conditions = read_initial_conditions(input0.modes,
                                                        input0.aerosols,
                                                        input0.gases,
                                                        root);
    input0.chemistry = read_chemistry_model(root);

    // TODO: We fetch ensemble parameters here.
    int ensemble_size = 1;

    // We fetch the simulation parameters, which include all fixed timestep
    // sizes.
    std::vector<SimulationParams> sim_params = read_simulation_params(root);
    int num_dts = sim_params.size();

    int num_members = ensemble_size * num_dts;
    if (num_members > 1)
    {
      std::vector<SimulationInput> input;
      for (int i = 0; i < ensemble_size; ++i)
      {
        // TODO: ensemble perturbations go here.
        for (int j = 0; j < sim_params.size(); ++j)
          input0.simulation = sim_params[j];
      }
      return input;
    }
    else
      return std::vector<SimulationInput>(1, input0);
  }
  catch (YAML::BadFile& e)
  {
    throw YamlException(e.what());
  }
  catch (YAML::ParserException& e)
  {
    throw YamlException(e.what());
  }
}

} // namespace haero
