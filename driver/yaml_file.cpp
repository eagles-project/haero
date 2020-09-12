#include "driver/yaml_file.hpp"

namespace haero {

YamlFile::YamlFile(const std::string& filename)
{
  // Try to load the thing.
  try
  {
    _root = YAML::LoadFile(filename);
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

YamlFile::~YamlFile() {
}

std::vector<Mode> YamlFile::read_modes() const
{
  std::vector<Mode> modes;
  if (_root["modes"] and _root["modes"].IsMap())
  {
    auto node = _root["modes"];
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
  return modes;
}

std::vector<Species> YamlFile::read_aerosol_species() const
{
  std::vector<Species> species;
  if (_root["aerosols"] and _root["aerosols"].IsMap())
  {
    auto node = _root["aerosols"];
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
  return species;
}

std::vector<Species> YamlFile::read_gas_species() const
{
  std::vector<Species> species;
  if (_root["gases"] and _root["gases"].IsMap())
  {
    auto node = _root["gases"];
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
  return species;
}

namespace {

std::vector<std::map<std::string, Real> >
read_vector_from_map(const YAML::Node& node, const std::string& map_name)
{
  // First read all the values directly into memory.
  std::map<std::string, std::vector<Real> > values;
  if (node[map_name] and node[map_name].IsMap())
  {
    YAML::Node n = node[map_name];
    for (auto iter = n.begin(); iter != n.end(); ++iter)
    {
      // Fetch the variable and its initial condition.
      std::string var = iter->first.as<std::string>();
      YAML::Node val = iter->second;
      if (val.IsSequence())
      {
        for (size_t i = 0; i < val.size(); ++i)
          values[var].push_back(val[i].as<Real>());
      }
      else
        values[var].push_back(val.as<Real>());
    }
  }

  // Count up the distinct sets of values and set index counters.
  size_t num_cases = 1;
  std::map<std::string, int> indices;
  for (auto iter = values.begin(); iter != values.end(); ++iter)
  {
    indices[iter->first] = 0;
    num_cases *= iter->second.size();
  }

  // Shuffle the values into a vector of complete cases.
  std::vector<std::map<std::string, Real> > vec(num_cases);
  for (auto iter = values.begin(); iter != values.end(); ++iter)
  {
    std::string var = iter->first;
    for (size_t i = 0; i < iter->second.size(); ++i)
    {
      vec[indices[iter->first]][var] = iter->second[i];
      ++indices[iter->first];
    }
  }
  return vec;
}

} // anonymous namespace

std::vector<std::map<std::string, Real> > YamlFile::read_initial_conditions() const
{
  return read_vector_from_map(_root, "initial_conditions");
}

std::map<std::string, bool> YamlFile::read_physics_settings() const
{
  std::map<std::string, bool> settings;
  if (_root["physics"] and _root["physics"].IsMap())
  {
    auto node = _root["physics"];
    for (auto iter = node.begin(); iter != node.end(); ++iter)
      settings[iter->first.as<std::string>()] = iter->second.as<bool>();
  }
  return settings;
}

std::pair<std::string, std::map<std::string, Real> > YamlFile::read_atmosphere() const
{
  std::pair<std::string, std::map<std::string, Real> > atm;
  if (_root["atmosphere"] and _root["atmosphere"].IsMap())
  {
    auto node = _root["atmosphere"];
    if (node["model"])
      atm.first = node["model"].as<std::string>();
    else
      throw YamlException("No model found in atmosphere section.");
    for (auto iter = node.begin(); iter != node.end(); ++iter)
    {
      auto key = iter->first.as<std::string>();
      if (key != "model")
        atm.second[key] = iter->second.as<Real>();
    }
  }
}

std::vector<Real> YamlFile::read_timesteps() const
{
  std::vector<Real> dts;
  if (_root["timestep"])
  {
    YAML::Node timestep = _root["timestep"];
    if (timestep.IsSequence())
    {
      for (size_t i = 0; i < timestep.size(); ++i)
        dts.push_back(timestep[i].as<Real>());
    }
    else
      dts.push_back(timestep.as<Real>());
  }
  return dts;
}

Real YamlFile::read_duration() const
{
  if (not _root["duration"])
    throw YamlException("No duration found in YAML file.");
  else
    return _root["duration"].as<Real>();
}

std::pair<std::string, std::string> YamlFile::read_output_params() const
{
  auto prefix = std::string("mam_output");
  auto dir = std::string(".");
  if (_root["output"])
  {
    auto node = _root["output"];
    if (node["prefix"])
      prefix = node["prefix"].as<std::string>();
    if (node["directory"])
      dir = node["directory"].as<std::string>();
  }
  return std::make_pair(prefix, dir);
}

int YamlFile::read_output_freq() const
{
  int freq = -1;
  if (_root["output"])
  {
    auto node = _root["output"];
    if (node["frequency"])
      freq = node["frequency"].as<int>();
  }
  return freq;
}

} // namespace mam
