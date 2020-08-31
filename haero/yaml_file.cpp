#include "mam_yaml_file.hpp"

namespace mam {

Yaml_file::Yaml_file(const std::string& filename)
{
  // Try to load the thing.
  try
  {
    _root = YAML::LoadFile(filename);
  }
  catch (YAML::BadFile& e)
  {
    throw Yaml_exception(e.what());
  }
  catch (YAML::ParserException& e)
  {
    throw Yaml_exception(e.what());
  }
}

Yaml_file::~Yaml_file() {
}

std::vector<Mode> Yaml_file::read_modes() const
{
  std::vector<Mode> modes;
  if (_root["modes"] and _root["modes"].IsSequence())
  {
    YAML::Node node = _root["modes"];
    for (size_t i = 0; i < node.size(); ++i)
    {
      if (not node[i]["name"])
        throw Yaml_exception("mode entry has no name.");
      else if (not node[i]["D_min"])
        throw Yaml_exception("mode entry has no minimum diameter (D_min).");
      else if (not node[i]["D_max"])
        throw Yaml_exception("mode entry has no maximum diameter (D_max).");
      else if (not node[i]["sigma"])
        throw Yaml_exception("mode entry has no geometric stddev (sigma).");
      else
      {
        modes.push_back(Mode(node[i]["name"].as<std::string>(),
                             node[i]["D_min"].as<Real>(),
                             node[i]["D_max"].as<Real>(),
                             node[i]["sigma"].as<Real>()));
      }
    }
  }
  return modes;
}

std::vector<Species> Yaml_file::read_species() const
{
  std::vector<Species> species;
  return species;
}

std::vector<Reaction> Yaml_file::read_reactions() const
{
  std::vector<Reaction> reactions;
  return reactions;
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

std::vector<std::map<std::string, Real> > Yaml_file::read_initial_conditions() const
{
  return read_vector_from_map(_root, "initial_conditions");
}

Real Yaml_file::read_perturb_factor() const
{
  Real pert_factor = 0.0;
  if (_root["perturbation"])
    pert_factor = _root["perturbation"].as<Real>();
  return pert_factor;
}

void Yaml_file::read_control_params(int& do_gaschem, int& do_cloudchem,
                                    int& do_gasaerexch, int& do_rename,
                                    int& do_newnuc, int& do_coag,
                                    int& do_calcsize, int& num_unit,
                                    int& frac_unit, int& gas_unit)
{
  do_gaschem = 0;
  do_cloudchem = 0;
  do_gasaerexch = 0;
  do_rename = 0;
  do_newnuc = 0;
  do_coag = 0;
  do_calcsize = 0;
  num_unit = 0;
  frac_unit = 0;
  gas_unit = 0;
  if (_root["parameters"] and _root["parameters"].IsMap())
  {
    auto node = _root["parameters"];
    if (node["do_gaschem"])
      do_gaschem = node["do_gaschem"].as<int>();
    if (node["do_cloudchem"])
      do_cloudchem = node["do_cloudchem"].as<int>();
    if (node["do_gasaerexch"])
      do_gasaerexch = node["do_gasaerexch"].as<int>();
    if (node["do_rename"])
      do_rename = node["do_rename"].as<int>();
    if (node["do_newnuc"])
      do_newnuc = node["do_newnuc"].as<int>();
    if (node["do_coag"])
      do_coag = node["do_coag"].as<int>();
    if (node["do_calcsize"])
      do_coag = node["do_calcsize"].as<int>();
    if (node["num_unit"])
      num_unit = node["num_unit"].as<int>();
    if (node["frac_unit"])
      frac_unit = node["frac_unit"].as<int>();
    if (node["gas_unit"])
      gas_unit = node["gas_unit"].as<int>();
  }
}

void Yaml_file::read_met_input(Real& temp, Real& press, Real& RH_CLEA,
                               Real& hgt, Real& cld_frac)
{
  temp = 300.0;
  press = 3.0e4;
  RH_CLEA = 0.5;
  hgt = 7500.0;
  cld_frac = 0.0;
  if (_root["met_input"] and _root["met_input"].IsMap())
  {
    auto node = _root["met_input"];
    if (node["temp"])
      temp = node["temp"].as<Real>();
    if (node["press"])
      press = node["press"].as<Real>();
    if (node["RH_CLEA"])
      RH_CLEA = node["RH_CLEA"].as<Real>();
    if (node["hgt"])
      hgt = node["hgt"].as<Real>();
    if (node["cld_frac"])
      cld_frac = node["cld_frac"].as<Real>();
  }
}

std::vector<Real> Yaml_file::read_timesteps() const
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

Real Yaml_file::read_duration() const
{
  if (not _root["duration"])
    throw Yaml_exception("No duration found in YAML file.");
  else
  {
    return _root["duration"].as<Real>();
  }
}

void Yaml_file::read_output_params(std::string& dir,
                                   std::string& prefix,
                                   int& freq) const
{
  prefix = std::string("mam_output");
  dir = std::string(".");
  freq = -1;
  if (_root["output"])
  {
    auto node = _root["output"];
    if (node["prefix"])
      prefix = node["prefix"].as<std::string>();
    if (node["directory"])
      dir = node["directory"].as<std::string>();
    if (node["frequency"])
      freq = node["frequency"].as<int>();
  }
}

} // namespace mam
