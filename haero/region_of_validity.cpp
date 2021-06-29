#include "haero/region_of_validity.hpp"

#include "numeric"

namespace haero {

RegionOfValidity::RegionOfValidity()
    : temp_bounds({0, 500}), rel_hum_bounds({0, 1}) {}

RegionOfValidity::~RegionOfValidity() {}

HostRegionOfValidity::HostRegionOfValidity() : RegionOfValidity() {}

HostRegionOfValidity::~HostRegionOfValidity() {}

const RegionOfValidity& HostRegionOfValidity::getRegionOfValidity() const {
  return *this;
}

RegionOfValidity::Token HostRegionOfValidity::find_gas_bounds(
    const std::string& symbol) const {
  return get_string_to_token_gases(symbol);
}

RegionOfValidity::Token HostRegionOfValidity::add_gas_bounds(
    const std::string& symbol) {
  Token return_val = find_gas_bounds(symbol);
  EKAT_REQUIRE_MSG(BOUNDS_NOT_FOUND == return_val,
                   "Bounds for gas species " << symbol << " already exist!");
  return_val = gas_bounds_.extent(0);
  set_string_to_token_gases(symbol, return_val);
  return return_val;
}

RegionOfValidity::Token HostRegionOfValidity::set_string_to_token(
    std::map<std::string, Token>& registered_strings, const std::string& name,
    const Token token) {
  Token return_val = get_string_to_token(registered_strings, name);

  if (BOUNDS_NOT_FOUND != token && BOUNDS_NOT_FOUND == return_val) {
    return_val = registered_strings.size();
    const std::pair<std::string, Token> value(name, return_val);
    registered_strings.insert(value);
  }
  return return_val;
}

RegionOfValidity::Token HostRegionOfValidity::get_string_to_token(
    const std::map<std::string, Token>& registered_strings,
    const std::string& name) {
  Token return_val = BOUNDS_NOT_FOUND;
  const auto iter = registered_strings.find(name);
  if (registered_strings.end() != iter) return_val = iter->second;
  return return_val;
}

RegionOfValidity::Token HostRegionOfValidity::set_string_to_token_gases(
    const std::string& symbol, const Token token) {
  return set_string_to_token(registered_gases_, symbol, token);
}

RegionOfValidity::Token HostRegionOfValidity::get_string_to_token_gases(
    const std::string& name) const {
  return get_string_to_token(registered_gases_, name);
}

}  // namespace haero
