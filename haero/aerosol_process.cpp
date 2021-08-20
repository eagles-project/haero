#include "haero/aerosol_process.hpp"
#include "haero/utils.hpp"

namespace haero {

void AerosolProcess::interpret_and_set_param(const std::string& name,
                                             const std::string& value) {
  // Is this an integer?
  try {
    int int_value = std::stoi(value);
    set_param(name, int_value);
  }
  catch (const std::invalid_argument&) {
    // Okay, it's not an integer. Is it a real number?
    try {
#if HAERO_DOUBLE_PRECISION
      Real real_value = std::stod(value);
#else
      Real real_value = std::stof(value);
#endif
      set_param(name, real_value);
    }
    catch (const std::invalid_argument&) {
      // Boolean?
      if (is_boolean(value)) {
        set_param(name, as_boolean(value));
      } else {
        // Okay, we can only interpret this value as a string. String parameters
        // aren't supported for Fortran processes, so make sure we're not one of
        // those.
#if HAERO_FORTRAN
        if (dynamic_cast<haero::FAerosolProcess*>(this) != nullptr) {
          fprintf(stderr, "Parameter '%s' with string value '%s' given for\n"
                  "Fortran aerosol process %s. Fortran aerosol processes cannot\n"
                  "accept string parameter values.", name.c_str(), value.c_str(),
                  name_.label().c_str());
          return;
        }
#endif
        set_param(name, value);
      }
    }
  }
}

} // haero
