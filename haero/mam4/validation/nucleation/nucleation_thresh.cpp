#include <haero/mam4/vehkamaki2002.hpp>

#include <iostream>
#include <skywalker.hpp>
#include <validation.hpp>

using namespace skywalker;
using namespace haero;
using namespace haero::mam4;

void nucleation_thresh(Ensemble* ensemble) {
  Settings settings = ensemble->settings();
  if (settings.has("nucleation_method")) {
    int nuc_method = std::stoi(settings.get("nucleation_method"));
    if ((nuc_method != 2)) {
      std::stringstream ss;
      ss << "Invalid nucleation method: " << nuc_method << std::endl;
      throw skywalker::Exception(ss.str());
    }
  }

  // Run the ensemble.
  ensemble->process([](const Input& input, Output& output) {
    // Ensemble parameters
    PackType rel_hum = input.get("relative_humidity");
    PackType temp = input.get("temperature");

    // Compute the threshold of H2SO4 above which nucleation occurs.
    PackType c_thresh =
        vehkamaki2002::h2so4_nucleation_threshold(temp, rel_hum);

    // Write the computed nucleation rate.
    output.set("nucleation_threshold", c_thresh[0]);
  });
}
