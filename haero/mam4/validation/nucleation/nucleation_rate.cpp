#include <haero/mam4/merikanto2007.hpp>
#include <haero/mam4/vehkamaki2002.hpp>
#include <haero/mam4/wang2008.hpp>
#include <iostream>
#include <skywalker.hpp>
#include <validation.hpp>

using namespace skywalker;
using namespace haero;
using namespace haero::mam4;

namespace {

void run_vehkamaki2002(Ensemble* ensemble, int pbl_method) {
  ensemble->process([](const Input& input, Output& output) {
    // Ensemble parameters
    PackType c_h2so4 = input.get("c_h2so4");
    PackType rel_hum = input.get("relative_humidity");
    PackType temp = input.get("temperature");

    // Compute the mole fraction of H2SO4 in a critical cluster, and from it
    // the nucleation rate.
    PackType x_crit =
        vehkamaki2002::h2so4_critical_mole_fraction(c_h2so4, temp, rel_hum);
    PackType J = vehkamaki2002::nucleation_rate(c_h2so4, temp, rel_hum, x_crit);

    // Write the computed nucleation rate.
    output.set("nucleation_rate", J[0]);
  });
}

void run_merikanto2007(Ensemble* ensemble, int pbl_method) {
  ensemble->process([](const Input& input, Output& output) {
    // Ensemble parameters
    PackType c_h2so4 = input.get("c_h2so4");
    PackType xi_nh3 = input.get("xi_nh3");
    PackType rel_hum = input.get("relative_humidity");
    PackType temp = input.get("temperature");

    // Compute the nucleation rate.
    PackType log_J =
        merikanto2007::log_nucleation_rate(temp, rel_hum, c_h2so4, xi_nh3);
    PackType J = exp(log_J);

    // Write the computed nucleation rate.
    output.set("nucleation_rate", J[0]);
  });
}

} // anonymous namespace

void nucleation_rate(Ensemble* ensemble) {
  // Figure out settings for binary/ternary nucleation and planetary boundary
  // layer treatment
  Settings settings = ensemble->settings();
  int nuc_method = std::stoi(settings.get("nucleation_method"));
  if ((nuc_method != 2) and (nuc_method != 3)) {
    std::stringstream ss;
    ss << "Invalid nucleation method: " << nuc_method << std::endl;
    throw skywalker::Exception(ss.str());
  }

  int pbl_method = 0;  // no PBL correction by default.
  if (settings.has("pbl_method")) {
    pbl_method = std::stoi(settings.get("pbl_method"));
  }
  if ((pbl_method < 0) or (pbl_method > 2)) {
    std::stringstream ss;
    ss << "Invalid planetary boundary layer method: " << pbl_method
       << std::endl;
    throw skywalker::Exception(ss.str());
  }

  // Run the ensemble.
  if (nuc_method == 2) {  // binary nucleation
    run_vehkamaki2002(ensemble, pbl_method);
  } else {  // ternary nucleation
    run_merikanto2007(ensemble, pbl_method);
  }
}
