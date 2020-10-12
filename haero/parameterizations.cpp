#include "haero/aero_process.hpp"
#include "haero/parameterizations.hpp"

namespace haero {

AeroProcess* select_parametrization(AeroProcessType type,
                                    const Parameterizations& selections)
{
  AeroProcess* process = nullptr;
  if (type == ActivationProcess) {
    if (selections.activation == Parameterizations::NoActivation) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == CloudBorneWetRemovalProcess) {
    if (selections.cloudborne_wet_removal == Parameterizations::NoCloudBorneWetRemoval) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == CoagulationProcess) {
    if (selections.coagulation == Parameterizations::NoCoagulation) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == CondensationProcess) {
    if (selections.condensation == Parameterizations::NoCondensation) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == DryDepositionProcess) {
    if (selections.cloudborne_wet_removal == Parameterizations::NoCloudBorneWetRemoval) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == EmissionsProcess) {
    if (selections.emissions == Parameterizations::NoEmissions) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == InterstitialWetRemovalProcess) {
    if (selections.interstitial_wet_removal == Parameterizations::NoInterstitialWetRemoval) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == NucleationProcess) {
    if (selections.nucleation == Parameterizations::NoNucleation) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == ResuspensionProcess) {
    if (selections.resuspension == Parameterizations::NoResuspension) {
      process = new NullAeroProcess(type);
    }
  }
  else if (type == WaterUptakeProcess) {
    if (selections.water_uptake == Parameterizations::NoWaterUptake) {
      process = new NullAeroProcess(type);
    }
  }
  EKAT_ASSERT_MSG(process != nullptr, "No process was selected!");
  return process;
}

}

