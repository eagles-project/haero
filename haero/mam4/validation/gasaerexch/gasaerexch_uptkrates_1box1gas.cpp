#include <haero/mam4.hpp>
#include <iostream>
#include <skywalker.hpp>
#include <vector>

using namespace haero;
using namespace skywalker;

void test_gasaerexch_uptkrates_1box1gas_process(
    const Input& input, Output& output, mam4::GasAerExchImpl& gasaerexch) {
  // Ensemble parameters
  if (!input.has("temp")) {
    std::cerr << "Required name: "
              << "temp" << std::endl;
    exit(0);
  }
  if (!input.has_array("dgncur_awet")) {
    std::cerr << "Required name: "
              << "dgncur_awet" << std::endl;
    exit(0);
  }
  if (!input.has_array("lnsg")) {
    std::cerr << "Required name: "
              << "lnsg" << std::endl;
    exit(0);
  }
  if (!input.has_array("aernum")) {
    std::cerr << "Required name: "
              << "aernum" << std::endl;
    exit(0);
  }
  const bool has_solution = input.has_array("uptkaer");

  //-------------------------------------------------------
  // Process input, do calculations, and prepare output
  //-------------------------------------------------------
  const int n_mode = 4;
  const int nghq = 2;
  const Real accom = 0.65000000000000002;
  const Real beta_inp = 0.0000000000000000;
  const Real pi = 3.1415926535897931;
  const Real r_universal = 8314.4675910000005;
  const Real mw_gas = 98.078400000000002;
  const Real mw_air = 28.966000000000001;
  const Real pmid = 100000.00000000000;
  const Real pstd = 101325.00000000000;
  const Real vol_molar_gas = 42.880000000000003;
  const Real vol_molar_air = 20.100000000000001;

  const Kokkos::Array<bool, n_mode> l_condense_to_mode = {true, true, true,
                                                          true};

  // Parse input
  Kokkos::Array<PackType, n_mode> dgncur_awet;
  Kokkos::Array<Real, n_mode> lnsg;
  {
    const std::vector<Real> array = input.get_array("dgncur_awet");
    for (size_t i = 0; i < array.size() && i < n_mode; ++i)
      dgncur_awet[i] = array[i];
  }
  {
    const std::vector<Real> array = input.get_array("lnsg");
    for (size_t i = 0; i < array.size() && i < n_mode; ++i) lnsg[i] = array[i];
  }
  const std::vector<Real> aernum = input.get_array("aernum");
  const Real temp = input.get("temp");
  std::vector<Real> test_uptkaer;
  if (has_solution) {
    test_uptkaer = input.get_array("uptkaer");
  }

  kokkos_device_type::view_1d<PackType> uptkaer_dev("uptkaer on device",
                                                    n_mode);
  Kokkos::parallel_for(
      "gasaerexch.gas_aer_uptkrates_1box1gas", 1, KOKKOS_LAMBDA(const int) {
        Kokkos::Array<PackType, n_mode> uptkaer;
        gasaerexch.gas_aer_uptkrates_1box1gas(
            l_condense_to_mode, temp, pmid, pstd, mw_gas, mw_air, vol_molar_gas,
            vol_molar_air, accom, r_universal, pi, beta_inp, nghq, dgncur_awet,
            lnsg, uptkaer);
        for (size_t i = 0; i < n_mode; ++i) uptkaer_dev(i) = uptkaer[i];
      });
  Kokkos::Array<PackType, n_mode> uptkaer;
  {
    auto host_view = Kokkos::create_mirror_view(uptkaer_dev);
    Kokkos::deep_copy(host_view, uptkaer_dev);
    for (size_t i = 0; i < n_mode; ++i) uptkaer[i] = host_view[i];
  }
  // Write the computed nucleation rate.
  {
    std::vector<Real> values(n_mode);
    for (size_t i = 0; i < values.size() && i < n_mode; ++i)
      values[i] = uptkaer[i][0];
    output.set("uptkaer", values);
  }
}

void test_gasaerexch_uptkrates_1box1gas(std::unique_ptr<Ensemble>& ensemble) {
  mam4::AeroConfig mam4_config;
  mam4::GasAerExchImpl gasaerexch;
  gasaerexch.init(mam4_config);

  ensemble->process([&](const Input& input, Output& output) {
    test_gasaerexch_uptkrates_1box1gas_process(input, output, gasaerexch);
  });
}
