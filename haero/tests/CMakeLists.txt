include(create_unit_test)
link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

# This simple unit test gets us off the ground.
CreateUnitTest(hello_world hello_world.cpp)

#----------------------------------------------------
# Unit tests for basic data structures and functions
#----------------------------------------------------

CreateUnitTest(view_pack_helpers view_helper_tests.cpp)
CreateUnitTest(floating_point_tests floating_point_tests.cpp)
CreateUnitTest(conversions_tests conversions_tests.cpp)
CreateUnitTest(mode_tests mode_tests.cpp)
CreateUnitTest(aerosol_species_tests aerosol_species_tests.cpp)
CreateUnitTest(gas_species_tests gas_species_tests.cpp)
CreateUnitTest(model_tests model_tests.cpp)
CreateUnitTest(prognostics_tests prognostics_tests.cpp)
CreateUnitTest(tendencies_tests tendencies_tests.cpp)
CreateUnitTest(region_of_validity_tests region_of_validity_tests.cpp)
CreateUnitTest(kokkos_team_execute kokkos_team_execute.cpp)
CreateUnitTest(process_tests process_tests.cpp)
CreateUnitTest(cuda_tests cuda_tests.cpp)
CreateUnitTest(math_tests math_tests.cpp)
CreateUnitTest(kohler_tests kohler_tests.cpp)
CreateUnitTest(mode_wet_radius mode_wet_radius_tests.cpp)
CreateUnitTest(microphysics_functions microphys_tests.cpp)

#----------------------------------------------------
# Unit tests for low-level aerosol parameterizations
#----------------------------------------------------

CreateUnitTest(vehkamaki2002_tests vehkamaki2002_tests.cpp)
CreateUnitTest(merikanto2007_tests merikanto2007_tests.cpp)
CreateUnitTest(wang2008_tests wang2008_tests.cpp)

#----------------------------------------------------
# Unit tests for aerosol processes
#----------------------------------------------------

CreateUnitTest(mam_nucleation_process_tests mam_nucleation_process_tests.cpp)
CreateUnitTest(mam_rename_process_tests mam_rename_process_tests.cpp)
CreateUnitTest(mam_nucleation_packed_tests mam_nucleation_packed_tests.cpp)
CreateUnitTest(aerosol_activation_tests aerosol_activation_tests.cpp)

if (HAERO_FORTRAN)
  # Fortran process stub tests -- checks that data can be passed between C++ and
  # Fortran
  add_library(haero_faerosol_process_stubs
              faerosol_process_stub.cpp faerosol_process_stub.F90
              mam_gasaerexch_test_bridge.F90
              mam_nucleation_test_bridge.F90)
  add_dependencies(haero_faerosol_process_stubs haero)
  EkatCreateUnitTest(faerosol_process_stub_tests faerosol_process_stub_tests.cpp
                     LIBS haero_faerosol_process_stubs;${HAERO_LIBRARIES})

  # TODO: the Fortran version of the nucleation process doesn't currently pass
  # TODO: tests in single-precision builds. We'll punt for now, but should
  # TODO: probably try to figure this out soon.
  if (HAERO_DOUBLE_PRECISION)
    # Aerosol process unit tests
    EkatCreateUnitTest(mam_gasaerexch_fprocess_tests mam_gasaerexch_fprocess_tests.cpp
                       LIBS haero_faerosol_process_stubs;${HAERO_LIBRARIES})

    EkatCreateUnitTest(mam_nucleation_fprocess_tests mam_nucleation_fprocess_tests.cpp
                       LIBS haero_faerosol_process_stubs;${HAERO_LIBRARIES})

    EkatCreateUnitTest(mam_calcsize_fprocess_tests mam_calcsize_fprocess_tests.cpp
                       LIBS haero_faerosol_process_stubs;${HAERO_LIBRARIES})

    EkatCreateUnitTest(mam_rename_fprocess_tests mam_rename_fprocess_tests.cpp
                       LIBS haero_faerosol_process_stubs;${HAERO_LIBRARIES})
  endif()
endif()
