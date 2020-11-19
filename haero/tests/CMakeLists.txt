include(EkatCreateUnitTest)

# This simple unit test gets us off the ground.
EkatCreateUnitTest(hello_world hello_world.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(view_pack_helpers view_helper_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(floating_point_tests floating_point_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(mode_tests mode_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(species_tests species_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(model_tests model_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(prognostics_tests prognostics_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(tendencies_tests tendencies_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

add_library(haero_process_stubs
            prog_process_stub.cpp prog_process_stub.f90
            diag_process_stub.cpp diag_process_stub.f90)

EkatCreateUnitTest(process_stub_tests process_stub_tests.cpp
                   LIBS haero_process_stubs;${HAERO_LIBRARIES})
