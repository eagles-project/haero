include(EkatCreateUnitTest)

# This simple unit test gets us off the ground.
EkatCreateUnitTest(hello_world hello_world.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(view_pack_helpers view_helper_tests.cpp
                   LIBS ${HAERO_LIBRARIES})

EkatCreateUnitTest(floating_point_tests floating_point_tests.cpp LIBS ${HAERO_LIBRARIES})
