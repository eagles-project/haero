include(create_unit_test)
#link_directories(${EKAT_LIBRARY_DIR}) # Needed to link unit tests.

CreateUnitTest(mam4_nucleation_unit_tests mam4_nucleation_unit_tests.cpp)
CreateSkywalkerTest(mam4_gasaerexch_unit_tests mam4_gasaerexch_unit_tests.cpp)

configure_file(${PROJECT_SOURCE_DIR}/haero/mam4/validation/mam_x_validation/gasaerexch/skywalker_gasaerexch_uptkrates_1box1gas.yaml 
    skywalker_gasaerexch_uptkrates_1box1gas.yaml COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/haero/mam4/validation/mam_x_validation/gasaerexch/skywalker_gasaerexch_uptkrates_1box1gas.py 
    skywalker_gasaerexch_uptkrates_1box1gas.py.ref COPYONLY)

add_test(NAME gasaerexch_compare_output
    COMMAND ${CMAKE_COMMAND} -E compare_files skywalker_gasaerexch_uptkrates_1box1gas.py skywalker_gasaerexch_uptkrates_1box1gas.py.ref
)

set_tests_properties(gasaerexch_compare_output  PROPERTIES DEPENDS mam4_gasaerexch_unit_tests)

