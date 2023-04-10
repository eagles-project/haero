include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Generate files containing our config parameters.
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/haero_config.hpp.in
  ${CMAKE_CURRENT_BINARY_DIR}/haero_config.hpp
  @ONLY
)

# Generate a blank haero_version.cpp file to get us off the ground.
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMAND ${CMAKE_COMMAND} -E touch ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMENT "Generating initial haero_version.cpp..."
)

# Update the haero version information if it's already been generated
add_custom_target(update_version_info ALL
  COMMAND ${BASH} ${PROJECT_SOURCE_DIR}/tools/update_version_info.sh ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
  COMMENT "Updating haero version information..."
)

# Generate a C++ file with physical constants defined, from our C++ header.
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/constants.cpp
  COMMAND ${CMAKE_COMMAND}
    -DCPP_HEADER=${CMAKE_CURRENT_SOURCE_DIR}/constants.hpp
    -DCPP_SOURCE=${CMAKE_CURRENT_BINARY_DIR}/constants.cpp
    -P ${PROJECT_SOURCE_DIR}/cmake/gen_consts_cpp.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/constants.hpp
          ${PROJECT_SOURCE_DIR}/cmake/gen_consts_cpp.cmake
  COMMENT "Generating constants.cpp..."
)

add_library(haero
            ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
            ${CMAKE_CURRENT_BINARY_DIR}/constants.cpp
            testing.cpp
            utils.cpp
            )
add_dependencies(haero ext_libraries update_version_info)
target_compile_options(haero PRIVATE  $<$<COMPILE_LANGUAGE:CXX>:${KOKKOS_CUDA_OPTIONS}>)
if(CMAKE_COMPILER_IS_GNUCXX)
  MESSAGE(STATUS "Using GNU compiler, turning on warnings are errors....")
  target_compile_options(haero PRIVATE -Werror)
endif()
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES})
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES} PARENT_SCOPE)

# unit tests
add_subdirectory(tests)

# Installation targets
install(TARGETS haero DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES aero_process.hpp
              aero_species.hpp
              atmosphere.hpp
              constants.hpp
              floating_point.hpp
              gas_species.hpp
              haero.hpp
              ${CMAKE_CURRENT_BINARY_DIR}/haero_config.hpp
              math.hpp
              testing.hpp
              utils.hpp
              root_finders.hpp
        DESTINATION include/haero)

