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
    -P ${PROJECT_SOURCE_DIR}/cmake/HaeroGenerateConstants.cmake
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/constants.hpp
          ${PROJECT_SOURCE_DIR}/cmake/HaeroGenerateConstants.cmake
  COMMENT "Generating constants.cpp..."
)

add_library(haero
            ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
            ${CMAKE_CURRENT_BINARY_DIR}/constants.cpp
            testing.cpp
            utils.cpp
            )
add_dependencies(haero update_version_info)
target_include_directories(haero PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(haero PUBLIC ekat::AllLibs)

if (HAERO_ENABLE_GPU)
  if (HAERO_AMD_GPU)
    # Nothing... yet
  else()
    target_compile_options(haero PRIVATE  $<$<COMPILE_LANGUAGE:CXX>:${KOKKOS_CUDA_OPTIONS}>)
  endif()
  if(CMAKE_COMPILER_IS_GNUCXX AND HAERO_STANDALONE)
    MESSAGE(STATUS "Using GNU compiler, turning on warnings are errors....")
    target_compile_options(haero PRIVATE -Werror)
  endif()
endif()
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES})
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES} PARENT_SCOPE)

# unit tests
if (HAERO_ENABLE_TESTS)
  add_subdirectory(tests)
endif()

# Installation targets
install(TARGETS haero
        EXPORT HaeroTargets
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES aero_process.hpp
              aero_species.hpp
              atmosphere.hpp
              surface.hpp
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
