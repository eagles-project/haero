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

if (HAERO_FORTRAN)
  # Generate a Fortran module with physical constants, from our C++ header.
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/haero_constants.F90
    COMMAND ${CMAKE_COMMAND}
      -DCPP_HEADER=${CMAKE_CURRENT_SOURCE_DIR}/physical_constants.hpp
      -DF90_MODULE=${CMAKE_CURRENT_BINARY_DIR}/haero_constants.F90
      -P ${PROJECT_SOURCE_DIR}/cmake/gen_consts_module.cmake
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/physical_constants.hpp
            ${PROJECT_SOURCE_DIR}/cmake/gen_consts_module.cmake
    COMMENT "Generating haero_constants.F90..."
  )

  set(HAERO_FORTRAN_SOURCES
      haero.F90
      haero_precision.F90
      ${CMAKE_CURRENT_BINARY_DIR}/haero_constants.F90)
endif()


# Gather Haero process implementations (placed in HAERO_PROCESSES).
add_subdirectory(processes)
# Generate a header file that brings in all C++ process headers.
set(available_processes_hpp ${CMAKE_CURRENT_BINARY_DIR}/available_processes.hpp)
file(WRITE ${available_processes_hpp} "#ifndef HAERO_AVAILABLE_PROCESSES_HPP\n#define HAERO_AVAILABLE_PROCESSES_HPP\n\n")
foreach(src ${HAERO_PROCESSES})
  if (${src} MATCHES ".cpp")
    message(STATUS "Adding Haero Process:" ${src})
    get_filename_component(filename ${src} NAME)
    string(REPLACE ".cpp" ".hpp" header ${filename})
    if (${header} MATCHES "fprocess")
      file(APPEND ${available_processes_hpp} "#if HAERO_FORTRAN\n")
      file(APPEND ${available_processes_hpp} "#include \"haero/processes/${header}\"\n")
      file(APPEND ${available_processes_hpp} "#endif\n")
    else()
      file(APPEND ${available_processes_hpp} "#include \"haero/processes/${header}\"\n")
    endif()
  endif()
endforeach()
file(APPEND ${available_processes_hpp} "#endif\n")

# Gather Haero process implementations (placed in HAERO_DIAGNOSTICS).
add_subdirectory(diagnostics)
# Generate a header file that brings in all C++ diagnostic headers.
set(available_diagnostics_hpp ${CMAKE_CURRENT_BINARY_DIR}/available_diagnostics.hpp)
file(WRITE ${available_diagnostics_hpp} "#ifndef HAERO_AVAILABLE_DIAGNOSTICS_HPP\n#define HAERO_AVAILABLE_DIAGNOSTICS_HPP\n\n")
foreach(src ${HAERO_DIAGNOSTICS})
  if (${src} MATCHES ".cpp")
    message(STATUS "Adding Haero diagnostic: " ${src})
    get_filename_component(filename ${src} NAME)
    string(REPLACE ".cpp" ".hpp" header ${filename})
    if (${header} MATCHES "fdiagnostic")
      file (APPEND ${available_diagnostics_hpp} "#if HAERO_FORTRAN\n")
      file (APPEND ${available_diagnostics_hpp} "#include \"haero/diagnostics/${header}\"\n")
      file (APPEND ${available_diagnostics_hpp} "#endif\n")
    else()
      file (APPEND ${available_diagnostics_hpp} "#include \"haero/diagnostics/${header}\"\n")
    endif()
  endif()
endforeach()
file(APPEND ${available_diagnostics_hpp} "#endif\n")

add_library(haero
            ${CMAKE_CURRENT_BINARY_DIR}/haero_version.cpp
            ${HAERO_FORTRAN_SOURCES}
            atmosphere.cpp
            diagnostics.cpp
            modal_aerosol_config.cpp
            model.cpp
            prognostics.cpp
            selected_processes.cpp
            tendencies.cpp
            utils.cpp
            view_pack_helpers.cpp
            ${HAERO_PROCESSES} # populated by the processes/ ѕubdirectory
            ${HAERO_DIAGNOSTICS} # populated by the diagnostics/ subdirectory
            )
add_dependencies(haero ekat update_version_info tchem)
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES})
set(HAERO_LIBRARIES haero;${HAERO_LIBRARIES} PARENT_SCOPE)

# Installation targets
install(TARGETS haero DESTINATION lib)
install(FILES haero.hpp
              ${CMAKE_CURRENT_BINARY_DIR}/haero_config.hpp
              diagnostics.hpp
              floating_point.hpp
              mode.hpp
              model.hpp
              selected_processes.hpp
              physical_constants.hpp
              prognostics.hpp
              aerosol_species.hpp
              utils.hpp
              view_pack_helpers.hpp
        DESTINATION include/haero)

# Library unit tests.
add_subdirectory(tests)

