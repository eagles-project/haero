# Generates haero.cmake (for downstream aerosol packages) and install it.
macro(HaeroGenerateConfig)
  set(HAERO_TPL_IMPORTED_LOCATIONS "") # <-- generated segment of haero.cmake
  foreach (lib kokkoscore kokkoscontainers kokkossimd spdlog yaml-cpp ekat ekat_test_main ekat_test_session haero)
    string(APPEND HAERO_TPL_IMPORTED_LOCATIONS "add_library(${lib} STATIC IMPORTED GLOBAL)\n")
    set(suffix "")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      # It seems to be a trend these days to append 'd' to your library in
      # Debug builds. :-/
      if (${lib} MATCHES "spdlog" OR ${lib} MATCHES "yaml-cpp")
        set(suffix "d")
      endif()
    endif()
    get_target_property(imported_loc ${lib} IMPORTED_LOCATION)
    if (imported_loc MATCHES "-NOTFOUND")
      set(imported_loc ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/lib${lib}${suffix}.a)
    endif()
    string(APPEND HAERO_TPL_IMPORTED_LOCATIONS "set_target_properties(${lib} PROPERTIES IMPORTED_LOCATION ${imported_loc})\n")
  endforeach()
  # generate the cmake file used by mam4xx
  configure_file(
    ${PROJECT_SOURCE_DIR}/cmake/haero.cmake.in
    ${PROJECT_BINARY_DIR}/share/haero.cmake
    @ONLY
  )
  install(DIRECTORY ${PROJECT_BINARY_DIR}/share DESTINATION .)
endmacro()
