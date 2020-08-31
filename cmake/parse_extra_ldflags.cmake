# This macro parses the given string of extra LDFLAGS into additional
# libraries and library paths.
macro(parse_extra_ldflags extra_ldflags)
  # Split the given string into a list.
  string(REPLACE " " ";" extra_ldflags_list ${extra_ldflags})
  foreach(item ${extra_ldflags_list})
    if (item MATCHES "-l")
      string(REPLACE "-l" "" library ${item})
      list(APPEND HAERO_LIBRARIES ${library})
    elseif(item MATCHES "-L")
      string(REPLACE "-L" "" libdir ${item})
      link_directories(${libdir})
    endif()
  endforeach()
endmacro()
