# This CMake function parses a C++ header containing physical constants for
# Haero and produces produces a C++ .cpp file that defines each constant, which
# is a silly requirement of pre-C++17 compilers. Gotta love this half-baked
# language.
#
# This brittle parser was copied from gen_consts_module.cmake, so we rely on
# some structure:
# 1. Constants are defined within a haero::constants C++ namespace
# 2. All constants are floating point numbers (of type Real or double).
function(gen_consts_cpp cxx_header cpp_file)
  set(in_haero_ns FALSE)
  set(in_constants_struct FALSE)

  # Generate C++ boilerplate.
  list(APPEND cxx_lines "// This file was generated from ${cxx_header}")
  list(APPEND cxx_lines "// using cmake/gen_consts_module.cmake. PLEASE DO NOT EDIT.")
  list(APPEND cxx_lines "#include \"haero/constants.hpp\"")
  list(APPEND cxx_lines "namespace haero {")

  # Parse the file into a list of strings containing lines (with newlines
  # removed).
  file(STRINGS ${cxx_header} lines)
  foreach(line ${lines})
    if (in_haero_ns)
      if (in_constants_struct)
        if (line MATCHES "};")
          break()
        endif()

        # Look for a Doxygen comment delimiter at the beginning of the line.
        string(FIND ${line} "///" pos)
        if (pos GREATER_EQUAL 0)
          math(EXPR begin "${pos} + 3")
          string(SUBSTRING ${line} ${begin} -1 comment)
        else()
          # This line isn't a Doxygen comment. Is it an assignment?
          string(FIND ${line} "=" equals)
          if (equals GREATER 0)
            # This is a constant. We want its type, its name, and its value.
            unset(type)
            if (line MATCHES "Real")
              set(type "Real")
            elseif (line MATCHES "double")
              set(type "double")
            endif()
            if (type)
              # Figure out the name in the assignment
              math(EXPR end "${equals} - 1")
              string(SUBSTRING ${line} 0 ${end} name)
              string(STRIP ${name} name)
              string(FIND ${name} " " last_space REVERSE)
              math(EXPR begin "${last_space} + 1")
              string(SUBSTRING ${name} ${begin} -1 name)
              string(STRIP ${name} name)

              list(APPEND cxx_lines "constexpr ${type} Constants::${name}_SEMICOLON_")
            endif()
          endif()
        endif()
      elseif(${line} MATCHES "struct Constants")
      set(in_constants_struct TRUE)
      endif()
    elseif(${line} MATCHES "namespace haero")
      set(in_haero_ns TRUE)
    endif()
  endforeach()
  list(APPEND cxx_lines "} // end namespace haero")

  # Write out the .cpp file.
  string(REPLACE ";" "\n" cxx_source "${cxx_lines}")
  string(REPLACE "_SEMICOLON_" "\;" cxx_source "${cxx_source}")
  file(WRITE ${cpp_file} ${cxx_source})
endfunction()

gen_consts_cpp(${CPP_HEADER} ${CPP_SOURCE})
