# This CMake function parses a C++ header containing physical constants for
# Haero and produces a Fortran 90 module containing those same constants.
# This is a pretty brittle parser, so we rely on some structure:
# 1. Constants are defined within a haero::constants C++ namespace
# 2. Every constant can be preceded by a single-line comment, to be converted
#    to a Fortran equivalent.
# 3. All constants are floating point numbers (of type Real or double).
function(gen_consts_module cxx_header fortran_module)
  set(in_haero_ns FALSE)
  set(in_constants_ns FALSE)

  # Generate Fortran boilerplate.
  list(APPEND fortran_lines "! This file was generated from ${cxx_header}")
  list(APPEND fortran_lines "! using cmake/gen_consts_module.cmake. PLEASE DO NOT EDIT.")
  list(APPEND fortran_lines "module haero_constants")
  list(APPEND fortran_lines "  use haero, only: wp")
  list(APPEND fortran_lines "  use iso_c_binding, only: c_double")
  list(APPEND fortran_lines "  implicit none")
  list(APPEND fortran_lines "  public")

  # Parse the file into a list of strings containing lines (with newlines
  # removed).
  file(STRINGS ${cxx_header} lines)
  foreach(line ${lines})
    if(in_haero_ns)
      if(in_constants_ns)
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
              set(type "real(wp)")
            elseif (line MATCHES "double")
              set(type "real(c_double)")
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

              # Figure out the value
              math(EXPR begin "${equals} + 1")
              string(SUBSTRING ${line} ${begin} -1 value)
              string(REPLACE ";" "" value ${value})
              string(STRIP ${value} value)

              # Now we write it to our Fortran code list.
              if (comment)
                list(APPEND fortran_lines "  !${comment}")
                unset(comment)
              endif()
              list(APPEND fortran_lines "  ${type}, parameter :: ${name} = ${value}")
            endif()
          endif()
        endif()
      elseif(${line} MATCHES "namespace constants")
      set(in_constants_ns TRUE)
      endif()
    elseif(${line} MATCHES "namespace haero")
      set(in_haero_ns TRUE)
    endif()
  endforeach()
  list(APPEND fortran_lines "end module")
  string(REPLACE ";" "\n" fortran_source "${fortran_lines}")
  file(WRITE ${fortran_module} ${fortran_source})
endfunction()

gen_consts_module(${CPP_HEADER} ${F90_MODULE})
