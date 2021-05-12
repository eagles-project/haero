#ifndef HAERO_SKYWALKER_WRITE_PY_MODULE_HPP
#define HAERO_SKYWALKER_WRITE_PY_MODULE_HPP

#include "skywalker.hpp"

namespace skywalker {

// Writes simulation output data to a Python module.
void write_py_module(const std::vector<InputData>& input,
                     const std::vector<OutputData>& output,
                     const char* py_module_name);

}

#endif
