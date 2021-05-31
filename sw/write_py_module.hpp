#ifndef HAERO_SKYWALKER_WRITE_PY_MODULE_HPP
#define HAERO_SKYWALKER_WRITE_PY_MODULE_HPP

#include "skywalker.hpp"

namespace skywalker {

// Writes simulation output data to a Python module.
void write_py_module(const std::vector<InputData>& inputs,
                     const std::vector<OutputData>& outputs,
                     const char* py_module_name);

}  // namespace skywalker

#endif
