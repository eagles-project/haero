#include "write_py_module.hpp"

namespace {

haero::Real fetch_input_var(const skywalker::InputData& input,
                            const std::string& name)
{
  return 0.0;
}

// Writes the given input variable to our Python module.
void write_input_var(FILE* file,
                     const std::vector<skywalker::InputData>& inputs,
                     const std::string& var_name) {
  fprintf(file, "%s = [", var_name.c_str());
  for (auto input: inputs) {
    auto var = fetch_input_var(input, var_name);
    fprintf(file, "%g, ", var);
  }
  fprintf(file, "]\n");
}

haero::Real fetch_output_var(const skywalker::OutputData& output,
                             const std::string& name)
{
  return 0.0;
}

// Writes the given output variable to our Python module.
void write_output_var(FILE* file,
                      const std::vector<skywalker::OutputData>& outputs,
                      const std::string& var_name) {
  fprintf(file, "%s = [", var_name.c_str());
  for (auto output: outputs) {
    auto var = fetch_output_var(output, var_name);
    fprintf(file, "%g, ", var);
  }
  fprintf(file, "]\n");
}

}

namespace skywalker {

void write_py_module(const std::vector<InputData>& inputs,
                     const std::vector<OutputData>& outputs,
                     const char* py_module_name) {
  FILE* file = fopen(py_module_name, "w");
  fprintf(file, "# This file was automatically generated by skywalker (HAERO edition).\n\n");
  fprintf(file, "# Object is just a dynamic container that stores input/output data.\n");
  fprintf(file, "class Object(object):\n");
  fprintf(file,  "    pass\n\n");

  // Write input data.
  fprintf(file, "# Input is stored here.\n");
  fprintf(file, "input = Object()\n");
  write_input_var(file, inputs, "temperature");

  // Write output data.
  fprintf(file, "# Output data is stored here.\n");
  write_output_var(file, outputs, "num_concs");

  fclose(file);
}

}
