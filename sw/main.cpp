#include "parse_yaml.hpp"
#include "haero/haero.hpp"
#include "haero/modal_aerosol_config.hpp"

#include "ekat/ekat_session.hpp"

namespace {

// Print driver usage information and exit.
void usage(const char* exe)
{
  fprintf(stderr, "%s: usage:\n", exe);
  fprintf(stderr, "%s <input.yml>\n", exe);
  exit(1);
}

} // anonymous namespace

int main(int argc, const char* argv[]) {
  ekat::initialize_ekat_session(argc, const_cast<char**>(argv), false);

  if (argc < 2)
    usage(argv[0]);

  // Read the input file and extract input.
  std::string input_file(argv[1]);
  auto param_walk = parse_yaml(input_file);

}
