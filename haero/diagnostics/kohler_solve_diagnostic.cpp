#include "kohler_solve_diagnostic.hpp"
#include <sstream>
#include <iomanip>

namespace haero {

template <typename T>
std::string KohlerPolynomial<T>::mathematica_verification_program(const int n, const std::string& output_dir) {
  std::ostringstream ss;
  ss << "kelvinCoeff = " << std::fixed << std::setprecision(16) << kelvin_droplet_effect_coeff << ";\n";
  ss << "rhMin = " << rel_humidity_min << ";\n";
  ss << "rhMax = " << rel_humidity_max << ";\n";
  ss << "hygMin = " << hygro_min << ";\n";
  ss << "hygMax = " << hygro_max << ";\n";
  ss << "radMin = " << dry_radius_min_microns << ";\n";
  ss << "radMax = " << dry_radius_max_microns << ";\n";
  ss << "nn = " << n << ";\n";
  ss << "drh = (rhMax - rhMin)/(nn-1);\n";
  ss << "dhyg = (hygMax - hygMin)/(nn-1);\n";
  ss << "drad = (radMax - radMin)/(nn-1);\n";
  ss << "kinputs = Flatten[Table[{rhMin + i*drh, hygMin + j*dhyg, radMin + k*drad}, {i, 0, nn - 1}, {j, 0, nn - 1}, {k, 0, nn - 1}], 2];\n";
  ss << "kroots = Flatten[Table[NSolve[Log[kinputs[[i]][[1]]] r^4 + kelvinCoeff*(kinputs[[i]][[3]]^3 - r^3) + (kinputs[[i]][[2]] - Log[kinputs[[i]][[1]]]) kinputs[[i]][[3]]^3 r == 0 && r > 0, r, Reals], {i, Length[kinputs]}]];\n";
  ss << "kout = Table[r /. kroots[[i]], {i, Length[kroots]}];\n";
  ss << "Export[\"" << output_dir << "/mm_kohler_roots.txt\", kout];\n";
  return ss.str();
}

// ETI
template struct KohlerPolynomial<PackType>;

} // namespace haero
