#include <cstdio>

#include "catch2/catch.hpp"
#include "haero/processes/merikanto2007.hpp"

using namespace haero;

// These tests aren't really unit tests in the usual sense--they are "smoke
// tests" in that they exercise the parameterizations of Merikanto et al
// (2007) and generate plot files that can be used to compare the results of
// the parameterizations with the calculations from classical nucleation
// theory (e.g. given in Chapter 11 of the Third Edition of Seinfeld and
// Pandis, 2016).

namespace {

// Create a vector of values between two endpoints v1 and v2 with a step size dv.
std::vector<PackType> create_inputs(Real v1, Real v2, Real dv) {
  size_t len = 1 + static_cast<size_t>(std::ceil((v2 - v1) / dv));
  std::vector<PackType> v(len);
  for (int i = 0; i < len; ++i) {
    v[i] = v1 + i * dv;
  }
  return v;
}

// Create a vector of output generated from traversing an input array and
// applying the given function to it.
std::vector<PackType> generate_output(
    const std::vector<PackType>& x,
    std::function<PackType(PackType)> f) {
  std::vector<PackType> y(x.size());
  for (int i = 0; i < x.size(); ++i) {
    y[i] = f(x[i]);
  }
  return y;
}

// Writes a Python script that generates a plot of a curve y(x).
void generate_curve_plot(const std::string& filename,
                         const std::string& title,
                         const std::string& x_label,
                         const std::string& x_scale,
                         const std::string& y_label,
                         const std::string& y_scale,
                         const std::vector<PackType>& x,
                         const std::vector<PackType>& y,
                         Real xmin = 1.0, Real xmax = -1.0,
                         Real ymin = 1.0, Real ymax = -1.0) {
  FILE* fp = fopen(filename.c_str(), "w");
  fprintf(fp, "import matplotlib.pyplot as plt\n\n");

  fprintf(fp, "# Plot data\n");
  fprintf(fp, "x = [");
  for (int i = 0; i < x.size(); ++i) {
    fprintf(fp, "%g, ", x[i][0]);
  }
  fprintf(fp, "]\n");
  fprintf(fp, "y = [");
  for (int i = 0; i < y.size(); ++i) {
    fprintf(fp, "%g, ", y[i][0]);
  }
  fprintf(fp, "]\n\n");

  fprintf(fp, "fig, ax = plt.subplots()\n");
  fprintf(fp, "ax.plot(x, y)\n\n");
  fprintf(fp, "ax.set(xlabel='%s', xscale='%s', ylabel='%s', yscale='%s', title='%s')\n",
          x_label.c_str(), x_scale.c_str(), y_label.c_str(), y_scale.c_str(), title.c_str());
  if (xmin < xmax) { // if x limits are specified
    fprintf(fp, "ax.set_xlim(%g, %g)\n", xmin, xmax);
  }
  if (ymin < ymax) { // if y limits are specified
    fprintf(fp, "ax.set_ylim(%g, %g)\n", ymin, ymax);
  }
  fprintf(fp, "ax.grid()\n\n");

  size_t dot = filename.find(".");
  std::string plotfile = filename.substr(0, dot) + ".png";
  fprintf(fp, "print('Generating plot: %s')\n", plotfile.c_str());
  fprintf(fp, "fig.savefig('%s')\n", plotfile.c_str());
  fclose(fp);
}

// Create a vector of output generated from traversing two input arrays and
// applying the given function to it.
std::vector<std::vector<PackType>> generate_output(
    const std::vector<PackType>& x,
    const std::vector<PackType>& y,
    std::function<PackType(PackType, PackType)> f) {
  std::vector<std::vector<PackType>> z(x.size());
  for (int i = 0; i < x.size(); ++i) {
    z[i].resize(y.size());
    for (int j = 0; j < y.size(); ++j) {
      z[i][j] = f(x[i], y[j]);
    }
  }
  return z;
}

// Writes a Python script that generates a plot of a family of curves y(x)
// parameterized by the values in p.
void generate_multicurve_plot(const std::string& filename,
                              const std::string& title,
                              const std::string& x_label,
                              const std::string& y_label,
                              const std::string& p_label,
                              const std::vector<PackType>& x,
                              const std::vector<std::vector<PackType>>& y,
                              const std::vector<PackType>& p,
                              Real xmin = 1.0, Real xmax = -1.0,
                              Real ymin = 1.0, Real ymax = -1.0) {
  FILE* fp = fopen(filename.c_str(), "w");
  fprintf(fp, "import matplotlib.pyplot as plt\n\n");

  fprintf(fp, "# Plot data\n");
  fprintf(fp, "x = [");
  for (int i = 0; i < x.size(); ++i) {
    fprintf(fp, "%g, ", x[i][0]);
  }
  fprintf(fp, "]\n");
  for (int i = 0; i < y.size(); ++i) {
    fprintf(fp, "y%d = [", i+1);
    for (int j = 0; j < y[i].size(); ++j) {
      fprintf(fp, "%g, ", y[i][j][0]);
    }
    fprintf(fp, "]\n");
  }
  fprintf(fp, "p = [");
  for (int i = 0; i < p.size(); ++i) {
    fprintf(fp, "%g, ", p[i][0]);
  }
  fprintf(fp, "]\n\n");

  for (int i = 0; i < y.size(); ++i) {
    fprintf(fp, "plt.plot(x, y%d, label='%s = %g')\n", i+1, p_label.c_str(), p[i][0]);
  }
  fprintf(fp, "plt.legend()\n");
  fprintf(fp, "plt.title('%s')\n", title.c_str());
  fprintf(fp, "plt.xlabel('%s')\n", x_label.c_str());
  fprintf(fp, "plt.ylabel('%s')\n\n", y_label.c_str());
  if (xmin < xmax) { // if x limits are specified
    fprintf(fp, "plt.xlim(%g, %g)\n", xmin, xmax);
  }
  if (ymin < ymax) { // if y limits are specified
    fprintf(fp, "plt.ylim(%g, %g)\n", ymin, ymax);
  }

  size_t dot = filename.find(".");
  std::string plotfile = filename.substr(0, dot) + ".png";
  fprintf(fp, "print('Generating plot: %s')\n", plotfile.c_str());
  fprintf(fp, "plt.savefig('%s')\n", plotfile.c_str());
  fclose(fp);
}

} // anonymous namespace

// Compare the output of this test with Merikanto et al (2007), figure 2.
// Because NH3 affects H2SO4 nucleation, the figure has a few different sets of
// curves
TEST_CASE("merikanto2007_figure_2") {
  PackType temps[2] = {235.15, 273.15};
  PackType xi_nh3s[2][2] = {
      // temperature dependent
      {0.1, 10},   // T = 235.15 K
      {10, 1000},  // T = 273.15 K
  };
  std::vector<PackType> c_h2so4s[2][2] = {
      // temperature and NH3 dependent
      {
          // T = 235.15 K
          {2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9},  // xi = 0.1
          {1e5, 2e5, 5e5, 1e6, 5e6, 1e7, 5e7, 1e8, 5e8, 1e9}   // xi = 10
      },
      {
          // T = 273.15 K
          {3e7, 1e8, 2e8, 3e8, 4e8, 5e8, 6e8, 7e8, 8e8, 9e8, 1e9},  // xi = 10
          {8e6, 5e7, 8e7, 1e8, 2e8, 5e8, 6e8, 7e8, 8e8, 9e8, 1e9}   // xi = 1000
      }};
  std::vector<PackType> Js[2][2];
  PackType rel_hum(0.5);
  for (int i = 0; i < 2; ++i) {
    PackType temp = temps[i];
    for (int j = 0; j < 2; ++j) {
      PackType xi_nh3 = xi_nh3s[i][j];
      Js[i][j] = generate_output(c_h2so4s[i][j], [=](PackType c) {
        // NOTE: there may be a glitch in the onset temperature calculation.
        // NOTE: when we use it here, the first data point for the T = 273.15 K
        // NOTE: nucleation rates is zero for both xi = 10 and xi = 1000.
        auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                        c, xi_nh3);
        return exp(log_J);
      });
    }
  }

  // Write the plot script.
  std::string filename = "merikanto2007_figure_2.py";
  FILE* fp = fopen(filename.c_str(), "w");
  fprintf(fp, "import matplotlib.pyplot as plt\n\n");

  for (int i = 0; i < 2; ++i) {
    PackType temp = temps[i];
    for (int j = 0; j < 2; ++j) {
      PackType xi_nh3 = xi_nh3s[i][j];
      fprintf(fp, "# T = %g, xi_nh3 = %g\n", temp[0], xi_nh3[0]);
      fprintf(fp, "x = [");
      for (int k = 0; k < c_h2so4s[i][j].size(); ++k) {
        fprintf(fp, "%g, ", c_h2so4s[i][j][k][0]);
      }
      fprintf(fp, "]\n");
      fprintf(fp, "y = [");
      for (int k = 0; k < Js[i][j].size(); ++k) {
        fprintf(fp, "%g, ", Js[i][j][k][0]);
      }
      fprintf(fp, "]\n");
      fprintf(fp, "label = 'T = %g, xi = %g'\n", temp[0], xi_nh3[0]);
      fprintf(fp, "plt.plot(x, y, label=label)\n");
    }
  }
  fprintf(fp, "plt.legend()\n");
  fprintf(fp, "plt.title('Ternary Nucleation rates')\n");
  fprintf(fp, "plt.xlabel('H2SO4 concentration [#/cc]')\n");
  fprintf(fp, "plt.xscale('log')\n");
  fprintf(fp, "plt.ylabel('Nucleation rate [#/cc/s]')\n\n");
  fprintf(fp, "plt.yscale('log')\n");
  fprintf(fp, "plt.xlim(1e5, 1e9)\n");
  fprintf(fp, "plt.ylim(1e-5, 1e8)\n");

  size_t dot = filename.find(".");
  std::string plotfile = filename.substr(0, dot) + ".png";
  fprintf(fp, "print('Generating plot: %s')\n", plotfile.c_str());
  fprintf(fp, "plt.savefig('%s')\n", plotfile.c_str());
  fclose(fp);
}

// Compare the output of this test with Merikanto et al (2007), figure 3.
TEST_CASE("merikanto2007_figure_3") {
  PackType temps[2] = {235.15, 273.15};
  PackType c_h2so4s[2] = {1e6, 1e9};
  PackType xi_nh3s[2][2] = {// T = 235.15
                            {0.1, 100},
                            // T = 273.15
                            {10, 1000}};
  FILE* fp = fopen("merikanto2007_figure_3.dat", "w");
  for (int i = 0; i < 2; ++i) {
    PackType temp = temps[i];
    PackType c_h2so4 = c_h2so4s[i];
    for (int j = 0; j < 2; ++j) {
      PackType xi_nh3 = xi_nh3s[i][j];
      for (int k = 0; k < 10; ++k) {
        PackType rel_hum(0.05 + k * 0.1);

        // If we're above the onset temperature, J = 0.
        auto onset_temp =
            merikanto2007::onset_temperature(rel_hum, c_h2so4, xi_nh3);
        PackType J(0);
        if (temp[0] < onset_temp[0]) {
          auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                          c_h2so4, xi_nh3);
          J = exp(log_J);
        }

        fprintf(fp, "%g %g %g %g\n", temp[0], xi_nh3[0], rel_hum[0], J[0]);
      }
    }
  }
  fclose(fp);
}

// Compare the output of this test with Merikanto et al (2007), figure 4.
TEST_CASE("merikanto2007_figure_4") {
  PackType temps[2] = {235.15, 273.15};
  PackType c_h2so4s[2][2] = {{1e7, 1e9},
                             {1e8, 1e9}};
  std::vector<PackType> xi_nh3s[2][2] = {
    // temperature and c_h2so4 dependent
    {
      // T = 235.15
      {0.1, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000}, // c_h2so4 = 1e7
      {0.1, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000}, // c_h2so4 = 1e9
    },
    {
      // T = 273.15
      {5, 10, 20, 50, 100, 200, 500, 1000}, // c_h2so4 = 1e8
      {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000}, // c_h2so4 = 1e9
    }
  };
  std::vector<PackType> Js[2][2];
  PackType rel_hum(0.5);
  for (int i = 0; i < 2; ++i) {
    PackType temp = temps[i];
    for (int j = 0; j < 2; ++j) {
      PackType c_h2so4 = c_h2so4s[i][j];
      Js[i][j] = generate_output(xi_nh3s[i][j], [=](PackType xi) {
        // NOTE: there may be a glitch in the onset temperature calculation.
        // NOTE: when we use it here, the first data point for the T = 273.15 K
        // NOTE: nucleation rates is zero for both xi = 10 and xi = 1000.
        auto log_J = merikanto2007::log_nucleation_rate(temp, rel_hum,
                                                        c_h2so4, xi);
        return exp(log_J);
      });
    }
  }

  // Write the plot script.
  std::string filename = "merikanto2007_figure_4.py";
  FILE* fp = fopen(filename.c_str(), "w");
  fprintf(fp, "import matplotlib.pyplot as plt\n\n");

  for (int i = 0; i < 2; ++i) {
    PackType temp = temps[i];
    for (int j = 0; j < 2; ++j) {
      PackType c_h2so4 = c_h2so4s[i][j];
      fprintf(fp, "# T = %g, c_h2so4 = %g\n", temp[0], c_h2so4[0]);
      fprintf(fp, "x = [");
      for (int k = 0; k < xi_nh3s[i][j].size(); ++k) {
        fprintf(fp, "%g, ", xi_nh3s[i][j][k][0]);
      }
      fprintf(fp, "]\n");
      fprintf(fp, "y = [");
      for (int k = 0; k < Js[i][j].size(); ++k) {
        fprintf(fp, "%g, ", Js[i][j][k][0]);
      }
      fprintf(fp, "]\n");
      fprintf(fp, "label = 'T = %g, c_h2so4 = %g'\n", temp[0], c_h2so4[0]);
      fprintf(fp, "plt.plot(x, y, label=label)\n");
    }
  }
  fprintf(fp, "plt.legend()\n");
  fprintf(fp, "plt.title('Ternary Nucleation rates')\n");
  fprintf(fp, "plt.xlabel('NH_3 mixing ratio [ppt]')\n");
  fprintf(fp, "plt.xscale('log')\n");
  fprintf(fp, "plt.ylabel('Nucleation rate [#/cc/s]')\n\n");
  fprintf(fp, "plt.yscale('log')\n");
  fprintf(fp, "plt.xlim(1e-1, 1e3)\n");
  fprintf(fp, "plt.ylim(1e-5, 1e8)\n");

  size_t dot = filename.find(".");
  std::string plotfile = filename.substr(0, dot) + ".png";
  fprintf(fp, "print('Generating plot: %s')\n", plotfile.c_str());
  fprintf(fp, "plt.savefig('%s')\n", plotfile.c_str());
  fclose(fp);
}
