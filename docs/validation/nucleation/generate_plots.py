# This script generates plots for cross-validating the particle nucleation
# process. To generate the data needed by the script, run skywalker in the
# current directory, using nucleation.yaml as input.

import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.tri as tri
import numpy as np
import haero_skywalker as haero_data

def plot_nucleation_rate_contours(filename):
    """Plot the contours of the rate of nucleation as a function of temperature
and relative humidity."""

    # Fetch scatter data from skywalker.
    RH, T, J = haero_data.input.atmosphere.relative_humidity, \
               haero_data.input.atmosphere.temperature, \
               haero_data.output.aerosols.interstitial.aitken.so4

    # Interpolate the data onto a triangulated grid.
    nx, ny = 100, 100
    RHi = np.linspace(min(RH), max(RH), nx)
    Ti = np.linspace(min(T), max(T), ny)
    triang = tri.Triangulation(RH, T)
    interpolator = tri.LinearTriInterpolator(triang, J)
#    interpolator = tri.CubicTriInterpolator(triang, J)
    Xi, Yi = np.meshgrid(RHi, Ti)
    Ji = interpolator(Xi, Yi)

    # Plot contours.
    plt.contour(RHi, Ti, Ji,
                locator=ticker.LogLocator(),
                colors='k')
    plt.contourf(RHi, Ti, Ji,
                 locator=ticker.LogLocator(),
                 cmap='jet')

    plt.xlabel('Relative humidity [-]')
    plt.ylabel('Temperature [K]')
    plt.title('Nucleation rate [#/cc]')
    plt.colorbar()
#    plt.show()
    plt.savefig(filename)

    # Side-by-side plots go here.
    #fig, axs = plt.subplots(ncols=2)
    #for ax in axs:

if __name__ == '__main__':
    plot_nucleation_rate_contours('nucleation_rate.png')

