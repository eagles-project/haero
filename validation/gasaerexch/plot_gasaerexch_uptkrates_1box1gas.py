# This script generates plots for cross-validating the gas aerosol exchange
# process. 
# After the skywalker_gasaerexch_uptkrates_1box1gas_plot test has been run
# it creates a skywalker_gasaerexch_uptkrates_1box1gas_plot.py file that
# this script reads.  This script must be run in the directory that the
# skywalker_gasaerexch_uptkrates_1box1gas_plot.py file is in.

import os, sys, importlib
import matplotlib.pyplot as plt
from   matplotlib import colors
import matplotlib.tri as tri
import numpy as np

data = importlib.import_module('skywalker_gasaerexch_uptkrates_1box1gas_plot')

dgncur_awet = np.array(data.input.dgncur_awet).T
uptkaer     = np.array(data.output.uptkaer).T

nplt = int(len(data.input.dgncur_awet[0]))
npts = int(len(data.input.temp)/nplt)
ntmp = int(len(data.input.temp)/npts)

fig, ax = plt.subplots(nrows=1, ncols=nplt, figsize=(15, 10))
fig_title = 'Num Mixing Ratio #/kmol-air: {:7.3g}    ln(sigmag): {:6.3f} \n X-axis is median wet diameter (m)'
fig.suptitle(fig_title.format(data.input.aernum[0][0], data.input.lnsg[0][0]), fontsize=20)
for j in range(ntmp) :
  for i in range(nplt) :
    lab = label='T={}'.format(data.input.temp[npts*j])
    ax[i].plot(dgncur_awet[i][npts*j:npts*(j+1)], uptkaer[i][npts*j:npts*(j+1)], label=lab)

ax[0].set_ylabel("Uptake Rate", fontsize=20)
for i in range(nplt) :
  ax[i].legend()

plt.show()
