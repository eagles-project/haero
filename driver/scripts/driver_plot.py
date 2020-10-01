import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from datetime import datetime

gravity_m_per_s2 = 9.80616
mSize = 8.0
mWidth = mSize/4.0

def make_parser():
  parser = argparse.ArgumentParser(description="Plot output from HAERO driver",
    usage="python %(prog)s <data_file.nc> [-t, --time_idx] [-c, --col_idx]")
  parser.add_argument("filename", help="HAERO Driver netCDF file",
    default=argparse.SUPPRESS)
  parser.add_argument("--col_idx", "-c", nargs='?', type=int, default=0,
    help="column index for requested data [default=0]")
  parser.add_argument("--time_idx", "-t", nargs='?', type=int, default=0,
    help="time index for requested data [default=0]")
#   parser.add_argument("--vars, -v", nargs='*', help="variables to plot")
  return parser

def process_input(args):
  ncd = Dataset(args.filename, "r")
  if args.time_idx >= len(ncd.dimensions["time"]):
    raise IndexError("time index out of bounds")
  if args.col_idx >= len(ncd.dimensions["ncol"]):
    raise IndexError("column index out of bounds.")
  time_idx = args.time_idx
  col_idx = args.col_idx
  return ncd, time_idx, col_idx

def haero_atts(dset):
  haero_version = dset.__dict__["HAERO_version"]
  haero_revision = dset.__dict__["HAERO_revision"]
  return "HAERO (version, revision) = (" + haero_version + ", "+ haero_revision + ")"

def plot_column(ax0,ax1,ax2, time_idx, col_idx, dset):
  tval = ncd.variables["time"][time_idx]
  # interface quantities
  phi = np.array(ncd.variables["geopotential"][time_idx][col_idx])
  z = phi/gravity_m_per_s2/1000 # convert to km
  w = np.array(ncd.variables["vertical_velocity"][time_idx][col_idx])
  mu = np.array(ncd.variables["mu"][time_idx][col_idx])
  pi = np.array(ncd.variables["hydrostatic_pressure"][time_idx][col_idx])/100 # convert to hPa

  # midpoint quantities
  qv = np.array(ncd.variables["water_vapor_mixing_ratio"][time_idx][col_idx])
  p = np.array(ncd.variables["pressure"][time_idx][col_idx])/100 # convert to hPa
  T = np.array(ncd.variables["temperature"][time_idx][col_idx])
  thetav = np.array(ncd.variables["virtual_potential_temperature"][time_idx][col_idx])
  theta = thetav/(1 + 0.61*qv)
  exner = np.array(ncd.variables["exner"][time_idx][col_idx])
  Tv = thetav * exner

  ### Plot 0 : w and qv
  wcolor = "tab:blue"
  ax0.plot(w, pi, '.', color=wcolor)
  ax0.set_xlabel("$w$ (m/s)", color=wcolor)
  ax0.tick_params(axis='x', labelcolor=wcolor)
  ax0z = ax0.twinx()
  ax0z.plot(w,z, '.', color=wcolor)
  ax0z.set_ylabel('$z$ (km)')
  ax0.set_ylabel('$p$ (hPa)')
  ax0.invert_yaxis()
  ax02 = ax0.twiny()
  qvcolor = "tab:red"
  ax02.plot(qv,p, 'x', color=qvcolor)
  ax02.set_xlabel("$q_v$ (kg/kg)", color=qvcolor)
  ax02.tick_params(axis='x', labelcolor=qvcolor)

  ### Plot 1 : temperature variables and mu
  thetavcolor = "tab:blue"
  thetacolor = "tab:cyan"
  Tcolor = "tab:purple"
  Tvcolor= "tab:green"
  ax1.plot(thetav, p, '+', label=r"$\theta_v$", markerfacecolor="None", markeredgecolor=thetavcolor)
  ax1.plot(theta, p, 'x', label=r"$\theta$", markerfacecolor="None", markeredgecolor=thetacolor)
  ax1.plot(T, p, 'v', label="$T$", markerfacecolor="None", markeredgecolor=Tcolor)
  ax1.plot(Tv, p, '^', label="$T_v$", markerfacecolor="None", markeredgecolor=Tvcolor)
  ax1.invert_yaxis()
  ax1.set_xlabel(r"$\theta, \theta_v, T$ (K)", color=thetavcolor)
  ax1.tick_params(axis='x', labelcolor=thetavcolor)
  ax1.legend()
  mucolor = "tab:red"
  ax12 = ax1.twiny()
  ax12.plot(mu, pi, color=mucolor)
  ax12.set_xlabel("$\mu$", color=mucolor)
  ax12.tick_params(axis='x', labelcolor=mucolor)

  ### Plot 2 : reserved for parameterization variables
  ax2.text(0.5,0.5, "Reserved for", fontsize=10, wrap=True, ha='center')
  ax2.text(0.5,0.25, "parameterizations", fontsize=10,ha='center')

if __name__ == "__main__":
  parser = make_parser()
  args = parser.parse_args()

  ncd, time_idx, col_idx = process_input(args)
  hatts = haero_atts(ncd)
  print(hatts)

  plt.rc('text', usetex=True)
  plt.rc('font', family='serif')
  plt.rc('ps', useafm=True)
  plt.rc('pdf', use14corefonts=True)

  fig = plt.figure(constrained_layout=True)
  gspc = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
  ax0 = fig.add_subplot(gspc[0,0])
  ax1 = fig.add_subplot(gspc[0,1])
  ax2 = fig.add_subplot(gspc[0,2])
  plot_column(ax0, ax1, ax2, time_idx, col_idx, ncd)
  fig.text(0, 0.4, "HAERO Standalone Driver", fontsize=6)
  fig.text(0, 0.37, hatts.replace('_', '\_'), fontsize=6)
  fig.text(0.015, 0.35, "data file: " + args.filename.replace('_','\_'), fontsize=6)
  fig.text(0.015, 0.33, "time\_idx = " + str(time_idx) + ", col\_idx = " + str(col_idx),fontsize=6)
  fig.text(0.015, 0.31, "plots made at " + datetime.utcnow().strftime("%H:%MZ%b-%d-%y"), fontsize=6)
  fig.savefig('column_profiles.pdf', bbox_inches='tight')
  plt.close(fig)
