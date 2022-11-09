# This file loads all necessary modules for building Haero on Compy (PNNL).

source /etc/profile.d/modules.sh

# Load relevant modules:
#(Intel 20.0)
module purge && module load gcc/10.2.0 openmpi/4.0.1 cmake/3.19.6

# Set relevant compilers.
export FC=gfortran
export CXX=c++
export CC=gcc

