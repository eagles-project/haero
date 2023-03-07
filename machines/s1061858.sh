# This file sets up the environment variables for building Haero on s1068158.srn.sandia.gov
#
# It's a 2020 macbook pro using Apple's clang-11.0, gfortran-11.2, and open-mpi-4.0.5
#

export FC=mpifort
export CC=mpicc
export CXX=mpicxx

export EXTRA_LDFLAGS="-L/usr/local/lib/gcc/11;-lgfortran"

