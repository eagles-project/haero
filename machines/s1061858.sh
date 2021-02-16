# This file sets up the environment variables for building MAM on s1068158.srn.sandia.gov
#
# It's a 2020 macbook pro using Apple's clang-11.0, gfortran-10.2, and open-mpi-4.0.5
#

export HDF5_INCLUDE_DIR=$HDF5_ROOT/include
export HDF5_LIBRARY_DIR=$HDF5_ROOT/lib
export HDF5_LIBRARY=libhdf5.a
export HDF5_HL_LIBRARY=libhdf5_hl.a

export NETCDF_INCLUDE_DIR=$NETCDF_ROOT/include
export NETCDF_LIBRARY_DIR=$NETCDF_ROOT/lib
export NETCDF_LIBRARY=libnetcdf.a

export NETCDFF_INCLUDE_DIR=$NETCDF_ROOT/include
export NETCDFF_LIBRARY_DIR=$NETCDF_ROOT/lib
export NETCDFF_LIBRARY=libnetcdff.a

export FC=mpifort
export CC=mpicc
export CXX=mpicxx

