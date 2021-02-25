# This file sets up the environment variables for building MAM on s1046231.srn.sandia.gov
#
# It's a 2019 macbook pro using macports' gcc6 compilers and openmpi 4.0.1
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

export OPENBLAS_INCLUDE_DIR=/usr/local/opt/openblas/include
export OPENBLAS_LIBRARY_DIR=/usr/local/opt/openblas/lib
export OPENBLAS_LIBRARY=libopenblas.a

export FC=mpifort
export CC=mpicc
export CXX=mpicxx

