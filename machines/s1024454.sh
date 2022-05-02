# This file sets up the environment variables for building MAM on s1024454.srn.sandia.gov
#
# It's a rhel8 linux box with an intel skylake gold (24 cores with 2 threads per core over 2 sockets)
# with an nvidia maxwell gpu.
#

export HDF5_INCLUDE_DIR=$HDF5_ROOT/include
export HDF5_LIBRARY_DIR=$HDF5_ROOT/lib
export HDF5_LIBRARY=libhdf5.a
export HDF5_HL_LIBRARY=libhdf5_hl.a

export NETCDF_INCLUDE_DIR=$NETCDF_C_ROOT/include
export NETCDF_LIBRARY_DIR=$NETCDF_C_ROOT/lib
export NETCDF_LIBRARY=libnetcdf.a

export NETCDFF_INCLUDE_DIR=$NETCDF_FORTRAN_ROOT/include
export NETCDFF_LIBRARY_DIR=$NETCDF_FORTRAN_ROOT/lib
export NETCDFF_LIBRARY=libnetcdff.a

# export OPENBLAS_INCLUDE_DIR=/home/mjschm/OpenBLAS/install/include
# export OPENBLAS_LIBRARY_DIR=/home/mjschm/OpenBLAS/install/lib
# export OPENBLAS_LIBRARY=libopenblas.a

export FC=mpifort
export CC=mpicc
export CXX=mpicxx

