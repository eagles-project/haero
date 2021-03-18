# This file loads all necessary modules for building MAM on PNNL's Compy Machine

source /etc/profile.d/modules.sh

# Load relevant modules (GNU 10.2.0)
module purge && module load cmake/3.19.6  gcc/10.2.0   mvapich2/2.3.5   netcdf/4.6.3


# Set third-party library paths
export HDF5_DIR=/share/apps/hdf5/1.10.5/serial
export HDF5_INCLUDE_DIR=$HDF5_DIR/include
export HDF5_LIBRARY_DIR=$HDF5_DIR/lib
export HDF5_LIBRARY=libhdf5.a
export HDF5_HL_LIBRARY=libhdf5_hl.a

export NETCDF_INCLUDE_DIR=$NETCDF_INCLUDE
export NETCDF_LIBRARY_DIR=$NETCDF_LIB
export NETCDF_LIBRARY=libnetcdf.a

export NETCDFF_INCLUDE_DIR=$NETCDF_INCLUDE_DIR
export NETCDFF_LIBRARY_DIR=$NETCDF_LIBRARY_DIR
export NETCDFF_LIBRARY=libnetcdff.a

# Set relevant compilers.
export FC=mpif90
export CXX=mpicxx

# Set extra LDFLAGS
export EXTRA_LDFLAGS=-lcurl
