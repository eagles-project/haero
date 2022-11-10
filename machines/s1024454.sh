# This file sets up the environment variables for building Haero on s1024454.srn.sandia.gov
#
# It's a rhel7 linux box with an intel skylake gold (24 cores with 2 threads per core over 2 sockets)
# with an nvidia maxwell gpu.
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


