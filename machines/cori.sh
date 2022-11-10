# This file loads all necessary modules for building Haero on NERSC Cori (KNL nodes).

# Load relevant modules
module purge
module load modules
module load cmake
module load PrgEnv-intel
module load impi
module load cray-hdf5
module load cray-netcdf

# Set third-party library paths
export HDF5_INCLUDE_DIR=$HDF5_DIR/include
export HDF5_LIBRARY_DIR=$HDF5_DIR/lib
export HDF5_LIBRARY=libhdf5.a
export HDF5_HL_LIBRARY=libhdf5_hl.a

export NETCDF_INCLUDE_DIR=$NETCDF_DIR/include
export NETCDF_LIBRARY_DIR=$NETCDF_DIR/lib
export NETCDF_LIBRARY=libnetcdf.a

export NETCDFF_INCLUDE_DIR=$NETCDF_DIR/include
export NETCDFF_LIBRARY_DIR=$NETCDF_DIR/lib
export NETCDFF_LIBRARY=libnetcdff.a

# Set relevant compilers and flags.
export CC=mpiicc
export CXX=mpiicpc
export FC=mpiifort
