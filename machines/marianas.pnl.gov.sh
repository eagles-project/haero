# This file loads all necessary modules for building MAM on PNNL's Marianas system

module purge
module load gcc/7.5.0 openmpi/4.1.0 cuda/10.2.89 cmake/3.19.6

module use -a /qfs/projects/eagles/modules/linux-centos7-broadwell
module load eagles-hdf5/1.10.7/openmpi-4.1.0/gcc-7.5.0/hash-oatrxcn
module load eagles-netcdf-c/4.8.0/openmpi-4.1.0/gcc-7.5.0/hash-ogxvk6f
module load eagles-netcdf-fortran/4.5.3/openmpi-4.1.0/gcc-7.5.0/hash-5qynzur
module load eagles-openblas/0.3.15/gcc-7.5.0/hash-jk2siiz
module load eagles-yaml-cpp/0.6.3/gcc-7.5.0/hash-okmj5ss
module load eagles-zlib/1.2.11/gcc-7.5.0/hash-gfg64sy

export HDF5_LIBRARY=libhdf5.a
export HDF5_HL_LIBRARY=libhdf5_hl.a

export NETCDF_INCLUDE_DIR=$NETCDF_C_INCLUDE_DIR
export NETCDF_LIBRARY_DIR=$NETCDF_C_LIBRARY_DIR
export NETCDF_LIBRARY=libnetcdf.a

export NETCDFF_INCLUDE_DIR=$NETCDF_FORTRAN_INCLUDE_DIR
export NETCDFF_LIBRARY_DIR=$NETCDF_FORTRAN_LIBRARY_DIR
export NETCDFF_LIBRARY=libnetcdff.a

export OPENBLAS_LIBRARY=libopenblas.so

export YAMLCPP_INCLUDE_DIR=$YAML_CPP_INCLUDE_DIR
export YAMLCPP_LIBRARY_DIR=$YAML_CPP_LIBRARY_DIR
export YAMLCPP_LIBRARY=libyaml-cpp.so

export CC=mpicc
export CXX=mpic++
export FC=mpifort

export OMPI_CC=gcc OMPI_FC=gfortran OMPI_CXX=g++
