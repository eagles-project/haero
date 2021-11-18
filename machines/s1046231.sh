# This file sets up the environment variables for building haero on s1046231.srn.sandia.gov
#
# It's a 2019 macbook pro using apple clang compilers and openmpi 4.0.5
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

export EKAT_INCLUDE_DIR=$HOME/EKAT/install/include
export EKAT_LIBRARY_DIR=$HOME/EKAT/install/lib
export EKAT_LIBRARY=libekat.a

export YAMLCPP_INCLUDE_DIR=$HOME/EKAT/install/include/yaml-cpp
export YAMLCPP_LIBRARY_DIR=$HOME/EKAT/install/lib
export YAMLCPP_LIBRARY=libyaml-cpp.a

export OPENBLAS_INCLUDE_DIR=/usr/local/opt/openblas/include
export OPENBLAS_LIBRARY_DIR=/usr/local/opt/openblas/lib
export OPENBLAS_LIBRARY=libopenblas.a

# this tines is built against local kokkos
# export TINES_INCLUDE_DIR=$HOME/Tines/build/install/include
# export TINES_LIBRARY_DIR=$HOME/Tines/build/install/lib
# export TINES_LIBRARY=libtines.a

# this tines is built against local ekat kokkos
export TINES_INCLUDE_DIR=$HOME/Tines/install_ekat/include
export TINES_LIBRARY_DIR=$HOME/Tines/install_ekat/lib
export TINES_LIBRARY=libtines.a

# this TChem is built against local ekat kokkos
export TCHEM_INCLUDE_DIR=$HOME/TChem/install_ekat/include
export TCHEM_LIBRARY_DIR=$HOME/TChem/install_ekat/lib
export TCHEM_LIBRARY=libtchem.a

export FC=mpifort
export CC=mpicc
export CXX=mpicxx

