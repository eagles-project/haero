# This file sets up the environment variables for building haero on s1046231.srn.sandia.gov
#
# It's a 2019 macbook pro using apple clang compilers and openmpi 4.0.5
#

export OPENBLAS_INCLUDE_DIR=/usr/local/opt/openblas/include
export OPENBLAS_LIBRARY_DIR=/usr/local/opt/openblas/lib
export OPENBLAS_LIBRARY=libopenblas.a

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

