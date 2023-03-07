# This file loads all necessary modules for building Haero on Marianas (PNNL).

module purge

export CI_CACHE_ROOT=/qfs/projects/eagles/manc568/ci-store

export GCC_VERSION=9.1.0
export OMPI_VERSION=4.1.1
export CUDA_VERSION=11.0

module load gcc/$GCC_VERSION
module load openmpi/$OMPI_VERSION
module load cuda/$CUDA_VERSION
module load cmake/3.19.6
module load git

export CI_CACHE="${CI_CACHE_ROOT}/gcc@${GCC_VERSION}_openmpi@${OMPI_VERSION}_cuda@${CUDA_VERSION}"

export LD_LIBRARY_PATH="${CI_CACHE}/lib:${LD_LIBRARY_PATH}"

export YAMLCPP_INCLUDE_DIR=$CI_CACHE/include
export YAMLCPP_LIBRARY_DIR=$CI_CACHE/lib
export YAMLCPP_LIBRARY=libyaml-cpp.a
 
export CC=mpicc
export CXX=mpic++
export FC=mpifort

# export OMPI_CC=gcc OMPI_FC=gfortran OMPI_CXX=g++
