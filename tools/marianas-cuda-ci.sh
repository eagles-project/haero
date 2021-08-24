#!/usr/bin/bash

set -ex
source /etc/profile.d/modules.sh
whoami
pwd
groups

# This will be the catch-all trap handler after arguments are parsed.
cleanup() {
  exit_code=$?
  echo
  echo Exit code $exit_code caught in build script.
  echo
  echo "BUILD_STATUS:$exit_code"
}

trap 'cleanup $?' EXIT

SOURCE_DIR=${WORKDIR:-$PWD}

# Use this ssh command wrapper which sets up permissions to have ssh access
# to github so submodules can be fetched
export GIT_SSH=/qfs/projects/eagles/ci-secrets/ssh_git_wrapper

# Ensure git submodules have been fetched
(
  source /etc/profile.d/modules.sh
  cd $SOURCE_DIR
  pwd
  module load git
  git --version
  git submodule update --init --recursive
  module purge
)

BUILD_DIR=${BUILD_DIR:-$SOURCE_DIR/build}
source $SOURCE_DIR/machines/marianas.pnl.gov.sh
PREFIX=$BUILD_DIR
MPI=ON
DEVICE=CUDA
DEVICE_ARCH=TURING75
PRECISION=double
PACK_SIZE=1
HAERO_COLUMN_VIEW_TYPE="DeviceType::view_1d<PackType>"
HAERO_SPECIES_COLUMN_VIEW_TYPE="DeviceType::view_2d<PackType>"
HAERO_MODE_COLUMN_VIEW_TYPE="DeviceType::view_2d<PackType>"
BUILD_TYPE=Debug
BUILD_DRIVER=ON
BUILD_CHEM_DRIVER=ON
GENERATOR="Unix Makefiles"
OPTIONS=""

# [ -d $BUILD_DIR ] && rm -rf $BUILD_DIR
[ -d $BUILD_DIR ] || mkdir -p $BUILD_DIR
cd $BUILD_DIR

if [ "$MPI" = "ON" ]; then
  OPTIONS="-DHAVE_MPI=ON"
fi
if [ "$VERBOSE" = "ON" ]; then
  OPTIONS="$OPTIONS -DCMAKE_VERBOSE_MAKEFILE=ON"
fi
if [ ! "$HDF5_LIBRARY_DIR" = "" -o ! "$HDF5_INCLUDE_DIR" = "" -o ! "$HDF5_LIBRARY" = "" -o ! "$HDF5_HL_LIBRARY" = "" ]; then
  if [ "$HDF5_LIBRARY" = "" ]; then
    echo "Error: HDF5_LIBRARY must be specified for a custom HDF5 library."
    exit 1
  fi
  if [ "$HDF5_HL_LIBRARY" = "" ]; then
    echo "Error: HDF5_HL_LIBRARY must be specified for a custom high-level HDF5 library."
    exit 1
  fi
  if [ "$HDF5_INCLUDE_DIR" = "" ]; then
    echo "Error: HDF5_INCLUDE_DIR must be specified for a custom HDF5 library."
    exit 1
  fi
  OPTIONS="$OPTIONS -DHDF5_LIBRARY=$HDF5_LIBRARY_DIR/$HDF5_LIBRARY"
  OPTIONS="$OPTIONS -DHDF5_HL_LIBRARY=$HDF5_LIBRARY_DIR/$HDF5_HL_LIBRARY"
  OPTIONS="$OPTIONS -DHDF5_INCLUDE_DIR=$HDF5_INCLUDE_DIR"
fi
if [ ! "$NETCDF_LIBRARY_DIR" = "" -o ! "$NETCDF_INCLUDE_DIR" = "" -o ! "$NETCDF_LIBRARY" = "" ]; then
  if [ "$NETCDF_LIBRARY" = "" ]; then
    echo "Error: NETCDF_LIBRARY must be specified for a custom NetCDF library."
    exit 1
  fi
  if [ "$NETCDF_INCLUDE_DIR" = "" ]; then
    echo "Error: NETCDF_INCLUDE_DIR must be specified for a custom NetCDF library."
    exit 1
  fi
  OPTIONS="$OPTIONS -DNETCDF_LIBRARY=$NETCDF_LIBRARY_DIR/$NETCDF_LIBRARY"
  OPTIONS="$OPTIONS -DNETCDF_INCLUDE_DIR=$NETCDF_INCLUDE_DIR"
fi
if [ ! "$YAMLCPP_LIBRARY_DIR" = "" -o ! "$YAMLCPP_INCLUDE_DIR" = "" -o ! "$YAMLCPP_LIBRARY" = "" ]; then
  if [ "$YAMLCPP_LIBRARY" = "" ]; then
    echo "Error: YAMLCPP_LIBRARY must be specified for a custom yaml-cpp library."
    exit 1
  fi
  if [ "$YAMLCPP_INCLUDE_DIR" = "" ]; then
    echo "Error: YAMLCPP_INCLUDE_DIR must be specified for a custom yaml-cpp library."
    exit 1
  fi
  OPTIONS="$OPTIONS -DYAMLCPP_LIBRARY=$YAMLCPP_LIBRARY_DIR/$YAMLCPP_LIBRARY"
  OPTIONS="$OPTIONS -DYAMLCPP_INCLUDE_DIR=$YAMLCPP_INCLUDE_DIR"
fi
if [ ! "$EKAT_LIBRARY_DIR" = "" -o ! "$EKAT_INCLUDE_DIR" = "" -o ! "$EKAT_LIBRARY" = "" ]; then
  if [ "$EKAT_LIBRARY" = "" ]; then
    echo "Error: EKAT_LIBRARY must be specified for a custom EKAT library."
    exit 1
  fi
  if [ "$EKAT_INCLUDE_DIR" = "" ]; then
    echo "Error: EKAT_INCLUDE_DIR must be specified for a custom EKAT library."
    exit 1
  fi
  OPTIONS="$OPTIONS -DEKAT_LIBRARY=$EKAT_LIBRARY_DIR/$EKAT_LIBRARY"
  OPTIONS="$OPTIONS -DEKAT_INCLUDE_DIR=$EKAT_INCLUDE_DIR"
fi
if [ ! "$OPENBLAS_LIBRARY_DIR" = "" -o ! "$OPENBLAS_INCLUDE_DIR" = "" -o ! "$OPENBLAS_LIBRARY" = "" ]; then
  if [ "$OPENBLAS_LIBRARY" = "" ]; then
    echo "Error: OPENBLAS_LIBRARY must be specified for a custom OpenBLAS library."
    exit 1
  fi
  if [ "$OPENBLAS_INCLUDE_DIR" = "" ]; then
    echo "Error: OPENBLAS_INCLUDE_DIR must be specified for a custom OpenBLAS library."
    exit 1
  fi
  OPTIONS="$OPTIONS -DOPENBLAS_LIBRARY=$OPENBLAS_LIBRARY_DIR/$OPENBLAS_LIBRARY"
  OPTIONS="$OPTIONS -DOPENBLAS_INCLUDE_DIR=$OPENBLAS_INCLUDE_DIR"
fi
if [ ! "$TINES_LIBRARY_DIR" = "" -o ! "$TINES_INCLUDE_DIR" = "" -o ! "$TINES_LIBRARY" = "" ]; then
  if [ "$TINES_LIBRARY" = "" ]; then
    echo "Error: TINES_LIBRARY must be specified for a custom Tines library."
    exit 1
  fi
  if [ "$TINES_INCLUDE_DIR" = "" ]; then
    echo "Error: TINES_INCLUDE_DIR must be specified for a custom Tines library."
    exit 1
  fi
  OPTIONS="$OPTIONS -DTINES_LIBRARY=$TINES_LIBRARY_DIR/$TINES_LIBRARY"
  OPTIONS="$OPTIONS -DTINES_INCLUDE_DIR=$TINES_INCLUDE_DIR"
fi
if [ ! "$TCHEM_LIBRARY_DIR" = "" -o ! "$TCHEM_INCLUDE_DIR" = "" -o ! "$TCHEM_LIBRARY" = "" ]; then
  if [ "$TCHEM_LIBRARY" = "" ]; then
    echo "Error: TCHEM_LIBRARY must be specified for a custom TChem library."
    exit 1
  fi
  if [ "$TCHEM_INCLUDE_DIR" = "" ]; then
    echo "Error: TCHEM_INCLUDE_DIR must be specified for a custom TChem library."
    exit 1
  fi
  OPTIONS="$OPTIONS -DTCHEM_LIBRARY=$TCHEM_LIBRARY_DIR/$TCHEM_LIBRARY"
  OPTIONS="$OPTIONS -DTCHEM_INCLUDE_DIR=$TCHEM_INCLUDE_DIR"
fi
if [ ! "$EXTRA_LDFLAGS" = "" ]; then
  OPTIONS="$OPTIONS -DHAERO_EXTRA_LDFLAGS=$EXTRA_LDFLAGS"
fi
rm -f CMakeCache.txt
if [ "$PACK_SIZE" = "1" -a "$DEVICE" = "CPU" ]; then
  OPTIONS="$OPTIONS -DCMAKE_Fortran_COMPILER=$FC"
fi

cmake \
  -DHAERO_DISABLE_SUBMODULE_CHECKS=ON \
  -DCMAKE_INSTALL_PREFIX:PATH=$PREFIX \
  -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
  -DCMAKE_C_COMPILER=$CC \
  -DCMAKE_CXX_COMPILER=$CXX \
  -DHAERO_PRECISION=$PRECISION \
  -DHAERO_ENABLE_DRIVER=$BUILD_DRIVER \
  -DHAERO_ENABLE_CHEM_DRIVER=$BUILD_CHEM_DRIVER \
  -DHAERO_DEVICE=$DEVICE \
  -DHAERO_DEVICE_ARCH=$DEVICE_ARCH \
  -DHAERO_PACK_SIZE=$PACK_SIZE \
  -DHAERO_COLUMN_VIEW_TYPE="$HAERO_COLUMN_VIEW_TYPE" \
  -DHAERO_SPECIES_COLUMN_VIEW_TYPE="$HAERO_SPECIES_COLUMN_VIEW_TYPE" \
  -DHAERO_MODE_COLUMN_VIEW_TYPE="$HAERO_MODE_COLUMN_VIEW_TYPE" \
  -DHAERO_DOC_INCLUDE_COMMENTS=$INCLUDE_DESIGN_DOC_COMMENTS \
  $OPTIONS \
  -G "$GENERATOR" \
  $SOURCE_DIR || exit 1

make VERBOSE=1 -j `nproc` || exit 1

# ctest -VV

echo
echo Build finished successfully
echo
