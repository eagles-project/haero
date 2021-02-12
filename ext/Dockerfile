# This builds a Docker image with Haero's third-party libraries pre-built
# and installed in /opt/haero.
FROM ubuntu:20.10

RUN apt-get update && apt-get install -y --no-install-recommends \
  autoconf \
  bash-completion \
  cmake \
  curl \
  g++ \
  gcc \
  gfortran \
  git \
  gfortran
  m4 \
  make \
  pkg-config \
  zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

# Process build-time arguments.
ARG BUILD_TYPE=Debug
ARG PRECISION=double
ARG PACK_SIZE=1

# We build Haero in /haero
WORKDIR /haero
COPY .. .
RUN \
  ./setup build && \
  cd build && \
  cmake \
    -DCMAKE_INSTALL_PREFIX=/opt/haero \
    -DCMAKE_BUILD_TYPE=$BUILD_TYPE
    -DCMAKE_CXX_COMPILER=mpicxx \
    -DCMAKE_C_COMPILER=mpicc \
    -DCMAKE_Fortran_COMPILER=mpif90 \
    -DHAERO_PRECISION=$PRECISION \
    -DHAERO_DEVICE=CPU \
    -DHAERO_DEVICE_ARCH=AMDAVX \
    -DHAERO_PACK_SIZE=$PACK_SIZE
    -DHAERO_ENABLE_DRIVER=ON \
    -G "Unix Makefiles" \
    .. && \
  make -j && \
  make install && \
  cd / && \
  rm -rf /haero

LABEL maintainer='Jeffrey N. Johnson <jeff@cohere-llc.com>'
LABEL description='Ubuntu with HAERO third-party libraries'
