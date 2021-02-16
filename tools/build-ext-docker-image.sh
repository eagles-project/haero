#!/usr/bin/env bash

# This script builds a Docker image that contains Haero's third-party libraries
# in /opt/haero, on top of a recent Ubuntu image. Run it like so:
#
# `./build-ext-docker-image.sh <build-type> <precision> <pack-size>`
#
# The arguments are:
# <build-type> - Debug or Release (determines build optimization level)
# <precision> - single or double (floating point precision)
# <pack-size> - a positive integer (size of SIMD packs)
#
# For this script to work, Docker must be installed on your machine.
BUILD_TYPE=$1
PRECISION=$2
PACK_SIZE=$3

if [[ "$1" == "" || "$2" == "" || "$3" == "" ]]; then
  echo "Usage: $0 <build-type> <precision> <pack-size>"
  exit
fi

TAG=$BUILD_TYPE-$PRECISION-pack-size-$PACK_SIZE

# Build the image locally.
cp Dockerfile.ext ../Dockerfile
docker build -t haero-tpl:$TAG --network=host \
  --build-arg BUILD_TYPE=$BUILD_TYPE \
  --build-arg PRECISION=$PRECISION \
  --build-arg PACK_SIZE=$PACK_SIZE \
  ..
rm ../Dockerfile

# Tag the image.
docker image tag haero-tpl:$TAG coherellc/haero-tpl:$TAG

echo "To upload this image to DockerHub, use the following:"
echo "docker login"
echo "docker image push coherellc/haero-tpl:$TAG"
echo "docker logout"
