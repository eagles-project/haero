# Tools

This directory contains the following helpful items:

* `build-ext-docker-image.sh`: this script is used to generate Docker images
  containing third-party libraries for Haero with different configurations.
  It uses `Dockerfile.ext` to build the images.
* `build-docker-images.sh`: this script runs build-ext-docker-image.sh in all
  configurations required by our GitHub Actions auto-tester. It can be called
  without arguments, but you must have push access to the cohere-llc DockerHub
  image repository to run it.
* `Dockerfile.ext`: this is a `Dockerfile` used to generate a Docker image that
  contains all of the third-party libraries needed by Haero. It's used by our
  automatic testing system to accelerate builds.
* `marianas-cuda-ci.sh`: this script is used by our continuous integration (CI)
  system on the CUDA-equipped Marianas machine at PNNL.
* `update_version_info.sh`: this shell script is used by the build system to
  fetch version and git revision information. This information is written to
  a C++ source file that provides this information to library clients.
