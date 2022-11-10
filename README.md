![auto_test](https://github.com/jeff-cohere/haero/actions/workflows/auto_test.yml/badge.svg)
![push_mirror](https://github.com/jeff-cohere/haero/actions/workflows/push_mirror.yml/badge.svg)

# Haero: A High-Performance Aerosol Library

Haero is a library that contains parameterizations that describe the dynamics of
aerosols in the atmosphere. Rather than providing an aerosol package to be
coupled in a specific way with a host models, it provides direct access to
individual aerosol parameterizations tied to specific governing equations. This
low-level approach allows an atmospheric "host model" to use its own coupling
and time integration logic with these parameterizations.

The short-term goal of Haero is to provide the capabilities of the MAM4 package
to [E3SM](https://github.com/E3SM-Project)'s state-of-the-science cloud-resolving
atmospheric model, [SCREAM](https://github.com/E3SM-Project/scream).

## Supported Platforms

You can build and run Haero on a Mac (**Intel, not M1**) or Linux laptop or
workstation. We also support a limited number of platforms on which the Haero
model is built and tested:

* NERSC Cori
* Compy, Constance, Deception at PNNL

## Required Software

To build Haero, you need:

* [CMake v3.12+](https://cmake.org/)
* GNU Make
* reliable C and C++ compilers
* a good Fortran compiler, such as GNU `gfortran` or Intel's `ifort` compiler
* a working MPI installation (like [OpenMPI](https://www.open-mpi.org/) or
  [Mpich](https://www.mpich.org/)), if you're interested in multi-node
  parallelism.

You can obtain all of these (except perhaps your favorite Fortran compiler)
freely on the Linux and Mac platforms. On Linux, just use your favorite package
manager. On a Mac, you can get the Clang C/C++ compiler by installing XCode, and
then use a package manager like [Homebrew](https://brew.sh/) or
[MacPorts](https://www.macports.org/) to get the rest.

For example, to download the relevant software on your Mac using Homebrew, type

```
brew install cmake gfortran openmpi
```

## Building the Model

To configure Haero:

1. Make sure you have the latest versions of all the required submodules:
   ```
   git submodule update --init --recursive
   ```
2. Create a build directory by running the `setup` script from the top-level
   source directory:
   ```
   ./setup build
   ```
3. Change to your build directory and edit the `config.sh` file to select
   configuration options. Then execute `config.sh` to configure the model.
   If you're on a machine that requires modules to get access to compilers, etc,
   use `source config.sh` to make sure your environment is updated.
4. From the build directory, type `make -j` to build the library. (If you've
   configured your build for a GPU, place a number after the `-j` flag, as in
   `make -j 8`).
5. To run tests for the library (and the driver, if configured), type
   `make test`.
6. To install the model to the location indicated by `PREFIX` in your
   `config.sh` script, type `make install`. By default, products are installed
   in `include`, `lib`, `bin`, and `share` ѕubdirectories within your build
   directory.

### Making code changes and rebuilding

This project uses **build trees** that are separate from source trees. This
is standard practice in CMake-based build systems, and it allows you to build
several different configurations without leaving generated and compiled files
all over your source directory. However, you might have to change the way you
work in order to be productive in this kind of environment.

When you make a code change, make sure you build from the build directory that
you created in step 1 above:

```
cd /path/to/haero/build
make -j
```

You can also run tests from this build directory with `make test`.

This is very different from how some people like to work. One method of making
this easier is to use an editor in a dedicated window, and have another window
open with a terminal, sitting in your `build` directory.

The build directory has a structure that mirrors the source directory, and you
can type `make` in any one of its subdirectories to do partial builds. In
practice, though, it's safest to always build from the top of the build tree.

## Generating Documentation

Documentation for Haero can be built using
[`mkdocs`](https://squidfunk.github.io/mkdocs-material/).
In order to build and view the
documentation, you must download `mkdocs` and its Material theme:

```pip3 install mkdocs mkdocs-material```

Then, run `mkdocs serve` from the root directory of your Haero repo,
and point your browser to [`http://localhost:8000`](http://localhost:8000).

At this time, Haero's documentation includes an extensive design document
describing the design approach used by Haero, including high-level descriptions
of its aerosol parameterizations.

# FAQ

## Building and Rebuilding

+ **When I run config.sh, I see an error complaining about a bad fd number!**
  You probably typed `sh config.sh` to run the configuration script. It's
  actually a `bash` script. Just type `./config.sh`.
+ **How do I "reconfigure my build"?** If you want to change a compile-time
  parameter in your model, you must reconfigure and rebuild it. To do this,
  edit your `config.sh` and change the parameter as needed. Then rerun it with
  `./config.sh`. After the script finishes, executing you can type `make -j` to
  rebuild the model.
+ **A pull request has the `reconfig required` label. What does this mean?**
  A pull request with the `reconfig required` label has made a change to the
  structure of the `config.sh` script, so you must rerun `setup <build_dir>`
  to regenerate your `config.sh` script. Once you've regenerated this script,
  you can reconfigure and build as usual.

## Testing

+ **Where are testing results stored when I type `make test`?** All testing
 results are logged to `Testing/Temporary/LastTest.log` within your build
 directory. A list of tests that failed is also written to `LastTestsFailed.log`
 in that same directory.

## Source Control and Repository

+ **Git thinks I have modifications in my submodules? Git submodules are
  annoying! Help!** We agree. These warnings about "dirty" modifications
  are irritating and useless. You can get rid of them by setting the
  following config parameter for Git:

  ```
  git config --global diff.ignoreSubmodules dirty
  ```
+ **Why must I clone the submodules for libraries that I already have installed
  locally? Git submodules are annoying! Help!** We agree. The submodule
  mechanism is a leaky abstraction and doesn't allow us to easily select which
  submodules should be cloned, so we just clone them all to keep things simple.
  We'll address this issue when a solution becomes more obvious.

