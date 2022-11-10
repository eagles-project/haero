# Installation

Haero builds and runs on the following platforms:

* Mac and Linux laptops and workstations
* NERSC Cori
* Compy and Constance at PNNL

## Required Software

To build Haero, you need:

* [CMake v3.12+](https://cmake.org/)
* GNU Make
* reliable C and C++ compilers
* a good Fortran compiler, such as GNU `gfortran` or Intel's `ifort` compiler
* a working MPI installation (like [OpenMPI](https://www.open-mpi.org/) or
  [Mpich](https://www.mpich.org/)).

You can obtain all of these (except perhaps your favorite Fortran compiler)
freely on the Linux and Mac platforms. On Linux, just use your favorite package
manager. On a Mac, you can get the Clang C/C++ compiler by installing XCode, and
then use a package manager like [Homebrew](https://brew.sh/) or
[MacPorts](https://www.macports.org/) to get the rest.

For example, to download the relevant software on your Mac using Homebrew, type

```
brew install cmake gfortran openmpi
```

## Clone the Repository

First, go get the [source code](https://github.com/eagles-project/haero)
at GitHub:

=== "SSH"
    ```
    git clone git@github.com:eagles-project/haero.git
    ```
=== "HTTPS"
    ```
    git clone https://github.com/eagles-project/haero.git
    ```

This places a `haero` folder into your current path.

## Configure Haero

Haero uses CMake, and accepts a number of options that specify how it should be
built. In order to simplify the build process, we've provided a simple `setup`
script that generates a shell script you can run to invoke CMake with the
appropriate options set.

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
   configuration options. Then run `source config.sh` to configure the model
   (this ensures that any settings by machine-specific scripts stick around in
   your environment).

If you prefer, you can fish the options out of the `setup` script (or your
generated `config.sh` file) and feed them directly to CMake.

## Build, Test, and Install Haero

After you've configured Haero, you can build it:

1. From the build directory, type `make -j` to build the library. (If you've
   configured your build for a GPU, place a number after the `-j` flag, as in
   `make -j 8`).
4. To run tests for the library (and the driver, if configured), type
   `make test`.
5. To install the model to the location indicated by `PREFIX` in your
   `config.sh` script (or `CMAKE_INSTALL_PREFIX`, if you specified it manually),
   type `make install`. By default, products are installed in `include`, `lib`,
   `bin`, and `share` ѕubdirectories within your build directory.

## Making code changes and rebuilding

Notice that you must build Haero in a  **build tree**, separate from its source
trees. This is standard practice in CMake-based build systems, and it allows you
to build several different configurations without leaving generated and compiled
files all over your source directory. However, you might have to change the way
you work in order to be productive in this kind of environment.

When you make a code change, make sure you build from the build directory that
you created in step 1 above:

```
cd /path/to/haero/build
make -j
```

You can also run tests from this build directory with `make test`.

This is very different from how some people like to work. One method of making
this easier is to use an editor in a dedicated window, and have another window
open with a terminal, sitting in your `build` directory. If you're using a fancy
modern editor, it might have a CMake-based workflow that handles all of this for
you.

The build directory has a structure that mirrors the source directory, and you
can type `make` in any one of its subdirectories to do partial builds. In
practice, though, it's safest to always build from the top of the build tree.

## Generating Documentation

Haero's documentation is built using [Material for Mkdocs](https://squidfunk.github.io/mkdocs-material/),
which is a static website generator with lots of features. Currently, though,
the Haero repository is private, so we don't publish the documentation to a web
site. Instead, if you've [installed Mkdocs](https://squidfunk.github.io/mkdocs-material/getting-started/),
you can [run a local server](https://squidfunk.github.io/mkdocs-material/creating-your-site/#previewing-as-you-write)
to view the documentation.

## FAQ

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

