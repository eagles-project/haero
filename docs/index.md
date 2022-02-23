# Overview

Haero is a software library that can be used to implement a high-performance
aerosol model. Haero defines a framework that allows one to create **aerosol
processes** representing various stages of the aerosol lifecycle. These
processes can be built in a way that allows them to run on CPUs or GPUs, using
a single set of source code that leverages the
[Kokkos C++ library](https://kokkos.org).

* The [Installation](installation.md) guide shows you how to build and install
  Haero on your own machine or on a supported high-performance platform.
* The [Glossary](glossary.md) introduces the terminology we use to describe
  Haero and the physical concepts it attempts to represent. Look here if you see
  a term or symbol you don't understand.
* The [Physics](physics.md) guide gives a brief description of the relevant
  governing equations that accommodate the particle size distributions inherent
  in aerosol modeling.
* The [Library](library.md) guide introduces the Haero library and its
  main abstraction, the **aerosol process**, which provides an elementaryu
  building block for implementing an aerosol model.
* The [Processes](processes.md) guide contains brief descriptions of the
  aerosol processes included with Haero itself, with references to the original
  models they implement.
* The [Driver](driver.md) guide describes Haero's standalone driver, which
  includes a one-dimensional hydrostatic dynamics package. The driver can be
  used to verify aerosol processes in a simple and self-contained context.
* In the [Testing](testing.md) guide, we describe in general terms the
  methodology we use to test Haero and its aerosol processes.

### Acknowledgements

Haero was developed by an interdisciplinary team consisting of aerosol and
atmospheric researchers, applied mathematicians, and software engineers. It was
created for the [EAGLES project](https://climatemodeling.science.energy.gov/projects/enabling-aerosol-cloud-interactions-global-convection-permitting-scales-eagles),
an effort to improve the treatment of aerosols in
[E3SM](https://climatemodeling.science.energy.gov/projects/energy-exascale-earth-system-model),
the Department of Energy's global climate model. The source code is available on
[GitHub](https://github.com/eagles-project/haero). This effort was funded by
the Office of Science's [Biological and Environmental
Research](https://science.osti.gov/ber) Program.
