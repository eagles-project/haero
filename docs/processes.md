# Aerosol Processes

Haero offers implementations for each of the important stages in the aerosol
life cycle. Here we describe the physics of each stage, the approximations
made by the associated parametrizations, and the implementations of the
underlying prognostic and diagnostic processes.

Each process is independent of all other processes. In other words, state
variables are taken directly from the prognostics and diagostics passed
to the process, and no assumptions are made about how these state variables
were computed.

!!! info "Aerosol process"
    An aerosol process computes a *tendency*, which provides the right-hand-side
    for a differential equation solved either by Haero or its model. In this
    sense, an aerosol process is *prognostic*.

Processes may need intermediate quantities using the available state data and
metadata to evaluate their tendencies. These intermediate quantities are
*diagnostic* in the sense they they can be computed, or diagnosed, from the
current state variables; they are not evolved in time by either Haero or its
host model. They are computed by diagnostic functions; storage of these
quantities depends on their use in various aerosol processes and/or the host
model.

!!! info "Diagnostic function"
    A diagnostic function updates any relevant diagnostic variables in place at
    simulation time $t$.

    Diagnostic functions are associated with the design of Kokkos kernel
    *functors*. In Haero, we apply a design priniciple: each diagnostic
    function updates exactly 1 diagnostic quantity. So each diagnostic function
    is launched as a separate kernel.

!!! info "Independent aerosol processes"
    Every aerosol process is independent of other aerosol processes. Specifically:
    a process takes a set of state variables and diagnostics, and *using only
    these variables*, it computes a set of tendencies for the prognostic state
    variables at simulation time $t$.

This process independence is a dramatic departure from prior implementations of
aerosol physics in MAM. Any sequential coupling of processes must be performed
by a host model.

With this assumption, it's possible to recover the legacy MAM algorithms by
constructing a coupling procedure that follows the relevant assumptions. But
it's also possible to construct more sophisticated coupling processes that take
advantage of advances in numerical analysis.

