# Glossary

This section contains definitions of important terms and mathematical quantities
used in Haero.

## Terminology

Here are terms you may encounter when reading about Haero's aerosol processes.

* **Cloud Condensation Nuclei (CCN)**:
* **Precursors**:
* **Secondary Aerosols**:

# Mathematical Notation

The following symbols are used in equations within Haero's aerosol processes.
Relevant units are given in square brackets next to their symbols. $[-]$
indicates that a quantity is unitless.

* $T$: air temperature [K]
* $p$: air pressure [Pa]
* $\mu_s$: molecular weight of chemical species $s$ [kg/mol]
* $\mathscr{R}$: universal gas constant [J/mol/K]
* $R_s = \mathscr{R}/\mu_s$: gas constant specific to species $s$ [J/kg/K]
* $c_{air}$: molar concentration of dry air [kmol air / m$^3$ air]:
  $$
    c\dsub{air} = \frac{p}{\mathscr{R}T} = \frac{\rho\dsub{air}}{M\dsub{air}}
  $$,
  where $\rho\dsub{air}$ and $M\dsub{air}$ are the density and molecular weight
  of dry air, respectively. Note that
  $M_{air}$ = 28.966$~kg~$\cdot$~kmol$^{-1}$ in E3SM.

  \item $\rho\dsub{\specidx}$ [kg/m$^3$]: density of species $\specidx$
        (gas, liquid or solid)
    \hwc{For liquid- and solid-phase species there seems to be two types of
      densities:
      \begin{enumerate}
        \item mass divided by the volume of air in which they are suspended
        \item mass divided by the volume the liquid or solid matter occupies
      \end{enumerate}
      $\rho\dsub{\specidx}$ seems to be fall in to the first type. Let's watch
      out for places of potential confusion.}

  \item $i$, $j$: indices of log-normal modes

  \item $I$: total number of log-normal modes

  \item $D\dsub{p}$: particle diameter, unit: $\rm m$

  \hwc{Which of the following notations need to be distinguished for
       interstitial and cloud-born aerosols?}
  \jsc{We at least need different notations for interstitial and cloud-borne
       aerosols for $\amass{\specidx}$ and $q\dsub{n,i}$ in the rename process.
       All the other parameters are related to $\amass{\specidx}$ and
       $q\dsub{n,i}$ so we may also need to use different notation.}

  \item $N\dsub{i}$ [\#/m$^{-3}$]: number concentration of particles in mode $i$

  \item $n\dsub{i} (D\dsub{p})$ [m$^{-4}$]: size distribution of aerosol
        particles in mode $i$, expressed as a function of $D\dsub{p}$
        (\refeq{n_Dp})

  \item $n\dsub{i} (\ln D\dsub{p})$ [m$^{-3}$]: size distribution of aerosol
        particles in mode $i$, expressed as a function of $\ln D\dsub{p}$
        (\refeq{n_lnDp})

  \item $n\dsub{norm,i} (\ln D\dsub{p})$ [-]: $n\dsub{i} (\ln D\dsub{p})$
        normalized by the mode's number concentration $N\dsub{i}$

  \item $D\dsub{gn,d,i}$ [m]: geometric mean of dry particle diameter
        $D\dsub{p}$ in mode $i$

  \item $D\dsub{gn,w,i}$ [m]: geometric mean of wet particle diameter
        $D\dsub{p}$ in mode $i$

  \item $\sigma\dsub{g,i}$ [-]: geometric standard deviation of $D\dsub{p}$ in
        mode $i$
        \hwc{Why is it unitless?}
        \jsc{The unit of $\sigma\dsub{g,i}$ should be the same as that of
             $D\dsub{p}$. $\ln D\dsub{p}$ and $\ln \sigma\dsub{g,i}$ are
             unitless. According to Seinfeld's book, when we write
             $\ln D\dsub{p}$ and $\ln \sigma\dsub{g,i}$, we really mean
             $\ln \frac{D\dsub{p}}{1}$ and $\ln \frac{\sigma\dsub{g,i}}{1}$,
             where 1 is the ``reference'' particle diameter or standard
             deviation and not explicitly indicated.}

  \item $V\dsub{d,i}$ [m$^3$ particles/m$^3$ dry air]: volume concentration of
        dry aerosol particles in mode $i$

  \item $V\dsub{w,i}$ [m$^3$ particles/m$^3$ dry air]: volume concentration of
        wet aerosol particles in mode $i$

  \item $V\dsub{\specidx,i}$ [m$^3$ species $\specidx$/m$^3$ dry air]: dry
        volume of species $\specidx$ in mode $i$

  \item $q\dsub{n,i}$ [\#/kmol]: total number mixing ratio of particles in mode
        $i$

  \item $\amass{\specidx}$ [kmol species $\specidx$/kmol dry air] in microphysics,
        [kg species $\specidx$/kg dry air] in dry/wet deposition: mass mixing
        ratio for aerosol species $\specidx$ in mode $i$
        % listed in Section~\ref{sec:MAM_procs_and_eqns})

* $\rho\dsub{d,i}$ [kg dry aerosol particles/m$^3$ dry aerosol particles]:
        density of dry aerosol particles in mode $i$

* $\rho_{w,i}$: density of wet aerosol particles in mode $i$
  [kg wet particles/m$^3$ wet particles]

* $\vmass{\specidx}$: mass mixing ratio of gas vapor species
  [kmol species $s$/kmol dry air in microphysics, kg species $s$/kg dry air]
   elsewhere]
* $\Delta t_{phys}$: time step size for aerosol physics [s]
