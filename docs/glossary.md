$$
\newcommand{\dsub}[1]{_{_{#1}}}
\newcommand{\specidx}{L}
\newcommand{\amass}[1]{q\dsub{m,#1,i}}
\newcommand{\vmass}[1]{q\dsub{v,#1}}
\newcommand{\mw}[1]{M\dsub{w,#1}}
\newcommand{\rsag}{\rm H_2SO_4}
\newcommand{\sag}{H_2SO_4}
\newcommand{\rasf}{\rm SO_4}
\newcommand{\asf}{SO_4}
$$

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
Relevant units are given in square brackets next to their symbols. [-]
indicates that a quantity is unitless.

* $T$: air temperature [K]
* $p$: air pressure [Pa]
* $M_s$: molecular weight of chemical species $s$ [kg/mol].
  Note that $M_{air}$ = 28.966$~kg~$\cdot$~kmol$^{-1}$ in E3SM.
* $\mathscr{R}$: universal gas constant [J/mol/K]
* $R_s = \mathscr{R}/M_s$: gas constant specific to species $s$ [J/kg/K]
* $c\dsub{air} = \frac{p}{\mathscr{R}T}$ [kmol air / m$^3$ air]: molar
  concentration of dry air
* $\rho\dsub{\specidx}$ [kg/m$^3$]: density of species $\specidx$ (gas, liquid
  or solid). The volume used to calculate the density is the volume of air in
  which the species is suspended.
* $i$, $j$ [-]: indices of log-normal modes
* $I$ [-]: total number of log-normal modes
* $D\dsub{p}$ [m]: particle diameter
* $N\dsub{i}$ [\#/m$^{-3}$]: number concentration of particles in mode $i$
* $n\dsub{i} (D\dsub{p})$ [m$^{-4}$]: size distribution of aerosol
  particles in mode $i$, expressed as a function of $D\dsub{p}$
  (\refeq{n_Dp})
* $n\dsub{i} (\ln D\dsub{p})$ [m$^{-3}$]: size distribution of aerosol
  particles in mode $i$, expressed as a function of $\ln D\dsub{p}$
  (\refeq{n_lnDp})
* $n\dsub{norm,i} (\ln D\dsub{p})$ [-]: $n\dsub{i} (\ln D\dsub{p})$
  normalized by the mode's number concentration $N\dsub{i}$
* $D\dsub{gn,d,i}$ [m]: geometric mean of dry particle diameter $D\dsub{p}$ in
  mode $i$
* $D\dsub{gn,w,i}$ [m]: geometric mean of wet particle diameter $D\dsub{p}$ in
  mode $i$
* $\sigma\dsub{g,i}$ [m]: geometric standard deviation of $D\dsub{p}$ in mode
  $i$
* $V\dsub{d,i}$ [m$^3$ particles/m$^3$ dry air]: volume concentration of dry
  aerosol particles in mode $i$
* $V\dsub{w,i}$ [m$^3$ particles/m$^3$ dry air]: volume concentration of wet
  aerosol particles in mode $i$
* $V\dsub{\specidx,i}$ [m$^3$ species $\specidx$/m$^3$ dry air]: dry volume of
  species $\specidx$ in mode $i$
* $q\dsub{n,i}$ [\#/kmol]: total number mixing ratio of particles in mode $i$
* $\amass{\specidx}$ [kmol species $\specidx$/kmol dry air in microphysics,
  kg species $\specidx$/kg dry air in dry/wet deposition]: mass mixing ratio for
  aerosol species $\specidx$ in mode $i$
* $\rho\dsub{d,i}$ [kg dry aerosol particles/m$^3$ dry aerosol particles]:
  density of dry aerosol particles in mode $i$
* $\rho_{w,i}$ [kg wet particles/m$^3$ wet particles]: density of wet aerosol
  particles in mode $i$
* $\vmass{\specidx}$ [kmol species $s$/kmol dry air in microphysics, elsewhere
  kg species $s$/kg dry air]: mass mixing ratio of gas vapor species
* $\Delta t_{phys}$ [s]: time step size for aerosol physics
