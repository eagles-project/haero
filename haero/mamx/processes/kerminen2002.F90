module kerminen2002

  use haero_precision, only: wp

  implicit none

  private

  public :: growth_rate, &
            condensation_sink, &
            growth_parameter, &
            growth_parameter_from_rate, &
            apparent_nucleation_factor, &
            apparent_nucleation_factor_from_eta

contains

  !> This function computes the growth rate @f$GR@f$ [nm/h] for particles with
  !> the given number concentration, mass density, and molecular weight at the
  !> given temperature, using KK2002 eq 21.
  !> @param [in] c The number concentration density of nucleated particles
  !> [kg/m3]
  !> @param [in] rho The mass density of nucleated particles [kg/m3]
  !> @param [in] mw The molecular weight of the nucleated particle species
  !> [kg/mol]
  !> @param [in] temp The atmospheric temperature [K]
  function growth_rate(c, rho, mw, temp) result(rate)

    real(wp), intent(in) :: c, rho, mw, temp

    real(wp) :: speed, rate

    speed = 14.7_wp * sqrt(temp) ! molecular speed [m/s]
    rate = 3.0e-9_wp * speed * mw * c / rho
  end function

  !> This function computes the condensation sink @f$CS'@f$ parameter [m-2] used
  !> to compute the @f$\eta@f$ parameter for the nucleated particle growth
  !> parameterization.
  !> @param [in] rho_air The mass density of dry air [kg/m3]
  !> @param [in] d_wet_grown The wet diameter of grown particles [nm]
  !> @param [in] c_tot The total number concentration of aerosol particles [#/cc]
  function condensation_sink(rho_air, d_wet_grown, c_tot) result(sink)
    use haero_constants, only: pi, molec_weight_dry_air

    real(wp), intent(in) :: rho_air, d_wet_grown, c_tot

    real(wp) :: alpha, Kn, beta
    real(wp) :: sink

    ! For the purposes of this calculation, we use alpha == 1 and we use the mean
    ! free path of air as computed from the air density in the calculation of the
    ! Knudsen number for the nucleation mode.
    ! NOTE: this differs from the MAM4 calculation, which uses an H2SO4
    ! NOTE: uptake rate that assumes a process ordering, which we're no
    ! NOTE: longer allowed to do.
    alpha = 1_wp ! accommodation coefficient

    ! The Knudsen number for the nucleated particles is Kn = 2 * lambda / d,
    ! where lambda is the mean free path of air, and d is the grown particle
    ! diameter. The mean free path is 1/(n * sigma), where n = rho_air/mw_air
    ! is the number density of air, and sigma = pi*d^2 is the cross section
    ! of a grown particle. Putting everything togther, we have
    !      2 * mw_air
    ! Kn = --------------- 3
    !      pi * rho_air * d
    ! TODO: should we attempt to estimate the wet number density?
    Kn = 2_wp * molec_weight_dry_air / (pi * rho_air * d_wet_grown**3)

    ! Compute the transitional correction for the condensational mass flux
    ! (Fuchs and Sutugin, 1971, or KK2002 eq 4).
    beta = &
      (1.0_wp + Kn) / (1.0_wp + 0.377_wp * Kn + 1.33_wp * Kn * (1_wp + Kn) / alpha)

    ! Compute the condensation sink from KK2002 eq 3.
    sink = 0.5_wp * d_wet_grown * beta * c_tot

  end function

  !> This function computes the growth parameter @f$\eta@f$ [nm] used in the
  !> conversion from the "real" (base) nucleation rate to the "apparent"
  !> nucleation rate:
  !> @f$J_{app} = J_{real} \exp\left[\frac{\eta}{d_f} -
  !>                                 \frac{\eta}{d_i}\right],@f$
  !> where @f$d_i@f$ and @f$d_f@f$ are the initial (nucleated) and final (grown)
  !> wet diameters of the particles in question.
  !> @param [in] c_so4 The number concentration of SO4 aerosol [#/cc]
  !> @param [in] c_nh4 The number concentration of NH4 aerosol [#/cc]
  !> @param [in] nh4_to_so4_molar_ratio The molar ratio of NH4 to SO4 [-]
  !> @param [in] temp The atmospheric temperature [K]
  !> @param [in] rel_hum The atmospheric relative humidity [-]
  !> @param [in] d_dry_crit The dry diameter of particles in a CC [nm]
  !> @param [in] d_wet_crit The wet diameter of particles in a CC [nm]
  !> @param [in] d_dry_grown The dry diameter of grown particles [nm]
  !> @param [in] rho_grown The mass density of grown particles [kg/m3]
  !> @param [in] rho_air The mass density of dry air [kg/m3]
  !> @param [in] mw_h2so4 The molecular weight of H2SO4 gas [kg/mol]
  function growth_parameter(c_so4, c_nh4, nh4_to_so4_molar_ratio, &
                            temp, rel_hum, d_dry_crit, d_wet_crit, &
                            d_dry_grown, rho_grown, rho_air, mw_h2so4) result(eta)

    real(wp), intent(in) :: c_so4, c_nh4, nh4_to_so4_molar_ratio, &
                            temp, rel_hum, d_dry_crit, d_wet_crit, &
                            d_dry_grown, rho_grown, rho_air, mw_h2so4

    real(wp) :: bounded_rel_hum, wet_dry_vol_ratio, V_frac_wet_so4, &
                cond_growth_rate, d_wet_grown, cond_sink
    real(wp) :: eta

    ! Compute the wet/dry volume ratio using the simple Kohler approximation
    ! for ammonium sulfate and bisulfate.
    bounded_rel_hum = max(0.10_wp, min(0.95_wp, rel_hum))
    wet_dry_vol_ratio = 1.0_wp - 0.56_wp / log(bounded_rel_hum)

    ! Compute the growth rate [nm/h] of new particles.

    ! Compute the fraction of the wet volume due to SO4 aerosol.
    V_frac_wet_so4 = &
      1.0_wp / (wet_dry_vol_ratio * (1.0_wp + nh4_to_so4_molar_ratio * 17.0_wp / 98.0_wp))

    ! Compute the condensation growth rate gr [nm/h] of new particles from
    ! KK2002 eq 21 for H2SO4 uptake and correct for NH3/H2O uptake.
    cond_growth_rate = growth_rate(c_so4, rho_grown, mw_h2so4, temp) / V_frac_wet_so4

    ! Wet diameter [nm] of grown particles with dry diameter d_dry_grown.
    d_wet_grown = 1e9_wp * d_dry_grown * wet_dry_vol_ratio**(1.0_wp/3.0_wp)

    ! Compute the condensation sink CS' from KK2002 eqs 3-4.
    cond_sink = condensation_sink(rho_air, d_wet_grown, c_so4 + c_nh4)

    ! Compute eta [nm] using KK2002 eq 11.
    eta = growth_parameter_from_rate(temp, d_dry_crit, d_wet_crit, d_dry_grown, &
                                     d_wet_grown, rho_grown, cond_growth_rate, &
                                     cond_sink)
  end function

  !> This function computes the growth parameter @f$\eta@f$ [nm] in terms of a
  !> given condensation growth rate GR and a given condensity sink CS'.
  !> @param [in] temp The atmospheric temperature [K]
  !> @param [in] d_dry_crit The dry diameter of particles in a CC [nm]
  !> @param [in] d_wet_crit The wet diameter of particles in a CC [nm]
  !> @param [in] d_dry_grown The dry diameter of grown particles [nm]
  !> @param [in] d_wet_grown The wet diameter of grown particles [nm]
  !> @param [in] rho_grown The mass density of grown particles [kg/m3]
  !> @param [in] cond_growth_rate The condensation growth rate GR [m/s]
  !> @param [in] cond_sink The condensation sink CS' [1/m2]
  function growth_parameter_from_rate(temp, d_dry_crit, d_wet_crit, &
                                      d_dry_grown, d_wet_grown, rho_grown, &
                                      cond_growth_rate, cond_sink) result(eta)

    real(wp), intent(in) :: temp, d_dry_crit, d_wet_crit, &
                            d_dry_grown, d_wet_grown, rho_grown, &
                            cond_growth_rate, cond_sink

    real(wp) :: gamma, eta

    ! Compute gamma from KK2002 eq 22 [nm2/m2/h], neglecting the
    ! (d_mean/150)^0.048 factor.
    gamma = 0.23_wp * d_wet_crit**0.2 * (d_wet_grown/3.0)**0.075 * &
            (1e-3*rho_grown)**-0.33 * (temp/293.0)**(-0.75)

    ! Compute eta [nm] using KK2002 eq 11.
    eta = gamma * cond_sink / cond_growth_rate
  end function

  !> Computes a conversion factor that transforms the "real" (base) nucleation
  !> rate @f$J@f$ to an "apparent" nucleation rate that accounts for the growth
  !> of nucleated particles in critical clusters (CC) needed to place them into
  !> an appropriate nucleation mode. The converstion from the "real" to
  !> "apparent" nucleation rate is
  !> @f$J_{app} = J_{real} \exp\left[\frac{\eta}{d_f} -
  !>                                 \frac{\eta}{d_i}\right],@f$
  !> where @f$\eta@f$ [nm]is a growth parameter computed by the parameterization,
  !> and @f$d_i@f$ and @f$d_f@f$ [nm] are the initial (nucleated) and final
  !> (grown) wet diameters of the particles in question.
  !> @param [in] c_so4 The number concentration of SO4 aerosol [#/cc]
  !> @param [in] c_nh4 The number concentration of NH4 aerosol [#/cc]
  !> @param [in] nh4_to_so4_molar_ratio The molar ratio of NH4 to SO4 [-]
  !> @param [in] temp The atmospheric temperature [K]
  !> @param [in] rel_hum The atmospheric relative humidity [-]
  !> @param [in] d_dry_crit The dry diameter of particles in a CC [nm]
  !> @param [in] d_wet_crit The wet diameter of particles in a CC [nm]
  !> @param [in] d_dry_grown The dry diameter of grown particles [nm]
  !> @param [in] rho_grown The mass density of grown particles [kg/m3]
  !> @param [in] rho_air The mass density of dry air [kg/m3]
  !> @param [in] mw_h2so4 The molecular weight of H2SO4 gas [kg/mol]
  function apparent_nucleation_factor(c_so4, c_nh4, nh4_to_so4_molar_ratio, &
                                      temp, rel_hum, d_dry_crit, &
                                      d_wet_crit, d_dry_grown, rho_grown, &
                                      rho_air, mw_h2so4) result(factor)

    real(wp), intent(in) :: c_so4, c_nh4, nh4_to_so4_molar_ratio, &
                            temp, rel_hum, d_dry_crit, &
                            d_wet_crit, d_dry_grown, rho_grown, &
                            rho_air, mw_h2so4

    real(wp) :: bounded_rel_hum, wet_dry_vol_ratio, d_wet_grown, eta
    real(wp) :: factor

    ! Compute the wet/dry volume ratio using the simple Kohler approximation
    ! for ammonium sulfate and bisulfate.
    bounded_rel_hum = max(0.10_wp, min(0.95_wp, rel_hum))
    wet_dry_vol_ratio = 1.0_wp - 0.56_wp / log(bounded_rel_hum)

    ! Wet diameter [nm] of grown particles with dry diameter d_dry_grown.
    d_wet_grown = 1e9_wp * d_dry_grown * wet_dry_vol_ratio**(1.0_wp / 3.0_wp)

    ! Growth parameter eta.
    eta = growth_parameter(c_so4, c_nh4, nh4_to_so4_molar_ratio, temp, &
                           rel_hum, d_dry_crit, d_wet_crit, d_dry_grown, &
                           rho_grown, rho_air, mw_h2so4)

    factor = apparent_nucleation_factor_from_eta(eta, d_wet_crit, d_wet_grown)
  end function


  !> Computes the conversion factor connecting the "real" (base) nucleation rate
  !> @f$J@f$ to the "apparent" nucleation rate using the given growth parameter
  !> @f$\eta@f$ and the initial and final wet particle diameters:
  !> @f$J_{app} = J_{real} \exp\left[\frac{\eta}{d_f} -
  !>                                 \frac{\eta}{d_i}\right],@f$
  !> where and @f$d_i@f$ and @f$d_f@f$ [nm] are the initial (nucleated) and final
  !> (grown) wet diameters of the particles in question.
  !> @param [in] eta The growth parameter @f$\eta@f$ [nm]
  !> @param [in] d_wet_crit The wet diameter of particles in a CC [nm]
  !> @param [in] d_wet_grown The wet diameter of grown particles [nm]
  function apparent_nucleation_factor_from_eta(eta, d_wet_crit, d_wet_grown) result(factor)
    real(wp), intent(in) :: eta, d_wet_crit, d_wet_grown

    real(wp) :: factor

    factor = exp(eta / d_wet_grown - eta / d_wet_crit)
  end function

end module
