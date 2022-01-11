module vehkamaki2002

  use haero_precision, only: wp
  implicit none

  private

  public :: h2so4_critical_mole_fraction, &
            nucleation_rate, &
            num_critical_molecules, &
            critical_radius, &
            h2so4_nucleation_threshold

contains

  !> Computes the mole fraction of sulfuric acid in a critical cluster as
  !> parameterized by Vehkmaki et al (2002), eq 11.
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] rel_hum The relative humidity [-]
  function h2so4_critical_mole_fraction(c_h2so4, temp, rel_hum) result(x_crit)
    real(wp), intent(in) :: c_h2so4, temp, rel_hum

    real(wp) :: x_crit

    x_crit = 0.740997_wp - 0.00266379_wp * temp   &
             - 0.00349998_wp * log (c_h2so4)   &
             + 0.0000504022_wp * temp * log (c_h2so4)   &
             + 0.00201048_wp * log (rel_hum)   &
             - 0.000183289_wp * temp * log (rel_hum)   &
             + 0.00157407_wp * (log (rel_hum)) ** 2.0_wp   &
             - 0.0000179059_wp * temp * (log (rel_hum)) ** 2.0_wp   &
             + 0.000184403_wp * (log (rel_hum)) ** 3.0_wp   &
             - 1.50345e-6_wp * temp * (log (rel_hum)) ** 3.0_wp
  end function

  !> Computes the binary nucleation rate [m-3 s-1] as parameterized by
  !> Vehkmaki et al (2002), eq 12.
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] rel_hum The relative humidity [-]
  !> @param [in] x_crit The mole fraction of H2SO4 in a critical cluster [-]
  function nucleation_rate(c_h2so4, temp, rel_hum, x_crit) result(rate)

    real(wp), intent(in) :: c_h2so4, temp, rel_hum, x_crit

    real(wp) :: a, b, c, d, e, f, g, h, i, j
    real(wp) :: log_rate, rate

    a = 0.14309_wp+2.21956_wp*temp   &
        - 0.0273911_wp * temp**2.0_wp   &
        + 0.0000722811_wp * temp**3.0_wp + 5.91822_wp/x_crit

    b = 0.117489_wp + 0.462532_wp *temp   &
        - 0.0118059_wp * temp**2.0_wp   &
        + 0.0000404196_wp * temp**3.0_wp + 15.7963_wp/x_crit

    c = -0.215554_wp-0.0810269_wp * temp   &
        + 0.00143581_wp * temp**2.0_wp   &
        - 4.7758e-6_wp * temp**3.0_wp   &
        - 2.91297_wp/x_crit

    d = -3.58856_wp+0.049508_wp * temp   &
        - 0.00021382_wp * temp**2.0_wp   &
        + 3.10801e-7_wp * temp**3.0_wp   &
        - 0.0293333_wp/x_crit

    e = 1.14598_wp - 0.600796_wp * temp   &
        + 0.00864245_wp * temp**2.0_wp   &
        - 0.0000228947_wp * temp**3.0_wp   &
        - 8.44985_wp/x_crit

    f = 2.15855_wp + 0.0808121_wp * temp   &
        - 0.000407382_wp * temp**2.0_wp   &
        - 4.01957e-7_wp * temp**3.0_wp   &
        + 0.721326_wp/x_crit

    g = 1.6241_wp - 0.0160106_wp * temp   &
        + 0.0000377124_wp * temp**2.0_wp   &
        + 3.21794e-8_wp * temp**3.0_wp   &
        - 0.0113255_wp/x_crit

    h = 9.71682_wp - 0.115048_wp * temp   &
        + 0.000157098_wp * temp**2.0_wp   &
        + 4.00914e-7_wp * temp**3.0_wp   &
        + 0.71186_wp/x_crit

    i = -1.05611_wp + 0.00903378_wp * temp   &
        - 0.0000198417_wp * temp**2.0_wp   &
        + 2.46048e-8_wp  * temp**3.0_wp   &
        - 0.0579087_wp/x_crit

    j = -0.148712_wp + 0.00283508_wp * temp   &
        - 9.24619e-6_wp  * temp**2.0_wp   &
        + 5.00427e-9_wp * temp**3.0_wp   &
        - 0.0127081_wp/x_crit

    log_rate = a &
               + b * log (rel_hum) &
               + c * ( log (rel_hum))**2.0_wp &
               + d * ( log (rel_hum))**3.0_wp &
               + e * log (c_h2so4) &
               + f * (log (rel_hum)) * (log (c_h2so4)) &
               + g * ((log (rel_hum) ) **2.0_wp) * (log (c_h2so4)) &
               + h * (log (c_h2so4)) **2.0_wp &
               + i * log (rel_hum) * ((log (c_h2so4)) **2.0_wp) &
               + j * (log (c_h2so4)) **3.0_wp
    log_rate = min(log_rate, log(1.0e38_wp))
    rate = exp(log_rate)
  end function

  !> Computes the total number of molecules in a critical cluster as
  !> parameterized in Vehkamaki et al (2002), eq 13.
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] rel_hum The relative humidity [-]
  !> @param [in] x_crit The mole fraction of H2SO4 in a critical cluster [-]
  function num_critical_molecules(c_h2so4, temp, rel_hum, x_crit) result(n_crit)

    real(wp), intent(in) :: c_h2so4, temp, rel_hum, x_crit

    real(wp) :: a, b, c, d, e, f, g, h, i, j
    real(wp) :: n_crit

    a = -0.00295413_wp - 0.0976834_wp*temp   &
        + 0.00102485_wp * temp**2.0_wp   &
        - 2.18646e-6_wp * temp**3.0_wp - 0.101717_wp/x_crit

    b = -0.00205064_wp - 0.00758504_wp*temp   &
        + 0.000192654_wp * temp**2.0_wp   &
        - 6.7043e-7_wp * temp**3.0_wp - 0.255774_wp/x_crit

    c = +0.00322308_wp + 0.000852637_wp * temp   &
        - 0.0000154757_wp * temp**2.0_wp   &
        + 5.66661e-8_wp * temp**3.0_wp   &
        + 0.0338444_wp/x_crit

    d = +0.0474323_wp - 0.000625104_wp * temp   &
        + 2.65066e-6_wp * temp**2.0_wp   &
        - 3.67471e-9_wp * temp**3.0_wp   &
        - 0.000267251_wp/x_crit

    e = -0.0125211_wp + 0.00580655_wp * temp   &
        - 0.000101674_wp * temp**2.0_wp   &
        + 2.88195e-7_wp * temp**3.0_wp   &
        + 0.0942243_wp/x_crit

    f = -0.038546_wp - 0.000672316_wp * temp   &
        + 2.60288e-6_wp * temp**2.0_wp   &
        + 1.19416e-8_wp * temp**3.0_wp   &
        - 0.00851515_wp/x_crit

    g = -0.0183749_wp + 0.000172072_wp * temp   &
        - 3.71766e-7_wp * temp**2.0_wp   &
        - 5.14875e-10_wp * temp**3.0_wp   &
        + 0.00026866_wp/x_crit

    h = -0.0619974_wp + 0.000906958_wp * temp   &
        - 9.11728e-7_wp * temp**2.0_wp   &
        - 5.36796e-9_wp * temp**3.0_wp   &
        - 0.00774234_wp/x_crit

    i = +0.0121827_wp - 0.00010665_wp * temp   &
        + 2.5346e-7_wp * temp**2.0_wp   &
        - 3.63519e-10_wp * temp**3.0_wp   &
        + 0.000610065_wp/x_crit

    j = +0.000320184_wp - 0.0000174762_wp * temp   &
        + 6.06504e-8_wp * temp**2.0_wp   &
        - 1.4177e-11_wp * temp**3.0_wp   &
        + 0.000135751_wp/x_crit

    n_crit = exp ( &
      a &
      + b * log (rel_hum)   &
      + c * ( log (rel_hum))**2.0_wp   &
      + d * ( log (rel_hum))**3.0_wp   &
      + e * log (c_h2so4)   &
      + f * (log (rel_hum)) * (log (c_h2so4))   &
      + g * ((log (rel_hum) ) **2.0_wp) * (log (c_h2so4))   &
      + h * (log (c_h2so4)) **2.0_wp   &
      + i * log (rel_hum) * ((log (c_h2so4)) **2.0_wp)   &
      + j * (log (c_h2so4)) **3.0_wp)

  end function

  !> Computes the radius [nm] of a critical cluster as parameterized in Vehkamaki
  !> et al (2002), eq 14.
  !> @param [in] x_crit The mole fraction of H2SO4 in a critical cluster [-]
  !> @param [in] n_tot The total number of molecules in the critical cluster [-]
  function critical_radius(x_crit, n_tot) result(r_crit)

    real(wp), intent(in) :: x_crit, n_tot

    real(wp) :: r_crit

    r_crit = exp( -1.6524245_wp + 0.42316402_wp*x_crit &
                 + 0.3346648_wp*log(n_tot) )
  end function

  !> Computes the threshold number concentration of H2SO4 [cm-3] that produces a
  !> nucleation rate of 1 cm-3 s-1 at the given temperature and relative
  !> humidity as parameterized by Vehkamaki et al (2002), eq 15.
  !> @param [in] temp The atmospherіc temperature [K]
  !> @param [in] rel_hum The relative humidity [-]
  function h2so4_nucleation_threshold(temp, rel_hum) result(thresh)

    real(wp), intent(in) :: temp, rel_hum

    real(wp) :: thresh

    thresh = exp(-279.243_wp &
      + 11.7344_wp * rel_hum &
      + 22700.9 / temp &
      - 1088.64 * rel_hum / temp &
      + 1.14436 * temp &
      - 0.0302331 * rel_hum * temp &
      - 0.00130254 * temp**2 &
      - 6.38697 * log(rel_hum) &
      + 854.98 * log(rel_hum) / temp &
      + 0.00879662 * temp * log(rel_hum))
  end function

end module
