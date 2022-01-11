module wang2008

  use haero_precision, only: wp

  implicit none

  private

  public :: first_order_pbl_nucleation_rate, &
            second_order_pbl_nucleation_rate

contains

  !> Computes the nucleation rate within the planetary boundary layer using a
  !> first-order reaction (Wang 2008 eq 1) adopted from the case studies in
  !> Shito et al (2006).
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  function first_order_pbl_nucleation_rate(c_h2so4) result(rate)
    real(wp), intent(in) :: c_h2so4

    real(wp) :: rate

    rate = 1e-6_wp * c_h2so4;
  end function

  !> Computes the nucleation rate within the planetary boundary layer using a
  !> second-order reaction (Wang 2008 eq 2) adopted from the case studies in
  !> Shito et al (2006).
  !> @param [in] c_h2so4 The number concentration of H2SO4 gas [cm-3]
  function second_order_pbl_nucleation_rate(c_h2so4) result(rate)
    real(wp), intent(in) :: c_h2so4

    real(wp) :: rate

    rate = 1e-12_wp * c_h2so4;
  end function
end module
