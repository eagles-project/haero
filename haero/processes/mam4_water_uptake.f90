module mam4_wateruptake

  use iso_c_binding
  use haero, only: wp => Real, model_num_modes

  implicit none
  private

  public :: wateruptake_init, &
            wateruptake_init_with_bisection, &
            wateruptake_update, &
            wateruptake_destroy

  !> This type stores parameters for the water uptake process.
  type :: waterupdate_process
    !> Do we use bisection, or the algebraic method?
    logical :: use_bisection
    !> How many modes are in this aerosol model?
    integer :: nmodes
  end type

  !> These parameters are used by the water uptake process.
  real(wp), parameter :: third  = 1._wp/3._wp
  real(wp), parameter :: pi43   = pi*4.0_wp/3.0_wp

contains

!> Creates and returns an initialized pointer for a water uptake process that
!> solves the Kohler equation using the algebraic method.
type(c_ptr) function wateruptake_init(model) bind(c)
  use iso_c_binding, only :: c_ptr, c_loc
  implicit none

  ! Arguments
  type(c_ptr),            intent(in) :: model

  type(water_uptake_process), pointer, target :: process

  ! Initialize the water uptake process.
  allocate(process)
  process%use_bisection = .false.
  process%nmodes = model_num_modes(model)

  ! Return a C pointer.
  wateruptake_init = c_loc(process)
end function

!> Creates and returns an initialized pointer for a water uptake process that
!> solves the Kohler equation using bisection.
type(c_ptr) function wateruptake_init_with_bisection(model) bind(c)
  use iso_c_binding, only :: c_ptr, c_loc
  implicit none

  ! Arguments
  type(c_ptr),            intent(in) :: model

  type(water_uptake_process), pointer, target :: process

  ! Initialize the water uptake process.
  allocate(process)
  process%use_bisection = .true.
  process%nmodes = model_num_modes(model)

  ! Return a C pointer.
  wateruptake_init = c_loc(process)
end function

subroutine wateruptake_destroy(p)
  use iso_c_binding, only :: c_ptr, c_f_pointer
  implicit none

  type(c_ptr), intent(inout) :: p

  type(wateruptake_process), pointer :: process

  call c_f_pointer(p, process)
  deallocate(process)
end subroutine

subroutine wateruptake_update(p, model, t, prognostics, diagnostics) bind(c)
  use iso_c_binding, only :: c_ptr, c_f_pointer
  implicit none

  ! Arguments
  type(c_ptr),       intent(in)    :: p
  type(c_ptr),       intent(in)    :: model
  real(wp), value,   intent(in)    :: t
  type(c_ptr),       intent(in)    :: prognostics
  type(c_ptr),       intent(inout) :: diagnostics

  real(wp), optional, pointer :: dgnumdry_m(:,:,:)
  real(wp), optional, pointer :: dgnumwet_m(:,:,:)
  real(wp), optional, pointer :: qaerwat_m(:,:,:)
  real(wp), optional, pointer :: wetdens_m(:,:,:)

  ! local variables

  type(wateruptake_process) :: process

  integer :: lchnk                         ! chunk index
  integer :: ncol                          ! number of columns
  integer :: list_idx                      ! radiative constituents list index
  integer :: stat

  integer :: i, k, l, m
  integer :: itim_old
  integer :: nmodes
  integer :: nspec

  real(wp), pointer :: h2ommr(:,:)          ! specific humidity
  real(wp), pointer :: t(:,:)               ! temperatures (K)
  real(wp), pointer :: pmid(:,:)            ! layer pressure (Pa)
  real(wp), pointer :: raer(:,:)            ! aerosol species MRs (kg/kg and #/kg)

  real(wp), pointer :: cldn(:,:)            ! layer cloud fraction (0-1)
  real(wp), pointer :: dgncur_a(:,:,:)
  real(wp), pointer :: dgncur_awet(:,:,:)
  real(wp), pointer :: wetdens(:,:,:)
  real(wp), pointer :: qaerwat(:,:,:)

  real(wp), allocatable :: maer(:,:,:)      ! aerosol wet mass MR (including water) (kg/kg-air)
  real(wp), allocatable :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
  real(wp), allocatable :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
  real(wp), allocatable :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
  real(wp), allocatable :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)
  real(wp), allocatable :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)

  real(wp), allocatable :: rhcrystal(:)
  real(wp), allocatable :: rhdeliques(:)
  real(wp), allocatable :: specdens_1(:)

  real(wp) :: dryvolmr(pcols,pver)          ! volume MR for aerosol mode (m3/kg)
  real(wp) :: specdens
  real(wp) :: spechygro, spechygro_1
  real(wp) :: duma, dumb
  real(wp) :: sigmag
  real(wp) :: alnsg
  real(wp) :: v2ncur_a
  real(wp) :: drydens                       ! dry particle density  (kg/m^3)
  real(wp) :: rh(pcols,pver)                ! relative humidity (0-1)

  real(wp) :: es(pcols)                     ! saturation vapor pressure
  real(wp) :: qs(pcols)                     ! saturation specific humidity
  real(wp) :: cldn_thresh
  real(wp) :: aerosol_water(pcols,pver)     ! sum of aerosol water (wat_a1 + wat_a2 + wat_a3 + wat_a4)
  logical :: history_aerosol                ! Output the MAM aerosol variables and tendencies
  logical :: history_verbose                ! produce verbose history output

  character(len=3) :: trnum       ! used to hold mode number (as characters)

  ! Get access to process data.
  call c_f_pointer(p, process)

  lchnk = state%lchnk
  ncol = state%ncol

  ! determine default variables
  call phys_getopts(history_aerosol_out = history_aerosol, &
                    history_verbose_out = history_verbose)

  list_idx = 0
  if (present(list_idx_in)) then
    list_idx = list_idx_in

    ! check that all optional args are present
    if (.not. present(dgnumdry_m) .or. .not. present(dgnumwet_m) .or. &
        .not. present(qaerwat_m)  .or. .not. present(wetdens_m)) then
      call endrun('modal_aero_wateruptake_dr called for'// &
                  'diagnostic list but required args not present')
    end if

    ! arrays for diagnostic calculations must be associated
    if (.not. associated(dgnumdry_m) .or. .not. associated(dgnumwet_m) .or. &
        .not. associated(qaerwat_m)  .or. .not. associated(wetdens_m)) then
      call endrun('modal_aero_wateruptake_dr called for'// &
                  'diagnostic list but required args not associated')
    end if
  end if

  ! loop over all aerosol modes
  call rad_cnst_get_info(list_idx, nmodes=nmodes)

  allocate( &
    maer(pcols,pver,nmodes),     &
    hygro(pcols,pver,nmodes),    &
    naer(pcols,pver,nmodes),     &
    dryvol(pcols,pver,nmodes),   &
    drymass(pcols,pver,nmodes),  &
    dryrad(pcols,pver,nmodes),   &
    rhcrystal(nmodes),           &
    rhdeliques(nmodes),          &
    specdens_1(nmodes)           )

  maer(:,:,:)     = 0._wp
  hygro(:,:,:)    = 0._wp

  if (list_idx == 0) then
    call pbuf_get_field(pbuf, dgnum_idx,      dgncur_a )
    call pbuf_get_field(pbuf, dgnumwet_idx,   dgncur_awet )
    call pbuf_get_field(pbuf, wetdens_ap_idx, wetdens)
    call pbuf_get_field(pbuf, qaerwat_idx,    qaerwat)
  else
    dgncur_a    => dgnumdry_m
    dgncur_awet => dgnumwet_m
    qaerwat     => qaerwat_m
    wetdens     => wetdens_m
  end if

 ! Step1: prepare the necessary input data (at chunk level) to solve Kohler equation
 !        the necessary input includes 1) clear-sky relative humidty; 2) volume-mean hygroscopicity; 3) dry radius/volume/mass
 do m = 1, nmodes

   dryvolmr(:,:) = 0._wp

   ! get mode properties: sigmag and critical RH for one mode
   call rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag,  &
     rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))

   ! get mode info: number of species in a mode
   call rad_cnst_get_info(list_idx, m, nspec=nspec)

   do l = 1, nspec

     ! get species interstitial mixing ratio ('a')
     call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', state, pbuf, raer)                       ! get species mass mixing ratio
     call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, hygro_aer=spechygro)  ! get species density and hygroscopicity

     if (l == 1) then
       ! save off these values to be used as defaults
       specdens_1(m)  = specdens
       spechygro_1    = spechygro
     end if

     do k = top_lev, pver
       do i = 1, ncol
         duma          = raer(i,k)
         maer(i,k,m)   = maer(i,k,m) + duma
         dumb          = duma/specdens
         dryvolmr(i,k) = dryvolmr(i,k) + dumb
         hygro(i,k,m)  = hygro(i,k,m) + dumb*spechygro
       end do
     end do
   end do

   alnsg = log(sigmag)

   do k = top_lev, pver
     do i = 1, ncol

       if (dryvolmr(i,k) > 1.0e-30_wp) then
         hygro(i,k,m) = hygro(i,k,m)/dryvolmr(i,k)     ! volume-mean hygroscopicity
       else
         hygro(i,k,m) = spechygro_1
       end if

       ! dry aerosol properties

       v2ncur_a = 1._wp / ( (pi/6._wp)*(dgncur_a(i,k,m)**3._wp)*exp(4.5_wp*alnsg**2._wp) )
       ! naer = aerosol number (#/kg)
       naer(i,k,m) = dryvolmr(i,k)*v2ncur_a

       ! compute mean (1 particle) dry volume and mass for each mode
       ! old coding is replaced because the new (1/v2ncur_a) is equal to
       ! the mean particle volume
       ! also moletomass forces maer >= 1.0e-30, so (maer/dryvolmr)
       ! should never cause problems (but check for maer < 1.0e-31 anyway)
       if (maer(i,k,m)  >  1.0e-31_wp) then
         drydens = maer(i,k,m)/dryvolmr(i,k)
       else
         drydens = 1.0_wp
       end if
       dryvol(i,k,m)   = 1.0_wp/v2ncur_a
       drymass(i,k,m)  = drydens*dryvol(i,k,m)
       dryrad(i,k,m)   = (dryvol(i,k,m)/pi43)**third

     end do
   end do
 end do ! modes

 ! relative humidity calc
 h2ommr => state%q(:,:,1)
 t      => state%t
 pmid   => state%pmid

 itim_old    =  pbuf_old_tim_idx()
 call pbuf_get_field(pbuf, cld_idx, cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

 do k = top_lev, pver
   call qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol))
   do i = 1, ncol
     if (qs(i) > h2ommr(i,k)) then
       rh(i,k) = h2ommr(i,k)/qs(i)
     else
       rh(i,k) = 0.98_wp
     endif
     rh(i,k) = max(rh(i,k), 0.0_wp)
     rh(i,k) = min(rh(i,k), 0.98_wp)
     if(pergro_mods) then
       cldn_thresh = 0.9998_wp
     else
       cldn_thresh = 1.0_wp !original code
     endif
     if (cldn(i,k)  <  cldn_thresh) then !BSINGH
       rh(i,k) = (rh(i,k) - cldn(i,k)) / (1.0_wp - cldn(i,k))  ! clear portion
     end if
     rh(i,k) = max(rh(i,k), 0.0_wp)
   end do
 end do

 ! Step2: solve Kohler equation with the array-based input generated above
 !        update dgn_wet, qaerwat and wetdens
 call wateruptake_sub( ncol, top_lev, pver, nmodes, rhcrystal, rhdeliques,   &       ! intent(in)
                       dryrad, naer, hygro, rh, dryvol, drymass, specdens_1, &       ! intent(in)
                       dgncur_a, dgncur_awet, qaerwat, wetdens )                     ! pointer, intent(inout)

 if (list_idx == 0) then

   aerosol_water(:ncol,:) = 0._wp
   do m = 1, nmodes
     ! output to history
     write( trnum, '(i3.3)' ) m
     call outfld( 'wat_a'//trnum(3:3),  qaerwat(:,:,m),     pcols, lchnk)
     call outfld( 'dgnd_a'//trnum(2:3), dgncur_a(:,:,m),    pcols, lchnk)
     call outfld( 'dgnw_a'//trnum(2:3), dgncur_awet(:,:,m), pcols, lchnk)
     if (history_aerosol .and. .not. history_verbose) &
       aerosol_water(:ncol,:) = aerosol_water(:ncol,:) + qaerwat(:ncol,:,m)
   end do

   if (history_aerosol .and. .not. history_verbose) &
     call outfld( 'aero_water',  aerosol_water(:ncol,:),    ncol, lchnk)

 end if

 deallocate( &
   maer,   hygro,     naer,       dryvol,    drymass,    &
   dryrad, rhcrystal, rhdeliques, specdens_1             )

end subroutine waterupdate_update

subroutine wateruptake_sub( &
   ncol, str_lev, end_lev, nmodes, rhcrystal, rhdeliques, dryrad, naer, hygro,
   rh, dryvol, drymass, specdens_1, dgncur_a, dgncur_awet, qaerwat, wetdens)

  ! Arguments
  integer, intent(in) :: ncol             ! number of columns
  integer, intent(in) :: str_lev, end_lev ! start and end vertical levels
  integer, intent(in) :: nmodes

  real(wp), dimension(:), intent(in)      :: rhcrystal, rhdeliques, specdens_1
  real(wp), dimension(:,:), intent(in)    :: rh
  real(wp), dimension(:,:,:), intent(in)  :: naer, dryvol, drymass, dryrad, hygro

  real(wp), pointer :: dgncur_a(:,:,:), dgncur_awet(:,:,:), qaerwat(:,:,:), wetdens(:,:,:)

  ! local variables

  integer :: i, k, m
  real(wp), dimension(:,:,:), allocatable :: wetrad, &  ! wet radius of aerosol (m)
                                             wetvol, &  ! single-particle-mean wet volume (m3)
                                             wtrvol     ! single-particle-mean water volume in wet aerosol (m3)

  real(wp) :: hystfac                                   ! working variable for hysteresis

  ! alocate local arrays
  allocate(wetrad(pcols,pver,nmodes), wetvol(pcols,pver,nmodes), wtrvol(pcols,pver,nmodes))

  ! loop over all aerosol modes
  do m = 1, nmodes

    hystfac = 1.0_wp / max(1.0e-5_wp, (rhdeliques(m) - rhcrystal(m)))

    do k = str_lev, end_lev            ! make this change so that the foloowing subroutine can be callable for a grid cell
      do i = 1, ncol

        ! compute wet radius for each mode
        call kohler(use_bisection, dryrad(i:i,k,m), &
                    hygro(i:i,k,m), rh(i:i,k), &
                    wetrad(i:i,k,m), 1)

        wetrad(i,k,m) = max(wetrad(i,k,m), dryrad(i,k,m))
        wetvol(i,k,m) = pi43*wetrad(i,k,m)**3
        wetvol(i,k,m) = max(wetvol(i,k,m), dryvol(i,k,m))
        wtrvol(i,k,m) = wetvol(i,k,m) - dryvol(i,k,m)
        wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_wp)

        ! apply simple treatment of deliquesence/crystallization hysteresis
        ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
        ! the "upper curve" value, and the fraction is a linear function of rh
        if (rh(i,k) < rhcrystal(m)) then
          wetrad(i,k,m) = dryrad(i,k,m)
          wetvol(i,k,m) = dryvol(i,k,m)
          wtrvol(i,k,m) = 0.0_wp
        else if (rh(i,k) < rhdeliques(m)) then
          wtrvol(i,k,m) = wtrvol(i,k,m)*hystfac*(rh(i,k) - rhcrystal(m))
          wtrvol(i,k,m) = max(wtrvol(i,k,m), 0.0_wp)
          wetvol(i,k,m) = dryvol(i,k,m) + wtrvol(i,k,m)
          wetrad(i,k,m) = (wetvol(i,k,m)/pi43)**third
        end if

        ! calculate new dgn_awet, qaerwat, wetdens
        dgncur_awet(i,k,m) = dgncur_a(i,k,m) * (wetrad(i,k,m)/dryrad(i,k,m))
        qaerwat(i,k,m)     = rhoh2o*naer(i,k,m)*wtrvol(i,k,m)
        ! compute aerosol wet density (kg/m3)
        if (wetvol(i,k,m) > 1.0e-30_wp) then
          wetdens(i,k,m) = (drymass(i,k,m) + rhoh2o*wtrvol(i,k,m))/wetvol(i,k,m)
        else
          wetdens(i,k,m) = specdens_1(m)
        end if

      end do ! columns
    end do ! levels
  end do ! modes

  deallocate( wetrad, wetvol, wtrvol )

end subroutine wateruptake_sub

!-----------------------------------------------------------------------
subroutine kohler(use_bisection, rdry_in, hygro, s, rwet_out, im )

  ! calculates equlibrium radius r of haze droplets as function of
  ! dry particle mass and relative humidity s using kohler solution
  ! given in pruppacher and klett (eqn 6-35)

  ! for multiple aerosol types, assumes an internal mixture of aerosols

  implicit none

  ! arguments
  logical, intent(in) :: use_bisection ! use bisection to solve
  real(wp) :: rdry_in(:)    ! aerosol dry radius (m)
  real(wp) :: hygro(:)      ! aerosol volume-mean hygroscopicity (--)
  real(wp) :: s(:)          ! relative humidity (1 = saturated)
  real(wp) :: rwet_out(:)   ! aerosol wet radius (m)
  integer  :: im            ! number of grid points to be processed

  ! local variables
  integer, parameter :: imax=200
  integer :: i, n, nsol

  real(wp) :: a, b
  real(wp) :: p40(imax),p41(imax),p42(imax),p43(imax) ! coefficients of polynomial
  real(wp) :: p30(imax),p31(imax),p32(imax) ! coefficients of polynomial
  real(wp) :: p
  real(wp) :: r3, r4
  real(wp) :: r(im)         ! wet radius (microns)
  real(wp) :: rdry(imax)    ! radius of dry particle (microns)
  real(wp) :: ss            ! relative humidity (1 = saturated)
  real(wp) :: slog(imax)    ! log relative humidity
  real(wp) :: vol(imax)     ! total volume of particle (microns**3)
  real(wp) :: xi, xr

  complex(wp) :: cx4(4,imax),cx3(3,imax)

  real(wp), parameter :: eps     = 1.e-4_wp
  ! real(wp), parameter :: mw      = mwh2o              ! original value: 18._wp
  ! real(wp), parameter :: pi      = 3.14159_wp        ! use the pi in the physconst
  ! real(wp), parameter :: rhow    = rhoh2o / 1.e3_wp   ! original value: 1._wp
  ! real(wp), parameter :: surften = 76._wp
  ! real(wp), parameter :: tair    = 273._wp
  ! real(wp), parameter :: third   = 1._wp/3._wp
  ! real(wp), parameter :: ugascon = r_universal*1.e4   ! 8.3e7_wp

  real(wp), parameter :: mw      = 18._wp
  real(wp), parameter :: pi      = 3.14159_wp
  real(wp), parameter :: rhow    = 1._wp
  real(wp), parameter :: surften = 76._wp
  real(wp), parameter :: tair    = 273._wp
  real(wp), parameter :: third   = 1._wp/3._wp
  real(wp), parameter :: ugascon = 8.3e7_wp
  real(wp)            :: q4_coeff(5)      ! coefficient for quartic equation: q45*x^4 + q44*x^3 + q43*x^2 + q42*x + q41 = 0

  ! effect of organics on surface tension is neglected
  a=2.e4_wp*mw*surften/(ugascon*tair*rhow)

  do i=1,im
    rdry(i) = rdry_in(i)*1.0e6_wp     ! convert (m) to (microns)
    vol(i) = rdry(i)**3               ! vol is r**3, not volume
    b = vol(i)*hygro(i)

    ! quartic
    ss=min(s(i),1._wp-eps)
    ss=max(ss,1.e-10_wp)
    slog(i)=log(ss)
    p43(i)=-a/slog(i)
    p42(i)=0._wp
    p41(i)=b/slog(i)-vol(i)
    p40(i)=a*vol(i)/slog(i)

    ! cubic for rh=1
    p32(i)=0._wp
    p31(i)=-b/a
    p30(i)=-vol(i)
  end do


  do 100 i=1,im
    if (vol(i).le.1.e-12_wp) then
      r(i)=rdry(i)
      go to 100
    endif

    p=abs(p31(i))/(rdry(i)*rdry(i))
    if (p < eps) then
      ! approximate solution for small particles
      r(i)=rdry(i)*(1._wp+p*third/(1._wp-slog(i)*rdry(i)/a))
    else
      if ( use_bisection ) then
        q4_coeff(1) = p40(i)
        q4_coeff(2) = p41(i)
        q4_coeff(3) = p42(i)
        q4_coeff(4) = p43(i)
        q4_coeff(5) = 1.0_wp
        call bisection ( rdry(i), q4_coeff(1:5), r(i) )
      else
        call makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1)
        ! find smallest real(wp) solution
        r(i)=1000._wp*rdry(i)
        nsol=0
        do n=1,4
          xr=real(cx4(n,i))
          xi=aimag(cx4(n,i))
          if (abs(xi) > abs(xr)*eps) cycle
          if (xr > r(i)) cycle
          if (xr < rdry(i)*(1._wp-eps)) cycle
          if (xr.ne.xr) cycle
          r(i)=xr
          nsol=n
        end do
        if(nsol == 0)then
          write(iulog,*) 'ccm kohlerc - no real(wp) solution found (quartic)'
          write(iulog,*)'roots =', (cx4(n,i),n=1,4)
          write(iulog,*)'p0-p3 =', p40(i), p41(i), p42(i), p43(i)
          write(iulog,*)'rh=',s(i)
          write(iulog,*)'setting radius to dry radius=',rdry(i)
          r(i)=rdry(i)
          ! stop
        endif
      endif
    endif

    if (s(i) > 1._wp-eps) then
      ! save quartic solution at s=1-eps
      r4=r(i)
      ! cubic for rh=1
      p=abs(p31(i))/(rdry(i)*rdry(i))

      if (p < eps) then
        r(i)=rdry(i)*(1._wp+p*third)
      else
        if ( use_bisection ) then
          q4_coeff(1) = p30(i)
          q4_coeff(2) = p31(i)
          q4_coeff(3) = p32(i)
          q4_coeff(4) = 1.0_wp
          q4_coeff(5) = 0.0_wp
          call bisection ( rdry(i), q4_coeff(1:5), r(i) )
        else
          call makoh_cubic(cx3,p32,p31,p30,im)
          ! find smallest real(wp) solution
          r(i)=1000._wp*rdry(i)
          nsol=0
          do n=1,3
            xr=real(cx3(n,i))
            xi=aimag(cx3(n,i))
            if(abs(xi) > abs(xr)*eps) cycle
            if(xr > r(i)) cycle
            if(xr < rdry(i)*(1._wp-eps)) cycle
            if (xr.ne.xr) cycle
            r(i)=xr
            nsol=n
          end do
          if (nsol == 0) then
            write (iulog,*) 'ccm kohlerc - no real(wp) solution found (cubic)'
            write (iulog,*) 'roots =', (cx3(n,i),n=1,3)
            write (iulog,*) 'p0-p2 =', p30(i), p31(i), p32(i)
            write (iulog,*) 'rh=',s(i)
            write (iulog,*) 'setting radius to dry radius=',rdry(i)
            r(i)=rdry(i)
            ! stop
          endif
        endif
      endif
      r3=r(i)
      ! now interpolate between quartic, cubic solutions
      r(i)=(r4*(1._wp-s(i))+r3*(s(i)-1._wp+eps))/eps
    endif

  100 continue

  ! bound and convert from microns to m
  do i=1,im
    r(i) = min(r(i),30._wp) ! upper bound based on 1 day lifetime
    rwet_out(i) = r(i)*1.e-6_wp
  end do
end subroutine kohler

subroutine makoh_cubic( cx, p2, p1, p0, im )
  ! solves  x**3 + p2 x**2 + p1 x + p0 = 0
  ! where p0, p1, p2 are real
  integer, parameter :: imx=200
  integer :: im
  real(wp) :: p0(imx), p1(imx), p2(imx)
  complex(wp) :: cx(3,imx)

  integer :: i
  real(wp) :: eps, q(imx), r(imx), sqrt3, third
  complex(wp) :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

  save eps
  data eps/1.e-20_wp/

  third=1._wp/3._wp
  ci=cmplx(0._wp,1._wp,dp)
  sqrt3=sqrt(3._wp)
  cw=0.5_wp*(-1+ci*sqrt3)
  cwsq=0.5_wp*(-1-ci*sqrt3)

  do i=1,im
    if(p1(i) == 0._wp)then
      ! completely insoluble particle
      cx(1,i)=(-p0(i))**third
      cx(2,i)=cx(1,i)
      cx(3,i)=cx(1,i)
    else
      q(i)=p1(i)/3._wp
      r(i)=p0(i)/2._wp
      crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
      crad(i)=sqrt(crad(i))

      cy(i)=r(i)-crad(i)
      if (abs(cy(i)) > eps) then
        cy(i)=cy(i)**third
      endif
      cq=q(i)
      cz(i)=-cq/cy(i)

      cx(1,i)=-cy(i)-cz(i)
      cx(2,i)=-cw*cy(i)-cwsq*cz(i)
      cx(3,i)=-cwsq*cy(i)-cw*cz(i)
    endif
  enddo
end subroutine makoh_cubic

subroutine makoh_quartic( cx, p3, p2, p1, p0, im )

  ! solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
  ! where p0, p1, p2, p3 are real
  integer, parameter :: imx=200
  integer :: im
  real(wp) :: p0(imx), p1(imx), p2(imx), p3(imx)
  complex(wp) :: cx(4,imx)

  integer :: i
  real(wp) :: third, q(imx), r(imx)
  complex(wp) :: cb(imx), cb0(imx), cb1(imx), crad(imx), cy(imx), czero


  czero=cmplx(0.0_wp,0.0_wp,wp)
  third=1._wp/3._wp

  do 10 i=1,im

    q(i)=-p2(i)*p2(i)/36._wp+(p3(i)*p1(i)-4*p0(i))/12._wp
    r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48._wp   &
     +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16

    crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
    crad(i)=sqrt(crad(i))

    cb(i)=r(i)-crad(i)
    if(cb(i) == czero)then
      ! insoluble particle
      cx(1,i)=(-p1(i))**third
      cx(2,i)=cx(1,i)
      cx(3,i)=cx(1,i)
      cx(4,i)=cx(1,i)
    else
      cb(i)=cb(i)**third

      cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6

      cb0(i)=sqrt(cy(i)*cy(i)-p0(i))
      cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i))

      cb(i)=p3(i)/2+cb1(i)
      crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i))
      crad(i)=sqrt(crad(i))
      cx(1,i)=(-cb(i)+crad(i))/2._wp
      cx(2,i)=(-cb(i)-crad(i))/2._wp

      cb(i)=p3(i)/2-cb1(i)
      crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i))
      crad(i)=sqrt(crad(i))
      cx(3,i)=(-cb(i)+crad(i))/2._wp
      cx(4,i)=(-cb(i)-crad(i))/2._wp
    endif
  10 continue
end subroutine makoh_quartic

! bisection method to find a root for the equation F(x) = 0 !
subroutine bisection( rd, c, rw )

  real(wp), intent(in) :: rd, c(5)
  real(wp), intent(out) :: rw

  ! Local variables
  real(wp) :: left, right, middle, f_left, f_right, f_middle
  integer  :: niter
  integer, parameter  :: niter_max = 1000
  real(wp), parameter :: toler = 10._wp * epsilon(1._wp) ! bisection tolerance

  left    = rd
  right   = rd * 1e3_wp
  f_left  = quartic_func(c(1:5), left)
  f_right = quartic_func(c(1:5), right)

  ! start the iteration
  niter = 0
  do
    middle   = (left + right) / 2.0_wp
    f_middle = quartic_func(c(1:5), middle)

    if (abs(left-right)/middle < toler) then
      rw = middle
      exit
    else if (f_middle  ==  0._wp) then
      rw = middle
      exit
    else if (f_middle*f_left < 0.0_wp) then
      right   = middle
      f_right = f_middle
    else
      left    = middle
      f_left  = f_middle
    end if
    niter = niter + 1
    if (niter > niter_max) then
      write(*, '(a)') "Too many iterations..."
      rw = rd
      exit
    end if
  end do
end subroutine bisection

real(wp) function quartic_func( c, x )
  implicit none
  real(wp), intent(in) :: c(5)
  real(wp), intent(in) :: x

  quartic_func = c(5)*(x**4) + c(4)*(x**3) + c(3)*(x**2) + c(2)*x + c(1)

end function quartic_func


end module mam4_wateruptake

