!< \section arg_table_held_suarez_1994
!! \htmlinclude held_suarez_1994.html
module held_suarez_1994
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Implement idealized Held-Suarez forcings
  !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
  !    intercomparison of the dynamical cores of atmospheric general
  !    circulation models.'
  !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
  ! 
  !-----------------------------------------------------------------------

  use ccpp_kinds, only: kind_phys

  implicit none
  private
  save

  public :: held_suarez_1994_init
  public :: held_suarez_1994_run

  !!
  !! Forcing parameters
  !!
  real(kind_phys), parameter :: efoldf  =  1._kind_phys  ! efolding time for wind dissipation
  real(kind_phys), parameter :: efolda  = 40._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: efolds  =  4._kind_phys  ! efolding time for T dissipation
  real(kind_phys), parameter :: sigmab  =  0.7_kind_phys ! threshold sigma level
  real(kind_phys), parameter :: t00     = 200._kind_phys ! minimum reference temperature
  real(kind_phys), parameter :: kf      = 1._kind_phys/(86400._kind_phys*efoldf) ! 1./efolding_time for wind dissipation

  real(kind_phys), parameter :: onemsig = 1._kind_phys - sigmab ! 1. - sigma_reference

  real(kind_phys), parameter :: ka      = 1._kind_phys/(86400._kind_phys * efolda) ! 1./efolding_time for temperature diss.
  real(kind_phys), parameter :: ks      = 1._kind_phys/(86400._kind_phys * efolds)

  !!
  !! Model constants, reset in init call
  !!
  real(kind_phys)              :: cappa = 2.0_kind_phys / 7.0_kind_phys  ! R/Cp
  real(kind_phys)              :: cpair = 1004.0_kind_phys        ! specific heat of dry air (J/K/kg)
  real(kind_phys)              :: psurf_ref = 0.0_kind_phys       ! Surface pressure
  ! pref_mid_norm are layer midpoints normalized by surface pressure ('eta' coordinate)
  real(kind_phys), allocatable :: pref_mid_norm(:)



!======================================================================= 
contains
!======================================================================= 

!> \section arg_table_held_suarez_1994_init Argument Table
!! \htmlinclude held_suarez_1994_init.html
  subroutine held_suarez_1994_init(pver_in, cappa_in, cpair_in, psurf_ref_in, pref_mid_norm_in, errmsg, errflg)
    !! Dummy arguments
    integer,           intent(in) :: pver_in
    real(kind_phys),   intent(in) :: cappa_in
    real(kind_phys),   intent(in) :: cpair_in
    real(kind_phys),   intent(in) :: psurf_ref_in
    real(kind_phys),   intent(in) :: pref_mid_norm_in(:)
    character(len=512),intent(out):: errmsg
    integer,           intent(out):: errflg

    integer               :: pver                     ! Num vertical levels

    errmsg = ' '
    errflg = 0

    pver = pver_in
    allocate(pref_mid_norm(pver))
    cappa         = cappa_in
    cpair         = cpair_in
    psurf_ref     = psurf_ref_in
    pref_mid_norm = pref_mid_norm_in

  end subroutine held_suarez_1994_init

!> \section arg_table_held_suarez_1994_run Argument Table
!! \htmlinclude held_suarez_1994_run.html
  subroutine held_suarez_1994_run(pver, ncol, clat, pmid, &
       u, v, t, du, dv, s, errmsg, errflg)

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pver      ! Num vertical levels
    integer,  intent(in)  :: ncol      ! Num active columns
    real(kind_phys), intent(in)  :: clat(:)   ! latitudes(radians) for columns
    real(kind_phys), intent(in)  :: pmid(:,:) ! mid-point pressure
    real(kind_phys), intent(in)  :: u(:,:)    ! Zonal wind (m/s)
    real(kind_phys), intent(in)  :: v(:,:)    ! Meridional wind (m/s)
    real(kind_phys), intent(in)  :: t(:,:)    ! Temperature (K)
    !
    ! Output arguments
    !
    real(kind_phys),   intent(out) :: du(:,:)   ! Zonal wind tend
    real(kind_phys),   intent(out) :: dv(:,:)   ! Meridional wind tend
    real(kind_phys),   intent(out) :: s(:,:)    ! Heating rate
    character(len=512),intent(out):: errmsg
    integer,           intent(out):: errflg
    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(kind_phys) :: kv            ! 1./efolding_time (normalized) for wind
    real(kind_phys) :: kt            ! 1./efolding_time for temperature diss.
    real(kind_phys) :: trefa         ! "radiative equilibrium" T
    real(kind_phys) :: trefc         ! used in calc of "radiative equilibrium" T
    real(kind_phys) :: cossq(ncol)   ! coslat**2
    real(kind_phys) :: cossqsq(ncol) ! coslat**4
    real(kind_phys) :: sinsq(ncol)   ! sinlat**2
    real(kind_phys) :: coslat(ncol)  ! cosine(latitude)

    ! CJ add
    ! FV3-type Rayleigh friction
    !
    real(kind_phys)              :: pih                         ! 0.5*pi
    real(kind_phys)              :: kr                          ! RF friction coefficient
    ! real(kind_phys)            :: coeff                       ! coefficient for time-implicit scheme, needed for time-implicit scheme
    ! real(kind_phys)            :: u_tmp, v_tmp                ! updated u and v (implicit time stepping)
    ! real(kind_phys), parameter :: fv3_rf_cutoff = 0.00001_kind_phys  ! 1 Pa in sigma coordinate, Rayleigh friction above it
    real(kind_phys), parameter   :: fv3_rf_cutoff = 0.001_kind_phys    ! 100 Pa in sigma coordinate, Rayleigh friction above it
    ! real(kind_phys), parameter :: fv3_tau_rev   = 1._kind_phys/(86400._kind_phys*1._kind_phys)  ! 1/(1 day) damping time scale 
    real(kind_phys), parameter   :: fv3_tau_rev   = 1._kind_phys/(86400._kind_phys*3._kind_phys)  ! 1/(3 days) damping time scale 
    !
    !-----------------------------------------------------------------------
    !
    pih = 2._kind_phys*atan(1._kind_phys)      ! Initialize pi/2

    errmsg = ' '
    errflg = 0

    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)
    end do

    !
    !-----------------------------------------------------------------------
    !
    ! Held/Suarez IDEALIZED physics algorithm:
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the
    !   intercomparison of the dynamical cores of atmospheric general
    !   circulation models.
    !   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
    !
    !-----------------------------------------------------------------------
    !
    ! Compute idealized radiative heating rates (as dry static energy)
    !
    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        do i = 1, ncol
          kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
          trefc   = 315._kind_phys - (60._kind_phys * sinsq(i))
          trefa = (trefc - 10._kind_phys*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*kt*cpair
        end do
      else
        do i = 1, ncol
          trefc   = 315._kind_phys - 60._kind_phys*sinsq(i)
          trefa = (trefc - 10._kind_phys*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*ka*cpair
        end do
      end if
    end do
    !
    ! Add diffusion near the surface for the wind fields
    !
    du(:,:) = 0._kind_phys
    dv(:,:) = 0._kind_phys

    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
        do i = 1, ncol
          du(i,k) = -kv*u(i,k)
          dv(i,k) = -kv*v(i,k)
        end do
      end if
      ! CJ, add Rayleigh friction in sponge layer
      !     pref_mid_norm(1) serves as the position of the model top (full level, actual half level has lower pressure)
      if (pref_mid_norm(k) < fv3_rf_cutoff) then   ! RF above the cutoff sigma level
         kr = fv3_tau_rev * (sin(pih*log(fv3_rf_cutoff/pref_mid_norm(k))/log(fv3_rf_cutoff/pref_mid_norm(1))))**2._kind_phys  ! FV3 RF coefficient
      !  coeff = 1._kind_phys/(1._kind_phys + kr*ztodt)     ! requires the time step ztodt
         do i = 1, ncol
      !    u_tmp = u(i,k)*coeff               ! updated u via implicit time stepping
      !    v_tmp = v(i,k)*coeff               ! updated v via implicit time stepping 
      !    du(i,k) = (u_tmp - u(i,k))/ztodt   ! tendency
      !    dv(i,k) = (v_tmp - v(i,k))/ztodt   ! tendency 
           du(i,k) = -kr*u(i,k)               ! tendency via explicit time stepping
           dv(i,k) = -kr*v(i,k)               ! tendency via explicit time stepping
         end do
       endif
 
    end do

  end subroutine held_suarez_1994_run

end module held_suarez_1994
