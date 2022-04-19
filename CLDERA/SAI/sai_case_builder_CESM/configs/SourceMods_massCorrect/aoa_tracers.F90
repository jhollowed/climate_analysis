!===============================================================================
! Age of air test tracers
! provides dissipation rate and surface fluxes for diagnostic constituents
!===============================================================================

module aoa_tracers

  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use ppgrid,       only: pcols, pver
  use constituents, only: pcnst, cnst_add, cnst_name, cnst_longname
  use cam_logfile,  only: iulog
  use ref_pres,     only: pref_mid_norm

  implicit none
  private
  save

  ! Public interfaces
  public :: aoa_tracers_register         ! register constituents
  public :: aoa_tracers_implements_cnst  ! true if named constituent is implemented by this package
  public :: aoa_tracers_init_cnst        ! initialize constituent field
  public :: aoa_tracers_init             ! initialize history fields, datasets
  public :: aoa_tracers_timestep_init    ! place to perform per timestep initialization
  public :: aoa_tracers_timestep_tend    ! calculate tendencies
  public :: aoa_tracers_readnl           ! read namelist options

  ! Private module data

  integer, parameter :: ncnst=4  ! number of constituents implemented by this module

  ! constituent names
  character(len=8), parameter :: c_names(ncnst) = (/'SAI_SO2', 'SAI_ASH', 'SAI_PT', 'SAI_AOA'/)
  !character(len=8), parameter :: c_names(ncnst) = (/'AOA1', 'AOA2', 'HORZ', 'VERT'/)

  ! constituent source/sink names
  !character(len=8), parameter :: src_names(ncnst) = (/'AOA1SRC', 'AOA2SRC', 'HORZSRC', 'VERTSRC'/)

  integer :: ifirst ! global index of first constituent
  integer :: ixsai1 ! global index for AOA1 tracer
  integer :: ixsai2 ! global index for AOA2 tracer
  integer :: ixsai3   ! global index for HORZ tracer
  integer :: ixsai4   ! global index for VERT tracer

  ! Data from namelist variables
  logical :: aoa_tracers_flag  = .false.    ! true => turn on test tracer code, namelist variable
  logical :: aoa_read_from_ic_file = .true. ! true => tracers initialized from IC file

!===============================================================================
contains
!===============================================================================

!================================================================================
  subroutine aoa_tracers_readnl(nlfile)

    use namelist_utils,     only: find_group_name
    use units,              only: getunit, freeunit
    use mpishorthand
    use cam_abortutils,     only: endrun

    implicit none

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aoa_tracers_readnl'


    namelist /aoa_tracers_nl/ aoa_tracers_flag, aoa_read_from_ic_file

    !-----------------------------------------------------------------------------

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'aoa_tracers_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aoa_tracers_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    call mpibcast(aoa_tracers_flag, 1, mpilog,  0, mpicom)
    call mpibcast(aoa_read_from_ic_file, 1, mpilog,  0, mpicom)
#endif

  endsubroutine aoa_tracers_readnl

!================================================================================

  subroutine aoa_tracers_register
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents
    !
    !-----------------------------------------------------------------------
    use physconst,  only: cpair, mwdry
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    call cnst_add(c_names(1), mwdry, cpair, 0._r8, ixsai1, readiv=aoa_read_from_ic_file, &
                  longname='Stratospheric aerosol injection plume SO2')
    ifirst = ixsai1
    call cnst_add(c_names(2), mwdry, cpair, 0._r8, ixsai2, readiv=aoa_read_from_ic_file, &
                  longname='Stratospheric aerosol injection plume ash')
    call cnst_add(c_names(3), mwdry, cpair, 0._r8, ixsai3,   readiv=aoa_read_from_ic_file, &
                  longname='potential temperature at initial time of stratospheric aerosol injection')
    call cnst_add(c_names(4), mwdry, cpair, 0._r8, ixsai4,   readiv=aoa_read_from_ic_file, &
                  longname='stratospheric aeosol injection plume clock tracer')

  end subroutine aoa_tracers_register

!===============================================================================

  function aoa_tracers_implements_cnst(name)
    !-----------------------------------------------------------------------
    !
    ! Purpose: return true if specified constituent is implemented by this package
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: aoa_tracers_implements_cnst        ! return value

    !---------------------------Local workspace-----------------------------
    integer :: m
    !-----------------------------------------------------------------------

    aoa_tracers_implements_cnst = .false.

    if (.not. aoa_tracers_flag) return

    do m = 1, ncnst
       if (name == c_names(m)) then
          aoa_tracers_implements_cnst = .true.
          return
       end if
    end do

  end function aoa_tracers_implements_cnst

!===============================================================================

  subroutine aoa_tracers_init_cnst(name, latvals, lonvals, mask, q)

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize test tracers mixing ratio fields
    !  This subroutine is called at the beginning of an initial run ONLY
    !
    !-----------------------------------------------------------------------

    character(len=*), intent(in)  :: name
    real(r8),         intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8),         intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,          intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8),         intent(out) :: q(:,:)   ! kg tracer/kg dry air (gcol, plev)

    integer :: m
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    do m = 1, ncnst
       if (name ==  c_names(m))  then
          ! pass global constituent index
          call init_cnst_3d(ifirst+m-1, latvals, lonvals, mask, q)
       endif
    end do

  end subroutine aoa_tracers_init_cnst

!===============================================================================

  subroutine aoa_tracers_init

    !-----------------------------------------------------------------------
    !
    ! Purpose: initialize age of air constituents
    !          (declare history variables)
    !-----------------------------------------------------------------------

    use cam_history,    only: addfld, add_default

    integer :: m, mm, k
    !-----------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    ! Set names of tendencies and declare them as history variables
   
    call addfld('SAI_MASS',  (/ 'lev' /), 'A', 'kg', 'mass of grid box' )
    call add_default('SAI_MASS', 1, ' ')

    do m = 1, ncnst
       mm = ifirst+m-1
       call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', cnst_longname(mm))
       
       call add_default(cnst_name(mm), 1, ' ')
       !call add_default (src_names(m),  1, ' ')
    end do

  end subroutine aoa_tracers_init

!===============================================================================

  subroutine aoa_tracers_timestep_init( phys_state )
    !-----------------------------------------------------------------------
    ! Provides a place to reinitialize diagnostic constituents HORZ and VERT
    !-----------------------------------------------------------------------

    use time_manager,   only: get_curr_date
    use ppgrid,         only: begchunk, endchunk
    use physics_types,  only: physics_state

    type(physics_state), intent(inout), dimension(begchunk:endchunk), optional :: phys_state


    integer c, i, k, ncol
    integer yr, mon, day, tod
    !--------------------------------------------------------------------------

    if (.not. aoa_tracers_flag) return

    call get_curr_date (yr,mon,day,tod)

    if ( day == 1 .and. tod == 0) then
       if (masterproc) then
         write(iulog,*) 'SAI_CONSTITUENTS: RE-INITIALIZING CONSTITUENTS'
       endif

       do c = begchunk, endchunk
          ncol = phys_state(c)%ncol
          do k = 1, pver
             do i = 1, ncol
                !phys_state(c)%q(i,k,ixsai3) = 2._r8 + sin(phys_state(c)%lat(i))
                !phys_state(c)%q(i,k,ixsai4) = qrel_vert(k)
                continue
             end do
          end do
       end do

    end if

  end subroutine aoa_tracers_timestep_init

!===============================================================================

  subroutine aoa_tracers_timestep_tend(state, ptend, cflx, landfrac, dt, ncol)

    use physics_types, only: physics_state, physics_ptend, physics_ptend_init
    use cam_history,   only: outfld
    use time_manager,  only: get_nstep
  
    !--JH--
    use ref_pres,     only: pref_mid_norm
    use time_manager, only: get_curr_time
    use physconst,    only: pi, rearth, cpair, rair 
    use phys_grid,    only : get_area_all_p

    ! Arguments
    !type(physics_state), intent(in)   :: state              ! state variables
    type(physics_state), intent(inout) :: state              ! --JH--
    type(physics_ptend), intent(out)   :: ptend              ! package tendencies
    real(r8),            intent(inout) :: cflx(pcols,pcnst)  ! Surface constituent flux (kg/m^2/s)
    real(r8),            intent(in)    :: landfrac(pcols)    ! Land fraction
    real(r8),            intent(in)    :: dt                 ! timestep
    integer,             intent(in)    :: ncol               ! no. of column in chunk

    !----------------- Local workspace-------------------------------

    integer  :: i, k
    integer  :: lchnk             ! chunk identifier
    integer  :: nstep             ! current timestep number
    logical  :: lq(pcnst)

    !--JH--
    real(r8), parameter :: deg2rad = pi/180._r8
    real(r8), parameter :: rad2deg = 180._r8/pi
    integer  :: day,sec
    real(r8) :: t                 
    real(r8) :: lat                 
    real(r8) :: lon                 
    real(r8) :: zz                 
    real(r8) :: rhoatm
    real(r8) :: rr             
    real(r8) :: lat0                 
    real(r8) :: lon0                 
    real(r8) :: z0                 
    real(r8) :: zs                 
    real(r8) :: rs                
    real(r8) :: dz                 
    real(r8) :: dr                 
    real(r8) :: tau                 
    real(r8) :: tf 
    real(r8) :: A_so2                 
    real(r8) :: A_ash                 
    real(r8) :: M_so2                
    real(r8) :: M_ash                 
    real(r8) :: k_so2                 
    real(r8) :: k_ash                 
    real(r8) :: P0             
    real(r8) :: alpha 
    real(r8) :: tmin 
    real(r8) :: tmmin 
    real(r8) :: Hint 
    real(r8) :: Vint 
    ! for mass correction
    real(r8) :: rhoi 
    real(r8) :: ksol_so2 
    real(r8) :: ksol_ash  
    ! for mass estiamte
    real(r8) :: area(ncol), mass(ncol,pver)
    

    !------------------------------------------------------------------

    if (.not. aoa_tracers_flag) then
       call physics_ptend_init(ptend,state%psetcols,'none') !Initialize an empty ptend for use with physics_update
       return
    end if

    lq(:)      = .FALSE.
    lq(ixsai1) = .TRUE.
    lq(ixsai2) = .TRUE.
    lq(ixsai3)   = .TRUE.
    lq(ixsai4)   = .TRUE.
    call physics_ptend_init(ptend,state%psetcols, 'aoa_tracers', lq=lq)
    
    nstep = get_nstep()
    lchnk = state%lchnk

    !--JH--
    ! get current time in seconds
    call get_curr_time(day,sec)
    t      = (day*24.0*60.0*60.0) + sec
    lat0   = 15.15 * deg2rad
    lon0   = 120.35 * deg2rad
    z0     = 25000.0
    dz     = 7500.0                
    dr     = 100000.0                
    rs     = 3.0*dr
    zs     = 3.0*dz                
    tf     = 172800.0        
    tau    = -LOG(0.05)/tf                
    M_so2  = 2.0e10
    M_ash  = 2.0e10
    k_so2  = 1.0/2592000.0
    k_ash  = 1.0/86400.0
    P0     = 100000.0
    ! select exponential time decay
    !alpha  = tf  ! constant T(t)
    alpha  = (1/tau) * (1 - EXP(-tau*tf))  ! exponential T(t)
    ! modified value of k to use for integral solution below
    ! ksol_so2 = k_so2  ! constant T(t)
    ! ksol_ash = k_so2  ! constant T(t)
    ksol_so2 = k_so2-tau  ! exponential T(t)
    ksol_ash = k_so2-tau  ! exponential T(t)

    Hint = 1.0 - EXP(-rs**2.0/(2.0*dr**2.0))
    Vint = ERF(z0/(SQRT(2.0)*dz)) - ERF((z0-zs)/(SQRT(2.0)*dz))
    
    A_so2 = M_so2 / (alpha * SQRT(2.0*pi**3.0) * dr**2.0 * dz) 
    A_so2 = A_so2 * 1.0/(Hint * Vint)     
    A_ash = M_ash / (alpha * SQRT(2.0*pi**3.0) * dr**2.0 * dz) 
    A_ash = A_ash * 1.0/(Hint * Vint)

    ! get area of column (assume height ~ a)
    call get_area_all_p(lchnk, ncol, area)
    area = area * rearth**2

    do k = 1, pver
       do i = 1, ncol
          
          lat = state%lat(i)
          lon = state%lon(i)
          zz  = state%zm(i, k)
          rr = (rearth+z0) * ACOS( SIN(lat)*SIN(lat0) + COS(lat)*COS(lat0) * COS(ABS(lon-lon0)))
          rhoatm = state%pmid(i, k) / (rair * state%t(i, k))
          
          if(t <= tf) then
              tmin = t
          else
              tmin = tf
          end if
          if((t-dt) <= tf) then
              tmmin = t-dt
          else
              tmmin = tf
          end if
          rhoi = state%q(i, k, ixsai1) * rhoatm ! current density (before injection)

          ! ---------- SAI_SO2 ----------
          ptend%q(i,k,ixsai1) = EXP(-k_so2*t) / ksol_so2 * &             ! rho(t_i)
                                (rhoi*ksol_so2 + &
                                 A_so2 * EXP(-(1.0/2.0) * (rr/dr)**2.0) * &
                                         EXP(-(1.0/2.0) * ((zz-z0)/dz)**2.0) * &
                                         (EXP(ksol_so2 * tmin) - 1))
          ptend%q(i,k,ixsai1) = ptend%q(i,k,ixsai1) - &           ! rho(t_i-1)
                                EXP(-k_so2*(t-dt)) / ksol_so2 * &
                                (rhoi*ksol_so2 + &
                                 A_so2 * EXP(-(1.0/2.0) * (rr/dr)**2.0) * &
                                         EXP(-(1.0/2.0) * ((zz-z0)/dz)**2.0) * &
                                         (EXP(ksol_so2 * tmmin) - 1))
          ptend%q(i,k,ixsai1) = ptend%q(i,k,ixsai1) / dt           ! divide Delta rho / Delta t
          
          ! scale density to concentration
          ptend%q(i,k,ixsai1) = ptend%q(i,k,ixsai1) / rhoatm


          ! ---------- SAI_ASH ----------
          ptend%q(i,k,ixsai2) = EXP(-k_ash*t) / ksol_ash * &             ! rho(t_i)
                                (rhoi*ksol_ash + &
                                 A_ash * EXP(-(1.0/2.0) * (rr/dr)**2.0) * &
                                         EXP(-(1.0/2.0) * ((zz-z0)/dz)**2.0) * &
                                         (EXP(ksol_ash * tmin) - 1))
          ptend%q(i,k,ixsai2) = ptend%q(i,k,ixsai2) - &           ! rho(t_i-1)
                                EXP(-k_ash*(t-dt)) / ksol_ash * &
                                (rhoi*ksol_ash + &
                                 A_ash * EXP(-(1.0/2.0) * (rr/dr)**2.0) * &
                                         EXP(-(1.0/2.0) * ((zz-z0)/dz)**2.0) * &
                                         (EXP(ksol_ash * tmmin) - 1))
          ptend%q(i,k,ixsai2) = ptend%q(i,k,ixsai2) / dt           ! divide Delta rho / Delta t
          
          ! scale density to concentration
          ptend%q(i,k,ixsai2) = ptend%q(i,k,ixsai2) / rhoatm


          ! ---------- SAI_PT ----------
          ! --JH--: Potential Temperature
          ! initialize within the first minute of the injection
          if (t < 60.0_r8) then
              state%q(i,k,ixsai3) = state%t(i,k) * (P0 / state%pmid(i, k))**(rair/cpair)
          end if
          ptend%q(i,k,ixsai3) = 0.0_r8


          ! ---------- SAI_AOA ----------
          ptend%q(i,k,ixsai4) = 0.0_r8
          
          ! ---------- MASS ----------
          ! mass(i,k) = state%pdel(i,k) * area(i) * 1/(rair * state%t(i,k))
          mass(i,k) = state%pdel(i,k)


       end do
    end do


    ! record tracer mass on history files
    call outfld('SAI_MASS', mass(:,:), ncol, lchnk)

    ! Set tracer fluxes
    do i = 1, ncol
       cflx(i,ixsai1) = 0._r8
       cflx(i,ixsai2) = 0._r8
       cflx(i,ixsai3) = 0._r8
       cflx(i,ixsai4) = 0._r8
    end do

  end subroutine aoa_tracers_timestep_tend

!===========================================================================

  subroutine init_cnst_3d(m, latvals, lonvals, mask, q)

    integer,  intent(in)  :: m          ! global constituent index
    real(r8), intent(in)  :: latvals(:) ! lat in degrees (ncol)
    real(r8), intent(in)  :: lonvals(:) ! lon in degrees (ncol)
    logical,  intent(in)  :: mask(:)    ! Only initialize where .true.
    real(r8), intent(out) :: q(:,:)     ! kg tracer/kg dry air (gcol,plev)

    integer :: j, k, gsize
    !-----------------------------------------------------------------------

    if (masterproc) then
      write(iulog,*) 'AGE-OF-AIR CONSTITUENTS: INITIALIZING ',cnst_name(m),m
    end if

    !--JH--: initialize everything to zero 
    if (m == ixsai1) then

       q(:,:) = 0.0_r8

    else if (m == ixsai2) then

       q(:,:) = 0.0_r8

    else if (m == ixsai3) then

       q(:,:) = 0.0_r8

    else if (m == ixsai4) then

       q(:,:) = 0.0_r8

    end if

  end subroutine init_cnst_3d

!=====================================================================


end module aoa_tracers
