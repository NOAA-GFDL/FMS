!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!> \brief ensemble_manager_mod
!!
module ensemble_manager_mod


  use fms_mod, only : open_namelist_file,close_file,check_nml_error
  use mpp_mod, only : mpp_npes, stdout, stdlog, mpp_error, FATAL
  use mpp_mod, only : mpp_pe, mpp_declare_pelist
  use mpp_mod, only : input_nml_file
  use fms_io_mod, only       : set_filename_appendix

  IMPLICIT NONE

  private

  integer, parameter :: MAX_ENSEMBLE_SIZE = 100


  integer, allocatable, dimension(:,:) :: ensemble_pelist
  integer, allocatable, dimension(:,:) :: ensemble_pelist_ocean
  integer, allocatable, dimension(:,:) :: ensemble_pelist_atmos
  integer, allocatable, dimension(:,:) :: ensemble_pelist_land
  integer, allocatable, dimension(:,:) :: ensemble_pelist_ice
  integer, allocatable, dimension(:)   :: ensemble_pelist_ocean_filter
  integer, allocatable, dimension(:)   :: ensemble_pelist_atmos_filter
  integer, allocatable, dimension(:)   :: ensemble_pelist_land_filter
  integer, allocatable, dimension(:)   :: ensemble_pelist_ice_filter

  integer :: ensemble_size = 1
  integer :: ensemble_id = 1
  integer :: pe, total_npes_pm=0,ocean_npes_pm=0,atmos_npes_pm=0
  integer :: land_npes_pm=0,ice_npes_pm=0

  public :: ensemble_manager_init, get_ensemble_id, get_ensemble_size, get_ensemble_pelist
  public :: ensemble_pelist_setup
  public :: get_ensemble_filter_pelist
contains

!> \brief ensemble_manager_init
!!
!! \throw FATAL, "ensemble_manager_mod: ensemble_nml variable ensemble_size must be a positive integer"
!! \throw FATAL, "ensemble_manager_mod: ensemble_nml variable ensemble_size should be no larger than MAX_ENSEMBLE_SIZE, change ensemble_size or increase MAX_ENSEMBLE_SIZE"
!! \throw FATAL, "ensemble_size must be divis by npes"
!! \throw FATAL, "get_ensemble_pelist: size of pelist 1st index < ensemble_size"
!! \throw FATAL, "get_ensemble_pelist: size of pelist 2nd index < ocean_npes_pm"
!! \throw FATAL, "get_ensemble_pelist: size of pelist 2nd index < atmos_npes_pm"
!! \throw FATAL, "get_ensemble_pelist: size of pelist 2nd index < land_npes_pm"
!! \throw FATAL, "get_ensemble_pelist: size of pelist 2nd index < ice_npes_pm"
!! \throw FATAL, "get_ensemble_pelist: unknown argument name=[name]"
!! \throw FATAL, "get_ensemble_pelist: size of pelist 2nd index < total_npes_pm"
  subroutine ensemble_manager_init()


    integer :: i, io_status, ioun, npes, ierr

    namelist /ensemble_nml/ ensemble_size

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, ensemble_nml, iostat=io_status)
#else
    ioun = open_namelist_file()
    read(ioun,nml=ensemble_nml,iostat = io_status)
    call close_file(ioun)
#endif
    ierr = check_nml_error(io_status, 'ensemble_nml')

    if(ensemble_size < 1) call mpp_error(FATAL, &
       'ensemble_manager_mod: ensemble_nml variable ensemble_size must be a positive integer')
    if(ensemble_size > max_ensemble_size)  call mpp_error(FATAL, &
       'ensemble_manager_mod: ensemble_nml variable ensemble_size should be no larger than MAX_ENSEMBLE_SIZE, '// &
       'change ensemble_size or increase MAX_ENSEMBLE_SIZE')

    pe = mpp_pe()
    npes = mpp_npes()
    if (npes < ensemble_size) then
      call mpp_error(FATAL,'npes must be >= ensemble_size')
    endif
    total_npes_pm = npes/ensemble_size
    if (mod(npes, total_npes_pm) /= 0) call mpp_error(FATAL,'ensemble_size must be divis by npes')

    call mpp_declare_pelist((/(i,i=0,npes-1)/),'_ens0') ! for ensemble driver

  end subroutine ensemble_manager_init

  function get_ensemble_id()
    integer :: get_ensemble_id
    get_ensemble_id = ensemble_id
  end function get_ensemble_id

  function get_ensemble_size()

    integer, dimension(6) :: get_ensemble_size

    get_ensemble_size(1) = ensemble_size
    get_ensemble_size(2) = total_npes_pm
    get_ensemble_size(3) = ocean_npes_pm
    get_ensemble_size(4) = atmos_npes_pm
    get_ensemble_size(5) = land_npes_pm
    get_ensemble_size(6) = ice_npes_pm

  end function get_ensemble_size


  subroutine get_ensemble_pelist(pelist, name)

    integer, intent(inout) :: pelist(:,:)
    character(len=*), intent(in), optional  :: name

    if (size(pelist,1) < ensemble_size) &
         call mpp_error(FATAL,'get_ensemble_pelist: size of pelist 1st index < ensemble_size')

    if(present(name)) then
       select case(name)
       case('ocean')
          if (size(pelist,2) < ocean_npes_pm)&
               call mpp_error(FATAL,'get_ensemble_pelist: size of pelist 2nd index < ocean_npes_pm')
          pelist = 0
          pelist(1:ensemble_size,1:ocean_npes_pm) = &
               ensemble_pelist_ocean(1:ensemble_size,1:ocean_npes_pm)

       case('atmos')
          if (size(pelist,2) < atmos_npes_pm)&
               call mpp_error(FATAL,'get_ensemble_pelist: size of pelist 2nd index < atmos_npes_pm')
          pelist = 0
          pelist(1:ensemble_size,1:atmos_npes_pm) = &
               ensemble_pelist_atmos(1:ensemble_size,1:atmos_npes_pm)

       case('land')
          if (size(pelist,2) < land_npes_pm)&
               call mpp_error(FATAL,'get_ensemble_pelist: size of pelist 2nd index < land_npes_pm')
          pelist = 0
          pelist(1:ensemble_size,1:land_npes_pm) = &
               ensemble_pelist_land(1:ensemble_size,1:land_npes_pm)

       case('ice')
          if (size(pelist,2) < ice_npes_pm)&
               call mpp_error(FATAL,'get_ensemble_pelist: size of pelist 2nd index < ice_npes_pm')
          pelist = 0
          pelist(1:ensemble_size,1:ice_npes_pm) = &
               ensemble_pelist_ice(1:ensemble_size,1:ice_npes_pm)

       case default
          call mpp_error(FATAL,'get_ensemble_pelist: unknown argument name='//name)
       end select
    else
       if (size(pelist,2) < total_npes_pm)&
            call mpp_error(FATAL,'get_ensemble_pelist: size of pelist 2nd index < total_npes_pm')
       pelist = 0
       pelist(1:ensemble_size,1:total_npes_pm) = &
            ensemble_pelist(1:ensemble_size,1:total_npes_pm)
    endif

    return
  end subroutine get_ensemble_pelist

!> \brief get_ensemble_filter_pelist
!!
!! \throw FATAL, "get_ensemble_filter_pelist: size of pelist argument < ensemble_size * ocean_npes_pm"
!! \throw FATAL, "get_ensemble_filter_pelist: size of pelist argument < ensemble_size * atmos_npes_pm"
!! \throw FATAL, "get_ensemble_filter_pelist: size of pelist argument < ensemble_size * land_npes_pm"
!! \throw FATAL, "get_ensemble_filter_pelist: size of pelist argument < ensemble_size * ice_npes_pm"
!! \throw FATAL, "get_ensemble_filter_pelist: unknown argument name=[name]"
  subroutine get_ensemble_filter_pelist(pelist, name)

    integer, intent(inout) :: pelist(:)
    character(len=*), intent(in)  :: name

    select case(name)
    case('ocean')
       if (size(pelist) < ensemble_size * ocean_npes_pm)&
            call mpp_error(FATAL,'get_ensemble_filter_pelist: size of pelist argument < ensemble_size * ocean_npes_pm')
       pelist = 0
       pelist(1:ensemble_size*ocean_npes_pm) = &
            ensemble_pelist_ocean_filter(1:ensemble_size*ocean_npes_pm)

    case('atmos')
       if (size(pelist) < ensemble_size * atmos_npes_pm)&
            call mpp_error(FATAL,'get_ensemble_filter_pelist: size of pelist argument < ensemble_size * atmos_npes_pm')
       pelist = 0
       pelist(1:ensemble_size*atmos_npes_pm) = &
            ensemble_pelist_atmos_filter(1:ensemble_size*atmos_npes_pm)

    case('land')
       if (size(pelist) < ensemble_size * land_npes_pm)&
            call mpp_error(FATAL,'get_ensemble_filter_pelist: size of pelist argument < ensemble_size * land_npes_pm')
       pelist = 0
       pelist(1:ensemble_size*land_npes_pm) = &
            ensemble_pelist_land_filter(1:ensemble_size*land_npes_pm)

    case('ice')
       if (size(pelist) < ensemble_size * ice_npes_pm)&
            call mpp_error(FATAL,'get_ensemble_filter_pelist: size of pelist argument < ensemble_size * ice_npes_pm')
       pelist = 0
       pelist(1:ensemble_size*ice_npes_pm) = &
            ensemble_pelist_ice_filter(1:ensemble_size*ice_npes_pm)

    case default
       call mpp_error(FATAL,'get_ensemble_filter_pelist: unknown argument name='//name)
    end select


    return
  end subroutine get_ensemble_filter_pelist

!nnz: I think the following block of code should be contained in a subroutine
!     to consolidate and ensure the consistency of declaring the various pelists.
!>\brief ensemble_pelist_setup
!!
!! \throw FATAL, "ensemble_manager_mod: land_npes > atmos_npes"
!! \throw FATAL, "ensemble_manager_mod: ice_npes > atmos_npes"
  subroutine ensemble_pelist_setup(concurrent, atmos_npes, ocean_npes, land_npes, ice_npes, &
                                   Atm_pelist, Ocean_pelist, Land_pelist, Ice_pelist)
    logical, intent(in)                  :: concurrent
    integer, intent(in)                  :: atmos_npes, ocean_npes
    integer, intent(in)                  :: land_npes, ice_npes
    integer, dimension(:), intent(inout) :: Atm_pelist, Ocean_pelist
    integer, dimension(:), intent(inout) :: Land_pelist, Ice_pelist
    integer           :: atmos_pe_start, atmos_pe_end, ocean_pe_start, ocean_pe_end
    integer           :: land_pe_start, land_pe_end, ice_pe_start, ice_pe_end
    character(len=10) :: pelist_name, text
    integer           :: npes, n, m ,i

    npes = total_npes_pm

    ! make sure land_npes and ice_npes are not greater than atmos_npes
    if(land_npes > atmos_npes) call mpp_error(FATAL, 'ensemble_manager_mod: land_npes > atmos_npes')
    if(ice_npes  > atmos_npes) call mpp_error(FATAL, 'ensemble_manager_mod: ice_npes > atmos_npes')

    allocate(ensemble_pelist(ensemble_size,total_npes_pm))
    allocate(ensemble_pelist_ocean(1:ensemble_size, 1:ocean_npes) )
    allocate(ensemble_pelist_atmos(1:ensemble_size, 1:atmos_npes) )
    allocate(ensemble_pelist_land (1:ensemble_size, 1:land_npes) )
    allocate(ensemble_pelist_ice  (1:ensemble_size, 1:ice_npes) )

    atmos_pe_start = 0
    ocean_pe_start = 0
    land_pe_start = 0
    ice_pe_start = 0
    if( concurrent .OR. atmos_npes+ocean_npes == npes )then
       ocean_pe_start = ensemble_size*atmos_npes
    endif
    do n=1,ensemble_size
       atmos_pe_end = atmos_pe_start + atmos_npes - 1
       ocean_pe_end = ocean_pe_start + ocean_npes - 1
       land_pe_end  = land_pe_start  + land_npes  - 1
       ice_pe_end   = ice_pe_start   + ice_npes   - 1
       ensemble_pelist_atmos(n, 1:atmos_npes) = (/(i,i=atmos_pe_start,atmos_pe_end)/)
       ensemble_pelist_ocean(n, 1:ocean_npes) = (/(i,i=ocean_pe_start,ocean_pe_end)/)
       ensemble_pelist_land (n, 1:land_npes) = (/(i,i=land_pe_start, land_pe_end)/)
       ensemble_pelist_ice  (n, 1:ice_npes) = (/(i,i=ice_pe_start,  ice_pe_end)/)
       ensemble_pelist(n, 1:atmos_npes)       = ensemble_pelist_atmos(n, 1:atmos_npes)
       if( concurrent .OR. atmos_npes+ocean_npes == npes ) &
            ensemble_pelist(n, atmos_npes+1:npes)  = ensemble_pelist_ocean(n, 1:ocean_npes)
       if(ANY(ensemble_pelist(n,:) == pe)) ensemble_id = n
       write(pelist_name,'(a,i2.2)')  '_ens',n
       call mpp_declare_pelist(ensemble_pelist(n,:), trim(pelist_name))
       atmos_pe_start = atmos_pe_end + 1
       ocean_pe_start = ocean_pe_end + 1
       land_pe_start  = atmos_pe_start
       ice_pe_start   = atmos_pe_start
    enddo

    Atm_pelist(:)   = ensemble_pelist_atmos(ensemble_id,:)
    Ocean_pelist(:) = ensemble_pelist_ocean(ensemble_id,:)
    Land_pelist(:)  = ensemble_pelist_land (ensemble_id,:)
    Ice_pelist(:)   = ensemble_pelist_ice  (ensemble_id,:)

    !    write(pelist_name,'(a,i2.2)')  '_ocn_ens',ensemble_id
    !    call mpp_declare_pelist(Ocean%pelist , trim(pelist_name) )

    !    write(pelist_name,'(a,i2.2)')  '_atm_ens',ensemble_id
    !    call mpp_declare_pelist(Atm%pelist , trim(pelist_name) )
    !
    !nnz: The above is sufficient for non-concurrent mode.
    !     BUT
    !     For atmosphere_init to work in ensemble, concurrent mode
    !     the following Atm_pelist should be declared (per ensemble member)
    !     instead of the above Atm%pelist!
    !
    !   allocate( Atm_pelist(1:ensemble_size, 1:atmos_npes) )
    !   do n=1,ensemble_size
    !       do i=1, atmos_npes
    !          Atm_pelist(n, i) = ensemble_pelist(n, i)
    !       enddo
    !       write(pelist_name,'(a,i2.2)')  '_atm_ens',n
    !       call mpp_declare_pelist(Atm_pelist(n,:) , trim(pelist_name) )
    !    enddo
    !
    !     The way I understand this with the help of Totalview is:
    !     With mpp_declare_pelist(Atm%pelist)
    !     When we are in fv_arrays_init when mp_init(comID) is called
    !     comID is the same for the atmospheric PE's for both ensemble members
    !     since peset(5)%id is the same (7) for those PE's, so the PE count is double what it should be inside
    !     mp_init().
    !     It is also true that for Ocean PE's, peset(4)%id is the same (6) for Ocean PE's in both ensemble members
    !     but for Ocean it is not a problem because Ocean is not trying to create new communicators
    !     from this peset whereas ATM does (vis mp_init).
    !
    !     Who sets peset(i)%id ? Can it be modified to assign different %id for the two subsets.
    !     peset(i)%id = 0 for Ocean PE's on ATM pesets and for ATM PE's on Ocean pesets.
    !
    !     With mpp_declare_pelist(Atm_pelist(n,:)) n=1,...,ensemble_size
    !     we get separate pesets for each ATM ensemble member and each with a different %id and mp_init is cured.
    !
    !     There is also a matter of precedence. If we have both calls
    !     call mpp_declare_pelist(Atm%pelist , trim(pelist_name) )
    !     and
    !     call mpp_declare_pelist(Atm_pelist(n,:) , trim(pelist_name) )
    !     then concurrent run fails because with call mpp_set_current_pelist( Atm%pelist   )
    !     peset(i) is searched for i=1,2,... and the first pelist that matches argument, its peset is set as current.
    !
    !     To be consistent with ATM and OCEAN we can do the following
    !     (eventhough mpp_declare_pelist(Ocean%pelist) is adequate right now.)

    if( concurrent )then
       do n=1,ensemble_size
          write(pelist_name,'(a,i2.2)')  'atm_ens',n
          call mpp_declare_pelist(ensemble_pelist_atmos(n,:) , trim(pelist_name) )
          write(pelist_name,'(a,i2.2)')  'ocn_ens',n
          call mpp_declare_pelist(ensemble_pelist_ocean(n,:) , trim(pelist_name) )
          write(pelist_name,'(a,i2.2)')  'lnd_ens',n
          call mpp_declare_pelist(ensemble_pelist_land(n,:) , trim(pelist_name) )
          write(pelist_name,'(a,i2.2)')  'ice_ens',n
          call mpp_declare_pelist(ensemble_pelist_ice(n,:) , trim(pelist_name) )
       enddo
    else
       write(pelist_name,'(a,i2.2)')  'atm_ens',ensemble_id
       call mpp_declare_pelist(Atm_pelist , trim(pelist_name) )
       write(pelist_name,'(a,i2.2)')  'ocn_ens',ensemble_id
       call mpp_declare_pelist(Ocean_pelist , trim(pelist_name) )
       write(pelist_name,'(a,i2.2)')  'lnd_ens',ensemble_id
       call mpp_declare_pelist(Land_pelist , trim(pelist_name) )
       write(pelist_name,'(a,i2.2)')  'ice_ens',ensemble_id
       call mpp_declare_pelist(Ice_pelist , trim(pelist_name) )
    endif

    ocean_npes_pm = ocean_npes
    atmos_npes_pm = atmos_npes
    land_npes_pm  = land_npes
    ice_npes_pm   = ice_npes

    !Declare pelist of all Ocean, Atmos, Land and Ice pes across all ensembles ( filters )
    allocate(ensemble_pelist_ocean_filter(ensemble_size * ocean_npes_pm))
    allocate(ensemble_pelist_atmos_filter(ensemble_size * atmos_npes_pm))
    allocate(ensemble_pelist_land_filter (ensemble_size * land_npes_pm))
    allocate(ensemble_pelist_ice_filter  (ensemble_size * ice_npes_pm))
    do n=1,ensemble_size
       do m=1,ocean_npes_pm
          i=(n-1)*ocean_npes_pm + m
          ensemble_pelist_ocean_filter(i) = ensemble_pelist_ocean(n,m)
       enddo
       do m=1,atmos_npes_pm
          i=(n-1)*atmos_npes_pm + m
          ensemble_pelist_atmos_filter(i) = ensemble_pelist_atmos(n,m)
       enddo
       do m=1,land_npes_pm
          i=(n-1)*land_npes_pm + m
          ensemble_pelist_land_filter(i)  = ensemble_pelist_land(n,m)
       enddo
       do m=1,ice_npes_pm
          i=(n-1)*ice_npes_pm + m
          ensemble_pelist_ice_filter(i)   = ensemble_pelist_ice(n,m)
       enddo
    enddo

    write(pelist_name,'(a)')  'ocn_filter'
    call mpp_declare_pelist(ensemble_pelist_ocean_filter, trim(pelist_name) )

    write(pelist_name,'(a)')  'atm_filter'
    call mpp_declare_pelist(ensemble_pelist_atmos_filter, trim(pelist_name) )

    write(pelist_name,'(a)')  'lnd_filter'
    call mpp_declare_pelist(ensemble_pelist_land_filter,  trim(pelist_name) )

    write(pelist_name,'(a)')  'ice_filter'
    call mpp_declare_pelist(ensemble_pelist_ice_filter,   trim(pelist_name) )

    !
    !Rename output files to identify the ensemble
    !If ensemble_size=1 do not rename files so that the same coupler
    !can be used for non-ensemble experiments
    !
    if (ensemble_size > 1) then
       write( text,'(a,i2.2)' ) 'ens_', ensemble_id
       !Append ensemble_id to the restart filenames
       call set_filename_appendix(trim(text))
    endif

  end subroutine ensemble_pelist_setup


end module ensemble_manager_mod
