! -*- f90 -*-
module oda_core_ecda_mod
  ! FMS Shared modules
  use fms_mod, only : file_exist, read_data
  use fms_mod, only : open_namelist_file, check_nml_error, close_file
  use fms_mod, only : error_mesg, FATAL, NOTE
#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#endif
  use mpp_mod, only : mpp_sum, stdout, stdlog, mpp_sync_self
  use mpp_mod, only : mpp_pe, mpp_root_pe
  use mpp_io_mod, only : mpp_open, mpp_close, MPP_ASCII, MPP_RDONLY, MPP_MULTI, MPP_SINGLE, MPP_NETCDF
  use mpp_io_mod, only : mpp_get_atts, mpp_get_info, mpp_get_fields, mpp_read, axistype, fieldtype, mpp_get_axes
  use mpp_io_mod, only : mpp_get_axis_data, mpp_get_field_name
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : domain2d, mpp_get_global_domain, mpp_update_domains
  use mpp_domains_mod, only : mpp_global_field
  use mpp_memutils_mod, only : mpp_print_memuse_stats
  use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time
  use time_manager_mod, only : operator( <= ), operator( - ), operator( > ), operator ( < )
  use get_cal_time_mod, only : get_cal_time
  use axis_utils_mod, only : frac_index  
  use horiz_interp_type_mod, only: horiz_interp_type
  use horiz_interp_bilinear_mod, only : horiz_interp_bilinear_new
  use constants_mod, only : DEG_TO_RAD

  ! ODA_tools modules
  use oda_types_mod, only : ocean_profile_type, ocn_obs_flag_type, grid_type, obs_clim_type
  use oda_types_mod, only : DROP_PROFILER, MOORING, SATELLITE, DRIFTER, SHIP, TEMP_ID, SALT_ID, MISSING_VALUE
  use oda_types_mod, only : UNKNOWN, TAO
  use xbt_adjust, only : xbt_drop_rate_adjust
  
  implicit none

  private

  public :: copy_obs, oda_core_init, open_profile_dataset, &
       get_obs, get_obs_woa05t, get_obs_woa05s, get_obs_sst, get_obs_suv, &
       get_obs_eta, open_profile_dataset_sst, ocn_obs, ssh_td, max_profiles

  ! Parameters
  integer, parameter :: PROFILE_FILE = 1
  integer, parameter :: SFC_FILE = 2

  ! oda_core_nml variables
  real :: max_misfit = 5.0 !< used to inflate observation errors where the difference from the first guess is large
  real :: ass_start_lat = -87.0 !< set obs domain
  real :: ass_end_lat = 87.0 !< set obs domain
  integer :: max_profiles = 50000
  namelist /oda_core_nml/ max_misfit, ass_start_lat, ass_end_lat, max_profiles

  ! Shared ocean_obs_nml namelist variables
  real :: eta_obs_start_lat = -80.0 !< set obs domain
  real :: eta_obs_end_lat = 85.0 !< set obs domain
  real :: sst_obs_start_lat = -82.0 !< set obs domain
  real :: sst_obs_end_lat = 89.0 !< set obs domain

  integer :: max_prflvs = 200 ! for vd test

  type(ocean_profile_type), target, dimension(:), allocatable :: profiles  

  integer :: num_profiles, no_sst, no_prf, no_temp, no_salt, no_suv, no_eta ! total number of observations 
  integer :: no_woa05

  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed  ! indices for local domain on model grid
  integer :: isg, ieg, jsg, jeg
  integer :: isd_filt, ied_filt, jsd_filt, jed_filt
  integer :: isd_flt0, ied_flt0, jsd_flt0, jed_flt0
  integer :: nk

  real, dimension(:,:), allocatable, save :: mask_tao

  ! sst obs grid information
  real, allocatable :: woa05_lon(:), woa05_lat(:), woa05_z(:)
  real, allocatable :: sst_lon(:), sst_lat(:), obs_sst(:,:)
  real, allocatable, save :: obs_woa05t(:,:,:), obs_woa05s(:,:,:)
  integer ::  nlon, nlat, nlev
  integer ::  nlon_woa, nlat_woa, nlev_woa
  
  ! time window for DROP, MOORING and SATELLITE data respectively

  type(time_type) , dimension(0:100), public :: time_window

  type(grid_type), pointer :: Grd

  type(horiz_interp_type), save :: Interp

  real, allocatable, dimension(:, :) :: x_grid, y_grid, x_grid_uv, y_grid_uv
  real :: lon_out(1, 1), lat_out(1, 1)

  type(ocn_obs_flag_type) :: ocn_obs

  integer :: ssh_td

    type obs_entry_type
       character(len=128) :: filename
       character(len=16)  :: file_type
    end type obs_entry_type

 
contains

  subroutine init_observations(time_s, time_e, filt_domain, localize)  
    type(time_type), intent(in) :: time_s, time_e
    type(domain2d), intent(in) :: filt_domain
    logical, intent(in), optional :: localize

    integer, parameter :: SUV_ID = 4, ETA_ID = 5, WOAT_ID = 11, WOAS_ID = 12

    ! ocean_obs_nml variables
    integer :: mooring_window = 5
    integer :: satellite_window = 10 
    integer :: drop_window = 30
    integer :: drifter_window = 30
    integer :: ship_window = 30
    integer :: unknown_window = 30

    logical :: prfs_obs, salt_obs, sst_obs, eta_obs, suv_obs
    logical :: temp_obs_argo, salt_obs_argo, temp_obs_gtspp
    logical :: temp_obs_woa05, salt_obs_woa05
    integer :: eta_obs_td = 10
    integer :: max_files = 30 
    integer :: max_files_argo = 10 
    integer :: max_files_gtspp = 10 
    namelist /ocean_obs_nml/ mooring_window, satellite_window, drop_window,&
         & drifter_window, ship_window, unknown_window,&
         & prfs_obs, salt_obs, sst_obs, eta_obs, suv_obs,&
         & temp_obs_argo, salt_obs_argo, temp_obs_gtspp,&
         & temp_obs_woa05, salt_obs_woa05, eta_obs_td,&
         & sst_obs_start_lat, sst_obs_end_lat, eta_obs_start_lat, eta_obs_end_lat,&
         & max_files, max_files_argo, max_files_gtspp

    integer :: i, j, n, obs_variable
    integer :: ioun, io_status, ierr
    integer :: stdout_unit, stdlog_unit
    integer :: nfiles, nrecs, unit
    integer :: nfiles_argo, nrecs_argo, unit_argo
    integer :: nfiles_gtspp, nrecs_gtspp, unit_gtspp
    integer, dimension(:), allocatable :: filetype
    integer, dimension(:), allocatable :: filetype_argo
    integer, dimension(:), allocatable :: filetype_gtspp

    character(len=128) :: input_files_woa05t, input_files_woa05s
    character(len=256) :: record
    character(len=128), dimension(:), allocatable :: input_files
    character(len=128), dimension(:), allocatable :: input_files_argo
    character(len=128), dimension(:), allocatable :: input_files_gtspp

    type(obs_entry_type) :: tbl_entry

    stdout_unit = stdout()
    stdlog_unit = stdlog()

#ifdef INTERNAL_FILE_NML
    read (input_nml_file, ocean_obs_nml, iostat=io_status)
#else
    ioun = open_namelist_file()
    read(UNIT=ioun, NML=ocean_obs_nml, IOSTAT=io_status)
    ierr = check_nml_error(io_status,'ocean_obs_nml')
    call close_file(ioun)    
#endif
    write (UNIT=stdlog_unit, NML=ocean_obs_nml)

    ! Allocate filetype* and input_files* variables
    allocate(filetype(max_files), input_files(max_files))
    allocate(filetype_argo(max_files_argo), input_files_argo(max_files_argo))
    allocate(filetype_gtspp(max_files_gtspp), input_files_gtspp(max_files_gtspp))

    filetype = -1
    filetype_argo = -1
    filetype_gtspp = -1
    input_files = ''
    input_files_argo = ''
    input_files_gtspp = ''

    if ( prfs_obs .or. salt_obs .or. temp_obs_argo .or. temp_obs_gtspp .or. salt_obs_argo ) then
       ocn_obs%use_prf_as_obs = .true.
    end if
    ocn_obs%use_sst_as_obs = sst_obs
    ocn_obs%use_ssh_as_obs = eta_obs
    ocn_obs%use_suv_as_obs = suv_obs
    ocn_obs%use_woa05_t = temp_obs_woa05
    ocn_obs%use_woa05_s = salt_obs_woa05
    ssh_td = eta_obs_td

    ! time window for DROP, MOORING and SATELLITE data respectively
    ! will be available from namelist
    time_window(:) = set_time(0,unknown_window)
    time_window(DROP_PROFILER:DROP_PROFILER+9) = set_time(0,drop_window)
    time_window(MOORING:MOORING+9) = set_time(0,mooring_window)
    time_window(SATELLITE:SATELLITE+9) = set_time(0,satellite_window)
    time_window(DRIFTER:DRIFTER+9) = set_time(0,drifter_window)
    time_window(SHIP:SHIP+9) = set_time(0,ship_window)
    
    nfiles = 0
    nrecs=0
    call mpp_open(unit, 'ocean_obs_table', action=MPP_RDONLY)
    read_obs: do while ( nfiles <= max_files )
       read (UNIT=unit, FMT='(A)', IOSTAT=io_status) record
       if ( io_status < 0 ) then
          exit read_obs
       else if ( io_status > 0 ) then
          cycle read_obs
       else
          nrecs = nrecs + 1
          if ( record(1:1) == '#' ) cycle read_obs
          read ( UNIT=record, FMT=*, IOSTAT=io_status ) tbl_entry
          if ( io_status < 0 ) then
             exit read_obs
          else if ( io_status > 0 ) then
             cycle read_obs
          else
             nfiles = nfiles + 1       
             input_files(nfiles) = tbl_entry%filename
             select case ( trim(tbl_entry%file_type) )
             case ('profiles')
                filetype(nfiles) = PROFILE_FILE
             case ('sfc')
                filetype(nfiles) = SFC_FILE
             case default
                call error_mesg('oda_core_mod::init_observations', 'error in obs_table entry format', FATAL)
             end select
          end if
       end if
    end do read_obs
    if ( nfiles > max_files ) then
       call error_mesg('oda_core_mod::init_observations', 'number of obs files exceeds max_files parameter', FATAL)
    end if
    CALL mpp_close(unit)
    
    nfiles_argo = 0
    nrecs_argo = 0
    call mpp_open(unit_argo, 'ocean_obs_argo_table', action=MPP_RDONLY)
    read_obs_argo: do while ( nfiles_argo <= max_files_argo )
       read (UNIT=unit_argo, FMT='(A)', IOSTAT=io_status) record
       if ( io_status < 0 ) then
          exit read_obs_argo
       else if ( io_status > 0 ) then
          cycle read_obs_argo
       else 
          nrecs_argo = nrecs_argo + 1
          if ( record(1:1) == '#' ) cycle read_obs_argo
          read (UNIT=record, FMT=*, IOSTAT=io_status) tbl_entry
          if ( io_status < 0 ) then
             exit read_obs_argo
          else if ( io_status > 0 ) then
             cycle read_obs_argo
          else 
             nfiles_argo = nfiles_argo + 1       
             input_files_argo(nfiles_argo) = tbl_entry%filename
             select case ( trim(tbl_entry%file_type) )
             case ('profiles')
                filetype_argo(nfiles_argo) = PROFILE_FILE
             case ('sfc')
                filetype_argo(nfiles_argo) = SFC_FILE
             case default
                call error_mesg('oda_core_mod::init_observations', 'error in obs_table entry format for argo', FATAL)
             end select
          end if
       end if
    end do read_obs_argo
    if ( nfiles_argo > max_files_argo ) then
       call error_mesg('oda_core_mod::init_observations', 'number of obs files exceeds max_files_argo parameter', FATAL)
    end if
    call mpp_close(unit_argo)

    nfiles_gtspp = 0
    nrecs_gtspp = 0
    call mpp_open(unit_gtspp, 'ocean_obs_gtspp_table', action=MPP_RDONLY)
    read_obs_gtspp: do while ( nfiles_gtspp <= max_files_gtspp )
       read (UNIT=unit_gtspp, FMT='(A)', IOSTAT=io_status) record
       if ( io_status < 0 ) then
          exit read_obs_gtspp
       else if ( io_status > 0 ) then
          cycle read_obs_gtspp
       else
          nrecs_gtspp = nrecs_gtspp + 1
          if ( record(1:1) == '#' ) cycle read_obs_gtspp
          read (UNIT=record, FMT=*, IOSTAT=io_status) tbl_entry
          if ( io_status < 0 ) then
             exit read_obs_gtspp
          else if ( io_status > 0 ) then
             cycle read_obs_gtspp
          else
             nfiles_gtspp = nfiles_gtspp + 1       
             input_files_gtspp(nfiles_gtspp) = tbl_entry%filename
             select case ( trim(tbl_entry%file_type) )
             case ('profiles')
                filetype_gtspp(nfiles_gtspp) = PROFILE_FILE
             case ('sfc')
                filetype_gtspp(nfiles_gtspp) = SFC_FILE
             case default
                call error_mesg('oda_core_mod::init_observations', 'error in obs_table entry format for gtspp', FATAL)
             end select
          end if
       end if
    end do read_obs_gtspp
    if ( nfiles_gtspp > max_files_gtspp ) then
       call error_mesg('oda_core_mod::init_observations', 'number of obs files exceeds max_files_gtspp parameter', FATAL)
    end if
    CALL mpp_close(unit_gtspp)

    num_profiles = 0
    no_prf = 0
    no_sst = 0
    no_temp = 0
    no_salt = 0
    no_suv = 0
    no_eta = 0
    no_woa05 = 0

    ! get local indices for Model grid
    allocate(x_grid(isg:ieg,jsg:jeg), x_grid_uv(isg:ieg,jsg:jeg))
    allocate(y_grid(isg:ieg,jsg:jeg), y_grid_uv(isg:ieg,jsg:jeg))

    call mpp_global_field(filt_domain, Grd%x(:,:), x_grid(:,:))
    call mpp_global_field(filt_domain, Grd%y(:,:), y_grid(:,:))

    ! Allocate profiles
    allocate(profiles(max_profiles))

    do j=jsg, jeg
       do i=isg, ieg
          if ( x_grid(i,j) .lt. 80.0 ) x_grid(i,j) = x_grid(i,j) + 360.0
       end do
    end do

    ! uv grid may not be precise, need to be carefully checked
    x_grid_uv(:,:) = x_grid(:,:) + 0.5
    do j=jsg, jeg-1
       do i=isg, ieg
          y_grid_uv(i,j) = y_grid(i,j) + 0.5*(y_grid(i,j+1)-y_grid(i,j))
       end do
    end do
    do i=isg, ieg
       y_grid_uv(i,jeg) = 90.0
    end do

    if ( prfs_obs ) then
       obs_variable = TEMP_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("TEMP_ID = ",I5)') TEMP_ID
       end if
       do n=1, nfiles
          select case ( filetype(n) )
          case (PROFILE_FILE)
             call open_profile_dataset(trim(input_files(n)), time_s, time_e, obs_variable, localize)           
          case default
             call error_mesg('oda_core_mod::init_observations', 'filetype not currently supported for prfs_obs', FATAL)
          end select
       end do
    end if

    if ( salt_obs ) then
       obs_variable = SALT_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("SALT_ID = ",I5)') SALT_ID
       end if
       do n=1, nfiles
          select case ( filetype(n) )
          case (PROFILE_FILE)
             call open_profile_dataset(trim(input_files(n)), time_s, time_e, obs_variable, localize)
          case default
             call error_mesg('oda_core_mod::init_observations', 'filetype not currently supported for salt_obs', FATAL)
          end select
       end do
    end if

    if ( temp_obs_gtspp ) then
       obs_variable = TEMP_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("TEMP_ID = ",I5)') TEMP_ID
       end if
       do n=1, nfiles_gtspp
          if ( mpp_pe() == mpp_root_pe() ) then
             write (UNIT=stdout_unit, FMT='("f_typ_gtspp = ",I8)') filetype_gtspp(n)
             write (UNIT=stdout_unit, FMT='("i_f_gtspp = ",A)') input_files_gtspp(n)
          end if
          select case ( filetype_gtspp(n) )
          case (PROFILE_FILE)
             call open_profile_dataset_gtspp()
          case default
             call error_mesg('oda_core_mod::init_observations', 'filetype_gtspp not currently supported', FATAL)
          end select
       end do
    end if

    if ( temp_obs_argo ) then
       obs_variable = TEMP_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("TEMP_ID = ",I5)') TEMP_ID
       end if
       do n=1, nfiles_argo
          if ( mpp_pe() == mpp_root_pe() ) then
             write (UNIT=stdout_unit, FMT='("f_typ_argo = ",I8)') filetype_argo(n)
             write (UNIT=stdout_unit, FMT='("i_f_argo = ",A)') input_files_argo(n)
          end if
          select case ( filetype_argo(n) )
          case (PROFILE_FILE)
             call open_profile_dataset_argo(trim(input_files_argo(n)), time_s, time_e, obs_variable, localize)
          case default
             call error_mesg('oda_core_mod::init_observations', 'filetype_argo not currently supported', FATAL)
          end select
       end do
    end if

    if ( salt_obs_argo ) then
       obs_variable = SALT_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("SALT_ID = ",I5)') SALT_ID
       end if
       do n=1, nfiles_argo
          select case ( filetype_argo(n) )
          case (PROFILE_FILE)
             call open_profile_dataset_argo(trim(input_files_argo(n)), time_s, time_e, obs_variable, localize)
          case default
             call error_mesg('oda_core_mod::init_observations', 'filetype_argo not currently supported', FATAL)
          end select
       end do
    end if

    if ( temp_obs_woa05 ) then
       obs_variable = WOAT_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("WOAT_ID = ",I5)') WOAT_ID
       end if
       input_files_woa05t = "INPUT/woa05_temp.nc"
       call open_profile_dataset_woa05t(trim(input_files_woa05t), obs_variable, localize)
    end if

    if ( salt_obs_woa05 ) then
       obs_variable = WOAS_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("WOAS_ID = ",I5)') WOAS_ID
       end if
       input_files_woa05s = "INPUT/woa05_salt.nc"
       call open_profile_dataset_woa05s(trim(input_files_woa05s), obs_variable, localize)
    end if

    if ( sst_obs ) then
       obs_variable = TEMP_ID
       if ( mpp_pe() == mpp_root_pe() ) then
          write (UNIT=stdout_unit, FMT='("TEMP_ID for sst = ",I5)') TEMP_ID
       end if
       nfiles = 1
       nrecs = 0
       input_files(nfiles) = "INPUT/sst_daily.nc"
       filetype(nfiles) = PROFILE_FILE
!!$       call open_profile_dataset_sst(trim(input_files(nfiles)), obs_variable, localize)
    end if

    if ( eta_obs ) then
       obs_variable = ETA_ID
       nfiles = 1
       nrecs = 0
       input_files(nfiles) = "INPUT/ocean.19760101-20001231.eta_t.nc"
       filetype(nfiles) = PROFILE_FILE
       call open_profile_dataset_eta(trim(input_files(nfiles)), obs_variable, localize)
    end if

    if ( suv_obs ) then
       obs_variable = SUV_ID
       nfiles = 1
       nrecs = 0
       input_files(nfiles) = "INPUT/sfc_current.197601-200012.nc"
       filetype(nfiles) = PROFILE_FILE
       call open_profile_dataset_suv(trim(input_files(nfiles)), obs_variable, localize)
    end if

    ! Deallocate before exiting routine
    deallocate(filetype, input_files)
    deallocate(filetype_argo, input_files_argo)
    deallocate(filetype_gtspp, input_files_gtspp)
  end subroutine init_observations

  subroutine open_profile_dataset(filename, time_start, time_end, obs_variable, localize)  
    character(len=*), intent(in) :: filename
    type(time_type), intent(in) :: time_start, time_end
    integer, intent(in) :: obs_variable
    logical, intent(in), optional :: localize

    integer, parameter :: MAX_LEVELS = 1000
    integer, parameter :: MAX_LNKS = 500

    real :: lon, lat, time, profile_error, rlink, flag_t, flag_s, fix_depth
    real :: ri0, rj0
    real, dimension(MAX_LEVELS) :: depth, data, t_flag, s_flag
    real, dimension(MAX_LNKS, MAX_LEVELS) :: data_bfr, depth_bfr, t_flag_bfr, s_flag_bfr

    integer :: unit, ndim, nvar, natt, nstation
    integer :: stdout_unit
    integer :: inst_type, var_id
    integer :: num_levs, k, kk, i, i0, j0, k0, nlevs, a, nn, ii, nlinks
    integer :: nprof_in_filt_domain
    integer :: bad_point, bad_point_g, out_bound_point

    logical :: data_is_local, localize_data, cont
    logical :: data_in_period
    logical :: prof_in_filt_domain
    logical, dimension(MAX_LEVELS) :: flag
    logical, dimension(MAX_LNKS, MAX_LEVELS) :: flag_bfr

    character(len=32) :: fldname, axisname, anal_fldname, time_units
    character(len=138) :: emsg_local

    type(time_type) :: profile_time
    type(axistype), pointer :: depth_axis, station_axis
    type(axistype), allocatable, dimension(:), target :: axes
    type(fieldtype), allocatable, dimension(:), target :: fields
    type(fieldtype), pointer :: field_lon, field_lat, field_flag, field_time, field_depth, field_data
    type(fieldtype), pointer :: field_error, field_link, field_t_flag, field_s_flag, field_fix_depth ! snz drop rate

    ! NOTE: fields are restricted to be in separate files

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    nprof_in_filt_domain = 0
    stdout_unit = stdout()

    anal_fldname = 'none'
    var_id=-1
    if ( obs_variable == TEMP_ID ) then
       anal_fldname = 'temp'
       var_id = TEMP_ID
    else if ( obs_variable == SALT_ID ) then
       anal_fldname = 'salt'
       var_id = SALT_ID
    end if

!    call mpp_print_memuse_stats('open_profile_dataset Start')

    call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_MULTI, action=MPP_RDONLY)
    call mpp_get_info(unit, ndim, nvar, natt, nstation)

    if ( mpp_pe() .eq. mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("Opened profile dataset: ",A)') trim(filename)
    end if

    ! get axis information
    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case ( trim(axisname) )
       case ('depth_index')
          depth_axis => axes(i)
       case ('station_index')
          station_axis => axes(i)
       end select
    end do

    ! get field information
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields)
    do i=1, nvar
       call mpp_get_atts(fields(i), name=fldname)
       if( var_id .eq. TEMP_ID ) then
          select case (trim(fldname))
          case ('longitude')
             field_lon => fields(i)
          case ('latitude')
             field_lat => fields(i)
          case ('profile_flag') 
             field_flag => fields(i)
          case ('time')
             field_time => fields(i)
          case ('temp')
             field_data => fields(i)
          case ('depth')
             field_depth => fields(i)
          case ('link')
             field_link => fields(i)
          case ('temp_error')
             field_error => fields(i)
          case ('temp_flag')
             field_t_flag => fields(i)
          case ('fix_depth') ! snz drop rate
             field_fix_depth => fields(i)
          end select
       else if( var_id .eq. SALT_ID ) then
          select case (trim(fldname))
          case ('longitude')
             field_lon => fields(i)
          case ('latitude')
             field_lat => fields(i)
          case ('profile_flag_s') 
             field_flag => fields(i)
          case ('time')
             field_time => fields(i)
          case ('salt')
             field_data => fields(i)
          case ('depth')
             field_depth => fields(i)
          case ('link')
             field_link => fields(i)
          case ('salt_error')
             field_error => fields(i)
          case ('salt_flag')
             field_s_flag => fields(i)
          case ('fix_depth') ! snz drop rate
             field_fix_depth => fields(i)
          end select
       end if
    end do

    call mpp_get_atts(depth_axis, len=nlevs)

    if ( nlevs > MAX_LEVELS ) then 
       call error_mesg('oda_core_mod::open_profile_dataset', 'increase parameter MAX_LEVELS', FATAL)
    else if (nlevs < 1) then
       call error_mesg('oda_core_mod::open_profile_dataset', 'Value of nlevs is less than 1.', FATAL)
    end if

    if ( .NOT.ASSOCIATED(field_data) ) then
       call error_mesg('oda_core_mod::open_profile_dataset',&
            & 'profile dataset not used because data not needed for Analysis', NOTE)
       return
    end if

    write(UNIT=stdout_unit, FMT='("There are ",I8," records in this dataset.")') nstation
    write(UNIT=stdout_unit, FMT='("Searching for profiles . . .")')
    
    call mpp_get_atts(field_time, units=time_units)

    bad_point = 0
    out_bound_point = 0

    i = 1
    cont = .true.
    
    do while ( cont )
       prof_in_filt_domain = .false.
       depth = missing_value  ! snz add
       data = missing_value   ! snz add

       call mpp_read(unit, field_lon, lon, tindex=i)
       call mpp_read(unit, field_lat, lat, tindex=i)
       call mpp_read(unit, field_time, time, tindex=i)
       call mpp_read(unit, field_depth, depth(1:nlevs), tindex=i)
       call mpp_read(unit, field_data, data(1:nlevs), tindex=i)
       call mpp_read(unit, field_error, profile_error, tindex=i)
       call mpp_read(unit, field_fix_depth, fix_depth, tindex=i) ! snz drop rate
       if ( var_id == TEMP_ID ) then
          call mpp_read(unit, field_t_flag, t_flag(1:nlevs), tindex=i)
          call mpp_read(unit, field_flag, flag_t, tindex=i)
       else if ( var_id == SALT_ID ) then
          call mpp_read(unit, field_s_flag, s_flag(1:nlevs), tindex=i)
          call mpp_read(unit, field_flag, flag_s, tindex=i)
       end if
       call mpp_read(unit, field_link, rlink, tindex=i)

       inst_type = 20 ! snz change one line
!!$       inst_type = DRIFTER + ARGO

       data_is_local = .false.
       data_in_period = .false.

       if ( lon .lt. 0.0 ) lon = lon + 360.0
       if ( lon .gt. 360.0 ) lon = lon - 360.0
       if ( lon .lt. 80.0 ) lon = lon + 360.0 

       if ( lat > ass_start_lat .and. lat < ass_end_lat ) data_is_local = .true.

       profile_time = get_cal_time(time, time_units, 'julian')
       if ( profile_time > time_start .and. profile_time < time_end ) data_in_period = .true.
       if ( (data_in_period .and. data_is_local) .and. (.NOT.localize_data) ) then ! localize

          if (isd_filt < 1 .and. ied_filt > ieg) then
             ! filter domain is a full x band
             if (lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  lat >= y_grid(1,jsd_flt0) .and. lat <= y_grid(ieg-1,jsd_flt0)) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt <= ieg) then
             ! Interior filter domain
            if (lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
               & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)) then
                prof_in_filt_domain = .true.
            end if
          else if (isd_filt < 1 .and. ied_filt <= ieg) then
             ! lhs filter domain
             isd_flt0 = isd_filt + ieg
             if ((lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
                  & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)).or.&
                  & (lon >= x_grid(isd_flt0,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_flt0,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1))) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt > ieg) then
             ! rhs filter domain
             ied_flt0 = ied_filt - ieg
             if ( lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1) ) then
                prof_in_filt_domain = .true.
             end if
             if (ied_flt0-1 > 1) then
                if ( lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_flt0-1,jsd_flt0) .and.&
                     & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_flt0-1,jed_flt0-1) ) then
                   prof_in_filt_domain = .true.
                end if
             end if
          end if

          if ( var_id == TEMP_ID .and. flag_t == 0.0 ) then
             num_profiles = num_profiles + 1
             no_temp = no_temp + 1
             no_prf = no_prf + 1
          end if
          if ( var_id == SALT_ID .and. flag_s == 0.0 ) then 
             num_profiles = num_profiles + 1
             no_salt = no_salt + 1
             no_prf = no_prf + 1
          end if

          if ( num_profiles > max_profiles ) then
             call error_mesg('oda_core_mod::open_profile_dataset',&
                  & 'maximum number of profiles exceeded, increase max_profiles in oda_core_nml', FATAL)
          end if

          num_levs = 0
          do k=1, MAX_LEVELS
             flag(k) = .true.

             if ( depth(k) > 2000.0 ) depth(k) = missing_value  ! snz add for rdat-hybn

             if ( var_id == TEMP_ID ) then
                if ( data(k) .eq. missing_value .or.&
                     & depth(k) .eq. missing_value .or. t_flag(k) .ne. 0.0 ) then
                   flag(k) = .false.
                else
                   num_levs = num_levs + 1
                end if
             else if ( var_id == SALT_ID ) then
                if ( data(k) .eq. missing_value .or.&
                     & depth(k) .eq. missing_value .or. s_flag(k) .ne. 0.0 ) then
                   flag(k) = .false.
                else
                   num_levs = num_levs+1
                end if
             end if
          end do

          ! large profile are stored externally in separate records
          ! read linked records and combine profile
          ii = i + 1
          nlinks = 0
          do while ( rlink > 0.0 )
             nlinks = nlinks + 1

             if ( nlinks > MAX_LNKS ) then
                write (emsg_local, '("nlinks (",I6,") > MAX_LNKS (",I6,")")')&
                     & nlinks, MAX_LNKS
                call error_mesg('oda_core_mod::open_profile_dataset',&
                     & trim(emsg_local)//' in file "'//trim(filename)//&
                     & '".  Increase parameter MAX_LNKS', FATAL)
             end if

             depth_bfr(nlinks,:) = missing_value
             data_bfr(nlinks,:) = missing_value
             call mpp_read(unit, field_depth, depth_bfr(nlinks,1:nlevs), tindex=ii)
             call mpp_read(unit, field_data, data_bfr(nlinks,1:nlevs), tindex=ii)
             if ( var_id == TEMP_ID ) then 
                call mpp_read(unit, field_t_flag, t_flag_bfr(nlinks,1:nlevs), tindex=ii)
             else if ( var_id == SALT_ID ) then
                call mpp_read(unit, field_s_flag, s_flag_bfr(nlinks,1:nlevs), tindex=ii)
             end if
             call mpp_read(unit, field_link, rlink, tindex=ii)
             ii = ii + 1
          end do
          i = ii ! set record counter to start of next profile

          if ( nlinks > 0 ) then
             do nn=1, nlinks
                do k=1, MAX_LEVELS
                   flag_bfr(nn,k) = .true.

                   if ( depth_bfr(nn,k) > 2000.0 ) depth_bfr(nn,k) = missing_value  ! snz add for rdat-hybn

                   if ( var_id == TEMP_ID ) then
                      if ( data_bfr(nn,k) .eq. missing_value .or.&
                           & depth_bfr(nn,k) .eq. missing_value .or.&
                           & t_flag_bfr(nn,k) .ne. 0.0 ) then
                         flag_bfr(nn,k) = .false.
                      else
                         num_levs = num_levs+1
                      end if
                   else if (var_id == SALT_ID) then
                      if ( data_bfr(nn,k) .eq. missing_value .or.&
                           & depth_bfr(nn,k) .eq. missing_value .or.&
                           & s_flag_bfr(nn,k) .ne. 0.0 ) then
                         flag_bfr(nn,k) = .false.
                      else
                         num_levs = num_levs+1
                      end if
                   end if
                end do
             end do
          end if

          ! mh2 asks to change from [if (num_levs == 0) cycle]
          if ( num_levs == 0 ) then
             if ( i .gt. nstation ) cont = .false.
             cycle
          end if

          if ( num_profiles > 0 .and. prof_in_filt_domain ) then ! snz - 05 Nov 2012

             allocate(profiles(num_profiles)%depth(num_levs))
             allocate(profiles(num_profiles)%data(num_levs))
             allocate(profiles(num_profiles)%flag(num_levs))
             profiles(num_profiles)%variable = var_id
             if ( inst_type < 1 ) inst_type = UNKNOWN
             profiles(num_profiles)%inst_type = inst_type
             profiles(num_profiles)%levels = num_levs
             profiles(num_profiles)%lat = lat
             profiles(num_profiles)%lon = lon
!             allocate(profiles(num_profiles)%ms(num_levs))
!             allocate(profiles(num_profiles)%ms_inv(num_levs))           
!             profiles(num_profiles)%ms(:) = 0.5
             kk = 1
             do k=1, MAX_LEVELS
                if ( flag(k) ) then
                   if ( kk > profiles(num_profiles)%levels ) then
                      call error_mesg('oda_core_mod::open_profile_dataset',&
                           & 'Loop value "kk" is greater than profile levels', FATAL)
                   end if
                   profiles(num_profiles)%depth(kk) = depth(k)
                   profiles(num_profiles)%data(kk) = data(k)
!                   profiles(num_profiles)%ms_inv(kk) = 1./profiles(num_profiles)%ms(kk)                  
                   kk = kk + 1
                end if
             end do

             do nn=1, nlinks
                do k=1, MAX_LEVELS
                   if ( flag_bfr(nn,k) ) then
                      if ( kk > profiles(num_profiles)%levels ) then
                         call error_mesg('oda_core_mod::open_profile_dataset',&
                              & 'Loop value "kk" is greater than profile levels (bfr loop)', FATAL)
                      end if
                      profiles(num_profiles)%depth(kk) = depth_bfr(nn,k)
                      profiles(num_profiles)%data(kk) = data_bfr(nn,k)
!                      profiles(num_profiles)%ms_inv(kk) = 1./profiles(num_profiles)%ms(kk)
                      kk = kk + 1
                   end if
                end do
             end do
           
             profiles(num_profiles)%time = profile_time
           
             ! calculate interpolation coefficients (make sure to account for grid offsets here!)
             if ( lat < 65.0 ) then ! regular grids
                ri0 = frac_index(lon, x_grid(:,1))
                rj0 = frac_index(lat, y_grid(90,:))
                i0 = floor(ri0)
                j0 = floor(rj0)
                if ( i0 > ieg .or. j0 > jeg ) then
                   write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') mpp_pe(), i0, j0
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
                end if
                if ( isd_filt >= 1 .and. ied_filt <= ieg ) then
                  if ( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
                     write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                          & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                     call error_mesg('oda_core_mod::open_profile_dataset',&
                          & 'i0,j0 out of bounds in prfs01. '//trim(emsg_local), FATAL)
                  end if
                end if
                if ( isd_filt < 1 .and. i0 > ied_filt-1 .and. i0 < isd_filt + ieg ) then
                   write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'i0,j0 out of bounds in prfs02. '//trim(emsg_local), FATAL)
                end if
                if ( ied_filt > ieg .and. i0 > ied_filt-ieg-1 .and. ied_filt < isd_filt ) then
                   write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'i0,j0 out of bounds in prfs03. '//trim(emsg_local), FATAL)
                end if
                Profiles(num_profiles)%i_index = ri0
                Profiles(num_profiles)%j_index = rj0
             else ! tripolar grids
                lon_out(1,1) = (lon-360.0)*DEG_TO_RAD
                lat_out(1,1) = lat*DEG_TO_RAD
                call horiz_interp_bilinear_new (Interp, (x_grid-360.0)*DEG_TO_RAD, y_grid*DEG_TO_RAD,&
                     & lon_out, lat_out, new_search=.true., no_crash_when_not_found=.true.)

                if ( Interp%i_lon(1,1,1) == -999. ) bad_point = bad_point + 1
                if ( Interp%wti(1,1,2) < 1.0 ) then
                   i0 = Interp%i_lon(1,1,1)
                else
                   i0 = Interp%i_lon(1,1,2)
                end if
                if ( Interp%wtj(1,1,2) < 1.0 ) then
                   j0 = Interp%j_lat(1,1,1)
                else
                   j0 = Interp%j_lat(1,1,2)
                end if
                if ( i0 > ieg .or. j0 > jeg ) then
                   write (UNIT=emsg_local, FMT='("i0 = ",I6,", j0 = ",I6)') mpp_pe(), i0, j0
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'For tripolar grids, either i0 > ieg or j0 > jeg', FATAL)
                end if
                if( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
!!$                   print*,'prfs.pe,i0,j0= ',mpp_pe(), i0, j0,&
!!$                        & 'isd_filt,ied_filt,jsd_filt,jed_filt= ',isd_filt,ied_filt,jsd_filt,jed_filt
!!$                   print*,'pe,lon,lat=',mpp_pe(),lon,lat,'x_grid(i0+-1)',x_grid(i0-1:i0+1,j0),&
!!$                        & 'y_grid(i0,j0+-1)=',y_grid(i0,j0-1:j0+1)
!!$                   print*,'lono11,lato11=',x_grid(i0,j0),y_grid(i0,j0),'lono21,lato21=',x_grid(i0+1,j0),y_grid(i0+1,j0)
!!$                   print*,'lono12,lato12=',x_grid(i0,j0+1),y_grid(i0,j0+1),'lono22,lato22=',x_grid(i0+1,j0+1),y_grid(i0+1,j0+1)
!!$                   print*,'lonm11,latm11=',x_grid(isd_filt,jsd_filt),y_grid(isd_filt,jsd_filt),&
!!$                        & 'lonm21,latm21=',x_grid(ied_filt,jsd_filt),y_grid(ied_filt,jsd_filt)
!!$                   print*,'lonm12,latm12=',x_grid(isd_filt,jed_filt),y_grid(isd_filt,jed_filt),&
!!$                        & 'lonm22,latm22=',x_grid(ied_filt,jed_filt),y_grid(ied_filt,jed_filt)
!!$                   print*,'wti(1:2)=',Interp%wti(1,1,:),'wtj(1:2)=',Interp%wtj(1,1,:)

                   out_bound_point = out_bound_point + 1
                end if
                if ( Interp%wti(1,1,2) < 1.0 ) then
                   Profiles(num_profiles)%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
                else
                   Profiles(num_profiles)%i_index =Interp%i_lon(1,1,2)
                end if
                if (Interp%wtj(1,1,2) < 1.0) then
                   Profiles(num_profiles)%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
                else
                   Profiles(num_profiles)%j_index =Interp%j_lat(1,1,2)
                end if
             end if ! grids

             Profiles(num_profiles)%accepted = .true.

             if ( var_id == TEMP_ID .and. flag_t /= 0.0 ) Profiles(num_profiles)%accepted = .false.
             if ( var_id == SALT_ID .and. flag_s /= 0.0 ) Profiles(num_profiles)%accepted = .false.

             if (i0 < 1 .or. j0 < 1) then
                Profiles(num_profiles)%accepted = .false.
             end if
             if( i0 < isd_filt .or. i0 >= ied_filt .or. j0 < jsd_filt .or. j0 >= jed_filt ) then
                Profiles(num_profiles)%accepted = .false.
             end if

             if ( Profiles(num_profiles)%accepted ) then ! here
                if ( i0 /= ieg .and. j0 /= jeg ) then
                   if (Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else if ( i0 == ieg .and. j0 /= jeg ) then
                   if (Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else if ( i0 /= ieg .and. j0 == jeg ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or. Grd%mask(i0+1,j0,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else
                   if ( Grd%mask(i0,j0,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                end if
             end if ! here

             if ( Profiles(num_profiles)%accepted .and.&
                  & Profiles(num_profiles)%inst_type == MOORING+TAO ) then
                if ( allocated(mask_tao) ) then
                   if ( mask_tao(i0,j0) < 1.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                      write (UNIT=stdout_unit,&
                           & FMT='("Rejecting tao mooring at (lat,lon) = (",F10.5,",",F10.5,") based on user-specified mask.")')&
                           & Profiles(num_profiles)%lat,&
                           & Profiles(num_profiles)%lon
                   end if
                end if
             end if

             if ( Profiles(num_profiles)%accepted ) then ! accepted
                Profiles(num_profiles)%flag(:) = .true.
                allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
                do k=1, Profiles(num_profiles)%levels
                   if (Profiles(num_profiles)%depth(k) < Grd%z(1)) then
                     Profiles(num_profiles)%k_index(k) = 1.0
                   else
                     Profiles(num_profiles)%k_index(k) = frac_index(Profiles(num_profiles)%depth(k), (/0.,Grd%z(:)/))! - 1 snz modify to v3.2 JAN3012
                   end if
                   if ( Profiles(num_profiles)%k_index(k) < 1.0 ) then
                      if ( Profiles(num_profiles)%depth(k) < 0.0 ) then
                         Profiles(num_profiles)%k_index(k) = 0.0
                      else if ( Profiles(num_profiles)%depth(k) > Grd%z(size(Grd%z,1)) ) then
                         Profiles(num_profiles)%k_index(k) = real(nk)
                      end if
                   else
                      Profiles(num_profiles)%k_index(k) = Profiles(num_profiles)%k_index(k) - 1.0
                   end if
                   if ( Profiles(num_profiles)%k_index(k) > real(nk) ) then
                      call error_mesg('oda_core_mod::open_profile_dataset', 'Profile k_index is greater than nk', FATAL)
                   else if ( Profiles(num_profiles)%k_index(k) < 0.0 ) then
                      call error_mesg('oda_core_mod::open_profile_dataset', 'Profile k_index is less than 0', FATAL)
                   end if
                   k0 = floor(Profiles(num_profiles)%k_index(k))

                   IF ( k0 >= 1 ) THEN ! snz add
                      if ( Profiles(num_profiles)%flag(k) ) then ! flag
                         if ( i0 /= ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0+1,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 == ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(1,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                                 & Grd%mask(1,j0+1,k0) == 0.0) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 /= ieg .and. j0 == jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0) == 0.0) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else
                            if ( Grd%mask(i0,j0,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         end if

                         if ( i0 /= ieg .and. j0 /= jeg) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0+1,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 == ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(1,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                                 & Grd%mask(1,j0+1,k0+1) == 0.0) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 /= ieg .and. j0 == jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0+1) == 0.0) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         end if
                        
                         if ( abs(Profiles(num_profiles)%data(k)) > 1.e4 &
                              & .or. abs(Profiles(num_profiles)%depth(k)) > 1.e4 ) then
                            Profiles(num_profiles)%flag(k) = .false.
                         end if
                      end if ! flag
                   end if ! snz add
                end do
             end if ! accepted         
          endif ! 05 Nov 2012
       else ! localize
          i = i+1
       end if ! localize

       if ( var_id == TEMP_ID .and. num_profiles > 0 ) call xbt_drop_rate_adjust(Profiles(num_profiles))

       if ( i .gt. nstation ) cont = .false.
    end do

    a = nprof_in_filt_domain
    bad_point_g = bad_point
    call mpp_sum(a)
    call mpp_sum(bad_point_g)
    call mpp_sum(out_bound_point)

    if ( no_prf /= num_profiles ) then
       write(UNIT=stdout_unit, FMT='("PE: ",I6," no_prf = ",I8,", num_profiles = ",I8)') mpp_pe(), no_prf, num_profiles
    end if
    if ( var_id == TEMP_ID ) then
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," temp prfs within global domain")') no_temp
       write(UNIT=stdout_unit, FMT='("The total of bad_point temp ",I8," within global domain")') bad_point_g
       write(UNIT=stdout_unit, FMT='("The total out_bound_point temp ",I8)') out_bound_point
    else if ( var_id == SALT_ID ) then
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," salt prfs within global domain")') no_salt
       write(UNIT=stdout_unit, FMT='("The total of bad_point salt",I8," within global domain")') bad_point_g
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," prfs within global domain")') no_prf
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," prfs within current PEs computer domain")') a
       write(UNIT=stdout_unit, FMT='("The total out_bound_point salt ",I8)') out_bound_point
    end if

    call mpp_sync_self()
    call mpp_close(unit)
    deallocate(axes)
    deallocate(fields)

!    call mpp_print_memuse_stats('open_profile_dataset End')

  end subroutine open_profile_dataset

  subroutine open_profile_dataset_gtspp()
    return
  end subroutine open_profile_dataset_gtspp

  subroutine open_profile_dataset_argo(filename, time_start, time_end, obs_variable, localize)  
    character(len=*), intent(in) :: filename
    type(time_type), intent(in) :: time_start, time_end
    integer, intent(in) :: obs_variable
    logical, intent(in), optional :: localize

    integer, parameter :: MAX_LEVELS = 1000
    integer, parameter :: MAX_LNKS = 500

    real :: lon, lat, time, rlink, prf_type
    real :: ri0, rj0
    real, dimension(MAX_LEVELS) :: depth, data
    real, dimension(MAX_LNKS, MAX_LEVELS) :: data_bfr, depth_bfr

    integer :: unit, ndim, nvar, natt, nstation
    integer :: stdout_unit
    integer :: inst_type, var_id
    integer :: num_levs, k, kk, i, i0, j0, k0, nlevs, a, nn, ii, nlinks
    integer :: nprof_in_filt_domain, out_bound_point

    character(len=32) :: fldname, axisname, anal_fldname, time_units
    character(len=128) :: emsg_local

    logical :: data_is_local, localize_data, cont
    logical :: data_in_period
    logical :: prof_in_filt_domain
    logical, dimension(MAX_LEVELS) :: flag 
    logical, dimension(MAX_LNKS, MAX_LEVELS) :: flag_bfr

    type(time_type) :: profile_time
    type(axistype), pointer :: depth_axis, station_axis
    type(axistype), allocatable, dimension(:), target :: axes
    type(fieldtype), allocatable, dimension(:), target :: fields
    type(fieldtype), pointer :: field_lon, field_lat, field_flag, field_time
    type(fieldtype), pointer :: field_depth, field_data, field_link, field_var_type

    ! NOTE: fields are restricted to be in separate files

    stdout_unit = stdout()

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    nprof_in_filt_domain = 0

    anal_fldname = 'none'
    var_id=-1
    if ( obs_variable == TEMP_ID ) then
       anal_fldname = 'temp'
       var_id = TEMP_ID
    else if ( obs_variable == SALT_ID ) then
       anal_fldname = 'salt'
       var_id = SALT_ID
    end if

!    call mpp_print_memuse_stats('open_profile_dataset_argo Start')

    call mpp_open(unit, filename, form=MPP_NETCDF, fileset=MPP_SINGLE, threading=MPP_MULTI, action=MPP_RDONLY)
    call mpp_get_info(unit, ndim, nvar, natt, nstation)

    write (UNIT=stdout_unit, FMT='("Opened profile dataset: ",A)') trim(filename)

    ! get axis information

    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case (trim(axisname))
       case ('depth_index')
          depth_axis => axes(i)
       case ('station_index')
          station_axis => axes(i)
       end select
    end do

    ! get field information
    allocate(fields(nvar))
    call mpp_get_fields(unit, fields)
    do i=1, nvar
       call mpp_get_atts(fields(i), name=fldname)
       if( var_id .eq. TEMP_ID ) then
          select case (trim(fldname))
          case ('longitude')
             field_lon => fields(i)
          case ('latitude')
             field_lat => fields(i)
          case ('dens_flag') 
             field_flag => fields(i)
          case ('time')
             field_time => fields(i)
          case ('temp')
             field_data => fields(i)
          case ('depth')
             field_depth => fields(i)
          case ('link')
             field_link => fields(i)
          case ('var_type')
             field_var_type => fields(i)
          end select
       else if( var_id .eq. SALT_ID ) then
          select case (trim(fldname))
          case ('longitude')
             field_lon => fields(i)
          case ('latitude')
             field_lat => fields(i)
          case ('dens_flag') 
             field_flag => fields(i)
          case ('time')
             field_time => fields(i)
          case ('salt')
             field_data => fields(i)
          case ('depth')
             field_depth => fields(i)
          case ('link')
             field_link => fields(i)
          case ('var_type')
             field_var_type => fields(i)
          end select
       end if
    end do

    call mpp_get_atts(depth_axis, len=nlevs)

    if ( nlevs > MAX_LEVELS ) then 
       call error_mesg('oda_core_mod::open_profile_dataset_argo', 'increase parameter MAX_LEVELS', FATAL)
    else if (nlevs < 1) then
       call error_mesg('oda_core_mod::open_profile_dataset_argo', 'nlevs less than 1.', FATAL)
    end if

    if ( .NOT.ASSOCIATED(field_data) ) then
       call error_mesg('oda_core_mod::open_profile_dataset_argo',&
            & 'profile dataset not used because data not needed for Analysis', NOTE)
       return
    end if

    write (UNIT=stdout_unit, FMT='("There are ",I8," records in this dataset")') nstation
    write (UNIT=stdout_unit, FMT='("Searching for profiles . . .")')
    
    call mpp_get_atts(field_time, units=time_units)

    out_bound_point = 0

    i=1
    cont=.true.
    
    do while (cont)
       prof_in_filt_domain = .false.
       depth = missing_value  ! snz add
       data = missing_value   ! snz add

       call mpp_read(unit, field_lon, lon, tindex=i)
       call mpp_read(unit, field_lat, lat, tindex=i)
       call mpp_read(unit, field_time, time, tindex=i)
       call mpp_read(unit, field_depth, depth(1:nlevs), tindex=i)
       call mpp_read(unit, field_data, data(1:nlevs), tindex=i)
       call mpp_read(unit, field_var_type, prf_type, tindex=i)
       call mpp_read(unit, field_link, rlink, tindex=i)

!!$       inst_type = DRIFTER + ARGO
       inst_type = 20 ! snz change one line
       data_is_local = .false.
       data_in_period = .false.

       if ( lon .lt. 0.0 ) lon = lon + 360.0
       if ( lon .gt. 360.0 ) lon = lon - 360.0
       if ( lon .lt. 80.0 ) lon = lon + 360.0 

       if ( lat > ass_start_lat .and. lat < ass_end_lat ) data_is_local = .true.

       profile_time = get_cal_time(time, time_units, 'NOLEAP')
       if ( profile_time > time_start .and. profile_time < time_end ) data_in_period = .true.
       if ( (data_in_period .and. data_is_local) .and. (.NOT.localize_data) ) then

          if (isd_filt < 1 .and. ied_filt > ieg) then
             ! filter domain is a full x band
             if (lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  lat >= y_grid(1,jsd_flt0) .and. lat <= y_grid(ieg-1,jsd_flt0)) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt <= ieg) then
             ! Interior filter domain
            if (lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
               & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)) then
                prof_in_filt_domain = .true.
            end if
          else if (isd_filt < 1 .and. ied_filt <= ieg) then
             ! lhs filter domain
             isd_flt0 = isd_filt + ieg
             if ((lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
                  & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)).or.&
                  & (lon >= x_grid(isd_flt0,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_flt0,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1))) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt > ieg) then
             ! rhs filter domain
             ied_flt0 = ied_filt - ieg
             if ( lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1) ) then
                prof_in_filt_domain = .true.
             end if
             if (ied_flt0-1 > 1) then
                if ( lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_flt0-1,jsd_flt0) .and.&
                     & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_flt0-1,jed_flt0-1) ) then
                   prof_in_filt_domain = .true.
                end if
             end if
          end if

          if ( var_id == TEMP_ID ) then
             num_profiles = num_profiles + 1
             no_temp = no_temp + 1
             no_prf = no_prf + 1
          else if ( var_id == SALT_ID .and. prf_type == 2.0 ) then
             num_profiles =num_profiles + 1
             no_salt = no_salt + 1
             no_prf = no_prf + 1
          end if

          if ( num_profiles > max_profiles ) then
             call error_mesg('oda_core_mod::open_profile_dataset_argo',&
                  & 'maximum number of profiles exceeded, increase max_profiles in oda_core_nml', FATAL)
          end if

          num_levs = 0
          do k=1, MAX_LEVELS
             flag(k) = .true.
             if ( depth(k) > 2000.0 ) depth(k) = missing_value  ! snz add for rdat-hybn
             if ( var_id == TEMP_ID ) then
                if ( data(k) .eq. missing_value .or. depth(k) .eq. missing_value ) then
                   flag(k) = .false.
                else
                   num_levs = num_levs+1
                end if
             else if ( var_id == SALT_ID ) then
                if ( data(k) .eq. missing_value .or. depth(k) .eq. missing_value ) then
                   flag(k) = .false.
                else
                   num_levs = num_levs+1
                end if
             end if
          end do

          ! large profile are stored externally in separate records
          ! read linked records and combine profile
          ii=i+1
          nlinks = 0
          do while ( rlink > 0.0 )
             nlinks = nlinks + 1

             if ( nlinks > MAX_LNKS ) then
                write (emsg_local, '("nlinks (",I6,") > MAX_LNKS (",I6,")")')&
                     & nlinks, MAX_LNKS
                call error_mesg('oda_core_mod::open_profile_dataset_argo',&
                     & trim(emsg_local)//' in file "'//trim(filename)//&
                     & '".  Increase parameter MAX_LNKS', FATAL)
             end if

             depth_bfr(nlinks,:) = missing_value
             data_bfr(nlinks,:) = missing_value
             call mpp_read(unit,field_depth,depth_bfr(nlinks,1:nlevs),tindex=ii)
             call mpp_read(unit,field_data,data_bfr(nlinks,1:nlevs),tindex=ii)
             call mpp_read(unit,field_link,rlink,tindex=ii)
             ii=ii+1
          end do
          i=ii ! set record counter to start of next profile

          if ( nlinks > 0 ) then
             do nn=1, nlinks
                do k=1, MAX_LEVELS
                   flag_bfr(nn,k) = .true.

                   if ( depth_bfr(nn,k) > 2000.0 ) depth_bfr(nn,k) = missing_value  ! snz add for rdat-hybn

                   if ( var_id == TEMP_ID ) then
                      if ( data_bfr(nn,k) .eq. missing_value .or. depth_bfr(nn,k) .eq. missing_value ) then
                         flag_bfr(nn,k) = .false.
                      else
                         num_levs = num_levs+1
                      end if
                   else if ( var_id == SALT_ID ) then
                      if ( data_bfr(nn,k) .eq. missing_value .or. depth_bfr(nn,k) .eq. missing_value ) then
                         flag_bfr(nn,k) = .false.
                      else
                         num_levs = num_levs+1
                      end if
                   end if
                end do
             end do
          end if

          ! mh2 asks to change from [if (num_levs == 0) cycle]
          if ( num_levs == 0 ) then
             if ( i .gt. nstation ) cont = .false.
             cycle
          end if

          if (nprof_in_filt_domain > 0 .and. prof_in_filt_domain) then ! snz 05 Nov 2012

          allocate(profiles(num_profiles)%depth(num_levs))
          allocate(profiles(num_profiles)%data(num_levs))
          allocate(profiles(num_profiles)%flag(num_levs))
          profiles(num_profiles)%variable = var_id
          if ( inst_type < 1 ) inst_type = UNKNOWN
          profiles(num_profiles)%inst_type = inst_type
          profiles(num_profiles)%levels = num_levs
          profiles(num_profiles)%lat = lat
          profiles(num_profiles)%lon = lon
!          allocate(profiles(num_profiles)%ms(num_levs))
!          allocate(profiles(num_profiles)%ms_inv(num_levs))           
!          profiles(num_profiles)%ms(:) = 0.5

          kk= 1
          do k=1, MAX_LEVELS
             if ( flag(k) ) then
                if ( kk > profiles(num_profiles)%levels ) then
                   call error_mesg('oda_core_mod::open_profile_dataset_argo',&
                        & 'Loop variable "kk" is greater than profile levels', FATAL)
                end if
                profiles(num_profiles)%depth(kk) = depth(k)
                profiles(num_profiles)%data(kk) = data(k)
!                profiles(num_profiles)%ms_inv(kk) = 1./profiles(num_profiles)%ms(kk)                  
                kk = kk + 1
             end if
          end do

          do nn=1, nlinks
             do k=1, MAX_LEVELS
                if ( flag_bfr(nn,k) ) then
                   if ( kk > profiles(num_profiles)%levels ) then
                      call error_mesg('oda_core_mod::open_profile_dataset_argo',&
                           & 'Loop variable "kk" is greater than profile levels (bfr loop)', FATAL)
                   end if
                   profiles(num_profiles)%depth(kk) = depth_bfr(nn,k)
                   profiles(num_profiles)%data(kk) = data_bfr(nn,k)
!                   profiles(num_profiles)%ms_inv(kk) = 1./profiles(num_profiles)%ms(kk)
                   kk = kk + 1
                end if
             end do
          end do
           
          profiles(num_profiles)%time = profile_time
           
! snz uses the following to test excluding the coast area salt profiles
!          if (profiles(num_profiles)%variable == SALT_ID .and. &
!              profiles(num_profiles)%depth(num_levs) < 900.0) profiles(num_profiles)%accepted = .false.

          ! calculate interpolation coefficients (make sure to account for grid offsets here!)
          ! note that this only works for lat/lon grids
          if ( lat < 65.0 ) then ! regular grids
             ri0 = frac_index(lon, x_grid(:,1))
             rj0 = frac_index(lat, y_grid(90,:))
             i0 = floor(ri0)
             j0 = floor(rj0)
             if ( i0 > ieg .or. j0 > jeg ) then
                write (UNIT=emsg_local, FMT='("i0 = ",I6,", j0 = ",I6)') mpp_pe(), i0, j0
                call error_mesg('oda_core_mod::open_profile_dataset_argo',&
                     & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
             end if
             if ( isd_filt >= 1 .and. ied_filt <= ieg ) then
               if ( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
                  write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                       & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                  call error_mesg('oda_core_mod::open_profile_dataset',&
                       & 'i0,j0 out of bounds in argo01. '//trim(emsg_local), FATAL)
               end if
             end if
             if ( isd_filt < 1 .and. i0 > ied_filt-1 .and. i0 < isd_filt + ieg ) then
                write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                     & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                call error_mesg('oda_core_mod::open_profile_dataset',&
                     & 'i0,j0 out of bounds in argo02. '//trim(emsg_local), FATAL)
             end if
             if ( ied_filt > ieg .and. i0 > ied_filt-ieg-1 .and. ied_filt < isd_filt ) then
                write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                     & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                call error_mesg('oda_core_mod::open_profile_dataset',&
                     & 'i0,j0 out of bounds in argo03. '//trim(emsg_local), FATAL)
             end if
             Profiles(num_profiles)%i_index = ri0
             Profiles(num_profiles)%j_index = rj0
          else ! tripolar grids
             lon_out(1,1) = lon*DEG_TO_RAD
             lat_out(1,1) = lat*DEG_TO_RAD
             call horiz_interp_bilinear_new (Interp, x_grid*DEG_TO_RAD, y_grid*DEG_TO_RAD, lon_out, lat_out)
             if(Interp%wti(1,1,2) < 1.0) then
                i0 = Interp%i_lon(1,1,1)
             else
                i0 = Interp%i_lon(1,1,2)
             end if
             if ( Interp%wtj(1,1,2) < 1.0 ) then
                j0 = Interp%j_lat(1,1,1)
             else
                j0 = Interp%j_lat(1,1,2)
             end if
             if ( i0 > ieg .or. j0 > jeg ) then
                write (UNIT=emsg_local, FMT='("i0 = ",I6,", j0 = ",I6)') mpp_pe(), i0, j0
                call error_mesg('oda_core_mod::open_profile_dataset_argo',&
                     & 'For tirpolar grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
             end if
             if ( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
!!$                print*,'argo.pe,i0,j0= ',mpp_pe(), i0, j0,&
!!$                     & 'isd_filt,ied_filt,jsd_filt,jed_filt= ',isd_filt,ied_filt,jsd_filt,jed_filt
!!$                print*,'pe,lon,lat=',mpp_pe(),lon,lat,'x_grid(i0+-1)',x_grid(i0-1:i0+1,j0),&
!!$                     & 'y_grid(i0,j0+-1)=',y_grid(i0,j0-1:j0+1)
!!$                print*,'lono11,lato11=',x_grid(i0,j0),y_grid(i0,j0),'lono21,lato21=',x_grid(i0+1,j0),y_grid(i0+1,j0)
!!$                print*,'lono12,lato12=',x_grid(i0,j0+1),y_grid(i0,j0+1),'lono22,lato22=',x_grid(i0+1,j0+1),y_grid(i0+1,j0+1)
!!$                print*,'lonm11,latm11=',x_grid(isd_filt,jsd_filt),y_grid(isd_filt,jsd_filt),&
!!$                     & 'lonm21,latm21=',x_grid(ied_filt,jsd_filt),y_grid(ied_filt,jsd_filt)
!!$                print*,'lonm12,latm12=',x_grid(isd_filt,jed_filt),y_grid(isd_filt,jed_filt),&
!!$                     & 'lonm22,latm22=',x_grid(ied_filt,jed_filt),y_grid(ied_filt,jed_filt)
!!$                print*,'wti(1:2)=',Interp%wti(1,1,:),'wtj(1:2)=',Interp%wtj(1,1,:)

                out_bound_point = out_bound_point + 1
             end if
             if ( Interp%wti(1,1,2) < 1.0 ) then
                Profiles(num_profiles)%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
             else
                Profiles(num_profiles)%i_index =Interp%i_lon(1,1,2)
             end if
             if ( Interp%wtj(1,1,2) < 1.0 ) then
                Profiles(num_profiles)%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
             else
                Profiles(num_profiles)%j_index =Interp%j_lat(1,1,2)
             end if
          end if ! grids

           Profiles(num_profiles)%accepted = .true.
           if ( i0 < 1 .or. j0 < 1 ) then
              Profiles(num_profiles)%accepted = .false.
           end if
           if ( i0 < isd_filt .or. i0 >= ied_filt .or. j0 < jsd_filt .or. j0 >= jed_filt ) then
              Profiles(num_profiles)%accepted = .false.
           end if

           if ( Profiles(num_profiles)%accepted ) then ! here
              if ( i0 /= ieg .and. j0 /= jeg ) then
                 if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                      & Grd%mask(i0+1,j0,1) == 0.0 .or.&
                      & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                      & Grd%mask(i0+1,j0+1,1) == 0.0) then
                    Profiles(num_profiles)%accepted = .false.
                 end if
              else if ( i0 == ieg .and. j0 /= jeg ) then
                 if (Grd%mask(i0,j0,1) == 0.0 .or.&
                      & Grd%mask(1,j0,1) == 0.0 .or.&
                      & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                      & Grd%mask(1,j0+1,1) == 0.0) then
                    Profiles(num_profiles)%accepted = .false.
                 end if
              else if ( i0 /= ieg .and. j0 == jeg ) then
                 if (Grd%mask(i0,j0,1) == 0.0 .or.&
                      & Grd%mask(i0+1,j0,1) == 0.0) then
                    Profiles(num_profiles)%accepted = .false.
                 end if
              else
                 if ( Grd%mask(i0,j0,1) == 0.0 ) then
                    Profiles(num_profiles)%accepted = .false.
                 end if
              end if
           end if ! here

           if ( Profiles(num_profiles)%accepted .and. Profiles(num_profiles)%inst_type == MOORING+TAO) then
              if ( allocated(mask_tao) ) then
                 if ( mask_tao(i0,j0) < 1.0 ) then
                    Profiles(num_profiles)%accepted = .false.
                    write (UNIT=stdout_unit,&
                         & FMT='("Rejecting tao mooring at (lat,lon) = (",F10.5,",",F10.5,") based on user-specified mask.")')&
                         & Profiles(num_profiles)%lat,&
                         & Profiles(num_profiles)%lon
                 end if
              end if
           end if

           if ( Profiles(num_profiles)%accepted ) then
              Profiles(num_profiles)%flag(:) = .true.
              allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
              do k=1, Profiles(num_profiles)%levels
                 if (Profiles(num_profiles)%depth(k) < Grd%z(1)) then
                   Profiles(num_profiles)%k_index(k) = 1.0
                 else
                   Profiles(num_profiles)%k_index(k) = frac_index(Profiles(num_profiles)%depth(k), (/0.,Grd%z(:)/))! - 1 snz modify to v3.2 JAN3012
                 end if
                 if ( Profiles(num_profiles)%k_index(k) < 1 ) then
                    if ( Profiles(num_profiles)%depth(k) < 0 ) then
                       Profiles(num_profiles)%k_index(k) = 0
                    else if ( Profiles(num_profiles)%depth(k) > Grd%z(size(Grd%z,1)) ) then
                       Profiles(num_profiles)%k_index(k) = nk
                    end if
                 else
                    Profiles(num_profiles)%k_index(k) = Profiles(num_profiles)%k_index(k) - 1
                 end if
                 if ( Profiles(num_profiles)%k_index(k) > nk ) then 
                    call error_mesg('oda_core_mod::open_profile_dataset_argo', 'Profile k_index is greater than nk', FATAL)
                 else if ( Profiles(num_profiles)%k_index(k) < 0 ) then
                    call error_mesg('oda_core_mod::open_profile_dataset_argo', 'Profile k_index is less than 0', FATAL)
                 end if
                 k0 = floor(Profiles(num_profiles)%k_index(k))

                 if ( k0 >= 1 ) then ! snz add
                    if ( Profiles(num_profiles)%flag(k) ) then ! flag
                       if ( i0 /= ieg .and. j0 /= jeg ) then
                          if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                               & Grd%mask(i0+1,j0,k0) == 0.0 .or.&
                               & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                               & Grd%mask(i0+1,j0+1,k0) == 0.0 ) then
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       else if ( i0 == ieg .and. j0 /= jeg ) then
                          if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                               & Grd%mask(1,j0,k0) == 0.0 .or.&
                               & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                               & Grd%mask(1,j0+1,k0) == 0.0 ) then 
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       else if ( i0 /= ieg .and. j0 == jeg ) then
                          if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                               & Grd%mask(i0+1,j0,k0) == 0.0 ) then 
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       else 
                          if ( Grd%mask(i0,j0,k0) == 0.0 ) then 
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       end if

                       if ( i0 /= ieg .and. j0 /= jeg ) then
                          if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                               & Grd%mask(i0+1,j0,k0+1) == 0.0 .or.&
                               & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                               & Grd%mask(i0+1,j0+1,k0+1) == 0.0 ) then
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       else if ( i0 == ieg .and. j0 /= jeg ) then
                          if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                               & Grd%mask(1,j0,k0+1) == 0.0 .or.&
                               & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                               & Grd%mask(1,j0+1,k0+1) == 0.0) then
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       else if ( i0 /= ieg .and. j0 == jeg ) then
                          if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                               & Grd%mask(i0+1,j0,k0+1) == 0.0 ) then
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       else
                          if ( Grd%mask(i0,j0,k0+1) == 0.0 ) then 
                             Profiles(num_profiles)%flag(k) = .false.
                          end if
                       end if

                       if ( abs(Profiles(num_profiles)%data(k)) > 1.e4 &
                            & .or. abs(Profiles(num_profiles)%depth(k)) > 1.e4 ) then
                          Profiles(num_profiles)%flag(k) = .false.
                       end if
                    end if ! flag
                 end if ! snz add
              end do
           end if  ! accepted 

           end if ! 05 Nov 2012

        else ! localize
           i = i+1
        end if ! localize

        if ( i .gt. nstation ) cont = .false.
     end do

     a = nprof_in_filt_domain
     call mpp_sum(a)
     call mpp_sum(out_bound_point)

    if ( no_prf /= num_profiles ) then
       write(UNIT=stdout_unit, FMT='("PE: ",I6," no_prf = ",I8,", num_profiles = ",I8)') mpp_pe(), no_prf, num_profiles
    end if
    if ( var_id == TEMP_ID ) then
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," argo temp prfs within global domain")') no_temp
       write(UNIT=stdout_unit, FMT='("A total out of bound points",I8," argo temp within global domain")') out_bound_point
    else if ( var_id == SALT_ID ) then
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," argo salt prfs within global domain")') no_salt
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," argo prfs within global domain")') no_prf
       write(UNIT=stdout_unit, FMT='("A grand total of ",I8," argo prfs within current PEs computer domain")') a
       write(UNIT=stdout_unit, FMT='("A total out of bound points",I8," argo salt within global domain")') out_bound_point
    end if

    call mpp_sync_self()
    call mpp_close(unit)
    deallocate(axes)
    deallocate(fields)

!    call mpp_print_memuse_stats('open_profile_dataset_argo End')

  end subroutine open_profile_dataset_argo

  ! get profiles and sfc
  ! obs relevant to current analysis interval
  subroutine get_obs(model_time, Prof, nprof)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), dimension(:), intent(inout) :: Prof
    integer, intent(inout) :: nprof

    integer :: i, k, kk, k_interval
    integer :: yr, mon, day, hr, min, sec
    integer :: stdout_unit

    type(time_type) :: tdiff

    nprof = 0
    stdout_unit = stdout()

    write (UNIT=stdout_unit, FMT='("Gathering profiles for current analysis time")')
    call get_date(model_time, yr, mon, day, hr, min, sec)
    write (UNIT=stdout_unit, FMT='("Current YYYY/MM/DD = ",I4,"/",I2,"/",I2)') yr, mon, day

    do i=1, no_prf
       if ( Profiles(i)%time <= model_time ) then
          tdiff = model_time - Profiles(i)%time
       else
          tdiff = Profiles(i)%time - model_time
       end if

       ! no tdiff criteria for monthly mean data like
       ! but tdiff criteria has to be set for daily data
       if ( tdiff <= time_window(Profiles(i)%inst_type) .and. Profiles(i)%accepted ) then
          ! for single profile test

          nprof = nprof + 1
          if ( nprof > size(Prof,1) ) then
             call error_mesg('oda_core_mod::get_obs',&
                  & 'Passed in array "Prof" is smaller than number of profiles, increase size of Prof before call.',&
                  & FATAL)
          end if
          call copy_obs(Profiles(i:i), Prof(nprof:nprof))

          Prof(nprof)%tdiff = tdiff
          
          ! snz add the following few lines for increasing deep water data
          if ( Prof(nprof)%levels > max_prflvs ) then
             k_interval = (Prof(nprof)%levels-max_prflvs+50)/50 + 1
             kk = max_prflvs - 50
             do k=max_prflvs-50+1, Prof(nprof)%levels, k_interval
                kk = kk + 1
                Prof(nprof)%depth(kk) = Prof(nprof)%depth(k)
                Prof(nprof)%k_index(kk) = Prof(nprof)%k_index(k)
                Prof(nprof)%data(kk) = Prof(nprof)%data(k)
                Prof(nprof)%flag(kk) = Prof(nprof)%flag(k)
             end do
             Prof(nprof)%levels = kk
          end if
          ! snz end the adding lines
       end if
    end do

    write (UNIT=stdout_unit,&
         & FMT='("A total of ",I8," profiles are being used for the current analysis step.")') nprof

    return
  end subroutine get_obs

  subroutine oda_core_init(Domain, Grid, time_s, time_e, filt_domain, localize)
    type(domain2d), intent(inout) :: Domain
    type(grid_type), target, intent(in) :: Grid
    logical, intent(in), optional :: localize
    type(time_type), intent(in) :: time_s, time_e
    type(domain2d), intent(in) :: filt_domain

    integer :: ioun, ierr, io_status
    integer :: stdlog_unit

    stdlog_unit = stdlog()

    ! Read in the namelist file
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, oda_core_nml, iostat=io_status)
#else
    ioun = open_namelist_file()
    read(ioun, NML=oda_core_nml, IOSTAT=io_status)
    ierr = check_nml_error(io_status, 'oda_core_nml')
    call close_file(ioun)
#endif

    write(stdlog_unit, NML=oda_core_nml)

    Grd => Grid
    
    call mpp_get_compute_domain(Domain, isc, iec, jsc, jec)
    call mpp_get_data_domain(Domain, isd, ied, jsd, jed)
    call mpp_get_global_domain(Domain, isg, ieg, jsg, jeg)
    call mpp_get_data_domain(filt_domain, isd_filt, ied_filt, jsd_filt, jed_filt)

    jsd_flt0 = jsd_filt
    jed_flt0 = jed_filt
    if (jsd_filt < 1) jsd_flt0 = 1
    if (jed_filt > jeg) jed_flt0 = jeg

    nk = size(Grid%z)

    call init_observations(time_s, time_e, filt_domain, localize)
  end subroutine oda_core_init
    
  subroutine copy_obs(obs_in, obs_out)
    type(ocean_profile_type), dimension(:), intent(in) :: obs_in
    type(ocean_profile_type), dimension(:), intent(inout) :: obs_out

    integer :: n

    if ( size(obs_in) .ne. size(obs_out) ) then
       call error_mesg('oda_core_mod::copy_obs', 'Size of in and out obs variables are not equal.', FATAL)
    end if

    do n=1, size(obs_in)
       Obs_out(n)%variable = Obs_in(n)%variable
       Obs_out(n)%inst_type = Obs_in(n)%inst_type
       Obs_out(n)%levels = Obs_in(n)%levels
       Obs_out(n)%lon = Obs_in(n)%lon
       Obs_out(n)%lat = Obs_in(n)%lat
       Obs_out(n)%accepted = Obs_in(n)%accepted
       if ( associated(Obs_out(n)%depth) ) then
          deallocate(Obs_out(n)%depth)
          nullify(Obs_out(n)%depth)
       end if
       allocate(Obs_out(n)%depth(Obs_in(n)%levels))
       Obs_out(n)%depth(:) = Obs_in(n)%depth(:)
       if ( associated(Obs_out(n)%data) ) then
          deallocate(Obs_out(n)%data)
          nullify(Obs_out(n)%data)
       end if
       allocate(Obs_out(n)%data(Obs_in(n)%levels))
       Obs_out(n)%data(:) = Obs_in(n)%data(:)
       if ( associated(Obs_out(n)%flag) ) then
          deallocate(Obs_out(n)%flag)
          nullify(Obs_out(n)%flag)
       end if
       allocate(Obs_out(n)%flag(Obs_in(n)%levels))
       Obs_out(n)%flag(:) = Obs_in(n)%flag(:)          
       Obs_out(n)%time = Obs_in(n)%time
       Obs_out(n)%yyyy = Obs_in(n)%yyyy
       Obs_out(n)%mmdd = Obs_in(n)%mmdd
       Obs_out(n)%i_index = Obs_in(n)%i_index
       Obs_out(n)%j_index = Obs_in(n)%j_index
       if ( associated(Obs_out(n)%k_index) ) then
          deallocate(Obs_out(n)%k_index)
          nullify(Obs_out(n)%k_index)
       end if
       allocate(Obs_out(n)%k_index(Obs_in(n)%levels))          
       Obs_out(n)%k_index = Obs_in(n)%k_index

!       if ( associated(Obs_out(n)%ms) ) then
!          deallocate(Obs_out(n)%ms)
!          nullify(Obs_out(n)%ms)
!       end if
!       allocate(Obs_out(n)%ms(Obs_in(n)%levels))          
!       Obs_out(n)%ms = Obs_in(n)%ms
!       if ( associated(Obs_out(n)%ms_inv) ) then
!          deallocate(Obs_out(n)%ms_inv)
!          nullify(Obs_out(n)%ms_inv)
!       end if
!       allocate(Obs_out(n)%ms_inv(Obs_in(n)%levels))          
!       Obs_out(n)%ms_inv = 1./Obs_in(n)%ms    

       Obs_out(n)%tdiff = Obs_in(n)%tdiff
       if ( associated(Obs_out(n)%Forward_model%wgt) ) then
          deallocate(Obs_out(n)%Forward_model%wgt)
          nullify(Obs_out(n)%Forward_model%wgt)
       end if
    end do
  end subroutine copy_obs

  subroutine open_profile_dataset_sst(filename, obs_variable, localize)  
    character(len=*), intent(in) :: filename
    integer, intent(in) :: obs_variable
    logical, intent(in), optional :: localize

    integer, parameter :: MAX_LEVELS = 1

    real :: lon, lat, rms_err
    real :: ri0, rj0
    real, dimension(MAX_LEVELS) :: depth, data

    integer :: unit, ndim, nvar, natt, ntime
    integer :: var_id, inst_type
    integer :: num_levs, k, kk, i, j, i0, j0
    integer :: stdout_unit

    logical :: data_is_local, localize_data
    logical, dimension(MAX_LEVELS) :: flag

    character(len=32) :: axisname, anal_fldname
    character(len=128) :: emsg_local

    type(axistype), dimension(:), allocatable, target :: axes
    type(axistype), pointer :: lon_axis, lat_axis

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    stdout_unit = stdout()

    anal_fldname = 'temp'
    var_id = obs_variable

    call mpp_open(unit, trim(filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)

    do i=1, ndim
       call mpp_get_atts(axes(i),name=axisname)
       select case (trim(axisname))
!!$       case ('GRIDLON_T') ! cm2 grids
!!$       case ('gridlon_t') ! cm2 grids
       case ('XT_OCEAN') ! cm2.5 grids
!!$       case ('lon') ! after 2008
          lon_axis => axes(i)
!!$       case ('GRIDLAT_T') ! cm2 grids
!!$       case ('gridlat_t') ! cm2 grids
       case ('YT_OCEAN') ! cm2 grids
!!$       case ('lat') ! for after 2008
          lat_axis => axes(i)
       end select
    end do

    call mpp_get_atts(lon_axis,len=nlon)
    call mpp_get_atts(lat_axis,len=nlat)
    if ( nlon /= 1440 .or. nlat /= 1070 ) then ! after 2008
       write (UNIT=emsg_local, FMT='("sst obs dim is not same as in file. nlon = ",I5,", nlat = ",I5)') nlon, nlat
       call error_mesg('oda_core_mod::open_profile_dataset_sst', trim(emsg_local), FATAL)
    end if

    ! idealized
    do j=1, nlat
       do i=1, nlon
          lon = x_grid(i,j)
          lat = y_grid(i,j)
          rms_err = 0.5
          inst_type = 20
          data_is_local = .true.

          if ( lon .lt. 0.0 ) lon = lon + 360.0
          if ( lon .gt. 360.0 ) lon = lon - 360.0
          if ( lon .lt. 80.0 ) lon = lon + 360.0

          if ( lat < sst_obs_start_lat .or. lat > sst_obs_end_lat ) data_is_local = .false. ! at the final test

          if ( Grd%mask(i,j,1) == 0 ) data_is_local = .false.

          if ( abs(lat) < 40.0 ) then
             if ( i/4*4 /= i .or. j/4*4 /= j ) data_is_local = .false.
          else if ( abs(lat) < 60.0 ) then
             if ( i/8*8 /= i .or. j/6*6 /= j ) data_is_local = .false.
          else
             if ( i/16*16 /= i .or. j/8*8 /= j ) data_is_local = .false.
          end if

          if ( data_is_local .and. (.NOT.localize_data) ) then
             if ( lat < 60.0 ) then ! regular grids
                ri0 = frac_index(lon, x_grid(:,nlat/2))
                rj0 = frac_index(lat, y_grid(nlon/4,:))
                i0 = floor(ri0)
                j0 = floor(rj0)
             else ! tripolar grids
                lon_out(1,1) = lon*DEG_TO_RAD
                lat_out(1,1) = lat*DEG_TO_RAD
                call horiz_interp_bilinear_new (Interp, x_grid*DEG_TO_RAD, y_grid*DEG_TO_RAD, lon_out, lat_out)
                if ( Interp%wti(1,1,2) < 1.0 ) then
                   i0 = Interp%i_lon(1,1,1)
                else
                   i0 = Interp%i_lon(1,1,2)
                end if
                if ( Interp%wtj(1,1,2) < 1.0 ) then
                   j0 = Interp%j_lat(1,1,1)
                else
                   j0 = Interp%j_lat(1,1,2)
                end if

                if ( i0 > ieg .or. j0 > jeg ) then
                   write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') i0, j0
                   call error_mesg('oda_core_mod::open_profile_dataset_sst',&
                        & 'For tripolar grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
                end if
             end if

             if ( i0 /= ieg .and. j0 /= jeg ) then ! exclude SSTs at ieg and jeg
                if ( Grd%mask(i0,j0,1) /= 0.0 .and. Grd%mask(i0+1,j0,1) /= 0.0 .and.&
                     & Grd%mask(i0,j0+1,1) /= 0.0 .and. Grd%mask(i0+1,j0+1,1) /= 0.0 ) then
                   no_sst = no_sst+1
                   num_profiles=num_profiles+1

                   if ( num_profiles > max_profiles ) then
                      call error_mesg('oda_core_mod::open_profile_dataset_sst',&
                           & 'Maximum number of profiles exceeded, increase max_profiles in oda_core_nml', FATAL)
                   end if

                   num_levs = 0
                   flag = .false.
                   do k=1, 1
                      flag(k) = .true.
                      data(k) = 0.0
                      depth(k) = 0.0
                      num_levs = num_levs + 1
                   end do
                   if ( num_levs == 0 ) cycle
                   allocate(profiles(num_profiles)%depth(num_levs))
                   allocate(profiles(num_profiles)%data(num_levs))
                   allocate(profiles(num_profiles)%flag(num_levs))
                   allocate(profiles(num_profiles)%ms(num_levs))
                   allocate(profiles(num_profiles)%ms_inv(num_levs))
                   profiles(num_profiles)%variable = var_id
                   profiles(num_profiles)%inst_type = inst_type
                   profiles(num_profiles)%levels = num_levs
                   profiles(num_profiles)%lat = lat
                   profiles(num_profiles)%lon = lon

                   kk = 1
                   do k=1, 1
                      if ( flag(k) ) then
                         profiles(num_profiles)%depth(kk) = depth(k)
                         profiles(num_profiles)%data(kk) = data(k)
                         profiles(num_profiles)%ms(kk) = 1.0
                         profiles(num_profiles)%ms_inv(kk) = 1.0
                         kk=kk+1
                      end if
                   end do

                   ! calculate interpolation coefficients (make sure to account for grid offsets here!)
                   if ( lat < 60.0 ) then ! for regular grids
                      Profiles(num_profiles)%i_index = ri0
                      Profiles(num_profiles)%j_index = rj0
                   else ! for tripolar grids
                      lon_out(1,1) = lon*DEG_TO_RAD
                      lat_out(1,1) = lat*DEG_TO_RAD
                      call horiz_interp_bilinear_new (Interp, x_grid*DEG_TO_RAD, y_grid*DEG_TO_RAD, lon_out, lat_out)
                      if ( Interp%wti(1,1,2) < 1.0 ) then
                         Profiles(num_profiles)%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
                      else
                         Profiles(num_profiles)%i_index =Interp%i_lon(1,1,2)
                      end if
                      if ( Interp%wtj(1,1,2) < 1.0 ) then
                         Profiles(num_profiles)%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
                      else
                         Profiles(num_profiles)%j_index =Interp%j_lat(1,1,2)
                      end if
                   end if

                   Profiles(num_profiles)%accepted = .true.
                   if ( i0 < 1 .or. j0 < 1 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                   if ( Profiles(num_profiles)%accepted ) then
                      if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                           & Grd%mask(i0+1,j0,1) == 0.0 .or.&
                           & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                           & Grd%mask(i0+1,j0+1,1) == 0.0 ) then
                         Profiles(num_profiles)%accepted = .false.
                      end if
                   end if
                   if ( Profiles(num_profiles)%accepted ) then
                      Profiles(num_profiles)%flag(:) = .true.
                      allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
                      do k=1, Profiles(num_profiles)%levels
                         Profiles(num_profiles)%k_index(k) = frac_index(depth(k), (/0.,Grd%z(:)/)) - 1
                         !::sdu:: Do we need the same out-of-range check here?
                      end do
                   end if
                end if ! exclude SSTs at ieg and jeg
             end if
          end if
       end do
    end do

    call mpp_close(unit)
    deallocate(axes)
    write (UNIT=stdout_unit, FMT='("A grand total of ",I8," sst points within global domain")') no_sst
    write (UNIT=stdout_unit, FMT='("A final total @sst of ",I8," prfs within global domain")') num_profiles
  end subroutine open_profile_dataset_sst

  subroutine open_profile_dataset_woa05t(filename, obs_variable, localize)  
    character(len=*), intent(in) :: filename
    integer, intent(in) :: obs_variable
    logical, intent(in), optional :: localize

    integer, parameter :: MAX_LEVELS = 24

    real :: lon, lat, rms_err
    real :: ri0, rj0
    real, dimension(MAX_LEVELS) :: depth, data

    integer :: unit, ndim, nvar, natt, ntime
    integer :: var_id, inst_type
    integer :: num_levs, k, kk, i, j, i0, j0, k0
    integer :: stdout_unit, istat
    integer :: out_bound_point

    logical :: data_is_local, localize_data
    logical :: prof_in_filt_domain
    logical, dimension(MAX_LEVELS) :: flag

    character(len=32) :: axisname, anal_fldname
    character(len=128) :: emsg_local

    type(axistype), dimension(:), allocatable, target :: axes
    type(axistype), pointer :: lon_axis, lat_axis, z_axis, t_axis

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    stdout_unit = stdout()
    anal_fldname = 'temp'
    var_id = obs_variable

!    call mpp_print_memuse_stats('open_profile_dataset_woa05t Start')

    call mpp_open(unit, trim(filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)

    write (UNIT=stdout_unit, FMT='("Opened profile woa05t dataset: ",A)') trim(filename)

    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(axes(ndim))
    call mpp_get_axes(unit,axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case (trim(axisname))
       case ('lon')
          lon_axis => axes(i)
       case ('lat')
          lat_axis => axes(i)
       case ('depth')
          z_axis => axes(i)
       end select
    end do

    call mpp_get_atts(lon_axis, len=nlon_woa)
    call mpp_get_atts(lat_axis, len=nlat_woa)
    call mpp_get_atts(z_axis, len=nlev_woa)

    if ( nlon_woa /= 360 .or. nlat_woa /= 180 ) then
       write (UNIT=emsg_local, FMT='("woa05 obs dim is not same as in file. nlon_woa = ",I8,", nlat_woa = ",I8)') nlon_woa, nlat_woa
       call error_mesg('oda_core_mod::open_profile_dataset_woa05t', trim(emsg_local), FATAL)
    end if

    allocate(woa05_lon(nlon_woa), woa05_lat(nlat_woa), woa05_z(nlev_woa))

    call mpp_get_axis_data(lon_axis, woa05_lon)
    call mpp_get_axis_data(lat_axis, woa05_lat)
    call mpp_get_axis_data(z_axis, woa05_z)

    out_bound_point = 0
    ! idealized
    do j=1, nlat_woa
       do i=1, nlon_woa
          lon = woa05_lon(i)
          lat = woa05_lat(j)
          rms_err = 5
          inst_type = 20
          data_is_local = .true.
          prof_in_filt_domain = .false.

          if ( lon .lt. 0.0 ) lon = lon + 360.0
          if ( lon .gt. 360.0 ) lon = lon - 360.0
          if ( lon .lt. 80.0 ) lon = lon + 360.0

          if ( lat < -80.0 .or. lat > 80.0 ) data_is_local = .false.
          if ( abs(lat) < 20.0 .and. (mod(i,2) /= 0 .or. mod(j,2) /= 0) ) data_is_local = .false.
          if ( abs(lat) >= 20.0 .and. (mod(i,4) /= 0 .or. mod(j,4) /= 0) ) data_is_local = .false.
          if ( abs(lat) >= 60.0 .and. (mod(i,6) /= 0 .or. mod(j,6) /= 0) ) data_is_local = .false.

          if (isd_filt < 1 .and. ied_filt > ieg) then
             ! filter domain is a full x band
             if (lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  lat >= y_grid(1,jsd_flt0) .and. lat <= y_grid(ieg-1,jsd_flt0)) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt <= ieg) then
             ! Interior filter domain
            if (lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
               & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)) then
                prof_in_filt_domain = .true.
            end if
          else if (isd_filt < 1 .and. ied_filt <= ieg) then
             ! lhs filter domain
             isd_flt0 = isd_filt + ieg
             if ((lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
                  & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)).or.&
                  & (lon >= x_grid(isd_flt0,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_flt0,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1))) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt > ieg) then
             ! rhs filter domain
             ied_flt0 = ied_filt - ieg
             if ( lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1) ) then
                prof_in_filt_domain = .true.
             end if
             if (ied_flt0-1 > 1) then
                if ( lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_flt0-1,jsd_flt0) .and.&
                     & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_flt0-1,jed_flt0-1) ) then
                   prof_in_filt_domain = .true.
                end if
             end if
          end if

          if ( data_is_local .and. (.NOT.localize_data) ) then ! global index

            no_woa05 = no_woa05 + 1
            num_profiles = num_profiles + 1
             if ( num_profiles > max_profiles ) then
                call error_mesg('oda_core_mod::open_profile_dataset_woa05t',&
                     & 'Maximum number of profiles exceeded, increase max_profiles in oda_core_nml', FATAL)
             end if

            num_levs = 0
            flag = .false.
            do k=1, nlev
              flag(k) = .true.
              data(k) = 0.0
              depth(k) = woa05_z(k)
              num_levs = num_levs+1
            end do
            if ( num_levs == 0 ) cycle

            if ( prof_in_filt_domain ) then ! localize

             allocate(profiles(num_profiles)%depth(num_levs))
             allocate(profiles(num_profiles)%data(num_levs))
             allocate(profiles(num_profiles)%flag(num_levs))
!             allocate(profiles(num_profiles)%ms(num_levs))
!             allocate(profiles(num_profiles)%ms_inv(num_levs))
             profiles(num_profiles)%variable = var_id
             profiles(num_profiles)%inst_type = inst_type
             profiles(num_profiles)%levels = num_levs
             profiles(num_profiles)%lat = lat
             profiles(num_profiles)%lon = lon
             
             kk = 1
             do k=1, nlev
                if ( flag(k) ) then
                   profiles(num_profiles)%depth(kk) = depth(k)
                   profiles(num_profiles)%data(kk) = data(k)
!                   profiles(num_profiles)%ms(kk) = 1.0
!                   profiles(num_profiles)%ms_inv(kk) = 1.0
                   kk = kk + 1
                end if
             end do

             ! calculate interpolation coefficients (make sure to account for grid offsets here!)
             if ( lat < 65.0 ) then ! regular grids
                ri0 = frac_index(lon, x_grid(:,1))
                rj0 = frac_index(lat, y_grid(90,:))
                i0 = floor(ri0)
                j0 = floor(rj0)
                if ( i0 > ieg .or. j0 > jeg ) then
                   write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') i0, j0
                   call error_mesg('oda_core_mod::open_profile_dataset_woa05t',&
                        & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
                end if
                if ( isd_filt >= 1 .and. ied_filt <= ieg ) then
                  if ( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
                     write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                          & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                     call error_mesg('oda_core_mod::open_profile_dataset',&
                          & 'i0,j0 out of bounds in woat01. '//trim(emsg_local), FATAL)
                  end if
                end if
                if ( isd_filt < 1 .and. i0 > ied_filt-1 .and. i0 < isd_filt + ieg ) then
                   write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'i0,j0 out of bounds in woat02. '//trim(emsg_local), FATAL)
                end if
                if ( ied_filt > ieg .and. i0 > ied_filt-ieg-1 .and. ied_filt < isd_filt ) then
                   write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'i0,j0 out of bounds in woat03. '//trim(emsg_local), FATAL)
                end if
                Profiles(num_profiles)%i_index = ri0
                Profiles(num_profiles)%j_index = rj0 
             else ! tripolar grids
                lon_out(1,1) = lon*DEG_TO_RAD
                lat_out(1,1) = lat*DEG_TO_RAD
                call horiz_interp_bilinear_new (Interp, x_grid*DEG_TO_RAD, y_grid*DEG_TO_RAD, lon_out, lat_out)
                if ( Interp%wti(1,1,2) < 1.0 ) then
                   i0 = Interp%i_lon(1,1,1)
                else
                   i0 = Interp%i_lon(1,1,2)
                end if
                if ( Interp%wtj(1,1,2) < 1.0 ) then
                   j0 = Interp%j_lat(1,1,1)
                else
                   j0 = Interp%j_lat(1,1,2)
                end if
                if ( i0 > ieg .or. j0 > jeg ) then
                   write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') i0, j0
                   call error_mesg('oda_core_mod::open_profile_dataset_woa05t',&
                        & 'For tripolar grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
                end if
                if ( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
                   write (UNIT=stdout_unit, FMT='("woat.pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(),i0, j0,isd_filt,ied_filt,jsd_filt,jed_filt
                   out_bound_point = out_bound_point + 1
                end if
                if ( Interp%wti(1,1,2) < 1.0 ) then
                   Profiles(num_profiles)%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
                else
                   Profiles(num_profiles)%i_index =Interp%i_lon(1,1,2)
                end if
                if ( Interp%wtj(1,1,2) < 1.0 ) then
                   Profiles(num_profiles)%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
                else
                   Profiles(num_profiles)%j_index =Interp%j_lat(1,1,2)
                end if
             end if ! grids
           
             Profiles(num_profiles)%accepted = .true.
             if ( i0 < 1 .or. j0 < 1 ) then
                Profiles(num_profiles)%accepted = .false.
             end if
             if ( i0 < isd_filt .or. i0 >= ied_filt .or. j0 < jsd_filt .or. j0 >= jed_filt ) then
                Profiles(num_profiles)%accepted = .false.
             end if

             if ( Profiles(num_profiles)%accepted ) then ! here
                if ( i0 /= ieg .and. j0 /= jeg ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else if ( i0 == ieg .and. j0 /= jeg ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else if ( i0 /= ieg .and. j0 == jeg ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else
                   if ( Grd%mask(i0,j0,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                end if
             end if ! here
             
             if ( Profiles(num_profiles)%accepted ) then ! accepted
                Profiles(num_profiles)%flag(:) = .true.
                allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
                do k=1, Profiles(num_profiles)%levels
                   if (depth(k) < Grd%z(1)) then
                     Profiles(num_profiles)%k_index(k) = 0.0
                   else
                     Profiles(num_profiles)%k_index(k) = frac_index(depth(k), (/0.,Grd%z(:)/)) - 1.0 ! snz modify to v3.2 JAN3012
                   end if
                   if ( Profiles(num_profiles)%k_index(k) > nk ) then
                      call error_mesg('oda_core_mod::open_profile_dataset_woa05t',&
                           & 'Profile k_index is greater than nk', FATAL)
                   end if
                   k0 = floor(Profiles(num_profiles)%k_index(k))
                   
                   if ( k0 >= 1 ) then ! snz add
                      if ( Profiles(num_profiles)%flag(k) ) then ! flag
                         if ( i0 /= ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0+1,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 == ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(1,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                                 & Grd%mask(1,j0+1,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 /= ieg .and. j0 == jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else
                            if ( Grd%mask(i0,j0,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         end if
                         
                         if ( i0 /= ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0+1,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 == ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(1,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                                 & Grd%mask(1,j0+1,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 /= ieg .and. j0 == jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0+1) == 0.0) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         end if
                         
                         if ( abs(Profiles(num_profiles)%data(k)) > 1.e4 &
                              & .or. abs(Profiles(num_profiles)%depth(k)) > 1.e4 ) then
                            Profiles(num_profiles)%flag(k) = .false.
                         end if
                      end if ! flag
                   end if ! snz add
                end do
             end if ! accepted
           end if ! localize
         end if ! global index
       end do
    end do

    call mpp_close(unit)
    deallocate(axes)
    write (UNIT=stdout_unit, FMT='("A grand total of ",I8," woa05t points within global domain")') no_woa05
    write (UNIT=stdout_unit, FMT='("A final total @woa05t of ",I8," prfs within global domain")') num_profiles

!    call mpp_print_memuse_stats('open_profile_dataset_woa05t End')

  end subroutine open_profile_dataset_woa05t

  subroutine open_profile_dataset_woa05s(filename, obs_variable, localize)  
    character(len=*), intent(in) :: filename
    integer, intent(in) :: obs_variable
    logical, intent(in), optional :: localize

    integer, parameter :: MAX_LEVELS = 24

    real :: lon, lat, rms_err
    real :: ri0, rj0
    real, dimension(MAX_LEVELS) :: depth, data

    integer :: unit, ndim, nvar, natt, ntime
    integer :: var_id, inst_type
    integer :: num_levs, k, kk, i, j, i0, j0, k0
    integer :: stdout_unit
    integer :: out_bound_point

    logical :: data_is_local, localize_data
    logical :: prof_in_filt_domain
    logical, dimension(MAX_LEVELS) :: flag

    character(len=32) :: axisname, anal_fldname
    character(len=128) :: emsg_local

    type(axistype), dimension(:), allocatable, target :: axes
    type(axistype), pointer :: lon_axis, lat_axis, z_axis

    if ( present(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    stdout_unit = stdout()

    anal_fldname = 'salt'
    var_id = obs_variable

!    call mpp_print_memuse_stats('open_profile_dataset_woa05s Start')

    call mpp_open(unit, trim(filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)

    write (UNIT=stdout_unit, FMT='("Opened profile woa05s dataset: ",A)') trim(filename)

    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    allocate(axes(ndim))
    call mpp_get_axes(unit, axes)
    do i=1, ndim
       call mpp_get_atts(axes(i), name=axisname)
       select case ( trim(axisname) )
       case ('lon')
          lon_axis => axes(i)
       case ('lat')
          lat_axis => axes(i)
       case ('depth')
          z_axis => axes(i)
       end select
    end do
    call mpp_get_atts(lon_axis, len=nlon_woa)
    call mpp_get_atts(lat_axis, len=nlat_woa)
    call mpp_get_atts(z_axis, len=nlev_woa)
    if ( nlon_woa /= 360 .or. nlat_woa /= 180 ) then
       write (UNIT=emsg_local, FMT='("woa05 obs dim is not same as in file nlon_woa = ",I8,", nlat_woa = ",I8)') nlon_woa, nlat_woa
       call error_mesg('oda_core_mod::open_profile_dataset_woa05s', trim(emsg_local), FATAL)
    end if

    out_bound_point = 0
    ! idealized
    do j=1, nlat_woa
       do i=1, nlon_woa
          lon = woa05_lon(i)
          lat = woa05_lat(j)
          rms_err = 5
          inst_type = 20
          data_is_local = .true.
          prof_in_filt_domain = .false.

          if ( lon .lt. 0.0 ) lon = lon + 360.0
          if ( lon .gt. 360.0 ) lon = lon - 360.0
          if ( lon .lt. 80.0 ) lon = lon + 360.0

          if ( lat < -80.0 .or. lat > 80.0 ) data_is_local = .false.
          if ( abs(lat) < 20.0 .and. (mod(i,2) /= 0 .or. mod(j,2) /= 0) ) data_is_local = .false.
          if ( abs(lat) >= 20.0 .and. (mod(i,4) /= 0 .or. mod(j,4) /= 0) ) data_is_local = .false.
          if ( abs(lat) >= 60.0 .and. (mod(i,6) /= 0 .or. mod(j,6) /= 0) ) data_is_local = .false.

          if (isd_filt < 1 .and. ied_filt > ieg) then
             ! filter domain is a full x band
             if (lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  lat >= y_grid(1,jsd_flt0) .and. lat <= y_grid(ieg-1,jsd_flt0)) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt <= ieg) then
             ! Interior filter domain
            if (lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
               & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)) then
                prof_in_filt_domain = .true.
            end if
          else if (isd_filt < 1 .and. ied_filt <= ieg) then
             ! lhs filter domain
             isd_flt0 = isd_filt + ieg
             if ((lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_filt-1,jsd_flt0) .and.&
                  & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_filt-1,jed_flt0-1)).or.&
                  & (lon >= x_grid(isd_flt0,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_flt0,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1))) then
                prof_in_filt_domain = .true.
             end if
          else if (isd_filt >= 1 .and. ied_filt > ieg) then
             ! rhs filter domain
             ied_flt0 = ied_filt - ieg
             if ( lon >= x_grid(isd_filt,jsd_flt0) .and. lon <= x_grid(ieg-1,jsd_flt0) .and.&
                  & lat >= y_grid(isd_filt,jsd_flt0) .and. lat <=  y_grid(ieg-1,jed_flt0-1) ) then
                prof_in_filt_domain = .true.
             end if
             if (ied_flt0-1 > 1) then
                if ( lon >= x_grid(1,jsd_flt0) .and. lon <= x_grid(ied_flt0-1,jsd_flt0) .and.&
                     & lat >= y_grid(1,jsd_flt0) .and. lat <=  y_grid(ied_flt0-1,jed_flt0-1) ) then
                   prof_in_filt_domain = .true.
                end if
             end if
          end if

          if ( data_is_local .and. (.NOT.localize_data) ) then ! global index

            no_woa05 = no_woa05 + 1
            num_profiles=num_profiles + 1
             if ( num_profiles > max_profiles ) then
                call error_mesg('oda_core_mod::open_profile_dataset_woa05s',&
                     & 'Maximum number of profiles exceeded, increase max_profiles in oda_core_nml.', FATAL)
             end if

            num_levs = 0
            flag = .false.
            do k=1, nlev_woa
              flag(k) = .true.
              data(k) = 0.0
              depth(k) = woa05_z(k)
              num_levs = num_levs + 1
            end do
            if ( num_levs == 0 ) cycle

           if ( prof_in_filt_domain ) then ! localize

             allocate(profiles(num_profiles)%depth(num_levs))
             allocate(profiles(num_profiles)%data(num_levs))
             allocate(profiles(num_profiles)%flag(num_levs))
!             allocate(profiles(num_profiles)%ms(num_levs))
!             allocate(profiles(num_profiles)%ms_inv(num_levs))
             profiles(num_profiles)%variable = var_id
             profiles(num_profiles)%inst_type = inst_type
             profiles(num_profiles)%levels = num_levs
             profiles(num_profiles)%lat = lat
             profiles(num_profiles)%lon = lon

             kk = 1
             do k=1, nlev_woa
                if ( flag(k) ) then
                   profiles(num_profiles)%depth(kk) = depth(k)
                   profiles(num_profiles)%data(kk) = data(k)
!                   profiles(num_profiles)%ms(kk) = 1.0
!                   profiles(num_profiles)%ms_inv(kk) = 1.0
                   kk = kk + 1
                end if
             end do

             ! calculate interpolation coefficients (make sure to account for grid offsets here!)
             if ( lat < 65.0 ) then ! regular grids
                ri0 = frac_index(lon, x_grid(:,1))
                rj0 = frac_index(lat, y_grid(90,:))
                i0 = floor(ri0)
                j0 = floor(rj0)
                if ( i0 > ieg .or. j0 > jeg ) then
                   write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') i0, j0
                   call error_mesg('oda_core_mod::open_profile_dataset_woa05s',&
                        & 'For regular grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
                end if
                if ( isd_filt >= 1 .and. ied_filt <= ieg ) then
                  if ( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
                     write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                          & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                     call error_mesg('oda_core_mod::open_profile_dataset',&
                          & 'i0,j0 out of bounds in woas01. '//trim(emsg_local), FATAL)
                  end if
                end if
                if ( isd_filt < 1 .and. i0 > ied_filt-1 .and. i0 < isd_filt + ieg ) then
                   write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'i0,j0 out of bounds in woas02. '//trim(emsg_local), FATAL)
                end if
                if ( ied_filt > ieg .and. i0 > ied_filt-ieg-1 .and. ied_filt < isd_filt ) then
                   write (UNIT=emsg_local, FMT='("pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(), i0, j0, isd_filt,ied_filt,jsd_filt,jed_filt
                   call error_mesg('oda_core_mod::open_profile_dataset',&
                        & 'i0,j0 out of bounds in woas03. '//trim(emsg_local), FATAL)
                end if
                Profiles(num_profiles)%i_index = ri0
                Profiles(num_profiles)%j_index = rj0 
             else ! tripolar grids
                lon_out(1,1) = lon*DEG_TO_RAD
                lat_out(1,1) = lat*DEG_TO_RAD
                call horiz_interp_bilinear_new (Interp, x_grid*DEG_TO_RAD, y_grid*DEG_TO_RAD, lon_out, lat_out)
                if ( Interp%wti(1,1,2) < 1.0 ) then
                   i0 = Interp%i_lon(1,1,1)
                else
                   i0 = Interp%i_lon(1,1,2)
                end if
                if ( Interp%wtj(1,1,2) < 1.0 ) then
                   j0 = Interp%j_lat(1,1,1)
                else
                   j0 = Interp%j_lat(1,1,2)
                end if
                if ( i0 > ieg .or. j0 > jeg ) then
                   write (UNIT=emsg_local, FMT='("i0 = ",I8,", j0 = ",I8)') i0, j0
                   call error_mesg('oda_core_mod::open_profile_dataset_woa05s',&
                        & 'For tripolar grids, either i0 > ieg or j0 > jeg.  '//trim(emsg_local), FATAL)
                end if
                if ( i0 < isd_filt .or. i0 > ied_filt .or. j0 < jsd_filt .or. j0 > jed_filt ) then
                   write (UNIT=stdout_unit, FMT='("woas.pe,i0,j0= ",3I8,"isd_filt,ied_filt,jsd_filt,jed_filt= ",4I8)')&
                        & mpp_pe(),i0, j0,isd_filt,ied_filt,jsd_filt,jed_filt
                   out_bound_point = out_bound_point + 1
                end if
                if ( Interp%wti(1,1,2) < 1.0 ) then
                   Profiles(num_profiles)%i_index =Interp%i_lon(1,1,1) + Interp%wti(1,1,2)
                else
                   Profiles(num_profiles)%i_index =Interp%i_lon(1,1,2)
                end if
                if ( Interp%wtj(1,1,2) < 1.0 ) then
                   Profiles(num_profiles)%j_index =Interp%j_lat(1,1,1) + Interp%wtj(1,1,2)
                else
                   Profiles(num_profiles)%j_index =Interp%j_lat(1,1,2)
                end if
             end if ! grids
           
             Profiles(num_profiles)%accepted = .true.
             if ( i0 < 1 .or. j0 < 1 ) then
                Profiles(num_profiles)%accepted = .false.
             end if
             if ( i0 < isd_filt .or. i0 >= ied_filt .or. j0 < jsd_filt .or. j0 >= jed_filt ) then
                Profiles(num_profiles)%accepted = .false.
             end if

             if ( Profiles(num_profiles)%accepted ) then ! here
                if ( i0 /= ieg .and. j0 /= jeg ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else if ( i0 == ieg .and. j0 /= jeg ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else if ( i0 /= ieg .and. j0 == jeg ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                else
                   if ( Grd%mask(i0,j0,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                end if
             end if ! here

             if ( Profiles(num_profiles)%accepted ) then ! accepted
                Profiles(num_profiles)%flag(:) = .true.
                allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
                do k=1, Profiles(num_profiles)%levels
                   if (depth(k) < Grd%z(1)) then
                     Profiles(num_profiles)%k_index(k) = 0.0
                   else
                     Profiles(num_profiles)%k_index(k) = frac_index(depth(k), (/0.,Grd%z(:)/)) - 1.0 ! snz modify to v3.2 JAN3012
                   end if
                   if ( Profiles(num_profiles)%k_index(k) > nk ) then
                      call error_mesg('oda_core_mod::open_profile_dataset_woa05s',&
                           & 'Profile k_index is greater than nk', FATAL)
                   end if
                   k0 = floor(Profiles(num_profiles)%k_index(k))

                   if ( k0 >= 1 ) then ! snz add
                      if ( Profiles(num_profiles)%flag(k) ) then ! flag
                         if ( i0 /= ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0+1,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 == ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(1,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0) == 0.0 .or.&
                                 & Grd%mask(1,j0+1,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 /= ieg .and. j0 == jeg ) then
                            if ( Grd%mask(i0,j0,k0) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else
                            if ( Grd%mask(i0,j0,k0) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         end if

                         if ( i0 /= ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0+1,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         elseif ( i0 == ieg .and. j0 /= jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(1,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0,j0+1,k0+1) == 0.0 .or.&
                                 & Grd%mask(1,j0+1,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else if ( i0 /= ieg .and. j0 == jeg ) then
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 .or.&
                                 & Grd%mask(i0+1,j0,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         else
                            if ( Grd%mask(i0,j0,k0+1) == 0.0 ) then
                               Profiles(num_profiles)%flag(k) = .false.
                            end if
                         end if

                         if ( abs(Profiles(num_profiles)%data(k)) > 1.e4 &
                              & .or. abs(Profiles(num_profiles)%depth(k)) > 1.e4 ) then
                            Profiles(num_profiles)%flag(k) = .false.
                         end if
                      end if ! flag
                   end if ! snz add
                end do
             end if ! accepted
           end if ! localize
         end if ! global index
       end do
    end do

    call mpp_close(unit)
    deallocate(axes)
    write (UNIT=stdout_unit, FMT='("A grand total of ",I8," woa05s points within global domain")') no_woa05
    write (UNIT=stdout_unit, FMT='("A final total @woa05s of ",I8," prfs within global domain")') num_profiles

!    call mpp_print_memuse_stats('open_profile_dataset_woa05s Ens')

  end subroutine open_profile_dataset_woa05s

  subroutine get_obs_sst(model_time, Prof, nprof, no_prf0, sst_climo, Filter_domain)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), dimension(:), intent(inout) :: Prof
    integer, intent(inout) :: nprof
    integer, intent(in) :: no_prf0
    type(obs_clim_type), intent(inout) :: sst_climo
    type(domain2d), intent(in) :: Filter_domain

    integer :: i0, j0, i, days, seconds, days1, seconds1
    integer :: unit, ndim, nvar, natt, ntime, time_idx
    integer :: iy0, in0, id0, ih0, im0, is0, i_m
    integer :: stdout_unit, istat
    integer, dimension(12) :: n_days
    integer, save :: year_on_first_read = 0 !< Year on first read of file during
                                            !! module run
    character(len=128) :: sst_filename, emsg_local

    type(fieldtype), dimension(:), allocatable :: fields
    type(time_type) :: tdiff, sst_time0, time1

    n_days  = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    nprof = 0

    stdout_unit = stdout()

    call get_time(model_time, seconds, days)
    call get_date(model_time, iy0, in0, id0, ih0, im0, is0)
    if ( year_on_first_read == 0 ) then
       year_on_first_read = iy0
    end if
    time1=set_date(year_on_first_read, 1, 1, 0, 0, 0)
    call get_time(time1, seconds1, days1)
    time_idx = days-days1+1

    ! daily data 

    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("time_idx = ",I8)') time_idx
    end if

    sst_filename = "INPUT/sst_daily.nc"

    call mpp_open(unit, trim(sst_filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    allocate(fields(nvar), STAT=istat)
    if ( istat .ne. 0 ) then 
       call error_mesg('oda_core_mod::get_obs_sst', 'Unable to allocate fields', FATAL)
    end if

    call mpp_get_fields(unit, fields)
    do i=1, nvar
       select case ( mpp_get_field_name(fields(i)) )
       case ('SST1') ! for AVHRR daily SST
          call mpp_read(unit, fields(i), Filter_domain, sst_climo%sst_obs, tindex=time_idx)
       end select
    end do

    ! get profiles and sst
    ! obs relevant to current analysis interval
    sst_time0 = set_date(iy0, in0, id0, ih0, im0, is0)

    if ( no_sst > 1 ) then

    do i=no_prf+no_woa05+1, no_prf+no_woa05+no_sst
       Profiles(i)%time = sst_time0

       tdiff = model_time - Profiles(i)%time

       i0 = floor(Profiles(i)%i_index)
       if ( i0 < 1 .or. i0 > 1440 ) then
          write (UNIT=emsg_local, FMT='("Profiles(",I8,")%lon = ",I4,", i0 = ",I8)') i, Profiles(i)%lon, i0
          call error_mesg('oda_core_mod::get_obs_sst',&
               & 'Profile longitude index outside range [1,1440].  '//trim(emsg_local), FATAL)
       end if

       j0 = floor(Profiles(i)%j_index)
       if ( j0 < 1 .or. j0 > 1070 ) then
          write (UNIT=emsg_local, FMT='("Profiles(",I8,")%lat = ",I4,", j0 = ",I8)') i, Profiles(i)%lon, j0
          call error_mesg('oda_core_mod::get_obs_sst',&
               & 'Profile latitude index outside range [1,1070].  '//trim(emsg_local), FATAL)
       end if

       if ( Profiles(i)%accepted ) then
          nprof = nprof + 1
          if ( nprof > size(Prof,1) ) then 
             call error_mesg('oda_core_mod::get_obs_sst',&
                  & 'Passed in array "Prof" is smaller than number of profiles, increase size of Prof before call.',&
                  & FATAL)
          end if

          call copy_obs(Profiles(i:i), Prof(no_prf0+nprof:no_prf0+nprof))

          Prof(nprof+no_prf0)%tdiff = tdiff
       end if
    end do

    end if ! for no_sst > 1

    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("no of sst records: ",I8)') nprof
    end if

    deallocate(fields)
    call mpp_close(unit)
  end subroutine get_obs_sst

  subroutine get_obs_woa05t(model_time, Prof, nprof, no_prf0)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), dimension(:), intent(inout) :: Prof
    integer, intent(inout) :: nprof
    integer, intent(in) :: no_prf0

    real :: ri0, rj0, lon_woa05

    integer :: i0, j0
    integer :: i, k, unit, time_idx
    integer :: iy0, in0, id0, ih0, im0, is0
    integer :: ndim, nvar, natt, ntime
    integer :: stdout_unit, istat

    character(len=32) :: axisname
    character(len=128) :: woa05t_filename

    type(fieldtype), dimension(:), allocatable :: fields
    type(time_type) :: tdiff, woa05_time0

    nprof = 0
    stdout_unit = stdout()

    call get_date(model_time, iy0, in0, id0, ih0, im0, is0)

    time_idx = in0

    ! daily data 
    woa05t_filename = "INPUT/woa05_temp.nc"
    
    call mpp_open(unit, trim(woa05t_filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    allocate(fields(nvar), STAT=istat)
    if ( istat .ne. 0 ) then 
       call error_mesg('oda_core_mod::get_obs_woa05t', 'Unable to allocate fields', FATAL)
    end if

    call mpp_get_fields(unit, fields)

    allocate(obs_woa05t(nlon_woa,nlat_woa,nlev_woa))

    do i=1, nvar
       select case ( mpp_get_field_name(fields(i)) )
       case ('t0112an1')
          call mpp_read(unit, fields(i), obs_woa05t, tindex=time_idx)
       end select
    end do

    woa05_time0 = set_date(iy0, in0, id0, ih0, im0, is0)

    do i=no_prf+1, no_prf+no_woa05/2
       Profiles(i)%time = woa05_time0

       tdiff = model_time - Profiles(i)%time

       lon_woa05 = Profiles(i)%lon
       if ( lon_woa05 < 0.0 ) lon_woa05 = lon_woa05 + 360.0
       if ( lon_woa05 > 360.0 ) lon_woa05 = lon_woa05 - 360.0
       ri0 = frac_index(lon_woa05, woa05_lon)
       i0 = floor(ri0)
       if ( i0 < 1 ) i0 = 1
       if ( i0 > nlon_woa ) i0 = nlon_woa

       rj0 = frac_index(Profiles(i)%lat, woa05_lat)
       j0 = floor(rj0)
       if(j0 < 1 ) j0 = 1
       if(j0 > nlat_woa) j0 = nlat_woa

       if ( Profiles(i)%accepted ) then
          nprof = nprof + 1
          if ( nprof > size(Prof,1) ) then
             call error_mesg('oda_core_mod::get_obs_woa05t',&
                  & 'Passed in array "Prof" is smaller than number of profiles, increase size of Prof before call.',&
                  & FATAL)
          end if
          Profiles(i)%data(1:nlev_woa) = obs_woa05t(i0,j0,1:nlev_woa)
          do k=1, nlev_woa
             if ( abs(Profiles(i)%data(k)) > 1.e3 .or.&
                  & abs(Profiles(i)%depth(k)) > 1.e5 ) then
                Profiles(i)%flag(k) = .false.
             end if
          end do
          call copy_obs(Profiles(i:i), Prof(no_prf0+nprof:no_prf0+nprof))
          Prof(no_prf0+nprof)%tdiff = tdiff
       end if
    end do

    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("no of woa05t records: ",I8)') nprof
    end if

    deallocate(fields, obs_woa05t)
    call mpp_close(unit)
  end subroutine get_obs_woa05t

  subroutine get_obs_woa05s(model_time, Prof, nprof, no_prf0)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), dimension(:), intent(inout) :: Prof
    integer, intent(inout) :: nprof
    integer, intent(in) :: no_prf0

    real :: ri0, rj0, lon_woa05

    integer :: i0, j0
    integer :: i, k, unit, time_idx
    integer :: iy0, in0, id0, ih0, im0, is0
    integer :: ndim, nvar, natt, ntime
    integer :: stdout_unit, istat

    character(len=128) :: woa05s_filename

    type(fieldtype), dimension(:), allocatable :: fields
    type(time_type) :: tdiff, woa05_time0

    nprof = 0
    stdout_unit = stdout()

    call get_date(model_time, iy0,in0,id0,ih0,im0,is0)

    time_idx = in0

    ! climatological data 
    woa05s_filename = "INPUT/woa05_salt.nc"

    call mpp_open(unit, trim(woa05s_filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)

    allocate(fields(nvar), STAT=istat)
    if ( istat .ne. 0 ) then 
       call error_mesg('oda_core_mod::get_obs_woa05s', 'Unable to allocate fields', FATAL)
    end if

    call mpp_get_fields(unit, fields)

    allocate(obs_woa05s(nlon_woa,nlat_woa,nlev_woa))

    do i=1, nvar
       select case ( mpp_get_field_name(fields(i)) )
       case ('s0112an1')
          call mpp_read(unit, fields(i), obs_woa05s, tindex=time_idx)
       end select
    end do

    woa05_time0 = set_date(iy0, in0, id0, ih0, im0, is0)

    do i=no_prf+no_woa05/2+1, no_prf+no_woa05
       Profiles(i)%time = woa05_time0

       tdiff = model_time - Profiles(i)%time

       lon_woa05 = Profiles(i)%lon
       if ( lon_woa05 < 0.0 ) lon_woa05 = lon_woa05 + 360.0
       if ( lon_woa05 > 360.0 ) lon_woa05 = lon_woa05 - 360.0
       ri0 = frac_index(lon_woa05, woa05_lon)
       i0 = floor(ri0)
       if ( i0 < 1 ) i0 = 1
       if ( i0 > nlon_woa ) i0 = nlon_woa

       rj0 = frac_index(Profiles(i)%lat, woa05_lat)
       j0 = floor(rj0)
       if ( j0 < 1 ) j0 = 1
       if ( j0 > nlat_woa ) j0 = nlat_woa

       if ( Profiles(i)%accepted ) then
          nprof = nprof + 1
          if ( nprof > size(Prof,1) ) then
             call error_mesg('oda_core_mod::get_obs_woa05s',&
                  & 'Passed in array "Prof" is smaller than number of profiles, increase size of Prof before call',&
                  & FATAL)
          end if
          Profiles(i)%data(1:nlev_woa) = obs_woa05s(i0,j0,1:nlev_woa)
          do k=1, nlev_woa
             if ( abs(Profiles(i)%data(k)) > 1.e3 .or.&
                  & abs(Profiles(i)%depth(k)) > 1.e5 ) then
                Profiles(i)%flag(k) = .false.
             end if
          end do
          call copy_obs(Profiles(i:i),Prof(no_prf0+nprof:no_prf0+nprof))
          Prof(no_prf0+nprof)%tdiff = tdiff
       end if
    end do

    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("no of woa05t records: ",I8)') nprof
    end if

    deallocate(fields, obs_woa05s)
    call mpp_close(unit)
  end subroutine get_obs_woa05s

  subroutine open_profile_dataset_eta(filename, obs_variable, localize)  
    character(len=*), intent(in) :: filename
    integer, intent(in) :: obs_variable
    logical, intent(in), optional :: localize

    integer, parameter :: MAX_LEVELS = 1000

    real :: lon, lat, rms_err
    real :: ri0, rj0
    real, dimension(MAX_LEVELS) :: depth, data

    integer :: inst_type, var_id
    integer :: num_levs, k, kk, i, j, i0, j0
    integer :: stdout_unit

    logical :: data_is_local, localize_data
    logical, dimension(MAX_LEVELS) :: flag

    character(len=32) :: anal_fldname

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    stdout_unit = stdout()

    anal_fldname = 'eta'
    var_id = obs_variable

    !snz idealized
    do j=1, size(x_grid, dim=1)
       do i=1, size(x_grid, dim=2)
          lon = x_grid(i,j)
          lat = y_grid(i,j)
          rms_err = 0.5

          inst_type = 20
          data_is_local = .true.

          if ( lon .lt. 0.0 ) lon = lon + 360.0
          if ( lon .gt. 360.0 ) lon = lon - 360.0
          if ( lon .lt. 80.0 ) lon = lon + 360.0 

          if ( lat < eta_obs_start_lat .or. lat > eta_obs_end_lat ) data_is_local = .false.
          if ( abs(lat) < 20.0 .and.&
               & (mod(floor(lon),2) /= 0 .or. mod(floor(lat),2) /= 0) ) data_is_local = .false.
          if ( (abs(lat) >= 20.0 .and. abs(lat) < 40.0) .and.&
               & (mod(floor(lon),4) /= 0 .or. mod(floor(lat),4) /= 0) ) data_is_local = .false.
          if ( (abs(lat) >= 40.0) .and.&
               & (mod(floor(lon),6) /= 0 .or. mod(floor(lat),6) /= 0) ) data_is_local = .false.

          if ( data_is_local .and. (.NOT.localize_data) ) then
             ri0 = frac_index(lon, x_grid(:,1))
             rj0 = frac_index(lat, y_grid(90,:))
             i0 = floor(ri0)
             j0 = floor(rj0)
             if ( Grd%mask(i0,j0,1) /= 0.0 .and. Grd%mask(i0+1,j0,1) /= 0.0 .and.&
                  & Grd%mask(i0,j0+1,1) /= 0.0 .and. Grd%mask(i0+1,j0+1,1) /= 0.0 ) then
                no_eta = no_eta+1
                num_profiles=num_profiles+1
                if ( num_profiles > max_profiles ) then
                   call error_mesg('oda_core_mod::open_profile_dataset_eta',&
                        & 'Maximum number of profiles exceeded, increase max_profiles in oda_core_nml.', FATAL)
                end if


                num_levs = 0
                flag = .false.
                do k=1, 1
                   flag(k) = .true.
                   data(k) = 0.0
                   depth(k) = 0.0
                   num_levs = num_levs+1
                end do
                if ( num_levs == 0 ) cycle
                allocate(profiles(num_profiles)%depth(num_levs))
                allocate(profiles(num_profiles)%data(num_levs))
                allocate(profiles(num_profiles)%flag(num_levs))
                allocate(profiles(num_profiles)%ms(num_levs))
                allocate(profiles(num_profiles)%ms_inv(num_levs))
                profiles(num_profiles)%variable = var_id
                profiles(num_profiles)%inst_type = inst_type
                profiles(num_profiles)%levels = num_levs
                profiles(num_profiles)%lat = lat
                profiles(num_profiles)%lon = lon
                kk = 1
                do k=1, 1
                   if ( flag(k) ) then
                      profiles(num_profiles)%depth(kk) = depth(k)
                      profiles(num_profiles)%data(kk) = data(k)
                      profiles(num_profiles)%ms(kk) = 1.0
                      profiles(num_profiles)%ms_inv(kk) = 1.0
                      kk = kk + 1
                   end if
                end do
           
                ! calculate interpolation coefficients (make sure to account for grid offsets here!)
                Profiles(num_profiles)%i_index = ri0
                Profiles(num_profiles)%j_index = rj0           
                Profiles(num_profiles)%accepted = .true.
                if ( i0 < 1 .or. j0 < 1 ) then
                   Profiles(num_profiles)%accepted = .false.
                end if
                if ( Profiles(num_profiles)%accepted ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                end if
                if ( Profiles(num_profiles)%accepted ) then
                   Profiles(num_profiles)%flag(:) = .true.
                   allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
                   do k=1, Profiles(num_profiles)%levels
                      if (depth(k) < Grd%z(1)) then
                        Profiles(num_profiles)%k_index(k) = 0.0
                      else
                        Profiles(num_profiles)%k_index(k) = frac_index(depth(k), (/0.,Grd%z(:)/)) - 1
                      end if
                      if ( Profiles(num_profiles)%k_index(k) < 1.0 ) Profiles(num_profiles)%flag(k) = .false.
                   end do
                end if
             end if
          end if
       end do
    end do

    write (UNIT=stdout_unit, FMT='("A grand total of ",I8," eta points within global domain")') no_eta
    write (UNIT=stdout_unit, FMT='("A final total @eta of ",I8," prfs within global domain")') num_profiles
  end subroutine open_profile_dataset_eta

  subroutine open_profile_dataset_suv(filename, obs_variable, localize)  
    character(len=*), intent(in) :: filename
    integer, intent(in) :: obs_variable
    logical, intent(in), optional :: localize

    integer, parameter :: MAX_LEVELS = 1000

    real :: ri0, rj0
    real :: lon, lat, rms_err
    real, dimension(MAX_LEVELS) :: depth, data

    integer :: inst_type, var_id
    integer :: num_levs, k, kk, i, j, i0, j0
    integer :: stdout_unit

    logical :: data_is_local, localize_data
    logical, dimension(MAX_LEVELS) :: flag

    character(len=32) :: anal_fldname

    if ( PRESENT(localize) ) then
       localize_data = localize
    else
       localize_data = .true.
    end if

    stdout_unit = stdout()

    anal_fldname = 'suv'
    var_id = obs_variable

    !snz idealized
    do j=1, size(x_grid_uv,dim=1)
       do i=1, size(x_grid_uv,dim=2)
          lon = x_grid_uv(i,j)
          lat = y_grid_uv(i,j)
          rms_err = 0.5

          inst_type = 20
          data_is_local = .true.

          if ( lon .lt. 0.0 ) lon = lon + 360.0
          if ( lon .gt. 360.0 ) lon = lon - 360.0
          if ( lon .lt. 80.0 ) lon = lon + 360.0 

          if ( lat < -40.0 .or. lat > 40.0 ) data_is_local = .false.
          if ( abs(lat) < 20.0 .and.&
               & (mod(floor(lon),2) /= 0 .or. mod(floor(lat),2) /= 0) ) data_is_local = .false.
          if ( (abs(lat) >= 20.0 .and. abs(lat) < 40.0) .and.&
               & (mod(floor(lon),4) /= 0 .or. mod(floor(lat),4) /= 0) ) data_is_local = .false.
          if ( (abs(lat) >= 40.0) .and.&
               & (mod(floor(lon),6) /= 0 .or. mod(floor(lat),6) /= 0) ) data_is_local = .false.

          if ( data_is_local .and. (.NOT.localize_data) ) then
             ri0 = frac_index(lon, x_grid_uv(:,1))
             rj0 = frac_index(lat, y_grid_uv(90,:))
             i0 = floor(ri0)
             j0 = floor(rj0)
             if ( Grd%mask(i0,j0,1) /= 0.0 .and. Grd%mask(i0+1,j0,1) /= 0.0 .and.&
                  & Grd%mask(i0,j0+1,1) /= 0.0 .and. Grd%mask(i0+1,j0+1,1) /= 0.0 ) then
                no_suv = no_suv+1
                num_profiles = num_profiles + 1
                if ( num_profiles > max_profiles ) then
                   call error_mesg('oda_core_mod::open_profile_dataset_suv',&
                        & 'Maximum number of profiles exceeded, increase max_profiles in oda_core_nml.', FATAL)
                end if


                num_levs = 0
                flag = .false.
                do k=1, 2
                   flag(k) = .true.
                   data(k) = 0.0
                   depth(k) = 0.0
                   num_levs = num_levs + 1
                end do
                if ( num_levs == 0 ) cycle
                allocate(profiles(num_profiles)%depth(num_levs))
                allocate(profiles(num_profiles)%data(num_levs))
                allocate(profiles(num_profiles)%flag(num_levs))
                allocate(profiles(num_profiles)%ms(num_levs))
                allocate(profiles(num_profiles)%ms_inv(num_levs))
                profiles(num_profiles)%variable = var_id
                profiles(num_profiles)%inst_type = inst_type
                profiles(num_profiles)%levels = num_levs
                profiles(num_profiles)%lat = lat
                profiles(num_profiles)%lon = lon
                kk = 1
                do k=1, 2
                   if ( flag(k) ) then
                      profiles(num_profiles)%depth(kk) = depth(k)
                      profiles(num_profiles)%data(kk) = data(k)
                      profiles(num_profiles)%ms(kk) = 1.0
                      profiles(num_profiles)%ms_inv(kk) = 1.0
                      kk = kk + 1
                   end if
                end do
           
                ! calculate interpolation coefficients (make sure to account for grid offsets here!)
                Profiles(num_profiles)%i_index = ri0
                Profiles(num_profiles)%j_index = rj0           
                Profiles(num_profiles)%accepted = .true.
                if ( i0 < 1 .or. j0 < 1 ) then
                   Profiles(num_profiles)%accepted = .false.
                end if
                if ( Profiles(num_profiles)%accepted ) then
                   if ( Grd%mask(i0,j0,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0,1) == 0.0 .or.&
                        & Grd%mask(i0,j0+1,1) == 0.0 .or.&
                        & Grd%mask(i0+1,j0+1,1) == 0.0 ) then
                      Profiles(num_profiles)%accepted = .false.
                   end if
                end if
                if ( Profiles(num_profiles)%accepted ) then
                   Profiles(num_profiles)%flag(:) = .true.
                   allocate(Profiles(num_profiles)%k_index(Profiles(num_profiles)%levels))
                   do k=1, Profiles(num_profiles)%levels
                      if (depth(k) < Grd%z(1)) then
                        Profiles(num_profiles)%k_index(k) = 0.0
                      else
                        Profiles(num_profiles)%k_index(k) = frac_index(depth(k), (/0.,Grd%z(:)/)) - 1.0 ! snz modify to v3.2 JAN3012
                      end if
                      if ( Profiles(num_profiles)%k_index(k) < 1.0 ) Profiles(num_profiles)%flag(k) = .false.
                   end do
                end if
             end if
          end if
       end do
    end do

    write (UNIT=stdout_unit, FMT='("A grand total of ",I8," suv points within global domain")') no_suv
    write (UNIT=stdout_unit, FMT='("A final total @suv of ",I8," prfs within global domain")') num_profiles
  end subroutine open_profile_dataset_suv

  subroutine get_obs_suv(model_time, Prof, nprof, no_prf0)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), dimension(:), intent(inout) :: Prof
    integer, intent(inout) :: nprof
    integer, intent(in) :: no_prf0   

    ! get sst data and put into profiles
    ! only current day                 

    real :: sfc_lon, sfc_lat, ri0, rj0
    real, dimension(1440,1070,1) :: sfc_u, sfc_v

    integer :: i, k, i_m, i0, j0
    integer :: unit, time_idx, ndim, nvar, natt, ntime 
    integer :: iy0, in0, id0, ih0, im0, is0
    integer :: stdout_unit, istat
    integer, dimension(12) :: n_days

    character(len=80) :: sfc_filename
    character(len=256) :: emsg_local

    type(fieldtype), dimension(:), allocatable :: fields
    type(time_type) :: tdiff, sfc_time0

    n_days  = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    nprof = 0
    stdout_unit = stdout()

    call get_date(model_time, iy0, in0, id0, ih0, im0, is0)
                                                                                 
    !monthly    
!!$    time_idx = (iy0-1984)*12+in0
    ! daily    
    if ( in0 == 1 ) then
       time_idx = (iy0-1984)*365 + id0
    else
       time_idx = 0
       do i_m=1, in0-1
          time_idx = time_idx + (iy0-1984)*365 + n_days(i_m)
       end do
       time_idx = time_idx + id0
    end if

    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("time_idx = ",I8)') time_idx
    end if

    sfc_time0 = set_date(iy0, in0, id0, ih0, im0, is0)

    sfc_filename = "INPUT/sfc_current.198401-198412.nc"

    call mpp_open(unit, trim(sfc_filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    
    allocate(fields(nvar), STAT=istat)
    if ( istat .ne. 0 ) then 
       call error_mesg('oda_core_mod::get_obs_sst', 'Unable to allocate fields', FATAL)
    end if

    call mpp_get_fields(unit, fields)
    do i=1, nvar
       select case ( mpp_get_field_name(fields(i)) )
       case ('U_SFC')
          call mpp_read(unit, fields(i), sfc_u, tindex=time_idx)
       case ('V_SFC')
          call mpp_read(unit, fields(i), sfc_v, tindex=time_idx)
       end select
    end do
    
    do i=no_prf+no_sst+no_eta+1, no_prf+no_sst+no_eta+no_suv
       Profiles(i)%time = sfc_time0

       tdiff = model_time - Profiles(i)%time

       sfc_lon = Profiles(i)%lon     
       ri0 = frac_index(sfc_lon, x_grid_uv(:,1))
       i0 = floor(ri0)
       if ( i0 < 1 .or. i0 > 1440 ) then  
          write (UNIT=emsg_local, FMT='("Profiles(",I8,")%lon = ",I4,", i0 = ",I8)') i, Profiles(i)%lon, i0
          call error_mesg('oda_core_mod::get_obs_suv',&
               & 'Profile longitude index outside range [1,1440].  '//trim(emsg_local), FATAL)
       end if
                 
       sfc_lat = Profiles(i)%lat
       rj0 = frac_index(sfc_lat, y_grid_uv(90,:))
       j0 = floor(rj0)
       if ( j0 < 1 .or. j0 > 1070 ) then
          write (UNIT=emsg_local, FMT='("Profiles(",I8,")%lat = ",I4,", j0 = ",I8)') i, Profiles(i)%lon, j0
          call error_mesg('oda_core_mod::get_obs_suv',&
               & 'Profile latitude index outside range [1,1070].  '//trim(emsg_local), FATAL)
       end if

       if ( Profiles(i)%accepted ) then 
          nprof = nprof + 1
          if ( nprof > size(Prof,1) ) then 
             call error_mesg('oda_core_mod::get_obs_suv',&
                  & 'Passed in array "Prof" is smaller than number of profiles, increase size of Prof before call',&
                  & FATAL)
          end if
          Profiles(i)%data(1) = sfc_u(i0,j0,1)
          Profiles(i)%data(2) = sfc_v(i0,j0,1)

          call copy_obs(Profiles(i:i),Prof(nprof+no_prf0:nprof+no_prf0))

          Prof(nprof+no_prf0)%tdiff = tdiff  
       end if
    end do

    deallocate(fields)
    call mpp_close(unit)
  end subroutine get_obs_suv

  subroutine get_obs_eta(model_time, Prof, nprof, no_prf0)
    type(time_type), intent(in) :: model_time
    type(ocean_profile_type), dimension(:), intent(inout) :: Prof
    integer, intent(inout) :: nprof
    integer, intent(in) :: no_prf0   

    ! get sst data and put into profiles
    ! only current day                 

    real :: eta_lon, eta_lat, ri0, rj0
    real, dimension(1440,1070) :: eta_t

    integer :: i, i0, j0, i_m
    integer :: iy0, in0, id0, ih0, im0, is0
    integer :: unit, time_idx, ndim, nvar, natt, ntime
    integer :: stdout_unit, istat
    integer, dimension(12) :: n_days

    character(len=80) :: eta_filename
    character(len=256) :: emsg_local
    
    type(fieldtype), dimension(:), allocatable :: fields
    type(time_type) :: tdiff, sfc_time0

    n_days  = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    nprof = 0
    stdout_unit = stdout()

    call get_date(model_time, iy0, in0, id0, ih0, im0, is0)

    !monthly    
!!$    time_idx = (iy0-1984)*12+in0
    ! daily    
    if ( in0 == 1 ) then
       time_idx = (iy0-1976)*365 + id0 - 1
    else
       time_idx = 0
       do i_m=1, in0-1
          time_idx = time_idx + n_days(i_m)
       end do
       time_idx = (iy0-1976)*365 + time_idx + id0 - 1
    end if

    if ( mpp_pe() == mpp_root_pe() ) then
       write (UNIT=stdout_unit, FMT='("time_idx = ",I8)') time_idx
    end if

    sfc_time0 = set_date(iy0, in0, id0, ih0, im0, is0)

    eta_filename='INPUT/ocean.19760101-20001231.eta_t.nc'

    call mpp_open(unit, trim(eta_filename), MPP_RDONLY, MPP_NETCDF, threading=MPP_MULTI, fileset=MPP_SINGLE)
    call mpp_get_info(unit, ndim, nvar, natt, ntime)
    
    allocate(fields(nvar), STAT=istat)
    if ( istat .ne. 0 ) then 
       call error_mesg('oda_core_mod::get_obs_sst', 'Unable to allocate fields', FATAL)
    end if

    call mpp_get_fields(unit, fields)
    do i=1, nvar
       select case ( mpp_get_field_name(fields(i)) )
       case ('eta_t')
          call mpp_read(unit, fields(i), eta_t, tindex=time_idx)
       end select
    end do

    do i=no_prf+no_sst+1, no_prf+no_sst+no_eta
       Profiles(i)%time = sfc_time0

       tdiff = model_time - Profiles(i)%time

       eta_lon = Profiles(i)%lon
       ri0 = frac_index(eta_lon, x_grid(:,1))
       i0 = floor(ri0)
       if ( i0 < 1 .or. i0 > 1440 ) then  
          write (UNIT=emsg_local, FMT='("Profiles(",I8,")%lon = ",I4,", i0 = ",I8)') i, Profiles(i)%lon, i0
          call error_mesg('oda_core_mod::get_obs_eta',&
               & 'Profile longitude index outside range [1,1440].  '//trim(emsg_local), FATAL)
       end if
                 
       eta_lat = Profiles(i)%lat
       rj0 = frac_index(eta_lat, y_grid(90,:))
       j0 = floor(rj0)
       if ( j0 < 1 .or. j0 > 1070 ) then
          write (UNIT=emsg_local, FMT='("Profiles(",I8,")%lat = ",I4,", j0 = ",I8)') i, Profiles(i)%lon, j0
          call error_mesg('oda_core_mod::get_obs_suv',&
               & 'Profile latitude index outside range [1,1070].  '//trim(emsg_local), FATAL)
       end if

       if ( Profiles(i)%accepted ) then 
          if ( eta_t(i0,j0) > -9.9 ) then !!! excluding missing values
             nprof = nprof + 1              
             if ( nprof > size(Prof,1) ) then 
                call error_mesg('oda_core_mod::get_obs_eta',&
                     & 'Passed in array "Prof" is smaller than number of profiles, increase size of Prof before call.',&
                     & FATAL)
             end if
             Profiles(i)%data(1) = eta_t(i0,j0)

             call copy_obs(Profiles(i:i),Prof(nprof+no_prf0:nprof+no_prf0))

             Prof(nprof+no_prf0)%tdiff = tdiff  
          end if
       end if
    end do

    deallocate(fields)
    call mpp_close(unit)
  end subroutine get_obs_eta
end module oda_core_ecda_mod
