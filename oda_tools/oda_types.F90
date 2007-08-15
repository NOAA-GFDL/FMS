module oda_types_mod
#define MAX_LEVS_FILE_ 50
  
!============================================================
! This module contains type declarations and default values
! for oda modules.  
!============================================================
  
! Contact: Matthew.Harrison@gfdl.noaa.gov

  use time_manager_mod, only : time_type, set_time
  use mpp_mod, only : stdout
  use mpp_domains_mod, only : domain2d
  
  implicit none

  private

! Controls record length for optimal storage  
  integer, parameter, public :: max_levels_file=MAX_LEVS_FILE_
! Maximum number of neighbors for QC or analysis
  integer, parameter, public :: max_neighbors=100 ! for profiles
! Maximum number of records per profile for storage   
  integer, parameter, public :: max_links=100 ! for profiles  
  
! List of variables for ODA   

  real, parameter, public :: missing_value=-1.e20
  
  type, public :: forward_model_type
     real, dimension(:,:,:,:), pointer :: wgt=>NULL() ! interpolation weights
  end type forward_model_type
  
  type, public :: ocean_profile_type
     integer :: nvar
     real    :: project ! e.g. FGGE, COARE, ACCE, ...
     real    :: probe ! MBT, XBT, drifting buoy
     real    :: ref_inst ! instrument (thermograph, hull sensor, ...)
     integer :: wod_cast_num
     real    :: fix_depth
     real    :: ocn_vehicle
     real    :: database_id
     integer :: levels
     integer :: profile_flag ! an overall flag for the profile
     integer :: profile_flag_s ! an overall flag for the profile salinity     
     real :: lat, lon
     logical :: accepted
     integer :: nlinks
     type(ocean_profile_type), pointer, dimension(:) :: next ! Large profiles are stored as linked list.
     integer, dimension(max_neighbors) :: nbr_index
     real, dimension(max_neighbors) :: nbr_dist ! distance in radians 
     real, dimension(:), pointer :: depth, data_t, data_s
     integer, dimension(:), pointer :: flag_t ! level-by-level flags for temp
     integer, dimension(:), pointer :: flag_s ! level-by-level flags for salinity
     real    :: temp_err, salt_err ! measurement error
     real, dimension(:), pointer :: ms_t      ! ms temperature by level
     real, dimension(:), pointer :: ms_s      ! ms salinity by level     
     type(time_type) :: time 
     integer         :: yyyy
     integer         :: mmdd
     type(time_type), pointer :: Model_time ! each profile can be associated with a first-guess field with an associated time and grid
     type(grid_type), pointer :: Model_grid 
     real :: i_index, j_index ! model longitude and latitude indices respectively
     real, dimension(:), pointer :: k_index     ! model depth indices
     type(forward_model_type) :: Forward_model  ! linear operation from model to observation
     type(time_type) :: tdiff      ! positive difference between model time and observation time
  end type ocean_profile_type

  type, public :: ocean_surface_type
     integer :: variable  ! variable ids are defined by the ocean_types module (e.g. TEMP_ID, SALT_ID, ...)
     integer :: inst_type  ! instrument types are defined by platform class (e.g. MOORING, DROP) and instrument type (XBT, CTD, ...)
     integer :: qc_flag, nobs
     logical :: is_gridded
     integer :: nlon, nlat
     real, pointer, dimension(:) :: lat, lon
     logical :: accepted     
     real, pointer, dimension(:) :: data
     real, dimension(:), pointer :: ms     =>NULL() 
     real, dimension(:), pointer :: i_index=>NULL() , j_index=>NULL()  ! model indices
     real, pointer, dimension(:,:) :: data2=>NULL() 
     real, dimension(:,:), pointer :: ms2     =>NULL() 
     real, dimension(:,:), pointer :: i_index2=>NULL() , j_index2=>NULL()  ! model indices
     real :: k_index          
     type(forward_model_type) :: Forward_model
     type(time_type) :: time
     integer :: yyyy
     integer :: mmdd
     character(len=8) :: wmo_id
     type(time_type), pointer :: Model_time=>NULL() 
     type(grid_type), pointer :: Model_grid=>NULL()    
     ! positive difference between current model time 
     ! and observation time
     type(time_type) :: tdiff
  end type ocean_surface_type


  type, public :: grid_type
     real, pointer, dimension(:,:) :: x=>NULL() , y=>NULL() 
     real, pointer, dimension(:,:) :: x_bound=>NULL() , y_bound=>NULL()      
     real, pointer, dimension(:,:) :: dx=>NULL() , dy=>NULL() 
     real, pointer, dimension(:) :: z=>NULL() , z_bound=>NULL() 
     real, pointer, dimension(:) :: dz=>NULL() 
     real, pointer, dimension(:,:,:) :: mask=>NULL()
     type(domain2d), pointer :: Dom ! FMS domain type
     logical :: cyclic
     integer :: ni, nj, nk
  end type grid_type

  type, public :: field_type
     type(grid_type) :: grid
     real, pointer, dimension(:,:,:) :: data=>NULL() 
  end type field_type


  type, public :: field_dist_type_3d
     integer :: error_model
     character(len=32) :: name
     type(grid_type), pointer :: grid=>NULL() 
     real, pointer, dimension(:,:,:) :: ex=>NULL() , vr=>NULL()
     real, pointer, dimension(:,:,:) :: obs_d=>NULL() ! obs minus expected value
  end type field_dist_type_3d

  type, public :: field_dist_type_2d
     integer :: error_model
     character(len=32) :: name
     type(grid_type), pointer :: grid=>NULL() 
     real, pointer, dimension(:,:) :: ex=>NULL() , vr=>NULL() 
  end type field_dist_type_2d
     
  type, public :: ocean_dist_type
     type(field_dist_type_3d) :: temp,salt,u,v
     type(field_dist_type_2d) :: eta
  end type ocean_dist_type

  public init_obs
  
  interface init_obs
     module procedure init_obs_profile
  end interface
  
  contains

    subroutine oda_types_init()

      use fms_mod, only : open_namelist_file, check_nml_error, close_file
      

    end subroutine oda_types_init


    subroutine init_obs_profile(profile)

      type(ocean_profile_type), intent(inout) :: profile

      profile%nvar = 0
      profile%project = -1.0
      profile%probe   = -1.0
      profile%wod_cast_num = -1
      profile%ref_inst = -1.0
      profile%fix_depth = -1.0
      profile%ocn_vehicle = -1.0
      profile%database_id = -1.0
      profile%levels = 0
      profile%profile_flag = 0
      profile%profile_flag_s = 0
      profile%lat = -1.e10
      profile%lon = -1.e10
      profile%accepted = .true.
      if (ASSOCIATED(profile%next)) deallocate(profile%next)
      profile%nlinks = 0
      profile%nbr_index(:) = -1
      profile%nbr_dist(:) = -1.0
      if (ASSOCIATED(profile%depth)) deallocate(profile%depth)
      if (ASSOCIATED(profile%data_t)) deallocate(profile%data_t)
      if (ASSOCIATED(profile%data_s)) deallocate(profile%data_s)
      if (ASSOCIATED(profile%flag_t)) deallocate(profile%flag_t)
      if (ASSOCIATED(profile%flag_s)) deallocate(profile%flag_s)
      if (ASSOCIATED(profile%ms_t)) deallocate(profile%ms_t)
      if (ASSOCIATED(profile%ms_s)) deallocate(profile%ms_s)
      profile%temp_err = -1.0
      profile%salt_err = -1.0
      profile%time = set_time(0,0)
      profile%yyyy = 0
      profile%mmdd = 0      
      if (ASSOCIATED(profile%model_time)) deallocate(profile%model_time)
      if (ASSOCIATED(profile%model_grid)) deallocate(profile%model_grid)
      profile%i_index = -1
      profile%j_index = -1
      if (ASSOCIATED(profile%k_index)) deallocate(profile%k_index)
      profile%tdiff = set_time(0,0)

      return
      
    end subroutine init_obs_profile
    
end module oda_types_mod
