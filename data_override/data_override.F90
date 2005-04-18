module data_override_mod
!
! <CONTACT EMAIL="Giang.Nong@noaa.gov">
! G.T. Nong
! </CONTACT>
!
! <CONTACT EMAIL="Matthew.Harrison@noaa.gov">
!  M.J. Harrison 
! </CONTACT>
!
! <CONTACT EMAIL="Michael.Winton@noaa.gov">
! M. Winton
! </CONTACT>

!<OVERVIEW>
! Given a gridname, fieldname and model time this routine will get data in a file whose
! path is described in a user-provided data_table, do spatial and temporal interpolation if 
! necessary to convert data to model's grid and time.
!
! Before using data_override a data_table must be created with the following entries:
! gridname, fieldname_code, fieldname_file, file_name, ongrid, factor.
!
! More explainations about data_table entries can be found in the source code (defining data_type)
!
! If user wants to override fieldname_code with a const, set fieldname_file in data_table = ""
! and factor = const
!
! If user wants to override fieldname_code with data from a file, set fieldname_file = name in
! the netCDF data file, factor then will be for unit conversion (=1 if no conversion required)
!
! A field can be overriden globally (by default) or users can specify a region in which
! data_override will take place, field values outside the region will not be affected. 
!</OVERVIEW>
use  platform_mod, only: r8_kind
use constants_mod, only: PI
use mpp_io_mod, only: axistype, mpp_close, mpp_open, mpp_get_axis_data, MPP_RDONLY
use mpp_mod, only : mpp_error,FATAL,mpp_pe, stdout, stdlog
use horiz_interp_mod, only : horiz_interp, horiz_interp_init, horiz_interp_type
use time_interp_external_mod, only:time_interp_external_init, time_interp_external, &
     init_external_field, get_external_field_size
use fms_io_mod, only: field_size, read_data, write_data,fms_io_init,nullify_domain,return_domain, &
     set_domain
use fms_mod, only: write_version_number, file_exist, field_exist
use axis_utils_mod, only: get_axis_bounds
use mpp_domains_mod, only : domain2d, mpp_get_compute_domain, NULL_DOMAIN2D,operator(.NE.),operator(.EQ.)
use time_manager_mod, only: time_type

implicit none
private

character(len=128) :: version = '$Id: data_override.F90,v 12.0 2005/04/14 17:55:22 fms Exp $'
character(len=128) :: tagname = '$Name: lima $'

type data_type
   character(len=3) :: gridname
   character(len=128) :: fieldname_code !fieldname used in user's code (model)
   character(len=128) :: fieldname_file ! fieldname used in the netcdf data file
   character(len=128) :: file_name   ! name of netCDF data file
   logical :: ongrid   ! true if data is on model's grid, false otherwise
   real :: factor ! For unit conversion, default=1, see OVERVIEW above
end type data_type

type override_type
   character(len=3) :: gridname  
   character(len=128) :: fieldname
   integer :: t_index !index for time interp
   type(horiz_interp_type) :: horz_interp ! index for horizontal spatial interp
   integer :: dims(4) ! dimensions(x,y,z,t) of the field stored in filename
   integer :: comp_domain(4) ! istart,iend,jstart,jend for compute domain
   integer :: region1(4) ! istart,iend,jstart,jend for region1
   integer :: region2(4) ! istart,iend,jstart,jend for region2
end type override_type

 integer, parameter :: max_table=100, max_array=100
 integer :: table_size ! actual size of data table
 integer, parameter :: ANNUAL=1, MONTHLY=2, DAILY=3, HOURLY=4, UNDEF=-1
 real :: deg_to_radian, radian_to_deg 
 logical:: module_is_initialized = .FALSE.

type(domain2D),save :: ocn_domain,atm_domain,lnd_domain, ice_domain 
real(r8_kind), dimension(:,:), allocatable :: glo_lat_ocn, glo_lon_ocn, glo_lat_atm, &
     glo_lon_atm, glo_lat_lnd, glo_lon_lnd
real, dimension(:,:), target, allocatable :: lon_local_ocn, lat_local_ocn,lon_local_atm,&
     lon_local_ice, lon_local_lnd, lat_local_atm, lat_local_ice, lat_local_lnd
integer:: num_fields = 0 ! number of fields in override_array already processed
type(data_type), dimension(max_table) :: data_table ! user-provided data table
type(data_type) :: default_table
type(override_type), dimension(max_array), save :: override_array ! to store processed fields
type(override_type), save :: default_array
logical :: atm_on, ocn_on, lnd_on, ice_on

interface data_override
     module procedure data_override_2d
     module procedure data_override_3d
end interface

public :: data_override_init, data_override

contains

! <SUBROUTINE NAME="data_override_init">
!   <DESCRIPTION>
! Assign default values for default_table, get domain of component models,
! get global grids of component models.
! Users should call data_override_init before calling data_override
!   </DESCRIPTION>
!   <TEMPLATE>
! call data_override_init
!   </TEMPLATE>
subroutine data_override_init(Atm_domain_in, Ocean_domain_in, Ice_domain_in, Land_domain_in)
  type (domain2d), intent(in), optional :: Atm_domain_in
  type (domain2d), intent(in), optional :: Ocean_domain_in, Ice_domain_in 
  type (domain2d), intent(in), optional :: Land_domain_in

! <NOTE>
! This subroutine should be called in coupler_init after
! (ocean/atmos/land/ice)_model_init have been called.
!
! data_override_init can be called more than once, in one call the user can pass 
! up to 4 domains of component models, at least one domain must be present in
! any call
!
! Data_table is initialized here with default values. Users should provide "real" values
! that will override the default values. Real values can be given using data_table, each
! line of data_table contains one data_entry. Items of data_entry are comma separated.
!
! </NOTE>
  integer :: is,ie,js,je
  integer :: i, iunit, io_status, ntable, ierr, ioun
  character(len=256) :: record
  type(data_type) :: data_entry
  type (domain2d) :: domain2 ! It should not be necessary to save and restore the
                             ! the current_domain of fms_io_mod because it should
                             ! not have a current_domain. domain should be a required
                             ! argument of read_data and write_data rather than an
                             ! optional argument.

!  if(module_is_initialized) return

  atm_on = PRESENT(Atm_domain_in)
  ocn_on = PRESENT(Ocean_domain_in)
  lnd_on = PRESENT(Land_domain_in)
  ice_on = PRESENT(Ice_domain_in)
   
  if(.not. module_is_initialized) then
    radian_to_deg = 180./PI
    deg_to_radian = PI/180.

    call write_version_number (version, tagname)

    domain2 = NULL_DOMAIN2D     ! See comment above
    call return_domain(domain2) ! See comment above
    call fms_io_init
    call nullify_domain
!    get global lat and lon of all three model grids
    call get_global_grid()
    if (domain2 .NE. NULL_DOMAIN2D) call set_domain(domain2) ! See comment above

    atm_domain = NULL_DOMAIN2D
    ocn_domain = NULL_DOMAIN2D
    lnd_domain = NULL_DOMAIN2D
    ice_domain = NULL_DOMAIN2D 

!2 Initialize user-provided data table  
    default_table%gridname = 'none'
    default_table%fieldname_code = 'none'
    default_table%fieldname_file = 'none'
    default_table%file_name = 'none'
    default_table%ongrid = .FALSE.
    default_table%factor = 1.
    do i = 1,max_table
       data_table(i) = default_table
    enddo

!3 Read coupler_table 
    call mpp_open(iunit, 'data_table', action=MPP_RDONLY)
    ntable = 0
    do while (ntable <= max_table)
       read(iunit,'(a)',end=100) record
       ntable=ntable+1     
       if (record(1:1) == '#') cycle
       if (record(1:10) == '          ') cycle
       read(record,*,err=99) data_entry         
       data_table(ntable) = data_entry
    enddo
    call mpp_error(FATAL,'too many enries in data_table')
99  call mpp_error(FATAL,'error in data_table format')
100 continue
    table_size = ntable
    call mpp_close(iunit)
!4 Initialize override array
    default_array%gridname = 'NONE'
    default_array%fieldname = 'NONE'
    default_array%t_index = -1
    default_array%dims = -1
    default_array%comp_domain = -1
    default_array%region1 = -1
    default_array%region2 = -1
    do i = 1, max_array
       override_array(i) = default_array
    enddo
    call time_interp_external_init
 endif

 if (atm_on) atm_domain = Atm_domain_in
 if (ocn_on) ocn_domain = Ocean_domain_in
 if (lnd_on) lnd_domain = Land_domain_in
 if (ice_on) ice_domain = Ice_domain_in 
! compute lon and lat local of 4 component model
  if (atm_on) then
      call mpp_get_compute_domain( atm_domain,is,ie,js,je) 
      allocate(lon_local_atm(is:ie,js:je), lat_local_atm(is:ie,js:je))
      lon_local_atm(is:ie,js:je) = glo_lon_atm(is:ie,js:je)
      lat_local_atm(is:ie,js:je) = glo_lat_atm(is:ie,js:je) 
  endif

  if (ocn_on) then
      call mpp_get_compute_domain( ocn_domain,is,ie,js,je) 
      allocate(lon_local_ocn(is:ie,js:je), lat_local_ocn(is:ie,js:je))
      lon_local_ocn(is:ie,js:je) = glo_lon_ocn(is:ie,js:je)
      lat_local_ocn(is:ie,js:je) = glo_lat_ocn(is:ie,js:je)
  endif

  if (lnd_on) then
      call mpp_get_compute_domain( lnd_domain,is,ie,js,je) 
      allocate(lon_local_lnd(is:ie,js:je), lat_local_lnd(is:ie,js:je))
      lon_local_lnd(is:ie,js:je) = glo_lon_lnd(is:ie,js:je)
      lat_local_lnd(is:ie,js:je) = glo_lat_lnd(is:ie,js:je)
  endif


  if (ice_on) then
      call mpp_get_compute_domain( ice_domain,is,ie,js,je) 
      allocate(lon_local_ice(is:ie,js:je), lat_local_ice(is:ie,js:je))
      lon_local_ice(is:ie,js:je) = glo_lon_ocn(is:ie,js:je)
      lat_local_ice(is:ie,js:je) = glo_lat_ocn(is:ie,js:je)
  endif    

  module_is_initialized = .TRUE.
 
end subroutine data_override_init
! </SUBROUTINE>
! #####
subroutine get_domain(gridname, domain, comp_domain)
! Given a gridname, this routine returns the working domain associated with this gridname
  character(len=3), intent(in) :: gridname
  type(domain2D), intent(inout) :: domain
  integer, intent(out), optional :: comp_domain(4) ! istart,iend,jstart,jend for compute domain

  domain = NULL_DOMAIN2D
  select case (gridname)
     case('OCN')
        domain = ocn_domain
     case('ATM')       
        domain = atm_domain
     case('LND')       
        domain = lnd_domain
     case('ICE')        
        domain = ice_domain
     case default
        call mpp_error(FATAL,'error in data_override get_domain')
  end select
  if(domain .EQ. NULL_DOMAIN2D) call mpp_error(FATAL,'data_override: failure in get_domain')
  if(present(comp_domain)) &
     call mpp_get_compute_domain(domain,comp_domain(1),comp_domain(2),comp_domain(3),comp_domain(4)) 
end subroutine get_domain
! #####

! <SUBROUTINE NAME="data_override_2d">
!   <DESCRIPTION>
! This routine performs data override for 2D fields; for usage, see data_override_3d.
!   </DESCRIPTION>
subroutine data_override_2d(gridname,fieldname,data_2D,time,override,region1,region2)
  character(len=3), intent(in) :: gridname ! model grid ID
  character(len=*), intent(in) :: fieldname ! field to override
  logical, intent(out), optional :: override ! true if the field has been overriden succesfully
  real, intent(in), optional :: region1(4),region2(4) !lat and lon of region where override is done
  type(time_type), intent(in) :: time !  model time
  real, dimension(:,:), intent(inout) :: data_2D !data returned by this call
!  real, dimension(size(data_2D,1),size(data_2D,2),1) :: data_3D
  real, dimension(:,:,:), allocatable ::  data_3D
  integer       :: index1
  integer       :: i

!1  Look  for the data file in data_table 
  if(PRESENT(override)) override = .false.
  index1 = -1
  do i = 1, table_size
     if( trim(gridname) /= trim(data_table(i)%gridname)) cycle
     if( trim(fieldname) /= trim(data_table(i)%fieldname_code)) cycle
     index1 = i                               ! field found        
     exit
  enddo
  if(index1 .eq. -1) return  ! NO override was performed

  allocate(data_3D(size(data_2D,1),size(data_2D,2),1))
  data_3D(:,:,1) = data_2D
  if(present(override) .and. PRESENT(region1).and. present(region2)) then
     call data_override_3d(gridname,fieldname,data_3D,time,override,region1,region2,data_index=index1)
  else if (PRESENT(region1).and. present(region2)) then
     call data_override_3d(gridname,fieldname,data_3D,time,region1=region1,region2=region2,data_index=index1)
  else if(present(override) .and. PRESENT(region1)) then
     call data_override_3d(gridname,fieldname,data_3D,time,override,region1,data_index=index1)
  else if(present(override)) then
     call data_override_3d(gridname,fieldname,data_3D,time,override,data_index=index1)
  else if (present(region1)) then
     call data_override_3d(gridname,fieldname,data_3D,time,region1=region1,data_index=index1)
  else
     call data_override_3d(gridname,fieldname,data_3D,time,data_index=index1)
  endif  
     
  data_2D(:,:) = data_3D(:,:,1)
  deallocate(data_3D)
end subroutine data_override_2d
! </SUBROUTINE>
! ####

! <SUBROUTINE NAME="data_override_3d">
!   <DESCRIPTION>
! This routine performs data override for 3D fields
!   <TEMPLATE>
! call data_override(gridname,fieldname,data,time,override)
!   </TEMPLATE>
!   </DESCRIPTION>

!   <IN NAME="gridname"  TYPE="character" DIM="(*)">
! Grid name (Ocean, Ice, Atmosphere, Land)
!   </IN>
!   <IN NAME="fieldname_code" TYPE="character" DIM="(*)">
!    Field name as used in the code (may be different from the name in NetCDF data file)
!   </IN>
!   <OUT NAME="data" TYPE="real" DIM="(:,:,:)">
!    array containing output data
!   </OUT>
!   <IN NAME="time" TYPE="time_type">
!    model time
!   </IN>
!   <OUT NAME="override" TYPE="logical">
!    TRUE if the field is overriden, FALSE otherwise
!   </OUT>
subroutine data_override_3d(gridname,fieldname_code,data1,time,override,region1,region2,data_index)
  character(len=3), intent(in) :: gridname ! model grid ID
  character(len=*), intent(in) :: fieldname_code ! field name as used in the model
  logical, intent(out), optional :: override ! true if the field has been overriden succesfully
  type(time_type), intent(in) :: time !(target) model time
  real, intent(in), optional :: region1(4),region2(4) !lat and lon of regions where override is done
!Note: region2 can not exist without region1. In other words, if only one region is specified, it
! should be region1
  integer, intent(in), optional :: data_index
  real, dimension(:,:,:), intent(out) :: data1 !data returned by this call
  real, dimension(:,:,:), allocatable :: data !temporary array for data
  character(len=128) :: filename !file containing source data
  character(len=128) :: fieldname ! fieldname used in the data file
  integer :: i,j
  integer :: dims(4)
  integer :: index1 ! field index in data_table
  integer :: id_time !index for time interp in override array
  integer :: axis_sizes(4)
  type(horiz_interp_type) :: id_horz_interp   !index for horizontal interp 
  real, dimension(:),allocatable :: lon_in, lat_in !of the input (source) grid
  real, dimension(:,:), pointer :: lon_local =>NULL(), &
                                   lat_local =>NULL() !of output (target) grid cells

  type(axistype) :: axis_centers(4), axis_bounds(4)
  logical :: bilinear_interp = .true. 
  logical :: data_file_is_2D = .false.  !data in netCDF file is 2D
  logical :: ongrid
  type(domain2D) :: domain
  integer :: curr_position ! position of the field currently processed in override_array
  real :: factor
  integer, dimension(4) :: comp_domain = 0  ! istart,iend,jstart,jend for compute domain
  integer, dimension(4) :: reg1 = -1         ! istart,iend,jstart,jend for region1
  integer, dimension(4) :: reg2 = -1         ! istart,iend,jstart,jend for region2
  if(.not.module_is_initialized) &
       call mpp_error(FATAL,'Error: need to call data_override_init first')

!1  Look  for the data file in data_table 
  if(PRESENT(override)) override = .false.
  if (present(data_index)) then
    index1 = data_index
  else
    index1 = -1
    do i = 1, table_size
       if( trim(gridname) /= trim(data_table(i)%gridname)) cycle
       if( trim(fieldname_code) /= trim(data_table(i)%fieldname_code)) cycle
       index1 = i                               ! field found        
       exit
    enddo
    if(index1 .eq. -1) return  ! NO override was performed
  endif
 
  if(present(region2) .and. .not. present(region1)) &
       call mpp_error(FATAL,'data_override: region2 is specified without region1')

  fieldname = data_table(index1)%fieldname_file ! fieldname in netCDF data file
  factor = data_table(index1)%factor



  if(fieldname == "") then
     data1 = factor
     if(PRESENT(override)) override = .true.
     return
  else
     filename = data_table(index1)%file_name
     if (filename == "") call mpp_error(FATAL,'data_override: filename not given in data_table')
  endif  
  ongrid = data_table(index1)%ongrid

!3 Check if fieldname has been previously processed
  curr_position = -1
  if(num_fields > 0 ) then
     do i = 1, num_fields
        if(trim(override_array(i)%gridname) /= trim(gridname))   cycle 
        if(trim(override_array(i)%fieldname) /= trim(fieldname_code)) cycle
        curr_position = i
        exit        
     enddo
  endif

  if(curr_position < 0) then ! the field has not been processed previously
! Get working domain from model's gridname
     if(present(region1)) call get_region_bounds(gridname,region1,reg1) 
     call get_domain(gridname,domain,comp_domain)                          
     if(present(region2)) call get_region_bounds(gridname,region1,reg2)
     num_fields = num_fields + 1
     curr_position = num_fields     
! record fieldname, gridname in override_array    
     override_array(curr_position)%fieldname = fieldname_code
     override_array(curr_position)%gridname = gridname
     override_array(curr_position)%comp_domain = comp_domain
     if(present(region1)) override_array(curr_position)%region1 = reg1          
     if(present(region2)) override_array(curr_position)%region2 = reg2
!4 get index for time interp   
     if(ongrid) then
        id_time = init_external_field(filename,fieldname,domain=domain,verbose=.false.)
        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 1') 
        override_array(curr_position)%t_index = id_time     
     else !ongrid=false
        id_time = init_external_field(filename,fieldname,domain=domain, axis_centers=axis_centers,&
             axis_sizes=axis_sizes, verbose=.false.,override=.true.)  
        dims = get_external_field_size(id_time)
        override_array(curr_position)%dims = dims
        if(id_time<0) call mpp_error(FATAL,'data_override:field not found in init_external_field 2')
        override_array(curr_position)%t_index = id_time
        
!5 Get local lon and lat of model grid
        select case(gridname)
        case('OCN')          
           lon_local => lon_local_ocn; lat_local => lat_local_ocn
        case('ICE')
           lon_local => lon_local_ice; lat_local => lat_local_ice
        case('ATM')
           lon_local => lon_local_atm; lat_local => lat_local_atm
        case('LND')
           lon_local => lon_local_lnd; lat_local => lat_local_lnd
        case default
           call mpp_error(FATAL,'error: gridname not recognized in data_override')
        end select

!7 get lon and lat of the input (source) grid, assuming that axis%data contains
!  lat and lon of the input grid (in degrees)
        call get_axis_bounds(axis_centers(1),axis_bounds(1), axis_centers)
        call get_axis_bounds(axis_centers(2),axis_bounds(2), axis_centers)
        allocate(lon_in(axis_sizes(1)+1), lat_in(axis_sizes(2)+1))
        call mpp_get_axis_data(axis_bounds(1),lon_in)
        call mpp_get_axis_data(axis_bounds(2),lat_in)
! convert lon_in and lat_in from deg to radian
        lon_in = lon_in * deg_to_radian
        lat_in = lat_in * deg_to_radian
!8 do horizontal_interp_init to get id_horz_interp 
        call horiz_interp_init (id_horz_interp,lon_in, lat_in, lon_local, lat_local,&
             interp_method="bilinear")
        override_array(curr_position)%horz_interp = id_horz_interp
        deallocate(lon_in)
        deallocate(lat_in)
     endif !(ongrid)
  else !curr_position >0
     dims = override_array(curr_position)%dims
     if(present(region1)) reg1 = override_array(curr_position)%region1     
     if(present(region2)) reg2 = override_array(curr_position)%region2 
     comp_domain = override_array(curr_position)%comp_domain
!9 Get id_time  previously stored in override_array
     id_time = override_array(curr_position)%t_index
! if ongrid == .false. need to get id_horz_interp
     if(.not. ongrid) &
          id_horz_interp = override_array(curr_position)%horz_interp
  endif !if curr_position < 0

  allocate(data(comp_domain(1):comp_domain(2),comp_domain(3):comp_domain(4),size(data1,3)))
  data = HUGE(1.0)
  ! Determine if  data in netCDF file is 2D or not  
  data_file_is_2D = .false.
  if((dims(3) == 1) .and. (size(data1,3)>1)) data_file_is_2D = .true. 

  if(ongrid) then
!10 do time interp to get data in compute_domain and return     
     if(data_file_is_2D) then        
        call time_interp_external(id_time,time,data(:,:,1),verbose=.false.)
        data(:,:,1) = data(:,:,1)*factor
        do i = 2, size(data,3)
           data(:,:,i) = data(:,:,1)
        enddo
     else
        call time_interp_external(id_time,time,data,verbose=.false.)
        data = data*factor
     endif    

     if(present(region1)) then
        do i = comp_domain(1), comp_domain(2)
           do j = comp_domain(3), comp_domain(4)
              if(i>=reg1(1) .and. i<=reg1(2) .and. j>=reg1(3) .and. j<=reg1(4)) &
                 data1(i,j,:) = data(i,j,:)
           enddo
        enddo
     else
        data1 = data
     endif

     if(present(region2)) then
        do i = comp_domain(1), comp_domain(2)
           do j = comp_domain(3), comp_domain(4)
              if(i>=reg2(1) .and. i<=reg2(2) .and. j>=reg2(3) .and. j<=reg2(4)) &
                 data1(i,j,:) = data(i,j,:)
           enddo
        enddo
     endif
     if(PRESENT(override)) override = .true.
     deallocate(data)
     return

  else  ! off grid case
! do time interp to get global data
     if(data_file_is_2D) then
        call time_interp_external(id_time,time,data(:,:,1),verbose=.false.,horz_interp=id_horz_interp) 
        data(:,:,1) = data(:,:,1)*factor
        do i = 2, size(data,3)
           data(:,:,i) = data(:,:,1)
        enddo
     else
        call time_interp_external(id_time,time,data,verbose=.false.,horz_interp=id_horz_interp)
        data = data*factor
     endif
  endif

  if(present(region1)) then
     do i = comp_domain(1), comp_domain(2)
        do j = comp_domain(3), comp_domain(4)
           if(i>=reg1(1) .and. i<=reg1(2) .and. j>=reg1(3) .and. j<=reg1(4)) &
                data1(i,j,:) = data(i,j,:)
        enddo
     enddo
  else
     data1 = data
  endif

  if(present(region2)) then
     do i = comp_domain(1), comp_domain(2)
        do j = comp_domain(3), comp_domain(4)
           if(i>=reg2(1) .and. i<=reg2(2) .and. j>=reg2(3) .and. j<=reg2(4)) &
                 data1(i,j,:) = data(i,j,:)
        enddo
     enddo
  endif
  if(PRESENT(override)) override = .true.
  deallocate(data)

end subroutine data_override_3d
! </SUBROUTINE>
!...

! Get global lon and lat of three model (target) grids from grid_spec.nc
subroutine get_global_grid()
  character(len=128) :: grid_file = 'INPUT/grid_spec.nc'
  integer :: i, j, siz(4)
  integer :: nlon_out, nlat_out ! size of global lon and lat
  logical :: file_open
  real(r8_kind), dimension(:,:,:), allocatable :: lon_vert_glo, lat_vert_glo !of OCN grid vertices
  real(r8_kind), dimension(:),   allocatable :: lon_atm, lat_atm ! lon and lat of ATM grid
  real(r8_kind), dimension(:),   allocatable :: lon_lnd, lat_lnd ! lon and lat of LND grid
  logical :: is_new_grid


! Test if grid_file is already opened
  inquire (file=trim(grid_file), opened=file_open)
  if(file_open) call mpp_error(FATAL, 'grid_spec.nc already opened')

!1 get global lon and lat of ocean grid vertices

  if (ocn_on .or. ice_on) then
    is_new_grid = .FALSE.
    if(field_exist(grid_file, 'x_T')) then
       is_new_grid = .true.
    else if(field_exist(grid_file, 'geolon_t')) then
       is_new_grid = .FALSE.
    else
       call mpp_error(FATAL,'data_override: both x_T and geolon_t is not in the grid file '//trim(grid_file) )
    endif
 
    if(is_new_grid) then
      call field_size(grid_file, 'x_T', siz)
      nlon_out = siz(1); nlat_out = siz(2)
      allocate(lon_vert_glo(nlon_out,nlat_out,4), lat_vert_glo(nlon_out,nlat_out,4) )
      call read_data(trim(grid_file), 'x_vert_T', lon_vert_glo)
      call read_data(trim(grid_file), 'y_vert_T', lat_vert_glo)
      
!2 Global lon and lat of ocean grid cell centers are determined from adjacent vertices
      if(.not. (allocated(glo_lat_ocn) .or. allocated(glo_lon_ocn))) &
      allocate(glo_lat_ocn(nlon_out, nlat_out), glo_lon_ocn(nlon_out,nlat_out))
      glo_lon_ocn(:,:) = (lon_vert_glo(:,:,1) + lon_vert_glo(:,:,2) + lon_vert_glo(:,:,3) + lon_vert_glo(:,:,4))*0.25
      glo_lat_ocn(:,:) = (lat_vert_glo(:,:,1) + lat_vert_glo(:,:,2) + lat_vert_glo(:,:,3) + lat_vert_glo(:,:,4))*0.25
      deallocate(lon_vert_glo)
      deallocate(lat_vert_glo)
! convert from degree to radian
      glo_lon_ocn = glo_lon_ocn * deg_to_radian
      glo_lat_ocn = glo_lat_ocn * deg_to_radian
    else      
      call field_size(grid_file, 'geolon_vert_t', siz)
      allocate(lon_vert_glo(siz(1),siz(2),1))
      call read_data(trim(grid_file), 'geolon_vert_t', lon_vert_glo)
      
      call field_size(grid_file, 'geolat_vert_t', siz)
      allocate(lat_vert_glo(siz(1),siz(2),1))
      call read_data(trim(grid_file), 'geolat_vert_t', lat_vert_glo)

!2 Global lon and lat of ocean grid cell centers are determined from adjacent vertices
      nlon_out = size(lon_vert_glo,1); nlat_out = size(lon_vert_glo,2)
      if(.not. (allocated(glo_lat_ocn) .or. allocated(glo_lon_ocn))) &
      allocate(glo_lat_ocn(nlon_out - 1, nlat_out - 1), &
           glo_lon_ocn(nlon_out - 1,nlat_out - 1))
      do i = 1, nlon_out - 1
         do j = 1, nlat_out - 1
            glo_lon_ocn(i,j) = (lon_vert_glo(i,j,1) + lon_vert_glo(i+1,j,1))/2.
            glo_lat_ocn(i,j) = (lat_vert_glo(i,j,1) + lat_vert_glo(i,j+1,1))/2.
         enddo
      enddo
      deallocate(lon_vert_glo)
      deallocate(lat_vert_glo)
! convert from degree to radian
      glo_lon_ocn = glo_lon_ocn * deg_to_radian
      glo_lat_ocn = glo_lat_ocn * deg_to_radian
    endif
  endif
!3 Get global lon and lat of ATM grid

  if (atm_on) then
     call field_size(grid_file, 'xta', siz)
     nlon_out = siz(1); allocate(lon_atm(nlon_out))
     call read_data(grid_file, 'xta', lon_atm)

     call field_size(grid_file, 'yta', siz)
     nlat_out = siz(1); allocate(lat_atm(nlat_out))
     call read_data(grid_file, 'yta', lat_atm)

!4 create 2D array of ATM lon and lat
     if(.not.(allocated(glo_lat_atm) .or. allocated(glo_lon_atm))) &
     allocate(glo_lat_atm(nlon_out,nlat_out),glo_lon_atm(nlon_out,nlat_out))
     do i = 1, nlon_out
        do j = 1, nlat_out
           glo_lon_atm(i,j) = lon_atm(i)
           glo_lat_atm(i,j) = lat_atm(j)
        enddo
     enddo
! convert to radian
     glo_lon_atm = glo_lon_atm * deg_to_radian
     glo_lat_atm = glo_lat_atm * deg_to_radian
     deallocate(lon_atm)
     deallocate(lat_atm)
 endif
 
!5 Get global lon and lat of LND grid

 if (lnd_on) then
     call field_size(grid_file, 'xtl', siz)
     nlon_out = siz(1); allocate(lon_lnd(nlon_out))
     call read_data(grid_file, 'xtl', lon_lnd)

     call field_size(grid_file, 'ytl', siz)
     nlat_out = siz(1); allocate(lat_lnd(nlat_out))
     call read_data(grid_file, 'ytl', lat_lnd)

!6 create 2D array of LND lon and lat
     if(.not.(allocated(glo_lat_lnd) .or. allocated(glo_lon_lnd))) &
     allocate(glo_lat_lnd(nlon_out,nlat_out),glo_lon_lnd(nlon_out,nlat_out))
     do i = 1, nlon_out
        do j = 1, nlat_out
           glo_lon_lnd(i,j) = lon_lnd(i)
           glo_lat_lnd(i,j) = lat_lnd(j)
        enddo
     enddo
! convert to radian
     glo_lon_lnd = glo_lon_lnd * deg_to_radian
     glo_lat_lnd = glo_lat_lnd * deg_to_radian
     deallocate(lon_lnd)
     deallocate(lat_lnd)
 endif
 
end subroutine get_global_grid

subroutine get_region_bounds(gridname, region_in, region_out)
! Given gridname and region limits (in lat and lon), this routine returns
! the region's indices (in i and j) determined on global array
! lat values are between (-90, +90), lon values:(0,360)     
! Do not give negative lon
  character(len=3), intent(in) :: gridname ! model grid ID
  real, intent(in) :: region_in(4) !(lat_start, lat_end, lon_start, lon_end)
  integer, intent(out) :: region_out(4) ! istart,iend,jstart,jend  
  real, dimension(:,:), allocatable :: lon_global, lat_global
  integer :: i,j
  real :: lat_start, lat_end, lon_start, lon_end
  integer :: size_lat, size_lon

! Get lon_global and lat_global
  select case(gridname)
  case('OCN')    
     allocate(lon_global(size(glo_lon_ocn,1),size(glo_lon_ocn,2)))
     allocate(lat_global(size(glo_lat_ocn,1),size(glo_lat_ocn,2)))
     lon_global = glo_lon_ocn; lat_global = glo_lat_ocn
  case('ICE')          
     allocate(lon_global(size(glo_lon_ocn,1),size(glo_lon_ocn,2)))
     allocate(lat_global(size(glo_lat_ocn,1),size(glo_lat_ocn,2)))
     lon_global = glo_lon_ocn; lat_global = glo_lat_ocn
  case('ATM')       
     allocate(lon_global(size(glo_lon_atm,1),size(glo_lon_atm,2)))
     allocate(lat_global(size(glo_lat_atm,1),size(glo_lat_atm,2)))
     lon_global = glo_lon_atm; lat_global = glo_lat_atm
  case('LND')          
     allocate(lon_global(size(glo_lon_lnd,1),size(glo_lon_lnd,2)))
     allocate(lat_global(size(glo_lat_lnd,1),size(glo_lat_lnd,2)))
     lon_global = glo_lon_lnd; lat_global = glo_lat_lnd
  case default
     call mpp_error(FATAL,'error: grid not recognized in data_override')
  end select

! convert lon_global and lat_global from radian to degree
  lon_global = lon_global * radian_to_deg
  lat_global = lat_global * radian_to_deg
  size_lat = size(lat_global,2); size_lon = size(lon_global,1)

! get indices in x and y axes
  lat_start = region_in(1)
  lat_end = region_in(2)
  lon_start = region_in(3)
  lon_end = region_in(4)

  if(lat_start < lat_global(1,1) .or. lat_end > lat_global(1,size_lat) .or.   &
       lon_start < lon_global(1,1) .or. lon_end > lon_global(size_lon,1))     &
     call mpp_error(FATAL,'data_override: error in region bounds')
  
  region_out(1) = get_index(lon_start,lon_global(:,1))
  region_out(2) = get_index(lon_end,lon_global(:,1))
  region_out(3) = get_index(lat_start,lat_global(1,:))
  region_out(4) = get_index(lat_end,lat_global(1,:))
  deallocate(lon_global,lat_global)
  if(region_out(1)==-1 .or. region_out(2)==-1 .or. region_out(3)==-1 .or. region_out(4) ==-1) &
        call mpp_error(FATAL,'ERROR:data_override: get_region_bounds ')
end subroutine get_region_bounds

function get_index(number, array)
  real :: number
  real, dimension(:) :: array
  integer get_index, i, n
!Find index i of array such that array(i) is closest to number
! array must be monotonically increasing

  n = size(array(:))
  do i = 2, n
     if(array(i) < array(i-1)) call mpp_error(FATAL,'ERROR:data_override: get_index') 
  enddo
  get_index = -1

     do i = 1, n-1         
           if ((array(i)<=number) .and. (array(i+1)>= number)) then
              if(number - array(i) <= array(i+1) - number) then
                 get_index = i                 
              else
                 get_index = i+1
              endif
              exit
           endif
     enddo
  end function get_index

end module data_override_mod
 
