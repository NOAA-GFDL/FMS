module station_data_mod 
! <CONTACT EMAIL="Giang.Nong@gfdl.noaa.gov">
!   Giang Nong
! </CONTACT>
! <OVERVIEW>
! This module is used for outputing model results in a list
! of stations (not gridded arrrays). The user needs to supply
! a list of stations with lat, lon values of each station.
! Data at a single point (I,J) that is closest to each station will
! be written to a file. No interpolation is made when a station
! is between two or more grid points.<BR/>
! In the output file, a 3D field will have a format of array(n1,n2) and
! a 2D field is array(n1) where n1 is number of stations and n2 is number
! of vertical levels or depths.
! </OVERVIEW>
! <DESCRIPTION>
! Here are some basic steps of how to use station_data_mod <BR/>
!1/Call <TT>data_station_init</TT>  <BR/>
! user needs to supply 2 tables: list_stations and station_data_table as follows:<BR/>
! example of  list of stations (# sign means comment)<BR/>
!               #  station_id          lat    lon <BR/>
!                   station_1          20.4   100.8 <BR/>
! example of station_data_table (# sign means comment) <BR/>
! # General descriptor <BR/>
! Am2p14 station data <BR/>

!#  start time (should be the same as model's initial time) <BR/>
! 19800101 <BR/>
!# file inforamtion <BR/>
!#   filename,    output_frequency, frequency_unit, time_axis_unit <BR/>
!   "ocean_day"         1              "days"           "hours" <BR/>
!# field information <BR/>
!# module     field_name    filename    time_method   pack <BR/>
!  Ice_mod    temperature   ocean_day     . true.       2 <BR/>
!  Ice_mod    pressure      ocean_day      .false.      2    <BR/>
! 2/
! Call register_station_field to register each field that needs to be written to a file, the call
! <TT>register_station_field</TT> returns a field_id that will be used later in send_station_data <BR/>
! 3/
! Call <TT> send_station_data</TT> will send data at each station in the list
! to a file <BR/>
! 4/ Finally, call <TT>station_data_end</TT> after the last time step.<BR/>
! </DESCRIPTION>
use axis_utils_mod, only: nearest_index
use mpp_io_mod,    only : mpp_open, MPP_RDONLY, MPP_ASCII, mpp_close,MPP_OVERWR,MPP_NETCDF, &
                          mpp_write_meta, MPP_SINGLE, mpp_write, fieldtype,mpp_flush
use fms_mod,       only : error_mesg, FATAL, WARNING, stdlog, write_version_number,&
                          mpp_pe, lowercase, stdout, close_file, open_namelist_file, check_nml_error
use mpp_mod,       only : mpp_npes,  mpp_sync, mpp_root_pe, mpp_send, mpp_recv, mpp_max, &
                          mpp_get_current_pelist, input_nml_file, &
                          COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4
use mpp_domains_mod,only: domain2d, mpp_get_compute_domain
use diag_axis_mod, only : diag_axis_init
use diag_output_mod,only:  write_axis_meta_data, write_field_meta_data,diag_fieldtype,done_meta_data
use diag_manager_mod,only : get_date_dif, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, &
                            DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS
use diag_util_mod,only    : diag_time_inc
use time_manager_mod, only: operator(>),operator(>=),time_type,get_calendar_type,NO_CALENDAR,set_time, &
                            set_date, increment_date, increment_time
implicit none
private
integer, parameter  :: max_fields_per_file = 150
integer, parameter  :: max_files = 31
integer             :: num_files = 0
integer             :: num_stations = 0
integer             :: max_stations = 20
integer             :: max_output_fields = 100
integer             :: num_output_fields = 0
real                :: EMPTY = 0.0
real                :: MISSING = 1.E20
logical             :: module_is_initialized = .false.
logical             :: need_write_axis = .true.
integer             :: base_year, base_month, base_day, base_hour, base_minute, base_second
type (time_type)    :: base_time
character (len=10)  :: time_unit_list(6) = (/'seconds   ', 'minutes   ', &
     'hours     ', 'days      ', 'months    ', 'years     '/)
integer, parameter  :: EVERY_TIME =  0
integer, parameter  :: END_OF_RUN = -1
! Include variable "version" to be written to log file.
#include<file_version.h>
character(len=256)  :: global_descriptor
character (len = 7) :: avg_name = 'average'
integer             :: total_pe
integer             :: lat_axis, lon_axis
integer, allocatable:: pelist(:)
character(len=32)   :: pelist_name
type file_type
   character(len=128)  :: name
   integer             :: output_freq
   integer             :: output_units
   integer             :: time_units
   integer             :: fields(max_fields_per_file)
   integer             :: num_fields
   integer             :: file_unit
   integer             :: time_axis_id, time_bounds_id
   type (time_type)    :: last_flush
   type(fieldtype)     :: f_avg_start, f_avg_end, f_avg_nitems, f_bounds
end type file_type
type station_type
   character(len=128)  :: name
   real                :: lat, lon
   integer             :: id
   integer             :: global_i, global_j  ! index of global grid
   integer             :: local_i, local_j    ! index on the current PE
   logical             :: need_compute        ! true if the station present in this PE
end type station_type

type group_field_type
   integer             :: output_file
   integer             :: num_station ! number of stations on this PE
   integer, pointer    :: station_id(:)  =>null() ! id of station on this PE
   character(len=128)  :: output_name, module_name,long_name, units
   logical             :: time_average,time_max,time_min, time_ops, register
   integer             :: pack, axes(2), num_axes
   character(len=8)    :: time_method   !could be: true, false, max, min, mean, ...
   real, pointer       :: buffer(:, :)=>null()
   integer             :: counter,nlevel
   type(time_type)     :: last_output, next_output
   type(fieldtype)     :: f_type
end type group_field_type

type global_field_type
   real, pointer       :: buffer(:,:)=>null()
   integer             :: counter
end type global_field_type   


type(global_field_type),save            :: global_field
type (file_type),save                   :: files(max_files)
type(group_field_type),allocatable,save :: output_fields(:)
type (station_type),allocatable         :: stations(:)
type(diag_fieldtype),save               :: diag_field
public register_station_field, send_station_data, station_data_init, station_data_end

interface register_station_field
    module procedure register_station_field2d
    module procedure register_station_field3d
end interface
interface send_station_data
    module procedure send_station_data_2d
    module procedure send_station_data_3d
end interface
contains

! <INTERFACE NAME="station_data_init">
! <TEMPLATE>
! station_data_init()
! </TEMPLATE>
!   <DESCRIPTION>
! read in lat. lon of each station<BR/>
! create station_id based on lat, lon<BR/>
! read station_data_table, initialize output_fields and output files<BR/>
!   </DESCRIPTION>
! </INTERFACE>
subroutine station_data_init()

character(len=128)    :: station_name
real                  :: lat, lon  
integer               :: iunit,nfiles,nfields,time_units,output_freq_units,j,station_id,io_status,logunit, ierr
logical               :: init_verbose
character(len=128)    :: record
type file_part_type
   character(len=128) :: name
   integer            :: output_freq
   character(len=10)  :: output_freq_units 
   integer            :: format    ! must always be 1 for netcdf files
   character(len=10)  :: time_unit
end type file_part_type
type field_part_type
   character(len=128) :: module_name,field_name,file_name
   character(len=8)   :: time_method   
   integer            :: pack
end type field_part_type

type(file_part_type)  :: record_files
type(field_part_type) :: record_fields

namelist /station_data_nml/ max_output_fields, max_stations,init_verbose

  if (module_is_initialized) return
  init_verbose = .false.
  total_pe = mpp_npes()
  allocate(pelist(total_pe))
  call mpp_get_current_pelist(pelist, pelist_name) 

! read namelist
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, station_data_nml, iostat=io_status)
  ierr = check_nml_error(io_status, 'station_data_nml')
#else
  iunit = open_namelist_file ()
  ierr=1; do while (ierr /= 0)
  read  (iunit, nml=station_data_nml, iostat=io_status, end=10)
  ierr = check_nml_error(io_status, 'station_data_nml')
  enddo
10 call close_file (iunit)

#endif
  logunit = stdlog()
  write(logunit, station_data_nml)

  allocate(output_fields(max_output_fields), stations(max_stations))
! read list of stations
  if(init_verbose) then
     logunit = stdout()
     write(logunit, *) ' '
     write(logunit, *) '****** Summary of STATION information from list_stations ********'
     write(logunit, *) ' '
     write(logunit, *) 'station name      ', '   latitude ', '   longitude '
     write(logunit, *) ' '
  endif
  call mpp_open(iunit, 'list_stations',form=MPP_ASCII,action=MPP_RDONLY)
  do while (num_stations<max_stations)
     read(iunit,'(a)',end=76,err=75) record
     if (record(1:1) == '#') cycle 
     if(len_trim(record) < 1) cycle
     read(record, *, end = 76, err = 75) station_name, lat, lon 
     station_id = get_station_id(lat, lon)
     if(station_id > 0) then       
        stations(station_id)%name = station_name
     else
        call error_mesg('station_data_init','station DUPLICATED in file list_stations', FATAL)
     endif
     logunit = stdout()
     if( init_verbose.and.  mpp_pe() == mpp_root_pe()) &
          write(logunit,1)stations(station_id)%name,stations(station_id)%lat,stations(station_id)%lon
1 format(1x,A18, 1x,F8.2,4x,F8.2)     
75   continue
  enddo
  call error_mesg('station_data_init','max_stations exceeded, increase it via namelist', FATAL)
76 continue
  call mpp_close (iunit)
  logunit = stdout()
  if(init_verbose)  write(logunit, *)'*****************************************************************'
     
! read station_data table
  call mpp_open(iunit, 'station_data_table',form=MPP_ASCII,action=MPP_RDONLY)
! Read in the global file labeling string
  read(iunit, *, end = 99, err=99) global_descriptor

! Read in the base date
  read(iunit, *, end = 99, err = 99) base_year, base_month, base_day, &
       base_hour, base_minute, base_second
  if (get_calendar_type() /= NO_CALENDAR) then
     base_time = set_date(base_year, base_month, base_day, base_hour, &
          base_minute, base_second)
  else
! No calendar - ignore year and month
     base_time = set_time(base_hour*3600+base_minute*60+base_second, base_day)
     base_year  = 0
     base_month = 0
  end if
  nfiles=0
  do while (nfiles <= max_files)
     read(iunit,'(a)',end=86,err=85) record
     if (record(1:1) == '#') cycle        
     read(record,*,err=85,end=85)record_files%name,record_files%output_freq, &
          record_files%output_freq_units,record_files%format,record_files%time_unit
     if(record_files%format /= 1) cycle   !avoid reading field part
     time_units = 0
     output_freq_units = 0
     do j = 1, size(time_unit_list(:))
        if(record_files%time_unit == time_unit_list(j)) time_units = j
        if(record_files%output_freq_units == time_unit_list(j)) output_freq_units = j     
     end do
     if(time_units == 0) &
          call error_mesg('station_data_init',' check time unit in station_data_table',FATAL)
     if(output_freq_units == 0) & 
          call error_mesg('station_data_init',', check output_freq in station_data_table',FATAL)
      call init_file(record_files%name,record_files%output_freq, output_freq_units,time_units)
85    continue
   enddo
   call error_mesg('station_data_init','max_files exceeded, increase max_files', FATAL)
86 continue
   rewind(iunit)
   nfields=0
   do while (nfields <= max_output_fields)
       read(iunit,'(a)',end=94,err=93) record
       if (record(1:1) == '#') cycle
       read(record,*,end=93,err=93) record_fields
       if (record_fields%pack .gt. 8 .or.record_fields%pack .lt. 1) cycle !avoid reading file part
       nfields=nfields+1
       call init_output_field(record_fields%module_name,record_fields%field_name, &
            record_fields%file_name,record_fields%time_method,record_fields%pack)
93     continue
    enddo
    call error_mesg('station_data_init','max_output_fields exceeded, increase it via nml ', FATAL)
94  continue
    call close_file(iunit)
    call check_duplicate_output_fields
    call write_version_number ("STATION_DATA_MOD", version)
    module_is_initialized = .true.
    return
99  continue
    call error_mesg('station_data_init','error reading station_datatable',FATAL)
end subroutine station_data_init
!----------------------------------------------------------------------
subroutine check_duplicate_output_fields()
! pair(output_name and output_file) should be unique in data_station_table, ERROR1
! pair(module_name and output_name) should be unique in data_station_table, ERROR2
integer            :: i, j, tmp_file
character(len=128) :: tmp_name, tmp_module

if(num_output_fields <= 1) return 
do i = 1, num_output_fields-1
   tmp_name = trim(output_fields(i)%output_name)
   tmp_file =  output_fields(i)%output_file
   tmp_module = trim(output_fields(i)%module_name)
   do j = i+1, num_output_fields
      if((tmp_name == trim(output_fields(j)%output_name)).and. &
           (tmp_file == output_fields(j)%output_file)) &
           call error_mesg (' ERROR1 in station_data_table:', &           
           &' module/field '//tmp_module//'/'//tmp_name//' duplicated', FATAL)
      if((tmp_name == trim(output_fields(j)%output_name)).and. &
           (tmp_module == trim(output_fields(j)%module_name))) &
           call error_mesg (' ERROR2 in station_data_table:', &           
           &' module/field '//tmp_module//'/'//tmp_name//' duplicated', FATAL)
   enddo
enddo
end subroutine check_duplicate_output_fields
!----------------------------------------------------------------------
function get_station_id(lat,lon)
  integer         :: get_station_id, i
  real, intent(in):: lat,lon
! each station should have distinct lat and lon
  get_station_id = -1
  do i = 1, num_stations
     if(stations(i)%lat == lat .and. stations(i)%lon == lon) return
  enddo
  num_stations = num_stations + 1
  stations(num_stations)%id = num_stations
  stations(num_stations)%lat = lat
  stations(num_stations)%lon = lon
  stations(num_stations)%need_compute = .false.
  stations(num_stations)%global_i = -1; stations(num_stations)%global_j = -1 
  stations(num_stations)%local_i = -1 ; stations(num_stations)%local_j = -1
  get_station_id = num_stations
end function get_station_id
!----------------------------------------------------------------------
subroutine init_file(filename, output_freq, output_units, time_units)
  character(len=*), intent(in) :: filename
  integer, intent(in)          :: output_freq, output_units, time_units
  character(len=128)           :: time_units_str
  real, dimension(1)           :: tdata

  num_files = num_files + 1
  if(num_files >= max_files) &
       call error_mesg('station_data, init_file', ' max_files exceeded, incease max_files', FATAL)
  files(num_files)%name = trim(filename)
  files(num_files)%output_freq = output_freq
  files(num_files)%output_units = output_units
  files(num_files)%time_units = time_units
  files(num_files)%num_fields = 0
  files(num_files)%last_flush = base_time
  files(num_files)%file_unit = -1
!---- register axis_id and time boundaries id
  write(time_units_str, 11) trim(time_unit_list(files(num_files)%time_units)), base_year, &
       base_month, base_day, base_hour, base_minute, base_second
11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
  files(num_files)%time_axis_id = diag_axis_init ('Time', tdata, time_units_str, 'T',  &
       'Time' , set_name=trim(filename))
  files(num_files)%time_bounds_id = diag_axis_init('nv',(/1.,2./),'none','N','vertex number',&
       set_name=trim(filename))
end subroutine init_file

!--------------------------------------------------------------------------
subroutine init_output_field(module_name,field_name,file_name,time_method,pack)
  character(len=*), intent(in)           :: module_name, field_name, file_name
  character(len=*), intent(in)           :: time_method
  integer, intent(in)                    :: pack
  integer                                :: out_num, file_num,num_fields, method_selected, l1 
  character(len=8)                       :: t_method
! Get a number for this output field
  num_output_fields = num_output_fields + 1
  if(num_output_fields > max_output_fields) &
       call error_mesg('station_data', 'max_output_fields exceeded, increase it via nml', FATAL)
  out_num = num_output_fields
  file_num = find_file(file_name)
   if(file_num < 0) &
        call error_mesg('station_data,init_output_field', 'file '//trim(file_name) &
        //' is NOT found in station_data_table', FATAL)
! Insert this field into list of fields of this file
   files(file_num)%num_fields = files(file_num)%num_fields + 1
   if(files(file_num)%num_fields > max_fields_per_file) &
        call error_mesg('station_data, init_output_field', 'max_fields_per_file exceeded ', FATAL)
   num_fields = files(file_num)%num_fields
   files(file_num)%fields(num_fields)  = out_num
   output_fields(out_num)%output_name  = trim(field_name)
   output_fields(out_num)%module_name  = trim(module_name)
   output_fields(out_num)%counter      = 0
   output_fields(out_num)%output_file  = file_num
   output_fields(out_num)%pack         = pack
   output_fields(out_num)%time_average = .false.
   output_fields(out_num)%time_min     = .false.
   output_fields(out_num)%time_max     = .false. 
   output_fields(out_num)%time_ops     = .false.
   output_fields(out_num)%register     = .false.
   t_method = lowercase(time_method)
 select case (trim(t_method))
 case('.true.')
    output_fields(out_num)%time_average = .true.
    output_fields(out_num)%time_method  = 'mean'
 case('mean')
    output_fields(out_num)%time_average = .true.
    output_fields(out_num)%time_method  = 'mean'
 case('average')
    output_fields(out_num)%time_average = .true.
    output_fields(out_num)%time_method  = 'mean'
 case('avg')
    output_fields(out_num)%time_average = .true.
    output_fields(out_num)%time_method  = 'mean'
 case('.false.')
    output_fields(out_num)%time_average = .false.
     output_fields(out_num)%time_method  = 'point'
 case ('max')
    call error_mesg('station_data, init_output_field','time_method MAX is not supported',&
         FATAL)
    output_fields(out_num)%time_max = .true.
    output_fields(out_num)%time_method  = 'max'
    l1 = len_trim(output_fields(out_num)%output_name)
    if(output_fields(out_num)%output_name(l1-2:l1) /= 'max') &
           output_fields(out_num)%output_name = trim(field_name)//'_max'      
 case ('min')
    call error_mesg('station_data, init_output_field','time_method MIN is not supported',&
         FATAL)
    output_fields(out_num)%time_min = .true.
    output_fields(out_num)%time_method  = 'min'
    l1 = len_trim(output_fields(out_num)%output_name)
    if(output_fields(out_num)%output_name(l1-2:l1) /= 'min') &
         output_fields(out_num)%output_name = trim(field_name)//'_min'
 case default
    call error_mesg('station_data, init_output_field', 'error in time_method of field '&
         //trim(field_name), FATAL)
 end select
 if (files(file_num)%output_freq == EVERY_TIME) &
      output_fields(out_num)%time_average = .false.
 output_fields(out_num)%time_ops = output_fields(out_num)%time_min.or.output_fields(out_num)%time_max &
      .or.output_fields(out_num)%time_average
 output_fields(out_num)%time_method = trim(time_method)
end subroutine init_output_field
!--------------------------------------------------------------------------
function find_file(name)
integer                      :: find_file
character(len=*), intent(in) :: name
integer                      :: i

find_file = -1
do i = 1, num_files
   if(trim(files(i)%name) == trim(name)) then
      find_file = i
      return
   end if
end do
end function find_file
! <INTERFACE NAME="register_station_field">
! <TEMPLATE>
! register_station_field (module_name,fieldname,glo_lat,glo_lon,levels,init_time, 
!     domain,longname,units) <BR/>
! </TEMPLATE>
!   <DESCRIPTION>
! This function is similar to register_diag_field of diag_manager_mod. All arguments
! are inputs that user needs to supply, some are optional. The names of input args are
! self-describing.<BR/> levels is absent for 2D fields. <BR/>
! Note that pair (module_name, fieldname) must be unique in the 
! station_data_table or a fatal error will occur. <BR/>
! A field id is returned from this call that will be used later in send_station_data. <BR/>
!   </DESCRIPTION>
! </INTERFACE>


!--------------------------------------------------------------------------
function register_station_field2d (module_name,fieldname,glo_lat,glo_lon,init_time, &
     domain,longname,units)
  integer                                :: register_station_field2d
  character(len=*), intent(in)           :: module_name, fieldname
  real,dimension(:), intent(in)          :: glo_lat,glo_lon
  type(domain2d), intent(in)             :: domain
  type(time_type), intent(in)            :: init_time
  character(len=*), optional, intent(in) :: longname, units
  real                                   :: levels(1:1)

  levels = 0.
  register_station_field2d = register_station_field3d (module_name,fieldname,glo_lat,glo_lon,&
       levels,init_time,domain,longname,units)
end function register_station_field2d
!--------------------------------------------------------------------------

function register_station_field3d (module_name,fieldname,glo_lat,glo_lon,levels,init_time, &
     domain,longname,units)

! write field meta data on ROOT PE only
! allocate buffer
  integer                                :: register_station_field3d
  character(len=*), intent(in)           :: module_name, fieldname
  real,dimension(:), intent(in)          :: glo_lat,glo_lon,levels !in X,Y,Z direction respectively
  type(domain2d), intent(in)             :: domain
  type(time_type), intent(in)            :: init_time
  character(len=*), optional, intent(in) :: longname,units
  integer                                :: i,ii, nlat, nlon,nlevel, isc, iec, jsc, jec
  character(len=128)                     :: error_msg
  integer                                :: local_num_stations ! number of stations on this PE
  integer                                :: out_num ! position of this field in array output_fields
  integer                                :: file_num, freq, output_units, outunit
  real, allocatable                      :: station_values(:), level_values(:)
  character(len=128)                     :: longname2,units2


  if(PRESENT(longname)) then
     longname2 = longname
  else
     longname2 = fieldname
  endif
  if(PRESENT(units)) then
     units2 = units
  else
     units2 = "none"
  endif

  nlat = size(glo_lat); nlon = size(glo_lon); nlevel=size(levels)
  allocate(station_values(num_stations), level_values(nlevel))
  do i = 1, nlevel
     level_values(i) = real(i)
  enddo
! determine global index of this field in all stations
  outunit = stdout()
  do i = 1,num_stations
     station_values(i) = real(i)
     if(stations(i)%lat<glo_lat(1) .or. stations(i)%lat>glo_lat(nlat)) then
        write(error_msg,'(F9.3)') stations(i)%lat
        write(outunit,*) 'Station with latitude '//trim(error_msg)//' outside global latitude values'
        call error_mesg ('register_station_field', 'latitude out of range', FATAL)
     endif
      if(stations(i)%lon<glo_lon(1) .or. stations(i)%lon>glo_lon(nlon)) then
        write(error_msg,'(F9.3)') stations(i)%lon
        write(outunit,*) 'Station with longitude '//trim(error_msg)//' outside global longitude values'
        call error_mesg ('register_station_field', 'longitude out of range', FATAL)
     endif
     stations(i)%global_i = nearest_index(stations(i)%lon, glo_lon)
     stations(i)%global_j = nearest_index(stations(i)%lat, glo_lat)
     if(stations(i)%global_i<0 .or. stations(i)%global_j<0) &
          call error_mesg ('register_station_field', 'Error in global index of station',FATAL)
  enddo
! determine local index of this field in all stations , local index starts from 1
  call mpp_get_compute_domain(domain, isc,iec,jsc,jec)
  local_num_stations = 0
  do i = 1,num_stations
     if(isc<=stations(i)%global_i .and. iec>= stations(i)%global_i .and. &
        jsc<=stations(i)%global_j .and. jec>= stations(i)%global_j) then
        stations(i)%need_compute = .true.
        stations(i)%local_i = stations(i)%global_i - isc + 1
        stations(i)%local_j = stations(i)%global_j - jsc + 1
        local_num_stations = local_num_stations +1    
     endif
  enddo
! get the position of this field in the array output_fields
  out_num = find_output_field(module_name, fieldname)  
  if(out_num < 0 .and. mpp_pe() == mpp_root_pe()) then 
     call error_mesg ('register_station_field', &
          'module/field_name '//trim(module_name)//'/'//&
          trim(fieldname)//' NOT found in station_data table', WARNING)
     register_station_field3d = out_num
     return
  endif
  if(local_num_stations>0) then
     allocate(output_fields(out_num)%station_id(local_num_stations))
     allocate(output_fields(out_num)%buffer(local_num_stations,nlevel))
     output_fields(out_num)%buffer = EMPTY
! fill out list of available stations in this PE 
     ii=0
     do i = 1,num_stations       
        if(stations(i)%need_compute) then
           ii = ii+ 1
           if(ii>local_num_stations) call error_mesg ('register_station_field', &
                'error in determining local_num_station', FATAL)
           output_fields(out_num)%station_id(ii)=stations(i)%id
        endif
     enddo
  endif
  output_fields(out_num)%num_station = local_num_stations
  if( mpp_pe() == mpp_root_pe()) then
     allocate(global_field%buffer(num_stations,nlevel))
     global_field%buffer = MISSING
  endif
  output_fields(out_num)%register = .true.
  output_fields(out_num)%output_name = fieldname
  file_num = output_fields(out_num)%output_file
  output_fields(out_num)%last_output = init_time
  freq = files(file_num)%output_freq
  output_units = files(file_num)%output_units
  output_fields(out_num)%next_output = diag_time_inc(init_time, freq, output_units)
  register_station_field3d = out_num
  output_fields(out_num)%long_name = longname2
  output_fields(out_num)%units = units2
  output_fields(out_num)%nlevel = nlevel
! deal with axes
 
  output_fields(out_num)%axes(1) = diag_axis_init('Stations',station_values,'station number', 'X')
  if(nlevel == 1) then
     output_fields(out_num)%num_axes = 1     
  else
     output_fields(out_num)%num_axes = 2     
     output_fields(out_num)%axes(2) = diag_axis_init('Levels',level_values,'level number', 'Y' )   
  endif
  if(need_write_axis) then
     lat_axis = diag_axis_init('Latitude', stations(1:num_stations)%lat,'station latitudes', 'n')
     lon_axis = diag_axis_init('Longitude',stations(1:num_stations)%lon,  'station longitudes', 'n')
  endif
  need_write_axis = .false.
 
!  call mpp_sync()

end function register_station_field3d

!-------------------------------------------------------------------------

function find_output_field(module_name, field_name)
  integer find_output_field
  character(len=*), intent(in) :: module_name, field_name
  integer                      :: i

  find_output_field = -1
  do i = 1, num_output_fields
     if(trim(output_fields(i)%module_name) == trim(module_name) .and. &
          lowercase(trim(output_fields(i)%output_name)) == &
          lowercase(trim(field_name))) then 
        find_output_field = i
        return
     endif
  end do
end function find_output_field

!-------------------------------------------------------------------------
subroutine opening_file(file)
! open file, write axis meta_data for all files (only on ROOT PE, 
!                        do nothing on other PEs)
 integer, intent(in)  :: file
 character(len=128)   :: time_units
 integer              :: j,field_num,num_axes,axes(5),k
 logical              :: time_ops
 integer              :: time_axis_id(1),time_bounds_id(1)

 write(time_units, 11) trim(time_unit_list(files(file)%time_units)), base_year, &
      base_month, base_day, base_hour, base_minute, base_second
11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
 call mpp_open(files(file)%file_unit, files(file)%name, action=MPP_OVERWR, &
       form=MPP_NETCDF, threading=MPP_SINGLE, fileset=MPP_SINGLE)
 call mpp_write_meta (files(file)%file_unit, 'title', cval=trim(global_descriptor))
 time_ops = .false.
 do j = 1, files(file)%num_fields
    field_num = files(file)%fields(j)
    if(output_fields(field_num)%time_ops) then
       time_ops = .true.
       exit
    endif
 enddo
!write axis meta data
 do j = 1, files(file)%num_fields
    field_num = files(file)%fields(j)
    num_axes = output_fields(field_num)%num_axes
    axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
    do k = 1,num_axes
       if(axes(k)<0) &
            call error_mesg ('station_data opening_file','output_name '// &
            trim(output_fields(field_num)%output_name)// &
            ' has axis_id = -1', FATAL)
    enddo
    axes(num_axes + 1) = lat_axis
    axes(num_axes + 2) = lon_axis
    axes(num_axes + 3) = files(file)%time_axis_id
! need in write_axis: name, unit,long_name
    call write_axis_meta_data(files(file)%file_unit, axes(1:num_axes + 3), time_ops)
    if(time_ops) then
       axes(num_axes + 4) = files(file)%time_bounds_id
       call write_axis_meta_data(files(file)%file_unit, axes(1:num_axes + 4))     
    endif
 end do
!write field meta data
 do j = 1, files(file)%num_fields
    field_num = files(file)%fields(j)
    num_axes = output_fields(field_num)%num_axes
    axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
    num_axes = num_axes + 1
    axes(num_axes) = files(file)%time_axis_id    
    diag_field = write_field_meta_data(files(file)%file_unit,output_fields(field_num)%output_name, &
         axes(1:num_axes),output_fields(field_num)%units, &
         output_fields(field_num)%long_name,time_method=output_fields(field_num)%time_method,&
         pack=output_fields(field_num)%pack)
    output_fields(field_num)%f_type = diag_field%Field
 end do
 if(time_ops) then
    time_axis_id(1) = files(file)%time_axis_id
    time_bounds_id(1) = files(file)%time_bounds_id
    diag_field=write_field_meta_data(files(file)%file_unit, avg_name // '_T1',time_axis_id, &
        time_units,"Start time for average period", pack=1)
    files(file)%f_avg_start = diag_field%Field

    diag_field=write_field_meta_data(files(file)%file_unit,avg_name // '_T2' ,time_axis_id, &
          time_units,"End time for average period", pack=1)
    files(file)%f_avg_end = diag_field%Field

    diag_field=write_field_meta_data(files(file)%file_unit,avg_name // '_DT' ,time_axis_id, &
          time_units,"Length of average period", pack=1)  
    files(file)%f_avg_nitems = diag_field%Field

    diag_field=write_field_meta_data(files(file)%file_unit, 'Time_bnds', (/time_bounds_id,time_axis_id/), &
           trim(time_unit_list(files(file)%time_units)), &
           'Time axis boundaries', pack=1) 
     files(file)%f_bounds =  diag_field%Field
 endif
 call done_meta_data(files(file)%file_unit)
end subroutine opening_file 
! <INTERFACE NAME="send_station_data">
! <TEMPLATE>
! send_station_data(field_id, data, time)
! </TEMPLATE>
!   <DESCRIPTION>
! data should have the size of compute domain(isc:iec,jsc:jec)<BR/>
! time is model's time<BR/>
! field_id is returned from <TT>register_station_field</TT><BR/>
! only data at stations will be be sent to root_pe which, in turn, sends to output file
!   </DESCRIPTION>
! </INTERFACE>
!-------------------------------------------------------------------------
subroutine send_station_data_2d(field_id, data, time)
  integer, intent(in)         :: field_id
  real,    intent(in)         :: data(:,:)
  type(time_type), intent(in) :: time
  real                        :: data3d(size(data,1),size(data,2),1)

  data3d(:,:,1) = data
  call send_station_data_3d(field_id, data3d, time)
end subroutine send_station_data_2d
!-------------------------------------------------------------------------
subroutine send_station_data_3d(field_id, data, time)
 
  integer, intent(in)         :: field_id
  real,    intent(in)         :: data(:,:,:)
  type(time_type), intent(in) :: time
  integer                     :: freq,units,file_num,local_num_stations,i,ii, max_counter
  integer                     :: index_x, index_y, station_id
  integer, allocatable        :: station_ids(:)  ! root_pe only,  to receive local station_ids
  real,    allocatable        :: tmp_buffer(:,:) ! root_pe only, to receive local buffer from each PE
  

  if (.not.module_is_initialized) &
     call error_mesg ('send_station_data_3d',' station_data NOT initialized', FATAL)

  if(field_id < 0) return  
  file_num = output_fields(field_id)%output_file
  if( mpp_pe() == mpp_root_pe() .and. files(file_num)%file_unit < 0) then
     call opening_file(file_num)
  endif
  freq = files(file_num)%output_freq
  units = files(file_num)%output_units
! compare time with next_output

  if (time > output_fields(field_id)%next_output .and. freq /= END_OF_RUN) then  ! time to write out     
! ALL PEs, including root PE, must send data to root PE        
     call mpp_send(output_fields(field_id)%num_station,plen=1,to_pe=mpp_root_pe(),tag=COMM_TAG_1)
     if(output_fields(field_id)%num_station > 0) then
        call mpp_send(output_fields(field_id)%station_id(1),plen=size(output_fields(field_id)%station_id),&
             to_pe=mpp_root_pe())
        call mpp_send(output_fields(field_id)%buffer(1,1),plen=size(output_fields(field_id)%buffer),&
             to_pe=mpp_root_pe())
     endif   
! get max_counter if the field is averaged
     if(output_fields(field_id)%time_average) then
        max_counter = output_fields(field_id)%counter
        call mpp_max(max_counter, pelist)
     endif
! receive local data from all PEs 
     if(mpp_pe() == mpp_root_pe()) then
        do i = 1,size(pelist)           
           call mpp_recv(local_num_stations,glen=1,from_pe=pelist(i),tag=COMM_TAG_1)
           if(local_num_stations> 0) then
              allocate(station_ids(local_num_stations))
              allocate(tmp_buffer(local_num_stations,output_fields(field_id)%nlevel))
              call mpp_recv(station_ids(1), glen=size(station_ids), from_pe=pelist(i))
              call mpp_recv(tmp_buffer(1,1),glen=size(tmp_buffer),  from_pe=pelist(i)) 
              do ii = 1,local_num_stations
                 global_field%buffer(station_ids(ii),:) = tmp_buffer(ii,:)
              enddo
              deallocate(station_ids, tmp_buffer)
           endif
        enddo
! send global_buffer content to file
        if(output_fields(field_id)%time_average) then  
           if(max_counter == 0 ) &
                call error_mesg ('send_station_data','counter=0 for averaged field '// &
                output_fields(field_id)%output_name, FATAL)
           global_field%buffer = global_field%buffer/real(max_counter)
        endif
! check if global_field contains any missing values
        if(any(global_field%buffer == MISSING)) &
             call error_mesg ('send_station_data','Global_field contains MISSING, field '// &
             output_fields(field_id)%output_name, FATAL)
        call station_data_out(file_num,field_id,global_field%buffer,output_fields(field_id)%next_output)
        global_field%buffer = MISSING
     endif
     call mpp_sync()
! clear buffer, increment next_output time and reset counter on ALL PEs
     if(output_fields(field_id)%num_station>0)  output_fields(field_id)%buffer = EMPTY
     output_fields(field_id)%last_output = output_fields(field_id)%next_output
     output_fields(field_id)%next_output =  diag_time_inc(output_fields(field_id)%next_output,&
          freq, units)
     output_fields(field_id)%counter = 0; max_counter = 0
    
  endif
! accumulate buffer only
  do i = 1 , output_fields(field_id)%num_station
     station_id = output_fields(field_id)%station_id(i)
     index_x = stations(station_id)%local_i; index_y = stations(station_id)%local_j
     if(index_x>size(data,1) .or. index_y>size(data,2)) &
          call error_mesg ('send_station_data','local index out of range for field '// &
          output_fields(field_id)%output_name, FATAL) 
     if(output_fields(field_id)%time_average) then
        output_fields(field_id)%buffer(i,:) = output_fields(field_id)%buffer(i,:) + &
             data(index_x,index_y,:)                                          ! accumulate buffer 
     else                                                                     ! not average
        output_fields(field_id)%buffer(i,:) = data(index_x,index_y,:)
     endif
  enddo
  if(output_fields(field_id)%time_average) &
       output_fields(field_id)%counter = output_fields(field_id)%counter + 1
end subroutine send_station_data_3d
!------------------------------------------------------------------------

subroutine station_data_out(file, field, data, time,final_call_in)

  integer, intent(in)          :: file, field
  real, intent(inout)          :: data(:, :)
  type(time_type), intent(in)  :: time
  logical, optional, intent(in):: final_call_in
  logical                      :: final_call
  integer                      :: i, num
  real :: dif, time_data(2, 1, 1), dt_time(1, 1, 1), start_dif, end_dif

  final_call = .false.
  if(present(final_call_in)) final_call = final_call_in
  dif = get_date_dif(time, base_time, files(file)%time_units)
  call mpp_write(files(file)%file_unit,output_fields(field)%f_type, data, dif)
  start_dif = get_date_dif(output_fields(field)%last_output, base_time,files(file)%time_units)
  end_dif = dif
  do i = 1, files(file)%num_fields
     num = files(file)%fields(i)     
      if(output_fields(num)%time_ops) then
         if(num == field) then
            time_data(1, 1, 1) = start_dif
            call mpp_write(files(file)%file_unit, files(file)%f_avg_start, &
                 time_data(1:1,:,:), dif)
            time_data(2, 1, 1) = end_dif
            call mpp_write(files(file)%file_unit, files(file)%f_avg_end, &
                 time_data(2:2,:,:), dif)
            dt_time(1, 1, 1) = end_dif - start_dif
            call mpp_write(files(file)%file_unit, files(file)%f_avg_nitems, &
                 dt_time(1:1,:,:), dif)
! Include boundary variable for CF compliance
            call mpp_write(files(file)%file_unit, files(file)%f_bounds, &
                 time_data(1:2,:,:), dif)
            exit
         endif
      end if
   end do
   if(final_call) then
      if(time >= files(file)%last_flush) then
         call mpp_flush(files(file)%file_unit)
         files(file)%last_flush = time
      endif
   else
      if(time > files(file)%last_flush) then
         call mpp_flush(files(file)%file_unit)
         files(file)%last_flush = time
      endif
   endif
end subroutine station_data_out
! <INTERFACE NAME="station_data_end">
! <TEMPLATE>
! station_data_end(time)
! </TEMPLATE>
!   <DESCRIPTION>
! Must be called <TT> after the last time step</TT> to write the buffer content
!   </DESCRIPTION>
! </INTERFACE>


!-----------------------------------------------------------------------------
subroutine station_data_end(time)

  type(time_type), intent(in) :: time            !model's time
  integer                     :: freq, max_counter, local_num_stations
  integer                     :: file, nfield, field, pe, col
  integer, allocatable        :: station_ids(:)  ! root_pe only,  to receive local station_ids
  real,    allocatable        :: tmp_buffer(:,:) ! root_pe only, to receive local buffer from each PE

  do file = 1, num_files
     freq = files(file)%output_freq
     do nfield = 1, files(file)%num_fields
        field = files(file)%fields(nfield)
        if(.not. output_fields(field)%register) cycle
        if(time >= output_fields(field)%next_output .or. freq == END_OF_RUN) then
! ALL PEs, including root PE, must send data to root PE        
           call mpp_send(output_fields(field)%num_station,plen=1,to_pe=mpp_root_pe(),tag=COMM_TAG_2)
           if(output_fields(field)%num_station > 0) then
              call mpp_send(output_fields(field)%station_id(1),plen=size(output_fields(field)%station_id),&
                   to_pe=mpp_root_pe(),tag=COMM_TAG_3)
              call mpp_send(output_fields(field)%buffer(1,1),plen=size(output_fields(field)%buffer),&
                   to_pe=mpp_root_pe(),tag=COMM_TAG_4)
           endif
! get max_counter if the field is averaged
           if(output_fields(field)%time_average) then
              max_counter = output_fields(field)%counter
              call mpp_max(max_counter, pelist)
           endif
! only root PE receives local data from all PEs 
           if(mpp_pe() == mpp_root_pe()) then
              do pe = 1,size(pelist)           
                 call mpp_recv(local_num_stations,glen=1,from_pe=pelist(pe),tag=COMM_TAG_2)
                 if(local_num_stations> 0) then
                    allocate(station_ids(local_num_stations))
                    allocate(tmp_buffer(local_num_stations,output_fields(field)%nlevel))
                    call mpp_recv(station_ids(1), glen=size(station_ids), from_pe=pelist(pe),tag=COMM_TAG_3)
                    call mpp_recv(tmp_buffer(1,1),glen=size(tmp_buffer),from_pe=pelist(pe),tag=COMM_TAG_4) 
                    do col = 1,local_num_stations
                       global_field%buffer(station_ids(col),:) = tmp_buffer(col,:)
                    enddo
                    deallocate(station_ids, tmp_buffer)
                 endif
              enddo
! send global_buffer content to file
              if(output_fields(field)%time_average)then           
                 if(max_counter == 0 )&
                      call error_mesg ('send_station_end','counter=0 for averaged field '// &
                      output_fields(field)%output_name, FATAL)
                 global_field%buffer = global_field%buffer/real(max_counter)
              endif
! check if global_field contains any missing values
              if(any(global_field%buffer == MISSING)) &
                   call error_mesg ('send_station_end','Global_field contains MISSING, field '// &
                   output_fields(field)%output_name, FATAL)              
              call station_data_out(file,field,global_field%buffer,output_fields(field)%next_output,.true.)
              global_field%buffer = MISSING
           endif
           call mpp_sync()
        endif
!deallocate field buffer
        if(output_fields(field)%num_station>0) &
             deallocate(output_fields(field)%buffer, output_fields(field)%station_id)
     enddo ! nfield
  enddo    ! file
  if(mpp_pe() == mpp_root_pe()) deallocate(global_field%buffer)
end subroutine station_data_end
!-----------------------------------------------------------------------------------

end module station_data_mod


