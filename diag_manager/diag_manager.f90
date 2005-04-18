module diag_manager_mod

! <CONTACT EMAIL="Matthew.Harrison@gfdl.noaa.gov">
!   Matt Harrison
! </CONTACT>
! <CONTACT EMAIL="Giang.Nong@noaa.gov">
!   Giang Nong
! </CONTACT>
! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

! <OVERVIEW>
!   <TT>diag_manager_mod</TT> is a set of simple calls for parallel diagnostics on
!   distributed systems. It is geared toward the writing of data in netCDF format.
! </OVERVIEW>

! <DESCRIPTION>
!   <TT>diag_manager_mod</TT> provides a convenient set of interfaces for
!   writing data to disk.  It is built upon the parallel I/O interface
!   <LINK SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/shared/mpp/mpp_io.html">mpp_io</LINK>.
!   A single group of calls to the <TT>diag_manager_mod</TT> interfaces provides data to disk
!   at any number of sampling and/or averaging intervals specified at run-time.
!   Run-time specification of diagnostics are input through the diagnostics table. <BR/>
! <B> Usage</B> of <TT> diag_manager</TT> includes the following steps:<BR/>
! 1. Create diag_table as described in 
!   <LINK SRC="diag_table_tk.html">diag_table_tk</LINK> documentation.<BR/>
! 2. Call <TT>diag_manager_init</TT> to initialize diag_manager_mod <BR/>
! 3. Call <TT> register_diag_field</TT> to register the field to be outputted. <BR/>
!    <B> NOTE</B> ALL fields in diag_table should be registered BEFORE the first send_data call
! 4. Call <TT>send_data</TT> to send data to output fields <BR/>
! 5. Call <TT>diag_manager_end</TT> to exit diag_manager <BR/><BR/>
!   
! <B> New Features </B> of <TT> diag_manager_mod </TT>: <BR/>
! 1. Ability to output scalars (see <TT>send_data </TT>and <TT>register_diag_field</TT>)<BR/>
! 2. Ability to output time average of fields that have time dependent mask. <BR/>
! 3. Give optional warning if <TT>register_diag_field </TT>fails due to misspelled module name
!    or field name. <BR/>
! 4. Check if a field is registered twice.<BR/>     
! 5. Check for duplicate lines in <LINK SRC="diag_table_tk.html">diag_table</LINK>. <BR/>
! 6. <LINK SRC="diag_table_tk.html">diag_table</LINK> can contain fields that are NOT written to any files. 
!    The file name in diag_table of these fields is <TT> null</TT> <BR/>
! 7. By default, a field is output in its global grid, it is now possible to output a field in 
!    a region specified by user, see <TT>send_data</TT> for more details <BR/>
! 8. To check if the diag table is set up correctly, user should set <TT>init_verbose=.true.</TT> in 
!    diag namelist, then the the content of diag_table is printed in stdout. <BR/>
! 9. New optional format of file information in <LINK SRC="diag_table_tk.html">diag_table</LINK><BR/>
!    It is possible to have just one file name and reuse it many times. A time string will be suffixed
!    to the base file name each time a new file is opened. The time string can be any combination from
!    year to second of current model time. Here is an example of file information: <BR/>
!    <TT>"file2_yr_dy%1yr%3dy",  5, "days",1,"hours","Time", 10, "days", "1 1 7 0 0 0"</TT>
!     <BR/>
!    From left to right we have: file name, output frequency, output frequency unit, Format (should always
!    be 1), time axis unit, time axis name, frequency for creating new file, unit for creating new file,
!    start time of the new file.<BR/>
!    In this example the string <TT>10, "days", "1 1 7 0 0 0"</TT> is optional.<BR/>
!    Keyword for the time string suffix is <TT>%xyr,%xmo,%xdy,%xhr,%xmi,%xsc</TT> where
!    <TT>x</TT> is mandatory 1 digit number specifying the width of field used in writing the string<BR/>
! 10. New time axis for time averaged fields <BR/>
!    users can use a namelist option to handle the time value written to time axis for time averaged fields.<BR/>
!    If  <TT>mix_snapshot_average_fields=.true.</TT> then a time averaged file will have time values 
!    corresponding to ending time_bound e.g. January monthly average is labeled Feb01. Users can have both 
!    snapshot and averaged fields in one file. <BR/>
!    If  <TT>mix_snapshot_average_fields=.false.</TT>
!    The time value written to time axis for time averaged fields is the middle on the averaging time.
!    For example, January monthly mean will be written at Jan.16 not Feb.01 as before. However, to use this 
!    new feature users should <B>separate</B> snapshot fields and time averaged fields in <B>different</B> files 
!    or a fatal error will occur.<BR/>
!    The namelist <B>default</B> value is <TT>mix_snapshot_average_fields=.false.</TT> <BR/><BR/>
! 11 Time average, Max and Min <BR/>
!    In addition to time average users can also get Max or Min value during the same interval of time as time
!    average. For this purpose, in the diag table users must replace <TT>.true.</TT> or <TT>.false.</TT> by 
!    <TT>"max"</TT> or <TT> "min"</TT>. <BR/>
!    Currently, Max and Min are not available for regional output. <BR/><BR/>

!   <B>Features of <TT>diag_manager_mod</TT> </B>include: <BR/>
!  
!   Integrated netCDF capability: <LINK
!   SRC="http://www.unidata.ucar.edu/packages/netcdf/">netCDF</LINK> is a
!   data format widely used in the climate/weather modeling
!   community. netCDF is considered the principal medium of data storage
!   for <TT>diag_manager_mod</TT>.
!
! </DESCRIPTION>

use time_manager_mod, only: get_time, set_time, get_date, set_date,    &
                            increment_date, operator(-), operator(>=), &
                            operator(>), operator(<), operator(==),    &
                            time_type, increment_time, month_name,     &
                            get_calendar_type, NO_CALENDAR, operator(/), operator(+), &
                            get_time, leap_year, GREGORIAN
use       mpp_io_mod, only: mpp_open, MPP_RDONLY, MPP_ASCII, mpp_close
use          fms_mod, only: error_mesg, FATAL, WARNING, NOTE,          &
                            close_file, stdlog, write_version_number,  &
                            file_exist, mpp_pe, open_namelist_file, &
                            check_nml_error, lowercase, stdout, mpp_error
use mpp_mod, only         : mpp_get_current_pelist, mpp_npes, mpp_sync, mpp_root_pe, stdout, mpp_sum    
use diag_axis_mod, only   : diag_axis_init, get_axis_length, get_diag_axis, get_domain1d, get_domain2d, &
     get_axis_global_length, diag_subaxes_init,get_diag_axis_cart, get_diag_axis_data, max_axes
use  diag_output_mod, only: diag_output_init, write_axis_meta_data, &
                            write_field_meta_data, done_meta_data,  &
                            diag_field_out, diag_output_end,        &
                            diag_flush, diag_fieldtype
use mpp_domains_mod, only : domain1d, domain2d, mpp_get_compute_domain, null_domain1d,&
     null_domain2d, operator(/=), mpp_modify_domain, mpp_get_domain_components
implicit none
private

public  diag_manager_init, send_data, send_tile_averaged_data, diag_manager_end,  &
        register_diag_field, register_static_field, &
        diag_axis_init, get_base_time, get_base_date, need_data, average_tiles
public :: DIAG_ALL,DIAG_OCEAN,DIAG_OTHER

! Specify storage limits for fixed size tables used for pointers, etc.
integer, parameter  :: max_fields_per_file = 150
integer, parameter  :: max_out_per_in_field = 10 !max number of output_fields per input_field
integer, parameter  :: max_files = 31
integer             :: num_files = 0
integer             :: max_input_fields = 300
integer             :: num_input_fields = 0
integer             :: max_output_fields = 300
integer             :: num_output_fields = 0
integer, parameter  :: DIAG_OTHER = 0
integer, parameter  :: DIAG_OCEAN = 1
integer, parameter  :: DIAG_ALL   = 2
integer, parameter  :: VERY_BIG_NUMBER = 10000
real                :: EMPTY = 0.0
integer             :: null_axis_id

! Global data for all files
type (time_type)    :: base_time
integer             :: base_year, base_month, base_day, base_hour, base_minute, base_second
character(len = 256):: global_descriptor
character(len = 64) :: iospec = "              "

type diag_grid
   real    :: start(3), end(3) ! coordinates (lat,lon,depth) of local domain to output   
   integer :: l_start_indx(3), l_end_indx(3) ! start and end indices at each LOCAL PE
   integer :: subaxes(3) ! id returned from diag_subaxes_init of 3 subaxes
end type diag_grid

type file_type
   character(len=128)  :: name
   integer             :: output_freq
   integer             :: output_units
   integer             :: format
   integer             :: time_units
   character(len=128)  :: long_name
   integer             :: fields(max_fields_per_file)
   integer             :: num_fields
   integer             :: file_unit
   integer             :: bytes_written
   integer             :: time_axis_id, time_bounds_id
   type (time_type)    :: last_flush
   type(diag_fieldtype):: f_avg_start, f_avg_end, f_avg_nitems, f_bounds
   logical             :: local ! true if fields are output in a region instead of global
   integer             :: new_file_freq ! frequency to create new file
   integer             :: new_file_freq_units ! time units of new_file_freq (days, hours, years, ...)
   type (time_type)    :: next_open ! time to open a new file
   type (time_type)    :: start_time ! time for opening the file for the first time 
end type file_type

type input_field_type
   character(len=128) :: module_name, field_name, long_name, units, standard_name
   integer            :: axes(3), num_axes 
   logical            :: missing_value_present, range_present
   real               :: missing_value, range(2)
   integer            :: output_fields(max_out_per_in_field)
   integer            :: num_output_fields, size(3) 
   logical            :: static, register, mask_variant, local
end type input_field_type

type output_field_type
   integer             :: input_field, output_file
   character(len=128)  :: output_name
   logical             :: time_average, static
   logical             :: time_max
   logical             :: time_min, time_ops
   integer             :: pack
   character(len=8)    :: time_method   
   real, pointer       :: buffer(:, :, :)=>null()
   real, pointer       :: counter(:, :, :)=>null()
   type(time_type)     :: last_output, next_output, next_next_output
   real                :: count_0d
   type(diag_fieldtype):: f_type
   integer             :: axes(3), num_axes, num_elements, total_elements, region_elements
   type(diag_grid)     :: output_grid
   logical             :: local_output, need_compute, phys_window
end type output_field_type

! local_output: true if this field is written out on a region, not global
! need_compute: true if this PE is involved in writing the field

type coord_type  ! define the region for outputting a field
   real :: xbegin
   real :: xend
   real :: ybegin
   real :: yend
   real :: zbegin
   real :: zend
end type coord_type

type (file_type), save :: files(max_files)
type (input_field_type),allocatable :: input_fields(:)
type (output_field_type), allocatable, save :: output_fields(:)
logical             :: append_pelist_name = .false.
logical             :: mix_snapshot_average_fields
logical             :: first_send_data_call = .true.
logical             :: first_send_data_local = .true.
logical             :: module_is_initialized = .false.
logical             :: do_diag_field_log = .false.
logical             :: write_bytes_in_file = .false.
integer             :: diag_log_unit
integer, parameter  :: EVERY_TIME =  0
integer, parameter  :: END_OF_RUN = -1
integer, parameter  :: DIAG_SECONDS = 1, DIAG_MINUTES = 2, DIAG_HOURS = 3
integer, parameter  :: DIAG_DAYS = 4, DIAG_MONTHS = 5, DIAG_YEARS = 6
character (len=10)  :: time_unit_list(6) = (/'seconds   ', 'minutes   ', &
   'hours     ', 'days      ', 'months    ', 'years     '/)
character (len = 7) :: avg_name = 'average'
character(len=32)   :: pelist_name

! version number of this module
character(len=128)  :: version = '$Id: diag_manager.f90,v 12.0 2005/04/14 17:55:35 fms Exp $'
character(len=128)  :: tagname = '$Name: lima $'  


! <INTERFACE NAME="send_data">
! <TEMPLATE>
!send_data(diag_field_id, field, time, is_in, js_in, ks_in,
!             mask, rmask, ie_in, je_in, ke_in, weight)
!</TEMPLATE>
!   <OVERVIEW>
!     Send data over to output fields. 
!   </OVERVIEW>
!   <DESCRIPTION>
! send_data is overloaded for fields having zero dimension (scalars) to 3 dimension.
! diag_field_id corresponds to the id returned from a previous call to
! register_diag_field.  The field array is restricted to the computational
! range of the array. Optional argument is_in can be used to update
! sub-arrays of the entire field.  Additionally, an optional logical or real
! mask can be used to apply missing values to the array.<BR/>
! If a field is declared to be mask_variant in <TT> register_diag_field</TT> logical mask 
! should be mandatory.<BR/>
!
! For the real  mask, the mask is applied if the mask value is less than 0.5.
!
! By default, a field will be written out entirely in its global grid. Users can also specify
! region in which the field will be output. The region is specified in diag-table just before
! the end of output_field replacing "none". For example:<BR/>
! by default:<BR/>
! "ocean_mod","Vorticity","vorticity","file1","all",.false.,"none",2 <BR/>
! for regional output:<BR/>
! "ocean_mod","Vorticity","vorticity_local","file2","all",.false.,"0.5 53.5 -89.5 -28.5 -1 -1",2<BR/>
! the format of region is "xbegin xend ybegin yend zbegin zend". If it is a 2D field use (-1 -1)
! for (zbegin zend) as in the example above. For a 3D field use (-1 -1) for (zbegin zend) when you want to
! write the whole vertical extent, otherwise specify real coordinates. The units used for region are the 
! actual units used in grid_spec.nc (for example degrees for lat, lon). Fatal error will occur if region's 
! boundaries are not found in grid_spec.nc<BR/>
! 
! <BR/>
! Special note when using regional output: result files containing regional outputs should be 
! different from files containing global (default) output. It is FATAL error to have one file
! containing both regional and global results. For maximum flexibility and independence from
! PE counts one file should contain just one region.<BR/>
! <BR/>
! Time averaging is supported in regional output.<BR/>
! <BR/>
! Physical fields (written in "physics windows" of atmospheric code) are 
! currently fully supported for regional outputs.<BR/><BR/>
! Note of dimension of field in send_data<BR/>
! Most fields are defined in data_domain but used in compute domain. In send_data users can pass EITHER
! field in data domain OR field in compute domain. If data domain is used, users should also pass the
! starting and ending indices of compute domain (isc, iec ...). If compute domain is used no indices
! are needed. These indices are for determining halo exclusively. If users want to ouput the field
! partially they should use regional output as mentioned above.<BR/><BR/>

! Weight in Time averaging is now supported, each time level may have a different weight. The default 
! of weight is 1.
! 
!   </DESCRIPTION>
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>

interface send_data
   module procedure send_data_0d
   module procedure send_data_1d
   module procedure send_data_2d
   module procedure send_data_3d
end interface

interface register_diag_field
   module procedure register_diag_field_scalar
   module procedure register_diag_field_array
end interface
! </INTERFACE>
! <INTERFACE NAME="send_tile_averaged_data">
!   <OVERVIEW>
!     Send tile-averaged data over to output fields. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     send_tile_averaged_data is overloaded for 3-d and 4-d arrays. 
!     diag_field_id corresponds to the id returned by previous call to 
!     register_diag_field. Logical mask can be used to mask out undefined
!     and/or unused values. Note that the dimension of output field is
!     smaller by one than the dimension of the data, since averaging over
!     tiles (3d dimension) is performed.
!   </DESCRIPTION>
!   <IN NAME="diag_field_id" TYPE="integer" >                </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:,:)" >  </IN>
!   <IN NAME="area" TYPE="real" DIM="(:,:,:)" >  </IN>
!   <IN NAME="time" TYPE="time_type" DIM="(:,:,:)" >  </IN>
!   <IN NAME="mask" TYPE="logical" DIM="(:,:,:)" >  </IN>
interface send_tile_averaged_data
   module procedure send_tile_averaged_data2d
   module procedure send_tile_averaged_data3d
end interface
!</INTERFACE>

contains

!-------------------------------------------------------------------------
function register_diag_field_scalar(module_name, field_name, init_time, &
   long_name, units, missing_value, range)

! Indicates the calling modules intent to supply data for this field.
integer                                ::  register_diag_field_scalar
character(len=*), intent(in)           :: module_name, field_name
type(time_type),  optional, intent(in) :: init_time
character(len=*), optional, intent(in) :: long_name, units
real, optional, intent(in)             :: missing_value, range(2)
 
if(present(init_time)) then
   register_diag_field_scalar = register_diag_field_array(module_name, field_name,&
        (/null_axis_id/), init_time,long_name, units, missing_value, range)
else
   register_diag_field_scalar = register_static_field(module_name, field_name,&
        (/null_axis_id/),long_name, units, missing_value, range)
endif

end function register_diag_field_scalar

! <FUNCTION NAME="register_diag_field">

!<OVERVIEW>
!     Register Diagnostic Field.
!</OVERVIEW>
!<DESCRIPTION>
!     Return field index for subsequent calls to <LINK SRC="#send_data"> send_data </LINK> <BR/>
!<TT> axes</TT> are axis ID returned from <TT>diag_axis_init</TT>, <TT>axes</TT>  are required
! for fields of 1-3 dimension and NOT required for scalars. <BR/>
! for a static scalar (constant) init_time is not needed. <BR/>
! optional <TT> mask_variant</TT> is for fields that have a time-dependent mask. If <TT>mask_variant</TT> is
! true then <TT>mask</TT> must be present in argument list of <TT>send_data</TT> <BR/>
! When optional <TT> verbose</TT> is true a warning will be given if <TT> register_diag_field</TT>
! fails due to misspelled field name or module name. The default <TT> verbose</TT> is false.
! <BR/>
! The pair (module_name, fieldname) should be registered only once or a FATAL error will occur
!
!</DESCRIPTION>
!<TEMPLATE>
!     register_diag_field(module_name,field_name,axes,init_time, &
!     long_name,units,missing_value,range,mask_variant,verbose)
!</TEMPLATE>
!
!   <IN NAME="module_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="field_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(:)"> </IN>
!   <IN NAME="init_time" TYPE="time_type"> </IN>
!   <IN NAME="long_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="units" TYPE="character(len=*)"> </IN>
!   <IN NAME="missing_value" TYPE="real"> </IN>
!   <IN NAME="range" TYPE="real" DIM="(2)"> </IN>
!   <IN NAME="mask_variant" TYPE="logical"> </IN> 
!   <IN NAME="verbose" TYPE="logical"> </IN>

function register_diag_field_array(module_name, field_name, axes, init_time, &
   long_name, units, missing_value, range, mask_variant,standard_name, verbose)

! Indicates the calling modules intent to supply data for this field.

integer                                :: register_diag_field_array
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
type(time_type), intent(in)            :: init_time
character(len=*), optional, intent(in) :: long_name, units, standard_name
real, optional, intent(in)             :: missing_value, range(2)
logical, optional, intent(in)          :: mask_variant, verbose
integer                                :: field, num_axes, i_size, j_size, k_size, j, ind, file_num, freq
integer                                :: output_units
logical                                :: mask_variant1, verbose1
character(len=128)                     :: msg
mask_variant1 = .false.
if(present(mask_variant)) mask_variant1 = mask_variant
verbose1= .false.
if(present(verbose)) verbose1 = verbose
! Call register static, then set static back to false

register_diag_field_array = register_static_field(module_name, field_name, axes, &
   long_name, units, missing_value, range, mask_variant1, dynamic =.true.)

if(verbose1.and.register_diag_field_array<0 ) &
     call error_mesg ('register_diag_field', &
     'module/output_field '//trim(module_name)//'/'//&
     &trim(field_name)//' NOT found in diag_table', WARNING) 

if(.not. first_send_data_call) call  error_mesg ('register_diag_field', &
     'module/output_field '//trim(module_name)//'/'//&
     &trim(field_name)//' registered AFTER first send_data call, TOO LATE', WARNING)  
if(register_diag_field_array >0) then
   input_fields(register_diag_field_array)%static = .false.

   field = register_diag_field_array
   if(present(standard_name))input_fields(field)%standard_name = standard_name  
   do j = 1, input_fields(field)%num_output_fields
      ind = input_fields(field)%output_fields(j)
      output_fields(ind)%static = .false.
! Set up times in output_fields
      output_fields(ind)%last_output = init_time
! Get output frequency from for the appropriate output file
      file_num = output_fields(ind)%output_file
      if(file_num == max_files) cycle
      if(output_fields(ind)%local_output ) then
         if(output_fields(ind)%need_compute) then         
            files(file_num)%local = .true.
         endif
         call mpp_sync()
      endif
! Need to sync start_time of file with init time of model
      if(files(file_num)%start_time < init_time) then
         files(file_num)%start_time = init_time
      endif
! Need to increase next_open until it is greater than init time
      do
         if(files(file_num)%next_open > init_time) exit
         files(file_num)%next_open = diag_time_inc(files(file_num)%next_open, &
              files(file_num)%new_file_freq, files(file_num)%new_file_freq_units)
      enddo
      freq = files(file_num)%output_freq
      output_units = files(file_num)%output_units
      output_fields(ind)%next_output = &
         diag_time_inc(init_time, freq, output_units)
      output_fields(ind)%next_next_output = &
         diag_time_inc(output_fields(ind)%next_output, freq, output_units)
      if(verbose1 .and. mpp_pe() == mpp_root_pe() .and. output_fields(ind)%local_output ) then
         write(msg,'(" lon(",F5.1,", ",F5.1,"), lat(",F5.1,", ",F5.1,")")')output_fields(ind)%output_grid%start(1),&
            output_fields(ind)%output_grid%end(1),output_fields(ind)%output_grid%start(2),&
            output_fields(ind)%output_grid%end(2)
         write(stdout(),* ) 'module/output_field '//trim(module_name)//'/'//trim(field_name)// &
              &' will be output in region:'//trim(msg)
      endif
   end do
endif
end function register_diag_field_array
! </FUNCTION>

!-------------------------------------------------------------------------
! <FUNCTION NAME="register_static_field">

!   <OVERVIEW>
!     Register Static Field.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Return field index for subsequent call to send_data.
!   </DESCRIPTION>
!   <TEMPLATE>
!     register_static_field(module_name, field_name, axes, &
!     long_name, units, missing_value, range, require)
!   </TEMPLATE>

!   <IN NAME="module_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="field_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(:)"> </IN>
!   <IN NAME="long_name" TYPE="character(len=*)"> </IN>
!   <IN NAME="units" TYPE="character(len=*)"> </IN>
!   <IN NAME="missing_value" TYPE="real"> </IN>
!   <IN NAME="range" TYPE="real" DIM="(2)"> </IN>

function register_static_field(module_name, field_name, axes, &
   long_name, units, missing_value, range, mask_variant, require, standard_name,dynamic)

integer                                :: register_static_field
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
character(len=*), optional, intent(in) :: long_name, units, standard_name
real, optional, intent(in)             :: missing_value, range(2)
logical, optional, intent(in)          :: mask_variant
logical, optional, intent(in)          :: require  ! require static field to be in every file, e.g. 2-d axes
logical, optional, intent(in)          :: dynamic
integer                                :: field, num_axes, j, out_num, file_num, freq, siz(3), local_siz(3), k
character(len=128)                     :: error_string
logical                                :: mask_variant1, dynamic1
integer                                :: local_start(3), local_end(3) ! indices of local domain of global axes
character(len=128)                     :: msg

mask_variant1 = .false.
if(present(mask_variant)) mask_variant1 = mask_variant
dynamic1 = .false.
if(present(dynamic)) dynamic1 = dynamic
if (.not.module_is_initialized) call error_mesg ('register_static_field',  &
                       'diag_manager has NOT been initialized', FATAL)

if ( do_diag_field_log) then
   call log_diag_field_info (module_name, field_name, axes, &
        long_name, units, missing_value=missing_value, range=range, &
        mask_variant=mask_variant1, require=require, dynamic=dynamic1)
endif

register_static_field = find_input_field(module_name, field_name)

if (PRESENT(require)) then
   if (require) then
      call init_input_field(module_name, field_name)
      register_static_field = find_input_field(module_name, field_name)
      do j=1, num_files
! need to think about what to do if the axes are not present, e.g. file only
! contains data slices
         call init_output_field(module_name, field_name,field_name, &
                               files(j)%name,time_method=".false.",pack=2)
      enddo
   endif
endif

! Negative index returned if this field is not used in table
if(register_static_field < 0) return

! Store information for this input field into input field table
field = register_static_field
! Set static to true, if called by register_diag_field this is flipped to false
input_fields(field)%static = .true.
! check if the field is registered twice
if (input_fields(field)%register .and. mpp_pe() == mpp_root_pe() ) then
    call error_mesg ('register_diag_field', &
     'module/output_field '//trim(module_name)//'/'//&
     &trim(field_name)//' ALREADY registered, should not register twice', FATAL) 
endif
! Set flag that this field was registered
input_fields(field)%register = .true.
! set flag for mask: does it change with time?
input_fields(field)%mask_variant = mask_variant1

! Store the axis info
num_axes = size(axes(:))
input_fields(field)%axes(1:num_axes) = axes
input_fields(field)%num_axes = num_axes
! Need to check for present, otherwise defaults
if(present(long_name)) then
   input_fields(field)%long_name = trim(long_name)
else
   input_fields(field)%long_name = input_fields(field)%field_name
endif
if(present(standard_name))input_fields(field)%standard_name = standard_name 
if(present(units)) then
   input_fields(field)%units = trim(units)
else
   input_fields(field)%units = 'none'
endif
if(present(missing_value)) then

   input_fields(field)%missing_value = missing_value
   input_fields(field)%missing_value_present = .true.
else
   input_fields(field)%missing_value_present = .false.
endif
if(present(range)) then
   input_fields(field)%range = range
   input_fields(field)%range_present = .true.
else
   input_fields(field)%range = (/ 1., 0. /)
   input_fields(field)%range_present = .false.
endif

siz = 1; local_siz = 1
local_start = 1;  local_end= 1
do j = 1, num_axes
   if(axes(j) .le. 0) then
      call error_mesg ('register_diag_field', &
           'module/output_field '//trim(module_name)//'/'//&
           &trim(field_name)//' has non-positive axis_id', FATAL) 
   endif
   siz(j) = get_axis_length(axes(j))
end do

! Default length for axes is 1
do j = 1, 3
   input_fields(field)%size(j) = siz(j)
end do

! Need to loop through all output_fields associated and allocate their buffers
do j = 1, input_fields(field)%num_output_fields
   out_num = input_fields(field)%output_fields(j)
! Range is required when pack >= 4 
   if(output_fields(out_num)%pack>=4 .and. .not.input_fields(field)%range_present) then
      if(mpp_pe() .eq. mpp_root_pe()) &
           call error_mesg ('register_diag_field ', 'output_field '//trim(field_name)// &
           ' has pack >=4, range is REQUIRED in register_diag_field', FATAL)
   endif
!  if local_output (size of output_fields does NOT equal size of input_fields)
   if(output_fields(out_num)%local_output) then
      if(size(axes(:)).le.1) &
           call error_mesg ('Register_diag_field', 'axes of '//trim(field_name)// &
           ' must >=2 for local output', FATAL)
      call get_subfield_size(axes,out_num)
      if(output_fields(out_num)%need_compute) then
         do k = 1,num_axes
            local_start(k) = output_fields(out_num)%output_grid%l_start_indx(k)
            local_end(k) = output_fields(out_num)%output_grid%l_end_indx(k)
            local_siz(k) = local_end(k) - local_start(k) +1                         
         enddo
         allocate(output_fields(out_num)%buffer(local_siz(1), local_siz(2), local_siz(3)))
         if(output_fields(out_num)%time_max) then
            output_fields(out_num)%buffer = -1*huge(1.)
         else if(output_fields(out_num)%time_min) then
            output_fields(out_num)%buffer = huge(1.) 
         else
            output_fields(out_num)%buffer = EMPTY
         endif
         output_fields(out_num)%region_elements = local_siz(1)*local_siz(2)*local_siz(3)
         output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
      endif
      call mpp_sync()     
   else ! the field is output globally
! size of output_fields equal size of input_fields 
      allocate(output_fields(out_num)%buffer(siz(1), siz(2), siz(3)))
      if(output_fields(out_num)%time_max) then
         output_fields(out_num)%buffer = -1*huge(1.)
      else if(output_fields(out_num)%time_min) then
         output_fields(out_num)%buffer = huge(1.) 
      else
         output_fields(out_num)%buffer = EMPTY
      endif
      output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
   endif
  
! Reset to false in register_field if this is not static
   output_fields(out_num)%static = .true.
! check if time average is true for static field
   if(.not.dynamic1 .and. output_fields(out_num)%time_ops) then      
      write(msg,'(a,"/",a)')trim(module_name), trim(field_name)
      if(mpp_pe() .eq. mpp_root_pe()) &
           call  error_mesg ('register_static_field', 'module/field '//trim(msg)//' is STATIC, can NOT time ops',WARNING)
      output_fields(out_num)%time_ops = .false.
      output_fields(out_num)%time_average = .false.
      output_fields(out_num)%time_method = 'point'
   endif
! Axes are copied from input_fields if output globally or from subaxes if output locally
   if(.not.output_fields(out_num)%local_output) then 
      output_fields(out_num)%axes = input_fields(field)%axes
   else      
      output_fields(out_num)%axes = output_fields(out_num)%output_grid%subaxes
   endif
! assume that the number of axes of output_fields = that of input_fields
   output_fields(out_num)%num_axes = input_fields(field)%num_axes
end do

if (input_fields(field)%mask_variant) then
   do j = 1, input_fields(field)%num_output_fields
      out_num = input_fields(field)%output_fields(j)
      if(output_fields(out_num)%time_average) then
         allocate(output_fields(out_num)%counter(siz(1), siz(2), siz(3)))
         output_fields(out_num)%counter = 0.0
      endif
   enddo
endif

end function register_static_field
! </FUNCTION>


! <SUBROUTINE NAME="log_diag_field_info">
!   <OVERVIEW>
!     Writes brief diagnostic field info to the log file.
!   </OVERVIEW>
!   <DESCRIPTION>
!     If <TT>init_verbose</TT> namelist parameter is true, adds a line briefly 
!     describing diagnostic field to the log file. Normally users should not call
!     this subroutine directly, since it is called by register_static_field and register_diag_field
!     if do_not_log is not set to true. It is used, however, in LM3 to avoid excessive
!     log due to a number of fields registered for each of the tile types. LM3 code uses
!     do_not_log parameter in the registration calls, and calls this subroutine to
!     log field information under generic name.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call log_diag_field_info ( module_name, field_name, axes, long_name, units,
!     missing_value, range, mask_variant, require, dynamic )
!   </TEMPLATE>
subroutine log_diag_field_info ( module_name, field_name, axes, long_name, units, &
                                 missing_value, range, mask_variant, require, dynamic )
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
character(len=*), optional, intent(in) :: long_name, units
real   , optional, intent(in)          :: missing_value, range(2)
logical, optional, intent(in)          :: mask_variant
logical, optional, intent(in)          :: require  ! require static field to be in every file, e.g. 2-d axes
logical, optional, intent(in)          :: dynamic

! ---- local vars
character(len=256) :: lmodule, lfield, lname, lunits
character(len=64) :: lmissval, lmin, lmax
character(len=8) :: numaxis, timeaxis
character(len=1) :: sep = '|'

if (.not.do_diag_field_log)    return
if (mpp_pe().ne.mpp_root_pe()) return

lmodule = trim(module_name)
lfield = trim(field_name)
lname  = ''; if(present(long_name)) lname  = trim(long_name)
lunits = ''; if(present(units))     lunits = trim(units)

write (numaxis,'(i1)') size(axes)

if (present(missing_value)) then
   write (lmissval,*) missing_value
else
   lmissval = ''
endif

if (present(range)) then
   write (lmin,*) range(1)
   write (lmax,*) range(2)
else
   lmin = ''
   lmax = ''
endif

if (present(dynamic)) then
   if (dynamic) then
      timeaxis = 'T'
   else
      timeaxis = 'F'
   endif
else
   timeaxis = ''
endif

!write (diag_log_unit,'(8(a,a),a)') &
write (diag_log_unit,'(17a)') &
             trim(lmodule),  sep, trim(lfield),  sep, trim(lname),    sep, &
             trim(lunits),   sep, trim(numaxis), sep, trim(timeaxis), sep, &
             trim(lmissval), sep, trim(lmin),    sep, trim(lmax)

end subroutine log_diag_field_info
! </SUBROUTINE>

subroutine get_subfield_size(axes,outnum)
! Get size, start and end indices for output_fields(outnum), fill in
! output_fields(outnum)%output_grid%(start_indx, end_indx)

integer, intent(in) :: axes(:) ! axes of the input_field
integer, intent(in) :: outnum  ! position in array output_fields
real, allocatable   :: global_lat(:), global_lon(:), global_depth(:)
integer             :: global_axis_size
integer             :: i,xbegin,xend,ybegin,yend,xbegin_l,xend_l,ybegin_l,yend_l 
character(len=1)    :: cart
type(domain2d)      :: Domain2, Domain2_new
type(domain1d)      :: Domain1,Domain1x,Domain1y,Domain1x_new,Domain1y_new
real                :: start(3), end(3) ! start and end coordinates in 3 axes
integer             :: gstart_indx(3), gend_indx(3) ! global start and end indices of output domain in 3 axes 
real, allocatable   :: subaxis_x(:), subaxis_y(:), subaxis_z(:) !containing local coordinates in x,y,z axes
character(len=128)  :: msg
!initilization for local output
start = -1.e10; end=-1.e10 ! initially out of (lat/lon/depth) range
gstart_indx = -1; gend_indx=-1

! get axis data (lat, lon, depth) and indices
   start= output_fields(outnum)%output_grid%start
   end = output_fields(outnum)%output_grid%end

do i = 1,size(axes(:))   
   global_axis_size = get_axis_global_length(axes(i))
   output_fields(outnum)%output_grid%subaxes(i) = -1
   call get_diag_axis_cart(axes(i), cart)
   select case(cart)
   case ('X')
      if(i.ne.1) &
           call error_mesg ('diag_manager, get subfield size', 'wrong order of axes, X should come first',FATAL)
      allocate(global_lon(global_axis_size))
      call get_diag_axis_data(axes(i),global_lon)
      gstart_indx(i) = get_index(start(i),global_lon)
      gend_indx(i) = get_index(end(i),global_lon)
      allocate(subaxis_x(gstart_indx(i):gend_indx(i)))
      subaxis_x=global_lon(gstart_indx(i):gend_indx(i))   
   case ('Y')
      if(i.ne.2) &
           call error_mesg ('diag_manager, get subfield size', 'wrong order of axes, Y should come second',FATAL)
      allocate(global_lat(global_axis_size))
      call get_diag_axis_data(axes(i),global_lat)
      gstart_indx(i) = get_index(start(i),global_lat)
      gend_indx(i) = get_index(end(i),global_lat)
      allocate(subaxis_y(gstart_indx(i):gend_indx(i)))
      subaxis_y=global_lat(gstart_indx(i):gend_indx(i))
   case ('Z')
      if(start(i)*end(i)<0) &
           call error_mesg ('diag_manager, get subfield size','wrong values in vertical axis of region',FATAL)
      if(start(i)>0 .and. end(i)>0) then 
         allocate(global_depth(global_axis_size))
         call get_diag_axis_data(axes(i),global_depth)
         gstart_indx(i) = get_index(start(i),global_depth)
         gend_indx(i) = get_index(end(i),global_depth)
         allocate(subaxis_z(gstart_indx(i):gend_indx(i)))
         subaxis_z=global_depth(gstart_indx(i):gend_indx(i))
         output_fields(outnum)%output_grid%subaxes(i) = &
              diag_subaxes_init(axes(i),subaxis_z, gstart_indx(i),gend_indx(i))
         deallocate(subaxis_z,global_depth)
      else ! regional vertical axis is the same as global vertical axis
         gstart_indx(i) = 1
         gend_indx(i) = global_axis_size
         output_fields(outnum)%output_grid%subaxes(i) = axes(i)
         if(i /= 3) &
              call error_mesg ('diag_manager, get subfield size','i should equal 3 for z axis',FATAL)
      endif      
   case default
       call error_mesg ('diag_manager, get_subfield_size', 'Wrong axis_cart', FATAL)
   end select
enddo
do i = 1,size(axes(:))
   if(gstart_indx(i)== -1 .or. gend_indx(i)== -1) &
        call error_mesg ('diag_manager, get_subfield_size', 'can not find gstart_indx/gend_indx for ' &
        //trim(output_fields(outnum)%output_name), FATAL)  
enddo

! get domain and compute_domain(xbegin,xend,ybegin,yend)
xbegin=-1; xend=-1
ybegin=-1; yend=-1

Domain2 = get_domain2d(axes)
if(Domain2 /= NULL_DOMAIN2D) then
   call mpp_get_compute_domain(Domain2,xbegin,xend,ybegin,yend)
   call mpp_get_domain_components(Domain2, Domain1x, Domain1y)
else
   do i = 1, MIN(size(axes(:)),2)    
      Domain1 = get_domain1d(axes(i))
      if(Domain1 /= NULL_DOMAIN1D) then
         call get_diag_axis_cart(axes(i),cart)
         select case(cart)
         case ('X')
            Domain1x = get_domain1d(axes(i))
            call mpp_get_compute_domain(Domain1x,xbegin,xend)   
         case ('Y')
            Domain1y = get_domain1d(axes(i))
            call mpp_get_compute_domain(Domain1y,ybegin,yend)
         case default ! do nothing here
         end select
      else
         call error_mesg ('diag_manager, get_subfield_size', 'NO domain available', FATAL)
      endif
   enddo
endif      
if(xbegin== -1 .or. xend==-1 .or. ybegin==-1 .or. yend==-1) &
   call error_mesg ('diag_manager, get_subfield_size', 'wrong compute domain indices',FATAL)  

! get the area containing BOTH compute domain AND local output area
if(gstart_indx(1)> xend .or. xbegin > gend_indx(1)) then
   output_fields(outnum)%output_grid%l_start_indx(1) = -1
   output_fields(outnum)%output_grid%l_end_indx(1) = -1
   output_fields(outnum)%need_compute = .false. ! not involved
elseif (gstart_indx(2)> yend .or. ybegin > gend_indx(2)) then
   output_fields(outnum)%output_grid%l_start_indx(2) = -1
   output_fields(outnum)%output_grid%l_end_indx(2) = -1
   output_fields(outnum)%need_compute = .false. ! not involved
else
   output_fields(outnum)%output_grid%l_start_indx(1) = MAX(xbegin, gstart_indx(1))
   output_fields(outnum)%output_grid%l_start_indx(2) = MAX(ybegin, gstart_indx(2))
   output_fields(outnum)%output_grid%l_end_indx(1) = MIN(xend, gend_indx(1))
   output_fields(outnum)%output_grid%l_end_indx(2) = MIN(yend, gend_indx(2))
   output_fields(outnum)%need_compute = .true.  ! involved in local output
endif

if(output_fields(outnum)%need_compute) then
! need to modify domain1d and domain2d for subaxes
   xbegin_l = output_fields(outnum)%output_grid%l_start_indx(1)
   xend_l = output_fields(outnum)%output_grid%l_end_indx(1)
   ybegin_l = output_fields(outnum)%output_grid%l_start_indx(2)
   yend_l = output_fields(outnum)%output_grid%l_end_indx(2)
   call mpp_modify_domain(Domain2, Domain2_new, xbegin_l,xend_l, ybegin_l,yend_l, &
                          gstart_indx(1),gend_indx(1), gstart_indx(2),gend_indx(2))
   call mpp_get_domain_components(Domain2_new, Domain1x_new, Domain1y_new)

   output_fields(outnum)%output_grid%subaxes(1) = &
        diag_subaxes_init(axes(1),subaxis_x, gstart_indx(1),gend_indx(1),Domain1x_new,Domain2_new)
   output_fields(outnum)%output_grid%subaxes(2) = &
        diag_subaxes_init(axes(2),subaxis_y, gstart_indx(2),gend_indx(2),Domain1y_new,Domain2_new)
   do i = 1, size(axes(:))
      if(output_fields(outnum)%output_grid%subaxes(i) == -1) then  
         write(msg,'(a,"/",I4)') 'at i = ',i
         call error_mesg ('diag_manager, get_subfield_size '//trim(output_fields(outnum)%output_name),&
              'error '//trim(msg), FATAL)   
      endif
   enddo
! local start index should start from 1
   output_fields(outnum)%output_grid%l_start_indx(1) = MAX(xbegin, gstart_indx(1)) - xbegin + 1   
   output_fields(outnum)%output_grid%l_start_indx(2) = MAX(ybegin, gstart_indx(2)) - ybegin + 1
   output_fields(outnum)%output_grid%l_end_indx(1) = MIN(xend, gend_indx(1)) - xbegin + 1 
   output_fields(outnum)%output_grid%l_end_indx(2) = MIN(yend, gend_indx(2)) - ybegin + 1
   if(size(axes(:))>2) then
      output_fields(outnum)%output_grid%l_start_indx(3) = gstart_indx(3)
      output_fields(outnum)%output_grid%l_end_indx(3) = gend_indx(3)
   else
      output_fields(outnum)%output_grid%l_start_indx(3) = 1
      output_fields(outnum)%output_grid%l_end_indx(3) = 1
   endif
endif
deallocate(subaxis_x, global_lon)
deallocate(subaxis_y, global_lat)

end subroutine get_subfield_size 

function get_index(number, array)
  real, intent(in) :: number
  real, intent(in), dimension(:) :: array
  integer get_index, i, n
  logical :: found
! Find index i of array such that array(i) is closest to number
! array must be  monotonouslly ordered

  n = size(array(:))
! check if array is monotonous
  do i = 2, n-1
     if((array(i-1)<array(i) .and. array(i)>array(i+1)).or.(array(i-1)>array(i) .and. array(i)<array(i+1))) &
          call error_mesg('diag_manager', 'get_index, array NOT monotonously ordered',FATAL) 
  enddo
  get_index = -1
  found = .false.
! search in increasing array 
  do i = 1, n-1                
     if ((array(i)<=number) .and. (array(i+1)>= number)) then
        if(number - array(i) <= array(i+1) - number) then
           get_index = i
           found=.true.
        else
           get_index = i+1
           found=.true.
        endif
        exit
     endif     
  enddo 
! if not found, search in decreasing array
  if(.not.found) then
      do i = 1, n-1
         if ((array(i)>=number) .and. (array(i+1)<= number)) then
            if(array(i)-number <= number-array(i+1)) then
               get_index = i              
            else
               get_index = i+1               
            endif
            exit
         endif
      enddo
   endif
  end function get_index


!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_0d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real"  > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_0d(diag_field_id, field, time)
    
logical                      :: send_data_0d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field
type (time_type), intent(in), optional :: time
real                         :: field_out(1, 1, 1)

! First copy the data to a three d array with last element 1
field_out(1, 1, 1) = field
send_data_0d = send_data_3d(diag_field_id, field_out, time)
end function send_data_0d

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_1d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight)

logical                      :: send_data_1d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:)
real, intent(in), optional   :: weight
type (time_type), intent(in), optional :: time
integer, optional            :: is_in, ie_in
logical, optional            :: mask(:)
real, optional               :: rmask(:)
real                         :: field_out(size(field(:)), 1, 1)
logical                      :: mask_out(size(field(:)), 1, 1)


! First copy the data to a three d array with last element 1
field_out(:, 1, 1) = field

! Default values for mask
mask_out = .true.
if(present(mask)) mask_out(:, 1, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, 1, 1) = .false.
if(present(mask) .or. present(rmask)) then
   send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in, 1, 1, mask_out, ie_in=ie_in,weight=weight)
else
   send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in, 1, 1, ie_in=ie_in, weight=weight)
endif
end function send_data_1d


!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_2d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_2d(diag_field_id, field, time, is_in, js_in, &
    mask, rmask, ie_in, je_in, weight)

logical                      :: send_data_2d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:, :)
real, intent(in), optional   :: weight
type (time_type), intent(in), optional :: time
integer, optional            :: is_in, js_in, ie_in, je_in
logical, optional            :: mask(:, :)
real, optional               :: rmask(:, :)
real                         :: field_out(size(field, 1), size(field, 2), 1)
logical                      :: mask_out(size(field, 1), size(field, 2), 1)


! First copy the data to a three d array with last element 1
field_out(:, :, 1) = field

! Default values for mask
mask_out = .true.
if(present(mask)) mask_out(:, :, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, :, 1) = .false.
if(present(mask) .or. present(rmask)) then
   send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in, js_in, 1, mask_out,&
        ie_in=ie_in, je_in=je_in,weight=weight)
else
   send_data_2d = send_data_3d(diag_field_id, field_out,time,is_in,js_in,1,ie_in=ie_in,je_in=je_in,weight=weight)
endif

end function send_data_2d

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_3d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             mask, rmask, ie_in,je_in, ke_in,weight)

logical                      :: send_data_3d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:,:,:)
real, intent(in), optional   :: weight
type (time_type), intent(in), optional :: time
integer, optional            :: is_in, js_in, ks_in,ie_in,je_in, ke_in 
logical, optional            :: mask(:, :, :)
real, optional               :: rmask(:, :, :)
real                         :: num, weight1
logical                      :: average, need_compute, local_output, phys_window
logical                      :: time_max, time_min
integer                      :: i, out_num, file_num, n1, n2, n3, number_of_outputs, ii,f1,f2,f3,f4
integer                      :: freq, units, is, js, ks, kount, ie, je, ke, i1, j1,k1, j, k
character(len=128)           :: error_string
integer                      :: l_start(3), l_end(3) ! local start and end indices on 3 axes for regional output
integer                      :: hi, hj  ! halo size in x and y direction
integer                      :: b1,b2,b3 ! size of buffer in x,y,z axes
type (time_type)             :: middle_time      
if (.not.module_is_initialized) call error_mesg ('send_data_3d',' diag_manager NOT initialized', FATAL)
! If is, js, or ks not present default them to 1
is = 1
if(present(is_in)) is = is_in
js = 1
if(present(js_in)) js = js_in
ks = 1
if(present(ks_in)) ks = ks_in
n1 = size(field, 1); n2 = size(field, 2); n3 = size(field, 3)
ie = is+n1-1; je = js+n2-1; ke = ks+n3-1
if (present(ie_in)) ie = ie_in
if (present(je_in)) je = je_in
if (present(ke_in)) ke = ke_in
if(present(is_in).and.present(ie_in).and.present(js_in).and.present(je_in)) then
   hi = (n1-(ie-is+1))/2; hj = (n2-(je-js+1))/2  ! compute halo in x and y direction
   is = 1 + hi; ie = n1-hi; js=1+hj; je=n2-hj
else
   hi = 0; hj = 0
endif
f1=1+hi; f2=n1-hi; f3=1+hj; f4=n2-hj ! used for field, mask and rmask bounds
! weight is for time averaging where each time level may has a different weight
weight1 = 1.
if(present(weight)) weight1 = weight

! If diag_field_id is < 0 it means that this field is not registered, simply return
if(diag_field_id < 0) then
   send_data_3d = .false.
   return
else
   send_data_3d = .true.
endif
if(input_fields(diag_field_id)%local) then
! need to increase number_of_outputs by 1 for mpp_sync() in case of local output
   number_of_outputs = input_fields(diag_field_id)%num_output_fields + 1  
else
   number_of_outputs = input_fields(diag_field_id)%num_output_fields
endif

! Loop through each output field that depends on this input field
do ii = 1, number_of_outputs
   if(ii == input_fields(diag_field_id)%num_output_fields + 1) then
      call mpp_sync()
      exit
   endif
! Get index to an output field
   out_num = input_fields(diag_field_id)%output_fields(ii)
! is this field output on a local domain only?
   local_output = output_fields(out_num)%local_output
! if local_output, does the current PE take part in send_data?
   need_compute = output_fields(out_num)%need_compute
! skip all PEs not participating in outputting this field
   if(local_output .and. (.not. need_compute)) cycle
! Get index to output file for this field
   file_num = output_fields(out_num)%output_file
   if(file_num == max_files) cycle
! Output frequency and units for this file is
   freq = files(file_num)%output_freq
   units = files(file_num)%output_units
! Is this output field being time averaged?
   average = output_fields(out_num)%time_average
! Looking for max and min value of this field over the sampling interval?
   time_max = output_fields(out_num)%time_max
   time_min = output_fields(out_num)%time_min   
   if(output_fields(out_num)%total_elements > size(field(f1:f2,f3:f4,ks:ke))) then
       output_fields(out_num)%phys_window = .true.
   else
      output_fields(out_num)%phys_window = .false.
   endif
   phys_window = output_fields(out_num)%phys_window
   if(need_compute) then        
      l_start = output_fields(out_num)%output_grid%l_start_indx
      l_end = output_fields(out_num)%output_grid%l_end_indx
   endif
   
! Initialize output time for fields output every time step
   if (freq == EVERY_TIME) then
     if (output_fields(out_num)%next_output == output_fields(out_num)%last_output) then
       if(present(time)) then
         output_fields(out_num)%next_output = time
       else
         write (error_string,'(a,"/",a)')  &
             trim(input_fields(diag_field_id)%module_name), &
             trim(output_fields(out_num)%output_name)
          call error_mesg('diag_manager send_data ', &
             'module/output_field '//trim(error_string)//&
             ', time must be present when output frequency = EVERY_TIME', FATAL)
       endif
     endif
   endif

   if(.not.output_fields(out_num)%static .and. .not.present(time)) then
         write (error_string,'(a,"/",a)')  &
             trim(input_fields(diag_field_id)%module_name), &
             trim(output_fields(out_num)%output_name)
          call error_mesg('diag_manager send_data ', &
             'module/output_field '//trim(error_string)//&
             ', time must be present for nonstatic field', FATAL)
   endif

! Is it time to output for this field; CAREFUL ABOUT > vs >= HERE
   if(.not.output_fields(out_num)%static .and. freq /= END_OF_RUN) then
      if(time > output_fields(out_num)%next_output ) then
! A non-static field that has skipped a time level is an error
      if(time >  output_fields(out_num)%next_next_output .and. freq > 0) then
         if(mpp_pe() .eq. mpp_root_pe()) then 
            write (error_string,'(a,"/",a)')  &
                 trim(input_fields(diag_field_id)%module_name), &
                 trim(output_fields(out_num)%output_name)
            call error_mesg('Warning1, diag_manager send_data ', &
                 'module/output_field '//trim(error_string)//&
                 ' is skipped one time level in output data', WARNING)
         endif
      end if
! If average get size: Average intervals are last_output, next_output
      if(average) then
         b1=size(output_fields(out_num)%buffer,1); b2=size(output_fields(out_num)%buffer,2) 
         b3=size(output_fields(out_num)%buffer,3)
         if (input_fields(diag_field_id)%mask_variant) then
            if (any(output_fields(out_num)%counter>0.)) then
               do i=1,b1; do j=1,b2; do k=1,b3
                  if(output_fields(out_num)%counter(i,j,k)>0.) output_fields(out_num)%buffer(i,j,k) = &
                       output_fields(out_num)%buffer(i,j,k)/output_fields(out_num)%counter(i,j,k)
               enddo;enddo;enddo
            else
               if(any(output_fields(out_num)%buffer /= input_fields(diag_field_id)%missing_value)) then
                  write (error_string,'(a,"/",a)')  &
                       trim(input_fields(diag_field_id)%module_name), &
                       trim(output_fields(out_num)%output_name)
                  call error_mesg ('error2 diag_manager send_data', &
                       'module/output_field '//trim(error_string)//&
                       &', write EMPTY buffer', FATAL) 
               endif
            endif
         else
            if(phys_window) then
               if(need_compute) then
                  num = real(output_fields(out_num)%num_elements/output_fields(out_num)%region_elements)
               else
                  num = real(output_fields(out_num)%num_elements/output_fields(out_num)%total_elements)
               endif
            else
               num = output_fields(out_num)%count_0d
            endif
            if (num>0.) then
               if (input_fields(diag_field_id)%missing_value_present) then
                  do i=1,b1; do j=1,b2; do k=1,b3
                     if(output_fields(out_num)%buffer(i,j,k)/= input_fields(diag_field_id)%missing_value) &
                          output_fields(out_num)%buffer(i,j,k) = output_fields(out_num)%buffer(i,j,k)/num  
                  enddo;enddo;enddo
               else
                  output_fields(out_num)%buffer = output_fields(out_num)%buffer/num
               endif
            else
               if (input_fields(diag_field_id)%missing_value_present) then
                  if(any(output_fields(out_num)%buffer /= input_fields(diag_field_id)%missing_value)) then
                     write (error_string,'(a,"/",a)')  &
                          trim(input_fields(diag_field_id)%module_name), &
                          trim(output_fields(out_num)%output_name)
                     if(mpp_pe() .eq. mpp_root_pe())  call error_mesg ('error3 diag_manager_mod', &
                          'module/output_field '//trim(error_string)//&
                          &', write EMPTY buffer', FATAL)
                  endif
               endif
            endif
         endif ! mask_variant
      endif !average
! Output field
      if((output_fields(out_num)%time_ops).and.(.not. mix_snapshot_average_fields)) then
         middle_time = (output_fields(out_num)%last_output+output_fields(out_num)%next_output)/2
         call diag_data_out(file_num, out_num, &
              output_fields(out_num)%buffer, middle_time)
      else
         call diag_data_out(file_num, out_num, &
              output_fields(out_num)%buffer, output_fields(out_num)%next_output)
      endif
! Take care of cleaning up the time counters and the storeage size
      output_fields(out_num)%last_output = output_fields(out_num)%next_output
      if (freq == EVERY_TIME) then
         output_fields(out_num)%next_output = time
      else
         output_fields(out_num)%next_output = output_fields(out_num)%next_next_output
         output_fields(out_num)%next_next_output = &
              diag_time_inc(output_fields(out_num)%next_next_output, freq, units)
      endif
      output_fields(out_num)%count_0d = 0.0
      output_fields(out_num)%num_elements = 0
      if (time_max) then 
          output_fields(out_num)%buffer = -1*HUGE(1.0)
      else if (time_min) then
          output_fields(out_num)%buffer = HUGE(1.0)
      else
          output_fields(out_num)%buffer = EMPTY
      endif      
      if(input_fields(diag_field_id)%mask_variant .and. average) output_fields(out_num)%counter = 0.0
   end if  !time > output_fields(out_num)%next_output
   end if  !.not.output_fields(out_num)%static .and. freq /= END_OF_RUN
! Finished output of previously buffered data, now deal with buffering new data   
 
! Take care of submitted field data
   if(average) then
      if (input_fields(diag_field_id)%mask_variant) then
         if(need_compute) then
            write (error_string,'(a,"/",a)')  &
                 trim(input_fields(diag_field_id)%module_name), &
                 trim(output_fields(out_num)%output_name)   
            call error_mesg ('error4 diag_manager_mod', 'module/output_field '//trim(error_string)//&
                 &', regional output NOT supported with mask_variant', FATAL)
         endif
         if (present(mask)) then
            if (input_fields(diag_field_id)%missing_value_present) then              
               do i=is,ie; do j=js,je; do k=ks,ke
                  if(mask(i-is+1+hi,j-js+1+hj,k)) then
                     output_fields(out_num)%buffer(i-hi,j-hj,k)=output_fields(out_num)%buffer(i-hi,j-hj,k) + &
                          field(i-is+1+hi,j-js+1+hj,k)*weight1  
                     output_fields(out_num)%counter(i-hi,j-hj,k)=output_fields(out_num)%counter(i-hi,j-hj,k) + weight1
                  else
                     output_fields(out_num)%buffer(i-hi,j-hj,k)= input_fields(diag_field_id)%missing_value
                  endif
               enddo; enddo; enddo

            else
               write (error_string,'(a,"/",a)')  &
                    trim(input_fields(diag_field_id)%module_name), &
                    trim(output_fields(out_num)%output_name)
               call error_mesg ('error5 diag_manager_mod', 'module/output_field '//trim(error_string)//&
                    &', variable mask but no missing value defined', FATAL)
            endif
         else
            write (error_string,'(a,"/",a)')  &
                 trim(input_fields(diag_field_id)%module_name), &
                 trim(output_fields(out_num)%output_name)
            call error_mesg ('error6 diag_manager_mod', 'module/output_field '//trim(error_string)//&
                 &', variable mask but no mask given', FATAL)
         endif         
      else ! mask_variant=false
         if (present(mask)) then
            if (input_fields(diag_field_id)%missing_value_present) then
               if(need_compute) then
                  do i = is, ie
                     do j = js,je                       
                        if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                           i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                           do k = l_start(3),l_end(3)
                              k1=k-l_start(3)+1
                              if(mask(i-is+1+hi,j-js+1+hj,k)) then
                                 output_fields(out_num)%buffer(i1,j1,k1) = output_fields(out_num)%buffer(i1,j1,k1)+ &
                                      field(i-is+1+hi,j-js+1+hj,k)*weight1                              
                              else
                                 output_fields(out_num)%buffer(i1,j1,k1)=input_fields(diag_field_id)%missing_value   
                              endif
                              output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + 1
                           enddo
                        endif
                     enddo
                  enddo
               else
                  do i=is,ie; do j=js,je; do k=ks,ke
                     if(mask(i-is+1+hi,j-js+1+hj,k)) then
                        output_fields(out_num)%buffer(i-hi,j-hj,k)=output_fields(out_num)%buffer(i-hi,j-hj,k) + &
                             field(i-is+1+hi,j-js+1+hj,k)*weight1  
                     else
                        output_fields(out_num)%buffer(i-hi,j-hj,k)= input_fields(diag_field_id)%missing_value
                     endif
                  enddo; enddo; enddo

               endif
               if(need_compute.and. .not.phys_window) then
                  if(any(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3)))) &
                       output_fields(out_num)%count_0d=output_fields(out_num)%count_0d + weight1
               else
                  if(any(mask(f1:f2,f3:f4,ks:ke))) output_fields(out_num)%count_0d=output_fields(out_num)%count_0d+weight1
               endif               
            else ! missing value NOT present
               if(.not.all(mask(f1:f2,f3:f4,ks:ke)) .and. mpp_pe() .eq. mpp_root_pe()) &
                    call error_mesg('warning2 send_data_3d', &
                    'Mask will be ignored since missing values were not specified',WARNING)
               if(need_compute) then
                   do i = is, ie
                     do j = js,je                       
                        if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                           i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1                                                
                           output_fields(out_num)%buffer(i1,j1,:)= output_fields(out_num)%buffer(i1,j1,:)+ &
                                field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           output_fields(out_num)%num_elements=output_fields(out_num)%num_elements+l_end(3)-l_start(3)+1
                        endif
                     enddo
                  enddo
               else                          
                  output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = &
                       output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) &
                       + field(f1:f2,f3:f4,ks:ke)*weight1   
               endif
               if(.not. phys_window) output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1
            endif
         else ! mask NOT present
            if (input_fields(diag_field_id)%missing_value_present) then
               if(need_compute) then 
                  do i = is, ie
                     do j = js,je                       
                        if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                           i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1 
                           do k = l_start(3),l_end(3)
                              k1=k-l_start(3)+1
                              if(field(i-is+1,j-js+1,k) /= input_fields(diag_field_id)%missing_value) then
                                 output_fields(out_num)%buffer(i1,j1,k1)= output_fields(out_num)%buffer(i1,j1,k1)+ &
                                   field(i-is+1+hi,j-js+1+hj,k)*weight1
                              else
                                 output_fields(out_num)%buffer(i1,j1,k1)= input_fields(diag_field_id)%missing_value
                              endif
                              output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + 1
                           enddo
                        endif
                     enddo
                  enddo
                  if(.not. phys_window) then
                     if(any(field(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3)) /= &
                          input_fields(diag_field_id)%missing_value)) &
                          output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1 
                  endif
               else
                  do i=is,ie; do j=js,je; do k=ks,ke
                     if(field(i-is+1+hi,j-js+1+hj,k) /= input_fields(diag_field_id)%missing_value)  then
                        output_fields(out_num)%buffer(i-hi,j-hj,k)=output_fields(out_num)%buffer(i-hi,j-hj,k) + &
                             field(i-is+1+hi,j-js+1+hj,k)*weight1  
                     else
                        output_fields(out_num)%buffer(i-hi,j-hj,k)= input_fields(diag_field_id)%missing_value
                     endif
                  enddo; enddo; enddo                  

                  if(any(field(f1:f2,f3:f4,ks:ke) /= input_fields(diag_field_id)%missing_value)) &
                       output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1    
               endif
            else       ! no missing value defined, No mask
               if(need_compute) then
                  do i = is, ie
                     do j = js,je                       
                        if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                           i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1                           
                           output_fields(out_num)%buffer(i1,j1,:)= output_fields(out_num)%buffer(i1,j1,:)+ &
                                field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           output_fields(out_num)%num_elements=output_fields(out_num)%num_elements+l_end(3)-l_start(3)+1 
                        endif
                     enddo
                  enddo
               else 
                  output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = &
                       output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) + field(f1:f2,f3:f4,ks:ke)*weight1
               endif
               if(.not. phys_window) output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1
            endif
         endif ! if mask present
      endif  !if mask_variant
      if(.not. need_compute) &
           output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + (ie-is+1)*(je-js+1)*(ke-ks+1)
! Add processing for Max and Min
   else if (time_max) then
      if (present(mask)) then
         where(mask(f1:f2,f3:f4,ks:ke).and.field(f1:f2,f3:f4,ks:ke)>output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke))&
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
      else
         where (field(f1:f2,f3:f4,ks:ke) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke)) &
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
      endif
      output_fields(out_num)%count_0d = 1
   else if (time_min) then
      if (present(mask)) then
         where(mask(f1:f2,f3:f4,ks:ke).and.field(f1:f2,f3:f4,ks:ke)<output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke)) &
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke) 
      else
         where (field(f1:f2,f3:f4,ks:ke) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke)) &
              output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
      endif
      output_fields(out_num)%count_0d = 1
   else  ! ( not average)
      output_fields(out_num)%count_0d = 1
      if(need_compute) then
         do i = is, ie
            do j = js,je 
               if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                  i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                  output_fields(out_num)%buffer(i1,j1,:) = field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))
               endif
            enddo
         enddo
      else
         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
      endif
               
      if (present(mask) .and. input_fields(diag_field_id)%missing_value_present) then
         if(need_compute) then
            do i = is, ie
               do j = js,je                       
                  if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                     i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                     do k = l_start(3),l_end(3)
                        k1=k-l_start(3)+1
                        if(.not. mask(i-is+1+hi,j-js+1+hj,k)) &
                             output_fields(out_num)%buffer(i1,j1,k1) = input_fields(diag_field_id)%missing_value
                     enddo
                  endif
               enddo
            enddo
         else
            do i=is,ie; do j=js,je; do k=ks,ke
               if(.not. mask(i-is+1+hi,j-js+1+hj,k)) &
                    output_fields(out_num)%buffer(i-hi,j-hj,k)= input_fields(diag_field_id)%missing_value
            enddo; enddo; enddo            

         endif
      endif
   endif !average
 
! If rmask and missing value present, then insert missing value     
   if (present(rmask) .and. input_fields(diag_field_id)%missing_value_present) then
      if(need_compute) then
         do i = is, ie
            do j = js,je                       
               if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                  i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                  do k = l_start(3),l_end(3)
                     k1=k-l_start(3)+1
                     if(rmask(i-is+1+hi,j-js+1+hj,k)<0.5) &
                          output_fields(out_num)%buffer(i1,j1,k1) = input_fields(diag_field_id)%missing_value                  
                  enddo
               endif
            enddo
         enddo
      else
         do i=is,ie; do j=js,je; do k=ks,ke
            if( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5) &
                 output_fields(out_num)%buffer(i-hi,j-hj,k)= input_fields(diag_field_id)%missing_value
         enddo; enddo; enddo  

      endif
   endif
end do

end function send_data_3d


!-------------------------------------------------------------------------
! <FUNCTION NAME="send_tile_averaged_data_2d" INTERFACE="send_tile_averaged_data">
!   <IN NAME="diag_field_id" TYPE="integer" >               </IN>
!   <IN NAME="field"         TYPE="real"      DIM="(:,:,:)">  </IN>
!   <IN NAME="area"          TYPE="real"      DIM="(:,:,:)">  </IN>
!   <IN NAME="time"          TYPE="time_type" DIM="(:,:,:)">  </IN>
!   <IN NAME="mask"          TYPE="logical"   DIM="(:,:,:)">  </IN>
! </FUNCTION>

function send_tile_averaged_data2d ( id, field, area, time, mask )

  ! --- return value ---------------------------------------------------------
  logical                      :: send_tile_averaged_data2d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id             ! id od the diagnostic field 
  real,    intent(in)          :: field(:,:,:)   ! field to average and send
  real,    intent(in)          :: area (:,:,:)   ! area of tiles (== averaging 
                                                 ! weights), arbitrary units
  type(time_type), intent(in)  :: time           ! current time
  logical, intent(in),optional :: mask (:,:,:)   ! land mask

  ! --- local vars -----------------------------------------------------------
  real  :: out(size(field,1), size(field,2))

  call average_tiles( field, area, mask, out )
  send_tile_averaged_data2d = send_data( id, out, time, mask=ANY(mask,DIM=3) )
end function send_tile_averaged_data2d


!-------------------------------------------------------------------------
! <FUNCTION NAME="send_tile_averaged_data3d" INTERFACE="send_tile_averaged_data">
!   <IN NAME="diag_field_id" TYPE="integer">                </IN>
!   <IN NAME="field"         TYPE="real"      DIM="(:,:,:,:)">  </IN>
!   <IN NAME="area"          TYPE="real"      DIM="(:,:,:)">  </IN>
!   <IN NAME="time"          TYPE="time_type" DIM="(:,:,:)">  </IN>
!   <IN NAME="mask"          TYPE="logical"   DIM="(:,:,:)">  </IN>
! </FUNCTION>

function send_tile_averaged_data3d( id, field, area, time, mask )

  ! --- return value ---------------------------------------------------------
  logical                      :: send_tile_averaged_data3d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id              ! id of the diagnostic field
  real,    intent(in)          :: field(:,:,:,:)  ! (lon, lat, tile, lev) field 
                                                  ! to average and send
  real,    intent(in)          :: area (:,:,:)    ! (lon, lat, tile) tile areas 
                                                  ! ( == averaging weights), 
                                                  ! arbitrary units
  type(time_type), intent(in)  :: time            ! current time
  logical, intent(in),optional :: mask (:,:,:)    ! (lon, lat, tile) land mask

  ! --- local vars -----------------------------------------------------------
  real    :: out(size(field,1), size(field,2), size(field,4))
  logical :: mask3(size(field,1), size(field,2), size(field,4))
  integer :: it

  do it=1,size(field,4)
     call average_tiles( field(:,:,:,it), area, mask, out(:,:,it) )
  enddo

  mask3(:,:,1) = ANY(mask,DIM=3)
  do it = 2, size(field,4)
     mask3(:,:,it) = mask3(:,:,1)
  enddo

  send_tile_averaged_data3d = send_data( id, out, time, mask=mask3 )
end function send_tile_averaged_data3d

! ============================================================================
subroutine average_tiles ( x, area, mask, out )
! ============================================================================
! average 2-dimensional field over tiles
  ! --- arguments ------------------------------------------------------------
  real,    intent(in)  :: x   (:,:,:) ! (lon, lat, tile) field to average
  real,    intent(in)  :: area(:,:,:) ! (lon, lat, tile) fractional area
  logical, intent(in)  :: mask(:,:,:) ! (lon, lat, tile) land mask
  real,    intent(out) :: out (:,:)   ! (lon, lat)       result of averaging

  ! --- local vars -----------------------------------------------------------------
  integer  :: it                      ! iterator over tile number
  real     :: s(size(x,1),size(x,2))  ! area accumulator

  s(:,:)   = 0.0
  out(:,:) = 0.0

  do it = 1,size(area,3)
     where (mask(:,:,it)) 
        out(:,:) = out(:,:) + x(:,:,it)*area(:,:,it)
        s(:,:)   = s(:,:) + area(:,:,it)
     endwhere
  enddo

  where( s(:,:) > 0 ) &
       out(:,:) = out(:,:)/s(:,:)

end subroutine average_tiles

!-------------------------------------------------------------------------
! <SUBROUTINE NAME="diag_manager_end">
!
!   <OVERVIEW>
!     Exit Diagnostics Manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Flushes diagnostic buffers where necessary. Close diagnostics files.<BR/>
!  A warning will be issued here if a field in diag_table is not registered
!   </DESCRIPTION>
!   <TEMPLATE>
!     call diag_manager_end (time)
!   </TEMPLATE>
!   <IN NAME="TIME" TYPE="time_type"></IN>
subroutine diag_manager_end (time)

type(time_type), intent(in) :: time
integer                     :: file

if (do_diag_field_log) then
   call mpp_close (diag_log_unit)
endif

do file = 1, num_files
      call closing_file(file, time)   
end do
end subroutine diag_manager_end
! </SUBROUTINE>
!-------------------------------------------------------------------------

subroutine closing_file(file, time)
! This subroutine replaces diag_manager_end
! close just one file: files(file)

integer, intent(in)         :: file
type(time_type), intent(in) :: time
integer                     :: j, i, input_num, freq
real                        ::  num 
character(len=128)          :: message
type(time_type)             :: middle_time
logical                     :: local_output, need_compute, phys_window

! Output all registered, non_static output_fields
do j = 1, files(file)%num_fields
   i = files(file)%fields(j) !this is position of output_field in array output_fields

! is this field output on a local domain only?
   local_output = output_fields(i)%local_output
! if local_output, does the current PE take part in send_data?
   need_compute = output_fields(i)%need_compute
   phys_window = output_fields(i)%phys_window
! skip all PEs not participating in outputting this field
   if(local_output .and. (.not. need_compute)) return
! skip fields that were not registered or non-static   
   input_num = output_fields(i)%input_field
   if (input_fields(input_num)%static) cycle
   if (.not.input_fields(input_num)%register) cycle 
   freq = files(file)%output_freq
   if(freq /= END_OF_RUN .and. files(file)%file_unit <0 &
        .and. output_fields(i)%num_elements == 0 .and. output_fields(i)%count_0d == 0 ) cycle
! Is it time to output for this field; CAREFUL ABOUT >= vs > HERE
! For end should be >= because no more data is coming 
   if(time >= output_fields(i)%next_output .or. freq == END_OF_RUN) then
      if(time >= output_fields(i)%next_next_output .and. freq > 0) then
         write (message,'(a,"/",a)') trim(input_fields(input_num)%module_name), &
              trim(output_fields(i)%output_name)
         if(mpp_pe() .eq. mpp_root_pe()) & 
              call error_mesg('diag_manager_end, closing_file', 'module/output_field ' //  &
              &trim(message)//', skip one time level, maybe send_data never called', WARNING)
      end if

! If average, get size: Need time average intervals here
      if(output_fields(i)%time_average) then
         if (input_fields(input_num)%mask_variant) then
            if (any(output_fields(i)%counter>0.)) then
               where(output_fields(i)%counter>0.) &
                    output_fields(i)%buffer = output_fields(i)%buffer/output_fields(i)%counter
            else
               if(any(output_fields(i)%buffer /= input_fields(input_num)%missing_value)) then
                  write (message,'(a,"/",a)')  &
                       trim(input_fields(input_num)%module_name), &
                       trim(output_fields(i)%output_name)
                  call error_mesg ('diag_manager_mod, closing_file','module/output_field '//trim(message)//&
                       &', write EMPTY buffer', FATAL) 
               endif
            endif
         else            
            if(phys_window) then
               if(need_compute) then
                  num = real(output_fields(i)%num_elements/output_fields(i)%region_elements)
               else
                  num = real(output_fields(i)%num_elements/output_fields(i)%total_elements)
               endif
            else
               num = output_fields(i)%count_0d
            endif
            if (num>0.) then
               if (input_fields(input_num)%missing_value_present) then
                  where (output_fields(i)%buffer /= input_fields(input_num)%missing_value) &
                       output_fields(i)%buffer = output_fields(i)%buffer/num
               else
                  output_fields(i)%buffer = output_fields(i)%buffer/num
               endif           
            endif
         endif
      endif
! Output field
      if (freq == END_OF_RUN) output_fields(i)%next_output = time
      if( (output_fields(i)%time_ops).and.(.not. mix_snapshot_average_fields)) then
         middle_time = (output_fields(i)%last_output+output_fields(i)%next_output)/2
         call diag_data_out(file, i, output_fields(i)%buffer, middle_time, .true.)
      else
         call diag_data_out(file, i, &
              output_fields(i)%buffer, output_fields(i)%next_output, .true.)
      endif
   end if 
end do
! Now it's time to output static fields
call  write_static(file)

! Write out the number of bytes of data saved to this file
if (write_bytes_in_file) then
   call mpp_sum (files(file)%bytes_written)
   write (stdout(), '(a,i12,a,a)') 'Diag_Manager: ',files(file)%bytes_written, &
        ' bytes of data written to file ',trim(files(file)%name)
endif
end subroutine closing_file
!---------------------------------------------------------------

subroutine write_static(file)
! Output all static fields in this file
integer, intent(in) :: file
integer             :: j, i, input_num

do j = 1, files(file)%num_fields
   i = files(file)%fields(j)
   input_num = output_fields(i)%input_field
! skip fields that were not registered
   if (.not.input_fields(input_num)%register) cycle
! only output static fields here
   if (.not.output_fields(i)%static) cycle
   call diag_data_out(file, i, output_fields(i)%buffer, files(file)%last_flush, .true., .true.)
end do
! Close up this file   
call mpp_close(files(file)%file_unit)
files(file)%file_unit = -1
end subroutine write_static

!-------------------------------------------------------------------------

!##################################################################################
subroutine init_file(name, output_freq, output_units, format, time_units, long_name, &
     new_file_freq, new_file_freq_units, start_time)

character(len=*), intent(in)          :: name, long_name
integer, intent(in)                   :: output_freq, output_units, format, time_units
integer, intent(in), optional         :: new_file_freq, new_file_freq_units
type(time_type), intent(in), optional :: start_time
integer                               :: new_file_freq1, new_file_freq_units1
real, dimension(1)                    :: tdata
character(len=128)                    :: time_units_str
! Get a number for this file
num_files = num_files + 1
if(num_files >= max_files) then
   call error_mesg('diag_manager, init_file', ' max_files exceeded, incease max_files', FATAL)
endif

new_file_freq1 = VERY_BIG_NUMBER
if(present(new_file_freq)) new_file_freq1 = new_file_freq
if (get_calendar_type() == NO_CALENDAR) then
  new_file_freq_units1 = DIAG_DAYS
else 
   new_file_freq_units1 = DIAG_YEARS
endif
if(present(new_file_freq_units)) new_file_freq_units1 = new_file_freq_units 
files(num_files)%name = trim(name)
files(num_files)%output_freq = output_freq
files(num_files)%output_units = output_units
files(num_files)%format = format
files(num_files)%time_units = time_units
files(num_files)%long_name = trim(long_name)
files(num_files)%num_fields = 0
files(num_files)%local = .false.
files(num_files)%last_flush = base_time
files(num_files)%file_unit = -1
files(num_files)%new_file_freq = new_file_freq1
files(num_files)%new_file_freq_units = new_file_freq_units1
if(present(start_time)) then 
   files(num_files)%start_time = start_time
else
   files(num_files)%start_time = base_time
endif
files(num_files)%next_open=diag_time_inc(files(num_files)%start_time,new_file_freq1,new_file_freq_units1)
! add time_axis_id and time_bounds_id here
write(time_units_str, 11) trim(time_unit_list(files(num_files)%time_units)), base_year, &
     base_month, base_day, base_hour, base_minute, base_second
11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
files(num_files)%time_axis_id = diag_axis_init (trim(long_name), tdata, time_units_str, 'T',  &
     trim(long_name) , set_name=trim(name) )
!---- register axis for storing time boundaries
files(num_files)%time_bounds_id = diag_axis_init( 'nv',(/1.,2./),'none','N','vertex number',&
                              set_name=trim(name))
end subroutine init_file

!--------------------------------------------------------------------------

subroutine init_input_field(module_name, field_name)

character(len=*), intent(in) :: module_name, field_name

! Get a number for this input_field if not already set up
if(find_input_field(module_name, field_name) < 0) then
   num_input_fields = num_input_fields + 1
   if(num_input_fields > max_input_fields) then
      call error_mesg('diag_manager,init_input_field', 'max_input_fields exceeded, increase it via diag_manager_nml',&
           FATAL)
   end if
else
! If this is already initialized don't need to do anything
   return
end if

input_fields(num_input_fields)%module_name = trim(module_name)
input_fields(num_input_fields)%field_name = trim(field_name)
input_fields(num_input_fields)%num_output_fields = 0
! Set flag that this field has not been registered
input_fields(num_input_fields)%register = .false.
input_fields(num_input_fields)%local = .false.
input_fields(num_input_fields)%standard_name = 'none'
end subroutine init_input_field

!---------------------------------------------------------------------------

subroutine init_output_field(module_name, field_name, output_name, output_file,&
   time_method, pack, local_coord)

character(len=*), intent(in)           :: module_name, field_name, output_name, output_file
character(len=*), intent(in)           :: time_method
integer, intent(in)                    :: pack
type(coord_type), intent(in), optional :: local_coord
integer                                :: out_num, in_num, file_num, &
     num_fields, i, method_selected, l1
character(len=128)                     :: error_msg
character(len=8)                       :: t_method

! Get a number for this output field
num_output_fields = num_output_fields + 1
if(num_output_fields > max_output_fields) then
   call error_mesg('diag_manager', 'max_output_fields exceeded, increase it via diag_manager_nml', FATAL)
endif
out_num = num_output_fields

! First, find the index to the associated input field
in_num = find_input_field(module_name, field_name)
if(in_num < 0) then
   write (error_msg,'(a,"/",a)') trim(module_name),trim(field_name)
   call error_mesg('diag_manager,init_output_field', &
      'module_name/field_name '//trim(error_msg)//' NOT registered', FATAL)
endif

! Add this output field into the list for this input field
input_fields(in_num)%num_output_fields = &
   input_fields(in_num)%num_output_fields + 1
if(input_fields(in_num)%num_output_fields > max_out_per_in_field) then
   call error_mesg('diag_manager,init_output_field', 'max_out_per_in_field exceeded, increase max_out_per_in_field', FATAL)
endif
input_fields(in_num)%output_fields(input_fields(in_num)%num_output_fields) &
   = out_num

! Also put pointer to input field in this output field
output_fields(out_num)%input_field = in_num

! Next, find the number for the corresponding file
if(trim(output_file).eq.'null') then
   file_num = max_files
else
   file_num = find_file(output_file)
   if(file_num < 0) then
      call error_mesg('diag_manager,init_output_field', 'file '//trim(output_file)//' is NOT found in diag_table', FATAL)
   end if
endif

! Insert this field into list for this file
files(file_num)%num_fields = files(file_num)%num_fields + 1
if(files(file_num)%num_fields > max_fields_per_file) then
   call error_mesg('diag_manager,init_output_field', 'max_fields_per_file exceeded, increase max_fields_per_file ', FATAL)
endif
num_fields = files(file_num)%num_fields
files(file_num)%fields(num_fields) = out_num
output_fields(out_num)%count_0d = 0

! Set the file for this output field
output_fields(out_num)%output_file = file_num

! Enter the other data for this output field
output_fields(out_num)%output_name = trim(output_name)
output_fields(out_num)%pack = pack
output_fields(out_num)%num_elements = 0
output_fields(out_num)%total_elements = 0
output_fields(out_num)%region_elements = 0

method_selected = 0
! Initialize all time method to false
output_fields(out_num)%time_average = .false.
output_fields(out_num)%time_min = .false.
output_fields(out_num)%time_max = .false. 
output_fields(out_num)%time_ops = .false.

! cannot time average fields output every time
if (files(file_num)%output_freq == EVERY_TIME) then
   output_fields(out_num)%time_average = .false.
   method_selected = method_selected+1
else
   t_method = lowercase(time_method)
   select case (trim(t_method))
   case('.true.')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('mean')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('average')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('avg')
      output_fields(out_num)%time_average = .true.
      method_selected = method_selected+1
      t_method = 'mean'
   case('.false.')
      output_fields(out_num)%time_average = .false.
      method_selected = method_selected+1
      t_method = 'point'
   case('none')
      output_fields(out_num)%time_average = .false.
      method_selected = method_selected+1
      t_method = 'point'
   case('point')
      output_fields(out_num)%time_average = .false.
      method_selected = method_selected+1
      t_method = 'point'   
   case ('maximum')
      output_fields(out_num)%time_max = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'max') &
           output_fields(out_num)%output_name = trim(output_name)//'_max'
      method_selected = method_selected+1
      t_method = 'max'        
   case ('max')
      output_fields(out_num)%time_max = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'max') &
           output_fields(out_num)%output_name = trim(output_name)//'_max'
      method_selected = method_selected+1
      t_method = 'max'
   case ('minimum')
      output_fields(out_num)%time_min = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'min') &
           output_fields(out_num)%output_name = trim(output_name)//'_min'
      method_selected = method_selected+1
      t_method = 'min'        
   case ('min')
      output_fields(out_num)%time_min = .true.
      l1 = len_trim(output_fields(out_num)%output_name)
      if(output_fields(out_num)%output_name(l1-2:l1) /= 'min') &
           output_fields(out_num)%output_name = trim(output_name)//'_min'
      method_selected = method_selected+1
      t_method = 'min'     
   end select
endif
output_fields(out_num)%time_ops = output_fields(out_num)%time_min.or.output_fields(out_num)%time_max &
     .or.output_fields(out_num)%time_average

output_fields(out_num)%phys_window = .false.
! need to initialize grid_type = -1(start, end, l_start_indx,l_end_indx etc...)
if(present(local_coord)) then
   input_fields(in_num)%local = .true.
   output_fields(out_num)%output_grid%start(1) = local_coord%xbegin
   output_fields(out_num)%output_grid%start(2) = local_coord%ybegin
   output_fields(out_num)%output_grid%start(3) = local_coord%zbegin
   output_fields(out_num)%output_grid%end(1) = local_coord%xend
   output_fields(out_num)%output_grid%end(2) = local_coord%yend
   output_fields(out_num)%output_grid%end(3) = local_coord%zend
   do i = 1,3
      output_fields(out_num)%output_grid%l_start_indx(i) = -1
      output_fields(out_num)%output_grid%l_end_indx(i) = -1
      output_fields(out_num)%output_grid%subaxes(i) = -1
   enddo
   output_fields(out_num)%local_output = .true.
   output_fields(out_num)%need_compute = .false.
else
   output_fields(out_num)%local_output = .false.
   output_fields(out_num)%need_compute = .false.
endif

if (method_selected /= 1) call error_mesg('init_output_field','improper &
     &time method in diag_table for output field:'//trim(output_name),FATAL)
output_fields(out_num)%time_method = trim(t_method)

end subroutine init_output_field

!-------------------------------------------------------------------------

! <SUBROUTINE NAME="diag_manager_init">

!   <OVERVIEW>
!     Initialize Diagnostics Manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Open and read diag_table. Select fields and files for diagnostic output.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call diag_manager_init()
!   </TEMPLATE>

subroutine diag_manager_init(diag_model_subset)
integer, optional, intent(IN) :: diag_model_subset
integer                       :: diag_subset_output

type tableB_type
   character(len=128) :: module_name,field_name,output_name,name
   character(len=50)  :: time_sampling
   character(len=8)   :: time_method   
   character(len=50)  :: spatial_ops
   integer            :: pack
end type tableB_type

type tableA_type
   character(len=128) :: name
   integer            :: output_freq
   character(len=10)  :: output_freq_units
   integer            :: format
   character(len=10)  :: time_units
   character(len=128) :: long_name
   integer            :: new_file_freq
   character(len=10)  :: new_file_freq_units
   character(len=25)  :: start_time_s  
end type tableA_type

character(len=256) :: record
character(len=9)   :: amonth
integer            :: iunit,n,m,num_fields,time_units, output_freq_units, nfiles,nfields
integer            :: j, log_unit, name_len, nrecs, ierr, io_status
integer, allocatable, dimension(:) :: pelist
integer            :: new_file_freq_units, record_len, time_pos
logical            :: new_file_type, start_time_present, init_verbose
integer            :: yr, mo, dy, hr, mi, sc
type(time_type)    :: start_time

namelist /diag_manager_nml/ append_pelist_name, mix_snapshot_average_fields, init_verbose,iospec, &
     max_output_fields, max_input_fields, max_axes, do_diag_field_log, write_bytes_in_file

type(tableB_type)  :: textB
type(tableA_type ) :: textA
type(coord_type)   :: local_coord !local coordinates used in local output

!  If the module was already initialized do nothing
if (module_is_initialized) return
init_verbose = .false.
diag_subset_output = DIAG_ALL
if (PRESENT(diag_model_subset))then
  if(diag_model_subset>=DIAG_OTHER .AND. diag_model_subset<=DIAG_ALL) then
    diag_subset_output = diag_model_subset
  else
    call error_mesg('diag_manager_init','file diag_table nonexistent',FATAL)
  endif
endif
mix_snapshot_average_fields = .false.
call mpp_open(iunit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
read(iunit,diag_manager_nml,iostat=io_status)
write(stdlog(), diag_manager_nml)
if (io_status > 0) then
   call error_mesg('diag_manager_init', 'Error reading diag_manager_nml',FATAL)
endif
call mpp_close (iunit)
if(mix_snapshot_average_fields) then
   if(mpp_pe() == mpp_root_pe()) call mpp_error(WARNING,'Namelist '// &
        'mix_snapshot_average_fields = true will cause ERROR in time coordinates '// &
        'of all time_averaged fields. Strongly recommend mix_snapshot_average_fields = false')
endif
allocate(output_fields(max_output_fields))
allocate(input_fields(max_input_fields))
if (.not.file_exist('diag_table') ) &
call error_mesg('diag_manager_init','file diag_table nonexistent',FATAL)

call mpp_open(iunit, 'diag_table', action=MPP_RDONLY)

! Read in the global file labeling string
read(iunit, *, end = 99, err=99) global_descriptor

! Read in the base date
read(iunit, *, end = 99, err = 99) base_year, base_month, base_day, &
   base_hour, base_minute, base_second

! Set up the time type for base time
if (get_calendar_type() /= NO_CALENDAR) then
   base_time = set_date(base_year, base_month, base_day, base_hour, &
                        base_minute, base_second)
   amonth = month_name(base_month)
else
! No calendar - ignore year and month
   base_time = set_time(base_hour*3600+base_minute*60+base_second, base_day)
   base_year  = 0
   base_month = 0
   amonth = 'day'
end if

allocate(pelist(mpp_npes()))
call mpp_get_current_pelist(pelist, pelist_name)

nrecs=0
nfiles=0
if(init_verbose) then
   write(stdout(), *) ' '
   write(stdout(), *) '******** Summary of output FILE information from diag_manager *********'
   write(stdout(), *) ' '
   write(stdout(), *) 'file name     ', '       saving frequency   ', '      saving frequency unit'
   write(stdout(), *) ' '
endif
do while (nfiles <= max_files)
   read(iunit,'(a)',end=86,err=85) record
   nrecs=nrecs+1
   if (record(1:1) == '#') cycle   
   record_len = len_trim(record)
   if(record_len == 0) cycle
   time_pos = MAX(index(record,'time'), index(record,'Time'))
   if (time_pos + 6 >= record_len) then
      new_file_type = .false.
      start_time_present = .false.
   else
      new_file_type = .true.
      if(INDEX(record,'s',.true.) +3 < record_len) then 
         start_time_present = .true.
      else
         start_time_present = .false.
      endif
   endif
! Start reading file information  
   if(.not. new_file_type) then
      read(record,*,err=85,end=85)textA%name,textA%output_freq,textA%output_freq_units, &
           textA%format,textA%time_units,textA%long_name   
   elseif(.not.start_time_present) then
      read(record,*,err=85,end=85)textA%name,textA%output_freq,textA%output_freq_units, &
           textA%format,textA%time_units,textA%long_name,textA%new_file_freq,textA%new_file_freq_units
   else
      read(record,*,err=85,end=85) textA
   endif

! does the record contain start time for the file ???   
   if(start_time_present .and. new_file_type) then
      read(textA%start_time_s,*,end=85,err=85) yr, mo, dy, hr, mi, sc
      start_time = set_date( yr, mo, dy, hr, mi, sc)
   else
      start_time = base_time
   endif

! test file format to make sure its OK
   if (textA%format .gt. 2 .or. textA%format .lt. 1) cycle
   if( diag_subset_output==DIAG_OTHER .AND. verify( 'ocean',lowercase(textA%name) )==0 )cycle
   if( diag_subset_output==DIAG_OCEAN .AND. verify( 'ocean',lowercase(textA%name) )/=0 )cycle
   nfiles=nfiles+1
   time_units = 0
   output_freq_units = 0
   new_file_freq_units = 0
   do j = 1, size(time_unit_list(:))
      if(textA%time_units == time_unit_list(j)) time_units = j
      if(textA%output_freq_units == time_unit_list(j)) output_freq_units = j
      if(new_file_type) then
         if(textA%new_file_freq_units == time_unit_list(j)) new_file_freq_units = j
      endif
   end do
   if(time_units == 0) &
        call error_mesg('diag_manager_init','invalid time units, check time unit in diag_table',FATAL)
   if(output_freq_units == 0) & 
        call error_mesg('diag_manager_init','invalid output frequency units, check diag table',FATAL)
   if (new_file_type .and. new_file_freq_units == 0 ) &
        call error_mesg('diag_manager_init','invalid new_file frequency units, check diag table',FATAL)
   ! remove trailing .nc extension from file name 
   name_len = len_trim(textA%name)
   if (textA%name(name_len-2:name_len) == '.nc') then
       textA%name = textA%name(1:name_len-3)
       name_len = name_len - 3
   endif
   ! add optional suffix based on pelist name
   if (append_pelist_name) then
      textA%name(name_len+1:) = trim(pelist_name)
   endif   
   ! assign values to file_types
   if(new_file_type) then
      call init_file(textA%name,textA%output_freq, output_freq_units, &
           textA%format, time_units,textA%long_name,textA%new_file_freq,new_file_freq_units, start_time)
   else
      call init_file(textA%name,textA%output_freq, output_freq_units, &
           textA%format, time_units,textA%long_name)
   endif
   if(init_verbose) write(stdout(), 1)textA%name, textA%output_freq, textA%output_freq_units 
85 continue
enddo
call error_mesg('diag_manager_init','max_files exceeded, increase max_files', FATAL)
86 continue
1 format(1x,A21, 4x,I3,4x,A30)
write(stdout(), *)' '
if(init_verbose) write(stdout(), *)'************************************************************************'
rewind(iunit)

if(init_verbose) then
   write(stdout(), *) ' '
   write(stdout(), *) '******* Summary of output FIELD information from diag_manager **********'
   write(stdout(), *) ' '
   write(stdout(), *) 'module         ',' field name      ', '   file name      ', '    time method '
   write(stdout(), *) ' '
endif


nfields=0;nrecs=0
do while (nfields <= max_output_fields)
   read(iunit,'(a)',end=94,err=93) record
   nrecs=nrecs+1
   if (record(1:1) == '#') cycle
   read(record,*,end=93,err=93) textB
   if (textB%pack .gt. 8 .or. textB%pack .lt. 1) cycle
   if( diag_subset_output==DIAG_OTHER .AND. verify( 'ocean',lowercase(textB%name) )==0 )cycle
   if( diag_subset_output==DIAG_OCEAN .AND. verify( 'ocean',lowercase(textB%name) )/=0 )cycle
   if(textB%spatial_ops /= 'none') read(textB%spatial_ops,*,end=93,err=93)local_coord
   nfields=nfields+1
   !   assign values to field_types
   call init_input_field(textB%module_name,textB%field_name)
   !   remove trailing .nc extension
   name_len= len_trim(textB%name)
   if (textB%name(name_len-2:name_len) == '.nc') then
       textB%name = textB%name(1:name_len-3)
       name_len = name_len-3
   endif
   ! add optional suffix based on pelist name
   if (append_pelist_name) then
       textB%name(name_len+1:) = trim(pelist_name)
   endif
   if(trim(textB%spatial_ops) == 'none') then
      call init_output_field(textB%module_name,textB%field_name,textB%output_name,&
           textB%name,textB%time_method,textB%pack)
   else
      call init_output_field(textB%module_name,textB%field_name,textB%output_name,&
           textB%name,textB%time_method,textB%pack, local_coord)
   endif
   if(init_verbose) &
        write(stdout(),2)textB%module_name,textB%field_name,textB%name,textB%time_method 
93 continue
enddo
call error_mesg('diag_manager_init','max_output_fields exceeded, increase it via diag_manager_nml', FATAL)
94 continue
call close_file(iunit)
2 format(1x,A15,1x,A15,4x,A21,A10)
if(init_verbose) write(stdout(), *)'************************************************************************'
! check duplicate output_fields in the diag_table
call check_duplicate_output_fields
!initialize files%bytes_written to zero
files(:)%bytes_written = 0

! open diag field log file
if (do_diag_field_log) then
   call mpp_open (diag_log_unit, 'diag_field_log.out', nohdrs=.TRUE.)
endif

! version number to logfile
call write_version_number (version, tagname)

log_unit = stdlog()
if ( mpp_pe() == 0 ) then
   write (log_unit,95) base_year, trim(amonth), base_day, &
        base_hour, base_minute, base_second
endif
call close_file (log_unit)
95 format ('base date used = ',i4,1x,a,2i3,2(':',i2.2),' gmt')
module_is_initialized = .true.
! create axis_id for scalars here
null_axis_id= diag_axis_init( 'scalar_axis', (/0./), 'none', 'X', 'none')
return
99 continue
call error_mesg('diag_manager_init','error reading table',FATAL)
end subroutine diag_manager_init
! </SUBROUTINE>

!-------------------------------------------------------------------------

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

!-------------------------------------------------------------------------

function find_input_field(module_name, field_name)
integer find_input_field
character(len=*), intent(in) :: module_name, field_name
integer                      :: i

find_input_field = -1
do i = 1, num_input_fields
   if(trim(input_fields(i)%module_name) == trim(module_name) .and. &
      lowercase(trim(input_fields(i)%field_name)) == &
      lowercase(trim(field_name))) then 
      find_input_field = i
      return
   endif
end do
end function find_input_field

!-------------------------------------------------------------------------

subroutine opening_file(file, time)

! WARNING: Assumes that all data structures are fully initialized
  integer, intent(in)           :: file
  type(time_type), intent(in)   :: time  
  integer                       :: i, j, field_num, n_fields, axes(5), input_field_num, num_axes, k
  character(len=128)            :: time_units, timeb_units, avg, error_string, filename
  logical                       :: time_ops
  integer                       :: time_axis_id(1), field_num1, time_bounds_id(1)
  character(len=7)              :: prefix
  character(len=128)            :: suffix, base_name
  integer                       :: yr, mo, dy, hr, mi, sc, position
  character(len=32)             :: time_name, timeb_name,time_longname, timeb_longname, cart_name
  type(domain1d)                :: domain
  integer                       :: dir, edges
  real, dimension(2)            :: data

! First, get a file_unit and a time axis identifier for each file
  i = file
! it's unlikely that a file starts with word "rregion", need to check anyway.
  if (len(files(i)%name) >=7 .and. .not. files(i)%local) then
     prefix = files(i)%name(1:7)
     if(lowercase(prefix) == 'rregion') &
          call error_mesg ('diag_manager opening_file', 'file name should not start with' &
          //' word "rregion"', WARNING)
  endif
  
! Here is where time_units string must be set up; time since base date
  write(time_units, 11) trim(time_unit_list(files(i)%time_units)), base_year, &
       base_month, base_day, base_hour, base_minute, base_second
11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2)
  base_name = files(i)%name
  if(files(i)%new_file_freq < VERY_BIG_NUMBER) then
     position = INDEX(files(i)%name, '%')
     if(position>0) then
        base_name = base_name(1:position-1)
     else
        call error_mesg ('diag_manager opening_file', 'file name '//TRIM(files(i)%name)// &
             ' does not contain % for time stamp string', FATAL) 
     endif
     suffix = get_time_string(files(i)%name, time)
  else
     suffix = ' '
  endif 
! Add CVS tag as prefix of filename  (currently not implemented)
!  i1 = INDEX(tagname,':') + 2
!  i2 = len_trim(tagname) - 2
!  if(i2 <=i1)  call error_mesg ('diag_manager opening_file','error in CVS tagname index',FATAL)
!  prefix2 = tagname(i1:i2)//'_'

  if(files(i)%local) then      
! prefix "rregion" to all local files for post processing, the prefix will be removed in postprocessing
     filename = 'rregion'//trim(base_name)//trim(suffix)
  else
!     filename = trim(prefix2)//trim(base_name)//trim(suffix)
     filename = trim(base_name)//trim(suffix)
  endif
  call diag_output_init(filename, files(i)%format, global_descriptor, &
       files(i)%long_name, time_units, files(i)%file_unit,iospec) 
  files(i)%bytes_written = 0 
! Does this file contain time_average fields?
  time_ops = .false.
  do j = 1, files(i)%num_fields
     field_num = files(i)%fields(j)
     if(output_fields(field_num)%time_ops) then
        time_ops = .true.
        exit
     endif
  enddo
! Loop through all fields with this file to output axes
  do j = 1, files(i)%num_fields
     field_num = files(i)%fields(j)
     input_field_num = output_fields(field_num)%input_field
     if (.not.input_fields(input_field_num)%register) then
        write (error_string,'(a,"/",a)')  &
             trim(input_fields(input_field_num)%module_name), &
             trim(input_fields(input_field_num)%field_name)
        if(mpp_pe() .eq. mpp_root_pe()) &
             call error_mesg ('diag_manager opening_file', &
             'module/field_name '//trim(error_string)//&
             &' NOT registered, ALL fields should be registered BEFORE the first send_data', WARNING)  
        cycle
     endif
! Put the time axis in the axis field
     num_axes = output_fields(field_num)%num_axes
     axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
! make sure that axis_id are not -1
     do k = 1,num_axes
        if(axes(k)<0) then
           write(error_string,'(a)') output_fields(field_num)%output_name
           call error_mesg ('diag_manager opening_file','output_name '//trim(error_string)// &
                ' has axis_id = -1', FATAL)
        endif
     enddo
     axes(num_axes + 1) = files(i)%time_axis_id
     call write_axis_meta_data(files(i)%file_unit, axes(1:num_axes + 1), time_ops)
     if(time_ops) then
        axes(num_axes + 2) = files(i)%time_bounds_id
        call write_axis_meta_data(files(i)%file_unit, axes(1:num_axes + 2))     
     endif
  end do

! Looking for the first NON-static field in a file
  field_num1 = files(i)%fields(1)
  do j = 1, files(i)%num_fields
     field_num = files(i)%fields(j)
     if(output_fields(field_num)%time_ops) then
        field_num1 = field_num
        exit
     endif
  enddo
  do j = 1, files(i)%num_fields
     field_num = files(i)%fields(j)
     input_field_num = output_fields(field_num)%input_field
     if (.not.input_fields(input_field_num)%register) cycle
! Make sure that 1 file contains either time_average or instantaneous fields
! cannot have both time_average and instantaneous in 1 file
     if(.not. mix_snapshot_average_fields) then
        if((output_fields(field_num)%time_ops.neqv.output_fields(field_num1)%time_ops) .and. &
             .not.output_fields(field_num1)%static .and. .not.output_fields(field_num)%static) then
           if(mpp_pe() == mpp_root_pe()) call error_mesg ('diag_manager opening_file','file '// &
                trim(files(i)%name)//' can NOT have BOTH time average AND instantaneous fields.'// &
                ' Create a new file or set mix_snapshot_average_fields=.true.' , FATAL)
        endif
     endif    
! Put the time axis in the axis field
     num_axes = output_fields(field_num)%num_axes
     axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
     if (.not. output_fields(field_num)%static) then
        num_axes=num_axes+1
        axes(num_axes) = files(i)%time_axis_id
     endif
     if(output_fields(field_num)%time_average) then
        avg = avg_name
     else if(output_fields(field_num)%time_max) then
        avg = avg_name
     else if(output_fields(field_num)%time_min) then
        avg = avg_name
     else
        avg = " "
     end if
     if(input_fields(input_field_num)%missing_value_present) then
        output_fields(field_num)%f_type = write_field_meta_data(files(i)%file_unit, &
             output_fields(field_num)%output_name, axes(1:num_axes), &
             input_fields(input_field_num)%units, &
             input_fields(input_field_num)%long_name, &
             input_fields(input_field_num)%range, output_fields(field_num)%pack,&
             input_fields(input_field_num)%missing_value, avg_name = avg,&
             time_method=output_fields(field_num)%time_method,&
             standard_name = input_fields(input_field_num)%standard_name)
        
! NEED TO TAKE CARE OF TIME AVERAGING INFO TOO BOTH CASES
     else
        output_fields(field_num)%f_type = write_field_meta_data(files(i)%file_unit, &
             output_fields(field_num)%output_name, axes(1:num_axes), &
             input_fields(input_field_num)%units, &
             input_fields(input_field_num)%long_name, &
             input_fields(input_field_num)%range, output_fields(field_num)%pack,&
             avg_name = avg,&
             time_method=output_fields(field_num)%time_method, &
             standard_name = input_fields(input_field_num)%standard_name)
     endif
  end do

! If any of the fields in the file are time averaged, need to output the axes
! Use double precision since time axis is double precision
   if(time_ops) then
      time_axis_id(1) = files(i)%time_axis_id
      files(i)%f_avg_start = write_field_meta_data(files(i)%file_unit, &
           avg_name // '_T1', time_axis_id, time_units, &
           "Start time for average period", pack=1)
      files(i)%f_avg_end = write_field_meta_data(files(i)%file_unit, &
           avg_name // '_T2', time_axis_id, time_units, &
           "End time for average period", pack=1)
      files(i)%f_avg_nitems = write_field_meta_data(files(i)%file_unit, &
           avg_name // '_DT', time_axis_id,     &
           trim(time_unit_list(files(i)%time_units)), &
           "Length of average period", pack=1)
  endif

  if (time_ops) then
      time_axis_id(1) = files(i)%time_axis_id
      time_bounds_id(1) = files(i)%time_bounds_id
      call get_diag_axis( time_axis_id(1), time_name, time_units, time_longname, &
           &cart_name, dir, edges, Domain, data)
      call get_diag_axis( time_bounds_id(1), timeb_name, timeb_units, timeb_longname, &
           &cart_name, dir, edges, Domain, data)     
      files(i)%f_bounds =  write_field_meta_data(files(i)%file_unit, &
           trim(time_name)//'_bounds', (/time_bounds_id,time_axis_id/),     &
           trim(time_unit_list(files(i)%time_units)), &
           trim(time_name)//' axis boundaries', pack=1)      
   end if
! Let lower levels know that all meta data has been sent
   call done_meta_data(files(i)%file_unit)
 end subroutine opening_file

!--------------------------------------------------------------------------

subroutine diag_data_out(file, field, dat, time, final_call_in, static_write_in)

integer, intent(in)          :: file, field
real, intent(inout)          :: dat(:, :, :)
type(time_type), intent(in)  :: time
logical, optional, intent(in):: final_call_in, static_write_in
logical                      :: final_call, too_early, static_write
integer                      :: i, num
real :: dif, time_data(2, 1, 1), dt_time(1, 1, 1), start_dif, end_dif

too_early = .false.
final_call = .false.
if(present(final_call_in)) final_call = final_call_in
static_write = .false.
if(present(static_write_in)) static_write = static_write_in
dif = get_date_dif(time, base_time, files(file)%time_units)
! get file_unit, open new file and close curent file if necessary
If(.not. static_write .or. files(file)%file_unit<0) call check_and_open(file, time, too_early)
if(too_early) return  ! still too early to write data
call diag_field_out(files(file)%file_unit,output_fields(field)%f_type, dat, dif)
! record number of bytes written to this file
files(file)%bytes_written = files(file)%bytes_written + (size(dat,1)*size(dat,2)*size(dat,3))*(8/output_fields(field)%pack)

! *** inserted this line because start_dif < 0 for static fields ***
if (.not. output_fields(field)%static) then 
   start_dif = get_date_dif(output_fields(field)%last_output, base_time,files(file)%time_units)
   if(.not. mix_snapshot_average_fields) then
      end_dif = get_date_dif(output_fields(field)%next_output, base_time, files(file)%time_units)
   else
      end_dif = dif
   endif
endif

! Need to write average axes out;
do i = 1, files(file)%num_fields
   num = files(file)%fields(i)
   if(output_fields(num)%time_ops .and. &
      input_fields(output_fields(num)%input_field)%register) then
      if(num == field) then
! Output the axes if this is first time-averaged field
         time_data(1, 1, 1) = start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_start, &
            time_data(1:1,:,:), dif)
         time_data(2, 1, 1) = end_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_end, &
            time_data(2:2,:,:), dif)
! Compute the length of the average
         dt_time(1, 1, 1) = end_dif - start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_nitems, &
            dt_time(1:1,:,:), dif)
!
! Include boundary variable for CF compliance
!
         call diag_field_out(files(file)%file_unit, files(file)%f_bounds, &
            time_data(1:2,:,:), dif)         
         exit
      endif
   end if
end do

! If write time is greater (equal for the last call) than last_flush for this file, flush it
if(final_call) then
   if(time >= files(file)%last_flush) then
      call diag_flush(files(file)%file_unit)
      files(file)%last_flush = time
   endif
else
   if(time > files(file)%last_flush) then
      call diag_flush(files(file)%file_unit)
      files(file)%last_flush = time
   endif
endif

end subroutine diag_data_out

!-------------------------------------------------------------------------

subroutine check_and_open(file, time, too_early)
! this routine checks if it is time to open a new file. If yes, it first closes the
! current file, open a new file and returns file_unit
! previous diag_manager_end is replaced by closing_file and output_setup by opening_file.

integer, intent(in)         :: file
type(time_type), intent(in) :: time
logical, intent(out)        :: too_early

if(time>=files(file)%start_time) then 
   if(files(file)%file_unit < 0)then ! need to open a new file
      call opening_file(file, time)
   else
         if(time>files(file)%next_open) then ! need to close current file and open a new one 
            call write_static(file)  ! write all static fields and close this file
            call opening_file(file, time)
            files(file)%next_open = diag_time_inc(files(file)%next_open, files(file)%new_file_freq, &
                 files(file)%new_file_freq_units)
         endif ! no need to open new file, simply return file_unit
   endif
   too_early=.false.
else
   too_early = .true.
endif
end subroutine check_and_open
!-------------------------------------------------------------------------

function get_date_dif(t2, t1, units)

real                        :: get_date_dif
type(time_type), intent(in) :: t2, t1
integer, intent(in)         :: units
integer                     :: year, month, day, hour, minute, second, dif_seconds, dif_days
type(time_type)             :: dif_time

! Compute time axis label value
if(t2 < t1)   call error_mesg('get_date_dif', &
                't2 is less than t1', FATAL)

dif_time = t2 - t1
!del call get_date(dif_time, year, month, day, hour, minute, second) ! not used
call get_time(dif_time, dif_seconds, dif_days)

if(units == DIAG_SECONDS) then
   get_date_dif = dif_seconds + 86400 * dif_days
else if(units == DIAG_MINUTES) then
   get_date_dif = 1440 * dif_days + dif_seconds / 60.
else if(units == DIAG_HOURS) then
   get_date_dif = 24 * dif_days + dif_seconds / 3600.
else if(units == DIAG_DAYS) then
   get_date_dif = dif_days + dif_seconds / 86400.
else if(units == DIAG_MONTHS) then
   call error_mesg('diag_data_out', 'months not supported as output units', FATAL)
else if(units == DIAG_YEARS) then
call error_mesg('diag_data_out', 'years not supported as output units', FATAL)
else
   call error_mesg('diag_data_out', 'illegal time units', FATAL)
end if

end function get_date_dif

! ---------------------------------------------------------------------------

function get_time_string(filename, current_time)
! This function determines a string based on current time.
! This string is used as suffix in output file name

type(time_type), intent(in)    :: current_time
character(len=128), intent(in) :: filename
character(len=128)             :: get_time_string
integer                        :: yr1, mo1, dy1, hr1, mi1, sc1  ! get from current time
integer                        :: yr2, dy2, hr2, mi2            ! for computing next_level time unit
integer                        :: yr1_s, mo1_s, dy1_s, hr1_s, mi1_s, sc1_s ! actual values to write string
character(len=20)              :: yr, mo, dy, hr, mi, sc        ! string of current time (output)
integer                        :: abs_sec, abs_day              ! component of current_time
integer                        :: days_per_month(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
integer                        :: julian_day, i, position, len, first_percent
character(len=10)              :: format
character(len=1)               :: width  ! width of the field in format write
character(len=128)             :: filetail

if (get_calendar_type() == GREGORIAN) &
     call error_mesg('diag_manager, get_time_string', 'GREGORIAN is NOT supported calendar', WARNING)
julian_day = 0
format = '("_",i*.*)'
call get_date(current_time, yr1, mo1, dy1, hr1, mi1, sc1)
len = len_trim(filename)
first_percent = index(filename, '%')
filetail = filename(first_percent:len)
! compute year string 
position = INDEX(filetail, 'yr')
if(position>0) then
   width = filetail(position-1:position-1)
   yr1_s = yr1
   format(7:9) = width//'.'//width
   write(yr, format) yr1_s   
   yr2 = 0
else  
   yr = ' '
   yr2 = yr1 - 1
endif
! compute month string 
position = INDEX(filetail, 'mo')
if(position>0) then   
   width = filetail(position-1:position-1)
   mo1_s = yr2*12 + mo1  
   format(7:9) = width//'.'//width
   write(mo, format) mo1_s
else
   mo = ' '
endif
! compute day string        
if(LEN_TRIM(mo) > 0) then       !  month present
   dy1_s = dy1 
   dy2 = dy1_s - 1
elseif(LEN_TRIM(yr) >0 )  then  ! no month, year present
! compute julian day
   if(mo1 == 1) then
      dy1_s = dy1
   else
      do i = 1, mo1-1
         julian_day = julian_day + days_per_month(i)
      enddo
      if(leap_year(current_time) .and. mo1>2) julian_day = julian_day + 1
      julian_day = julian_day + dy1
      dy1_s = julian_day
   endif
   dy2 = dy1_s - 1
else                            ! no month, no year
   call get_time(current_time, abs_sec, abs_day)
   dy1_s = abs_day  
   dy2 = dy1_s 
endif
position = INDEX(filetail, 'dy')
if(position>0) then 
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(dy, format) dy1_s
else
   dy = ' '
endif
! compute hour string
if(LEN_TRIM(dy) > 0) then
   hr1_s = hr1
else
   hr1_s = dy2*24 + hr1
endif
hr2 = hr1_s
position = INDEX(filetail, 'hr')
if(position>0) then
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(hr, format) hr1_s
else
   hr = ' '
endif
! compute minute string
if(LEN_TRIM(hr) > 0) then
   mi1_s = mi1
else
   mi1_s = hr2*60 + mi1
endif
mi2 = mi1_s
position = INDEX(filetail, 'mi')
if(position>0) then
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(mi, format) mi1_s
else
   mi = ' '
endif
! compute second string
if(LEN_TRIM(mi) > 0) then
   sc1_s = sc1
else
   sc1_s = mi2*60 + sc1
endif
position = INDEX(filetail, 'sc')
if(position>0) then
   width = filetail(position-1:position-1)
   format(7:9) = width//'.'//width
   write(sc, format) sc1_s
else
   sc = ' '
endif
get_time_string = trim(yr)//trim(mo)//trim(dy)//trim(hr)//trim(mi)//trim(sc)
end function get_time_string
      
!----------------------------------------------------------------------
subroutine check_duplicate_output_fields()
! pair(output_name and output_file) should be unique in output_fields
integer :: i, j, tmp_file
character(len=128) :: tmp_name
! Do the checking when more than 1 output_fileds present
if(num_output_fields <= 1) return 
do i = 1, num_output_fields-1
   tmp_name = trim(output_fields(i)%output_name)
   tmp_file =  output_fields(i)%output_file
   do j = i+1, num_output_fields
      if((tmp_name == trim(output_fields(j)%output_name)).and. &
           (tmp_file == output_fields(j)%output_file)) &
           call error_mesg (' ERROR in diag_table', &           
           &' output_field '//tmp_name//' duplicated', FATAL)
   enddo
enddo
end subroutine check_duplicate_output_fields

!-------------------------------------------------------------------------

function diag_time_inc(time, output_freq, output_units)

type (time_type)             :: diag_time_inc
type (time_type), intent(in) :: time
integer, intent(in)          :: output_freq, output_units

! special values for output frequency are -1 for output at end of run
! and 0 for every timestep.  Need to check for these here?
! Return zero time increment, hopefully this value is never used

if (output_freq == END_OF_RUN .or. output_freq == EVERY_TIME) then
    diag_time_inc = time
    return
endif

! Make sure calendar was not set after initialization
if (base_year == 0 .and. get_calendar_type() /= NO_CALENDAR) &
  call error_mesg('diag_time_inc',  &
        'calendar_type was set after diag_manager_init', FATAL)
   

if(output_units == DIAG_SECONDS) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, output_freq, 0)
   else
      diag_time_inc = increment_date(time, 0, 0, 0, 0, 0, output_freq)
   endif
else if(output_units == DIAG_MINUTES) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, output_freq*60, 0)
   else
      diag_time_inc = increment_date(time, 0, 0, 0, 0, output_freq, 0)
   endif
else if(output_units == DIAG_HOURS) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, output_freq*3600, 0)
   else
      diag_time_inc = increment_date(time, 0, 0, 0, output_freq, 0, 0)
   endif
else if(output_units == DIAG_DAYS) then
   if (get_calendar_type() == NO_CALENDAR) then
      diag_time_inc = increment_time(time, 0, output_freq)
   else
      diag_time_inc = increment_date(time, 0, 0, output_freq, 0, 0, 0)
   endif
else if(output_units == DIAG_MONTHS) then
   if (get_calendar_type() == NO_CALENDAR) then
      call error_mesg('diag_time_inc', &
         'output units of months NOT allowed with no calendar', FATAL)
   else
      diag_time_inc = increment_date(time, 0, output_freq, 0, 0, 0, 0)
   endif
else if(output_units == DIAG_YEARS) then
   if (get_calendar_type() == NO_CALENDAR) then
      call error_mesg('diag_time_inc', &
         'output units of years NOT allowed with no calendar', FATAL)
   else
      diag_time_inc = increment_date(time, output_freq, 0, 0, 0, 0, 0)
   endif
else 
   call error_mesg('diag_time_inc','illegal output units',FATAL)
endif

end function diag_time_inc

!-------------------------------------------------------------------------

! <FUNCTION NAME="get_base_time">

!   <OVERVIEW>
!     Return base time for diagnostics. 
!   </OVERVIEW>
!   <DESCRIPTION>
!     Return base time for diagnostics (note: base time must be >= model time).
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_base_time ()
!   </TEMPLATE>

 function get_base_time ()
 type(time_type) :: get_base_time

   if (.not.module_is_initialized) call error_mesg (  &
                        'get_base_time in diag_manager_mod', &
                        'module has not been initialized', FATAL)

   get_base_time = base_time

 end function get_base_time
! </FUNCTION>

!-------------------------------------------------------------------------
! <SUBROUTINE NAME="get_base_date">
!   <OVERVIEW>
!     Return base date for diagnostics.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Return date information for diagnostic reference time.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call get_base_date (year, month, day, hour, minute, second)
!   </TEMPLATE>

 subroutine get_base_date (year, month, day, hour, minute, second)
   integer, intent(out) :: year, month, day, hour, minute, second

   if (.not.module_is_initialized) call error_mesg (  &
                        'get_base_date in diag_manager_mod', &
                        'module has not been initialized', FATAL)
   year   = base_year
   month  = base_month
   day    = base_day
   hour   = base_hour
   minute = base_minute
   second = base_second

 end subroutine get_base_date
! </SUBROUTINE>

!-------------------------------------------------------------------------


!-------------------------------------------------------------------------

! <FUNCTION NAME="need_data">

!   <OVERVIEW>
!     Determine whether data is needed for the current model time step.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Determine whether data is needed for the current model time step.
!     Since diagnostic data are buffered, the "next" model time is passed
!     instead of the current model time. This call can be used to minimize
!     overhead for complicated diagnostics.
!   </DESCRIPTION>
!   <TEMPLATE>
!     need_data(diag_field_id,next_model_time)
!   </TEMPLATE>

!   <IN NAME="inext_model_time" TYPE="time_type"  >
!     next_model_time = current model time + model time_step
!   </IN>
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>

function need_data(diag_field_id,next_model_time)
!
! next_model_time = current model time + model time_step
!
type (time_type), intent(in) :: next_model_time
integer, intent(in) :: diag_field_id
logical :: need_data
integer :: i, out_num 
! loop through output fields

need_data=.false.

! If diag_field_id is < 0 it means that this field is unused, return
if (diag_field_id < 0 ) return
do i= 1,input_fields(diag_field_id)%num_output_fields
! Get index to an output field
   out_num = input_fields(diag_field_id)%output_fields(i)
   if (.not.output_fields(out_num)%static) then
      if (next_model_time > output_fields(out_num)%next_output) need_data=.true.
! Is this output field being time averaged?
! assume average data based on every timestep
! needs to be changed when different forms of averaging are implemented 
      if (output_fields(out_num)%time_average) need_data=.true. 
   endif
enddo
return

end function need_data
! </FUNCTION>

!-------------------------------------------------------------------------
 
end module diag_manager_mod

! <INFO>

!   <COMPILER NAME="COMPILING AND LINKING SOURCE">
!     Any module or program unit using <TT>diag_manager_mod</TT> must contain the line

!   <PRE>
!   use diag_manager_mod
!   </PRE>

!   If netCDF output is desired, the cpp flag <TT>-Duse_netCDF</TT>
!   must be turned on. The loader step requires an explicit link to the
!   netCDF library (typically something like <TT>-L/usr/local/lib
!   -lnetcdf</TT>, depending on the path to the netCDF library).
!   <LINK SRC="http://www.unidata.ucar.edu/packages/netcdf/guidef">netCDF
!   release 3 for fortran</LINK> is required.
!   </COMPILER>
!   <PRECOMP FLAG="PORTABILITY"> 
!     <TT>diag_manager_mod</TT> uses standard f90.
!   </PRECOMP>
!   <LOADER FLAG="ACQUIRING SOURCE">
!     GFDL users can checkout diag_manager_mod using the cvs command
!     <TT>setenv CVSROOT '/home/fms/cvs';cvs co diag_manager</TT>.  
!   </LOADER>

! </INFO>
