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
!   of FMS code /shared/mpp/mpp_io.F90.
!   A single group of calls to the <TT>diag_manager_mod</TT> interfaces provides data to disk
!   at any number of sampling and/or averaging intervals specified at run-time.
!   Run-time specification of diagnostics are input through the diagnostics table. <BR/>
! <B> Usage</B> of <TT> diag_manager</TT> includes the following steps:<BR/>
! 1. Create diag_table as described in 
!   <LINK SRC="diag_table_tk.html">diag_table_tk</LINK> documentation.<BR/>
! 2. Call <TT>diag_manager_init</TT> to initialize diag_manager_mod <BR/>
! 3. Call <TT> register_diag_field</TT> to register the field to be output. <BR/>
!    <B> NOTE</B> ALL fields in diag_table should be registered BEFORE the first send_data call
! 4. Call <TT>send_data</TT> to send data to output fields <BR/>
! 5. Call <TT>diag_manager_end</TT> to exit diag_manager <BR/><BR/>
!   
! <B> Features </B> of <TT> diag_manager_mod </TT>: <BR/>
! 1. Ability to output from 0-D array (scalars) to 3-D arrays. <BR/>
! 2. Ability to output time average of fields that have time dependent mask. <BR/>
! 3. Give optional warning if <TT>register_diag_field </TT>fails due to misspelled module name
!    or field name. <BR/>
! 4. Check if a field is registered twice.<BR/>     
! 5. Check for duplicate lines in diag_table. <BR/>
! 6. <LINK SRC="diag_table_tk.html">diag_table</LINK> can contain fields that are NOT written to any files. 
!    The file name in diag_table of these fields is <TT> null</TT> <BR/>
! 7. By default, a field is output in its global grid, it is now possible to output a field in 
!    a region specified by user, see <TT>send_data</TT> for more details <BR/>
! 8. To check if the diag table is set up correctly, user should set <TT>debug_diag_manager=.true.</TT> in 
!    diag_manager namelist, then the the content of diag_table is printed in stdout. <BR/>
! 9. New optional format of file information in <LINK SRC="diag_table_tk.html">diag_table</LINK><BR/>
!    It is possible to have just one file name and reuse it many times. A time string will be appended
!    to the base file name each time a new file is opened. The time string can be any combination from
!    year to second of current model time. Here is an example of file information: <BR/>
!    <TT>"file2_yr_dy%1yr%3dy",2,"hours",1,"hours","Time", 10, "days", "1 1 7 0 0 0", 6, "hours"</TT>
!     <BR/>
!    From left to right we have: file name, output frequency, output frequency unit, Format (should always
!    be 1), time axis unit, time axis name, frequency for creating new file, unit for creating new file,
!    start time of the new file, file duration, file duration unit.<BR/>
!    file duration, if absent, will equal to frequency for creating new file <BR/>
!    the above file means: create a new file every 10 days, each file will last 6 hours from creation time,
!    no files will be created before time "1 1 7 0 0 0" <BR/>
!    In this example the string <TT>10, "days", "1 1 7 0 0 0", 6, "hours"</TT> is optional.<BR/>
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
!    new feature users should <B>separate</B> snapshot fields and time averaged fields in 
!    <B>different</B> files or a fatal error will occur.<BR/>
!    The namelist <B>default</B> value is <TT>mix_snapshot_average_fields=.false.</TT> <BR/><BR/>
! 11 Time average, Max and Min <BR/>
!    In addition to time average userscan also get Max or Min value during the same interval of time as time
!    average. For this purpose, in the diag table users must replace <TT>.true.</TT> or <TT>.false.</TT> by 
!    <TT>"max"</TT> or <TT> "min"</TT>. <BR/>
!    Currently, Max and Min are not available for regional output. <BR/><BR/>
! 12 Standard_name is added as optional in register_diag_field. <BR/><BR/>
! 13 When namelist debug_diag_manager = true array bounds are checked in send_data. <BR/><BR/>
! 14 Coordinate attribute can be written in output file if argument "aux" is given in diag_axis_init.
!    The corresponding fields (geolat/geolon) should also be written in the same file.<BR/><BR/>

!
! </DESCRIPTION>

use time_manager_mod, only: set_time, set_date, operator(>=), operator(>), operator(<), operator(==), &
                            operator(/=), time_type, month_name, get_calendar_type, NO_CALENDAR, &
                            operator(/), operator(+)                            
use mpp_io_mod, only      : mpp_open, MPP_RDONLY, MPP_ASCII, mpp_close
use fms_mod, only         : error_mesg, FATAL, WARNING, NOTE, close_file, stdlog, write_version_number,  &
                            file_exist, mpp_pe, open_namelist_file, check_nml_error, lowercase, stdout, &
                            mpp_error, fms_error_handler
use mpp_mod, only         : mpp_get_current_pelist, mpp_npes, mpp_sync, mpp_root_pe, mpp_sum    
use diag_axis_mod, only   : diag_axis_init, get_axis_length, max_axes                            
use diag_util_mod, only   : get_subfield_size, log_diag_field_info, update_bounds, check_out_of_bounds, &
                            check_bounds_are_exact_dynamic, check_bounds_are_exact_static, init_file, &
                            diag_time_inc, find_input_field, init_input_field, init_output_field, &
                            diag_data_out, write_static, check_duplicate_output_fields, get_date_dif

use diag_data_mod, only   : max_files, DIAG_OTHER, DIAG_OCEAN, DIAG_ALL, EVERY_TIME, END_OF_RUN, &
                            DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS, &
                            num_files, max_input_fields, max_output_fields, num_output_fields, EMPTY, &
                            FILL_VALUE, null_axis_id, MAX_VALUE, MIN_VALUE, single, base_time, base_year, &
                            base_month, base_day, base_hour, base_minute, base_second, global_descriptor, &
                            coord_type, files, input_fields, output_fields, Time_zero, append_pelist_name, &
                            mix_snapshot_average_fields, first_send_data_call, do_diag_field_log, &
                            write_bytes_in_file, debug_diag_manager, diag_log_unit, time_unit_list, &
                            pelist_name, max_axes, module_is_initialized, max_num_axis_sets

use diag_output_mod, only : get_diag_global_att, set_diag_global_att

use constants_mod, only   : SECONDS_PER_HOUR, SECONDS_PER_MINUTE

implicit none
private

public  diag_manager_init, send_data, send_tile_averaged_data, diag_manager_end,  &
        register_diag_field, register_static_field, diag_axis_init, get_base_time, &
        get_base_date, need_data, average_tiles, DIAG_ALL, DIAG_OCEAN, DIAG_OTHER, &
        get_date_dif, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, &
        DIAG_YEARS, get_diag_global_att, set_diag_global_att


! version number of this module
character(len=128)  :: version = '$Id: diag_manager.F90,v 15.0.4.1 2007/12/04 17:12:54 slm Exp $'
character(len=128)  :: tagname = '$Name: omsk_2008_03 $'  


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
   long_name, units, missing_value, range, do_not_log, err_msg)

! Indicates the calling modules intent to supply data for this field.
integer                                :: register_diag_field_scalar
character(len=*), intent(in)           :: module_name, field_name
type(time_type),  optional, intent(in) :: init_time
character(len=*), optional, intent(in) :: long_name, units
real, optional, intent(in)             :: missing_value, range(2)
logical,          optional, intent(in) :: do_not_log ! if TRUE, field information is not logged
character(len=*), optional, intent(out):: err_msg
 
if(present(init_time)) then
   register_diag_field_scalar = register_diag_field_array(module_name, field_name,&
        (/null_axis_id/), init_time,long_name, units, missing_value, range, &
        do_not_log=do_not_log, err_msg=err_msg)
else
   if(present(err_msg)) err_msg = ''
   register_diag_field_scalar = register_static_field(module_name, field_name,&
        (/null_axis_id/),long_name, units, missing_value, range, do_not_log=do_not_log)
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

function register_diag_field_array(module_name, field_name, axes, init_time, &
   long_name, units, missing_value, range, mask_variant,standard_name,verbose,&
   do_not_log,err_msg)

! Indicates the calling modules intent to supply data for this field.

integer                                :: register_diag_field_array
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
type(time_type), intent(in)            :: init_time
character(len=*), optional, intent(in) :: long_name, units, standard_name
real, optional, intent(in)             :: missing_value, range(2)
logical, optional, intent(in)          :: mask_variant,verbose
logical, optional, intent(in)          :: do_not_log ! if TRUE, field info is not logged
character(len=*), optional, intent(out):: err_msg
integer                                :: field, j, ind, file_num, freq
integer                                :: output_units
logical                                :: mask_variant1, verbose1
character(len=128)                     :: msg
mask_variant1 = .false.
verbose1 = .false.
if(present(mask_variant)) mask_variant1 = mask_variant
if(present(verbose)) verbose1 = verbose
if(present(err_msg)) err_msg = ''

! Call register static, then set static back to false

register_diag_field_array = register_static_field(module_name, field_name, axes, &
   long_name, units, missing_value, range, mask_variant1, dynamic =.true., &
   do_not_log=do_not_log)

if((debug_diag_manager.or.verbose1) .and. register_diag_field_array<0 ) &
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
! Need to sync  close_time too
         files(file_num)%close_time = diag_time_inc(files(file_num)%start_time,  &
                                                    files(file_num)%duration,    &
                                                    files(file_num)%duration_units)
!-- NZ
!-- Why is this block of code here in this function? 
!-- This kind of resetting the member variables should be done by one function call
!-- in the module that first sets them (diag_util_mod).

      endif
! Need to increase next_open until it is greater than init time
      do
         if(files(file_num)%next_open > init_time) exit
         files(file_num)%next_open = diag_time_inc(files(file_num)%next_open, &
              files(file_num)%new_file_freq, files(file_num)%new_file_freq_units, err_msg=msg)
         if(msg /= '') then
           if(fms_error_handler('register_diag_field',' file='//trim(files(file_num)%name)//': '//trim(msg),err_msg)) return
         endif
      enddo
      freq = files(file_num)%output_freq
      output_units = files(file_num)%output_units
      output_fields(ind)%next_output = diag_time_inc(init_time, freq, output_units, err_msg=msg)
      if(msg /= '') then
        if(fms_error_handler('register_diag_field',' file='//trim(files(file_num)%name)//': '//trim(msg),err_msg)) return
      endif
      output_fields(ind)%next_next_output = &
         diag_time_inc(output_fields(ind)%next_output, freq, output_units, err_msg=msg)
      if(msg /= '') then
        if(fms_error_handler('register_diag_field',' file='//trim(files(file_num)%name)//': '//trim(msg),err_msg)) return
      endif
      if(debug_diag_manager .and. mpp_pe() == mpp_root_pe() .and. output_fields(ind)%local_output ) then
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
   long_name, units, missing_value, range, mask_variant, require, standard_name, dynamic, &
   do_not_log)

integer                                :: register_static_field
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
character(len=*), optional, intent(in) :: long_name, units, standard_name
real, optional, intent(in)             :: missing_value, range(2)
logical, optional, intent(in)          :: mask_variant
logical, optional, intent(in)          :: require  ! require static field to be in every file, e.g. 2-d axes
logical, optional, intent(in)          :: dynamic
logical, optional, intent(in)          :: do_not_log ! if TRUE, field information is not logged
integer                                :: field, num_axes, j, out_num, siz(3), local_siz(3), k
logical                                :: mask_variant1, dynamic1, allow_log
integer                                :: local_start(3), local_end(3) ! indices of local domain of global axes
character(len=128)                     :: msg

mask_variant1 = .false.
if(present(mask_variant)) mask_variant1 = mask_variant
dynamic1 = .false.
if(present(dynamic)) dynamic1 = dynamic
if (.not.module_is_initialized) call error_mesg ('register_static_field',  &
                       'diag_manager has NOT been initialized', FATAL)

allow_log = .TRUE. ; if (present(do_not_log)) allow_log = .not.do_not_log
if ( do_diag_field_log.and.allow_log ) then
   call log_diag_field_info (module_name, field_name, axes, &
        long_name, units, missing_value=missing_value, range=range, &
        dynamic=dynamic1)
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
            output_fields(out_num)%buffer = MAX_VALUE
         else if(output_fields(out_num)%time_min) then
            output_fields(out_num)%buffer = MIN_VALUE
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
         output_fields(out_num)%buffer = MAX_VALUE
      else if(output_fields(out_num)%time_min) then
         output_fields(out_num)%buffer = MIN_VALUE
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

! Initialize a time variable used in an error check
   output_fields(out_num)%Time_of_prev_field_data = Time_zero
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

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_0d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real"  > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_0d(diag_field_id, field, time, err_msg)
    
logical                      :: send_data_0d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field
type (time_type), intent(in),  optional :: time
character(len=*), intent(out), optional :: err_msg
real :: field_out(1, 1, 1)

! First copy the data to a three d array with last element 1
field_out(1, 1, 1) = field
send_data_0d = send_data_3d(diag_field_id, field_out, time, err_msg=err_msg)
end function send_data_0d

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_1d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg)

logical                      :: send_data_1d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:)
real, intent(in), optional   :: weight
type (time_type), intent(in), optional :: time
integer, intent(in), optional          :: is_in, ie_in
logical, intent(in), optional          :: mask(:)
real,    intent(in), optional          :: rmask(:)
character(len=*), intent(out), optional :: err_msg
real    :: field_out(size(field(:)), 1, 1)
logical ::  mask_out(size(field(:)), 1, 1)


! First copy the data to a three d array with last element 1
field_out(:, 1, 1) = field

! Default values for mask
mask_out = .true.
if(present(mask)) mask_out(:, 1, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, 1, 1) = .false.
if(present(mask) .or. present(rmask)) then
   send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in, 1, 1, mask_out, &
        ie_in=ie_in,weight=weight, err_msg=err_msg)
else
   send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in, 1, 1, ie_in=ie_in, weight=weight, err_msg=err_msg)
endif
end function send_data_1d


!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_2d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_2d(diag_field_id, field, time, is_in, js_in, &
    mask, rmask, ie_in, je_in, weight, err_msg)

logical                      :: send_data_2d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:, :)
real, intent(in), optional   :: weight
type (time_type), intent(in),  optional :: time
integer, intent(in), optional           :: is_in, js_in, ie_in, je_in
logical, intent(in), optional           :: mask(:, :)
real,    intent(in), optional           :: rmask(:, :)
character(len=*), intent(out), optional :: err_msg
real    :: field_out(size(field, 1), size(field, 2), 1)
logical ::  mask_out(size(field, 1), size(field, 2), 1)


! First copy the data to a three d array with last element 1
field_out(:, :, 1) = field

! Default values for mask
mask_out = .true.
if(present(mask)) mask_out(:, :, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, :, 1) = .false.
if(present(mask) .or. present(rmask)) then
   send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in, js_in, 1, mask_out,&
        ie_in=ie_in, je_in=je_in,weight=weight, err_msg=err_msg)
else
   send_data_2d = send_data_3d(diag_field_id, field_out,time,is_in,js_in,1,ie_in=ie_in,je_in=je_in,weight=weight, err_msg=err_msg)
endif

end function send_data_2d

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_3d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             mask, rmask, ie_in,je_in, ke_in,weight, err_msg)

logical                      :: send_data_3d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:,:,:)
real, intent(in), optional   :: weight
type (time_type), intent(in), optional :: time
integer, intent(in),optional :: is_in, js_in, ks_in,ie_in,je_in, ke_in 
logical, intent(in),optional :: mask(:, :, :)
real,    intent(in),optional :: rmask(:, :, :)
character(len=*), intent(out), optional :: err_msg
character(len=256)           :: err_msg_local
real                         :: num, weight1
logical                      :: average, need_compute, local_output, phys_window
logical                      :: time_max, time_min
integer                      :: i, out_num, file_num, n1, n2, n3, number_of_outputs, ii,f1,f2,f3,f4
integer                      :: freq, units, is, js, ks, ie, je, ke, i1, j1,k1, j, k
character(len=128)           :: error_string
integer                      :: l_start(3), l_end(3) ! local start and end indices on 3 axes for regional output
integer                      :: hi, hj, twohi, twohj  ! halo size in x and y direction
integer                      :: b1,b2,b3 ! size of buffer in x,y,z axes
type (time_type)             :: middle_time      
logical                      :: missvalue_present
real                         :: missvalue

if(present(err_msg)) err_msg = ''
if (.not.module_is_initialized) then
  if(fms_error_handler('send_data_3d','diag_manager NOT initialized',err_msg)) return
endif
err_msg_local = ''

! send_data works in either one or another of two modes.
! 1. Input field is a window (e.g. FMS physics)
! 2. Input field includes halo data
! It cannot handle a window of data that has halos.
! (A field with no windows or halos can be thought of as a special case of either mode.)
! The logic for indexing is quite different for these two modes, but is not clearly separated.
! If both the beggining and ending indices are present, then field is assumed to have halos.
! If only beggining indices are present, then field is assumed to be a window.

! There are a number of ways a user could mess up this logic, depending on the combination
! of presence/absence of is,ie,js,je. The checks below should catch improper combinations.
if(present(ie_in)) then
  if(.not.present(is_in)) then
    if(fms_error_handler('send_data_3d','ie_in present without is_in',err_msg)) return
  endif
  if(present(js_in) .and. .not.present(je_in)) then
    if(fms_error_handler('send_data_3d','is_in and ie_in present, but js_in present without je_in',err_msg)) return
  endif
endif
if(present(je_in)) then
  if(.not.present(js_in)) then
    if(fms_error_handler('send_data_3d','je_in present without js_in',err_msg)) return
  endif
  if(present(is_in) .and. .not.present(ie_in)) then
    if(fms_error_handler('send_data_3d','js_in and je_in present, but is_in present without ie_in',err_msg)) return
  endif
endif

! If is, js, or ks not present default them to 1
is = 1; js = 1; ks = 1
if(present(is_in)) is = is_in
if(present(js_in)) js = js_in
if(present(ks_in)) ks = ks_in
n1 = size(field, 1); n2 = size(field, 2); n3 = size(field, 3)
ie = is+n1-1; je = js+n2-1; ke = ks+n3-1
if (present(ie_in)) ie = ie_in
if (present(je_in)) je = je_in
if (present(ke_in)) ke = ke_in
twohi = n1-(ie-is+1)
if(mod(twohi,2) /= 0) then
  if(fms_error_handler('send_data_3d','non-symmetric halos in first dimension',err_msg)) return
endif
twohj = n2-(je-js+1)
if(mod(twohj,2) /= 0) then
  if(fms_error_handler('send_data_3d','non-symmetric halos in second dimension',err_msg)) return
endif
hi = twohi/2; hj = twohj/2

! The next line is necessary to ensure that is,ie,js,ie are relative to field(1:,1:)
! But this works only when there is no windowing.
if(present(ie_in) .and. present(je_in)) then
  is=1+hi; ie=n1-hi; js=1+hj; je=n2-hj
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

missvalue_present = input_fields(diag_field_id)%missing_value_present
if(missvalue_present) missvalue = input_fields(diag_field_id)%missing_value

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
          if(fms_error_handler('send_data_3d','module/output_field '//trim(error_string)//&
             ', time must be present when output frequency = EVERY_TIME', err_msg)) return
       endif
     endif
   endif
   if(.not.output_fields(out_num)%static .and. .not.present(time)) then
         write (error_string,'(a,"/",a)')  &
             trim(input_fields(diag_field_id)%module_name), &
             trim(output_fields(out_num)%output_name)
          if(fms_error_handler('send_data_3d','module/output_field '//trim(error_string)//&
             ', time must be present for nonstatic field', err_msg)) return
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
               if(fms_error_handler('send_data_3d','module/output_field '//trim(error_string)//&
                    ' is skipped one time level in output data', err_msg)) return
            endif
         end if
! If average get size: Average intervals are last_output, next_output
         if(average) then
            b1=size(output_fields(out_num)%buffer,1); b2=size(output_fields(out_num)%buffer,2) 
            b3=size(output_fields(out_num)%buffer,3)
            if (input_fields(diag_field_id)%mask_variant) then           
               do k=1,b3; do j=1,b2;  do i=1,b1 
                  if(output_fields(out_num)%counter(i,j,k)>0.)then
                     output_fields(out_num)%buffer(i,j,k) = &
                          output_fields(out_num)%buffer(i,j,k)/output_fields(out_num)%counter(i,j,k)
                  else
                     output_fields(out_num)%buffer(i,j,k) =  missvalue
                  endif
               enddo;enddo;enddo                               
            else  !not mask variant
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
                  if (missvalue_present) then
                     do k=1,b3; do j=1,b2;  do i=1,b1
                        if(output_fields(out_num)%buffer(i,j,k)/= missvalue) &
                             output_fields(out_num)%buffer(i,j,k) = output_fields(out_num)%buffer(i,j,k)/num  
                     enddo;enddo;enddo
                  else
                     output_fields(out_num)%buffer = output_fields(out_num)%buffer/num
                  endif
               else
                  if (missvalue_present) then
                     if(any(output_fields(out_num)%buffer /= missvalue)) then
                        write (error_string,'(a,"/",a)')  &
                             trim(input_fields(diag_field_id)%module_name), &
                             trim(output_fields(out_num)%output_name)
                        if(mpp_pe() .eq. mpp_root_pe()) then
                          if(fms_error_handler('send_data_3d','module/output_field '//trim(error_string)//&
                             &', write EMPTY buffer', err_msg)) return
                        endif
                     endif
                  endif
               endif
            endif ! mask_variant
         elseif(time_min .or. time_max) then
            if(missvalue_present) then
               where(abs(output_fields(out_num)%buffer) == MIN_VALUE) 
                  output_fields(out_num)%buffer = missvalue
               end where
            endif ! if missvalue is NOT present buffer retains max_value or min_value
         endif !average
! Output field
         if((output_fields(out_num)%time_ops).and.(.not. mix_snapshot_average_fields)) then
            middle_time = (output_fields(out_num)%last_output+output_fields(out_num)%next_output)/2
            call diag_data_out(file_num, out_num, output_fields(out_num)%buffer, middle_time)
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
            output_fields(out_num)%buffer = MAX_VALUE
         else if (time_min) then
            output_fields(out_num)%buffer = MIN_VALUE
         else
            output_fields(out_num)%buffer = EMPTY
         endif
         if(input_fields(diag_field_id)%mask_variant .and. average) output_fields(out_num)%counter = 0.0
      end if  !time > output_fields(out_num)%next_output
   end if  !.not.output_fields(out_num)%static .and. freq /= END_OF_RUN
! Finished output of previously buffered data, now deal with buffering new data   

   if(.not.output_fields(out_num)%static .and. .not.need_compute .and. debug_diag_manager ) then
      call check_bounds_are_exact_dynamic(out_num, diag_field_id, Time, err_msg=err_msg_local)
      if(err_msg_local /= '') then
        if(fms_error_handler('send_data',err_msg_local,err_msg)) return
      endif
   endif
 
! Take care of submitted field data
   if(average) then
      if (input_fields(diag_field_id)%mask_variant) then
         if(need_compute) then
            write (error_string,'(a,"/",a)')  &
                 trim(input_fields(diag_field_id)%module_name), &
                 trim(output_fields(out_num)%output_name)   
            if(fms_error_handler('send_data_3d','module/output_field '//trim(error_string)//&
                 &', regional output NOT supported with mask_variant', err_msg)) return
         endif
         if (present(mask)) then
            if (missvalue_present) then              
               if(debug_diag_manager) then
                  call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                  call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                  if(err_msg_local /= '') then
                    if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
                  endif
               endif
!               do i=is,ie; do j=js,je; do k=ks,ke
               do k=ks,ke; do j=js,je;
!CDIR NODEP
                                      do i=is,ie
                  if(mask(i-is+1+hi,j-js+1+hj,k)) then
                     output_fields(out_num)%buffer(i-hi,j-hj,k)=output_fields(out_num)%buffer(i-hi,j-hj,k) + &
                          field(i-is+1+hi,j-js+1+hj,k)*weight1  
                     output_fields(out_num)%counter(i-hi,j-hj,k)=output_fields(out_num)%counter(i-hi,j-hj,k) + weight1                 
                  endif
               enddo; enddo; enddo

            else
               write (error_string,'(a,"/",a)')  &
                    trim(input_fields(diag_field_id)%module_name), &
                    trim(output_fields(out_num)%output_name)
               if(fms_error_handler('send_data_3d','module/output_field '//trim(error_string)//&
                    &', variable mask but no missing value defined', err_msg)) return
            endif
         else  ! no mask present
            write (error_string,'(a,"/",a)')  &
                 trim(input_fields(diag_field_id)%module_name), &
                 trim(output_fields(out_num)%output_name)
            if(fms_error_handler('send_data_3d','module/output_field '//trim(error_string)//&
                 &', variable mask but no mask given', err_msg)) return
         endif         
      else ! mask_variant=false
         if (present(mask)) then
            if (missvalue_present) then
               if(need_compute) then
                  do k = l_start(3),l_end(3)
                     k1=k-l_start(3)+1
                     do j = js,je 
                        do i = is, ie
                           if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                              i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
!CDIR NODEP                          
                              if(mask(i-is+1+hi,j-js+1+hj,k)) then
                                 output_fields(out_num)%buffer(i1,j1,k1) = output_fields(out_num)%buffer(i1,j1,k1)+ &
                                      field(i-is+1+hi,j-js+1+hj,k)*weight1                              
                              else
                                 output_fields(out_num)%buffer(i1,j1,k1) = missvalue   
                              endif
                              output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + 1
!                              output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + &
!                                                                  l_end(3)-l_start(3)+1
                           endif
                        enddo
                     enddo
                  enddo
               else
                  if(debug_diag_manager) then
                     call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     if(err_msg_local /= '') then
                       if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
                     endif
                  endif
!                  do i=is,ie; do j=js,je; do k=ks,ke
                  do k=ks,ke; do j=js,je;
!CDIR NODEP
                                         do i=is,ie
                     if(mask(i-is+1+hi,j-js+1+hj,k)) then
                        output_fields(out_num)%buffer(i-hi,j-hj,k)=output_fields(out_num)%buffer(i-hi,j-hj,k)+&
                             field(i-is+1+hi,j-js+1+hj,k)*weight1  
                     else
                        output_fields(out_num)%buffer(i-hi,j-hj,k)= missvalue
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
               if(.not.all(mask(f1:f2,f3:f4,ks:ke)) .and. mpp_pe() .eq. mpp_root_pe()) then
                  call error_mesg('warning2 send_data_3d', &
                  'Mask will be ignored since missing values were not specified',WARNING)
               endif
               if(need_compute) then                 
                  do j = js,je 
                     do i = is, ie
                        if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                           i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                           output_fields(out_num)%buffer(i1,j1,:)= output_fields(out_num)%buffer(i1,j1,:)+ &
                                field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           output_fields(out_num)%num_elements=output_fields(out_num)%num_elements+l_end(3)-l_start(3)+1
                        endif
                     enddo
                  enddo
               else                          
                  if(debug_diag_manager) then
                     call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     if(err_msg_local /= '') then
                       if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
                     endif
                  endif
                  output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = &
                       output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) &
                       + field(f1:f2,f3:f4,ks:ke)*weight1   
               endif
               if(.not. phys_window) output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1
            endif
         else ! mask NOT present
            if (missvalue_present) then
               if(need_compute) then 
                  do k = l_start(3),l_end(3)
                     k1=k-l_start(3)+1                    
                     do j = js,je 
                        do i = is, ie
                           if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                              i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1 
!CDIR NODEP                          
                              if(field(i-is+1,j-js+1,k) /= missvalue) then
                                 output_fields(out_num)%buffer(i1,j1,k1)= output_fields(out_num)%buffer(i1,j1,k1)+ &
                                   field(i-is+1+hi,j-js+1+hj,k)*weight1
                              else
                                 output_fields(out_num)%buffer(i1,j1,k1) = missvalue
                              endif
                              output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + 1
!                              output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + &
!                                                                    l_end(3)-l_start(3)+1
                           endif
                        enddo
                     enddo
                  enddo
                  if(.not. phys_window) then
                     if(any(field(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3)) /= &
                          missvalue)) &
                          output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1 
                  endif
               else
                  if(debug_diag_manager) then
                     call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     if(err_msg_local /= '') then
                       if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
                     endif
                  endif
!                  do i=is,ie; do j=js,je; do k=ks,ke
                  do k=ks,ke; do j=js,je;
!CDIR NODEP
                                         do i=is,ie
                     if(field(i-is+1+hi,j-js+1+hj,k) /= missvalue)  then
                        output_fields(out_num)%buffer(i-hi,j-hj,k)=output_fields(out_num)%buffer(i-hi,j-hj,k) + &
                             field(i-is+1+hi,j-js+1+hj,k)*weight1  
                     else
                        output_fields(out_num)%buffer(i-hi,j-hj,k)= missvalue
                     endif
                  enddo; enddo; enddo                  

                  if(any(field(f1:f2,f3:f4,ks:ke) /= missvalue)) &
                       output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + weight1    
               endif
            else       ! no missing value defined, No mask
               if(need_compute) then
                  do j = js,je  
                     do i = is, ie                                         
                        if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                           i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1                           
                           output_fields(out_num)%buffer(i1,j1,:)= output_fields(out_num)%buffer(i1,j1,:)+ &
                                field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                           output_fields(out_num)%num_elements=output_fields(out_num)%num_elements+l_end(3)-l_start(3)+1
                        endif
                     enddo
                  enddo
               else 
                  if(debug_diag_manager) then
                     call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                     call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                     if(err_msg_local /= '') then
                       if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
                     endif
                  endif
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
         if (need_compute) then
            do k = l_start(3),l_end(3)
               k1=k-l_start(3)+1
               do j = js,je
                  do i = is, ie
                     if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                        i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                        if(mask(i-is+1+hi,j-js+1+hj,k).and. &
                           field(i-is+1+hi,j-js+1+hj,k)>output_fields(out_num)%buffer(i1,j1,k1)) then
                           output_fields(out_num)%buffer(i1,j1,k1) = field(i-is+1+hi,j-js+1+hj,k)
                        endif
                     endif
                  enddo
               enddo
            enddo
         else
            if(debug_diag_manager)then
               call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
               call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
               if(err_msg_local /= '') then
                  if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
               endif
            endif
            where(mask(f1:f2,f3:f4,ks:ke).and. &
                  field(f1:f2,f3:f4,ks:ke)>output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke))&
                 output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
         endif
      else
         if (need_compute) then
            do k = l_start(3),l_end(3)
               k1=k-l_start(3)+1
               do j = js,je
                  do i = is, ie
                     if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                        i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                        if(field(i-is+1+hi,j-js+1+hj,k) > output_fields(out_num)%buffer(i1,j1,k1)) then
                           output_fields(out_num)%buffer(i1,j1,k1) = field(i-is+1+hi,j-js+1+hj,k)
                        endif
                     endif
                  enddo
               enddo
            enddo
         else
            if(debug_diag_manager)then
               call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
               call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
               if(err_msg_local /= '') then
                  if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
               endif
            endif            
            where (field(f1:f2,f3:f4,ks:ke) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke)) &
                 output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
         endif
      endif
      output_fields(out_num)%count_0d = 1
   else if (time_min) then
      if (present(mask)) then
         if (need_compute) then
            do k = l_start(3),l_end(3)
               k1=k-l_start(3)+1
               do j = js,je
                  do i = is, ie
                     if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                        i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                        if(mask(i-is+1+hi,j-js+1+hj,k).and. &
                           field(i-is+1+hi,j-js+1+hj,k)<output_fields(out_num)%buffer(i1,j1,k1)) then
                           output_fields(out_num)%buffer(i1,j1,k1) = field(i-is+1+hi,j-js+1+hj,k)
                        endif
                     endif
                  enddo
               enddo
            enddo
         else
            if(debug_diag_manager) then
               call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
               call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
               if(err_msg_local /= '') then
                  if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
               endif
            endif
            where(mask(f1:f2,f3:f4,ks:ke).and. &
                  field(f1:f2,f3:f4,ks:ke)<output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke)) &
                 output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke) 
         endif
      else
         if (need_compute) then
            do k = l_start(3),l_end(3)
               k1=k-l_start(3)+1
               do j = js,je
                  do i = is, ie
                     if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                        i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                        if(field(i-is+1+hi,j-js+1+hj,k) < output_fields(out_num)%buffer(i1,j1,k1)) then
                           output_fields(out_num)%buffer(i1,j1,k1) = field(i-is+1+hi,j-js+1+hj,k)
                        endif
                     endif
                  enddo
               enddo
            enddo
         else
            if(debug_diag_manager) then
               call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
               call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
               if(err_msg_local /= '') then
                  if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
               endif
            endif
            where (field(f1:f2,f3:f4,ks:ke) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke)) &
                 output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
         endif
      endif
      output_fields(out_num)%count_0d = 1
   else  ! ( not average, not min, max)
      output_fields(out_num)%count_0d = 1
      if(need_compute) then
         do j = js,je 
            do i = is, ie           
               if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                  i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1
                  output_fields(out_num)%buffer(i1,j1,:) = field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))
               endif
            enddo
         enddo
      else
         if(debug_diag_manager) then
            call update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
            call check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
            if(err_msg_local /= '') then
              if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
            endif
         endif
         output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke) = field(f1:f2,f3:f4,ks:ke)
      endif
               
      if (present(mask) .and. missvalue_present) then
         if(need_compute) then
            do k = l_start(3),l_end(3)
               k1=k-l_start(3)+1
               do j = js,je ;do i = is, ie                                                       
                  if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                     i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1                    
                     if(.not. mask(i-is+1+hi,j-js+1+hj,k)) &
                          output_fields(out_num)%buffer(i1,j1,k1) = missvalue                     
                  endif
               enddo;  enddo              
            enddo
         else
!            do i=is,ie; do j=js,je; do k=ks,ke
            do k=ks,ke; do j=js,je; do i=is,ie
               if(.not. mask(i-is+1+hi,j-js+1+hj,k)) &
                    output_fields(out_num)%buffer(i-hi,j-hj,k)= missvalue
            enddo; enddo; enddo                        
         endif
      endif
   endif !average

   if(output_fields(out_num)%static .and. .not.need_compute .and. debug_diag_manager) then
      call check_bounds_are_exact_static(out_num, diag_field_id, err_msg=err_msg_local)
      if(err_msg_local /= '') then
        if(fms_error_handler('send_data in diag_manager_mod',err_msg_local,err_msg)) return
      endif
   endif
 
! If rmask and missing value present, then insert missing value     
   if (present(rmask) .and. missvalue_present) then
      if(need_compute) then
         do k = l_start(3),l_end(3)
            k1=k-l_start(3)+1
            do j = js,je ;  do i = is, ie                                               
               if(l_start(1)+hi<=i.and.i<=l_end(1)+hi.and.l_start(2)+hj<=j.and.j<=l_end(2)+hj) then
                  i1 = i-l_start(1)-hi+1 ; j1=  j-l_start(2)-hj+1                 
                  if(rmask(i-is+1+hi,j-js+1+hj,k)<0.5) &
                       output_fields(out_num)%buffer(i1,j1,k1) = missvalue                                   
               endif
            enddo; enddo         
         enddo
      else
!         do i=is,ie; do j=js,je; do k=ks,ke
         do k=ks,ke; do j=js,je; do i=is,ie
            if( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5) &
                 output_fields(out_num)%buffer(i-hi,j-hj,k)= missvalue
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
real                        :: num 
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
            where (output_fields(i)%counter>0.) 
               output_fields(i)%buffer = output_fields(i)%buffer/output_fields(i)%counter
            elsewhere
               output_fields(i)%buffer = input_fields(input_num)%missing_value                            
            end where
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
      elseif( output_fields(i)%time_max .or. output_fields(i)%time_min) then
         if(input_fields(input_num)%missing_value_present) then
            where(abs(output_fields(i)%buffer) == MIN_VALUE) 
               output_fields(i)%buffer = input_fields(input_num)%missing_value
            end where
         endif ! if missvalue is NOT present buffer retains max_value or min_value 
      endif
! Output field
      if (freq == END_OF_RUN) output_fields(i)%next_output = time
      if( (output_fields(i)%time_ops).and.(.not. mix_snapshot_average_fields)) then
         middle_time = (output_fields(i)%last_output+output_fields(i)%next_output)/2
         call diag_data_out(file, i, output_fields(i)%buffer, middle_time, .true.)
      else
         call diag_data_out(file, i, output_fields(i)%buffer, output_fields(i)%next_output, .true.)
      endif
   elseif (.not. output_fields(i)%written_once) then
      call error_mesg('Potential error in diag_manager_end ',trim(output_fields(i)%output_name)//' NOT available,'//&
           ' check if output interval > runlength. Netcdf fill_values are written', NOTE)
      output_fields(i)%buffer = FILL_VALUE
      call diag_data_out(file, i, output_fields(i)%buffer, time, .true.)   
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

subroutine diag_manager_init(diag_model_subset, err_msg)
integer, optional, intent(IN) :: diag_model_subset
character(len=*), intent(out), optional :: err_msg
integer                       :: diag_subset_output
character(len=256) :: err_msg_local

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
   integer            :: file_duration
   character(len=25)  :: file_duration_units
end type tableA_type

character(len=256) :: record
character(len=128) :: record_tail
character(len=9)   :: amonth
integer            :: iunit,time_units, output_freq_units, nfiles,nfields
integer            :: j, log_unit, name_len, nrecs, io_status,file_duration_units
integer            :: file_duration
integer, allocatable, dimension(:) :: pelist
integer            :: record_len, record_tail_len,time_pos
logical            :: file_freq_present, file_start_present, file_duration_present
integer            :: yr, mo, dy, hr, mi, sc, new_file_freq_units
type(time_type)    :: start_time
real(kind=single)  :: foo

namelist /diag_manager_nml/ append_pelist_name, mix_snapshot_average_fields, max_output_fields, &
     max_input_fields, max_axes, do_diag_field_log, write_bytes_in_file, debug_diag_manager, max_num_axis_sets

type(tableB_type)  :: textB
type(tableA_type ) :: textA
type(coord_type)   :: local_coord !local coordinates used in local output

if (present(err_msg)) err_msg=''
!  If the module was already initialized do nothing
if (module_is_initialized) return
min_value = huge(foo)
max_value = -min_value

Time_zero = set_time(0,0)
diag_subset_output = DIAG_ALL
if (PRESENT(diag_model_subset))then
  if(diag_model_subset>=DIAG_OTHER .AND. diag_model_subset<=DIAG_ALL) then
    diag_subset_output = diag_model_subset
  else
    if(fms_error_handler('diag_manager_init','invalid value of diag_model_subset',err_msg)) return
  endif
endif
call mpp_open(iunit, 'input.nml',form=MPP_ASCII,action=MPP_RDONLY)
read(iunit,diag_manager_nml,iostat=io_status)
write(stdlog(), diag_manager_nml)
if (io_status > 0) then
   if(fms_error_handler('diag_manager_init', 'Error reading diag_manager_nml',err_msg)) return
endif
call mpp_close (iunit)
if(mix_snapshot_average_fields) then
   if(mpp_pe() == mpp_root_pe()) call mpp_error(WARNING,'Namelist '// &
        'mix_snapshot_average_fields = true will cause ERROR in time coordinates '// &
        'of all time_averaged fields. Strongly recommend mix_snapshot_average_fields = false')
endif
allocate(output_fields(max_output_fields))
allocate(input_fields(max_input_fields))
if (.not.file_exist('diag_table') ) then
  if(fms_error_handler('diag_manager_init','file diag_table nonexistent', err_msg)) return
endif
call mpp_open(iunit, 'diag_table', action=MPP_RDONLY)
! Read in the global file labeling string
read(iunit, *, end = 99, err=99) global_descriptor

! Read in the base date
read(iunit, *, end = 99, err = 99) base_year, base_month, base_day, &
   base_hour, base_minute, base_second

! Set up the time type for base time
if (get_calendar_type() /= NO_CALENDAR) then
   if(base_year==0 .or. base_month==0 .or. base_day==0) then
     if(fms_error_handler('diag_manager_init','base_year/month/day can not equal zero', err_msg)) return
   endif
   base_time = set_date(base_year, base_month, base_day, base_hour, base_minute, base_second)
   amonth = month_name(base_month)
else
! No calendar - ignore year and month
   base_time = set_time(nint(base_hour*SECONDS_PER_HOUR)+nint(base_minute*SECONDS_PER_MINUTE)+base_second, base_day)
   base_year  = 0
   base_month = 0
   amonth = 'day'
end if

allocate(pelist(mpp_npes()))
call mpp_get_current_pelist(pelist, pelist_name)

nrecs=0
nfiles=0
if(debug_diag_manager) then
   write(stdout(), *) ' '
   write(stdout(), *) '******** Summary of output FILES information from diag_manager *********'
   write(stdout(), *) ' '
   write(stdout(), *) 'file name     ', '       saving frequency   ', '      saving frequency unit'
   write(stdout(), *) ' '
endif
do while (nfiles <= max_files)
   file_freq_present=.false.; file_start_present=.false.; file_duration_present=.false.
   read(iunit,'(a)',end=86,err=85) record
   nrecs=nrecs+1
   if (record(1:1) == '#') cycle   
   record_len = len_trim(record)
   if(record_len == 0) cycle
   time_pos = MAX(index(record,'time'), index(record,'Time'))
! We are only looking for time axis here.
   if(time_pos == 0) cycle

!-- NZ   
!-- Even this can be problematic if one of the fields in table contains the string "time".
!-- This part of the code (sparcing the diag_table) needs enhancements.
!-- How about designing a .xml format for diag_table?
!-- Fortran is not good for regexp. How about a Perl preprocessing of the table?

   record_tail = record(time_pos+6:record_len)
   record_tail_len = len_trim(record_tail)
   if (time_pos + 7 < record_len) then
      file_freq_present = .true.
      if(INDEX(record_tail,'s',back=.true.)== INDEX(record_tail,'s')) then         
         if(INDEX(record_tail,'s') + 3 < record_tail_len) file_start_present=.true.
      else
         file_start_present=.true.
         file_duration_present = .true.
      endif
   endif
! Start reading file information  
   if(.not. file_freq_present) then
      read(record,*,err=85,end=85)textA%name,textA%output_freq,textA%output_freq_units, &
           textA%format,textA%time_units,textA%long_name   
   elseif(.not.file_start_present) then
      read(record,*,err=85,end=85)textA%name,textA%output_freq,textA%output_freq_units, &
           textA%format,textA%time_units,textA%long_name,textA%new_file_freq,textA%new_file_freq_units
   elseif(.not. file_duration_present) then
      read(record,*,err=85,end=85) textA%name,textA%output_freq,textA%output_freq_units, &
           textA%format,textA%time_units,textA%long_name,textA%new_file_freq,&
           textA%new_file_freq_units, textA%start_time_s
   else
      read(record,*,err=85,end=85) textA
   endif
! does the record contain start time for the file ???   
   if(file_start_present) then
      read(textA%start_time_s,*,end=85,err=85) yr, mo, dy, hr, mi, sc
      start_time = set_date( yr, mo, dy, hr, mi, sc, err_msg=err_msg_local)
      if(err_msg_local /= '') then
        if(fms_error_handler('diag_manager_init',err_msg_local,err_msg)) return
      endif
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
   file_duration_units = 0
   do j = 1, size(time_unit_list(:))
      if(textA%time_units == time_unit_list(j)) time_units = j
      if(textA%output_freq_units == time_unit_list(j)) output_freq_units = j
      if(file_freq_present) then
         if(textA%new_file_freq_units == time_unit_list(j)) new_file_freq_units = j
      endif
      if(file_duration_present) then
        if(textA%file_duration_units == time_unit_list(j))  file_duration_units = j
     endif
   end do
   if(time_units == 0) then
     if(fms_error_handler('diag_manager_init','invalid time units, check time unit in diag_table',err_msg)) return
   endif
   if(output_freq_units == 0) then
     if(fms_error_handler('diag_manager_init','invalid output frequency units, check diag table',err_msg)) return
   endif
   if(file_freq_present .and. new_file_freq_units == 0 ) then
     if(fms_error_handler('diag_manager_init','invalid new_file frequency units, check diag table',err_msg)) return
   endif
   if(file_duration_present .and. file_duration_units == 0 ) then
     if(fms_error_handler('diag_manager_init','invalid new_file frequency units, check diag table',err_msg)) return
   endif

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
   if(file_freq_present) then
      if(file_duration_present) then
         file_duration = textA%file_duration
      else
         file_duration = textA%new_file_freq
         file_duration_units = new_file_freq_units
      endif
      call init_file(textA%name,textA%output_freq, output_freq_units, &
           textA%format, time_units,textA%long_name,textA%new_file_freq,new_file_freq_units, start_time, &
           file_duration, file_duration_units )
   else
      call init_file(textA%name,textA%output_freq, output_freq_units, &
           textA%format, time_units,textA%long_name)
   endif
   if(debug_diag_manager) write(stdout(), 1)textA%name, textA%output_freq, textA%output_freq_units 

85 continue
enddo
if(fms_error_handler('diag_manager_init','max_files exceeded, increase max_files',err_msg)) return
86 continue
1 format(1x,A21, 4x,I5,4x,A30)
write(stdout(), *)' '
if(debug_diag_manager) write(stdout(), *)'************************************************************************'
rewind(iunit)

if(debug_diag_manager) then
   write(stdout(), *) ' '
   write(stdout(), *) '******* Summary of output FIELDS information from diag_manager *********'
   write(stdout(), *) ' '
   write(stdout(), *) 'module         ',' field name      ', '   file name      ', '    time method '
   write(stdout(), *) ' '
endif


nfields=0; nrecs=0
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
   if(debug_diag_manager) &
        write(stdout(),2)textB%module_name,textB%field_name,textB%name,textB%time_method 
   if(trim(textB%spatial_ops) == 'none') then
      call init_output_field(textB%module_name,textB%field_name,textB%output_name,&
           textB%name,textB%time_method,textB%pack)
   else
      call init_output_field(textB%module_name,textB%field_name,textB%output_name,&
           textB%name,textB%time_method,textB%pack, local_coord)
   endif
93 continue
enddo
if(fms_error_handler('diag_manager_init','max_output_fields exceeded, increase it via diag_manager_nml',err_msg)) return
94 continue
call close_file(iunit)
2 format(1x,A15,1x,A15,4x,A21,A10)
if(debug_diag_manager) write(stdout(), *)'************************************************************************'
! check duplicate output_fields in the diag_table
call check_duplicate_output_fields(err_msg=err_msg_local)
if(err_msg_local /= '') then
  if(fms_error_handler('diag_manager_init',err_msg_local,err_msg)) return
endif
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
if(fms_error_handler('diag_manager_init','error reading table',err_msg)) return
end subroutine diag_manager_init
! </SUBROUTINE>

!----------------------------------------------------------------------
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

function need_data(diag_field_id, next_model_time)
!
! next_model_time = current model time + model time_step
!
type (time_type), intent(in) :: next_model_time
integer, intent(in)          :: diag_field_id
logical                      :: need_data
integer                      :: i, out_num 

need_data=.false.
if (diag_field_id < 0 ) return   !this field is unused
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
!###################################################################################################


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
!###################################################################################################

#ifdef test_diag_manager

 ! This program runs only one of many possible tests with each execution.
 ! Each test ends with an intentional fatal error.
 ! diag_manager_mod is not a stateless module, and there are situations
 ! where a fatal error leaves the module in a state that does not allow
 ! it to function properly if used again. Therefore, the program must
 ! be terminated after each intentional fatal error.

 ! Each test is dependent on the diag_table, and different diag_tables
 ! exist for each test. Depending on the test, an intentional fatal error
 ! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
 ! Because of this, the calls to all of those routines differ depending on the test.

 ! The diag_table for each test is included below.

 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 1

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 2

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 3

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 4

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 !  "diag_test2", 1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 !  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 5

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 !  "diag_test2", 1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 !  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 6

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 !  "diag_test2", 1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 !  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 7

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 8

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 !  "diag_test2", 1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 !  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 9

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 10

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 11

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 12

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 ! # Test of the error check that duplicate field names do not appear in same file,
 !  "test_mod",              "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 13

 ! test_diag_manager
 ! 1 3 1 0 0 0
 ! #output files
 !  "diag_test",  1, "days", 1, "days", "time",
 !  "diag_test2", 1, "months", 1, "days", "time",
 ! #output variables
 !  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
 ! # Test of WARNING message that no data is written when run length is less than output interval  
 !  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------
 ! diag_table for test 14

 ! test_diag_manager
 ! 1990 1 29 0 0 0
 ! #output files
 !  "diag_test2", 1, "months", 1, "days", "time",
 ! #output variables
 ! # Test of check for invalid date. (Jan 29 1990 + one month = Feb 29 1990)
 !  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
 !--------------------------------------------------------------------------------------------------

 program test

 ! This program runs only one of many possible tests with each execution.
 ! Each test ends with an intentional fatal error.
 ! diag_manager_mod is not a stateless module, and there are situations
 ! where a fatal error leaves the module in a state that does not allow
 ! it to function properly if used again. Therefore, the program must
 ! be terminated after each intentional fatal error.

 ! Each test is dependent on the diag_table, and different diag_tables
 ! exist for each test. Depending on the test, an intentional fatal error
 ! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
 ! Because of this, the calls to all of those routines differ depending on the test.

 use          mpp_mod, only: mpp_pe
 use  mpp_domains_mod, only: domain2d, mpp_define_domains, mpp_get_compute_domain
 use          fms_mod, only: fms_init, fms_end, mpp_npes, file_exist, open_namelist_file, check_nml_error, close_file, open_file
 use          fms_mod, only: error_mesg, FATAL, stdlog
 use       fms_io_mod, only: fms_io_exit
 use    constants_mod, only: constants_init, pi

 use time_manager_mod, only: time_type, set_calendar_type, set_date, decrement_date, operator(+), set_time
 use time_manager_mod, only: NOLEAP, JULIAN, GREGORIAN, THIRTY_DAY_MONTHS

 use diag_manager_mod, only: diag_manager_init, send_data, diag_axis_init, diag_manager_end
 use diag_manager_mod, only: register_static_field, register_diag_field

 implicit none

 type(domain2d) :: Domain1
 type(domain2d) :: Domain2
 integer, parameter :: nlon1=18, nlat1=18, nlev=2
 integer, parameter :: nlon2=36, nlat2=36
 real, dimension(nlon1  ) :: lon_global1
 real, dimension(nlon1+1) :: lonb_global1
 real, dimension(nlat1)   :: lat_global1
 real, dimension(nlat1+1) :: latb_global1
 real, dimension(nlon2  ) :: lon_global2
 real, dimension(nlon2+1) :: lonb_global2
 real, dimension(nlat2)   :: lat_global2
 real, dimension(nlat2+1) :: latb_global2
 real, dimension(nlev  )  :: pfull, bk
 real, dimension(nlev+1)  :: phalf
 real, allocatable, dimension(:) :: lon1, lat1, lonb1, latb1
 real, allocatable, dimension(:) :: lon2, lat2, lonb2, latb2
 real, allocatable, dimension(:,:,:) :: dat1, dat1h
 real, allocatable, dimension(:,:,:) :: dat2, dat2h
 real, allocatable, dimension(:,:) :: dat2_2d
 integer :: id_phalf, id_pfull, id_bk
 integer :: id_lon1, id_lonb1, id_latb1, id_lat1, id_dat1
 integer :: id_lon2, id_lat2, id_dat2, id_dat2_2d
 integer :: i, j, k, is1, ie1, js1, je1, nml_unit, io, ierr, log_unit, out_unit
 integer :: is_in, ie_in, js_in, je_in
 integer :: is2, ie2, js2, je2, hi=1, hj=1
 real :: rad_to_deg, dp, surf_press=1.e5
 type(time_type) :: Time
 logical :: used, test_successful
 integer, dimension(2) :: layout = (/0,0/)
 integer :: test_number=1
 character(len=256) :: err_msg

 namelist / test_diag_manager_nml / layout, test_number

 call fms_init
 nml_unit = open_namelist_file()
 log_unit = stdlog()
 out_unit = open_file(file='test_diag_manager.out', form='formatted', threading='multi', action='write')
 call constants_init
 call set_calendar_type(JULIAN)

 rad_to_deg = 180./pi

 if (file_exist('input.nml')) then
   ierr=1
   do while (ierr /= 0)
     read(nml_unit, nml=test_diag_manager_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'test_diag_manager_nml')
   enddo
10 call close_file(nml_unit)
 endif
 write(log_unit,test_diag_manager_nml)

 if(test_number == 12) then
   call diag_manager_init(err_msg=err_msg)
   if(err_msg /= '') then
     write(out_unit,'(a)') 'test12 successful: err_msg='//trim(err_msg)
     call error_mesg('test_diag_manager','test12 successful.',FATAL)
   else
     write(out_unit,'(a)') 'test12 fails'
     call error_mesg('test_diag_manager','test12 fails',FATAL)
   endif
 else
   call diag_manager_init
 endif

 if(any(layout == (/0,0/))) then
   layout = (/1,mpp_npes()/)
 endif
 call mpp_define_domains( (/1,nlon1,1,nlat1/), layout, Domain1, name='test_diag_manager')
 call mpp_get_compute_domain(Domain1, is1, ie1, js1, je1)
 allocate(lon1(is1:ie1), lat1(js1:je1), lonb1(is1:ie1+1), latb1(js1:je1+1))
 call compute_grid(nlon1, nlat1, is1,ie1,js1,je1, lon_global1, lat_global1, lonb_global1, latb_global1, lon1, lat1, lonb1, latb1)
 call mpp_define_domains( (/1,nlon2,1,nlat2/), layout, Domain2, name='test_diag_manager')
 call mpp_get_compute_domain(Domain2, is2, ie2, js2, je2)
 allocate(lon2(is2:ie2), lat2(js2:je2), lonb2(is2:ie2+1), latb2(js2:je2+1))
 call compute_grid(nlon2, nlat2, is2,ie2,js2,je2, lon_global2, lat_global2, lonb_global2, latb_global2, lon2, lat2, lonb2, latb2)
 dp = surf_press/nlev
 do k=1,nlev+1
   phalf(k) = dp*(k-1)
 enddo
 do k=1,nlev
   pfull(k) = .5*(phalf(k) + phalf(k+1))
   bk(k) = pfull(k)/surf_press
 enddo

 allocate(dat1 (is1   :ie1   , js1   :je1   , nlev))
 allocate(dat1h(is1-hi:ie1+hi, js1-hj:je1+hj, nlev))
 dat1h = 0.
 do j=js1,je1
   do i=is1,ie1
     dat1 (i,j,1) = sin(lon1(i))*cos(lat1(j))
   enddo
 enddo
 dat1h(is1:ie1,js1:je1,1) = dat1(:,:,1)
 dat1 (:,:,2) = -dat1 (:,:,1)
 dat1h(:,:,2) = -dat1h(:,:,1)

 allocate(dat2    (is2   :ie2   , js2   : je2, nlev))
 allocate(dat2_2d (is2   :ie2   , js2   : je2 ))
 allocate(dat2h(is2-hi:ie2+hi, js2-hj: je2+hj, nlev))
 dat2h = 0.
 do j=js2,je2
   do i=is2,ie2
     dat2(i,j,1) = sin(lon2(i))*cos(lat2(j))
   enddo
 enddo
 dat2h(is2:ie2,js2:je2,1) = dat2(:,:,1)
 dat2 (:,:,2) = -dat2(:,:,1)
 dat2h(:,:,2) = -dat2h(:,:,1)
 dat2_2d = dat2(:,:,1)

 id_lonb1 = diag_axis_init('lonb1', rad_to_deg*lonb_global1, 'degrees_E', 'x', long_name='longitude edges', Domain2=Domain1)
 id_latb1 = diag_axis_init('latb1', rad_to_deg*latb_global1, 'degrees_N', 'y', long_name='latitude edges',  Domain2=Domain1)

 id_lon1  = diag_axis_init('lon1',  rad_to_deg*lon_global1, 'degrees_E','x',long_name='longitude',Domain2=Domain1,edges=id_lonb1)
 id_lat1  = diag_axis_init('lat1',  rad_to_deg*lat_global1, 'degrees_N','y',long_name='latitude', Domain2=Domain1,edges=id_latb1)

 id_phalf= diag_axis_init('phalf', phalf, 'Pa', 'z', long_name='half pressure level', direction=-1)
 id_pfull= diag_axis_init('pfull', pfull, 'Pa', 'z', long_name='full pressure level', direction=-1, edges=id_phalf)

 id_lon2 = diag_axis_init('lon2',  rad_to_deg*lon_global2,  'degrees_E', 'x', long_name='longitude', Domain2=Domain2)
 id_lat2 = diag_axis_init('lat2',  rad_to_deg*lat_global2,  'degrees_N', 'y', long_name='latitude',  Domain2=Domain2)

 if(test_number == 14) then
   Time = set_date(1990,1,29,0,0,0)
 else
   Time = set_date(1990,1,1,0,0,0)
 endif

 id_dat1 = register_diag_field('test_diag_manager_mod', 'dat1', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K')
 id_dat2 = register_diag_field('test_diag_manager_mod', 'dat2', (/id_lon2,id_lat2,id_pfull/), Time, 'sample data', 'K')

 if(test_number == 14) then
   id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data', 'K', err_msg=err_msg)
   if(err_msg /= '') then
     write(out_unit,'(a)') 'test14 successful. err_msg='//trim(err_msg)
   else
     write(out_unit,'(a)') 'test14 fails.'
   endif
 else
   id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data', 'K')
 endif

 id_bk = register_static_field('test_diag_manager_mod', 'bk', (/id_pfull/), 'half level sigma', 'none')

 if(test_number == 13) then
   if(id_dat2_2d > 0) used=send_data(id_dat2_2d, dat2(:,:,1), Time, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test13: successful if a WARNING message appears that refers to output interval greater than runlength'
   else
     write(out_unit,'(a)') 'test13 fails: err_msg='//trim(err_msg)
   endif
 endif

! Note: test12 involves diag_manager_init, it does not require a call to send_data.
!       See call to diag_manager_init above.

 if(test_number == 11) then
   is_in=1+hi
   js_in=1+hj
   ie_in=ie2-is2+1+hi
   je_in=je2-js2+1+hj

   if(id_dat2 > 0) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test11.1 successful.'
   else
     write(out_unit,'(a)') 'test11.1 fails. err_msg='//trim(err_msg)
   endif

!  intentional_error: je_in is missing
   if(id_dat2 > 0) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test11.2 fails.'
   else
     write(out_unit,'(a)') 'test11.2 successful. err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 10) then
!  1 window, no halos, static, 1 dimension, global data.

   if(id_bk > 0) used = send_data(id_bk, bk, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test10.1 successful.'
   else
     write(out_unit,'(a)') 'test10.1 fails: err_msg='//trim(err_msg)
   endif

!  intentional_error: data array too large.
   if(id_bk > 0) used = send_data(id_bk, phalf, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test10.2 fails.'
   else
     write(out_unit,'(a)') 'test10.2 successful: err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 9) then
!  1 window, no halos, static, 1 dimension, global data

   if(id_bk > 0) used = send_data(id_bk, bk, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test9.1 successful.'
   else
     write(out_unit,'(a)') 'test9.1 fails: err_msg='//trim(err_msg)
   endif

!  intentional_error: data array too small
   if(id_bk > 0) used = send_data(id_bk, bk(1:nlev-1), err_msg=err_msg) ! intentional_error
   if(err_msg == '') then
     write(out_unit,'(a)') 'test9.2 fails.'
   else
     write(out_unit,'(a)') 'test9.2 successful: err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 8) then
!  1 window with halos
   is_in=1+hi; js_in=1+hj

   ie_in=ie2-is2+1+hi
   je_in=je2-js2+1+hj
   if(id_dat2 > 0) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, &
                   ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test8.1 successful.'
   else
     write(out_unit,'(a)') 'test8.1 fails: err_msg='//trim(err_msg)
   endif

!  intentional_error: data array too small in both x and y directions
!  Error check is done on second call to send_data. Change in value of Time triggers the check.
   Time = Time + set_time(0,1)
   ie_in=ie1-is1+1+hi
   je_in=je1-js1+1+hj
   if(id_dat2 > 0) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in, &
                   ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
   Time = Time + set_time(0,1)
   if(id_dat2 > 0) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in, &
                   ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test8.2 fails.'
   else
     write(out_unit,'(a)') 'test8.2 successful: err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 7) then
!  1 window with halos
   is_in=1+hi; js_in=1+hj

   ie_in=ie1-is1+1+hi
   je_in=je1-js1+1+hj
   if(id_dat1 > 0) used=send_data(id_dat1, dat1h, Time, is_in=is_in, js_in=js_in, &
                   ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test7.1 successful.'
   else
     write(out_unit,'(a)') 'test7.1 fails: err_msg='//trim(err_msg)
   endif

!  intentional_error: data array too large in both x and y directions
   ie_in=ie2-is2+1+hi
   je_in=je2-js2+1+hj
   if(id_dat1 > 0) used=send_data(id_dat1, dat2h, Time, is_in=is_in, js_in=js_in, &
                   ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test7.2 fails.'
   else
     write(out_unit,'(a)') 'test7.2 successful: err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 6) then
!  multiple windows, no halos

!  No error messages should appear at any point within either do loop for test6.1
   test_successful = .true.
   do i=is2,ie2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
     if(err_msg /= '') then
       write(out_unit,'(a)') 'test6.1 fails: err_msg='//trim(err_msg)
       test_successful = .false.
     endif
   enddo
   Time = Time + set_time(0,1)
   do i=is2,ie2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
     if(err_msg /= '') then
       write(out_unit,'(a)') 'test6.1 fails: err_msg='//trim(err_msg)
       test_successful = .false.
     endif
   enddo
   if(test_successful) then
     write(out_unit,'(a)') 'test6.1 successful.'
   else
     write(out_unit,'(a)') 'test6.1 fails.'
   endif

!  intentional_error: data array too small in y direction
!  Error check is done on second call to send_data. Change in value of Time triggers the check.
   Time = Time + set_time(0,1)
   do i=is2,ie2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1)
   enddo
   Time = Time + set_time(0,1)
   do i=is2,ie2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1, err_msg=err_msg)
     if(err_msg /= '') exit ! exit immediately after error is detected. No need to continue.
   enddo
   if(err_msg == '') then
     write(out_unit,'(a)') 'test6.2 fails.'
   else
     write(out_unit,'(a)') 'test6.2 successful: err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 5) then
!  multiple windows, no halos

!  No error messages should appear at any point within either do loop for test5.1
   test_successful = .true.
   do j=js2,je2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
     if(err_msg /= '') then
       write(out_unit,'(a)') 'test5.1 fails: err_msg='//trim(err_msg)
       test_successful = .false.
     endif
   enddo
   Time = Time + set_time(0,1)
   do j=js2,je2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
     if(err_msg /= '') then
       write(out_unit,'(a)') 'test5.1 fails: err_msg='//trim(err_msg)
       test_successful = .false.
     endif
   enddo
   if(test_successful) then
     write(out_unit,'(a)') 'test5.1 successful.'
   else
     write(out_unit,'(a)') 'test5.1 fails.'
   endif

!  intentional_error: data array too small in x direction.
!  Error check is done on second call to send_data. Change in value of Time triggers the check.
   Time = Time + set_time(0,1)
   do j=js2,je2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1)
   enddo
   Time = Time + set_time(0,1)
   do j=js2,je2
     if(id_dat2 > 0) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
     if(err_msg /= '') exit ! exit immediately after error is detected. No need to continue.
   enddo
   if(err_msg == '') then
     write(out_unit,'(a)') 'test5.2 fails.'
   else
     write(out_unit,'(a)') 'test5.2 successful: err_msg='//trim(err_msg)
   endif
 endif

 if(test_number == 4) then
!  1 window, no halos

   if(id_dat2 > 0) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
   Time = Time + set_time(0,1)
   if(id_dat2 > 0) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test4.1 successful.'
   else
     write(out_unit,'(a)') 'test4.1 fails: err_msg='//trim(err_msg)
   endif

!  intentional_error: data array too small in both x and y directions
!  Error check is done on second call to send_data. Change in value of Time triggers the check.
   Time = Time + set_time(0,1)
   if(id_dat2 > 0) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
   Time = Time + set_time(0,1)
   if(id_dat2 > 0) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test4.2 fails.'
   else
     write(out_unit,'(a)') 'test4.2 successful: err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 3) then
!  multiple windows, no halos

!  No error messages should appear at any point within do loop for test3.1
   test_successful = .true.
   do i=is1,ie1
     if(id_dat1 > 0) used = send_data(id_dat1, dat1(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
     if(err_msg /= '') then
       write(out_unit,'(a)') 'test3.1 fails: err_msg='//trim(err_msg)
       test_successful = .false.
     endif
   enddo
   if(test_successful) then
     write(out_unit,'(a)') 'test3.1 successful.'
   else
     write(out_unit,'(a)') 'test3.1 fails.'
   endif

!  intentional_error: data array too large in y direction
   do i=is1,ie1
     if(id_dat1 > 0) used = send_data(id_dat1, dat2(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
     if(err_msg /= '') exit ! exit immediately after error is detected. No need to continue.
   enddo
   if(err_msg == '') then
     write(out_unit,'(a)') 'test3.2 fails.'
   else
     write(out_unit,'(a)') 'test3.2 successful: err_msg='//trim(err_msg)
   endif

 endif

 if(test_number == 2) then
!  multiple windows, no halos

!  No error messages should appear at any point within do loop for test2.1
   test_successful = .true.
   do j=js1,je1
     if(id_dat1 > 0) used = send_data(id_dat1, dat1(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
     if(err_msg /= '') then
       write(out_unit,'(a)') 'test2.1 fails: err_msg='//trim(err_msg)
       test_successful = .false.
     endif
   enddo
   if(test_successful) then
     write(out_unit,'(a)') 'test2.1 successful.'
   else
     write(out_unit,'(a)') 'test2.1 fails.'
   endif

!  intentional_error: data array too large in x direction
   do j=js1,je1
     if(id_dat1 > 0) used = send_data(id_dat1, dat2(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
     if(err_msg /= '') exit ! exit immediately after error is detected. No need to continue.
   enddo
   if(err_msg == '') then
     write(out_unit,'(a)') 'test2.2 fails.'
   else
     write(out_unit,'(a)') 'test2.2 successful: err_msg='//trim(err_msg)
   endif

 endif
 if(test_number == 1) then
!  1 window, no halos
   if(id_dat1 > 0) used = send_data(id_dat1, dat2, Time, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test1.1 fails: Intentional error not detected'
   else
     write(out_unit,'(a)') 'test1.1 successful: '//trim(err_msg)
   endif

!  intentional_error: data array too large in both x and y directions
   if(id_dat1 > 0) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)
   if(err_msg == '') then
     write(out_unit,'(a)') 'test1.2 successful'
   else
     write(out_unit,'(a)') 'test1.2 fails: '//trim(err_msg)
   endif

 endif

 call diag_manager_end(Time)
 call fms_io_exit
 call fms_end

 contains
!=================================================================================================================================
 subroutine compute_grid(nlon, nlat, is, ie, js, je, lon_global, lat_global, lonb_global, latb_global, lon, lat, lonb, latb)
 integer, intent(in) :: nlon, nlat, is, ie, js, je
 real, intent(out), dimension(:) :: lon_global, lat_global, lonb_global, latb_global, lon, lat, lonb, latb
 real :: dlon, dlat
 integer :: i, j

 dlon = 2*pi/nlon
 dlat =   pi/nlat

 do i=1,nlon+1
   lonb_global(i) = dlon*(i-1)
 enddo
 do j=1,nlat+1
   latb_global(j) = dlat*(j-1) - .5*pi
 enddo
 do i=1,nlon
   lon_global(i) = .5*(lonb_global(i) + lonb_global(i+1))
 enddo
 do j=1,nlat
   lat_global(j) = .5*(latb_global(j) + latb_global(j+1))
 enddo
 lon  = lon_global(is:ie)
 lat  = lat_global(js:je)
 lonb = lonb_global(is:ie+1)
 latb = latb_global(js:je+1)

 end subroutine compute_grid
!=================================================================================================================================
 end program test
#endif
