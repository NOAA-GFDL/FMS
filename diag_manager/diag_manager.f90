module diag_manager_mod

! <CONTACT EMAIL="Matthew.Harrison@noaa.gov">
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
!   Run-time specification of diagnostics are input through the diagnostics table.
!<BR/>
!
!
! <B> Usage</B> of <TT> diag_manager</TT> includes the following steps:<BR/>
! 1. Create diag_table as described in 
!   <LINK SRC="diag_table_tk.html">diag_table_tk</LINK> documentation.<BR/>
! 2. Call <TT>diag_manager_init</TT> to initialize diag_manager_mod <BR/>
! 3. Call <TT> register_diag_field</TT> to register the field to be outputted. <BR/>
!    <B> NOTE</B> ALL fields in diag_table should be registered BEFORE the first send_data call
! 4. Call <TT>send_data</TT> to send data to output fields <BR/>
! 5. Call <TT>diag_manager_end</TT> to exit diag_manager <BR/>
!   
!
! <B> New Features </B> of <TT> diag_manager_mod </TT>: <BR/>
! 1. Ability to output scalars (see <TT>send_data </TT>and <TT> register_diag_field</TT>)<BR/>
! 2. Ability to output time average of fields that have time dependent mask. <BR/>
! 3. Give optional warning if <TT>register_diag_field </TT>fails due to misspelled module name
!    or field name. <BR/>
! 4. Give warning if a field in diag_table is never registered by the end of the program. <BR/>
! 5. Check for duplicate lines in diag_table<BR/>
! 6. diag_table can contain fields that are NOT written to any files. The file name in diag_table
!    of these fields is <TT> null</TT> <BR/>
! 7. By default, a field is output in its global grid, it is now possible to output a field in 
!    a region specified by user, see send_data for more details <BR/>
!
!   <B>Features of <TT>diag_manager_mod</TT> </B>include: <BR/>
!   Simple, minimal API.<BR/>
!   Run-time choice of diagnostics.<BR/>
!   Self-describing files: comprehensive header information
!   (metadata) in the file itself.<BR/>
!   Strong parallel write performance.<BR/>
!   Integrated netCDF capability: <LINK
!   SRC="http://www.unidata.ucar.edu/packages/netcdf/">netCDF</LINK> is a
!   data format widely used in the climate/weather modeling
!   community. netCDF is considered the principal medium of data storage
!   for <TT>diag_manager_mod</TT>. Raw unformatted
!   fortran I/O capability is also available.<BR/>
!   Requires off-line post-processing: a tool for this purpose,
!   <TT>mppnccombine</TT>, is available. GFDL users may use
!   <TT>~hnv/pub/mppnccombine</TT>. Outside users may obtain the
!   source <LINK
!   SRC="ftp://ftp.gfdl.gov/perm/hnv/mpp/mppnccombine.c">here</LINK>.  It
!   can be compiled on any C compiler and linked with the netCDF
!   library. The program is free and is covered by the <LINK
!   SRC="ftp://ftp.gfdl.gov/perm/hnv/mpp/LICENSE">GPL license</LINK>.
! </DESCRIPTION>

use time_manager_mod, only: get_time, set_time, get_date, set_date,    &
                            increment_date, operator(-), operator(>=), &
                            operator(>), operator(<), operator(==),    &
                            time_type, increment_time, month_name,     &
                            get_calendar_type, NO_CALENDAR
use       mpp_io_mod, only: mpp_open, MPP_RDONLY
use          fms_mod, only: error_mesg, FATAL, WARNING, NOTE,          &
                            close_file, stdlog, write_version_number,  &
                            file_exist, mpp_pe, open_namelist_file, &
                            check_nml_error, lowercase
use mpp_mod, only         : mpp_get_current_pelist, mpp_npes, mpp_sync, mpp_pe, mpp_root_pe,stdout
use diag_axis_mod, only   : diag_axis_init, get_axis_length, get_diag_axis, get_domain1d, get_domain2d, &
     get_axis_global_length, diag_subaxes_init,get_diag_axis_cart, get_diag_axis_data
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
integer, parameter :: max_fields_per_file = 150
integer, parameter :: max_out_per_in_field = 10 !max number of output_fields per input_field
integer, parameter :: max_files = 21
integer            :: num_files = 0
integer, parameter :: max_input_fields = 300
integer            :: num_input_fields = 0
integer, parameter :: max_output_fields = 300
integer            :: num_output_fields = 0
integer, parameter :: DIAG_OTHER = 0
integer, parameter :: DIAG_OCEAN = 1
integer, parameter :: DIAG_ALL   = 2
real               :: EMPTY = 0.0
integer            :: null_axis_id

! Global data for all files
type (time_type)    :: base_time
integer             :: base_year, base_month, base_day, base_hour, base_minute, base_second
character(len = 256):: global_descriptor

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
   integer             :: time_axis_id
   type (time_type)    :: last_flush
   type(diag_fieldtype):: f_avg_start, f_avg_end, f_avg_nitems
   logical             :: local ! true if fields are outputted in a region instead of global
end type file_type

type input_field_type
   character(len=128) :: module_name, field_name, long_name, units
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
   integer             :: pack
   real, pointer       :: buffer(:, :, :)=>null()
   real, pointer       :: counter(:, :, :)=>null()
   type(time_type)     :: last_output, next_output, next_next_output
   real                :: count_0d
   type(diag_fieldtype):: f_type
   integer             :: axes(3), num_axes, num_elements, total_elements
   type(diag_grid)     :: output_grid
   logical             :: local_output, need_compute
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
type (input_field_type) :: input_fields(max_input_fields)
type (output_field_type), save :: output_fields(max_output_fields)

logical             :: first_send_data_call = .true.
logical             :: first_send_data_local = .true.
logical             :: module_is_initialized = .false.
integer, parameter  :: EVERY_TIME =  0
integer, parameter  :: END_OF_RUN = -1
integer, parameter  :: DIAG_SECONDS = 1, DIAG_MINUTES = 2, DIAG_HOURS = 3
integer, parameter  :: DIAG_DAYS = 4, DIAG_MONTHS = 5, DIAG_YEARS = 6
character (len=10)  :: time_unit_list(6) = (/'seconds   ', 'minutes   ', &
   'hours     ', 'days      ', 'months    ', 'years     '/)
character (len = 7) :: avg_name = 'average'
character(len=32)   :: pelist_name

! version number of this module
character(len=128) :: version = '$Id: diag_manager.f90,v 10.0 2003/10/24 22:01:27 fms Exp $'
character(len=128) :: tagname = '$Name: jakarta $'  


! <INTERFACE NAME="send_data">
! <TEMPLATE>
!send_data(diag_field_id, field, time, is_in, js_in, ks_in, &
!             mask, rmask, weight)
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
! For the real  mask, the mask is applied if the mask value is less than 0.5.  The
! weight array is currently not implemented.
!
! By default, a field will be written out entirely in its global grid. Users can also specify
! region in which the field will be outputted. The region is specified in diag-table just before
! the end of output_field replacing "none". For example:<BR/>
! by default:<BR/>
! "ocean_mod","Vorticity","vorticity","file1","all",.false.,"none",2 <BR/>
! for regional output:<BR/>
! "ocean_mod","Vorticity","vorticity_local","file2","all",.false.,"0.5 53.5 -89.5 -28.5 1 1",2<BR/>
! the format of region is "xbegin xend ybegin yend zbegin zend". If it is a 2D field use (1 1)
! for (zbegin zend) as in the example above. The units used for region are the actual units used in 
! grid_spec.nc (for example lat, lon, meter). Fatal error will occur if region's boundaries
! are not found in grid_spec.nc<BR/>
! <BR/>
! Special note when using regional output: result files containing regional outputs should be 
! different from files containing global (default) output. It is FATAL error to have one file
! containing both regional and global results. For maximum flexibility and independence from
! PE counts one file should contain just one region
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
type(time_type), intent(in)            :: init_time
character(len=*), optional, intent(in) :: long_name, units
real, optional, intent(in)             :: missing_value, range(2)
 
register_diag_field_scalar = register_diag_field_array(module_name, field_name,&
   (/null_axis_id/), init_time,long_name, units, missing_value, range)

end function register_diag_field_scalar

! <FUNCTION NAME="register_diag_field">

!<OVERVIEW>
!     Register Diagnostic Field.
!</OVERVIEW>
!<DESCRIPTION>
!     Return field index for subsequent calls to <LINK SRC="#send_data"> send_data </LINK> <BR/>
!<TT> axes</TT> are axis ID returned from <TT>diag_axis_init</TT>, <TT>axes</TT>  are required
! for fields of 1-3 dimension and NOT required for scalars. <BR/>
! optional <TT> mask_variant</TT> is for fields that have a time-dependent mask. If <TT>mask_variant</TT> is
! true then <TT>mask</TT> must be present in argument list of <TT>send_data</TT> <BR/>
! When optional <TT> verbose</TT> is true a warning will be given if <TT> register_diag_field</TT>
! fails due to misspelled field name or module name. The default <TT> verbose</TT> is false.
!</DESCRIPTION>
!<TEMPLATE>
!     register_diag_field(module_name,field_name,axes,init_time, &
!     long_name,units,missing_value,range,mask_variant,verbose)
!</TEMPLATE>

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
   long_name, units, missing_value, range, mask_variant,verbose)

! Indicates the calling modules intent to supply data for this field.

integer                                :: register_diag_field_array
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
type(time_type), intent(in)            :: init_time
character(len=*), optional, intent(in) :: long_name, units
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
! Need set output fields to not static, too
! Need to loop through all output_fields associated and allocate their buffers
! Would have to consider horizontal operations at this point
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
              &' will be outputtted in region:'//trim(msg)
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
   long_name, units, missing_value, range, mask_variant, require, dynamic)

integer                                :: register_static_field
character(len=*), intent(in)           :: module_name, field_name
integer, intent(in)                    :: axes(:)
character(len=*), optional, intent(in) :: long_name, units
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

register_static_field = find_input_field(module_name, field_name)

if (PRESENT(require)) then
   if (require) then
      call init_input_field(module_name, field_name)
      register_static_field = find_input_field(module_name, field_name)
      do j=1, num_files
! need to think about what to do if the axes are not present, e.g. file only
! contains data slices
         call init_output_field(module_name, field_name,field_name, &
                               files(j)%name,.false.,2)
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
num_axes = size(axes)
input_fields(field)%axes(1:num_axes) = axes
input_fields(field)%num_axes = num_axes

! Need to check for present, otherwise defaults
if(present(long_name)) then
   input_fields(field)%long_name = trim(long_name)
else
   input_fields(field)%long_name = input_fields(field)%field_name
endif
if(present(units)) then
   input_fields(field)%units = trim(units)
else
   input_fields(field)%units = 'none'
endif
if(present(missing_value)) then
! ??? What about passing missing value through to the output_fields???
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

! Next need to compute the size of the domain for this set of axes
! Comes from axis handler; assume result is i_size, j_size, k_size
! for three-d problem
siz = 1; local_siz = 1
local_start = 1;  local_end= 1
do j = 1, num_axes
   siz(j) = get_axis_length(axes(j))
end do

! Default length for axes is 1
do j = 1, 3
   input_fields(field)%size(j) = siz(j)
end do

! Need to loop through all output_fields associated and allocate their buffers
! Would have to consider horizontal operations at this point
do j = 1, input_fields(field)%num_output_fields
   out_num = input_fields(field)%output_fields(j)
!  if local_output (size of output_fields does NOT equal size of input_fields)
   if(output_fields(out_num)%local_output) then
      if(size(axes).le.1) &
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
         output_fields(out_num)%total_elements = local_siz(1)*local_siz(2)*local_siz(3)
      endif
      call mpp_sync()     
   else ! the field is outputted globally
! size of output_fields equal size of input_fields 
      allocate(output_fields(out_num)%buffer(siz(1), siz(2), siz(3)))
      output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
   endif
   output_fields(out_num)%buffer = EMPTY
! Reset to false in register_field if this is not static
   output_fields(out_num)%static = .true.
! check if time average is true for static field
   if(.not.dynamic1 .and. output_fields(out_num)%time_average) then      
      write(msg,'(a,"/",a)')trim(module_name), trim(field_name)
      if(mpp_pe() .eq. mpp_root_pe()) &
           call  error_mesg ('register_static_field', 'module/field '//trim(msg)//' is STATIC, can NOT time average',WARNING)
      output_fields(out_num)%time_average = .false.
   endif
! Axes are copied from input_fields if outputted globally or from subaxes if outputted locally
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

!initilization for local output
start = -1.e10; end=-1.e10 ! initially out of (lat/lon/depth) range
gstart_indx = -1; gend_indx=-1

! get axis data (lat, lon, depth) and indices
   start= output_fields(outnum)%output_grid%start
   end = output_fields(outnum)%output_grid%end

do i = 1,size(axes)   
   global_axis_size = get_axis_global_length(axes(i))
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
      allocate(global_depth(global_axis_size))
      call get_diag_axis_data(axes(i),global_depth)
      gstart_indx(i) = get_index(start(i),global_depth)
      gend_indx(i) = get_index(end(i),global_depth)
      allocate(subaxis_z(gstart_indx(i):gend_indx(i)))
      subaxis_z=global_depth(gstart_indx(i):gend_indx(i))
      output_fields(outnum)%output_grid%subaxes(i) = &
           diag_subaxes_init(axes(i),subaxis_z, gstart_indx(i),gend_indx(i))
      deallocate(subaxis_z,global_depth)
   case default
       call error_mesg ('diag_manager, get_subfield_size', 'Wrong axis_cart', FATAL)
   end select
enddo
do i = 1,size(axes)
   if(gstart_indx(i)== -1 .or. gend_indx(i)== -1) &
        call error_mesg ('diag_manager, get_subfield_size', 'can not find gstart_indx/gend_indx for ' &
        //trim(output_fields(outnum)%output_name), FATAL)  
enddo
if(size(axes)>2 .and. (gstart_indx(3)== -1 .or. gend_indx(3) ==-1)) &
     call error_mesg('diag_manager, get_subfield_size', 'can not find local depth in Z axis',FATAL)

! get domain and compute_domain(xbegin,xend,ybegin,yend)
xbegin=-1; xend=-1
ybegin=-1; yend=-1

Domain2 = get_domain2d(axes)
if(Domain2 /= NULL_DOMAIN2D) then
   call mpp_get_compute_domain(Domain2,xbegin,xend,ybegin,yend)
   call mpp_get_domain_components(Domain2, Domain1x, Domain1y)
else
   do i = 1, MIN(size(axes),2)    
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
   call mpp_modify_domain(Domain1x,Domain1x_new, xbegin_l,xend_l,gstart_indx(1),gend_indx(1))
   call mpp_modify_domain(Domain1y,Domain1y_new, ybegin_l,yend_l,gstart_indx(2),gend_indx(2))
   call mpp_modify_domain(Domain2_new, Domain1x_new, Domain1y_new)
   output_fields(outnum)%output_grid%subaxes(1) = &
        diag_subaxes_init(axes(1),subaxis_x, gstart_indx(1),gend_indx(1),Domain1x_new,Domain2_new)
   output_fields(outnum)%output_grid%subaxes(2) = &
        diag_subaxes_init(axes(2),subaxis_y, gstart_indx(2),gend_indx(2),Domain1y_new,Domain2_new)
   if(output_fields(outnum)%output_grid%subaxes(1) == -1 .or. output_fields(outnum)%output_grid%subaxes(2) == -1)&
        call error_mesg ('diag_manager, get_subfield_size', 'wrong subaxis id', FATAL) 
! local start index should start from 1
   output_fields(outnum)%output_grid%l_start_indx(1) = MAX(xbegin, gstart_indx(1)) - xbegin + 1   
   output_fields(outnum)%output_grid%l_start_indx(2) = MAX(ybegin, gstart_indx(2)) - ybegin + 1
   output_fields(outnum)%output_grid%l_end_indx(1) = MIN(xend, gend_indx(1)) - xbegin + 1 
   output_fields(outnum)%output_grid%l_end_indx(2) = MIN(yend, gend_indx(2)) - ybegin + 1
   if(size(axes)>2) then
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

  n = size(array)
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
type (time_type), intent(in) :: time
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

function send_data_1d(diag_field_id, field, time, is_in, mask, rmask, weight)

logical                      :: send_data_1d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:)
type (time_type), intent(in) :: time
integer, optional            :: is_in
logical, optional            :: mask(:)
real, optional               :: rmask(:), weight(:)
real                         :: field_out(size(field), 1, 1)
real                         :: weight_out(size(field), 1, 1)
logical                      :: mask_out(size(field), 1, 1)


! First copy the data to a three d array with last element 1
field_out(:, 1, 1) = field

! Default values for mask and weight
mask_out = .true.
weight_out = 1.0

if(present(mask)) mask_out(:, 1, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, 1, 1) = .false.
if(present(weight)) weight_out(:, 1, 1) = weight
if(present(mask) .or. present(rmask)) then
   send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in, 1, 1, mask_out, weight_out)
else
   send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in, 1, 1)
endif
end function send_data_1d


!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_2d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_2d(diag_field_id, field, time, is_in, js_in, &
    mask, rmask, weight)

logical                      :: send_data_2d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:, :)
type (time_type), intent(in) :: time
integer, optional            :: is_in, js_in
logical, optional            :: mask(:, :)
real, optional               :: rmask(:, :), weight(:, :)
real                         :: field_out(size(field, 1), size(field, 2), 1)
real                         :: weight_out(size(field, 1), size(field, 2), 1)
logical                      :: mask_out(size(field, 1), size(field, 2), 1)


! First copy the data to a three d array with last element 1
field_out(:, :, 1) = field

! Default values for mask and weight
mask_out = .true.
weight_out = 1.0

if(present(mask)) mask_out(:, :, 1) = mask
if(present(rmask)) where (rmask < 0.5) mask_out(:, :, 1) = .false.
if(present(weight)) weight_out(:, :, 1) = weight

if(present(mask) .or. present(rmask)) then
   send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in, js_in, 1, mask_out, weight_out)
else
   send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in, js_in, 1)
endif

end function send_data_2d

!-------------------------------------------------------------------------
! <FUNCTION NAME="send_data_3d" INTERFACE="send_data">
!   <IN NAME="diag_field_id" TYPE="integer"  > </IN>
!   <IN NAME="field" TYPE="real" DIM="(:,:,:)" > </IN>
!   <IN NAME="time" TYPE="time_type"  > </IN>
! </FUNCTION>

function send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             mask, rmask, weight)

logical                      :: send_data_3d
integer, intent(in)          :: diag_field_id
real, intent(in)             :: field(:, :, :)
type (time_type), intent(in) :: time
integer, optional            :: is_in, js_in, ks_in
logical, target, optional    :: mask(:, :, :)
real, target, optional       :: rmask(:, :, :)
real, optional               :: weight(:, :, :) ! not implemented at this time
real                         :: num
logical                      :: average, need_compute, local_output
integer                      :: i, out_num, file_num, n1, n2, n3
integer                      :: freq, units, is, js, ks, kount
character(len=128)           :: error_string
integer                      :: local_start(3), local_end(3) ! start and end indices on 3 axes
logical, pointer             :: mask_local(:,:,:)
real, pointer                :: rmask_local(:,:,:)
integer                      :: number_of_outputs

if (.not.module_is_initialized) call error_mesg ('send_data_3d',' diag_manager NOT initialized', FATAL)
! If is, js, or ks not present default them to 1
is = 1
if(present(is_in)) is = is_in
js = 1
if(present(js_in)) js = js_in
ks = 1
if(present(ks_in)) ks = ks_in

if(first_send_data_call .and. .not. input_fields(diag_field_id)%static) then
   call output_setup()
   first_send_data_call = .false.
! Also set the time last flushed for all files to this time
   do i = 1, num_files
      files(i)%last_flush = time
   end do
endif
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
   if(first_send_data_local) then
      call output_setup(.true.)
      first_send_data_local = .false.
   endif
else
   number_of_outputs = input_fields(diag_field_id)%num_output_fields
endif

! Loop through each output field that depends on this input field
do i = 1, number_of_outputs
   if(i == input_fields(diag_field_id)%num_output_fields + 1) then
      call mpp_sync()
      exit
   endif
! Get index to an output field
   out_num = input_fields(diag_field_id)%output_fields(i)

! is this field outputted on a local domain only?
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
   
   if(need_compute) then        
      local_start = output_fields(out_num)%output_grid%l_start_indx
      local_end = output_fields(out_num)%output_grid%l_end_indx
      if(present(mask)) then
         mask_local => mask(local_start(1):local_end(1),local_start(2):local_end(2),local_start(3):local_end(3))
      endif
      if(present(rmask)) then
         rmask_local => rmask(local_start(1):local_end(1),local_start(2):local_end(2),local_start(3):local_end(3))
      endif
   endif
   
! Initialize output time for fields output every time step
! this will guarantee that first time here field will not be written
   if (freq == EVERY_TIME) then
     if (output_fields(out_num)%next_output == output_fields(out_num)%last_output) &
         output_fields(out_num)%next_output = time
   endif

! Is it time to output for this field; CAREFUL ABOUT > vs >= HERE
   if(.not.output_fields(out_num)%static .and. freq /= END_OF_RUN .and. &
        time > output_fields(out_num)%next_output ) then
! A non-static field that has skipped a time level is an error
      if(time >  output_fields(out_num)%next_next_output .and. freq > 0) then
         if(mpp_pe() .eq. mpp_root_pe()) &
              call error_mesg('warning1, diag_manager', &
              'skip one time level in output data', WARNING)
      end if
! If average get size: Average intervals are last_output, next_output
      if(average) then
         if (input_fields(diag_field_id)%mask_variant) then
            if (any(output_fields(out_num)%counter>0.)) then
               where(output_fields(out_num)%counter>0.) &
                    output_fields(out_num)%buffer = output_fields(out_num)%buffer/output_fields(out_num)%counter 
            else
               if(any(output_fields(out_num)%buffer /= input_fields(diag_field_id)%missing_value)) then
                  write (error_string,'(a,"/",a)')  &
                       trim(input_fields(diag_field_id)%module_name), &
                       trim(output_fields(out_num)%output_name)
                  call error_mesg ('error2 diag_manager_mod', &
                       'module/output_field '//trim(error_string)//&
                       &', write EMPTY buffer', FATAL) 
               endif
            endif
         else
            num = output_fields(out_num)%count_0d
            if(output_fields(out_num)%total_elements > size(field)) then
               num = real(output_fields(out_num)%num_elements/output_fields(out_num)%total_elements)
               if(local_output) then
                  write (error_string,'(a,"/",a)')  &
                       trim(input_fields(diag_field_id)%module_name), &
                       trim(output_fields(out_num)%output_name)
                  call error_mesg ('diag_manager_mod, send_data', &
                       'module/output_field '//trim(error_string)//&
                       &', does NOT work with regional output', FATAL)
               endif
            endif
            if (num>0.) then
               if (input_fields(diag_field_id)%missing_value_present) then
                  where (output_fields(out_num)%buffer /= input_fields(diag_field_id)%missing_value) &
                       output_fields(out_num)%buffer = output_fields(out_num)%buffer/num
               else
                  output_fields(out_num)%buffer = output_fields(out_num)%buffer/num
               endif
            else
               if(any(output_fields(out_num)%buffer /= input_fields(diag_field_id)%missing_value)) then
                  write (error_string,'(a,"/",a)')  &
                       trim(input_fields(diag_field_id)%module_name), &
                       trim(output_fields(out_num)%output_name)
                  call error_mesg ('error3 diag_manager_mod', &
                       'module/output_field '//trim(error_string)//&
                       &', write EMPTY buffer', FATAL)
               endif
            endif
         endif
      endif !average
! Output field
      call diag_data_out(file_num, out_num, &
           output_fields(out_num)%buffer, output_fields(out_num)%next_output)
      
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
      output_fields(out_num)%buffer = EMPTY
      output_fields(out_num)%num_elements = 0
      if(input_fields(diag_field_id)%mask_variant .and. average) output_fields(out_num)%counter = 0.0
   end if  !time > output_fields(out_num)%next_output

! Finished output of previously buffered data, now deal with buffering new data   
   n1 = size(field, 1); n2 = size(field, 2); n3 = size(field, 3) ! need the size of buffer,not field
! Take care of submitted field data
   if(average) then
      if (input_fields(diag_field_id)%mask_variant) then
         if(need_compute) then
            write (error_string,'(a,"/",a)')  &
                 trim(input_fields(diag_field_id)%module_name), &
                 trim(output_fields(out_num)%output_name)   
            call error_mesg ('error4 diag_manager_mod', &
                 'module/output_field '//trim(error_string)//&
                 &', local output NOT supported with mask_variant', FATAL)
         endif
         if (present(mask)) then
            if (input_fields(diag_field_id)%missing_value_present) then            
               where (mask)                  
                  output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = &
                       output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) + field
                  output_fields(out_num)%counter = output_fields(out_num)%counter +1.                                      
               elsewhere
                  output_fields(out_num)%buffer(is:is+n1-1, js:js+n2-1, ks:ks+n3-1) = &
                       input_fields(diag_field_id)%missing_value
               endwhere
            else
               write (error_string,'(a,"/",a)')  &
                    trim(input_fields(diag_field_id)%module_name), &
                    trim(output_fields(out_num)%output_name)
               call error_mesg ('error5 diag_manager_mod', &
                    'module/output_field '//trim(error_string)//&
                    &', variable mask but no missing value defined', FATAL)
            endif
         else
            write (error_string,'(a,"/",a)')  &
                 trim(input_fields(diag_field_id)%module_name), &
                 trim(output_fields(out_num)%output_name)
            call error_mesg ('error6 diag_manager_mod', &
                 'module/output_field '//trim(error_string)//&
                 &', variable mask but no mask given', FATAL)
         endif         
      else ! mask_variant=false
         if (present(mask)) then
            if (input_fields(diag_field_id)%missing_value_present) then
               if(need_compute) then
                  where(mask_local)
                     output_fields(out_num)%buffer = output_fields(out_num)%buffer + &
                          field(local_start(1):local_end(1),local_start(2):local_end(2),local_start(3):local_end(3))
                  elsewhere
                     output_fields(out_num)%buffer = input_fields(diag_field_id)%missing_value
                  endwhere
               else   
                  where(mask)                                           
                     output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = &
                          output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) + field                  
                  elsewhere                  
                     output_fields(out_num)%buffer(is:is+n1-1, js:js+n2-1, ks:ks+n3-1) = &
                          input_fields(diag_field_id)%missing_value                    
                  endwhere
               endif
               if(need_compute) then
                  if(any(mask_local)) output_fields(out_num)%count_0d=output_fields(out_num)%count_0d + 1.
               else
                  if(any(mask)) output_fields(out_num)%count_0d=output_fields(out_num)%count_0d + 1.
               endif               
            else ! missing value NOT present
               if(.not.all(mask) .and. mpp_pe() .eq. mpp_root_pe()) &
                    call error_mesg('warning2 send_data_3d', &
                    'Mask will be ignored since missing values were not specified',WARNING)
               if(need_compute) then
                  output_fields(out_num)%buffer = output_fields(out_num)%buffer + &
                       field(local_start(1):local_end(1),local_start(2):local_end(2),local_start(3):local_end(3))
               else                          
                  output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = &
                       output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) + field   
               endif
               output_fields(out_num)%count_0d = output_fields(out_num)%count_0d +1.
            endif
         else ! mask NOT present
            if (input_fields(diag_field_id)%missing_value_present) then
               if(need_compute) then                  
                  where(field /= input_fields(diag_field_id)%missing_value)                  
                     output_fields(out_num)%buffer = output_fields(out_num)%buffer + &
                          field(local_start(1):local_end(1),local_start(2):local_end(2),local_start(3):local_end(3))
                  elsewhere
                     output_fields(out_num)%buffer = input_fields(diag_field_id)%missing_value
                  endwhere
               else
                  where(field /= input_fields(diag_field_id)%missing_value)
                     output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = &
                          output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) + field                  
                  elsewhere                   
                     output_fields(out_num)%buffer(is:is+n1-1, js:js+n2-1, ks:ks+n3-1) = &
                          input_fields(diag_field_id)%missing_value                  
                  endwhere
               endif
               if(any(field /= input_fields(diag_field_id)%missing_value)) &
                    output_fields(out_num)%count_0d = output_fields(out_num)%count_0d + 1 
            else
               if(need_compute) then
                  output_fields(out_num)%buffer = output_fields(out_num)%buffer + &
                       field(local_start(1):local_end(1),local_start(2):local_end(2),local_start(3):local_end(3))
               else 
                  output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = &
                       output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) + field
               endif
               output_fields(out_num)%count_0d = output_fields(out_num)%count_0d +1.
            endif
         endif
      endif
      output_fields(out_num)%num_elements = output_fields(out_num)%num_elements + n1*n2*n3
   else  ! ( not average)
      if(need_compute) then
         output_fields(out_num)%buffer =  &
              field(local_start(1):local_end(1),local_start(2):local_end(2),local_start(3):local_end(3))
      else
         output_fields(out_num)%buffer(is:is+n1-1,js: js+n2-1, ks:ks+n3-1) = field
      endif
               
      if (present(mask) .and. input_fields(diag_field_id)%missing_value_present) then
         if(need_compute) then
            where(.not. mask_local)
               output_fields(out_num)%buffer = input_fields(diag_field_id)%missing_value
            endwhere
         else
            where (.not. mask)             
               output_fields(out_num)%buffer(is:is+n1-1, js:js+n2-1, ks:ks+n3-1) = &
                    input_fields(diag_field_id)%missing_value           
            endwhere
         endif
      endif
   endif !average
 
! If rmask and missing value present, then insert missing value     
   if (present(rmask) .and. input_fields(diag_field_id)%missing_value_present) then
      if(need_compute) then
         where (rmask_local < 0.5)
            output_fields(out_num)%buffer = input_fields(diag_field_id)%missing_value
         endwhere
      else
         where (rmask < 0.5)         
            output_fields(out_num)%buffer(is:is+n1-1, js:js+n2-1, ks:ks+n3-1) = &
                 input_fields(diag_field_id)%missing_value        
         endwhere
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

!   <OVERVIEW>
!     Exit Diagnostics Manager.
!   </OVERVIEW>
!   <DESCRIPTION>
!     Flushes diagnostic buffers where necessary. Close diagnostics files.<BR/>
! A warning will be issued here if a field in diag_table is not registered
!   </DESCRIPTION>
!   <TEMPLATE>
!     call diag_manager_end (time)
!   </TEMPLATE>
!   <IN NAME="TIME" TYPE="time_type"></IN>

subroutine diag_manager_end (time)

! PROBABLY JUST WRITE STATIC FIELDS HERE FOR NOW???

type(time_type), intent(in) :: time

integer            :: i, file_num, freq, j, input_num
integer            :: out_num
logical            :: average, local_output, need_compute
character(len=128) :: message
real               :: num, num1, num2

if(.not.module_is_initialized) return

do i = 1, num_output_fields

! is this field outputted on a local domain only?
   local_output = output_fields(i)%local_output
! if local_output, does the current PE take part in send_data?
   need_compute = output_fields(i)%need_compute
! skip all PEs not participating in outputting this field
   if(local_output .and. (.not. need_compute)) cycle

! only output non-static fields here
   input_num = output_fields(i)%input_field
   if (input_fields(input_num)%static) cycle
! skip fields that were not registered
   if (.not.input_fields(input_num)%register) cycle

! Get index to output file for this field
   file_num = output_fields(i)%output_file
   if(file_num == max_files) cycle
! Output frequency for this file is
   freq = files(file_num)%output_freq

! Is it time to output for this field; CAREFUL ABOUT >= vs > HERE
! For end should be >= because no more data is coming 
   if(time >= output_fields(i)%next_output .or. freq == END_OF_RUN) then
      if(time >= output_fields(i)%next_next_output .and. freq > 0) then
         write (message,'(a,"/",a)') trim(input_fields(input_num)%module_name), &
              trim(output_fields(i)%output_name)
         if(mpp_pe() .eq. mpp_root_pe()) & 
              call error_mesg('warning4 diag_manager_end', 'module/output_field ' //  &
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
                  call error_mesg ('error7 diag_manager_mod', &
                       'module/output_field '//trim(message)//&
                       &', write EMPTY buffer', FATAL) 
               endif
            endif
         else
            
            num1 = output_fields(i)%count_0d
            num2 = real(output_fields(i)%num_elements/output_fields(i)%total_elements)
            num = MIN(num1,num2)
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
      call diag_data_out(file_num, i, &
            output_fields(i)%buffer, output_fields(i)%next_output, .true.)
   end if 
end do

! If there are no non-static fields then need to setup output data
if(first_send_data_call) then
   call output_setup()
   first_send_data_call = .false.
endif

! Output static fields and signal end of output for each file
do i = 1, num_files
! Loop to look for static fields
   do j = 1, files(i)%num_fields
      out_num = files(i)%fields(j)
      input_num = output_fields(out_num)%input_field

! skip fields that were not registered
      if (.not.input_fields(input_num)%register) cycle
! only output static fields here
      if (.not.output_fields(out_num)%static) cycle

      call diag_data_out(i, out_num, output_fields(out_num)%buffer, time, .true.)
   end do   
! Close up this file   
   call diag_output_end(files(i)%file_unit)
end do

module_is_initialized = .FALSE.

end subroutine diag_manager_end
! </SUBROUTINE>

!-------------------------------------------------------------------------

subroutine init_file(name, output_freq, output_units, format, time_units, long_name)

character(len=*), intent(in) :: name, long_name
integer, intent(in)          :: output_freq, output_units, format, time_units

! Get a number for this file
num_files = num_files + 1
if(num_files >= max_files) then
   call error_mesg('diag_manager, init_file', ' max_files exceeded, incease max_files', FATAL)
endif

files(num_files)%name = trim(name)
files(num_files)%output_freq = output_freq
files(num_files)%output_units = output_units
files(num_files)%format = format
files(num_files)%time_units = time_units
files(num_files)%long_name = trim(long_name)
files(num_files)%num_fields = 0
files(num_files)%local = .false.
! This value should be updated with first call to send_data for non-static field
files(num_files)%last_flush = set_time(0, 0)

end subroutine init_file

!--------------------------------------------------------------------------

subroutine init_input_field(module_name, field_name)

character(len=*), intent(in) :: module_name, field_name

! Get a number for this input_field if not already set up
if(find_input_field(module_name, field_name) < 0) then
   num_input_fields = num_input_fields + 1
   if(num_input_fields > max_input_fields) then
      call error_mesg('diag_manager,init_input_field', 'max_input_fields exceeded, increase max_input_fields', FATAL)
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
end subroutine init_input_field

!---------------------------------------------------------------------------

subroutine init_output_field(module_name, field_name, output_name, output_file,&
   time_average, pack, local_coord)

character(len=*), intent(in)           :: module_name, field_name, output_name, output_file
logical, intent(in)                    :: time_average
integer, intent(in)                    :: pack
type(coord_type), intent(in), optional :: local_coord
integer                                :: out_num, in_num, file_num, num_fields, i
character(len=128)                     :: error_msg

! Get a number for this output field
num_output_fields = num_output_fields + 1
if(num_output_fields > max_output_fields) then
   call error_mesg('diag_manager,init_output_field', 'max_output_fields exceeded, increase max_output_fields', FATAL)
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

! Set the file for this output field
output_fields(out_num)%output_file = file_num

! Enter the other data for this output field
output_fields(out_num)%output_name = trim(output_name)
output_fields(out_num)%pack = pack
output_fields(out_num)%num_elements = 0
output_fields(out_num)%total_elements = 0

! cannot time average fields output every time
if (files(file_num)%output_freq == EVERY_TIME) then
  output_fields(out_num)%time_average = .false.
else
  output_fields(out_num)%time_average = time_average
endif

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
   logical            :: time_avg
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
end type tableA_type

character(len=256) :: record
character(len=9)   :: amonth

integer :: iunit,n,m,num_fields,time_units, output_freq_units, nfiles,nfields
integer :: j, log_unit, name_len, nrecs, ierr, io_status
integer, allocatable, dimension(:) :: pelist
logical :: append_pelist_name = .false.

namelist /diag_manager_nml/ append_pelist_name

type(tableB_type) :: textB
type(tableA_type) :: textA
type(coord_type) :: local_coord !local coordinates used in local output

!  If the module was already initialized do nothing
if (module_is_initialized) return

diag_subset_output = DIAG_ALL
if (PRESENT(diag_model_subset))then
  if(diag_model_subset>=DIAG_OTHER .AND. diag_model_subset<=DIAG_ALL) then
    diag_subset_output = diag_model_subset
  else
    call error_mesg('diag_manager_init','file diag_table nonexistent',FATAL)
  endif
endif

iunit = open_namelist_file()
read  (iunit, diag_manager_nml,iostat=io_status)
write (stdlog(), diag_manager_nml)
!ierr = check_nml_error(io_status,'diag_manager_nml')
call close_file (iunit)
    
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
do while (nfiles <= max_files)
   read(iunit,'(a)',end=86,err=85) record
   nrecs=nrecs+1
   if (record(1:1) == '#') cycle
   read(record,*,err=85,end=85) textA
   ! test file format to make sure its OK
   if (textA%format .gt. 2 .or. textA%format .lt. 1) cycle
   if( diag_subset_output==DIAG_OTHER .AND. verify( 'ocean',lcase(textA%name) )==0 )cycle
   if( diag_subset_output==DIAG_OCEAN .AND. verify( 'ocean',lcase(textA%name) )/=0 )cycle
   nfiles=nfiles+1
   time_units = 0
   output_freq_units = 0
   do j = 1, size(time_unit_list)
      if(textA%time_units == time_unit_list(j)) time_units = j
      if(textA%output_freq_units == time_unit_list(j)) output_freq_units = j
   end do
   if(time_units == 0) &
        call error_mesg('diag_manager_init','invalid time units, check time unit in diag_table',FATAL)
   if(output_freq_units == 0) & 
        call error_mesg('diag_manager_init','invalid output frequency units, check diag table',FATAL)
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
   call init_file(textA%name,textA%output_freq, output_freq_units, &
        textA%format, time_units,textA%long_name)
85 continue
enddo
call error_mesg('diag_manager_init','too many files in diag_table, increase max_files', FATAL)
86 continue



rewind(iunit)
!if (nfiles .lt. 1) call error_mesg('diag_manager_init','error reading file records',FATAL)

nfields=0;nrecs=0
do while (nfields <= max_output_fields)
   read(iunit,'(a)',end=94,err=93) record
   nrecs=nrecs+1
   if (record(1:1) == '#') cycle
   read(record,*,end=93,err=93) textB
   if (textB%pack .gt. 8 .or. textB%pack .lt. 1) cycle
   if( diag_subset_output==DIAG_OTHER .AND. verify( 'ocean',lcase(textB%name) )==0 )cycle
   if( diag_subset_output==DIAG_OCEAN .AND. verify( 'ocean',lcase(textB%name) )/=0 )cycle
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
           textB%name,textB%time_avg,textB%pack)
   else
      call init_output_field(textB%module_name,textB%field_name,textB%output_name,&
           textB%name,textB%time_avg,textB%pack, local_coord)
   endif
93 continue
enddo
call error_mesg('diag_manager_init','too many fields in table, increase max_output_fields', FATAL)
94 continue
call close_file(iunit)

! check duplicate output_fields in the diag_table
call check_duplicate_output_fields

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

integer :: find_file
character(len=*), intent(in) :: name

integer :: i

find_file = -1
do i = 1, num_files
   if(files(i)%name == name) then
      find_file = i
      return
   end if
end do

end function find_file

!-------------------------------------------------------------------------

function find_input_field(module_name, field_name)

integer find_input_field
character(len=*), intent(in) :: module_name, field_name

integer :: i

find_input_field = -1
do i = 1, num_input_fields
   if(trim(input_fields(i)%module_name) == trim(module_name) .and. &
      lcase(trim(input_fields(i)%field_name)) == &
      lcase(trim(field_name))) then 
      find_input_field = i
      return
   endif
end do

end function find_input_field

!-------------------------------------------------------------------------

subroutine output_setup(local1)

! WARNING: Assumes that all data structures are fully initialized

integer                       ::i, j, field_num, n_fields, axes(4), input_field_num, num_axes, k
character(len=128)            ::time_units, avg, error_string, filename
logical                       :: file_time_avg, local
integer                       :: time_axis_id(1)
logical, intent(in), optional :: local1
character(len=7)              :: prefix

local = .false.
if(present(local1)) local = local1
! Set up for output, focused on netcdf

! First, get a file_unit and a time axis identifier for each file
outer:do i = 1, num_files
   if(.not.local) then
      do j = 1, files(i)%num_fields
         field_num = files(i)%fields(j)
         if(output_fields(field_num)%local_output)cycle outer
      enddo
! it's unlikely that a file starts with word "rregion", need to check anyway.
      if (len(files(i)%name) >=7) then
         prefix = files(i)%name(1:7)
         if(lowercase(prefix) == 'rregion') &
              call error_mesg ('diag_manager output_setup', 'file name should not start with' &
              //' word "rregion"', WARNING)
      endif
   else      
      if(.not.files(i)%local) cycle outer     
   endif   
   file_time_avg = .false.
! Skip this file if no fields are to be output
   if (files(i)%num_fields == 0) cycle
! Here is where time_units string must be set up; time since base date
   write(time_units, 11) trim(time_unit_list(files(i)%time_units)), base_year, &
      base_month, base_day, base_hour, base_minute, base_second
 11 format(a, ' since ', i4.4, '-', i2.2, '-', i2.2, ' ', i2.2, ':', i2.2, ':', i2.2) 
   if(files(i)%local) then      
! prefix "rregion" to all local files for post processing, the prefix will be removed in postprocessing
         filename = 'rregion'//trim(files(i)%name)
   else
      filename = trim(files(i)%name)
   endif
   call diag_output_init(filename, files(i)%format, global_descriptor, &
      files(i)%long_name, time_units, files(i)%file_unit, files(i)%time_axis_id)

! Loop through all fields with this file to output axes
   do j = 1, files(i)%num_fields
      field_num = files(i)%fields(j)
      input_field_num = output_fields(field_num)%input_field
      if (.not.input_fields(input_field_num)%register) then
         write (error_string,'(a,"/",a)')  &
              trim(input_fields(input_field_num)%module_name), &
              trim(input_fields(input_field_num)%field_name)
         if(mpp_pe() .eq. mpp_root_pe()) &
              call error_mesg ('diag_manager output_setup', &
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
            call error_mesg ('diag_manager output_setup','output_name '//trim(error_string)// &
                 ' has axis_id = -1', FATAL)
         endif
      enddo
      axes(num_axes + 1) = files(i)%time_axis_id
      call write_axis_meta_data(files(i)%file_unit, axes(1:num_axes + 1))
   end do

! Now output metadata for each of the fields associated with this file
   do j = 1, files(i)%num_fields
      field_num = files(i)%fields(j)
      input_field_num = output_fields(field_num)%input_field
      if (.not.input_fields(input_field_num)%register) cycle
! Put the time axis in the axis field
      num_axes = output_fields(field_num)%num_axes
      axes(1:num_axes) = output_fields(field_num)%axes(1:num_axes)
      if (.not. output_fields(field_num)%static) then
        num_axes=num_axes+1
        axes(num_axes) = files(i)%time_axis_id
      endif
      if(output_fields(field_num)%time_average) then
         avg = avg_name
         file_time_avg = .true.
      else
         avg = ""
      end if
      if(input_fields(input_field_num)%missing_value_present) then
         output_fields(field_num)%f_type = write_field_meta_data(files(i)%file_unit, &
            output_fields(field_num)%output_name, axes(1:num_axes), &
            input_fields(input_field_num)%units, &
            input_fields(input_field_num)%long_name, &
            input_fields(input_field_num)%range, output_fields(field_num)%pack,&
            input_fields(input_field_num)%missing_value, avg_name = avg)

! NEED TO TAKE CARE OF TIME AVERAGING INFO TOO BOTH CASES
      else
         output_fields(field_num)%f_type = write_field_meta_data(files(i)%file_unit, &
            output_fields(field_num)%output_name, axes(1:num_axes), &
            input_fields(input_field_num)%units, &
            input_fields(input_field_num)%long_name, &
            input_fields(input_field_num)%range, output_fields(field_num)%pack,&
               avg_name = avg)
      endif
   end do

! If any of the fields in the file are time averaged, need to output the axes
! Use double precision since time axis is double precision
   if(file_time_avg) then
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
   end if

! Let lower levels know that all meta data has been sent
   call done_meta_data(files(i)%file_unit)

end do outer

! ALSO NEED TO SEND DATA FOR ANY STATIC FIELD THAT HAS ALREADY BEEN WRITTEN

end subroutine output_setup

!-------------------------------------------------------------------------

subroutine diag_data_out(file, field, dat, time, final_call_in)

integer, intent(in) :: file, field
real, intent(inout) :: dat(:, :, :)
type(time_type), intent(in) :: time
logical, optional, intent(in) :: final_call_in
logical :: final_call
integer :: i, num
real :: dif, time_data(1, 1, 1), dt_time(1, 1, 1), start_dif

final_call = .false.
if(present(final_call_in)) final_call = final_call_in

dif = get_date_dif(time, base_time, files(file)%time_units)
call diag_field_out(files(file)%file_unit,output_fields(field)%f_type, dat, dif)

! *** inserted this line because start_dif < 0 for static fields ***
if (output_fields(field)%static) return

start_dif = get_date_dif(output_fields(field)%last_output, base_time, &
   files(file)%time_units)

! Need to write average axes out;
do i = 1, files(file)%num_fields
   num = files(file)%fields(i)
   if(output_fields(num)%time_average .and. .not.output_fields(num)%static .and. &
      input_fields(output_fields(num)%input_field)%register) then
      if(num == field) then
! Output the axes if this is first time-averaged field
         time_data(1, 1, 1) = start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_start, &
            time_data, dif)
         time_data(1, 1, 1) = dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_end, &
            time_data, dif)
! Compute the length of the average
         dt_time(1, 1, 1) = dif - start_dif
         call diag_field_out(files(file)%file_unit, files(file)%f_avg_nitems, &
            dt_time, dif)
      endif
      goto 10
   end if
end do

! If write time is greater (equal for the last call) than last_flush for this file, flush it
10 if(final_call) then
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

function get_date_dif(t2, t1, units)

real :: get_date_dif
type(time_type), intent(in) :: t2, t1
integer, intent(in) :: units

integer :: year, month, day, hour, minute, second, dif_seconds, dif_days
type(time_type) :: dif_time

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

subroutine print_file(i)

integer, intent(in) :: i
integer :: j

write(*, *)'----------------------------------------------------'
write(*, *) 'contents of file entry ', i
write(*, *) 'name'
write(*, *) trim(files(i)%name)
write(*, *) 'output_freq ', files(i)%output_freq
write(*, *) 'format ', files(i)%format
write(*, *) 'time units '
write(*, *) time_unit_list(files(i)%time_units)
write(*, *) 'long_name '
write(*, *) trim(files(i)%long_name)
write(*, *) 'num fields ', files(i)%num_fields
do j = 1, files(i)%num_fields
   write(*, *) 'field ', j, ' is ', files(i)%fields(j)
end do
write(*, *)'----------------------------------------------------'

end subroutine print_file

!------------------------------------------------------------------------

subroutine print_input_field(i)

integer, intent(in) :: i
integer :: j

write(*, *)'----------------------------------------------------'
write(*, *) 'contents of field entry ', i
write(*, *) 'module and field name'
write(*, *) trim(input_fields(i)%module_name)
write(*, *) trim(input_fields(i)%field_name)
write(*, *) 'long_name '
write(*, *) trim(input_fields(i)%long_name)
write(*, *) 'units '
write(*, *) trim(input_fields(i)%units)
write(*, *) 'axes ', input_fields(i)%axes(1), input_fields(i)%axes(2), &
   input_fields(i)%axes(3)
write(*, *) 'num_axes ', input_fields(i)%num_axes
write(*, *) 'missing value and range present ', &
   input_fields(i)%missing_value_present, input_fields(i)%range_present
write(*, *) 'missing value and range ', input_fields(i)%missing_value, &
   input_fields(i)%range(1), input_fields(i)%range(2)
write(*, *) 'sizes and total ', input_fields(i)%size(1), &
   input_fields(i)%size(2), input_fields(i)%size(3)
!   input_fields(i)%total_elements
write(*, *) 'static ', input_fields(i)%static
write(*, *) 'num_output_fields ', input_fields(i)%num_output_fields
do j = 1, input_fields(i)%num_output_fields
   write(*, *) 'output field ', j, ' is ', input_fields(i)%output_fields(j)
end do
write(*, *)'----------------------------------------------------'

end subroutine print_input_field

!-------------------------------------------------------------------------

subroutine print_output_field(i)

integer, intent(in) :: i
integer :: seconds, days

write(*, *)'----------------------------------------------------'
write(*, *) 'contents of output field ', i
write(*, *) 'input field and output file ', output_fields(i)%input_field, &
   output_fields(i)%output_file
write(*, *) 'output name ', trim(output_fields(i)%output_name)
write(*, *) 'time average ', output_fields(i)%time_average
write(*, *) 'pack ', output_fields(i)%pack
!write(*, *) 'total_elements and num_elements ', &
!   output_fields(i)%total_elements, output_fields(i)%num_elements
write(*, *) 'static ', output_fields(i)%static
write(*, *) 'num_axes ', output_fields(i)%num_axes
write(*, *) 'axes ', output_fields(i)%axes
call get_time(output_fields(i)%last_output, seconds, days)
write(*, *) 'last output ', seconds, ' seconds ', days, ' days '
call get_time(output_fields(i)%next_output, seconds, days)
write(*, *) 'next output ', seconds, ' seconds ', days, ' days '
call get_time(output_fields(i)%next_next_output, seconds, days)
write(*, *) 'next_next_output ', seconds, ' seconds ', days, ' days '
write(*, *)'----------------------------------------------------'

end subroutine print_output_field

!-------------------------------------------------------------------------

function diag_time_inc(time, output_freq, output_units)

type (time_type) :: diag_time_inc
type (time_type), intent(in) :: time
integer, intent(in) :: output_freq, output_units

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
         'output units of months not allowed with no calendar', FATAL)
   else
      diag_time_inc = increment_date(time, 0, output_freq, 0, 0, 0, 0)
   endif
else if(output_units == DIAG_YEARS) then
   if (get_calendar_type() == NO_CALENDAR) then
      call error_mesg('diag_time_inc', &
         'output units of years not allowed with no calendar', FATAL)
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

 function lcase (cs) 
!
    integer, parameter :: co=iachar('a')-iachar('A') ! case offset
    
    character(len=*), intent(in) :: cs ! character string 
    character(len=len(cs)) :: lcase 
    character :: ca(len(cs)) ! character array
    
    ca=transfer(cs,"x",len(cs)) 
    where (ca >= "A" .and. ca <= "Z") ca=achar(iachar(ca)+co) 
    lcase=transfer(ca,cs) 
    return
    
 end function lcase 

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
