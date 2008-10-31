#include <fms_platform.h>

module diag_data_mod


use time_manager_mod, only: time_type
use  mpp_domains_mod, only: domain1d, domain2d
use  mpp_io_mod,      only: fieldtype
public


#ifdef use_netCDF
#include <netcdf.inc>
#endif

! Specify storage limits for fixed size tables used for pointers, etc.
integer, parameter  :: max_fields_per_file = 300
integer, parameter  :: max_out_per_in_field = 10 !max number of output_fields per input_field
integer, parameter  :: max_files = 31
integer, parameter  :: DIAG_OTHER = 0
integer, parameter  :: DIAG_OCEAN = 1
integer, parameter  :: DIAG_ALL   = 2
integer, parameter  :: VERY_LARGE_FILE_FREQ = 10000
integer, parameter  :: VERY_LARGE_AXIS_LENGTH = 10000
integer, parameter  :: EVERY_TIME =  0
integer, parameter  :: END_OF_RUN = -1
integer, parameter  :: DIAG_SECONDS = 1, DIAG_MINUTES = 2, DIAG_HOURS = 3
integer, parameter  :: DIAG_DAYS = 4, DIAG_MONTHS = 5, DIAG_YEARS = 6
integer             :: num_files = 0
integer             :: max_input_fields = 300
integer             :: num_input_fields = 0
integer             :: max_output_fields = 300
integer             :: num_output_fields = 0
real                :: EMPTY = 0.0
real                :: FILL_VALUE = nf_fill_real  ! from file /usr/local/include/netcdf.inc
integer             :: null_axis_id
real                :: MAX_VALUE, MIN_VALUE
integer, parameter  :: single = selected_real_kind(p=6,r=37)
! Global data for all files
type (time_type)    :: base_time
integer             :: base_year, base_month, base_day, base_hour, base_minute, base_second
character(len = 256):: global_descriptor
integer             :: max_num_axis_sets = 25

type diag_grid
   real    :: start(3), end(3) ! coordinates (lat,lon,depth) of local domain to output   
   integer :: l_start_indx(3), l_end_indx(3) ! start and end indices at each LOCAL PE
   integer :: subaxes(3) ! id returned from diag_subaxes_init of 3 subaxes
end type diag_grid

type diag_fieldtype
   type(fieldtype)      :: Field
   type(domain2d)       :: Domain
   real                 :: miss, miss_pack
   logical              :: miss_present, miss_pack_present
   integer              :: tile_count
end type

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
   integer             :: duration
   integer             :: duration_units
   type (time_type)    :: next_open ! time to open a new file
   type (time_type)    :: start_time ! time for opening the file for the first time  
   type (time_type)    :: close_time  !file does not allow data after close time
end type file_type

type input_field_type
   character(len=128) :: module_name, field_name, long_name, units, standard_name
   character(len=64)  :: interp_method
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
   real, _ALLOCATABLE  :: buffer(:, :, :) _NULL
   real, _ALLOCATABLE  :: counter(:, :, :) _NULL
   type(time_type)     :: last_output, next_output, next_next_output
   real                :: count_0d
   type(diag_fieldtype):: f_type
   integer             :: axes(3), num_axes, num_elements, total_elements, region_elements
   type(diag_grid)     :: output_grid
   logical             :: local_output, need_compute, phys_window, written_once
   integer             :: imin, imax, jmin, jmax, kmin, kmax
   type(time_type)     :: Time_of_prev_field_data
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
type (input_field_type), allocatable :: input_fields(:)
type (output_field_type), allocatable, save :: output_fields(:)
type (time_type)    :: Time_zero
logical             :: append_pelist_name = .false.
logical             :: mix_snapshot_average_fields =.false.
logical             :: first_send_data_call = .true.
logical             :: module_is_initialized = .false.
logical             :: do_diag_field_log = .false.
logical             :: write_bytes_in_file = .false.
logical             :: debug_diag_manager = .false.
integer             :: diag_log_unit

character (len=10)  :: time_unit_list(6) = (/'seconds   ', 'minutes   ', &
   'hours     ', 'days      ', 'months    ', 'years     '/)
character(len=32)   :: pelist_name


character(len=128),private  :: version = '$Id: diag_data.F90,v 16.0.4.1 2008/09/08 17:42:27 z1l Exp $'
character(len=128),private  :: tagname = '$Name: perth_2008_10 $'


! definitions for diag_axis_mod
!maximum number of independent axes
integer, parameter  :: max_subaxes = 10
integer             :: max_axes = 60
type diag_axis_type
   character(len=128)   :: name
   character(len=256)   :: units, long_name
   character(len=1)     :: cart_name
   real, pointer        :: data(:)
   integer              :: start(max_subaxes)
   integer              :: end(max_subaxes)
   character(len=128)   :: subaxis_name(max_subaxes)
   integer              :: length, direction, edges, set, shift
   type(domain1d)       :: Domain
   type(domain2d)       :: Domain2
   character(len=128)   :: aux
   integer              :: tile_count
end type diag_axis_type

type diag_global_att_type
   character(len=128)   :: grid_type='regular'
   character(len=128)   :: tile_name='N/A'
end type diag_global_att_type


end module diag_data_mod
