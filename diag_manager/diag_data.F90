#include <fms_platform.h>

MODULE diag_data_mod


  USE time_manager_mod, ONLY: time_type
  USE mpp_domains_mod,  ONLY: domain1d, domain2d
  USE mpp_io_mod,       ONLY: fieldtype

#ifdef use_netCDF
  USE netcdf
#endif

  PUBLIC

  ! Specify storage limits for fixed size tables used for pointers, etc.
  INTEGER, PARAMETER  :: MAX_FIELDS_PER_FILE = 300
  INTEGER, PARAMETER  :: MAX_OUT_PER_IN_FIELD = 10 !max number of output_fields per input_field
  INTEGER, PARAMETER  :: DIAG_OTHER = 0
  INTEGER, PARAMETER  :: DIAG_OCEAN = 1
  INTEGER, PARAMETER  :: DIAG_ALL   = 2
  INTEGER, PARAMETER  :: VERY_LARGE_FILE_FREQ = 100000
  INTEGER, PARAMETER  :: VERY_LARGE_AXIS_LENGTH = 10000
  INTEGER, PARAMETER  :: EVERY_TIME =  0
  INTEGER, PARAMETER  :: END_OF_RUN = -1
  INTEGER, PARAMETER  :: DIAG_SECONDS = 1, DIAG_MINUTES = 2, DIAG_HOURS = 3
  INTEGER, PARAMETER  :: DIAG_DAYS = 4, DIAG_MONTHS = 5, DIAG_YEARS = 6
  INTEGER :: MAX_FILES = 31
  INTEGER :: num_files = 0
  INTEGER :: max_input_fields = 300
  INTEGER :: num_input_fields = 0
  INTEGER :: max_output_fields = 300
  INTEGER :: num_output_fields = 0
  REAL :: EMPTY = 0.0
#ifdef use_netCDF
  REAL                :: FILL_VALUE = NF90_FILL_REAL  ! from file /usr/local/include/netcdf.inc
#else
  REAL                :: FILL_VALUE = 9.9692099683868690e+36 
#endif

  INTEGER             :: null_axis_id
  REAL                :: MAX_VALUE, MIN_VALUE
  INTEGER, PARAMETER  :: SINGLE = SELECTED_REAL_KIND(p=6,r=37)
  ! Global data for all files
  TYPE (time_type)    :: base_time
  INTEGER             :: base_year, base_month, base_day, base_hour, base_minute, base_second
  CHARACTER(len = 256):: global_descriptor
  INTEGER             :: max_num_axis_sets = 25
  
  TYPE diag_grid
     REAL    :: start(3), END(3) ! coordinates (lat,lon,depth) of local domain to output   
     INTEGER :: l_start_indx(3), l_end_indx(3) ! start and end indices at each LOCAL PE
     INTEGER :: subaxes(3) ! id returned from diag_subaxes_init of 3 subaxes
  END TYPE diag_grid
  
  TYPE diag_fieldtype
     TYPE(fieldtype)      :: Field
     TYPE(domain2d)       :: Domain
     REAL                 :: miss, miss_pack
     LOGICAL              :: miss_present, miss_pack_present
     INTEGER              :: tile_count
  END TYPE diag_fieldtype
  
  TYPE coord_type  ! define the region for outputting a field
     REAL :: xbegin
     REAL :: xend
     REAL :: ybegin
     REAL :: yend
     REAL :: zbegin
     REAL :: zend
  END TYPE coord_type
  
  
  TYPE file_type
     CHARACTER(len=128)  :: name
     INTEGER             :: output_freq
     INTEGER             :: output_units
     INTEGER             :: FORMAT
     INTEGER             :: time_units
     CHARACTER(len=128)  :: long_name
     INTEGER             :: fields(max_fields_per_file)
     INTEGER             :: num_fields
     INTEGER             :: file_unit
     INTEGER             :: bytes_written
     INTEGER             :: time_axis_id, time_bounds_id
     TYPE (time_type)    :: last_flush
     TYPE(diag_fieldtype):: f_avg_start, f_avg_end, f_avg_nitems, f_bounds
     LOGICAL             :: local ! true if fields are output in a region instead of global
     INTEGER             :: new_file_freq ! frequency to create new file
     INTEGER             :: new_file_freq_units ! time units of new_file_freq (days, hours, years, ...)
     INTEGER             :: duration
     INTEGER             :: duration_units
     INTEGER             :: tile_count
     TYPE (time_type)    :: next_open ! time to open a new file
     TYPE (time_type)    :: start_time ! time for opening the file for the first time  
     TYPE (time_type)    :: close_time  !file does not allow data after close time
  END TYPE file_type
  
  
  ! Notation: 
  !   "input field": the data structure describing the field as registered by the 
  !        model code
  !   "output field": the data structure describing the actual diagnostic output
  !        with requested frequency and other options 
  !
  ! Input fields, output fields, and output files are gathered in arrays called
  ! "input_fields", "output_fields", and "files", respectively. Indices in these
  ! arrays are used as pointers to create associations between various data
  ! structures.
  !
  ! Each input field associated with one or several output fields via array of
  ! indices output_fields; each output field points to the single "parent" input
  ! field with the input_field index, and to the output file with the output_file 
  ! index
  TYPE input_field_type
     CHARACTER(len=128) :: module_name, field_name, long_name, units, standard_name
     CHARACTER(len=64)  :: interp_method
     INTEGER            :: axes(3), num_axes 
     LOGICAL            :: missing_value_present, range_present
     REAL               :: missing_value, RANGE(2)
     INTEGER            :: output_fields(max_out_per_in_field)
     INTEGER            :: num_output_fields, SIZE(3) 
     LOGICAL            :: static, register, mask_variant, local
     INTEGER            :: tile_count
     TYPE(coord_type)   :: local_coord
  END TYPE input_field_type
  
  TYPE output_field_type
     INTEGER             :: input_field ! index of the corresponding input field in the table
     INTEGER             :: output_file ! index of the output file in the table
     CHARACTER(len=128)  :: output_name
     LOGICAL :: time_average ! true if the output field is averaged over time interval
     LOGICAL :: static
     LOGICAL :: time_max ! true if the output field is maximum over time interval
     LOGICAL :: time_min ! true if the output field is minimum over time interval
     LOGICAL :: time_ops ! true if any of time_min, time_max, or time_average is true
     INTEGER             :: pack
     CHARACTER(len=50) :: time_method ! time method field from the input file 
     ! coordianes of the buffer and counter are (x, y, z, time-of-day)
     REAL, _ALLOCATABLE  :: buffer (:, :, :, :) _NULL
     REAL, _ALLOCATABLE  :: counter(:, :, :, :) _NULL
     ! the following two counters are used in time-averaging for some 
     ! combination of the field options. Their size is the length of the 
     ! diurnal axis; the counters must be tracked separately for each of
     ! the diurnal interval, becaus the number of time slices accumulated
     ! in each can be different, depending on time step and the number of
     ! diurnal samples.
     REAL, _ALLOCATABLE  :: count_0d(:)
     INTEGER, _ALLOCATABLE :: num_elements(:)
     
     TYPE(time_type)     :: last_output, next_output, next_next_output
     TYPE(diag_fieldtype):: f_type
     INTEGER             :: axes(4), num_axes, total_elements, region_elements
     INTEGER :: n_diurnal_samples ! number of diurnal sample intervals, 1 or more
     TYPE(diag_grid)     :: output_grid
     LOGICAL             :: local_output, need_compute, phys_window, written_once
     LOGICAL             :: reduced_k_range
     INTEGER             :: imin, imax, jmin, jmax, kmin, kmax
     TYPE(time_type)     :: Time_of_prev_field_data
  END TYPE output_field_type
  
  ! local_output: true if this field is written out on a region, not global
  ! need_compute: true if this PE is involved in writing the field
  
  TYPE(file_type), SAVE, ALLOCATABLE :: files(:)
  TYPE(input_field_type), ALLOCATABLE :: input_fields(:)
  TYPE(output_field_type), ALLOCATABLE :: output_fields(:)
  TYPE(time_type)    :: Time_zero
  LOGICAL             :: append_pelist_name = .FALSE.
  LOGICAL             :: mix_snapshot_average_fields =.FALSE.
  LOGICAL             :: first_send_data_call = .TRUE.
  LOGICAL             :: module_is_initialized = .FALSE.
  LOGICAL             :: do_diag_field_log = .FALSE.
  LOGICAL             :: write_bytes_in_file = .FALSE.
  LOGICAL             :: debug_diag_manager = .FALSE.
  INTEGER             :: diag_log_unit
 
  CHARACTER (len=10)  :: time_unit_list(6) = (/'seconds   ', 'minutes   ',&
       & 'hours     ', 'days      ', 'months    ', 'years     '/)
  CHARACTER(len=32)   :: pelist_name
  
  
  CHARACTER(len=128),PRIVATE  :: version = '$Id: diag_data.F90,v 17.0 2009/07/21 03:18:43 fms Exp $'
  CHARACTER(len=128),PRIVATE  :: tagname = '$Name: quebec_200910 $'
  
  
  ! definitions for diag_axis_mod
  !maximum number of independent axes
  INTEGER, PARAMETER  :: max_subaxes = 10
  INTEGER             :: max_axes = 60
  TYPE diag_axis_type
     CHARACTER(len=128)   :: name
     CHARACTER(len=256)   :: units, long_name
     CHARACTER(len=1)     :: cart_name
     REAL, POINTER        :: data(:)
     INTEGER              :: start(max_subaxes)
     INTEGER              :: END(max_subaxes)
     CHARACTER(len=128)   :: subaxis_name(max_subaxes)
     INTEGER              :: length, direction, edges, set, shift
     TYPE(domain1d)       :: Domain
     TYPE(domain2d)       :: Domain2
     CHARACTER(len=128)   :: aux
     INTEGER              :: tile_count
  END TYPE diag_axis_type

  TYPE diag_global_att_type
     CHARACTER(len=128)   :: grid_type='regular'
     CHARACTER(len=128)   :: tile_name='N/A'
  END TYPE diag_global_att_type
  
  
END MODULE diag_data_mod
