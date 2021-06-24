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
!> @defgroup diag_data_mod diag_data_mod
!> @ingroup diag_manager
!! @brief Type descriptions and global variables for the diag_manager modules.
!! @author Seth Underwood <seth.underwood@noaa.gov>
!!
!!   Notation:
!!   <DL>
!!     <DT>input field</DT>
!!     <DD>The data structure describing the field as
!!       registered by the model code.</DD>
!!
!!     <DT>output field</DT>
!!     <DD>The data structure describing the actual
!!       diagnostic output with requested frequency and
!!       other options.</DD>
!!   </DL>
!!
!!   Input fields, output fields, and output files are gathered in arrays called
!!   "input_fields", "output_fields", and "files", respectively. Indices in these
!!   arrays are used as pointers to create associations between various data
!!   structures.
!!
!!   Each input field associated with one or several output fields via array of
!!   indices output_fields; each output field points to the single "parent" input
!!   field with the input_field index, and to the output file with the output_file
!!   index.

!> @file
!> @brief File for @ref diag_data_mod

!> @addtogroup diag_data_mod
!> @{
MODULE diag_data_mod
use platform_mod

  USE time_manager_mod, ONLY: time_type
  USE mpp_domains_mod, ONLY: domain1d, domain2d, domainUG
  USE mpp_io_mod, ONLY: fieldtype
  USE fms_mod, ONLY: WARNING, write_version_number
#ifdef use_netCDF
  ! NF90_FILL_REAL has value of 9.9692099683868690e+36.
  USE netcdf, ONLY: NF_FILL_REAL => NF90_FILL_REAL
#endif
  use fms2_io_mod

  IMPLICIT NONE

  PUBLIC

  ! Specify storage limits for fixed size tables used for pointers, etc.
  INTEGER, PARAMETER :: MAX_FIELDS_PER_FILE = 300 !< Maximum number of fields per file.
  INTEGER, PARAMETER :: DIAG_OTHER = 0
  INTEGER, PARAMETER :: DIAG_OCEAN = 1
  INTEGER, PARAMETER :: DIAG_ALL   = 2
  INTEGER, PARAMETER :: VERY_LARGE_FILE_FREQ = 100000
  INTEGER, PARAMETER :: VERY_LARGE_AXIS_LENGTH = 10000
  INTEGER, PARAMETER :: EVERY_TIME =  0
  INTEGER, PARAMETER :: END_OF_RUN = -1
  INTEGER, PARAMETER :: DIAG_SECONDS = 1, DIAG_MINUTES = 2, DIAG_HOURS = 3
  INTEGER, PARAMETER :: DIAG_DAYS = 4, DIAG_MONTHS = 5, DIAG_YEARS = 6
  INTEGER, PARAMETER :: MAX_SUBAXES = 10
  INTEGER, PARAMETER :: GLO_REG_VAL = -999 !< Value used in the region specification of the diag_table
                                           !! to indicate to use the full axis instead of a sub-axis
  INTEGER, PARAMETER :: GLO_REG_VAL_ALT = -1 !< Alternate value used in the region specification of the
                                             !! diag_table to indicate to use the full axis instead of a sub-axis
  REAL, PARAMETER :: CMOR_MISSING_VALUE = 1.0e20 !< CMOR standard missing value
  INTEGER, PARAMETER :: DIAG_FIELD_NOT_FOUND = -1 !< Return value for a diag_field that isn't found in the diag_table

  !> @}

  !> @brief Contains the coordinates of the local domain to output.
  !> @ingroup diag_data_mod
  TYPE diag_grid
     REAL, DIMENSION(3) :: start !< start coordinates (lat,lon,depth) of local domain to output
     REAL, DIMENSION(3) :: END !< end coordinates (lat,lon,depth) of local domain to output
     INTEGER, DIMENSION(3) :: l_start_indx !< start indices at each LOCAL PE
     INTEGER, DIMENSION(3) :: l_end_indx !< end indices at each LOCAL PE
     INTEGER, DIMENSION(3) :: subaxes !< id returned from diag_subaxes_init of 3 subaxes
  END TYPE diag_grid

  !> @brief Diagnostic field type
  !> @ingroup diag_data_mod
  TYPE diag_fieldtype
     TYPE(fieldtype) :: Field
     TYPE(domain2d) :: Domain
     TYPE(domainUG) :: DomainU
     REAL :: miss, miss_pack
     LOGICAL :: miss_present, miss_pack_present
     INTEGER :: tile_count
  END TYPE diag_fieldtype

  !> @brief Attribute type for diagnostic fields
  !> @ingroup diag_data_mod
  TYPE :: diag_atttype
     INTEGER             :: type !< Data type of attribute values (NF_INT, NF_FLOAT, NF_CHAR)
     INTEGER             :: len !< Number of values in attribute, or if a character string then
                                !! length of the string.
     CHARACTER(len=128)  :: name !< Name of the attribute
     CHARACTER(len=1280) :: catt !< Character string to hold character value of attribute
     REAL, allocatable, DIMENSION(:)    :: fatt !< REAL array to hold value of REAL attributes
     INTEGER, allocatable, DIMENSION(:) :: iatt !< INTEGER array to hold value of INTEGER attributes
  END TYPE diag_atttype

  !> @brief Define the region for field output
  !> @ingroup diag_data_mod
  TYPE coord_type
     REAL :: xbegin
     REAL :: xend
     REAL :: ybegin
     REAL :: yend
     REAL :: zbegin
     REAL :: zend
  END TYPE coord_type

  !> @brief Type to define the diagnostic files that will be written as defined by the diagnostic table.
  !> @ingroup diag_data_mod
  TYPE file_type
     CHARACTER(len=128) :: name !< Name of the output file.
     CHARACTER(len=128) :: long_name
     INTEGER, DIMENSION(max_fields_per_file) :: fields
     INTEGER :: num_fields
     INTEGER :: output_freq
     INTEGER :: output_units
     INTEGER :: FORMAT
     INTEGER :: time_units
     INTEGER :: file_unit
     INTEGER :: bytes_written
     INTEGER :: time_axis_id, time_bounds_id
     INTEGER :: new_file_freq !< frequency to create new file
     INTEGER :: new_file_freq_units !< time units of new_file_freq (days, hours, years, ...)
     INTEGER :: duration
     INTEGER :: duration_units
     INTEGER :: tile_count
     LOGICAL :: local !< .TRUE. if fields are output in a region instead of global.
     TYPE(time_type) :: last_flush
     TYPE(time_type) :: next_open !< Time to open a new file.
     TYPE(time_type) :: start_time !< Time file opened.
     TYPE(time_type) :: close_time !< Time file closed.  File does not allow data after close time
     TYPE(diag_fieldtype):: f_avg_start, f_avg_end, f_avg_nitems, f_bounds
     TYPE(diag_atttype), allocatable, dimension(:) :: attributes !< Array to hold user definable attributes
     INTEGER :: num_attributes !< Number of defined attibutes
!----------
!ug support
     logical(I4_KIND) :: use_domainUG = .false.
     logical(I4_KIND) :: use_domain2D = .false.
!----------
!Check if time axis was already registered
     logical, allocatable :: is_time_axis_registered
!Support for fms2_io time
     real :: rtime_current
     integer :: time_index
     CHARACTER(len=10) :: filename_time_bounds
  END TYPE file_type

  !> @brief Type to hold the input field description
  !> @ingroup diag_data_mod
  TYPE input_field_type
     CHARACTER(len=128) :: module_name, field_name, long_name, units
     CHARACTER(len=256) :: standard_name
     CHARACTER(len=64) :: interp_method
     INTEGER, DIMENSION(3) :: axes
     INTEGER :: num_axes
     LOGICAL :: missing_value_present, range_present
     REAL :: missing_value
     REAL, DIMENSION(2) :: range
     INTEGER, allocatable, dimension(:) :: output_fields
     INTEGER :: num_output_fields
     INTEGER, DIMENSION(3) :: size
     LOGICAL :: static, register, mask_variant, local
     INTEGER :: numthreads
     INTEGER :: active_omp_level !< The current level of OpenMP nesting
     INTEGER :: tile_count
     TYPE(coord_type) :: local_coord
     TYPE(time_type)  :: time
     LOGICAL :: issued_mask_ignore_warning !< Indicates if the mask_ignore_warning
                                           !! has been issued for this input
                                           !! field.  Once .TRUE. the warning message
                                           !! is suppressed on all subsequent
                                           !! send_data calls.
  END TYPE input_field_type

  !> @brief Type to hold the output field description.
  !> @ingroup diag_data_mod
  TYPE output_field_type
     INTEGER :: input_field !< index of the corresponding input field in the table
     INTEGER :: output_file !< index of the output file in the table
     CHARACTER(len=128) :: output_name
     LOGICAL :: time_average !< true if the output field is averaged over time interval
     LOGICAL :: time_rms !< true if the output field is the rms.  If true, then time_average is also
     LOGICAL :: static
     LOGICAL :: time_max !< true if the output field is maximum over time interval
     LOGICAL :: time_min !< true if the output field is minimum over time interval
     LOGICAL :: time_sum !< true if the output field is summed over time interval
     LOGICAL :: time_ops !< true if any of time_min, time_max, time_rms or time_average is true
     INTEGER  :: pack
     INTEGER :: pow_value !< Power value to use for mean_pow(n) calculations
     CHARACTER(len=50) :: time_method !< time method field from the input file
     ! coordinates of the buffer and counter are (x, y, z, time-of-day)
     REAL, allocatable, DIMENSION(:,:,:,:) :: buffer !< coordinates of the buffer and counter are (x, y, z, time-of-day)
     REAL, allocatable, DIMENSION(:,:,:,:) :: counter !< coordinates of the buffer and counter are (x, y, z, time-of-day)
     ! the following two counters are used in time-averaging for some
     ! combination of the field options. Their size is the length of the
     ! diurnal axis; the counters must be tracked separately for each of
     ! the diurnal interval, because the number of time slices accumulated
     ! in each can be different, depending on time step and the number of
     ! diurnal samples.
     REAL, allocatable, DIMENSION(:)  :: count_0d !< the following two counters are used in time-averaging for some
     !! combination of the field options. Their size is the length of the
     !! diurnal axis; the counters must be tracked separately for each of
     !! the diurnal interval, because the number of time slices accumulated
     !! in each can be different, depending on time step and the number of
     !! diurnal samples.
     INTEGER, allocatable, dimension(:) :: num_elements !< the following two counters are used in time-averaging for some
     !! combination of the field options. Their size is the length of the
     !! diurnal axis; the counters must be tracked separately for each of
     !! the diurnal interval, because the number of time slices accumulated
     !! in each can be different, depending on time step and the number of
     !! diurnal samples.

     TYPE(time_type) :: last_output, next_output, next_next_output
     TYPE(diag_fieldtype) :: f_type
     INTEGER, DIMENSION(4) :: axes
     INTEGER :: num_axes, total_elements, region_elements
     INTEGER :: n_diurnal_samples !< number of diurnal sample intervals, 1 or more
     TYPE(diag_grid) :: output_grid
     LOGICAL :: local_output, need_compute, phys_window, written_once
     LOGICAL :: reduced_k_range
     INTEGER :: imin, imax, jmin, jmax, kmin, kmax
     TYPE(time_type) :: Time_of_prev_field_data
     TYPE(diag_atttype), allocatable, dimension(:) :: attributes
     INTEGER :: num_attributes
!----------
!ug support
     logical :: reduced_k_unstruct = .false.
!----------
  END TYPE output_field_type

  !> @brief Type to hold the diagnostic axis description.
  !> @ingroup diag_data_mod
  TYPE diag_axis_type
     CHARACTER(len=128) :: name
     CHARACTER(len=256) :: units, long_name
     CHARACTER(len=1) :: cart_name
     REAL, DIMENSION(:), POINTER :: data
     INTEGER, DIMENSION(MAX_SUBAXES) :: start
     INTEGER, DIMENSION(MAX_SUBAXES) :: end
     CHARACTER(len=128), DIMENSION(MAX_SUBAXES) :: subaxis_name
     INTEGER :: length, direction, edges, set, shift
     TYPE(domain1d) :: Domain
     TYPE(domain2d) :: Domain2
     TYPE(domain2d), dimension(MAX_SUBAXES) :: subaxis_domain2
     type(domainUG) :: DomainUG
     CHARACTER(len=128) :: aux, req
     INTEGER :: tile_count
     TYPE(diag_atttype), allocatable, dimension(:) :: attributes !< Array to hold user definable attributes
     INTEGER :: num_attributes !< Number of defined attibutes
     INTEGER :: domain_position !< The position in the doman (NORTH or EAST or CENTER)
  END TYPE diag_axis_type

  !> @ingroup diag_data_mod
  TYPE diag_global_att_type
     CHARACTER(len=128)   :: grid_type='regular'
     CHARACTER(len=128)   :: tile_name='N/A'
  END TYPE diag_global_att_type

! Include variable "version" to be written to log file.
#include<file_version.h>

  !> @addtogroup diag_data_mod
  !> @{

  ! <!-- Other public variables -->
  INTEGER :: num_files = 0 !< Number of output files currenly in use by the diag_manager.
  INTEGER :: num_input_fields = 0 !< Number of input fields in use.
  INTEGER :: num_output_fields = 0 !< Number of output fields in use.
  INTEGER :: null_axis_id

  ! <!-- Namelist variables -->
  LOGICAL :: append_pelist_name = .FALSE.
  LOGICAL :: mix_snapshot_average_fields =.FALSE.
  INTEGER :: max_files = 31 !< Maximum number of output files allowed.  Increase via diag_manager_nml.
  INTEGER :: max_output_fields = 300 !< Maximum number of output fields.  Increase via diag_manager_nml.
  INTEGER :: max_input_fields = 600 !< Maximum number of input fields.  Increase via diag_manager_nml.
  INTEGER :: max_out_per_in_field = 150 !< Maximum number of output_fields per input_field.  Increase via diag_manager_nml.
  INTEGER :: max_axes = 60 !< Maximum number of independent axes.
  LOGICAL :: do_diag_field_log = .FALSE.
  LOGICAL :: write_bytes_in_file = .FALSE.
  LOGICAL :: debug_diag_manager = .FALSE.
  LOGICAL :: flush_nc_files = .FALSE. !< Control if diag_manager will force a
                                      !! flush of the netCDF file on each write.
                                      !! Note: changing this to .TRUE. can greatly
                                      !! reduce the performance of the model, as the
                                      !! model must wait until the flush to disk has
                                      !! completed.
  INTEGER :: max_num_axis_sets = 25
  LOGICAL :: use_cmor = .FALSE. !< Indicates if we should overwrite the MISSING_VALUE to use the CMOR missing value.
  LOGICAL :: issue_oor_warnings = .TRUE. !< Issue warnings if the output field has values outside the given
                                         !! range for a variable.
  LOGICAL :: oor_warnings_fatal = .FALSE. !< Cause a fatal error if the output field has a value outside the
                                          !! given range for a variable.
  LOGICAL :: region_out_use_alt_value = .TRUE. !< Will determine which value to use when checking a regional
                                               !! output if the region is the full axis or a sub-axis.
                                               !! The values are defined as <TT>GLO_REG_VAL</TT>
                                               !! (-999) and <TT>GLO_REG_VAL_ALT</TT> (-1) in <TT>diag_data_mod</TT>.

  INTEGER :: max_field_attributes = 4 !< Maximum number of user definable attributes per field. Liptak: Changed from 2 to 4 20170718
  INTEGER :: max_file_attributes = 2 !< Maximum number of user definable global attributes per file.
  INTEGER :: max_axis_attributes = 4 !< Maximum number of user definable attributes per axis.
  LOGICAL :: prepend_date = .TRUE. !< Should the history file have the start date prepended to the file name.
                                   !! <TT>.TRUE.</TT> is only supported if the diag_manager_init
                                   !! routine is called with the optional time_init parameter.
  LOGICAL :: use_mpp_io = .false. !< false is fms2_io (default); true is mpp_io

  ! <!-- netCDF variable -->

#ifdef use_netCDF
  REAL :: FILL_VALUE = NF_FILL_REAL !< Fill value used.  Value will be <TT>NF90_FILL_REAL</TT> if using the
                                    !! netCDF module, otherwise will be 9.9692099683868690e+36.
                                    ! from file /usr/local/include/netcdf.inc
#else
  REAL :: FILL_VALUE = 9.9692099683868690e+36
#endif

  INTEGER :: pack_size = 1 !< 1 for double and 2 for float

  ! <!-- REAL public variables -->
  REAL :: EMPTY = 0.0
  REAL :: MAX_VALUE, MIN_VALUE

  ! <!-- Global data for all files -->
  TYPE(time_type) :: diag_init_time !< Time diag_manager_init called.  If init_time not included in
                                    !! diag_manager_init call, then same as base_time
  TYPE(time_type) :: base_time
  INTEGER :: base_year, base_month, base_day, base_hour, base_minute, base_second
  CHARACTER(len = 256):: global_descriptor

  ! <!-- ALLOCATABLE variables -->
  TYPE(file_type), SAVE, ALLOCATABLE :: files(:)
  TYPE(input_field_type), ALLOCATABLE :: input_fields(:)
  TYPE(output_field_type), ALLOCATABLE :: output_fields(:)
  ! used if use_mpp_io = .false.
  type(FmsNetcdfUnstructuredDomainFile_t),allocatable, target :: fileobjU(:)
  type(FmsNetcdfDomainFile_t),allocatable, target :: fileobj(:)
  type(FmsNetcdfFile_t),allocatable, target :: fileobjND(:)
  character(len=2),allocatable :: fnum_for_domain(:) !< If this file number in the array is for the "unstructured" or "2d" domain

  ! <!-- Even More Variables -->
  TYPE(time_type) :: time_zero
  LOGICAL :: first_send_data_call = .TRUE.
  LOGICAL :: module_is_initialized = .FALSE. !< Indicate if diag_manager has been initialized
  INTEGER :: diag_log_unit
  CHARACTER(len=10), DIMENSION(6) :: time_unit_list = (/'seconds   ', 'minutes   ',&
       & 'hours     ', 'days      ', 'months    ', 'years     '/)
  CHARACTER(len=32) :: pelist_name
  INTEGER :: oor_warning = WARNING

CONTAINS

  !> @brief Initialize and write the version number of this file to the log file.
  SUBROUTINE diag_data_init()
    IF (module_is_initialized) THEN
       RETURN
    END IF

    ! Write version number out to log file
    call write_version_number("DIAG_DATA_MOD", version)
  END SUBROUTINE diag_data_init

END MODULE diag_data_mod
!> @}
! close documentation grouping
