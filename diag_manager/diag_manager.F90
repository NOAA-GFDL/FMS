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
!> @defgroup diag_manager_mod diag_manager_mod
!> @ingroup diag_manager
!! @brief diag_manager_mod is a set of simple calls for parallel diagnostics
!!   on distributed systems. It is geared toward the writing of data in netCDF
!!   format. See @ref diag_manager for diag table information.
!! @author Matt Harrison, Giang Nong, Seth Underwood
!!
!!   <TT>diag_manager_mod</TT> provides a convenient set of interfaces for
!!   writing data to disk.  It is built upon the parallel I/O interface of FMS
!!   code <TT>/shared/mpp/mpp_io.F90</TT>.
!!
!!   A single group of calls to the <TT>diag_manager_mod</TT> interfaces
!!   provides data to disk at any number of sampling and/or averaging intervals
!!   specified at run-time. Run-time specification of diagnostics are input
!!   through the diagnostics table.
!!
!!   <H4>Usage</H4>
!!   Use of <TT>diag_manager</TT> includes the following steps:
!!   <OL>
!!     <LI> Create diag_table as described in the @ref diag_table_mod
!!          documentation.</LI>
!!     <LI> Call @ref diag_manager_init to initialize
!!          diag_manager_mod.</LI>
!!     <LI> Call @ref register_diag_field to register the field to be
!!          output.
!!          <B>NOTE:</B> ALL fields in diag_table should be registered <I>BEFORE</I>
!!          the first send_data call</LI>
!!     <LI> Call @ref send_data to send data to output fields </LI>
!!     <LI> Call @ref diag_manager_end to exit diag_manager </LI>
!!   </OL>
!!
!!   <H4>Features</H4>
!!   Features of <TT>diag_manager_mod</TT>:
!!   <OL>
!!     <LI> Ability to output from 0D arrays (scalars) to 3D arrays.</LI>
!!     <LI> Ability to output time average of fields that have time dependent
!!          mask.</LI>
!!     <LI> Give optional warning if <TT>register_diag_field</TT> fails due to
!!          misspelled module name or field name.</LI>
!!     <LI> Check if a field is registered twice.</LI>
!!     <LI> Check for duplicate lines in diag_table. </LI>
!!     <LI> @ref diag_table_mod can contain fields
!!          that are NOT written to any files. The file name in diag_table of
!!          these fields is <TT>null</TT>.</LI>
!!     <LI> By default, a field is output in its global grid.  The user can now
!!          output a field in a specified region.  See
!!          @ref send_data for more details.</LI>
!!     <LI> To check if the diag table is set up correctly, user should set
!!          <TT>debug_diag_manager=.true.</TT> in diag_manager namelist, then
!!          the the content of diag_table is printed in stdout.</LI>
!!     <LI> New optional format of file information in @ref diag_table_mod.
!!          It is possible to have just
!!          one file name and reuse it many times. A time string will be appended to
!!          the base file name each time a new file is opened. The time string can be
!!          any combination from year to second of current model time.
!!
!!          Here is an example file line: <BR />
!!          <PRE>"file2_yr_dy%1yr%3dy",2,"hours",1,"hours","Time", 10, "days", "1 1 7 0 0 0", 6, "hours"</PRE>
!!          <BR />
!!
!!          From left to right we have:
!!          <UL>
!!            <LI>file name</LI>
!!            <LI>output frequency</LI>
!!            <LI>output frequency unit</LI>
!!            <LI>Format (should always be 1)</LI>
!!            <LI>time axis unit</LI>
!!            <LI>time axis name</LI>
!!            <LI>frequency for creating new file</LI>
!!            <LI>unit for creating new file</LI>
!!            <LI>start time of the new file</LI>
!!            <LI>file duration</LI>
!!            <LI>file duration unit.</LI>
!!          </UL>
!!          The 'file duration', if absent, will be equal to frequency for creating a new file.
!!
!!          Thus, the above means: create a new file every 10 days, each file will last 6 hours from creation time, no files will
!!          be created before time "1 1 7 0 0 0".
!!
!!          In this example the string
!!          <TT>10, "days", "1 1 7 0 0 0", 6, "hours"</TT> is optional.
!!
!!          Keywords for the time string suffix is
!!          <TT>%xyr,%xmo,%xdy,%xhr,%xmi,%xsc</TT> where <TT>x</TT> is a
!!          mandatory 1 digit number specifying the width of field used in
!!          writing the string</LI>
!!     <LI> New time axis for time averaged fields.  Users can use a namelist option to handle the time value written
!!          to time axis for time averaged fields.
!!
!!          If <TT>mix_snapshot_average_fields=.true.</TT> then a time averaged file will have time values corresponding to
!!          ending time_bound e.g. January monthly average is labeled Feb01. Users can have both snapshot and averaged fields in
!!          one file.
!!
!!          If <TT>mix_snapshot_average_fields=.false.</TT> The time value written to time axis for time averaged fields is the
!!          middle on the averaging time. For example, January monthly mean will be written at Jan 16 not Feb 01 as
!!          before. However, to use this new feature users should <B>separate</B> snapshot fields and time averaged fields in
!!          <B>different</B> files or a fatal error will occur.
!!
!!          The namelist <B>default</B> value is <TT>mix_snapshot_average_fields=.false.</TT></LI>
!!     <LI> Time average, Root Mean Square, Max and Min, and diurnal. In addition to time average users can also get then Root Mean Square, Max or Min value
!!          during the same interval of time as time average. For this purpose, in the diag table users must replace
!!          <TT>.true.</TT> or <TT>.false.</TT> by <TT>rms</TT>, <TT>max</TT> or <TT>min</TT>.  <B><I>Note:</I></B> Currently, max
!!          and min are not available for regional output.
!!
!!          A diurnal average or the average of an integer power can also be requested using <TT>diurnal##</TT> or <TT>pow##</TT> where
!!          <TT>##</TT> are the number of diurnal sections or integer power to average.</LI>
!!     <LI> <TT>standard_name</TT> is added as optional argument in @ref register_diag_field. </LI>
!!     <LI>When namelist variable <TT>debug_diag_manager = .true.</TT> array
!!         bounds are checked in @ref send_data.</LI>
!!     <LI>Coordinate attributes can be written in the output file if the
!!         argument "<TT>aux</TT>" is given in @ref diag_axis_mod#diag_axis_init . The
!!         corresponding fields (geolat/geolon) should also be written to the
!!         same file.</LI>
!!   </OL>

!> @file
!> @brief File for @ref diag_manager_mod

MODULE diag_manager_mod
use platform_mod

  ! <NAMELIST NAME="diag_manager_nml">
  !   <DATA NAME="append_pelist_name" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="mix_snapshot_average_fields" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !     Set to .TRUE. to allow both time average and instantaneous fields in the same output file.
  !   </DATA>
  !   <DATA NAME="max_files" TYPE="INTEGER" DEFULT="31">
  !   </DATA>
  !   <DATA NAME="max_output_fields" TYPE="INTEGER" DEFAULT="300">
  !   </DATA>
  !   <DATA NAME="max_input_fields" TYPE="INTEGER" DEFAULT="300">
  !   </DATA>
  !   <DATA NAME="max_axes" TYPE="INTEGER" DEFAULT="60">
  !   </DATA>
  !   <DATA NAME="do_diag_field_log" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="write_bytes_in_files" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="debug_diag_manager" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !   </DATA>
  !   <DATA NAME="max_num_axis_sets" TYPE="INTEGER" DEFAULT="25">
  !   </DATA>
  !   <DATA NAME="use_cmor" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !     Let the <TT>diag_manager</TT> know if the missing value (if supplied) should be overridden to be the
  !     CMOR standard value of -1.0e20.
  !   </DATA>
  !   <DATA NAME="issue_oor_warnings" TYPE="LOGICAL" DEFAULT=".TRUE.">
  !     If <TT>.TRUE.</TT>, then the <TT>diag_manager</TT> will check for values outside the valid range.  This range is defined in
  !     the model, and passed to the <TT>diag_manager_mod</TT> via the OPTIONAL variable range in the <TT>register_diag_field</TT>
  !     function.
  !   </DATA>
  !   <DATA NAME="oor_warnings_fatal" TYPE="LOGICAL" DEFAULT=".FALSE.">
  !     If <TT>.TRUE.</TT> then <TT>diag_manager_mod</TT> will issue a <TT>FATAL</TT> error if any values for the output field are
  !     outside the given range.
  !   </DATA>
  !   <DATA NAME="max_field_attributes" TYPE="INTEGER" DEFAULT="4">
  !     Maximum number of user definable attributes per field.
  !   </DATA>
  !   <DATA NAME="max_file_attributes" TYPE="INTEGER" DEFAULT="2">
  !     Maximum number of user definable global attributes per file.
  !   </DATA>
  !   <DATA NAME="prepend_date" TYPE="LOGICAL" DEFAULT=".TRUE.">
  !     If <TT>.TRUE.</TT> then prepend the file start date to the output file.  <TT>.TRUE.</TT> is only supported if the
  !      diag_manager_init routine is called with the optional time_init parameter.  Note: This was usually done by FRE after the
  !     model run.
  !   </DATA>
  !   <DATA NAME="region_out_use_alt_value" TYPE="LOGICAL" DEFAULT=".TRUE.">
  !     Will determine which value to use when checking a regional output if the region is the full axis or a sub-axis.
  !     The values are defined as <TT>GLO_REG_VAL</TT> (-999) and <TT>GLO_REG_VAL_ALT</TT> (-1) in <TT>diag_data_mod</TT>.
  !   </DATA>
  !   <DATA NAME="use_mpp_io" TYPE="LOGICAL" DEFAULT=".false.">
  !    Set to true, diag_manager uses mpp_io.  Default is fms2_io.
  !   </DATA>
  ! </NAMELIST>

  USE time_manager_mod, ONLY: set_time, set_date, get_time, time_type, OPERATOR(>=), OPERATOR(>),&
       & OPERATOR(<), OPERATOR(==), OPERATOR(/=), OPERATOR(/), OPERATOR(+), ASSIGNMENT(=), get_date, &
       & get_ticks_per_second
  USE mpp_io_mod, ONLY: mpp_open, mpp_close, mpp_get_maxunits
  USE mpp_mod, ONLY: mpp_get_current_pelist, mpp_pe, mpp_npes, mpp_root_pe, mpp_sum

  USE mpp_mod, ONLY: input_nml_file

  USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdout, stdlog, write_version_number,&
       & file_exist, fms_error_handler, check_nml_error, get_mosaic_tile_file, lowercase
  USE fms_io_mod, ONLY: get_instance_filename
  USE diag_axis_mod, ONLY: diag_axis_init, get_axis_length, get_axis_num, get_domain2d, get_tile_count,&
       & diag_axis_add_attribute, axis_compatible_check, CENTER, NORTH, EAST
  USE diag_util_mod, ONLY: get_subfield_size, log_diag_field_info, update_bounds,&
       & check_out_of_bounds, check_bounds_are_exact_dynamic, check_bounds_are_exact_static,&
       & diag_time_inc, find_input_field, init_input_field, init_output_field,&
       & diag_data_out, write_static, get_date_dif, get_subfield_vert_size, sync_file_times,&
       & prepend_attribute, attribute_init, diag_util_init
  USE diag_data_mod, ONLY: max_files, CMOR_MISSING_VALUE, DIAG_OTHER, DIAG_OCEAN, DIAG_ALL, EVERY_TIME,&
       & END_OF_RUN, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS, num_files,&
       & max_input_fields, max_output_fields, num_output_fields, EMPTY, FILL_VALUE, null_axis_id,&
       & MAX_VALUE, MIN_VALUE, base_time, base_year, base_month, base_day,&
       & base_hour, base_minute, base_second, global_descriptor, coord_type, files, input_fields,&
       & output_fields, Time_zero, append_pelist_name, mix_snapshot_average_fields,&
       & first_send_data_call, do_diag_field_log, write_bytes_in_file, debug_diag_manager,&
       & diag_log_unit, time_unit_list, pelist_name, max_axes, module_is_initialized, max_num_axis_sets,&
       & use_cmor, issue_oor_warnings, oor_warnings_fatal, oor_warning, pack_size,&
       & max_out_per_in_field, flush_nc_files, region_out_use_alt_value, max_field_attributes, output_field_type,&
       & max_file_attributes, max_axis_attributes, prepend_date, DIAG_FIELD_NOT_FOUND, diag_init_time, diag_data_init, &
       & use_mpp_io
  USE diag_data_mod, ONLY:  fileobj, fileobjU, fnum_for_domain, fileobjND
  USE diag_table_mod, ONLY: parse_diag_table
  USE diag_output_mod, ONLY: get_diag_global_att, set_diag_global_att
  USE diag_grid_mod, ONLY: diag_grid_init, diag_grid_end
  USE constants_mod, ONLY: SECONDS_PER_DAY

#ifdef use_netCDF
  USE netcdf, ONLY: NF90_INT, NF90_FLOAT, NF90_CHAR
#endif

!----------
!ug support
  use diag_axis_mod, only: DIAG_AXIS_2DDOMAIN
  use diag_axis_mod, only: DIAG_AXIS_UGDOMAIN
!----------

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diag_manager_init, send_data, send_tile_averaged_data, diag_manager_end,&
       & register_diag_field, register_static_field, diag_axis_init, get_base_time, get_base_date,&
       & need_data, DIAG_ALL, DIAG_OCEAN, DIAG_OTHER, get_date_dif, DIAG_SECONDS,&
       & DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS, get_diag_global_att,&
       & set_diag_global_att, diag_field_add_attribute, diag_field_add_cell_measures,&
       & get_diag_field_id, diag_axis_add_attribute, CMOR_MISSING_VALUE, null_axis_id
  PUBLIC :: CENTER, NORTH, EAST !< Used for diag_axis_init
  ! Public interfaces from diag_grid_mod
  PUBLIC :: diag_grid_init, diag_grid_end
  PUBLIC :: diag_manager_set_time_end, diag_send_complete
  PUBLIC :: diag_send_complete_instant
  ! Public interfaces from diag_data_mod
  PUBLIC :: DIAG_FIELD_NOT_FOUND

  ! version number of this module
  ! Include variable "version" to be written to log file.
#include<file_version.h>

  type(time_type) :: Time_end

  !> @brief Send data over to output fields.
  !!
  !> <TT>send_data</TT> is overloaded for fields having zero dimension
  !! (scalars) to 3 dimension.  <TT>diag_field_id</TT> corresponds to the id
  !! returned from a previous call to <TT>register_diag_field</TT>. The field
  !! array is restricted to the computational range of the array. Optional
  !! argument <TT>is_in</TT> can be used to update sub-arrays of the entire
  !! field. Additionally, an optional logical or real mask can be used to
  !! apply missing values to the array.
  !!
  !! If a field is declared to be <TT>mask_variant</TT> in
  !! <TT>register_diag_field</TT> logical mask should be mandatory.
  !!
  !! For the real  mask, the mask is applied if the mask value is less than
  !! 0.5.
  !!
  !! By default, a field will be written out entirely in its global grid.
  !! Users can also specify regions in which the field will be output. The
  !! region is specified in diag-table just before the end of output_field
  !! replacing "none".
  !!
  !! For example, by default:
  !!
  !! "ocean_mod","Vorticity","vorticity","file1","all",.false.,"none",2
  !!
  !! for regional output:
  !!
  !! "ocean_mod","Vorticity","vorticity_local","file2","all",.false.,"0.5 53.5 -89.5 -28.5 -1 -1",2
  !!
  !! The format of a region is "<TT>xbegin xend ybegin yend zbegin zend</TT>".
  !! If it is a 2D field use (-1 -1) for (zbegin zend) as in the example above.
  !! For a 3D field use (-1 -1) for (zbegin zend) when you want to write the
  !! entire vertical extent, otherwise specify real coordinates.  The units
  !! used for region are the actual units used in grid_spec.nc (for example
  !! degrees for lat, lon).  <B><I>NOTE:</I></B> A FATAL error will occur if
  !! the region's boundaries are not found in grid_spec.nc.
  !!
  !! Regional output on the cubed sphere grid is also supported.  To use regional
  !! output on the cubed sphere grid, first the grid information needs to be sent to
  !! <TT>diag_manager_mod</TT> using the @ref diag_grid#diag_grid_init subroutine.
  !!
  !! @note When using regional output the files containing regional
  !! outputs should be different from files containing global (default) output.
  !! It is a FATAL error to have one file containing both regional and global
  !! results. For maximum flexibility and independence from PE counts one file
  !! should contain just one region.
  !!
  !!
  !! Time averaging is supported in regional output.
  !!
  !! Physical fields (written in "physics windows" of atmospheric code) are
  !! fully supported for regional outputs.
  !!
  !! <B><I>NOTE:</I></B> Most fields are defined in the data domain but use the
  !! compute domain. In <TT>send_data</TT> the field can be passed in either
  !! the data domain or in the compute domain.  If the data domain is used, the
  !! start and end indicies of the compute domain (isc, iec, . . .) should be
  !! passed.  If the compute domain is used no indices are needed.  The indices
  !! are for determining halo exclusively.  If users want to output the field
  !! partially they should use regional output as mentioned above.
  !!
  !! Weight in Time averaging is now supported, each time level may have a
  !! different weight. The default of weight is 1.
  !> @ingroup diag_manager_mod
  INTERFACE send_data
     MODULE PROCEDURE send_data_0d
     MODULE PROCEDURE send_data_1d
     MODULE PROCEDURE send_data_2d
     MODULE PROCEDURE send_data_3d
#ifdef OVERLOAD_R8
     MODULE PROCEDURE send_data_2d_r8
     MODULE PROCEDURE send_data_3d_r8
#endif
  END INTERFACE

  !> @brief Register a diagnostic field for a given module
  !> @ingroup diag_manager_mod
  INTERFACE register_diag_field
     MODULE PROCEDURE register_diag_field_scalar
     MODULE PROCEDURE register_diag_field_array
  END INTERFACE

  !> @brief Send tile-averaged data over to output fields.
  !> @ingroup diag_manager_mod
  INTERFACE send_tile_averaged_data
     MODULE PROCEDURE send_tile_averaged_data1d
     MODULE PROCEDURE send_tile_averaged_data2d
     MODULE PROCEDURE send_tile_averaged_data3d
  END INTERFACE

  !> @brief Add a attribute to the output field
  !> @ingroup diag_manager_mod
  INTERFACE diag_field_add_attribute
     MODULE PROCEDURE diag_field_add_attribute_scalar_r
     MODULE PROCEDURE diag_field_add_attribute_scalar_i
     MODULE PROCEDURE diag_field_add_attribute_scalar_c
     MODULE PROCEDURE diag_field_add_attribute_r1d
     MODULE PROCEDURE diag_field_add_attribute_i1d
  END INTERFACE diag_field_add_attribute

!> @addtogroup diag_manager_mod
!> @{
CONTAINS

  !> @brief Registers a scalar field
  !! @return field index for subsequent call to send_data.
  INTEGER FUNCTION register_diag_field_scalar(module_name, field_name, init_time, &
       & long_name, units, missing_value, range, standard_name, do_not_log, err_msg,&
       & area, volume, realm)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    TYPE(time_type), OPTIONAL, INTENT(in) :: init_time
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name, units, standard_name
    REAL, OPTIONAL, INTENT(in) :: missing_value
    REAL,  DIMENSION(2), OPTIONAL, INTENT(in) :: RANGE
    LOGICAL, OPTIONAL, INTENT(in) :: do_not_log !< if TRUE, field information is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(out):: err_msg
    INTEGER, OPTIONAL, INTENT(in) :: area, volume
    CHARACTER(len=*), OPTIONAL, INTENT(in):: realm !< String to set as the value to the modeling_realm attribute

    IF ( PRESENT(err_msg) ) err_msg = ''

    IF ( PRESENT(init_time) ) THEN
       register_diag_field_scalar = register_diag_field_array(module_name, field_name,&
            & (/null_axis_id/), init_time,long_name, units, missing_value, range, &
            & standard_name=standard_name, do_not_log=do_not_log, err_msg=err_msg,&
            & area=area, volume=volume, realm=realm)
    ELSE
       register_diag_field_scalar = register_static_field(module_name, field_name,&
            & (/null_axis_id/),long_name, units, missing_value, range,&
            & standard_name=standard_name, do_not_log=do_not_log, realm=realm)
    END IF
  END FUNCTION register_diag_field_scalar

  !> @brief Registers an array field
  !> @return field index for subsequent call to send_data.
  INTEGER FUNCTION register_diag_field_array(module_name, field_name, axes, init_time, &
       & long_name, units, missing_value, range, mask_variant, standard_name, verbose,&
       & do_not_log, err_msg, interp_method, tile_count, area, volume, realm)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    INTEGER, INTENT(in) :: axes(:)
    TYPE(time_type), INTENT(in) :: init_time
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name, units, standard_name
    REAL, OPTIONAL, INTENT(in) :: missing_value, RANGE(2)
    LOGICAL, OPTIONAL, INTENT(in) :: mask_variant,verbose
    LOGICAL, OPTIONAL, INTENT(in) :: do_not_log !< if TRUE, field info is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(out):: err_msg
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method !< The interp method to be used when
                                                            !! regridding the field in post-processing.
                                                            !! Valid options are "conserve_order1",
                                                            !! "conserve_order2", and "none".
    INTEGER, OPTIONAL, INTENT(in) :: tile_count
    INTEGER, OPTIONAL, INTENT(in) :: area !< diag_field_id containing the cell area field
    INTEGER, OPTIONAL, INTENT(in) :: volume !< diag_field_id containing the cell volume field
    CHARACTER(len=*), OPTIONAL, INTENT(in):: realm !< String to set as the value to the modeling_realm attribute

    INTEGER :: field, j, ind, file_num, freq
    INTEGER :: i, cm_ind, cm_file_num
    INTEGER :: output_units
    INTEGER :: stdout_unit
    LOGICAL :: mask_variant1, verbose1
    LOGICAL :: cm_found
    CHARACTER(len=128) :: msg

    ! get stdout unit number
    stdout_unit = stdout()

    IF ( PRESENT(mask_variant) ) THEN
       mask_variant1 = mask_variant
    ELSE
       mask_variant1 = .FALSE.
    END IF

    IF ( PRESENT(verbose) ) THEN
       verbose1 = verbose
    ELSE
       verbose1 = .FALSE.
    END IF

    IF ( PRESENT(err_msg) ) err_msg = ''

    ! Call register static, then set static back to false
    register_diag_field_array = register_static_field(module_name, field_name, axes,&
         & long_name, units, missing_value, range, mask_variant1, standard_name=standard_name,&
         & DYNAMIC=.TRUE., do_not_log=do_not_log, interp_method=interp_method, tile_count=tile_count, realm=realm)

    IF ( .NOT.first_send_data_call ) THEN
       ! <ERROR STATUS="WARNING">
       !   module/output_field <module_name>/<field_name> registered AFTER first
       !   send_data call, TOO LATE
       ! </ERROR>
       IF ( mpp_pe() == mpp_root_pe() ) &
            & CALL  error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '&
            &//TRIM(module_name)//'/'// TRIM(field_name)//&
            &' registered AFTER first send_data call, TOO LATE', WARNING)
    END IF

    IF ( register_diag_field_array < 0 ) THEN
       ! <ERROR STATUS="WARNING">
       !   module/output_field <modul_name>/<field_name> NOT found in diag_table
       ! </ERROR>
       IF ( debug_diag_manager .OR. verbose1 ) THEN
          IF ( mpp_pe() == mpp_root_pe() ) &
               & CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '&
               &//TRIM(module_name)//'/'// TRIM(field_name)//' NOT found in diag_table',&
               & WARNING)
       END IF
    ELSE
       input_fields(register_diag_field_array)%static = .FALSE.
       field = register_diag_field_array


       ! Verify that area and volume do not point to the same variable
       IF ( PRESENT(volume).AND.PRESENT(area) ) THEN
          IF ( area.EQ.volume ) THEN
             IF (PRESENT(err_msg)) THEN
                err_msg = 'diag_manager_mod::register_diag_field: module/output_field '&
                  &//TRIM(module_name)//'/'// TRIM(field_name)//' AREA and VOLUME CANNOT be the same variable.&
                  & Contact the developers.'
                register_diag_field_array = -1
                RETURN
             ELSE
                CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '&
                  &//TRIM(module_name)//'/'// TRIM(field_name)//' AREA and VOLUME CANNOT be the same variable.&
                  & Contact the developers.',&
                  & FATAL)
             ENDIF
          END IF
       END IF

       ! Check for the existence of the area/volume field(s)
       IF ( PRESENT(area) ) THEN
          IF ( area < 0 ) THEN
             IF (PRESENT(err_msg)) THEN
                err_msg = 'diag_manager_mod::register_diag_field: module/output_field '&
                  &//TRIM(module_name)//'/'// TRIM(field_name)//' AREA measures field NOT found in diag_table.&
                  & Contact the model liaison.'
                register_diag_field_array = -1
                RETURN
             ELSE
                CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '&
                  &//TRIM(module_name)//'/'// TRIM(field_name)//' AREA measures field NOT found in diag_table.&
                  & Contact the model liaison.',&
                  & FATAL)
             ENDIF
          END IF
       END IF
       IF ( PRESENT(volume) ) THEN
          IF ( volume < 0 ) THEN
             IF (PRESENT(err_msg)) THEN
                err_msg = 'diag_manager_mod::register_diag_field: module/output_field '&
                  &//TRIM(module_name)//'/'// TRIM(field_name)//' VOLUME measures field NOT found in diag_table.&
                  & Contact the model liaison.'
                register_diag_field_array = -1
                RETURN
             ELSE
                CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '&
                  &//TRIM(module_name)//'/'// TRIM(field_name)//' VOLUME measures field NOT found in diag_table.&
                  & Contact the model liaison.',&
                  & FATAL)
             ENDIF
          END IF
       END IF

       IF ( PRESENT(standard_name) ) input_fields(field)%standard_name = standard_name

       DO j = 1, input_fields(field)%num_output_fields
          ind = input_fields(field)%output_fields(j)
          output_fields(ind)%static = .FALSE.
          ! Set up times in output_fields
          output_fields(ind)%last_output = init_time
          ! Get output frequency from for the appropriate output file
          file_num = output_fields(ind)%output_file
          IF ( file_num == max_files ) CYCLE
          IF ( output_fields(ind)%local_output ) THEN
             IF ( output_fields(ind)%need_compute) THEN
                files(file_num)%local = .TRUE.
             END IF
          END IF

          ! Need to sync start_time of file with init time of model
          ! and close_time calculated with the duration of the file.
          ! Also, increase next_open until it is greater than init_time.
          CALL sync_file_times(file_num, init_time, err_msg=msg)
          IF ( msg /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::register_diag_field', TRIM(msg), err_msg) ) RETURN
          END IF

          freq = files(file_num)%output_freq
          output_units = files(file_num)%output_units
          output_fields(ind)%next_output = diag_time_inc(init_time, freq, output_units, err_msg=msg)
          IF ( msg /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::register_diag_field',&
                  & ' file='//TRIM(files(file_num)%name)//': '//TRIM(msg),err_msg)) RETURN
          END IF
          output_fields(ind)%next_next_output = &
               & diag_time_inc(output_fields(ind)%next_output, freq, output_units, err_msg=msg)
          IF ( msg /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::register_diag_field',&
                  &' file='//TRIM(files(file_num)%name)//': '//TRIM(msg),err_msg) ) RETURN
          END IF
          IF ( debug_diag_manager .AND. mpp_pe() == mpp_root_pe() .AND. output_fields(ind)%local_output ) THEN
             WRITE (msg,'(" lon(",F5.1,", ",F5.1,"), lat(",F5.1,", ",F5.1,"), dep(",F5.1,", ",F5.1,")")') &
                  & output_fields(ind)%output_grid%start(1),output_fields(ind)%output_grid%end(1),&
                  & output_fields(ind)%output_grid%start(2),output_fields(ind)%output_grid%end(2),&
                  & output_fields(ind)%output_grid%start(3),output_fields(ind)%output_grid%end(3)
             WRITE(stdout_unit,* ) 'module/output_field '//TRIM(module_name)//'/'//TRIM(field_name)// &
                  & ' will be output in region:'//TRIM(msg)
          END IF

          ! Set the cell_measures attribute in the out file
          CALL init_field_cell_measures(output_fields(ind), area=area, volume=volume, err_msg=err_msg)
          IF ( LEN_TRIM(err_msg).GT.0 ) THEN
             CALL error_mesg ('diag_manager_mod::register_diag_field',&
                  & TRIM(err_msg)//' for module/field '//TRIM(module_name)//'/'//TRIM(field_name),&
                  & FATAL)
          END IF

       END DO
    END IF
  END FUNCTION register_diag_field_array

  !> @brief Return field index for subsequent call to send_data.
  !! @return field index for subsequent call to send_data.
  INTEGER FUNCTION register_static_field(module_name, field_name, axes, long_name, units,&
       & missing_value, range, mask_variant, standard_name, DYNAMIC, do_not_log, interp_method,&
       & tile_count, area, volume, realm)
    CHARACTER(len=*), INTENT(in) :: module_name, field_name
    INTEGER, DIMENSION(:), INTENT(in) :: axes
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: long_name, units, standard_name
    REAL, OPTIONAL, INTENT(in) :: missing_value
    REAL, DIMENSION(2), OPTIONAL, INTENT(in) :: range
    LOGICAL, OPTIONAL, INTENT(in) :: mask_variant
    LOGICAL, OPTIONAL, INTENT(in) :: DYNAMIC
    LOGICAL, OPTIONAL, INTENT(in) :: do_not_log !< if TRUE, field information is not logged
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method !< The interp method to be used when
                                                            !! regridding the field in post-processing.
                                                            !! Valid options are "conserve_order1",
                                                            !! "conserve_order2", and "none".
    INTEGER,          OPTIONAL, INTENT(in) :: tile_count
    INTEGER,          OPTIONAL, INTENT(in) :: area !< Field ID for the area field associated with this field
    INTEGER,          OPTIONAL, INTENT(in) :: volume !< Field ID for the volume field associated with this field
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: realm !< String to set as the value to the modeling_realm attribute

    REAL :: missing_value_use
    INTEGER :: field, num_axes, j, out_num, k
    INTEGER, DIMENSION(3) :: siz, local_siz, local_start, local_end ! indices of local domain of global axes
    INTEGER :: tile, file_num
    LOGICAL :: mask_variant1, dynamic1, allow_log
    CHARACTER(len=128) :: msg
    INTEGER :: domain_type

    ! Fatal error if the module has not been initialized.
    IF ( .NOT.module_is_initialized ) THEN
       ! <ERROR STATUS="FATAL">diag_manager has NOT been initialized</ERROR>
       CALL error_mesg ('diag_manager_mod::register_static_field', 'diag_manager has NOT been initialized', FATAL)
    END IF

    ! Check if OPTIONAL parameters were passed in.
    IF ( PRESENT(missing_value) ) THEN
       IF ( use_cmor ) THEN
          missing_value_use = CMOR_MISSING_VALUE
       ELSE
          missing_value_use = missing_value
       END IF
    END IF

    IF ( PRESENT(mask_variant) ) THEN
       mask_variant1 = mask_variant
    ELSE
       mask_variant1 = .FALSE.
    END IF

    IF ( PRESENT(DYNAMIC) ) THEN
       dynamic1 = DYNAMIC
    ELSE
       dynamic1 = .FALSE.
    END IF

    IF ( PRESENT(tile_count) ) THEN
       tile = tile_count
    ELSE
       tile = 1
    END IF

    IF ( PRESENT(do_not_log) ) THEN
       allow_log = .NOT.do_not_log
    ELSE
       allow_log = .TRUE.
    END IF

    ! Namelist do_diag_field_log is by default false.  Thus to log the
    ! registration of the data field, but the OPTIONAL parameter
    ! do_not_log == .FALSE. and the namelist variable
    ! do_diag_field_log == .TRUE..
    IF ( do_diag_field_log.AND.allow_log ) THEN
       CALL log_diag_field_info (module_name, field_name, axes, &
            & long_name, units, missing_value=missing_value, range=range, &
            & DYNAMIC=dynamic1)
    END IF

    register_static_field = find_input_field(module_name, field_name, 1)
    field = register_static_field
    ! Negative index returned if this field was not found in the diag_table.
    IF ( register_static_field < 0 ) RETURN

    ! Check that the axes are compatible with each other
    domain_type = axis_compatible_check(axes,field_name)

    IF ( tile > 1 ) THEN
       IF ( .NOT.input_fields(field)%register ) THEN
          ! <ERROR STATUS="FATAL">
          !   module/output_field <module_name>/<field_name> is not registered for tile_count = 1,
          !   should not register for tile_count > 1
          ! </ERROR>
          CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '//trim(module_name)//'/'//&
           & TRIM(field_name)//' is not registered for tile_count = 1, should not register for tile_count > 1',&
           & FATAL)
       END IF

       CALL init_input_field(module_name, field_name, tile)
       register_static_field = find_input_field(module_name, field_name, tile)
       DO j = 1, input_fields(field)%num_output_fields
          out_num = input_fields(field)%output_fields(j)
          file_num = output_fields(out_num)%output_file
          IF(input_fields(field)%local) THEN
             CALL init_output_field(module_name, field_name,output_fields(out_num)%output_name,&
                  & files(file_num)%name,output_fields(out_num)%time_method, output_fields(out_num)%pack,&
                  & tile, input_fields(field)%local_coord)
          ELSE
             CALL init_output_field(module_name, field_name,output_fields(out_num)%output_name,&
                  & files(file_num)%name,output_fields(out_num)%time_method, output_fields(out_num)%pack, tile)
          END IF
       END DO
       field = register_static_field
    END IF

    ! Store information for this input field into input field table

    ! Set static to true, if called by register_diag_field this is
    ! flipped back to false
    input_fields(field)%static = .TRUE.
    ! check if the field is registered twice
    IF ( input_fields(field)%register .AND. mpp_pe() == mpp_root_pe() ) THEN
       ! <ERROR STATUS="FATAL">
       !   module/output_field <module_name>/<field_name> ALREADY Registered, should
       !   not register twice
       ! </ERROR>
       CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '//trim(module_name)//'/'//&
            & TRIM(field_name)//' ALREADY registered, should not register twice', FATAL)
    END IF

    ! Verify that area and volume do not point to the same variable
    IF ( PRESENT(volume).AND.PRESENT(area) ) THEN
       IF ( area.EQ.volume ) THEN
          CALL error_mesg ('diag_manager_mod::register_static_field', 'module/output_field '&
               &//TRIM(module_name)//'/'// TRIM(field_name)//' AREA and VOLUME CANNOT be the same variable.&
               & Contact the developers.',&
               & FATAL)
       END IF
    END IF

    ! Check for the existence of the area/volume field(s)
    IF ( PRESENT(area) ) THEN
       IF ( area < 0 ) THEN
          CALL error_mesg ('diag_manager_mod::register_static_field', 'module/output_field '&
               &//TRIM(module_name)//'/'// TRIM(field_name)//' AREA measures field NOT found in diag_table.&
               & Contact the model liaison.n',&
               & FATAL)
       END IF
    END IF
    IF ( PRESENT(volume) ) THEN
       IF ( volume < 0 ) THEN
          CALL error_mesg ('diag_manager_mod::register_static_field', 'module/output_field '&
               &//TRIM(module_name)//'/'// TRIM(field_name)//' VOLUME measures field NOT found in diag_table&
               & Contact the model liaison.',&
               & FATAL)
       END IF
    END IF

    ! Set flag that this field was registered
    input_fields(field)%register = .TRUE.
    ! set flag for mask: does it change with time?
    input_fields(field)%mask_variant = mask_variant1
    ! Set flag for mask warning
    input_fields(field)%issued_mask_ignore_warning = .FALSE.

    ! Check for more OPTIONAL parameters.
    IF ( PRESENT(long_name) ) THEN
       input_fields(field)%long_name = TRIM(long_name)
    ELSE
       input_fields(field)%long_name = input_fields(field)%field_name
    END IF

    IF ( PRESENT(standard_name) ) input_fields(field)%standard_name = standard_name

    IF ( PRESENT(units) ) THEN
       input_fields(field)%units = TRIM(units)
    ELSE
       input_fields(field)%units = 'none'
    END IF

    IF ( PRESENT(missing_value) ) THEN
       input_fields(field)%missing_value = missing_value_use
       input_fields(field)%missing_value_present = .TRUE.
    ELSE
       input_fields(field)%missing_value_present = .FALSE.
    END IF

    IF ( PRESENT(range) ) THEN
       input_fields(field)%range = range
       ! don't use the range if it is not a valid range
       input_fields(field)%range_present = range(2) .gt. range(1)
    ELSE
       input_fields(field)%range = (/ 1., 0. /)
       input_fields(field)%range_present = .FALSE.
    END IF

    IF ( PRESENT(interp_method) ) THEN
       IF ( TRIM(interp_method) .NE. 'conserve_order1' .AND.&
            & TRIM(interp_method) .NE. 'conserve_order2' .AND.&
            & TRIM(interp_method) .NE. 'none' ) THEN
          ! <ERROR STATUS="FATAL">
          !   when registering module/output_field <module_name>/<field_name> then optional
          !   argument interp_method = <interp_method>, but it should be "conserve_order1",
          !   "conserve_order2", or "none"
          ! </ERROR>
          CALL error_mesg ('diag_manager_mod::register_diag_field',&
               & 'when registering module/output_field '//TRIM(module_name)//'/'//&
               & TRIM(field_name)//', the optional argument interp_method = '//TRIM(interp_method)//&
               & ', but it should be "conserve_order1", "conserve_order2", or "none"', FATAL)
       END IF
       input_fields(field)%interp_method = TRIM(interp_method)
    ELSE
       input_fields(field)%interp_method = ''
    END IF

    ! Store the axis info
    num_axes = SIZE(axes(:)) ! num_axes should be <= 3.
    input_fields(field)%axes(1:num_axes) = axes
    input_fields(field)%num_axes = num_axes

    siz = 1
    DO j = 1, num_axes
       IF ( axes(j) .LE. 0 ) THEN
          ! <ERROR STATUS="FATAL">
          !   module/output_field <module_name>/<field_name> has non-positive axis_id
          ! </ERROR>
          CALL error_mesg ('diag_manager_mod::register_diag_field', 'module/output_field '//trim(module_name)//'/'//&
               & TRIM(field_name)//' has non-positive axis_id', FATAL)
       END IF
       siz(j) = get_axis_length(axes(j))
    END DO

    ! Default length for axes is 1
    DO j = 1, 3
       input_fields(field)%size(j) = siz(j)
    END DO

    local_siz = 1
    local_start = 1
    local_end= 1
    ! Need to loop through all output_fields associated and allocate their buffers
    DO j = 1, input_fields(field)%num_output_fields
       out_num = input_fields(field)%output_fields(j)
       ! Range is required when pack >= 4
       IF ( output_fields(out_num)%pack>=4 .AND. .NOT.input_fields(field)%range_present ) THEN
          IF(mpp_pe() .EQ. mpp_root_pe()) THEN
             ! <ERROR STATUS="FATAL">
             !   output_field <field_name> has pack >= 4, range is REQUIRED in register_diag_field
             ! </ERROR>
             CALL error_mesg ('diag_manager_mod::register_diag_field ', 'output_field '//TRIM(field_name)// &
                  ' has pack >=4, range is REQUIRED in register_diag_field', FATAL)
          END IF
       END IF
       ! reset the number of diurnal samples to 1 if the field is static (and, therefore,
       ! doesn't vary diurnally)
       IF ( .NOT.dynamic1 ) output_fields(out_num)%n_diurnal_samples = 1

       !Check that the domain associated with the inputted field matches
       !the domain associated output files to which it will be written.
       file_num = output_fields(out_num)%output_file
       if (domain_type .eq. DIAG_AXIS_2DDOMAIN) then
           if (files(file_num)%use_domainUG) then
               call error_mesg("diag_manager_mod::register_static_field", &
                               "Diagnostics living on a structured grid" &
                               //" and an unstructured grid cannot exist" &
                               //" in the same file (" &
                               //trim(files(file_num)%name)//")", &
                               FATAL)
           elseif (.not. files(file_num)%use_domain2D) then
               files(file_num)%use_domain2D = .true.
           endif
       elseif (domain_type .eq. DIAG_AXIS_UGDOMAIN) then
           if (files(file_num)%use_domain2D) then
               call error_mesg("diag_manager_mod::register_static_field", &
                               "Diagnostics living on a structured grid" &
                               //" and an unstructured grid cannot exist" &
                               //" in the same file (" &
                               //trim(files(file_num)%name)//")", &
                               FATAL)
           elseif (.not. files(file_num)%use_domainUG) then
               files(file_num)%use_domainUG = .true.
           endif
       endif


       !  if local_output (size of output_fields does NOT equal size of input_fields)
       IF ( output_fields(out_num)%reduced_k_range ) THEN
          CALL get_subfield_vert_size(axes, out_num)

!----------
!ug support
         !Send_data requires that the reduced k dimension be the 3rd dimension
         !of the buffer, so set it to be the correct size.  If the diagnostic
         !is unstructured, set the second dimension of the buffer to be 1.
          if (domain_type .eq. DIAG_AXIS_UGDOMAIN) then
              local_start(2) = output_fields(out_num)%output_grid%l_start_indx(2)
              local_end(2) = output_fields(out_num)%output_grid%l_end_indx(2)
              local_siz(2) = local_end(2) - local_start(2) + 1
              allocate(output_fields(out_num)%buffer(siz(1),local_siz(2),siz(3), &
                                                     output_fields(out_num)%n_diurnal_samples))
              output_fields(out_num)%region_elements = siz(1)*local_siz(2)*siz(3)
              output_fields(out_num)%reduced_k_unstruct = .true.
          else
              local_start(3) = output_fields(out_num)%output_grid%l_start_indx(3)
              local_end(3) = output_fields(out_num)%output_grid%l_end_indx(3)
              local_siz(3) = local_end(3) - local_start(3) + 1
              allocate(output_fields(out_num)%buffer(siz(1),siz(2),local_siz(3), &
                                                     output_fields(out_num)%n_diurnal_samples))
              output_fields(out_num)%region_elements = siz(1)*siz(2)*local_siz(3)
              output_fields(out_num)%reduced_k_unstruct = .false.
          endif
          output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
!----------

          IF ( output_fields(out_num)%time_max ) THEN
             output_fields(out_num)%buffer = MAX_VALUE
          ELSE IF ( output_fields(out_num)%time_min ) THEN
             output_fields(out_num)%buffer = MIN_VALUE
          ELSE
             output_fields(out_num)%buffer = EMPTY
          END IF
       ELSE IF ( output_fields(out_num)%local_output ) THEN
          IF ( SIZE(axes(:)) .LE. 1 ) THEN
             ! <ERROR STATUS="FATAL">axes of <field_name> must >= 2 for local output</ERROR>
             CALL error_mesg ('diag_manager_mod::register_diag_field', 'axes of '//TRIM(field_name)//&
                  & ' must >= 2 for local output', FATAL)
          END IF
          CALL get_subfield_size(axes, out_num)
          IF ( output_fields(out_num)%need_compute ) THEN
             DO k = 1, num_axes
                local_start(k) = output_fields(out_num)%output_grid%l_start_indx(k)
                local_end(k) = output_fields(out_num)%output_grid%l_end_indx(k)
                local_siz(k) = local_end(k) - local_start(k) +1
             END DO
             ALLOCATE(output_fields(out_num)%buffer(local_siz(1), local_siz(2), local_siz(3),&
                  & output_fields(out_num)%n_diurnal_samples))
             IF(output_fields(out_num)%time_max) THEN
                output_fields(out_num)%buffer = MAX_VALUE
             ELSE IF(output_fields(out_num)%time_min) THEN
                output_fields(out_num)%buffer = MIN_VALUE
             ELSE
                output_fields(out_num)%buffer = EMPTY
             END IF
             output_fields(out_num)%region_elements = local_siz(1)*local_siz(2)*local_siz(3)
             output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
             files(output_fields(out_num)%output_file)%local = .true.
          END IF
       ELSE ! the field is output globally
          ! size of output_fields equal size of input_fields
          ALLOCATE(output_fields(out_num)%buffer(siz(1), siz(2), siz(3),&
               & output_fields(out_num)%n_diurnal_samples))
          IF(output_fields(out_num)%time_max) THEN
             output_fields(out_num)%buffer = MAX_VALUE
          ELSE IF(output_fields(out_num)%time_min) THEN
             output_fields(out_num)%buffer = MIN_VALUE
          ELSE
             output_fields(out_num)%buffer = EMPTY
          END IF
          output_fields(out_num)%total_elements = siz(1)*siz(2)*siz(3)
       END IF

       ! Reset to false in register_field if this is not static
       output_fields(out_num)%static = .TRUE.
       ! check if time average is true for static field
       IF ( .NOT.dynamic1 .AND. output_fields(out_num)%time_ops ) THEN
          WRITE (msg,'(a,"/",a)') TRIM(module_name), TRIM(field_name)
          IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
             ! <ERROR STATUS="WARNING">
             !   module/field <module_name>/<field_name> is STATIC.
             !   Cannot perform time operations average, maximum or
             !   minimum on static fields.  Setting the time operation to 'NONE'
             !   for this field.
             ! </ERROR>
             CALL error_mesg ('diag_manager_mod::register_static_field',&
                  & 'module/field '//TRIM(msg)//' is STATIC.  Cannot perform time operations&
                  & average, maximum, or minimum on static fields.  Setting the time operation&
                  & to "NONE" for this field.', WARNING)
          END IF
          output_fields(out_num)%time_ops = .FALSE.
          output_fields(out_num)%time_average = .FALSE.
          output_fields(out_num)%time_method = 'point'
       END IF

       ! assume that the number of axes of output_fields = that of input_fields
       ! this should be changed later to take into account time-of-day axis
       output_fields(out_num)%num_axes = input_fields(field)%num_axes
       ! Axes are copied from input_fields if output globally or from subaxes if output locally
       IF ( .NOT.output_fields(out_num)%local_output ) THEN
          output_fields(out_num)%axes(1:input_fields(field)%num_axes) =&
               & input_fields(field)%axes(1:input_fields(field)%num_axes)
       ELSE
          output_fields(out_num)%axes(1:input_fields(field)%num_axes) =&
               & output_fields(out_num)%output_grid%subaxes(1:input_fields(field)%num_axes)
       END IF

       ! if necessary, initialize the diurnal time axis and append its index in the
       ! output field axes array
       IF ( output_fields(out_num)%n_diurnal_samples > 1 ) THEN
          output_fields(out_num)%axes(output_fields(out_num)%num_axes+1) =&
               & init_diurnal_axis(output_fields(out_num)%n_diurnal_samples)
          output_fields(out_num)%num_axes = output_fields(out_num)%num_axes+1
       END IF

       IF ( output_fields(out_num)%reduced_k_range ) THEN
!----------
!ug support
          if (domain_type .eq. DIAG_AXIS_UGDOMAIN) then
              output_fields(out_num)%axes(2) = output_fields(out_num)%output_grid%subaxes(2)
          else
              output_fields(out_num)%axes(3) = output_fields(out_num)%output_grid%subaxes(3)
          endif
!----------
       END IF

       ! Initialize a time variable used in an error check
       output_fields(out_num)%Time_of_prev_field_data = Time_zero

       ! Set the cell_measures attribute in the out file
       CALL init_field_cell_measures(output_fields(out_num), area=area, volume=volume, err_msg=msg)
       IF ( LEN_TRIM(msg).GT.0 ) THEN
          CALL error_mesg ('diag_manager_mod::register_static_field',&
               & TRIM(msg)//' for module/field '//TRIM(module_name)//'/'//TRIM(field_name),&
               & FATAL)
       END IF

       ! Add the modeling_realm attribute
       IF ( PRESENT(realm) ) THEN
          CALL prepend_attribute(output_fields(out_num), 'modeling_realm', lowercase(TRIM(realm)))
       END IF
    END DO

    IF ( input_fields(field)%mask_variant ) THEN
       DO j = 1, input_fields(field)%num_output_fields
          out_num = input_fields(field)%output_fields(j)
          IF(output_fields(out_num)%time_average) THEN
!----------
!ug support
             !Send_data requires that the reduced k dimension be the 3rd dimension
             !of the counter array, so set it to be the correct size.  If the diagnostic
             !is unstructured, set the second dimension of the counter array to be 1.
              if (output_fields(out_num)%reduced_k_range .and. &
                  domain_type .eq. DIAG_AXIS_UGDOMAIN) then
                  allocate(output_fields(out_num)%counter(siz(1),local_siz(2),siz(3), &
                                                          output_fields(out_num)%n_diurnal_samples))
              else
                  allocate(output_fields(out_num)%counter(siz(1),siz(2),siz(3), &
                                                          output_fields(out_num)%n_diurnal_samples))
              endif
!----------
             output_fields(out_num)%counter = 0.0
          END IF
       END DO
    END IF
  END FUNCTION register_static_field

  !> @brief Return the diagnostic field ID of a given variable.
  !! @return get_diag_field_id will return the ID returned during the register_diag_field call.
  !!   If the variable is not in the diag_table, then the value "DIAG_FIELD_NOT_FOUND" will be
  !!   returned.
  INTEGER FUNCTION get_diag_field_id(module_name, field_name)
    CHARACTER(len=*), INTENT(in) :: module_name !< Module name that registered the variable
    CHARACTER(len=*), INTENT(in) :: field_name !< Variable name

    ! find_input_field will return DIAG_FIELD_NOT_FOUND if the field is not
    ! included in the diag_table
    get_diag_field_id = find_input_field(module_name, field_name, tile_count=1)
  END FUNCTION get_diag_field_id

  !> @brief Finds the corresponding related output field and file for a given input field
  !! @return Logical get_related_field
  LOGICAL FUNCTION get_related_field(field, rel_field, out_field_id, out_file_id)
    INTEGER, INTENT(in) :: field !< input field ID to find the corresponding
    TYPE(output_field_type), INTENT(in) :: rel_field !< Output field that field must correspond to
    INTEGER, INTENT(out) :: out_field_id !< output_field index of related output field
    INTEGER, INTENT(out) :: out_file_id !< file index of the out_field_id output field

    INTEGER :: i, cm_ind, cm_file_num
    INTEGER :: rel_file

    ! Output file index of field to compare to
    rel_file = rel_field%output_file

    ! Default return values
    out_field_id = -1
    out_file_id = -1
    get_related_field = .FALSE.

    ! First check if any fields are in the same file as rel_field
    DO i = 1, input_fields(field)%num_output_fields
       cm_ind = input_fields(field)%output_fields(i)
       cm_file_num = output_fields(cm_ind)%output_file

       IF ( cm_file_num.EQ.rel_file.AND.&
            & (( (output_fields(cm_ind)%time_ops.EQV.rel_field%time_ops) .AND.&
            & (output_fields(cm_ind)%next_output.EQ.rel_field%next_output) .AND.&
            & (output_fields(cm_ind)%last_output.EQ.rel_field%last_output) ).OR.&
            & (output_fields(cm_ind)%static.OR.rel_field%static) ) ) THEN
          get_related_field = .TRUE.
          out_field_id = cm_ind
          out_file_id = cm_file_num
          EXIT
       END IF
    END DO

    ! Now look for the field in a different file
    IF ( .NOT.get_related_field ) THEN
       DO i = 1, input_fields(field)%num_output_fields
          cm_ind = input_fields(field)%output_fields(i)
          cm_file_num = output_fields(cm_ind)%output_file

          ! If time_method, freq, output_units, next_output, and last_output the same, or
          ! the output_field is static then valid for cell_measures
!!$ For now, only static fields can be in an external file
!!$          IF ( ( (files(cm_file_num)%output_freq.EQ.files(rel_file)%output_freq) .AND.&
!!$               & (files(cm_file_num)%output_units.EQ.files(rel_file)%output_units) .AND.&
!!$               & (output_fields(cm_ind)%time_ops.EQV.rel_field%time_ops) .AND.&
!!$               & (output_fields(cm_ind)%next_output.EQ.rel_field%next_output) .AND.&
!!$               & (output_fields(cm_ind)%last_output.EQ.rel_field%last_output) ).OR.&
!!$               & ( output_fields(cm_ind)%static.OR.rel_field%static ) ) THEN
          IF ( output_fields(cm_ind)%static.OR.rel_field%static ) THEN
             get_related_field = .TRUE.
             out_field_id = cm_ind
             out_file_id = cm_file_num
             EXIT
          END IF
       END DO
    END IF
  END FUNCTION get_related_field

  !> @brief If needed, add cell_measures and associated_file attribute to out field/file
  SUBROUTINE init_field_cell_measures(output_field, area, volume, err_msg)
    TYPE(output_field_type), INTENT(inout) :: output_field !< Output field that needs the cell_measures
    INTEGER, INTENT(in), OPTIONAL :: area !< Field ID for area
    INTEGER, INTENT(in), OPTIONAL :: volume !< Field ID for volume
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    INTEGER :: cm_ind, cm_file_num, file_num

    IF ( PRESENT(err_msg) ) THEN
       err_msg = ''
    END IF

    ! Verify that area/volume are defined (.gt.0
    IF ( PRESENT(area) ) THEN
       IF ( area.LE.0 ) THEN
          IF ( fms_error_handler('diag_manager_mod::init_field_cell_measure',&
               & 'AREA field not in diag_table for field '//TRIM(input_fields(output_field%input_field)%module_name)//&
               & '/'//TRIM(input_fields(output_field%input_field)%field_name), err_msg) ) RETURN
       END IF
    END IF

    IF ( PRESENT(volume) ) THEN
       IF ( volume.LE.0 ) THEN
          IF ( fms_error_handler('diag_manager_mod::init_field_cell_measure',&
               & 'VOLUME field not in diag_table for field '//TRIM(input_fields(output_field%input_field)%module_name)//&
               & '/'//TRIM(input_fields(output_field%input_field)%field_name), err_msg) ) RETURN
       END IF
    END IF

    ! Get the file number that the output_field will be written to
    file_num = output_field%output_file

    ! Take care of the cell_measures attribute
    IF ( PRESENT(area) ) THEN
       IF ( get_related_field(area, output_field, cm_ind, cm_file_num) ) THEN
          CALL prepend_attribute(output_field, 'cell_measures',&
               & 'area: '//TRIM(output_fields(cm_ind)%output_name))
          IF ( cm_file_num.NE.file_num ) THEN
             ! Not in the same file, set the global attribute associated_files
             CALL add_associated_files(file_num, cm_file_num, cm_ind)
          END IF
       ELSE
          IF ( fms_error_handler('diag_manager_mod::init_field_cell_measures',&
               & 'AREA measures field "'//TRIM(input_fields(area)%module_name)//'/'//&
               & TRIM(input_fields(area)%field_name)//&
               & '" NOT in diag_table with correct output frequency for field '//&
               & TRIM(input_fields(output_field%input_field)%module_name)//&
               & '/'//TRIM(input_fields(output_field%input_field)%field_name), err_msg) ) RETURN
       END IF
    END IF


    IF ( PRESENT(volume) ) THEN
       IF ( get_related_field(volume, output_field, cm_ind, cm_file_num) ) THEN
          CALL prepend_attribute(output_field, 'cell_measures',&
               & 'volume: '//TRIM(output_fields(cm_ind)%output_name))
          IF ( cm_file_num.NE.file_num ) THEN
             ! Not in the same file, set the global attribute associated_files
             CALL add_associated_files(file_num, cm_file_num, cm_ind)
          END IF
       ELSE
          IF ( fms_error_handler('diag_manager_mod::init_field_cell_measures',&
               & 'VOLUME measures field "'//TRIM(input_fields(volume)%module_name)//'/'//&
               & TRIM(input_fields(volume)%field_name)//&
               & '" NOT in diag_table with correct output frequency for field '//&
               & TRIM(input_fields(output_field%input_field)%module_name)//&
               & '/'//TRIM(input_fields(output_field%input_field)%field_name), err_msg) ) RETURN
       END IF
    END IF
  END SUBROUTINE init_field_cell_measures

  !> @brief Add to the associated files attribute
  !!
  !! @throw FATAL, "Length of asso_file_name is not long enough to hold the associated file name."
  !!     The length of character array asso_file_name is not long enough to hold the full file name
  !!     of the associated_file.  Please contact the developer to increase the length of the  variable.
  SUBROUTINE add_associated_files(file_num, cm_file_num, cm_ind)
    INTEGER, intent(in) :: file_num !< File number that needs the associated_files attribute
    INTEGER, intent(in) :: cm_file_num !< file number that contains the associated field
    INTEGER, intent(in) :: cm_ind !< index of the output_field in the associated file

    INTEGER :: year, month, day, hour, minute, second
    INTEGER :: n
    CHARACTER(len=25) :: date_prefix
    CHARACTER(len=256) :: asso_file_name

    ! Create the date_string
    IF ( prepend_date ) THEN
       CALL get_date(diag_init_time, year, month, day, hour, minute, second)
       WRITE (date_prefix, '(1I20.4, 2I2.2,".")') year, month, day
       date_prefix=ADJUSTL(date_prefix)
    ELSE
       date_prefix=''
    END IF

    ! Get the base file name
    ! Verify asso_file_name is long enough to hold the file name,
    ! plus 17 for the additional '.ens_??.tile?.nc' (and a null character)
    IF ( LEN_TRIM(files(cm_file_num)%name)+17 > LEN(asso_file_name) ) THEN
       CALL error_mesg ('diag_manager_mod::add_associated_files',&
            & 'Length of asso_file_name is not long enough to hold the associated file name. '&
            & //'Contact the developer', FATAL)
    ELSE
       asso_file_name = TRIM(files(cm_file_num)%name)
    END IF

    ! Add the ensemble number string into the file name
    ! As frepp does not have native support for multiple ensemble runs
    ! this will not be done.  However, the code is left here for the time
    ! frepp does.
    !CALL get_instance_filename(TRIM(asso_file_name), asso_file_name)

    ! Append .nc suffix, if needed. Note that we no longer try to append cubic sphere tile
    ! number to the name of the associated file.
    n = max(len_trim(asso_file_name),3)
    if (asso_file_name(n-2:n).NE.'.nc') asso_file_name = trim(asso_file_name)//'.nc'

    ! Should look like :associated_files = " output_name: output_file_name " ;
    CALL prepend_attribute(files(file_num), 'associated_files',&
         & TRIM(output_fields(cm_ind)%output_name)//': '//&
         & TRIM(date_prefix)//TRIM(asso_file_name))
  END SUBROUTINE add_associated_files

  !> @return true if send is successful
  LOGICAL FUNCTION send_data_0d(diag_field_id, field, time, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, INTENT(in) :: field
    TYPE(time_type), INTENT(in), OPTIONAL :: time
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL :: field_out(1, 1, 1)

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_0d = .FALSE.
       RETURN
    END IF
    ! First copy the data to a three d array with last element 1
    field_out(1, 1, 1) = field
    send_data_0d = send_data_3d(diag_field_id, field_out, time, err_msg=err_msg)
  END FUNCTION send_data_0d

  !> @return true if send is successful
  LOGICAL FUNCTION send_data_1d(diag_field_id, field, time, is_in, mask, rmask, ie_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:), INTENT(in) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    REAL, INTENT(in), DIMENSION(:), OPTIONAL :: rmask
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, ie_in
    LOGICAL, INTENT(in), DIMENSION(:), OPTIONAL :: mask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field(:)), 1, 1) :: field_out
    LOGICAL, DIMENSION(SIZE(field(:)), 1, 1) ::  mask_out

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_1d = .FALSE.
       RETURN
    END IF

    ! First copy the data to a three d array with last element 1
    field_out(:, 1, 1) = field

    ! Default values for mask
    IF ( PRESENT(mask) ) THEN
       mask_out(:, 1, 1) = mask
    ELSE
       mask_out = .TRUE.
    END IF

    IF ( PRESENT(rmask) ) WHERE (rmask < 0.5) mask_out(:, 1, 1) = .FALSE.
    IF ( PRESENT(mask) .OR. PRESENT(rmask) ) THEN
       IF ( PRESENT(is_in) .OR. PRESENT(ie_in) ) THEN
          send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=1, ks_in=1,&
               & mask=mask_out, ie_in=ie_in, je_in=1, ke_in=1, weight=weight, err_msg=err_msg)
       ELSE
          send_data_1d = send_data_3d(diag_field_id, field_out, time, mask=mask_out,&
               & weight=weight, err_msg=err_msg)
       END IF
    ELSE
       IF ( PRESENT(is_in) .OR. PRESENT(ie_in) ) THEN
          send_data_1d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=1, ks_in=1,&
               & ie_in=ie_in, je_in=1, ke_in=1, weight=weight, err_msg=err_msg)
       ELSE
          send_data_1d = send_data_3d(diag_field_id, field_out, time, weight=weight, err_msg=err_msg)
       END IF
    END IF
  END FUNCTION send_data_1d

  !> @return true if send is successful
  LOGICAL FUNCTION send_data_2d(diag_field_id, field, time, is_in, js_in, &
       & mask, rmask, ie_in, je_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, INTENT(in), DIMENSION(:,:) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ie_in, je_in
    LOGICAL, INTENT(in), DIMENSION(:,:), OPTIONAL :: mask
    REAL, INTENT(in), DIMENSION(:,:),OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) :: field_out
    LOGICAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) ::  mask_out

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_2d = .FALSE.
       RETURN
    END IF

    ! First copy the data to a three d array with last element 1
    field_out(:, :, 1) = field

    ! Default values for mask
    IF ( PRESENT(mask) ) THEN
       mask_out(:, :, 1) = mask
    ELSE
       mask_out = .TRUE.
    END IF

    IF ( PRESENT(rmask) ) WHERE ( rmask < 0.5 ) mask_out(:, :, 1) = .FALSE.
    IF ( PRESENT(mask) .OR. PRESENT(rmask) ) THEN
       send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=1, mask=mask_out,&
            & ie_in=ie_in, je_in=je_in, ke_in=1, weight=weight, err_msg=err_msg)
    ELSE
       send_data_2d = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=1,&
            & ie_in=ie_in, je_in=je_in, ke_in=1, weight=weight, err_msg=err_msg)
    END IF
  END FUNCTION send_data_2d

#ifdef OVERLOAD_R8
  !> @return true if send is successful
  LOGICAL FUNCTION send_data_2d_r8(diag_field_id, field, time, is_in, js_in, &
       & mask, rmask, ie_in, je_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL(kind=8), INTENT(in), DIMENSION(:,:) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ie_in, je_in
    LOGICAL, INTENT(in), DIMENSION(:,:), OPTIONAL :: mask
    REAL, INTENT(in), DIMENSION(:,:),OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) :: field_out
    LOGICAL, DIMENSION(SIZE(field,1),SIZE(field,2),1) ::  mask_out

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_2d_r8 = .FALSE.
       RETURN
    END IF

    ! First copy the data to a three d array with last element 1
    field_out(:, :, 1) = field

    ! Default values for mask
    IF ( PRESENT(mask) ) THEN
       mask_out(:, :, 1) = mask
    ELSE
       mask_out = .TRUE.
    END IF

    IF ( PRESENT(rmask) ) WHERE ( rmask < 0.5 ) mask_out(:, :, 1) = .FALSE.
    IF ( PRESENT(mask) .OR. PRESENT(rmask) ) THEN
       send_data_2d_r8 = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=1, mask=mask_out,&
            & ie_in=ie_in, je_in=je_in, ke_in=1, weight=weight, err_msg=err_msg)
    ELSE
       send_data_2d_r8 = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=1,&
            & ie_in=ie_in, je_in=je_in, ke_in=1, weight=weight, err_msg=err_msg)
    END IF
  END FUNCTION send_data_2d_r8

  !> @return true if send is successful
  LOGICAL FUNCTION send_data_3d_r8(diag_field_id, field, time, is_in, js_in, ks_in, &
             & mask, rmask, ie_in, je_in, ke_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL(kind=8), INTENT(in), DIMENSION(:,:,:) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ks_in,ie_in,je_in, ke_in
    LOGICAL, INTENT(in), DIMENSION(:,:,:), OPTIONAL :: mask
    REAL, INTENT(in), DIMENSION(:,:,:),OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2),size(field,3)) :: field_out
    LOGICAL, DIMENSION(SIZE(field,1),SIZE(field,2),size(field,3)) ::  mask_out

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_3d_r8 = .FALSE.
       RETURN
    END IF

    ! First copy the data to a three d array with last element 1
    field_out = field

    ! Default values for mask
    IF ( PRESENT(mask) ) THEN
       mask_out = mask
    ELSE
       mask_out = .TRUE.
    END IF

    IF ( PRESENT(rmask) ) WHERE ( rmask < 0.5 ) mask_out = .FALSE.
    IF ( PRESENT(mask) .OR. PRESENT(rmask) ) THEN
       send_data_3d_r8 = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=ks_in, mask=mask_out,&
            & ie_in=ie_in, je_in=je_in, ke_in=ke_in, weight=weight, err_msg=err_msg)
    ELSE
       send_data_3d_r8 = send_data_3d(diag_field_id, field_out, time, is_in=is_in, js_in=js_in, ks_in=ks_in,&
            & ie_in=ie_in, je_in=je_in, ke_in=ke_in, weight=weight, err_msg=err_msg)
    END IF
  END FUNCTION send_data_3d_r8
#endif

  !> @return true if send is successful
  LOGICAL FUNCTION send_data_3d(diag_field_id, field, time, is_in, js_in, ks_in, &
             & mask, rmask, ie_in, je_in, ke_in, weight, err_msg)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:,:,:), INTENT(in) :: field
    REAL, INTENT(in), OPTIONAL :: weight
    TYPE (time_type), INTENT(in), OPTIONAL :: time
    INTEGER, INTENT(in), OPTIONAL :: is_in, js_in, ks_in,ie_in,je_in, ke_in
    LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
    REAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: rmask
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    REAL :: weight1
    REAL :: missvalue
    INTEGER :: pow_value
    INTEGER :: ksr, ker
    INTEGER :: i, out_num, file_num, n1, n2, n3, number_of_outputs, ii,f1,f2,f3,f4
    INTEGER :: freq, units, is, js, ks, ie, je, ke, i1, j1,k1, j, k
    INTEGER, DIMENSION(3) :: l_start !< local start indices on 3 axes for regional output
    INTEGER, DIMENSION(3) :: l_end !< local end indices on 3 axes for regional output
    INTEGER   :: hi !< halo size in x direction
    INTEGER   :: hj !< halo size in y direction
    INTEGER   :: twohi !< halo size in x direction
    INTEGER   :: twohj !< halo size in y direction
    INTEGER :: sample !< index along the diurnal time axis
    INTEGER :: day !< components of the current date
    INTEGER :: second !< components of the current date
    INTEGER :: tick !< components of the current date
    INTEGER :: status
    INTEGER :: numthreads
    INTEGER :: active_omp_level
#if defined(_OPENMP)
    INTEGER :: omp_get_num_threads !< OMP function
    INTEGER :: omp_get_level !< OMP function
#endif
    LOGICAL :: average, phys_window, need_compute
    LOGICAL :: reduced_k_range, local_output
    LOGICAL :: time_max, time_min, time_rms, time_sum
    LOGICAL :: missvalue_present
    LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: oor_mask
    CHARACTER(len=256) :: err_msg_local
    CHARACTER(len=128) :: error_string, error_string1

    ! If diag_field_id is < 0 it means that this field is not registered, simply return
    IF ( diag_field_id <= 0 ) THEN
       send_data_3d = .FALSE.
       RETURN
    ELSE
       send_data_3d = .TRUE.
    END IF

    IF ( PRESENT(err_msg) ) err_msg = ''
    IF ( .NOT.module_is_initialized ) THEN
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'diag_manager NOT initialized', err_msg) ) RETURN
    END IF
    err_msg_local = ''
    ! The following lines are commented out as they have not been included in the code prior to now,
    ! and there are a lot of send_data calls before register_diag_field calls.  A method to do this safely
    ! needs to be developed.
    !
    ! Set first_send_data_call to .FALSE. on first non-static field.
!!$    IF ( .NOT.input_fields(diag_field_id)%static .AND. first_send_data_call ) THEN
!!$       first_send_data_call = .FALSE.
!!$    END IF

    ! oor_mask is only used for checking out of range values.
    ALLOCATE(oor_mask(SIZE(field,1),SIZE(field,2),SIZE(field,3)), STAT=status)
    IF ( status .NE. 0 ) THEN
       WRITE (err_msg_local, FMT='("Unable to allocate oor_mask(",I5,",",I5,",",I5,"). (STAT: ",I5,")")')&
            & SIZE(field,1), SIZE(field,2), SIZE(field,3), status
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) RETURN
    END IF

    IF ( PRESENT(mask) ) THEN
       oor_mask = mask
    ELSE
       oor_mask = .TRUE.
    END IF
    IF ( PRESENT(rmask) ) WHERE ( rmask < 0.5 ) oor_mask = .FALSE.

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
    IF ( PRESENT(ie_in) ) THEN
       IF ( .NOT.PRESENT(is_in) ) THEN
          IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'ie_in present without is_in', err_msg) ) THEN
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
       IF ( PRESENT(js_in) .AND. .NOT.PRESENT(je_in) ) THEN
          IF ( fms_error_handler('diag_manager_modsend_data_3d',&
               & 'is_in and ie_in present, but js_in present without je_in', err_msg) ) THEN
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
    END IF
    IF ( PRESENT(je_in) ) THEN
       IF ( .NOT.PRESENT(js_in) ) THEN
          IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'je_in present without js_in', err_msg) ) THEN
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
       IF ( PRESENT(is_in) .AND. .NOT.PRESENT(ie_in) ) THEN
          IF ( fms_error_handler('diag_manager_mod::send_data_3d',&
               & 'js_in and je_in present, but is_in present without ie_in', err_msg)) THEN
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF
    END IF

    ! If is, js, or ks not present default them to 1
    is = 1
    js = 1
    ks = 1
    IF ( PRESENT(is_in) ) is = is_in
    IF ( PRESENT(js_in) ) js = js_in
    IF ( PRESENT(ks_in) ) ks = ks_in
    n1 = SIZE(field, 1)
    n2 = SIZE(field, 2)
    n3 = SIZE(field, 3)
    ie = is+n1-1
    je = js+n2-1
    ke = ks+n3-1
    IF ( PRESENT(ie_in) ) ie = ie_in
    IF ( PRESENT(je_in) ) je = je_in
    IF ( PRESENT(ke_in) ) ke = ke_in
    twohi = n1-(ie-is+1)
    IF ( MOD(twohi,2) /= 0 ) THEN
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'non-symmetric halos in first dimension', err_msg) ) THEN
          DEALLOCATE(oor_mask)
          RETURN
       END IF
    END IF
    twohj = n2-(je-js+1)
    IF ( MOD(twohj,2) /= 0 ) THEN
       IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'non-symmetric halos in second dimension', err_msg) ) THEN
          DEALLOCATE(oor_mask)
          RETURN
       END IF
    END IF
    hi = twohi/2
    hj = twohj/2

    ! The next line is necessary to ensure that is,ie,js,ie are relative to field(1:,1:)
    ! But this works only when there is no windowing.
    IF ( PRESENT(ie_in) .AND. PRESENT(je_in) ) THEN
       is=1+hi
       ie=n1-hi
       js=1+hj
       je=n2-hj
    END IF

    ! used for field, mask and rmask bounds
    f1=1+hi
    f2=n1-hi
    f3=1+hj
    f4=n2-hj

    ! weight is for time averaging where each time level may has a different weight
    IF ( PRESENT(weight) ) THEN
       weight1 = weight
    ELSE
       weight1 = 1.
    END IF

    ! Is there a missing_value?
    missvalue_present = input_fields(diag_field_id)%missing_value_present
    IF ( missvalue_present ) missvalue = input_fields(diag_field_id)%missing_value

    number_of_outputs = input_fields(diag_field_id)%num_output_fields
!$OMP CRITICAL
    input_fields(diag_field_id)%numthreads = 1
    active_omp_level=0
#if defined(_OPENMP)
    input_fields(diag_field_id)%numthreads = omp_get_num_threads()
    input_fields(diag_field_id)%active_omp_level = omp_get_level()
#endif
    numthreads = input_fields(diag_field_id)%numthreads
    active_omp_level = input_fields(diag_field_id)%active_omp_level
!$OMP END CRITICAL

    if(present(time)) input_fields(diag_field_id)%time = time

    ! Issue a warning if any value in field is outside the valid range
    IF ( input_fields(diag_field_id)%range_present ) THEN
       IF ( ISSUE_OOR_WARNINGS .OR. OOR_WARNINGS_FATAL ) THEN
          WRITE (error_string, '("[",ES14.5E3,",",ES14.5E3,"]")')&
               & input_fields(diag_field_id)%range(1:2)
          WRITE (error_string1, '("(Min: ",ES14.5E3,", Max: ",ES14.5E3, ")")')&
                  & MINVAL(field(f1:f2,f3:f4,ks:ke),MASK=oor_mask(f1:f2,f3:f4,ks:ke)),&
                  & MAXVAL(field(f1:f2,f3:f4,ks:ke),MASK=oor_mask(f1:f2,f3:f4,ks:ke))
          IF ( missvalue_present ) THEN
             IF ( ANY(oor_mask(f1:f2,f3:f4,ks:ke) .AND.&
                  &   ((field(f1:f2,f3:f4,ks:ke) < input_fields(diag_field_id)%range(1) .OR.&
                  &     field(f1:f2,f3:f4,ks:ke) > input_fields(diag_field_id)%range(2)).AND.&
                  &     field(f1:f2,f3:f4,ks:ke) .NE. missvalue)) ) THEN
                ! <ERROR STATUS="WARNING/FATAL">
                !   A value for <module_name> in field <field_name> (Min: <min_val>, Max: <max_val>)
                !   is outside the range [<lower_val>,<upper_val>] and not equal to the missing
                !   value.
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::send_data_3d',&
                     & 'A value for '//&
                     &TRIM(input_fields(diag_field_id)%module_name)//' in field '//&
                     &TRIM(input_fields(diag_field_id)%field_name)//' '&
                     &//TRIM(error_string1)//&
                     &' is outside the range '//TRIM(error_string)//',&
                     & and not equal to the missing value.',&
                     &OOR_WARNING)
             END IF
          ELSE
             IF ( ANY(oor_mask(f1:f2,f3:f4,ks:ke) .AND.&
                  &   (field(f1:f2,f3:f4,ks:ke) < input_fields(diag_field_id)%range(1) .OR.&
                  &    field(f1:f2,f3:f4,ks:ke) > input_fields(diag_field_id)%range(2))) ) THEN
                ! <ERROR STATUS="WARNING/FATAL">
                !   A value for <module_name> in field <field_name> (Min: <min_val>, Max: <max_val>)
                !   is outside the range [<lower_val>,<upper_val>].
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::send_data_3d',&
                     & 'A value for '//&
                     &TRIM(input_fields(diag_field_id)%module_name)//' in field '//&
                     &TRIM(input_fields(diag_field_id)%field_name)//' '&
                     &//TRIM(error_string1)//&
                     &' is outside the range '//TRIM(error_string)//'.',&
                     &OOR_WARNING)
             END IF
          END IF
       END IF
    END IF

    ! Loop through each output field that depends on this input field
    num_out_fields: DO ii = 1, number_of_outputs
       ! Get index to an output field
       out_num = input_fields(diag_field_id)%output_fields(ii)

      ! is this field output on a local domain only?
       local_output = output_fields(out_num)%local_output
       ! if local_output, does the current PE take part in send_data?
       need_compute = output_fields(out_num)%need_compute

       reduced_k_range = output_fields(out_num)%reduced_k_range

      ! skip all PEs not participating in outputting this field
       IF ( local_output .AND. (.NOT.need_compute) ) CYCLE

       ! Get index to output file for this field
       file_num = output_fields(out_num)%output_file
       IF(file_num == max_files) CYCLE
       ! Output frequency and units for this file is
       freq = files(file_num)%output_freq
       units = files(file_num)%output_units
       ! Is this output field being time averaged?
       average = output_fields(out_num)%time_average
       ! Is this output field the rms?
       ! If so, then average is also .TRUE.
       time_rms = output_fields(out_num)%time_rms
       ! Power value for rms or pow(x) calculations
       pow_value = output_fields(out_num)%pow_value
       ! Looking for max and min value of this field over the sampling interval?
       time_max = output_fields(out_num)%time_max
       time_min = output_fields(out_num)%time_min
       ! Sum output over time interval
       time_sum = output_fields(out_num)%time_sum
       IF ( output_fields(out_num)%total_elements > SIZE(field(f1:f2,f3:f4,ks:ke)) ) THEN
          output_fields(out_num)%phys_window = .TRUE.
       ELSE
          output_fields(out_num)%phys_window = .FALSE.
       END IF
       phys_window = output_fields(out_num)%phys_window
       IF ( need_compute ) THEN
          l_start = output_fields(out_num)%output_grid%l_start_indx
          l_end = output_fields(out_num)%output_grid%l_end_indx
       END IF

       ! compute the diurnal index
       sample = 1
       IF ( PRESENT(time) ) THEN
          CALL get_time(time,second,day,tick) ! current date
          sample = floor((second+real(tick)/get_ticks_per_second())*output_fields(out_num)%n_diurnal_samples/SECONDS_PER_DAY) + 1
       END IF

       ! Get the vertical layer start and end index.
       IF ( reduced_k_range ) THEN
!----------
!ug support
           if (output_fields(out_num)%reduced_k_unstruct) then
               js = output_fields(out_num)%output_grid%l_start_indx(2)
               je = output_fields(out_num)%output_grid%l_end_indx(2)
           endif
           l_start(3) = output_fields(out_num)%output_grid%l_start_indx(3)
           l_end(3) = output_fields(out_num)%output_grid%l_end_indx(3)
!----------
       END IF
       ksr= l_start(3)
       ker= l_end(3)

       ! Initialize output time for fields output every time step
       IF ( freq == EVERY_TIME .AND. .NOT.output_fields(out_num)%static ) THEN
          IF (output_fields(out_num)%next_output == output_fields(out_num)%last_output) THEN
             IF(PRESENT(time)) THEN
                output_fields(out_num)%next_output = time
             ELSE
                WRITE (error_string,'(a,"/",a)')&
                     & TRIM(input_fields(diag_field_id)%module_name),&
                     & TRIM(output_fields(out_num)%output_name)
                IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                     & ', time must be present when output frequency = EVERY_TIME', err_msg)) THEN
                   DEALLOCATE(oor_mask)
                   RETURN
                END IF
             END IF
          END IF
       END IF
       IF ( .NOT.output_fields(out_num)%static .AND. .NOT.PRESENT(time) ) THEN
          WRITE (error_string,'(a,"/",a)')&
               & TRIM(input_fields(diag_field_id)%module_name), &
               & TRIM(output_fields(out_num)%output_name)
          IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
               & ', time must be present for nonstatic field', err_msg)) THEN
             DEALLOCATE(oor_mask)
             RETURN
          END IF
       END IF

       ! Is it time to output for this field; CAREFUL ABOUT > vs >= HERE
       !--- The fields send out within openmp parallel region will be written out in
       !--- diag_send_complete.
       IF ( (numthreads == 1) .AND. (active_omp_level.LE.1) ) then
          IF ( .NOT.output_fields(out_num)%static .AND. freq /= END_OF_RUN ) THEN
             IF ( time > output_fields(out_num)%next_output ) THEN
                ! A non-static field that has skipped a time level is an error
                IF ( time > output_fields(out_num)%next_next_output .AND. freq > 0 ) THEN
                   IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
                      WRITE (error_string,'(a,"/",a)')&
                           & TRIM(input_fields(diag_field_id)%module_name), &
                           & TRIM(output_fields(out_num)%output_name)
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                           & ' is skipped one time level in output data', err_msg)) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF

                status = writing_field(out_num, .FALSE., error_string, time)
                IF(status == -1) THEN
                   IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
                      IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
                           & ', write EMPTY buffer', err_msg)) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
             END IF  !time > output_fields(out_num)%next_output
          END IF  !.not.output_fields(out_num)%static .and. freq /= END_OF_RUN
          ! Finished output of previously buffered data, now deal with buffering new data
       END IF

       IF ( .NOT.output_fields(out_num)%static .AND. .NOT.need_compute .AND. debug_diag_manager ) THEN
          CALL check_bounds_are_exact_dynamic(out_num, diag_field_id, Time, err_msg=err_msg_local)
          IF ( err_msg_local /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                DEALLOCATE(oor_mask)
                RETURN
             END IF
          END IF
       END IF

       ! Take care of submitted field data
       IF ( average ) THEN
          IF ( input_fields(diag_field_id)%mask_variant ) THEN
             IF ( need_compute ) THEN
                WRITE (error_string,'(a,"/",a)')  &
                     & TRIM(input_fields(diag_field_id)%module_name), &
                     & TRIM(output_fields(out_num)%output_name)
                IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                     & ', regional output NOT supported with mask_variant', err_msg)) THEN
                   DEALLOCATE(oor_mask)
                   RETURN
                END IF
             END IF

             ! Should reduced_k_range data be supported with the mask_variant option   ?????
             ! If not, error message should be produced and the reduced_k_range loop below eliminated
             IF ( PRESENT(mask) ) THEN
                IF ( missvalue_present ) THEN
                   IF ( debug_diag_manager ) THEN
                      CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                      CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                      IF ( err_msg_local /= '' ) THEN
                         IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                            DEALLOCATE(oor_mask)
                            RETURN
                         END IF
                      END IF
                   END IF
                   IF( numthreads>1 .AND. phys_window ) then
                      IF ( reduced_k_range ) THEN
                         DO k= ksr, ker
                            k1= k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & (field(i-is+1+hi, j-js+1+hj, k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & field(i-is+1+hi, j-js+1+hj, k) * weight1
                                     END IF
                                     output_fields(out_num)%counter(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%counter(i-hi,j-hj,k1,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k)*weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k)*weight1
                                     END IF
                                     output_fields(out_num)%counter(i-hi,j-hj,k,sample) =&
                                          &output_fields(out_num)%counter(i-hi,j-hj,k,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      END IF
                   ELSE
!$OMP CRITICAL
                      IF ( reduced_k_range ) THEN
                         DO k= ksr, ker
                            k1= k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & (field(i-is+1+hi, j-js+1+hj, k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & field(i-is+1+hi, j-js+1+hj, k) * weight1
                                     END IF
                                     output_fields(out_num)%counter(i-hi,j-hj,k1,sample) =&
                                          & output_fields(out_num)%counter(i-hi,j-hj,k1,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k)*weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k)*weight1
                                     END IF
                                     output_fields(out_num)%counter(i-hi,j-hj,k,sample) =&
                                          &output_fields(out_num)%counter(i-hi,j-hj,k,sample) + weight1
                                  END IF
                               END DO
                            END DO
                         END DO
                      END IF
!$OMP END CRITICAL
                   END IF
                ELSE
                   WRITE (error_string,'(a,"/",a)')&
                        & TRIM(input_fields(diag_field_id)%module_name), &
                        & TRIM(output_fields(out_num)%output_name)
                   IF ( fms_error_handler('diag_manager_mod::send_data_3d', 'module/output_field '//TRIM(error_string)//&
                        & ', variable mask but no missing value defined', err_msg)) THEN
                      DEALLOCATE(oor_mask)
                      RETURN
                   END IF
                END IF
             ELSE  ! no mask present
                WRITE (error_string,'(a,"/",a)')&
                     & TRIM(input_fields(diag_field_id)%module_name), &
                     & TRIM(output_fields(out_num)%output_name)
                IF(fms_error_handler('diag_manager_mod::send_data_3d','module/output_field '//TRIM(error_string)//&
                     & ', variable mask but no mask given', err_msg)) THEN
                   DEALLOCATE(oor_mask)
                   RETURN
                END IF
             END IF
          ELSE ! mask_variant=false
             IF ( PRESENT(mask) ) THEN
                IF ( missvalue_present ) THEN
                   IF ( need_compute ) THEN
                      IF (numthreads>1 .AND. phys_window) then
                         DO k = l_start(3), l_end(3)
                            k1 = k-l_start(3)+1
                            DO j = js, je
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                     i1 = i-l_start(1)-hi+1
                                     j1=  j-l_start(2)-hj+1
                                     IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                        IF ( pow_value /= 1 ) THEN
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                        ELSE
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                        END IF
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k = l_start(3), l_end(3)
                            k1 = k-l_start(3)+1
                            DO j = js, je
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                     i1 = i-l_start(1)-hi+1
                                     j1=  j-l_start(2)-hj+1
                                     IF ( mask(i-is+1+hi, j-js+1+hj, k) ) THEN
                                        IF ( pow_value /= 1 ) THEN
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                        ELSE
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                        END IF
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      ENDIF
!$OMP CRITICAL
                      DO j = js, je
                         DO i = is, ie
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                               output_fields(out_num)%num_elements(sample) = &
                                    output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                            END IF
                         END DO
                      END DO
!$OMP END CRITICAL
                   ELSE IF ( reduced_k_range ) THEN
                      IF (numthreads>1 .AND. phys_window) then
                         DO k=ksr, ker
                            k1 = k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k=ksr, ker
                            k1 = k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
                   ELSE
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '' ) THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF (numthreads>1 .AND. phys_window) then
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
                   END IF
!$OMP CRITICAL
                   IF ( need_compute .AND. .NOT.phys_window ) THEN
                      IF ( ANY(mask(l_start(1)+hi:l_end(1)+hi,l_start(2)+hj:l_end(2)+hj,l_start(3):l_end(3))) ) &
                           & output_fields(out_num)%count_0d(sample) =&
                           & output_fields(out_num)%count_0d(sample) + weight1
                   ELSE
                      IF ( ANY(mask(f1:f2,f3:f4,ks:ke)) ) output_fields(out_num)%count_0d(sample) =&
                           & output_fields(out_num)%count_0d(sample)+weight1
                   END IF
!$OMP END CRITICAL

                ELSE ! missing value NOT present
                   IF (   (.NOT.ALL(mask(f1:f2,f3:f4,ks:ke)) .AND. mpp_pe() .EQ. mpp_root_pe()).AND.&
                        &  .NOT.input_fields(diag_field_id)%issued_mask_ignore_warning ) THEN
                      ! <ERROR STATUS="WARNING">
                      !   Mask will be ignored since missing values were not specified for field <field_name>
                      !   in module <module_name>
                      ! </ERROR>
                      CALL error_mesg('diag_manager_mod::send_data_3d',&
                           & 'Mask will be ignored since missing values were not specified for field '//&
                           & trim(input_fields(diag_field_id)%field_name)//' in module '//&
                           & trim(input_fields(diag_field_id)%module_name), WARNING)
                      input_fields(diag_field_id)%issued_mask_ignore_warning = .TRUE.
                   END IF
                   IF ( need_compute ) THEN
                      IF (numthreads>1 .AND. phys_window) then
                         DO j = js, je
                            DO i = is, ie
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1
                                  j1 =  j-l_start(2)-hj+1
                                  IF ( pow_value /= 1 ) THEN
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                          & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                                  ELSE
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                          & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                               END IF
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO j = js, je
                            DO i = is, ie
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1
                                  j1 =  j-l_start(2)-hj+1
                                  IF ( pow_value /= 1 ) THEN
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                          & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                                  ELSE
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample)+ &
                                          & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                               END IF
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
!$OMP CRITICAL
                      DO j = js, je
                         DO i = is, ie
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                               output_fields(out_num)%num_elements(sample)=&
                                    & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1

                            END IF
                         END DO
                      END DO
!$OMP END CRITICAL
                   ELSE IF ( reduced_k_range ) THEN
                      IF (numthreads>1 .AND. phys_window) then
                         ksr= l_start(3)
                         ker= l_end(3)
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                                 & (field(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                                 & field(f1:f2,f3:f4,ksr:ker)*weight1
                         END IF
                      ELSE
!$OMP CRITICAL
                         ksr= l_start(3)
                         ker= l_end(3)
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                                 & (field(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) +&
                                 & field(f1:f2,f3:f4,ksr:ker)*weight1
                         END IF
!$OMP END CRITICAL
                      END IF
                   ELSE
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '') THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF (numthreads>1 .AND. phys_window) then
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & (field(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & field(f1:f2,f3:f4,ks:ke)*weight1
                         END IF
                      ELSE
!$OMP CRITICAL
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & (field(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & field(f1:f2,f3:f4,ks:ke)*weight1
                         END IF
!$OMP END CRITICAL
                      END IF
                   END IF
!$OMP CRITICAL
                   IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
                        & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL
                END IF
             ELSE ! mask NOT present
                IF ( missvalue_present ) THEN
                   IF ( need_compute ) THEN
                      if( numthreads>1 .AND. phys_window ) then
                         DO k = l_start(3), l_end(3)
                            k1 = k - l_start(3) + 1
                            DO j = js, je
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj) THEN
                                     i1 = i-l_start(1)-hi+1
                                     j1=  j-l_start(2)-hj+1
                                     IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                        IF ( pow_value /= 1 ) THEN
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                        ELSE
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                        END IF
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k = l_start(3), l_end(3)
                            k1 = k - l_start(3) + 1
                            DO j = js, je
                               DO i = is, ie
                                  IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj) THEN
                                     i1 = i-l_start(1)-hi+1
                                     j1=  j-l_start(2)-hj+1
                                     IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                        IF ( pow_value /= 1 ) THEN
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                        ELSE
                                           output_fields(out_num)%buffer(i1,j1,k1,sample) =&
                                                & output_fields(out_num)%buffer(i1,j1,k1,sample) +&
                                                & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                        END IF
                                     ELSE
                                        output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                                     END IF
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
!$OMP CRITICAL
                      DO j = js, je
                         DO i = is, ie
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj) THEN
                               output_fields(out_num)%num_elements(sample) =&
                                    & output_fields(out_num)%num_elements(sample) + l_end(3) - l_start(3) + 1
                            END IF
                         END DO
                      END DO
                      IF ( .NOT.phys_window ) THEN
                         outer0: DO k = l_start(3), l_end(3)
                            DO j=l_start(2)+hj, l_end(2)+hj
                               DO i=l_start(1)+hi, l_end(1)+hi
                                  IF ( field(i,j,k) /= missvalue ) THEN
                                     output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) + weight1
                                     EXIT outer0
                                  END IF
                               END DO
                            END DO
                         END DO outer0
                      END IF
!$OMP END CRITICAL
                   ELSE IF ( reduced_k_range ) THEN
                      if( numthreads>1 .AND. phys_window ) then
                         ksr= l_start(3)
                         ker= l_end(3)
                         DO k = ksr, ker
                            k1 = k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      else
!$OMP CRITICAL
                         ksr= l_start(3)
                         ker= l_end(3)
                         DO k = ksr, ker
                            k1 = k - ksr + 1
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue ) THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k1,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
!$OMP CRITICAL
                      outer3: DO k = ksr, ker
                         k1=k-ksr+1
                         DO j=f3, f4
                            DO i=f1, f2
                               IF ( field(i,j,k) /= missvalue ) THEN
                                  output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) + weight1
                                  EXIT outer3
                               END IF
                            END DO
                         END DO
                      END DO outer3
!$OMP END CRITICAL
                   ELSE
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '' ) THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF( numthreads > 1 .AND. phys_window ) then
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO k=ks, ke
                            DO j=js, je
                               DO i=is, ie
                                  IF ( field(i-is+1+hi,j-js+1+hj,k) /= missvalue )  THEN
                                     IF ( pow_value /= 1 ) THEN
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & (field(i-is+1+hi,j-js+1+hj,k) * weight1)**(pow_value)
                                     ELSE
                                        output_fields(out_num)%buffer(i-hi,j-hj,k,sample) =&
                                             & output_fields(out_num)%buffer(i-hi,j-hj,k,sample) +&
                                             & field(i-is+1+hi,j-js+1+hj,k) * weight1
                                     END IF
                                  ELSE
                                     output_fields(out_num)%buffer(i-hi,j-hj,k,sample) = missvalue
                                  END IF
                               END DO
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF
!$OMP CRITICAL
                      outer1: DO k=ks, ke
                         DO j=f3, f4
                            DO i=f1, f2
                               IF ( field(i,j,k) /= missvalue ) THEN
                                  output_fields(out_num)%count_0d(sample) = output_fields(out_num)%count_0d(sample) + weight1
                                  EXIT outer1
                               END IF
                            END DO
                         END DO
                      END DO outer1
!$OMP END CRITICAL
                   END IF
                ELSE ! no missing value defined, No mask
                   IF ( need_compute ) THEN
                      IF( numthreads > 1 .AND. phys_window ) then
                         DO j = js, je
                            DO i = is, ie
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1
                                  j1=  j-l_start(2)-hj+1
                                  IF ( pow_value /= 1 ) THEN
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                          & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                                  ELSE
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                               END IF
                            END DO
                         END DO
                      ELSE
!$OMP CRITICAL
                         DO j = js, je
                            DO i = is, ie
                               IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                                  i1 = i-l_start(1)-hi+1
                                  j1=  j-l_start(2)-hj+1
                                  IF ( pow_value /= 1 ) THEN
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                          & (field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1)**(pow_value)
                                  ELSE
                                     output_fields(out_num)%buffer(i1,j1,:,sample)= output_fields(out_num)%buffer(i1,j1,:,sample) +&
                                          & field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))*weight1
                               END IF
                               END IF
                            END DO
                         END DO
!$OMP END CRITICAL
                      END IF

!$OMP CRITICAL
                      DO j = js, je
                         DO i = is, ie
                            IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                               output_fields(out_num)%num_elements(sample) =&
                                    & output_fields(out_num)%num_elements(sample)+l_end(3)-l_start(3)+1
                            END IF
                         END DO
                      END DO
!$OMP END CRITICAL
                      ! Accumulate time average
                   ELSE IF ( reduced_k_range ) THEN
                      ksr= l_start(3)
                      ker= l_end(3)
                      IF( numthreads > 1 .AND. phys_window ) then
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                                 & (field(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                                 & field(f1:f2,f3:f4,ksr:ker)*weight1
                         END IF
                      ELSE
!$OMP CRITICAL
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                                 & (field(f1:f2,f3:f4,ksr:ker)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                                 & field(f1:f2,f3:f4,ksr:ker)*weight1
                         END IF
!$OMP END CRITICAL
                      END IF
                   ELSE
                      IF ( debug_diag_manager ) THEN
                         CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                         CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                         IF ( err_msg_local /= '' ) THEN
                            IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                               DEALLOCATE(oor_mask)
                               RETURN
                            END IF
                         END IF
                      END IF
                      IF( numthreads > 1 .AND. phys_window ) then
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & (field(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & field(f1:f2,f3:f4,ks:ke)*weight1
                         END IF
                      ELSE
!$OMP CRITICAL
                         IF ( pow_value /= 1 ) THEN
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & (field(f1:f2,f3:f4,ks:ke)*weight1)**(pow_value)
                         ELSE
                            output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) =&
                                 & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) +&
                                 & field(f1:f2,f3:f4,ks:ke)*weight1
                         END IF
!$OMP END CRITICAL
                      END IF
                   END IF
!$OMP CRITICAL
                   IF ( .NOT.phys_window ) output_fields(out_num)%count_0d(sample) =&
                        & output_fields(out_num)%count_0d(sample) + weight1
!$OMP END CRITICAL
                END IF
             END IF ! if mask present
          END IF  !if mask_variant
!$OMP CRITICAL
          IF ( .NOT.need_compute .AND. .NOT.reduced_k_range )&
               & output_fields(out_num)%num_elements(sample) =&
               & output_fields(out_num)%num_elements(sample) + (ie-is+1)*(je-js+1)*(ke-ks+1)
          IF ( reduced_k_range ) &
               & output_fields(out_num)%num_elements(sample) = output_fields(out_num)%num_elements(sample) +&
               & (ie-is+1)*(je-js+1)*(ker-ksr+1)
!$OMP END CRITICAL
             ! Add processing for Max and Min
       ELSE IF ( time_max ) THEN
          IF ( PRESENT(mask) ) THEN
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1
                            j1=  j-l_start(2)-hj+1
                            IF ( mask(i-is+1+hi,j-js+1+hj,k) .AND.&
                                 & field(i-is+1+hi,j-js+1+hj,k)>output_fields(out_num)%buffer(i1,j1,k1,sample)) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Maximum time value with masking
             ELSE IF ( reduced_k_range ) THEN
                ksr = l_start(3)
                ker = l_end(3)
                WHERE ( mask(f1:f2,f3:f4,ksr:ker) .AND. &
                     & field(f1:f2,f3:f4,ksr:ker) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample))&
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( mask(f1:f2,f3:f4,ks:ke) .AND.&
                     & field(f1:f2,f3:f4,ks:ke)>output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample))&
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
             END IF
          ELSE
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF(l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1
                            j1 =  j-l_start(2)-hj+1
                            IF ( field(i-is+1+hi,j-js+1+hj,k) > output_fields(out_num)%buffer(i1,j1,k1,sample) ) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Maximum time value
             ELSE IF ( reduced_k_range ) THEN
                ksr = l_start(3)
                ker = l_end(3)
                WHERE ( field(f1:f2,f3:f4,ksr:ker) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) ) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( field(f1:f2,f3:f4,ks:ke) > output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) ) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
             END IF
          END IF
          output_fields(out_num)%count_0d(sample) = 1
       ELSE IF ( time_min ) THEN
          IF ( PRESENT(mask) ) THEN
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1
                            j1 =  j-l_start(2)-hj+1
                            IF ( mask(i-is+1+hi,j-js+1+hj,k) .AND.&
                                 & field(i-is+1+hi,j-js+1+hj,k) < output_fields(out_num)%buffer(i1,j1,k1,sample) ) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Minimum time value with masking
             ELSE IF ( reduced_k_range ) THEN
                ksr= l_start(3)
                ker= l_end(3)
                WHERE ( mask(f1:f2,f3:f4,ksr:ker) .AND.&
                     & field(f1:f2,f3:f4,ksr:ker) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample)) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( mask(f1:f2,f3:f4,ks:ke) .AND.&
                     & field(f1:f2,f3:f4,ks:ke) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) ) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
             END IF
          ELSE
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                            i1 = i-l_start(1)-hi+1
                            j1=  j-l_start(2)-hj+1
                            IF ( field(i-is+1+hi,j-js+1+hj,k) < output_fields(out_num)%buffer(i1,j1,k1,sample) ) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Minimum time value
             ELSE IF ( reduced_k_range ) THEN
                ksr= l_start(3)
                ker= l_end(3)
                WHERE ( field(f1:f2,f3:f4,ksr:ker) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) ) &
                     output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( field(f1:f2,f3:f4,ks:ke) < output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) )&
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
             END IF
          END IF
          output_fields(out_num)%count_0d(sample) = 1
       ELSE IF ( time_sum ) THEN
          IF ( PRESENT(mask) ) THEN
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1
                            j1 =  j-l_start(2)-hj+1
                            IF ( mask(i-is+1+hi,j-js+1+hj,k) ) THEN
                               output_fields(out_num)%buffer(i1,j1,k1,sample) = &
                                    output_fields(out_num)%buffer(i1,j1,k1,sample) + &
                                    field(i-is+1+hi,j-js+1+hj,k)
                            END IF
                         END IF
                      END DO
                   END DO
                END DO
                ! Minimum time value with masking
             ELSE IF ( reduced_k_range ) THEN
                ksr= l_start(3)
                ker= l_end(3)
                output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = &
                     &   output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                     &   field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                WHERE ( mask(f1:f2,f3:f4,ks:ke) ) &
                     & output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = &
                     &  output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) + &
                     &  field(f1:f2,f3:f4,ks:ke)
             END IF
          ELSE
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <=i.AND.i<=l_end(1)+hi.AND.l_start(2)+hj<=j.AND.j<=l_end(2)+hj) THEN
                            i1 = i-l_start(1)-hi+1
                            j1=  j-l_start(2)-hj+1
                            output_fields(out_num)%buffer(i1,j1,k1,sample) = &
                               &    output_fields(out_num)%buffer(i1,j1,k1,sample) + &
                               &    field(i-is+1+hi,j-js+1+hj,k)
                         END IF
                      END DO
                   END DO
                END DO
             ELSE IF ( reduced_k_range ) THEN
                ksr= l_start(3)
                ker= l_end(3)
                output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = &
                     &  output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) + &
                     &  field(f1:f2,f3:f4,ksr:ker)
             ELSE
                IF ( debug_diag_manager ) THEN
                   CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                   CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                   IF ( err_msg_local /= '' ) THEN
                      IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                         DEALLOCATE(oor_mask)
                         RETURN
                      END IF
                   END IF
                END IF
                output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = &
                &    output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) + &
                &    field(f1:f2,f3:f4,ks:ke)
             END IF
          END IF
          output_fields(out_num)%count_0d(sample) = 1
       ELSE  ! ( not average, not min, not max, not sum )
          output_fields(out_num)%count_0d(sample) = 1
          IF ( need_compute ) THEN
             DO j = js, je
                DO i = is, ie
                   IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                      i1 = i-l_start(1)-hi+1
                      j1 = j-l_start(2)-hj+1
                      output_fields(out_num)%buffer(i1,j1,:,sample) = field(i-is+1+hi,j-js+1+hj,l_start(3):l_end(3))
                   END IF
                END DO
             END DO
             ! instantaneous output
          ELSE IF ( reduced_k_range ) THEN
             ksr = l_start(3)
             ker = l_end(3)
             output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,:,sample) = field(f1:f2,f3:f4,ksr:ker)
          ELSE
             IF ( debug_diag_manager ) THEN
                CALL update_bounds(out_num, is-hi, ie-hi, js-hj, je-hj, ks, ke)
                CALL check_out_of_bounds(out_num, diag_field_id, err_msg=err_msg_local)
                IF ( err_msg_local /= '' ) THEN
                   IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg) ) THEN
                      DEALLOCATE(oor_mask)
                      RETURN
                   END IF
                END IF
             END IF
             output_fields(out_num)%buffer(is-hi:ie-hi,js-hj:je-hj,ks:ke,sample) = field(f1:f2,f3:f4,ks:ke)
          END IF

          IF ( PRESENT(mask) .AND. missvalue_present ) THEN
             IF ( need_compute ) THEN
                DO k = l_start(3), l_end(3)
                   k1 = k - l_start(3) + 1
                   DO j = js, je
                      DO i = is, ie
                         IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                            i1 = i-l_start(1)-hi+1
                            j1 =  j-l_start(2)-hj+1
                            IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) )&
                                 & output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                         END IF
                      END DO
                   END DO
                END DO
             ELSE IF ( reduced_k_range ) THEN
                ksr= l_start(3)
                ker= l_end(3)
                DO k=ksr, ker
                   k1= k - ksr + 1
                   DO j=js, je
                      DO i=is, ie
                         IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) ) &
                              & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                      END DO
                   END DO
                END DO
             ELSE
                DO k=ks, ke
                   DO j=js, je
                      DO i=is, ie
                         IF ( .NOT.mask(i-is+1+hi,j-js+1+hj,k) )&
                              & output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                      END DO
                   END DO
                END DO
             END IF
          END IF
       END IF !average

       IF ( output_fields(out_num)%static .AND. .NOT.need_compute .AND. debug_diag_manager ) THEN
          CALL check_bounds_are_exact_static(out_num, diag_field_id, err_msg=err_msg_local)
          IF ( err_msg_local /= '' ) THEN
             IF ( fms_error_handler('diag_manager_mod::send_data_3d', err_msg_local, err_msg)) THEN
                DEALLOCATE(oor_mask)
                RETURN
             END IF
          END IF
       END IF

       ! If rmask and missing value present, then insert missing value
       IF ( PRESENT(rmask) .AND. missvalue_present ) THEN
          IF ( need_compute ) THEN
             DO k = l_start(3), l_end(3)
                k1 = k - l_start(3) + 1
                DO j = js, je
                   DO i = is, ie
                      IF ( l_start(1)+hi <= i .AND. i <= l_end(1)+hi .AND. l_start(2)+hj <= j .AND. j <= l_end(2)+hj ) THEN
                         i1 = i-l_start(1)-hi+1
                         j1 =  j-l_start(2)-hj+1
                         IF ( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5 ) &
                              & output_fields(out_num)%buffer(i1,j1,k1,sample) = missvalue
                      END IF
                   END DO
                END DO
             END DO
          ELSE IF ( reduced_k_range ) THEN
             ksr= l_start(3)
             ker= l_end(3)
             DO k= ksr, ker
                k1 = k - ksr + 1
                DO j=js, je
                   DO i=is, ie
                      IF ( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5 ) &
                           & output_fields(out_num)%buffer(i-hi,j-hj,k1,sample)= missvalue
                   END DO
                END DO
             END DO
          ELSE
             DO k=ks, ke
                DO j=js, je
                   DO i=is, ie
                      IF ( rmask(i-is+1+hi,j-js+1+hj,k) < 0.5 ) &
                           & output_fields(out_num)%buffer(i-hi,j-hj,k,sample)= missvalue
                   END DO
                END DO
             END DO
          END IF
       END IF

    END DO num_out_fields

    DEALLOCATE(oor_mask)
  END FUNCTION send_data_3d

  !> @return true if send is successful
  LOGICAL FUNCTION send_tile_averaged_data1d ( id, field, area, time, mask )
    INTEGER, INTENT(in) :: id  !< id od the diagnostic field
    REAL, INTENT(in) :: field(:,:) !< field to average and send
    REAL, INTENT(in) :: area (:,:) !< area of tiles (== averaging weights), arbitrary units
    TYPE(time_type), INTENT(in)  :: time !< current time
    LOGICAL, INTENT(in),OPTIONAL :: mask (:,:) !< land mask

    REAL, DIMENSION(SIZE(field,1)) :: out(SIZE(field,1))

    ! If id is < 0 it means that this field is not registered, simply return
    IF ( id <= 0 ) THEN
       send_tile_averaged_data1d = .FALSE.
       RETURN
    END IF

    CALL average_tiles1d (id, field, area, mask, out)
    send_tile_averaged_data1d = send_data(id, out, time=time, mask=ANY(mask,DIM=2))
  END FUNCTION send_tile_averaged_data1d

  !> @brief Calculates average for a field with the given area and land mask
  SUBROUTINE average_tiles1d(diag_field_id, x, area, mask, out)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:,:), INTENT(in) :: x !< (ug_index, tile) field to average
    REAL, DIMENSION(:,:), INTENT(in) :: area !< (ug_index, tile) fractional area
    LOGICAL, DIMENSION(:,:), INTENT(in) :: mask !< (ug_index, tile) land mask
    REAL, DIMENSION(:), INTENT(out) :: out !< (ug_index) result of averaging

    INTEGER  :: it !< iterator over tile number
    REAL, DIMENSION(SIZE(x,1)) :: s !< area accumulator
    REAL :: local_missing_value

    ! # FATAL if diag_field_id is less than 0, indicates field was not in diag_table.
    ! The calling functions should not have passed in an invalid diag_field_id
    IF ( diag_field_id <= 0 ) THEN
       ! <ERROR STATUS="FATAL">
       !   diag_field_id less than 0.  Contact developers.
       ! </ERROR>
       CALL error_mesg('diag_manager_mod::average_tiles1d',&
            & "diag_field_id less than 0.  Contact developers.", FATAL)
    END IF

    ! Initialize local_missing_value
    IF ( input_fields(diag_field_id)%missing_value_present ) THEN
       local_missing_value = input_fields(diag_field_id)%missing_value
    ELSE
       local_missing_value = 0.0
    END IF

    ! Initialize s and out to zero.
    s(:) = 0.0
    out(:) = 0.0

    DO it = 1, SIZE(area,dim=2)
       WHERE ( mask(:,it) )
          out(:) = out(:) + x(:,it)*area(:,it)
          s(:) = s(:) + area(:,it)
       END WHERE
    END DO

    WHERE ( s(:) > 0 )
       out(:) = out(:)/s(:)
    ELSEWHERE
       out(:) = local_missing_value
    END WHERE
  END SUBROUTINE average_tiles1d

  !> @return true if send is successful
  LOGICAL FUNCTION send_tile_averaged_data2d ( id, field, area, time, mask )
    INTEGER, INTENT(in) :: id  !< id od the diagnostic field
    REAL, INTENT(in) :: field(:,:,:) !< field to average and send
    REAL, INTENT(in) :: area (:,:,:) !< area of tiles (== averaging weights), arbitrary units
    TYPE(time_type), INTENT(in)  :: time !< current time
    LOGICAL, INTENT(in),OPTIONAL :: mask (:,:,:) !< land mask

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2)) :: out(SIZE(field,1), SIZE(field,2))

    ! If id is < 0 it means that this field is not registered, simply return
    IF ( id <= 0 ) THEN
       send_tile_averaged_data2d = .FALSE.
       RETURN
    END IF

    CALL average_tiles(id, field, area, mask, out)
    send_tile_averaged_data2d = send_data(id, out, time, mask=ANY(mask,DIM=3))
  END FUNCTION send_tile_averaged_data2d

  !> @return true if send is successful
  LOGICAL FUNCTION send_tile_averaged_data3d( id, field, area, time, mask )
    INTEGER, INTENT(in) :: id !< id of the diagnostic field
    REAL, DIMENSION(:,:,:,:), INTENT(in) :: field !< (lon, lat, tile, lev) field to average and send
    REAL, DIMENSION(:,:,:), INTENT(in) :: area (:,:,:) !< (lon, lat, tile) tile areas ( == averaging weights), arbitrary units
    TYPE(time_type), INTENT(in)  :: time !< current time
    LOGICAL, DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask !< (lon, lat, tile) land mask

    REAL, DIMENSION(SIZE(field,1),SIZE(field,2),SIZE(field,4)) :: out
    LOGICAL, DIMENSION(SIZE(field,1),SIZE(field,2),SIZE(field,4)) :: mask3
    INTEGER :: it

    ! If id is < 0 it means that this field is not registered, simply return
    IF ( id <= 0 ) THEN
       send_tile_averaged_data3d = .FALSE.
       RETURN
    END IF

    DO it=1, SIZE(field,4)
       CALL average_tiles(id, field(:,:,:,it), area, mask, out(:,:,it) )
    END DO

    mask3(:,:,1) = ANY(mask,DIM=3)
    DO it = 2, SIZE(field,4)
       mask3(:,:,it) = mask3(:,:,1)
    END DO

    send_tile_averaged_data3d = send_data( id, out, time, mask=mask3 )
  END FUNCTION send_tile_averaged_data3d

  !> @brief Calculates tile average of a field
  SUBROUTINE average_tiles(diag_field_id, x, area, mask, out)
    INTEGER, INTENT(in) :: diag_field_id
    REAL, DIMENSION(:,:,:), INTENT(in) :: x !< (lon, lat, tile) field to average
    REAL, DIMENSION(:,:,:), INTENT(in) :: area !< (lon, lat, tile) fractional area
    LOGICAL, DIMENSION(:,:,:), INTENT(in) :: mask !< (lon, lat, tile) land mask
    REAL, DIMENSION(:,:), INTENT(out) :: out !< (lon, lat) result of averaging

    INTEGER  :: it !< iterator over tile number
    REAL, DIMENSION(SIZE(x,1),SIZE(x,2)) :: s !< area accumulator
    REAL :: local_missing_value

    ! # FATAL if diag_field_id is less than 0, indicates field was not in diag_table.
    ! The calling functions should not have passed in an invalid diag_field_id
    IF ( diag_field_id <= 0 ) THEN
       ! <ERROR STATUS="FATAL">
       !   diag_field_id less than 0.  Contact developers.
       ! </ERROR>
       CALL error_mesg('diag_manager_mod::average_tiles',&
            & "diag_field_id less than 0.  Contact developers.", FATAL)
    END IF

    ! Initialize local_missing_value
    IF ( input_fields(diag_field_id)%missing_value_present ) THEN
       local_missing_value = input_fields(diag_field_id)%missing_value
    ELSE
       local_missing_value = 0.0
    END IF

    ! Initialize s and out to zero.
    s(:,:) = 0.0
    out(:,:) = 0.0

    DO it = 1, SIZE(area,3)
       WHERE ( mask(:,:,it) )
          out(:,:) = out(:,:) + x(:,:,it)*area(:,:,it)
          s(:,:) = s(:,:) + area(:,:,it)
       END WHERE
    END DO

    WHERE ( s(:,:) > 0 )
       out(:,:) = out(:,:)/s(:,:)
    ELSEWHERE
       out(:,:) = local_missing_value
    END WHERE
  END SUBROUTINE average_tiles

  !> @return Integer writing_field
  INTEGER FUNCTION writing_field(out_num, at_diag_end, error_string, time)
    INTEGER, INTENT(in) :: out_num
    LOGICAL, INTENT(in) :: at_diag_end
    CHARACTER(len=*), INTENT(out) :: error_string
    TYPE(time_type), INTENT(in) :: time

    TYPE(time_type) :: middle_time
    TYPE(time_type) :: filename_time
    LOGICAL :: time_max, time_min, reduced_k_range, missvalue_present
    LOGICAL :: average, time_rms, need_compute, phys_window
    INTEGER :: in_num, file_num, freq, units
    INTEGER :: b1,b2,b3,b4 !< size of buffer along x,y,z,and diurnal axes
    INTEGER :: i, j, k, m
    REAL    :: missvalue, num
    real,allocatable,dimension(:,:,:,:) :: diurnal_buffer
    writing_field = 0

    need_compute = output_fields(out_num)%need_compute

    in_num = output_fields(out_num)%input_field
    IF ( input_fields(in_num)%static ) RETURN

    missvalue = input_fields(in_num)%missing_value
    missvalue_present = input_fields(in_num)%missing_value_present
    reduced_k_range = output_fields(out_num)%reduced_k_range
    phys_window = output_fields(out_num)%phys_window
    ! Is this output field being time averaged?
    average = output_fields(out_num)%time_average
    ! Are we taking the rms of the field?
    ! If so, then average is also .TRUE.
    time_rms = output_fields(out_num)%time_rms
    ! Looking for max and min value of this field over the sampling interval?
    time_max = output_fields(out_num)%time_max
    time_min = output_fields(out_num)%time_min
    file_num = output_fields(out_num)%output_file
    freq = files(file_num)%output_freq
    units = files(file_num)%output_units

    ! If average get size: Average intervals are last_output, next_output
    IF ( average ) THEN
       b1=SIZE(output_fields(out_num)%buffer,1)
       b2=SIZE(output_fields(out_num)%buffer,2)
       b3=SIZE(output_fields(out_num)%buffer,3)
       b4=SIZE(output_fields(out_num)%buffer,4)
       IF ( input_fields(in_num)%mask_variant ) THEN
          DO m=1, b4
             DO k=1, b3
                DO j=1, b2
                   DO i=1, b1
                      IF ( output_fields(out_num)%counter(i,j,k,m) > 0. )THEN
                         output_fields(out_num)%buffer(i,j,k,m) = &
                              & output_fields(out_num)%buffer(i,j,k,m)/output_fields(out_num)%counter(i,j,k,m)
                         IF ( time_rms ) output_fields(out_num)%buffer(i,j,k,m) = &
                              SQRT(output_fields(out_num)%buffer(i,j,k,m))
                      ELSE
                         output_fields(out_num)%buffer(i,j,k,m) =  missvalue
                      END IF
                   END DO
                END DO
             END DO
          END DO
       ELSE  !not mask variant
          DO m = 1, b4
             IF ( phys_window ) THEN
                IF ( need_compute .OR. reduced_k_range ) THEN
                   num = REAL(output_fields(out_num)%num_elements(m)/output_fields(out_num)%region_elements)
                ELSE
                   num = REAL(output_fields(out_num)%num_elements(m)/output_fields(out_num)%total_elements)
                END IF
             ELSE
                num = output_fields(out_num)%count_0d(m)
             END IF
             IF ( num > 0. ) THEN
                IF ( missvalue_present ) THEN
                   DO k=1, b3
                      DO j=1, b2
                         DO i=1, b1
                            IF ( output_fields(out_num)%buffer(i,j,k,m) /= missvalue ) THEN
                               output_fields(out_num)%buffer(i,j,k,m) = output_fields(out_num)%buffer(i,j,k,m)/num
                               IF ( time_rms ) output_fields(out_num)%buffer(i,j,k,m) =&
                                    & SQRT(output_fields(out_num)%buffer(i,j,k,m))
                            END IF
                         END DO
                      END DO
                   END DO
                ELSE
                   output_fields(out_num)%buffer(:,:,:,m) = output_fields(out_num)%buffer(:,:,:,m)/num
                   IF ( time_rms ) output_fields(out_num)%buffer(:,:,:,m) =&
                        & SQRT(output_fields(out_num)%buffer(:,:,:,m))
                END IF
             ELSE IF ( .NOT. at_diag_end ) THEN
                IF ( missvalue_present ) THEN
                   IF(ANY(output_fields(out_num)%buffer /= missvalue)) THEN
                      WRITE (error_string,'(a,"/",a)')&
                           & TRIM(input_fields(in_num)%module_name), &
                           & TRIM(output_fields(out_num)%output_name)
                      writing_field = -1
                      RETURN
                   END IF
                END IF
             END IF
          END DO
       END IF ! mask_variant
    ELSE IF ( time_min .OR. time_max ) THEN
       IF ( missvalue_present ) THEN
          WHERE ( ABS(output_fields(out_num)%buffer) == MIN_VALUE )
             output_fields(out_num)%buffer = missvalue
          END WHERE
       END IF ! if missvalue is NOT present buffer retains max_value or min_value
    END IF !average

    ! Output field
    IF ( at_diag_end .AND. freq == END_OF_RUN ) output_fields(out_num)%next_output = time
! if (time .eq. output_fields(out_num)%next_output) then
    if (output_fields(out_num)%n_diurnal_samples > 1) then
          !> allocate the buffer for diurnal data
          allocate(diurnal_buffer(size(output_fields(out_num)%buffer,1),size(output_fields(out_num)%buffer,2),&
                    size(output_fields(out_num)%buffer,4),size(output_fields(out_num)%buffer,3)))
          !> swap the last 2 axes in the data buffer to match the netcdf output order
          do i = 1,size(output_fields(out_num)%buffer,4)
          do j = 1,size(output_fields(out_num)%buffer,3)
               diurnal_buffer(:,:,i,j) = output_fields(out_num)%buffer(:,:,j,i)
          enddo
          enddo
    endif
    IF ( (output_fields(out_num)%time_ops) .AND. (.NOT. mix_snapshot_average_fields) ) THEN
       middle_time = (output_fields(out_num)%last_output+output_fields(out_num)%next_output)/2
       if (trim(files(file_num)%filename_time_bounds) == "begin") then
           filename_time = output_fields(out_num)%last_output
       elseif (trim(files(file_num)%filename_time_bounds) == "middle") then
           filename_time = middle_time
       elseif (trim(files(file_num)%filename_time_bounds) == "end") then
           filename_time = output_fields(out_num)%next_output
       endif

       if (output_fields(out_num)%n_diurnal_samples > 1) then
          CALL diag_data_out(file_num, out_num, diurnal_buffer, middle_time, &
                 & filename_time=filename_time)
       else
          CALL diag_data_out(file_num, out_num, output_fields(out_num)%buffer, middle_time, &
                 & filename_time=filename_time)
       endif
    ELSE
       if (output_fields(out_num)%n_diurnal_samples > 1) then
            CALL diag_data_out(file_num, out_num, &
                 & diurnal_buffer, output_fields(out_num)%next_output)
       else
            CALL diag_data_out(file_num, out_num, &
                 & output_fields(out_num)%buffer, output_fields(out_num)%next_output)
       endif
    END IF
!output_fields(out_num)%last_output = output_fields(out_num)%next_output
! endif
    IF ( at_diag_end ) RETURN

    ! Take care of cleaning up the time counters and the storeage size
    output_fields(out_num)%last_output = output_fields(out_num)%next_output
    IF ( freq == END_OF_RUN ) THEN
       output_fields(out_num)%next_output = time
    ELSE
       IF ( freq == EVERY_TIME ) THEN
          output_fields(out_num)%next_output = time
       ELSE
          output_fields(out_num)%next_output = output_fields(out_num)%next_next_output
          output_fields(out_num)%next_next_output = &
               & diag_time_inc(output_fields(out_num)%next_next_output, freq, units)
       END IF
       output_fields(out_num)%count_0d(:) = 0.0
       output_fields(out_num)%num_elements(:) = 0
       IF ( time_max ) THEN
          output_fields(out_num)%buffer = MAX_VALUE
       ELSE IF ( time_min ) THEN
          output_fields(out_num)%buffer = MIN_VALUE
       ELSE
          output_fields(out_num)%buffer = EMPTY
       END IF
       IF ( input_fields(in_num)%mask_variant .AND. average ) output_fields(out_num)%counter = 0.0
    END IF
    if (allocated(diurnal_buffer))deallocate(diurnal_buffer)
  END FUNCTION writing_field

  SUBROUTINE diag_manager_set_time_end(Time_end_in)
    TYPE (time_type), INTENT(in) :: Time_end_in

    Time_end = Time_end_in

  END SUBROUTINE diag_manager_set_time_end

  !-----------------------------------------------------------------------
  !> @brief The subroutine 'diag_send_complete_instant' allows the user to
  !! save diagnostic data on variable intervals (user defined in code logic)
  !! to the same file.  The argument (time_type) will be written to the
  !! time axis correspondingly.
  !!
  !> The user is responsible for any averaging of accumulated data
  !! as this routine is not designed for instantaneous values.  This routine
  !! works only for send_data calls within OpenMP regions as they are buffered
  !! until the complete signal is given.
  SUBROUTINE diag_send_complete_instant(time)
    TYPE (time_type), INTENT(in) :: time
    !--- local variables
    integer :: file, j, freq, in_num, file_num, out_num

    DO file = 1, num_files
      freq = files(file)%output_freq
      IF (freq == 0) then
        DO j = 1, files(file)%num_fields
          out_num = files(file)%fields(j)
          in_num = output_fields(out_num)%input_field
          IF ( (input_fields(in_num)%numthreads == 1) .AND.&
               & (input_fields(in_num)%active_omp_level.LE.1) ) CYCLE
          file_num = output_fields(out_num)%output_file
          CALL diag_data_out(file_num, out_num, &
               & output_fields(out_num)%buffer, time)
        END DO
      END IF
    END DO
  END SUBROUTINE diag_send_complete_instant

  !-----------------------------------------------------------------------
  !> @brief Saves diagnostic data for the given time value.
  SUBROUTINE diag_send_complete(time_step, err_msg)
    TYPE (time_type), INTENT(in)           :: time_step
    character(len=*), INTENT(out), optional :: err_msg

    type(time_type)    :: next_time, time
    integer            :: file, j, out_num, in_num, freq, status
    logical            :: local_output, need_compute
    CHARACTER(len=128) :: error_string

    IF ( Time_end == Time_zero ) THEN
       ! <ERROR STATUS="FATAL">
       !   diag_manager_set_time_end must be called before diag_send_complete
       ! </ERROR>
       CALL error_mesg('diag_manager_mod::diag_send_complete',&
            & "diag_manager_set_time_end must be called before diag_send_complete", FATAL)
    END IF

    DO file = 1, num_files
       freq = files(file)%output_freq
       DO j = 1, files(file)%num_fields
          out_num = files(file)%fields(j) !this is position of output_field in array output_fields
          in_num = output_fields(out_num)%input_field

          IF ( (input_fields(in_num)%numthreads == 1) .AND. (input_fields(in_num)%active_omp_level.LE.1) ) CYCLE
          IF ( output_fields(out_num)%static .OR. freq == END_OF_RUN ) CYCLE
          time = input_fields(in_num)%time
          IF ( time >= time_end ) CYCLE

          ! is this field output on a local domain only?
          local_output = output_fields(out_num)%local_output
          ! if local_output, does the current PE take part in send_data?
          need_compute = output_fields(out_num)%need_compute
          ! skip all PEs not participating in outputting this field
          IF ( local_output .AND. (.NOT.need_compute) ) CYCLE
          next_time = time + time_step

          IF ( next_time > output_fields(out_num)%next_output ) THEN
             ! A non-static field that has skipped a time level is an error
             IF ( next_time > output_fields(out_num)%next_next_output .AND. freq > 0 ) THEN
                IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
                   WRITE (error_string,'(a,"/",a)')&
                        & TRIM(input_fields(in_num)%module_name), &
                        & TRIM(output_fields(out_num)%output_name)
                   IF ( fms_error_handler('diag_send_complete',&
                        & 'module/output_field '//TRIM(error_string)//&
                        & ' is skipped one time level in output data', err_msg)) RETURN
                END IF
             END IF

             status = writing_field(out_num, .FALSE., error_string, next_time)
             IF ( status == -1 ) THEN
                IF ( mpp_pe() .EQ. mpp_root_pe() ) THEN
                   IF(fms_error_handler('diag_manager_mod::diag_send_complete','module/output_field '//TRIM(error_string)//&
                        & ', write EMPTY buffer', err_msg)) RETURN
                END IF
             END IF
          END IF  !time > output_fields(out_num)%next_output
       END DO
    END DO

  END SUBROUTINE diag_send_complete

  !> @brief Flushes diagnostic buffers where necessary. Close diagnostics files.
  !!   A warning will be issued here if a field in diag_table is not registered
  SUBROUTINE diag_manager_end(time)
    TYPE(time_type), INTENT(in) :: time

    INTEGER :: file

    IF ( do_diag_field_log ) THEN
       CALL mpp_close (diag_log_unit)
    END IF
    DO file = 1, num_files
       CALL closing_file(file, time)
    END DO
    if (allocated(fileobjU)) deallocate(fileobjU)
    if (allocated(fileobj)) deallocate(fileobj)
    if (allocated(fileobjND)) deallocate(fileobjND)
    if (allocated(fnum_for_domain)) deallocate(fnum_for_domain)
  END SUBROUTINE diag_manager_end

  !> @brief Replaces diag_manager_end; close just one file: files(file)
  SUBROUTINE closing_file(file, time)
    INTEGER, INTENT(in) :: file
    TYPE(time_type), INTENT(in) :: time

    INTEGER :: j, i, input_num, freq, status, loop1, loop2
    INTEGER :: stdout_unit
    LOGICAL :: reduced_k_range, need_compute, local_output
    CHARACTER(len=128) :: message
    real,allocatable,dimension(:,:,:,:) :: diurnal_buffer

    stdout_unit = stdout()

    ! Output all registered, non_static output_fields
    DO j = 1, files(file)%num_fields
       i = files(file)%fields(j) !this is position of output_field in array output_fields

       ! is this field output on a local domain only?
       local_output = output_fields(i)%local_output
       ! if local_output, does the current PE take part in send_data?
       need_compute = output_fields(i)%need_compute

       reduced_k_range = output_fields(i)%reduced_k_range

       ! skip all PEs not participating in outputting this field
       IF ( local_output .AND. (.NOT. need_compute) ) CYCLE
       ! skip fields that were not registered or non-static
       input_num = output_fields(i)%input_field
       IF ( input_fields(input_num)%static ) CYCLE
       IF ( .NOT.input_fields(input_num)%register ) CYCLE
       freq = files(file)%output_freq
       IF ( freq /= END_OF_RUN .AND. files(file)%file_unit < 0 &
            & .AND. ALL(output_fields(i)%num_elements(:) == 0)&
            & .AND. ALL(output_fields(i)%count_0d(:) == 0) ) CYCLE
       ! Is it time to output for this field; CAREFUL ABOUT >= vs > HERE
       ! For end should be >= because no more data is coming
       IF ( time >= output_fields(i)%next_output .OR. freq == END_OF_RUN ) THEN
          IF ( time >= output_fields(i)%next_next_output .AND. freq > 0 ) THEN
             WRITE (message,'(a,"/",a)') TRIM(input_fields(input_num)%module_name), &
                  & TRIM(output_fields(i)%output_name)
             ! <ERROR STATUS="WARNING">
             !   <input_fields(input_num)%module_name>/<output_fields(i)%output_name> skip one time
             !   level, maybe send_data never called
             ! </ERROR>
             IF ( mpp_pe() .EQ. mpp_root_pe() ) &
                  & CALL error_mesg('diag_manager_mod::closing_file', 'module/output_field ' //&
                  & TRIM(message)//', skip one time level, maybe send_data never called', WARNING)
             status = writing_field(i, .TRUE.,message,time)
          ELSE
             status = writing_field(i, .TRUE., message, time)
          END IF
       ELSEIF ( .NOT.output_fields(i)%written_once ) THEN
          ! <ERROR STATUS="NOTE">
          !   <output_fields(i)%output_name) NOT available, check if output interval > runlength.
          !   NetCDF fill_values are written
          ! </ERROR>
          CALL error_mesg('Potential error in diag_manager_end ',&
               & TRIM(output_fields(i)%output_name)//' NOT available,'//&
               & ' check if output interval > runlength. Netcdf fill_values are written', NOTE)
          output_fields(i)%buffer = FILL_VALUE
          if (output_fields(i)%n_diurnal_samples > 1) then
               !> allocate the buffer for diurnal data
               if (.not. allocated(diurnal_buffer)) &
                    allocate(diurnal_buffer(size(output_fields(i)%buffer,1),size(output_fields(i)%buffer,2),&
                    size(output_fields(i)%buffer,4),size(output_fields(i)%buffer,3)))
               !> swap the last 2 axes in the data buffer to match the netcdf output order
               do loop1 = 1,size(output_fields(i)%buffer,4)
               do loop2 = 1,size(output_fields(i)%buffer,3)
                    diurnal_buffer(:,:,loop1,loop2) = output_fields(i)%buffer(:,:,loop2,loop1)
               enddo
               enddo
               CALL diag_data_out(file, i, diurnal_buffer, time, .TRUE.)
          else
          CALL diag_data_out(file, i, output_fields(i)%buffer, time, .TRUE.)
          endif
       END IF
    END DO
    ! Now it's time to output static fields
    CALL write_static(file)

    ! Write out the number of bytes of data saved to this file
    IF ( write_bytes_in_file ) THEN
       CALL mpp_sum (files(file)%bytes_written)
       IF ( mpp_pe() == mpp_root_pe() )&
            & WRITE (stdout_unit,'(a,i12,a,a)') 'Diag_Manager: ',files(file)%bytes_written, &
            & ' bytes of data written to file ',TRIM(files(file)%name)
    END IF
     if (allocated(diurnal_buffer)) deallocate(diurnal_buffer)
  END SUBROUTINE closing_file

  !> @brief Initialize Diagnostics Manager.
  !! @details Open and read diag_table. Select fields and files for diagnostic output.
  SUBROUTINE diag_manager_init(diag_model_subset, time_init, err_msg)
    INTEGER, OPTIONAL, INTENT(IN) :: diag_model_subset
    INTEGER, DIMENSION(6), OPTIONAL, INTENT(IN) :: time_init !< Model time diag_manager initialized
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg

    CHARACTER(len=*), PARAMETER :: SEP = '|'

    INTEGER, PARAMETER :: FltKind = R4_KIND
    INTEGER, PARAMETER :: DblKind = R8_KIND
    INTEGER :: diag_subset_output
    INTEGER :: mystat
    INTEGER, ALLOCATABLE, DIMENSION(:) :: pelist
    INTEGER :: stdlog_unit, stdout_unit
    integer :: j
    CHARACTER(len=256) :: err_msg_local

    NAMELIST /diag_manager_nml/ append_pelist_name, mix_snapshot_average_fields, max_output_fields, &
         & max_input_fields, max_axes, do_diag_field_log, write_bytes_in_file, debug_diag_manager,&
         & max_num_axis_sets, max_files, use_cmor, issue_oor_warnings,&
         & oor_warnings_fatal, max_out_per_in_field, flush_nc_files, region_out_use_alt_value, max_field_attributes,&
         & max_file_attributes, max_axis_attributes, prepend_date, use_mpp_io

    ! If the module was already initialized do nothing
    IF ( module_is_initialized ) RETURN

    ! Clear the err_msg variable if contains any residual information
    IF ( PRESENT(err_msg) ) err_msg = ''

    ! Initialize diag_util_mod and diag_data_mod
    ! These init routine only write out the version number to the log file
    call diag_util_init()
    call diag_data_init()

    ! Determine pack_size from how many bytes a real value has (how compiled)
    pack_size = SIZE(TRANSFER(0.0_DblKind, (/0.0, 0.0, 0.0, 0.0/)))
    IF ( pack_size.NE.1 .AND. pack_size.NE.2 ) THEN
       IF ( fms_error_handler('diag_manager_mod::diag_manager_init', 'unknown pack_size.  Must be 1, or 2.', err_msg) ) RETURN
    END IF

    ! Get min and max values for real(kind=R4_KIND)
    min_value = HUGE(0.0_FltKind)
    max_value = -min_value

    ! get stdlog and stdout unit number
    stdlog_unit = stdlog()
    stdout_unit = stdout()

    ! version number to logfile
    CALL write_version_number("DIAG_MANAGER_MOD", version)

    Time_zero = set_time(0,0)
    !--- initialize time_end to time_zero
    Time_end  = Time_zero
    diag_subset_output = DIAG_ALL
    IF ( PRESENT(diag_model_subset) ) THEN
       IF ( diag_model_subset >= DIAG_OTHER .AND. diag_model_subset <= DIAG_ALL ) THEN
          diag_subset_output = diag_model_subset
       ELSE
          IF ( fms_error_handler('diag_manager_mod::diag_manager_init', 'invalid value of diag_model_subset',err_msg) ) RETURN
       END IF
    END IF

    READ (input_nml_file, NML=diag_manager_nml, IOSTAT=mystat)
    ! Check the status of reading the diag_manager_nml

    IF ( check_nml_error(IOSTAT=mystat, NML_NAME='DIAG_MANAGER_NML') < 0 ) THEN
       IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_manager_mod::diag_manager_init', 'DIAG_MANAGER_NML not found in input nml file.  Using defaults.',&
               & WARNING)
       END IF
    END IF

    IF ( mpp_pe() == mpp_root_pe() ) THEN
       WRITE (stdlog_unit, diag_manager_nml)
    END IF

    ! Issue note about using the CMOR missing value.
    IF ( use_cmor ) THEN
       err_msg_local = ''
       WRITE (err_msg_local,'(ES8.1E2)') CMOR_MISSING_VALUE
       CALL error_mesg('diag_manager_mod::diag_manager_init', 'Using CMOR missing value ('//TRIM(err_msg_local)//').', NOTE)
    END IF

    ! Issue note if attempting to set diag_manager_nml::max_files larger than
    ! mpp_get_maxunits() -- Default is 1024 set in mpp_io.F90
    IF ( max_files .GT. mpp_get_maxunits() ) THEN
       err_msg_local = ''
       WRITE (err_msg_local,'(A,I6,A,I6,A,I6,A)') "DIAG_MANAGER_NML variable 'max_files' (",max_files,") is larger than '",&
            & mpp_get_maxunits(),"'.  Forcing 'max_files' to be ",mpp_get_maxunits(),"."
       CALL error_mesg('diag_manager_mod::diag_managet_init', TRIM(err_msg_local), NOTE)
       max_files = mpp_get_maxunits()
    END IF

    ! How to handle Out of Range Warnings.
    IF ( oor_warnings_fatal ) THEN
       oor_warning = FATAL
       CALL error_mesg('diag_manager_mod::diag_manager_init', 'Out &
            &of Range warnings are fatal.', NOTE)
    ELSEIF ( .NOT.issue_oor_warnings ) THEN
       CALL error_mesg('diag_manager_mod::diag_manager_init', 'Out &
            &of Range warnings will be ignored.', NOTE)
    END IF

    IF ( mix_snapshot_average_fields ) THEN
       IF ( mpp_pe() == mpp_root_pe() ) THEN
          CALL error_mesg('diag_manager_mod::diag_manager_init', 'Setting diag_manager_nml variable '//&
               & 'mix_snapshot_average_fields = .TRUE. will cause ERRORS in the time coordinates '//&
               & 'of all time averaged fields.  Strongly recommend setting mix_snapshot_average_fields '//&
               & '= .FALSE.', WARNING)
       END IF
    END IF
    ALLOCATE(output_fields(max_output_fields))
    ALLOCATE(input_fields(max_input_fields))
    DO j = 1, max_input_fields
      ALLOCATE(input_fields(j)%output_fields(MAX_OUT_PER_IN_FIELD))
    END DO
    ALLOCATE(files(max_files))
    if (.not.use_mpp_io) then
      ALLOCATE(fileobjU(max_files))
      ALLOCATE(fileobj(max_files))
      ALLOCATE(fileobjND(max_files))
      ALLOCATE(fnum_for_domain(max_files))
    !> Initialize fnum_for_domain with "dn" which stands for done
      fnum_for_domain(:) = "dn"
       CALL error_mesg('diag_manager_mod::diag_manager_init',&
               & 'diag_manager is using fms2_io', NOTE)
    else
       CALL error_mesg('diag_manager_mod::diag_manager_init',&
             &'MPP_IO is no longer supported.  Please remove use_mpp_io from diag_manager_nml namelist',&
              &FATAL)
    endif
    ALLOCATE(pelist(mpp_npes()))
    CALL mpp_get_current_pelist(pelist, pelist_name)

    ! set the diag_init_time if time_init present.  Otherwise, set it to base_time
    IF ( PRESENT(time_init) ) THEN
       diag_init_time = set_date(time_init(1), time_init(2), time_init(3), time_init(4),&
            & time_init(5), time_init(6))
    ELSE
       diag_init_time = base_time
       IF ( prepend_date .EQV. .TRUE. ) THEN
          CALL error_mesg('diag_manager_mod::diag_manager_init',&
               & 'prepend_date only supported when diag_manager_init is called with time_init present.', NOTE)
          prepend_date = .FALSE.
       END IF
    END IF

    CALL parse_diag_table(DIAG_SUBSET=diag_subset_output, ISTAT=mystat, ERR_MSG=err_msg_local)
    IF ( mystat /= 0 ) THEN
       IF ( fms_error_handler('diag_manager_mod::diag_manager_init',&
            & 'Error parsing diag_table. '//TRIM(err_msg_local), err_msg) ) RETURN
    END IF

    !initialize files%bytes_written to zero
    files(:)%bytes_written = 0

    ! open diag field log file
    IF ( do_diag_field_log.AND.mpp_pe().EQ.mpp_root_pe() ) THEN
       CALL mpp_open(diag_log_unit, 'diag_field_log.out', nohdrs=.TRUE.)
       WRITE (diag_log_unit,'(777a)') &
            & 'Module',        SEP, 'Field',          SEP, 'Long Name',    SEP,&
            & 'Units',         SEP, 'Number of Axis', SEP, 'Time Axis',    SEP,&
            & 'Missing Value', SEP, 'Min Value',      SEP, 'Max Value',    SEP,&
            & 'AXES LIST'
    END IF

    module_is_initialized = .TRUE.
    ! create axis_id for scalars here
    null_axis_id = diag_axis_init('scalar_axis', (/0./), 'none', 'N', 'none')
    RETURN
  END SUBROUTINE diag_manager_init

  !> @brief Return base time for diagnostics.
  !! @return time_type get_base_time
  !! @details Return base time for diagnostics (note: base time must be >= model time).
  TYPE(time_type) FUNCTION get_base_time ()
    ! <ERROR STATUS="FATAL">
    !   MODULE has not been initialized
    ! </ERROR>
    IF ( .NOT.module_is_initialized ) CALL error_mesg('diag_manager_mod::get_base_time', &
         & 'module has not been initialized', FATAL)
    get_base_time = base_time
  END FUNCTION get_base_time

  !> @brief Return base date for diagnostics.
  !! @details Return date information for diagnostic reference time.
  SUBROUTINE get_base_date(year, month, day, hour, minute, second)
    INTEGER, INTENT(out) :: year, month, day, hour, minute, second

    ! <ERROR STATUS="FATAL">module has not been initialized</ERROR>
    IF (.NOT.module_is_initialized) CALL error_mesg ('diag_manager_mod::get_base_date', &
         & 'module has not been initialized', FATAL)
    year   = base_year
    month  = base_month
    day    = base_day
    hour   = base_hour
    minute = base_minute
    second = base_second
  END SUBROUTINE get_base_date

  !> @brief Determine whether data is needed for the current model time step.
  !! @return Logical need_data
  !! @details Determine whether data is needed for the current model time step.
  !!     Since diagnostic data are buffered, the "next" model time is passed
  !!     instead of the current model time. This call can be used to minimize
  !!     overhead for complicated diagnostics.
  LOGICAL FUNCTION need_data(diag_field_id, next_model_time)
    TYPE(time_type), INTENT(in) :: next_model_time !< next_model_time = current model time + model time_step
    INTEGER, INTENT(in) :: diag_field_id

    INTEGER :: i, out_num

    need_data = .FALSE.
    IF ( diag_field_id < 0 ) RETURN ! this field is unused
    DO i = 1, input_fields(diag_field_id)%num_output_fields
       ! Get index to an output field
       out_num = input_fields(diag_field_id)%output_fields(i)
       IF ( .NOT.output_fields(out_num)%static ) THEN
          IF ( next_model_time > output_fields(out_num)%next_output ) need_data=.TRUE.
          ! Is this output field being time averaged?
          ! assume average data based on every timestep
          ! needs to be changed when different forms of averaging are implemented
          IF ( output_fields(out_num)%time_average) need_data = .TRUE.
       END IF
    END DO
    RETURN
  END FUNCTION need_data

  !> @brief Finds or initializes a diurnal time axis and returns its' ID.
  !! @return Integer init_diurnal_axis
  !! @details Given number of time intervals in the day, finds or initializes a diurnal time axis
  !!     and returns its ID. It uses get_base_date, so should be in the file where it's accessible.
  !!     The units are 'days since BASE_DATE', all diurnal axes belong to the set 'diurnal'
  INTEGER FUNCTION init_diurnal_axis(n_samples)
    INTEGER, INTENT(in) :: n_samples !< number of intervals during the day

    REAL :: DATA  (n_samples)   !< central points of time intervals
    REAL :: edges (n_samples+1) !< boundaries of time intervals
    INTEGER :: edges_id !< id of the corresponding edges
    INTEGER :: i
    INTEGER :: year !< components of the base date
    INTEGER :: month !< components of the base date
    INTEGER :: day !< components of the base date
    INTEGER :: hour !< components of the base date
    INTEGER :: minute !< components of the base date
    INTEGER :: second !< components of the base date
    CHARACTER(32)  :: name  !< name of the axis
    CHARACTER(128) :: units !< units of time

    CALL get_base_date(year, month, day, hour, minute, second)
    WRITE (units,11) 'hours', year, month, day, hour, minute, second
11  FORMAT(a,' since ',i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2)
    ! compute central points and units
    edges(1) = 0.0
    DO i = 1, n_samples
       DATA (i) = 24.0*(REAL(i)-0.5)/n_samples
       edges(i+1) = 24.0* REAL(i)/n_samples
    END DO

    ! define edges
    name = ''
    WRITE (name,'(a,i2.2)') 'time_of_day_edges_', n_samples
    edges_id = get_axis_num(name, 'diurnal')
    IF ( edges_id <= 0 ) THEN
       edges_id =  diag_axis_init(name,edges,units,'N','time of day edges', set_name='diurnal')
    END IF

    ! define axis itself
    name = ''
    WRITE (name,'(a,i2.2)') 'time_of_day_', n_samples
    init_diurnal_axis = get_axis_num(name, 'diurnal')
    IF ( init_diurnal_axis <= 0 ) THEN
       init_diurnal_axis = diag_axis_init(name, DATA, units, 'N', 'time of day', set_name='diurnal', edges=edges_id)
    END IF
  END FUNCTION init_diurnal_axis

  SUBROUTINE diag_field_attribute_init(diag_field_id, name, type, cval, ival, rval)
    INTEGER, INTENT(in) :: diag_field_id !< input field ID, obtained from diag_manager_mod::register_diag_field.
    CHARACTER(len=*), INTENT(in) :: name !< Name of the attribute
    INTEGER, INTENT(in) :: type !< NetCDF type (NF90_FLOAT, NF90_INT, NF90_CHAR)
    CHARACTER(len=*), INTENT(in), OPTIONAL :: cval !< Character string attribute value
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: ival !< Integer attribute value(s)
    REAL, DIMENSION(:), INTENT(in), OPTIONAL :: rval !< Real attribute value(s)

    INTEGER :: istat, length, i, j, this_attribute, out_field
    CHARACTER(len=1024) :: err_msg

    IF ( .NOT.first_send_data_call ) THEN
       ! Call error due to unable to add attribute after send_data called
       ! <ERROR STATUS="FATAL">
       !   Attempting to add attribute <name> to module/input_field <module_name>/<field_name>
       !   after first send_data call.  Too late.
       ! </ERROR>
       CALL error_mesg('diag_manager_mod::diag_field_add_attribute', 'Attempting to add attribute "'&
            &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
            &//TRIM(input_fields(diag_field_id)%field_name)//'" after first send_data call.  Too late.', FATAL)
    END IF

    ! Simply return if diag_field_id <= 0 --- not in diag_table
    IF ( diag_field_id .LE. 0 ) THEN
       RETURN
    ELSE
       DO j=1,input_fields(diag_field_id)%num_output_fields
          out_field = input_fields(diag_field_id)%output_fields(j)

          ! Allocate memory for the attributes
          CALL attribute_init(output_fields(out_field))

          ! Check if attribute already exists
          this_attribute = 0
          DO i=1, output_fields(out_field)%num_attributes
             IF ( TRIM(output_fields(out_field)%attributes(i)%name) .EQ. TRIM(name) ) THEN
                this_attribute = i
                EXIT
             END IF
          END DO

          IF ( this_attribute.NE.0 .AND. (type.EQ.NF90_INT .OR. type.EQ.NF90_FLOAT) ) THEN
             ! <ERROR STATUS="FATAL">
             !   Attribute <name> already defined for module/input_field <module_name>/<field_name>.
             !   Contact the developers
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_field_add_attribute',&
                  & 'Attribute "'//TRIM(name)//'" already defined for module/input_field "'&
                  &//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                  &//TRIM(input_fields(diag_field_id)%field_name)//'".  Contact the developers.', FATAL)
          ELSE IF ( this_attribute.NE.0 .AND. type.EQ.NF90_CHAR .AND. debug_diag_manager ) THEN
             ! <ERROR STATUS="NOTE">
             !   Attribute <name> already defined for module/input_field <module_name>/<field_name>.
             !   Prepending.
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_field_add_attribute',&
                  & 'Attribute "'//TRIM(name)//'" already defined for module/input_field "'&
                  &//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                  &//TRIM(input_fields(diag_field_id)%field_name)//'".  Prepending.', NOTE)
          ELSE IF ( this_attribute.EQ.0 ) THEN
             ! Defining a new attribute
             ! Increase the number of field attributes
             this_attribute = output_fields(out_field)%num_attributes + 1
             ! Checking to see if num_attributes == max_field_attributes, and return error message
             IF ( this_attribute .GT. max_field_attributes ) THEN
                ! <ERROR STATUS="FATAL">
                !   Number of attributes exceeds max_field_attributes for attribute <name> to module/input_field <module_name>/<field_name>.
                !   Increase diag_manager_nml:max_field_attributes.
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::diag_field_add_attribute',&
                     & 'Number of attributes exceeds max_field_attributes for attribute "'&
                     &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                     &//TRIM(input_fields(diag_field_id)%field_name)//'".  Increase diag_manager_nml:max_field_attributes.',&
                     & FATAL)
             ELSE
                output_fields(out_field)%num_attributes = this_attribute
                ! Set name and type
                output_fields(out_field)%attributes(this_attribute)%name = name
                output_fields(out_field)%attributes(this_attribute)%type = type
                ! Initialize catt to a blank string, as len_trim doesn't always work on an uninitialized string
                output_fields(out_field)%attributes(this_attribute)%catt = ''
             END IF
          END IF

          SELECT CASE (type)
          CASE (NF90_INT)
             IF ( .NOT.PRESENT(ival) ) THEN
                ! <ERROR STATUS="FATAL">
                !   Number type claims INTEGER, but ival not present for attribute <name> to module/input_field <module_name>/<field_name>.
                !   Contact the developers.
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::diag_field_add_attribute',&
                     & 'Attribute type claims INTEGER, but ival not present for attribute "'&
                     &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                     &//TRIM(input_fields(diag_field_id)%field_name)//'". Contact then developers.', FATAL)
             END IF
             length = SIZE(ival)
             ! Allocate iatt(:) to size of ival
             ALLOCATE(output_fields(out_field)%attributes(this_attribute)%iatt(length), STAT=istat)
             IF ( istat.NE.0 ) THEN
                ! <ERROR STATUS="FATAL">
                !   Unable to allocate iatt for attribute <name> to module/input_field <module_name>/<field_name>
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::diag_field_add_attribute', 'Unable to allocate iatt for attribute "'&
                     &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                     &//TRIM(input_fields(diag_field_id)%field_name)//'"', FATAL)
             END IF
             ! Set remaining fields
             output_fields(out_field)%attributes(this_attribute)%len = length
             output_fields(out_field)%attributes(this_attribute)%iatt = ival
          CASE (NF90_FLOAT)
             IF ( .NOT.PRESENT(rval) ) THEN
                ! <ERROR STATUS="FATAL">
                !   Attribute type claims READ, but rval not present for attribute <name> to module/input_field <module_name>/<field_name>.
                !   Contact the developers.
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::diag_field_add_attribute',&
                     & 'Attribute type claims REAL, but rval not present for attribute "'&
                     &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                     &//TRIM(input_fields(diag_field_id)%field_name)//'". Contact the developers.', FATAL)
             END IF
             length = SIZE(rval)
             ! Allocate iatt(:) to size of rval
             ALLOCATE(output_fields(out_field)%attributes(this_attribute)%fatt(length), STAT=istat)
             IF ( istat.NE.0 ) THEN
                ! <ERROR STATUS="FATAL">
                !   Unable to allocate fatt for attribute <name> to module/input_field <module_name>/<field_name>
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::diag_field_add_attribute', 'Unable to allocate fatt for attribute "'&
                     &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                     &//TRIM(input_fields(diag_field_id)%field_name)//'"', FATAL)
             END IF
             ! Set remaining fields
             output_fields(out_field)%attributes(this_attribute)%len = length
             output_fields(out_field)%attributes(this_attribute)%fatt = rval
          CASE (NF90_CHAR)
             IF ( .NOT.PRESENT(cval) ) THEN
                ! <ERROR STATUS="FATAL">
                !   Attribute type claims CHARACTER, but cval not present for attribute <name> to module/input_field <module_name>/<field_name>.
                !   Contact the developers.
                ! </ERROR>
                CALL error_mesg('diag_manager_mod::diag_field_add_attribute',&
                     & 'Attribute type claims CHARACTER, but cval not present for attribute "'&
                     &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                     &//TRIM(input_fields(diag_field_id)%field_name)//'". Contact the developers.', FATAL)
             END IF
             CALL prepend_attribute(output_fields(out_field), TRIM(name), TRIM(cval))
          CASE default
             ! <ERROR STATUS="FATAL">
             !   Unknown attribute type for attribute <name> to module/input_field <module_name>/<field_name>.
             !   Contact the developers.
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_field_add_attribute', 'Unknown attribute type for attribute "'&
                  &//TRIM(name)//'" to module/input_field "'//TRIM(input_fields(diag_field_id)%module_name)//'/'&
                  &//TRIM(input_fields(diag_field_id)%field_name)//'". Contact the developers.', FATAL)
          END SELECT
       END DO
    END IF
  END SUBROUTINE diag_field_attribute_init

  !> @brief Add a scalar real attribute to the diag field corresponding to a given id
  SUBROUTINE diag_field_add_attribute_scalar_r(diag_field_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_field_id !< ID number for field to add attribute to
    CHARACTER(len=*), INTENT(in) :: att_name !< new attribute name
    REAL, INTENT(in) :: att_value !< new attribute value

    CALL diag_field_add_attribute_r1d(diag_field_id, att_name, (/ att_value /))
  END SUBROUTINE diag_field_add_attribute_scalar_r

  !> @brief Add a scalar integer attribute to the diag field corresponding to a given id
  SUBROUTINE diag_field_add_attribute_scalar_i(diag_field_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_field_id !< ID number for field to add attribute to
    CHARACTER(len=*), INTENT(in) :: att_name !< new attribute name
    INTEGER, INTENT(in) :: att_value !< new attribute value

    CALL diag_field_add_attribute_i1d(diag_field_id, att_name, (/ att_value /))
  END SUBROUTINE diag_field_add_attribute_scalar_i

  !> @brief Add a scalar character attribute to the diag field corresponding to a given id
  SUBROUTINE diag_field_add_attribute_scalar_c(diag_field_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_field_id !< ID number for field to add attribute to
    CHARACTER(len=*), INTENT(in) :: att_name !< new attribute name
    CHARACTER(len=*), INTENT(in) :: att_value !< new attribute value

    CALL diag_field_attribute_init(diag_field_id, att_name, NF90_CHAR, cval=att_value)
  END SUBROUTINE diag_field_add_attribute_scalar_c

  !> @brief Add a real 1D array attribute to the diag field corresponding to a given id
  SUBROUTINE diag_field_add_attribute_r1d(diag_field_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_field_id !< ID number for field to add attribute to
    CHARACTER(len=*), INTENT(in) :: att_name !< new attribute name
    REAL, DIMENSION(:), INTENT(in) :: att_value !< new attribute value

    INTEGER :: num_attributes, len
    CHARACTER(len=512) :: err_msg

    CALL diag_field_attribute_init(diag_field_id, att_name, NF90_FLOAT, rval=att_value)
  END SUBROUTINE diag_field_add_attribute_r1d

  !> @brief Add an integer 1D array attribute to the diag field corresponding to a given id
  SUBROUTINE diag_field_add_attribute_i1d(diag_field_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_field_id !< ID number for field to add attribute to
    CHARACTER(len=*), INTENT(in) :: att_name !< new attribute name
    INTEGER, DIMENSION(:), INTENT(in) :: att_value !< new attribute value

    CALL diag_field_attribute_init(diag_field_id, att_name, NF90_INT, ival=att_value)
  END SUBROUTINE diag_field_add_attribute_i1d

  !> @brief Add the cell_measures attribute to a diag out field
  !!
  !> Add the cell_measures attribute to a give diag field.  This is useful if the
  !! area/volume fields for the diagnostic field are defined in another module after
  !! the diag_field.
  SUBROUTINE diag_field_add_cell_measures(diag_field_id, area, volume)
    INTEGER, INTENT(in) :: diag_field_id
    INTEGER, INTENT(in), OPTIONAL :: area !< diag ids of area
    INTEGER, INTENT(in), OPTIONAL :: volume !< diag ids of volume

    integer :: j, ind

    IF ( diag_field_id.GT.0 ) THEN
       IF ( .NOT.PRESENT(area) .AND. .NOT.present(volume) ) THEN
          CALL ERROR_MESG('diag_manager_mod::diag_field_add_cell_measures', &
               & 'either area or volume arguments must be present', FATAL )
       END IF

       DO j=1, input_fields(diag_field_id)%num_output_fields
          ind = input_fields(diag_field_id)%output_fields(j)
          CALL init_field_cell_measures(output_fields(ind), area=area, volume=volume)
       END DO
    END IF
  END SUBROUTINE diag_field_add_cell_measures

END MODULE diag_manager_mod
!> @}
! close documentation grouping
