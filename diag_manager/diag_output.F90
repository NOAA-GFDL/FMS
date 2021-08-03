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
!> @defgroup diag_output_mod diag_output_mod
!> @ingroup diag_manager
!! @brief diag_output_mod is an integral part of
!!   diag_manager_mod. Its function is to write axis-meta-data,
!!   field-meta-data and field data.
!! @author Seth Underwood

!> @file
!> @brief File for @ref diag_output_mod

!> @addtogroup diag_output_mod
!> @{
MODULE diag_output_mod

use platform_mod
use,intrinsic :: iso_fortran_env, only: real128
use,intrinsic :: iso_c_binding, only: c_double,c_float,c_int64_t, &
                                      c_int32_t,c_int16_t,c_intptr_t
! use_mpp_io = .false.
  USE mpp_io_mod, ONLY: axistype, fieldtype, mpp_io_init, &
       & mpp_get_id, MPP_WRONLY, MPP_OVERWR,&
       & MPP_NETCDF, MPP_MULTI, MPP_SINGLE, mpp_get_field_name, &
       & fillin_fieldtype
  USE mpp_domains_mod, ONLY: domain1d, domain2d, mpp_define_domains, mpp_get_pelist,&
       &  mpp_get_global_domain, mpp_get_compute_domains, null_domain1d, null_domain2d,&
       & domainUG, null_domainUG, CENTER, EAST, NORTH, mpp_get_compute_domain,&
       & OPERATOR(.NE.), mpp_get_layout, OPERATOR(.EQ.), mpp_get_io_domain, &
       & mpp_get_compute_domain, mpp_get_global_domain
  USE mpp_mod, ONLY: mpp_npes, mpp_pe, mpp_root_pe, mpp_get_current_pelist
  USE diag_axis_mod, ONLY: diag_axis_init, get_diag_axis, get_axis_length,&
       & get_axis_global_length, get_domain1d, get_domain2d, get_axis_aux, get_tile_count,&
       & get_domainUG, get_diag_axis_name
  USE diag_data_mod, ONLY: pack_size, diag_fieldtype, diag_global_att_type, CMOR_MISSING_VALUE, diag_atttype, files
  USE time_manager_mod, ONLY: get_calendar_type, valid_calendar_types
  USE fms_mod, ONLY: error_mesg, mpp_pe, write_version_number, fms_error_handler, FATAL, note

#ifdef use_netCDF
  USE netcdf, ONLY: NF90_INT, NF90_FLOAT, NF90_CHAR
#endif

  use mpp_domains_mod, only: mpp_get_UG_io_domain
  use mpp_domains_mod, only: mpp_get_UG_domain_npes
  use mpp_domains_mod, only: mpp_get_UG_domain_pelist
  use mpp_mod,         only: mpp_gather
  use mpp_mod,         only: uppercase,lowercase
  use fms2_io_mod
  use axis_utils2_mod,   only: axis_edges


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: diag_output_init, write_axis_meta_data, write_field_meta_data, done_meta_data,&
       & diag_fieldtype, get_diag_global_att, set_diag_global_att
  PUBLIC :: diag_field_write, diag_write_time !< use_mpp_io = .false.
  TYPE(diag_global_att_type), SAVE :: diag_global_att

  INTEGER, PARAMETER      :: NETCDF1 = 1
  INTEGER, PARAMETER      :: mxch  = 128
  INTEGER, PARAMETER      :: mxchl = 256
  INTEGER                 :: current_file_unit = -1
  INTEGER, DIMENSION(2,2) :: max_range = RESHAPE((/ -32767, 32767, -127,   127 /),(/2,2/))
  INTEGER, DIMENSION(2)   :: missval = (/ -32768, -128 /)

  INTEGER, PARAMETER      :: max_axis_num = 20
  INTEGER                 :: num_axis_in_file = 0
  INTEGER, DIMENSION(max_axis_num) :: axis_in_file
  LOGICAL, DIMENSION(max_axis_num) :: time_axis_flag, edge_axis_flag
  TYPE(axistype), DIMENSION(max_axis_num), SAVE :: Axis_types

  LOGICAL :: module_is_initialized = .FALSE.

  ! Include variable "version" to be written to log file.
  character(len=*), parameter :: version = '2020.03'
  !> @}

  !> Write diag field using @ref fms2_io
  !> @ingroup diag_output_mod
  interface diag_field_write
     module procedure diag_field_write_field
     module procedure diag_field_write_varname
  end interface

  !> Initialize output for writing.
  !> @ingroup diag_output_mod
  interface diag_output_init
     module procedure diag_output_init_fms2_io
  end interface

  !> Writes axis metadata to a file.
  !> @ingroup diag_output_mod
  interface write_axis_meta_data
     module procedure write_axis_meta_data_fms2_io
  end interface

  !> Writes field metadata to a file.
  !> @ingroup diag_output_mod
  interface write_field_meta_data
     module procedure write_field_meta_data_fms2_io
  end interface

  !> Private interface to write metadata for an attribute to a file.
  !!
  !> @note Added for mpp_io support
  !> @ingroup diag_output_mod
  interface write_attribute_meta
     module procedure write_attribute_meta_fms2_io
  end interface

!> @addtogroup diag_output_mod
!> @{
CONTAINS

  !> @brief Registers the time axis and opens the output file.
  SUBROUTINE diag_output_init_fms2_io (file_name, FORMAT, file_title, file_unit,&
       & all_scalar_or_1d, domain, domainU, fileobj, fileobjU, fileobjND, fnum_domain, &
       & attributes)
    CHARACTER(len=*), INTENT(in)  :: file_name !< Output file name
    CHARACTER(len=*), INTENT(in)  :: file_title !< Descriptive title for the file
    INTEGER         , INTENT(in)  :: FORMAT !< File format (Currently only 'NETCDF' is valid)
    INTEGER         , INTENT(out) :: file_unit !< File unit number assigned to the output file.
                                               !! Needed for subsuquent calls to
                                               !! diag_output_mod
    LOGICAL         , INTENT(in)  :: all_scalar_or_1d
    TYPE(domain2d)  , INTENT(in)  :: domain
    TYPE(diag_atttype), INTENT(in), DIMENSION(:), OPTIONAL :: attributes
    TYPE(domainUG), INTENT(in)    :: domainU !< The unstructure domain
    type(FmsNetcdfUnstructuredDomainFile_t),intent(inout),target :: fileobjU
    type(FmsNetcdfDomainFile_t),intent(inout),target :: fileobj
    type(FmsNetcdfFile_t),intent(inout),target :: fileobjND
    class(FmsNetcdfFile_t), pointer :: fileob => NULL()
    character(*),intent(out) :: fnum_domain
    INTEGER :: form, threading, fileset, i
    TYPE(diag_global_att_type) :: gAtt
    character(len=:),allocatable :: fname_no_tile
    integer :: len_file_name
    integer, allocatable, dimension(:) :: current_pelist
    integer :: mype  !< The pe you are on
    character(len=9) :: mype_string !< a string to store the pe
    !---- initialize mpp_io ----
    IF ( .NOT.module_is_initialized ) THEN
       CALL mpp_io_init ()
       module_is_initialized = .TRUE.
       CALL write_version_number("DIAG_OUTPUT_MOD", version)
    END IF
    !---- set up output file ----
    SELECT CASE (FORMAT)
    CASE (NETCDF1)
       form      = MPP_NETCDF
       threading = MPP_MULTI
       fileset   = MPP_MULTI
    CASE default
       ! <ERROR STATUS="FATAL">invalid format</ERROR>
       CALL error_mesg('diag_output_init', 'invalid format', FATAL)
    END SELECT

    IF(all_scalar_or_1d) THEN
       threading = MPP_SINGLE
       fileset   = MPP_SINGLE
    END IF

    len_file_name = len(trim(file_name))
!> If the file name has .tileX or .tileX.nc where X is a one or two digit tile number, removes
!! that suffix from the time name because fms2_io will add it
!! \note If mpp_domains accepts more than 99 tiles, this will need to be updated
    allocate(character(len=len_file_name) :: fname_no_tile)
    if (len_file_name < 6) then
       if (trim(file_name) == "tile") then
          call error_mesg('diag_output_init', 'You can not name your history file "tile"',FATAL)
       else
          fname_no_tile = trim(file_name)
       endif
    !> One-digit tile numbers example
    !! \verbatim
    !! filename.tile1.nc
    !!       09876543210
    !!          ^  ^
    !! filename.tile1
    !!    09876543210
    !!          ^  ^
    !! \endverbatim
    elseif (lowercase(file_name(len_file_name-4:len_file_name-1)) .eq. "tile") then
       fname_no_tile = file_name(1:len_file_name-6)
    elseif (len_file_name < 9) then
       fname_no_tile = trim(file_name)
    elseif (lowercase(file_name(len_file_name-7:len_file_name-4)) .eq. "tile") then
       fname_no_tile = file_name(1:len_file_name-9)
    !> Two-digit tile numbers example
    !! \verbatim
    !! filename.tile10.nc
    !!        09876543210
    !!          ^  ^
    !! filename.tile10
    !!     09876543210
    !!          ^  ^
    !! \endverbatim
    elseif (lowercase(file_name(len_file_name-5:len_file_name-2)) .eq. "tile") then
       fname_no_tile = file_name(1:len_file_name-7)

    elseif (lowercase(file_name(len_file_name-5:len_file_name-8)) .eq. "tile") then
       fname_no_tile = file_name(1:len_file_name-10)
    else
       fname_no_tile = trim(file_name)
    endif
!> If there is a .nc suffix on the file name, removes the .nc
    if (len(trim(fname_no_tile)) > 3 ) then
       checkNC: do i = 3,len(trim(fname_no_tile))
         if (fname_no_tile(i-2:i) == ".nc") then
            fname_no_tile(i-2:i) = "   "
            exit checkNC
         endif
       enddo checkNC
    endif

!> Checks to make sure that only domain2D or domainUG is used.  If both are not null, then FATAL
    if (domain .NE. NULL_DOMAIN2D .AND. domainU .NE. NULL_DOMAINUG)&
          & CALL error_mesg('diag_output_init', "Domain2D and DomainUG can not be used at the same time in "//&
          & trim(file_name), FATAL)

    !---- open output file (return file_unit id) -----
    IF ( domain .NE. NULL_DOMAIN2D ) THEN
     !> Check if there is an io_domain
     iF ( associated(mpp_get_io_domain(domain)) ) then
       fileob => fileobj
       if (.not.check_if_open(fileob)) call open_check(open_file(fileobj, trim(fname_no_tile)//".nc", "overwrite", &
                            domain, is_restart=.false.))
       fnum_domain = "2d" ! 2d domain
       file_unit = 2
     elSE !< No io domain, so every core is going to write its own file.
       fileob => fileobjND
       mype = mpp_pe()
       write(mype_string,'(I0.4)') mype
        if (.not.check_if_open(fileob)) then
               call open_check(open_file(fileobjND, trim(fname_no_tile)//".nc."//trim(mype_string), "overwrite", &
                            is_restart=.false.))
               !< For regional subaxis add the NumFilesInSet attribute, which is added by fms2_io for (other)
               !< domains with sufficient decomposition info. Note mppnccombine will work with an entry of zero.
               call register_global_attribute(fileobjND, "NumFilesInSet", 0)
       endif
       fnum_domain = "nd" ! no domain
       if (file_unit < 0) file_unit = 10
     endiF
    ELSE IF (domainU .NE. NULL_DOMAINUG) THEN
       fileob => fileobjU
       if (.not.check_if_open(fileob)) call open_check(open_file(fileobjU, trim(fname_no_tile)//".nc", "overwrite", &
                            domainU, is_restart=.false.))
       fnum_domain = "ug" ! unstructured grid
       file_unit=3
    ELSE
       fileob => fileobjND
!        if (.not.check_if_open(fileob) .and. mpp_pe() == mpp_root_pe()) then
        allocate(current_pelist(mpp_npes()))
        call mpp_get_current_pelist(current_pelist)
        if (.not.check_if_open(fileob)) then
               call open_check(open_file(fileobjND, trim(fname_no_tile)//".nc", "overwrite", &
                            pelist=current_pelist, is_restart=.false.))
        endif
       fnum_domain = "nd" ! no domain
       if (file_unit < 0) file_unit = 10
       deallocate(current_pelist)
    END IF

    !---- write global attributes ----
    IF ( file_title(1:1) /= ' ' ) THEN
       call register_global_attribute(fileob, 'title', TRIM(file_title), str_len=len_trim(file_title))
    END IF

    IF ( PRESENT(attributes) ) THEN
       DO i=1, SIZE(attributes)
          SELECT CASE (attributes(i)%type)
          CASE (NF90_INT)
             call register_global_attribute(fileob, TRIM(attributes(i)%name), attributes(i)%iatt)
          CASE (NF90_FLOAT)

             call register_global_attribute(fileob, TRIM(attributes(i)%name), attributes(i)%fatt)
          CASE (NF90_CHAR)

             call register_global_attribute(fileob, TRIM(attributes(i)%name), attributes(i)%catt, str_len=len_trim(attributes(i)%catt))
          CASE default
             ! <ERROR STATUS="FATAL">
             !   Unknown attribute type for attribute <name> to module/input_field <module_name>/<field_name>.
             !   Contact the developers.
             ! </ERROR>
             CALL error_mesg('diag_output_mod::diag_output_init', 'Unknown attribute type for global attribute "'&
                  &//TRIM(attributes(i)%name)//'" in file "'//TRIM(file_name)//'". Contact the developers.', FATAL)
          END SELECT
       END DO
    END IF
    !---- write grid type (mosaic or regular)
    CALL get_diag_global_att(gAtt)

    call register_global_attribute(fileob, 'grid_type', TRIM(gAtt%grid_type), str_len=len_trim(gAtt%grid_type))

    call register_global_attribute(fileob, 'grid_tile', TRIM(gAtt%tile_name), str_len=len_trim(gAtt%tile_name))

  END SUBROUTINE diag_output_init_fms2_io

  !> @brief Write the axis meta data to file.
  SUBROUTINE write_axis_meta_data_fms2_io(file_unit, axes, fileob, time_ops, time_axis_registered)
    INTEGER, INTENT(in) :: file_unit !< File unit number
    INTEGER, INTENT(in) :: axes(:) !< Array of axis ID's, including the time axis
    class(FmsNetcdfFile_t) , intent(inout),target :: fileob
    class(FmsNetcdfFile_t) ,pointer                        :: fptr
    LOGICAL, INTENT(in), OPTIONAL :: time_ops !< .TRUE. if this file contains any min, max, time_rms, or time_average
    logical, intent(inout) , optional :: time_axis_registered
    TYPE(domain1d)       :: Domain

    TYPE(domainUG)       :: domainU

    CHARACTER(len=mxch)  :: axis_name, axis_units, axis_name_current
    CHARACTER(len=mxchl) :: axis_long_name
    CHARACTER(len=1)     :: axis_cart_name
    INTEGER              :: axis_direction, axis_edges
    REAL, ALLOCATABLE    :: axis_data(:)
    INTEGER, ALLOCATABLE :: axis_extent(:), pelist(:)
integer :: domain_size, axis_length, axis_pos
    INTEGER              :: num_attributes
    TYPE(diag_atttype), DIMENSION(:), ALLOCATABLE :: attributes
    INTEGER              :: calendar, id_axis, id_time_axis
    INTEGER              :: i, j, index, num, length, edges_index
    INTEGER              :: gbegin, gend, gsize, ndivs
    LOGICAL              :: time_ops1
    CHARACTER(len=2048)  :: err_msg
    type(domainUG),pointer                     :: io_domain
    integer(I4_KIND)                          :: io_domain_npes
    integer(I4_KIND),dimension(:),allocatable :: io_pelist
    integer(I4_KIND),dimension(:),allocatable :: unstruct_axis_sizes
    real,dimension(:),allocatable              :: unstruct_axis_data
    integer                                    :: id_axis_current
    logical :: is_time_axis_registered
    integer :: istart, iend
    integer :: gstart, cstart, cend !< Start and end of global and compute domains
    integer :: clength !< Length of compute domain
    integer :: data_size
    integer, allocatable, dimension(:) :: all_indicies
    character(len=32) :: type_str !< Str indicating the type of the axis data

    ! Make sure err_msg is initialized
    err_msg = ''
    fptr => fileob !Use for selecting a type
    IF ( PRESENT(time_ops) ) THEN
       time_ops1 = time_ops
    ELSE
       time_ops1 = .FALSE.
    END IF
    if (present(time_axis_registered)) then
     is_time_axis_registered = time_axis_registered
    else
     is_time_axis_registered = .false.
    endif
    !---- save the current file_unit ----
    IF ( num_axis_in_file == 0 ) current_file_unit = file_unit

    !---- dummy checks ----
    num = SIZE(axes(:))
    ! <ERROR STATUS="FATAL">number of axes < 1 </ERROR>
    IF ( num < 1 ) CALL error_mesg('write_axis_meta_data', 'number of axes < 1.', FATAL)

    ! <ERROR STATUS="FATAL">writing meta data out-of-order to different files.</ERROR>
    IF ( file_unit /= current_file_unit ) CALL error_mesg('write_axis_meta_data',&
         & 'writing meta data out-of-order to different files.', FATAL)

    IF (pack_size .eq. 1) then
       type_str = "double"
    ELSE IF (pack_size .eq. 2) then
       type_str = "float"
    ENDIF

    !---- check all axes ----
    !---- write axis meta data for new axes ----
    DO i = 1, num
       id_axis = axes(i)
       index = get_axis_index ( id_axis )

       !---- skip axes already written -----
       IF ( index > 0 ) CYCLE

       !---- create new axistype (then point to) -----
       num_axis_in_file = num_axis_in_file + 1
       axis_in_file(num_axis_in_file) = id_axis
       edge_axis_flag(num_axis_in_file) = .FALSE.
       length = get_axis_global_length(id_axis)
       ALLOCATE(axis_data(length))

       CALL get_diag_axis(id_axis, axis_name, axis_units, axis_long_name,&
            & axis_cart_name, axis_direction, axis_edges, Domain, DomainU, axis_data,&
            & num_attributes, attributes, domain_position=axis_pos)

       IF ( Domain .NE. null_domain1d ) THEN
          IF ( length > 0 ) THEN
             if (trim(uppercase(trim(axis_cart_name))) .eq. "X" .or. trim(uppercase(trim(axis_cart_name))) .eq. "Y") then
                  select type (fptr)
                    type is (FmsNetcdfDomainFile_t)
                         call register_axis(fptr, axis_name, lowercase(trim(axis_cart_name)), domain_position=axis_pos )
                      if (allocated(fptr%pelist)) then
                         call get_global_io_domain_indices(fptr, trim(axis_name), istart, iend)
                         call register_field(fptr, axis_name, type_str, (/axis_name/) )
                         if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                         call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                         call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                         select case (axis_direction)
                              case (1)
                                   call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                              case (-1)
                                   call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                         end select
                         call write_data(fptr, axis_name, axis_data(istart:iend) )
                      endif
                    type is (FmsNetcdfFile_t) !< For regional X and Y axes, treat as any other axis
                         call mpp_get_global_domain(domain, begin=gstart, end=gend)  !< Get the global indicies
                         call mpp_get_compute_domain(domain, begin=cstart, end=cend, size=clength) !< Get the compute indicies
                         iend =  cend - gstart + 1     !< Get the array indicies for the axis data
                         istart = cstart - gstart + 1
                         call register_axis(fptr, axis_name, dimension_length=clength)
                         call register_field(fptr, axis_name, type_str, (/axis_name/) )
                         call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                         call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                         call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                         select case (axis_direction)
                              case (1)
                                   call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                              case (-1)
                                   call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                         end select
                         !< For regional subaxis add the "domain_decomposition" attribute, which is added
                         !< fms2_io for (other) domains with sufficient decomposition info.
                         call register_variable_attribute(fptr, axis_name, "domain_decomposition", &
                              (/gstart, gend, cstart, cend/))
                         call write_data(fptr, axis_name, axis_data(istart:iend) )
                    class default
                         call error_mesg("diag_output_mod::write_axis_meta_data", &
                              "The file object is not the right type. It must be FmsNetcdfDomainFile_t or "//&
                                "FmsNetcdfFile_t for a X or Y axis, ", FATAL)
                  end select
             endif

          ELSE
               select type (fptr)
                    type is (FmsNetcdfDomainFile_t)
                         call register_axis(fptr, axis_name, lowercase(trim(axis_cart_name)), domain_position=axis_pos )
                      if (allocated(fptr%pelist)) then
                         call get_global_io_domain_indices(fptr, trim(axis_name), istart, iend)
                         call register_field(fptr, axis_name, type_str, (/axis_name/) )
                      endif
                    type is (FmsNetcdfUnstructuredDomainFile_t)
                        call register_axis(fptr, axis_name )
                    type is (FmsNetcdfFile_t)
                         call register_axis(fptr, axis_name, dimension_length=size(axis_data))
                      if (allocated(fptr%pelist)) then
!                         call get_global_io_domain_indices(fptr, trim(axis_name), istart, iend)
                         istart = lbound(axis_data,1)
                         iend = ubound(axis_data,1)
                         call register_field(fptr, axis_name, type_str, (/axis_name/) )
                      endif
                    class default
                         call error_mesg("diag_output_mod::write_axis_meta_data", &
                              "The FmsNetcdfDomain file object is not the right type.", FATAL)
                end select
                    call register_field(fileob, axis_name, type_str, (/axis_name/) )
                    call register_variable_attribute(fileob, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                    call register_variable_attribute(fileob, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                    call register_variable_attribute(fileob, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                    select case (axis_direction)
                         case (1)
                              call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                         case (-1)
                              call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                    end select
                    call write_data(fileob, axis_name, axis_data(istart:iend) )
          END IF
       ELSE
          IF ( length > 0 ) THEN

            !For an unstructured dimension, only the root rank of the io_domain
            !pelist will perform the wirte, so a gather of the unstructured
            !axis size and axis data is required.
             if (uppercase(trim(axis_cart_name)) .eq. "U") then
                 if (DomainU .eq. null_domainUG) then
                     call error_mesg("diag_output_mod::write_axis_meta_data", &
                                     "A non-nul domainUG is required to" &
                                     //" write unstructured axis metadata.", &
                                     FATAL)
                 endif
                 io_domain => null()
                 io_domain => mpp_get_UG_io_domain(DomainU)
                 io_domain_npes = mpp_get_UG_domain_npes(io_domain)
                 allocate(io_pelist(io_domain_npes))
                 call mpp_get_UG_domain_pelist(io_domain, &
                                               io_pelist)
                 allocate(unstruct_axis_sizes(io_domain_npes))
                 unstruct_axis_sizes = 0
                 call mpp_gather((/size(axis_data)/), &
                                 unstruct_axis_sizes, &
                                 io_pelist)
                 if (mpp_pe() .eq. io_pelist(1)) then
                     allocate(unstruct_axis_data(sum(unstruct_axis_sizes)))
                 else
                     allocate(unstruct_axis_data(1))
                 endif
                 unstruct_axis_data = 0.0
                 call mpp_gather(axis_data, &
                                 size(axis_data), &
                                 unstruct_axis_data, &
                                 unstruct_axis_sizes, &
                                 io_pelist)
                  select type (fptr)
                   type is (FmsNetcdfUnstructuredDomainFile_t)
                        call register_axis(fptr, axis_name )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        call write_data(fptr, axis_name, axis_data)
                   class default
                        call error_mesg("diag_output_mod::write_axis_meta_data", &
                             "The file unstructred 1 object is not the right type.", NOTE)
                  end select
                 deallocate(io_pelist)
                 deallocate(unstruct_axis_sizes)
                 deallocate(unstruct_axis_data)
                 io_domain => null()

             else
                 select type (fptr)
                   type is (FmsNetcdfUnstructuredDomainFile_t)
                        call register_axis(fptr, axis_name, size(axis_data) )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        select case (axis_direction)
                             case (1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                             case (-1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                        end select
                        call write_data(fptr, axis_name, axis_data)
                   type is (FmsNetcdfDomainFile_t)
                    if (.not.variable_exists(fptr, axis_name)) then
                        call register_axis(fptr, axis_name, size(axis_data) )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        select case (axis_direction)
                             case (1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                             case (-1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                        end select
                        call write_data(fptr, axis_name, axis_data)
                    endif
                   type is (FmsNetcdfFile_t)
                    if (.not.variable_exists(fptr, axis_name)) then
                        call register_axis(fptr, axis_name, size(axis_data) )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        select case (axis_direction)
                             case (1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                             case (-1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                        end select
                        call write_data(fptr, axis_name, axis_data)
                    endif
                   class default
                        call error_mesg("diag_output_mod::write_axis_meta_data", &
                             "The file object is not the right type.", FATAL)
                 end select
             endif

          ELSE
!> @note Check if the time variable is registered.  It's possible that is_time_axis_registered is set to true if using
!! time-templated files because they aren't closed when done writing.  An alternative to this set up would be to put
!! variable_exists into the if statement with an .or. so that it gets registered.
                is_time_axis_registered = variable_exists(fptr,trim(axis_name),.true.)
                if (allocated(fptr%pelist) .and. .not. is_time_axis_registered) then
                 select type (fptr)
                   type is (FmsNetcdfDomainFile_t)
                        call register_axis(fptr, trim(axis_name), unlimited )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))

                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        is_time_axis_registered = .true.
                        if (present(time_axis_registered)) time_axis_registered = is_time_axis_registered
                   type is (FmsNetcdfUnstructuredDomainFile_t)
                        call register_axis(fptr, axis_name, size(axis_data) )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        is_time_axis_registered = .true.
                   type is (FmsNetcdfFile_t)
                        call register_axis(fptr, trim(axis_name), unlimited )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        is_time_axis_registered = .true.
                        if (present(time_axis_registered)) time_axis_registered = is_time_axis_registered
                   class default
                        call error_mesg("diag_output_mod::write_axis_meta_data", &
                             "The file object is not the right type.", FATAL)
                 end select
                endif
          END IF
       END IF

       ! Write axis attributes
       id_axis = mpp_get_id(Axis_types(num_axis_in_file))
       CALL write_attribute_meta(file_unit, id_axis, num_attributes, attributes, err_msg, varname=axis_name, fileob=fileob)
       IF ( LEN_TRIM(err_msg) .GT. 0 ) THEN
          CALL error_mesg('diag_output_mod::write_axis_meta_data', TRIM(err_msg), FATAL)
       END IF

       !---- write additional attribute (calendar_type) for time axis ----
       !---- NOTE: calendar attribute is compliant with CF convention
       !---- http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-current.htm#cal
       IF ( axis_cart_name == 'T' ) THEN
          time_axis_flag(num_axis_in_file) = .TRUE.
          id_time_axis = mpp_get_id(Axis_types(num_axis_in_file))
          calendar = get_calendar_type()


          call register_variable_attribute(fileob, axis_name, "calendar_type", &
                                    UPPERCASE(TRIM(valid_calendar_types(calendar))), str_len=len_trim(valid_calendar_types(calendar)) )
          call register_variable_attribute(fileob, axis_name, "calendar", &
                                    lowercase(TRIM(valid_calendar_types(calendar))), str_len=len_trim(valid_calendar_types(calendar)) )
          IF ( time_ops1 ) THEN

             call register_variable_attribute(fileob, axis_name, 'bounds', TRIM(axis_name)//'_bnds', str_len=len_trim(TRIM(axis_name)//'_bnds'))
          END IF
          call set_fileobj_time_name(fileob, axis_name)
       ELSE
          time_axis_flag(num_axis_in_file) = .FALSE.
       END IF

       DEALLOCATE(axis_data)

       ! Deallocate attributes
       IF ( ALLOCATED(attributes) ) THEN
          DO j=1, num_attributes
             IF ( allocated(attributes(j)%fatt ) ) THEN
                DEALLOCATE(attributes(j)%fatt)
             END IF
             IF ( allocated(attributes(j)%iatt ) ) THEN
                DEALLOCATE(attributes(j)%iatt)
             END IF
          END DO
          DEALLOCATE(attributes)
       END IF

       !------------- write axis containing edge information ---------------

       !  --- this axis has no edges -----
       IF ( axis_edges <= 0 ) CYCLE

       !  --- was this axis edge previously defined? ---
       id_axis_current = id_axis
       axis_name_current = axis_name
       id_axis = axis_edges
       edges_index = get_axis_index(id_axis)
       IF ( edges_index > 0 ) CYCLE

       !  ---- get data for axis edges ----
       length = get_axis_global_length ( id_axis )
       ALLOCATE(axis_data(length))
       CALL get_diag_axis(id_axis, axis_name, axis_units, axis_long_name, axis_cart_name,&
            & axis_direction, axis_edges, Domain, DomainU, axis_data, num_attributes, attributes)

       !  ---- write edges attribute to original axis ----
       call register_variable_attribute(fileob, axis_name_current, "edges",trim(axis_name), str_len=len_trim(axis_name))
       !  ---- add edges index to axis list ----
       !  ---- assume this is not a time axis ----
       num_axis_in_file = num_axis_in_file + 1
       axis_in_file(num_axis_in_file) = id_axis
       edge_axis_flag(num_axis_in_file) = .TRUE.
       time_axis_flag (num_axis_in_file) = .FALSE.

       !  ---- write edges axis to file ----
       IF ( Domain .NE. null_domain1d ) THEN
          ! assume domain decomposition is irregular and loop through all prev and next
          ! domain pointers extracting domain extents.  Assume all pes are used in
          ! decomposition
          CALL mpp_get_global_domain(Domain, begin=gbegin, END=gend, size=gsize)
          CALL mpp_get_layout(Domain, ndivs)
          IF ( ndivs .NE. 1 ) THEN
             IF ( ALLOCATED(axis_extent) ) DEALLOCATE(axis_extent)
             ALLOCATE(axis_extent(0:ndivs-1))
             CALL mpp_get_compute_domains(Domain,size=axis_extent(0:ndivs-1))
             gend=gend+1
             axis_extent(ndivs-1)= axis_extent(ndivs-1)+1
             IF ( ALLOCATED(pelist) ) DEALLOCATE(pelist)
             ALLOCATE(pelist(0:ndivs-1))
             CALL mpp_get_pelist(Domain,pelist)
          END IF
       END IF

!> Add edges axis with fms2_io
                 select type (fptr)
                   type is (FmsNetcdfUnstructuredDomainFile_t)
                        call register_axis(fptr, axis_name, size(axis_data) )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        select case (axis_direction)
                             case (1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                             case (-1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                        end select
                        call write_data(fptr, axis_name, axis_data)
                   type is (FmsNetcdfDomainFile_t)
                    if (.not.variable_exists(fptr, axis_name)) then
                        call register_axis(fptr, axis_name, size(axis_data) )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        select case (axis_direction)
                             case (1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                             case (-1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                        end select
                        call write_data(fptr, axis_name, axis_data)
                    endif
                   type is (FmsNetcdfFile_t)
                    if (.not.variable_exists(fptr, axis_name)) then
                        call register_axis(fptr, axis_name, size(axis_data) )
                        call register_field(fptr, axis_name, type_str, (/axis_name/) )
                        if(trim(axis_units) .ne. "none") call register_variable_attribute(fptr, axis_name, "units", trim(axis_units), str_len=len_trim(axis_units))
                        call register_variable_attribute(fptr, axis_name, "long_name", trim(axis_long_name), str_len=len_trim(axis_long_name))
                        if(trim(axis_cart_name).ne."N") call register_variable_attribute(fptr, axis_name, "axis",trim(axis_cart_name), str_len=len_trim(axis_cart_name))
                        select case (axis_direction)
                             case (1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "up", str_len=len_trim("up"))
                             case (-1)
                                  call register_variable_attribute(fptr, axis_name, "positive", "down", str_len=len_trim("down"))
                        end select
                        call write_data(fptr, axis_name, axis_data)
                    endif
                   class default
                        call error_mesg("diag_output_mod::write_axis_meta_data", &
                             "The file object unstructured 2 is not the right type.", FATAL)
                 end select
       ! Write edge axis attributes
       id_axis = mpp_get_id(Axis_types(num_axis_in_file))
!       CALL write_attribute_meta(file_unit, id_axis, num_attributes, attributes, err_msg)
       IF ( LEN_TRIM(err_msg) .GT. 0 ) THEN
          CALL error_mesg('diag_output_mod::write_axis_meta_data', TRIM(err_msg), FATAL)
       END IF

       DEALLOCATE (axis_data)
       ! Deallocate attributes
       IF ( ALLOCATED(attributes) ) THEN
          DO j=1, num_attributes
             IF ( allocated(attributes(j)%fatt ) ) THEN
                DEALLOCATE(attributes(j)%fatt)
             END IF
             IF ( allocated(attributes(j)%iatt ) ) THEN
                DEALLOCATE(attributes(j)%iatt)
             END IF
          END DO
          DEALLOCATE(attributes)
       END IF
    END DO
  END SUBROUTINE write_axis_meta_data_fms2_io

  !> @brief Write the field meta data to file.
  !! @return diag_fieldtype Field
  !! @details The meta data for the field is written to the file indicated by file_unit
  FUNCTION write_field_meta_data_fms2_io ( file_unit, name, axes, units, long_name, range, pack, mval,&
       & avg_name, time_method, standard_name, interp_method, attributes, num_attributes,     &
       & use_UGdomain, fileob) result ( Field )
    INTEGER, INTENT(in) :: file_unit !< Output file unit number
    INTEGER, INTENT(in) :: axes(:) !< Array of axis IDs
    CHARACTER(len=*), INTENT(in) :: name !< Field name
    CHARACTER(len=*), INTENT(in) :: units !< Field units
    CHARACTER(len=*), INTENT(in) :: long_name !< Field's long name
    REAL, OPTIONAL, INTENT(in) :: RANGE(2) !< Valid range (min, max).  If min > max, the range will be ignored
    REAL, OPTIONAL, INTENT(in) :: mval !< Missing value, must be within valid range
    INTEGER, OPTIONAL, INTENT(in) :: pack !< Packing flag.  Only valid when range specified.  Valid values:
                                          !! Flag | Size
                                          !! --- | ---
                                          !! 1 | 64bit
                                          !! 2 | 32bit
                                          !! 4 | 16bit
                                          !! 8 | 8bit
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: avg_name !< Name of variable containing time averaging info
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: time_method !< Name of transformation applied to the time-varying data,
                                                          !! i.e. "avg", "min", "max"
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: standard_name !< Standard name of field
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: interp_method
    TYPE(diag_atttype), DIMENSION(:), allocatable, OPTIONAL, INTENT(in) :: attributes
    INTEGER, OPTIONAL, INTENT(in) :: num_attributes
    LOGICAL, OPTIONAL, INTENT(in) :: use_UGdomain
class(FmsNetcdfFile_t), intent(inout)     :: fileob

    logical :: is_time_bounds !< Flag indicating if the variable is time_bounds
    CHARACTER(len=256) :: standard_name2
    CHARACTER(len=1280) :: att_str
    TYPE(diag_fieldtype) :: Field
    LOGICAL :: coord_present
    CHARACTER(len=40) :: aux_axes(SIZE(axes))
    CHARACTER(len=160) :: coord_att
    CHARACTER(len=1024) :: err_msg

character(len=128),dimension(size(axes)) :: axis_names
    REAL :: scale, add
    INTEGER :: i, indexx, num, ipack, np, att_len
    LOGICAL :: use_range
    INTEGER :: axis_indices(SIZE(axes))
    logical :: use_UGdomain_local
    !---- Initialize err_msg to bank ----
    err_msg = ''

    !---- dummy checks ----
    coord_present = .FALSE.
    IF( PRESENT(standard_name) ) THEN
       standard_name2 = standard_name
    ELSE
       standard_name2 = 'none'
    END IF

    use_UGdomain_local = .false.
    if(present(use_UGdomain)) use_UGdomain_local = use_UGdomain

    num = SIZE(axes(:))
    ! <ERROR STATUS="FATAL">number of axes < 1</ERROR>
    IF ( num < 1 ) CALL error_mesg ( 'write_meta_data', 'number of axes < 1', FATAL)
    ! <ERROR STATUS="FATAL">writing meta data out-of-order to different files</ERROR>
    IF ( file_unit /= current_file_unit ) CALL error_mesg ( 'write_meta_data',  &
         & 'writing meta data out-of-order to different files', FATAL)

    IF (trim(name) .eq. "time_bnds") then
       is_time_bounds = .true.
    ELSE
       is_time_bounds = .false.
    ENDIF

    !---- check all axes for this field ----
    !---- set up indexing to axistypes ----
    DO i = 1, num
       indexx = get_axis_index(axes(i))
       !---- point to existing axistype -----
       IF ( indexx > 0 ) THEN
          axis_indices(i) = indexx
       ELSE
          ! <ERROR STATUS="FATAL">axis data not written for field</ERROR>
          CALL error_mesg ('write_field_meta_data',&
               & 'axis data not written for field '//TRIM(name), FATAL)
       END IF
       !Get the axes names
          call get_diag_axis_name(axes(i),axis_names(i))
    END DO

    !  Create coordinate attribute
    IF ( num >= 2 .OR. (num==1 .and. use_UGdomain_local) ) THEN
       coord_att = ' '
       DO i = 1, num
          aux_axes(i) = get_axis_aux(axes(i))
          IF( TRIM(aux_axes(i)) /= 'none' ) THEN
             IF(LEN_TRIM(coord_att) == 0) THEN
                coord_att = TRIM(aux_axes(i))
             ELSE
                coord_att = TRIM(coord_att)// ' '//TRIM(aux_axes(i))
             ENDIF
             coord_present = .TRUE.
          END IF
       END DO
    END IF

    !--------------------- write field meta data ---------------------------

    !---- select packing? ----
    !(packing option only valid with range option)
    IF ( PRESENT(pack) ) THEN
       ipack = pack
    ELSE
       ipack = 2
    END IF

    !---- check range ----
    use_range = .FALSE.
    add = 0.0
    scale = 1.0
    IF ( PRESENT(range) ) THEN
       IF ( RANGE(2) > RANGE(1) ) THEN
          use_range = .TRUE.
          !---- set packing parameters ----
          IF ( ipack > 2 ) THEN
             np = ipack/4
             add = 0.5*(RANGE(1)+RANGE(2))
             scale = (RANGE(2)-RANGE(1)) / real(max_range(2,np)-max_range(1,np))
          END IF
       END IF
    END IF

    !---- select packing? ----
    IF ( PRESENT(mval) ) THEN
       Field%miss = mval
       Field%miss_present = .TRUE.
       IF ( ipack > 2 ) THEN
          np = ipack/4
          Field%miss_pack = REAL(missval(np))*scale+add
          Field%miss_pack_present = .TRUE.
       ELSE
          Field%miss_pack = mval
          Field%miss_pack_present = .FALSE.
       END IF
    ELSE
       Field%miss_present = .FALSE.
       Field%miss_pack_present = .FALSE.
    END IF

    !------ write meta data and return fieldtype -------
!!! Fill in mpp fieldtype for field%field
    IF ( use_range ) THEN
       IF ( Field%miss_present ) THEN
          CALL fillin_fieldtype( Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units, long_name,&
               & RANGE(1), RANGE(2),&
               & missing=Field%miss_pack,&
               & fill=Field%miss_pack,&
               & scale=scale, add=add, pack=ipack,&
               & time_method=time_method)
       ELSE
          CALL fillin_fieldtype( Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units,  long_name,&
               & RANGE(1), RANGE(2),&
               & missing=CMOR_MISSING_VALUE,&
               & fill=CMOR_MISSING_VALUE,&
               & scale=scale, add=add, pack=ipack,&
               & time_method=time_method)
       END IF
    ELSE
       IF ( Field%miss_present ) THEN
          CALL fillin_fieldtype( Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units, long_name,&
               & missing=Field%miss_pack,&
               & fill=Field%miss_pack,&
               & pack=ipack, time_method=time_method)
       ELSE
          CALL fillin_fieldtype( Field%Field,&
               & Axis_types(axis_indices(1:num)),&
               & name, units, long_name,&
               & missing=CMOR_MISSING_VALUE,&
               & fill=CMOR_MISSING_VALUE,&
               & pack=ipack, time_method=time_method)
       END IF
    END IF
  if (.not. variable_exists(fileob,name)) then
  ! ipack Valid values:
  !        1 = 64bit </LI>
  !        2 = 32bit </LI>
  !        4 = 16bit </LI>
  !        8 =  8bit </LI>
     select case (ipack)
     case (1)
          call register_field(fileob,name,"double",axis_names)
          !< Don't write the _FillValue, missing_value if the variable is
          !time_bounds to be cf compliant
          if (.not. is_time_bounds) then
          IF ( Field%miss_present ) THEN
               call register_variable_attribute(fileob,name,"_FillValue",real(Field%miss_pack,8))
               call register_variable_attribute(fileob,name,"missing_value",real(Field%miss_pack,8))
          ELSE
               call register_variable_attribute(fileob,name,"_FillValue",real(CMOR_MISSING_VALUE,8))
               call register_variable_attribute(fileob,name,"missing_value",real(CMOR_MISSING_VALUE,8))
          ENDIF
          IF ( use_range ) then
               call register_variable_attribute(fileob,name,"valid_range", real(RANGE,8))
          ENDIF
          endif !< if (.not. is_time_bounds)
     case (2) !default
          call register_field(fileob,name,"float",axis_names)
          !< Don't write the _FillValue, missing_value if the variable is
          !time_bounds to be cf compliant
          if (.not. is_time_bounds) then
          IF ( Field%miss_present ) THEN
               call register_variable_attribute(fileob,name,"_FillValue",real(Field%miss_pack,4))
               call register_variable_attribute(fileob,name,"missing_value",real(Field%miss_pack,4))
          ELSE
               call register_variable_attribute(fileob,name,"_FillValue",real(CMOR_MISSING_VALUE,4))
               call register_variable_attribute(fileob,name,"missing_value",real(CMOR_MISSING_VALUE,4))
          ENDIF
          IF ( use_range ) then
               call register_variable_attribute(fileob,name,"valid_range", real(RANGE,4))
          ENDIF
          endif !< if (.not. is_time_bounds)
     case default
          CALL error_mesg('diag_output_mod::write_field_meta_data',&
               &"Pack values must be 1 or 2. Contact the developers.", FATAL)
     end select
     if (trim(units) .ne. "none") call register_variable_attribute(fileob,name,"units",trim(units), str_len=len_trim(units))
     call register_variable_attribute(fileob,name,"long_name",long_name, str_len=len_trim(long_name))
     IF (present(time_method) ) then
          call register_variable_attribute(fileob,name,'cell_methods','time: '//trim(time_method), str_len=len_trim('time: '//trim(time_method)))
     ENDIF
  endif
    !---- write user defined attributes -----
    IF ( PRESENT(num_attributes) ) THEN
       IF ( PRESENT(attributes) ) THEN
          IF ( num_attributes .GT. 0 .AND. allocated(attributes) ) THEN
             CALL write_attribute_meta(file_unit, mpp_get_id(Field%Field), num_attributes, attributes, time_method, err_msg, fileob=fileob, varname=name)
             IF ( LEN_TRIM(err_msg) .GT. 0 ) THEN
                CALL error_mesg('diag_output_mod::write_field_meta_data',&
                     & TRIM(err_msg)//" Contact the developers.", FATAL)
             END IF
          ELSE
             ! Catch some bad cases
             IF ( num_attributes .GT. 0 .AND. .NOT.allocated(attributes) ) THEN
                CALL error_mesg('diag_output_mod::write_field_meta_data',&
                     & 'num_attributes > 0 but attributes is not allocated for attribute '&
                     &//TRIM(attributes(i)%name)//' for field '//TRIM(name)//'. Contact the developers.', FATAL)
             ELSE IF ( num_attributes .EQ. 0 .AND. allocated(attributes) ) THEN
                CALL error_mesg('diag_output_mod::write_field_meta_data',&
                     & 'num_attributes == 0 but attributes is allocated for attribute '&
                     &//TRIM(attributes(i)%name)//' for field '//TRIM(name)//'. Contact the developers.', FATAL)
             END IF
          END IF
       ELSE
          ! More edge error cases
          CALL error_mesg('diag_output_mod::write_field_meta_data',&
               & 'num_attributes present but attributes missing for attribute '&
               &//TRIM(attributes(i)%name)//' for field '//TRIM(name)//'. Contact the developers.', FATAL)
       END IF
    ELSE IF ( PRESENT(attributes) ) THEN
       CALL error_mesg('diag_output_mod::write_field_meta_data',&
            & 'attributes present but num_attributes missing for attribute '&
            &//TRIM(attributes(i)%name)//' for field '//TRIM(name)//'. Contact the developers.', FATAL)
    END IF

    !---- write additional attribute for time averaging -----
    IF ( PRESENT(avg_name) ) THEN
       IF ( avg_name(1:1) /= ' ' ) THEN
          call register_variable_attribute(fileob,name,'time_avg_info',&
             & trim(avg_name)//'_T1,'//trim(avg_name)//'_T2,'//trim(avg_name)//'_DT', &
             & str_len=len_trim(trim(avg_name)//'_T1,'//trim(avg_name)//'_T2,'//trim(avg_name)//'_DT'))
       END IF
    END IF

    ! write coordinates attribute for CF compliance
    IF ( coord_present ) then
         call register_variable_attribute(fileob,name,'coordinates',TRIM(coord_att), str_len=len_trim(coord_att))
    ENDIF
    IF ( TRIM(standard_name2) /= 'none' ) then
         call register_variable_attribute(fileob,name,'standard_name',TRIM(standard_name2), str_len=len_trim(standard_name2))
    ENDIF
    !---- write attribute for interp_method ----
    IF( PRESENT(interp_method) ) THEN
       call register_variable_attribute(fileob,name,'interp_method', TRIM(interp_method), str_len=len_trim(interp_method))
    END IF

    !---- get axis domain ----
    Field%Domain = get_domain2d ( axes )
    Field%tile_count = get_tile_count ( axes )
    Field%DomainU = get_domainUG ( axes(1) )

  END FUNCTION write_field_meta_data_fms2_io

  !> \brief Write out attribute meta data to file
  !!
  !! Write out the attribute meta data to file, for field and axes
  SUBROUTINE write_attribute_meta_fms2_io(file_unit, id, num_attributes, attributes, time_method, err_msg, varname, fileob)
    INTEGER, INTENT(in) :: file_unit !< File unit number
    INTEGER, INTENT(in) :: id !< ID of field, file, axis to get attribute meta data
    INTEGER, INTENT(in) :: num_attributes !< Number of attributes to write
    TYPE(diag_atttype), DIMENSION(:), INTENT(in) :: attributes !< Array of attributes
    CHARACTER(len=*), INTENT(in), OPTIONAL :: time_method !< To include in cell_methods attribute if present
    CHARACTER(len=*), INTENT(out), OPTIONAL :: err_msg !< Return error message
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: varname !< The name of the variable
class(FmsNetcdfFile_t), intent(inout)     :: fileob

    INTEGER :: i, att_len
    CHARACTER(len=1280) :: att_str

    ! Clear err_msg if present
    IF ( PRESENT(err_msg) ) err_msg = ''

    DO i = 1, num_attributes
       SELECT CASE (attributes(i)%type)
       CASE (NF90_INT)
          IF ( .NOT.allocated(attributes(i)%iatt) ) THEN
             IF ( fms_error_handler('diag_output_mod::write_attribute_meta',&
                  & 'Integer attribute type indicated, but array not allocated for attribute '&
                  &//TRIM(attributes(i)%name)//'.', err_msg) ) THEN
                RETURN
             END IF
          END IF
          if (present(varname))call register_variable_attribute(fileob, varname,TRIM(attributes(i)%name)  , attributes(i)%iatt)
       CASE (NF90_FLOAT)
          IF ( .NOT.allocated(attributes(i)%fatt) ) THEN
             IF ( fms_error_handler('diag_output_mod::write_attribute_meta',&
                  & 'Real attribute type indicated, but array not allocated for attribute '&
                  &//TRIM(attributes(i)%name)//'.', err_msg) ) THEN
                RETURN
             END IF
          END IF
          if (present(varname))call register_variable_attribute(fileob, varname,TRIM(attributes(i)%name)  , real(attributes(i)%fatt,4) )
       CASE (NF90_CHAR)
          att_str = attributes(i)%catt
          att_len = attributes(i)%len
          IF ( TRIM(attributes(i)%name).EQ.'cell_methods' .AND. PRESENT(time_method) ) THEN
             ! Append ",time: time_method" if time_method present
             att_str = attributes(i)%catt(1:attributes(i)%len)//' time: '//time_method
             att_len = LEN_TRIM(att_str)
          END IF
          if (present(varname))&
               call register_variable_attribute(fileob, varname,TRIM(attributes(i)%name)  , att_str(1:att_len), str_len=att_len)

       CASE default
          IF ( fms_error_handler('diag_output_mod::write_attribute_meta', 'Invalid type for attribute '&
               &//TRIM(attributes(i)%name)//'.', err_msg) ) THEN
             RETURN
          END IF
       END SELECT
    END DO
  END SUBROUTINE write_attribute_meta_fms2_io

  !> @brief Writes axis data to file.
  !! @details Writes axis data to file.  This subroutine is to be called once per file
  !!     after all <TT>write_meta_data</TT> calls, and before the first
  !!     <TT>diag_field_out</TT> call.
  SUBROUTINE done_meta_data(file_unit)
    INTEGER,  INTENT(in)  :: file_unit !< Output file unit number

    INTEGER               :: i

    !---- write data for all non-time axes ----
    num_axis_in_file = 0
  END SUBROUTINE done_meta_data

  !> @brief Outputs the diagnostic data to a file using fms2_io taking a field object as input
  subroutine diag_field_write_field (field, buffer, static, fileob, file_num, fileobjU, fileobj, fileobjND, fnum_for_domain, time_in)
    TYPE(diag_fieldtype), INTENT(inout) :: Field !<
    REAL , INTENT(inout) :: buffer(:,:,:,:)
    logical, intent(in), optional :: static
    class(FmsNetcdfFile_t), optional, intent(inout),target :: fileob
    class(FmsNetcdfFile_t), pointer :: fptr => null()
    integer, intent(in), optional  :: file_num
    type(FmsNetcdfUnstructuredDomainFile_t),intent(inout), optional :: fileobjU(:)
    type(FmsNetcdfDomainFile_t),intent(inout), optional:: fileobj(:)
    type(FmsNetcdfFile_t),intent(inout), optional:: fileobjND(:)
    character(len=2), intent(in), optional :: fnum_for_domain
    INTEGER, OPTIONAL, INTENT(in) :: time_in
    integer :: time
    real(kind=4),allocatable :: local_buffer(:,:,:,:)
     if (present(static)) then
          if (static) time = 0
     elseif (present(time_in)) then
          time = time_in
     else
          time = 0
     endif

     if (present(fileob)) then !> Write output to the fileob file
          fptr => fileob
          select type (fptr)
          type is (FmsNetcdfFile_t)
               call write_data (fptr,trim(mpp_get_field_name(field%field)),buffer)
          type is (FmsNetcdfDomainFile_t)
               call write_data (fptr,trim(mpp_get_field_name(field%field)),buffer)
          type is (FmsNetcdfUnstructuredDomainFile_t)
               call write_data (fptr,trim(mpp_get_field_name(field%field)),buffer)
          class default
               call error_mesg("diag_field_write","fileob passed in is not one of the FmsNetcdfFile_t types",fatal)
          end select
     elseif (present(file_num) .and. present(fileobjU) .and. present(fileobjND) .and. present(fileobj) .and. present(fnum_for_domain)) then
          allocate(local_buffer(size(buffer,1),size(buffer,2),size(buffer,3),size(buffer,4)))
          local_buffer = real(buffer,4)
     !> Figure out which file object to write output to
!          if (fnum_for_domain == "2d" .or. fnum_for_domain == "nd") then
          if (fnum_for_domain == "2d" ) then
               if (check_if_open(fileobj(file_num))) then
                    if (time == 0) then
                         call write_data (fileobj (file_num), trim(mpp_get_field_name(field%field)), local_buffer)
                    else
                         call write_data (fileobj (file_num), trim(mpp_get_field_name(field%field)), local_buffer, unlim_dim_level=time)
                    endif
               endif
          elseif (fnum_for_domain == "nd") then
               if (check_if_open(fileobjND (file_num)) ) then
                    if (time == 0) then
                         call write_data (fileobjND (file_num), trim(mpp_get_field_name(field%field)), local_buffer)
                    else
                         call write_data (fileobjND (file_num), trim(mpp_get_field_name(field%field)), local_buffer, unlim_dim_level=time)
                    endif
               endif
          elseif (fnum_for_domain == "ug") then
                    if (time == 0) then
                         call write_data (fileobjU(file_num), trim(mpp_get_field_name(field%field)), local_buffer)
                    else
                         call write_data (fileobjU(file_num), trim(mpp_get_field_name(field%field)), local_buffer, unlim_dim_level=time)
                    endif
          else
               call error_mesg("diag_field_write","No file object is associated with this file number",fatal)
          endif
     elseif (present(file_num) ) then
          write (6,*) present(file_num) ,present(fileobjU) , present(fileobjND) , present(fileobj) , present(fnum_for_domain)
          call error_mesg("diag_field_write","When FILE_NUM is used to determine which file object to use,"&
           //" You must also include fileobjU, fileobj, fileonjND, and fnum_for_domain",fatal)
     else
          call error_mesg("diag_field_write","You must include a fileob or a file_num.",fatal)
     endif
     if (allocated(local_buffer)) deallocate(local_buffer)
  end subroutine diag_field_write_field

  !> \brief Writes diagnostic data out using fms2_io routine.
  subroutine diag_field_write_varname (varname, buffer, static, fileob, file_num, fileobjU, fileobj, fileobjND, fnum_for_domain, time_in)
    CHARACTER(len=*), INTENT(in) :: varname !<
    REAL , INTENT(inout) :: buffer(:,:,:,:)
    logical, intent(in), optional :: static
    class(FmsNetcdfFile_t), intent(inout), optional, target :: fileob
    class(FmsNetcdfFile_t), pointer :: fptr => null()
    integer, intent(in), optional  :: file_num
    type(FmsNetcdfUnstructuredDomainFile_t),intent(inout), optional :: fileobjU(:)
    type(FmsNetcdfDomainFile_t),intent(inout), optional:: fileobj(:)
    type(FmsNetcdfFile_t),intent(inout), optional:: fileobjND(:)
    character(len=2), intent(in), optional :: fnum_for_domain
    INTEGER, OPTIONAL, INTENT(in) :: time_in
    integer :: time
    real(kind=4),allocatable :: local_buffer(:,:,:,:)
!> Set up the time.  Static field and default time is 0
     if (present(static) .and. static) then
          time = 0
     elseif (present(time_in)) then
          time = time_in
     else
          time = 0
     endif

     if (present(fileob)) then !> Write output to the fileob file
          fptr => fileob
          select type (fptr)
          type is (FmsNetcdfFile_t)
               call write_data (fptr,trim(varname),buffer)
          type is (FmsNetcdfDomainFile_t)
               call write_data (fptr,trim(varname),buffer)
          type is (FmsNetcdfUnstructuredDomainFile_t)
               call write_data (fptr,trim(varname),buffer)
          class default
               call error_mesg("diag_field_write","fileob passed in is not one of the FmsNetcdfFile_t types",fatal)
          end select
          call write_data (fileob,trim(varname),buffer)
     elseif (present(file_num) .and. present(fileobjU) .and. present(fileobj) .and. present(fileobjND) .and. present(fnum_for_domain)) then
     !> Figure out which file object to write output to
          if (fnum_for_domain == "2d" ) then
               if (check_if_open(fileobj(file_num))) then
                    call write_data (fileobj (file_num), trim(varname), buffer, unlim_dim_level=time )
               endif
          elseif (fnum_for_domain == "nd") then
               if (check_if_open(fileobjND (file_num)) ) then
                    call write_data (fileobjND (file_num), trim(varname), buffer, unlim_dim_level=time)
               endif
          elseif (fnum_for_domain == "ug") then
               call write_data (fileobjU(file_num), trim(varname), buffer, unlim_dim_level=time)
          else
               call error_mesg("diag_field_write","No file object is associated with this file number",fatal)
          endif
     elseif (present(file_num) ) then
          call error_mesg("diag_field_write","When FILE_NUM is used to determine which file object to use,"&
           //" You must also include fileobjU, fileobj, and fnum_for_domain",fatal)
     else
          call error_mesg("diag_field_write","You must include a fileob or a file_num.",fatal)
     endif
  end subroutine diag_field_write_varname
!> \brief Writes the time data to the history file
  subroutine diag_write_time (fileob,rtime_value,time_index,time_name)
     class(FmsNetcdfFile_t), intent(inout),target  :: fileob      !< fms2_io file object
     class(FmsNetcdfFile_t), pointer                        :: fptr => null()
     real, intent(in)                                       :: rtime_value !< The value of time to be written
     integer, intent(in)                                    :: time_index  !< The index of the time variable
     character(len=*),intent(in),optional                   :: time_name   !< The name of the time variable
     character(len=:),allocatable                           :: name_time   !< The name of the time variable
!> Get the name of the time variable
     if (present(time_name)) then
          allocate(character(len=len(time_name)) :: name_time)
          name_time = time_name
     else
          allocate(character(len=4) :: name_time)
          name_time = "time"
     endif
!> Write the time data
     call write_data (fileob, trim(name_time), rtime_value, unlim_dim_level=time_index)
!> Cleanup
     if (allocated(name_time)) deallocate(name_time)
     if (associated(fptr)) nullify(fptr)
  end subroutine diag_write_time

  !> @brief Return the axis index number.
  !! @return Integer index
  FUNCTION get_axis_index(num) RESULT ( index )
    INTEGER, INTENT(in) :: num

    INTEGER :: index
    INTEGER :: i

    !---- get the array index for this axis type ----
    !---- set up pointers to axistypes ----
    !---- write axis meta data for new axes ----
    index = 0
    DO i = 1, num_axis_in_file
       IF ( num == axis_in_file(i) ) THEN
          index = i
          EXIT
       END IF
    END DO
  END FUNCTION get_axis_index

  !> @brief Return the global attribute type.
  SUBROUTINE get_diag_global_att(gAtt)
    TYPE(diag_global_att_type), INTENT(out) :: gAtt

    gAtt=diag_global_att
  END SUBROUTINE get_diag_global_att

  !> @brief Set the global attribute type.
  SUBROUTINE set_diag_global_att(component, gridType, tileName)
    CHARACTER(len=*),INTENT(in) :: component, gridType, tileName

    ! The following two lines are set to remove compile time warnings
    ! about 'only used once'.
    CHARACTER(len=64) :: component_tmp
    component_tmp = component
    ! Don't know how to set these for specific component
    ! Want to be able to say
    ! if(output_file has component) then
    diag_global_att%grid_type = gridType
    diag_global_att%tile_name = tileName
    ! endif
  END SUBROUTINE set_diag_global_att
END MODULE diag_output_mod
!> @}
! close documentation grouping
