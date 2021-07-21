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
!> @defgroup diag_axis_mod diag_axis_mod
!> @ingroup diag_manager
!> @brief An integral part of @ref diag_manager_mod. It helps to create axis IDs
!! that are used in @ref register_diag_field.
!!
!> @author Seth Underwood
!!
!! Users first create axis ID by calling diag_axis_init, then use this axis ID in
!! register_diag_field.

!> @file
!> @brief File for @ref diag_axis_mod

!> @addtogroup diag_axis_mod
!> @{
MODULE diag_axis_mod
use platform_mod

  USE mpp_domains_mod, ONLY: domainUG, domain1d, domain2d, mpp_get_compute_domain,&
       & mpp_get_domain_components, null_domain1d, null_domain2d, null_domainUG,&
       & NORTH, EAST, CENTER, &
       & OPERATOR(.NE.), mpp_get_global_domain, mpp_get_domain_name
  USE fms_mod, ONLY: error_mesg, write_version_number, lowercase, uppercase,&
       & fms_error_handler, FATAL, NOTE
  USE diag_data_mod, ONLY: diag_axis_type, max_subaxes, max_axes,&
       & max_num_axis_sets, max_axis_attributes, debug_diag_manager,&
       & first_send_data_call, diag_atttype
#ifdef use_netCDF
  USE netcdf, ONLY: NF90_INT, NF90_FLOAT, NF90_CHAR
#endif

  IMPLICIT NONE

  PRIVATE
  PUBLIC  diag_axis_init, get_diag_axis, get_domain1d, get_domain2d,&
       & get_axis_length, get_axis_global_length, diag_subaxes_init,&
       & get_diag_axis_cart, get_diag_axis_data, max_axes, get_axis_aux,&
       & get_tile_count, get_axes_shift, get_diag_axis_name,&
       & get_axis_num, get_diag_axis_domain_name, diag_axis_add_attribute,&
       & get_domainUG, axis_compatible_check, axis_is_compressed, &
       & get_compressed_axes_ids, get_axis_reqfld, &
       & NORTH, EAST, CENTER

  ! Include variable "version" to be written to log file
#include<file_version.h>

!----------
  integer(I4_KIND),parameter,public :: DIAG_AXIS_NODOMAIN = 0 !< For unstructured grid support
  integer(I4_KIND),parameter,public :: DIAG_AXIS_2DDOMAIN = 1 !< For unstructured grid support
  integer(I4_KIND),parameter,public :: DIAG_AXIS_UGDOMAIN = 2 !< For unstructured grid support
!----------

  INTEGER, DIMENSION(:), ALLOCATABLE :: num_subaxes !< counter of number of axes defined
  INTEGER :: num_def_axes = 0

  CHARACTER(len=128), DIMENSION(:), ALLOCATABLE, SAVE :: Axis_sets !< storage for axis set names
  INTEGER :: num_axis_sets = 0

  TYPE(diag_axis_type), ALLOCATABLE, SAVE :: Axes(:) !< global storage for all defined axes
  LOGICAL :: module_is_initialized = .FALSE.

  !> @}

  !> @brief Add an arbitrary attribute and value to the diagnostic axis.
  !!
  !> Any number of attributes can be added to a given axis.  All attribute addition must
  !! be done before first <TT>send_data</TT> call.<br>
  !!
  !! If a real or integer attribute is already defined, a FATAL error will be called.
  !! If a character attribute is already defined, then it will be prepended to the
  !! existing attribute value.
  !! <br>Example usage:
  !! @code{.F90} call diag_axis_add_attribute(diag_axis_id, att_name, att_value) @endcode
  !> @ingroup diag_axis_mod
  INTERFACE diag_axis_add_attribute
     MODULE PROCEDURE diag_axis_add_attribute_scalar_r
     MODULE PROCEDURE diag_axis_add_attribute_scalar_i
     MODULE PROCEDURE diag_axis_add_attribute_scalar_c
     MODULE PROCEDURE diag_axis_add_attribute_r1d
     MODULE PROCEDURE diag_axis_add_attribute_i1d
  END INTERFACE diag_axis_add_attribute

  !> @addtogroup diag_axis_mod
  !> @{

CONTAINS

  !> @brief Initialize the axis, and return the axis ID.
  !!
  !> <TT>diag_axis_init</TT> initializes an axis and returns the axis ID that
  !! is to be used with <TT>register_diag_field</TT>.  This function also
  !! increments the axis counter and fills in the axes
  !!
  !! @return integer axis ID
  INTEGER FUNCTION diag_axis_init(name, DATA, units, cart_name, long_name, direction,&
       & set_name, edges, Domain, Domain2, DomainU, aux, req, tile_count, domain_position )
    CHARACTER(len=*), INTENT(in) :: name !< Short name for axis
    REAL, DIMENSION(:), INTENT(in) :: DATA !< Array of coordinate values
    CHARACTER(len=*), INTENT(in) :: units !< Units for the axis
    CHARACTER(len=*), INTENT(in) :: cart_name !< Cartesian axis ("X", "Y", "Z", "T")
    CHARACTER(len=*), INTENT(in), OPTIONAL :: long_name !< Long name for the axis.
    CHARACTER(len=*), INTENT(in), OPTIONAL :: set_name
    INTEGER, INTENT(in), OPTIONAL :: direction !< Indicates the direction of the axis
    INTEGER, INTENT(in), OPTIONAL :: edges !< Axis ID for the previously defined "edges axis"
    TYPE(domain1d), INTENT(in), OPTIONAL :: Domain
    TYPE(domain2d), INTENT(in), OPTIONAL :: Domain2
    TYPE(domainUG), INTENT(in), OPTIONAL :: DomainU
    CHARACTER(len=*), INTENT(in), OPTIONAL :: aux !< Auxiliary name, can only be <TT>geolon_t</TT> or <TT>geolat_t</TT>
    CHARACTER(len=*), INTENT(in), OPTIONAL :: req !< Required field names.
    INTEGER, INTENT(in), OPTIONAL :: tile_count
    INTEGER, INTENT(in), OPTIONAL :: domain_position


    TYPE(domain1d) :: domain_x, domain_y
    INTEGER :: ierr, axlen
    INTEGER :: i, set, tile
    INTEGER :: isc, iec, isg, ieg
    CHARACTER(len=128) :: emsg

    IF ( .NOT.module_is_initialized ) THEN
       CALL write_version_number("DIAG_AXIS_MOD", version)
    ENDIF

    IF ( PRESENT(tile_count)) THEN
       tile = tile_count
    ELSE
       tile = 1
    END IF

    ! Allocate the axes
    IF (.NOT. ALLOCATED(Axis_sets)) ALLOCATE(Axis_sets(max_num_axis_sets))
    IF (.NOT. ALLOCATED(Axes)) ALLOCATE(Axes(max_axes))
    IF (.NOT. ALLOCATED(num_subaxes)) THEN
       ALLOCATE(num_subaxes(max_axes))
       num_subaxes = 0
    END IF

    !---- is there an axis set? ----
    IF ( PRESENT(set_name) ) THEN
       set = get_axis_set_num (set_name)
       !---- add new set name ----
       IF (set == 0) THEN
          num_axis_sets = num_axis_sets + 1
          IF ( num_axis_sets > max_num_axis_sets ) THEN
             WRITE (emsg, FMT='("num_axis_sets (",I2,") exceeds max_num_axis_sets (",I2,"). ")')&
                  & num_axis_sets, max_num_axis_sets
             ! <ERROR STATUS="FATAL">
             !   num_axis_sets (<num_axis_sets>) exceeds max_num_axis_sets(<num_axis_sets>).
             !   Increase max_num_axis_sets via diag_manager_nml.
             ! </ERROR>
             CALL error_mesg('diag_axis_mod::diag_axis_init',&
                  & TRIM(emsg)//'  Increase max_num_axis_sets via diag_manager_nml.', FATAL)
          END IF
          set = num_axis_sets
          Axis_sets(set) = set_name
       END IF
    ELSE
       set = 0
    END IF

    !---- see if axis already exists --
    ! if this is time axis, return the ID of a previously defined
    ! if this is spatial axis, FATAL error
    DO i = 1, num_def_axes
       IF ( TRIM(name) == Axes(i)%name ) THEN
          IF ( TRIM(name) == 'Stations' .OR. TRIM(name) == 'Levels') THEN
             diag_axis_init = i
             RETURN
          ELSE IF ( set == Axes(i)%set ) THEN
             IF ( TRIM(lowercase(name)) == 'time' .OR.&
                  & TRIM(lowercase(cart_name)) == 't' .OR.&
                  & TRIM(lowercase(name)) == 'nv' .OR.&
                  & TRIM(lowercase(cart_name)) == 'n' ) THEN
                diag_axis_init = i
                RETURN
             ELSE IF ( (lowercase(cart_name) /= 'x' .AND. lowercase(cart_name) /= 'y')&
                  & .OR. tile /= Axes(i)%tile_count) THEN
                ! <ERROR STATUS="FATAL">axis_name <NAME> and axis_set already exist.</ERROR>
                CALL error_mesg('diag_axis_mod::diag_axis_init',&
                     & 'axis_name '//TRIM(name)//' and axis_set already exist.', FATAL)
             END IF
          END IF
       END IF
    END DO

    !---- register axis ----
    num_def_axes = num_def_axes + 1
    ! <ERROR STATUS="FATAL">max_axes exceeded, increase it via diag_manager_nml</ERROR>
    IF (num_def_axes > max_axes) CALL error_mesg ('diag_axis_mod::diag_axis_init',&
         & 'max_axes exceeded, increase via diag_manager_nml', FATAL)
    diag_axis_init = num_def_axes

    !---- check and then save cart_name name ----
    IF ( TRIM(uppercase(cart_name)) == 'X' .OR.&
         & TRIM(uppercase(cart_name)) == 'Y' .OR.&
         & TRIM(uppercase(cart_name)) == 'Z' .OR.&
         & TRIM(uppercase(cart_name)) == 'T' .OR.&
         & TRIM(uppercase(cart_name)) == 'U' .OR.&
         & TRIM(uppercase(cart_name)) == 'N' ) THEN
       Axes(diag_axis_init)%cart_name = TRIM(uppercase(cart_name))
    ELSE
       ! <ERROR STATUS="FATAL">Invalid cart_name name.</ERROR>
       CALL error_mesg('diag_axis_mod::diag_axis_init', 'Invalid cart_name name. '//TRIM(uppercase(cart_name)), FATAL)
    END IF

    !---- allocate storage for coordinate values of axis ----
    IF ( Axes(diag_axis_init)%cart_name == 'T' ) THEN
       axlen = 0
    ELSE
       axlen = SIZE(DATA(:))
    END IF
    ALLOCATE ( Axes(diag_axis_init)%data(1:axlen) )

    ! Initialize Axes(diag_axis_init)
    Axes(diag_axis_init)%name   = TRIM(name)
    Axes(diag_axis_init)%data   = DATA(1:axlen)
    Axes(diag_axis_init)%units  = units
    Axes(diag_axis_init)%length = axlen
    Axes(diag_axis_init)%set    = set
    ! start and end are used in subaxes information only
    Axes(diag_axis_init)%start = -1
    Axes(diag_axis_init)%end = -1
    Axes(diag_axis_init)%subaxis_name = ""
    Axes(diag_axis_init)%shift = 0
    Axes(diag_axis_init)%num_attributes = 0

    IF ( PRESENT(long_name) ) THEN
       Axes(diag_axis_init)%long_name = long_name
    ELSE
       Axes(diag_axis_init)%long_name = name
    END IF

    IF ( PRESENT(aux) ) THEN
       Axes(diag_axis_init)%aux = TRIM(aux)
    ELSE
       Axes(diag_axis_init)%aux = 'none'
    END IF

    IF ( PRESENT(req) ) THEN
       Axes(diag_axis_init)%req = TRIM(req)
    ELSE
       Axes(diag_axis_init)%req = 'none'
    END IF
    IF ( PRESENT(domain_position) ) THEN
       if (domain_position == NORTH .or. domain_position == EAST .or. domain_position == CENTER) then
          Axes(diag_axis_init)%domain_position = domain_position
       else
          CALL error_mesg('diag_axis_mod::diag_axis_init', "Position must be NORTH, EAST, or CENTER" ,&
                         FATAL)
       endif
    ELSE
       Axes(diag_axis_init)%domain_position = CENTER
    END IF

    !---- axis direction (-1, 0, or +1) ----
    IF ( PRESENT(direction) )THEN
       IF ( ABS(direction) /= 1 .AND. direction /= 0 )&
            ! <ERROR STATUS="FATAL">direction must be 0, +1, or -1</ERROR>
            & CALL error_mesg('diag_axis_mod::diag_axis_init', 'direction must be 0, +1 or -1', FATAL)
       Axes(diag_axis_init)%direction = direction
    ELSE
       Axes(diag_axis_init)%direction = 0
    END IF

    !---- Handle the DomainU check
    IF (present(DomainU) .AND. (PRESENT(Domain2) .OR. PRESENT(Domain)) ) THEN
       ! <ERROR STATUS="FATAL">Presence of DomainU and another Domain at the same time is prohibited</ERROR>
       CALL error_mesg('diag_axis_mod::diag_axis_init',&
            & 'Presence of DomainU and another Domain at the same time is prohibited', FATAL)
    !---- domain2d type ----
    ELSE IF ( PRESENT(Domain2) .AND. PRESENT(Domain)) THEN
       ! <ERROR STATUS="FATAL">Presence of both Domain and Domain2 at the same time is prohibited</ERROR>
       CALL error_mesg('diag_axis_mod::diag_axis_init',&
            & 'Presence of both Domain and Domain2 at the same time is prohibited', FATAL)
    ELSE IF ( PRESENT(Domain2) .OR. PRESENT(Domain)) THEN
       IF ( Axes(diag_axis_init)%cart_name /= 'X' .AND. Axes(diag_axis_init)%cart_name /= 'Y') THEN
          ! <ERROR STATUS="FATAL">Domain must not be present for an axis which is not in the X or Y direction.</ERROR>
          CALL error_mesg('diag_axis_mod::diag_axis_init',&
               & 'A Structured Domain must not be present for an axis which is not in the X or Y direction', FATAL)
       END IF
    ELSE IF (present(DomainU) .AND. Axes(diag_axis_init)%cart_name /= 'U') THEN
          CALL error_mesg('diag_axis_mod::diag_axis_init',&
               & 'In the unstructured domain, the axis cart_name must be U', FATAL)
    END IF

    Axes(diag_axis_init)%tile_count = tile

    IF ( PRESENT(Domain2) ) THEN
       Axes(diag_axis_init)%Domain2 = Domain2
       CALL mpp_get_domain_components(Domain2, domain_x, domain_y, tile_count=tile_count)
       IF ( Axes(diag_axis_init)%cart_name == 'X' ) Axes(diag_axis_init)%Domain = domain_x
       IF ( Axes(diag_axis_init)%cart_name == 'Y' ) Axes(diag_axis_init)%Domain = domain_y
       Axes(diag_axis_init)%DomainUG = null_DomainUG
    ELSE IF ( PRESENT(Domain)) THEN
       !---- domain1d type ----
       Axes(diag_axis_init)%Domain2 = null_domain2d ! needed since not 2-D domain
       Axes(diag_axis_init)%Domain = Domain
       Axes(diag_axis_init)%DomainUG = null_DomainUG
    ELSE IF (present(DomainU)) THEN
       Axes(diag_axis_init)%Domain2 = null_domain2d
       Axes(diag_axis_init)%Domain = null_domain1d
       Axes(diag_axis_init)%DomainUG = DomainU
    ELSE
       Axes(diag_axis_init)%Domain2 = null_domain2d
       Axes(diag_axis_init)%Domain = null_domain1d
       Axes(diag_axis_init)%DomainUG = null_domainUG
    END IF

    !--- set up the shift value for x-y axis
    IF ( Axes(diag_axis_init)%Domain .NE. null_domain1d ) THEN
       CALL mpp_get_compute_domain(Axes(diag_axis_init)%Domain, isc, iec)
       CALL mpp_get_global_domain(Axes(diag_axis_init)%Domain, isg, ieg)
       IF ( Axes(diag_axis_init)%length == ieg - isg + 2 ) THEN
          Axes(diag_axis_init)%shift = 1
       END IF
    END IF

    !---- have axis edges been defined ? ----
    Axes(diag_axis_init)%edges = 0
    IF (PRESENT(edges) ) THEN
       IF ( edges > 0 .AND. edges < num_def_axes ) THEN
          ierr=0
          IF ( Axes(edges)%cart_name /= Axes(diag_axis_init)%cart_name) ierr=1
          IF ( Axes(edges)%length    /= Axes(diag_axis_init)%length+1 ) ierr=ierr+2
          IF ( Axes(edges)%set       /= Axes(diag_axis_init)%set      ) ierr=ierr+4
          IF ( ierr > 0 )   THEN
             ! <ERROR STATUS="FATAL">Edges axis does not match axis (code <CODE>).</ERROR>
             WRITE (emsg,'("Edges axis does not match axis (code ",I1,").")') ierr
             CALL error_mesg('diag_axis_mod::diag_axis_init', emsg, FATAL)
          END IF
          Axes(diag_axis_init)%edges = edges
       ELSE
          ! <ERROR STATUS="FATAL">Edges axis is not defined.</ERROR>
          CALL error_mesg('diag_axis_mod::diag_axis_init', 'Edges axis is not defined', FATAL)
       END IF
    END IF

    ! Module is now initialized
    module_is_initialized = .TRUE.

  END FUNCTION diag_axis_init

  !> @brief Create a subaxis on a parent axis.
  !!
  !> Given the ID of a parent axis, create a subaxis and fill it with data,
  !!     and return the ID of the corresponding subaxis.
  !!
  !!     The subaxis is defined on the parent axis from <TT>start_indx</TT>
  !!     to <TT>end_indx</TT>.
  !!
  !! @return Integer ID of the corresponding subaxis.
  INTEGER FUNCTION diag_subaxes_init(axis, subdata, start_indx, end_indx, domain_2d)
    INTEGER, INTENT(in) :: axis !< ID of the parent axis
    REAL, DIMENSION(:), INTENT(in) :: subdata !< Data of the subaxis
    INTEGER, INTENT(in) :: start_indx !< Start index of the subaxis
    INTEGER, INTENT(in) :: end_indx !< End index of the subaxis
    TYPE(domain2d), INTENT(in), OPTIONAL  :: domain_2d

    INTEGER :: i, nsub_axis, direction
    INTEGER :: xbegin, xend, ybegin, yend
    INTEGER :: ad_xbegin, ad_xend, ad_ybegin, ad_yend
    CHARACTER(len=128) :: name, nsub_name
    CHARACTER(len=128) :: units
    CHARACTER(len=128) :: cart_name
    CHARACTER(len=128) :: long_name
    CHARACTER(len=128) :: emsg
    LOGICAL :: subaxis_set, hasDomain

    ! there may be more than 1 subaxis on a parent axis, check for redundancy
    nsub_axis = 0
    subaxis_set = .FALSE.

    IF ( PRESENT(domain_2d) ) THEN
       hasDomain = .TRUE.
       CALL mpp_get_compute_domain(domain_2d, xbegin, xend, ybegin, yend)
    ELSE
       hasDomain = .FALSE.
    END IF
    sa_search: DO i = 1, num_subaxes(axis)
       IF ( start_indx == Axes(axis)%start(i) .AND. end_indx == Axes(axis)%end(i) ) THEN
          IF ( hasDomain ) THEN
             CALL mpp_get_compute_domain(Axes(axis)%subaxis_domain2(i), ad_xbegin, ad_xend, ad_ybegin, ad_yend)
             IF ( .NOT.((xbegin == ad_xbegin .AND. xend == ad_xend) .AND.&
                  & (ybegin == ad_ybegin .AND. yend == ad_yend)) ) THEN
                CYCLE sa_search
             END IF
          END IF
          nsub_axis = i
          subaxis_set = .TRUE.    !subaxis already exists
          name = TRIM(Axes(axis)%subaxis_name(nsub_axis))
          EXIT sa_search
       END IF
    END DO sa_search

    IF ( nsub_axis == 0 ) THEN  ! create new subaxis
       num_subaxes(axis) = num_subaxes(axis) + 1
       IF (num_subaxes(axis) > max_subaxes) THEN
          ! <ERROR STATUS="FATAL">max_subaxes (value <VALUE>) is too small.  Consider increasing max_subaxes.</ERROR>
          WRITE (emsg,'("max_subaxes (value ",I4,") is too small.  Consider increasing max_subaxes.")') max_subaxes
          CALL error_mesg('diag_axis_mod::diag_subaxes_init', emsg, FATAL)
       END IF
       nsub_axis = num_subaxes(axis)
       Axes(axis)%start(nsub_axis) = start_indx
       Axes(axis)%end(nsub_axis)   = end_indx
       if ( hasDomain ) Axes(axis)%subaxis_domain2(nsub_axis) = domain_2d
    END IF

    ! Create new name for the subaxis from name of parent axis
    ! If subaxis already exists, get the index and return
    IF(subaxis_set) THEN
       IF ( Axes(axis)%set > 0 ) THEN
          diag_subaxes_init = get_axis_num(name, set_name=TRIM(Axis_sets(Axes(axis)%set)))
       ELSE
          diag_subaxes_init = get_axis_num(name)
       END IF
    ELSE
       ! get a new index for subaxis
       !::sdu:: Need a check to allow larger numbers in the index number.
       WRITE (nsub_name,'(I2.2)') nsub_axis
       name = TRIM(Axes(axis)%name)//'_sub'//TRIM(nsub_name)
       Axes(axis)%subaxis_name(nsub_axis) = name
       long_name = TRIM(Axes(axis)%long_name)
       units = TRIM(Axes(axis)%units)
       cart_name = TRIM(Axes(axis)%cart_name)
       direction = Axes(axis)%direction
       IF (Axes(axis)%set > 0) THEN
          diag_subaxes_init =  diag_axis_init (TRIM(name), subdata, TRIM(units), TRIM(cart_name), TRIM(long_name),&
               & set_name=TRIM(Axis_sets(Axes(axis)%set)), direction=direction, Domain2=domain_2d)
       ELSE
          diag_subaxes_init =  diag_axis_init (TRIM(name), subdata, TRIM(units), TRIM(cart_name), TRIM(long_name),&
               & direction=direction, Domain2=domain_2d)
       END IF
    END IF
  END FUNCTION diag_subaxes_init
  !> @brief Return information about the axis with index ID
  SUBROUTINE get_diag_axis(id, name, units, long_name, cart_name,&
       & direction, edges, Domain, DomainU, DATA, num_attributes, attributes, domain_position)
    CHARACTER(len=*), INTENT(out) :: name, units, long_name, cart_name
    INTEGER, INTENT(in) :: id !< Axis ID
    TYPE(domain1d), INTENT(out) :: Domain
    TYPE(domainUG), INTENT(out) :: DomainU
    INTEGER, INTENT(out) :: direction !< Direction of data. (See <TT>@ref diag_axis_init</TT> for a description of
                                      !! allowed values)
    INTEGER, INTENT(out) :: edges !< Axis ID for the previously defined "edges axis".
    REAL, DIMENSION(:), INTENT(out) :: DATA !< Array of coordinate values for this axis.
    INTEGER, INTENT(out), OPTIONAL :: num_attributes
    TYPE(diag_atttype), ALLOCATABLE, DIMENSION(:), INTENT(out), OPTIONAL :: attributes
    INTEGER, INTENT(out), OPTIONAL :: domain_position

    INTEGER :: i, j, istat

    CALL valid_id_check(id, 'get_diag_axis')
    name      = Axes(id)%name
    units     = Axes(id)%units
    long_name = Axes(id)%long_name
    cart_name = Axes(id)%cart_name
    direction = Axes(id)%direction
    edges     = Axes(id)%edges
    Domain    = Axes(id)%Domain
    DomainU   = Axes(id)%DomainUG
    if (present(domain_position)) domain_position = Axes(id)%domain_position
    IF ( Axes(id)%length > SIZE(DATA(:)) ) THEN
       ! <ERROR STATUS="FATAL">array data is too small.</ERROR>
       CALL error_mesg('diag_axis_mod::get_diag_axis', 'array data is too small', FATAL)
    ELSE
       DATA(1:Axes(id)%length) = Axes(id)%data(1:Axes(id)%length)
    END IF
    IF ( PRESENT(num_attributes) ) THEN
       num_attributes = Axes(id)%num_attributes
    END IF
    IF ( PRESENT(attributes) ) THEN
       IF ( allocated(Axes(id)%attributes) ) THEN
          IF ( ALLOCATED(attributes) ) THEN
             ! If allocate, make sure attributes is large enough to hold Axis(id)%attributes
             IF ( Axes(id)%num_attributes .GT. SIZE(attributes(:)) ) THEN
                CALL error_mesg('diag_axis_mod::get_diag_axis', 'array attribute is too small', FATAL)
             END IF
          ELSE
             ! Allocate attributes
             ALLOCATE(attributes(Axes(id)%num_attributes), STAT=istat)
             IF ( istat .NE. 0 ) THEN
                CALL error_mesg('diag_axis_mod::get_diag_axis', 'Unable to allocate memory for attribute', FATAL)
             END IF
          END IF
          DO i=1, Axes(id)%num_attributes
             ! Unallocate all att arrays in preparation for new data
             IF ( allocated(attributes(i)%fatt) ) THEN
                DEALLOCATE(attributes(i)%fatt)
             END IF
             IF ( allocated(attributes(i)%iatt) ) THEN
                DEALLOCATE(attributes(i)%iatt)
             END IF

             ! Copy in attribute data
             attributes(i)%type = Axes(id)%attributes(i)%type
             attributes(i)%len = Axes(id)%attributes(i)%len
             attributes(i)%name = Axes(id)%attributes(i)%name
             attributes(i)%catt = Axes(id)%attributes(i)%catt
             ! Allocate fatt arrays (if needed), and copy in data
             IF ( allocated(Axes(id)%attributes(i)%fatt) ) THEN
                ALLOCATE(attributes(i)%fatt(SIZE(Axes(id)%attributes(i)%fatt(:))), STAT=istat)
                IF ( istat .NE. 0 ) THEN
                   CALL error_mesg('diag_axis_mod::get_diag_axis', 'Unable to allocate memory for attribute%fatt', FATAL)
                END IF
                DO j=1, SIZE(attributes(i)%fatt(:))
                   attributes(i)%fatt(j) = Axes(id)%attributes(i)%fatt(j)
                END DO
             END IF
             ! Allocate iatt arrays (if needed), and copy in data
             IF ( allocated(Axes(id)%attributes(i)%iatt) ) THEN
                ALLOCATE(attributes(i)%iatt(SIZE(Axes(id)%attributes(i)%iatt(:))), STAT=istat)
                IF ( istat .NE. 0 ) THEN
                   CALL error_mesg('diag_axis_mod::get_diag_axis', 'Unable to allocate memory for attribute%iatt', FATAL)
                END IF
                DO j=1, SIZE(attributes(i)%iatt(:))
                   attributes(i)%iatt(j) = Axes(id)%attributes(i)%iatt(j)
                END DO
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE get_diag_axis

  !> @brief Return the axis cartesian.
  SUBROUTINE get_diag_axis_cart(id, cart_name)
    INTEGER, INTENT(in)           :: id !< Axis ID
    CHARACTER(len=*), INTENT(out) :: cart_name !< Cartesian axis

    CALL valid_id_check(id, 'get_diag_axis_cart')
    cart_name = Axes(id)%cart_name
  END SUBROUTINE get_diag_axis_cart

  !> @brief Return the axis data.
  SUBROUTINE get_diag_axis_data(id, DATA)
    INTEGER, INTENT(in) :: id !< Axis ID
    REAL, DIMENSION(:), INTENT(out) :: DATA !< Axis data

    CALL valid_id_check(id, 'get_diag_axis_data')
    IF (Axes(id)%length > SIZE(DATA(:))) THEN
       ! <ERROR STATUS="FATAL">array data is too small</ERROR>
       CALL error_mesg('diag_axis_mod::get_diag_axis_data', 'array data is too small', FATAL)
    ELSE
       DATA(1:Axes(id)%length) = Axes(id)%data
    END IF
  END SUBROUTINE get_diag_axis_data

  !> @brief Return the short name of the axis.
  SUBROUTINE get_diag_axis_name(id, name)
    INTEGER         , INTENT(in)  :: id !< Axis ID
    CHARACTER(len=*), INTENT(out) :: name !< Axis short name

    CALL valid_id_check(id, 'get_diag_axis_name')
    name = Axes(id)%name
  END SUBROUTINE get_diag_axis_name

  !> @brief Return the name of the axis' domain
  SUBROUTINE get_diag_axis_domain_name(id, name)
    INTEGER, INTENT(in) :: id !< Axis ID
    CHARACTER(len=*), INTENT(out) :: name !< Axis' domain name

    CALL valid_id_check(id, 'get_diag_axis_domain_name')
    name = mpp_get_domain_name(Axes(id)%domain2)
  END SUBROUTINE get_diag_axis_domain_name

  !> @brief Return the length of the axis.
  !> @return length of axis as an integer
  INTEGER FUNCTION get_axis_length(id)
    INTEGER, INTENT(in) :: id !< Axis ID
    INTEGER :: length

    CALL valid_id_check(id, 'get_axis_length')
    IF ( Axes(id)%Domain .NE. null_domain1d ) THEN
       CALL mpp_get_compute_domain(Axes(id)%Domain,size=length)
       !---one extra point is needed for some case. ( like symmetry domain )
       get_axis_length = length + Axes(id)%shift
    ELSE
       get_axis_length = Axes(id)%length
    END IF
  END FUNCTION get_axis_length

  !> @brief Return the auxiliary name for the axis.
  !! @return auxiliary name for the axis
  CHARACTER(len=128) FUNCTION get_axis_aux(id)
    INTEGER, INTENT(in) :: id !< Axis ID

    CALL valid_id_check(id, 'get_axis_aux')
    get_axis_aux =  Axes(id)%aux
  END FUNCTION get_axis_aux

  !> @brief Return the required field names for the axis.
  !! @return required field names for the axis
  CHARACTER(len=128) FUNCTION get_axis_reqfld(id)
    INTEGER, INTENT(in) :: id !< Axis ID

    CALL valid_id_check(id, 'get_axis_reqfld')
    get_axis_reqfld =  Axes(id)%req
  END FUNCTION get_axis_reqfld

  !> @brief Return the global length of the axis.
  !! @return global length of the axis
  INTEGER FUNCTION get_axis_global_length(id)
    INTEGER, INTENT(in) :: id !< Axis ID

    CALL valid_id_check(id, 'get_axis_global_length')
    get_axis_global_length = Axes(id)%length
  END FUNCTION get_axis_global_length

  !> @brief Return the tile count for the axis.
  !! @return tile count for the axis
  INTEGER FUNCTION get_tile_count(ids)
    INTEGER, DIMENSION(:), INTENT(in) :: ids !< Axis IDs.
                                             !! Possible dimensions: 1 <= <TT>size(ids(:))</TT> <= 4.

    INTEGER :: i, id, flag

    IF ( SIZE(ids(:)) < 1 ) THEN
       ! <ERROR STATUS="FATAL">input argument has incorrect size.</ERROR>
       CALL error_mesg('diag_axis_mod::get_tile_count', 'input argument has incorrect size', FATAL)
    END IF
    get_tile_count = 1
    flag = 0
    DO i = 1, SIZE(ids(:))
       id = ids(i)
       CALL valid_id_check(id, 'get_tile_count')
       IF ( Axes(id)%cart_name == 'X' .OR.  &
            Axes(id)%cart_name == 'Y' ) flag = flag + 1
       !     --- both x/y axes found ---
       IF ( flag == 2 ) THEN
          get_tile_count = Axes(id)%tile_count
          EXIT
       END IF
    END DO
  END FUNCTION get_tile_count

  !> @brief Retrun the 1D domain for the axis ID given.
  !! @return 1D domain for the axis ID given
  TYPE(domain1d) FUNCTION get_domain1d(id)
    INTEGER, INTENT(in) :: id !< Axis ID

    CALL valid_id_check(id, 'get_domain1d')
    IF (Axes(id)%Domain .NE. NULL_DOMAIN1D) THEN
       get_domain1d = Axes(id)%Domain
    ELSE
       get_domain1d = NULL_DOMAIN1D
    ENDIF
  END FUNCTION get_domain1d

  !> @brief Return the 2D domain for the axis IDs given.
  !! @return 2D domain for the axis IDs given
  TYPE(domain2d) FUNCTION get_domain2d(ids)
    INTEGER, DIMENSION(:), INTENT(in) :: ids !< Axis IDs.
                                             !! Possible dimensions: 1 <= <TT>size(ids(:))</TT> <= 4.

    INTEGER :: i, id, flag

    IF ( SIZE(ids(:)) < 1 ) THEN
       ! <ERROR STATUS="FATAL">input argument has incorrect size.</ERROR>
       CALL error_mesg('diag_axis_mod::get_domain2d', 'input argument has incorrect size', FATAL)
    END IF
    get_domain2d = null_domain2d
    flag = 0
    DO i = 1, SIZE(ids(:))
       id = ids(i)
       CALL valid_id_check(id, 'get_domain2d')
       IF ( Axes(id)%cart_name == 'X' .OR. Axes(id)%cart_name == 'Y' ) flag = flag + 1
       !     --- both x/y axes found ---
       IF ( flag == 2 ) THEN
          IF (Axes(id)%Domain2 .NE. NULL_DOMAIN2D) get_domain2d = Axes(id)%Domain2
          EXIT
       END IF
    END DO
  END FUNCTION get_domain2d

  !> @brief Retrun the 1D domain for the axis ID given.
  !! @return 1D domain for the axis ID given
  TYPE(domainUG) FUNCTION get_domainUG(id)
    INTEGER, INTENT(in) :: id !< Axis ID

    CALL valid_id_check(id, 'get_domainUG')
    IF (Axes(id)%DomainUG .NE. NULL_DOMAINUG) THEN
       get_domainUG = Axes(id)%DomainUG
    ELSE
       get_domainUG = NULL_DOMAINUG
    ENDIF
  END FUNCTION get_domainUG

!ug support
  !> @brief Checks if the axes are compatible
  !! @return integer domain_type
  function axis_compatible_check(id,varname) result(domain_type)

   !Inputs/Outputs
    integer,dimension(:),intent(in)  :: id          !<The array of axis IDs
    character(*),intent(in),optional :: varname     !<The name of the variable
    integer(I4_KIND)                :: domain_type !<DIAG_AXIS_NODOMAIN = no domain.
                                                    !<DIAG_AXIS_2DDOMAIN = structured domain.
                                                    !<DIAG_AXIS_UGDOMAIN = unstructured domain.

   !Local variables
    logical :: XorY          !<XorY set to true if X or Y is found as a cart_name.
    logical :: UG            !<UG set to true if U is found as a cart_name.
    integer :: n             !<Looping index.
    logical :: uses_domain2D !<True if an axis is associated with a 2D domain.
    logical :: uses_domainUG !<True if an axis is associated with an unstructured domain.

   !Initialize flags.
    XorY = .false.
    UG = .false.
    uses_domain2D = .false.
    uses_domainUG = .false.

   !Make sure that valid set of axes was passed, and determine the domain-type
   !associated with the axes.
    do n = 1,size(id)
        call valid_id_check(id(n), &
                            "axis_compatible_check")
        if (Axes(id(n))%cart_name .eq. "X" .or. &
            Axes(id(n))%cart_name .eq. "Y") then
            XorY = .true.
        elseif (Axes(id(n))%cart_name .eq. "U") then
            UG = .true.
        endif
        if (Axes(id(n))%Domain2 .ne. null_domain2d) then
            uses_domain2D = .true.
        elseif (Axes(id(n))%DomainUG .ne. null_domainUG) then
            uses_domainUG = .true.
        endif
    enddo
    if (UG .and. XorY) then
        if (present(varname)) then
            call error_mesg("axis_compatible_check", &
                            "Can not use an unstructured grid with a "// &
                            "horizontal cartesian coordinate for the field " &
                            //trim(varname), &
                            FATAL)
        else
            call error_mesg("axis_compatible_check", &
                            "Can not use an unstructured grid with a horizontal "// &
                            "cartesian coordinate", &
                            FATAL)
        endif
    endif
    if (uses_domain2D .and. uses_domainUG) then
        if (present(varname)) then
            call error_mesg("axis_compatible_check", &
                            "Can not use an unstructured grid with a"// &
                            "structured grid for the field "//trim(varname), &
                            FATAL)
        else
            call error_mesg("axis_compatible_check", &
                            "Can not use an unstructured grid with a"// &
                            "structured grid.", &
                            FATAL)
        endif
    endif
    if (uses_domain2D) then
        domain_type = DIAG_AXIS_2DDOMAIN
    elseif (uses_domainUG) then
        domain_type = DIAG_AXIS_UGDOMAIN
    else
        domain_type = DIAG_AXIS_NODOMAIN
    endif

    return
  end function axis_compatible_check

  !> @brief Return the value of the shift for the axis IDs given.
  SUBROUTINE get_axes_shift(ids, ishift, jshift)
    INTEGER, DIMENSION(:), INTENT(in) :: ids
    INTEGER, INTENT(out) :: ishift !< X shift value.
    INTEGER, INTENT(out) :: jshift !< Y shift value.

    INTEGER :: i, id

    !-- get the value of the shift.
    ishift = 0
    jshift = 0
    DO i = 1, SIZE(ids(:))
       id = ids(i)
       CALL valid_id_check(id, 'get_axes_shift')
       SELECT CASE (Axes(id)%cart_name)
       CASE ( 'X' )
          ishift = Axes(id)%shift
       CASE ( 'Y' )
          jshift = Axes(id)%shift
       END SELECT
    END DO
  END SUBROUTINE get_axes_shift

  !> @brief Returns index into axis table corresponding to a given axis name.
  !! @return Returns index into axis table corresponding to a given axis name.
  INTEGER FUNCTION get_axis_num(axis_name, set_name)
    CHARACTER(len=*), INTENT(in) :: axis_name !< Axis name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: set_name !< Set name

    INTEGER :: set, n

    IF ( PRESENT(set_name) ) THEN
       set = get_axis_set_num (TRIM(set_name))
    ELSE
       set = 0
    END IF
    get_axis_num = 0
    DO n = 1, num_def_axes
       IF ( TRIM(axis_name) == TRIM(Axes(n)%name) .AND. Axes(n)%set == set ) THEN
          get_axis_num = n
          RETURN
       END IF
    END DO
  END FUNCTION get_axis_num

  !> @brief Returns index in axis set table corresponding to a given axis set name
  !! @return Returns index in axis set table corresponding to a given axis set name
  INTEGER FUNCTION get_axis_set_num (set_name)
    CHARACTER(len=*), INTENT(in) :: set_name !< Set name

    INTEGER :: iset

    get_axis_set_num = 0
    DO iset = 1, num_axis_sets
       IF ( set_name == Axis_sets(iset) ) THEN
          get_axis_set_num = iset
          RETURN
       END IF
    END DO
  END FUNCTION get_axis_set_num

  !> @brief Check to see if the given axis id is a valid id.  If the axis id is invalid,
  !!     call a FATAL error.  If the ID is valid, just return.
  SUBROUTINE valid_id_check(id, routine_name)
    INTEGER, INTENT(in) :: id !< Axis is to check for validity
    CHARACTER(len=*), INTENT(in) :: routine_name !< Name of the subroutine checking for a valid axis id

    CHARACTER(len=5) :: emsg

    IF ( id < 1 .OR. id > num_def_axes) THEN
       ! <ERROR STATUS="FATAL">
       !   Illegal value for axis used (value <VALUE>).
       ! </ERROR>
       WRITE (emsg, '(I2)') id
       CALL error_mesg('diag_axis_mod::'//TRIM(routine_name),&
            & 'Illegal value for axis_id used (value '//TRIM(emsg)//').', FATAL)
    END IF
  END SUBROUTINE valid_id_check

  SUBROUTINE diag_axis_attribute_init(diag_axis_id, name, type, cval, ival, rval)
    INTEGER, INTENT(in) :: diag_axis_id !< input field ID, obtained from diag_axis_mod::diag_axis_init.
    CHARACTER(len=*) :: name !< Name of the attribute
    INTEGER, INTENT(in) :: type !< NetCDF type (NF_FLOAT, NF_INT, NF_CHAR)
    CHARACTER(len=*), INTENT(in), OPTIONAL :: cval !< Character string attribute value
    INTEGER, DIMENSION(:), INTENT(in), OPTIONAL :: ival !< Integer attribute value(s)
    REAL, DIMENSION(:), INTENT(in), OPTIONAL :: rval !< Real attribute value(s)

    INTEGER :: istat, length, i, j, this_attribute, out_field
    CHARACTER(len=1024) :: err_msg

    IF ( .NOT.first_send_data_call ) THEN
       ! Call error due to unable to add attribute after send_data called
       ! <ERROR STATUS="FATAL">
       !   Attempting to add attribute <name> to axis <axis_name>
       !   after first send_data call.  Too late.
       ! </ERROR>
       CALL error_mesg('diag_manager_mod::diag_axis_add_attribute', 'Attempting to add attribute "'&
            &//TRIM(name)//'" to axis after first send_data call.  Too late.', FATAL)
    END IF

    ! Simply return if diag_axis_id <= 0 --- not registered
    IF ( diag_axis_id .LE. 0 ) THEN
       RETURN
    ELSE IF ( diag_axis_id .GT. num_def_axes ) THEN
       ! Call error axis_id greater than num_def_axes.  Axis is unknown
       ! <ERROR STATUS="FATAL">
       !   Attempting to add attribute <name> to axis ID <axis_ID>, however ID unknown.
       ! </ERROR>
       WRITE(err_msg, '(I5)') diag_axis_id
       CALL error_mesg('diag_manager_mod::diag_axis_add_attribute', 'Attempting to add attribute "'&
            &//TRIM(name)//'" to axis ID "'//TRIM(err_msg)//'", however ID unknown.', FATAL)

    ELSE
       ! Allocate memory for the attributes
       CALL attribute_init_axis(Axes(diag_axis_id))

       ! Check if attribute already exists
       this_attribute = 0
       DO i=1, Axes(diag_axis_id)%num_attributes
          IF ( TRIM(Axes(diag_axis_id)%attributes(i)%name) .EQ. TRIM(name) ) THEN
             this_attribute = i
             EXIT
          END IF
       END DO

       IF ( this_attribute.NE.0 .AND. (type.EQ.NF90_INT .OR. type.EQ.NF90_FLOAT) ) THEN
          ! <ERROR STATUS="FATAL">
          !   Attribute <name> already defined for axis <axis_name>.
          !   Contact the developers
          ! </ERROR>
          CALL error_mesg('diag_manager_mod::diag_axis_add_attribute',&
               & 'Attribute "'//TRIM(name)//'" already defined for axis "'&
               &//TRIM(Axes(diag_axis_id)%name)//'".  Contact the developers.', FATAL)
       ELSE IF ( this_attribute.NE.0 .AND. type.EQ.NF90_CHAR .AND. debug_diag_manager ) THEN
          ! <ERROR STATUS="NOTE">
          !   Attribute <name> already defined for axis <axis_name>.
          !   Prepending.
          ! </ERROR>
          CALL error_mesg('diag_manager_mod::diag_axis_add_attribute',&
               & 'Attribute "'//TRIM(name)//'" already defined for axis"'&
               &//TRIM(Axes(diag_axis_id)%name)//'".  Prepending.', NOTE)
       ELSE
          ! Defining a new attribute
          ! Increase the number of field attributes
          this_attribute = Axes(diag_axis_id)%num_attributes + 1
          ! Checking to see if num_attributes == max_axis_attributes, and return error message
          IF ( this_attribute .GT. max_axis_attributes ) THEN
             ! <ERROR STATUS="FATAL">
             !   Number of attributes exceeds max_axis_attributes for attribute <name> for axis <axis_name>.
             !   Increase diag_manager_nml:max_axis_attributes.
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_axis_add_attribute',&
                  & 'Number of attributes exceeds max_axis_attributes for attribute "'&
                  & //TRIM(name)//'" for axis "'//TRIM(Axes(diag_axis_id)%name)&
                  & //'".  Increase diag_manager_nml:max_axis_attributes.',&
                  & FATAL)
          ELSE
             Axes(diag_axis_id)%num_attributes = this_attribute
             ! Set name and type
             Axes(diag_axis_id)%attributes(this_attribute)%name = name
             Axes(diag_axis_id)%attributes(this_attribute)%type = type
             ! Initialize catt to a blank string, as len_trim doesn't always work on an uninitialized string
             Axes(diag_axis_id)%attributes(this_attribute)%catt = ''
          END IF
       END IF

       SELECT CASE (type)
       CASE (NF90_INT)
          IF ( .NOT.PRESENT(ival) ) THEN
             ! <ERROR STATUS="FATAL">
             !   Number type claims INTEGER, but ival not present for attribute <name> for axis <axis_name>.
             !   Contact the developers.
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_axis_add_attribute',&
                  & 'Attribute type claims INTEGER, but ival not present for attribute "'&
                  & //TRIM(name)//'" for axis "'//TRIM(Axes(diag_axis_id)%name)&
                  & //'". Contact the developers.', FATAL)
          END IF
          length = SIZE(ival)
          ! Allocate iatt(:) to size of ival
          ALLOCATE(Axes(diag_axis_id)%attributes(this_attribute)%iatt(length), STAT=istat)
          IF ( istat.NE.0 ) THEN
             ! <ERROR STATUS="FATAL">
             !   Unable to allocate iatt for attribute <name> for axis <axis_name>
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_axis_add_attribute', 'Unable to allocate iatt for attribute "'&
                  & //TRIM(name)//'" for axis "'//TRIM(Axes(diag_axis_id)%name)//'"', FATAL)
          END IF
          ! Set remaining fields
          Axes(diag_axis_id)%attributes(this_attribute)%len = length
          Axes(diag_axis_id)%attributes(this_attribute)%iatt = ival
       CASE (NF90_FLOAT)
          IF ( .NOT.PRESENT(rval) ) THEN
             ! <ERROR STATUS="FATAL">
             !   Attribute type claims REAL, but rval not present for attribute <name> for axis <axis_name>.
             !   Contact the developers.
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_axis_add_attribute',&
                  & 'Attribute type claims REAL, but rval not present for attribute "'&
                  & //TRIM(name)//'" for axis "'//TRIM(Axes(diag_axis_id)%name)&
                  & //'". Contact the developers.', FATAL)
          END IF
          length = SIZE(rval)
          ! Allocate iatt(:) to size of rval
          ALLOCATE(Axes(diag_axis_id)%attributes(this_attribute)%fatt(length), STAT=istat)
          IF ( istat.NE.0 ) THEN
             ! <ERROR STATUS="FATAL">
             !   Unable to allocate fatt for attribute <name> for axis <axis_name>
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_axis_add_attribute', 'Unable to allocate fatt for attribute "'&
                  & //TRIM(name)//'" for axis "'//TRIM(Axes(diag_axis_id)%name)&
                  & //'"', FATAL)
          END IF
          ! Set remaining fields
          Axes(diag_axis_id)%attributes(this_attribute)%len = length
          Axes(diag_axis_id)%attributes(this_attribute)%fatt = rval
       CASE (NF90_CHAR)
          IF ( .NOT.PRESENT(cval) ) THEN
             ! <ERROR STATUS="FATAL">
             !   Attribute type claims CHARACTER, but cval not present for attribute <name> for axis <axis_name>.
             !   Contact the developers.
             ! </ERROR>
             CALL error_mesg('diag_manager_mod::diag_axis_add_attribute',&
                  & 'Attribute type claims CHARACTER, but cval not present for attribute "'&
                  & //TRIM(name)//'" for axis "'//TRIM(Axes(diag_axis_id)%name)&
                  & //'". Contact the developers.', FATAL)
          END IF
          CALL prepend_attribute_axis(Axes(diag_axis_id), TRIM(name), TRIM(cval))
       CASE default
          ! <ERROR STATUS="FATAL">
          !   Unknown attribute type for attribute <name> for axis <axis_name>.
          !   Contact the developers.
          ! </ERROR>
          CALL error_mesg('diag_manager_mod::diag_axis_add_attribute', 'Unknown attribute type for attribute "'&
               & //TRIM(name)//'" for axis "'//TRIM(Axes(diag_axis_id)%name)&
               & //'". Contact the developers.', FATAL)
       END SELECT
    END IF
  END SUBROUTINE diag_axis_attribute_init

  SUBROUTINE diag_axis_add_attribute_scalar_r(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    REAL, INTENT(in) :: att_value

    CALL diag_axis_add_attribute_r1d(diag_axis_id, att_name, (/ att_value /))
  END SUBROUTINE diag_axis_add_attribute_scalar_r

  SUBROUTINE diag_axis_add_attribute_scalar_i(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    INTEGER, INTENT(in) :: att_value

    CALL diag_axis_add_attribute_i1d(diag_axis_id, att_name, (/ att_value /))
  END SUBROUTINE diag_axis_add_attribute_scalar_i

  SUBROUTINE diag_axis_add_attribute_scalar_c(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    CHARACTER(len=*), INTENT(in) :: att_value

    CALL diag_axis_attribute_init(diag_axis_id, att_name, NF90_CHAR, cval=att_value)
  END SUBROUTINE diag_axis_add_attribute_scalar_c

  SUBROUTINE diag_axis_add_attribute_r1d(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    REAL, DIMENSION(:), INTENT(in) :: att_value

    INTEGER :: num_attributes, len
    CHARACTER(len=512) :: err_msg

    CALL diag_axis_attribute_init(diag_axis_id, att_name, NF90_FLOAT, rval=att_value)
  END SUBROUTINE diag_axis_add_attribute_r1d

  SUBROUTINE diag_axis_add_attribute_i1d(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    INTEGER, DIMENSION(:), INTENT(in) :: att_value

    CALL diag_axis_attribute_init(diag_axis_id, att_name, NF90_INT, ival=att_value)
  END SUBROUTINE diag_axis_add_attribute_i1d

  !> @brief Allocates memory in out_file for the attributes.  Will <TT>FATAL</TT> if err_msg is not included
  !!   in the subroutine call.
  SUBROUTINE attribute_init_axis(out_axis, err_msg)
    TYPE(diag_axis_type), INTENT(inout) :: out_axis !< output file to allocate memory for attribute
    CHARACTER(LEN=*), INTENT(out), OPTIONAL :: err_msg !< Error message, passed back to calling function

    INTEGER :: istat

    ! Initialize err_msg
    IF ( PRESENT(err_msg) ) err_msg = ''

    ! Allocate memory for the attributes
    IF ( .NOT.allocated(out_axis%attributes) ) THEN
       ALLOCATE(out_axis%attributes(max_axis_attributes), STAT=istat)
       IF ( istat.NE.0 ) THEN
          ! <ERROR STATUS="FATAL">
          !   Unable to allocate memory for diag axis attributes
          ! </ERROR>
          IF ( fms_error_handler('diag_util_mod::attribute_init_axis', 'Unable to allocate memory for diag axis attributes', err_msg) ) THEN
             RETURN
          END IF
       ELSE
          ! Set equal to 0.  It will be increased when attributes added
          out_axis%num_attributes = 0
       END IF
    END IF
  END SUBROUTINE attribute_init_axis

  !> @brief Prepends the attribute value to an already existing attribute.  If the
  !!    attribute isn't yet defined, then creates a new attribute
  SUBROUTINE prepend_attribute_axis(out_axis, att_name, prepend_value, err_msg)
    TYPE(diag_axis_type), INTENT(inout) :: out_axis !< diagnostic axis that will get the attribute
    CHARACTER(len=*), INTENT(in) :: att_name !< Name of the attribute
    CHARACTER(len=*), INTENT(in) :: prepend_value !< Value to prepend
    CHARACTER(len=*), INTENT(out) , OPTIONAL :: err_msg !< Error message, passed back to calling routine

    INTEGER :: length, i, this_attribute
    CHARACTER(len=512) :: err_msg_local

    ! Initialize string variables
    err_msg_local = ''
    IF ( PRESENT(err_msg) ) err_msg = ''

    ! Make sure the attributes for this out file have been initialized
    CALL attribute_init_axis(out_axis, err_msg_local)
    IF ( TRIM(err_msg_local) .NE. '' ) THEN
       IF ( fms_error_handler('diag_util_mod::prepend_attribute_axis', TRIM(err_msg_local), err_msg) ) THEN
          RETURN
       END IF
    END IF

    ! Find if attribute exists
    this_attribute = 0
    DO i=1, out_axis%num_attributes
       IF ( TRIM(out_axis%attributes(i)%name) .EQ. TRIM(att_name) ) THEN
          this_attribute = i
          EXIT
       END IF
    END DO

    IF ( this_attribute > 0 ) THEN
       IF ( out_axis%attributes(this_attribute)%type .NE. NF90_CHAR ) THEN
          ! <ERROR STATUS="FATAL">
          !   Attribute <name> is not a character attribute.
          ! </ERROR>
          IF ( fms_error_handler('diag_util_mod::prepend_attribute_axis',&
               & 'Attribute "'//TRIM(att_name)//'" is not a character attribute.',&
               & err_msg) ) THEN
             RETURN
          END IF
       END IF
    ELSE
       ! Defining a new attribute
       ! Increase the number of file attributes
       this_attribute = out_axis%num_attributes + 1
       IF ( this_attribute .GT. max_axis_attributes ) THEN
          ! <ERROR STATUS="FATAL">
          !   Number of attributes exceeds max_axis_attributes for attribute <name>.
          !   Increase diag_manager_nml:max_axis_attributes.
          ! </ERROR>
          IF ( fms_error_handler('diag_util_mod::prepend_attribute_axis',&
               & 'Number of attributes exceeds max_axis_attributes for attribute "'&
               &//TRIM(att_name)//'".  Increase diag_manager_nml:max_axis_attributes.',&
               & err_msg) ) THEN
             RETURN
          END IF
       ELSE
          out_axis%num_attributes = this_attribute
          ! Set name and type
          out_axis%attributes(this_attribute)%name = att_name
          out_axis%attributes(this_attribute)%type = NF90_CHAR
          ! Initialize catt to a blank string, as len_trim doesn't always work on an uninitialized string
          out_axis%attributes(this_attribute)%catt = ''
       END IF
    END IF

    ! Only add string only if not already defined
    IF ( INDEX(TRIM(out_axis%attributes(this_attribute)%catt), TRIM(prepend_value)).EQ.0 ) THEN
       ! Check if new string length goes beyond the length of catt
       length = LEN_TRIM(TRIM(prepend_value)//" "//TRIM(out_axis%attributes(this_attribute)%catt))
       IF ( length.GT.LEN(out_axis%attributes(this_attribute)%catt) ) THEN
          ! <ERROR STATUS="FATAL">
          !   Prepend length for attribute <name> is longer than allowed.
          ! </ERROR>
          IF ( fms_error_handler('diag_util_mod::prepend_attribute_file',&
               & 'Prepend length for attribute "'//TRIM(att_name)//'" is longer than allowed.',&
               & err_msg) ) THEN
             RETURN
          END IF
       END IF
       ! Set files
       out_axis%attributes(this_attribute)%catt =&
            & TRIM(prepend_value)//' '//TRIM(out_axis%attributes(this_attribute)%catt)
       out_axis%attributes(this_attribute)%len = length
    END IF
  END SUBROUTINE prepend_attribute_axis

  !> @brief given an axis, returns TRUE if the axis uses compression-by-gathering: that is, if
  !!   this is an axis for fields on unstructured grid
  !! @return logical whether or not the axis uses compression-by-gathering
  logical function axis_is_compressed(id)
    integer, intent(in) :: id

    integer :: i

    CALL valid_id_check(id, 'axis_is_compressed')

    axis_is_compressed = .FALSE.
    if (.not.allocated(Axes(id)%attributes)) return
    do i = 1, Axes(id)%num_attributes
       if (trim(Axes(id)%attributes(i)%name)=='compress') then
          axis_is_compressed = .TRUE.
          return
       endif
    enddo
  end function axis_is_compressed


  !> @brief given an index of compressed-by-gathering axis, return an array of axes used in
  !!   compression. It is a fatal error to call it on axis that is not compressed
  subroutine get_compressed_axes_ids(id, r)
    integer, intent(in)  :: id
    integer, intent(out), allocatable :: r(:)

    integer iatt, k, k1, k2, n
    logical :: space

    character(*), parameter :: tag = 'get_compressed_axes_ids'

    CALL valid_id_check(id, tag)

    associate (axis=>Axes(id))
    if (.not.allocated(axis%attributes)) call error_mesg(tag, &
       'attempt to get compression dimensions from axis "'//trim(axis%name)//'" which is not compressed (does not have any attributes)', FATAL)

    iatt = 0
    do k = 1,axis%num_attributes
       if (trim(axis%attributes(k)%name)=='compress') then
          iatt = k; exit ! from loop
       endif
    enddo

    if (iatt == 0) call error_mesg(tag, &
       'attempt to get compression dimensions from axis "'//trim(axis%name)//&
       '" which is not compressed (does not have "compress" attributes).', FATAL)
    if (axis%attributes(iatt)%type/=NF90_CHAR) call error_mesg(tag, &
       'attempt to get compression dimensions from axis "'//trim(axis%name)//&
       '" but the axis attribute "compress" has incorrect type.', FATAL)

    ! parse the "compress" attribute
    ! calculate the number of compression axes
    space = .TRUE.; n=0
    do k = 1, len(axis%attributes(iatt)%catt)
       if (space.and.(axis%attributes(iatt)%catt(k:k)/=' ')) then
          n = n+1
       endif
       space = (axis%attributes(iatt)%catt(k:k)==' ')
    enddo

    allocate(r(n))
    ! make array of compression axes indices. Go from the last to the first to get the
    ! array in FORTRAN order: they are listed in "compress" attribute  C order (fastest
    ! dimension last)
    k2 = 0
    do k = n, 1, -1
       do k1 = k2+1, len(axis%attributes(iatt)%catt)
          if (axis%attributes(iatt)%catt(k1:k1)/=' ') exit
       enddo
       do k2 = k1+1, len(axis%attributes(iatt)%catt)
          if (axis%attributes(iatt)%catt(k2:k2)==' ') exit
       enddo
       r(k) = get_axis_num(axis%attributes(iatt)%catt(k1:k2),Axis_sets(axis%set))
       if (r(k)<=0) call error_mesg(tag, &
           'compression dimension "'//trim(axis%attributes(iatt)%catt(k1:k2))//&
           '" not found among the axes of set "'//trim(Axis_sets(axis%set))//'".', FATAL)
    enddo
    end associate ! axis
  end subroutine get_compressed_axes_ids
END MODULE diag_axis_mod
!> @}
! close documentation grouping
