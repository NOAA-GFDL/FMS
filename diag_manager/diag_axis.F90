#include <fms_platform.h>

MODULE diag_axis_mod
  ! <CONTACT EMAIL="seth.underwood@noaa.gov">
  !   Seth Underwood
  ! </CONTACT>

  ! <OVERVIEW> <TT>diag_axis_mod</TT> is an integral part
  !   of diag_manager_mod. It helps to create axis IDs
  !   that are used in register_diag_field.
  ! </OVERVIEW>

  ! <DESCRIPTION> Users first create axis ID by calling
  !   diag_axis_init, then use this axis ID in
  !   register_diag_field.
  ! </DESCRIPTION>

  USE mpp_domains_mod, ONLY: domain1d, domain2d, mpp_get_compute_domain,&
       & mpp_get_domain_components, null_domain1d, null_domain2d,&
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
       & get_axis_num, get_diag_axis_domain_name, diag_axis_add_attribute


  ! Module variables
  ! Parameters
  ! Include variable "version" to be written to log file.
#include<file_version.h>

  ! counter of number of axes defined
  INTEGER, DIMENSION(:), ALLOCATABLE :: num_subaxes
  INTEGER :: num_def_axes = 0

  ! storage for axis set names
  CHARACTER(len=128), DIMENSION(:), ALLOCATABLE, SAVE :: Axis_sets
  INTEGER :: num_axis_sets = 0

  ! ---- global storage for all defined axes ----
  TYPE(diag_axis_type), ALLOCATABLE, SAVE :: Axes(:)
  LOGICAL :: module_is_initialized = .FALSE.

  ! <INTERFACE NAME="diag_axis_add_attribute">
  !   <OVERVIEW>
  !     Add a attribute to the diag axis
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE diag_axis_add_attribute(diag_axis_id, att_name, att_value)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Add an arbitrary attribute and value to the diagnostic axis.  Any number
  !     of attributes can be added to a given axis.  All attribute addition must
  !     be done before first <TT>send_data</TT> call.
  !
  !     If a real or integer attribute is already defined, a FATAL error will be called.
  !     If a character attribute is already defined, then it will be prepended to the
  !     existing attribute value.
  !   </DESCRIPTION>
  !   <IN NAME="diag_axis_id" TYPE="INTEGER" />
  !   <IN NAME="att_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="att_value" TYPE="REAL|INTEGER|CHARACTER(len=*)" />
  INTERFACE diag_axis_add_attribute
     MODULE PROCEDURE diag_axis_add_attribute_scalar_r
     MODULE PROCEDURE diag_axis_add_attribute_scalar_i
     MODULE PROCEDURE diag_axis_add_attribute_scalar_c
     MODULE PROCEDURE diag_axis_add_attribute_r1d
     MODULE PROCEDURE diag_axis_add_attribute_i1d
  END INTERFACE diag_axis_add_attribute
  ! </INTERFACE>

CONTAINS

  ! <FUNCTION NAME="diag_axis_init">
  !   <OVERVIEW>
  !     Initialize the axis, and return the axis ID.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION diag_axis_init(name, data, units, cart_name, long_name,
  !           direction, set_name, edges, Domain, Domain2, aux, tile_count)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     <TT>diag_axis_init</TT> initializes an axis and returns the axis ID that
  !     is to be used with <TT>register_diag_field</TT>.  This function also
  !     increments the axis counter and fills in the axes
  !   </DESCRIPTION>
  !   <IN NAME="name" TYPE="CHARACTER(len=*)">Short name for axis</IN>
  !   <IN NAME="data" TYPE="REAL, DIMENSION(:)">Array of coordinate values</IN>
  !   <IN NAME="units" TYPE="CHARACTER(len=*)">Units for the axis</IN>
  !   <IN NAME="cart_name" TYPE="CHARACTER(len=*)">
  !     Cartesian axis ("X", "Y", "Z", "T")
  !   </IN>
  !   <IN NAME="direction" TYPE="INTEGER, OPTIONAL" DEFAULT="0">
  !     Indicates the direction of the axis:
  !     <UL>
  !       <LI>Up if +1</LI>
  !       <LI>Down if -1</LI>
  !       <LI>Neither up or down if 0</LI>
  !     </UL>
  !   </IN>
  !   <IN NAME="long_name" TYPE="CHARACTER(len=*), OPTIONAL" DEFAULT="name">
  !     Long name for the axis.
  !   </IN>
  !   <IN NAME="edges" TYPE="INTEGER, OPTIONAL">
  !     Axis ID for the previously defined "edges axis"
  !   </IN>
  !   <IN NAME="Domain" TYPE="TYPE(domain1d), OPTIONAL" />
  !   <IN NAME="Domain2" TYPE="TYPE(domain2d), OPTIONAL" />
  !   <IN NAME="aux" TYPE="CHARACTER(len=*), OPTIONAL">
  !     Auxiliary name, can only be <TT>geolon_t</TT> or <TT>geolat_t</TT>
  !   </IN>
  !   <IN NAME="tile_count" TYPE="INTEGER, OPTIONAL" />
  INTEGER FUNCTION diag_axis_init(name, DATA, units, cart_name, long_name, direction,&
       & set_name, edges, Domain, Domain2, aux, tile_count)
    CHARACTER(len=*), INTENT(in) :: name
    REAL, DIMENSION(:), INTENT(in) :: DATA
    CHARACTER(len=*), INTENT(in) :: units
    CHARACTER(len=*), INTENT(in) :: cart_name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: long_name, set_name
    INTEGER, INTENT(in), OPTIONAL :: direction, edges
    TYPE(domain1d), INTENT(in), OPTIONAL :: Domain
    TYPE(domain2d), INTENT(in), OPTIONAL :: Domain2
    CHARACTER(len=*), INTENT(in), OPTIONAL :: aux
    INTEGER, INTENT(in), OPTIONAL :: tile_count

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
         & TRIM(uppercase(cart_name)) == 'N' ) THEN
       Axes(diag_axis_init)%cart_name = TRIM(uppercase(cart_name))
    ELSE
       ! <ERROR STATUS="FATAL">Invalid cart_name name.</ERROR>
       CALL error_mesg('diag_axis_mod::diag_axis_init', 'Invalid cart_name name.', FATAL)
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

    !---- axis direction (-1, 0, or +1) ----
    IF ( PRESENT(direction) )THEN
       IF ( ABS(direction) /= 1 .AND. direction /= 0 )&
            ! <ERROR STATUS="FATAL">direction must be 0, +1, or -1</ERROR>
            & CALL error_mesg('diag_axis_mod::diag_axis_init', 'direction must be 0, +1 or -1', FATAL)
       Axes(diag_axis_init)%direction = direction
    ELSE
       Axes(diag_axis_init)%direction = 0
    END IF

    !---- domain2d type ----
    IF ( PRESENT(Domain2) .AND. PRESENT(Domain)) THEN
       ! <ERROR STATUS="FATAL">Presence of both Domain and Domain2 at the same time is prohibited</ERROR>
       CALL error_mesg('diag_axis_mod::diag_axis_init',&
            & 'Presence of both Domain and Domain2 at the same time is prohibited', FATAL)
    ELSE IF ( PRESENT(Domain2) .OR. PRESENT(Domain)) THEN
       IF ( Axes(diag_axis_init)%cart_name /= 'X' .AND. Axes(diag_axis_init)%cart_name /= 'Y') THEN
          ! <ERROR STATUS="FATAL">Domain must not be present for an axis which is not in the X or Y direction.</ERROR>
          CALL error_mesg('diag_axis_mod::diag_axis_init',&
               & 'Domain must not be present for an axis which is not in the X or Y direction', FATAL)
       END IF
    END IF

    Axes(diag_axis_init)%tile_count = tile

    IF ( PRESENT(Domain2) ) THEN
       Axes(diag_axis_init)%Domain2 = Domain2
       CALL mpp_get_domain_components(Domain2, domain_x, domain_y, tile_count=tile_count)
       IF ( Axes(diag_axis_init)%cart_name == 'X' ) Axes(diag_axis_init)%Domain = domain_x
       IF ( Axes(diag_axis_init)%cart_name == 'Y' ) Axes(diag_axis_init)%Domain = domain_y
    ELSE IF ( PRESENT(Domain)) THEN
       !---- domain1d type ----
       Axes(diag_axis_init)%Domain2 = null_domain2d ! needed since not 2-D domain
       Axes(diag_axis_init)%Domain = Domain
    ELSE
       Axes(diag_axis_init)%Domain2 = null_domain2d
       Axes(diag_axis_init)%Domain = null_domain1d
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
  ! </FUNCTION>

  ! <FUNCTION NAME="diag_subaxes_init">
  !   <OVERVIEW>
  !     Create a subaxis on a parent axis.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION diag_subaxes_init(axis, subdata, start_indx, end_indx,
  !           domain_1d, domain_2d)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Given the ID of a parent axis, create a subaxis and fill it with data,
  !     and return the ID of the corresponding subaxis.
  !
  !     The subaxis is defined on the parent axis from <TT>start_indx</TT>
  !     to <TT>end_indx</TT>.
  !   </DESCRIPTION>
  !   <IN NAME="axis" TYPE="INTEGER">ID of the parent axis</IN>
  !   <IN NAME="subdata" TYPE="REAL, DIMENSION(:)">Data of the subaxis</IN>
  !   <IN NAME="start_indx" TYPE="INTEGER">Start index of the subaxis</IN>
  !   <IN NAME="end_indx" TYPE="INTEGER">End index of the subaxis</IN>
  !   <IN NAME="domain_1d" TYPE="TYPE(domain1d), OPTIONAL" />
  !   <IN NAME="domain_2d" TYPE="TYPE(domain2d), OPTIONAL" />
  INTEGER FUNCTION diag_subaxes_init(axis, subdata, start_indx, end_indx, domain_2d)
    INTEGER, INTENT(in) :: axis
    REAL, DIMENSION(:), INTENT(in) :: subdata
    INTEGER, INTENT(in) :: start_indx
    INTEGER, INTENT(in) :: end_indx
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
  ! </FUNCTION>

  ! <SUBROUTINE NAME="get_diag_axis">
  !   <OVERVIEW>
  !     Return information about the axis with index ID
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_diag_axis(id, name, units, long_name, cart_name,
  !          direction, edges, Domain, data, num_attributes, attributes)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return information about the axis with index ID
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  !   <OUT NAME="name" TYPE="CHARACTER(len=*)">Short name for axis</OUT>
  !   <OUT NAME="units" TYPE="CHARACTER(len=*)">Units for axis</OUT>
  !   <OUT NAME="long_name" TYPE="CHARACTER(len=*)">Long name for axis</OUT>
  !   <OUT NAME="cart_name" TYPE="CHARACTER(len=*)">
  !     Cartesian axis ("x", "y", "z", "t").
  !   </OUT>
  !   <OUT NAME="direction" TYPE="INTEGER">
  !     Direction of data. (See <TT>diag_axis_init</TT> for a description of
  !     allowed values)
  !   </OUT>
  !   <OUT NAME="edges" TYPE="INTEGER">
  !     Axis ID for the previously defined "edges axis".
  !   </OUT>
  !   <OUT NAME="Domain" TYPE="TYPE(domain1d)" />
  !   <OUT NAME="data" TYPE="REAL, DIMENSION(:)">
  !     Array of coordinate values for this axis.
  !   </OUT>
  SUBROUTINE get_diag_axis(id, name, units, long_name, cart_name,&
       & direction, edges, Domain, DATA, num_attributes, attributes)
    CHARACTER(len=*), INTENT(out) :: name, units, long_name, cart_name
    INTEGER, INTENT(in) :: id
    TYPE(domain1d), INTENT(out) :: Domain
    INTEGER, INTENT(out) :: direction, edges
    REAL, DIMENSION(:), INTENT(out) :: DATA
    INTEGER, INTENT(out), OPTIONAL :: num_attributes
    TYPE(diag_atttype), ALLOCATABLE, DIMENSION(:), INTENT(out), OPTIONAL :: attributes

    INTEGER :: i, j, istat

    CALL valid_id_check(id, 'get_diag_axis')
    name      = Axes(id)%name
    units     = Axes(id)%units
    long_name = Axes(id)%long_name
    cart_name = Axes(id)%cart_name
    direction = Axes(id)%direction
    edges     = Axes(id)%edges
    Domain    = Axes(id)%Domain
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
       IF ( _ALLOCATED(Axes(id)%attributes) ) THEN
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
             IF ( _ALLOCATED(attributes(i)%fatt) ) THEN
                DEALLOCATE(attributes(i)%fatt)
             END IF
             IF ( _ALLOCATED(attributes(i)%iatt) ) THEN
                DEALLOCATE(attributes(i)%iatt)
             END IF

             ! Copy in attribute data
             attributes(i)%type = Axes(id)%attributes(i)%type
             attributes(i)%len = Axes(id)%attributes(i)%len
             attributes(i)%name = Axes(id)%attributes(i)%name
             attributes(i)%catt = Axes(id)%attributes(i)%catt
             ! Allocate fatt arrays (if needed), and copy in data
             IF ( _ALLOCATED(Axes(id)%attributes(i)%fatt) ) THEN
                ALLOCATE(attributes(i)%fatt(SIZE(Axes(id)%attributes(i)%fatt(:))), STAT=istat)
                IF ( istat .NE. 0 ) THEN
                   CALL error_mesg('diag_axis_mod::get_diag_axis', 'Unable to allocate memory for attribute%fatt', FATAL)
                END IF
                DO j=1, SIZE(attributes(i)%fatt(:))
                   attributes(i)%fatt(j) = Axes(id)%attributes(i)%fatt(j)
                END DO
             END IF
             ! Allocate iatt arrays (if needed), and copy in data
             IF ( _ALLOCATED(Axes(id)%attributes(i)%iatt) ) THEN
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
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="get_diag_axis_cart">
  !   <OVERVIEW>
  !     Return the axis cartesian.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_diag_axis_cart(id, cart_name)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the axis cartesian ('X', 'Y', 'Z' or 'T') for the axis ID given.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  !   <OUT NAME="cart_name" TYPE="CHARACTER(len=*)">Cartesian axis</OUT>
  SUBROUTINE get_diag_axis_cart(id, cart_name)
    INTEGER, INTENT(in)           :: id
    CHARACTER(len=*), INTENT(out) :: cart_name

    CALL valid_id_check(id, 'get_diag_axis_cart')
    cart_name = Axes(id)%cart_name
  END SUBROUTINE get_diag_axis_cart
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="get_diag_axis_data">
  !   <OVERVIEW>
  !     Return the axis data.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_diag_axis_data(id, data)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the axis data for the axis ID given.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  !   <OUT NAME="data" TYPE="REAL, DIMENSION(:)">Axis data</OUT>
  SUBROUTINE get_diag_axis_data(id, DATA)
    INTEGER, INTENT(in) :: id
    REAL, DIMENSION(:), INTENT(out) :: DATA

    CALL valid_id_check(id, 'get_diag_axis_data')
    IF (Axes(id)%length > SIZE(DATA(:))) THEN
       ! <ERROR STATUS="FATAL">array data is too small</ERROR>
       CALL error_mesg('diag_axis_mod::get_diag_axis_data', 'array data is too small', FATAL)
    ELSE
       DATA(1:Axes(id)%length) = Axes(id)%data
    END IF
  END SUBROUTINE get_diag_axis_data
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="get_diag_axis_name">
  !   <OVERVIEW>
  !     Return the short name of the axis.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_diag_axis_name (id, name)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the short name for the axis ID given.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  !   <OUT NAME="name" TYPE="CHARACTER(len=*)">Axis short name</OUT>
  SUBROUTINE get_diag_axis_name(id, name)
    INTEGER         , INTENT(in)  :: id
    CHARACTER(len=*), INTENT(out) :: name

    CALL valid_id_check(id, 'get_diag_axis_name')
    name = Axes(id)%name
  END SUBROUTINE get_diag_axis_name
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="get_diag_axis_domain_name">
  !   <OVERVIEW>
  !     Return the name of the axis' domain
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_diag_axis_domain_name(id, name)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Retruns the name of the axis' domain.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  !   <OUT NAME="name" TYPE="CHARACTER(len=*)">Axis' domain name</OUT>
  SUBROUTINE get_diag_axis_domain_name(id, name)
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(out) :: name

    CALL valid_id_check(id, 'get_diag_axis_domain_name')
    name = mpp_get_domain_name(Axes(id)%domain2)
  END SUBROUTINE get_diag_axis_domain_name
  ! </SUBROUTINE>

  ! <FUNCTION NAME="get_axis_length">
  !   <OVERVIEW>
  !     Return the length of the axis.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION get_axis_length(id)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the length of the axis ID given.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  INTEGER FUNCTION get_axis_length(id)
    INTEGER, INTENT(in) :: id

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
  ! </FUNCTION>

  ! <FUNCTION NAME="get_axis_aux">
  !   <OVERVIEW>
  !     Return the auxiliary name for the axis.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     CHARACTER(len=128) FUNCTION get_axis_aux(id)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Returns the auxiliary name for the axis.  The only possible values for
  !     the auxiliary names is <TT>geolon_t</TT> or <TT>geolat_t</TT>.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  CHARACTER(len=138) FUNCTION get_axis_aux(id)
    INTEGER, INTENT(in) :: id

    CALL valid_id_check(id, 'get_axis_aux')
    get_axis_aux =  Axes(id)%aux
  END FUNCTION get_axis_aux
  ! </FUNCTION>

  ! <FUNCTION NAME="get_axis_global_length">
  !   <OVERVIEW>
  !     Return the global length of the axis.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION get_axis_global_length (id)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Returns the global length of the axis ID given.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  INTEGER FUNCTION get_axis_global_length(id)
    INTEGER, INTENT(in) :: id

    CALL valid_id_check(id, 'get_axis_global_length')
    get_axis_global_length = Axes(id)%length
  END FUNCTION get_axis_global_length
  ! </FUNCTION>

  ! <FUNCTION NAME="get_tile_count">
  !   <OVERVIEW>
  !     Return the tile count for the axis.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION get_tile_count (ids)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the tile count for the axis IDs given.
  !   </DESCRIPTION>
  !   <IN NAME="ids" TYPE="INTEGER, DIMENSION(:)">
  !     Axis IDs.  Possible dimensions: 1 <= <TT>size(ids(:))</TT> <= 4.
  !   </IN>
  INTEGER FUNCTION get_tile_count(ids)
    INTEGER, DIMENSION(:), INTENT(in) :: ids

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
  ! </FUNCTION>

  ! <FUNCTION NAME="get_domain1d">
  !   <OVERVIEW>
  !     Return the 1D domain.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     TYPE(domain1d) FUNCTION get_domain1d(id)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Retrun the 1D domain for the axis ID given.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis ID</IN>
  TYPE(domain1d) FUNCTION get_domain1d(id)
    INTEGER, INTENT(in) :: id

    CALL valid_id_check(id, 'get_domain1d')
    IF (Axes(id)%Domain .NE. NULL_DOMAIN1D) THEN
       get_domain1d = Axes(id)%Domain
    ELSE
       get_domain1d = NULL_DOMAIN1D
    ENDIF
  END FUNCTION get_domain1d
  ! </FUNCTION>

  ! <FUNCTION NAME="get_domain2d">
  !   <OVERVIEW>
  !     Return the 2D domain.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     TYPE(domain2d) FUNCTION get_domain2d(ids)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the 2D domain for the axis IDs given.
  !   </DESCRIPTION>
  !   <IN NAME="ids" TYPE="INTEGER, DIMENSION(:)">
  !     Axis IDs.  Possible dimensions: 1 <= <TT>size(ids(:))</TT> <= 4.
  !   </IN>
  TYPE(domain2d) FUNCTION get_domain2d(ids)
    INTEGER, DIMENSION(:), INTENT(in) :: ids

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
  ! </FUNCTION>

  ! <SUBROUTINE NAME="get_axes_shift">
  !   <OVERVIEW>
  !     Return the value of the shift.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE get_axes_shift(ids, ishift, jshift)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Return the value of the shift for the axis IDs given.
  !   </DESCRIPTION>
  !   <IN NAME="ids" TYPE="INTEGER, DIMENSION(:)">
  !     Axis IDs.  Possible dimensions: 1 <= <TT>size(ids(:))</TT> <= 4
  !   </IN>
  !   <OUT NAME="ishift" TYPE="INTEGER">X shift value.</OUT>
  !   <OUT NAME="jshift" TYPE="INTEGER">Y shift value.</OUT>
  SUBROUTINE get_axes_shift(ids, ishift, jshift)
    INTEGER, DIMENSION(:), INTENT(in) :: ids
    INTEGER, INTENT(out) :: ishift, jshift

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
  ! </SUBROUTINE>

  ! <PRIVATE>
  ! <FUNCTION NAME="get_axis_num">
  !   <OVERVIEW>
  !     Returns index into axis table corresponding to a given axis name.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION get_axis_num(axis_name, set_name)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Returns index into axis table corresponding to a given axis name.
  !   </DESCRIPTION>
  !   <IN NAME="axis_name" TYPE="CHARACTER(len=*)">Axis name.</IN>
  !   <IN NAME="set_name" TYPE="CHARACTER(len=*), OPTIONAL">Set name.</IN>
  INTEGER FUNCTION get_axis_num(axis_name, set_name)
    CHARACTER(len=*), INTENT(in) :: axis_name
    CHARACTER(len=*), INTENT(in), OPTIONAL :: set_name

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
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <FUNCTION NAME="get_axis_set_num">
  !   <OVERVIEW>
  !     Returns index in axis set table corresponding to a given axis set name.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     INTEGER FUNCTION get_axis_set_num(set_name)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Returns index in axis set table corresponding to a given axis set name.
  !   </DESCRIPTION>
  !   <IN NAME="set_name" TYPE="CHARACTER(len=*)">Set name.</IN>
  INTEGER FUNCTION get_axis_set_num (set_name)
    CHARACTER(len=*), INTENT(in) :: set_name

    INTEGER :: iset

    get_axis_set_num = 0
    DO iset = 1, num_axis_sets
       IF ( set_name == Axis_sets(iset) ) THEN
          get_axis_set_num = iset
          RETURN
       END IF
    END DO
  END FUNCTION get_axis_set_num
  ! </FUNCTION>
  ! </PRIVATE>

  ! <PRIVATE>
  ! <SUBROUTINE NAME="valid_id_check">
  !   <OVERVIEW>
  !     Check to see if the axis id is a vaild id.
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE valid_id_check(id, routine_name)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Check to see if the given axis id is a valid id.  If the axis id is invalid,
  !     call a FATAL error.  If the ID is valid, just return.
  !   </DESCRIPTION>
  !   <IN NAME="id" TYPE="INTEGER">Axis id to check for validity</IN>
  !   <IN NAME="routine_name" TYPE="CHARACTER(len=*)">Name of the subroutine checking for a valid axis id.</IN>
  SUBROUTINE valid_id_check(id, routine_name)
    INTEGER, INTENT(in) :: id
    CHARACTER(len=*), INTENT(in) :: routine_name

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
  ! </SUBROUTINE>
  ! </PRIVATE>

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

  ! <SUBROUTINE NAME="diag_axis_add_attribute_scalar_r" INTERFACE="diag_axis_add_attribute">
  !   <IN NAME="diag_axis_id" TYPE="INTEGER" />
  !   <IN NAME="att_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="att_value" TYPE="REAL" />
  SUBROUTINE diag_axis_add_attribute_scalar_r(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    REAL, INTENT(in) :: att_value

    CALL diag_axis_add_attribute_r1d(diag_axis_id, att_name, (/ att_value /))
  END SUBROUTINE diag_axis_add_attribute_scalar_r
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_axis_add_attribute_scalar_i" INTERFACE="diag_axis_add_attribute">
  !   <IN NAME="diag_axis_id" TYPE="INTEGER" />
  !   <IN NAME="att_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="att_value" TYPE="INTEGER" />
  SUBROUTINE diag_axis_add_attribute_scalar_i(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    INTEGER, INTENT(in) :: att_value

    CALL diag_axis_add_attribute_i1d(diag_axis_id, att_name, (/ att_value /))
  END SUBROUTINE diag_axis_add_attribute_scalar_i
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_axis_add_attribute_scalar_c" INTERFACE="diag_axis_add_attribute">
  !   <IN NAME="diag_axis_id" TYPE="INTEGER" />
  !   <IN NAME="att_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="att_value" TYPE="CHARACTER(len=*)" />
  SUBROUTINE diag_axis_add_attribute_scalar_c(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    CHARACTER(len=*), INTENT(in) :: att_value

    CALL diag_axis_attribute_init(diag_axis_id, att_name, NF90_CHAR, cval=att_value)
  END SUBROUTINE diag_axis_add_attribute_scalar_c
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_axis_add_attribute_r1d" INTERFACE="diag_axis_add_attribute">
  !   <IN NAME="diag_axis_id" TYPE="INTEGER" />
  !   <IN NAME="att_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="att_value" TYPE="REAL, DIMENSION(:)" />
  SUBROUTINE diag_axis_add_attribute_r1d(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    REAL, DIMENSION(:), INTENT(in) :: att_value

    INTEGER :: num_attributes, len
    CHARACTER(len=512) :: err_msg

    CALL diag_axis_attribute_init(diag_axis_id, att_name, NF90_FLOAT, rval=att_value)
  END SUBROUTINE diag_axis_add_attribute_r1d
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="diag_axis_add_attribute_i1d" INTERFACE="diag_axis_add_attribute">
  !   <IN NAME="diag_axis_id" TYPE="INTEGER" />
  !   <IN NAME="att_name" TYPE="CHARACTER(len=*)" />
  !   <IN NAME="att_value" TYPE="INTEGER, DIMENSION(:)" />
  SUBROUTINE diag_axis_add_attribute_i1d(diag_axis_id, att_name, att_value)
    INTEGER, INTENT(in) :: diag_axis_id
    CHARACTER(len=*), INTENT(in) :: att_name
    INTEGER, DIMENSION(:), INTENT(in) :: att_value

    CALL diag_axis_attribute_init(diag_axis_id, att_name, NF90_INT, ival=att_value)
  END SUBROUTINE diag_axis_add_attribute_i1d
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="attribute_init_axis" INTERFACE="attribute_init">
  !   <OVERVIEW>
  !     Allocates the atttype in out_axis
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE attribute_init(out_axis, err_msg)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Allocates memory in out_file for the attributes.  Will <TT>FATAL</TT> if err_msg is not included
  !     in the subroutine call.
  !   </DESCRIPTION>
  !   <INOUT NAME="out_file" TYPE="TYPE(file_type)">output file to allocate memory for attribute</INOUT>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">Error message, passed back to calling function</OUT>
  SUBROUTINE attribute_init_axis(out_axis, err_msg)
    TYPE(diag_axis_type), INTENT(inout) :: out_axis
    CHARACTER(LEN=*), INTENT(out), OPTIONAL :: err_msg

    INTEGER :: istat

    ! Initialize err_msg
    IF ( PRESENT(err_msg) ) err_msg = ''

    ! Allocate memory for the attributes
    IF ( .NOT._ALLOCATED(out_axis%attributes) ) THEN
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
  ! </SUBROUTINE>

  ! <SUBROUTINE NAME="prepend_attribute_axis">
  !   <OVERVIEW>
  !     Prepends the attribute value to an already existing attribute.  If the
  !     attribute isn't yet defined, then creates a new attribute
  !   </OVERVIEW>
  !   <TEMPLATE>
  !     SUBROUTINE prepend_attribute(out_file, attribute, prepend_value)
  !   </TEMPLATE>
  !   <DESCRIPTION>
  !     Checks if the attribute already exists for <TT>out_axis</TT>.  If it exists,
  !     then prepend the <TT>prepend_value</TT>, otherwise, create the attribute
  !     with the <TT>prepend_value</TT>. <TT>err_msg</TT> indicates no duplicates found.
  !   </DESCRIPTION>
  !   <INOUT NAME="out_axis" TYPE="TYPE(diag_axis_type)">diagnostic axis that will get the attribute</INOUT>
  !   <IN NAME="att_name" TYPE="CHARACTER(len=*)">Name of the attribute</IN>
  !   <IN NAME="prepend_value" TYPE="CHARACTER(len=*)">Value to prepend</IN>
  !   <OUT NAME="err_msg" TYPE="CHARACTER(len=*), OPTIONAL">Error message, passed back to calling routine</OUT>
  SUBROUTINE prepend_attribute_axis(out_axis, att_name, prepend_value, err_msg)
    TYPE(diag_axis_type), INTENT(inout) :: out_axis
    CHARACTER(len=*), INTENT(in) :: att_name, prepend_value
    CHARACTER(len=*), INTENT(out) , OPTIONAL :: err_msg

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
  ! </SUBROUTINE>
END MODULE diag_axis_mod
