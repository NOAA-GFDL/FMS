module diag_axis_mod
! <CONTACT EMAIL="Giang.Nong@noaa.gov">
!   Giang Nong
! </CONTACT>

! <OVERVIEW>
!   <TT>diag_axis_mod</TT> is an integral part of diag_manager_mod. It helps to create
! axis IDs that are used in register_diag_field.
! </OVERVIEW>

! <DESCRIPTION>
! Users first create axis ID by calling diag_axis_init, then use this axis ID in 
! register_diag_field.
!
!</DESCRIPTION>

use mpp_domains_mod, only: domain1d, domain2d, mpp_get_compute_domain, &
                           mpp_get_domain_components, null_domain1d, &
                           null_domain2d, operator(/=), mpp_get_global_domain
use         fms_mod, only: error_mesg, write_version_number, lowercase, FATAL
use diag_data_mod, only  : diag_axis_type, max_subaxes, max_axes, max_num_axis_sets
implicit none

private
public  diag_axis_init, get_diag_axis, get_domain1d, get_domain2d, &
        get_axis_length, get_axis_global_length, diag_subaxes_init, &
        get_diag_axis_cart, get_diag_axis_data, max_axes, get_axis_aux, &
        get_tile_count, get_axes_shift, get_diag_axis_name



!counter of number of axes defined
integer,allocatable :: num_subaxes(:)
integer             :: num_def_axes = 0

!-----------------------------------------------------------------------
!
!  name        = short name for axis
!  units       = units for axis
!  long_name   = long name for axis
!  cart_name   = cartesian axis ('x','y','z', or 't')
!  length      = number of coordinates for axis
!  direction   = if +1, data are in a positive direction (default:+1)
!                if -1, data are in a negative direction
!  data        = array of coordinate values for this axis
!
!-----------------------------------------------------------------------

! ---- storage for axis set names ----

integer            :: num_axis_sets = 0
character(len=128), allocatable, save :: Axis_sets(:)
! ---- global storage for all defined axes ----
type (diag_axis_type), allocatable, save :: Axes(:)
logical            :: module_is_initialized = .FALSE.
character(len=128) :: &
     version='$Id: diag_axis.F90,v 16.0 2008/07/30 22:44:59 fms Exp $'
character(len=128) :: tagname='$Name: perth $'

contains
!#######################################################################

function diag_axis_init (name, data, units, cart_name, long_name,     &
       direction, set_name, edges, Domain, Domain2, aux, tile_count) &
       result (indexx)

! increment axis counter and fill in axes     
!-----------------------------------------------------------------------
!  name                = short name for axis
!  data                = array of coordinate values for this axis
!  units               = units for axis
!  cart_name           = cartesian axis ("x",'y','z','t')
!  direction(optional) = if +1, data are in a up   direction
!                      = if -1, data are in a down direction
!                      = if  0, neither up or down (default)
!  long_name(optional) = long name for axis (default: name)
!  edges    (optional) = axis id for the previously defined "edges axis"
!  aux      (optional) = auxiliary name, can be only either geolon_t or geolat_t
!-----------------------------------------------------------------------

  character(len=*), intent(in)           :: name
  real            , intent(in)           :: data(:)
  character(len=*), intent(in)           :: units
  character(len=*), intent(in)           :: cart_name  
  character(len=*), intent(in), optional :: long_name, set_name
  integer         , intent(in), optional :: direction, edges
  type(domain1d)  , intent(in), optional :: Domain
  type(domain2d)  , intent(in), optional :: Domain2
  character(len=*), intent(in), optional :: aux
  integer         , intent(in), optional :: tile_count
  type(domain1d)                         :: domain_x, domain_y
  integer                                :: indexx, ierr, axlen
  integer                                :: i, set
  integer                                :: isc, iec, isg, ieg

  if ( .not.module_is_initialized ) then
     call write_version_number( version, tagname )
  endif
  if (.not. allocated(Axis_sets)) allocate(Axis_sets(max_num_axis_sets))
  if (.not. allocated(Axes)) allocate(Axes(max_axes))
  if (.not. allocated(num_subaxes)) then
     allocate(num_subaxes(max_axes))
     num_subaxes = 0
  endif
!---- is there an axis set? ----
  if ( present(set_name) ) then
     set = get_axis_set_num (set_name)
!    ---- add new set name ----
     if (set == 0) then
        num_axis_sets = num_axis_sets + 1
        if (num_axis_sets > max_num_axis_sets)   &
             call error_mesg('diag_axis_init in diag_axis_mod',  &
             'exceeded max_num_axis_sets', FATAL)
        set = num_axis_sets
        Axis_sets(set) = set_name
     endif
  else
     set = 0
  endif

!---- see if axis already exists --
! if this is time axis, return the ID of a previously defined
! if this is spatial axis, FATAL error
  do i = 1, num_def_axes
     if (trim(name) == Axes(i)%name .and. trim(name) == 'Stations') then
        indexx = i
        return
     endif
     if (trim(name) == Axes(i)%name .and. trim(name) == 'Levels') then
        indexx = i
        return
     endif
     if (name == Axes(i)%name .and. set == Axes(i)%set) then
        if(trim(lowercase(name))=='time'.or.trim(lowercase(cart_name))=='t') then
           indexx = i
           return
        else
           call error_mesg ('diag_axis_init in diag_axis_mod', &
                'axis_name '//trim(name)//' and axis_set already exists', FATAL)
        endif
     endif
  enddo

!---- register axis ----
  num_def_axes = num_def_axes + 1
  if (num_def_axes > max_axes) call error_mesg ('diag_axis_init in diag_axis_mod', &
       'max_axes exceeded, increase it via diag_manager_nml', FATAL)
  indexx = num_def_axes

!---- check and then save cart_name name ----
  if (cart_name == 'x' .or.cart_name == 'X' ) then
     Axes(indexx)%cart_name = 'X'
  else if (cart_name == 'y' .or.cart_name == 'Y' ) then
     Axes(indexx)%cart_name = 'Y'
  else if (cart_name == 'z' .or.cart_name == 'Z' ) then
     Axes(indexx)%cart_name = 'Z'
  else if (cart_name == 't' .or.cart_name == 'T' ) then
      Axes(indexx)%cart_name = 'T'
  else if (cart_name == 'n' .or.cart_name == 'N' ) then
      Axes(indexx)%cart_name = 'N' 
  else     
      call error_mesg('diag_axis_init in diag_axis_mod', 'Invalid cart_name name.', FATAL)
  endif

!---- allocate storage for coordinate values of axis ----
  axlen = size(data(:))
  if ( Axes(indexx)%cart_name == 'T' ) axlen = 0
  allocate ( Axes(indexx)%data(1:axlen) )
  Axes(indexx)%name   = trim(name)
  Axes(indexx)%data   = data(1:axlen)
  Axes(indexx)%units  = units  
  Axes(indexx)%length = axlen
  Axes(indexx)%set    = set
! start and end are used in subaxes information only
  Axes(indexx)%start = -1
  Axes(indexx)%end = -1
  Axes(indexx)%subaxis_name = ""
  Axes(indexx)%shift = 0

  if (present(long_name))then
     Axes(indexx)%long_name = long_name
  else
     Axes(indexx)%long_name = name
  endif
  if (present(aux))then
     Axes(indexx)%aux = trim(aux)
  else
     Axes(indexx)%aux = 'none'
  endif
 
!---- axis direction (-1, 0, or +1) ----
  if (present(direction))then
     if(abs(direction) /= 1 .AND. direction /= 0) call error_mesg('diag_axis_init in diag_axis_mod', &
          'direction must be 0, +1 or -1',FATAL)
     Axes(indexx)%direction = direction
  else
     Axes(indexx)%direction = 0
  endif
!---- domain2d type ----
  if ( present(Domain2) .or. present(Domain)) then
    if ( Axes(indexx)%cart_name /= 'X' .and. Axes(indexx)%cart_name /= 'Y') then
      call error_mesg('diag_axis_init in diag_axis_mod', &
           'Domain must not be present for an axis which is not in the X or Y direction', &
            FATAL)
    endif
  endif
  if ( present(Domain2) .and. present(Domain)) then
    call error_mesg('diag_axis_init in diag_axis_mod', &
           'Presence of both Domain and Domain2 at the same time is prohibited', &
            FATAL)
  endif

  Axes(indexx)%tile_count = 1
  if(present(tile_count)) Axes(indexx)%tile_count = tile_count

  if ( present(Domain2) ) then
     Axes(indexx)%Domain2 = Domain2
     call mpp_get_domain_components(Domain2, domain_x, domain_y, tile_count=tile_count)
     if ( Axes(indexx)%cart_name == 'X' ) Axes(indexx)%Domain = domain_x
     if ( Axes(indexx)%cart_name == 'Y' ) Axes(indexx)%Domain = domain_y
  else
     Axes(indexx)%Domain2 = null_domain2d 
!---- domain1d type ----     
     if ( present(Domain)) then
        Axes(indexx)%Domain = Domain
     else
        Axes(indexx)%Domain = null_domain1d
     endif
  endif

  !--- set up the shift value for x-y axis
  if(Axes(indexx)%Domain .ne. null_domain1d ) then
     call mpp_get_compute_domain(Axes(indexx)%Domain, isc, iec)
     call mpp_get_global_domain(Axes(indexx)%Domain, isg, ieg)
     if(  Axes(indexx)%length == ieg - isg + 2 ) then
        Axes(indexx)%shift = 1 
     endif
  endif

!---- have axis edges been defined ? ----
  Axes(indexx)%edges = 0
  if (present(edges))then
     if ( edges > 0 .and. edges < num_def_axes ) then
        ierr=0
        if ( Axes(edges)%cart_name /= Axes(indexx)%cart_name) ierr=1
        if ( Axes(edges)%length    /= Axes(indexx)%length+1 ) ierr=ierr+2
        if ( Axes(edges)%set       /= Axes(indexx)%set      ) ierr=ierr+4
        if ( ierr > 0 )   call error_mesg ('diag_axis_init in diag_axis_mod', &
             'Edges axis does not match axis', FATAL)
        Axes(indexx)%edges = edges
     else
       call error_mesg ('diag_axis_init in diag_axis_mod', &
                        'Edges axis is not defined', FATAL)
     endif
  endif
  module_is_initialized = .TRUE.
!-----------------------------------------------------------------------
end function diag_axis_init

function diag_subaxes_init(axis,subdata,start_indx,end_indx,domain_1d,domain_2d)  result(index)
! Given ID of parent axis and data of subaxis, this function returns ID of corresponding subaxis 
! subaxis is defined on parent axis from start_indx to end_indx
  integer, intent(in)            :: axis ! ID of parent axis
  real, dimension(:), intent(in) :: subdata ! data of subaxis
  integer, intent(in)            :: start_indx ! start index of subaxis
  integer, intent(in)            :: end_indx ! end index of subaxis
  type(domain1d), intent(in), optional  :: domain_1d
  type(domain2d), intent(in), optional  :: domain_2d

  integer                        :: index
  integer                        :: i,nsub_axis, direction
  character(len=128)             :: name, nsub_name   
  character(len=128)             :: units
  character(len=128)             :: cart_name
  character(len=128)             :: long_name
  logical                        :: subaxis_set 

! there may be more than 1 subaxis on a parent axis, check for redundancy
  nsub_axis = 0; subaxis_set = .false.
  do i = 1, num_subaxes(axis)
     if(start_indx == Axes(axis)%start(i) .and. end_indx == Axes(axis)%end(i)) then
        nsub_axis = i
        subaxis_set = .true.    !subaxis already exists
        name = trim(Axes(axis)%subaxis_name(nsub_axis))
        exit
     endif
  enddo
  if (nsub_axis == 0) then  ! create new subaxis
     num_subaxes(axis) = num_subaxes(axis) + 1
     if (num_subaxes(axis) > max_subaxes) &
          call error_mesg ('diag_subaxes_init in diag_axis_mod',' max_subaxes too small, increase max_subaxes', FATAL)
     nsub_axis = num_subaxes(axis)
     Axes(axis)%start(nsub_axis) = start_indx
     Axes(axis)%end(nsub_axis)   = end_indx
  endif
  
  ! Create new name for the subaxis from name of parent axis

! If subaxis already exists, get the index and return       
  if(subaxis_set) then
     if (Axes(axis)%set > 0) then
        index = get_axis_num(name,set_name=trim(Axis_sets(Axes(axis)%set)))     
     else
        index = get_axis_num(name)     
     endif
  else
! get a new index for subaxis
     write(nsub_name,'(I1)') nsub_axis
     name = trim(Axes(axis)%name)//'_sub'//trim(nsub_name)
     Axes(axis)%subaxis_name(nsub_axis) = name
     long_name = trim(Axes(axis)%long_name)
     units = trim(Axes(axis)%units)
     cart_name = trim(Axes(axis)%cart_name)
     direction = Axes(axis)%direction
     if (Axes(axis)%set > 0) then
        index =  diag_axis_init (trim(name), subdata, trim(units), trim(cart_name), trim(long_name), &
             set_name=trim(Axis_sets(Axes(axis)%set)), Domain2=domain_2d)
     else
        index =  diag_axis_init (trim(name), subdata, trim(units), trim(cart_name), trim(long_name), &
             Domain2=domain_2d)
     endif
  endif        
end function diag_subaxes_init
!#######################################################################
         
subroutine get_diag_axis (id, name, units, long_name, cart_name, &
                           direction, edges, Domain, data)

! Return information about the axis with index id

!-----------------------------------------------------------------------
!  id         =  axis number
!  name       = short name for axis
!  units      = units for axis
!  long_name  = long name for axis
!  cart_name  = cartesian axis ("x",'y','z', or 't')
!  direction  = if +1, data are in a positive direction (default:+1)
!               if -1, data are in a negative direction
!  data       = array of coordinate values for this axis
!-----------------------------------------------------------------------

  character(len=*), intent(out) :: name, units, long_name, cart_name
  integer, intent(in)           :: id
  type(domain1d), intent(out)   :: Domain
  integer, intent(out)          :: direction, edges
  real, intent(out)             :: data(:)
  character(len=128)            :: error_msg

  if (id < 1 .or. id > num_def_axes) then
     write(error_msg,'(i2)')id
     call error_mesg('get_diag_axis in diag_axis_mod', &
          trim(error_msg)//' is illegal value for axis_id', FATAL)
  endif
  name      = Axes(id)%name
  units     = Axes(id)%units
  long_name = Axes(id)%long_name
  cart_name = Axes(id)%cart_name
  direction = Axes(id)%direction
  edges     = Axes(id)%edges
  Domain    = Axes(id)%Domain
  if (Axes(id)%length > size(data(:))) call error_mesg ('get_diag_axis in diag_axis_mod', &
       'array data is too small', FATAL)
  data(1:Axes(id)%length) = Axes(id)%data

end subroutine get_diag_axis
!#####################################################################

subroutine get_diag_axis_cart(id, cart_name)
! Return the axis cartesian from  index id
!  id         =  axis number
!  cart_name  = cartesian axis ("x",'y','z', or 't')
  character(len=*), intent(out) :: cart_name
  integer, intent(in)           :: id
  cart_name = Axes(id)%cart_name
end subroutine get_diag_axis_cart

!######################################################################
subroutine get_diag_axis_data(id,data)
! Return the axis data from  index id
!  id         =  axis number
  integer, intent(in) :: id
  real, intent(out)   :: data(:)
  if (Axes(id)%length > size(data(:))) call error_mesg ('get_diag_axis_data in diag_axis_mod', &
       'array data is too small', FATAL)
  data(1:Axes(id)%length) = Axes(id)%data
end subroutine get_diag_axis_data


!#######################################################################
subroutine get_diag_axis_name (id, name)
  integer         , intent(in)  :: id
  character(len=*), intent(out) :: name

  character(len=128) :: error_msg

  if (id < 1 .or. id > num_def_axes) then
     write(error_msg,'(i2)')id
     call error_mesg('get_diag_axis_name in diag_axis_mod', &
          trim(error_msg)//' is illegal value for axis_id', FATAL)
  endif
  name      = Axes(id)%name
end subroutine get_diag_axis_name

!######################################################################

function get_axis_length (id) result (length)
  integer, intent(in) :: id
  integer             :: length   
  if ( Axes(id)%Domain /= null_domain1d ) then
     call mpp_get_compute_domain(Axes(id)%Domain,size=length)
     !---one extra point is needed for some case. ( like symmetry domain )
     length = length + Axes(id)%shift
  else
     length = Axes(id) % length
  endif
end function get_axis_length
!######################################################################
function get_axis_aux (id) result (aux)

  integer, intent(in)   :: id
  character(len=128)    :: aux

  aux =  Axes(id)%aux
end function get_axis_aux

!#######################################################################

function get_axis_global_length (id) result (length)
  integer, intent(in) :: id
  integer             :: length
  length = Axes(id) % length
end function get_axis_global_length
!#######################################################################

function get_tile_count (ids) result (tile_count)
  integer, intent(in) :: ids(:)
  integer             :: tile_count
  integer             :: i, id, flag

  if ( size(ids(:)) < 1 .or. size(ids(:)) > 4 ) call error_mesg  &
       ('get_tile_count in diag_axis_mod', &
       'input argument has incorrect size', FATAL)
  tile_count = 1
  flag = 0
  do i = 1, size(ids(:))
     id = ids(i)
     if ( Axes(id)%cart_name == 'X' .or.  &
          Axes(id)%cart_name == 'Y' ) flag = flag + 1
!     --- both x/y axes found ---
     if ( flag == 2 ) then
        tile_count = Axes(id)%tile_count
        exit
     endif
  enddo

end function get_tile_count
!#######################################################################

function get_domain1d (id) result (Domain1)
  integer, intent(in) :: id
  type(domain1d)      :: Domain1   
  if (Axes(id)%Domain .NE. NULL_DOMAIN1D) then
     Domain1 = Axes(id)%Domain
  else
     Domain1 = NULL_DOMAIN1D
  endif
end function get_domain1d
!#######################################################################

function get_domain2d (ids) result (Domain2)
  integer, intent(in) :: ids(:)
  type(domain2d)      :: Domain2
  integer             :: i, id, flag

  if ( size(ids(:)) < 1 .or. size(ids(:)) > 4 ) call error_mesg  &
       ('get_domain2d in diag_axis_mod', &
       'input argument has incorrect size', FATAL)
  Domain2 = null_domain2d
  flag = 0
  do i = 1, size(ids(:))
     id = ids(i)
     if ( Axes(id)%cart_name == 'X' .or.  &
          Axes(id)%cart_name == 'Y' ) flag = flag + 1
!     --- both x/y axes found ---
     if ( flag == 2 ) then
        if (Axes(id)%Domain2 .NE. NULL_DOMAIN2D) then
           Domain2 = Axes(id)%Domain2
        endif
        exit
     endif
  enddo
end function get_domain2d

!#######################################################################
subroutine get_axes_shift ( ids, ishift, jshift ) 
  integer, intent(in)  :: ids(:)
  integer, intent(out) :: ishift, jshift
  integer              :: i, id

  !-- get the value of the shift.
  ishift = 0; jshift = 0
  do i = 1, size(ids(:))
     id = ids(i)
     select case (Axes(id)%cart_name)
     case ( 'X' )
        ishift = Axes(id)%shift
     case ( 'Y' )
        jshift = Axes(id)%shift
     end select
  end do

end subroutine get_axes_shift

!#######################################################################
function get_axis_num (axis_name, set_name) result (num)

! Returns index into axis table corresponding to a given axis name
  character(len=*), intent(in)           :: axis_name
  character(len=*), intent(in), optional ::  set_name
  integer                                :: num, set, n

  if (present(set_name)) then
     set = get_axis_set_num (trim(set_name))
  else
     set = 0
  endif
  num = 0
  do n = 1, num_def_axes
     if ( trim(axis_name) == trim(Axes(n)%name) .and. Axes(n)%set == set ) then
        num = n
        return
     endif
  enddo 
end function get_axis_num

!#######################################################################
function get_axis_set_num (set_name) result (num)

! Returns index in axis set table corresponding to a given axis set name

  character(len=*), intent(in) :: set_name
  integer                      :: num, iset

  num = 0
  do iset = 1, num_axis_sets
     if (set_name == Axis_sets(iset))then
        num = iset
        return
     endif
  enddo
end function get_axis_set_num

!#######################################################################
end module diag_axis_mod
