!> \author Miguel Zuniga
!> \email miguel.zuniga@noaa.gov
!! \brief A module for a double linked-list.
!!
!! \description A double linked list module that is rougly a fortran translation of the C++
!! boubly linked list from the book  ``Data Structures And Algorithm Analysis in C++", 3rd Edition,
!!  by Mark Allen Weiss
module fms_doubly_linked_list_mod
  implicit none

  !> Linked list node
  type, private :: node
    private
    class(*), pointer :: data => null()
    type(node), pointer :: next => null()
    type(node), pointer :: prev => null()
  end type node

  !> Linked list literator
  type, public :: literator
    type(node), pointer :: current
    type(node), pointer :: end
  contains
    procedure :: has_data => literator_has_data
    procedure :: next => literator_next
    procedure :: get => literator_data
  end type literator

  type, public :: linked_list
    !! Note we are overriding the default constructor with an
    ! interface of the same name
    type(node), pointer :: head
    type(node), pointer :: tail
    integer :: theSize
  contains
    procedure :: push_front => push_at_front
    procedure :: push_back => push_at_back
    procedure :: pop_front => pop_at_front
    procedure :: pop_back => pop_at_back
    procedure :: insert => insert_data
    procedure :: remove => remove_node
    procedure :: get_literator => get_forward_literator
    procedure :: size => get_size
    procedure :: is_empty => is_size_zero
    procedure :: clear => clear_all
    final :: destructor
  end type linked_list

  interface node
    module procedure :: node_constructor
  end interface node

  interface linked_list
    module procedure :: linked_list_constructor
  end interface linked_list

  interface literator
    module procedure :: literator_constructor
  end interface literator

contains

  !> Insert data d in a new node to be placed in front of the
  ! target node t_nd. Return an iterator that starts with the
  ! newly inserted node.
  function  insert_data( this, t_nd,  d ) result(liter)
    class(linked_list), intent(in out) :: this
    class(node), pointer, intent (in)  :: t_nd
    class(*), target, intent(in)       :: d
    class(literator), allocatable      :: liter
    !!
    class(node), pointer :: nd
    allocate(nd)
    nd%data => d
    nd%prev => t_nd%prev
    nd%next => t_nd
    t_nd%prev%next => nd
    t_nd%prev => nd
    this%theSize = this%theSize + 1
    liter = literator(nd, this%tail)
  end function insert_data

  !> Remove Node nd from the linked tree. Return the iterator
  ! that begins with the nex node after nd, and ends with the
  ! list end node. Return the list iterator if the node cannot
  ! be removed.
  function  remove_node( this, nd ) result( litr)
    class(linked_list), intent(in out) :: this
    type(node), pointer, intent(in out) :: nd
    class(literator), ALLOCATABLE  :: litr
    !Dont even try to remove the head and tail nodes!
    if (.not. ( associated (this%head , nd ) .or. &
      associated (this%tail , nd ) )) THEN
      litr = literator( nd%next, this%tail )
      nd%prev%next => nd%next
      nd%next%prev => nd%prev
      deallocate(nd)
      this%theSize =  this%theSize - 1
    else
      litr = this%get_literator()
    endif
  end function remove_node


  !>  Remove the data at the head of the list. Return an
  !  Iterator starting with the next (new front) node, and
  ! Return the list iterator if the list is already empty
  function pop_at_front (this ) result( liter )
    class(linked_list), intent(in out) :: this
    class(literator), allocatable :: liter
    !!
    class(node), pointer :: nd
    if(this%theSize /= 0) then
      nd => this%head%next !!  TODO: verify
      liter =  this%remove(nd)
    else
      liter = this%get_literator()
    endif
  end function pop_at_front

  function pop_at_back (this ) result( liter )
    class(linked_list), intent(in out) :: this
    class(literator) , allocatable :: liter
    !!
    class(node), pointer :: nd
    if(this%theSize /= 0) then
      nd => this%tail%prev !! TODO: verify
      liter = this%remove( nd )
    else
      liter = this%get_literator()
    endif
  end function pop_at_back



  !> push (insert) data at the head of the list
  function push_at_front( this, d ) result(litr)
    class(linked_list), intent(in out) :: this
    class(*), target, intent(in out) :: d
    class(literator), allocatable :: litr
    litr = this%insert (this%head%next, d)
  end function push_at_front

  !> push (insert) data at the end of the list
  function push_at_back( this, d ) result(litr)
    class(linked_list), intent(in out) :: this
    class(*), target, intent(in out) :: d
    class(literator), allocatable :: litr
    litr =  this%insert (this%tail, d)
  end function push_at_back

  function node_constructor () result (nd)
    type(node), allocatable :: nd
    allocate(nd)
    nd%data => null()
    nd%prev => null()
    nd%next => null()
    print * , "In node constructor"
  end function node_constructor

  function linked_list_constructor () result (ll)
    type(linked_list), allocatable :: ll
    allocate(ll)
    allocate(ll%head)
    allocate(ll%tail)
    print *, 'associated(ll%head) :' , associated(ll%head), &
      ' associated(ll%head) :' , associated(ll%head)
    ll%head%next => ll%tail
    ll%tail%prev => ll%head
    print * , "In list constructor"
  end function linked_list_constructor


  function literator_constructor ( fnd, tnd ) result (litr)
    type (node), pointer  :: fnd
    type (node), pointer  :: tnd
    type (literator), allocatable :: litr
    !type (literator), pointer :: litr
    allocate(litr)
    litr%current => fnd
    litr%end => tnd
  end function literator_constructor

  function get_size (this) result (sz)
    class(linked_list), intent(in out) :: this
    integer  :: sz
    sz = this%theSize
  end function get_size


  function is_size_zero (this) result (r)
    class(linked_list), intent(in out) :: this
    logical :: r
    if (this%theSize == 0) then
      r = .true.
    else
      r = .false.
    end if
  end function is_size_zero

  function get_forward_literator(this) result (litr)
    class(linked_list), intent(in) :: this
    class(literator), ALLOCATABLE :: litr
    litr = literator( this%head%next, this%tail )
  end function get_forward_literator


  !! Should return true unless the currrent is pointing to the next
  function literator_has_data( this ) result( r )
    class(literator), intent(in) :: this
    logical :: r
    if( associated (this%current, this%end )) then
      r = .false.
    else
      r = .true.
    end if
  end function literator_has_data

  !> Move to the next iff there is a next
  function literator_next( this ) result( status )
    class(literator), intent(in out ) :: this
    integer :: status
    if(this%has_data() .eqv. .true.) then
      this%current => this%current%next
    endif
    status = 0
  end function literator_next



  !> Get the current data object pointed to by the iterator.
  !! function does not allocate or assign the result if
  !! the user mistakenly called it without data present.
  function  literator_data( this ) result( rd )
    class(literator), intent(in) :: this
    class(*),  pointer  :: rd
    if (this%has_data() .eqv. .true.) then
      rd => this%current%data
    endif
  end function literator_data

  !! Iterate over all the nodes, remove them and deallocate the client data
  !! that the node was holding.
  subroutine clear_all( this  )
    class(linked_list), intent(inout) :: this
    class(node), pointer :: nd
    class(literator), allocatable :: iter
    class(*),  pointer  :: pdata
    !
    do while( this% theSize /= 0)
      nd => this%head%next
      iter =  this%remove(nd)
      pdata => iter%get()
      if (associated(pdata) .eqv. .false.) then
        print *, "clear all. pdata is not associated"
      else
       deallocate(pdata)
      endif
    end do
  end subroutine clear_all

  ! A destructor that deallocates every node and each nodes data element.
    subroutine destructor(this)
    type(linked_list ) :: this  !Note for destructors its needs to be type and not class!
    print *, "linked list . Starting destructor size=", this%theSize
    call this%clear()
    if (this%theSize /= 0) then
      print *, "linked list . destructor Error size=", this%theSize
    endif
    deallocate(this%head)
    deallocate(this%tail)
  end subroutine destructor


end module fms_doubly_linked_list_mod
