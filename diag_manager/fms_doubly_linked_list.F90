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

!> @defgroup fms_doubly_linked_list_mod fms_doubly_linked_list_mod
!> @ingroup diag_manager
!> @brief fms_doubly_linked_list_mod defines a generic doubly linked
!! list class and an iterator class for traversing the list.
!!
!> @author Miguel Zuniga Miguel.Zuniga@noaa.gov
!!
!! <TT>fms_doubly_linked_list_mod</TT>  defines a generic doubly linked
!! list class and an iterator class for traversing the list. It is
!! generic in the sense that the elements or objects it contains are
!! "class(*)" objects.  If additional typecheking or psossibly a
!! slightly different user interface is desired, consider creating
!! a wrapper or another class with this one for a memeber element and
!! procedures that are trivially implemeted by using this class.
!!
!! This versio is roughly a fortran translation of the C++ doubly linked list
!! class in the book ``Data Structures And Algorithm Analysis in C++", 3rd Edition,
!! by Mark Allen Weiss.

!> @file
!> @brief File for @ref fms_doubly_linked_list_mod
MODULE fms_doubly_linked_list_mod
  USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE
  implicit none
  !!TODO: COnsider setting the access (public,private) to functions, etc.
  !> Linked doubly-linked list node type.
  type, private :: node_type
    private
    class(*), pointer :: data => null()        !> The data pointed to by the node.
    type(node_type), pointer :: next => null() !> Pointer to the previous node.
    type(node_type), pointer :: prev => null() !> Pointer to the next node.
  end type node_type

  !> Linked list iterator
  type, public :: literator_type
    type(node_type), pointer :: current  !> Pointer to the current node.
    type(node_type), pointer :: end      !> A sentinel (non-data) node.
  contains
    procedure :: has_data => literator_has_data !< function returns true is there is data in the iterator.
    procedure :: next => literator_next !< function moves the iterator to the next data element.
    procedure :: get => literator_data !< function return a pointer to the current data.
  end type literator_type


  type, public :: linked_list_type
    !! Note we are overriding the default constructor with an
    ! interface of the same name
    type(node_type), pointer :: head !> The sentinal (non-data) head node of the linked list. .
    type(node_type), pointer :: tail !> The sentinel (non-data) tail node of the linked list.
    integer :: the_size               !> The number of data elements in the linked list.
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
  end type linked_list_type

  interface node_type
    module procedure :: node_constructor
  end interface node_type

  interface linked_list_type
    module procedure :: linked_list_constructor
  end interface linked_list_type

  interface literator_type
    module procedure :: literator_constructor
  end interface literator_type

contains

  !> Insert data d in a new node to be placed in front of the
  !! target node t_nd. Return an iterator that starts with the
  !! newly inserted node.
  function  insert_data( this, t_nd,  d ) result(liter)
    class(linked_list_type), intent(in out) :: this
    class(node_type), pointer, intent (in)  :: t_nd !> The target node.
    class(*), target, intent(in)            :: d    !> The data to insert.
    class(literator_type), allocatable      :: liter !> A linked list iterator
    !!
    class(node_type), pointer :: nd                  !> The new node that is to "hold" the data.
    allocate(nd)
    nd%data => d
    !!  Insert nd into list so that list section [prev node <--> target node ] looks like
    !!  [prev node <--> new nd <--> target node]. The four pointers pointing to and/or
    !! from "new nd" need to be set. Therefore :
    !!  a) The new nd's prev needs to be whatever was the targets prev:
    nd%prev => t_nd%prev
    !!  b) New node nd's next is obviously the target node:
    nd%next => t_nd
    !!  c) the next of the prev node needs to point to the new node nd:
    t_nd%prev%next => nd
    !!  d) target node's prev needs to point to the new node :
    t_nd%prev => nd
    this%the_size = this%the_size + 1
    liter = literator_type(nd, this%tail)
  end function insert_data

  !> Remove Node nd from the linked tree. Return the iterator
  !! that begins with the nex node after nd, and ends with the
  !! list end node. Return the list iterator if the node cannot
  !! be removed.
  function  remove_node( this, nd ) result( litr)
    class(linked_list_type), intent(in out) :: this
    type(node_type), pointer, intent(in out) :: nd  !> The node to remove from the list.
    class(literator_type), ALLOCATABLE  :: litr     !> The iterator starting from whthe node that was
                                                    !>   following the removed node.
    !Dont even try to remove the head and tail nodes!
    if (.not. ( associated (this%head , nd ) .or. &
      associated (this%tail , nd ) )) THEN
      litr = literator_type( nd%next, this%tail )
      nd%prev%next => nd%next
      nd%next%prev => nd%prev
      deallocate(nd)
      this%the_size =  this%the_size - 1
    else
      litr = this%get_literator()
    endif
  end function remove_node


  !>  Remove the head (first data node) of the list. Return an
  !! iterator to the remaining list.
  function pop_at_front (this ) result( liter )
    class(linked_list_type), intent(in out) :: this
    class(literator_type), allocatable :: liter  !> The iterator for the remaining list.
    !!
    class(node_type), pointer :: nd
    if(this%the_size /= 0) then
      nd => this%head%next
      liter =  this%remove(nd)
    else
      liter = this%get_literator()
    endif
  end function pop_at_front

  !>  Remove the tail (last data node) of the list. Return an
  !! iterator to the remaining list.
  function pop_at_back (this ) result( liter )
    class(linked_list_type), intent(in out) :: this
    class(literator_type) , allocatable :: liter  !> The iterator for the remaining list.
    !!
    class(node_type), pointer :: nd
    if(this%the_size /= 0) then
      nd => this%tail%prev
      liter = this%remove( nd )
    else
      liter = this%get_literator()
    endif
  end function pop_at_back



  !> push (insert) data at the head of the list
  function push_at_front( this, d ) result(litr)
    class(linked_list_type), intent(in out) :: this
    class(*), target, intent(in out) :: d      !> The data to insert.
    class(literator_type), allocatable :: litr !> The iterator for the resultant list.
    litr = this%insert (this%head%next, d)
  end function push_at_front

  !> push (insert) data at the end of the list
  function push_at_back( this, d ) result(litr)
    class(linked_list_type), intent(in out) :: this
    class(*), target, intent(in out) :: d              !> The data to insert.
    class(literator_type), allocatable :: litr         !> The iterator for the resultant list.
    litr =  this%insert (this%tail, d)
  end function push_at_back

  !> Constructor implementation for the node_type
  function node_constructor () result (nd)
    type(node_type), allocatable :: nd  !> The allocated node.
    allocate(nd)
    nd%data => null()
    nd%prev => null()
    nd%next => null()
  end function node_constructor

  !> The linked list construcotr.
  function linked_list_constructor () result (ll)
    type(linked_list_type), allocatable :: ll !> The resultant linked list to be reutrned.
    allocate(ll)
    allocate(ll%head)
    allocate(ll%tail)
    !!print *, 'associated(ll%head) :' , associated(ll%head), &
    !! ' associated(ll%head) :' , associated(ll%head)
    ll%head%next => ll%tail
    ll%tail%prev => ll%head
    ll%the_size = 0
  end function linked_list_constructor


  !> The  list iterator constructor.
  function literator_constructor ( fnd, tnd ) result (litr)
    type (node_type), pointer  :: fnd    !> What will be the first (and current) data node of the iterator.
    type (node_type), pointer  :: tnd    !> What will be the last (and a non-data) node for the iterator.
    !! node pointed to be the iterator.
    type (literator_type), allocatable :: litr  !> The resultant linked list to be reutrned.
    allocate(litr)
    litr%current => fnd
    litr%end => tnd
  end function literator_constructor

  !> Return the size (the number of data elements held) of the lined list.
  function get_size (this) result (sz)
    class(linked_list_type), intent(in out) :: this
    integer  :: sz           !> The size (number of data elements)
    sz = this%the_size
  end function get_size

!> Returns true if there are zero (0) elements in the container; false otherwise.
  function is_size_zero (this) result (r)
    class(linked_list_type), intent(in out) :: this
    logical :: r !> True iff the size (number of data elements) is zero.
    if (this%the_size == 0) then
      r = .true.
    else
      r = .false.
    end if
  end function is_size_zero

  !> Return a forward iterator for the linked list.
  function get_forward_literator(this) result (litr)
    class(linked_list_type), intent(in) :: this
    class(literator_type), ALLOCATABLE :: litr !> The iterator to be returned
    litr = literator_type( this%head%next, this%tail )
  end function get_forward_literator


  !> Returns true iff the iterator has data.
  function literator_has_data( this ) result( r )
    class(literator_type), intent(in) :: this
    logical :: r !> The result true/false.
    if( associated (this%current, this%end )) then
      r = .false.
    else
      r = .true.
    end if
  end function literator_has_data

  !> Move the iterators current data node pointer to the next data node iff there is a next available.
  function literator_next( this ) result( status )
    class(literator_type), intent(in out ) :: this
    integer :: status !> Zero iff success. Failure possible if iterator does not have data.
    status = -1
    if(this%has_data() .eqv. .true.) then
      this%current => this%current%next
      status = 0
    endif
  end function literator_next



  !> Get the current data object pointed to by the iterator.
  !! function does not allocate or assign the result if
  !! the user mistakenly called it without data present.
  function  literator_data( this ) result( rd )
    class(literator_type), intent(in) :: this
    class(*),  pointer  :: rd !> The current data element of the iterator.
    if (this%has_data() .eqv. .true.) then
      rd => this%current%data
    endif
  end function literator_data

  !> Iterate over all the nodes, remove them and deallocate the client data
  !! that the node was holding.
  subroutine clear_all( this  )
    class(linked_list_type), intent(inout) :: this
    class(node_type), pointer :: nd                      !> A pointer to linked list node
    class(literator_type), allocatable :: iter           !> A linked list iterator
    class(*),  pointer  :: pdata                         !> A pointer to the data.
    !
    do while( this% the_size /= 0)
      nd => this%head%next
      iter =  this%remove(nd)
      pdata => iter%get()
      if (associated(pdata) .eqv. .false.) then
        call  error_mesg ('doubly_linked_list:clear_all', &
          'linked list destructor containes unassociated data pointer', &
          WARNING)
      else
        deallocate(pdata)
      endif
    end do
  end subroutine clear_all

  !> A destructor that deallocates every node and each nodes data element.
    subroutine destructor(this)
    type(linked_list_type ) :: this  !Note for destructors its needs to be type and not class!
    call this%clear()
    deallocate(this%head)
    deallocate(this%tail)
  end subroutine destructor


end module fms_doubly_linked_list_mod
!> @}
! close documentation grouping
