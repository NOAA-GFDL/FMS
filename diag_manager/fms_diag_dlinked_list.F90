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

!> @defgroup fms_diag_dlinked_list_mod fms_diag_dlinked_list_mod
!> @ingroup diag_manager
!> @brief fms_diag_dlinked_list_mod defines a generic doubly linked
!! list class and an iterator class for traversing the list.
!!
!> @author Miguel Zuniga
!!
!! <TT>fms_diag_dlinked_list_mod</TT>  defines a generic doubly linked
!! list class and an iterator class for traversing the list. It is
!! generic in the sense that the elements or objects it contains are
!! "class(*)" objects.  If additional typecheking or psossibly a
!! slightly different user interface is desired, consider creating
!! a wrapper or another class with this one for a memeber element and
!! procedures that are trivially implemeted by using this class.
!!
!! This version is roughly a fortran translation of the C++ doubly linked list
!! class in the book ``Data Structures And Algorithm Analysis in C++", 3rd Edition,
!! by Mark Allen Weiss.

!> @file
!> @brief File for @ref fms_diag_dlinked_list_mod
!> @addtogroup fms_diag_dlinked_list_mod
!> @{
MODULE fms_diag_dlinked_list_mod
   USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE
   implicit none
   !> The doubly-linked list node type.
   type, public:: FmsDlListNode_t
      private
      class(*), pointer :: data_ptr => null()        !< The data pointed to by the node.
      type(FmsDlListNode_t), pointer :: next => null() !< A pointer to the previous node.
      type(FmsDlListNode_t), pointer :: prev => null() !< A pointer to the next node.
   end type FmsDlListNode_t

   !> Linked list iterator
   type, public :: FmsDllIterator_t
      private
      type(FmsDlListNode_t), pointer :: current=>null()  !< A pointer to the current node.
      type(FmsDlListNode_t), pointer :: end =>null()     !< A sentinel (non-data) node.
   contains
      procedure :: has_data => literator_has_data !< Function returns true if there is data in the iterator.
      procedure :: next => literator_next !< Function moves the iterator to the next data element. Used in
                                          !< conjunction with function has_data().
      procedure :: get => literator_data !< Function return a pointer to the current data. Used in conjunction
                                          !< with function has_data().
      procedure :: get_current_node_pointer => get_current_node_ptr !< Return  the current node pointer.
   end type FmsDllIterator_t

   !> The doubly-linked list type. Besides the member functions, see the
   !! associated iterator class ( FmsDllIterator_t) for traversal, and note that
   !! the default constructor is overriden with an interface of the same name.
   type, public :: FmsDlList_t
      private
      type(FmsDlListNode_t), pointer :: head=>null() !< The sentinal (non-data) head node of the linked list. .
      type(FmsDlListNode_t), pointer :: tail=>null() !< The sentinel (non-data) tail node of the linked list.
      integer :: the_size        !< The number of data elements in the linked list.
   contains
      procedure :: push_back => push_at_back
      procedure :: pop_back => pop_at_back
      procedure :: remove => remove_node
      procedure :: get_literator => get_forward_literator
      procedure :: size => get_size
      procedure :: is_empty => is_size_zero
      procedure :: clear => clear_all
      procedure :: initialize => linked_list_initializer
      final :: destructor
      procedure  :: insert => insert_data

   end type FmsDlList_t

   interface FmsDlList_t
      module procedure :: linked_list_constructor
   end interface FmsDlList_t

   interface FmsDllIterator_t
      module procedure :: literator_constructor
   end interface FmsDllIterator_t

contains

   !> @brief Insert data d in a new node to be placed in front of the
   !! target node t_nd.
   !! @return Returns an iterator that starts with the newly inserted node.
   function  insert_data( this, t_nd,  d ) result(liter)
      class(FmsDlList_t), intent(in out) :: this !<The instance of the class that this function is bound to.
      type(FmsDlListNode_t), pointer, intent (in)  :: t_nd  !< The target node.
      class(*), target, intent(in)                 :: d     !< The data to insert.
      class(FmsDllIterator_t), allocatable         :: liter !< A linked list iterator.
      type(FmsDlListNode_t), pointer :: nd              !< The new node that is to "hold" the data.
      allocate(nd)
      nd%data_ptr => d
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
      liter = FmsDllIterator_t(nd, this%tail)
   end function insert_data

   !> @brief Remove Node nd from the linked tree.
   !! @return Return the iterator that begins with the next node after nd, and ends with
   !! the list end node. Returns the list iterator if the node cannot be removed.
   function  remove_node( this, nd ) result( litr)
      class(FmsDlList_t), intent(in out) :: this !<The instance of the class that this function is bound to.
      type(FmsDlListNode_t), pointer, intent(in out) :: nd  !< The node to remove from the list.
      class(FmsDllIterator_t), ALLOCATABLE  :: litr         !< The iterator starting from the node that was
      !<  following the removed node.
      !Dont even try to remove the head and tail nodes!
      if (.not. ( associated (this%head , nd ) .or. &
         associated (this%tail , nd ) )) THEN
         litr = FmsDllIterator_t( nd%next, this%tail )
         nd%prev%next => nd%next
         nd%next%prev => nd%prev
         deallocate(nd)
         this%the_size =  this%the_size - 1
      else
         litr = this%get_literator()
      endif
   end function remove_node


   !> @brief Remove the tail (last data node) of the list.
   !! @return  Returns an iterator to the remaining list.
   function pop_at_back (this ) result( liter )
      class(FmsDlList_t), intent(in out) :: this !<The instance of the class that this function is bound to.
      class(FmsDllIterator_t) , allocatable :: liter  !< The iterator for the remaining list.
      !!
      type(FmsDlListNode_t), pointer :: nd !< The node being removed.
      if(this%the_size /= 0) then
         nd => this%tail%prev
         liter = this%remove( nd )
      else
         liter = this%get_literator()
      endif
   end function pop_at_back

   !> @brief Push (insert) data at the end of the list
   !> @return Returns an iterator that starts at the tail of the list.
   function push_at_back( this, d ) result(litr)
      class(FmsDlList_t), intent(in out) :: this   !<The instance of the class that this function is bound to.
      class(*), target, intent(in out) :: d        !< The data to insert.
      class(FmsDllIterator_t), allocatable :: litr  !< The iterator for the resultant list.
      litr =  this%insert (this%tail, d)
   end function push_at_back

   !> @brief Constructor for the linked list.
   !! @return Returns a newly allocated linked list instance.
   !! TODO: This function is not used since (observed on Intel compilers) with
   !! a finalize keyword on the destructor, when this function returns and ll
   !! goes out of scope, th allocations in initialized are undome
   !! whether ot not ll is declared a pointer or allocatable
   function linked_list_constructor () result (ll)
      type(FmsDlList_t), pointer :: ll !< The resultant linked list to be reutrned.
      allocate(ll)
      call ll%initialize()
   end function linked_list_constructor

   !> @brief Initializer for the linked list.
   !! @return Returns a newly allocated linked list instance.
   subroutine linked_list_initializer( this )
      class(FmsDlList_t), intent(inout) :: this  !<The instance of the class that this function is bound to
      if( associated(this%head) .or. associated(this%tail)) then
         call error_mesg('fms_diag_dlinked_list','linked list is already initalized', WARNING)
      else
         allocate(this%head)
         allocate(this%tail)
         this%head%next => this%tail
         this%tail%prev => this%head
         this%the_size = 0
      endif
   end subroutine linked_list_initializer


   !> @brief The list iterator constructor.
   !! @return Returns a newly allocated list iterator.
   function literator_constructor ( fnd, tnd ) result (litr)
      type (FmsDlListNode_t), pointer  :: fnd
      !< The sentinal (non-data) "first node" of the iterator will be fnd
      type (FmsDlListNode_t), pointer  :: tnd
      !< The sentinal (non-data) "last node" of the iterator will be tnd.
      type (FmsDllIterator_t), allocatable :: litr  !< The resultant linked list to be reutrned.
      allocate(litr)
      litr%current => fnd
      litr%end => tnd
   end function literator_constructor

   !> @brief Getter for the size (the number of data elements) of the linked list.
   !! @return Returns the size of the lined list.
   function get_size (this) result (sz)
      class(FmsDlList_t), intent(in out) :: this
      !<The instance of the class that this function is bound to.
      integer  :: sz           !< The size (number of data elements)
      sz = this%the_size
   end function get_size

!> @brief Determines if the size (number of data elements) of the list is zero.
!! @return Returns true if there are zero (0) data elements in the list; false otherwise.
   function is_size_zero (this) result (r)
      class(FmsDlList_t), intent(in out) :: this
      !<The instance of the class that this function is bound to.
      logical :: r !< True iff the size (number of data elements) is zero.
      if (this%the_size == 0) then
         r = .true.
      else
         r = .false.
      end if
   end function is_size_zero

   !> @brief Create and return a new forward iterator for the list.
   !> @return Returns a forward iterator for the linked list.
   function get_forward_literator(this) result (litr)
      class(FmsDlList_t), intent(in) :: this !<The instance of the class that this function is bound to.
      class(FmsDllIterator_t), ALLOCATABLE :: litr !< The iterator to be returned
      litr = FmsDllIterator_t( this%head%next, this%tail )
   end function get_forward_literator

   !> @brief Determine if the iterator has data.
   !> @return Returns true iff the iterator has data.
   function literator_has_data( this ) result( r )
      class(FmsDllIterator_t), intent(in) :: this
      !<The instance of the class that this function is bound to.
      logical :: r !< The result true/false.
      if( associated (this%current, this%end )) then
         r = .false.
      else
         r = .true.
      end if
   end function literator_has_data

   !> @brief Move the iterators current data node pointer to the next data node.
   !! @return Returns a status of 0 if succesful, -1 otherwise.
   function literator_next( this ) result( status )
      class(FmsDllIterator_t), intent(in out ) :: this
      integer :: status !< The returned status. Failure possible is if iterator does not have data.
      status = -1
      if(this%has_data() .eqv. .true.) then
         this%current => this%current%next
         status = 0
      endif
   end function literator_next

   !> @brief Get the current data object pointed to by the iterator.
   !! function does not allocate or assign the result if
   !! the user mistakenly called it without data present.
   !! @return Returns a pointer to the current data.
   function  literator_data( this ) result( rd )
      class(FmsDllIterator_t), intent(in) :: this !<The instance of the class that this function is bound to.
      class(*),  pointer  :: rd !< The current data element of the iterator.
      rd => null()
      if (this%has_data() .eqv. .true.) then
         rd => this%current%data_ptr
      endif
   end function literator_data

  !> @brief Get the current data object pointed to by the iterator.
   !! function does not allocate or assign the result if
   !! the user mistakenly called it without data present.
   !! @return Returns a pointer to the current data.
   function  get_current_node_ptr( this ) result( pn )
      class(FmsDllIterator_t), intent(in) :: this !<The instance of the class that this function is bound to.
      type(FmsDlListNode_t), pointer :: pn        !< The sentinel (non-data) tail node of the linked list.
      pn  => this%current
   end function get_current_node_ptr

   !> @brief Iterate over all the nodes and remove them. Also (by overridable default), it deallocates the
   !! client data associated with the nodes.
   subroutine clear_all( this, data_dealloc_flag)
      class(FmsDlList_t), intent(inout) :: this !<The instance of the class that this function is bound to.
      logical, optional, intent(in) :: data_dealloc_flag    !< If not present or .true., client data is deallocated.
      type(FmsDlListNode_t), pointer :: nd           !< A pointer to linked list node.
      class(FmsDllIterator_t), allocatable :: iter   !< A linked list iterator.
      class(*),  pointer  :: pdata =>null()          !< A pointer to the data.
      logical :: data_dealloc_f    !< Set to data_dealloc_flag if present, otherwise its .true.
      !
      data_dealloc_f = .true.
      if( PRESENT(data_dealloc_flag) ) then
         data_dealloc_f = data_dealloc_flag
      endif
      do while( this% the_size /= 0)
         nd => this%head%next
         pdata => nd%data_ptr
         iter =  this%remove(nd)
         if(data_dealloc_f .eqv. .true.) then
            if (associated(pdata) .eqv. .false.) then
               call  error_mesg ('fms_diag_dlinked_list', &
                  'In clear_all; linked node contains node with unassociated data pointer', &
                  WARNING)
            else
               deallocate(pdata)
            endif
         endif
      end do
   end subroutine clear_all

   !>  @brief A destructor that deallocates every node and each nodes data element. !Note
   !! that for the data elements to not be de-allocated, function clear() (or clear_all() )
   !! with the appropriate arguments must be called.
   subroutine destructor(this)
      type(FmsDlList_t) :: this  !<The instance of the type that this function is bound to.
      !! Note in the line above we use "type' and not "class" - needed for destructor definitions.
      !! TODO: This NOTE message below may be inaproppriate for this class.
      call error_mesg('fms_diag_dlinked_list', 'Entered destructor.',NOTE)
      call this%clear()
      deallocate(this%head)
      this%head=>null()
      deallocate(this%tail)
      this%tail=>null()
   end subroutine destructor

end module fms_diag_dlinked_list_mod
!> @}
! close documentation grouping
