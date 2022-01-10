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

!> @defgroup fms_diag_object_container_mod fms_diag_object_container_mod
!> @ingroup diag_manager
!> @brief fms_diag_object_container_mod defines a container class and iterator class
!! for inserting, removing and searching for <TT>fmsDiagObject_type</TT> instances
!!
!> @author Miguel Zuniga
!!
!! <TT>fms_diag_object_container_mod</TT> defines a container for inserting, removing and
!! searching for <TT>fmsDiagObject_type</TT> instances. It also defined an iterator for
!! the data in the container. The value returned by the fms_diag_object function get_id()
!! is used for search key comparison.
!!
!! Most of the functions in class FmsDiagObjectContainer_t are simple wrappers over
!! those of the underlying  <TT>fms_doubly_linked_list_mod</TT> class. The find/search
!! are a little more than that, and what FmsDiagObjectContainer_t provides over the
!! underlying liked  list is the search function, type checking, convenience, and a
!! fixed user interface defined for the intended use.
!!
!> @file
!> @brief File for @ref fms_diag_object_container_mod
!> @addtogroup fms_diag_object_container_mod
!> @{
MODULE fms_diag_object_container_mod
   use fms_diag_object_mod, only: fmsDiagObject_type
   USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE

   !! Since this version is based on the FDS linked list:
   use fms_diag_dlinked_list_mod, only : FmsDlList_t, FmsDllIterator_t, FmsDlListNode_t

   implicit none

   !> @brief A container of fmsDiagObject_type instances providing insert, remove ,
   !!  find/search, and size public member functions. Iterator is provided by
   !!  the associated iterator class (see  dig_obj_iterator class).
   !!
   !!  This version does not enforce uniqueness of ID keys (I.e. it is not a set).
   !!
   type, public:: FmsDiagObjectContainer_t
   private
      TYPE (FmsDlList_t), ALLOCATABLE :: the_linked_list !< This version based on the FDS linked_list.
   contains
      procedure :: insert => insert_diag_object
      procedure :: remove => remove_diag_object
      procedure :: find   => find_diag_object
      procedure :: size   => get_num_objects
      procedure :: iterator  => get_iterator
      final :: destructor
   end type FmsDiagObjectContainer_t


   !> @brief Iterator used to traverse the objects of the container.
   type, public :: FmsDiagObjIterator_t
   private
      type(FmsDllIterator_t) :: liter !< This version based on the FDS linked_list (and its iterator).
   contains
      procedure :: has_data => literator_has_data
      procedure :: next => literator_next
      procedure :: get => literator_data
   end type FmsDiagObjIterator_t

   interface FmsDiagObjectContainer_t
      module procedure :: diag_object_container_constructor
   end interface FmsDiagObjectContainer_t

   interface FmsDiagObjIterator_t
      module procedure :: diag_obj_iterator_constructor
   end interface FmsDiagObjIterator_t


contains

   !> @brief Returns an empty iterator if a diag object with this ID was not found.
   !! If the diag object was found, return an iterator with the current object being
   !! the found object, ad the last/anchor being the last/anchor of the container.
   !! Note that this routine can accept an optional iterator as input, which
   !! is useful for chaining searches, which may be needed if there are key duplicates.
   !! @return In iterator that starts from the inserted object.
   function find_diag_object (this, id , iiter) result (riter)
     class (FmsDiagObjectContainer_t), intent (in out) :: this
     !<The instance of the class that this function is bound to.
      integer , intent (in) :: id   !< The id of the object to find.
     class(FmsDiagObjIterator_t), intent (in), optional :: iiter
     !< An (optional) iterator over the searchable set.
      class(FmsDiagObjIterator_t) , allocatable :: riter !< The resultant iterator to the object.
      class(fmsDiagObject_type),  pointer:: ptdo  !< A pointer to temporaty diagnostic object
      integer :: status                        !< A status from iterator operations.
      !!
      if(present (iiter)) then
         riter = iiter
      else
         riter = this%iterator() !!An iterator over the entire container
      endif
      do while( riter%has_data() .eqv. .true.)
         ptdo => riter%get()
         if(id == ptdo%get_id() ) then
            EXIT
         end if
         status = riter%next()
      end do
   end function find_diag_object

   !> @brief insert diagnostic object obj with given id.
   !! Objects are inserted at the back / end of the list
   !! This version of the container also enforces that the
   !! objects ID is equal the input id.
   !! @return A status of -1 if there was an error, and 0 otherwise.
   function insert_diag_object (this, id, obj) result (status)
      class (FmsDiagObjectContainer_t), intent (in out) :: this
      integer,  intent (in) :: id                     !< The id of the object to insert.
      class(fmsDiagObject_type) , intent (in out) :: obj !< The object to insert
      integer :: status                               !< The returned status. 0 for success.
      class(FmsDllIterator_t), allocatable ::  tliter   !< A temporary iterator.

      status = -1
      if ( id .ne. obj%get_id() ) then
         !!TODO: log error
      endif
      tliter =  this%the_linked_list%push_back( obj )
      if(tliter%has_data() .eqv. .true. ) then
         status = 0
      endif
   end function

   !> @brief Remove and return the first object in the container with the corresponding id .
   !! Note that if the client code does not already have a reference to the object being
   !! removed, then  the client may want to to use procedure find before using procedure remove.
   !! If procedure find is used, consider calling remove with the iterator returned from find.
   !! @return In iterator starting from the node that was following the removed node.
   function remove_diag_object (this, id, iiter ) result (riter)
     class (FmsDiagObjectContainer_t), intent (in out) :: this
     !<The instance of the class that this function is bound to.
     integer , intent (in) :: id !< The id of the object to remove.
     class(FmsDiagObjIterator_t), intent (in), optional :: iiter
      !< An (optional) iterator over the searchable set.
      class(FmsDiagObjIterator_t), allocatable:: riter !< The resultant iterator
      class(FmsDllIterator_t), allocatable :: temp_liter !< A temporary iterator
      type(FmsDlListNode_t), pointer :: pn
      if(present (iiter)) then
         riter = iiter
      else
         riter = this%iterator() !!An iterator over the entire container
      endif
      !Find the object in the container.
      riter = this%find ( id , riter)
      pn => riter%liter%get_current_node_pointer()
      temp_liter = this%the_linked_list%remove( pn )
      riter = FmsDiagObjIterator_t(temp_liter)
   end function

   !> @brief  Getter for the number of objects help in the container.
   !! @return Return the number of objects..
   function get_num_objects (this ) result (sz)
     class (FmsDiagObjectContainer_t), intent (in out) :: this
     !< The instance of the class that this function is bound to.
      integer :: sz                     !< The returned result - the number of objects in container.
      sz = this%the_linked_list%size()
   end function


   !> @brief Return an iterator for the objects in the container.
   !! @return An iterator for the objects in the container.
   function get_iterator (this) result (oliter)
     class (FmsDiagObjectContainer_t), intent (in) :: this
     !<The instance of the class that this function is bound to.
      class(FmsDiagObjIterator_t) , allocatable :: oliter !< The returned iterator.
      class(FmsDllIterator_t), allocatable ::  temp_iter !< A temporary linked list iterator
      temp_iter = this%the_linked_list%get_literator()
      oliter = FmsDiagObjIterator_t( temp_iter )
   end function

   !> @brief A consructor for a container's iterator.
   !! @return  An for a container's iterator.
   function diag_obj_iterator_constructor( iliter ) result (diag_itr)
     class (FmsDllIterator_t), allocatable :: iliter
     !< An iterator. Normally the one that the container is based on.
      class (FmsDiagObjIterator_t), allocatable :: diag_itr !< The returned diag object iterator.
      allocate(diag_itr)
      diag_itr%liter = iliter;
   end function diag_obj_iterator_constructor

   !> @brief The default consructor for the container.
   !! @return Returns a container.
   function diag_object_container_constructor () result (doc)
      type(FmsDiagObjectContainer_t), allocatable :: doc !< The resultant container.
      allocate(doc)
      doc%the_linked_list = FmsDlList_t()
      !! print * , "In DOC constructor"
   end function diag_object_container_constructor

   !> @brief Determines if there is more data that can be accessed via the iterator.
   !> @return Returns true iff more data can be accessed via the iterator.
   function literator_has_data( this ) result( r )
     class(FmsDiagObjIterator_t), intent(in) :: this
     !<The instance of the class that this function is bound to.
      logical :: r                  !< True if this iterator has data.
      r = this%liter%has_data()
   end function literator_has_data

   !> @brief Move the iterator to the next object.
   !! @return Returns a status 0 if sucessful, or -1 if failed.
   function literator_next( this ) result( status )
     class(FmsDiagObjIterator_t), intent(in out ) :: this
     !<The instance of the class that this function is bound to.
      integer :: status  !< The returned status.
      status = -1
      status = this%liter%next()
   end function literator_next

   !> @brief Get the current data the iterator is pointing to.
   !! Note the common use case is to call function has_data to decide if
   !! this function should be called (again).
   !! @return Returns a pointer to the current data.
   function  literator_data( this ) result( rdo )
     class(FmsDiagObjIterator_t), intent(in) :: this
     !<The instance of the class that this function is bound to.
      class(fmsDiagObject_type),  pointer :: rdo !< The resultant diagnostic object.
      class(*),  pointer :: gp !< A eneric typed object in the container.

      rdo => null()
      gp => this%liter%get()
      select type(gp)
       type is (fmsDiagObject_type)  !! "type is", not the (polymorphic) "class is"
         rdo => gp
       class default
         CALL  error_mesg ('diag_object_container:literator_data', &
            'Data to be accessed via iterator is not of expected type.',FATAL)
      end select
   end function literator_data

   !> @brief The destructor for the container.
   subroutine destructor(this)
     type(FmsDiagObjectContainer_t) :: this
     !<The instance of the class that this function is bound to.
     !Note in the line above we use "type' and not "class" - needed for destructor definitions.
      deallocate(this%the_linked_list)
   end subroutine destructor


end module fms_diag_object_container_mod
!> @}
! close documentation grouping

