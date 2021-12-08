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
!> @brief fms_diag_object_container_mod defines a container for inserting, removing and
!! searching for <TT>fms_diag_object</TT> instances
!!
!> @author Miguel Zuniga Miguel.Zuniga@noaa.gov
!!
!! <TT>fms_diag_object_container_mod</TT> defines a container for inserting, removing and
!! searching for <TT>fms_diag_object</TT> instances.
!!
!! This version uses the <TT>fms_doubly_linked_list_mod</TT> as the underlying
!! implementation. Future vesions may change the underlying implementation but
!! the user interface will not change.
!!
!> @file
!> @brief File for @ref fms_diag_object_container_mod
MODULE fms_diag_object_container_mod
   use fms_diag_object_mod, only: fms_diag_object
   USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE

   !! Since this version is based on the FDS linked list:
   use fms_doubly_linked_list_mod, only : linked_list_type, literator_type

   implicit none

   !> \brief A container of fms_diag_object instances providing insert, remove ,
   !!  find/search, and size (functions). Iterator is provided by the associated
   !! iterator class (see  dig_obj_iterator class).
   !!
   !!  This version does not enforce uniqueness of ID keys (I.e. it is not a set).
   !!
   type, public :: fms_diag_object_container_type
      TYPE (linked_list_type), ALLOCATABLE :: the_linked_list !!This version based on the FDS linked_list.
   contains
      procedure :: insert => insert_diag_object
      procedure :: remove => remove_diag_object
      procedure :: find   => find_diag_object
      procedure :: size   => get_num_objects
      procedure :: iterator  => get_iterator
      final :: destructor
   end type fms_diag_object_container_type


   !> \brief Iterator used to traverse the objects of the container.
   type, public :: diag_obj_iterator_type
      type(literator_type) :: liter !!this version based on the FDS linked_list
   contains
      procedure :: has_data => literator_has_data
      procedure :: next => literator_next
      procedure :: get => literator_data
   end type diag_obj_iterator_type

   interface fms_diag_object_container_type
      module procedure :: diag_object_container_constructor
   end interface fms_diag_object_container_type

   interface diag_obj_iterator_type
      module procedure :: diag_obj_iterator_constructor
   end interface diag_obj_iterator_type


contains

   !> \brief Returns an empty itetator if a diag object with this ID was not found.
   !! If the diag object was found, return an iterator with the current object being
   !! the found object, ad the last/anchor being the last/anchor of the container.
   !! Note that this routine can accept an optional iterator as input, which
   !! is useful for chaining searches, which may be needed if there are key duplicates.
   function find_diag_object (this, id , iiter) result (riter)
      class (fms_diag_object_container_type), intent (in out) :: this
      integer , intent (in) :: id   !> The id of the object to find.
      class(diag_obj_iterator_type), intent (in), optional :: iiter !> Iterator over the searchable set, if provided.
      class(diag_obj_iterator_type) , allocatable :: riter !> The resultant iterator to the object.
      class(fms_diag_object),  allocatable , target:: tdo  !> A temporaty diagnostic object
      integer :: status                                    !> Status from iterator operations.
      !!
      if(present (iiter)) then
         riter = iiter
      else
         riter = this%iterator() !!An iterator over the entire container
      endif
      do while( riter%has_data() .eqv. .true.)
         tdo = riter%get()
         if(id == tdo%get_id() ) then
            EXIT
         end if
         status = riter%next()
      end do
   end function find_diag_object

   !> \brief insert diagnostic object obj with given id.
   !! Objects are inserted at the back / end of the list
   !! This version of the container also enforces that the
   !! objects ID is equal the input id. status of -1
   !! is returned from errors
   function insert_diag_object (this, id, obj) result (status)
      class (fms_diag_object_container_type), intent (in out) :: this
      integer,  intent (in) :: id                     !> The id of the object to insert.
      class(fms_diag_object) , intent (in out) :: obj !> The object to insert
      integer :: status                               !> The returned status. 0 for success.
      class(literator_type), allocatable ::  tliter   !> A temporary iterator.

      status = -1
      if ( id .ne. obj%get_id() ) then
         !!TODO: log error
      endif
      tliter =  this%the_linked_list%push_back( obj )
      if(tliter%has_data() .eqv. .true. ) then
         status = 0
      endif
   end function

   !> \brief Remove and return the first object in the container with the corresponding id .
   !! Return an iterator that starts with the next object on the container. (Note, if the client
   !! code does not already have a referenc to the object being removed, then  the client may
   !! want to to use subroutine find before using subroutine remove. If subroutine find is used,
   !! consider then calling remove with the iterator returned from find.)
   function remove_diag_object (this, id, iiter ) result (riter)
      class (fms_diag_object_container_type), intent (in out) :: this
      integer , intent (in) :: id !> The id of the object to remove.
      class(diag_obj_iterator_type), intent (in), optional :: iiter !> Iterator over the searchable set, if provided.
      class(diag_obj_iterator_type), allocatable:: riter !> The resultant iterator
      class(literator_type), allocatable :: temp_liter !> A temporary iterator
      integer :: status
      if(present (iiter)) then
        riter = iiter
      else
        riter = this%iterator() !!An iterator over the entire container
      endif
      !Find the object in the container.
      riter = this%find ( id , riter)
      temp_liter = this%the_linked_list%remove( riter%liter%current )
      riter = diag_obj_iterator_type(temp_liter)
   end function

   !> return the number of objects held in the container
   function get_num_objects (this )result (sz)
      class (fms_diag_object_container_type), intent (in out) :: this
      integer :: sz                     !> The returned result - the number of objects in container.
      sz = this%the_linked_list%size()
   end function


   !> return a pointer to an iterator for the objects in the container.
   function get_iterator (this) result (oliter)
      class (fms_diag_object_container_type), intent (in) :: this
      class(diag_obj_iterator_type) , allocatable :: oliter !> The reurned iterator to the objects in the container.
      oliter = diag_obj_iterator_type( this%the_linked_list%get_literator() )
   end function

  !> A consructor for a container's iterator.
   function diag_obj_iterator_constructor( iliter ) result (diag_itr)
    class (literator_type), allocatable :: iliter !> An iterator. Normally the one that the container is based on.
    class (diag_obj_iterator_type), allocatable :: diag_itr !> The returned diag object iterator.
    allocate(diag_itr)
    diag_itr%liter = iliter;
 end function diag_obj_iterator_constructor

 !> A the default consructor for the container.
   function diag_object_container_constructor () result (doc)
      type(fms_diag_object_container_type), allocatable :: doc !> The resultant container.
      allocate(doc)
      doc%the_linked_list = linked_list_type()
      !! print * , "In DOC constructor"
   end function diag_object_container_constructor


   !> Return true if there is more data this iterator can access
   function literator_has_data( this ) result( r )
      class(diag_obj_iterator_type), intent(in) :: this
      logical :: r                  !> True if this iterator has data,
      r = this%liter%has_data()
   end function literator_has_data

   !> Move the iterator to the next data
   function literator_next( this ) result( status )
      class(diag_obj_iterator_type), intent(in out ) :: this
      integer :: status !> !> Returned status of operation. 0 if the iterator can be moved to the next object.
      status = 0
      status = this%liter%next()
   end function literator_next

   !> Get the current data the iterator is poiting to
   function  literator_data( this ) result( rdo )
      class(diag_obj_iterator_type), intent(in) :: this
      class(fms_diag_object),  pointer :: rdo !> The resultant diagnostic object.
      class(*),  pointer :: gp !> A eneric typed object in the container.

      gp => this%liter%get()
      select type(gp)
       type is (fms_diag_object)  !! "type is", not the (polymorphic) "class is"
         rdo => gp
       class default
          CALL  error_mesg ('diag_object_container:literator_data', &
            'Data to be accessed via iterator is not of expected type.',FATAL)
      end select
   end function literator_data

   !> The destroctor for the container
   subroutine destructor(this)
      type(fms_diag_object_container_type) :: this  !Note for destructors its needs to be type and not class!
      deallocate(this%the_linked_list)
   end subroutine destructor


end module fms_diag_object_container_mod
!> @}
! close documentation grouping

