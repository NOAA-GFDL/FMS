module fms_diag_object_container_mod
   use fms_diag_object_mod, only: fms_diag_object

   !This version is based on the FDS linked list:
   use fms_doubly_linked_list_mod, only : linked_list, literator

   !If the FTL is the underlying container.
   !use ftlHashMapIntDiagoModule
   !type(ftlHashMapIntDiago) :: object_container


   implicit none

   !> \brief A container of diag objects providing insert, remove , find/search, and size (functions).
   !  Iterator is provided by the associated iterator class (see  dig_obj_iterator class).
   !
   !  This version does not enforce uniqueness of ID keys (I.e. it is not a set).
   !
   !  This class is just a wrapper for the underlying implementation, but it provides
   !  simple and easy to use type safety for the end user.
   !! The underlying container implementation is TBD. Candidates: FTL, gFTL, FDS list, FDS AA-tree ?
   type, public :: fms_diag_object_container
      type (linked_list), allocatable :: the_ll !!This version based on the FDS linked_list.
   contains
      procedure :: insert => insert_diag_object
      procedure :: remove => remove_diag_object
      procedure :: find   => find_diag_object
      procedure :: size   => get_num_objects
      procedure :: iterator  => get_iterator
      final :: destructor
   end type fms_diag_object_container


   !> \brief Iterator used to traverse the objects of the container.
   type, public :: diag_obj_iterator
      type(literator) :: liter !!this version based on the FDS linked_list
   contains
      procedure :: has_data => literator_has_data
      procedure :: next => literator_next
      procedure :: get => literator_data
   end type diag_obj_iterator

   interface fms_diag_object_container
      module procedure :: diag_object_container_constructor
   end interface fms_diag_object_container

   interface diag_obj_iterator
      module procedure :: diag_obj_iterator_constructor
   end interface diag_obj_iterator


   !! A module variable instance (a singleton) for the entire executable:
   type (fms_diag_object_container) :: the_object_container


contains

   !> \brief Returns an empty itetator if a diag object with this ID was not found.
   ! If the diag object was found, return an iterator with the current object being
   ! the found object, ad the last/anchor being the last/anchor of the container.
   ! Note that this routine can accept an optional iterator as input, which
   ! is useful for chaining searches, which may be needed if there are key duplicates.
   function find_diag_object (this, id , iiter) result (riter)
      class (fms_diag_object_container), intent (in out) :: this
      integer , intent (in) :: id
      class(diag_obj_iterator), intent (in), optional :: iiter
      class(diag_obj_iterator) , allocatable :: riter
      !!
      class(fms_diag_object),  allocatable , target:: tdo !!a temporaty diagnostic object
      integer :: status
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
   ! Objects are inserted at the back / end of the list
   ! This version of the container also enforces that the
   ! objects ID is equal the input key. status of -1
   ! is returned fro errors
   !! TODO: what should be the reutrn type
   !!TODO: WE ASSUME id == obj%diagID
   function insert_diag_object (this, id, obj) result (status)
      class (fms_diag_object_container), intent (in out) :: this
      integer,  intent (in) :: id
      class(fms_diag_object) , intent (in out) :: obj
      integer :: status
      !!
      class(literator), allocatable ::  tliter

      status = -1
      if ( id .ne. obj%get_id() ) then
         !!TODO: log error
      endif
      tliter =  this%the_ll%push_back( obj )
      if(tliter%has_data() .eqv. .true. ) then
         status = 0
      endif
   end function

   !> \brief Remove and return the first object in the container with the corresponding key .
   ! Return an iterator that starts with the next object on the container. (Note, if the client
   ! code does not already have a referenc to the object being removed, then  the client may
   ! want to to use subroutine find before using subroutine remove. If subroutine find is used,
   ! consider then calling remove with the iterator returned from find.)
   function remove_diag_object (this, key, iiter ) result (riter)
      class (fms_diag_object_container), intent (in out) :: this
      integer , intent (in) :: key
      class(diag_obj_iterator), intent (in), optional :: iiter
      class(diag_obj_iterator), allocatable:: riter
      integer :: status
      !!
      !Find the object in the container.
      riter = this%find ( key , iiter)
      riter = this%the_ll%remove( riter%liter%current )
   end function

   !> return the number of objects held in the container
   function get_num_objects (this )result (sz)
      class (fms_diag_object_container), intent (in out) :: this
      integer :: sz
      sz = this%the_ll%size()
   end function


   !> return a pointer to an iterator for the objects in the continer.
   function get_iterator (this) result (oliter)
      class (fms_diag_object_container), intent (in) :: this
      class(diag_obj_iterator) , allocatable :: oliter
      oliter = diag_obj_iterator()
      oliter%liter = this%the_ll%get_literator()
   end function

   function diag_obj_iterator_constructor ( ) result (litr)
      class (diag_obj_iterator), allocatable :: litr
      allocate(litr)
   end function diag_obj_iterator_constructor

   function diag_object_container_constructor () result (doc)
      type(fms_diag_object_container), allocatable :: doc
      allocate(doc)
      doc%the_ll = linked_list()
      print * , "In DOC constructor"
   end function diag_object_container_constructor



   !> Return true if there is more data this iterator can access
   function literator_has_data( this ) result( r )
      class(diag_obj_iterator), intent(in) :: this
      logical :: r
      r = this%liter%has_data()
   end function literator_has_data

   !> Move the iterator to the next data
   function literator_next( this ) result( status )
      class(diag_obj_iterator), intent(in out ) :: this
      integer :: status
      status = 0
      status = this%liter%next()
   end function literator_next

   !> Get the current data the iterator is poiting to
   function  literator_data( this ) result( rdo )
      class(diag_obj_iterator), intent(in) :: this
      class(fms_diag_object),  pointer :: rdo !!the resultant diagnostic object
      class(*),  pointer :: gp !! generic typed object in the container

      gp => this%liter%get()
      select type(gp)
       type is (fms_diag_object)  !! "type is", not the (polymorphic) "class is"
         rdo => gp
       class default
         print * , "error iter data type error"
      end select
   end function literator_data

   subroutine destructor(this)
      type(fms_diag_object_container) :: this  !Note for destructors its needs to be type and not class!
      print *, "Diag_obj_container . Starting destructor"
      deallocate(this%the_ll)
   end subroutine destructor


end module fms_diag_object_container_mod

