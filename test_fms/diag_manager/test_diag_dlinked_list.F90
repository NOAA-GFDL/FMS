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

!! <TT>fms_diag_dlinked_list_mod</TT>  defines a generic doubly linked
!! list class  and an associated iterator class  for traversing the list. It
!! is generic in the sense that the elements or objects it contains are
!! "class(*)" objects. Note the public interface functions and the lack
!! of a search (or find) function as per the definition of a linked list.
!! If a search function, additional type cheeking,  or possibly a
!! slightly different user interface is desired, then consider creating
!! another iterator and another wrapper, or another class with this one for
!! a member element and procedures that are trivially implemented by using
!! this class. (See, for example, class FmsDiagObjectContainer_t and its
!! associated iterator.
!!
!! This version is roughly a Fortran translation of the C++ doubly linked list
!! class in the book ``Data Structures And Algorithm Analysis in C++",
!! 3rd Edition, by Mark Allen Weiss.
program test_diag_dlinked_list
   use mpp_mod, only: mpp_init, mpp_set_stack_size, mpp_init_test_requests_allocated
   use fms_mod, ONLY: error_mesg, FATAL,NOTE
   use fms_diag_object_mod, only : fmsDiagObject_type
   use fms_diag_dlinked_list_mod, only : FmsDlList_t, FmsDllIterator_t

   implicit  none

   !> @brief  This class is the type for the data to insert in the linked list.
   type TestDummy_t
      integer :: id  = 0
      real    :: weight = 1000
   end type TestDummy_t

   !!
   type (FmsDlList_t), allocatable :: list !< Instance of the linked list
   class(FmsDllIterator_t), allocatable :: iter !< An iterator for the list
   type (TestDummy_t), pointer::   p_td_obj     !< A pointer to a test_dummy object
   class(*), pointer :: p_obj                   !< A pointer to a class(*) object
   integer, parameter :: num_objs = 40          !< Total number of objects tested
   integer ::  full_id_sum                      !< Sum of all the possible object id values
   integer :: sum                               !< Temp sum of vaalues of id sets
   !!
   integer :: ierr                              !< An error flag
   logical :: test_passed                       !< Flag indicating if the test_passed
   !! These fields below used to initialize diag object data. TBD
   integer :: id
   !!

   call mpp_init(mpp_init_test_requests_allocated)
   call mpp_set_stack_size(145746)

   call error_mesg("test_diag_linked_list", "Starting tests",NOTE)

   test_passed = .true.  !! will be set to false if there are any issues.

   !! Ids will initially be from 1 to num_objs, so :
   full_id_sum = (num_objs * (num_objs + 1)) / 2

   !! Create the list
   allocate(list)
   call list%initialize()

   if( list%size() /= 0) then
     test_passed = .false.
     call error_mesg("test_diag_linked_list", "list incorrect size. Expected 0 at start",FATAL)
   endif

  !! Initialize num_objs objects and insert into list one at a time.
  !! The loop iterator is same as id - created in order to facilitate
   !! some tests.
   do id = 1, num_objs
      !!Allocate on heap another test dummy object :
      allocate (p_td_obj)
      !! And set some of its dummy data :
      p_td_obj%id = id
      p_td_obj%weight = id + 1000
      !! And have the "Char(*) pointer also point to it:
      p_obj => p_td_obj

      !! Test insertion the common way :
      iter = list%push_back( p_obj)
      if(iter%has_data() .eqv. .false. ) then
         test_passed = .false.
         call error_mesg("test_diag_dlinked_list", "List push_back error.",FATAL)
      endif

   enddo

   if( list%size() /= num_objs) then
      test_passed = .false.
      call error_mesg("test_diag_dlinked_list", "List has incorrect size after inserts.",FATAL)
   endif


   !! Test iteration over the entire list :
   sum = 0
   sum = sum_ids_in_list ( list )

   if( sum  /=  full_id_sum) then
      test_passed = .false.
      call error_mesg("test_diag_dlinked_list", &
        &"Id sums via iteration over the list objects is not as expected",FATAL)
   endif

   if( list%size() /= num_objs) then
      test_passed = .false.
      call error_mesg("test_diag_dlinked_list", &
        &"The list size is not as expected post inserts.",FATAL)
   endif

   !! Test a removal from the back (id should be num_objs)
   p_obj => find_back_of_list( list)
   iter  = list%pop_back()
   !! Note the client is resposible for managing memory of anything he explicitly
   !! removes from the list:
   deallocate(p_obj)
   sum = sum_ids_in_list ( list )
   if( sum  /=  full_id_sum - num_objs ) then
      test_passed = .false.
      call error_mesg("test_diag_dlinked_list", &
        &"Id sums via iteration over the list objects is not as expected",FATAL)
   endif

   !! Repeat - test removal from the back of list (should be (num_objs -1)).
   p_obj => find_back_of_list( list)
   iter = list%pop_back()
   !! Note the client is resposible for managing memory of anything he explicitly
   !! removes from the list:
   deallocate(p_obj)
   sum = sum_ids_in_list ( list )
   if( sum  /=  (full_id_sum - num_objs - (num_objs -1) )) then
      test_passed = .false.
      call error_mesg("test_diag_dlinked_list", &
        & "Id sums via iteration over the list objects is not as expected",FATAL)
   endif

   !! List.clear() is called by the destructor automatically, but for further testing
   !! we will use it to renove (and deallocate) the data nodes and associated data
   !! of the list.
   call list%clear()
   if( list%size() /= 0) then
      test_passed = .false.
      call error_mesg("test_diag_dlinked_list", &
        "List is incorrect size after clearing.",FATAL)
   endif

   !! Allocated objects are deallocated automatically, but one can aslo make the call.
   deallocate(list)

   call error_mesg('test_diag_dlinked_list', 'Test has finished',NOTE)

   call MPI_finalize(ierr)

CONTAINS

   !> @brief Cast the "class(*) input data to the expected type.
   function  get_typed_data( pci ) result( pdo )
      class(*), intent(in), pointer :: pci !< An input pointer to the class(*) data object.
      class(TestDummy_t),  pointer :: pdo !< The resultant pointer to the expected underlying object type.
      !
      pdo => null()
      select type(pci)
       type is (TestDummy_t)  !! "type is", not the (polymorphic) "class is"
         pdo => pci
       class default
         call error_mesg("test_diag_dlinked_list", &
           & "Data to access is not of expected type.",FATAL)
      end select
   end function get_typed_data

   !> Calcualte the sum of the ids.
   !! Exercises iteration over the list.
   function sum_ids_in_list (the_list) result (rsum)
      type (FmsDlList_t), intent(inout) , allocatable :: the_list !< The linked list instance
      integer  :: rsum                        !< The resultant sum of ids
      class(FmsDllIterator_t), allocatable :: iter !< An iterator over the list
      type (TestDummy_t), pointer:: p_td_obj => null() !< A pointer to a test_dummy object
      class(*), pointer :: p_obj  => null()            !< A pointer to a class(*) object
      integer :: ic_status                             !< A list insertion status.
      !!
      rsum = 0
      iter = the_list%get_literator()
      do while( iter%has_data() .eqv. .true.)
         p_obj => iter%get()
         p_td_obj => get_typed_data (p_obj )
         rsum = rsum + p_td_obj%id
         ic_status = iter%next()
      end do
   end function sum_ids_in_list

  !> Find the past object in list. This also is a kind of search function,
   !! so if the provided wrapper is not used, you have to write your own.
   !! @return a pointer the object at the end of the list, or null if none
   function find_back_of_list (the_list) result (pdo)
    type (FmsDlList_t), intent(inout) , allocatable ::the_list  !< The linked list instance
    class(TestDummy_t),  pointer :: pdo !< The resultant back of list,
    class(FmsDllIterator_t), allocatable :: iter !< An iterator over the list
    class(*), pointer :: p_obj => null()         !< A pointer to a class(*) object
    integer :: ic_status                         !< A list insertion status.
    !!
    pdo=>null()
    iter = the_list%get_literator()
    do while( iter%has_data() .eqv. .true.)
       p_obj => iter%get()
       pdo => get_typed_data (p_obj )
       ic_status = iter%next()
    end do
 end function find_back_of_list

end program test_diag_dlinked_list
