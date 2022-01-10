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
   use mpp_mod, only: mpp_init, mpp_exit, mpp_error, FATAL, WARNING
   use mpp_mod, only : mpp_set_stack_size, mpp_init_test_requests_allocated
   use mpp_io_mod, only: mpp_io_init

   use fms_diag_object_mod, only : fmsDiagObject_type
   use fms_diag_dlinked_list_mod, only : FmsDlList_t, FmsDllIterator_t

   implicit  none

   !> @brief  This class is the type for the data to insert in the linked list.
   type TestDummy_t
      integer :: id  = 0
      character(len=20) :: name
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
   character(:), allocatable :: mname, mname_pre
   !!


   test_passed = .true.  !! will be set to false if there are any issues.

   call mpp_init(mpp_init_test_requests_allocated)
   call mpp_io_init()
   call mpp_set_stack_size(145746)

   !! Ids will initially be from 1 to num_objs, so :
   full_id_sum = (num_objs * (num_objs + 1)) / 2

   !!Create the list
   list = FmsDlList_t()

   if( list%size() /= 0) then
      test_passed = .false.
      call mpp_error(FATAL, "list incorrect size. Expected 0 at start")
   endif
   mname_pre = "ATM"

  !! Initialize num_objs objects and insert into list one at a time.
  !! The loop iterator is same as id - created in order to facilitate
   !! some tests.
   do id = 1, num_objs
      !!Allocate on heap another test dummy object :
      allocate (p_td_obj)
      !! And set some of its dummy data :
      call combine_str_int(mname_pre, id, mname)
      p_td_obj%id = id
      p_td_obj%name = mname
      !! And have the "Char(*) pointer also point to it:
      p_obj => p_td_obj

      !! Test insertion the common way :
      iter = list%push_back( p_obj)
      if(iter%has_data() .eqv. .false. ) then
         test_passed = .false.
         call mpp_error(FATAL, "List push_back error.")
      endif

   enddo

   if( list%size() /= num_objs) then
      test_passed = .false.
      call mpp_error(FATAL, "List has incorrect size after inserts.")
   endif


   !! Test iteration over the entire list :
   sum = 0
   sum = sum_ids_in_list ( list )

   if( sum  /=  full_id_sum) then
      test_passed = .false.
      call mpp_error(FATAL, "Id sums via iteration over the list objects is not as expected")
   endif

   if( list%size() /= num_objs) then
      test_passed = .false.
      call mpp_error(FATAL, "The list size is not as expected post inserts.")
   endif

   !! Test a removal from the back (id should be num_objs)
   p_obj => find_back_of_list( list)
   iter = list%pop_back()
   !! Note the client is resposible for managing memory of anything he explicitly
   !! removes from the list:
   deallocate(p_obj)
   sum = sum_ids_in_list ( list )
   if( sum  /=  full_id_sum - num_objs ) then
      test_passed = .false.
      call mpp_error(FATAL, "Id sums via iteration over the list objects is not as expected")
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
      call mpp_error(FATAL, "Id sums via iteration over the list objects is not as expected")
   endif

   call list%clear()
   if( list%size() /= 0) then
      test_passed = .false.
      call mpp_error(FATAL, "List is incorrect size after clearing.")
   endif

   write (6,*) "Finishing diag_dlinked_list tests."

   !! the list has a finalize/destructor which will deallocate data that is still it list.
   !! equivalent to calling list%clear() as above.
   deallocate(list)

   call MPI_finalize(ierr)

CONTAINS



   !> @brief Cast the "class(*) input data to the expected type.
   function  get_typed_data( data_in  ) result( rdo )
      class(*), intent(in), pointer :: data_in !< An input pointer to the class(*) object.
      class(TestDummy_t),  pointer :: rdo !< The resultant pointer to the expected underlying object type.
      rdo => null()

      select type(data_in)
       type is (TestDummy_t)  !! "type is", not the (polymorphic) "class is"
         rdo => data_in
       class default
         call mpp_error(FATAL, "Data to access is not of expected type.",FATAL)
      end select
   end function get_typed_data

   !> Calcualte the sum of the ids.
   !! Exercises iteration over the list.
   function sum_ids_in_list (list) result (rsum)
      type (FmsDlList_t), allocatable :: list !< The linked list instance
      integer  :: rsum                        !< The resultant sum of ids
      class(FmsDllIterator_t), allocatable :: iter !< An iterator over the list
      type (TestDummy_t), pointer::   p_td_obj     !< A pointer to a test_dummy object
      class(*), pointer :: p_obj                   !< A pointer to a class(*) object
      integer :: ic_status                         !< A list insertion status.
      !!
      rsum = 0
      iter = list%get_literator()
      do while( iter%has_data() .eqv. .true.)
         p_obj => iter%get()
         p_td_obj => get_typed_data (p_obj )
         id =  p_td_obj%id
         rsum = rsum + id
         ic_status = iter%next()
      end do
   end function sum_ids_in_list

  !> Calcualate the sum of the ids. This also is a kind of search function,
   !! so if the provided wrapper is not used, you have to write your own.
   !! @return a pointer the object at the end of the list, or null if none
   function find_back_of_list (list) result (p_rdo)
    type (FmsDlList_t), allocatable :: list !< The linked list instance
    class(TestDummy_t),  pointer :: p_rdo !< The resultant back of list,
    class(FmsDllIterator_t), allocatable :: iter !< An iterator over the list
    class(*), pointer :: p_obj                   !< A pointer to a class(*) object
    integer :: ic_status                         !< A list insertion status.
    !!
    p_rdo => null()
    iter = list%get_literator()
    do while( iter%has_data() .eqv. .true.)
       p_obj => iter%get()
       p_rdo => get_typed_data (p_obj )
       ic_status = iter%next()
    end do
 end function find_back_of_list

 subroutine combine_str_int (str, num, rs)
    character(:), allocatable, intent (in):: str
    integer ,    intent (in) :: num
    character(:), allocatable, intent (out) :: rs
    character(len_trim(str) + 8) :: tmp

    write (tmp, "(A4,I4)") str,num
    tmp = trim(tmp)
    rs = tmp
 end subroutine combine_str_int


end program test_diag_dlinked_list


