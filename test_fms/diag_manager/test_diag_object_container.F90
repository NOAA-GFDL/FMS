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

!> @brief  This programs tests  public member functions of the
!!  FmsDiagObjectContainer_t and FmsDiagObjIterator_t. As these two classes
!!  are largely wrappers to their underlying classes, it is also
!!  testing the underlying container and iterator classes. The container
!!  functions being tested are insert, remove, and size. The use of the iterators
!!  is also being tested.
program test_diag_obj_container
  use mpp_mod, only: mpp_init, mpp_exit, mpp_error, FATAL, WARNING
  use mpp_mod, only : mpp_set_stack_size, mpp_init_test_requests_allocated
  use mpp_io_mod, only: mpp_io_init

  use fms_diag_object_mod, only : fmsDiagObject_type
  use fms_diag_object_container_mod, only : FmsDiagObjectContainer_t, FmsDiagObjIterator_t
  USE time_manager_mod, ONLY: time_type

  implicit  none
  !!
  type (FmsDiagObjectContainer_t), allocatable :: container !< Instance of the container
  class(FmsDiagObjIterator_t), allocatable :: iter          !< An iterator for the container
  type (fmsDiagObject_type), allocatable , target ::  obj_vec(:)     !< A vector of objects
  type (fmsDiagObject_type), pointer::   pobj                        !< A pointer to an object
  integer, parameter :: num_objs = 10                             !< Total number of objects tested
  integer ::  full_id_sum                                         !< Sum of all the possible object id values
  integer :: sum                                                  !< Temp sum of vaalues of id sets
  !!
  integer :: ic_status                                            !< A status flag returned from container functions
  integer :: ierr                                                 !< An error flag
  !!
  logical :: test_passed                                          !< Flag indicating if the test_passed
  !! These fields below used to initialize diag object data. TBD
  integer :: id
  integer, dimension(2) :: axes
  TYPE(time_type)  :: init_time
  !!type (diag_fields_type)  :: diag_field
  character(:), allocatable :: mname, vname, mname_pre, vname_pre
  !!


  test_passed = .true.  !! will be set to false if there are any issues.

  call mpp_init(mpp_init_test_requests_allocated)
  call mpp_io_init()
  call mpp_set_stack_size(145746)

  !! Ids will initially be from 1 to num_objs, so :
  full_id_sum = (num_objs * (num_objs + 1)) / 2

  !!Create the container
  container = FmsDiagObjectContainer_t()
  !!In diag_manager, one module level container may be used instead of a local one like above.


  !! Allocate some test objects.
  !! NOTE: normally objects will be allocated one at a time with a stament like:
  !!   allocate(pobj, source = fms_diag_object(argument list ))
  !! or via constructor like :
  !!   pobj => fms_diag_object(argument list )
  !! Once the object ID is set, it should be inserted into the container and then the
  !!   container will be considered the manager of that object and its memory (unless the object is removed).
  !! Since type fms_diag_obj doesn't have a proper constructor yet, well be lazy by making array of objects
  !! ( normal fixed size array the thing whose use we are replacing to begin with ) and consider these particular
  !! objects to not be managed by the container.
  allocate(obj_vec(num_objs))

  !! Initialize each object and isnert into container one at a time.

  if( container%size() /= 0) then
    test_passed = .false.
    call mpp_error(FATAL, "Container incorrect size. Expected 0 at start")
  endif
  mname_pre = "ATM"
  vname_pre = "xvar"
  do id = 1, num_objs
    call combine_str_int(mname_pre, id, mname)
    call combine_str_int(vname_pre, id, vname )

    pobj => obj_vec( id ) !!Note use of pointer to obj.
    call pobj%setID(id)

    call pobj%register ("test_mod", vname, axes, init_time, "a_long_name")

    !!Insert object into the container.
    ic_status = container%insert(pobj%get_id(), pobj)
    if(ic_status .ne. 0)then
      test_passed = .false.
      call mpp_error(FATAL, "Container Insertion error.")
    endif
  enddo

  if( container%size() /= num_objs) then
    test_passed = .false.
    call mpp_error(FATAL, "Container has incorrect size after inserts.")
  endif

  !!Search the container for a an object of specified key
  iter =  container%find(123)
  if ( iter%has_data() .eqv. .true. ) then
    test_passed = .false.
    call mpp_error(FATAL, "Found in container unexpected object of id=123")
  endif

  !!Again, search the container for a an object of specified key
  iter = container%find(4)
  if (iter%has_data() .neqv. .true. ) then
    test_passed = .false.
    call mpp_error(FATAL, "Did not find expected container object of id=4")
  endif

  !! Iterate over all the objects in the container;
  sum = 0
  iter = container%iterator()
  do while( iter%has_data() .eqv. .true.)
      pobj => iter%get()  !!Note use of pointer and pointer assignment is preferred.
      id =  pobj%get_id( )
      !! vname =  pobj%get_varname()  !! print ...
      sum = sum + id
      ic_status = iter%next()
  end do

  if( sum  /=  full_id_sum) then
    test_passed = .false.
    call mpp_error(FATAL, "Id sums via iteration over the container objects is not as expected")
  endif

  if( container%size() /= num_objs) then
    test_passed = .false.
    call mpp_error(FATAL, "The container size is not as expected post inserts.")
  endif


  !! Test a removal ****
  iter = container%iterator()
  iter  = container%remove( 4, iter )
  iter  = container%find(4)
  !! Verify  the removal , part 1:
  if (  iter%has_data() .eqv. .true.) then
    test_passed = .false.
    call mpp_error(FATAL, "Found object of id = 4 after removing it")
  endif
   !! Verify  the removal , part 2 :
  if (container%size() /= (num_objs - 1)) then
     test_passed = .false.
    call mpp_error(FATAL,"The_container%size() \= num_obj -1 after a removal ")
  endif

   !! Verify  the removal , part 3 :
   !! Iterate over all the objects in the container AFTER the removal of id=4 object;
  sum = 0
  iter = container%iterator()
  do while( iter%has_data() .eqv. .true.)
      pobj => iter%get()  !!Note use of pointer and pointer assignment is preferred.
      id =  pobj%get_id( )
      !! vname =  pobj%get_varname()  !! print ...
      sum = sum + id
      ic_status = iter%next()
  end do
  if( sum  /=  full_id_sum - 4) then
    test_passed = .false.
    call mpp_error(FATAL, "Container incorrect id sums post removal of 4")
  endif
  !! End test a removal ****

  !! Test find and access object in the container
  iter = container%find(7)
  if (iter%has_data() .neqv. .true. ) then
    test_passed = .false.
    call mpp_error(FATAL, "Container did not find object of id=7")
  endif
  !! Check the find results more :
  pobj => iter%get()
  if(pobj%get_id() /=  7) then
    test_passed = .false.
    call mpp_error(FATAL," Id of returned object was not 7 ")
  endif
  !!TODO further access tests.


  !! Manually clear out the container.
  !! NOTE: In normal use this is NOT PERFORMED since with its finalize function, the  container
  !! deallocates all pointers and data it manages. However, the client needs to take care of
  !! the diag objects the client has decided that the container should not manage.
  !! In this wierd test case, all the diag objects were originally from a vector (a container itself!)
  !! and not allocated on the heap one at a time, so this step is needed before program completion.
  do id = 1, num_objs
    iter  = container%find(id)
    if (  iter%has_data() .eqv. .true.) then
      iter  = container%remove( id, iter )
    endif
  end do

  if( container%size() /= 0) then
    test_passed = .false.
    call mpp_error(FATAL, "Container is incorrect size after clearing.")
  endif

  write (6,*) "Finishing diag_obj_container tests."

  !! the container has a finalize/destructor which will
deallocate(container)

call MPI_finalize(ierr)

CONTAINS

subroutine combine_str_int (str, num, rs)
  character(:), allocatable, intent (in):: str
  integer ,    intent (in) :: num
  character(:), allocatable, intent (out) :: rs
  character(len_trim(str) + 8) :: tmp

  write (tmp, "(A4,I4)") str,num
  tmp = trim(tmp)
  rs = tmp
end subroutine combine_str_int

end program test_diag_obj_container


