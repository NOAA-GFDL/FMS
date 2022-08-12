!> @author Ryan Mulhall
!> @email ryan.mulhall@noaa.gov
!! @brief Contains buffer types and routines for the diag manager
!!
!! @description Holds buffered data for fmsDiagVars_type objects
!! buffer0-5d types extend fmsDiagBuffer_class, and upon allocation
!! are added to the module's buffer_lists depending on it's dimension
module fms_diag_buffer_mod

use platform_mod
use iso_c_binding
use fms_diag_axis_object_mod, only: diagDomain_t
use time_manager_mod, only: time_type
use mpp_mod, only: mpp_error, FATAL
use diag_data_mod, only: DIAG_NULL

implicit none

private

!> @brief Object that holds buffered data and other diagnostics
!! Abstract to ensure use through its extensions(buffer0-5d types)
type, abstract :: fmsDiagBuffer_class
    class(*), dimension(:,:,:,:,:), allocatable :: remap_buffer !< remapped buffer data
    integer, allocatable, private               :: buffer_id !< index in buffer list
    type(time_type), private                    :: init_time !< initialization time
    character(len=128), allocatable, private    :: interp_method
    integer                                     :: frequency
    integer, allocatable, private               :: tile_count
    integer, pointer, dimension(:), private     :: axis_ids
    class(diagDomain_t), pointer, private       :: domain
    integer, allocatable, private               :: area, volume
    class(*), allocatable, private              :: missing_value
    class(*), allocatable, private              :: data_RANGE(:)

    contains

    procedure :: get_remapped_buffer_pointer
    procedure :: get_area
    procedure :: get_volume
    procedure :: get_missing_value
    procedure :: get_data_RANGE
    procedure :: flush_buffer
    procedure :: remap_buffer_data => remap_buffer
    ! TODO deferred routines, will require some interfaces 
    !procedure(allocate_buffer), deferred :: allocate_buffer
    !procedure, deferred :: get_buffer 
    !procedure, deferred :: initialize_buffer

end type fmsDiagBuffer_class

!> holds pointer to allocated buffer0-5d objects
type :: fmsDiagBufferContainer_type
    class(fmsDiagBuffer_class), pointer :: buffer_obj
end type

!> Scalar buffer type to extend fmsDiagBufferContainer_type
type, extends(fmsDiagBuffer_class) :: buffer0d
    class(*), allocatable :: buffer
    contains
    procedure :: allocate_buffer => allocate_buffer_0d
    procedure :: get_buffer_data => get_buffer_0d
    procedure :: initialize_buffer => initialize_buffer_0d
    procedure :: add_to_buffer => add_to_buffer_0d

end type buffer0d

!> 1D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer1d
    class(*), allocatable :: buffer(:)
    contains
    ! TODO add
    procedure :: allocate_buffer => allocate_buffer_1d
    procedure :: get_buffer_data => get_buffer_1d
    procedure :: initialize_buffer => initialize_buffer_1d
    procedure :: add_to_buffer => add_to_buffer_1d
end type buffer1d

!> 2D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer2d
    class(*), allocatable :: buffer(:,:)
    contains
    ! TODO add
    procedure :: allocate_buffer => allocate_buffer_2d
    procedure :: get_buffer_data => get_buffer_2d
    procedure :: initialize_buffer => initialize_buffer_2d
    procedure :: add_to_buffer => add_to_buffer_2d
end type buffer2d

!> 3D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer3d
    class(*), allocatable :: buffer(:,:,:)
    contains
    ! TODO add 
    procedure :: allocate_buffer => allocate_buffer_3d
    procedure :: get_buffer_data => get_buffer_3d
    procedure :: initialize_buffer => initialize_buffer_3d
    procedure :: add_to_buffer => add_to_buffer_3d
end type buffer3d

!> 4D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer4d
    class(*), allocatable :: buffer(:,:,:,:)
    contains
    !! TODO add
    procedure :: allocate_buffer => allocate_buffer_4d
    procedure :: get_buffer_data => get_buffer_4d
    procedure :: initialize_buffer => initialize_buffer_4d
    procedure :: add_to_buffer => add_to_buffer_4d
end type buffer4d

!> 5D buffer type to extend fmsDiagBuffer_class
type, extends(fmsDiagBuffer_class) :: buffer5d
    class(*), allocatable :: buffer(:,:,:,:,:)
    contains
    ! TODO  add
    procedure :: allocate_buffer => allocate_buffer_5d
    procedure :: get_buffer_data => get_buffer_5d
    procedure :: initialize_buffer => initialize_buffer_5d
    procedure :: add_to_buffer => add_to_buffer_5d
end type buffer5d

! public types
public :: buffer0d
public :: buffer1d
public :: buffer2d
public :: buffer3d
public :: buffer4d
public :: buffer5d
public :: fmsDiagBuffer_class
public :: fmsDiagBufferContainer_type
! public routines
public :: fms_diag_buffer_init 
! Module variables
logical,private :: module_is_initialized = .false. !< Flag indicating if the module is initialized
integer, parameter :: DEFAULT_BUFFER_LIST_SIZE = 64
integer         :: num_buffers = 0 !< number of buffers allocated, might be better in diag_object


logical, parameter, private :: DEBUG = .true. !< debugging output

! TODO interface with function declarations needed for deferred routines
!abstract interface
    !integer function allocate_buffer(buffobj, mold) &
    !result(rslt)
        !import buffer0d 
        !!import time_type 
        !class(fmsDiag), intent(inout) :: buffobj !< scalar buffer object
        !class(*),intent(in) :: mold !< allocates to the type of mold
        !!type(time_type), intent(in), optional :: init_time
    !end function
!end interface
!
contains


!!--------module routines

logical function fms_diag_buffer_init(buffobjs)
    class(fmsDiagBufferContainer_type), allocatable, intent(out) :: buffobjs(:)
    allocate(buffobjs(DEFAULT_BUFFER_LIST_SIZE))
    fms_diag_buffer_init = allocated(buffobjs)
end function

!> creates a container type with the given buffer object
function fms_diag_buffer_create_container(buffobj) &
result(rslt)
    class(fmsDiagBuffer_class), allocatable, target :: buffobj
    class(fmsDiagBufferContainer_type), allocatable :: rslt
    allocate(rslt)
    rslt%buffer_obj => buffobj
end function 

!!--------generic routines for any fmsDiagBuffer_class objects 

!> Remaps 0-5d data buffer from the given object onto a 5d array and sets it in the type
subroutine remap_buffer(buffobj)
    class(fmsDiagBuffer_class), intent(inout) :: buffobj
    integer                                   :: buff_bounds(5)

    if( DEBUG) print *, 'remapping buffer'
    ! get num dimensions from type extension
    select type (buffobj)
        type is (buffer0d)
            if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
            if (DEBUG) print *, '0d buffer'
            ! get buffer data type to allocate and remap
            select type (buff => buffobj%buffer)
                type is(integer(i4_kind))
                    allocate (integer(kind=i4_kind) :: buffobj%remap_buffer(1,1,1,1,1))
                    buffobj%remap_buffer = buff
                type is(integer(i8_kind))
                    allocate (integer(kind=i8_kind) :: buffobj%remap_buffer(1,1,1,1,1))
                    buffobj%remap_buffer = buff
                type is(real(r4_kind))
                    allocate (real(kind=r4_kind) :: buffobj%remap_buffer(1,1,1,1,1))
                    buffobj%remap_buffer = buff
                type is(real(r8_kind))
                    allocate (real(kind=r8_kind) :: buffobj%remap_buffer(1,1,1,1,1))
                    buffobj%remap_buffer = buff 
            end select
        type is (buffer1d)
            if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
            if( DEBUG) print *, '1d buffer'
            buff_bounds(1) = size(buffobj%buffer, 1)
            select type (buff => buffobj%buffer)
                type is(integer(i4_kind))
                    allocate (integer(kind=i4_kind) :: buffobj%remap_buffer(buff_bounds(1),1,1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
                type is(integer(i8_kind))
                    allocate (integer(kind=i8_kind) :: buffobj%remap_buffer(buff_bounds(1),1,1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
                type is(real(r4_kind))
                    allocate (real(kind=r4_kind) :: buffobj%remap_buffer(buff_bounds(1),1,1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
                type is(real(r8_kind))
                    allocate (real(kind=r8_kind) :: buffobj%remap_buffer(buff_bounds(1),1,1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),1,1,1,1 /))
            end select
        type is (buffer2d)
            if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
            buff_bounds(1) = size(buffobj%buffer, 1)
            buff_bounds(2) = size(buffobj%buffer, 2)
            select type (buff => buffobj%buffer)
                type is(integer(i4_kind))
                    allocate (integer(kind=i4_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
                type is(integer(i8_kind))
                    allocate (integer(kind=i8_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
                type is(real(r4_kind))
                    allocate (real(kind=r4_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
                type is(real(r8_kind))
                    allocate (real(kind=r8_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2),1,1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),1,1,1 /))
            end select
        type is (buffer3d)
            if (.not. allocated(buffobj%buffer)) call mpp_error(FATAL, "remap_buffer: buffer data not yet allocated")
            buff_bounds(1) = size(buffobj%buffer, 1)
            buff_bounds(2) = size(buffobj%buffer, 2)
            buff_bounds(3) = size(buffobj%buffer, 3)
            select type (buff => buffobj%buffer)
                type is(integer(i4_kind))
                    allocate (integer(kind=i4_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2), buff_bounds(3),1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
                type is(integer(i8_kind))
                    allocate (integer(kind=i8_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
                type is(real(r4_kind))
                    allocate (real(kind=r4_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
                type is(real(r8_kind))
                    allocate (real(kind=r8_kind) :: buffobj%remap_buffer(buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1))
                    buffobj%remap_buffer = RESHAPE( buff, (/ buff_bounds(1),buff_bounds(2),buff_bounds(3),1,1 /))
            end select
        type is (buffer4d)
        type is (buffer5d)
        class default
            call mpp_error( FATAL, 'remap_buffer_pointer: invalid buffer type for remapping')
    end select
   
end subroutine 

!> @brief Gets the remapped buffer pointer
!! Will do the remapping if not set already
function get_remapped_buffer_pointer (obj) &
result(rslt)
    class(fmsDiagBuffer_class), target, intent(inout) :: obj
    class(fmsDiagBuffer_class), pointer       :: objptr
    class(*), pointer, dimension(:,:,:,:,:)  :: rslt
    if(.not. allocated(obj%remap_buffer)) then
        call remap_buffer(obj) 
    endif
    objptr => obj
    rslt => objptr%remap_buffer
end function

!> @brief Gets area
!! @return copy of the area or diag_null if not allocated
pure function get_area (obj) &
result(rslt)
     class (fmsDiagBuffer_class), intent(in) :: obj !< diag object
     integer :: rslt
     if (allocated(obj%area)) then
       rslt = obj%area
     else
       rslt = diag_null
     endif
end function get_area
!> @brief Gets volume
!! @return copy of the volume or diag_null if volume is not allocated
pure function get_volume (obj) &
result(rslt)
     class (fmsDiagBuffer_class), intent(in) :: obj !< diag object
     integer :: rslt
     if (allocated(obj%volume)) then
       rslt = obj%volume
     else
       rslt = diag_null
     endif
end function get_volume
!> @brief Gets missing_value
!! @return copy of The missing value
function get_missing_value (obj) &
result(rslt)
     class (fmsDiagBuffer_class), intent(in) :: obj !< diag object
     class(*),allocatable :: rslt
     if (allocated(obj%missing_value)) then
       select type (miss => obj%missing_value)
         type is (integer(kind=i4_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = miss
         type is (integer(kind=i8_kind))
             allocate (integer(kind=i8_kind) :: rslt)
             rslt = miss
         type is (real(kind=r4_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = miss
         type is (real(kind=r8_kind))
             allocate (integer(kind=i4_kind) :: rslt)
             rslt = miss
         class default
             call mpp_error ("get_missing_value", &
                     "The missing value is not a r8, r4, i8, or i4",&
                     FATAL)
         end select
       else
         call mpp_error ("get_missing_value", &
                 "The missing value is not allocated", FATAL)
       endif
end function get_missing_value
!> @brief Gets data_range
!! @return copy of the data range
function get_data_RANGE (obj) &
result(rslt)
     class (fmsDiagBuffer_class), intent(in) :: obj !< diag object
     class(*),allocatable :: rslt(:)
     if (allocated(obj%data_RANGE)) then
       select type (r => obj%data_RANGE)
         type is (integer(kind=i4_kind))
             allocate (integer(kind=i4_kind) :: rslt(2))
             rslt = r
         type is (integer(kind=i8_kind))
             allocate (integer(kind=i8_kind) :: rslt(2))
             rslt = r
         type is (real(kind=r4_kind))
             allocate (integer(kind=i4_kind) :: rslt(2))
             rslt = r
         type is (real(kind=r8_kind))
             allocate (integer(kind=i4_kind) :: rslt(2))
             rslt = r
         class default
             call mpp_error ("get_data_RANGE", &
                     "The data_RANGE value is not a r8, r4, i8, or i4",&
                     FATAL)
         end select
       else
         call mpp_error ("get_data_RANGE", &
                 "The data_RANGE value is not allocated", FATAL)
       endif
end function get_data_RANGE

subroutine flush_buffer(buffobj)
    class(fmsDiagBuffer_class), intent(inout) :: buffobj
    select type (buffobj)
        type is (buffer0d)
            deallocate(buffobj%buffer)
        type is (buffer1d)
            deallocate(buffobj%buffer)
        type is (buffer2d)
            deallocate(buffobj%buffer)
        type is (buffer3d)
            deallocate(buffobj%buffer)
        type is (buffer4d)
            deallocate(buffobj%buffer)
        type is (buffer5d)
            deallocate(buffobj%buffer)
    end select
    deallocate(buffobj%missing_value)
    deallocate(buffobj%buffer_id)
    deallocate(buffobj%interp_method)
    deallocate(buffobj%tile_count)
    deallocate(buffobj%area)
    deallocate(buffobj%volume)
end subroutine

!! -----------Type-specific routines for buffer0-5d 

!! allocations could be done in one routine for 0-5d if buffobj is changed to fmsDiagBuffer_class
!! not sure which approach would be better

!> allocates a scalar buffer to given mold type
!> @returns buffer id for 
integer function allocate_buffer_0d(buffobj, mold, init_time) &
result(rslt)
    class(buffer0d), intent(inout), target :: buffobj !< scalar buffer object
    class(*),intent(in) :: mold !< allocates to the type of mold
    type(time_type), intent(in), optional :: init_time
    type(fmsDiagBufferContainer_type), allocatable :: buff_type

    select type (mold)
        type is (integer(kind=i4_kind))
            allocate(integer(kind=i4_kind) :: buffobj%buffer)
            allocate(integer(kind=i4_kind) :: buffobj%missing_value)
            allocate(integer(kind=i4_kind) :: buffobj%data_RANGE(2))
        type is (integer(kind=i8_kind))
            allocate(integer(kind=i8_kind) :: buffobj%buffer)
            allocate(integer(kind=i8_kind) :: buffobj%missing_value)
            allocate(integer(kind=i8_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r4_kind))
            allocate(real(kind=r4_kind) :: buffobj%buffer)
            allocate(real(kind=r4_kind) :: buffobj%missing_value)
            allocate(real(kind=r4_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r8_kind))
            allocate(real(kind=r8_kind) :: buffobj%buffer)
            allocate(real(kind=r8_kind) :: buffobj%missing_value)
            allocate(real(kind=r8_kind) :: buffobj%data_RANGE(2))
        class default
             call mpp_error("allocate_buffer_0d", &
                     "The mold value passed to allocate a buffer is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
    allocate(buffobj%buffer_id)
    allocate(buffobj%interp_method)
    allocate(buffobj%tile_count)
    allocate(buffobj%area)
    allocate(buffobj%volume)

    num_buffers = num_buffers + 1
    rslt = num_buffers
    buffobj%buffer_id = num_buffers

end function

!> allocates a 1D buffer to given mold type
integer function allocate_buffer_1d(buffobj, mold, size) &
result(rslt)
    class(buffer1d), intent(inout), target :: buffobj !< scalar buffer object
    class(*),intent(in) :: mold !< allocates to the type of mold
    integer, intent(in) :: size !< dimension bounds
    type(fmsDiagBufferContainer_type), allocatable :: buff_type

    !if(.not. allocated(buffer_list_1d) .or. num_buffers(1)+1 .gt. BUFFER_LIST_SIZE) call allocate_buffer_list(1)
    !num_buffers(1) = num_buffers(1) + 1

    select type (mold)
        type is (integer(kind=i4_kind))
            allocate(integer(kind=i4_kind) :: buffobj%buffer(size))
            allocate(integer(kind=i4_kind) :: buffobj%missing_value)
            allocate(integer(kind=i4_kind) :: buffobj%data_RANGE(2))
        type is (integer(kind=i8_kind))
            allocate(integer(kind=i8_kind) :: buffobj%buffer(size))
            allocate(integer(kind=i8_kind) :: buffobj%missing_value)
            allocate(integer(kind=i8_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r4_kind))
            allocate(real(kind=r4_kind) :: buffobj%buffer(size))
            allocate(real(kind=r4_kind) :: buffobj%missing_value)
            allocate(real(kind=r4_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r8_kind))
            allocate(real(kind=r8_kind) :: buffobj%buffer(size))
            allocate(real(kind=r8_kind) :: buffobj%missing_value)
            allocate(real(kind=r8_kind) :: buffobj%data_RANGE(2))
        class default
             call mpp_error("allocate_buffer_1d", &
                     "The mold value passed to allocate a buffer is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
    allocate(buffobj%buffer_id)
    allocate(buffobj%interp_method)
    allocate(buffobj%tile_count)
    allocate(buffobj%area)
    allocate(buffobj%volume)
    num_buffers = num_buffers + 1
    rslt = num_buffers
    buffobj%buffer_id = num_buffers

end function
!> allocates a 2D buffer to given mold type
integer function allocate_buffer_2d(buffobj, mold, sizes) &
result(rslt)
    class(buffer2d), intent(inout), target :: buffobj !< 2D buffer object
    class(*),intent(in) :: mold !< allocates to the type of mold
    integer, intent(in) :: sizes(2) !< dimension bounds
    type(fmsDiagBufferContainer_type), allocatable :: buff_type

    !if(.not. allocated(buffer_list_2d) .or. num_buffers(2)+1 .gt. BUFFER_LIST_SIZE) call allocate_buffer_list(2)
    !num_buffers(2) = num_buffers(2) + 1
    select type (mold)
        type is (integer(kind=i4_kind))
            allocate(integer(kind=i4_kind) :: buffobj%buffer(sizes(1), sizes(2)))
            allocate(integer(kind=i4_kind) :: buffobj%missing_value)
            allocate(integer(kind=i4_kind) :: buffobj%data_RANGE(2))
        type is (integer(kind=i8_kind))
            allocate(integer(kind=i8_kind) :: buffobj%buffer(sizes(1), sizes(2)))
            allocate(integer(kind=i8_kind) :: buffobj%missing_value)
            allocate(integer(kind=i8_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r4_kind))
            allocate(real(kind=r4_kind) :: buffobj%buffer(sizes(1), sizes(2)))
            allocate(real(kind=r4_kind) :: buffobj%missing_value)
            allocate(real(kind=r4_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r8_kind))
            allocate(real(kind=r8_kind) :: buffobj%buffer(sizes(1), sizes(2)))
            allocate(real(kind=r8_kind) :: buffobj%missing_value)
            allocate(real(kind=r8_kind) :: buffobj%data_RANGE(2))
        class default
             call mpp_error("allocate_buffer_1d", &
                     "The mold value passed to allocate a buffer is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
    allocate(buffobj%buffer_id)
    allocate(buffobj%interp_method)
    allocate(buffobj%tile_count)
    allocate(buffobj%area)
    allocate(buffobj%volume)

    num_buffers = num_buffers + 1
    rslt = num_buffers
    buffobj%buffer_id = num_buffers
end function

!> allocates a 3D buffer to given mold type
integer function allocate_buffer_3d(buffobj, mold, sizes) &
result(rslt)
    class(buffer3d), intent(inout), target :: buffobj !< 3D buffer object
    class(*),intent(in) :: mold !< allocates to the type of mold
    integer, intent(in) :: sizes(3) !< dimension sizes
    type(fmsDiagBufferContainer_type), allocatable :: buff_type

    !if(.not. allocated(buffer_list_3d) .or. num_buffers(3)+1 .gt. BUFFER_LIST_SIZE) call allocate_buffer_list(3)
    !num_buffers(3) = num_buffers(3) + 1
    select type (mold)
        type is (integer(kind=i4_kind))
            allocate(integer(kind=i4_kind) :: buffobj%buffer( sizes(1),sizes(2), sizes(3)))
            allocate(integer(kind=i4_kind) :: buffobj%missing_value)
            allocate(integer(kind=i4_kind) :: buffobj%data_RANGE(2))
        type is (integer(kind=i8_kind))
            allocate(integer(kind=i8_kind) :: buffobj%buffer( sizes(1),sizes(2), sizes(3)))
            allocate(integer(kind=i8_kind) :: buffobj%missing_value)
            allocate(integer(kind=i8_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r4_kind))
            allocate(real(kind=r4_kind) :: buffobj%buffer( sizes(1),sizes(2), sizes(3)))
            allocate(real(kind=r4_kind) :: buffobj%missing_value)
            allocate(real(kind=r4_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r8_kind))
            allocate(real(kind=r8_kind) :: buffobj%buffer( sizes(1),sizes(2), sizes(3)))
            allocate(real(kind=r8_kind) :: buffobj%missing_value)
            allocate(real(kind=r8_kind) :: buffobj%data_RANGE(2))
        class default
             call mpp_error("allocate_buffer_3d", &
                     "The mold value passed to allocate a buffer is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
    allocate(buffobj%buffer_id)
    allocate(buffobj%interp_method)
    allocate(buffobj%tile_count)
    allocate(buffobj%area)
    allocate(buffobj%volume)
    num_buffers = num_buffers + 1
    rslt = num_buffers
    buffobj%buffer_id = num_buffers
end function

!> allocates a 4D buffer to given mold type
integer function allocate_buffer_4d(buffobj, mold, sizes) &
result(rslt)
    class(buffer4d), intent(inout), target :: buffobj !< 4D buffer object
    class(*),intent(in) :: mold !< allocates to the type of mold
    integer, intent(in) :: sizes(4) !< dimension sizes
    type(fmsDiagBufferContainer_type), allocatable :: buff_type

    !if(.not. allocated(buffer_list_4d) .or. num_buffers(4)+1 .gt. BUFFER_LIST_SIZE) call allocate_buffer_list(4)
    !num_buffers(4) = num_buffers(4) + 1
    select type (mold)
        type is (integer(kind=i4_kind))
            allocate(integer(kind=i4_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4)))
            allocate(integer(kind=i4_kind) :: buffobj%missing_value)
            allocate(integer(kind=i4_kind) :: buffobj%data_RANGE(2))
        type is (integer(kind=i8_kind))
            allocate(integer(kind=i8_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4)))
            allocate(integer(kind=i8_kind) :: buffobj%missing_value)
            allocate(integer(kind=i8_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r4_kind))
            allocate(real(kind=r4_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4)))
            allocate(real(kind=r4_kind) :: buffobj%missing_value)
            allocate(real(kind=r4_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r8_kind))
            allocate(real(kind=r8_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4)))
            allocate(real(kind=r8_kind) :: buffobj%missing_value)
            allocate(real(kind=r8_kind) :: buffobj%data_RANGE(2))
        class default
             call mpp_error("allocate_buffer_4d", &
                     "The mold value passed to allocate a buffer is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
    allocate(buffobj%buffer_id)
    allocate(buffobj%interp_method)
    allocate(buffobj%tile_count)
    allocate(buffobj%area)
    allocate(buffobj%volume)
    num_buffers = num_buffers + 1
    rslt = num_buffers
    buffobj%buffer_id = num_buffers

end function

!> allocates a 5D buffer to given mold type
integer function allocate_buffer_5d(buffobj, mold, sizes) &
result(rslt)
    class(buffer5d), intent(inout), target :: buffobj !< 5D buffer object
    class(*),intent(in) :: mold !< allocates to the type of mold
    integer, intent(in) :: sizes(5) !< dimension sizes
    !type(fmsDiagBufferContainer_type), allocatable :: buff_type

    !if(.not. allocated(buffer_list_5d) .or. num_buffers(5)+1 .gt. BUFFER_LIST_SIZE) call allocate_buffer_list(5)
    !num_buffers(5) = num_buffers(5) + 1
    select type (mold)
        type is (integer(kind=i4_kind))
            allocate(integer(kind=i4_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
            allocate(integer(kind=i4_kind) :: buffobj%missing_value)
            allocate(integer(kind=i4_kind) :: buffobj%data_RANGE(2))
        type is (integer(kind=i8_kind))
            allocate(integer(kind=i8_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
            allocate(integer(kind=i8_kind) :: buffobj%missing_value)
            allocate(integer(kind=i8_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r4_kind))
            allocate(real(kind=r4_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
            allocate(real(kind=r4_kind) :: buffobj%missing_value)
            allocate(real(kind=r4_kind) :: buffobj%data_RANGE(2))
        type is (real(kind=r8_kind))
            allocate(real(kind=r8_kind) :: buffobj%buffer(sizes(1),sizes(2),sizes(3),sizes(4),sizes(5)))
            allocate(real(kind=r8_kind) :: buffobj%missing_value)
            allocate(real(kind=r8_kind) :: buffobj%data_RANGE(2))
        class default
             call mpp_error("allocate_buffer_5d", &
                     "The mold value passed to allocate a buffer is not a r8, r4, i8, or i4",&
                     FATAL)
    end select
    allocate(buffobj%buffer_id)
    allocate(buffobj%interp_method)
    allocate(buffobj%tile_count)
    allocate(buffobj%area)
    allocate(buffobj%volume)
    num_buffers = num_buffers + 1
    rslt = num_buffers
    buffobj%buffer_id = num_buffers
end function

!> @brief Gets buffer data from buffer0d type
!! @return copy of the buffer data
function get_buffer_0d (obj) &
result(rslt)
    class (buffer0d), intent(in) :: obj
    class(*), allocatable :: rslt
    if (allocated(obj%buffer)) then
        rslt = obj%buffer
    else
        call mpp_error(FATAL, 'get_buffer_0d: buffer not allocated')
    endif
end function
!> @brief Gets buffer data from buffer1d type
!! @return copy of the buffer data
function get_buffer_1d (obj) &
result(rslt)
    class (buffer1d), intent(in) :: obj
    class(*), allocatable :: rslt(:)
    if (allocated(obj%buffer)) then
        rslt = obj%buffer
    else
        call mpp_error(FATAL, 'get_buffer_1d: buffer not allocated')
    endif
end function
!> @brief Gets buffer data from buffer2d type
!! @return copy of the buffer data
function get_buffer_2d (obj) &
result(rslt)
    class (buffer2d), intent(in) :: obj
    class(*), allocatable :: rslt(:,:)
    if (allocated(obj%buffer)) then
        rslt = obj%buffer
    else
        call mpp_error(FATAL, 'get_buffer_2d: buffer not allocated')
    endif
end function
!> @brief Gets buffer data from buffer3d type
!! @return copy of the buffer data
function get_buffer_3d (obj) &
result(rslt)
    class (buffer3d), intent(in) :: obj
    class(*), allocatable :: rslt(:,:,:)
    if (allocated(obj%buffer)) then
        rslt = obj%buffer
    else
        call mpp_error(FATAL, 'get_buffer_3d: buffer not allocated')
    endif
end function
!> @brief Gets buffer data from buffer4d type
!! @return copy of the buffer data
function get_buffer_4d (obj) &
result(rslt)
    class (buffer4d), intent(in) :: obj
    class(*), allocatable :: rslt(:,:,:,:)
    if (allocated(obj%buffer)) then
        rslt = obj%buffer
    else
        call mpp_error(FATAL, 'get_buffer_4d: buffer not allocated')
    endif
end function
!> @brief Gets buffer data from buffer5d type
!! @return copy of the buffer data
function get_buffer_5d (obj) &
result(rslt)
    class (buffer5d), intent(in) :: obj
    class(*), allocatable :: rslt(:,:,:,:,:)
    if (allocated(obj%buffer)) then
        rslt = obj%buffer
    else
        call mpp_error(FATAL, 'get_buffer_5d: buffer not allocated')
    endif
end function

!> @brief Initializes a buffer to a given fill value
!! TODO default to 0?
subroutine initialize_buffer_0d (buffobj, fillval)
    class(buffer0d), intent(inout) :: buffobj !< any buffer object
    class(*) :: fillval

    if(.not. allocated(buffobj%buffer)) call mpp_error(FATAL, 'initialize_buffer_0d:' // &
            'buffer not yet allocated, allocate_buffer() must be called on this object first.')
    select type(buff => buffobj%buffer)
    type is(real(r8_kind))
        select type(fillval)
        type is(real(r8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
        end select
    type is(real(r4_kind))
        select type(fillval)
        type is(real(r4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
        end select
    type is(integer(i8_kind))
        select type(fillval)
        type is(integer(i8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
        end select
    type is(integer(i4_kind))
        select type(fillval)
        type is(integer(i4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_0d: mismatch between fill and buffer values')
        end select
    class default
        call mpp_error(FATAL, 'initialize buffer_0d: buffer allocated to invalid data type, this shouldnt happen')
    end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_1d (buffobj, fillval)
    class(buffer1d), intent(inout) :: buffobj !< any 1d buffer object
    class(*) :: fillval !< fill value, must be same type as the allocated buffer in buffobj

    if(.not. allocated(buffobj%buffer)) call mpp_error(FATAL, 'initialize_buffer_1d:' // &
            'buffer not yet allocated, allocate_buffer() must be called on this object first.')
    ! have to check fill value and buffer types match
    select type(buff => buffobj%buffer)
    type is(real(r8_kind))
        select type(fillval)
        type is(real(r8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
        end select
    type is(real(r4_kind))
        select type(fillval)
        type is(real(r4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
        end select
    type is(integer(i8_kind))
        select type(fillval)
        type is(integer(i8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
        end select
    type is(integer(i4_kind))
        select type(fillval)
        type is(integer(i4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_1d: mismatch between fill and buffer values')
        end select
    class default
        call mpp_error(FATAL, 'initialize buffer_1d: buffer allocated to invalid data type, this shouldnt happen')
    end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_2d (buffobj, fillval)
    class(buffer2d), intent(inout) :: buffobj !< any buffer object
    class(*) :: fillval

    if(.not. allocated(buffobj%buffer)) call mpp_error(FATAL, 'initialize_buffer_2d:' // &
            'buffer not yet allocated, allocate_buffer() must be called on this object first.')
    ! have to check fill value and buffer types match
    select type(buff => buffobj%buffer)
    type is(real(r8_kind))
        select type(fillval)
        type is(real(r8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
        end select
    type is(real(r4_kind))
        select type(fillval)
        type is(real(r4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
        end select
    type is(integer(i8_kind))
        select type(fillval)
        type is(integer(i8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
        end select
    type is(integer(i4_kind))
        select type(fillval)
        type is(integer(i4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_2d: mismatch between fill and buffer values')
        end select
    class default
        call mpp_error(FATAL, 'initialize buffer_2d: buffer allocated to invalid data type, this shouldnt happen')
    end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_3d (buffobj, fillval)
    class(buffer3d), intent(inout) :: buffobj !< any buffer object
    class(*) :: fillval

    if(.not. allocated(buffobj%buffer)) call mpp_error(FATAL, 'initialize_buffer_3d:' // &
            'buffer not yet allocated, allocate_buffer() must be called on this object first.')
    ! have to check fill value and buffer types match
    select type(buff => buffobj%buffer)
    type is(real(r8_kind))
        select type(fillval)
        type is(real(r8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
        end select
    type is(real(r4_kind))
        select type(fillval)
        type is(real(r4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
        end select
    type is(integer(i8_kind))
        select type(fillval)
        type is(integer(i8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
        end select
    type is(integer(i4_kind))
        select type(fillval)
        type is(integer(i4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_3d: mismatch between fill and buffer values')
        end select
    class default
        call mpp_error(FATAL, 'initialize buffer_3d: buffer allocated to invalid data type, this shouldnt happen')
    end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_4d (buffobj, fillval)
    class(buffer4d), intent(inout) :: buffobj !< any buffer object
    class(*) :: fillval

    if(.not. allocated(buffobj%buffer)) call mpp_error(FATAL, 'initialize_buffer_4d:' // &
            'buffer not yet allocated, allocate_buffer() must be called on this object first.')
    ! have to check fill value and buffer types match
    select type(buff => buffobj%buffer)
    type is(real(r8_kind))
        select type(fillval)
        type is(real(r8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
        end select
    type is(real(r4_kind))
        select type(fillval)
        type is(real(r4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
        end select
    type is(integer(i8_kind))
        select type(fillval)
        type is(integer(i8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
        end select
    type is(integer(i4_kind))
        select type(fillval)
        type is(integer(i4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_4d: mismatch between fill and buffer values')
        end select
    class default
        call mpp_error(FATAL, 'initialize buffer_4d: buffer allocated to invalid data type, this shouldnt happen')
    end select

end subroutine

!> @brief Initializes a buffer to a given fill value
subroutine initialize_buffer_5d (buffobj, fillval)
    class(buffer5d), intent(inout) :: buffobj !< any buffer object
    class(*) :: fillval

    if(.not. allocated(buffobj%buffer)) call mpp_error(FATAL, 'initialize_buffer_5d:' // &
            'buffer not yet allocated, allocate_buffer() must be called on this object first.')
    ! have to check fill value and buffer types match
    select type(buff => buffobj%buffer)
    type is(real(r8_kind))
        select type(fillval)
        type is(real(r8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
        end select
    type is(real(r4_kind))
        select type(fillval)
        type is(real(r4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
        end select
    type is(integer(i8_kind))
        select type(fillval)
        type is(integer(i8_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
        end select
    type is(integer(i4_kind))
        select type(fillval)
        type is(integer(i4_kind))
            buff = fillval
        class default
            call mpp_error(FATAL, 'initialize_buffer_5d: mismatch between fill and buffer values')
        end select
    class default
        call mpp_error(FATAL, 'initialize buffer_5d: buffer allocated to invalid data type, this shouldnt happen')
    end select

end subroutine

!> @brief Add values to 0d buffer
!! this will just call the init routine since there's only one value
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_0d(buffobj, input_data)
    class(buffer0d), intent(inout) :: buffobj
    class(*)            :: input_data
    if( .not. allocated(buffobj%buffer)) call mpp_error (FATAL, 'add_to_buffer_1d: buffer not yet allocated')
    call buffobj%initialize_buffer(input_data)
end subroutine

!> @brief Add values to 1d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_1d(buffobj, input_data)
    class(buffer1d), intent(inout) :: buffobj
    class(*)                       :: input_data(:)
    integer                        :: n !< number of elements in input data
    logical                        :: type_error = .false.
    if( .not. allocated(buffobj%buffer)) call mpp_error (FATAL, 'add_to_buffer_1d: buffer not yet allocated')
    n = size(input_data)
    if( n .gt. size(buffobj%buffer)) call mpp_error( FATAL,"add_to_buffer_1d: input data larger than allocated buffer")
    ! have to check both types for assignment
    select type( buffer => buffobj%buffer )
    type is(integer(i4_kind))
        select type(input_data)
        type is(integer(i4_kind))
            buffer(1:n) = input_data(1:n)
        class default
            type_error = .true.
        end select
    type is(integer(i8_kind))
        select type(input_data)
        type is(integer(i8_kind))
            buffer(1:n) = input_data(1:n)
        class default
            type_error = .true.
        end select
    type is(real(r4_kind))
        select type(input_data)
        type is(real(r4_kind))
            buffer(1:n) = input_data(1:n)
        class default
            type_error = .true.
        end select
    type is(real(r8_kind))
        select type(input_data)
        type is(real(r8_kind))
            buffer(1:n) = input_data(1:n)
        class default
            type_error = .true.
        end select
    end select
    if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 2d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_2d(buffobj, input_data)
    class(buffer2d), intent(inout) :: buffobj
    class(*)                       :: input_data(:,:)
    integer                        :: n1, n2 !< number of elements per dimension
    logical                        :: type_error
    if( .not. allocated(buffobj%buffer)) call mpp_error (FATAL, 'add_to_buffer_2d: buffer not yet allocated')
    n1 = size(input_data, 1)
    n2 = size(input_data, 2)
    if( n1 .gt. size(buffobj%buffer, 1) .or. n2 .gt. size(buffobj%buffer, 2)) then
        call mpp_error( FATAL,"add_to_buffer_2d: input data larger than allocated buffer")
    endif
    !buffobj%buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
    ! have to check both types for assignment
    select type( buffer => buffobj%buffer )
    type is(integer(i4_kind))
        select type(input_data)
        type is(integer(i4_kind))
            buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
        class default
            type_error = .true.
        end select
    type is(integer(i8_kind))
        select type(input_data)
        type is(integer(i8_kind))
            buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
        class default
            type_error = .true.
        end select
    type is(real(r4_kind))
        select type(input_data)
        type is(real(r4_kind))
            buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
        class default
            type_error = .true.
        end select
    type is(real(r8_kind))
        select type(input_data)
        type is(real(r8_kind))
            buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
        class default
            type_error = .true.
        end select
    end select
    if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 3d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_3d(buffobj, input_data)
    class(buffer3d), intent(inout) :: buffobj
    class(*)                       :: input_data(:,:,:)
    integer                        :: n1, n2, n3 !< number of elements per dimension
    logical                        :: type_error
    if( .not. allocated(buffobj%buffer)) call mpp_error (FATAL, 'add_to_buffer_3d: buffer not yet allocated')
    n1 = size(input_data, 1)
    n2 = size(input_data, 2)
    n3 = size(input_data, 3)
    if( n1 .gt. size(buffobj%buffer, 1) .or. n2 .gt. size(buffobj%buffer, 2) .or. &
        n3 .gt. size(buffobj%buffer, 3)) then
        call mpp_error( FATAL,"add_to_buffer_2d: input data larger than allocated buffer")
    endif
    !buffobj%buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
    ! have to check both types for assignment
    select type( buffer => buffobj%buffer )
    type is(integer(i4_kind))
        select type(input_data)
        type is(integer(i4_kind))
            buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
        class default
            type_error = .true.
        end select
    type is(integer(i8_kind))
        select type(input_data)
        type is(integer(i8_kind))
            buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
        class default
            type_error = .true.
        end select
    type is(real(r4_kind))
        select type(input_data)
        type is(real(r4_kind))
            buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
        class default
            type_error = .true.
        end select
    type is(real(r8_kind))
        select type(input_data)
        type is(real(r8_kind))
            buffer(1:n1, 1:n2, 1:n3) = input_data(1:n1, 1:n2, 1:n3)
        class default
            type_error = .true.
        end select
    end select
    if( type_error ) call mpp_error (FATAL,'add_to_buffer_1d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 4d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_4d(buffobj, input_data)
    class(buffer4d), intent(inout) :: buffobj
    class(*)                       :: input_data(:,:,:,:)
    integer                        :: n1, n2, n3, n4!< number of elements per dimension
    logical                        :: type_error
    if( .not. allocated(buffobj%buffer)) call mpp_error (FATAL, 'add_to_buffer_4d: buffer not yet allocated')
    n1 = size(input_data, 1)
    n2 = size(input_data, 2)
    n3 = size(input_data, 3)
    n4 = size(input_data, 4)
    if( n1 .gt. size(buffobj%buffer, 1) .or. n2 .gt. size(buffobj%buffer, 2) .or. &
        n3 .gt. size(buffobj%buffer, 3) .or. n4 .gt. size(buffobj%buffer, 4)) then
        call mpp_error( FATAL,"add_to_buffer_4d: input data larger than allocated buffer")
    endif
    ! have to check both types for assignment
    select type( buffer => buffobj%buffer )
    type is(integer(i4_kind))
        select type(input_data)
        type is(integer(i4_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
        class default
            type_error = .true.
        end select
    type is(integer(i8_kind))
        select type(input_data)
        type is(integer(i8_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
        class default
            type_error = .true.
        end select
    type is(real(r4_kind))
        select type(input_data)
        type is(real(r4_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
        class default
            type_error = .true.
        end select
    type is(real(r8_kind))
        select type(input_data)
        type is(real(r8_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4) = input_data(1:n1, 1:n2, 1:n3, 1:n4)
        class default
            type_error = .true.
        end select
    end select
    if( type_error ) call mpp_error (FATAL,'add_to_buffer_4d: mismatch between allocated buffer and input data types')
end subroutine

!> @brief Add values to 5d buffer
!! input_data must match allocated type of buffer object
subroutine add_to_buffer_5d(buffobj, input_data)
    class(buffer5d), intent(inout) :: buffobj
    class(*)                       :: input_data(:,:,:,:,:)
    integer                        :: n1, n2, n3, n4, n5 !< number of elements per dimension
    logical                        :: type_error
    if( .not. allocated(buffobj%buffer)) call mpp_error (FATAL, 'add_to_buffer_5d: buffer not yet allocated')
    n1 = size(input_data, 1)
    n2 = size(input_data, 2)
    n3 = size(input_data, 3)
    n4 = size(input_data, 4)
    n5 = size(input_data, 5)
    if( n1 .gt. size(buffobj%buffer, 1) .or. n2 .gt. size(buffobj%buffer, 2) .or. &
        n3 .gt. size(buffobj%buffer, 3) .or. n4 .gt. size(buffobj%buffer, 4) .or. &
        n5 .gt. size(buffobj%buffer, 5)) then
        call mpp_error( FATAL,"add_to_buffer_4d: input data larger than allocated buffer")
    endif
    !buffobj%buffer(1:n1, 1:n2) = input_data(1:n1, 1:n2)
    ! have to check both types for assignment
    select type( buffer => buffobj%buffer )
    type is(integer(i4_kind))
        select type(input_data)
        type is(integer(i4_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
        class default
            type_error = .true.
        end select
    type is(integer(i8_kind))
        select type(input_data)
        type is(integer(i8_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
        class default
            type_error = .true.
        end select
    type is(real(r4_kind))
        select type(input_data)
        type is(real(r4_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
        class default
            type_error = .true.
        end select
    type is(real(r8_kind))
        select type(input_data)
        type is(real(r8_kind))
            buffer(1:n1, 1:n2, 1:n3, 1:n4, 1:n5) = input_data(1:n1, 1:n2, 1:n3, 1:n4, 1:n5)
        class default
            type_error = .true.
        end select
    end select
    if( type_error ) call mpp_error (FATAL,'add_to_buffer_5d: mismatch between allocated buffer and input data types')
end subroutine
end module fms_diag_buffer_mod
