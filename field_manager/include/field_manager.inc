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
!> @returns The number of values that have been decoded. This allows
!! a user to define a large array and fill it partially with
!! values from a list. This should be the size of the value array.
function PARSE_REALS_ ( text, label, values ) result (parse)
character(len=*), intent(in)  :: text !< The text string from which the values will be parsed.
character(len=*), intent(in)  :: label !< A label which describes the values being decoded.
real(FMS_FM_KIND_), intent(out) :: values(:) !< The value or values that have been decoded.

include 'parse.inc'
end function PARSE_REALS_

function PARSE_REAL_ ( text, label, parse_rval ) result (parse)
character(len=*), intent(in)  :: text !< The text string from which the values will be parsed.
character(len=*), intent(in)  :: label !< A label which describes the values being decoded.
real(FMS_FM_KIND_), intent(out) :: parse_rval !< The value or values that have been decoded.
integer :: parse

real(FMS_FM_KIND_) :: values(1)

   parse = PARSE_REALS_( text, label, values )
   if (parse > 0) parse_rval = values(1)
 end function PARSE_REAL_

!> @returns A flag to indicate whether the function operated with (false) or without
!! (true) errors.
function  FM_GET_VALUE_REAL_(name, get_rval, index)                 &
          result (success)
logical                                :: success
character(len=*), intent(in)           :: name  !< The name of a field that the user wishes to get a value for.
real(FMS_FM_KIND_), intent(out)        :: get_rval !< The value associated with the named field
integer,          intent(in), optional :: index !< An optional index to retrieve a single value from an array.

integer                         :: index_t
type (field_def), pointer, save :: temp_field_p
integer                         :: out_unit

integer, parameter :: lkind=FMS_FM_KIND_

out_unit = stdout()
!        Initialize the field manager if needed
if (.not. module_is_initialized) then
  call initialize_module_variables
endif
!        Must supply a field field name
if (name .eq. ' ') then
  get_rval = 0.0_lkind
  success = .false.
  return
endif
!        Set index to retrieve
if (present(index)) then
  index_t = index
else
  index_t = 1
endif
!        Get a pointer to the field
temp_field_p => get_field(name, current_list_p)

if (associated(temp_field_p)) then
!        check that the field is the correct type
  if (temp_field_p%field_type .eq. real_type) then
    if (index_t .lt. 1 .or. index_t .gt. temp_field_p%max_index) then
!        Index is not positive or is too large
       get_rval = 0.0_lkind
      success = .false.
    else
!        extract the value; the value is stored as r8_kind
      get_rval = real(temp_field_p%r_value(index_t),lkind)
      success = .true.
    endif
  else
    get_rval = 0.0_lkind
    success = .false.
  endif
else
  get_rval = 0.0_lkind
  success = .false.
endif

end function  FM_GET_VALUE_REAL_

!> @brief Assigns a given value to a given field
!> @returns An index for the named field
function  FM_NEW_VALUE_REAL_(name, new_rval, create, index, append) &
          result (field_index)
integer                                :: field_index
character(len=*), intent(in)           :: name !< The name of a field that the user wishes to create
                                               !! a value for.
real(FMS_FM_KIND_),    intent(in)      :: new_rval !< The value that the user wishes to apply to the
                                                !! named field.
logical,          intent(in), optional :: create !< If present and .true., then a value for this
                                                 !! field will be created.
integer,          intent(in), optional :: index !< The index to an array of values that the user
                                                !! wishes to apply a new value.
logical,          intent(in), optional :: append !< If present and .true., then append the value to
      !! an array of the present values. If present and .true., then index cannot be greater than 0.

logical                          :: create_t
integer                          :: i
integer                          :: index_t
real(r8_kind), allocatable, dimension(:) :: temp_r_value
character(len=fm_path_name_len)  :: path
character(len=fm_field_name_len) :: base
type (field_def), pointer, save  :: temp_list_p
type (field_def), pointer, save  :: temp_field_p
integer                          :: out_unit

out_unit = stdout()
!        Initialize the field manager if needed
if (.not. module_is_initialized) then
  call initialize_module_variables
endif
!        Must supply a field name
if (name .eq. ' ') then
  field_index = NO_FIELD
  return
endif
!        Check for optional arguments
if (present(create)) then
  create_t = create
else
  create_t = .false.
endif
!        Check that append is not true and index greater than 0
if (present(index) .and. present(append)) then
  if (append .and. index .gt. 0) then
    field_index = NO_FIELD
    return
  endif
endif
!        Set index to define
if (present(index)) then
  index_t = index
  if (index_t .lt. 0) then
!        Index is negative
    field_index = NO_FIELD
    return
  endif
else
  index_t = 1
endif

!        Get a pointer to the parent list
call find_base(name, path, base)
temp_list_p => find_list(path, current_list_p, create_t)

if (associated(temp_list_p)) then
  temp_field_p => find_field(base, temp_list_p)
  if (.not. associated(temp_field_p)) then
!        Create the field if it doesn't exist
    temp_field_p => create_field(temp_list_p, base)
  endif
  if (associated(temp_field_p)) then
!        Check if the field_type is the same as previously
!        If not then reset max_index to 0
    if (temp_field_p%field_type == integer_type) then
       ! promote integer field to real
       allocate(temp_field_p%r_value(size(temp_field_p%i_value)))
       do i = 1, size(temp_field_p%i_value)
          ! all real field values are stored as r8_kind
          temp_field_p%r_value(i) = real(temp_field_p%i_value(i),r8_kind)
       enddo
       temp_field_p%field_type = real_type
       deallocate(temp_field_p%i_value)
    else if (temp_field_p%field_type /= real_type ) then
      ! slm: why reset index to 0? does it make any sense? It sounds like this is the
      ! case where the values in the array have different types, so is it not an error?
      ! Or, alternatively, if string follows a real value, should not be the entire
      ! array converted to string type?
      temp_field_p%max_index = 0
    endif
!        Assign the type
    temp_field_p%field_type = real_type
!        Set the index if appending
    if (present(append)) then
      if (append) then
        index_t = temp_field_p%max_index + 1
      endif
    endif
    if (index_t .gt. temp_field_p%max_index + 1) then
!        Index too large
      field_index = NO_FIELD
      return
    elseif (index_t .eq. 0 .and.                                &
            temp_field_p%max_index .gt. 0) then
!        Can't set non-null field to null
      field_index = NO_FIELD
      return
    elseif (.not. allocated(temp_field_p%r_value) .and.        &
            index_t .gt. 0) then
!        Array undefined, so allocate the array
      allocate(temp_field_p%r_value(1))
      temp_field_p%max_index = 1
      temp_field_p%array_dim = 1
    elseif (index_t .gt. temp_field_p%array_dim) then
!        Array is too small, so allocate new array and copy over
!        old values
      temp_field_p%array_dim = temp_field_p%array_dim + array_increment
      allocate (temp_r_value(temp_field_p%array_dim))
      do i = 1, temp_field_p%max_index
        temp_r_value(i) = temp_field_p%r_value(i)
      enddo
      if (allocated(temp_field_p%r_value)) deallocate(temp_field_p%r_value)
      temp_field_p%r_value = temp_r_value
      temp_field_p%max_index = index_t
    endif
!        Assign the value and set the field_index for return
!        for non-null fields (index_t > 0)
    if (index_t .gt. 0) then
       ! all real field values are stored as r8_kind
      temp_field_p%r_value(index_t) = real(new_rval,r8_kind)
      if (index_t .gt. temp_field_p%max_index) then
        temp_field_p%max_index = index_t
      endif
    endif
    field_index = temp_field_p%index
  else
!        Error in making the field
    field_index = NO_FIELD
  endif
else
!        Error following the path
  field_index = NO_FIELD
endif

end function  FM_NEW_VALUE_REAL_
!> @}
! close documentation grouping
