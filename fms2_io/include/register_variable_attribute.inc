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

!> @brief Add an attribute to a variable.
subroutine register_variable_attribute_0d(fileobj, variable_name, attribute_name, &
                                          attribute_value, str_len)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  character(len=*), intent(in) :: attribute_name !< Attribute name.
  class(*), intent(in) :: attribute_value !< Attribute value
  integer, intent(in), optional :: str_len !< Length of the string

  integer :: varid
  integer :: err
  integer :: axtype
  integer :: xtype
  integer(kind=int32), dimension(2) :: i32_range
  integer(kind=int64), dimension(2) :: i64_range
  real(kind=real32), dimension(2) :: r32_range
  real(kind=real64), dimension(2) :: r64_range

  if (fileobj%is_root) then
    call set_netcdf_mode(fileobj%ncid, define_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(attribute_value)
      type is (character(len=*))
        if (.not. present(str_len)) call error("register_variable_attribute_0d: Need to include str length")
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           trim(attribute_value(1:str_len)))
      type is (integer(kind=int32))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (integer(kind=int64))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=real32))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=real64))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)

    !The missing_value attribute is not NUG compliant, but is included here to support
    !legacy model component code.
    axtype = get_attribute_type(fileobj%ncid, varid, attribute_name)
    xtype = get_variable_type(fileobj%ncid, varid)
    if (string_compare(attribute_name, "_FillValue") .or. &
        string_compare(attribute_name, "valid_min") .or. &
        string_compare(attribute_name, "valid_max") .or. &
        string_compare(attribute_name, "scale_factor") .or. &
        string_compare(attribute_name, "add_offset") .or. &
        string_compare(attribute_name, "missing_value")) then
      if (axtype .ne. xtype) then
        call error("type mismatch for "//trim(attribute_name))
      endif
    endif

    if (string_compare(attribute_name, "_FillValue")) then
      if (attribute_exists(fileobj%ncid, varid, "valid_range")) then
        select type(attribute_value)
          type is (integer(kind=int32))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        i32_range, .false.)
            if (attribute_value .lt. i32_range(1) .or. &
                attribute_value .gt. i32_range(2)) then
              call error("_FillValue inside valid_range.")
            endif
          type is (integer(kind=int64))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        i64_range, .false.)
            if (attribute_value .lt. i64_range(1) .or. &
                attribute_value .gt. i64_range(2)) then
              call error("_FillValue inside valid_range.")
            endif
          type is (real(kind=real32))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        r32_range, .false.)
            if (attribute_value .lt. r32_range(1) .or. &
                attribute_value .gt. r32_range(2)) then
              call error("_FillValue inside valid_range.")
            endif
          type is (real(kind=real64))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        r64_range, .false.)
            if (attribute_value .lt. r64_range(1) .or. &
                attribute_value .gt. r64_range(2)) then
              call error("_FillValue inside valid_range.")
            endif
          class default
            call error("unsupported type.")
        end select
      endif
    elseif (string_compare(attribute_name, "valid_min") .or. &
            string_compare(attribute_name, "valid_max")) then
      if (attribute_exists(fileobj%ncid, varid, "valid_range")) then
        call error("cannot have valid_range and valid_min/max.")
      endif
    endif
  endif
end subroutine register_variable_attribute_0d


!> @brief Add an attribute to a variable.
subroutine register_variable_attribute_1d(fileobj, variable_name, attribute_name, &
                                          attribute_value, str_len)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  character(len=*), intent(in) :: attribute_name !< Attribute name.
  class(*), dimension(:), intent(in) :: attribute_value !< Attribute value
  integer, intent(in), optional :: str_len !< Length of the string

  integer :: varid
  integer :: err
  integer :: axtype
  integer :: xtype

  if (fileobj%is_root) then
    call set_netcdf_mode(fileobj%ncid, define_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name))
    select type(attribute_value)
      type is (integer(kind=int32))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (integer(kind=int64))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=real32))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=real64))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      class default
        call error("unsupported type.")
    end select
    call check_netcdf_code(err)

    if (string_compare(attribute_name, "_FillValue") .or. &
        string_compare(attribute_name, "valid_min") .or. &
        string_compare(attribute_name, "valid_max") .or. &
        string_compare(attribute_name, "scale_factor") .or. &
        string_compare(attribute_name, "add_offset") .or. &
        string_compare(attribute_name, "missing_value")) then
      call error(trim(attribute_name)//" must be a scalar.")
    elseif (string_compare(attribute_name, "valid_range")) then
      if (attribute_exists(fileobj%ncid, varid, "valid_min") .or. &
          attribute_exists(fileobj%ncid, varid, "valid_max")) then
        call error("cannot have valid_range and valid_min/max.")
      endif
      axtype = get_attribute_type(fileobj%ncid, varid, attribute_name)
      xtype = get_variable_type(fileobj%ncid, varid)
      if (axtype .ne. xtype) then
        call error("type mismatch for valid_range.")
      endif
      if (size(attribute_value) .ne. 2) then
        call error("valid_range must be a vector with two values.")
      endif
    endif
  endif
end subroutine register_variable_attribute_1d
