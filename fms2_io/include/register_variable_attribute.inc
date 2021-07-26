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
!> @file
!> @ingroup fms2_io

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
  integer(kind=i4_kind), dimension(2) :: i32_range
  integer(kind=i8_kind), dimension(2) :: i64_range
  real(kind=r4_kind), dimension(2) :: r32_range
  real(kind=r8_kind), dimension(2) :: r64_range
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  append_error_msg = "register_variable_attribute_0d: file:"//trim(fileobj%path)//"- variable:"//&
                     &trim(variable_name)//" attribute: "//trim(attribute_name)

  if (fileobj%is_root) then
    call set_netcdf_mode(fileobj%ncid, define_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(attribute_value)
      type is (character(len=*))
        if (.not. present(str_len)) call error("Need to include str length:"//trim(append_error_msg))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           trim(attribute_value(1:str_len)))
      type is (integer(kind=i4_kind))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (integer(kind=i8_kind))
        if ( .not. fileobj%allow_int8) call error(trim(fileobj%path)//": 64 bit integers are only supported with 'netcdf4' file format"//&
                                                 &". Set netcdf_default_format='netcdf4' in the fms2_io namelist OR "//&
                                                 &"add nc_format='netcdf4' to your open_file call")
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=r4_kind))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=r8_kind))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      class default
        call error("Unsupported attribute type:"//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)

    !The missing_value attribute is not NUG compliant, but is included here to support
    !legacy model component code.
    axtype = get_attribute_type(fileobj%ncid, varid, attribute_name, append_error_msg)
    xtype = get_variable_type(fileobj%ncid, varid, append_error_msg)
    if (string_compare(attribute_name, "_FillValue") .or. &
        string_compare(attribute_name, "valid_min") .or. &
        string_compare(attribute_name, "valid_max") .or. &
        string_compare(attribute_name, "scale_factor") .or. &
        string_compare(attribute_name, "add_offset") .or. &
        string_compare(attribute_name, "missing_value")) then
      if (axtype .ne. xtype) then
        call error("The variable type does not match the attribute type: "//trim(append_error_msg))
      endif
    endif

    if (string_compare(attribute_name, "_FillValue")) then
      if (attribute_exists(fileobj%ncid, varid, "valid_range", msg=append_error_msg)) then
        select type(attribute_value)
          type is (integer(kind=i4_kind))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        i32_range, .false.)
            if (attribute_value .lt. i32_range(1) .or. &
                attribute_value .gt. i32_range(2)) then
              call error("_FillValue inside valid_range: "//trim(append_error_msg))
            endif
          type is (integer(kind=i8_kind))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        i64_range, .false.)
            if (attribute_value .lt. i64_range(1) .or. &
                attribute_value .gt. i64_range(2)) then
              call error("_FillValue inside valid_range: "//trim(append_error_msg))
            endif
          type is (real(kind=r4_kind))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        r32_range, .false.)
            if (attribute_value .lt. r32_range(1) .or. &
                attribute_value .gt. r32_range(2)) then
              call error("_FillValue inside valid_range: "//trim(append_error_msg))
            endif
          type is (real(kind=r8_kind))
            call get_variable_attribute(fileobj, variable_name, "valid_range", &
                                        r64_range, .false.)
            if (attribute_value .lt. r64_range(1) .or. &
                attribute_value .gt. r64_range(2)) then
              call error("_FillValue inside valid_range: "//trim(append_error_msg))
            endif
          class default
            call error("unsupported attribute type: "//trim(append_error_msg))
        end select
      endif
    elseif (string_compare(attribute_name, "valid_min") .or. &
            string_compare(attribute_name, "valid_max")) then
      if (attribute_exists(fileobj%ncid, varid, "valid_range", msg=append_error_msg)) then
        call error("cannot have attributes valid_range and valid_min/max: "//trim(append_error_msg))
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
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  append_error_msg = "register_variable_attribute_0d: file:"//trim(fileobj%path)//"- variable:"//&
                     &trim(variable_name)//" attribute: "//trim(attribute_name)

  if (fileobj%is_root) then
    call set_netcdf_mode(fileobj%ncid, define_mode)
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(attribute_value)
      type is (integer(kind=i4_kind))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (integer(kind=i8_kind))
        if ( .not. fileobj%allow_int8) call error(trim(fileobj%path)//": 64 bit integers are only supported with 'netcdf4' file format"//&
                                                 &". Set netcdf_default_format='netcdf4' in the fms2_io namelist OR "//&
                                                 &"add nc_format='netcdf4' to your open_file call")
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=r4_kind))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      type is (real(kind=r8_kind))
        err = nf90_put_att(fileobj%ncid, varid, trim(attribute_name), &
                           attribute_value)
      class default
        call error("unsupported attribute type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)

    if (string_compare(attribute_name, "_FillValue") .or. &
        string_compare(attribute_name, "valid_min") .or. &
        string_compare(attribute_name, "valid_max") .or. &
        string_compare(attribute_name, "scale_factor") .or. &
        string_compare(attribute_name, "add_offset") .or. &
        string_compare(attribute_name, "missing_value")) then
      call error(trim(attribute_name)//" must be a scalar.")
    elseif (string_compare(attribute_name, "valid_range")) then
      if (attribute_exists(fileobj%ncid, varid, "valid_min", msg=append_error_msg) .or. &
          attribute_exists(fileobj%ncid, varid, "valid_max", msg=append_error_msg)) then
        call error("cannot have valid_range and valid_min/max: "//trim(append_error_msg))
      endif
      axtype = get_attribute_type(fileobj%ncid, varid, attribute_name, append_error_msg)
      xtype = get_variable_type(fileobj%ncid, varid, append_error_msg)
      if (axtype .ne. xtype) then
        call error("The type of the variable is not the same as the type of the variable: "//trim(append_error_msg))
      endif
      if (size(attribute_value) .ne. 2) then
        call error("valid_range must be a vector with two values: "//trim(append_error_msg))
      endif
    endif
  endif
end subroutine register_variable_attribute_1d
