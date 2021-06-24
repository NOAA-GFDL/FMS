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

!> @brief Get the value of a variable's attribute.
subroutine get_variable_attribute_0d(fileobj, variable_name, attribute_name, &
                                     attribute_value, broadcast, reproduce_null_char_bug_flag)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  character(len=*), intent(in) :: attribute_name !< Attribute name.
  class(*), intent(inout) :: attribute_value !< Attribute value
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
  logical, intent(in), optional :: reproduce_null_char_bug_flag !< Flag indicating to reproduce the mpp_io bug where
                                                                !! the null characters were not removed after reading
                                                                !! a string attribute

  integer :: varid
  integer :: err
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message
  integer :: j
  character(len=1024), dimension(1) :: charbuf !< 1D Character buffer
  logical :: reproduce_null_char_bug !< Local flag indicating to reproduce the mpp_io bug where
                                     !! the null characters were not removed after reading a string attribute

  !! If the reproduce_null_char_bug flag is turned on fms2io will not remove the null character
  reproduce_null_char_bug = .false.
  if (present(reproduce_null_char_bug_flag)) reproduce_null_char_bug=reproduce_null_char_bug_flag

  append_error_msg = "get_variable_attribute_0d: file:"//trim(fileobj%path)//"- variable:"//&
                     &trim(variable_name)//" attribute: "//trim(attribute_name)

  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(attribute_value)
      type is (character(len=*))
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), charbuf(1))
        if (reproduce_null_char_bug) then
           call string_copy(attribute_value, charbuf(1), check_for_null=.false.)
        else
           call string_copy(attribute_value, charbuf(1), check_for_null=.true.)
        endif
      type is (integer(kind=i4_kind))
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      type is (integer(kind=i8_kind))
        if ( .not. fileobj%allow_int8) call error(trim(fileobj%path)//": 64 bit integers are only supported with 'netcdf4' file format"//&
                                                 &". Set netcdf_default_format='netcdf4' in the fms2_io namelist OR "//&
                                                 &"add nc_format='netcdf4' to your open_file call")
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      type is (real(kind=r4_kind))
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      type is (real(kind=r8_kind))
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      class default
        call error("unsupported attribute type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)
  endif
  if (present(broadcast)) then
    if (.not. broadcast) then
      return
    endif
  endif
  select type(attribute_value)
    type is (character(len=*))
      call mpp_broadcast(charbuf, len(charbuf), fileobj%io_root, pelist=fileobj%pelist)
      if (reproduce_null_char_bug) then
         call string_copy(attribute_value, charbuf(1), check_for_null=.false.)
      else
         call string_copy(attribute_value, charbuf(1), check_for_null=.true.)
      endif
    type is (integer(kind=i4_kind))
      call mpp_broadcast(attribute_value, fileobj%io_root, pelist=fileobj%pelist)
    type is (integer(kind=i8_kind))
      call mpp_broadcast(attribute_value, fileobj%io_root, pelist=fileobj%pelist)
    type is (real(kind=r4_kind))
      call mpp_broadcast(attribute_value, fileobj%io_root, pelist=fileobj%pelist)
    type is (real(kind=r8_kind))
      call mpp_broadcast(attribute_value, fileobj%io_root, pelist=fileobj%pelist)
    class default
      call error("unsupported attribute type: "//trim(append_error_msg))
  end select
end subroutine get_variable_attribute_0d


!> @brief Get the value of a variable's attribute.
subroutine get_variable_attribute_1d(fileobj, variable_name, attribute_name, &
                                     attribute_value, broadcast)

  class(FmsNetcdfFile_t), intent(in) :: fileobj !< File object.
  character(len=*), intent(in) :: variable_name !< Variable name.
  character(len=*), intent(in) :: attribute_name !< Attribute name.
  class(*), dimension(:), intent(inout) :: attribute_value !< Attribute value
  logical, intent(in), optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.

  integer :: varid
  integer :: err
  character(len=200) :: append_error_msg !< Msg to be appended to FATAL error message

  append_error_msg = "get_variable_attribute_1d: file:"//trim(fileobj%path)//"- variable:"//&
                     &trim(variable_name)//" attribute: "//trim(attribute_name)

  if (fileobj%is_root) then
    varid = get_variable_id(fileobj%ncid, trim(variable_name), msg=append_error_msg)
    select type(attribute_value)
      type is (integer(kind=i4_kind))
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      type is (integer(kind=i8_kind))
        if ( .not. fileobj%allow_int8) call error(trim(fileobj%path)//": 64 bit integers are only supported with 'netcdf4' file format"//&
                                                 &". Set netcdf_default_format='netcdf4' in the fms2_io namelist OR "//&
                                                 &"add nc_format='netcdf4' to your open_file call")
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      type is (real(kind=r4_kind))
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      type is (real(kind=r8_kind))
        err = nf90_get_att(fileobj%ncid, varid, trim(attribute_name), attribute_value)
      class default
        call error("unsupported attribute type: "//trim(append_error_msg))
    end select
    call check_netcdf_code(err, append_error_msg)
  endif
  if (present(broadcast)) then
    if (.not. broadcast) then
      return
    endif
  endif
  select type(attribute_value)
    type is (integer(kind=i4_kind))
      call mpp_broadcast(attribute_value, size(attribute_value), fileobj%io_root, &
                         pelist=fileobj%pelist)
    type is (integer(kind=i8_kind))
      call mpp_broadcast(attribute_value, size(attribute_value), fileobj%io_root, &
                         pelist=fileobj%pelist)
    type is (real(kind=r4_kind))
      call mpp_broadcast(attribute_value, size(attribute_value), fileobj%io_root, &
                         pelist=fileobj%pelist)
    type is (real(kind=r8_kind))
      call mpp_broadcast(attribute_value, size(attribute_value), fileobj%io_root, &
                         pelist=fileobj%pelist)
    class default
      call error("unsupported attribute type: "//trim(append_error_msg))
  end select
end subroutine get_variable_attribute_1d
