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

!> @brief Add a global attribute.
subroutine register_global_attribute_0d(fileobj, &
                                                  attribute_name, &
                                                  attribute_value, str_len)
    class(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    class(*), intent(in) :: attribute_value !< Attribute value
    integer, intent(in), optional :: str_len !< Length of the string
    integer :: err
    if (fileobj%is_root) then
        call set_netcdf_mode(fileobj%ncid, &
                             define_mode)
        select type(attribute_value)
            type is (character(len=*))
                if (.not. present(str_len)) call error("register_global_attribute_0d: Need to include str length")
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   trim(attribute_value(1:str_len)))
            type is (integer(kind=i4_kind))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (integer(kind=i8_kind))
                if ( .not. fileobj%allow_int8) call error(trim(fileobj%path)//": 64 bit integers are only supported with 'netcdf4' file format"//&
                                                       &". Set netcdf_default_format='netcdf4' in the fms2_io namelist OR "//&
                                                       &"add nc_format='netcdf4' to your open_file call")
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=r4_kind))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=r8_kind))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            class default
                call error("register_global_attribute_0d: unsupported type for "//&
                trim(attribute_name)//" for file: "//trim(fileobj%path)//"")
        end select
        call check_netcdf_code(err, "register_global_attribute_0d: file:"//trim(fileobj%path)//"- attribute:"//trim(attribute_name))
    endif
end subroutine register_global_attribute_0d
!> @brief Add a global attribute.
subroutine register_global_attribute_1d(fileobj, &
                                                  attribute_name, &
                                                  attribute_value, str_len)
    class(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    class(*),dimension(:), intent(in) :: attribute_value !< Attribute value
    integer, intent(in), optional :: str_len !< Length of the string
    integer :: err
    if (fileobj%is_root) then
        call set_netcdf_mode(fileobj%ncid, &
                             define_mode)
        select type(attribute_value)

            type is (integer(kind=i4_kind))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (integer(kind=i8_kind))
                if ( .not. fileobj%allow_int8) call error(trim(fileobj%path)//": 64 bit integers are only supported with 'netcdf4' file format"//&
                                                       &". Set netcdf_default_format='netcdf4' in the fms2_io namelist OR "//&
                                                       &"add nc_format='netcdf4' to your open_file call")
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=r4_kind))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=r8_kind))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            class default
                call error("register_global_attribute_1d: unsupported type for "//&
                trim(attribute_name)//" for file: "//trim(fileobj%path)//"")
        end select
        call check_netcdf_code(err, "register_global_attribute_1d: file:"//trim(fileobj%path)//"- attribute:"//trim(attribute_name))
    endif
end subroutine register_global_attribute_1d
