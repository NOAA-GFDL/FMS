!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @file
!> @brief Routines for different dimensions to add global attributes in netcdf files for the
!! @ref register_global_attribute interface

!> @addtogroup netcdf_io_mod
!> @{

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
                if ( .not. fileobj%is_netcdf4) call error(trim(fileobj%path)// &
                                               & ": 64 bit integers are only supported with 'netcdf4' file format"//  &
                                               & ". Set netcdf_default_format='netcdf4' in the fms2_io namelist OR "//&
                                               & "add nc_format='netcdf4' to your open_file call")
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
        call check_netcdf_code(err, "register_global_attribute_0d: file:"//trim(fileobj%path)//"- attribute:"// &
                             & trim(attribute_name))
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
                if ( .not. fileobj%is_netcdf4) call error(trim(fileobj%path)// &
                                                &": 64 bit integers are only supported with 'netcdf4' file format"//&
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
        call check_netcdf_code(err, "register_global_attribute_1d: file:"//trim(fileobj%path)//"- attribute:"// &
                             & trim(attribute_name))
    endif
end subroutine register_global_attribute_1d
!> @}
