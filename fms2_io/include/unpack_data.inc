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
!> @brief Routines for the @ref gather_data_bc interface

!> @addtogroup netcdf_io_mod
!> @{
subroutine unpack_data_0d(fileobj, varid, varname, var_data)
  class(FmsNetcdfFile_t), intent(in)    :: fileobj  !< File object.
  integer,                intent(in)    :: varid    !< Netcdf variable ID
  character(len=*),       intent(in)    :: varname  !< Name of the variable (for error messages)
  class(*),               intent(inout) :: var_data !< Array that the data
                                                    !! will be read into.

  character(len=128)  :: msg         !< Message to append in error message
  real(kind=r4_kind)  :: buf_r4_kind !< r4_kind buffer to read the scale_factor/add_offset to
  real(kind=r8_kind)  :: buf_r8_kind !< r8_kind buffer to read the scale_factor/add_offset to
  integer             :: err         !< netcdf error code

  msg = "Check your read_data call for the variable:"//trim(varname)//&
    " in file:"//trim(fileobj%path)

  if (attribute_exists(fileobj%ncid, varid, "scale_factor") .or. &
    attribute_exists(fileobj%ncid, varid, "add_offset")) then

    select type(var_data)
    type is (real(kind=r4_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r4_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r4_kind

    type is (real(kind=r8_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r8_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r8_kind

    class default
      call error("If using the scale_factor and add_offset variable, the buffer reading the data to needs to be &
        &r4_kind or r8_kind."//trim(msg))
    end select
  end if
end subroutine

subroutine unpack_data_1d(fileobj, varid, varname, var_data)
  class(FmsNetcdfFile_t), intent(in)    :: fileobj     !< File object.
  integer,                intent(in)    :: varid       !< Netcdf variable ID
  character(len=*),       intent(in)    :: varname     !< Name of the variable (for error messages)
  class(*),               intent(inout) :: var_data(:) !< Array that the data
                                                       !! will be read into.

  character(len=128)  :: msg         !< Message to append in error message
  real(kind=r4_kind)  :: buf_r4_kind !< r4_kind buffer to read the scale_factor/add_offset to
  real(kind=r8_kind)  :: buf_r8_kind !< r8_kind buffer to read the scale_factor/add_offset to
  integer             :: err         !< netcdf error code

  msg = "Check your read_data call for the variable:"//trim(varname)//&
    " in file:"//trim(fileobj%path)

  if (attribute_exists(fileobj%ncid, varid, "scale_factor") .or. &
    attribute_exists(fileobj%ncid, varid, "add_offset")) then

    select type(var_data)
    type is (real(kind=r4_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r4_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r4_kind

    type is (real(kind=r8_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r8_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r8_kind

    class default
      call error("If using the scale_factor and add_offset variable, the buffer reading the data to needs to be &
        &r4_kind or r8_kind."//trim(msg))
    end select
  end if
end subroutine

subroutine unpack_data_2d(fileobj, varid, varname, var_data)
  class(FmsNetcdfFile_t), intent(in)    :: fileobj       !< File object.
  integer,                intent(in)    :: varid         !< Netcdf variable ID
  character(len=*),       intent(in)    :: varname       !< Name of the variable (for error messages)
  class(*),               intent(inout) :: var_data(:,:) !< Array that the data
                                                         !! will be read into.

  character(len=128)  :: msg         !< Message to append in error message
  real(kind=r4_kind)  :: buf_r4_kind !< r4_kind buffer to read the scale_factor/add_offset to
  real(kind=r8_kind)  :: buf_r8_kind !< r8_kind buffer to read the scale_factor/add_offset to
  integer             :: err         !< netcdf error code

  msg = "Check your read_data call for the variable:"//trim(varname)//&
    " in file:"//trim(fileobj%path)

  if (attribute_exists(fileobj%ncid, varid, "scale_factor") .or. &
    attribute_exists(fileobj%ncid, varid, "add_offset")) then

    select type(var_data)
    type is (real(kind=r4_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r4_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r4_kind

    type is (real(kind=r8_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r8_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r8_kind

    class default
      call error("If using the scale_factor and add_offset variable, the buffer reading the data to needs to be &
        &r4_kind or r8_kind."//trim(msg))
    end select
  end if
end subroutine

subroutine unpack_data_3d(fileobj, varid, varname, var_data)
  class(FmsNetcdfFile_t), intent(in)    :: fileobj         !< File object.
  integer,                intent(in)    :: varid           !< Netcdf variable ID
  character(len=*),       intent(in)    :: varname         !< Name of the variable (for error messages)
  class(*),               intent(inout) :: var_data(:,:,:) !< Array that the data
                                                           !! will be read into.

  character(len=128)  :: msg         !< Message to append in error message
  real(kind=r4_kind)  :: buf_r4_kind !< r4_kind buffer to read the scale_factor/add_offset to
  real(kind=r8_kind)  :: buf_r8_kind !< r8_kind buffer to read the scale_factor/add_offset to
  integer             :: err         !< netcdf error code

  msg = "Check your read_data call for the variable:"//trim(varname)//&
    " in file:"//trim(fileobj%path)

  if (attribute_exists(fileobj%ncid, varid, "scale_factor") .or. &
    attribute_exists(fileobj%ncid, varid, "add_offset")) then

    select type(var_data)
    type is (real(kind=r4_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r4_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r4_kind

    type is (real(kind=r8_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r8_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r8_kind

    class default
      call error("If using the scale_factor and add_offset variable, the buffer reading the data to needs to be &
        &r4_kind or r8_kind."//trim(msg))
    end select
  end if
end subroutine

subroutine unpack_data_4d(fileobj, varid, varname, var_data)
  class(FmsNetcdfFile_t), intent(in)    :: fileobj           !< File object.
  integer,                intent(in)    :: varid             !< Netcdf variable ID
  character(len=*),       intent(in)    :: varname           !< Name of the variable (for error messages)
  class(*),               intent(inout) :: var_data(:,:,:,:) !< Array that the data
                                                             !! will be read into.

  character(len=128)  :: msg         !< Message to append in error message
  real(kind=r4_kind)  :: buf_r4_kind !< r4_kind buffer to read the scale_factor/add_offset to
  real(kind=r8_kind)  :: buf_r8_kind !< r8_kind buffer to read the scale_factor/add_offset to
  integer             :: err         !< netcdf error code

  msg = "Check your read_data call for the variable:"//trim(varname)//&
    " in file:"//trim(fileobj%path)

  if (attribute_exists(fileobj%ncid, varid, "scale_factor") .or. &
    attribute_exists(fileobj%ncid, varid, "add_offset")) then

    select type(var_data)
    type is (real(kind=r4_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r4_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r4_kind

    type is (real(kind=r8_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r8_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r8_kind

    class default
      call error("If using the scale_factor and add_offset variable, the buffer reading the data to needs to be &
        &r4_kind or r8_kind."//trim(msg))
    end select
  end if
end subroutine

subroutine unpack_data_5d(fileobj, varid, varname, var_data)
  class(FmsNetcdfFile_t), intent(in)    :: fileobj             !< File object.
  integer,                intent(in)    :: varid               !< Netcdf variable ID
  character(len=*),       intent(in)    :: varname             !< Name of the variable (for error messages)
  class(*),               intent(inout) :: var_data(:,:,:,:,:) !< Array that the data
                                                               !! will be read into.

  character(len=128)  :: msg         !< Message to append in error message
  real(kind=r4_kind)  :: buf_r4_kind !< r4_kind buffer to read the scale_factor/add_offset to
  real(kind=r8_kind)  :: buf_r8_kind !< r8_kind buffer to read the scale_factor/add_offset to
  integer             :: err         !< netcdf error code

  msg = "Check your read_data call for the variable:"//trim(varname)//&
    " in file:"//trim(fileobj%path)

  if (attribute_exists(fileobj%ncid, varid, "scale_factor") .or. &
    attribute_exists(fileobj%ncid, varid, "add_offset")) then

    select type(var_data)
    type is (real(kind=r4_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r4_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r4_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r4_kind

    type is (real(kind=r8_kind))
      err = nf90_get_att(fileobj%ncid, varid, "scale_factor", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data * buf_r8_kind

      err = nf90_get_att(fileobj%ncid, varid, "add_offset", buf_r8_kind )
      call check_netcdf_code(err, msg)
      var_data = var_data + buf_r8_kind

    class default
      call error("If using the scale_factor and add_offset variable, the buffer reading the data to needs to be &
        &r4_kind or r8_kind."//trim(msg))
    end select
  end if
end subroutine
