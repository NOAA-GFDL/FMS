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
!> @brief Routines to return a given arguments data type as a string.
!! Provides multiple dimensions for @ref get_data_type_string interface.

!> @addtogroup fms_io_utils_mod
!> @{

!> @brief Return a string describing the input data's type.
subroutine get_data_type_string_0d(sdata, type_string)

  class(*),  intent(in) :: sdata !> Data.
  character(len=*), intent(inout) :: type_string !> Data type.

  select type(sdata)
    type is (integer(kind=i4_kind))
      call string_copy(type_string, "int")
    type is (integer(kind=i8_kind))
      call string_copy(type_string, "i8_kind")
    type is (real(kind=r4_kind))
      call string_copy(type_string, "float")
    type is (real(kind=r8_kind))
      call string_copy(type_string, "double")
    type is (character(len=*))
      call string_copy(type_string, "char")
    class default
      call error("unsupported type.")
    end select
end subroutine get_data_type_string_0d



!> @brief Return a string describing the input data's type.
subroutine get_data_type_string_1d(sdata, type_string)

  class(*), dimension(:), intent(in) :: sdata !> Data.
  character(len=*), intent(inout) :: type_string !> Data type.

  select type(sdata)
    type is (integer(kind=i4_kind))
      call string_copy(type_string, "int")
    type is (integer(kind=i8_kind))
      call string_copy(type_string, "i8_kind")
    type is (real(kind=r4_kind))
      call string_copy(type_string, "float")
    type is (real(kind=r8_kind))
      call string_copy(type_string, "double")
    type is (character(len=*))
      call string_copy(type_string, "char")
    class default
      call error("unsupported type.")
    end select
end subroutine get_data_type_string_1d



!> @brief Return a string describing the input data's type.
subroutine get_data_type_string_2d(sdata, type_string)

  class(*), dimension(:,:), intent(in) :: sdata !> Data.
  character(len=*), intent(inout) :: type_string !> Data type.

  select type(sdata)
    type is (integer(kind=i4_kind))
      call string_copy(type_string, "int")
    type is (integer(kind=i8_kind))
      call string_copy(type_string, "i8_kind")
    type is (real(kind=r4_kind))
      call string_copy(type_string, "float")
    type is (real(kind=r8_kind))
      call string_copy(type_string, "double")
    type is (character(len=*))
      call string_copy(type_string, "char")
    class default
      call error("unsupported type.")
    end select
end subroutine get_data_type_string_2d



!> @brief Return a string describing the input data's type.
subroutine get_data_type_string_3d(sdata, type_string)

  class(*), dimension(:,:,:), intent(in) :: sdata !> Data.
  character(len=*), intent(inout) :: type_string !> Data type.

  select type(sdata)
    type is (integer(kind=i4_kind))
      call string_copy(type_string, "int")
    type is (integer(kind=i8_kind))
      call string_copy(type_string, "i8_kind")
    type is (real(kind=r4_kind))
      call string_copy(type_string, "float")
    type is (real(kind=r8_kind))
      call string_copy(type_string, "double")
    type is (character(len=*))
      call string_copy(type_string, "char")
    class default
      call error("unsupported type.")
    end select
end subroutine get_data_type_string_3d



!> @brief Return a string describing the input data's type.
subroutine get_data_type_string_4d(sdata, type_string)

  class(*), dimension(:,:,:,:), intent(in) :: sdata !> Data.
  character(len=*), intent(inout) :: type_string !> Data type.

  select type(sdata)
    type is (integer(kind=i4_kind))
      call string_copy(type_string, "int")
    type is (integer(kind=i8_kind))
      call string_copy(type_string, "i8_kind")
    type is (real(kind=r4_kind))
      call string_copy(type_string, "float")
    type is (real(kind=r8_kind))
      call string_copy(type_string, "double")
    type is (character(len=*))
      call string_copy(type_string, "char")
    class default
      call error("unsupported type.")
    end select
end subroutine get_data_type_string_4d



!> @brief Return a string describing the input data's type.
subroutine get_data_type_string_5d(sdata, type_string)

  class(*), dimension(:,:,:,:,:), intent(in) :: sdata !> Data.
  character(len=*), intent(inout) :: type_string !> Data type.

  select type(sdata)
    type is (integer(kind=i4_kind))
      call string_copy(type_string, "int")
    type is (integer(kind=i8_kind))
      call string_copy(type_string, "i8_kind")
    type is (real(kind=r4_kind))
      call string_copy(type_string, "float")
    type is (real(kind=r8_kind))
      call string_copy(type_string, "double")
    type is (character(len=*))
      call string_copy(type_string, "char")
    class default
      call error("unsupported type.")
    end select
end subroutine get_data_type_string_5d
!> @}
