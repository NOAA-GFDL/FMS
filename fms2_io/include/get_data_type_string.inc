!> @file


!> @brief Return a string describing the input data's type.
subroutine get_data_type_string_0d(sdata, type_string)

  class(*),  intent(in) :: sdata !> Data.
  character(len=*), intent(inout) :: type_string !> Data type.

  select type(sdata)
    type is (integer(kind=int32))
      call string_copy(type_string, "int")
    type is (integer(kind=int64))
      call string_copy(type_string, "int64")
    type is (real(kind=real32))
      call string_copy(type_string, "float")
    type is (real(kind=real64))
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
    type is (integer(kind=int32))
      call string_copy(type_string, "int")
    type is (integer(kind=int64))
      call string_copy(type_string, "int64")
    type is (real(kind=real32))
      call string_copy(type_string, "float")
    type is (real(kind=real64))
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
    type is (integer(kind=int32))
      call string_copy(type_string, "int")
    type is (integer(kind=int64))
      call string_copy(type_string, "int64")
    type is (real(kind=real32))
      call string_copy(type_string, "float")
    type is (real(kind=real64))
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
    type is (integer(kind=int32))
      call string_copy(type_string, "int")
    type is (integer(kind=int64))
      call string_copy(type_string, "int64")
    type is (real(kind=real32))
      call string_copy(type_string, "float")
    type is (real(kind=real64))
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
    type is (integer(kind=int32))
      call string_copy(type_string, "int")
    type is (integer(kind=int64))
      call string_copy(type_string, "int64")
    type is (real(kind=real32))
      call string_copy(type_string, "float")
    type is (real(kind=real64))
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
    type is (integer(kind=int32))
      call string_copy(type_string, "int")
    type is (integer(kind=int64))
      call string_copy(type_string, "int64")
    type is (real(kind=real32))
      call string_copy(type_string, "float")
    type is (real(kind=real64))
      call string_copy(type_string, "double")
    type is (character(len=*))
      call string_copy(type_string, "char")
    class default
      call error("unsupported type.")
    end select
end subroutine get_data_type_string_5d


