include(`macros.m4')dnl
!> @brief Return a string describing the input data's type.
`subroutine get_data_type_string_'NUM_DIMS`d(sdata, &'
                                             type_string)
    class(*),dim_declare(NUM_DIMS) intent(in) :: sdata
    character(len=*),intent(inout) :: type_string
    select type(sdata)
        type is (integer(kind=int32))
            call string_copy(type_string, &
                             "int")
        type is (integer(kind=int64))
            call string_copy(type_string, &
                             "int64")
        type is (real(kind=real32))
            call string_copy(type_string, &
                             "float")
        type is (real(kind=real64))
            call string_copy(type_string, &
                             "double")
        class default
            call error("unsupported type.")
    end select
`end subroutine get_data_type_string_'NUM_DIMS`d'
