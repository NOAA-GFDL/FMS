include(`macros.m4')dnl
!> @brief Get the value of a variable's attribute.
`subroutine get_variable_attribute_'NUM_DIMS`d(fileobj, &'
                                               variable_name, &
                                               attribute_name, &
                                               attribute_value, &
                                               broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    class(*),dim_declare(NUM_DIMS) intent(inout) :: attribute_value !< Attribute value
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: varid
    integer :: err
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(attribute_value)
            type is (integer(kind=int32))
                err = nf90_get_att(fileobj%ncid, &
                                   varid, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (integer(kind=int64))
                err = nf90_get_att(fileobj%ncid, &
                                   varid, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real32))
                err = nf90_get_att(fileobj%ncid, &
                                   varid, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real64))
                err = nf90_get_att(fileobj%ncid, &
                                   varid, &
                                   trim(attribute_name), &
                                   attribute_value)
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (present(broadcast)) then
        if (.not. broadcast) then
            return
        endif
    endif
    select type(attribute_value)
        type is (integer(kind=int32))
            call mpp_broadcast(attribute_value, &
                               var_size(NUM_DIMS,attribute_value)
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (integer(kind=int64))
            call mpp_broadcast(attribute_value, &
                               var_size(NUM_DIMS,attribute_value)
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (real(kind=real32))
            call mpp_broadcast(attribute_value, &
                               var_size(NUM_DIMS,attribute_value)
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (real(kind=real64))
            call mpp_broadcast(attribute_value, &
                               var_size(NUM_DIMS,attribute_value)
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        class default
            call error("unsupported type.")
    end select
`end subroutine get_variable_attribute_'NUM_DIMS`d'
