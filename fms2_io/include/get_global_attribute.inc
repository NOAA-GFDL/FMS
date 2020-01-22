!> @brief Get the value of a global attribute.
subroutine get_global_attribute_0d(fileobj, &
                                             attribute_name, &
                                             attribute_value, &
                                             broadcast)
    class(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    class(*), intent(inout) :: attribute_value !< Attribute value
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: err
    if (fileobj%is_root) then
        select type(attribute_value)
            type is (integer(kind=int32))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (integer(kind=int64))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real32))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real64))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
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

                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (integer(kind=int64))
            call mpp_broadcast(attribute_value, &

                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (real(kind=real32))
            call mpp_broadcast(attribute_value, &

                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (real(kind=real64))
            call mpp_broadcast(attribute_value, &

                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        class default
            call error("unsupported type.")
    end select
end subroutine get_global_attribute_0d
!> @brief Get the value of a global attribute.
subroutine get_global_attribute_1d(fileobj, &
                                             attribute_name, &
                                             attribute_value, &
                                             broadcast)
    class(FmsNetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    class(*),dimension(:), intent(inout) :: attribute_value !< Attribute value
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    integer :: err
    if (fileobj%is_root) then
        select type(attribute_value)
            type is (integer(kind=int32))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (integer(kind=int64))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real32))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real64))
                err = nf90_get_att(fileobj%ncid, &
                                   nf90_global, &
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
                               size(attribute_value), &
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (integer(kind=int64))
            call mpp_broadcast(attribute_value, &
                               size(attribute_value), &
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (real(kind=real32))
            call mpp_broadcast(attribute_value, &
                               size(attribute_value), &
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        type is (real(kind=real64))
            call mpp_broadcast(attribute_value, &
                               size(attribute_value), &
                               fileobj%io_root, &
                               pelist=fileobj%pelist)
        class default
            call error("unsupported type.")
    end select
end subroutine get_global_attribute_1d
