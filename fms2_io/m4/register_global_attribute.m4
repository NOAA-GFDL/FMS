include(`macros.m4')dnl
!> @brief Add a global attribute.
`subroutine register_global_attribute_'NUM_DIMS`d(fileobj, &'
                                                  attribute_name, &
                                                  attribute_value)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: attribute_name !< Attribute name.
    class(*),dim_declare(NUM_DIMS) intent(in) :: attribute_value !< Attribute value
    integer :: err
    if (fileobj%is_root) then
        call set_netcdf_mode(fileobj%ncid, &
                             define_mode)
        select type(attribute_value)
ifelse(NUM_DIMS,0,
           `type is (character(len=*))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   trim(attribute_value))')
            type is (integer(kind=int32))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (integer(kind=int64))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real32))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            type is (real(kind=real64))
                err = nf90_put_att(fileobj%ncid, &
                                   nf90_global, &
                                   trim(attribute_name), &
                                   attribute_value)
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
`end subroutine register_global_attribute_'NUM_DIMS`d'
