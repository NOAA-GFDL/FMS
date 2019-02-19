include(`macros.m4')dnl
!> @brief I/O domain root gathers data (variable must not have an unlimited
!!        dimension) from the rest of the ranks and then writes the data to the
!!        netcdf file.  This routine may only be used with variables that
!!        are "compressed".
`subroutine compressed_write_'NUM_DIMS`d(fileobj, &'
                                         variable_name, &
                                         cdata, &
                                         unlim_dim_level)

    !Inputs/outputs
    class(FmsNetcdfCompressedFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(in) :: cdata !< Compressed data that
                                                       !! will be gathered and
                                                       !! written to the 
                                                       !! netcdf file.
    integer,intent(in),optional :: unlim_dim_level !< Level for the unlimited
                                                   !! dimension.

    !Local variables.
ifelse(NUM_DIMS,0,,
   `integer :: compressed_dim_index
    integer,dimension(:),allocatable :: dimension_sizes
    character(len=nf90_max_name),dimension(:),allocatable :: dimension_names
    integer :: cindex
    integer :: i
    integer(kind=int32),dim_declare(NUM_DIMS) allocatable :: buf_int32
    integer(kind=int64),dim_declare(NUM_DIMS) allocatable :: buf_int64
    real(kind=real32),dim_declare(NUM_DIMS) allocatable :: buf_real32
    real(kind=real64),dim_declare(NUM_DIMS) allocatable :: buf_real64
    integer(kind=int32),dim_declare(NUM_DIMS) allocatable :: lbuf_int32
    integer(kind=int64),dim_declare(NUM_DIMS) allocatable :: lbuf_int64
    real(kind=real32),dim_declare(NUM_DIMS) allocatable :: lbuf_real32
    real(kind=real64),dim_declare(NUM_DIMS) allocatable :: lbuf_real64
    integer,dimension(NUM_DIMS) :: e
    integer,dimension(NUM_DIMS) :: c')

    if (.not. is_in_list(fileobj%compressed_vars,variable_name,.false.)) then
        call netcdf_write_data(fileobj, &
                               variable_name, &
                               cdata, &
                               unlim_dim_level=unlim_dim_level)
        return
ifelse(NUM_DIMS,0,
   `else
        call error("this branch should never be reached.")')
    endif
ifelse(NUM_DIMS,0,,
   `!Gather the data onto the I/O root and write it out.
    if (fileobj%is_root) then
        compressed_dim_index = get_variable_compressed_dimension_index(fileobj, &
                                                                       variable_name, &
                                                                       broadcast=.false.)
        call get_variable_size(fileobj, &
                               variable_name, &
                               dimension_sizes, &
                               broadcast=.false.)
        call get_dimension_names(fileobj, &
                                 dimension_names, &
                                 broadcast=.false.)
        do i = 1,size(dimension_names)
            cindex = get_compressed_dimension_index(fileobj,dimension_names(i))
            if (cindex .ne. dimension_not_found) then
                exit
            endif
        enddo
        c(:) = 1
        c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(1)
        e(:) = shape(cdata)
        e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(1)
        select type(cdata)
            type is (integer(kind=int32))
                call allocate_array(buf_int32, &
                                    dimension_sizes)
                call put_array_section(cdata, &
                                       buf_int32, &
                                       c, &
                                       e)
            type is (integer(kind=int64))
                call allocate_array(buf_int64, &
                                    dimension_sizes)
                call put_array_section(cdata, &
                                       buf_int64, &
                                       c, &
                                       e)
            type is (real(kind=real32))
                call allocate_array(buf_real32, &
                                    dimension_sizes)
                call put_array_section(cdata, &
                                       buf_real32, &
                                       c, &
                                       e)
            type is (real(kind=real64))
                call allocate_array(buf_real64, &
                                    dimension_sizes)
                call put_array_section(cdata, &
                                       buf_real64, &
                                       c, &
                                       e)
            class default
                call error("unsupported type.")
        end select
        do i = 2,size(fileobj%pelist)
            c(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_corner(i)
            e(compressed_dim_index) = fileobj%compressed_dims(cindex)%npes_nelems(i)
            select type(cdata)
                type is (integer(kind=int32))
                    call allocate_array(lbuf_int32, &
                                        e)
                    call mpp_recv(lbuf_int32, &
                                  size(lbuf_int32), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_RECV)
                    call put_array_section(lbuf_int32, &
                                           buf_int32, &
                                           c, &
                                           e)
                    deallocate(lbuf_int32)
                type is (integer(kind=int64))
                    call allocate_array(lbuf_int64, &
                                        e)
                    call mpp_recv(lbuf_int64, &
                                  size(lbuf_int64), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_RECV)
                    call put_array_section(lbuf_int64, &
                                           buf_int64, &
                                           c, &
                                           e)
                    deallocate(lbuf_int64)
                type is (real(kind=real32))
                    call allocate_array(lbuf_real32, &
                                        e)
                    call mpp_recv(lbuf_real32, &
                                  size(lbuf_real32), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_RECV)
                    call put_array_section(lbuf_real32, &
                                           buf_real32, &
                                           c, &
                                           e)
                    deallocate(lbuf_real32)
                type is (real(kind=real64))
                    call allocate_array(lbuf_real64, &
                                        e)
                    call mpp_recv(lbuf_real64, &
                                  size(lbuf_real64), &
                                  fileobj%pelist(i))
                    call mpp_sync_self(check=EVENT_RECV)
                    call put_array_section(lbuf_real64, &
                                           buf_real64, &
                                           c, &
                                           e)
                    deallocate(lbuf_real64)
                class default
                    call error("unsupported type.")
            end select
        enddo
        select type(cdata)
            type is (integer(kind=int32))
                call netcdf_write_data(fileobj, &
                                       variable_name, &
                                       buf_int32, &
                                       unlim_dim_level=unlim_dim_level)
                deallocate(buf_int32)
            type is (integer(kind=int64))
                call netcdf_write_data(fileobj, &
                                       variable_name, &
                                       buf_int64, &
                                       unlim_dim_level=unlim_dim_level)
                deallocate(buf_int64)
            type is (real(kind=real32))
                call netcdf_write_data(fileobj, &
                                       variable_name, &
                                       buf_real32, &
                                       unlim_dim_level=unlim_dim_level)
                deallocate(buf_real32)
            type is (real(kind=real64))
                call netcdf_write_data(fileobj, &
                                       variable_name, &
                                       buf_real64, &
                                       unlim_dim_level=unlim_dim_level)
                deallocate(buf_real64)
            class default
                call error("unsupported type.")
        end select
        deallocate(dimension_sizes)
        deallocate(dimension_names)
    else
        select type(cdata)
            type is (integer(kind=int32))
                call mpp_send(cdata, &
                              size(cdata), &
                              fileobj%io_root)
            type is (integer(kind=int64))
                call mpp_send(cdata, &
                              size(cdata), &
                              fileobj%io_root)
            type is (real(kind=real32))
                call mpp_send(cdata, &
                              size(cdata), &
                              fileobj%io_root)
            type is (real(kind=real64))
                call mpp_send(cdata, &
                              size(cdata), &
                              fileobj%io_root)
            class default
                call error("unsupported type.")
        end select
        call mpp_sync_self(check=EVENT_SEND)
    endif')
`end subroutine compressed_write_'NUM_DIMS`d'
