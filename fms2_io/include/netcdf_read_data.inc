!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_0d(fileobj, &
                                         variable_name, &
                                         buf, &
                                         unlim_dim_level, &
                                         corner, &

                                         broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.

    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: bcast
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(1) :: c

    if (present(broadcast)) then
        bcast = broadcast
    else
        bcast = .true.
    endif
    c(:) = 1
    if (present(corner)) then
        c(1:0) = corner(:)
    endif

    if (present(unlim_dim_level)) then
        unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                 variable_name, &
                                                                 broadcast=bcast)
        if (unlim_dim_index .ne. 1) then
            call error("unlimited dimension must be the slowest varying" &
                       //" dimension.")
        endif
        c(unlim_dim_index) = unlim_dim_level
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(buf)
            type is (integer(kind=int32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &

                                  )
            type is (integer(kind=int64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &

                                  )
            type is (real(kind=real32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &

                                  )
            type is (real(kind=real64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &

                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (bcast) then
        select type(buf)
            type is (integer(kind=int32))
                call mpp_broadcast(buf, &

                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (integer(kind=int64))
                call mpp_broadcast(buf, &

                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real32))
                call mpp_broadcast(buf, &

                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real64))
                call mpp_broadcast(buf, &

                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            class default
                call error("unsupported type.")
        end select
    endif
end subroutine netcdf_read_data_0d
!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_1d(fileobj, &
                                         variable_name, &
                                         buf, &
                                         unlim_dim_level, &
                                         corner, &
edge_lengths, &
                                         broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: bcast
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(2) :: c
integer,dimension(2) :: e
    if (present(broadcast)) then
        bcast = broadcast
    else
        bcast = .true.
    endif
    c(:) = 1
    if (present(corner)) then
        c(1:1) = corner(:)
    endif
e(:) = 1
    if (present(edge_lengths)) then
        e(1:1) = edge_lengths(:)
    else
        e(1:1) = shape(buf)
    endif
    if (present(unlim_dim_level)) then
        unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                 variable_name, &
                                                                 broadcast=bcast)
        if (unlim_dim_index .ne. 2) then
            call error("unlimited dimension must be the slowest varying" &
                       //" dimension.")
        endif
        c(unlim_dim_index) = unlim_dim_level
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(buf)
            type is (integer(kind=int32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (integer(kind=int64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (bcast) then
        select type(buf)
            type is (integer(kind=int32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (integer(kind=int64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            class default
                call error("unsupported type.")
        end select
    endif
end subroutine netcdf_read_data_1d
!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_2d(fileobj, &
                                         variable_name, &
                                         buf, &
                                         unlim_dim_level, &
                                         corner, &
edge_lengths, &
                                         broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: bcast
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(3) :: c
integer,dimension(3) :: e
    if (present(broadcast)) then
        bcast = broadcast
    else
        bcast = .true.
    endif
    c(:) = 1
    if (present(corner)) then
        c(1:2) = corner(:)
    endif
e(:) = 1
    if (present(edge_lengths)) then
        e(1:2) = edge_lengths(:)
    else
        e(1:2) = shape(buf)
    endif
    if (present(unlim_dim_level)) then
        unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                 variable_name, &
                                                                 broadcast=bcast)
        if (unlim_dim_index .ne. 3) then
            call error("unlimited dimension must be the slowest varying" &
                       //" dimension.")
        endif
        c(unlim_dim_index) = unlim_dim_level
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(buf)
            type is (integer(kind=int32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (integer(kind=int64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (bcast) then
        select type(buf)
            type is (integer(kind=int32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (integer(kind=int64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            class default
                call error("unsupported type.")
        end select
    endif
end subroutine netcdf_read_data_2d
!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_3d(fileobj, &
                                         variable_name, &
                                         buf, &
                                         unlim_dim_level, &
                                         corner, &
edge_lengths, &
                                         broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: bcast
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(4) :: c
integer,dimension(4) :: e
    if (present(broadcast)) then
        bcast = broadcast
    else
        bcast = .true.
    endif
    c(:) = 1
    if (present(corner)) then
        c(1:3) = corner(:)
    endif
e(:) = 1
    if (present(edge_lengths)) then
        e(1:3) = edge_lengths(:)
    else
        e(1:3) = shape(buf)
    endif
    if (present(unlim_dim_level)) then
        unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                 variable_name, &
                                                                 broadcast=bcast)
        if (unlim_dim_index .ne. 4) then
            call error("unlimited dimension must be the slowest varying" &
                       //" dimension.")
        endif
        c(unlim_dim_index) = unlim_dim_level
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(buf)
            type is (integer(kind=int32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (integer(kind=int64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (bcast) then
        select type(buf)
            type is (integer(kind=int32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (integer(kind=int64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            class default
                call error("unsupported type.")
        end select
    endif
end subroutine netcdf_read_data_3d
!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_4d(fileobj, &
                                         variable_name, &
                                         buf, &
                                         unlim_dim_level, &
                                         corner, &
edge_lengths, &
                                         broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: bcast
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(5) :: c
integer,dimension(5) :: e
    if (present(broadcast)) then
        bcast = broadcast
    else
        bcast = .true.
    endif
    c(:) = 1
    if (present(corner)) then
        c(1:4) = corner(:)
    endif
e(:) = 1
    if (present(edge_lengths)) then
        e(1:4) = edge_lengths(:)
    else
        e(1:4) = shape(buf)
    endif
    if (present(unlim_dim_level)) then
        unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                 variable_name, &
                                                                 broadcast=bcast)
        if (unlim_dim_index .ne. 5) then
            call error("unlimited dimension must be the slowest varying" &
                       //" dimension.")
        endif
        c(unlim_dim_index) = unlim_dim_level
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(buf)
            type is (integer(kind=int32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (integer(kind=int64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (bcast) then
        select type(buf)
            type is (integer(kind=int32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (integer(kind=int64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            class default
                call error("unsupported type.")
        end select
    endif
end subroutine netcdf_read_data_4d
!> @brief Read in data from a variable in a netcdf file.
subroutine netcdf_read_data_5d(fileobj, &
                                         variable_name, &
                                         buf, &
                                         unlim_dim_level, &
                                         corner, &
edge_lengths, &
                                         broadcast)
    class(NetcdfFile_t),intent(in) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dimension(:,:,:,:,:), intent(inout) :: buf !< Array that the data
                                                        !! will be read into.
    integer,intent(in),optional :: unlim_dim_level !< Unlimited dimension
                                                   !! level.
    integer,dimension(:),intent(in),optional :: corner !< Array of starting
                                                       !! indices describing
                                                       !! where the data
                                                       !! will be read from.
integer,dimension(:),intent(in),optional :: edge_lengths !< The number of
                                                             !! elements that
                                                             !! will be read
                                                             !! in each dimension.
    logical,intent(in),optional :: broadcast !< Flag controlling whether or
                                             !! not the data will be
                                             !! broadcasted to non
                                             !! "I/O root" ranks.
                                             !! The broadcast will be done
                                             !! by default.
    logical :: bcast
    integer :: err
    integer :: varid
    integer :: unlim_dim_index
    integer,dimension(6) :: c
integer,dimension(6) :: e
    if (present(broadcast)) then
        bcast = broadcast
    else
        bcast = .true.
    endif
    c(:) = 1
    if (present(corner)) then
        c(1:5) = corner(:)
    endif
e(:) = 1
    if (present(edge_lengths)) then
        e(1:5) = edge_lengths(:)
    else
        e(1:5) = shape(buf)
    endif
    if (present(unlim_dim_level)) then
        unlim_dim_index = get_variable_unlimited_dimension_index(fileobj, &
                                                                 variable_name, &
                                                                 broadcast=bcast)
        if (unlim_dim_index .ne. 6) then
            call error("unlimited dimension must be the slowest varying" &
                       //" dimension.")
        endif
        c(unlim_dim_index) = unlim_dim_level
    endif
    if (fileobj%is_root) then
        varid = get_variable_id(fileobj%ncid, &
                                trim(variable_name))
        select type(buf)
            type is (integer(kind=int32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (integer(kind=int64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real32))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            type is (real(kind=real64))
                err = nf90_get_var(fileobj%ncid, &
                                   varid, &
                                   buf, &
                                   start=c &
,count=e &
                                  )
            class default
                call error("unsupported type.")
        end select
        call check_netcdf_code(err)
    endif
    if (bcast) then
        select type(buf)
            type is (integer(kind=int32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (integer(kind=int64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real32))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            type is (real(kind=real64))
                call mpp_broadcast(buf, &
size(buf), &
                                   fileobj%io_root, &
                                   pelist=fileobj%pelist)
            class default
                call error("unsupported type.")
        end select
    endif
end subroutine netcdf_read_data_5d
