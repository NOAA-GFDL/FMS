function mpp_nc_create(name, mode, pelist, ncid) result(error_code)
    character(len=*), intent(in) :: name
    integer, intent(in) :: mode
    integer, intent(in) :: pelist(:)

    integer, intent(out) :: ncid

    integer :: idx
    integer :: io_info
    integer :: error_code

    ! NOTE: nf_create never takes NF_NOWRITE
    call MPI_Info_create(io_info, error_code)
    call set_hint(io_info, 'WR')

    idx = get_peset(pelist)
    error_code = nf_create_par(name, mode, peset(idx)%id, io_info, ncid)
end function mpp_nc_create


function mpp_nc_open(name, mode, pelist, ncid) result(error_code)
    character(len=*), intent(in) :: name
    integer, intent(in) :: mode
    integer, intent(in) :: pelist(:)

    integer, intent(out) :: ncid

    integer :: idx
    integer :: io_info
    integer :: error_code

    if (IAND(mode, NF_NOWRITE) /= 0) then
        call MPI_Info_create(io_info, error_code)
        call set_hint(io_info, 'RD')
    else
        call MPI_Info_create(io_info, error_code)
        call set_hint(io_info, 'WR')
    end if

    idx = get_peset(pelist)
    error_code = nf_open_par(name, mode, peset(idx)%id, io_info, ncid)
end function mpp_nc_open


subroutine set_hint(info, rwtype)
    integer, intent (inout) :: info
    character (len=*), intent(in) :: rwtype

    integer :: ierror, ihint, lenvar, ierr
    character (len=12) :: ncblk, stpcnt, stpsze, cbnodes
    integer, parameter :: nhints = 23
    character (len=25) :: hint_name(nhints), hint_val(nhints), hint_in

    hint_name(1) = "direct_read"
    hint_name(2) = "direct_write"
    hint_name(3) = "romio_lustre_co_ratio"
    hint_name(4) = "romio_lustre_coll_threshold"
    hint_name(5) = "romio_lustre_ds_in_coll"
    hint_name(6) = "cb_buffer_size"
    hint_name(7) = "romio_cb_read"
    hint_name(8) = "romio_cb_write"
    hint_name(9) = "cb_nodes"
    hint_name(10) = "romio_no_indep_rw"
    hint_name(11) = "romio_cb_pfr"
    hint_name(12) = "romio_cb_fr_types"
    hint_name(13) = "romio_cb_fr_alignment"
    hint_name(14) = "romio_cb_ds_threshold"
    hint_name(15) = "romio_cb_alltoall"
    hint_name(16) = "ind_rd_buffer_size"
    hint_name(17) = "ind_wr_buffer_size"
    hint_name(18) = "romio_ds_read"
    hint_name(19) = "romio_ds_write"
    hint_name(20) = "cb_config_list"
    hint_name(21) = "striping_unit"
    hint_name(22) = "striping_factor"
    hint_name(23) = "romio_lustre_start_iodevice"

    do ihint = 1, nhints
        call get_environment_variable( &
            trim(hint_name(ihint))//'_'//trim(rwtype), &
            hint_in, &
            lenvar, &
            ierror &
        )
        if (ierror == 0) then
            call MPI_Info_set(info, HINT_NAME(ihint),HINT_IN,ierr)
            if (ierr == 0) then
                if(mpp_pe() == mpp_root_pe()) &
                    ! TODO: Replace unit 6 with variable?
                    write (6, '(A,A,A,A,A)') rwtype, &
                        ':SET MPIIO_HINT for ', trim(HINT_NAME(ihint)), &
                        ' AS ', trim(HINT_IN)
            end if
        endif
    end do
end subroutine set_hint
