include(`macros.m4')dnl
!> @brief Add a restart variable to a netcdf file.
`subroutine netcdf_add_restart_variable_'NUM_DIMS`d(fileobj, &'
                                                    variable_name, &
                                                    vdata, &
                                                    dimensions)
    class(NetcdfFile_t),intent(inout) :: fileobj !< File object.
    character(len=*),intent(in) :: variable_name !< Variable name.
    class(*),dim_declare(NUM_DIMS) intent(in),target :: vdata !< Pointer to
                                                              !! variable data.
    character(len=*),dimension(ifelse(NUM_DIMS,0,1,:)),intent(in),optional :: dimensions !< Dimension names.
    character(len=8) :: buf
ifelse(NUM_DIMS,0,,
   `integer :: ndims
    integer :: vdata_rank')
    call add_restart_var_to_array(fileobj, &
                                  variable_name)
   `fileobj%restart_vars(fileobj%num_restart_vars)%data'NUM_DIMS`d => vdata'
    if (.not. fileobj%is_readonly) then
        call get_data_type_string(vdata, &
                                  buf)
ifelse(NUM_DIMS,0,
       `if (present(dimensions)) then
            if (.not. is_dimension_unlimited(fileobj,dimensions(1),.true.)) then
                call error("a scalar input variable can only have an unlimited" &
                           //" dimension.")
            endif
        endif',
       `if (.not. present(dimensions)) then
            call error("dimension names required if the file is" &
                       //" not read-only.")
        endif
        ndims = size(dimensions)
        vdata_rank = size(shape(vdata))
        if (ndims .eq. vdata_rank+1) then
            if (.not. is_dimension_unlimited(fileobj,dimensions(ndims),.true.)) then
                call error("the slowest dimension must be unlimited.")
            endif
        elseif (ndims .ne. vdata_rank) then
            call error("rank mismatch between vdata and dimensions arrays.")
        endif')
        call netcdf_add_variable(fileobj, &
                                 variable_name, &
                                 buf, &
                                 dimensions)
    endif
`end subroutine netcdf_add_restart_variable_'NUM_DIMS`d'
