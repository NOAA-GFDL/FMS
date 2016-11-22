!----------
!ug support

!>Add a field to a restart object (restart_file_type).  Return the index of the
!!inputted field in the fileObj%var array.
subroutine fms_io_unstructured_setup_one_field(fileObj, &
                                               filename, &
                                               fieldname, &
                                               field_dimension_sizes, &
                                               index_field, &
                                               domain, &
                                               mandatory, &
                                               is_scalar_field, &
                                               data_default, &
                                               longname, &
                                               units, &
                                               read_only, &
                                               owns_data)

   !Inputs/Outputs
    type(restart_file_type),intent(inout)     :: fileObj                 !<A restart object.
    character(len=*),intent(in)               :: filename                !<The name of the restart file.
    character(len=*),intent(in)               :: fieldname               !<The name of a field.
    integer(INT_KIND),dimension(:),intent(in) :: field_dimension_sizes   !<Array of sizes of the dimensions of the inputted field.
    integer(INT_KIND),intent(out)             :: index_field             !<Index of the inputted field in the fileObj%var array.
    type(domainUG),intent(in)                 :: domain                  !<An unstructured mpp domain.
    logical(INT_KIND),intent(in),optional     :: mandatory               !<Flag telling if the field is mandatory for the restart.
    logical(INT_KIND),intent(in),optional     :: is_scalar_field         !<Flag teling if the inputed field is a scalar field.
    real,intent(in),optional                  :: data_default            !<A default value for the data.
    character(len=*),intent(in),optional      :: longname                !<A more descriptive name of the field.
    character(len=*),intent(in),optional      :: units                   !<Units for the field.
    logical(INT_KIND),intent(in),optional     :: read_only               !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional     :: owns_data               !<Tells if the data will be deallocated when the restart object is deallocated.

   !Local variables
    real(DOUBLE_KIND)      :: default_data    !<The "default" data value.  This defaults to MPP_FILL_DOUBLE. Shouldn't this be a real(DOUBLE_KIND)?
    character(len=256)     :: filename2       !<A string used to manipulate the inputted filename.
    integer(INT_KIND)      :: length          !<the length of the (trimmed) inputted file name.
    character(len=256)     :: append_string   !<A string used to append the filename_appendix module variable string to the inputted filename.
    character(len=256)     :: fname           !<A string to hold a file name.
    type(var_type),pointer :: cur_var         !<A convenience pointer.
    integer(INT_KIND)      :: i               !<Loop variable.
    character(len=256)     :: error_msg       !<An error message string.
    integer(INT_KIND)      :: n_field_siz     !<Equal to size(field_dimension_sizes)
    integer(INT_KIND)      :: domain_index    !<Index of the inputted unstructured domain in the domain_UG_array module array.
    integer(INT_KIND)      :: global_size     !<Size of the global domain for the inputted unstructured domain.
    integer(INT_KIND)      :: compute_size    !<Size of the compute domain for the inputted unstructured domain.

   !Make sure that valid dimension sizes were passed in. The
   !array should be structured as such:
   !field_dimension_sizes(1) = unsturctured dimension size of the data domain.
   !field_dimension_sizes(2) = size of an extra dimension of the data domain.
   !field_dimension_sizes(3) = z-dimension size of the data domain.
   !field_dimension_sizes(4) = t-dimension size of the data domain.
   !field_dimension_sizes(5) = size of another extra dimension of the data domain.
    if (size(field_dimension_sizes) .gt. 5) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_setup_one_field:" &
                       //" the field_dimension_sizes array cannot have" &
                       //" more than 5 elements.")
    endif
    if (any(field_dimension_sizes .lt. 0)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_setup_one_field:" &
                       //" all dimensions must have a size that is a non-" &
                       //" negative integer.")
    endif

   !Set the "default" data value for the field.
    if (present(data_default)) then
        default_data = data_default
    else
        default_data = MPP_FILL_DOUBLE
    endif

   !Remove the ".nc" from file name.
    length = len_trim(filename)
    if (filename(length-2:length) .eq. ".nc") then
        filename2 = filename(1:length-3)
    else
        filename2 = filename(1:length)
    endif

   !Append the filename_appendix string to the file name.
   !filename_appendix is a module variable.
    append_string = ""
    if (len_trim(filename_appendix) .gt. 0) then
        append_string = filename_appendix
    endif
    if (len_trim(append_string) .gt. 0) then
        filename2 = trim(filename2)//'.'//trim(append_string)
    endif

   !If necessary, add the correct domain ".tilexxxx" string to the inputted
   !file name.  For a file named foo.nc, this would become foo.tilexxxx.nc.
    call get_mosaic_tile_file_ug(filename2, &
                                 fname, &
                                 domain)

    if (associated(fileObj%var)) then

       !Make sure that the filename stored in fileObj matches the filename
       !returned from get_mosaic_tile_file_ug.
        if (trim(fileObj%name) .ne. trim(fname)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_setup_one_field:" &
                           //" filename = "//trim(fname)//" is not" &
                           //" consistent with the filename of the" &
                           //" restart object = "//trim(fileObj%name))
        endif
    else

       !If any axis has already been registered, then make sure that the
       !filename returned from get_mosaic_tile_file_ug matches the filename
       !stored in the fileObj restart object.  If this is the first axis/
       !field registered to the restart object, then store the filename
       !returned from get_mosaic_tile_file_ug in the restart object.
        if (allocated(fileObj%axes)) then
            if (trim(fileObj%name) .ne. trim(fname)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_setup_one_field:" &
                               //" filename = "//trim(fname)//" is not" &
                               //" consistent with the filename of the" &
                               //" restart object = "//trim(fileObj%name))
            endif
        else
            fileObj%name = trim(fname)
        endif

       !Allocate necessary space in hte restart object.
        allocate(fileObj%var(max_fields))
        allocate(fileObj%p0dr(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p1dr(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p2dr(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p3dr(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p4dr(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p0di(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p1di(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p2di(MAX_TIME_LEVEL_REGISTER,max_fields))
        allocate(fileObj%p3di(MAX_TIME_LEVEL_REGISTER,max_fields))

       !Make sure that the restart file name is not currently being used by
       !an other restart objects.  Should this be fatal?
       !num_registered files is a module variable.
        do i = 1,num_registered_files
            if (trim(fname) .eq. trim(registered_file(i))) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_setup_one_field: " &
                               //trim(fname)//" is already registered with" &
                               //" another restart_file_type data")
                exit
            endif
        enddo

       !Iterate the number of registered restart files, and add the inputted
       !file to the array.  Should this be fatal?
       !max_files_w is a module variable.
        num_registered_files = num_registered_files + 1
        if (num_registered_files .gt. max_files_w) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_setup_one_field:" &
                           //" the number of registered files is greater" &
                           //" than max_files_w.  Please increase" &
                           //" max_files_w in the fms_io_nml namelist")
        endif
        registered_file(num_registered_files) = trim(fname)

       !Set values for the restart object.
       !max_fields is a module variable.
        fileObj%register_id = num_registered_files
        fileObj%max_ntime = field_dimension_sizes(4)
        fileObj%is_root_pe = mpp_domain_UG_is_tile_root_pe(domain)
        fileObj%nvar = 0
        do i = 1,max_fields
            fileObj%var(i)%name = "none"
            fileObj%var(i)%longname = "";
            fileObj%var(i)%units = "none";
            fileObj%var(i)%domain_present = .false.
            fileObj%var(i)%domain_idx = -1
            fileObj%var(i)%is_dimvar = .false.
            fileObj%var(i)%read_only = .false.
            fileObj%var(i)%owns_data = .false.
            fileObj%var(i)%position = CENTER
            fileObj%var(i)%ndim = -1
            fileObj%var(i)%siz(:) = -1
            fileObj%var(i)%gsiz(:) = -1
            fileObj%var(i)%id_axes(:) = -1
            fileObj%var(i)%initialized = .false.
            fileObj%var(i)%mandatory = .true.
            fileObj%var(i)%is = -1
            fileObj%var(i)%ie = -1
            fileObj%var(i)%js = -1
            fileObj%var(i)%je = -1
            fileObj%var(i)%default_data = -1
            fileObj%var(i)%compressed_axis = ""
            fileObj%var(i)%ishift = -1
            fileObj%var(i)%jshift = -1
            fileObj%var(i)%x_halo = -1
            fileObj%var(i)%y_halo = -1
            fileObj%var(i)%domain_UG_index = -1
        enddo
    endif

   !Get the index of the field in the fileObj%var array, if it exists.  If
   !it doesn't exist, set the index to be -1.
    index_field = -1
    do i = 1,fileObj%nvar
        if (trim(fileObj%var(i)%name) .eq. trim(fieldname)) then
            index_field = i
            exit
        endif
    enddo

    if (index_field > 0) then

       !If the field already exists in the fileObj%var array, then update its
       !time level.
        cur_var => null()
        cur_var => fileObj%var(index_field)

       !Make sure that the array of field dimension sizes matches the
       !dimension sizes of the found field for all dimensions except the
       !time level.
        if (cur_var%siz(1) .ne. field_dimension_sizes(1) .or. &
            cur_var%siz(2) .ne. field_dimension_sizes(2) .or. &
            cur_var%siz(3) .ne. field_dimension_sizes(3)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_setup_one_field:" &
                           //" field dimension size mismatch for field " &
                           //trim(fieldname)//" of file "//trim(filename))
        endif

       !Update the time level.
        cur_var%siz(4) = cur_var%siz(4) + field_dimension_sizes(4)
        if (fileObj%max_ntime .lt. cur_var%siz(4)) then
            fileObj%max_ntime = cur_var%siz(4)
        endif
        if (cur_var%siz(4) .gt. MAX_TIME_LEVEL_REGISTER) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_setup_one_field:" &
                           //" the time level of field "//trim(cur_var%name) &
                           //" in file "//trim(fileObj%name)//" is greater" &
                           //" than MAX_TIME_LEVEL_REGISTER(=2), increase" &
                           //" MAX_TIME_LEVEL_REGISTER or check your code.")
        endif
    else

       !If this is a new field, then add it the restart object.
        fileObj%nvar = fileObj%nvar + 1
        if (fileObj%nvar .gt. max_fields) then
            write(error_msg,'(I3,"/",I3)') fileObj%nvar,max_fields
            call mpp_error(FATAL, &
                           "fms_io_unstructured_setup_one_field:" &
                           //" max_fields exceeded, needs increasing," &
                           //" nvar/max_fields = "//trim(error_msg))
        endif
        index_field = fileObj%nvar
        cur_var => null()
        cur_var => fileObj%var(index_field)

       !Copy in the dimension sizes of the data domain (siz, used for
       !writes), and of the global domain (gsiz, used for reads).
        cur_var%registered_with_domain_UG = .true.
        n_field_siz = size(field_dimension_sizes)
        cur_var%siz(1:n_field_siz) = field_dimension_sizes(1:n_field_siz)
        cur_var%gsiz(2) = field_dimension_sizes(2)
        cur_var%gsiz(3) = field_dimension_sizes(3)
        if (n_field_siz .eq. 5) then
            cur_var%gsiz(4) = field_dimension_sizes(5)
        endif

       !Copy in the rest of the data.
        cur_var%name = fieldname
        cur_var%default_data = real(default_data)
        if (present(mandatory)) then
            cur_var%mandatory = mandatory
        endif
        if (present(read_only)) then
            cur_var%read_only = read_only
        endif
        if (present(owns_data)) then
            cur_var%owns_data = owns_data
        endif
        if (present(longname)) then
            cur_var%longname = longname
        else
            cur_var%longname = fieldname
        endif
        if (present(units)) then
            cur_var%units = units
        endif

       !Set the starting and ending indices for the unstructured axis compute
       !domain.  The js and je indices are not used with the unstructured grid,
       !so they are set to -1.
        cur_var%is = 1
        cur_var%ie = cur_var%siz(1)
        cur_var%js = -1
        cur_var%je = -1

       !Set the value of the domain_UG_scalar flag.  This defaults to false.
        cur_var%domain_UG_scalar = .false.
        if (present(is_scalar_field)) then
            cur_var%domain_UG_scalar = is_scalar_field
        endif

       !Modify the vales of cur_var%is and cur_var%ie if necessary.
        if (.not. cur_var%domain_UG_scalar) then
            cur_var%domain_present = .true.

           !Find the index of the domain in the domain_UG_array module array.
           !If it does not exist in the array, then add it.
            domain_index = fms_io_unstructured_lookup_domain(domain)
            if (domain_index .eq. -1) then
                domain_UG_array_size = domain_UG_array_size + 1
                if (domain_UG_array_size .gt. max_domains) then
                    call mpp_error(FATAL, &
                                   "fms_io_unstructured_setup_one_field:" &
                                   //" the number of unstructured grids has" &
                                   //" exceeded max_domains.")
                endif
                domain_index = domain_UG_array_size
                domain_UG_array(domain_index) = domain
            endif
            cur_var%domain_UG_index = domain_index

           !Adjust the starting and ending indices of cur_var to match
           !the global and compute domains for the unstructured domain.
            call mpp_get_UG_global_domain(domain_UG_array(domain_index), &
                                          size=global_size)
            call mpp_get_UG_compute_domain(domain_UG_array(domain_index), &
                                           size=compute_size)
            if (cur_var%siz(1) .ne. global_size .and. cur_var%siz(1) .ne. &
                compute_size) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_setup_one_field:" &
                               //" the field should exist on either the" &
                               //" compute or global domain.")
            endif
            cur_var%is = 1 + (cur_var%siz(1) - compute_size)/2
            cur_var%ie = cur_var%is + compute_size - 1;
            cur_var%gsiz(1) = global_size
        else
            cur_var%domain_present = .false.
            cur_var%gsiz(1) = field_dimension_sizes(1)
        endif
    endif

   !Nullify local pointer.
    cur_var => null()

    return
end subroutine fms_io_unstructured_setup_one_field

!----------
