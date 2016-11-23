!----------
!ug support

!>Read in field data from a restart file associated with an unstructured mpp
!!domain.
subroutine fms_io_unstructured_restore_state(fileObj, &
                                             directory)

   !Inputs/Outputs
    type(restart_file_type),intent(inout) :: fileObj   !<A restart object.
    character(len=*),intent(in),optional  :: directory !<Directory where input files live.

   !Local variables
    type(domainUG),pointer                      :: domain              !<A pointer to an unstructured mpp domain.
    character(len=128)                          :: dir                 !<Directory where the input files live.
    character(len=256)                          :: restartpath         !<The input file path (dir/file).
    logical(INT_KIND)                           :: fexist              !<Flag telling if a file exists.
    type(domainUG),pointer                      :: io_domain           !<Pointer to an unstructured domain I/O domain.
    character(len=256)                          :: filename            !<String used to check the existence of various file names.
    integer(INT_KIND)                           :: io_tile_id          !<I/O domain tile id.
    integer(INT_KIND)                           :: funit               !<A file unit returned from mpp_open.
    integer(INT_KIND)                           :: ndim                !<The number of dimensions in an input file.
    integer(INT_KIND)                           :: nvar                !<The number of fields in an input file.
    integer(INT_KIND)                           :: natt                !<The number of attributes in an input file.
    integer(INT_KIND)                           :: ntime               !<The number of time levels in an input file.
    character(len=8)                            :: axis_compressed_id  !<The string that identifies an unstructured axis in the input restart file.
    character(len=8)                            :: axis_cartesian_id   !<The string that identifies a cartesian axis in the input restart file.
    type(axistype),dimension(:),allocatable     :: axes                !<An array of axes from the input restart file.
    integer(INT_KIND)                           :: axis_length         !<The size of the data for an axis
    real,dimension(:),allocatable               :: axis_data           !<An array that holds axis data read in from the input restart file.
    type(fieldtype),dimension(:),allocatable    :: fields              !<An array of fields from the input restart file.
    type(var_type),pointer                      :: cur_var             !<Pointer to a field.
    character(len=80)                           :: varname             !<The name of a field.
    logical(INT_KIND)                           :: check_exist         !<Flag that tells whether or not a checksum attribute exists in the file for a field.
    integer(LONG_KIND),dimension(3)             :: checksum_file       !<Check-sum values read in from the file.
    logical(INT_KIND)                           :: is_there_a_checksum !<Flag telling if a check-sum should be calculated for a read-in field.
    integer(INT_KIND)                           :: tlev                !<A time value.
    integer(LONG_KIND)                          :: checksum_data       !<Check-sum value calculated for the read-in data.
    real                                        :: r0d                 !<Used to read-in an integer scalar field.
    real,dimension(:),allocatable               :: r1d                 !<Used to read-in an integer 1D field.
    real,dimension(:,:),allocatable             :: r2d                 !<Used to read-in an integer 2D field.
    character(len=256)                          :: mesg                !<An error message.
    integer(INT_KIND)                           :: j                   !<Loop variable.
    integer(INT_KIND)                           :: n                   !<Loop variable.
    integer(INT_KIND)                           :: l                   !<Loop variable.
    integer(INT_KIND)                           :: k                   !<Loop variable.


   !Make sure at least one field was registered to the restart object.
    if (.not. associated(fileObj%var)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_restore_state:" &
                       //" the restart object does not contain any fields.")
    endif

   !Make sure that an unstructured axis has been registered to this restart
   !object.  Also make sure that a x- and y-axis are registered, since these
   !will be needed to convert from the unstructured domain back to the
   !structured domain later.  Make sure that no "compressed" axes were
   !registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_restore_state:" &
                       //" there are no registered axes for the restart" &
                       //" object "//trim(fileObj%name))
    endif
    if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_restore_state:" &
                       //" an unstructured axis was not registered for the" &
                       //" restart object "//trim(fileObj%name))
    endif
    if (.not. associated(fileObj%axes(XIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_restore_state:" &
                       //" a x-axis was not registered for the restart" &
                       //" object "//trim(fileObj%name))
    endif
    if (.not. associated(fileObj%axes(YIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_restore_state:" &
                       //" a y-axis was not registered for the restart" &
                       //" object "//trim(fileObj%name))
    endif
    if (allocated(fileObj%axes(CIDX)%idx) .or. &
        allocated(fileObj%axes(HIDX)%idx)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_restore_state: the restart object " &
                       //trim(fileObj%name)//" cannot contain both a" &
                       //" compressed and an unstructured axis.")
    endif

   !Point to the unstructured domain associated with the unstructured axis.
    domain => null()
    domain => fileObj%axes(UNSTRUCTURED_GRID_INDEX)%domain_ug

   !Make sure that all fields are associated with the correct domain.
    do j = 1,fileObj%nvar
       !Make sure that all fields in the restart object have been registered
       !with an unstructured mpp domain.
        if (.not. fileObj%var(j)%registered_with_domain_UG) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_restore_state:" &
                           //" the field "//trim(fileObj%var(j)%name) &
                           //" in the restart object "//trim(fileObj%name) &
                           //" was not registered with an unstructured" &
                           //" mpp domain.")
        endif

       !For non-scalar fields, make sure that the unstructured domain with
       !which they were registered matches the domain associated with the
       !unstructured axis.
        if (.not. fileObj%var(j)%domain_UG_scalar) then
            if (.not. mpp_compare_UG_domains(domain, &
                                             domain_UG_array(fileObj%var(j)%domain_UG_index))) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the unstructured domain associated with" &
                               //" field "//trim(fileObj%var(j)%name) &
                               //" in the restart object "//trim(fileObj%name) &
                               //" does not match the unstructured domain" &
                               //" associated with the registered" &
                               //" unstructured axis.")
            endif
        endif
    enddo

   !Set the directory where the input files live.
    dir = "INPUT"
    if (present(directory)) then
        dir = trim(directory)
    endif

   !Set the name of the input file (including the path to it).
    if (len_trim(dir) .gt. 0) then
        restartpath = trim(dir)//"/"//trim(fileObj%name)
    else
        restartpath = trim(fileObj%name)
    endif

   !Check if the input file has the I/O domain tile id appended onto its
   !end.
    fexist = .false.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)
    io_tile_id = mpp_get_UG_domain_tile_id(io_domain)
    write(filename,'(a,i4.4)') trim(restartpath)//'.',io_tile_id
    inquire(file=trim(filename),exist=fexist)
    if (.not. fexist) then
        write(filename,'(a,i6.6)') trim(restartpath)//'.',io_tile_id
        inquire(file=trim(filename),exist=fexist)
    endif
    io_domain => null()

   !If the file exists, open it.  If not, then throw a fatal error.  In
   !addition to opening the file, mpp_open also reads the metadata for each
   !axis and field.
    if (fexist) then
        call mpp_open(funit, &
                      trim(restartpath), &
                      form=form, &
                      action=MPP_RDONLY, &
                      threading=MPP_MULTI, &
                      domain_ug=domain)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_restore_state:" &
                       //" unable to find any restart files " &
                       //"specified by "//trim(restartpath))
    endif

   !Get the number of dimensions, fields, attributes, and time levels
   !for the input file.
    call mpp_get_info(funit, &
                      ndim, &
                      nvar, &
                      natt, &
                      ntime)

   !Get the axes from the file.
    allocate(axes(ndim))
    call mpp_get_axes(funit, &
                      axes)

   !Make sure the axis data from the input restart file matches the axis
   !data registered to the restart object.
    do j = 1,ndim
        call mpp_get_atts(axes(j), &
                          cartesian=axis_cartesian_id)
        call mpp_get_atts(axes(j), &
                          compressed=axis_compressed_id)
        if (trim(axis_compressed_id) .eq. "U") then

           !Make sure that the unstructured axis defined in the input restart
           !file contains the same data as the unstructured axis registered
           !to the restart object.
            axis_length = mpp_get_axis_length(axes(j))
            if (axis_length .ne. &
                size(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered unstructured axis has a" &
                               //" different size than the unstructured axis" &
                               //" in the input restart file.")
            endif
            allocate(axis_data(axis_length))
            call mpp_get_axis_data(axes(j), &
                                   axis_data)
            if (maxval(abs(int(axis_data) - &
                fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) .ne. 0) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered unstructured axis has" &
                               //" different data than the unstructured" &
                               //" axis in the input restart file")
            endif
            deallocate(axis_data)
        elseif (trim(axis_cartesian_id) .eq. "X") then

           !Make sure that x-axis defined in the input restart
           !file contains the same data as the x-axis registered
           !to the restart object.
            axis_length = mpp_get_axis_length(axes(j))
            if (axis_length .ne. &
                size(fileObj%axes(XIDX)%data)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered x-axis has a" &
                               //" different size than the x-axis" &
                               //" in the input restart file.")
            endif
            allocate(axis_data(axis_length))
            call mpp_get_axis_data(axes(j), &
                                   axis_data)
            if (maxval(abs(axis_data - fileObj%axes(XIDX)%data)) .ne. 0.0) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered x-axis has" &
                               //" different data than the" &
                               //" x-axis in the input restart file")
            endif
            deallocate(axis_data)
        elseif (trim(axis_cartesian_id) .eq. "Y") then

           !Make sure that y-axis defined in the input restart
           !file contains the same data as the y-axis registered
           !to the restart object.
            axis_length = mpp_get_axis_length(axes(j))
            if (axis_length .ne. &
                size(fileObj%axes(YIDX)%data)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered y-axis has a" &
                               //" different size than the y-axis" &
                               //" in the input restart file.")
            endif
            allocate(axis_data(axis_length))
            call mpp_get_axis_data(axes(j), &
                                   axis_data)
            if (maxval(abs(axis_data - fileObj%axes(YIDX)%data)) .ne. 0.0) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered y-axis has" &
                               //" different data than the" &
                               //" y-axis in the input restart file")
            endif
            deallocate(axis_data)
        elseif (trim(axis_cartesian_id) .eq. "Z") then

           !Make sure that z-axis defined in the input restart (if it exists)
           !file contains the same data as the z-axis registered
           !to the restart object.
            if (.not. associated(fileObj%axes(ZIDX)%data)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the z-axis exists in the input restart" &
                               //" file, but was not registered in the" &
                               //" inputted restart object.")
            endif
            axis_length = mpp_get_axis_length(axes(j))
            if (axis_length .ne. &
                size(fileObj%axes(ZIDX)%data)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered z-axis has a" &
                               //" different size than the z-axis" &
                               //" in the input restart file.")
            endif
            allocate(axis_data(axis_length))
            call mpp_get_axis_data(axes(j), &
                                   axis_data)
            if (maxval(abs(axis_data - fileObj%axes(ZIDX)%data)) .ne. 0.0) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered z-axis has" &
                               //" different data than the" &
                               //" z-axis in the input restart file")
            endif
            deallocate(axis_data)
        elseif (trim(axis_cartesian_id) .eq. "CC") then

           !Make sure that cc-axis defined in the input restart (if it exists)
           !file contains the same data as the cc-axis registered
           !to the restart object.
            if (.not. associated(fileObj%axes(CCIDX)%data)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the cc-axis exists in the input restart" &
                               //" file, but was not registered in the" &
                               //" inputted restart object.")
            endif
            axis_length = mpp_get_axis_length(axes(j))
            if (axis_length .ne. &
                size(fileObj%axes(CCIDX)%data)) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered cc-axis has a" &
                               //" different size than the cc-axis" &
                               //" in the input restart file.")
            endif
            allocate(axis_data(axis_length))
            call mpp_get_axis_data(axes(j), &
                                   axis_data)
            if (maxval(abs(axis_data - fileObj%axes(CCIDX)%data)) .ne. 0.0) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" the registered cc-axis has" &
                               //" different data than the" &
                               //" cc-axis in the input restart file")
            endif
            deallocate(axis_data)
        elseif (trim(axis_cartesian_id) .ne. "T") then

           !Throw a fatal error if an unexpected axis exists in the input
           !restart file.
            call mpp_error(FATAL, &
                           "fms_io_unstructured_restore_state:" &
                           //" unexpected axis contained in the input" &
                           //" restart file.")
        endif
    enddo
    deallocate(axes)

   !Get all of the fields from the input restart file.
    allocate(fields(nvar))
    call mpp_get_fields(funit, &
                        fields(1:nvar))
    do j = 1,fileObj%nvar

       !Set the convenience pointer to point to a field.
        cur_var => null()
        cur_var => fileObj%var(j)

       !Set the starting and ending indices of the compute domain.
        do l = 1,nvar
           !Get the name of the field.
            call mpp_get_atts(fields(l), &
                              name=varname)

           !Check if the name of the field in the file matches the
           !name of the field in the restart object.
            if (lowercase(trim(varname)) .eq. &
                lowercase(trim(cur_var%name))) then

               !If the names match, then set the "initialized" flag to
               !true.  This flag tells if a field has been read in or
               !not.
                cur_var%initialized = .true.

               !Determine whether or not a check-sum is required for this
               !field.  If so, then get the check-sum from the file.  If
               !not, then set the check-sum equal to zero.
                check_exist = mpp_attribute_exist(fields(l), &
                                                  "checksum")
                if (check_exist) then
                    call mpp_get_atts(fields(l), &
                                      checksum=checksum_file)
                    is_there_a_checksum = .true.
                else
                    checksum_file = 0
                    is_there_a_checksum = .false.
                endif

               !Use the "checksum_required" module varaiable to override
               !the use of check-sums if desired.
                if (.not. checksum_required) then
                    is_there_a_checksum = .false.
                endif

                do k = 1,cur_var%siz(4)
                    tlev = k
                    if (associated(fileObj%p0dr(k,j)%p)) then
                        call mpp_read(funit, &
                                      fields(l), &
                                      fileObj%p0dr(k,j)%p, &
                                      tlev)
                        if (is_there_a_checksum) then
                            checksum_data = mpp_chksum(fileObj%p0dr(k,j)%p, &
                                                       (/mpp_pe()/), &
                                                       mask_val=cur_var%default_data)
                        endif
                    elseif (associated(fileObj%p1dr(k,j)%p)) then
                        call mpp_io_unstructured_read(funit, &
                                                      fields(l), &
                                                      domain, &
                                                      fileObj%p1dr(k,j)%p, &
                                                      tindex=tlev)
                        if (is_there_a_checksum) then
                            checksum_data = mpp_chksum(fileObj%p1dr(k,j)%p, &
                                                       mask_val=cur_var%default_data)
                        endif
                    elseif (associated(fileObj%p2dr(k,j)%p)) then
                        call mpp_io_unstructured_read(funit, &
                                                      fields(l), &
                                                      domain, &
                                                      fileObj%p2dr(k,j)%p, &
                                                      tindex=tlev)
                        if (is_there_a_checksum) then
                            checksum_data = mpp_chksum(fileObj%p2dr(k,j)%p, &
                                                       mask_val=cur_var%default_data)
                        endif
                    elseif (associated(fileObj%p3dr(k,j)%p)) then
                        call mpp_io_unstructured_read(funit, &
                                                      fields(l), &
                                                      domain, &
                                                      fileObj%p3dr(k,j)%p, &
                                                      tindex=tlev)
                        if (is_there_a_checksum) then
                            checksum_data = mpp_chksum(fileObj%p3dr(k,j)%p, &
                                                       mask_val=cur_var%default_data)
                        endif
                    elseif (associated(fileObj%p0di(k,j)%p)) then
                        r0d = 0
                        call mpp_read(funit, &
                                      fields(l), &
                                      r0d, &
                                      tlev)
                        fileObj%p0di(k,j)%p = int(r0d)
                        if (is_there_a_checksum) then
                            checksum_data = int(fileObj%p0di(k,j)%p,kind=LONG_KIND)
                        endif
                    elseif (associated(fileObj%p1di(k,j)%p)) then
                        allocate(r1d(cur_var%siz(1)))
                        r1d = 0
                        call mpp_io_unstructured_read(funit, &
                                                      fields(l), &
                                                      domain, &
                                                      r1d, &
                                                      tindex=tlev)
                        fileObj%p1di(k,j)%p = int(r1d)
                        if (is_there_a_checksum) then
                            checksum_data = mpp_chksum(fileObj%p1di(k,j)%p, &
                                                       mask_val=cur_var%default_data)
                        endif
                        deallocate(r1d)
                    elseif (associated(fileObj%p2di(k,j)%p)) then
                        allocate(r2d(cur_var%siz(1),cur_var%siz(2)))
                        r2d = 0
                        call mpp_io_unstructured_read(funit, &
                                                      fields(l), &
                                                      domain, &
                                                      r2d, &
                                                      tindex=tlev)
                        fileObj%p2di(k,j)%p = int(r2d)
                        if (is_there_a_checksum) then
                            checksum_data = mpp_chksum(fileObj%p2di(k,j)%p, &
                                                       mask_val=cur_var%default_data)
                        endif
                        deallocate(r2d)
                    else
                        call mpp_error(FATAL, &
                                       "fms_io_unstructured_restore_state:" &
                                       //" domain is present for the field " &
                                       //trim(varname)//" of file " &
                                       //trim(fileObj%name)//", but" &
                                       //" none of p1dr, p2dr, p3dr, p1di" &
                                       //" and p2di are associated.")
                    endif

                   !Make sure that the calculated checksum matches the
                   !checksum from the file.
                    if (is_there_a_checksum .and. (checksum_file(k) &
                        .ne. checksum_data)) then
                        write (mesg,'(a,Z16,a,Z16,a)') "Checksum of input field " &
                                                       //uppercase(trim(varname)) &
                                                       //" ",checksum_data, &
                                                       " does not match value ", &
                                                       checksum_file(k),  &
                                                       " stored in " &
                                                       //uppercase(trim(fileObj%name) &
                                                       //"." )
                        call mpp_error(FATAL, &
                                       "fms_io_unstructured_restore_state: " &
                                       //trim(mesg))
                    endif
                enddo
               !Start the search for the next restart variable.
                exit
            endif
        enddo
    enddo
    deallocate(fields)

   !Close the input file.
    call close_file(funit)

   !Make sure that all mandatory fields have been read in from the input files.
    do j = 1,fileObj%nvar
        if (.not. fileObj%var(j)%initialized) then
            if (fileObj%var(j)%mandatory) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_restore_state:" &
                               //" unable to find mandatory variable " &
                               //trim(fileObj%var(j)%name)//" in restart" &
                               //" file "//trim(fileObj%name))
            endif
        endif
    enddo

   !Nullify local pointers.
    domain => null()
    cur_var => null()

   !Print out the check-sum values.
    if (print_chksum) then
        call write_chksum(fileObj, &
                          MPP_RDONLY)
    endif

    return
end subroutine fms_io_unstructured_restore_state

!----------
