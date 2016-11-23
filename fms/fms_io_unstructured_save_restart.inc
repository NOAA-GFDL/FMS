!----------
!ug support

!>Write out metadata and data for axes and fields to a restart file
!!associated with an unstructured mpp domain.
subroutine fms_io_unstructured_save_restart(fileObj, &
                                            time_stamp, &
                                            directory, &
                                            append, &
                                            time_level)

   !Inputs/Outputs
    type(restart_file_type),intent(inout),target :: fileObj     !<A restart object.
    character(len=*),intent(in),optional         :: time_stamp  !<A time stamp for the file.
    character(len=*),intent(in),optional         :: directory   !<The directory where the restart file lives.
    logical(INT_KIND),intent(in),optional        :: append      !<Flag telling whether to append to or overwrite the restart file.
    real,intent(in),optional                     :: time_level  !<A time level value (do not specify a kind value).

   !Optional arguments:

   !If neither append or time_level is present:
   !  routine writes both meta data and field data.

   !If append is present and append=.true.:
   !   Only field data is written.
   !   The field data is appended to a new time level.
   !   time_level must also be present and it must be >= 0.0
   !   The value of time_level is written as a new value of the time axis data.

   !If time_level is present and time_level < 0.0:
   !   A new file is opened and only the meta data is written.

   !If append is present and append=.false.:
   !   Behaves the same was as if it were not present. That is, meta data is
   !   written and whether or not field data is written is determined by time_level.

   !Local variables
    type(domainUG),pointer                      :: domain            !<A pointer to an unstructured mpp domain.
    integer(INT_KIND)                           :: mpp_action        !<Parameter specifying how the file will be acted on (overwritten or appended to).
    logical(INT_KIND)                           :: write_meta_data   !<Flag telling whether or not metadata will be written to the restart file.
    logical(INT_KIND)                           :: write_field_data  !<Flag telling whether or not field data will be written to the restart file.
    character(len=128)                          :: dir               !<Directory where the restart file lives.
    character(len=80)                           :: restartname       !<The name of the restart file.
    character(len=256)                          :: restartpath       !<The restart file path (dir/file).
    integer(INT_KIND)                           :: funit             !<The file unit returned by mpp_open.
    type(ax_type),pointer                       :: axis              !<A pointer to an fms_io_axis_type.
    type(axistype)                              :: x_axis            !<An mpp_io_axis_type, used to write the x-axis to the restart file.
    type(axistype)                              :: y_axis            !<An mpp_io_axis_type, used to write the y-axis to the restart file.
    type(axistype)                              :: z_axis            !<An mpp_io_axis_type, used to write the z-axis to the restart file.
    logical(INT_KIND)                           :: z_axis_defined    !<Flag telling whether or not a z-axis has been define for the inputted restart object.
    type(axistype)                              :: cc_axis           !<An mpp_io_axis_type, used to write the cc-axis (???) to the restart file.
    logical(INT_KIND)                           :: cc_axis_defined   !<Flag telling whether or not a cc-axis (???) has been define for the inputted restart object.
    type(axistype)                              :: unstructured_axis !<An mpp_io_axis_type, used to write the unstructured axis to the restart file.
    type(axistype)                              :: t_axis            !<An mpp_io_axis_type, used to write the t-axis to the restart file.
    type(var_type),pointer                      :: cur_var           !<A pointer to an fms_io_field_type.
    integer(INT_KIND)                           :: num_var_axes      !<Number of dimensions for a field.
    type(axistype),dimension(4)                 :: var_axes          !<Array of axis for each field.
    integer(INT_KIND)                           :: cpack             !<(Number of bits in a real(8))/(Number of bits in a real)
    integer(LONG_KIND),dimension(:),allocatable :: check_val         !<An array of check-sums of a field at each time level.
    real                                        :: tlev              !<Time value for a time level (do not specify a kind value).
    real                                        :: r0d               !<Used to convert a scalar integer field into a scalar real field.
    real,dimension(:),allocatable               :: r1d               !<Used to convert a 1D integer field into a 1D real field.
    real,dimension(:,:),allocatable             :: r2d               !<Used to convert a 2D integer field into a 2D real field.
    integer(INT_KIND)                           :: j                 !<Loop variable.
    integer(INT_KIND)                           :: k                 !<Loop variable.

   !Make sure at least one field was registered to the restart object.
    if (.not. associated(fileObj%var)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_save_restart:" &
                       //" the restart object does not conatin any fields.")
    endif

   !If all fields in the file are read only, then simply return without
   !writing any data to the restart file.  If the restart file does not yet
   !exist, it is not created.
    if (all_field_read_only(fileObj)) then
        return
    endif

   !Make sure that an unstructured axis has been registered to this restart
   !object.  Also make sure that a x- and y-axis are registered, since these
   !will be needed to convert from the unstructured domain back to the
   !structured domain later.  Make sure that no "compressed" axes were
   !registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_save_restart: there are no" &
                       //" registered axes for the file "//trim(fileObj%name))
    endif
    if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_save_restart: an unstructured" &
                       //" axis was not registered for the file " &
                       //trim(fileObj%name))
    endif
    if (.not. associated(fileObj%axes(XIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_save_restart: a x-axis was not" &
                       //" registered for the file "//trim(fileObj%name))
    endif
    if (.not. associated(fileObj%axes(YIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_save_restart: a y-axis was not" &
                       //" registered for the file "//trim(fileObj%name))
    endif
    if (allocated(fileObj%axes(CIDX)%idx) .or. &
        allocated(fileObj%axes(HIDX)%idx)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_save_restart: the restart file " &
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
                            "fms_io_unstructured_save_restart:" &
                            //" the field "//trim(fileObj%var(j)%name) &
                            //" in the file "//trim(fileObj%name) &
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
                               "fms_io_unstructured_save_restart:" &
                               //" the unstructured domain associated with" &
                               //" field "//trim(fileObj%var(j)%name) &
                               //" in the file "//trim(fileObj%name) &
                               //" does not match the unstructured domain" &
                               //" associated with the registered" &
                               //" unstructured axis.")
            endif
        endif
    enddo

   !If necessary, make sure a valid set of optional arguments was provided.
    if (present(append)) then
        if (append .and. .not. present(time_level)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_save_compressed_restart:" &
                           //" a time_level must be present when" &
                           //" append=.true. for file "//trim(fileObj%name))
        endif
    endif

   !Determine whether or not metadata will be written to the restart file.  If
   !no optional arguments are specified, metadata will be written to the file,
   !with any old data overwritten.  If the optional append flag is true, then
   !it is assumed that the metadata already exists in the file, and thus
   !metadata will not be written to the file.
    mpp_action = MPP_OVERWR
    write_meta_data  = .true.
    if (present(append)) then
        if (append) then
            mpp_action = MPP_APPEND
            write_meta_data = .false.
            if (time_level .lt. 0.0) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_save_restart:" &
                               //" the inputted time_level cannot be" &
                               //" negative when append is .true." &
                               //" for file "//trim(fileObj%name))
            endif
        endif
    endif

   !Determine whether or not field data will be written to the restart file.
   !Field data will be written to the restart file unless a negative
   !time_level value is passed in.
    write_field_data = .true.
    if (present(time_level)) then
        if (time_level .lt. 0) then
            write_field_data = .false.
        endif
    endif

   !Set the directory where the restart file lives.  This defaults to 
   !"./RESTART".
    dir = "RESTART"
    if (present(directory)) then
        dir = trim(directory)
    endif

   !Set the name of the restart file excluding its path.
   !time_stamp_restart is a module variable.
    restartname = trim(fileObj%name)
    if (time_stamp_restart) then
        if (present(time_stamp)) then
            if (len_trim(restartname) + len_trim(time_stamp) .gt. 79) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_save_restart:" &
                               //" length of restart file name including" &
                               //" time stamp is greater than allowed" &
                               //" restart file name length.")
            endif
            restartname = trim(time_stamp)//"."//trim(restartname)
        endif
    endif

   !Set the name of the restart file including the path to it.
    if (len_trim(dir) .gt. 0) then
        restartpath = trim(dir)//"/"//trim(restartname)
    else
        restartpath = trim(restartname)
    endif

   !Open the restart file.
    call mpp_open(funit, &
                  trim(restartpath), &
                  action=mpp_action, &
                  form=form, &
                  is_root_pe=fileObj%is_root_pe, &
                  domain_ug=domain)

   !Write out the metadata for the axes and fields.
    axis => null()
    cur_var => null()
    if (write_meta_data) then

       !Write out the metadata for the x-axis to the restart file.
        axis => fileObj%axes(XIDX)
        call mpp_write_meta(funit, &
                            x_axis, &
                            axis%name, &
                            axis%units, &
                            axis%longname, &
                            data=axis%data, &
                            cartesian="X")
        axis => null()

       !Write out the metadata for the y-axis to the restart file.
        axis => fileObj%axes(YIDX)
        call mpp_write_meta(funit, &
                            y_axis, &
                            axis%name, &
                            axis%units, &
                            axis%longname, &
                            data=axis%data, &
                            cartesian="Y")
        axis => null()

       !If it is registered, then write out the metadata for the z-axis
       !to the restart file.
        if (associated(fileObj%axes(ZIDX)%data)) then
            axis => fileObj%axes(ZIDX)
            call mpp_write_meta(funit, &
                                z_axis, &
                                axis%name, &
                                axis%units, &
                                axis%longname, &
                                data=axis%data, &
                                cartesian="Z")
            axis => null()
            z_axis_defined = .true.
        else
            z_axis_defined = .false.
        endif

       !If it is registered, then write out the metadata for the cc-axis (???)
       !to the restart file.
        if (associated(fileObj%axes(CCIDX)%data)) then
            axis => fileObj%axes(CCIDX)
            call mpp_write_meta(funit, &
                                cc_axis, &
                                axis%name, &
                                axis%units, &
                                axis%longname, &
                                data=axis%data, &
                                cartesian="CC")
            axis => null()
            cc_axis_defined = .true.
        else
            cc_axis_defined = .false.
        endif

       !Write out the unstructured axis to the restart file.
        axis => fileObj%axes(UNSTRUCTURED_GRID_INDEX)
        call mpp_write_meta(funit, &
                            unstructured_axis, &
                            axis%name, &
                            axis%units, &
                            axis%longname, &
                            axis%unstructured_axis_data, &
                            compressed="U")
        axis => null()

       !Write out the time axis to the restart file.
        if (associated(fileObj%axes(TIDX)%data)) then
            axis => fileObj%axes(TIDX)
            call mpp_write_meta(funit, &
                                t_axis, &
                                axis%name, &
                                units=axis%units, &
                                longname=axis%longname, &
                                cartesian="T", &
                                calendar=axis%calendar)
            axis => null()
        else
            call mpp_write_meta(funit, &
                                t_axis, &
                                "Time", &
                                "time level", &
                                "Time", &
                                cartesian="T")
        endif

       !Loop through the fields and write out the metadata.
        do j = 1,fileObj%nvar

           !Point to the current field.
            cur_var => fileObj%var(j)

           !Cycle to the next field if the current field is read only.
            if (cur_var%read_only) then
                cur_var => null()
                cycle
            endif

           !Make sure the field has a valid number of time levels.
            if (cur_var%siz(4) .gt. 1 .and. cur_var%siz(4) .ne. &
                fileObj%max_ntime) then
                call mpp_error(FATAL, &
                               "fms_io_unstructured_save_restart: " &
                               //trim(cur_var%name)//" in file " &
                               //trim(fileObj%name)//" has more than one" &
                               //" time level, but the number of time levels" &
                               //" is not equal to max_ntime")
            endif

           !Determine the dimensions for the field.
            if (cur_var%ndim .eq. 0) then
               !For a scalar field foo, it is assumed foo = foo(t)
                num_var_axes = 1
                var_axes(1) = t_axis
            elseif (cur_var%ndim .eq. 1) then
               !For a 1D field foo, it is assumed foo = foo(unstructured) or
               !foo = foo(unstructured,t)
                num_var_axes = 1
                var_axes(1) = unstructured_axis
                if (cur_var%siz(4) .eq. fileObj%max_ntime) then
                    num_var_axes = 2
                    var_axes(2) = t_axis
                endif
            elseif (cur_var%ndim .eq. 2) then
               !For a 2D field foo, it is assumed foo = foo(unstructured,z) or
               !foo = foo(unstructured,z,t).  Make sure that a z-axis has
               !been registered to the restart object.
                if (.not. z_axis_defined) then
                    call mpp_error(FATAL, &
                                   "fms_io_unstructured_save_restart:" &
                                   //" two-dimensional field " &
                                   //trim(cur_var%name)//" requires that" &
                                   //"a z-axis must be registered to the" &
                                   //" restart file "//trim(fileObj%name))
                endif
                num_var_axes = 2
                var_axes(1) = unstructured_axis
                var_axes(2) = z_axis
                if (cur_var%siz(4) .eq. fileObj%max_ntime) then
                    num_var_axes = 3
                    var_axes(3) = t_axis
                endif
            elseif(cur_var%ndim .eq. 3) then
               !For a 3D field foo, it is assumed foo = foo(unstructured,z,cc)
               !or foo = foo(unstructured,z,cc,t)
                if (.not. z_axis_defined) then
                    call mpp_error(FATAL, &
                                   "fms_io_unstructured_save_restart:" &
                                   //" three-dimensional field " &
                                   //trim(cur_var%name)//" requires that" &
                                   //"a z-axis must be registered to the" &
                                   //" restart file "//trim(fileObj%name))
                endif
                if (.not. cc_axis_defined) then
                    call mpp_error(FATAL, &
                                   "fms_io_unstructured_save_restart:" &
                                   //" three-dimensional field " &
                                   //trim(cur_var%name)//" requires that" &
                                   //"a cc-axis must be registered to the" &
                                   //" restart file "//trim(fileObj%name))
                endif
                num_var_axes = 3
                var_axes(1) = unstructured_axis
                var_axes(2) = z_axis
                var_axes(3) = cc_axis
                if (cur_var%siz(4) .eq. fileObj%max_ntime) then
                    num_var_axes = 4
                    var_axes(4) = t_axis
                endif
            else
               !For dimension and higher fields are currently not allowed.
                call mpp_error(FATAL, &
                               "fms_io_unstructured_save_restart): " &
                               //trim(cur_var%name)//" in file " &
                               //trim(fileObj%name)//" has more than three" &
                               //" dimensions (not including time level)")
            endif

           !Get the "pack size" for default real types, where
           !pack_size = (Number of bits in a real(8))/(Number of bits in a real).
            cpack = pack_size

           !For each time level, calculate a check-sum of the field data.
           !Fields with integer(4) data are handled differently then real
           !fields.  To signify an integer(4) field, set cpack = 0.
            allocate(check_val(max(1,cur_var%siz(4))))
            do k = 1,cur_var%siz(4)
                if (associated(fileObj%p0dr(k,j)%p)) then
                    check_val(k) = mpp_chksum(fileObj%p0dr(k,j)%p, &
                                              (/mpp_pe()/), &
                                              mask_val=cur_var%default_data)
                elseif (associated(fileObj%p1dr(k,j)%p)) then
                    check_val(k) = mpp_chksum(fileObj%p1dr(k,j)%p, &
                                              mask_val=cur_var%default_data)
                elseif (associated(fileObj%p2dr(k,j)%p)) then
                    check_val(k) = mpp_chksum(fileObj%p2dr(k,j)%p, &
                                              mask_val=cur_var%default_data)
                elseif (associated(fileObj%p3dr(k,j)%p)) then
                    check_val(k) = mpp_chksum(fileObj%p3dr(k,j)%p, &
                                              mask_val=cur_var%default_data)
                elseif (associated(fileObj%p0di(k,j)%p)) then
                    check_val(k) = int(fileObj%p0di(k,j)%p,kind=LONG_KIND)
                    cpack = 0
                elseif (associated(fileObj%p1di(k,j)%p)) then
                    check_val(k) = mpp_chksum(fileObj%p1di(k,j)%p, &
                                              mask_val=cur_var%default_data)
                    cpack = 0
                elseif (associated(fileObj%p2di(k,j)%p)) then
                    check_val(k) = mpp_chksum(fileObj%p2di(k,j)%p, &
                                              mask_val=cur_var%default_data)
                    cpack = 0
                elseif (associated(fileObj%p3di(k,j)%p)) then
                    call mpp_error(FATAL, &
                                   "fms_io_unstructured_save_restart:" &
                                   //" 3D integer restart fields are not" &
                                   //" currently supported. (" &
                                   //trim(cur_var%name)//" of file " &
                                   //trim(fileObj%name)//")")
                else
                    call mpp_error(FATAL, &
                                   "fms_io_unstructured_save_restart:" &
                                   //" there is no pointer associated with " &
                                   //" the data of field " &
                                   //trim(cur_var%name)//" of file " &
                                   //trim(fileObj%name))
                endif
            enddo

           !Write out the metadata from a field.  Check-sums are only written
           !if field data is written to the restart file.
            if (write_field_data) then ! Write checksums only if valid field data exists
                call mpp_write_meta(funit, &
                                    cur_var%field, &
                                    var_axes(1:num_var_axes), &
                                    cur_var%name, &
                                    cur_var%units, &
                                    cur_var%longname, &
                                    pack=cpack, &
                                    checksum=check_val, &
                                    fill=cur_var%default_data)
            else
                call mpp_write_meta(funit, &
                                    cur_var%field, &
                                    var_axes(1:num_var_axes), &
                                    cur_var%name, &
                                    cur_var%units, &
                                    cur_var%longname, &
                                    pack=cpack, &
                                    fill=cur_var%default_data)
            endif
            deallocate(check_val)
            cur_var => null()
        enddo

       !Write the axis data to the restart file.
        call mpp_write(funit, &
                       x_axis)
        call mpp_write(funit, &
                       y_axis)
        call mpp_write(funit, &
                       unstructured_axis)
        if (z_axis_defined) then
            call mpp_write(funit, &
                           z_axis)
        endif
        if (cc_axis_defined) then
            call mpp_write(funit, &
                           cc_axis)
        endif
    endif

   !Write out field data to the restart file.
    if (write_field_data) then

       !Loop through all time levels.
        do k = 1,fileObj%max_ntime

           !Get the time value for the time level.
            if (present(time_level)) then
                tlev = time_level
            else
                tlev = real(k)
            endif

           !Loop through the fields.
            do j = 1,fileObj%nvar

               !Point to the current field.
                cur_var => fileObj%var(j)

               !Cycle to the next field if the current field is read only.
                if (cur_var%read_only) then
                    cur_var => null()
                    cycle
                endif

               !Write out the field data to the file.
                if (k .le. cur_var%siz(4)) then
                    if (associated(fileObj%p0dr(k,j)%p)) then
                        call mpp_write(funit, &
                                       cur_var%field, &
                                       fileObj%p0dr(k,j)%p, &
                                       tlev)
                    elseif (associated(fileObj%p1dr(k,j)%p)) then
                        call mpp_io_unstructured_write(funit, &
                                                       cur_var%field, &
                                                       domain, &
                                                       fileObj%p1dr(k,j)%p, &
                                                       nelems_io_in=fileObj%axes(UNSTRUCTURED_GRID_INDEX)%nelems, &
                                                       tstamp=tlev, &
                                                       default_data=cur_var%default_data)
                    elseif (associated(fileObj%p2dr(k,j)%p)) then
                        if (.not. z_axis_defined) then
                            call mpp_error(FATAL, &
                                           "fms_io_unstructured_save_restart:" &
                                           //" two-dimensional field " &
                                           //trim(cur_var%name)//" requires that" &
                                           //"a z-axis must be registered to the" &
                                           //" restart file "//trim(fileObj%name))
                        endif
                        call mpp_io_unstructured_write(funit, &
                                                       cur_var%field, &
                                                       domain, &
                                                       fileObj%p2dr(k,j)%p, &
                                                       nelems_io_in=fileObj%axes(UNSTRUCTURED_GRID_INDEX)%nelems, &
                                                       tstamp=tlev, &
                                                       default_data=cur_var%default_data)
                    elseif (associated(fileObj%p3dr(k,j)%p)) then
                        if (.not. z_axis_defined) then
                            call mpp_error(FATAL, &
                                           "fms_io_unstructured_save_restart:" &
                                           //" three-dimensional field " &
                                           //trim(cur_var%name)//" requires that" &
                                           //"a z-axis must be registered to the" &
                                           //" restart file "//trim(fileObj%name))
                        endif
                        if (.not. cc_axis_defined) then
                            call mpp_error(FATAL, &
                                           "fms_io_unstructured_save_restart:" &
                                           //" three-dimensional field " &
                                           //trim(cur_var%name)//" requires that" &
                                           //"a cc-axis must be registered to the" &
                                           //" restart file "//trim(fileObj%name))
                        endif
                        call mpp_io_unstructured_write(funit, &
                                                       cur_var%field, &
                                                       domain, &
                                                       fileObj%p3dr(k,j)%p, &
                                                       nelems_io_in=fileObj%axes(UNSTRUCTURED_GRID_INDEX)%nelems, &
                                                       tstamp=tlev, &
                                                       default_data=cur_var%default_data)
                    elseif (associated(fileObj%p0di(k,j)%p)) then
                        r0d =  real(fileObj%p0di(k,j)%p)
                        call mpp_write(funit, &
                                       cur_var%field, &
                                       r0d, &
                                       tlev)
                    elseif (associated(fileObj%p1di(k,j)%p)) then
                        allocate(r1d(cur_var%siz(1)))
                        r1d = real(fileObj%p1di(k,j)%p)
                        call mpp_io_unstructured_write(funit, &
                                                       cur_var%field, &
                                                       domain, &
                                                       r1d, &
                                                       nelems_io_in=fileObj%axes(UNSTRUCTURED_GRID_INDEX)%nelems, &
                                                       tstamp=tlev, &
                                                       default_data=cur_var%default_data)
                        deallocate(r1d)
                    elseif (associated(fileObj%p2di(k,j)%p)) then
                        if (.not. z_axis_defined) then
                            call mpp_error(FATAL, &
                                           "fms_io_unstructured_save_restart:" &
                                           //" two-dimensional field " &
                                           //trim(cur_var%name)//" requires that" &
                                           //"a z-axis must be registered to the" &
                                           //" restart file "//trim(fileObj%name))
                        endif
                        allocate(r2d(cur_var%siz(1),cur_var%siz(2)))
                        r2d = real(fileObj%p2di(k,j)%p)
                        call mpp_io_unstructured_write(funit, &
                                                       cur_var%field, &
                                                       domain, &
                                                       r2d, &
                                                       nelems_io_in=fileObj%axes(UNSTRUCTURED_GRID_INDEX)%nelems, &
                                                       tstamp=tlev, &
                                                       default_data=cur_var%default_data)
                        deallocate(r2d)
                    else
                        call mpp_error(FATAL, &
                                       "fms_io_unstructured_save_restart:" &
                                       //" there is no pointer associated" &
                                       //" with the data of field " &
                                       //trim(cur_var%name)//" of file " &
                                       //trim(fileObj%name))
                    endif
                endif
            enddo
        enddo
        cur_var => null()
    endif

   !Close the restart file.
    call mpp_close(funit)

   !Nullify local pointers.
    domain => null()
    axis => null()
    cur_var => null()

    return
end subroutine fms_io_unstructured_save_restart

!----------
