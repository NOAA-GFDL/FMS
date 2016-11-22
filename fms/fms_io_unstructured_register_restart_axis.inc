!----------
!ug support

!>Store an unstructured axis in a restart object.
subroutine fms_io_unstructured_register_restart_axis(fileObj, &
                                                     filename, &
                                                     fieldname, &
                                                     fdata, &
                                                     domain, &
                                                     units, &
                                                     longname)
   !Inputs/Outputs
    type(restart_file_type),intent(inout)            :: fileObj   !<A restart object.
    character(len=*),intent(in)                      :: filename  !<A name of a file.
    character(len=*),intent(in)                      :: fieldname !<A name for the axis field.
    integer(INT_KIND),dimension(:),intent(in),target :: fdata     !<Data for the axis.
    type(domainUG),intent(in),target                 :: domain    !<An unustructured mpp domain.
    character(len=*),intent(in),optional             :: units     !<Units for the axis.
    character(len=*),intent(in),optional             :: longname  !<A more descriptive name for the axis.

   !Local variables
    integer(INT_KIND)                          :: ug_index              !<Index of the unstructured axis in the fileObj%axes array.
    type(domainUG),pointer                     :: io_domain             !<Pointer to an unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes        !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist                !<A pelist.
    integer(INT_KIND)                          :: input_filename_length !<The length of the trimmed input filename.
    character(len=256)                         :: tmp_filename          !<A character buffer used to store various file names.
    character(len=256)                         :: filename_suffix       !<A string appended to the end of the inputted file name.
    character(len=256)                         :: mosaic_filename       !<The filename returned by the get_mosaic_tile_file_ug routine.

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis:" &
                       //" you must first to call fms_io_init")
    endif

   !All axes must be registered before any fields.  Make sure that no
   !fields have been registered to the restart object.
    if (associated(fileObj%var)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis: " &
                       //" you cannot register any fields before an axis.")
    endif

   !Use this code to make the filename consistent with the
   !fms_io_unstructured_setup_one_field routine.
    input_filename_length = len_trim(filename)
    if (input_filename_length .gt. 128) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis:" &
                       //" the inputted file name is longer than 128" &
                       //" characters.")
    endif
    if (filename(input_filename_length-2:input_filename_length) .eq. ".nc") then
        tmp_filename = filename(1:input_filename_length-3)
    else
        tmp_filename = filename(1:input_filename_length)
    endif
    filename_suffix = ""
    if (len_trim(filename_appendix) .gt. 0) then
        filename_suffix = filename_appendix
    endif
    if (len_trim(filename_suffix) .gt. 0) then
        tmp_filename = trim(tmp_filename)//"."//trim(filename_suffix)
    endif
    call get_mosaic_tile_file_ug(tmp_filename, &
                                 mosaic_filename, &
                                 domain)

   !Make sure that the correct file name was passed in, or set the filename
   !if this is the first axis/field registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        fileObj%name = trim(mosaic_filename)
    else
        if (trim(mosaic_filename) .ne. trim(fileObj%name)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_axis:" &
                           //" the inputted file name does not match the" &
                           //" existing file name for this restart object.")
        endif
    endif

   !If this is the first axis registered for the restart object, then
   !allocate the fileObj%axes array.  The size of the fileObj%axes array
   !is determined by the NIDX module parameter.
    if (.not. allocated(fileObj%axes)) then
        allocate(fileObj%axes(NIDX))
    endif

   !The unstructured axis is specified by the index UNSTRUCTURED_GRID_INDEX
   !in the fileObj%axes array.  Make sure that data has not already been
   !placed there.
    ug_index = UNSTRUCTURED_GRID_INDEX
    if (allocated(fileObj%axes(ug_index)%unstructured_axis_data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis: " &
                       //" the unstructured axis for this restart object" &
                       //" has already been defined.")
    endif

   !Set the name of the axis.
    fileObj%axes(ug_index)%name = trim(fieldname)

   !Point to the inputted unstructured domain.
    fileObj%axes(ug_index)%domain_ug => domain

   !Initialize the number of unstructured data elements each rank in
   !an unstructured I/O domain is responsible for.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(fileObj%axes(ug_index)%nelems(io_domain_npes))
    fileObj%axes(ug_index)%nelems = 0

   !Gather the sizes of the inputted unstructured data arrays for each rank
   !onto the root rank of the I/O domain pelist.
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)
    call mpp_gather((/size(fdata)/), &
                    fileObj%axes(ug_index)%nelems, &
                    pelist)

   !Gather the inputted unstructured data from each rank onto the root
   !rank of the I/O domain pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(fileObj%axes(ug_index)%unstructured_axis_data(sum(fileObj%axes(ug_index)%nelems)))
    else
       !This array for a non-root rank on the I/O domain pelist should never
       !be used, but is allocated to signify that this axis is defined for
       !this restart object.
        allocate(fileObj%axes(ug_index)%unstructured_axis_data(1))
        fileObj%axes(ug_index)%unstructured_axis_data = 0
    endif
    call mpp_gather(fdata, &
                    size(fdata), &
                    fileObj%axes(ug_index)%unstructured_axis_data, &
                    fileObj%axes(ug_index)%nelems, &
                    pelist)

   !Nullify local pointers and deallocate local allocatables.
    io_domain => null()
    deallocate(pelist)

   !Set the units and longname for the axis.
    if (present(units)) then
        fileObj%axes(ug_index)%units = units
    endif
    if (present(longname)) then
        fileObj%axes(ug_index)%longname = longname
    endif

   !Set the "cartesian" string for the axis.
    fileObj%axes(ug_index)%cartesian = "U"

    return
end subroutine fms_io_unstructured_register_restart_axis

!----------
