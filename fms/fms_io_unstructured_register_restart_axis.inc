!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* FMS is free software: you can redistribute it and/or modify it under
!* the terms of the GNU Lesser General Public License as published by
!* the Free Software Foundation, either version 3 of the License, or (at
!* your option) any later version.
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!* for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
!----------
!ug support

!------------------------------------------------------------------------------
!>Store a real axis (x,y,z,...) in a restart object assoicated with an
!!unstructured mpp domain.
subroutine fms_io_unstructured_register_restart_axis_r1D(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata, &
                                                         cartesian, &
                                                         domain, &
                                                         units, &
                                                         longname, &
                                                         sense, &
                                                         fmin, &
                                                         calendar)
   !Inputs/Outputs
    type(restart_file_type),intent(inout) :: fileObj   !<A restart object.
    character(len=*),intent(in)           :: filename  !<A name of a file.
    character(len=*),intent(in)           :: fieldname !<A name for the axis field.
    real,dimension(:),intent(in),target   :: fdata     !<Data for the axis.
    character(len=*),intent(in)           :: cartesian !<String indicating which cartesian axis this is (i.e. X, Y, Z).
    type(domainUG),intent(in),target      :: domain    !<An unustructured mpp domain.
    character(len=*),intent(in),optional  :: units     !<Units for the axis.
    character(len=*),intent(in),optional  :: longname  !<A more descriptive name for the axis.
    integer(INT_KIND),intent(in),optional :: sense     !<Positive direction.
    real,intent(in),optional              :: fmin      !<Minimum value for this axis.
    character(len=*),intent(in),optional  :: calendar  !<Type of calendar (only for time axis.)

   !Local variables
    integer(INT_KIND)                          :: input_filename_length !<The length of the trimmed input filename.
    character(len=256)                         :: tmp_filename          !<A character buffer used to store various file names.
    character(len=256)                         :: filename_suffix       !<A string appended to the end of the inputted file name.
    character(len=256)                         :: mosaic_filename       !<The filename returned by the get_mosaic_tile_file_ug routine.
    integer(INT_KIND)                          :: axis_index            !<Index of the inputted axis in the fileObj%axes array.
    type(domainUG),pointer                     :: io_domain             !<Pointer to an unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes        !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist                !<A pelist.
    integer(INT_KIND),dimension(:),allocatable :: fdata_sizes           !<Size of the axis data for each rank in the I/O domain pelist.

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_r1D:" &
                       //" you must first to call fms_io_init.")
    endif

   !All axes must be registered before any fields.  Make sure that no
   !fields have been registered to the restart object.
    if (associated(fileObj%var)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_r1D: " &
                       //" you cannot register any fields before an axis.")
    endif

   !Use this code to make the filename consistent with the
   !fms_io_unstructured_setup_one_field routine.
    input_filename_length = len_trim(filename)
    if (input_filename_length .gt. 128) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_r1D:" &
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
        filename_suffix = trim(filename_appendix)
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
                           "fms_io_unstructured_register_restart_axis_r1D:" &
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

   !Determine the index of the inputted axis in the fileObj%axes array from
   !the inputted cartesian string.
    select case (trim(cartesian))
        case ("X")
            axis_index = XIDX
        case ("Y")
            axis_index = YIDX
        case ("Z")
            axis_index = ZIDX
        case ("T")
            axis_index = TIDX
        case ("CC")
            axis_index = CCIDX
        case default
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_axis_r1D:" &
                           //" an invalid cartesian string was passed in.")
    end select

   !Make sure that data has not already been registered for the inputted
   !axis.
    if (associated(fileObj%axes(axis_index)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_r1D:" &
                       //" the "//trim(cartesian)//" axis for this restart" &
                       //" object has already been defined.")
    endif

   !Make sure that the axis size is consistent for all ranks on the
   !unstructured I/O domain pelist.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)
    allocate(fdata_sizes(io_domain_npes))
    fdata_sizes = 0
    call mpp_gather((/size(fdata)/), &
                    fdata_sizes, &
                    pelist)
    if (mpp_pe() .eq. pelist(1)) then
        if (maxval(fdata_sizes) .ne. size(fdata) .or. &
            minval(fdata_sizes) .ne. size(fdata)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_axis_r1D:" &
                           //" the "//trim(cartesian)//" axis must be the" &
                           //" the same size for all ranks in the" &
                           //" unstructured I/O domain pelist.")
        endif
    endif
    io_domain => null()
    deallocate(pelist)
    deallocate(fdata_sizes)

   !Set the name of the axis.
    fileObj%axes(axis_index)%name = trim(fieldname)

   !Point to the inputted unstructured domain for the axis.
    fileObj%axes(axis_index)%domain_ug => domain

   !Point to the inputted axis data.
    fileObj%axes(axis_index)%data => fdata

   !Store the inputted cartesian string. (Why?)
    fileObj%axes(axis_index)%cartesian = trim(cartesian)

   !Set the dimension length for the axis to -1 to signify that this is
   !not a "compressed" axis.
    fileObj%axes(axis_index)%dimlen = -1

   !Store the units for the axis.
    if (present(units)) then
        fileObj%axes(axis_index)%units = trim(units)
    else
        fileObj%axes(axis_index)%units = ""
    endif

   !Store the longname for the axis.
    if (present(longname)) then
        fileObj%axes(axis_index)%longname = trim(longname)
    else
        fileObj%axes(axis_index)%longname = ""
    endif

   !Store the "sense" for the axis.  Inputs must be for the z-dimension.
    if (present(sense)) then
        if (axis_index .ne. ZIDX) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_axis_r1D:" &
                           //" sense may only be defined for the z-axis.")
        endif
        if (abs(sense) .ne. 1) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_axis_r1D:" &
                           //" sense may only have the values +/- 1")
        endif
        fileObj%axes(axis_index)%sense = sense
    else
        fileObj%axes(axis_index)%sense = 0
    endif

   !Store the minimum value allowed for the axis.
    if (present(fmin)) then
        fileObj%axes(axis_index)%min = fmin
    else
        fileObj%axes(axis_index)%min = 0
    endif

   !Store the calendar for the axis.  This is only done for the time dimension.
    if (axis_index .eq. TIDX) then
        fileObj%axes(axis_index)%calendar = trim(calendar)
    endif

    return
end subroutine fms_io_unstructured_register_restart_axis_r1D

!------------------------------------------------------------------------------
!>Store an integer "compressed" axis in a restart object assoicated with an
!!unstructured mpp domain.
subroutine fms_io_unstructured_register_restart_axis_i1D(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata, &
                                                         compressed, &
                                                         compressed_axis, &
                                                         dimlen, &
                                                         domain, &
                                                         dimlen_name, &
                                                         dimlen_lname, &
                                                         units, &
                                                         longname, &
                                                         imin)
   !Inputs/Outputs
    type(restart_file_type),intent(inout)            :: fileObj         !<A restart object.
    character(len=*),intent(in)                      :: filename        !<A name of a file.
    character(len=*),intent(in)                      :: fieldname       !<A name for the axis field.
    integer(INT_KIND),dimension(:),intent(in),target :: fdata           !<Data for the axis.
    character(len=*),intent(in)                      :: compressed      !<"Compressed" string (???)
    character(len=*),intent(in)                      :: compressed_axis !<"Compressed" axis string.
    integer(INT_KIND),intent(in)                     :: dimlen          !<Length of the compressed dimension.
    type(domainUG),intent(in),target                 :: domain          !<An unustructured mpp domain.
    character(len=*),intent(in),optional             :: dimlen_name     !<(???)
    character(len=*),intent(in),optional             :: dimlen_lname    !<(???)
    character(len=*),intent(in),optional             :: units           !<Units for the axis.
    character(len=*),intent(in),optional             :: longname        !<A more descriptive name for the axis.
    integer(INT_KIND),intent(in),optional            :: imin            !<Minium value for the dimension.

   !Local variables
    integer(INT_KIND)                          :: input_filename_length !<The length of the trimmed input filename.
    character(len=256)                         :: tmp_filename          !<A character buffer used to store various file names.
    character(len=256)                         :: filename_suffix       !<A string appended to the end of the inputted file name.
    character(len=256)                         :: mosaic_filename       !<The filename returned by the get_mosaic_tile_file_ug routine.
    integer(INT_KIND)                          :: axis_index            !<Index of the inputted axis in the fileObj%axes array.
    type(domainUG),pointer                     :: io_domain             !<Pointer to an unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes        !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist                !<A pelist.

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_i1D:" &
                       //" you must first to call fms_io_init.")
    endif

   !All axes must be registered before any fields.  Make sure that no
   !fields have been registered to the restart object.
    if (associated(fileObj%var)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_i1D:" &
                       //" you cannot register any fields before an axis.")
    endif

   !Use this code to make the filename consistent with the
   !fms_io_unstructured_setup_one_field routine.
    input_filename_length = len_trim(filename)
    if (input_filename_length .gt. 128) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_i1D:" &
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
        filename_suffix = trim(filename_appendix)
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
                           "fms_io_unstructured_register_restart_axis_i1D:" &
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

   !Get the index of the inputted axis in the fileObj%axes array from the
   !inputted compressed_axis string.
    select case (trim(compressed_axis))
        case ("C")
            axis_index = CIDX
        case ("H")
            axis_index = HIDX
        case default
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_axis_i1D:" &
                           //" invalid compressed_axis string was passed in.")
    end select

   !Make sure that data has not already been registered for the inputted axis.
    if (allocated(fileObj%axes(axis_index)%idx)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_i1D:" &
                       //" the "//trim(compressed_axis)//" axis for this" &
                       //" restart object has already been defined.")
    endif

   !Set the name of the axis.
    fileObj%axes(axis_index)%name = trim(fieldname)

   !Point to the inputted unstructured domain.
    fileObj%axes(axis_index)%domain_ug => domain

   !Initialize the number of data elements each rank in an unstructured I/O
   !domain is responsible for.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(fileObj%axes(axis_index)%nelems(io_domain_npes))
    fileObj%axes(axis_index)%nelems = 0
    fileObj%axes(axis_index)%nelems_for_current_rank = size(fdata)

   !Gather the sizes of the inputted data arrays for each rank onto the root
   !rank of the I/O domain pelist.
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)
    call mpp_gather((/size(fdata)/), &
                    fileObj%axes(axis_index)%nelems, &
                    pelist)

   !Gather the inputted data from each rank onto the root rank of the I/O
   !domain pelist.
    if (mpp_pe() .eq. pelist(1)) then
        allocate(fileObj%axes(axis_index)%idx(sum(fileObj%axes(axis_index)%nelems)))
    else
       !This array for a non-root rank on the I/O domain pelist should never
       !be used, but is allocated to signify that this axis is defined for
       !this restart object.
        allocate(fileObj%axes(axis_index)%idx(1))
        fileObj%axes(axis_index)%idx = 0
    endif
    call mpp_gather(fdata, &
                    size(fdata), &
                    fileObj%axes(axis_index)%idx, &
                    fileObj%axes(axis_index)%nelems, &
                    pelist)

   !Nullify local pointers and deallocate local allocatables.
    io_domain => null()
    deallocate(pelist)

   !Set the "compressed" string for the axis.
    fileObj%axes(axis_index)%compressed = trim(compressed)

   !Set the dimension length for the axis.
    fileObj%axes(axis_index)%dimlen = dimlen

   !Set the dimlen_name (???) for the axis.
    if (present(dimlen_name)) then
        fileObj%axes(axis_index)%dimlen_name = trim(dimlen_name)
    else
        fileObj%axes(axis_index)%dimlen_name = ""
    endif

   !Set the dimlen_lname (???) for the axis.
    if (present(dimlen_lname)) then
        fileObj%axes(axis_index)%dimlen_lname = trim(dimlen_lname)
    else
        fileObj%axes(axis_index)%dimlen_lname = ""
    endif

   !Set the units for the axis.
    if (present(units)) then
        fileObj%axes(axis_index)%units = trim(units)
    else
        fileObj%axes(axis_index)%units = ""
    endif

   !Set the longname for the axis.
    if (present(longname)) then
        fileObj%axes(axis_index)%longname = trim(longname)
    else
        fileObj%axes(axis_index)%longname = ""
    endif

   !Set the minimum value for the axis.
    if (present(imin)) then
        fileObj%axes(axis_index)%imin = imin
    else
        fileObj%axes(axis_index)%imin = 0
    endif

    return
end subroutine fms_io_unstructured_register_restart_axis_i1D

!------------------------------------------------------------------------------
!>Store an unlimited axis in a restart object assoicated with an unstructured
!!mpp domain.
subroutine fms_io_unstructured_register_restart_axis_u(fileObj, &
                                                       filename, &
                                                       fieldname, &
                                                       nelems, &
                                                       domain, &
                                                       units, &
                                                       longname)
   !Inputs/Outputs
    type(restart_file_type),intent(inout) :: fileObj   !<A restart object.
    character(len=*),intent(in)           :: filename  !<A name of a file.
    character(len=*),intent(in)           :: fieldname !<A name for the axis field.
    integer(INT_KIND),intent(in)          :: nelems    !<Number of elements on the axis for the current rank.
    type(domainUG),intent(in),target      :: domain    !<An unustructured mpp domain.
    character(len=*),intent(in),optional  :: units     !<Units for the axis.
    character(len=*),intent(in),optional  :: longname  !<A more descriptive name for the axis.

   !Local variables
    integer(INT_KIND)                          :: input_filename_length !<The length of the trimmed input filename.
    character(len=256)                         :: tmp_filename          !<A character buffer used to store various file names.
    character(len=256)                         :: filename_suffix       !<A string appended to the end of the inputted file name.
    character(len=256)                         :: mosaic_filename       !<The filename returned by the get_mosaic_tile_file_ug routine.
    integer(INT_KIND)                          :: axis_index            !<Index of the inputted axis in the fileObj%axes array.
    type(domainUG),pointer                     :: io_domain             !<Pointer to an unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes        !<The total number of ranks in an I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist                !<A pelist.

   !Make sure that the module is initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_u:" &
                       //" you must first to call fms_io_init.")
    endif

   !All axes must be registered before any fields.  Make sure that no
   !fields have been registered to the restart object.
    if (associated(fileObj%var)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_u:" &
                       //" you cannot register any fields before an axis.")
    endif

   !Use this code to make the filename consistent with the
   !fms_io_unstructured_setup_one_field routine.
    input_filename_length = len_trim(filename)
    if (input_filename_length .gt. 128) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_u:" &
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
        filename_suffix = trim(filename_appendix)
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
                           "fms_io_unstructured_register_restart_axis_u:" &
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

   !Get the index of the inputted axis in the fileObj%axes array.
    axis_index = UIDX

   !Make sure that data has not already been registered for the inputted axis.
    if (allocated(fileObj%axes(axis_index)%idx)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_axis_u:" &
                       //" the unlimited axis for this restart object" &
                       //" has already been defined.")
    endif

   !Set the name of the axis.
    fileObj%axes(axis_index)%name = trim(fieldname)

   !Point to the inputted unstructured domain.
    fileObj%axes(axis_index)%domain_ug => domain

   !Initialize the number of data elements each rank in an unstructured I/O
   !domain is responsible for.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(fileObj%axes(axis_index)%nelems(io_domain_npes))
    fileObj%axes(axis_index)%nelems = 0

   !Gather the inputted number of elements each rank is responsible for onto
   !the root rank of the I/O domain pelist.
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                  pelist)
    call mpp_gather((/nelems/), &
                    fileObj%axes(axis_index)%nelems, &
                    pelist)

   !Nullify local pointers and deallocate local allocatables.
    io_domain => null()
    deallocate(pelist)

   !Set the units for the axis.
    if (present(units)) then
        fileObj%axes(axis_index)%units = trim(units)
    else
        fileObj%axes(axis_index)%units = ""
    endif

   !Set the longname for the axis.
    if (present(longname)) then
        fileObj%axes(axis_index)%longname = trim(longname)
    else
        fileObj%axes(axis_index)%longname = ""
    endif

    return
end subroutine fms_io_unstructured_register_restart_axis_u

!------------------------------------------------------------------------------

!----------
