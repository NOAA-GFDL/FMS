!----------
!ug support

!------------------------------------------------------------------------------
!>Add a real scalar field to a restart object (restart_file_type).  Return
!!the index of the inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_r_0d(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata_0d, &
                                                         domain, &
                                                         mandatory, &
                                                         data_default, &
                                                         longname, &
                                                         units, &
                                                         read_only, &
                                                         restart_owns_data) &
                                                         result(restart_index)

   !Inputs/Outputs
    type(restart_file_type),intent(inout) :: fileObj           !<A restart object.
    character(len=*),intent(in)           :: filename          !<The name of a file.
    character(len=*),intent(in)           :: fieldname         !<The name of a field.
    real,intent(in),target                :: fdata_0d          !<Some data.
    type(domainUG),intent(in),target      :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional           :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional              :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional  :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional  :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                     :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND) :: index_field !<Index of the inputted field in the fileObj%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_0d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that an unstructured axis has been registered to the restart
   !object.  This is needed to determine the name of the restart file.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_0d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    else
        if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_0d:" &
                           //" no unstructured axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             (/1,1,1,1/), &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             is_scalar_field=.true., &
                                             data_default=data_default, &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p0dr(fileObj%var(index_field)%siz(4),index_field)%p => fdata_0d
    fileObj%var(index_field)%ndim = 0
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_r_0d

!------------------------------------------------------------------------------
!>Add a real 1D field to a restart object (restart_file_type), where the
!!field is assumed to be along the unstructured axis.  Return
!!the index of the inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_r_1d(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata_1d, &
                                                         domain, &
                                                         mandatory, &
                                                         data_default, &
                                                         longname, &
                                                         units, &
                                                         read_only, &
                                                         restart_owns_data) &
                                                         result(restart_index)

   !Inputs/Outputs
    type(restart_file_type),intent(inout) :: fileObj           !<A restart object.
    character(len=*),intent(in)           :: filename          !<The name of a file.
    character(len=*),intent(in)           :: fieldname         !<The name of a field.
    real,dimension(:),intent(in),target   :: fdata_1d          !<Some data.
    type(domainUG),intent(in),target      :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional           :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional              :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional  :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional  :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                     :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND) :: index_field !<Index of the inputted field in the fileObj%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_1d:" &
                       //" you must first call fms_io_init")
    endif

   !One-dimensional fields (foo) are assumed to be functions of the
   !unstructured axis (foo = foo(unstructured)).  Make sure that an
   !unstructured axis has been registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_1d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    else
        if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_1d:" &
                           //" no unstructured axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
    endif

   !Make sure that the inputted data is the same size as the registered axis.
    if (size(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data) &
        .ne. size(fdata_1d)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_1d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered unstructured axis.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             (/size(fdata_1d,1),1,1,1/), &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             is_scalar_field=.false., &
                                             data_default=data_default, &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p1dr(fileObj%var(index_field)%siz(4),index_field)%p => fdata_1d
    fileObj%var(index_field)%ndim = 1
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_r_1d

!------------------------------------------------------------------------------
!>Add a real 2D field to a restart object (restart_file_type), where the
!!field's 1st axis assumed to be along the unstructured axis and the field's
!!2nd axis is assumed to be along the z-axis.  Return the index of the
!!inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_r_2d(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata_2d, &
                                                         domain, &
                                                         mandatory, &
                                                         data_default, &
                                                         longname, &
                                                         units, &
                                                         read_only, &
                                                         restart_owns_data) &
                                                         result(restart_index)

   !Inputs/Outputs
    type(restart_file_type),intent(inout) :: fileObj           !<A restart object.
    character(len=*),intent(in)           :: filename          !<The name of a file.
    character(len=*),intent(in)           :: fieldname         !<The name of a field.
    real,dimension(:,:),intent(in),target :: fdata_2d          !<Some data.
    type(domainUG),intent(in),target      :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional           :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional              :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional  :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional  :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                     :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND) :: index_field !<Index of the inputted field in the fileObj%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" you must first call fms_io_init")
    endif

   !Two-dimensional fields (foo) are assumed to be functions of the
   !unstructured axis and the z-axis (foo = foo(unstructured,z)).  Make sure
   !that an unstructured axis and a z-axis have been registered to the
   !restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    else
        if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" no unstructured axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
        if (.not. associated(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" no z-axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
    endif

   !Make sure that the inputted data is the same size as the registered axes.
    if (size(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data) &
        .ne. size(fdata_2d,1)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered unstructured axis.")
    endif
    if (size(fileObj%axes(ZIDX)%data) .ne. size(fdata_2d,2)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered z-axis.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             (/size(fdata_2d,1),size(fdata_2d,2),1,1/), &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             is_scalar_field=.false., &
                                             data_default=data_default, &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p2dr(fileObj%var(index_field)%siz(4),index_field)%p => fdata_2d
    fileObj%var(index_field)%ndim = 2
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_r_2d

!------------------------------------------------------------------------------
!>Add a real 3D field to a restart object (restart_file_type), where the
!!field's 1st axis assumed to be along the unstructured axis, the fields's
!!second axis is assumed to be along the z-axis, and the field's third axis
!!is assumed to be along the cc-axis (???).  Return the index of the
!!inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_r_3d(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata_3d, &
                                                         domain, &
                                                         mandatory, &
                                                         data_default, &
                                                         longname, &
                                                         units, &
                                                         read_only, &
                                                         restart_owns_data) &
                                                         result(restart_index)

   !Inputs/Outputs
    type(restart_file_type),intent(inout)   :: fileObj           !<A restart object.
    character(len=*),intent(in)             :: filename          !<The name of a file.
    character(len=*),intent(in)             :: fieldname         !<The name of a field.
    real,dimension(:,:,:),intent(in),target :: fdata_3d          !<Some data.
    type(domainUG),intent(in),target        :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional             :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional                :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional    :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional    :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional   :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional   :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                       :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND) :: index_field !<Index of the inputted field in the fileObj%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" you must first call fms_io_init")
    endif

   !Three-dimensional fields (foo) are assumed to be functions of the
   !unstructured axis, z-axis, and cc-axis (foo = foo(unstructured,z,cc)).
   !Make sure that an unstructured axis, z-axis, and cc-axis have been
   !registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    else
        if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" no unstructured axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
        if (.not. associated(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" no z-axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
        if (.not. associated(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" no cc-axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
    endif

   !Make sure that the inputted data is the same size as the registered axes.
    if (size(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data) &
        .ne. size(fdata_3d,1)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered unstructured axis.")
    endif
    if (size(fileObj%axes(ZIDX)%data) .ne. size(fdata_3d,2)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered z-axis.")
    endif
    if (size(fileObj%axes(CCIDX)%data) .ne. size(fdata_3d,3)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered cc-axis.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             (/size(fdata_3d,1),size(fdata_3d,2),size(fdata_3d,3),1/), &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             is_scalar_field=.false., &
                                             data_default=data_default, &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p3dr(fileObj%var(index_field)%siz(4),index_field)%p => fdata_3d
    fileObj%var(index_field)%ndim = 3
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_r_3d

!------------------------------------------------------------------------------
!>Add an integer scalar field to a restart object (restart_file_type).  Return
!!the index of the inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_i_0d(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata_0d, &
                                                         domain, &
                                                         mandatory, &
                                                         data_default, &
                                                         longname, &
                                                         units, &
                                                         read_only, &
                                                         restart_owns_data) &
                                                         result(restart_index)

   !Inputs/Outputs
    type(restart_file_type),intent(inout) :: fileObj           !<A restart object.
    character(len=*),intent(in)           :: filename          !<The name of a file.
    character(len=*),intent(in)           :: fieldname         !<The name of a field.
    integer,intent(in),target             :: fdata_0d          !<Some data.
    type(domainUG),intent(in),target      :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional           :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional              :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional  :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional  :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                     :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND) :: index_field !<Index of the inputted field in the fileObj%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_0d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that an unstructured axis has been registered to the restart
   !object.  This is needed to determine the name of the restart file.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_0d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    else
        if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_0d:" &
                           //" no unstructured axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             (/1,1,1,1/), &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             is_scalar_field=.true., &
                                             data_default=data_default, &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p0di(fileObj%var(index_field)%siz(4),index_field)%p => fdata_0d
    fileObj%var(index_field)%ndim = 0
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_i_0d

!------------------------------------------------------------------------------
!>Add an integer 1D field to a restart object (restart_file_type), where the
!!field is assumed to be along the unstructured axis.  Return
!!the index of the inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_i_1d(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata_1d, &
                                                         domain, &
                                                         mandatory, &
                                                         data_default, &
                                                         longname, &
                                                         units, &
                                                         read_only, &
                                                         restart_owns_data) &
                                                         result(restart_index)

   !Inputs/Outputs
    type(restart_file_type),intent(inout)  :: fileObj           !<A restart object.
    character(len=*),intent(in)            :: filename          !<The name of a file.
    character(len=*),intent(in)            :: fieldname         !<The name of a field.
    integer,dimension(:),intent(in),target :: fdata_1d          !<Some data.
    type(domainUG),intent(in),target       :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional            :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional               :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional   :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional   :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional  :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional  :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                      :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND) :: index_field !<Index of the inputted field in the fileObj%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_1d:" &
                       //" you must first call fms_io_init")
    endif

   !One-dimensional fields (foo) are assumed to be functions of the
   !unstructured axis (foo = foo(unstructured)).  Make sure that an
   !unstructured axis has been registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_1d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    else
        if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_1d:" &
                           //" no unstructured axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
    endif

   !Make sure that the inputted data is the same size as the registered axis.
    if (size(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data) &
        .ne. size(fdata_1d)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_1d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered unstructured axis.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             (/size(fdata_1d,1),1,1,1/), &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             is_scalar_field=.false., &
                                             data_default=data_default, &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p1di(fileObj%var(index_field)%siz(4),index_field)%p => fdata_1d
    fileObj%var(index_field)%ndim = 1
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_i_1d

!------------------------------------------------------------------------------
!>Add an integer 2D field to a restart object (restart_file_type), where the
!!field's 1st axis assumed to be along the unstructured axis and the field's
!!2nd axis is assumed to be along the z-axis.  Return the index of the
!!inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_i_2d(fileObj, &
                                                         filename, &
                                                         fieldname, &
                                                         fdata_2d, &
                                                         domain, &
                                                         mandatory, &
                                                         data_default, &
                                                         longname, &
                                                         units, &
                                                         read_only, &
                                                         restart_owns_data) &
                                                         result(restart_index)

   !Inputs/Outputs
    type(restart_file_type),intent(inout)    :: fileObj           !<A restart object.
    character(len=*),intent(in)              :: filename          !<The name of a file.
    character(len=*),intent(in)              :: fieldname         !<The name of a field.
    integer,dimension(:,:),intent(in),target :: fdata_2d          !<Some data.
    type(domainUG),intent(in),target         :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional              :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional                 :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional     :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional     :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional    :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional    :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                        :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND) :: index_field !<Index of the inputted field in the fileObj%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" you must first call fms_io_init")
    endif

   !Two-dimensional fields (foo) are assumed to be functions of the
   !unstructured axis and the z-axis (foo = foo(unstructured,z)).  Make sure
   !that an unstructured axis and a z-axis have been registered to the
   !restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    else
        if (.not. allocated(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" no unstructured axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
        if (.not. associated(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" no z-axis has been registered for" &
                           //" the restart object "//trim(fileObj%name))
        endif
    endif

   !Make sure that the inputted data is the same size as the registered axes.
    if (size(fileObj%axes(UNSTRUCTURED_GRID_INDEX)%unstructured_axis_data) &
        .ne. size(fdata_2d,1)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered unstructured axis.")
    endif
    if (size(fileObj%axes(ZIDX)%data) .ne. size(fdata_2d,2)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" the inputted data size does not match the size" &
                       //" of the registered z-axis.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             (/size(fdata_2d,1),size(fdata_2d,2),1,1/), &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             is_scalar_field=.false., &
                                             data_default=data_default, &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p2di(fileObj%var(index_field)%siz(4),index_field)%p => fdata_2d
    fileObj%var(index_field)%ndim = 2
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_i_2d

!------------------------------------------------------------------------------

!----------
