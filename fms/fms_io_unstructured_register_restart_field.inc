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
    type(domainUG),pointer                     :: io_domain             !<Pointer to an unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes        !<The number of ranks in the unstructured I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist                !<A pelist.
    real,dimension(:),allocatable              :: fdata_per_rank        !<Array used to gather the scalar field values.
    integer(INT_KIND)                          :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX)          :: field_dimension_sizes !<Array of dimension sizes for the field.
    integer(INT_KIND),dimension(1)             :: field_dimension_order !<Array telling the ordering of the dimensions for the field.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_0d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that the value of the scalar field is same across all ranks
   !in an I/O domain pelist.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                   pelist)
    allocate(fdata_per_rank(io_domain_npes))
    fdata_per_rank = 0.0
    call mpp_gather((/fdata_0d/), &
                    fdata_per_rank, &
                    pelist)
    if (mpp_pe() .eq. pelist(1)) then
        if (maxval(fdata_per_rank) .ne. fdata_0d .or. &
            minval(fdata_per_rank) .ne. fdata_0d) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_0d:" &
                           //" the scalar field data is not consistent across" &
                           //" all ranks in the I/O domain pelist.")
        endif
    endif
    io_domain => null()
    deallocate(pelist)
    deallocate(fdata_per_rank)

   !Set the dimension sizes for the field.  These correspond to:
   !field_dimension_sizes(XIDX) = x-dimension size
   !field_dimension_sizes(YIDX) = y-dimension size
   !field_dimension_sizes(CIDX) = c-dimension size
   !field_dimension_sizes(ZIDX) = z-dimension size
   !field_dimension_sizes(HIDX) = h-dimension size
   !field_dimension_sizes(TIDX) = t-dimension size
   !field_dimension_sizes(UIDX) = u-dimension size
   !field_dimension_sizes(CCIDX) = cc-dimension size
    field_dimension_sizes = 1

   !Set the ordering of the dimensions for the field.
    field_dimension_order(1) = TIDX

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             field_dimension_order, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
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
                                                         fdata_1d_axes, &
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
    integer(INT_KIND),dimension(1)        :: fdata_1d_axes     !<An array describing the axes for the data.
    type(domainUG),intent(in),target      :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional           :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional              :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional  :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional  :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                     :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND)                 :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !<Array of dimension sizes for the field.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_1d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that at least one axis was registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_1d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    endif

   !Make sure that the first dimension of the field is a "compressed" axis,
   !and that it corresponds to an axis that has been registered to the
   !restart object.
    field_dimension_sizes = 1
    if (fdata_1d_axes(1) .eq. CIDX) then
        if (.not. allocated(fileObj%axes(CIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_1d:" &
                           //" a compressed c-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_1d,1) .ne. fileObj%axes(CIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_1d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed c-axis.")
        endif
        field_dimension_sizes(CIDX) = size(fdata_1d,1)
    elseif (fdata_1d_axes(1) .eq. HIDX) then
        if (.not. allocated(fileObj%axes(HIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_1d:" &
                           //" a compressed h-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_1d,1) .ne. fileObj%axes(HIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_1d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed h-axis.")
        endif
        field_dimension_sizes(HIDX) = size(fdata_1d,1)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_1d:" &
                       //" One dimensional fields must be compressed.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             fdata_1d_axes, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
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
                                                         fdata_2d_axes, &
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
    integer(INT_KIND),dimension(2)        :: fdata_2d_axes     !<An array describing the axes for the data.
    type(domainUG),intent(in),target      :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional           :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional              :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional  :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional  :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                     :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND)                 :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !<Array of dimension sizes for the field.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that at least one axis was registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    endif

   !Make sure that the first dimension of the field is a "compressed" axis,
   !and that it corresponds to an axis that has been registered to the
   !restart object.
    field_dimension_sizes = 1
    if (fdata_2d_axes(1) .eq. CIDX) then
        if (.not. allocated(fileObj%axes(CIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" a compressed c-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_2d,1) .ne. fileObj%axes(CIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed c-axis.")
        endif
        field_dimension_sizes(CIDX) = size(fdata_2d,1)
    elseif (fdata_2d_axes(1) .eq. HIDX) then
        if (.not. allocated(fileObj%axes(HIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" a compressed h-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_2d,1) .ne. fileObj%axes(HIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed h-axis.")
        endif
        field_dimension_sizes(HIDX) = size(fdata_2d,1)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" The first dimension of the field must be a" &
                       //" compressed dimension.")
    endif

   !Make sure that the second dimension of the inputted field corresponds to
   !either a registered z- or cc-axis.
    if (fdata_2d_axes(2) .eq. ZIDX) then
        if (.not. associated(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" a z-axis was not registered to the" &
                           //" restart object.")
        endif
        if (size(fdata_2d,2) .ne. size(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" z-axis.")
        endif
        field_dimension_sizes(ZIDX) = size(fdata_2d,2)
    elseif (fdata_2d_axes(2) .eq. CCIDX) then
        if (.not. associated(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" a cc-axis was not registered to the" &
                           //" restart object.")
        endif
        if (size(fdata_2d,2) .ne. size(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" cc-axis.")
        endif
        field_dimension_sizes(CCIDX) = size(fdata_2d,2)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_2d:" &
                       //" unsupported axis parameter for the second" &
                       //" dimension of the field.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             fdata_2d_axes, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
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
                                                         fdata_3d_axes, &
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
    integer(INT_KIND),dimension(3)          :: fdata_3d_axes     !<An array describing the axes for the data.
    type(domainUG),intent(in),target        :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional             :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional                :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional    :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional    :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional   :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional   :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                       :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND)                 :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !<Array of dimension sizes for the field.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that at least one axis was registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    endif

   !Make sure that the first dimension of the field is a "compressed" axis,
   !and that it corresponds to an axis that has been registered to the
   !restart object.
    field_dimension_sizes = 1
    if (fdata_3d_axes(1) .eq. CIDX) then
        if (.not. allocated(fileObj%axes(CIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" a compressed c-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_3d,1) .ne. fileObj%axes(CIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed c-axis.")
        endif
        field_dimension_sizes(CIDX) = size(fdata_3d,1)
    elseif (fdata_3d_axes(1) .eq. HIDX) then
        if (.not. allocated(fileObj%axes(HIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" a compressed h-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_3d,1) .ne. fileObj%axes(HIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed h-axis.")
        endif
        field_dimension_sizes(HIDX) = size(fdata_3d,1)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" The first dimension of the field must be a" &
                       //" compressed dimension.")
    endif

   !Make sure that the second and third dimensions of the inputted field
   !corresponds to some combination of registered z- and cc-axes.
    if (.not. associated(fileObj%axes(ZIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" a z-axis was not registered to the" &
                       //" restart object.")
    endif
    if (.not. associated(fileObj%axes(CCIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" a cc-axis was not registered to the" &
                       //" restart object.")
    endif
    if (fdata_3d_axes(2) .eq. ZIDX) then
        if (size(fdata_3d,2) .ne. size(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" z-axis.")
        endif
        field_dimension_sizes(ZIDX) = size(fdata_3d,2)
        if (fdata_3d_axes(3) .ne. CCIDX) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" unsupported axis parameter for the third" &
                           //" dimension of the field.")
        elseif (size(fdata_3d,3) .ne. size(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" cc-axis.")

        else
            field_dimension_sizes(CCIDX) = size(fdata_3d,3)
        endif
    elseif (fdata_3d_axes(2) .eq. CCIDX) then
        if (size(fdata_3d,2) .ne. size(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" cc-axis.")
        endif
        field_dimension_sizes(CCIDX) = size(fdata_3d,2)
        if (fdata_3d_axes(3) .ne. ZIDX) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" unsupported axis parameter for the third" &
                           //" dimension of the field.")
        elseif (size(fdata_3d,3) .ne. size(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" z-axis.")
        else
            field_dimension_sizes(ZIDX) = size(fdata_3d,3)
        endif
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r_3d:" &
                       //" unsupported axis parameter for the second" &
                       //" dimension of the field.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             fdata_3d_axes, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
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

#ifdef OVERLOAD_R8
!------------------------------------------------------------------------------
!>Add a double_kind 2D field to a restart object (restart_file_type), where the
!!field's 1st axis assumed to be along the unstructured axis and the field's
!!2nd axis is assumed to be along the z-axis.  Return the index of the
!!inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_r8_2d(fileObj, &
                                                          filename, &
                                                          fieldname, &
                                                          fdata_2d, &
                                                          fdata_2d_axes, &
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
    real(DOUBLE_KIND),dimension(:,:),intent(in),target :: fdata_2d          !<Some data.
    integer(INT_KIND),dimension(2)        :: fdata_2d_axes     !<An array describing the axes for the data.
    type(domainUG),intent(in),target      :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional           :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real(DOUBLE_KIND),intent(in),optional              :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional  :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional  :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                     :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND)                 :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !<Array of dimension sizes for the field.

   !QUICK ERROR OUT AS SUPPORT NOT YET FULLY IMPLEMENTED
    call mpp_error(FATAL, &
                  "fms_io_unstructured_register_restart_field_r8_2d:" &
                  //" support has not yet been fully implemented")

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_2d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that at least one axis was registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_2d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    endif

   !Make sure that the first dimension of the field is a "compressed" axis,
   !and that it corresponds to an axis that has been registered to the
   !restart object.
    field_dimension_sizes = 1
    if (fdata_2d_axes(1) .eq. CIDX) then
        if (.not. allocated(fileObj%axes(CIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" a compressed c-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_2d,1) .ne. fileObj%axes(CIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed c-axis.")
        endif
        field_dimension_sizes(CIDX) = size(fdata_2d,1)
    elseif (fdata_2d_axes(1) .eq. HIDX) then
        if (.not. allocated(fileObj%axes(HIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" a compressed h-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_2d,1) .ne. fileObj%axes(HIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed h-axis.")
        endif
        field_dimension_sizes(HIDX) = size(fdata_2d,1)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_2d:" &
                       //" The first dimension of the field must be a" &
                       //" compressed dimension.")
    endif

   !Make sure that the second dimension of the inputted field corresponds to
   !either a registered z- or cc-axis.
    if (fdata_2d_axes(2) .eq. ZIDX) then
        if (.not. associated(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" a z-axis was not registered to the" &
                           //" restart object.")
        endif
        if (size(fdata_2d,2) .ne. size(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" z-axis.")
        endif
        field_dimension_sizes(ZIDX) = size(fdata_2d,2)
    elseif (fdata_2d_axes(2) .eq. CCIDX) then
        if (.not. associated(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" a cc-axis was not registered to the" &
                           //" restart object.")
        endif
        if (size(fdata_2d,2) .ne. size(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" cc-axis.")
        endif
        field_dimension_sizes(CCIDX) = size(fdata_2d,2)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_2d:" &
                       //" unsupported axis parameter for the second" &
                       //" dimension of the field.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             fdata_2d_axes, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             data_default=real(data_default), &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p2dr8(fileObj%var(index_field)%siz(4),index_field)%p => fdata_2d
    fileObj%var(index_field)%ndim = 2
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_r8_2d

!------------------------------------------------------------------------------
!>Add a double_kind 3D field to a restart object (restart_file_type), where the
!!field's 1st axis assumed to be along the unstructured axis, the fields's
!!second axis is assumed to be along the z-axis, and the field's third axis
!!is assumed to be along the cc-axis (???).  Return the index of the
!!inputted field in the fileObj%var array.
function fms_io_unstructured_register_restart_field_r8_3d(fileObj, &
                                                          filename, &
                                                          fieldname, &
                                                          fdata_3d, &
                                                          fdata_3d_axes, &
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
    real(DOUBLE_KIND),dimension(:,:,:),intent(in),target :: fdata_3d          !<Some data.
    integer(INT_KIND),dimension(3)          :: fdata_3d_axes     !<An array describing the axes for the data.
    type(domainUG),intent(in),target        :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional             :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real(DOUBLE_KIND),intent(in),optional                :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional    :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional    :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional   :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional   :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                       :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND)                 :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !<Array of dimension sizes for the field.

   !QUICK ERROR OUT AS SUPPORT NOT YET FULLY IMPLEMENTED
    call mpp_error(FATAL, &
                  "fms_io_unstructured_register_restart_field_r8_3d:" &
                  //" support has not yet been fully implemented")

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_3d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that at least one axis was registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_3d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    endif

   !Make sure that the first dimension of the field is a "compressed" axis,
   !and that it corresponds to an axis that has been registered to the
   !restart object.
    field_dimension_sizes = 1
    if (fdata_3d_axes(1) .eq. CIDX) then
        if (.not. allocated(fileObj%axes(CIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" a compressed c-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_3d,1) .ne. fileObj%axes(CIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed c-axis.")
        endif
        field_dimension_sizes(CIDX) = size(fdata_3d,1)
    elseif (fdata_3d_axes(1) .eq. HIDX) then
        if (.not. allocated(fileObj%axes(HIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" a compressed h-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_3d,1) .ne. fileObj%axes(HIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed h-axis.")
        endif
        field_dimension_sizes(HIDX) = size(fdata_3d,1)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_3d:" &
                       //" The first dimension of the field must be a" &
                       //" compressed dimension.")
    endif

   !Make sure that the second and third dimensions of the inputted field
   !corresponds to some combination of registered z- and cc-axes.
    if (.not. associated(fileObj%axes(ZIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_3d:" &
                       //" a z-axis was not registered to the" &
                       //" restart object.")
    endif
    if (.not. associated(fileObj%axes(CCIDX)%data)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_3d:" &
                       //" a cc-axis was not registered to the" &
                       //" restart object.")
    endif
    if (fdata_3d_axes(2) .eq. ZIDX) then
        if (size(fdata_3d,2) .ne. size(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" z-axis.")
        endif
        field_dimension_sizes(ZIDX) = size(fdata_3d,2)
        if (fdata_3d_axes(3) .ne. CCIDX) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" unsupported axis parameter for the third" &
                           //" dimension of the field.")
        elseif (size(fdata_3d,3) .ne. size(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" cc-axis.")

        else
            field_dimension_sizes(CCIDX) = size(fdata_3d,3)
        endif
    elseif (fdata_3d_axes(2) .eq. CCIDX) then
        if (size(fdata_3d,2) .ne. size(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" cc-axis.")
        endif
        field_dimension_sizes(CCIDX) = size(fdata_3d,2)
        if (fdata_3d_axes(3) .ne. ZIDX) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" unsupported axis parameter for the third" &
                           //" dimension of the field.")
        elseif (size(fdata_3d,3) .ne. size(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_r8_3d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" z-axis.")
        else
            field_dimension_sizes(ZIDX) = size(fdata_3d,3)
        endif
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_r8_3d:" &
                       //" unsupported axis parameter for the second" &
                       //" dimension of the field.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             fdata_3d_axes, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
                                             data_default=real(data_default), &
                                             longname=longname, &
                                             units=units, &
                                             read_only=read_only, &
                                             owns_data=restart_owns_data)

   !Point to the inputted data and return the "index_field" for the field.
    fileObj%p3dr8(fileObj%var(index_field)%siz(4),index_field)%p => fdata_3d
    fileObj%var(index_field)%ndim = 3
    restart_index = index_field

    return
end function fms_io_unstructured_register_restart_field_r8_3d
#endif

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
    type(domainUG),pointer                     :: io_domain             !<Pointer to an unstructured I/O domain.
    integer(INT_KIND)                          :: io_domain_npes        !<The number of ranks in the unstructured I/O domain pelist.
    integer(INT_KIND),dimension(:),allocatable :: pelist                !<A pelist.
    integer,dimension(:),allocatable           :: fdata_per_rank        !<Array used to gather the scalar field values.
    integer(INT_KIND)                          :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX)          :: field_dimension_sizes !<Array of dimension sizes for the field.
    integer(INT_KIND),dimension(1)             :: field_dimension_order !<Array telling the ordering of the dimensions for the field.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_0d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that the value of the scalar field is same across all ranks
   !in an I/O domain pelist.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)
    io_domain_npes = mpp_get_UG_domain_npes(io_domain)
    allocate(pelist(io_domain_npes))
    call mpp_get_UG_domain_pelist(io_domain, &
                                   pelist)
    allocate(fdata_per_rank(io_domain_npes))
    fdata_per_rank = 0.0
    call mpp_gather((/fdata_0d/), &
                    fdata_per_rank, &
                    pelist)
    if (mpp_pe() .eq. pelist(1)) then
        if (maxval(fdata_per_rank) .ne. fdata_0d .or. &
            minval(fdata_per_rank) .ne. fdata_0d) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_0d:" &
                           //" the scalar field data is not consistent across" &
                           //" all ranks in the I/O domain pelist.")
        endif
    endif
    io_domain => null()
    deallocate(pelist)
    deallocate(fdata_per_rank)

   !Set the dimension sizes for the field.  These correspond to:
   !field_dimension_sizes(XIDX) = x-dimension size
   !field_dimension_sizes(YIDX) = y-dimension size
   !field_dimension_sizes(CIDX) = c-dimension size
   !field_dimension_sizes(ZIDX) = z-dimension size
   !field_dimension_sizes(HIDX) = h-dimension size
   !field_dimension_sizes(TIDX) = t-dimension size
   !field_dimension_sizes(UIDX) = u-dimension size
   !field_dimension_sizes(CCIDX) = cc-dimension size
    field_dimension_sizes = 1

   !Set the ordering of the dimensions for the field.
    field_dimension_order(1) = TIDX

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             field_dimension_order, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
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
                                                         fdata_1d_axes, &
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
    integer(INT_KIND),dimension(1)         :: fdata_1d_axes     !<An array describing the axes for the data.
    type(domainUG),intent(in),target       :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional            :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional               :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional   :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional   :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional  :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional  :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                      :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND)                 :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !<Array of dimension sizes for the field.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_1d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that at least one axis was registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_1d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    endif

   !Make sure that the first dimension of the field is a "compressed" axis,
   !and that it corresponds to an axis that has been registered to the
   !restart object.
    field_dimension_sizes = 1
    if (fdata_1d_axes(1) .eq. CIDX) then
        if (.not. allocated(fileObj%axes(CIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_1d:" &
                           //" a compressed c-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_1d,1) .ne. fileObj%axes(CIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_1d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed c-axis.")
        endif
        field_dimension_sizes(CIDX) = size(fdata_1d,1)
    elseif (fdata_1d_axes(1) .eq. HIDX) then
        if (.not. allocated(fileObj%axes(HIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_1d:" &
                           //" a compressed h-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_1d,1) .ne. fileObj%axes(HIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_1d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed h-axis.")
        endif
        field_dimension_sizes(HIDX) = size(fdata_1d,1)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_1d:" &
                       //" One dimensional fields must be compressed.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             fdata_1d_axes, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
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
                                                         fdata_2d_axes, &
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
    integer(INT_KIND),dimension(2)           :: fdata_2d_axes     !<An array describing the axes for the data.
    type(domainUG),intent(in),target         :: domain            !<An unstructured mpp_domain.
    logical,intent(in),optional              :: mandatory         !<Flag telling if the field is mandatory for the restart.
    real,intent(in),optional                 :: data_default      !<A default value for the data.
    character(len=*),intent(in),optional     :: longname          !<A more descriptive name of the field.
    character(len=*),intent(in),optional     :: units             !<Units for the field.
    logical(INT_KIND),intent(in),optional    :: read_only         !<Tells whether or not the variable will be written to the restart file.
    logical(INT_KIND),intent(in),optional    :: restart_owns_data !<Tells if the data will be deallocated when the restart object is deallocated.
    integer(INT_KIND)                        :: restart_index     !<Index of the inputted field in the fileObj%var array.

   !Local variables
    integer(INT_KIND)                 :: index_field           !<Index of the inputted field in the fileObj%var array.
    integer(INT_KIND),dimension(NIDX) :: field_dimension_sizes !<Array of dimension sizes for the field.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" you must first call fms_io_init")
    endif

   !Make sure that at least one axis was registered to the restart object.
    if (.not. allocated(fileObj%axes)) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" no axes have been registered for the restart" &
                       //" object.")
    endif

   !Make sure that the first dimension of the field is a "compressed" axis,
   !and that it corresponds to an axis that has been registered to the
   !restart object.
    field_dimension_sizes = 1
    if (fdata_2d_axes(1) .eq. CIDX) then
        if (.not. allocated(fileObj%axes(CIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" a compressed c-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_2d,1) .ne. fileObj%axes(CIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed c-axis.")
        endif
        field_dimension_sizes(CIDX) = size(fdata_2d,1)
    elseif (fdata_2d_axes(1) .eq. HIDX) then
        if (.not. allocated(fileObj%axes(HIDX)%idx)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" a compressed h-axis was not registered" &
                           //" to the restart object.")
        endif
        if (size(fdata_2d,1) .ne. fileObj%axes(HIDX)%nelems_for_current_rank) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" compressed h-axis.")
        endif
        field_dimension_sizes(HIDX) = size(fdata_2d,1)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" The first dimension of the field must be a" &
                       //" compressed dimension.")
    endif

   !Make sure that the second dimension of the inputted field corresponds to
   !either a registered z- or cc-axis.
    if (fdata_2d_axes(2) .eq. ZIDX) then
        if (.not. associated(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" a z-axis was not registered to the" &
                           //" restart object.")
        endif
        if (size(fdata_2d,2) .ne. size(fileObj%axes(ZIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" z-axis.")
        endif
        field_dimension_sizes(ZIDX) = size(fdata_2d,2)
    elseif (fdata_2d_axes(2) .eq. CCIDX) then
        if (.not. associated(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" a cc-axis was not registered to the" &
                           //" restart object.")
        endif
        if (size(fdata_2d,2) .ne. size(fileObj%axes(CCIDX)%data)) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_register_restart_field_i_2d:" &
                           //" the size of the input data does not" &
                           //" match the size of the registered" &
                           //" cc-axis.")
        endif
        field_dimension_sizes(CCIDX) = size(fdata_2d,2)
    else
        call mpp_error(FATAL, &
                       "fms_io_unstructured_register_restart_field_i_2d:" &
                       //" unsupported axis parameter for the second" &
                       //" dimension of the field.")
    endif

   !Add a field to a restart object (restart_file_type).  Get the index of the
   !inputted field in the fileObj%var array.
    call fms_io_unstructured_setup_one_field(fileObj, &
                                             filename, &
                                             fieldname, &
                                             fdata_2d_axes, &
                                             field_dimension_sizes, &
                                             index_field, &
                                             domain, &
                                             mandatory=mandatory, &
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
