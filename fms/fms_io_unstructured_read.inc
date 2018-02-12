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
!>Read in a scalar field from a file associated with an unstructured mpp
!!domain.
subroutine fms_io_unstructured_read_r_scalar(filename, &
                                             fieldname, &
                                             fdata, &
                                             domain, &
                                             timelevel, &
                                             start, &
                                             nread, &
                                             threading)

   !Inputs/Outputs
    character(len=*),intent(in)              :: filename  !<The name of a file.
    character(len=*),intent(in)              :: fieldname !<The name of field in the input file.
    real,intent(inout)                       :: fdata     !<Data to be read in from the file.
    type(domainUG),intent(in)                :: domain    !<Unstructured mpp domain associated with the input file.
    integer,intent(in),optional              :: timelevel !<Time level at which the data is read in.
    integer,dimension(:),intent(in),optional :: start     !<NetCDF "Corner" indices for the data.
    integer,dimension(:),intent(in),optional :: nread     !<NetCDF "edge lengths" for the data.
    integer,intent(in),optional              :: threading !<Threading flag.

   !Local variables
    real,dimension(1) :: tmp    !<Dummy variable.
    integer           :: tlevel !<Dummy variable.

   !Use the 1D case.
    if (present(timelevel)) then
        if (tlevel .le. 0) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_read_r_scalar:" &
                           //" the inputted time level must be at" &
                           //" least one.")
        endif
        tlevel = timelevel
    else
        tlevel = 1
    endif
    call fms_io_unstructured_read_r_1D(filename, &
                                       fieldname, &
                                       tmp, &
                                       domain, &
                                       tlevel, &
                                       start, &
                                       nread, &
                                       threading)
    fdata = tmp(1)

    return
end subroutine fms_io_unstructured_read_r_scalar

!------------------------------------------------------------------------------
!>Read in a one dimensional "compressed" field from a file associated with
!!an unstructured mpp domain.
subroutine fms_io_unstructured_read_r_1D(filename, &
                                         fieldname, &
                                         fdata, &
                                         domain, &
                                         timelevel, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/Outputs
    character(len=*),intent(in)              :: filename  !<The name of a file.
    character(len=*),intent(in)              :: fieldname !<The name of field in the input file.
    real,dimension(:),intent(inout)          :: fdata     !<Data to be read in from the file.
    type(domainUG),intent(in)                :: domain    !<Unstructured mpp domain associated with the input file.
    integer,intent(in),optional              :: timelevel !<Time level at which the data is read in.
    integer,dimension(:),intent(in),optional :: start     !<NetCDF "Corner" indices for the data.
    integer,dimension(:),intent(in),optional :: nread     !<NetCDF "edge lengths" for the data.
    integer,intent(in),optional              :: threading !<Threading flag.

   !Local variables
    logical(INT_KIND)  :: found_file  !<Flag telling if the input file or any of its variants exist.
    character(len=256) :: fname       !<Name of file that is actually found.
    logical(INT_KIND)  :: read_dist   !<Flag telling if the file is "distributed" (has I/O domain tile id appended onto the end).
    integer(INT_KIND)  :: funit       !<File unit for the inputted file.
    integer(INT_KIND)  :: file_index  !<Index of the inputted file in the "files_read" module array.
    integer(INT_KIND)  :: index_field !<Index of the inputted field in the files_read(file_index)%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_read_r_1D:" &
                       //" you must first call fms_io_init.")
    endif

   !Get the full name of the input file.
    found_file = fms_io_unstructured_get_file_name(filename, &
                                                   fname, &
                                                   read_dist, &
                                                   domain)

   !If the file does not exit, then throw a fatal error.
    if (.not. found_file) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_read_r_1D:" &
                       //" file "//trim(filename) &
                       //" (with the consideration of the domain tile ids)" &
                       //" was not found.")
    endif

   !Get the file unit and in the index of the file in the "files_read" module
   !array.
    call fms_io_unstructured_get_file_unit(fname, &
                                           funit, &
                                           file_index, &
                                           read_dist, &
                                           domain)

   !Get the index of the inputted field in the files_read(file_index)%var
   !array.
    call get_field_id(funit, &
                      file_index, &
                      fieldname, &
                      index_field, &
                      .false., &
                      .false.)

   !Read in the field data.
    if (files_read(file_index)%var(index_field)%is_dimvar) then
        call mpp_get_axis_data(files_read(file_index)%var(index_field)%axis, &
                               fdata)
    else
        call mpp_io_unstructured_read(funit, &
                                      files_read(file_index)%var(index_field)%field, &
                                      domain, &
                                      fdata, &
                                      timelevel, &
                                      start, &
                                      nread, &
                                      threading)
    endif

    return
end subroutine fms_io_unstructured_read_r_1D

!------------------------------------------------------------------------------
!>Read in a two dimensional "compressed" field from a file associated with
!!an unstructured mpp domain.
subroutine fms_io_unstructured_read_r_2D(filename, &
                                         fieldname, &
                                         fdata, &
                                         domain, &
                                         timelevel, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/Outputs
    character(len=*),intent(in)              :: filename  !<The name of a file.
    character(len=*),intent(in)              :: fieldname !<The name of field in the input file.
    real,dimension(:,:),intent(inout)        :: fdata     !<Data to be read in from the file.
    type(domainUG),intent(in)                :: domain    !<Unstructured mpp domain associated with the input file.
    integer,intent(in),optional              :: timelevel !<Time level at which the data is read in.
    integer,dimension(:),intent(in),optional :: start     !<NetCDF "Corner" indices for the data.
    integer,dimension(:),intent(in),optional :: nread     !<NetCDF "edge lengths" for the data.
    integer,intent(in),optional              :: threading !<Threading flag.

   !Local variables
    logical(INT_KIND)  :: found_file  !<Flag telling if the input file or any of its variants exist.
    character(len=256) :: fname       !<Name of file that is actually found.
    logical(INT_KIND)  :: read_dist   !<Flag telling if the file is "distributed" (has I/O domain tile id appended onto the end).
    integer(INT_KIND)  :: funit       !<File unit for the inputted file.
    integer(INT_KIND)  :: file_index  !<Index of the inputted file in the "files_read" module array.
    integer(INT_KIND)  :: index_field !<Index of the inputted field in the files_read(file_index)%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_read_r_2D:" &
                       //" you must first call fms_io_init.")
    endif

   !Get the full name of the input file.
    found_file = fms_io_unstructured_get_file_name(filename, &
                                                   fname, &
                                                   read_dist, &
                                                   domain)

   !If the file does not exit, then throw a fatal error.
    if (.not. found_file) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_read_r_2D:" &
                       //" file "//trim(filename) &
                       //" (with the consideration of the domain tile ids)" &
                       //" was not found.")
    endif

   !Get the file unit and in the index of the file in the "files_read" module
   !array.
    call fms_io_unstructured_get_file_unit(fname, &
                                           funit, &
                                           file_index, &
                                           read_dist, &
                                           domain)

   !Get the index of the inputted field in the files_read(file_index)%var
   !array.
    call get_field_id(funit, &
                      file_index, &
                      fieldname, &
                      index_field, &
                      .false., &
                      .false.)

   !Read in the field data.
    if (files_read(file_index)%var(index_field)%is_dimvar) then
        call mpp_get_axis_data(files_read(file_index)%var(index_field)%axis, &
                               fdata(:,1))
    else
        call mpp_io_unstructured_read(funit, &
                                      files_read(file_index)%var(index_field)%field, &
                                      domain, &
                                      fdata, &
                                      timelevel, &
                                      start, &
                                      nread, &
                                      threading)
    endif

    return
end subroutine fms_io_unstructured_read_r_2D

!------------------------------------------------------------------------------
!>Read in a three dimensional "compressed" field from a file associated with
!!an unstructured mpp domain.
subroutine fms_io_unstructured_read_r_3D(filename, &
                                         fieldname, &
                                         fdata, &
                                         domain, &
                                         timelevel, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/Outputs
    character(len=*),intent(in)              :: filename  !<The name of a file.
    character(len=*),intent(in)              :: fieldname !<The name of field in the input file.
    real,dimension(:,:,:),intent(inout)      :: fdata     !<Data to be read in from the file.
    type(domainUG),intent(in)                :: domain    !<Unstructured mpp domain associated with the input file.
    integer,intent(in),optional              :: timelevel !<Time level at which the data is read in.
    integer,dimension(:),intent(in),optional :: start     !<NetCDF "Corner" indices for the data.
    integer,dimension(:),intent(in),optional :: nread     !<NetCDF "edge lengths" for the data.
    integer,intent(in),optional              :: threading !<Threading flag.

   !Local variables
    logical(INT_KIND)  :: found_file  !<Flag telling if the input file or any of its variants exist.
    character(len=256) :: fname       !<Name of file that is actually found.
    logical(INT_KIND)  :: read_dist   !<Flag telling if the file is "distributed" (has I/O domain tile id appended onto the end).
    integer(INT_KIND)  :: funit       !<File unit for the inputted file.
    integer(INT_KIND)  :: file_index  !<Index of the inputted file in the "files_read" module array.
    integer(INT_KIND)  :: index_field !<Index of the inputted field in the files_read(file_index)%var array.

   !Make sure that the module has been initialized.
    if (.not. module_is_initialized) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_read_r_3D:" &
                       //" you must first call fms_io_init.")
    endif

   !Get the full name of the input file.
    found_file = fms_io_unstructured_get_file_name(filename, &
                                                   fname, &
                                                   read_dist, &
                                                   domain)

   !If the file does not exit, then throw a fatal error.
    if (.not. found_file) then
        call mpp_error(FATAL, &
                       "fms_io_unstructured_read_r_3D:" &
                       //" file "//trim(filename) &
                       //" (with the consideration of the domain tile ids)" &
                       //" was not found.")
    endif

   !Get the file unit and in the index of the file in the "files_read" module
   !array.
    call fms_io_unstructured_get_file_unit(fname, &
                                           funit, &
                                           file_index, &
                                           read_dist, &
                                           domain)

   !Get the index of the inputted field in the files_read(file_index)%var
   !array.
    call get_field_id(funit, &
                      file_index, &
                      fieldname, &
                      index_field, &
                      .false., &
                      .false.)

   !Read in the field data.
    if (files_read(file_index)%var(index_field)%is_dimvar) then
        call mpp_get_axis_data(files_read(file_index)%var(index_field)%axis, &
                               fdata(:,1,1))
    else
        call mpp_io_unstructured_read(funit, &
                                      files_read(file_index)%var(index_field)%field, &
                                      domain, &
                                      fdata, &
                                      timelevel, &
                                      start, &
                                      nread, &
                                      threading)
    endif

    return
end subroutine fms_io_unstructured_read_r_3D

!------------------------------------------------------------------------------
!>Read in a scalar field from a file associated with an unstructured mpp
!!domain.
subroutine fms_io_unstructured_read_i_scalar(filename, &
                                             fieldname, &
                                             fdata, &
                                             domain, &
                                             timelevel, &
                                             start, &
                                             nread, &
                                             threading)

   !Inputs/Outputs
    character(len=*),intent(in)              :: filename  !<The name of a file.
    character(len=*),intent(in)              :: fieldname !<The name of field in the input file.
    integer,intent(inout)                    :: fdata     !<Data to be read in from the file.
    type(domainUG),intent(in)                :: domain    !<Unstructured mpp domain associated with the input file.
    integer,intent(in),optional              :: timelevel !<Time level at which the data is read in.
    integer,dimension(:),intent(in),optional :: start     !<NetCDF "Corner" indices for the data.
    integer,dimension(:),intent(in),optional :: nread     !<NetCDF "edge lengths" for the data.
    integer,intent(in),optional              :: threading !<Threading flag.

   !Local variables
    real,dimension(1) :: tmp    !<Dummy variable.
    integer           :: tlevel !<Dummy variable.

   !Read in the data.
    if (present(timelevel)) then
        if (tlevel .le. 0) then
            call mpp_error(FATAL, &
                           "fms_io_unstructured_read_i_scalar:" &
                           //" the inputted time level must be at" &
                           //" least one.")
        endif
        tlevel = timelevel
    else
        tlevel = 1
    endif
    call fms_io_unstructured_read_r_1D(filename, &
                                       fieldname, &
                                       tmp, &
                                       domain, &
                                       tlevel, &
                                       start, &
                                       nread, &
                                       threading)
    fdata = ceiling(tmp(1))

    return
end subroutine fms_io_unstructured_read_i_scalar

!------------------------------------------------------------------------------
!>Read in a one dimensional "compressed" field from a file associated with
!!an unstructured mpp domain.
subroutine fms_io_unstructured_read_i_1D(filename, &
                                         fieldname, &
                                         fdata, &
                                         domain, &
                                         timelevel, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/Outputs
    character(len=*),intent(in)              :: filename  !<The name of a file.
    character(len=*),intent(in)              :: fieldname !<The name of field in the input file.
    integer,dimension(:),intent(inout)       :: fdata     !<Data to be read in from the file.
    type(domainUG),intent(in)                :: domain    !<Unstructured mpp domain associated with the input file.
    integer,intent(in),optional              :: timelevel !<Time level at which the data is read in.
    integer,dimension(:),intent(in),optional :: start     !<NetCDF "Corner" indices for the data.
    integer,dimension(:),intent(in),optional :: nread     !<NetCDF "edge lengths" for the data.
    integer,intent(in),optional              :: threading !<Threading flag.

   !Local variables
    real,dimension(size(fdata)) :: tmp !<Dummy variable.
    integer(INT_KIND)           :: i   !<Loop variable.

   !Read in the data.
    call fms_io_unstructured_read_r_1D(filename, &
                                       fieldname, &
                                       tmp, &
                                       domain, &
                                       timelevel, &
                                       start, &
                                       nread, &
                                       threading)
    do i = 1,size(fdata)
        fdata(i) = ceiling(tmp(i))
    enddo

    return
end subroutine fms_io_unstructured_read_i_1D

!------------------------------------------------------------------------------
!>Read in a two dimensional "compressed" field from a file associated with
!!an unstructured mpp domain.
subroutine fms_io_unstructured_read_i_2D(filename, &
                                         fieldname, &
                                         fdata, &
                                         domain, &
                                         timelevel, &
                                         start, &
                                         nread, &
                                         threading)

   !Inputs/Outputs
    character(len=*),intent(in)              :: filename  !<The name of a file.
    character(len=*),intent(in)              :: fieldname !<The name of field in the input file.
    integer,dimension(:,:),intent(inout)     :: fdata     !<Data to be read in from the file.
    type(domainUG),intent(in)                :: domain    !<Unstructured mpp domain associated with the input file.
    integer,intent(in),optional              :: timelevel !<Time level at which the data is read in.
    integer,dimension(:),intent(in),optional :: start     !<NetCDF "Corner" indices for the data.
    integer,dimension(:),intent(in),optional :: nread     !<NetCDF "edge lengths" for the data.
    integer,intent(in),optional              :: threading !<Threading flag.

   !Local variables
    real,dimension(size(fdata,1),size(fdata,2)) :: tmp !<Dummy variable.
    integer(INT_KIND)                           :: i   !<Loop variable.
    integer(INT_KIND)                           :: j   !<Loop variable.

   !Read in the data.
    call fms_io_unstructured_read_r_2D(filename, &
                                       fieldname, &
                                       tmp, &
                                       domain, &
                                       timelevel, &
                                       start, &
                                       nread, &
                                       threading)
    do i = 1,size(fdata,2)
        do j = 1,size(fdata,1)
            fdata(j,i) = ceiling(tmp(j,i))
        enddo
    enddo

    return
end subroutine fms_io_unstructured_read_i_2D

!------------------------------------------------------------------------------

!----------
