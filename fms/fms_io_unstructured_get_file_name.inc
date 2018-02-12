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

!>For an inputted file name, check if it or any of its variants exist.
!!For a file named "foo", variants checked (in order) include:
!!
!!    foo
!!    foo.nc
!!    foo.<domain_tile_id_string>.nc
!!    foo.nc.<IO_domain_tile_id_string>
!!    foo.<domain_tile_id_string>.nc.<IO_domain_tile_id_string>
!!    foo.<ensemble_id>
!!    foo.<ensemble_id>.nc
!!    foo.<ensemble_id>.<domain_tile_id_string>.nc
!!    foo.<ensemble_id>.nc.<IO_domain_tile_id_string>
!!    foo.<ensemble_id>.<domain_tile_id_string>.nc.<IO_domain_tile_id_string>
!!
!!If a match is found, the value true is returned for the "does_file_exist"
!!flag.  In addition, the actual file name is returned and the "read_dist"
!!flag, which tells whether or not the filename contains the
!!IO_domain_tile_id_string appended.
!!
!!Should this be a subroutine instead of a funtion for clarity since it
!!returns more than one value?
function fms_io_unstructured_get_file_name(orig_file, &
                                           actual_file, &
                                           read_dist, &
                                           domain) &
                                           result(does_file_exist)

   !Inputs/Outputs
    character(len=*),intent(in)   :: orig_file       !<The name of file we're looking for.
    character(len=*),intent(out)  :: actual_file     !<Name of the file we found.
    logical(INT_KIND),intent(out) :: read_dist       !<Flag telling if the file is "distributed" (has IO domain tile id appended onto the end).
    type(domainUG),intent(in)     :: domain          !<Unstructured mpp domain.
    logical(INT_KIND)             :: does_file_exist !<Flag telling if the inputted file exists or one its variants.

   !Local variables
    logical(INT_KIND)      :: fexist          !<Flag that tells if a file exists.
    type(domainUG),pointer :: io_domain       !<Pointer to an unstructured I/O domain.
    integer(INT_KIND)      :: io_tile_id      !<Tile id for the I/O domain.
    character(len=256)     :: fname           !<A character buffer used to test different file names.
    character(len=512)     :: actual_file_tmp !<A character buffer used to test different file names.

   !Set the default return values for the function.
    actual_file = ""
    does_file_exist = .false.
    read_dist = .false.

   !Check if the file name does not contain ".nc".
    fexist = .false.
    if (index(orig_file,".nc",back=.true.) .eq. 0) then
        inquire(file=trim(orig_file),exist=fexist)
        if (fexist) then
            actual_file = orig_file
            does_file_exist = .true.
            return
        endif
    endif

   !If necessary, add the correct domain ".tilexxxx" string to the inputted
   !file name.  For a file named foo.nc, this would become foo.tilexxxx.nc.
   !Check if the new file name exists.
    call get_mosaic_tile_file_ug(orig_file, &
                                 actual_file, &
                                 domain)
    inquire(file=trim(actual_file),exist=fexist)
    if (fexist) then
        does_file_exist = .true.
        return
    endif

   !Point to the I/O domain for the unstructured grid.  This function call
   !will throw a fatal error if the I/O domain does not exist.
    io_domain => null()
    io_domain => mpp_get_UG_io_domain(domain)

   !Get the tile id for the I/O domain.
    io_tile_id = mpp_get_UG_domain_tile_id(io_domain)
    io_domain => null()

   !Check if the file has the I/O domain's tile id appended to the end of its
   !name.  For a file named foo.nc, this would become foo.nc.yyyy, where
   !"yyyy" would in reality be the I/O domain's tile id.  If the file exists,
   !then set the read_dist and does_file_exist flags to true and return.
    write(fname,'(a,i4.4)') trim(actual_file)//'.',io_tile_id
    inquire(file=trim(fname),exist=fexist)
    if (.not. fexist) then
        write(fname,'(a,i6.6)') trim(actual_file)//'.',io_tile_id
        inquire(file=trim(fname),exist=fexist)
    endif
    if (fexist) then
        read_dist = .true.
        does_file_exist = .true.
        return
    endif

   !Check if the file is part of an ensemble.
   !filename_appendix is a module variable.
    if (len_trim(filename_appendix) .gt. 0) then
        call get_instance_filename(orig_file, &
                                   actual_file)
        if (index(orig_file,'.nc',back=.true.) .eq. 0) then
            inquire(file=trim(actual_file),exist=fexist)
            if (fexist) then
                does_file_exist = .true.
                return
            endif
        endif

       !Make a local copy of "actual_file", and the use the local copy to
       !add the domain ".tilexxxx" string to "actual_file".
        actual_file_tmp = actual_file
        call get_mosaic_tile_file_ug(actual_file_tmp, &
                                     actual_file, &
                                     domain)
        inquire(file=trim(actual_file),exist=fexist)
        if (fexist) then
            does_file_exist = .true.
            return
        endif

       !Point to the I/O domain for the unstructured grid.  This function call
       !will throw a fatal error if the I/O domain does not exist.
        io_domain => mpp_get_UG_io_domain(domain)

       !Get the tile id for the I/O domain.
        io_tile_id = mpp_get_UG_domain_tile_id(io_domain)
        io_domain => null()

       !Check if the file has the I/O domain's tile id appended to the end of
       !its name.  If it does then set the read_dist and does_file_exist flags
       !to true and return.
        write(fname,'(a,i4.4)') trim(actual_file)//'.',io_tile_id
        inquire(file=trim(fname),exist=fexist)
        if (.not. fexist) then
            write(fname,'(a,i6.6)') trim(actual_file)//'.',io_tile_id
            inquire(file=trim(fname),exist=fexist)
        endif
        if (fexist) then
            read_dist = .true.
            does_file_exist = .true.
            return
        endif
    endif

    return
end function fms_io_unstructured_get_file_name

!------------------------------------------------------------------------------
