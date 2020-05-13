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

!> @file
!! @brief Tests the read_input_nml subroutine
!! @author Colin Gladue
!! @email gfdl.climate.model.info@noaa.gov

program test_read_input_nml

  use mpp_mod, only : mpp_init, mpp_exit
  use mpp_mod, only : mpp_error, FATAL, NOTE
  use mpp_mod, only : read_input_nml, mpp_get_current_pelist_name

#include<file_version.h>

character(len=200) :: line !< Current line being read
character(len=200) :: linelog, linenml !< Current log and nml lines being read
character(len=128) :: filename !< Name of file being read
integer :: n !< Looping variable
logical :: version_bool, filename_bool !< Booleans that tell whether or not the
                                       !! lines where version and filename should be 
                                       !! written in logfile have been found
                                       !! yet. Default value is false

namelist /test_read_input_nml_nml/ test_numb

open(10, file="test_numb.nml", form="formatted", status="old")
read(10, nml = test_read_input_nml_nml)
close(10)

if (test_numb.eq.1) then
  ! Test 1: Tests the subroutine on a valid input nml full of data, 
  ! with no arguments passed to read_input_nml()

  filename = "input.nml"
  call mpp_init() ! Initialize mpp
  call read_input_nml()
  call mpp_exit()
  open(110, file='logfile.000000.out', iostat=ioslog) ! Open logfile
  open(111, file='input.nml', iostat=iosnml) ! Open of input nml
  do while (ioslog.eq.0) ! Check for the first two written lines, then stay at this
                         ! position so we can compare from here on to the namelist
    read(110, '(A)', iostat=ioslog) line
    ! Check if we have found the version line
    if (index(line, "READ_INPUT_NML: "//trim(version)).ne.0) then
      version_bool = .true.
    end if
    ! Check if we have found the filename line
    if (index(line, "READ_INPUT_NML: "//trim(filename)).ne.0) then
      filename_bool = .true.
    end if
    ! Check if we have found all we are looking for 
    if (version_bool.and.filename_bool) then
      write(*,*) "SUCCESS: Found the first 2 written lines, version&
                       & and filename"
      exit ! Successful test portion
    end if
    ! If we have reached the end of the file, was anything missed?
    if (index(line, "Total runtime").ne.0) then
      if (.not.version_bool) then
        call mpp_error(FATAL, "Version not written to &
                                                    &logfile")
      else if (.not.filename_bool) then
        call mpp_error(FATAL, "Filename not written to &
                                                    &logfile")
      else
        call mpp_error(FATAL, "Logfile not written to by &
                                            &read_input_nml correctly.")
      end if
    end if
  end do
  ! Make sure the read pointers for logfile and input nml are in the same
  ! location below.
  read(110, '(A)', iostat=ioslog) linelog
  read(111, '(A)', iostat=iosnml) linenml
  do while (TRIM(linelog).ne.TRIM(" "//linenml))
    read(110, '(A)', iostat=ioslog) linelog
  end do
  ! Compare contents of logfile and the input nml
  do while (iosnml.eq.0)
    if (TRIM(linelog).ne.TRIM(" "//linenml)) then
      call mpp_error(FATAL, linelog//" -Does not equal- "//&
                          &linenml//". Namelist not written to logfile &
                                                          &correctly")
    end if
    read(110, '(A)', iostat=ioslog) linelog
    read(111, '(A)', iostat=iosnml) linenml
  end do
  write(*,*) "SUCCESS: Matched all lines from input nml to the logfile"
  close(110)
  close(111)

else if (test_numb.eq.2) then
  ! Test 2: Tests the same valid input nml, but with a pelist_name_in passed as an 
  ! argument to read_input_nml().  

  filename = "input_alternative.nml"
  call mpp_init() ! Initialize mpp
  call read_input_nml("alternative")
  call mpp_exit()
  open(110, file='logfile.000000.out', iostat=ioslog) ! Open logfile
  open(111, file='input_alternative.nml', iostat=iosnml) ! Open of input nml
  do while (ioslog.eq.0) ! Check for the first two written lines, then stay at this
                         ! position so we can compare from here on to the namelist
    read(110, '(A)', iostat=ioslog) line
    ! Check if we have found the version line
    if (index(line, "READ_INPUT_NML: "//trim(version)).ne.0) then
      version_bool = .true.
    end if
    ! Check if we have found the filename line
    if (index(line, "READ_INPUT_NML: "//trim(filename)).ne.0) then
      filename_bool = .true.
    end if
    ! Check if we have found all we are looking for 
    if (version_bool.and.filename_bool) then
      write(*,*) "SUCCESS: Found the first 2 written lines, version&
                       & and filename"
      exit ! Successful test portion
    end if
    ! If we have reached the end of the file, was anything missed?
    if (index(line, "Total runtime").ne.0) then
      if (.not.version_bool) then
        call mpp_error(FATAL, "Version not written to &
                                                    &logfile")
      else if (.not.filename_bool) then
        call mpp_error(FATAL, "Filename not written to &
                                                    &logfile")
      else
        call mpp_error(FATAL, "Logfile not written to by &
                                            &read_input_nml correctly.")
      end if
    end if
  end do
  ! Make sure the read pointers for logfile and input nml are in the same
  ! location below.
  read(110, '(A)', iostat=ioslog) linelog
  read(111, '(A)', iostat=iosnml) linenml
  do while (TRIM(linelog).ne.TRIM(" "//linenml))
    read(110, '(A)', iostat=ioslog) linelog
  end do
  ! Compare contents of logfile and the input nml
  do while (iosnml.eq.0)
    if (TRIM(linelog).ne.TRIM(" "//linenml)) then
      call mpp_error(FATAL, linelog//" -Does not equal- "//&
                          &linenml//". Namelist not written to logfile &
                                                          &correctly")
    end if
    read(110, '(A)', iostat=ioslog) linelog
    read(111, '(A)', iostat=iosnml) linenml
  end do
  write(*,*) "SUCCESS: Matched all lines from input nml to the logfile"
  close(110)
  close(111)

else if (test_numb.eq.3) then
  ! Test 3: Tests with an invalid pelist_name_in pass as an argument. An invalid
  ! pelist_name_in would be one who's size is greater than local pelist_name

  call mpp_init ! Initialize mpp        
  call read_input_nml(mpp_get_current_pelist_name()//"e")
                                                          ! Call read_input_nml
                                                          ! with the local
                                                          ! pelist_name plus an
                                                          ! extra character "e"
  call mpp_exit() ! Exit mpp

else if (test_numb.eq.4) then
  ! Test 4: Tests an empty input nml. No arguments are passed. 

  filename = "input.nml"
  call mpp_init() ! Initialize mpp
  call read_input_nml()
  call mpp_exit()
  open(44, file='logfile.000000.out', iostat=ios)
  do while (ios.eq.0) ! Check for the first two written lines
    read(44, '(A)', iostat=ios) line
    ! Check if we have found the version line
    if (index(line, "READ_INPUT_NML: "//trim(version)).ne.0) then
      version_bool = .true.
    end if
    ! Check if we have found the filename line
    if (index(line, "READ_INPUT_NML: "//trim(filename)).ne.0) then
      filename_bool = .true.
    end if
    ! Check if we have found all we are looking for 
    if (version_bool.and.filename_bool) then
      write(*,*) "SUCCESS: Found the first 2 written lines, version&
                       & and filename"
      exit ! Successful test
    end if
    ! If we have reached the end of the file, was anything missed?
    if (index(line, "Total runtime").ne.0) then
      if (.not.version_bool) then
        call mpp_error(FATAL, "Version not written to &
                                                    &logfile")
      else if (.not.filename_bool) then
        call mpp_error(FATAL, "Filename not written to &
                                                    &logfile")
      else
        call mpp_error(FATAL, "Logfile not written to by &
                                            &read_input_nml correctly.")
      end if
    end if
  end do
  close(44)
end if

end program test_read_input_nml
