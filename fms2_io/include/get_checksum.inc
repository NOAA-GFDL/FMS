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

!> @brief Given a data array, return a string containing the mpp_checksum
!!        in hex.
function get_checksum_0d(data) result(chksum)

  class(*),  intent(in) :: data !< Data to be checksummed.
  character(len=16) :: chksum

  integer,dimension(1) :: myrank

  myrank(1) = mpp_pe()
  chksum = ""
  select type(data)
    type is (integer(int32))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (integer(int64))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real32))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
  end select

end function get_checksum_0d
!> @brief Given a data array, return a string containing the mpp_checksum
!!        in hex.
function get_checksum_1d(data) result(chksum)

  class(*), dimension(:), intent(in) :: data !< Data to be checksummed.
  character(len=16) :: chksum

  integer,dimension(1) :: myrank

  myrank(1) = mpp_pe()
  chksum = ""
  select type(data)
    type is (integer(int32))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (integer(int64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real32))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
  end select

end function get_checksum_1d
!> @brief Given a data array, return a string containing the mpp_checksum
!!        in hex.
function get_checksum_2d(data) result(chksum)

  class(*), dimension(:,:), intent(in) :: data !< Data to be checksummed.
  character(len=16) :: chksum

  integer,dimension(1) :: myrank

  myrank(1) = mpp_pe()
  chksum = ""
  select type(data)
    type is (integer(int32))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (integer(int64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real32))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
  end select

end function get_checksum_2d
!> @brief Given a data array, return a string containing the mpp_checksum
!!        in hex.
function get_checksum_3d(data) result(chksum)

  class(*), dimension(:,:,:), intent(in) :: data !< Data to be checksummed.
  character(len=16) :: chksum

  integer,dimension(1) :: myrank

  myrank(1) = mpp_pe()
  chksum = ""
  select type(data)
    type is (integer(int32))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (integer(int64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real32))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
  end select

end function get_checksum_3d
!> @brief Given a data array, return a string containing the mpp_checksum
!!        in hex.
function get_checksum_4d(data) result(chksum)

  class(*), dimension(:,:,:,:), intent(in) :: data !< Data to be checksummed.
  character(len=16) :: chksum

  integer,dimension(1) :: myrank

  myrank(1) = mpp_pe()
  chksum = ""
  select type(data)
    type is (integer(int32))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (integer(int64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real32))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
  end select

end function get_checksum_4d
!> @brief Given a data array, return a string containing the mpp_checksum
!!        in hex.
function get_checksum_5d(data) result(chksum)

  class(*), dimension(:,:,:,:,:), intent(in) :: data !< Data to be checksummed.
  character(len=16) :: chksum

  integer,dimension(1) :: myrank

  myrank(1) = mpp_pe()
  chksum = ""
  select type(data)
    type is (integer(int32))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (integer(int64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real32))
!     write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
  end select

end function get_checksum_5d
