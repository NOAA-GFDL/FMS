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
