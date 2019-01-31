include(`macros.m4')dnl
!> @brief Given a data array, return a string containing the mpp_checksum
!!        in hex.
`function get_checksum_'NUM_DIMS`d(data) result(chksum)'

  class(*), dim_declare(NUM_DIMS) intent(in) :: data !< Data to be checksummed.
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
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
    type is (real(real64))
      write(chksum, "(Z16)") mpp_chksum(data, myrank)
  end select

`end function get_checksum_'NUM_DIMS`d'
