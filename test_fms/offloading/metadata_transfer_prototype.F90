module metadata_offload
  use platform_mod

  implicit none
  public

  integer, parameter :: real8_type = 1
  integer, parameter :: real4_type = 2
  integer, parameter :: int8_type = 3
  integer, parameter :: int4_type = 4
  integer, parameter :: str_type = 5

  type :: metadata_type
    integer            :: attribute_type
    real(kind=r8_kind) :: attribute_value(2)
    character(len=50)  :: attribute_name
  end type

end module metadata_offload

program test
  use fms_mod, only: fms_init, fms_end, string
  use mpp_mod
  use metadata_offload
  use platform_mod
  use mpi
  use, intrinsic :: iso_c_binding

  implicit none

  type(metadata_type) :: file_metadata(2)
  INTEGER             :: dummy_int
  REAL                :: dummy_real
  INTEGER :: lengths(0:2)
  INTEGER :: types(0:2)
  INTEGER(KIND=MPI_ADDRESS_KIND) :: displacements(0:2)
  INTEGER :: metadata_t
  integer :: sender, receiver
  integer :: ierror

  call fms_init()

  lengths = (/1, 2, 50/)
  displacements = (/0, sizeof(dummy_int), sizeof(dummy_int) + sizeof(dummy_real)/)
  types = (/MPI_INTEGER, MPI_REAL, MPI_CHARACTER/)
  call MPI_Type_create_struct(3, lengths, displacements, types, metadata_t, ierror)
  call MPI_Type_commit(metadata_t, ierror)

  receiver = 1
  sender = 0
  if (mpp_pe() .eq. mpp_root_pe()) then
    !A root pe keeps track of all of the metadata that has been sent, one metadata type per attribute type
    file_metadata(1)%attribute_name  = "_FillValue"//c_null_char
    file_metadata(1)%attribute_type  = real8_type
    file_metadata(1)%attribute_value = -666_r8_kind

    file_metadata(2)%attribute_name  = "missing_value"//c_null_char
    file_metadata(2)%attribute_type  = real8_type
    file_metadata(2)%attribute_value = -666_r8_kind
    CALL MPI_Send(file_metadata, 2, metadata_t, receiver, 0, MPI_COMM_WORLD, ierror)
  else
    CALL MPI_Recv(file_metadata, 2, metadata_t, sender, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
    call dump_metadata(file_metadata)
  endif

  call fms_end()

  contains

  subroutine dump_metadata(this)
    type(metadata_type), intent(in) :: this(:)

    integer :: i
    do i = 1, size(this)
      print *, mpp_pe(), " knows that the attribute_name is ", trim(ADJUSTL(this(i)%attribute_name)), &
        " and the attribute_type is ", string(this(i)%attribute_type), &
        " and the value is ", this(i)%attribute_value
    enddo
  end subroutine
end program