module metadata_transfer_mod
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

end module metadata_transfer_mod
