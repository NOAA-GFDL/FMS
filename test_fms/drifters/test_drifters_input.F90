
program test_drifters_input
  use drifters_input_mod
  implicit none
  character(len=128) :: ermesg
  integer :: i
  
  type(drifters_input_type) :: obj
  
  call drifters_input_new(obj, 'input.nc', ermesg)
  if(ermesg/='') print *,'ERROR: ', ermesg

  print *,'field_names:'
  do i = 1, size(obj%field_names)
     print *,trim(obj%field_names(i))
  enddo

  print *,'velocity_names:'
  do i = 1, size(obj%velocity_names)
     print *,trim(obj%velocity_names(i))
  enddo

  print *,'ids = ', obj%ids

  print *,'positions: '
  do i = 1, size(obj%positions, 2)
     print *,obj%positions(:,i)
  enddo

  call drifters_input_del(obj, ermesg)
end program test_drifters_input
