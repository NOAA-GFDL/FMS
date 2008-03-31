subroutine _SUBNAME_(errortype, errormsg1, scalar, errormsg2, array, errormsg3)
  integer,            intent(in) :: errortype
  _ARRAY1TYPE_,               intent(in) :: scalar
  _ARRAY2TYPE_, dimension(:), intent(in) :: array
  character(len=*),   intent(in) :: errormsg1, errormsg2
  character(len=*),   intent(in), optional :: errormsg3

  call mpp_error( errortype, errormsg1, (/scalar/), errormsg2, array, errormsg3)

end subroutine _SUBNAME_
