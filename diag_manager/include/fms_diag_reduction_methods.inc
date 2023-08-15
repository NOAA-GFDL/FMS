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

function DO_TIME_NONE_ (data_out, data_in, mask, bounds) &
result(error_msg)
  real(FMS_TRM_KIND_), intent(inout) :: data_out(:,:,:,:,:)
  real(FMS_TRM_KIND_), intent(in)    :: data_in(:,:,:,:)
  logical, intent(in) :: mask(:,:,:,:)
  type(fmsDiagIbounds_type), intent(in) :: bounds

  character(len=50) :: error_msg
  integer :: is, ie, js, je, ks, ke

  is = bounds%get_imin()
  ie = bounds%get_imax()
  js = bounds%get_jmin()
  je = bounds%get_jmax()
  ks = bounds%get_kmin()
  ke = bounds%get_kmax()

  data_out(is:ie, js:je, ks:ke, :, 1) = data_in(:,:,:,:)

  error_msg = ""
end function DO_TIME_NONE_