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

function DO_TIME_NONE_ (data_out, data_in, mask, bounds_in, bounds_out) &
result(error_msg)
  real(FMS_TRM_KIND_), intent(inout) :: data_out(:,:,:,:,:)
  real(FMS_TRM_KIND_), intent(in)    :: data_in(:,:,:,:)
  logical, intent(in), target :: mask(:,:,:,:)
  type(fmsDiagIbounds_type), intent(in) :: bounds_in
  type(fmsDiagIbounds_type), intent(in) :: bounds_out

  character(len=50) :: error_msg
  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out

  is_out = bounds_out%get_imin()
  ie_out = bounds_out%get_imax()
  js_out = bounds_out%get_jmin()
  je_out = bounds_out%get_jmax()
  ks_out = bounds_out%get_kmin()
  ke_out = bounds_out%get_kmax()

  is_in = bounds_in%get_imin()
  ie_in = bounds_in%get_imax()
  js_in = bounds_in%get_jmin()
  je_in = bounds_in%get_jmax()
  ks_in = bounds_in%get_kmin()
  ke_in = bounds_in%get_kmax()

!  write ((mpp_pe()+1)*100, *) "BOUNDS", lbound(mask), ubound(mask)
  write ((mpp_pe()+1)*100, *) "indices_in:", is_in,ie_in, js_in,je_in, ks_in,ke_in
  write ((mpp_pe()+1)*100, *) "indices_out:", is_out,ie_out, js_out,je_out, ks_out,ke_out
  write ((mpp_pe()+1)*100, *) "bounds_in", lbound(data_in), ubound(data_in)
  write ((mpp_pe()+1)*100, *) "bounds_in", lbound(data_out), ubound(data_out)

  write ((mpp_pe()+1)*100, *) "WUTT:", mask(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
  where (mask(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)) &
    data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
    data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
    
    write ((mpp_pe()+1)*100, *) "data_in:", data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
    write ((mpp_pe()+1)*100, *) "data_out:", data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1)
  error_msg = ""
end function DO_TIME_NONE_