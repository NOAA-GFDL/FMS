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

!> @brief Do the time_none reduction method (i.e copy the correct portion of the input data)
subroutine DO_TIME_NONE_ (data_out, data_in, mask, bounds_in, bounds_out, missing_value)
  real(FMS_TRM_KIND_),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(FMS_TRM_KIND_),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(FMS_TRM_KIND_),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

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

  where (mask(is_in:ie_in, js_in:je_in, ks_in:ke_in, :))
    data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
    data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
  elsewhere
    data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = missing_value
  end where

end subroutine DO_TIME_NONE_

!> @brief Do the time_min reduction method (i.e maintain the minimum value of the averaging time)
subroutine DO_TIME_MIN_ (data_out, data_in, mask, bounds_in, bounds_out, missing_value)
  real(FMS_TRM_KIND_),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(FMS_TRM_KIND_),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(FMS_TRM_KIND_),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  integer :: i, j, k, l !< For looping

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

  do i = 0, ie_out - is_out
    do j = 0, je_out - js_out
      do k = 0, ke_out - ks_out
        do l = 0, size(data_out, 4) - 1
          if (mask(is_in + i, js_in + j, ks_in + k, l + 1)) then
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .gt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          else
            data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = missing_value
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine DO_TIME_MIN_

!> @brief Do the time_max reduction method (i.e maintain the maximum value of the averaging time)
subroutine DO_TIME_MAX_ (data_out, data_in, mask, bounds_in, bounds_out, missing_value)
  real(FMS_TRM_KIND_),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(FMS_TRM_KIND_),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(FMS_TRM_KIND_),       intent(in)    :: missing_value       !< Missing_value for data points that are masked

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer

  integer :: i, j, k, l !< For looping

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

  do i = 0, ie_out - is_out
    do j = 0, je_out - js_out
      do k = 0, ke_out - ks_out
        do l = 0, size(data_out, 4) - 1
          if (mask(is_in + i, js_in + j, ks_in + k, l + 1)) then
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .lt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          else
            data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = missing_value
          endif
        enddo
      enddo
    enddo
  enddo
end subroutine DO_TIME_MAX_