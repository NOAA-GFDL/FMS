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

! for any debug prints
#ifndef DEBUG_REDUCT
#define DEBUG_REDUCT .false.
#endif

!> @brief Do the time_none reduction method (i.e copy the correct portion of the input data)
subroutine DO_TIME_NONE_ (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(FMS_TRM_KIND_),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(FMS_TRM_KIND_),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
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

  if (is_masked) then
    where (mask(is_in:ie_in, js_in:je_in, ks_in:ke_in, :))
      data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
      data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
    elsewhere
      data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = missing_value
    end where
  else
    data_out(is_out:ie_out, js_out:je_out, ks_out:ke_out, :, 1) = &
      data_in(is_in:ie_in, js_in:je_in, ks_in:ke_in, :)
  endif

end subroutine DO_TIME_NONE_

!> @brief Do the time_min reduction method (i.e maintain the minimum value of the averaging time)
subroutine DO_TIME_MIN_ (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(FMS_TRM_KIND_),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(FMS_TRM_KIND_),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
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

  !> Seperated this loops for performance. If is_masked = .false. (i.e "mask" and "rmask" were never passed in)
  !! then mask will always be .True. so the if (mask) is redudant.
  if (is_masked) then
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
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
  else
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .gt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          enddo
        enddo
      enddo
    enddo
  endif

end subroutine DO_TIME_MIN_

!> @brief Do the time_max reduction method (i.e maintain the maximum value of the averaging time)
subroutine DO_TIME_MAX_ (data_out, data_in, mask, is_masked, bounds_in, bounds_out, missing_value)
  real(FMS_TRM_KIND_),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(FMS_TRM_KIND_),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
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

  !> Seperated this loops for performance. If is_masked = .false. (i.e "mask" and "rmask" were never passed in)
  !! then mask will always be .True. so the if (mask) is redudant.
  if (is_masked) then
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
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
  else
    do l = 0, size(data_out, 4) - 1
      do k = 0, ke_out - ks_out
        do j = 0, je_out - js_out
          do i = 0, ie_out - is_out
            if (data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) .lt. &
              data_in(is_in + i, js_in + j, ks_in + k, l + 1) ) then
                data_out(is_out + i, js_out + j, ks_out + k, l + 1, 1) = &
                  data_in(is_in +i, js_in + j, ks_in + k, l + 1)
            endif
          enddo
        enddo
      enddo
    enddo
  endif
end subroutine DO_TIME_MAX_

!> Update the output buffer for reductions that involve summation (sum, avg, rms, pow).
!! Elements of the running field output buffer (data_out) are set with the following:
!!
!!    buffer(l) = buffer(l) + (weight * field(l)) ^ pow
!!
!! Where l are the indices passed in through the bounds_in/out
subroutine DO_TIME_SUM_UPDATE_(data_out, weight_sum, data_in, mask, is_masked, bounds_in, bounds_out, &
                               missing_value, increase_counter, weight, pow)
  real(FMS_TRM_KIND_),       intent(inout) :: data_out(:,:,:,:,:) !< output data
  real(r8_kind),       intent(inout) :: weight_sum          !< Sum of weights from the output buffer object
  real(FMS_TRM_KIND_),       intent(in)    :: data_in(:,:,:,:)    !< data to update the buffer with
  logical,                   intent(in)    :: mask(:,:,:,:)       !< mask
  logical,                   intent(in)    :: is_masked           !< .True. if the field is using a mask
  type(fmsDiagIbounds_type), intent(in)    :: bounds_in           !< indices indicating the correct portion
                                                                  !! of the input buffer
  type(fmsDiagIbounds_type), intent(in)    :: bounds_out          !< indices indicating the correct portion
                                                                  !! of the output buffer
  real(FMS_TRM_KIND_),       intent(in)    :: missing_value       !< Missing_value for data points that are masked
  logical,                   intent(in)    :: increase_counter    !< .True. if data has not been received for
                                                                  !! time, so the counter needs to be increased
  real(r8_kind),optional, intent(in) :: weight              !< Weight applied to data_in before added to data_out
                                                                  !! used for weighted averages, default 1.0
  real(FMS_TRM_KIND_),optional, intent(in) :: pow                 !< Used for pow reduction, adds field^pow to buffer

  integer :: is_in, ie_in, js_in, je_in, ks_in, ke_in       !< Starting and ending indices of each dimention for
                                                            !! the input buffer
  integer :: is_out, ie_out, js_out, je_out, ks_out, ke_out !< Starting and ending indices of each dimention for
                                                            !! the output buffer
  integer :: i, j, k, l !< For looping
  real(FMS_TRM_KIND_) :: weight_loc, pow_loc !< local copies of optional arguments
  integer, parameter  :: kindl = FMS_TRM_KIND_ !< real kind size as set by macro

  if(present(weight)) then
    weight_loc = weight
  else
    weight_loc = 1.0_kindl
  endif

  if(present(pow)) then
    pow_loc = weight
  else
    pow_loc = 1.0_kindl
  endif

  ! update with given weight for average before write
  if (increase_counter) weight_sum = weight_sum + weight_loc

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

  !> Seperated this loops for performance. If is_masked = .false. (i.e "mask" and "rmask" were never passed in)
  !! then mask will always be .True. so the if (mask) is redudant.
  ! TODO check if performance gain by not doing weight and pow if not needed
  if (is_masked) then
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          where (mask(is_in + i, js_in + j, ks_in + k, :))
            data_out(is_out + i, js_out + j, ks_out + k, :, 1) =           &
                       data_out(is_out + i, js_out + j, ks_out + k, :, 1)  &
                     + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_loc) ** pow_loc
          elsewhere
            data_out(is_out + i, js_out + j, ks_out + k, :, 1) = missing_value
          endwhere
        enddo
      enddo
    enddo
  else
    ! doesn't need to loop through l if no mask, just sums the 1d slices
    do k = 0, ke_out - ks_out
      do j = 0, je_out - js_out
        do i = 0, ie_out - is_out
          data_out(is_out + i, js_out + j, ks_out + k,  :, 1) =  &
             data_out(is_out + i, js_out + j, ks_out + k, :, 1)  &
            + (data_in(is_in +i, js_in + j, ks_in + k, :) * weight_loc) ** pow_loc
        enddo
      enddo
    enddo
  endif
end subroutine DO_TIME_SUM_UPDATE_

!> To be called with diag_send_complete, finishes reductions
!! Just divides the buffer by the counter array(which is just the sum of the weights used in the buffer's reduction)
!! TODO: change has_mask to an actual logical mask so we don't have to check for missing values
subroutine SUM_UPDATE_DONE_(out_buffer_data, weight_sum, reduction_method, missing_val, has_mask)
  real(FMS_TRM_KIND_), intent(inout) :: out_buffer_data(:,:,:,:,:) !< data buffer previously updated with
                                                                   !! do_time_sum_update
  real(r8_kind), intent(in)          :: weight_sum !< sum of weights for averaging, provided via argument to send data
  integer, intent(in)                :: reduction_method !< which reduction method to use, should be time_avg
  real(FMS_TRM_KIND_), intent(in)    :: missing_val !< missing value for masked elements
  logical, intent(in)                :: has_mask !< indicates if mask is used so missing values can be skipped
  !! TODO replace conditional in the `where` with passed in and ajusted mask from the original call
  !logical, optional, intent(in)      :: mask(:,:,:,:) !< logical mask from accept data call, if using one.
  !logical                            :: has_mask !< whether or not mask is present

  if ( has_mask ) then
    where(out_buffer_data(:,:,:,:,1) .ne. missing_val)
      out_buffer_data(:,:,:,:,1) = out_buffer_data(:,:,:,:,1) &
                                   / weight_sum
    endwhere
  else  !not mask variant
    out_buffer_data(:,:,:,:,1) = out_buffer_data(:,:,:,:,1) &
                               / weight_sum
  endif

end subroutine

