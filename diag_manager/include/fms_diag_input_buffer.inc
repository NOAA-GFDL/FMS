!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************

!> @brief Appends the input_data_buffer and the mask (only when the mask is set to .True.)
subroutine APPEND_DATA_BUFFER_(mask_out, mask_in, data_out, data_in)
  logical,             intent(inout) :: mask_out(:,:,:,:) !< Mask currently in the input_data_buffer
  logical,             intent(in)    :: mask_in(:,:,:,:)  !< Mask passed in to send_data
  real(FMS_TRM_KIND_), intent(inout) :: data_out(:,:,:,:) !< Data currently in the input_data_buffer
  real(FMS_TRM_KIND_), intent(in)    :: data_in(:,:,:,:)  !< Data passed in to send_data

  integer :: i, j, k, l !< For looping through the input_data_buffer

  do l = 1, size(data_out, 4)
    do k = 1, size(data_out, 3)
      do j = 1, size(data_out, 2)
        do i = 1, size(data_out, 1)
          if (mask_in(i,j,k,l)) then
            mask_out(i,j,k,l) = .True.
            data_out(i,j,k,l) = data_in(i,j,k,l)
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine APPEND_DATA_BUFFER_

!> @brief Sums the data in the input_data_buffer
subroutine SUM_DATA_BUFFER_(mask, data_out, data_in, counter, var_is_masked)
  logical,             intent(in)    :: mask(:,:,:,:)     !< Mask passed into send_data
  real(FMS_TRM_KIND_), intent(inout) :: data_out(:,:,:,:) !< Data currently saved in the input_data_buffer
  real(FMS_TRM_KIND_), intent(in)    :: data_in(:,:,:,:)  !< Data passed into send_data
  integer,             intent(inout) :: counter(:,:,:,:)  !< Number of times data has been summed
  logical,             intent(in)    :: var_is_masked     !< .True. if the variable is masked

  if (var_is_masked) then
    where (mask)
      data_out = data_out + data_in
    endwhere
  else
    data_out = data_out + data_in
  endif

  counter = counter + 1
end subroutine SUM_DATA_BUFFER_