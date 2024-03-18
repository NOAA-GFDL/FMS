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