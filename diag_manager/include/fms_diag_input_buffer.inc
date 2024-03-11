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

subroutine APPEND_DATA_BUFFER_(mask_out, mask_in, data_out, data_in)
  logical, intent(inout) :: mask_out(:,:,:,:)
  logical, intent(in)    :: mask_in(:,:,:,:)
  real(FMS_TRM_KIND_), intent(inout) :: data_out(:,:,:,:)
  real(FMS_TRM_KIND_), intent(in) :: data_in(:,:,:,:)

  integer :: i, j, k, l

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

end subroutine

subroutine SUM_DATA_BUFFER_(mask, data_out, data_in, counter, var_is_masked)
  logical,  intent(in)    :: mask(:,:,:,:)
  real(FMS_TRM_KIND_), intent(inout) :: data_out(:,:,:,:)
  real(FMS_TRM_KIND_), intent(in)    :: data_in(:,:,:,:)
  integer,  intent(inout)    :: counter(:,:,:,:)
  logical,  intent(in)    :: var_is_masked

  if (var_is_masked) then
    where (mask)
      data_out = data_out + data_in
    endwhere
  else
    data_out = data_out + data_in
  endif

  counter = counter + 1
end subroutine SUM_DATA_BUFFER_