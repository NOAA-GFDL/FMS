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
!> @addtogroup horiz_interp_type_mod
!> @{
!> @brief This statistics is for bilinear interpolation and spherical regrid.
 subroutine STATS_ ( dat, low, high, avg, miss, missing_value, mask )
 real(FMS_HI_KIND_),    intent(in)  :: dat(:,:)
 real(FMS_HI_KIND_),    intent(out) :: low, high, avg
 integer, intent(out) :: miss
 real(FMS_HI_KIND_), intent(in), optional :: missing_value
 real(FMS_HI_KIND_),    intent(in), optional :: mask(:,:)

 real(FMS_HI_KIND_) :: dsum, buffer_real(3)
 integer :: pe, root_pe, npes, p, buffer_int(2), npts
 integer, parameter :: kindl = FMS_HI_KIND_ !< compiled kind size

   pe = mpp_pe()
   root_pe = mpp_root_pe()
   npes = mpp_npes()

   dsum = 0.0_kindl
   miss = 0

   if (present(missing_value)) then
      miss = count(dat(:,:) == missing_value)
      low  = minval(dat(:,:), dat(:,:) /= missing_value)
      high = maxval(dat(:,:), dat(:,:) /= missing_value)
      dsum = sum(dat(:,:), dat(:,:) /= missing_value)
   else if(present(mask)) then
      miss = count(mask(:,:) <= 0.5_kindl )
      low  = minval(dat(:,:),mask=mask(:,:) > 0.5_kindl)
      high = maxval(dat(:,:),mask=mask(:,:) > 0.5_kindl)
      dsum = sum(dat(:,:), mask=mask(:,:) > 0.5_kindl)
   else
      miss = 0
      low  = minval(dat(:,:))
      high = maxval(dat(:,:))
      dsum = sum(dat(:,:))
   endif
   avg = 0.0_kindl

   npts = size(dat(:,:)) - miss
   if(pe == root_pe) then
      do p = 1, npes - 1  ! root_pe receive data from other pe
      ! Force use of "scalar", integer pointer mpp interface
         call mpp_recv(buffer_real(1),glen=3, from_pe=p+root_pe, tag=COMM_TAG_1)
         dsum = dsum + buffer_real(1)
         low  = min(low, buffer_real(2))
         high = max(high, buffer_real(3))
         call mpp_recv(buffer_int(1), glen=2, from_pe=p+root_pe, tag=COMM_TAG_2)
         miss = miss + buffer_int(1)
         npts = npts + buffer_int(2)
      enddo
      if(npts == 0) then
         print*, 'Warning: no points is valid'
      else
         avg = dsum/real(npts, FMS_HI_KIND_)
      endif
    else   ! other pe send data to the root_pe.
      buffer_real(1) = dsum
      buffer_real(2) = low
      buffer_real(3) = high
      ! Force use of "scalar", integer pointer mpp interface
      call mpp_send(buffer_real(1),plen=3,to_pe=root_pe, tag=COMM_TAG_1)
      buffer_int(1) = miss
      buffer_int(2) = npts
      call mpp_send(buffer_int(1), plen=2, to_pe=root_pe, tag=COMM_TAG_2)
    endif

    call mpp_sync_self()

    return

 end subroutine STATS_
