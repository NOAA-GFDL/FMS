! -*-f90-*-

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
!> @file
!> @brief Routines for data transmission between PE's

!> @addtogroup mpp_mod
!> @{

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!                                MPP_BROADCAST                                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MPP_BROADCAST_SCALAR_( broadcast_data, from_pe, pelist )
      MPP_TYPE_, intent(inout) :: broadcast_data
      integer, intent(in) :: from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(1)

      pointer( ptr, data1D )

      ptr = LOC(broadcast_data)
      call MPP_BROADCAST_( data1D, 1, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_SCALAR_

    subroutine MPP_BROADCAST_2D_( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: broadcast_data(:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_2D_

    subroutine MPP_BROADCAST_3D_( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: broadcast_data(:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
   end subroutine MPP_BROADCAST_3D_

    subroutine MPP_BROADCAST_4D_( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: broadcast_data(:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_4D_

    subroutine MPP_BROADCAST_5D_( broadcast_data, length, from_pe, pelist )
!this call was originally bundled in with mpp_transmit, but that doesn't allow
!broadcast to a subset of PEs. This version will, and mpp_transmit will remain
!backward compatible.
      MPP_TYPE_, intent(inout) :: broadcast_data(:,:,:,:,:)
      integer, intent(in) :: length, from_pe
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: data1D(length)

      pointer( ptr, data1D )
      ptr = LOC(broadcast_data)
      call mpp_broadcast( data1D, length, from_pe, pelist )

      return
    end subroutine MPP_BROADCAST_5D_
!> @}
