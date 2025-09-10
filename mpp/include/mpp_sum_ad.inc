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
!> @brief Adjoint summation routines for @ref mpp_mod

!> @addtogroup mpp_mod
!> @{

!#######################################################################

    !> Sums array a.
    !! when only first element is passed: this routine just converts to a call to MPP_SUM_
    subroutine MPP_SUM_SCALAR_AD_( a, pelist )
      MPP_TYPE_, intent(inout) :: a !< scalar value to sum over pelist
      integer, intent(in), optional :: pelist(:)
      MPP_TYPE_ :: b(1)

      b(1) = a
      if( debug )call mpp_error( NOTE, 'MPP_SUM_SCALAR_: calling MPP_SUM_ ...' )
      call MPP_SUM_AD_( b, 1, pelist )
      a = b(1)
      return
    end subroutine MPP_SUM_SCALAR_AD_

!#######################################################################
    !> Sums 2d array across pes
    subroutine MPP_SUM_2D_AD_( a, length, pelist )
      MPP_TYPE_, target, contiguous, intent(inout) :: a(:,:) !< 2d array to sum
      integer, intent(in) :: length !< amount of indices in given 2d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      MPP_TYPE_, pointer :: a1D(:)

      a1D(1:length) => a(:,:)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine MPP_SUM_2D_AD_

!#######################################################################
    !> Sums 3d array across pes
    subroutine MPP_SUM_3D_AD_( a, length, pelist )
      MPP_TYPE_, target, contiguous, intent(inout) :: a(:,:,:) !< 3d array to sum
      integer, intent(in) :: length !< amount of indices in given 3d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      MPP_TYPE_, pointer :: a1D(:)

      a1D(1:length) => a(:,:,:)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine MPP_SUM_3D_AD_

!#######################################################################
    !> Sums 4d array across pes
    subroutine MPP_SUM_4D_AD_( a, length, pelist )
      MPP_TYPE_, target, contiguous, intent(inout) :: a(:,:,:,:) !< 4d array to sum
      integer, intent(in) :: length !< amount of indices in given 4d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      MPP_TYPE_, pointer :: a1D(:)

      a1D(1:length) => a(:,:,:,:)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine MPP_SUM_4D_AD_

!#######################################################################
    !> Sums 5d array across pes
    subroutine MPP_SUM_5D_AD_( a, length, pelist )
      MPP_TYPE_, target, contiguous, intent(inout) :: a(:,:,:,:,:) !< 5d array to sum
      integer, intent(in) :: length !< amount of indices in given 5d array
      integer, intent(in), optional :: pelist(:) !< pelist to calculate sum across
      MPP_TYPE_, pointer :: a1D(:)

      a1D(1:length) => a(:,:,:,:,:)
      call mpp_sum_ad( a1D, length, pelist )

      return
    end subroutine MPP_SUM_5D_AD_
