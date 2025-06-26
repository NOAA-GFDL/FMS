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
!> @defgroup platform_mod platform_mod
!> @ingroup platform
!> @brief Uses @ref fms_platform.h to define byte sizes for variable kinds
!! to be used in fms.

!> @addtogroup platform_mod
!> @{
module platform_mod
!platform-dependent settings
#include <fms_platform.h>
  public
  integer, parameter :: r16_kind=QUAD_KIND, r8_kind=DOUBLE_KIND, r4_kind=FLOAT_KIND, &
                        c8_kind=DOUBLE_KIND, c4_kind=FLOAT_KIND, &
                        l8_kind=LONG_KIND, l4_kind=INT_KIND, &
                        i8_kind=LONG_KIND, i4_kind=INT_KIND, i2_kind=SHORT_KIND, &
                        ptr_kind=POINTER_KIND
  integer, parameter :: FMS_PATH_LEN = FMS_MAX_PATH_LEN
  integer, parameter :: FMS_FILE_LEN = FMS_MAX_FILE_LEN
!could additionally define things like OS, compiler...: useful?
end module platform_mod
!> @}
! close documentation grouping
