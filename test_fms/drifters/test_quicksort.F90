
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

      program test_quicksort
#ifdef use_drifters
        implicit none
        integer :: list(16) = (/6, 2, 3, 4, 1, 45, 3432, 3245, 32545, 66555, 32, 1,3, -43254, 324, 54/)
        print *,'before list=', list
        call qksrt_quicksort(size(list), list, 1, size(list))
        print *,'after  list=', list
#endif
      end program test_quicksort
