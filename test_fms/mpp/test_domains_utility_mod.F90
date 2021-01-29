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
!> @author Miguel Zuniga
!> @brief A module with utility auxiliary interface supporting test_mpp_domains.
!> @note Note: the source code of this module is largely originally in test_mpp_domains.F90.
module test_domains_utility_mod
    use mpp_mod,         only : FATAL, WARNING, MPP_DEBUG, NOTE
    use mpp_mod, only : mpp_error
    use mpp_domains_mod, only : ZERO, NINETY, MINUS_NINETY
    use platform_mod, only: r4_kind, r8_kind

  interface fill_coarse_data
    module procedure fill_coarse_data_r8
    module procedure fill_coarse_data_r4
  end interface fill_coarse_data

  interface fill_nest_data
    module procedure fill_nest_data_r8
    module procedure fill_nest_data_r4
  end interface fill_nest_data

  contains

subroutine fill_coarse_data_r8(data, rotate, iadd, jadd, is_c, ie_c, js_c, je_c, nz, isd, jsd, nx, ny, &
                              ishift, jshift, x_add, y_add, sign1, sign2, x_cyclic, y_cyclic, ieg, jeg)
    integer, intent(in)    :: rotate, is_c, ie_c, js_c, je_c, nz, isd, jsd, iadd, jadd, nx, ny, ishift, jshift
    integer, intent(in)    :: sign1, sign2
    real(kind=r8_kind),    intent(inout) :: data(isd:, jsd:, :)
    real(kind=r8_kind),    intent(in)    :: x_add, y_add
    logical, intent(in)    :: x_cyclic, y_cyclic
    integer, intent(in)    :: ieg, jeg
    integer :: i, j, k

    select case (rotate)
    case (ZERO)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = (i+iadd)*1.e+6 + (j+jadd)*1.e+3 + k + x_add
             enddo
          enddo
       enddo
    case (NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign1*((nx-j+1+iadd+jshift)*1.e+6 + (i+jadd)*1.e+3 + k + y_add)
             enddo
          enddo
       enddo
    case (MINUS_NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign2*((j+iadd)*1.e+6 + (ny-i+1+jadd+ishift)*1.e+3 + k + y_add)
             enddo
          enddo
       enddo
    case default
       call mpp_error(FATAL,"fill_coarse_data: rotate_coarse must be ZERO, NINETY, MINUS_NINETY")
    end select

    !---handle cyclic condition
    if(x_cyclic) then
       if(ie_c+ishift+iadd == ieg) then
          i = ie_c+ishift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = i*1.e+6 + (j+jadd)*1.e+3 + k + x_add
             enddo
          enddo
       endif
    endif


    if(y_cyclic) then
       if(je_c+jshift+jadd == jeg) then
          j = je_c+jshift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = (i+iadd)*1.e+6 + j*1.e+3 + k + x_add
             enddo
          enddo
       endif
    endif

  end subroutine fill_coarse_data_r8


subroutine fill_coarse_data_r4(data, rotate, iadd, jadd, is_c, ie_c, js_c, je_c, nz, isd, jsd, nx, ny, &
                              ishift, jshift, x_add, y_add, sign1, sign2, x_cyclic, y_cyclic, ieg, jeg)
    integer, intent(in)    :: rotate, is_c, ie_c, js_c, je_c, nz, isd, jsd, iadd, jadd, nx, ny, ishift, jshift
    integer, intent(in)    :: sign1, sign2
    real(kind=r4_kind),    intent(inout) :: data(isd:, jsd:, :)
    real(kind=r4_kind),    intent(in)    :: x_add, y_add
    logical, intent(in)    :: x_cyclic, y_cyclic
    integer, intent(in)    :: ieg, jeg
    integer :: i, j, k

    select case (rotate)
    case (ZERO)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = (i+iadd)*1.e+6 + (j+jadd)*1.e+3 + k + x_add
             enddo
          enddo
       enddo
    case (NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign1*((nx-j+1+iadd+jshift)*1.e+6 + (i+jadd)*1.e+3 + k + y_add)
             enddo
          enddo
       enddo
    case (MINUS_NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign2*((j+iadd)*1.e+6 + (ny-i+1+jadd+ishift)*1.e+3 + k + y_add)
             enddo
          enddo
       enddo
    case default
       call mpp_error(FATAL,"fill_coarse_data: rotate_coarse must be ZERO, NINETY, MINUS_NINETY")
    end select

    !---handle cyclic condition
    if(x_cyclic) then
       if(ie_c+ishift+iadd == ieg) then
          i = ie_c+ishift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = i*1.e+6 + (j+jadd)*1.e+3 + k + x_add
             enddo
          enddo
       endif
    endif


    if(y_cyclic) then
       if(je_c+jshift+jadd == jeg) then
          j = je_c+jshift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = (i+iadd)*1.e+6 + j*1.e+3 + k + x_add
             enddo
          enddo
       endif
    endif

  end subroutine fill_coarse_data_r4

  !###########################################################################################

    subroutine fill_nest_data_r8(buffer, is, ie, js, je, nnest, tile, ishift, jshift, iadd, jadd, rotate, &
                            isl, iel, jsl, jel, xadd, yadd, sign1, sign2, nx, ny)
     real(kind=r8_kind), dimension(is:,js:,:), intent(inout) :: buffer
     integer,                       intent(in) :: is, ie, js, je, nnest
     integer,                       intent(in) :: ishift, jshift
     integer, dimension(:),         intent(in) :: tile, iadd, jadd, rotate, isl, iel, jsl, jel
     real(kind=r8_kind),                          intent(in) :: xadd, yadd
     integer,                       intent(in) :: sign1, sign2
     integer,                       intent(in) :: nx, ny
     integer :: i, j, k, n, nk
     integer :: ioff, joff

     ioff = 0
     joff = 0
     nk = size(buffer,3)
     do k = 1, nk
        do n = 1, nnest
           if(iel(n) == ie) ioff = ishift
           if(jel(n) == je) joff = jshift

           select case (rotate(n))
           case(ZERO)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = xadd + tile(n) + (i-iadd(n))*1.e-3 + (j-jadd(n))*1.e-6 + k*1.e-9
                 enddo
              enddo
           case (NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign2*(yadd + tile(n) + (j-jadd(n))*1.e-3 + (nx-i+iadd(n)+1+ioff)*1.e-6 + k*1.e-9)
                 enddo
              enddo
           case (MINUS_NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign1*(yadd + tile(n) + (ny-j+jadd(n)+1+joff)*1.e-3 + (i-iadd(n))*1.e-6 + k*1.e-9)
                 enddo
              enddo
           case default
              call mpp_error(FATAL,"fill_nest_data: rotate must be ZERO, NINETY, MINUS_NINETY")
           end select
        enddo
     enddo

   end subroutine fill_nest_data_r8

   !###########################################################################################

    subroutine fill_nest_data_r4(buffer, is, ie, js, je, nnest, tile, ishift, jshift, iadd, jadd, rotate, &
                            isl, iel, jsl, jel, xadd, yadd, sign1, sign2, nx, ny)
     real(kind=r4_kind), dimension(is:,js:,:), intent(inout) :: buffer
     integer,                       intent(in) :: is, ie, js, je, nnest
     integer,                       intent(in) :: ishift, jshift
     integer, dimension(:),         intent(in) :: tile, iadd, jadd, rotate, isl, iel, jsl, jel
     real(kind=r4_kind),                          intent(in) :: xadd, yadd
     integer,                       intent(in) :: sign1, sign2
     integer,                       intent(in) :: nx, ny
     integer :: i, j, k, n, nk
     integer :: ioff, joff

     ioff = 0
     joff = 0
     nk = size(buffer,3)
     do k = 1, nk
        do n = 1, nnest
           if(iel(n) == ie) ioff = ishift
           if(jel(n) == je) joff = jshift

           select case (rotate(n))
           case(ZERO)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = xadd + tile(n) + (i-iadd(n))*1.e-3 + (j-jadd(n))*1.e-6 + k*1.e-9
                 enddo
              enddo
           case (NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign2*(yadd + tile(n) + (j-jadd(n))*1.e-3 + (nx-i+iadd(n)+1+ioff)*1.e-6 + k*1.e-9)
                 enddo
              enddo
           case (MINUS_NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign1*(yadd + tile(n) + (ny-j+jadd(n)+1+joff)*1.e-3 + (i-iadd(n))*1.e-6 + k*1.e-9)
                 enddo
              enddo
           case default
              call mpp_error(FATAL,"fill_nest_data: rotate must be ZERO, NINETY, MINUS_NINETY")
           end select
        enddo
     enddo

  end subroutine fill_nest_data_r4


end module test_domains_utility_mod
