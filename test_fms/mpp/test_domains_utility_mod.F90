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
    use mpp_domains_mod, only : ZERO, NINETY, MINUS_NINETY, &
                                domain2d, mpp_define_mosaic
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
                data(i,j,k) = dble(i+iadd)*1.d+6 + dble(j+jadd)*1.d+3 + dble(k) + x_add
             enddo
          enddo
       enddo
    case (NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign1*( dble(nx-j+1+iadd+jshift)*1.d+6 + dble(i+jadd)*1.d+3 + dble(k) + y_add)
             enddo
          enddo
       enddo
    case (MINUS_NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign2*( dble(j+iadd)*1.d+6 + dble(ny-i+1+jadd+ishift)*1.d+3 + dble(k) + y_add)
             enddo
          enddo
       enddo
    case default
       call mpp_error(FATAL,"fill_coarse_data_r8: rotate_coarse must be ZERO, NINETY, MINUS_NINETY")
    end select

    !---handle cyclic condition
    if(x_cyclic) then
       if(ie_c+ishift+iadd == ieg) then
          i = ie_c+ishift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = dble(i)*1.d+6 + dble(j+jadd)*1.d+3 + dble(k) + x_add
             enddo
          enddo
       endif
    endif


    if(y_cyclic) then
       if(je_c+jshift+jadd == jeg) then
          j = je_c+jshift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = dble(i+iadd)*1.d+6 + j*1.d+3 + dble(k) + x_add
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
       call mpp_error(FATAL,"fill_coarse_data_r4: rotate_coarse must be ZERO, NINETY, MINUS_NINETY")
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
                    buffer(i,j,k) = xadd + dble(tile(n)) + dble(i-iadd(n))*1.d-3 &
                                  & + dble(j-jadd(n))*1.d-6 + dble(k)*1.d-9
                 enddo
              enddo
           case (NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign2*(yadd + dble(tile(n)) + dble(j-jadd(n))*1.d-3 &
                                  & + dble(nx-i+iadd(n)+1+ioff)*1.d-6 + dble(k)*1.d-9)
                 enddo
              enddo
           case (MINUS_NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign1*(yadd + dble(tile(n)) + dble(ny-j+jadd(n)+1+joff)*1.d-3 &
                                  & + dble(i-iadd(n))*1.d-6 + dble(k)*1.d-9)
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

  !#######################################################################################
  !--- define mosaic domain for cubic grid
  subroutine define_cubic_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end, &
                                 use_memsize, halo)
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    logical, optional, intent(in) :: use_memsize
    integer, optional, intent(in) :: halo
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)
    logical                       :: use_memsize_local
    integer :: whalo = 2, shalo = 2, ehalo = 2, nhalo = 2

    use_memsize_local = .true.
    if(present(use_memsize)) use_memsize_local = use_memsize
    if(present(halo)) then
      nhalo = halo; shalo = halo; whalo = halo; shalo = halo
    endif

    ntiles = 6
    num_contact = 12
    if(size(pe_start(:)) .NE. 6 .OR. size(pe_end(:)) .NE. 6 ) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of pe_start and pe_end should be 6")
    if(size(global_indices,1) .NE. 4) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of first dimension of global_indices should be 4")
    if(size(global_indices,2) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of second dimension of global_indices should be 6")
    if(size(layout,1) .NE. 2) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of first dimension of layout should be 2")
    if(size(layout,2) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of second dimension of layout should be 6")
    if(size(ni(:)) .NE. 6 .OR. size(nj(:)) .NE. 6) call mpp_error(FATAL, &
         "define_cubic_mosaic: size of ni and nj should be 6")

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1; tile2(1) = 2
    istart1(1) = ni(1);  iend1(1) = ni(1);  jstart1(1) = 1;      jend1(1) = nj(1)
    istart2(1) = 1;      iend2(1) = 1;      jstart2(1) = 1;      jend2(1) = nj(2)
    !--- Contact line 2, between tile 1 (NORTH) and tile 3 (WEST)
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;      iend1(2) = ni(1);  jstart1(2) = nj(1);  jend1(2) = nj(1)
    istart2(2) = 1;      iend2(2) = 1;      jstart2(2) = nj(3);  jend2(2) = 1
    !--- Contact line 3, between tile 1 (WEST) and tile 5 (NORTH)
    tile1(3) = 1; tile2(3) = 5
    istart1(3) = 1;      iend1(3) = 1;      jstart1(3) = 1;      jend1(3) = nj(1)
    istart2(3) = ni(5);  iend2(3) = 1;      jstart2(3) = nj(5);  jend2(3) = nj(5)
    !--- Contact line 4, between tile 1 (SOUTH) and tile 6 (NORTH)
    tile1(4) = 1; tile2(4) = 6
    istart1(4) = 1;      iend1(4) = ni(1);  jstart1(4) = 1;      jend1(4) = 1
    istart2(4) = 1;      iend2(4) = ni(6);  jstart2(4) = nj(6);  jend2(4) = nj(6)
    !--- Contact line 5, between tile 2 (NORTH) and tile 3 (SOUTH)
    tile1(5) = 2; tile2(5) = 3
    istart1(5) = 1;      iend1(5) = ni(2);  jstart1(5) = nj(2);  jend1(5) = nj(2)
    istart2(5) = 1;      iend2(5) = ni(3);  jstart2(5) = 1;      jend2(5) = 1
    !--- Contact line 6, between tile 2 (EAST) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = ni(2);  iend1(6) = ni(2);  jstart1(6) = 1;      jend1(6) = nj(2)
    istart2(6) = ni(4);  iend2(6) = 1;      jstart2(6) = 1;      jend2(6) = 1
    !--- Contact line 7, between tile 2 (SOUTH) and tile 6 (EAST)
    tile1(7) = 2; tile2(7) = 6
    istart1(7) = 1;      iend1(7) = ni(2);  jstart1(7) = 1;      jend1(7) = 1
    istart2(7) = ni(6);  iend2(7) = ni(6);  jstart2(7) = nj(6);  jend2(7) = 1
    !--- Contact line 8, between tile 3 (EAST) and tile 4 (WEST)
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = ni(3);  iend1(8) = ni(3);  jstart1(8) = 1;      jend1(8) = nj(3)
    istart2(8) = 1;      iend2(8) = 1;      jstart2(8) = 1;      jend2(8) = nj(4)
    !--- Contact line 9, between tile 3 (NORTH) and tile 5 (WEST)
    tile1(9) = 3; tile2(9) = 5
    istart1(9) = 1;      iend1(9) = ni(3);  jstart1(9) = nj(3);  jend1(9) = nj(3)
    istart2(9) = 1;      iend2(9) = 1;      jstart2(9) = nj(5);  jend2(9) = 1
    !--- Contact line 10, between tile 4 (NORTH) and tile 5 (SOUTH)
    tile1(10) = 4; tile2(10) = 5
    istart1(10) = 1;     iend1(10) = ni(4); jstart1(10) = nj(4); jend1(10) = nj(4)
    istart2(10) = 1;     iend2(10) = ni(5); jstart2(10) = 1;     jend2(10) = 1
    !--- Contact line 11, between tile 4 (EAST) and tile 6 (SOUTH)
    tile1(11) = 4; tile2(11) = 6
    istart1(11) = ni(4); iend1(11) = ni(4); jstart1(11) = 1;     jend1(11) = nj(4)
    istart2(11) = ni(6); iend2(11) = 1;     jstart2(11) = 1;     jend2(11) = 1
    !--- Contact line 12, between tile 5 (EAST) and tile 6 (WEST)
    tile1(12) = 5; tile2(12) = 6
    istart1(12) = ni(5); iend1(12) = ni(5); jstart1(12) = 1;     jend1(12) = nj(5)
    istart2(12) = 1;     iend2(12) = 1;     jstart2(12) = 1;     jend2(12) = nj(6)
    msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1 ! make sure memory domain size is no smaller than
    msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1 ! data domain size

    if(use_memsize_local) then
       call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2, &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
         pe_start, pe_end, symmetry = .true., whalo=whalo, ehalo=ehalo,   &
         shalo=shalo, nhalo=nhalo, name = trim(type), memory_size = msize  )
    else
       call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2, &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
         pe_start, pe_end, symmetry = .true., whalo=whalo, ehalo=ehalo,   &
         shalo=shalo, nhalo=nhalo, name = trim(type) )
    endif

    return

  end subroutine define_cubic_mosaic

  !######################################################################################
  subroutine define_fourtile_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end, symmetry, halo)
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    logical,        intent(in)    :: symmetry
    integer, dimension(8)         :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(8)         :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)
    integer, optional, intent(in) :: halo
    integer :: whalo = 2, shalo = 2, ehalo = 2, nhalo = 2

    if(present(halo)) then
      nhalo = halo; shalo = halo; whalo = halo; shalo = halo
    endif

    ntiles = 4
    num_contact = 8
    if(size(pe_start(:)) .NE. 4 .OR. size(pe_end(:)) .NE. 4 ) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of pe_start and pe_end should be 4")
    if(size(global_indices,1) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of first dimension of global_indices should be 4")
    if(size(global_indices,2) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of second dimension of global_indices should be 4")
    if(size(layout,1) .NE. 2) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of first dimension of layout should be 2")
    if(size(layout,2) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of second dimension of layout should be 4")
    if(size(ni(:)) .NE. 4 .OR. size(nj(:)) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of ni and nj should be 4")

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1; tile2(1) = 2
    istart1(1) = ni(1); iend1(1) = ni(1); jstart1(1) = 1;     jend1(1) = nj(1)
    istart2(1) = 1;     iend2(1) = 1;     jstart2(1) = 1;     jend2(1) = nj(2)
    !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;     iend1(2) = ni(1); jstart1(2) = 1;     jend1(2) = 1
    istart2(2) = 1;     iend2(2) = ni(3); jstart2(2) = nj(3); jend2(2) = nj(3)
    !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
    tile1(3) = 1; tile2(3) = 2
    istart1(3) = 1;     iend1(3) = 1;     jstart1(3) = 1;     jend1(3) = nj(1)
    istart2(3) = ni(2); iend2(3) = ni(2); jstart2(3) = 1;     jend2(3) = nj(2)
    !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
    tile1(4) = 1; tile2(4) = 3
    istart1(4) = 1;     iend1(4) = ni(1); jstart1(4) = nj(1); jend1(4) = nj(1)
    istart2(4) = 1;     iend2(4) = ni(3); jstart2(4) = 1;     jend2(4) = 1
    !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
    tile1(5) = 2; tile2(5) = 4
    istart1(5) = 1;     iend1(5) = ni(2); jstart1(5) = 1;     jend1(5) = 1
    istart2(5) = 1;     iend2(5) = ni(4); jstart2(5) = nj(4); jend2(5) = nj(4)
    !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = 1;     iend1(6) = ni(2); jstart1(6) = nj(2); jend1(6) = nj(2)
    istart2(6) = 1;     iend2(6) = ni(4); jstart2(6) = 1;     jend2(6) = 1
    !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
    tile1(7) = 3; tile2(7) = 4
    istart1(7) = ni(3); iend1(7) = ni(3); jstart1(7) = 1;     jend1(7) = nj(3)
    istart2(7) = 1;     iend2(7) = 1;     jstart2(7) = 1;     jend2(7) = nj(4)
    !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = 1;     iend1(8) = 1;     jstart1(8) = 1;     jend1(8) = nj(3)
    istart2(8) = ni(4); iend2(8) = ni(4); jstart2(8) = 1;     jend2(8) = nj(4)
    msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1 ! make sure memory domain size is no smaller than
    msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1 ! data domain size
    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2,       &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,          &
         pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,    &
         name = type, memory_size = msize, symmetry = symmetry )

    return

  end subroutine define_fourtile_mosaic

end module test_domains_utility_mod
