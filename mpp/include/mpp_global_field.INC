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
    subroutine MPP_GLOBAL_FIELD_2D_( domain, local, global, flags, position,tile_count, default_data)
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:)
      MPP_TYPE_, intent(out) :: global(:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count
      MPP_TYPE_, intent(in), optional :: default_data

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),1)
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),1)
      pointer( lptr,  local3D )
      pointer( gptr, global3D )
      lptr = LOC( local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, position,tile_count, default_data )

    end subroutine MPP_GLOBAL_FIELD_2D_

    subroutine MPP_GLOBAL_FIELD_3D_( domain, local, global, flags, position, tile_count, default_data)
!get a global field from a local field
!local field may be on compute OR data domain
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count
      MPP_TYPE_, intent(in), optional :: default_data

      integer :: ishift, jshift
      integer :: tile
      integer :: isize, jsize

      tile = 1; if(PRESENT(tile_count)) tile = tile_count

      call mpp_get_domain_shift(domain, ishift, jshift, position)

      ! The alltoallw method requires that local and global be contiguous.
      ! We presume that `local` is contiguous if it matches the data domain;
      ! `global` is presumed to always be contiguous.
      ! Ideally we would use the F2015 function IS_CONTIGUOUS() to validate
      !     contiguity, but it is not yet suppored in many compilers.

      ! Also worth noting that many of the nD->3D conversion also assumes
      ! contiguity, so there many be other issues here.

      isize = domain%x(tile)%data%size + ishift
      jsize = domain%y(tile)%data%size + jshift
      if ((size(local, 1) .eq. isize) .and. (size(local, 2) .eq. jsize) &
              .and. use_alltoallw) then
          call mpp_do_global_field_a2a(domain, local, global, tile, &
                                       ishift, jshift, flags, default_data)
      else
          call mpp_do_global_field(domain, local, global, tile, &
                                   ishift, jshift, flags, default_data)
      end if
    end subroutine MPP_GLOBAL_FIELD_3D_

    subroutine MPP_GLOBAL_FIELD_4D_( domain, local, global, flags, position,tile_count, default_data )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count
      MPP_TYPE_, intent(in), optional :: default_data

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size(local,4))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(local,4))
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, position,tile_count, default_data )
    end subroutine MPP_GLOBAL_FIELD_4D_

    subroutine MPP_GLOBAL_FIELD_5D_( domain, local, global, flags, position,tile_count, default_data )
      type(domain2D), intent(in) :: domain
      MPP_TYPE_, intent(in)  ::  local(:,:,:,:,:)
      MPP_TYPE_, intent(out) :: global(:,:,:,:,:)
      integer, intent(in), optional :: flags
      integer, intent(in), optional :: position
      integer, intent(in), optional :: tile_count
      MPP_TYPE_, intent(in), optional :: default_data

      MPP_TYPE_ :: local3D (size( local,1),size( local,2),size( local,3)*size( local,4)*size(local,5))
      MPP_TYPE_ :: global3D(size(global,1),size(global,2),size(global,3)*size(global,4)*size(local,5))
      pointer( lptr, local3D  )
      pointer( gptr, global3D )
      lptr = LOC(local)
      gptr = LOC(global)
      call mpp_global_field( domain, local3D, global3D, flags, position,tile_count, default_data )
    end subroutine MPP_GLOBAL_FIELD_5D_
