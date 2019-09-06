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

program test_drifters_io

  use drifters_io_mod
  use mpp_mod, only : mpp_error, FATAL, stdout, mpp_init, mpp_exit

  implicit none
  type(drifters_io_type) :: drfts_io
  character(len=128) :: ermesg
  character(len=31) :: filename
  integer :: np, nd, nf, nt, i, j, k, npmax
  real :: dt, time, xmin, xmax, ymin, ymax, u, v, dr, x, y
  integer, allocatable :: ids(:)
  real, allocatable :: positions(:,:), fields(:,:)

  call mpp_init()

  ! number of dimensions
  nd = 3
  ! number of fields
  nf = 2
  ! max number of dirfters
  npmax = 20
  ! number of time steps
  nt = 50
  ! starting time
  time = 0.

  ! domain boundary. (drifters outside domain will not be written to file.)
  xmin = 0.
  ymin = 0.
  xmax = 1.
  ymax = 1.

  ! constant velocity
  u = (xmax-xmin)*sqrt(2.)
  v = (ymax-ymin)*sqrt(2.)
  dt = 1/real(nt)

  ! open file

  filename = 'test_drifter_io.nc'
  call drifters_io_new(drfts_io, filename, nd, nf, ermesg)
  if(ermesg/='') then
     print *,'ERROR after drifters_io_new: ', ermesg
     call mpp_error(FATAL, ermesg)
  endif
  ! set attributes

  call drifters_io_set_position_names(drfts_io, (/'x','y','z'/), ermesg)
  if(ermesg/='') then
     print *,'ERROR after drifters_io_position_names: ', ermesg
     call mpp_error(FATAL, ermesg)
  endif

  ! note the trailing blanks in the first field, which are added here to
  ! ensure that "salinity" will not be truncated (all names must have the
  ! same length)
  call drifters_io_set_field_names(drfts_io, (/'temp    ','salinity'/), ermesg)
  if(ermesg/='') then
      print *,'ERROR after drifters_io_field_names: ', ermesg
      call mpp_error(FATAL, ermesg)
  endif

  call drifters_io_set_position_units(drfts_io, (/'deg east ','deg north','meters   '/), ermesg)
  if(ermesg/='') then
      print *,'ERROR after drifters_io_position_units: ', ermesg
      call mpp_error(FATAL, ermesg)
  endif

  call drifters_io_set_field_units(drfts_io, (/'deg K','ppm  '/), ermesg)
  if(ermesg/='') then
      print *,'ERROR after drifters_io_field_units: ', ermesg
      call mpp_error(FATAL, ermesg)
  endif

  allocate(positions(nd, npmax), ids(npmax), fields(nf, npmax))
  dr = sqrt( (xmax-xmin)**2 + (ymax-ymin)**2 )/real(npmax)


  ! x
  positions(1, :) = +(/ (i*dr,i=0,npmax-1) /)/sqrt(2.)
  ! y
  positions(2, :) = -(/ (i*dr,i=0,npmax-1) /)/sqrt(2.)
  ! z
  positions(3, :) = 0.

  ! drifters' identity array (can be any integer number)
  ids = (/ (i, i=1, npmax) /)

  ! set fields as a function of space time
  fields(1, :) = sqrt( (positions(1,:)-xmin)**2 + (positions(2,:)-ymin)**2 )
  fields(2, :) = positions(1,:)-u*time + positions(2,:)-v*time ! invariant

  ! write to disk only drifters inside domain
  do i = 1, npmax
     x = positions(1,i)
     y = positions(2,i)
     if(x>=xmin .and. x<=xmax .and. y>=ymin .and. y<=ymax) then
        call drifters_io_write(drfts_io, time, np=1, nd=nd, nf=nf, &
             & ids=ids(i), positions=positions(:,i), fields=fields(:,i), ermesg=ermesg)
        if(ermesg/='') then
            print *,'ERROR after drifters_io_write: ', ermesg
            call mpp_error(FATAL, ermesg)
        endif
     endif
  enddo

  ! advect

  do j = 1, nt
     time = time + dt
     positions(1, :) = positions(1, :) + u*dt
     positions(2, :) = positions(2, :) + v*dt
     fields(1, :) = sqrt( (positions(1,:)-xmin)**2 + (positions(2,:)-ymin)**2 )
     fields(2, :) = positions(1,:)-u*time + positions(2,:)-v*time ! invariant

     do i = 1, npmax
        x = positions(1,i)
        y = positions(2,i)
        if(x>=xmin .and. x<=xmax .and. y>=ymin .and. y<=ymax) then
           call drifters_io_write(drfts_io, time, np=1, nd=nd, nf=nf, &
                & ids=ids(i), positions=positions(:,i), fields=fields(:,i), ermesg=ermesg)
           if(ermesg/='') then
              print *,'ERROR after drifters_io_write: ', ermesg
              call mpp_error(FATAL, ermesg)
           endif
        endif
     enddo

  enddo

  deallocate(positions, ids, fields)

  call drifters_io_del(drfts_io, ermesg)
  if(ermesg/='') then
      print *,'ERROR after drifters_io_del: ', ermesg
      call mpp_error(FATAL, ermesg)
  endif
  call mpp_exit()
end program test_drifters_io
