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
!> @author Jessica Liptak
!> @brief This module contains routines to fill halos in different domain configurations
!! It is required by test_mpp_update_domains_real and test_mpp_update_domains_int.
module fill_halo

use :: platform_mod

implicit none
private
integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
integer :: nx=64, ny=64, nz=10

public :: fill_halo_zero, fill_regular_refinement_halo, fill_regular_mosaic_halo
public :: fill_folded_north_halo, fill_folded_south_halo, fill_folded_east_halo, fill_folded_west_halo

!> Routines to fill halo regions of 64-bit and 32-bit real arrays on a regular grid
interface fill_regular_refinement_halo
  module procedure fill_regular_refinement_halo_r8
  module procedure fill_regular_refinement_halo_r4
  module procedure fill_regular_refinement_halo_i8
  module procedure fill_regular_refinement_halo_i4
end interface

!> Routines to fill halo regions of 64-bit and 32-bit real arrays with zeros
interface fill_halo_zero
  module procedure fill_halo_zero_r8
  module procedure fill_halo_zero_r4
  module procedure fill_halo_zero_i8
  module procedure fill_halo_zero_i4
end interface

!> Routines to fill halo regions of 64-bit and 32-bit real arrays on a mosaic grid
interface fill_regular_mosaic_halo
  module procedure fill_regular_mosaic_halo_r8
  module procedure fill_regular_mosaic_halo_r4
  module procedure fill_regular_mosaic_halo_i8
  module procedure fill_regular_mosaic_halo_i4
end interface fill_regular_mosaic_halo

!> Routines to fill halo regions of 64-bit and 32-bit real arrays on a domain with a folded north edge
interface fill_folded_north_halo
  module procedure fill_folded_north_halo_r8
  module procedure fill_folded_north_halo_r4
  module procedure fill_folded_north_halo_i8
  module procedure fill_folded_north_halo_i4
end interface fill_folded_north_halo

!> Routines to fill halo regions of 64-bit and 32-bit real arrays on a domain with a folded south edge
interface fill_folded_south_halo
  module procedure fill_folded_south_halo_r8
  module procedure fill_folded_south_halo_r4
  module procedure fill_folded_south_halo_i8
  module procedure fill_folded_south_halo_i4
end interface fill_folded_south_halo

!> Routines to fill halo regions of 64-bit and 32-bit real arrays on a domain with a folded east edge
interface fill_folded_east_halo
  module procedure fill_folded_east_halo_r8
  module procedure fill_folded_east_halo_r4
  module procedure fill_folded_east_halo_i8
  module procedure fill_folded_east_halo_i4
end interface fill_folded_east_halo

!> Routines to fill halo regions of 64-bit and 32-bit real arrays on a domain with a folded west edge
interface fill_folded_west_halo
  module procedure fill_folded_west_halo_r8
  module procedure fill_folded_west_halo_r4
  module procedure fill_folded_west_halo_i8
  module procedure fill_folded_west_halo_i4
end interface fill_folded_west_halo

contains

  !> fill the halo region of a 64-bit real array with zeros
  subroutine fill_halo_zero_r8(data, whalo, ehalo, shalo, nhalo, xshift, yshift, isc, iec, jsc, jec, isd, ied, jsd, jed)
    real(kind=r8_kind), dimension(isd:,jsd:,:), intent(inout) :: data
    integer,                         intent(in) :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer,                         intent(in) :: whalo, ehalo, shalo, nhalo, xshift, yshift

    if(whalo >=0) then
       data(iec+ehalo+1+xshift:ied+xshift,jsd:jed+yshift,:) = 0
       data(isd:isc-whalo-1,jsd:jed+yshift,:) = 0
    else
       data(iec+1+xshift:iec-ehalo+xshift,jsc+shalo:jec-nhalo+yshift,:) = 0
       data(isc+whalo:isc-1,jsc+shalo:jec-nhalo+yshift,:) = 0
    end if

    if(shalo>=0) then
       data(isd:ied+xshift, jec+nhalo+1+yshift:jed+yshift,:) = 0
       data(isd:ied+xshift, jsd:jsc-shalo-1,:) = 0
    else
       data(isc+whalo:iec-ehalo+xshift,jec+1+yshift:jec-nhalo+yshift,:) = 0
       data(isc+whalo:iec-ehalo+xshift,jsc+shalo:jsc-1,:) = 0
    end if
  end subroutine fill_halo_zero_r8

  !> fill the halo region of a 32-bit real array with zeros
  subroutine fill_halo_zero_r4(data, whalo, ehalo, shalo, nhalo, xshift, yshift, isc, iec, jsc, jec, isd, ied, jsd, jed)
    real(kind=r4_kind), dimension(isd:,jsd:,:), intent(inout) :: data
    integer,                         intent(in) :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer,                         intent(in) :: whalo, ehalo, shalo, nhalo, xshift, yshift

    if(whalo >=0) then
       data(iec+ehalo+1+xshift:ied+xshift,jsd:jed+yshift,:) = 0
       data(isd:isc-whalo-1,jsd:jed+yshift,:) = 0
    else
       data(iec+1+xshift:iec-ehalo+xshift,jsc+shalo:jec-nhalo+yshift,:) = 0
       data(isc+whalo:isc-1,jsc+shalo:jec-nhalo+yshift,:) = 0
    end if

    if(shalo>=0) then
       data(isd:ied+xshift, jec+nhalo+1+yshift:jed+yshift,:) = 0
       data(isd:ied+xshift, jsd:jsc-shalo-1,:) = 0
    else
       data(isc+whalo:iec-ehalo+xshift,jec+1+yshift:jec-nhalo+yshift,:) = 0
       data(isc+whalo:iec-ehalo+xshift,jsc+shalo:jsc-1,:) = 0
    end if
  end subroutine fill_halo_zero_r4

!> fill the halo region of a 64-bit integer array with zeros
  subroutine fill_halo_zero_i8(data, whalo, ehalo, shalo, nhalo, xshift, yshift, isc, iec, jsc, jec, isd, ied, jsd, jed)
    integer(kind=i8_kind), dimension(isd:,jsd:,:), intent(inout) :: data
    integer,                         intent(in) :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer,                         intent(in) :: whalo, ehalo, shalo, nhalo, xshift, yshift

    if(whalo >=0) then
       data(iec+ehalo+1+xshift:ied+xshift,jsd:jed+yshift,:) = 0
       data(isd:isc-whalo-1,jsd:jed+yshift,:) = 0
    else
       data(iec+1+xshift:iec-ehalo+xshift,jsc+shalo:jec-nhalo+yshift,:) = 0
       data(isc+whalo:isc-1,jsc+shalo:jec-nhalo+yshift,:) = 0
    end if

    if(shalo>=0) then
       data(isd:ied+xshift, jec+nhalo+1+yshift:jed+yshift,:) = 0
       data(isd:ied+xshift, jsd:jsc-shalo-1,:) = 0
    else
       data(isc+whalo:iec-ehalo+xshift,jec+1+yshift:jec-nhalo+yshift,:) = 0
       data(isc+whalo:iec-ehalo+xshift,jsc+shalo:jsc-1,:) = 0
    end if
  end subroutine fill_halo_zero_i8

!> fill the halo region of a 32-bit integer array with zeros
  subroutine fill_halo_zero_i4(data, whalo, ehalo, shalo, nhalo, xshift, yshift, isc, iec, jsc, jec, isd, ied, jsd, jed)
    integer(kind=i4_kind), dimension(isd:,jsd:,:), intent(inout) :: data
    integer,                         intent(in) :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer,                         intent(in) :: whalo, ehalo, shalo, nhalo, xshift, yshift

    if(whalo >=0) then
       data(iec+ehalo+1+xshift:ied+xshift,jsd:jed+yshift,:) = 0
       data(isd:isc-whalo-1,jsd:jed+yshift,:) = 0
    else
       data(iec+1+xshift:iec-ehalo+xshift,jsc+shalo:jec-nhalo+yshift,:) = 0
       data(isc+whalo:isc-1,jsc+shalo:jec-nhalo+yshift,:) = 0
    end if

    if(shalo>=0) then
       data(isd:ied+xshift, jec+nhalo+1+yshift:jed+yshift,:) = 0
       data(isd:ied+xshift, jsd:jsc-shalo-1,:) = 0
    else
       data(isc+whalo:iec-ehalo+xshift,jec+1+yshift:jec-nhalo+yshift,:) = 0
       data(isc+whalo:iec-ehalo+xshift,jsc+shalo:jsc-1,:) = 0
    end if
  end subroutine fill_halo_zero_i4


  !> fill the halo region of 64-bit array on a regular grid
  subroutine fill_regular_refinement_halo_r8( data, data_all, ni, nj, tm, te, tse, ts, &
                                             tsw, tw, tnw, tn, tne, ioff, joff )
    real(kind=r8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real(kind=r8_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, dimension(:),                intent(in)    :: ni, nj
    integer,                              intent(in)    :: tm, te, tse, ts, tsw, tw, tnw, tn, tne
    integer,                              intent(in)    :: ioff, joff


    if(te>0) data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1:nj(tm)+joff,                   :) = &
             data_all(1+ioff:ehalo+ioff,               1:nj(te)+joff,                   :,te)  ! east
    if(ts>0) data    (1:ni(tm)+ioff,                   1-shalo:0,                       :) = &
             data_all(1:ni(ts)+ioff,                   nj(ts)-shalo+1:nj(ts),           :,ts)  ! south
    if(tw>0) data    (1-whalo:0,                       1:nj(tm)+joff,                   :) = &
             data_all(ni(tw)-whalo+1:ni(tw),           1:nj(tw)+joff,                   :,tw)  ! west
    if(tn>0) data    (1:ni(tm)+ioff,                   nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1:ni(tn)+ioff,                   1+joff:nhalo+joff,               :,tn)  ! north
    if(tse>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1-shalo:0,                       :) = &
             data_all(1+ioff:ehalo+ioff,               nj(tse)-shalo+1:nj(tse),         :,tse) ! southeast
    if(tsw>0)data    (1-whalo:0,                       1-shalo:0,                       :) = &
             data_all(ni(tsw)-whalo+1:ni(tsw),         nj(tsw)-shalo+1:nj(tsw),         :,tsw) ! southwest
    if(tne>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1+ioff:ehalo+ioff,               1+joff:nhalo+joff,               :,tnw) ! northeast
    if(tnw>0)data    (1-whalo:0,                       nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(ni(tnw)-whalo+1:ni(tnw),         1+joff:nhalo+joff,               :,tne) ! northwest

  end subroutine fill_regular_refinement_halo_r8

  !> fill the halo region of 32-bit array on a regular grid
  subroutine fill_regular_refinement_halo_r4( data, data_all, ni, nj, tm, te, tse, ts, tsw, tw, tnw, tn, tne, ioff, joff )
    real(kind=r4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real(kind=r4_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, dimension(:),                intent(in)    :: ni, nj
    integer,                              intent(in)    :: tm, te, tse, ts, tsw, tw, tnw, tn, tne
    integer,                              intent(in)    :: ioff, joff


    if(te>0) data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1:nj(tm)+joff,                   :) = &
             data_all(1+ioff:ehalo+ioff,               1:nj(te)+joff,                   :,te)  ! east
    if(ts>0) data    (1:ni(tm)+ioff,                   1-shalo:0,                       :) = &
             data_all(1:ni(ts)+ioff,                   nj(ts)-shalo+1:nj(ts),           :,ts)  ! south
    if(tw>0) data    (1-whalo:0,                       1:nj(tm)+joff,                   :) = &
             data_all(ni(tw)-whalo+1:ni(tw),           1:nj(tw)+joff,                   :,tw)  ! west
    if(tn>0) data    (1:ni(tm)+ioff,                   nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1:ni(tn)+ioff,                   1+joff:nhalo+joff,               :,tn)  ! north
    if(tse>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1-shalo:0,                       :) = &
             data_all(1+ioff:ehalo+ioff,               nj(tse)-shalo+1:nj(tse),         :,tse) ! southeast
    if(tsw>0)data    (1-whalo:0,                       1-shalo:0,                       :) = &
             data_all(ni(tsw)-whalo+1:ni(tsw),         nj(tsw)-shalo+1:nj(tsw),         :,tsw) ! southwest
    if(tne>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1+ioff:ehalo+ioff,               1+joff:nhalo+joff,               :,tnw) ! northeast
    if(tnw>0)data    (1-whalo:0,                       nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(ni(tnw)-whalo+1:ni(tnw),         1+joff:nhalo+joff,               :,tne) ! northwest

  end subroutine fill_regular_refinement_halo_r4

!> fill the halo region of 64-bit integer array on a regular grid
  subroutine fill_regular_refinement_halo_i8( data, data_all, ni, nj, tm, te, tse, ts, tsw, &
                                             tw, tnw, tn, tne, ioff, joff )
    integer(kind=i8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer(kind=i8_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, dimension(:),                intent(in)    :: ni, nj
    integer,                              intent(in)    :: tm, te, tse, ts, tsw, tw, tnw, tn, tne
    integer,                              intent(in)    :: ioff, joff


    if(te>0) data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1:nj(tm)+joff,                   :) = &
             data_all(1+ioff:ehalo+ioff,               1:nj(te)+joff,                   :,te)  ! east
    if(ts>0) data    (1:ni(tm)+ioff,                   1-shalo:0,                       :) = &
             data_all(1:ni(ts)+ioff,                   nj(ts)-shalo+1:nj(ts),           :,ts)  ! south
    if(tw>0) data    (1-whalo:0,                       1:nj(tm)+joff,                   :) = &
             data_all(ni(tw)-whalo+1:ni(tw),           1:nj(tw)+joff,                   :,tw)  ! west
    if(tn>0) data    (1:ni(tm)+ioff,                   nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1:ni(tn)+ioff,                   1+joff:nhalo+joff,               :,tn)  ! north
    if(tse>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1-shalo:0,                       :) = &
             data_all(1+ioff:ehalo+ioff,               nj(tse)-shalo+1:nj(tse),         :,tse) ! southeast
    if(tsw>0)data    (1-whalo:0,                       1-shalo:0,                       :) = &
             data_all(ni(tsw)-whalo+1:ni(tsw),         nj(tsw)-shalo+1:nj(tsw),         :,tsw) ! southwest
    if(tne>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1+ioff:ehalo+ioff,               1+joff:nhalo+joff,               :,tnw) ! northeast
    if(tnw>0)data    (1-whalo:0,                       nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(ni(tnw)-whalo+1:ni(tnw),         1+joff:nhalo+joff,               :,tne) ! northwest

  end subroutine fill_regular_refinement_halo_i8

!> fill the halo region of 32-bit integer array on a regular grid
  subroutine fill_regular_refinement_halo_i4( data, data_all, ni, nj, tm, te, tse, ts, tsw, tw, tnw, tn, tne, ioff, joff )
    integer(kind=i4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer(kind=i4_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, dimension(:),                intent(in)    :: ni, nj
    integer,                              intent(in)    :: tm, te, tse, ts, tsw, tw, tnw, tn, tne
    integer,                              intent(in)    :: ioff, joff


    if(te>0) data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1:nj(tm)+joff,                   :) = &
             data_all(1+ioff:ehalo+ioff,               1:nj(te)+joff,                   :,te)  ! east
    if(ts>0) data    (1:ni(tm)+ioff,                   1-shalo:0,                       :) = &
             data_all(1:ni(ts)+ioff,                   nj(ts)-shalo+1:nj(ts),           :,ts)  ! south
    if(tw>0) data    (1-whalo:0,                       1:nj(tm)+joff,                   :) = &
             data_all(ni(tw)-whalo+1:ni(tw),           1:nj(tw)+joff,                   :,tw)  ! west
    if(tn>0) data    (1:ni(tm)+ioff,                   nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1:ni(tn)+ioff,                   1+joff:nhalo+joff,               :,tn)  ! north
    if(tse>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, 1-shalo:0,                       :) = &
             data_all(1+ioff:ehalo+ioff,               nj(tse)-shalo+1:nj(tse),         :,tse) ! southeast
    if(tsw>0)data    (1-whalo:0,                       1-shalo:0,                       :) = &
             data_all(ni(tsw)-whalo+1:ni(tsw),         nj(tsw)-shalo+1:nj(tsw),         :,tsw) ! southwest
    if(tne>0)data    (ni(tm)+1+ioff:ni(tm)+ehalo+ioff, nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(1+ioff:ehalo+ioff,               1+joff:nhalo+joff,               :,tnw) ! northeast
    if(tnw>0)data    (1-whalo:0,                       nj(tm)+1+joff:nj(tm)+nhalo+joff, :) = &
             data_all(ni(tnw)-whalo+1:ni(tnw),         1+joff:nhalo+joff,               :,tne) ! northwest

  end subroutine fill_regular_refinement_halo_i4

  ! Fill the halo points of a 64-bit real array on the regular mosaic grid
  subroutine fill_regular_mosaic_halo_r8(data, data_all, te, tse, ts, tsw, tw, tnw, tn, tne)
    real(kind=r8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real(kind=r8_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, intent(in)    :: te, tse, ts, tsw, tw, tnw, tn, tne

    data(nx+1:nx+ehalo, 1:ny,          :) = data_all(1:ehalo,       1:ny,          :, te) ! east
    data(1:nx,          1-shalo:0,     :) = data_all(1:nx,          ny-shalo+1:ny, :, ts) ! south
    data(1-whalo:0,     1:ny,          :) = data_all(nx-whalo+1:nx, 1:ny,          :, tw) ! west
    data(1:nx,          ny+1:ny+nhalo, :) = data_all(1:nx,          1:nhalo,       :, tn) ! north
    data(nx+1:nx+ehalo, 1-shalo:0,     :) = data_all(1:ehalo,       ny-shalo+1:ny, :,tse) ! southeast
    data(1-whalo:0,     1-shalo:0,     :) = data_all(nx-whalo+1:nx, ny-shalo+1:ny, :,tsw) ! southwest
    data(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = data_all(1:ehalo,       1:nhalo,       :,tnw) ! northeast
    data(1-whalo:0,     ny+1:ny+nhalo, :) = data_all(nx-whalo+1:nx, 1:nhalo,       :,tne) ! northwest
  end subroutine fill_regular_mosaic_halo_r8

  !> Fill the halo points of a 32-bit real array on the regular mosaic grid
  subroutine fill_regular_mosaic_halo_r4(data, data_all, te, tse, ts, tsw, tw, tnw, tn, tne)
    real(kind=r4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real(kind=r4_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, intent(in)    :: te, tse, ts, tsw, tw, tnw, tn, tne

    data(nx+1:nx+ehalo, 1:ny,          :) = data_all(1:ehalo,       1:ny,          :, te) ! east
    data(1:nx,          1-shalo:0,     :) = data_all(1:nx,          ny-shalo+1:ny, :, ts) ! south
    data(1-whalo:0,     1:ny,          :) = data_all(nx-whalo+1:nx, 1:ny,          :, tw) ! west
    data(1:nx,          ny+1:ny+nhalo, :) = data_all(1:nx,          1:nhalo,       :, tn) ! north
    data(nx+1:nx+ehalo, 1-shalo:0,     :) = data_all(1:ehalo,       ny-shalo+1:ny, :,tse) ! southeast
    data(1-whalo:0,     1-shalo:0,     :) = data_all(nx-whalo+1:nx, ny-shalo+1:ny, :,tsw) ! southwest
    data(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = data_all(1:ehalo,       1:nhalo,       :,tnw) ! northeast
    data(1-whalo:0,     ny+1:ny+nhalo, :) = data_all(nx-whalo+1:nx, 1:nhalo,       :,tne) ! northwest
  end subroutine fill_regular_mosaic_halo_r4

  ! Fill the halo points of a 64-bit integer array on the regular mosaic grid
  subroutine fill_regular_mosaic_halo_i8(data, data_all, te, tse, ts, tsw, tw, tnw, tn, tne)
    integer(kind=i8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer(kind=i8_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, intent(in)    :: te, tse, ts, tsw, tw, tnw, tn, tne

    data(nx+1:nx+ehalo, 1:ny,          :) = data_all(1:ehalo,       1:ny,          :, te) ! east
    data(1:nx,          1-shalo:0,     :) = data_all(1:nx,          ny-shalo+1:ny, :, ts) ! south
    data(1-whalo:0,     1:ny,          :) = data_all(nx-whalo+1:nx, 1:ny,          :, tw) ! west
    data(1:nx,          ny+1:ny+nhalo, :) = data_all(1:nx,          1:nhalo,       :, tn) ! north
    data(nx+1:nx+ehalo, 1-shalo:0,     :) = data_all(1:ehalo,       ny-shalo+1:ny, :,tse) ! southeast
    data(1-whalo:0,     1-shalo:0,     :) = data_all(nx-whalo+1:nx, ny-shalo+1:ny, :,tsw) ! southwest
    data(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = data_all(1:ehalo,       1:nhalo,       :,tnw) ! northeast
    data(1-whalo:0,     ny+1:ny+nhalo, :) = data_all(nx-whalo+1:nx, 1:nhalo,       :,tne) ! northwest
  end subroutine fill_regular_mosaic_halo_i8

  !> Fill the halo points of a 64-bit integer array on the regular mosaic grid
  subroutine fill_regular_mosaic_halo_i4(data, data_all, te, tse, ts, tsw, tw, tnw, tn, tne)
    integer(kind=i4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer(kind=i4_kind), dimension(:,:,:,:),             intent(in)    :: data_all
    integer, intent(in)    :: te, tse, ts, tsw, tw, tnw, tn, tne

    data(nx+1:nx+ehalo, 1:ny,          :) = data_all(1:ehalo,       1:ny,          :, te) ! east
    data(1:nx,          1-shalo:0,     :) = data_all(1:nx,          ny-shalo+1:ny, :, ts) ! south
    data(1-whalo:0,     1:ny,          :) = data_all(nx-whalo+1:nx, 1:ny,          :, tw) ! west
    data(1:nx,          ny+1:ny+nhalo, :) = data_all(1:nx,          1:nhalo,       :, tn) ! north
    data(nx+1:nx+ehalo, 1-shalo:0,     :) = data_all(1:ehalo,       ny-shalo+1:ny, :,tse) ! southeast
    data(1-whalo:0,     1-shalo:0,     :) = data_all(nx-whalo+1:nx, ny-shalo+1:ny, :,tsw) ! southwest
    data(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = data_all(1:ehalo,       1:nhalo,       :,tnw) ! northeast
    data(1-whalo:0,     ny+1:ny+nhalo, :) = data_all(nx-whalo+1:nx, 1:nhalo,       :,tne) ! northwest
  end subroutine fill_regular_mosaic_halo_i4

  !> Fill the halo region of a 64-bit array real on a domain with a folded north edge
  subroutine fill_folded_north_halo_r8(data, ioff, joff, ishift, jshift, sign)
    real(kind=r8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:ny+jshift,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift,1:ny+jshift,:) ! east
    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,nyp+1:nyp+nhalo,:) = sign*data(whalo+m2:1+ishift:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(m1+1:nx+m2,nyp+1:nyp+nhalo,:) = sign*data(nx+ishift:1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(nx+m2+1:nxp+ehalo,nyp+1:nyp+nhalo,:) = sign*data(nx:nx-ehalo+m1+1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)

  end subroutine fill_folded_north_halo_r8

  !> Fill the halo region of a 32-bit real array on a domain with a folded north edge
  subroutine fill_folded_north_halo_r4(data, ioff, joff, ishift, jshift, sign)
    real(kind=r4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:ny+jshift,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift,1:ny+jshift,:) ! east

    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,nyp+1:nyp+nhalo,:) = sign*data(whalo+m2:1+ishift:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(m1+1:nx+m2,nyp+1:nyp+nhalo,:) = sign*data(nx+ishift:1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(nx+m2+1:nxp+ehalo,nyp+1:nyp+nhalo,:) = sign*data(nx:nx-ehalo+m1+1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)

  end subroutine fill_folded_north_halo_r4

  !> Fill the halo region of a 64-bit integer array on a domain with a folded north edge
  subroutine fill_folded_north_halo_i8(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:ny+jshift,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift,1:ny+jshift,:) ! east
    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,nyp+1:nyp+nhalo,:) = sign*data(whalo+m2:1+ishift:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(m1+1:nx+m2,nyp+1:nyp+nhalo,:) = sign*data(nx+ishift:1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(nx+m2+1:nxp+ehalo,nyp+1:nyp+nhalo,:) = sign*data(nx:nx-ehalo+m1+1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)

  end subroutine fill_folded_north_halo_i8

  !> Fill the halo region of a 32-bit integer array on a domain with a folded north edge
  subroutine fill_folded_north_halo_i4(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:ny+jshift,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift,1:ny+jshift,:) ! east

    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,nyp+1:nyp+nhalo,:) = sign*data(whalo+m2:1+ishift:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(m1+1:nx+m2,nyp+1:nyp+nhalo,:) = sign*data(nx+ishift:1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(nx+m2+1:nxp+ehalo,nyp+1:nyp+nhalo,:) = sign*data(nx:nx-ehalo+m1+1:-1,nyp-joff:nyp-nhalo-joff+1:-1,:)

  end subroutine fill_folded_north_halo_i4

  !> Fill the halo region of a 64-bit real array on a domain with a folded south edge
  subroutine fill_folded_south_halo_r8(data, ioff, joff, ishift, jshift, sign)
    real(kind=r8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:nyp,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift, 1:nyp,:) ! east
    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,1-shalo:0,:) = sign*data(whalo+m2:1+ishift:-1, shalo+jshift:1+jshift:-1,:)

    data(m1+1:nx+m2,1-shalo:0,:) = sign*data(nxp:1:-1,shalo+jshift:1+jshift:-1,:)
    data(nx+m2+1:nxp+ehalo,1-shalo:0,:) = sign*data(nx:nx-ehalo+m1+1:-1,shalo+jshift:1+jshift:-1,:)

  end subroutine fill_folded_south_halo_r8

  !> Fill the halo region of a 32-bit real array on a domain with a folded south edge
  subroutine fill_folded_south_halo_r4(data, ioff, joff, ishift, jshift, sign)
    real(kind=r4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:nyp,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift, 1:nyp,:) ! east
    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,1-shalo:0,:) = sign*data(whalo+m2:1+ishift:-1, shalo+jshift:1+jshift:-1,:)

    data(m1+1:nx+m2,1-shalo:0,:) = sign*data(nxp:1:-1,shalo+jshift:1+jshift:-1,:)
    data(nx+m2+1:nxp+ehalo,1-shalo:0,:) = sign*data(nx:nx-ehalo+m1+1:-1,shalo+jshift:1+jshift:-1,:)

  end subroutine fill_folded_south_halo_r4

  !> Fill the halo region of a 64-bit intger array on a domain with a folded south edge
  subroutine fill_folded_south_halo_i8(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:nyp,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift, 1:nyp,:) ! east
    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,1-shalo:0,:) = sign*data(whalo+m2:1+ishift:-1, shalo+jshift:1+jshift:-1,:)

    data(m1+1:nx+m2,1-shalo:0,:) = sign*data(nxp:1:-1,shalo+jshift:1+jshift:-1,:)
    data(nx+m2+1:nxp+ehalo,1-shalo:0,:) = sign*data(nx:nx-ehalo+m1+1:-1,shalo+jshift:1+jshift:-1,:)

  end subroutine fill_folded_south_halo_i8

  !> Fill the halo region of a 32-bit integer array on a domain with a folded south edge
  subroutine fill_folded_south_halo_i4(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,1:nyp,:) = data(nx-whalo+1:nx,1:nyp,:) ! west
    data(nx+1:nx+ehalo+ishift,1:nyp,:) = data(1:ehalo+ishift, 1:nyp,:) ! east
    if(m1 .GE. 1-whalo) &
      data(1-whalo:m1,1-shalo:0,:) = sign*data(whalo+m2:1+ishift:-1, shalo+jshift:1+jshift:-1,:)

    data(m1+1:nx+m2,1-shalo:0,:) = sign*data(nxp:1:-1,shalo+jshift:1+jshift:-1,:)
    data(nx+m2+1:nxp+ehalo,1-shalo:0,:) = sign*data(nx:nx-ehalo+m1+1:-1,shalo+jshift:1+jshift:-1,:)

  end subroutine fill_folded_south_halo_i4

  !> Fill the halo region of a 64-bit real array on a domain with a folded west edge
  subroutine fill_folded_west_halo_r8(data, ioff, joff, ishift, jshift, sign)
    real(kind=r8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0,:) = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo,:) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(1-whalo:0, 1-shalo:m1,:) = sign*data(whalo+ishift:1+ishift:-1,shalo+m2:1+jshift:-1,:)
    data(1-whalo:0, m1+1:ny+m2,:) = sign*data(whalo+ishift:1+ishift:-1, nyp:1:-1, :)
    data(1-whalo:0, ny+m2+1:nyp+nhalo,:) = sign*data(whalo+ishift:1+ishift:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_west_halo_r8

  !> Fill the halo region of a 32-bit real array on a domain with a folded west edge
  subroutine fill_folded_west_halo_r4(data, ioff, joff, ishift, jshift, sign)
    real(kind=r4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0,:) = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo,:) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(1-whalo:0, 1-shalo:m1,:) = sign*data(whalo+ishift:1+ishift:-1,shalo+m2:1+jshift:-1,:)
    data(1-whalo:0, m1+1:ny+m2,:) = sign*data(whalo+ishift:1+ishift:-1, nyp:1:-1, :)
    data(1-whalo:0, ny+m2+1:nyp+nhalo,:) = sign*data(whalo+ishift:1+ishift:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_west_halo_r4

  !> Fill the halo region of a 64-bit integer array on a domain with a folded west edge
  subroutine fill_folded_west_halo_i8(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0,:) = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo,:) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(1-whalo:0, 1-shalo:m1,:) = sign*data(whalo+ishift:1+ishift:-1,shalo+m2:1+jshift:-1,:)
    data(1-whalo:0, m1+1:ny+m2,:) = sign*data(whalo+ishift:1+ishift:-1, nyp:1:-1, :)
    data(1-whalo:0, ny+m2+1:nyp+nhalo,:) = sign*data(whalo+ishift:1+ishift:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_west_halo_i8

  !> Fill the halo region of a 32-bit integer array on a domain with a folded west edge
  subroutine fill_folded_west_halo_i4(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0,:) = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo,:) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(1-whalo:0, 1-shalo:m1,:) = sign*data(whalo+ishift:1+ishift:-1,shalo+m2:1+jshift:-1,:)
    data(1-whalo:0, m1+1:ny+m2,:) = sign*data(whalo+ishift:1+ishift:-1, nyp:1:-1, :)
    data(1-whalo:0, ny+m2+1:nyp+nhalo,:) = sign*data(whalo+ishift:1+ishift:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_west_halo_i4

  !> Fill the halo region of a 64-bit real array on a domain with a folded east edge
  subroutine fill_folded_east_halo_r8(data, ioff, joff, ishift, jshift, sign)
    real(kind=r8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)  = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(nxp+1:nxp+ehalo, 1-shalo:m1, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, shalo+m2:1+jshift:-1,:)

    data(nxp+1:nxp+ehalo, m1+1:ny+m2, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, nyp:1:-1, :)
    data(nxp+1:nxp+ehalo, ny+m2+1:nyp+nhalo,:) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_east_halo_r8

  !> Fill the halo region of a 32-bit real array on a domain with a folded east edge
  subroutine fill_folded_east_halo_r4(data, ioff, joff, ishift, jshift, sign)
    real(kind=r4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)  = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(nxp+1:nxp+ehalo, 1-shalo:m1, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, shalo+m2:1+jshift:-1,:)

    data(nxp+1:nxp+ehalo, m1+1:ny+m2, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, nyp:1:-1, :)
    data(nxp+1:nxp+ehalo, ny+m2+1:nyp+nhalo,:) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_east_halo_r4

  !> Fill the halo region of a 64-bit integer array on a domain with a folded east edge
  subroutine fill_folded_east_halo_i8(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i8_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)  = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(nxp+1:nxp+ehalo, 1-shalo:m1, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, shalo+m2:1+jshift:-1,:)

    data(nxp+1:nxp+ehalo, m1+1:ny+m2, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, nyp:1:-1, :)
    data(nxp+1:nxp+ehalo, ny+m2+1:nyp+nhalo,:) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_east_halo_i8

  !> Fill the halo region of a 32-bit integer array on a domain with a folded east edge
  subroutine fill_folded_east_halo_i4(data, ioff, joff, ishift, jshift, sign)
    integer(kind=i4_kind), dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer, intent(in) :: ioff, joff, ishift, jshift, sign
    ! local
    integer :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)  = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) &
      data(nxp+1:nxp+ehalo, 1-shalo:m1, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, shalo+m2:1+jshift:-1,:)

    data(nxp+1:nxp+ehalo, m1+1:ny+m2, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, nyp:1:-1, :)
    data(nxp+1:nxp+ehalo, ny+m2+1:nyp+nhalo,:) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_east_halo_i4

end module fill_halo

