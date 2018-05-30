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

!> \file
!! \brief Contains the \ref block_control_mod module

module block_control_mod
#include <fms_platform.h>

use mpp_mod,         only: mpp_error, NOTE, WARNING, FATAL
use mpp_domains_mod, only: mpp_compute_extent

 public block_control_type

 type ix_type
   integer, dimension(:,:), _ALLOCATABLE :: ix _NULL
 end type ix_type

 type pk_type
   integer, dimension(:), _ALLOCATABLE :: ii _NULL
   integer, dimension(:), _ALLOCATABLE :: jj _NULL
 end type pk_type

 type block_control_type
   integer :: nx_block, ny_block  !< blocking factor using mpp-style decomposition
   integer :: nblks               !< number of blocks cover MPI domain
   integer :: isc, iec, jsc, jec  !< MPI domain global extents
   integer :: npz                 !< vertical extent
   integer, dimension(:),        _ALLOCATABLE :: ibs   _NULL, &  !< block extents for mpp-style
                                                 ibe   _NULL, &  !! decompositions
                                                 jbs   _NULL, &
                                                 jbe   _NULL
   type(ix_type), dimension(:),  _ALLOCATABLE :: ix    _NULL !< dereference packed index from global index
   !--- packed blocking fields
   integer, dimension(:),        _ALLOCATABLE :: blksz _NULL !< number of points in each individual block
                                                             !! blocks are not required to be uniforom in size
   integer, dimension(:,:),      _ALLOCATABLE :: blkno _NULL !< dereference block number using global indices
   integer, dimension(:,:),      _ALLOCATABLE :: ixp   _NULL !< dereference packed index from global indices
                                                             !! must be used in conjuction with blkno
   type(pk_type), dimension(:),  _ALLOCATABLE :: index _NULL !< dereference global indices from
                                                             !! block/ixp combo
 end type block_control_type

public :: define_blocks, define_blocks_packed

contains

!###############################################################################
!> \fn define_blocks
!!
!! \brief Sets up "blocks" used for OpenMP threading of column-based
!!        calculations using rad_n[x/y]xblock from coupler_nml
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),         intent(in)    :: component
!! type(block_control_type), intent(inout) :: Block
!! integer,                  intent(in)    :: isc, iec, jsc, jec, kpts
!! integer,                  intent(in)    :: nx_block, ny_block
!! logical,                  intent(inout) :: message
!! \endcode
!!
!! \param [in]    <component>
!! \param [inout] <Block>
!! \param [in]    <isc>
!! \param [in]    <iec>
!! \param [in]    <jsc>
!! \param [in]    <jec>
!! \param [in]    <kpts>
!! \param [in]    <nx_block>
!! \param [in]    <ny_block>
!! \param [inout] <message>
!!
  subroutine define_blocks (component, Block, isc, iec, jsc, jec, kpts, &
                            nx_block, ny_block, message)
    character(len=*),         intent(in)    :: component
    type(block_control_type), intent(inout) :: Block
    integer,                  intent(in)    :: isc, iec, jsc, jec, kpts
    integer,                  intent(in)    :: nx_block, ny_block
    logical,                  intent(inout) :: message

!-------------------------------------------------------------------------------
! Local variables:
!       blocks
!       i1
!       i2
!       j1
!       j2
!       text
!       i
!       j
!       nblks
!       ii
!       jj
!-------------------------------------------------------------------------------

    integer :: blocks
    integer, dimension(nx_block) :: i1, i2
    integer, dimension(ny_block) :: j1, j2
    character(len=256) :: text
    integer :: i, j, nblks, ii, jj

    if (message) then
      if ((mod(iec-isc+1,nx_block) .ne. 0) .or. (mod(jec-jsc+1,ny_block) .ne. 0)) then
        write( text,'(a,a,2i4,a,2i4,a)' ) trim(component),'define_blocks: domain (',&
             (iec-isc+1), (jec-jsc+1),') is not an even divisor with definition (',&
             nx_block, ny_block,') - blocks will not be uniform'
        call mpp_error (WARNING, trim(text))
      endif
      message = .false.
    endif

!--- set up blocks
    if (iec-isc+1 .lt. nx_block) &
        call mpp_error(FATAL, 'block_control: number of '//trim(component)//' nxblocks .gt. &
                             &number of elements in MPI-domain size')
    if (jec-jsc+1 .lt. ny_block) &
        call mpp_error(FATAL, 'block_control: number of '//trim(component)//' nyblocks .gt. &
                             &number of elements in MPI-domain size')
    call mpp_compute_extent(isc,iec,nx_block,i1,i2)
    call mpp_compute_extent(jsc,jec,ny_block,j1,j2)

    nblks = nx_block*ny_block
    Block%isc = isc
    Block%iec = iec
    Block%jsc = jsc
    Block%jec = jec
    Block%npz = kpts
    Block%nx_block = nx_block
    Block%ny_block = ny_block
    Block%nblks = nblks

    if (.not._ALLOCATED(Block%ibs)) &
         allocate (Block%ibs(nblks), &
                   Block%ibe(nblks), &
                   Block%jbs(nblks), &
                   Block%jbe(nblks), &
                   Block%ix(nblks) )

    blocks=0
    do j = 1, ny_block
      do i = 1, nx_block
        blocks = blocks + 1
        Block%ibs(blocks) = i1(i)
        Block%jbs(blocks) = j1(j)
        Block%ibe(blocks) = i2(i)
        Block%jbe(blocks) = j2(j)
        allocate(Block%ix(blocks)%ix(i1(i):i2(i),j1(j):j2(j)) )
        ix = 0
        do jj = j1(j), j2(j)
          do ii = i1(i), i2(i)
            ix = ix+1
            Block%ix(blocks)%ix(ii,jj) = ix
          enddo
        enddo
      enddo
    enddo

  end subroutine define_blocks



!###############################################################################
!> \fn define_blocks_packed
!!
!! \brief Creates and populates a data type which is used for defining the
!!        sub-blocks of the MPI-domain to enhance OpenMP and memory performance.
!!        Uses a packed concept
!!
!! <b> Parameters: </b>
!!
!! \code{.f90}
!! character(len=*),         intent(in)    :: component
!! type(block_control_type), intent(inout) :: Block
!! integer,                  intent(in)    :: isc, iec, jsc, jec, kpts
!! integer,                  intent(inout) :: blksz
!! logical,                  intent(inout) :: message
!! \endcode
!!
!! \param [in]    <component>
!! \param [inout] <Block>
!! \param [in]    <isc>
!! \param [in]    <iec>
!! \param [in]    <jsc>
!! \param [in]    <jec>
!! \param [in]    <kpts>
!! \param [inout] <blksz>
!! \param [inout] <message>
!!
  subroutine define_blocks_packed (component, Block, isc, iec, jsc, jec, &
                                   kpts, blksz, message)
    character(len=*),         intent(in)    :: component
    type(block_control_type), intent(inout) :: Block
    integer,                  intent(in)    :: isc, iec, jsc, jec, kpts
    integer,                  intent(inout) :: blksz
    logical,                  intent(inout) :: message

!-------------------------------------------------------------------------------
! Local variables:
!       nblks
!       lblksz
!       tot_pts
!       ii
!       jj
!       nb
!       ix
!       text
!-------------------------------------------------------------------------------

    integer :: nblks, lblksz, tot_pts, ii, jj,  nb, ix
    character(len=256) :: text

    tot_pts = (iec - isc + 1) * (jec - jsc + 1)
    if (blksz < 0) then
      nblks = 1
      blksz = tot_pts
    else
      if (mod(tot_pts,blksz) .eq. 0) then
        nblks = tot_pts/blksz
      else
        nblks = ceiling(real(tot_pts)/real(blksz))
      endif
    endif

    if (message) then
      if (mod(tot_pts,blksz) .ne. 0) then
        write( text,'(a,a,2i4,a,i4,a,i4)' ) trim(component),'define_blocks_packed: domain (',&
             (iec-isc+1), (jec-jsc+1),') is not an even divisor with definition (',&
             blksz,') - blocks will not be uniform with a remainder of ',mod(tot_pts,blksz)
        call mpp_error (WARNING, trim(text))
      endif
      message = .false.
    endif

    Block%isc   = isc
    Block%iec   = iec
    Block%jsc   = jsc
    Block%jec   = jec
    Block%npz   = kpts
    Block%nblks = nblks
    if (.not. _ALLOCATED(Block%blksz)) &
      allocate (Block%blksz(nblks), &
                Block%index(nblks), &
                Block%blkno(isc:iec,jsc:jec), &
                Block%ixp(isc:iec,jsc:jec))

!--- set up blocks
    do nb = 1, nblks
      lblksz = blksz
      if (nb .EQ. nblks) lblksz = tot_pts - (nb-1) * blksz
      Block%blksz(nb) = lblksz
      allocate (Block%index(nb)%ii(lblksz), &
                Block%index(nb)%jj(lblksz))
    enddo

!--- set up packed indices
    nb = 1
    ix = 0
    do jj = jsc, jec
      do ii = isc, iec
        ix = ix + 1
        if (ix .GT. blksz) then
          ix = 1
          nb = nb + 1
        endif
        Block%ixp(ii,jj) = ix
        Block%blkno(ii,jj) = nb
        Block%index(nb)%ii(ix) = ii
        Block%index(nb)%jj(ix) = jj
      enddo
    enddo

  end subroutine define_blocks_packed

end module block_control_mod
