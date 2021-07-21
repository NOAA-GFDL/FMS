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
!> @file
!> @brief Test mpp_update_domains on arrays of real numbers using different layouts and data precision
!!
!> @author Jessica Liptak
!!
!> @note This test is an extension of the routine test_halo_upate in test_mpp_domains.

module test_mpp_update_domains_real

  use compare_data_checksums, only : compare_checksums
  use fill_halo
  use mpp_mod, only : FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod, only : mpp_pe, mpp_npes, mpp_root_pe, mpp_error
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_broadcast, mpp_transmit, mpp_sum, mpp_max, mpp_chksum, ALL_PES
  use mpp_mod, only : mpp_gather, mpp_sync_self
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE, DGRID_NE, AGRID
  use mpp_domains_mod, only : FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE, FOLD_WEST_EDGE, FOLD_EAST_EDGE
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D, BITWISE_EFP_SUM
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod, only : mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST, mpp_define_mosaic_pelist
  use mpp_domains_mod, only : mpp_get_global_domain, ZERO, NINETY, MINUS_NINETY
  use mpp_domains_mod, only : mpp_deallocate_domain
  use platform_mod

  implicit none
  private
  integer :: id
  integer :: nx=64, ny=64, nz=10, stackmax=10000000
  integer :: i, j, k, n
  integer :: layout(2)
  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  integer :: x_cyclic_offset = 3   ! to be used in test_cyclic_offset
  integer :: y_cyclic_offset = -4  ! to be used in test_cyclic_offset
  character(len=32) :: warn_level = "fatal"
  integer :: wide_halo_x = 0, wide_halo_y = 0
  integer :: nx_cubic = 0, ny_cubic = 0
  integer :: ensemble_size = 1
  integer :: layout_cubic(2) = (/0,0/)
  integer :: layout_tripolar(2) = (/0,0/)
  integer :: layout_ensemble(2) = (/0,0/)
  public :: test_halo_update_r8, test_halo_update_r4, test_subset_update_r8, test_subset_update_r4

  contains

  !> Perform simple addition on 64-bit real arrays in different domain configurations and update the domains
  subroutine test_halo_update_r8( domain_type )
    character(len=*), intent(in) :: domain_type !< the domain type that will be tested
    real(kind=r8_kind), allocatable, dimension(:,:,:) :: xr8, x1r8, x2r8, x3r8, x4r8
    real(kind=r8_kind), allocatable, dimension(:,:,:) :: yr8, y1r8, y2r8, y3r8, y4r8
    type(domain2D) :: domain
    real(kind=r8_kind),    allocatable :: global1r8(:,:,:), global2r8(:,:,:), globalr8(:,:,:)
    logical, allocatable :: maskmap(:,:)
    integer              :: shift, i, xhalo, yhalo
    logical              :: is_symmetry, folded_south, folded_west, folded_east
    integer              :: is, ie, js, je, isd, ied, jsd, jed
    integer :: pe, npes

    pe = mpp_pe()
    npes = mpp_npes()
    ! when testing maskmap option, nx*ny should be able to be divided by both npes and npes+1
    if ((domain_type == 'Masked') .OR. (domain_type == 'Masked symmetry')) then
      if((mod(nx*ny, npes) .NE. 0) .OR. (mod(nx*ny, npes+1) .NE. 0)) then
        call mpp_error(NOTE,'test_halo_update_r8: nx*ny can not be divided by both npes and npes+1, '//&
                       'Masked test_halo_update_r8 will not be tested')
        return
      end if
    end if

    if(trim(domain_type) == 'Folded xy_halo' ) then
      xhalo = max(whalo, ehalo); yhalo = max(shalo, nhalo)
      allocate(globalr8(1-xhalo:nx+xhalo,1-yhalo:ny+yhalo,nz) )
    else
      allocate(globalr8(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )
    end if
    ! populate the global array
    globalr8 = 0.0
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          globalr8(i,j,k) = k + i*1e-3 + j*1e-6
        end do
      end do
    end do

    if(index(domain_type, 'symmetry') == 0) then
       is_symmetry = .false.
    else
       is_symmetry = .true.
    end if
    select case(trim(domain_type))
      case( 'Simple', 'Simple symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, name=trim(domain_type), symmetry = is_symmetry )
      case( 'Cyclic', 'Cyclic symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,        &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=trim(domain_type), symmetry = is_symmetry )
        globalr8(1-whalo:0,                 1:ny,:) = globalr8(nx-whalo+1:nx,             1:ny,:)
        globalr8(nx+1:nx+ehalo,             1:ny,:) = globalr8(1:ehalo,                   1:ny,:)
        globalr8(1-whalo:nx+ehalo,     1-shalo:0,:) = globalr8(1-whalo:nx+ehalo, ny-shalo+1:ny,:)
        globalr8(1-whalo:nx+ehalo, ny+1:ny+nhalo,:) = globalr8(1-whalo:nx+ehalo,       1:nhalo,:)
      case( 'Folded-north', 'Folded-north symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
             name=domain_type, symmetry = is_symmetry  )
        call fill_folded_north_halo(globalr8, 0, 0, 0, 0, 1)
      case( 'Folded-south symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_SOUTH_EDGE, &
             name=domain_type, symmetry = is_symmetry  )
        call fill_folded_south_halo(globalr8, 0, 0, 0, 0, 1)
      case( 'Folded-west symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=FOLD_WEST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=domain_type, symmetry = is_symmetry  )
        call fill_folded_west_halo(globalr8, 0, 0, 0, 0, 1)
      case( 'Folded-east symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=FOLD_EAST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=domain_type, symmetry = is_symmetry  )
        call fill_folded_east_halo(globalr8, 0, 0, 0, 0, 1)
      case( 'Folded xy_halo' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=xhalo, yhalo=yhalo,   &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=domain_type, symmetry = is_symmetry  )
        globalr8(1-xhalo:0,                1:ny,:) = globalr8(nx-xhalo+1:nx,                   1:ny,:)
        globalr8(nx+1:nx+xhalo,            1:ny,:) = globalr8(1:xhalo,                         1:ny,:)
        globalr8(1-xhalo:nx+xhalo,ny+1:ny+yhalo,:) = globalr8(nx+xhalo:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
      case( 'Masked', 'Masked symmetry' )
      !with fold and cyclic, assign to npes+1 and mask out the top-rightdomain
        call mpp_define_layout( (/1,nx,1,ny/), npes+1, layout )
        allocate( maskmap(layout(1),layout(2)) )
        maskmap(:,:) = .TRUE.; maskmap(layout(1),layout(2)) = .FALSE.
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
          shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
          maskmap=maskmap, name=domain_type, symmetry = is_symmetry  )
        deallocate(maskmap)
       !we need to zero out the globalr8 data on the missing domain.
       !this logic assumes top-right, in an even division
        if( mod(nx,layout(1)).NE.0 .OR. mod(ny,layout(2)).NE.0 )call mpp_error( FATAL, &
             'TEST_MPP_DOMAINS: test for masked domains needs (nx,ny) to divide evenly on npes+1 PEs.' )
        globalr8(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny,:) = 0
        call fill_folded_north_halo(globalr8, 0, 0, 0, 0, 1)
      case default
        call mpp_error( FATAL, 'test_halo_update_r8: no such test: '//domain_type )
    end select

    ! define the arrays
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain( domain, isd, ied, jsd, jed )
    allocate(xr8(isd:ied,jsd:jed,nz) )
    allocate(x1r8(isd:ied,jsd:jed,nz) )
    allocate(x2r8(isd:ied,jsd:jed,nz) )
    allocate(x3r8(isd:ied,jsd:jed,nz) )
    allocate(x4r8(isd:ied,jsd:jed,nz) )
    xr8(:,:,:) = 0.0
    xr8 (is:ie,js:je,:) = globalr8(is:ie,js:je,:)
    x1r8 = xr8; x2r8 = xr8; x3r8 = xr8; x4r8 = xr8

    ! update the halo region
    id = mpp_clock_id( domain_type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( xr8, domain )
    call mpp_clock_end  (id)
    call compare_checksums( xr8, globalr8(isd:ied,jsd:jed,:), domain_type )

    ! update part of the halo region
    id = mpp_clock_id( domain_type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1r8, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2r8, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3r8, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4r8, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x1r8(is:ied,js:jed,:), globalr8(is:ied,js:jed,:), domain_type//' partial x1r8' )
    call compare_checksums( x2r8(is:ied,js:jed,:), globalr8(is:ied,js:jed,:), domain_type//' partial x2r8' )
    call compare_checksums( x3r8(is:ied,js:jed,:), globalr8(is:ied,js:jed,:), domain_type//' partial x3r8' )
    call compare_checksums( x4r8(is:ied,js:jed,:), globalr8(is:ied,js:jed,:), domain_type//' partial x4r8' )

    !--- test vector update for FOLDED and MASKED case.
    if(domain_type == 'Simple' .or. domain_type == 'Simple symmetry' .or. domain_type == 'Cyclic' .or. domain_type == 'Cyclic symmetry') then
      deallocate(xr8,x1r8,x2r8,x3r8,x4r8)
      return
    end if

    !------------------------------------------------------------------
    !              vector update : BGRID_NE
    !------------------------------------------------------------------
    shift = 0
    if(is_symmetry) then
      shift = 1
      deallocate(globalr8)
      allocate(globalr8(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
      globalr8(:,:,:) = 0.0
      do k = 1,nz
        do j = 1,ny+1
          do i = 1,nx+1
            globalr8(i,j,k) = k + i*1e-3 + j*1e-6
          end do
        end do
      end do
      if(domain_type == 'Masked symmetry') then
        globalr8(nx-nx/layout(1)+1:nx+1,ny-ny/layout(2)+1:ny+1,:) = 0.0
      endif
      deallocate(xr8, x1r8, x2r8, x3r8, x4r8)
      allocate( xr8(isd:ied+1,jsd:jed+1,nz) )
      allocate( x1r8(isd:ied+1,jsd:jed+1,nz) )
      allocate( x2r8(isd:ied+1,jsd:jed+1,nz) )
      allocate( x3r8(isd:ied+1,jsd:jed+1,nz) )
      allocate( x4r8(isd:ied+1,jsd:jed+1,nz) )
    endif

    folded_south = .false.
    folded_west  = .false.
    folded_east  = .false.
    select case (domain_type)
    case ('Folded-north', 'Masked')
      ! fill in folded north edge, cyclic east and west edge
      call fill_folded_north_halo(globalr8, 1, 1, 0, 0, -1)
    case ('Folded xy_halo')
      ! fill in folded north edge, cyclic east and west edge
      globalr8(1-xhalo:0,                  1:ny,:) =  globalr8(nx-xhalo+1:nx,                     1:ny,:)
      globalr8(nx+1:nx+xhalo,              1:ny,:) =  globalr8(1:xhalo,                           1:ny,:)
      globalr8(1-xhalo:nx+xhalo-1,ny+1:ny+yhalo,:) = -globalr8(nx+xhalo-1:1-xhalo:-1,ny-1:ny-yhalo:-1,:)
      globalr8(nx+xhalo,          ny+1:ny+yhalo,:) = -globalr8(nx-xhalo,             ny-1:ny-yhalo:-1,:)
    case ('Folded-north symmetry', 'Masked symmetry' )
      call fill_folded_north_halo(globalr8, 1, 1, 1, 1, -1)
    case ('Folded-south symmetry' )
      folded_south = .true.
      call fill_folded_south_halo(globalr8, 1, 1, 1, 1, -1)
    case ('Folded-west symmetry' )
      folded_west = .true.
      call fill_folded_west_halo(globalr8, 1, 1, 1, 1, -1)
    case ('Folded-east symmetry' )
      folded_east = .true.
      call fill_folded_east_halo(globalr8, 1, 1, 1, 1, -1)
    case default
      call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//domain_type )
    end select

    xr8(:,:,:) = 0.0
    xr8(is:ie+shift,js:je+shift,:) = globalr8(is:ie+shift,js:je+shift,:)
    !set up yr8 array
    allocate( yr8 (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y1r8(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y2r8(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y3r8(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y4r8(isd:ied+shift,jsd:jed+shift,nz) )
    yr8 = xr8; x1r8 = xr8; x2r8 = xr8; x3r8 = xr8; x4r8 = xr8
    yr8 = xr8; y1r8 = xr8; y2r8 = xr8; y3r8 = xr8; y4r8 = xr8

    id = mpp_clock_id( domain_type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( xr8,  yr8,  domain, gridtype=BGRID_NE)
    call mpp_update_domains( x1r8, y1r8, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x2r8, y2r8, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x3r8, y3r8, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x4r8, y4r8, domain, gridtype=BGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    ! redundant points must be equal and opposite

    if(folded_south) then
      globalr8(nx/2+shift,                1,:) = 0.0  !pole points must have 0 velocity
      globalr8(nx+shift  ,                1,:) = 0.0  !pole points must have 0 velocity
      globalr8(nx/2+1+shift:nx-1+shift,   1,:) = -globalr8(nx/2-1+shift:1+shift:-1, 1,:)
      globalr8(1-whalo:shift,             1,:) = -globalr8(nx-whalo+1:nx+shift,     1,:)
      globalr8(nx+1+shift:nx+ehalo+shift, 1,:) = -globalr8(1+shift:ehalo+shift,     1,:)
      !--- the following will fix the +0/-0 problem on altix
      if(shalo >0) globalr8(shift,1,:) = 0.0  !pole points must have 0 velocity
    else if(folded_west) then
      globalr8(1, ny/2+shift, :) = 0.0 !pole points must have 0 velocity
      globalr8(1, ny+shift,   :) = 0.0 !pole points must have 0 velocity
      globalr8(1, ny/2+1+shift:ny-1+shift,   :) = -globalr8(1, ny/2-1+shift:1+shift:-1, :)
      globalr8(1, 1-shalo:shift,             :) = -globalr8(1, ny-shalo+1:ny+shift,     :)
      globalr8(1, ny+1+shift:ny+nhalo+shift, :) = -globalr8(1, 1+shift:nhalo+shift,     :)
      !--- the following will fix the +0/-0 problem on altix
      if(whalo>0) globalr8(1, shift, :) = 0.0  !pole points must have 0 velocity
    else if(folded_east) then
      globalr8(nx+shift, ny/2+shift, :) = 0.0 !pole points must have 0 velocity
      globalr8(nx+shift, ny+shift,   :) = 0.0 !pole points must have 0 velocity
      globalr8(nx+shift, ny/2+1+shift:ny-1+shift,   :) = -globalr8(nx+shift, ny/2-1+shift:1+shift:-1, :)
      globalr8(nx+shift, 1-shalo:shift,             :) = -globalr8(nx+shift, ny-shalo+1:ny+shift,     :)
      globalr8(nx+shift, ny+1+shift:ny+nhalo+shift, :) = -globalr8(nx+shift, 1+shift:nhalo+shift,     :)
      if(ehalo >0) globalr8(nx+shift, shift, :) = 0.0  !pole points must have 0 velocity
    else
      globalr8(nx/2+shift,                ny+shift,:) = 0.0  !pole points must have 0 velocity
      globalr8(nx+shift  ,                ny+shift,:) = 0.0  !pole points must have 0 velocity
      globalr8(nx/2+1+shift:nx-1+shift,   ny+shift,:) = -globalr8(nx/2-1+shift:1+shift:-1, ny+shift,:)
      if (domain_type == 'Folded xy_halo') then
        globalr8(1-xhalo:shift,             ny+shift,:) = -globalr8(nx-xhalo+1:nx+shift,     ny+shift,:)
        globalr8(nx+1+shift:nx+xhalo+shift, ny+shift,:) = -globalr8(1+shift:xhalo+shift,     ny+shift,:)
      else
        globalr8(1-whalo:shift,             ny+shift,:) = -globalr8(nx-whalo+1:nx+shift,     ny+shift,:)
        globalr8(nx+1+shift:nx+ehalo+shift, ny+shift,:) = -globalr8(1+shift:ehalo+shift,     ny+shift,:)
      end if
      !--- the following will fix the +0/-0 problem on altix
      if(nhalo >0) globalr8(shift,ny+shift,:) = 0.0  !pole points must have 0 velocity
    endif

    call compare_checksums( xr8,  globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE xr8' )
    call compare_checksums( yr8,  globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE yr8' )
    call compare_checksums( x1r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x1r8' )
    call compare_checksums( x2r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x2r8' )
    call compare_checksums( x3r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x3r8' )
    call compare_checksums( x4r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x4r8' )
    call compare_checksums( y1r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y1r8' )
    call compare_checksums( y2r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y2r8' )
    call compare_checksums( y3r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y3r8' )
    call compare_checksums( y4r8, globalr8(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y4r8' )

    deallocate(globalr8, xr8, x1r8, x2r8, x3r8, x4r8, yr8, y1r8, y2r8, y3r8, y4r8)

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- global1r8 is x-component and global2r8 is y-component
    if(domain_type == 'Folded xy_halo') then
      allocate(global1r8(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
      allocate(global2r8(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
    else
      allocate(global1r8(1-whalo:nx+ehalo+shift, 1-shalo:ny+nhalo, nz))
      allocate(global2r8(1-whalo:nx+ehalo, 1-shalo:ny+nhalo+shift, nz))
    end if
    allocate(xr8 (isd:ied+shift,jsd:jed,nz), yr8 (isd:ied,jsd:jed+shift,nz) )
    allocate(x1r8(isd:ied+shift,jsd:jed,nz), y1r8(isd:ied,jsd:jed+shift,nz) )
    allocate(x2r8(isd:ied+shift,jsd:jed,nz), y2r8(isd:ied,jsd:jed+shift,nz) )
    allocate(x3r8(isd:ied+shift,jsd:jed,nz), y3r8(isd:ied,jsd:jed+shift,nz) )
    allocate(x4r8(isd:ied+shift,jsd:jed,nz), y4r8(isd:ied,jsd:jed+shift,nz) )

    global1r8(:,:,:) = 0.0
    global2r8(:,:,:) = 0.0
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx+shift
          global1r8(i,j,k) = k + i*1e-3 + j*1e-6
        end do
      end do
      do j = 1,ny+shift
        do i = 1,nx
          global2r8(i,j,k) = k + i*1e-3 + j*1e-6
        end do
      end do
    end do

    if(domain_type == 'Masked' .or. domain_type == 'Masked symmetry') then
      global1r8(nx-nx/layout(1)+1:nx+shift,ny-ny/layout(2)+1:ny,:) = 0
      global2r8(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny+shift,:) = 0
    end if

    select case (domain_type)
    case ('Folded-north', 'Masked')
      !fill in folded north edge, cyclic east and west edge
      call fill_folded_north_halo(global1r8, 1, 0, 0, 0, -1)
      call fill_folded_north_halo(global2r8, 0, 1, 0, 0, -1)
    case ('Folded xy_halo')
      global1r8(1-xhalo:0,                   1:ny,:) =  global1r8(nx-xhalo+1:nx,                     1:ny,:)
      global1r8(nx+1:nx+xhalo,               1:ny,:) =  global1r8(1:xhalo,                           1:ny,:)
     global2r8(1-xhalo:0,                   1:ny,:) =  global2r8(nx-xhalo+1:nx,                     1:ny,:)
      global2r8(nx+1:nx+xhalo,               1:ny,:) =  global2r8(1:xhalo,                           1:ny,:)
      global1r8(1-xhalo:nx+xhalo-1, ny+1:ny+yhalo,:) = -global1r8(nx+xhalo-1:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
      global1r8(nx+xhalo,           ny+1:ny+yhalo,:) = -global1r8(nx-xhalo,              ny:ny-yhalo+1:-1,:)
      global2r8(1-xhalo:nx+xhalo,   ny+1:ny+yhalo,:) = -global2r8(nx+xhalo:1-xhalo:-1,   ny-1:ny-yhalo:-1,:)
    case ('Folded-north symmetry')
      call fill_folded_north_halo(global1r8, 1, 0, 1, 0, -1)
      call fill_folded_north_halo(global2r8, 0, 1, 0, 1, -1)
    case ('Folded-south symmetry')
      call fill_folded_south_halo(global1r8, 1, 0, 1, 0, -1)
      call fill_folded_south_halo(global2r8, 0, 1, 0, 1, -1)
    case ('Folded-west symmetry')
      call fill_folded_west_halo(global1r8, 1, 0, 1, 0, -1)
      call fill_folded_west_halo(global2r8, 0, 1, 0, 1, -1)
    case ('Folded-east symmetry')
      call fill_folded_east_halo(global1r8, 1, 0, 1, 0, -1)
      call fill_folded_east_halo(global2r8, 0, 1, 0, 1, -1)
    case default
       call mpp_error( FATAL, 'test_halo_update_r8: invalid test name: '//domain_type )
    end select

    xr8(:,:,:) = 0.; yr8(:,:,:) = 0.
    xr8(is:ie+shift,js:je,      :) = global1r8(is:ie+shift,js:je,      :)
    yr8(is:ie      ,js:je+shift,:) = global2r8(is:ie,      js:je+shift,:)
    x1r8 = xr8; x2r8 = xr8; x3r8 = xr8; x4r8 = xr8
    y1r8 = yr8; y2r8 = yr8; y3r8 = yr8; y4r8 = yr8

    id = mpp_clock_id( domain_type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( xr8,  yr8,  domain, gridtype=CGRID_NE)
    call mpp_update_domains( x1r8, y1r8, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x2r8, y2r8, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x3r8, y3r8, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x4r8, y4r8, domain, gridtype=CGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    ! redundant points must be equal and opposite
    if(folded_south) then
      global2r8(nx/2+1:nx,     1,:) = -global2r8(nx/2:1:-1, 1,:)
      global2r8(1-whalo:0,     1,:) = -global2r8(nx-whalo+1:nx, 1, :)
      global2r8(nx+1:nx+ehalo, 1,:) = -global2r8(1:ehalo,       1, :)
    else if(folded_west) then
      global1r8(1, ny/2+1:ny,     :) = -global1r8(1, ny/2:1:-1,     :)
      global1r8(1, 1-shalo:0,     :) = -global1r8(1, ny-shalo+1:ny, :)
      global1r8(1, ny+1:ny+nhalo, :) = -global1r8(1, 1:nhalo,       :)
    else if(folded_east) then
      global1r8(nx+shift, ny/2+1:ny,     :) = -global1r8(nx+shift, ny/2:1:-1,     :)
      global1r8(nx+shift, 1-shalo:0,     :) = -global1r8(nx+shift, ny-shalo+1:ny, :)
      global1r8(nx+shift, ny+1:ny+nhalo, :) = -global1r8(nx+shift, 1:nhalo,       :)
    else
      global2r8(nx/2+1:nx,     ny+shift,:) = -global2r8(nx/2:1:-1, ny+shift,:)
      if(domain_type == 'Folded xy_halo') then
        global2r8(1-xhalo:0,     ny+shift,:) = -global2r8(nx-xhalo+1:nx, ny+shift,:)
        global2r8(nx+1:nx+xhalo, ny+shift,:) = -global2r8(1:xhalo,       ny+shift,:)
      else
        global2r8(1-whalo:0,     ny+shift,:) = -global2r8(nx-whalo+1:nx, ny+shift,:)
        global2r8(nx+1:nx+ehalo, ny+shift,:) = -global2r8(1:ehalo,       ny+shift,:)
      end if
    endif

    call compare_checksums( xr8,  global1r8(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE xr8' )
    call compare_checksums( yr8,  global2r8(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE yr8' )
    call compare_checksums( x1r8, global1r8(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x1r8' )
    call compare_checksums( x2r8, global1r8(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x2r8' )
    call compare_checksums( x3r8, global1r8(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x3r8' )
    call compare_checksums( x4r8, global1r8(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x4r8' )
    call compare_checksums( y1r8, global2r8(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y1r8' )
    call compare_checksums( y2r8, global2r8(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y2r8' )
    call compare_checksums( y3r8, global2r8(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y3r8' )
    call compare_checksums( y4r8, global2r8(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y4r8' )

    deallocate(global1r8, global2r8, xr8, x1r8, x2r8, x3r8, x4r8, yr8, y1r8, y2r8, y3r8, y4r8)

  end subroutine test_halo_update_r8

  !> Perform simple addition on 32-bit real arrays in different domain configurations and update the domains
  subroutine test_halo_update_r4( domain_type )
   character(len=*), intent(in) :: domain_type !< the domain type that will be tested
   real(kind=r4_kind), allocatable, dimension(:,:,:) :: xr4, x1r4, x2r4, x3r4, x4r4
   real(kind=r4_kind), allocatable, dimension(:,:,:) :: yr4, y1r4, y2r4, y3r4, y4r4
   type(domain2D) :: domain
   real(kind=r4_kind),    allocatable :: global1r4(:,:,:), global2r4(:,:,:), globalr4(:,:,:)
   logical, allocatable :: maskmap(:,:)
   integer              :: shift, i, xhalo, yhalo
   logical              :: is_symmetry, folded_south, folded_west, folded_east
   integer              :: is, ie, js, je, isd, ied, jsd, jed
   integer :: pe, npes

   pe = mpp_pe()
   npes = mpp_npes()

   ! when testing maskmap option, nx*ny should be able to be divided by both npes and npes+1
   if((domain_type == 'Masked') .OR. (domain_type == 'Masked symmetry')) then
      if((mod(nx*ny, npes) .NE. 0) .OR. (mod(nx*ny, npes+1) .NE. 0)) then
        call mpp_error(NOTE,'test_halo_update_r4: nx*ny can not be divided by both npes and npes+1, '//&
              'Masked test_halo_update_r4 will not be tested')
         return
      end if
   end if

   if(trim(domain_type) == 'Folded xy_halo' ) then
     xhalo = max(whalo, ehalo); yhalo = max(shalo, nhalo)
     allocate(globalr4(1-xhalo:nx+xhalo,1-yhalo:ny+yhalo,nz) )
   else
     allocate(globalr4(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )
   end if

   globalr4 = 0
   do k = 1,nz
     do j = 1,ny
       do i = 1,nx
         globalr4(i,j,k) = k + i*1e-3 + j*1e-6
       end do
     end do
   end do

   if(index(domain_type, 'symmetry') == 0) then
     is_symmetry = .false.
   else
     is_symmetry = .true.
   end if
   select case(trim(domain_type))
   case( 'Simple', 'Simple symmetry' )
     call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                             shalo=shalo, nhalo=nhalo, name=trim(domain_type), symmetry = is_symmetry )
   case( 'Cyclic', 'Cyclic symmetry' )
     call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,        &
          shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, &
            name=trim(domain_type), symmetry = is_symmetry )
     globalr4(1-whalo:0,                 1:ny,:) = globalr4(nx-whalo+1:nx,             1:ny,:)
     globalr4(nx+1:nx+ehalo,             1:ny,:) = globalr4(1:ehalo,                   1:ny,:)
     globalr4(1-whalo:nx+ehalo,     1-shalo:0,:) = globalr4(1-whalo:nx+ehalo, ny-shalo+1:ny,:)
     globalr4(1-whalo:nx+ehalo, ny+1:ny+nhalo,:) = globalr4(1-whalo:nx+ehalo,       1:nhalo,:)
   case( 'Folded-north', 'Folded-north symmetry' )
     call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
            shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
            name=domain_type, symmetry = is_symmetry  )
     call fill_folded_north_halo(globalr4, 0, 0, 0, 0, 1)
   case( 'Folded-south symmetry' )
     call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
            shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_SOUTH_EDGE, &
            name=domain_type, symmetry = is_symmetry  )
     call fill_folded_south_halo(globalr4, 0, 0, 0, 0, 1)
   case( 'Folded-west symmetry' )
     call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
            shalo=shalo, nhalo=nhalo, xflags=FOLD_WEST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
            name=domain_type, symmetry = is_symmetry  )
     call fill_folded_west_halo(globalr4, 0, 0, 0, 0, 1)
   case( 'Folded-east symmetry' )
     call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
            shalo=shalo, nhalo=nhalo, xflags=FOLD_EAST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
            name=domain_type, symmetry = is_symmetry  )
     call fill_folded_east_halo(globalr4, 0, 0, 0, 0, 1)
   case( 'Folded xy_halo' )
     call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=xhalo, yhalo=yhalo,   &
            xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=domain_type, symmetry = is_symmetry)
       globalr4(1-xhalo:0,                1:ny,:) = globalr4(nx-xhalo+1:nx,                   1:ny,:)
       globalr4(nx+1:nx+xhalo,            1:ny,:) = globalr4(1:xhalo,                         1:ny,:)
       globalr4(1-xhalo:nx+xhalo,ny+1:ny+yhalo,:) = globalr4(nx+xhalo:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
   case( 'Masked', 'Masked symmetry' )
     ! with fold and cyclic, assign to npes+1 and mask out the top-rightdomain
     call mpp_define_layout( (/1,nx,1,ny/), npes+1, layout )
     allocate( maskmap(layout(1),layout(2)) )
     maskmap(:,:) = .TRUE.; maskmap(layout(1),layout(2)) = .FALSE.
     call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
            shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
            maskmap=maskmap, name=domain_type, symmetry = is_symmetry  )
     deallocate(maskmap)
     ! we need to zero out the globalr4 data on the missing domain.
     ! this logic assumes top-right, in an even division
     if ( mod(nx,layout(1)).NE.0 .OR. mod(ny,layout(2)).NE.0) call mpp_error( FATAL, &
      'test_halo_update_r4: test for masked domains needs (nx,ny) to divide evenly on npes+1 PEs.')
       globalr4(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny,:) = 0
     call fill_folded_north_halo(globalr4, 0, 0, 0, 0, 1)
   case default
     call mpp_error( FATAL, 'test_halo_update_r4: '//domain_type//' is not a valid test.')
   end select
   ! define the arrays
   call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
   call mpp_get_data_domain( domain, isd, ied, jsd, jed )
   allocate(xr4(isd:ied,jsd:jed,nz) )
   allocate(x1r4(isd:ied,jsd:jed,nz) )
   allocate(x2r4(isd:ied,jsd:jed,nz) )
   allocate(x3r4(isd:ied,jsd:jed,nz) )
   allocate(x4r4(isd:ied,jsd:jed,nz) )
   xr4 = 0.
   xr4 (is:ie,js:je,:) = globalr4(is:ie,js:je,:)
   x1r4 = xr4; x2r4 = xr4; x3r4 = xr4; x4r4 = xr4
   ! update the halo region
   id = mpp_clock_id( domain_type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
   call mpp_clock_begin(id)
   call mpp_update_domains( xr4, domain )
   call mpp_clock_end  (id)
   call compare_checksums( xr4, globalr4(isd:ied,jsd:jed,:), domain_type )
   ! update part of the halo region
   id = mpp_clock_id( domain_type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
   call mpp_clock_begin(id)
   call mpp_update_domains( x1r4, domain, NUPDATE+EUPDATE, complete=.false. )
   call mpp_update_domains( x2r4, domain, NUPDATE+EUPDATE, complete=.false. )
   call mpp_update_domains( x3r4, domain, NUPDATE+EUPDATE, complete=.false. )
   call mpp_update_domains( x4r4, domain, NUPDATE+EUPDATE, complete=.true. )
   call mpp_clock_end  (id)
   call compare_checksums( x1r4(is:ied,js:jed,:), globalr4(is:ied,js:jed,:), domain_type//' partial x1r4' )
   call compare_checksums( x2r4(is:ied,js:jed,:), globalr4(is:ied,js:jed,:), domain_type//' partial x2r4' )
   call compare_checksums( x3r4(is:ied,js:jed,:), globalr4(is:ied,js:jed,:), domain_type//' partial x3r4' )
   call compare_checksums( x4r4(is:ied,js:jed,:), globalr4(is:ied,js:jed,:), domain_type//' partial x4r4' )

   !--- test vector update for FOLDED and MASKED case.
   if(domain_type == 'Simple' .or. domain_type == 'Simple symmetry' .or. domain_type == 'Cyclic' .or. domain_type == 'Cyclic symmetry') then
      deallocate(xr4,x1r4,x2r4,x3r4,x4r4)
      return
   end if

   !------------------------------------------------------------------
   !              vector update : BGRID_NE
   !------------------------------------------------------------------
   shift = 0
   if (is_symmetry) then
     shift = 1
     deallocate(globalr4)
     allocate(globalr4(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
     globalr4 = 0.0
     do k = 1,nz
       do j = 1,ny+1
         do i = 1,nx+1
           globalr4(i,j,k) = k + i*1e-3 + j*1e-6
         end do
       end do
     end do
     if (domain_type == 'Masked symmetry') then
       globalr4(nx-nx/layout(1)+1:nx+1,ny-ny/layout(2)+1:ny+1,:) = 0
     endif
     deallocate(xr4, x1r4, x2r4, x3r4, x4r4)
     allocate( xr4(isd:ied+1,jsd:jed+1,nz) )
     allocate( x1r4(isd:ied+1,jsd:jed+1,nz) )
     allocate( x2r4(isd:ied+1,jsd:jed+1,nz) )
     allocate( x3r4(isd:ied+1,jsd:jed+1,nz) )
     allocate( x4r4(isd:ied+1,jsd:jed+1,nz) )
   endif

   folded_south = .false.
   folded_west  = .false.
   folded_east  = .false.
   select case (domain_type)
   case ('Folded-north', 'Masked')
     ! fill in folded north edge, cyclic east and west edge
     call fill_folded_north_halo(globalr4, 1, 1, 0, 0, -1)
   case ('Folded xy_halo')
     ! fill in folded north edge, cyclic east and west edge
     globalr4(1-xhalo:0,                  1:ny,:) =  globalr4(nx-xhalo+1:nx,                     1:ny,:)
     globalr4(nx+1:nx+xhalo,              1:ny,:) =  globalr4(1:xhalo,                           1:ny,:)
     globalr4(1-xhalo:nx+xhalo-1,ny+1:ny+yhalo,:) = -globalr4(nx+xhalo-1:1-xhalo:-1,ny-1:ny-yhalo:-1,:)
     globalr4(nx+xhalo,          ny+1:ny+yhalo,:) = -globalr4(nx-xhalo,             ny-1:ny-yhalo:-1,:)
   case ('Folded-north symmetry', 'Masked symmetry' )
     call fill_folded_north_halo(globalr4, 1, 1, 1, 1, -1)
   case ('Folded-south symmetry' )
     folded_south = .true.
     call fill_folded_south_halo(globalr4, 1, 1, 1, 1, -1)
   case ('Folded-west symmetry' )
     folded_west = .true.
     call fill_folded_west_halo(globalr4, 1, 1, 1, 1, -1)
   case ('Folded-east symmetry' )
     folded_east = .true.
     call fill_folded_east_halo(globalr4, 1, 1, 1, 1, -1)
   case default
     call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//domain_type )
   end select

   xr4 = 0.0
   xr4(is:ie+shift,js:je+shift,:) = globalr4(is:ie+shift,js:je+shift,:)
   !set up yr4 array
   allocate( yr4 (isd:ied+shift,jsd:jed+shift,nz) )
   allocate( y1r4(isd:ied+shift,jsd:jed+shift,nz) )
   allocate( y2r4(isd:ied+shift,jsd:jed+shift,nz) )
   allocate( y3r4(isd:ied+shift,jsd:jed+shift,nz) )
   allocate( y4r4(isd:ied+shift,jsd:jed+shift,nz) )
   yr4 = xr4; x1r4 = xr4; x2r4 = xr4; x3r4 = xr4; x4r4 = xr4
   yr4 = xr4; y1r4 = xr4; y2r4 = xr4; y3r4 = xr4; y4r4 = xr4

   id = mpp_clock_id( domain_type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
   call mpp_clock_begin(id)
   call mpp_update_domains( xr4,  yr4,  domain, gridtype=BGRID_NE)
   call mpp_update_domains( x1r4, y1r4, domain, gridtype=BGRID_NE, complete=.false. )
   call mpp_update_domains( x2r4, y2r4, domain, gridtype=BGRID_NE, complete=.false. )
   call mpp_update_domains( x3r4, y3r4, domain, gridtype=BGRID_NE, complete=.false. )
   call mpp_update_domains( x4r4, y4r4, domain, gridtype=BGRID_NE, complete=.true.  )
   call mpp_clock_end  (id)

   ! redundant points must be equal and opposite
   if(folded_south) then
     globalr4(nx/2+shift,                1,:) = 0.  !pole points must have 0 velocity
     globalr4(nx+shift  ,                1,:) = 0.  !pole points must have 0 velocity
     globalr4(nx/2+1+shift:nx-1+shift,   1,:) = -globalr4(nx/2-1+shift:1+shift:-1, 1,:)
     globalr4(1-whalo:shift,             1,:) = -globalr4(nx-whalo+1:nx+shift,     1,:)
     globalr4(nx+1+shift:nx+ehalo+shift, 1,:) = -globalr4(1+shift:ehalo+shift,     1,:)

     if(shalo >0) globalr4(shift,1,:) = 0.0  !pole points must have 0 velocity
   else if(folded_west) then
     globalr4(1, ny/2+shift, :) = 0.0 !pole points must have 0 velocity
     globalr4(1, ny+shift,   :) = 0.0 !pole points must have 0 velocity
     globalr4(1, ny/2+1+shift:ny-1+shift,   :) = -globalr4(1, ny/2-1+shift:1+shift:-1, :)
     globalr4(1, 1-shalo:shift,             :) = -globalr4(1, ny-shalo+1:ny+shift,     :)
     globalr4(1, ny+1+shift:ny+nhalo+shift, :) = -globalr4(1, 1+shift:nhalo+shift,     :)
     !--- the following will fix the +0/-0 problem on altix
     if(whalo>0) globalr4(1, shift, :) = 0.0  !pole points must have 0 velocity
   else if(folded_east) then
     globalr4(nx+shift, ny/2+shift, :) = 0.0 !pole points must have 0 velocity
     globalr4(nx+shift, ny+shift,   :) = 0.0 !pole points must have 0 velocity
     globalr4(nx+shift, ny/2+1+shift:ny-1+shift,   :) = -globalr4(nx+shift, ny/2-1+shift:1+shift:-1, :)
     globalr4(nx+shift, 1-shalo:shift,             :) = -globalr4(nx+shift, ny-shalo+1:ny+shift,     :)
     globalr4(nx+shift, ny+1+shift:ny+nhalo+shift, :) = -globalr4(nx+shift, 1+shift:nhalo+shift,     :)
     if(ehalo >0) globalr4(nx+shift, shift, :) = 0.0  !pole points must have 0 velocity
   else
      globalr4(nx/2+shift,                ny+shift,:) = 0.0  !pole points must have 0 velocity
      globalr4(nx+shift  ,                ny+shift,:) = 0.0  !pole points must have 0 velocity
      globalr4(nx/2+1+shift:nx-1+shift,   ny+shift,:) = -globalr4(nx/2-1+shift:1+shift:-1, ny+shift,:)
      if(domain_type == 'Folded xy_halo') then
        globalr4(1-xhalo:shift,             ny+shift,:) = -globalr4(nx-xhalo+1:nx+shift,     ny+shift,:)
        globalr4(nx+1+shift:nx+xhalo+shift, ny+shift,:) = -globalr4(1+shift:xhalo+shift,     ny+shift,:)
      else
        globalr4(1-whalo:shift,             ny+shift,:) = -globalr4(nx-whalo+1:nx+shift,     ny+shift,:)
        globalr4(nx+1+shift:nx+ehalo+shift, ny+shift,:) = -globalr4(1+shift:ehalo+shift,     ny+shift,:)
      end if
      !--- the following will fix the +0/-0 problem on altix
      if(nhalo >0) globalr4(shift,ny+shift,:) = 0.0  !pole points must have 0 velocity
   endif

   call compare_checksums( xr4,  globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE xr4' )
   call compare_checksums( yr4,  globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE yr4' )
   call compare_checksums( x1r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x1r4' )
   call compare_checksums( x2r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x2r4' )
   call compare_checksums( x3r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x3r4' )
   call compare_checksums( x4r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE x4r4' )
   call compare_checksums( y1r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y1r4' )
   call compare_checksums( y2r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y2r4' )
   call compare_checksums( y3r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y3r4' )
   call compare_checksums( y4r4, globalr4(isd:ied+shift,jsd:jed+shift,:), domain_type//' BGRID_NE y4r4' )

   deallocate(globalr4, xr4, x1r4, x2r4, x3r4, x4r4, yr4, y1r4, y2r4, y3r4, y4r4)
   !------------------------------------------------------------------
   !              vector update : CGRID_NE
   !------------------------------------------------------------------
   !--- global1r4 is x-component and global2r4 is y-component
   if(domain_type == 'Folded xy_halo') then
     allocate(global1r4(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
     allocate(global2r4(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
   else
     allocate(global1r4(1-whalo:nx+ehalo+shift, 1-shalo:ny+nhalo, nz))
     allocate(global2r4(1-whalo:nx+ehalo, 1-shalo:ny+nhalo+shift, nz))
   end if
   allocate(xr4 (isd:ied+shift,jsd:jed,nz), yr4 (isd:ied,jsd:jed+shift,nz) )
   allocate(x1r4(isd:ied+shift,jsd:jed,nz), y1r4(isd:ied,jsd:jed+shift,nz) )
   allocate(x2r4(isd:ied+shift,jsd:jed,nz), y2r4(isd:ied,jsd:jed+shift,nz) )
   allocate(x3r4(isd:ied+shift,jsd:jed,nz), y3r4(isd:ied,jsd:jed+shift,nz) )
   allocate(x4r4(isd:ied+shift,jsd:jed,nz), y4r4(isd:ied,jsd:jed+shift,nz) )

   global1r4 = 0.0
   global2r4 = 0.0
   do k = 1,nz
     do j = 1,ny
       do i = 1,nx+shift
         global1r4(i,j,k) = k + i*1e-3 + j*1e-6
       end do
     end do
     do j = 1,ny+shift
       do i = 1,nx
         global2r4(i,j,k) = k + i*1e-3 + j*1e-6
       end do
     end do
   end do

   if (domain_type == 'Masked' .or. domain_type == 'Masked symmetry') then
     global1r4(nx-nx/layout(1)+1:nx+shift,ny-ny/layout(2)+1:ny,:) = 0
     global2r4(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny+shift,:) = 0
   end if

   select case (domain_type)
   case ('Folded-north', 'Masked')
     !fill in folded north edge, cyclic east and west edge
     call fill_folded_north_halo(global1r4, 1, 0, 0, 0, -1)
     call fill_folded_north_halo(global2r4, 0, 1, 0, 0, -1)
   case ('Folded xy_halo')
     global1r4(1-xhalo:0,                   1:ny,:) =  global1r4(nx-xhalo+1:nx,                     1:ny,:)
     global1r4(nx+1:nx+xhalo,               1:ny,:) =  global1r4(1:xhalo,                           1:ny,:)
     global2r4(1-xhalo:0,                   1:ny,:) =  global2r4(nx-xhalo+1:nx,                     1:ny,:)
     global2r4(nx+1:nx+xhalo,               1:ny,:) =  global2r4(1:xhalo,                           1:ny,:)
     global1r4(1-xhalo:nx+xhalo-1, ny+1:ny+yhalo,:) = -global1r4(nx+xhalo-1:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
     global1r4(nx+xhalo,           ny+1:ny+yhalo,:) = -global1r4(nx-xhalo,              ny:ny-yhalo+1:-1,:)
     global2r4(1-xhalo:nx+xhalo,   ny+1:ny+yhalo,:) = -global2r4(nx+xhalo:1-xhalo:-1,   ny-1:ny-yhalo:-1,:)
   case ('Folded-north symmetry')
     call fill_folded_north_halo(global1r4, 1, 0, 1, 0, -1)
     call fill_folded_north_halo(global2r4, 0, 1, 0, 1, -1)
   case ('Folded-south symmetry')
     call fill_folded_south_halo(global1r4, 1, 0, 1, 0, -1)
     call fill_folded_south_halo(global2r4, 0, 1, 0, 1, -1)
   case ('Folded-west symmetry')
     call fill_folded_west_halo(global1r4, 1, 0, 1, 0, -1)
     call fill_folded_west_halo(global2r4, 0, 1, 0, 1, -1)
   case ('Folded-east symmetry')
     call fill_folded_east_halo(global1r4, 1, 0, 1, 0, -1)
     call fill_folded_east_halo(global2r4, 0, 1, 0, 1, -1)
   case default
     call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//domain_type )
   end select

   xr4 = 0.0; yr4 = 0.0
   xr4(is:ie+shift,js:je,      :) = global1r4(is:ie+shift,js:je,      :)
   yr4(is:ie      ,js:je+shift,:) = global2r4(is:ie,      js:je+shift,:)
   x1r4 = xr4; x2r4 = xr4; x3r4 = xr4; x4r4 = xr4
   y1r4 = yr4; y2r4 = yr4; y3r4 = yr4; y4r4 = yr4

   id = mpp_clock_id( domain_type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
   call mpp_clock_begin(id)
   call mpp_update_domains( xr4,  yr4,  domain, gridtype=CGRID_NE)
   call mpp_update_domains( x1r4, y1r4, domain, gridtype=CGRID_NE, complete=.false. )
   call mpp_update_domains( x2r4, y2r4, domain, gridtype=CGRID_NE, complete=.false. )
   call mpp_update_domains( x3r4, y3r4, domain, gridtype=CGRID_NE, complete=.false. )
   call mpp_update_domains( x4r4, y4r4, domain, gridtype=CGRID_NE, complete=.true.  )
   call mpp_clock_end  (id)

   ! redundant points must be equal and opposite
   if(folded_south) then
     global2r4(nx/2+1:nx,     1,:) = -global2r4(nx/2:1:-1, 1,:)
     global2r4(1-whalo:0,     1,:) = -global2r4(nx-whalo+1:nx, 1, :)
     global2r4(nx+1:nx+ehalo, 1,:) = -global2r4(1:ehalo,       1, :)
   else if(folded_west) then
     global1r4(1, ny/2+1:ny,     :) = -global1r4(1, ny/2:1:-1,     :)
     global1r4(1, 1-shalo:0,     :) = -global1r4(1, ny-shalo+1:ny, :)
     global1r4(1, ny+1:ny+nhalo, :) = -global1r4(1, 1:nhalo,       :)
   else if(folded_east) then
     global1r4(nx+shift, ny/2+1:ny,     :) = -global1r4(nx+shift, ny/2:1:-1,     :)
     global1r4(nx+shift, 1-shalo:0,     :) = -global1r4(nx+shift, ny-shalo+1:ny, :)
     global1r4(nx+shift, ny+1:ny+nhalo, :) = -global1r4(nx+shift, 1:nhalo,       :)
   else
     global2r4(nx/2+1:nx,     ny+shift,:) = -global2r4(nx/2:1:-1, ny+shift,:)
     if (domain_type == 'Folded xy_halo') then
       global2r4(1-xhalo:0,     ny+shift,:) = -global2r4(nx-xhalo+1:nx, ny+shift,:)
       global2r4(nx+1:nx+xhalo, ny+shift,:) = -global2r4(1:xhalo,       ny+shift,:)
     else
       global2r4(1-whalo:0,     ny+shift,:) = -global2r4(nx-whalo+1:nx, ny+shift,:)
       global2r4(nx+1:nx+ehalo, ny+shift,:) = -global2r4(1:ehalo,       ny+shift,:)
     end if
   endif

   call compare_checksums( xr4,  global1r4(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE xr4' )
   call compare_checksums( yr4,  global2r4(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE yr4' )
   call compare_checksums( x1r4, global1r4(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x1r4' )
   call compare_checksums( x2r4, global1r4(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x2r4' )
   call compare_checksums( x3r4, global1r4(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x3r4' )
   call compare_checksums( x4r4, global1r4(isd:ied+shift,jsd:jed,      :), domain_type//' CGRID_NE x4r4' )
   call compare_checksums( y1r4, global2r4(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y1r4' )
   call compare_checksums( y2r4, global2r4(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y2r4' )
   call compare_checksums( y3r4, global2r4(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y3r4' )
   call compare_checksums( y4r4, global2r4(isd:ied,      jsd:jed+shift,:), domain_type//' CGRID_NE y4r4' )

   deallocate(global1r4, global2r4, xr4, x1r4, x2r4, x3r4, x4r4, yr4, y1r4, y2r4, y3r4, y4r4)

 end subroutine test_halo_update_r4

 !> test a domain update of a 64-bit 3D array on a 9-pe subset of total allotted pes
 !> @note requires at least 16 pes
 subroutine test_subset_update_r8( )
   real(kind=r8_kind), allocatable, dimension(:,:,:) :: x
   type(domain2D) :: domain
   real(kind=r8_kind), allocatable :: global(:,:,:)
   integer              :: i, xhalo, yhalo
   integer              :: is, ie, js, je, isd, ied, jsd, jed
   integer :: pes9(9)=(/1,2,4,6,8,10,12,13,15/)
   integer :: ni, nj
   integer :: pe, npes

   pe = mpp_pe()
   npes = mpp_npes()

   call mpp_declare_pelist(pes9)
   if(any(mpp_pe()==pes9)) then
     call mpp_set_current_pelist(pes9)
     layout = (/3,3/)
     ni = 3; nj =3
     call mpp_define_domains((/1,ni,1,nj/), layout, domain, xhalo=1 &
       &, yhalo=1, xflags=CYCLIC_GLOBAL_DOMAIN, yflags&
       &=CYCLIC_GLOBAL_DOMAIN, name='subset domain')
     call mpp_get_compute_domain(domain, is, ie, js, je)
     print*, "pe=", mpp_pe(), is, ie, js, je

     allocate(global(0:ni+1,0:nj+1,nz) )

     global = 0
     do k = 1,nz
       do j = 1,nj
         do i = 1,ni
           global(i,j,k) = k + i*1e-3 + j*1e-6
         end do
       end do
     end do

     global(0,      1:nj,:) = global(ni,     1:nj,:)
     global(ni+1,   1:nj,:) = global(1,      1:nj,:)
     global(0:ni+1, 0,   :) = global(0:ni+1, nj,  :)
     global(0:ni+1, nj+1,:) = global(0:ni+1, 1,   :)

     ! set up x array
     call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
     call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
     allocate( x (isd:ied,jsd:jed,nz) )

     x(:,:,:) = 0.0
     x (is:ie,js:je,:) = global(is:ie,js:je,:)

     ! full update
     call mpp_update_domains(x, domain)
     call compare_checksums(x, global(isd:ied,jsd:jed,:), '64-bit array 9 pe subset')

     deallocate(x, global)
     call mpp_deallocate_domain(domain)
   endif

   call mpp_set_current_pelist()

  end subroutine test_subset_update_r8

 !> test a domain update of a 32-bit 3D array on a 9-pe subset of total allotted pes
 !> @note requires at least 16 pes
 subroutine test_subset_update_r4( )
   real(kind=r4_kind), allocatable, dimension(:,:,:) :: x
   type(domain2D) :: domain
   real(kind=r4_kind), allocatable :: global(:,:,:)
   integer              :: i, xhalo, yhalo
   integer              :: is, ie, js, je, isd, ied, jsd, jed
   integer :: pes9(9)=(/1,2,4,6,8,10,12,13,15/)
   integer :: ni, nj
   integer :: pe, npes

   pe = mpp_pe()
   npes = mpp_npes()

   call mpp_declare_pelist(pes9)
   if(any(mpp_pe()==pes9)) then
     call mpp_set_current_pelist(pes9)
     layout = (/3,3/)
     ni = 3; nj =3
     call mpp_define_domains((/1,ni,1,nj/), layout, domain, xhalo=1 &
       &, yhalo=1, xflags=CYCLIC_GLOBAL_DOMAIN, yflags&
       &=CYCLIC_GLOBAL_DOMAIN, name='subset domain')
     call mpp_get_compute_domain(domain, is, ie, js, je)
     print*, "pe=", mpp_pe(), is, ie, js, je

     allocate(global(0:ni+1,0:nj+1,nz) )

     global = 0
     do k = 1,nz
       do j = 1,nj
         do i = 1,ni
           global(i,j,k) = k + i*1e-3 + j*1e-6
         end do
       end do
     end do

     global(0,      1:nj,:) = global(ni,     1:nj,:)
     global(ni+1,   1:nj,:) = global(1,      1:nj,:)
     global(0:ni+1, 0,   :) = global(0:ni+1, nj,  :)
     global(0:ni+1, nj+1,:) = global(0:ni+1, 1,   :)

     ! set up x array
     call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
     call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
     allocate( x (isd:ied,jsd:jed,nz) )

     x(:,:,:) = 0.0
     x (is:ie,js:je,:) = global(is:ie,js:je,:)

     ! full update
     call mpp_update_domains(x, domain)
     call compare_checksums(x, global(isd:ied,jsd:jed,:), '32-bit array on 9 pe subset')

     deallocate(x, global)
     call mpp_deallocate_domain(domain)
   endif

   call mpp_set_current_pelist()

  end subroutine test_subset_update_r4

end module test_mpp_update_domains_real
