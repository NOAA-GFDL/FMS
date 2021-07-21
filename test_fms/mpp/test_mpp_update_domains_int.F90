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
!> @brief Test mpp_update_domains on arrays of integers using different layouts and data precision
!!
!> @author Jessica Liptak
!!
!> @note This test is an extension of the routine test_halo_upate in test_mpp_domains.

module test_mpp_update_domains_int

!  use compare_data_checksums, only : compare_checksums
  use fill_halo
  use compare_data_checksums_int, only : compare_checksums => compare_checksums_int
  use mpp_mod, only : FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod, only : mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_sync_self
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end
  use mpp_mod, only : mpp_declare_pelist, mpp_set_current_pelist, mpp_set_stack_size
  use mpp_mod, only : mpp_broadcast, mpp_transmit, mpp_sum, mpp_max, mpp_chksum, ALL_PES
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
  use mpp_io_mod, only: mpp_io_init
  use platform_mod, only: i4_kind, i8_kind

  implicit none
  private
  integer :: id
  integer :: nx=64, ny=64, nz=10, stackmax=10000000
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
  integer :: n
  integer :: stdunit = 6

  public :: test_halo_update_i8, test_halo_update_i4, test_subset_update_i8, test_subset_update_i4

  contains

  !> Perform simple addition on arrays of 64-bit integersin different domain configurations and update the domains
  subroutine test_halo_update_i8( domain_type )
    character(len=*), intent(in) :: domain_type !< the domain type that will be tested
    integer(kind=i8_kind), allocatable, dimension(:,:,:) :: xi8, x1i8, x2i8, x3i8, x4i8
    type(domain2D) :: domain
    integer(kind=i8_kind),    allocatable :: globali8(:,:,:)
    integer              :: shift, xhalo, yhalo
    logical              :: is_symmetry
    integer              :: is, ie, js, je, isd, ied, jsd, jed
    integer(kind=i8_kind) :: i, j, k ! kind specified because i,j,k define the x#i8 and globali8 values
    integer :: pe, npes

    pe = mpp_pe()
    npes = mpp_npes()

    allocate(globali8(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )
    ! populate the global array
    globali8 = 0.0
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          globali8(i,j,k) = k + i*1e-3 + j*1e-6
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
        globali8(1-whalo:0,                 1:ny,:) = globali8(nx-whalo+1:nx,             1:ny,:)
        globali8(nx+1:nx+ehalo,             1:ny,:) = globali8(1:ehalo,                   1:ny,:)
        globali8(1-whalo:nx+ehalo,     1-shalo:0,:) = globali8(1-whalo:nx+ehalo, ny-shalo+1:ny,:)
        globali8(1-whalo:nx+ehalo, ny+1:ny+nhalo,:) = globali8(1-whalo:nx+ehalo,       1:nhalo,:)
      case default
        call mpp_error( FATAL, 'test_halo_update_i8: no such test: '//domain_type )
    end select

    ! define the arrays
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain( domain, isd, ied, jsd, jed )
    allocate(xi8(isd:ied,jsd:jed,nz) )
    allocate(x1i8(isd:ied,jsd:jed,nz) )
    allocate(x2i8(isd:ied,jsd:jed,nz) )
    allocate(x3i8(isd:ied,jsd:jed,nz) )
    allocate(x4i8(isd:ied,jsd:jed,nz) )
    xi8(:,:,:) = 0.0
    xi8 (is:ie,js:je,:) = globali8(is:ie,js:je,:)
    x1i8 = xi8; x2i8 = xi8; x3i8 = xi8; x4i8 = xi8

    ! update the halo region
    id = mpp_clock_id( domain_type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( xi8, domain )
    call mpp_clock_end  (id)
    call compare_checksums( xi8, globali8(isd:ied,jsd:jed,:), domain_type )

    ! update part of the halo region
    id = mpp_clock_id( domain_type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1i8, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2i8, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3i8, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4i8, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x1i8(is:ied,js:jed,:), globali8(is:ied,js:jed,:), domain_type//' partial x1i8' )
    call compare_checksums( x2i8(is:ied,js:jed,:), globali8(is:ied,js:jed,:), domain_type//' partial x2i8' )
    call compare_checksums( x3i8(is:ied,js:jed,:), globali8(is:ied,js:jed,:), domain_type//' partial x3i8' )
    call compare_checksums( x4i8(is:ied,js:jed,:), globali8(is:ied,js:jed,:), domain_type//' partial x4i8' )

    deallocate(globali8, xi8, x1i8, x2i8, x3i8, x4i8)

  end subroutine test_halo_update_i8

  !> Perform simple addition on 32-bit real arrays in different domain configurations and update the domains
  subroutine test_halo_update_i4( domain_type )
   character(len=*), intent(in) :: domain_type !< the domain type that will be tested
   integer(kind=i4_kind), allocatable, dimension(:,:,:) :: xi4, x1i4, x2i4, x3i4, x4i4
   type(domain2D) :: domain
   integer(kind=i4_kind),    allocatable :: globali4(:,:,:)
   integer              :: shift, xhalo, yhalo
   logical              :: is_symmetry
   integer              :: is, ie, js, je, isd, ied, jsd, jed
   integer(kind=i4_kind) :: i, j, k ! kind specified because i,j,k define the x#i4 and globali4 values
   integer :: pe, npes

   pe = mpp_pe()
   npes = mpp_npes()

   allocate(globali4(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )

   globali4 = 0
   do k = 1,nz
     do j = 1,ny
       do i = 1,nx
         globali4(i,j,k) = k + i*1e-3 + j*1e-6
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
     globali4(1-whalo:0,                 1:ny,:) = globali4(nx-whalo+1:nx,             1:ny,:)
     globali4(nx+1:nx+ehalo,             1:ny,:) = globali4(1:ehalo,                   1:ny,:)
     globali4(1-whalo:nx+ehalo,     1-shalo:0,:) = globali4(1-whalo:nx+ehalo, ny-shalo+1:ny,:)
     globali4(1-whalo:nx+ehalo, ny+1:ny+nhalo,:) = globali4(1-whalo:nx+ehalo,       1:nhalo,:)
   case default
     call mpp_error( FATAL, 'test_halo_update_i4: '//domain_type//' is not a valid test.')
   end select
   ! define the arrays
   call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
   call mpp_get_data_domain( domain, isd, ied, jsd, jed )
   allocate(xi4(isd:ied,jsd:jed,nz) )
   allocate(x1i4(isd:ied,jsd:jed,nz) )
   allocate(x2i4(isd:ied,jsd:jed,nz) )
   allocate(x3i4(isd:ied,jsd:jed,nz) )
   allocate(x4i4(isd:ied,jsd:jed,nz) )
   xi4 = 0.0
   xi4 (is:ie,js:je,:) = globali4(is:ie,js:je,:)
   x1i4 = xi4; x2i4 = xi4; x3i4 = xi4; x4i4 = xi4
   ! update the halo region
   id = mpp_clock_id( domain_type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
   call mpp_clock_begin(id)
   call mpp_update_domains( xi4, domain )
   call mpp_clock_end  (id)
   call compare_checksums( xi4, globali4(isd:ied,jsd:jed,:), domain_type )
   ! update part of the halo region
   id = mpp_clock_id( domain_type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
   call mpp_clock_begin(id)
   call mpp_update_domains( x1i4, domain, NUPDATE+EUPDATE, complete=.false. )
   call mpp_update_domains( x2i4, domain, NUPDATE+EUPDATE, complete=.false. )
   call mpp_update_domains( x3i4, domain, NUPDATE+EUPDATE, complete=.false. )
   call mpp_update_domains( x4i4, domain, NUPDATE+EUPDATE, complete=.true. )
   call mpp_clock_end  (id)
   call compare_checksums( x1i4(is:ied,js:jed,:), globali4(is:ied,js:jed,:), domain_type//' partial x1i4' )
   call compare_checksums( x2i4(is:ied,js:jed,:), globali4(is:ied,js:jed,:), domain_type//' partial x2i4' )
   call compare_checksums( x3i4(is:ied,js:jed,:), globali4(is:ied,js:jed,:), domain_type//' partial x3i4' )
   call compare_checksums( x4i4(is:ied,js:jed,:), globali4(is:ied,js:jed,:), domain_type//' partial x4i4' )

   deallocate(globali4, xi4, x1i4, x2i4, x3i4, x4i4)

 end subroutine test_halo_update_i4

 !> test a domain update of a 3D array of 64-bit integers on a 9-pe subset of total allotted pes
 !> @note requires at least 16 pes
 subroutine test_subset_update_i8( )
   integer(kind=i8_kind), allocatable, dimension(:,:,:) :: x
   type(domain2D) :: domain
   integer(kind=i8_kind), allocatable :: global(:,:,:)
   integer              :: xhalo, yhalo
   integer              :: is, ie, js, je, isd, ied, jsd, jed
   integer :: pes9(9)=(/1,2,4,6,8,10,12,13,15/)
   integer :: ni, nj
   integer(kind=i8_kind) :: i, j, k ! kind specified because i,j,k define the x and global values
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

  end subroutine test_subset_update_i8

 !> test a domain update of a 3D array of 32-bit integers on a 9-pe subset of total allotted pes
 !> @note requires at least 16 pes
 subroutine test_subset_update_i4( )
   integer(kind=i4_kind), allocatable, dimension(:,:,:) :: x
   type(domain2D) :: domain
   integer(kind=i4_kind), allocatable :: global(:,:,:)
   integer              :: xhalo, yhalo
   integer              :: is, ie, js, je, isd, ied, jsd, jed
   integer :: pes9(9)=(/1,2,4,6,8,10,12,13,15/)
   integer :: ni, nj
   integer(kind=i4_kind) :: i, j, k ! kind specified because i,j,k define the xand globalvalues
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

  end subroutine test_subset_update_i4

end module test_mpp_update_domains_int
