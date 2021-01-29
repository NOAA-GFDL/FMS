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
program test_mpp_global_field

  use platform_mod
  use compare_data_checksums
  use compare_data_checksums_int
  use mpp_mod,         only : mpp_init, mpp_error, FATAL, mpp_init_test_requests_allocated
  use mpp_mod,         only : mpp_declare_pelist, mpp_pe, mpp_npes, mpp_root_pe
  !use mpp_mod,        only : mpp_clock_begin, mpp_clock_end, mpp_clock_id, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED
  use mpp_domains_mod, only : domain2D
  use mpp_domains_mod, only : CENTER, EAST, NORTH, CORNER, XUPDATE, YUPDATE
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field

  implicit none

  integer, parameter :: nx=20, ny=20, nz=40
  integer, parameter :: whalo=2, ehalo=2, shalo=2, nhalo=2
  integer, parameter :: stackmax=4000000

  integer :: pe, npes, ierr
  integer :: layout(2)


  !> call mpp_init
  call mpp_init(test_level=mpp_init_test_requests_allocated)

  !> get pe info
  pe = mpp_pe()
  npes = mpp_npes()

  !> initialize mpp domain(s)
  call mpp_domains_init()
  call mpp_domains_set_stack_size(stackmax)

  !> call test_global_field_r4_2d
  call test_global_field_r4_2d( 'Non-symmetry' )
  call test_global_field_r4_2d( 'Symmetry center' )
  call test_global_field_r4_2d( 'Symmetry corner' )
  call test_global_field_r4_2d( 'Symmetry east' )
  call test_global_field_r4_2d( 'Symmetry north' )
  !> call test_global_field_r8_2d
  call test_global_field_r8_2d( 'Non-symmetry' )
  call test_global_field_r8_2d( 'Symmetry center' )
  call test_global_field_r8_2d( 'Symmetry corner' )
  call test_global_field_r8_2d( 'Symmetry east' )
  call test_global_field_r8_2d( 'Symmetry north' )
  !> call test_global_field_i4_2d
  call test_global_field_i4_2d( 'Non-symmetry' )
  call test_global_field_i4_2d( 'Symmetry center' )
  call test_global_field_i4_2d( 'Symmetry corner' )
  call test_global_field_i4_2d( 'Symmetry east' )
  call test_global_field_i4_2d( 'Symmetry north' )
  !> call test_global_field_i8_2d
  call test_global_field_i8_2d( 'Non-symmetry' )
  call test_global_field_i8_2d( 'Symmetry center' )
  call test_global_field_i8_2d( 'Symmetry corner' )
  call test_global_field_i8_2d( 'Symmetry east' )
  call test_global_field_i8_2d( 'Symmetry north' )
  !> call test_global_field_r4_3d tests
  call test_global_field_r4_3d( 'Non-symmetry' )
  call test_global_field_r4_3d( 'Symmetry center' )
  call test_global_field_r4_3d( 'Symmetry corner' )
  call test_global_field_r4_3d( 'Symmetry east' )
  call test_global_field_r4_3d( 'Symmetry north' )
  !> call test_global_field_r8_3d tests
  call test_global_field_r8_3d( 'Non-symmetry' )
  call test_global_field_r8_3d( 'Symmetry center' )
  call test_global_field_r8_3d( 'Symmetry corner' )
  call test_global_field_r8_3d( 'Symmetry east' )
  call test_global_field_r8_3d( 'Symmetry north' )
  !> call test_global_field_i4_3d tests
  call test_global_field_i4_3d( 'Non-symmetry' )
  call test_global_field_i4_3d( 'Symmetry center' )
  call test_global_field_i4_3d( 'Symmetry corner' )
  call test_global_field_i4_3d( 'Symmetry east' )
  call test_global_field_i4_3d( 'Symmetry north' )
  !> call test_global_field_i8_3d tests
  call test_global_field_i8_3d( 'Non-symmetry' )
  call test_global_field_i8_3d( 'Symmetry center' )
  call test_global_field_i8_3d( 'Symmetry corner' )
  call test_global_field_i8_3d( 'Symmetry east' )
  call test_global_field_i8_3d( 'Symmetry north' )

  !> exit
  call mpp_domains_exit()
  call MPI_finalize(ierr)

contains
  !>
  !> test_global_field_r4_2d
  !>
  subroutine test_global_field_r4_2d( type )

    implicit none

    character(len=*), intent(in) :: type

    real(kind=r4_kind), parameter :: zero = 0.0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    real(kind=r4_kind), allocatable  :: global1(:,:), x(:,:), gcheck(:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position=CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo) )
    global1 = zero
    do j=1, nj
       do i=1, ni
          global1(i,j) = real( i*1e-3+j*1e-6, kind=r4_kind )
       end do
    enddo

    allocate( gcheck(ni,nj) )

    !> allocate for global domain
    allocate( x(isd:ied,jsd:jed) )
    x(:,:) = global1(isd:ied,jsd:jed)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on r4 data domain' )


    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on r4 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on r4 data domain' )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on r4 data domain' )

    !> test the data on compute domain

    deallocate(x)
    allocate( x(is:ie,js:je) )
    x(is:ie,js:je) = global1(is:ie,js:je)

    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on r4 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on r4 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on r4 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_r4_2d
  !>
  !> test_global_field_r8_2d
  !>
  subroutine test_global_field_r8_2d( type )

    implicit none

    character(len=*), intent(in) :: type

    real(kind=r8_kind), parameter :: zero = 0.0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    real(kind=r8_kind), allocatable  :: global1(:,:), x(:,:), gcheck(:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position=CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo) )
    global1 = zero
    do j=1, nj
       do i=1, ni
          global1(i,j) = real( i*1e-3+j*1e-6, kind=r8_kind )
       end do
    enddo

    allocate( gcheck(ni,nj) )

    !> allocate for global domain
    allocate( x(isd:ied,jsd:jed) )
    x(:,:) = global1(isd:ied,jsd:jed)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on r8 data domain' )


    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on r8 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on r8 data domain' )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on r8 data domain' )

    !> test the data on compute domain

    deallocate(x)
    allocate( x(is:ie,js:je) )
    x(is:ie,js:je) = global1(is:ie,js:je)

    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on r8 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on r8 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on r8 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_r8_2d
  !>
  !> test_global_field_i4_2d
  !>
  subroutine test_global_field_i4_2d( type )

    implicit none

    character(len=*), intent(in) :: type

    integer(kind=i4_kind), parameter :: zero = 0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    integer(kind=i4_kind), allocatable  :: global1(:,:), x(:,:), gcheck(:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position=CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo) )
    global1 = zero
    do j=1, nj
       do i=1, ni
          global1(i,j) = int( i*1e3+j*1e6, kind=i4_kind )
       end do
    enddo

    allocate( gcheck(ni,nj) )

    !> allocate for global domain
    allocate( x(isd:ied,jsd:jed) )
    x(:,:) = global1(isd:ied,jsd:jed)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on i4 data domain' )


    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on i4 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on i4 data domain' )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on i4 data domain' )

    !> test the data on compute domain

    deallocate(x)
    allocate( x(is:ie,js:je) )
    x(is:ie,js:je) = global1(is:ie,js:je)

    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on i4 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on i4 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on i4 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_i4_2d
  !>
  !> test_global_field_i8_2d
  !>
  subroutine test_global_field_i8_2d( type )

    implicit none

    character(len=*), intent(in) :: type

    integer(kind=i8_kind), parameter :: zero = 0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    integer(kind=i8_kind), allocatable  :: global1(:,:), x(:,:), gcheck(:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position=CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo) )
    global1 = zero
    do j=1, nj
       do i=1, ni
          global1(i,j) = int( i*1e3+j*1e6, kind=i8_kind )
       end do
    enddo

    allocate( gcheck(ni,nj) )

    !> allocate for global domain
    allocate( x(isd:ied,jsd:jed) )
    x(:,:) = global1(isd:ied,jsd:jed)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on i8 data domain' )


    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on i8 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on i8 data domain' )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on i8 data domain' )

    !> test the data on compute domain

    deallocate(x)
    allocate( x(is:ie,js:je) )
    x(is:ie,js:je) = global1(is:ie,js:je)

    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj), gcheck, type//' mpp_global_field on i8 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je), gcheck(1:ni,js:je), type//' mpp_global_field xupdate only on i8 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj), gcheck(is:ie,1:nj), type//' mpp_global_field yupdate only on i8 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_i8_2d
  !>
  !> test_global_field_r4_3d
  !>
  subroutine test_global_field_r4_3d( type )

    implicit none

    character(len=*), intent(in) :: type

    real(kind=r4_kind)   :: zero = 0.0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    real(kind=r4_kind), allocatable  :: global1(:,:,:), x(:,:,:), gcheck(:,:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position = CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global1
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo,nz) )
    global1 = zero
    do k=1, nz
       do j=1, nj
          do i=1, ni
             global1(i,j,k) = real( k+i*1e-3+j*1e-6, kind=r4_kind )
          end do
       end do
    enddo

    allocate( gcheck(ni,nj,nz) )

    !> for data domain
    allocate( x(isd:ied,jsd:jed, nz) )
    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on r4 data domain' )

    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:),type//' mpp_global_field xupdate only on r4 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:),type//' mpp_global_field yupdate only on r4 data domain' )

    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck,type//' mpp_global_field on r4 data domain' )

    !> test the data on compute domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on r4 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
         type//' mpp_global_field xupdate only on r4 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
         type//' mpp_global_field yupdate only on r4 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_r4_3d
  !>
  !> test_global_field_r8_3d
  !>
  subroutine test_global_field_r8_3d( type )

    implicit none

    character(len=*), intent(in) :: type

    real(kind=r8_kind) :: zero = 0.0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    real(kind=r8_kind), allocatable  :: global1(:,:,:), x(:,:,:), gcheck(:,:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position = CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global1
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo,nz) )
    global1 = zero
    do k=1, nz
       do j=1, nj
          do i=1, ni
             global1(i,j,k) = real( k+i*1e-3+j*1e-6, kind=r8_kind )
          end do
       end do
    enddo

    allocate( gcheck(ni,nj,nz) )

    !> for data domain
    allocate( x(isd:ied,jsd:jed, nz) )
    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on r8 data domain' )

    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:),type//' mpp_global_field xupdate only on r8 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:),type//' mpp_global_field yupdate only on r8 data domain' )

    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck,type//' mpp_global_field on r8 data domain' )

    !> test the data on compute domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on r8 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
         type//' mpp_global_field xupdate only on r8 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
         type//' mpp_global_field yupdate only on r8 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_r8_3d
  !>
  !> test_global_field_i4_3d
  !>
  subroutine test_global_field_i4_3d( type )

    implicit none

    character(len=*), intent(in) :: type

    integer(kind=i4_kind)   :: zero = 0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    integer(kind=i4_kind), allocatable  :: global1(:,:,:), x(:,:,:), gcheck(:,:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position = CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global1
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo,nz) )
    global1 = zero
    do k=1, nz
       do j=1, nj
          do i=1, ni
             global1(i,j,k) = int( k+i*1e3+j*1e6, kind=i4_kind )
          end do
       end do
    enddo

    allocate( gcheck(ni,nj,nz) )

    !> for data domain
    allocate( x(isd:ied,jsd:jed, nz) )
    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on i4 data domain' )

    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:),type//' mpp_global_field xupdate only on i4 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:),type//' mpp_global_field yupdate only on i4 data domain' )

    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj,:), gcheck,type//' mpp_global_field on i4 data domain' )

    !> test the data on compute domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on i4 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
         type//' mpp_global_field xupdate only on i4 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
         type//' mpp_global_field yupdate only on i4 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_i4_3d
  !>
  !> test_global_field_i8_3d
  !>
  subroutine test_global_field_i8_3d( type )

    implicit none

    character(len=*), intent(in) :: type

    integer(kind=i8_kind)   :: zero = 0

    type(domain2D)  :: domain
    integer         :: position, ishift, jshift, ni, nj, i, j, k
    integer         :: is, ie, js, je, isd, ied, jsd, jed
    !integer        :: id
    integer, allocatable :: pelist(:)
    integer(kind=i8_kind), allocatable  :: global1(:,:,:), x(:,:,:), gcheck(:,:,:)


    !> set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Non-symmetry' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type )
    case( 'Symmetry center', 'Symmetry corner', 'Symmetry east', 'Symmetry north' )
       call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
            shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case default
       call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    !> get compute domain
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    !> get data domain
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !> determine if an extra point is needed
    ishift = 0 ; jshift = 0 ; position = CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1 ; jshift = 1 ; position=CORNER
    case ('Symmetry east')
       ishift = 1 ; jshift = 0 ; position=EAST
    case ('Symmetry north')
       ishift = 0 ; jshift = 1 ; position=NORTH
    end select

    ie  = ie+ishift  ; je  = je+jshift
    ied = ied+ishift ; jed = jed+jshift
    ni  = nx+ishift  ; nj  = ny+jshift

    !> assign global1
    allocate( global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo,nz) )
    global1 = zero
    do k=1, nz
       do j=1, nj
          do i=1, ni
             global1(i,j,k) = int( k+i*1e3+j*1e6, kind=i8_kind )
          end do
       end do
    enddo

    allocate( gcheck(ni,nj,nz) )

    !> for data domain
    allocate( x(isd:ied,jsd:jed, nz) )
    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !> test the data on data domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on i8 data domain' )

    !> Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !! will be declared. But for the x-direction global field, mpp_sync_self will
    !! be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !! in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !! deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !! will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !! on all pe is needed for those partial pelist. But for y-update, it is ok.
    !! because the pelist in y-update is not continous.
    allocate( pelist(0:layout(1)-1) )
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:),type//' mpp_global_field xupdate only on i8 data domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:),type//' mpp_global_field yupdate only on i8 data domain' )

    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj,:), gcheck,type//' mpp_global_field on i8 data domain' )

    !> test the data on compute domain
    gcheck = zero
    !id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, position=position )
    !call mpp_clock_end(id)
    !>compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on i8 compute domain' )

    !> xupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=XUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
         type//' mpp_global_field xupdate only on i8 compute domain' )

    !> yupdate
    gcheck = zero
    !call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie,js:je,:), gcheck, flags=YUPDATE, position=position )
    !call mpp_clock_end(id)
    !> compare checksums between global and x arrays
    call compare_checksums_int( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
         type//' mpp_global_field yupdate only on i8 compute domain' )

    deallocate(global1, gcheck, x)

  end subroutine test_global_field_i8_3d

end program test_mpp_global_field
