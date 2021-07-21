!***********************************************************************
!                   GNU Lesser General Public License
!
! This file is part of the GFDL Flexible Modeling System (FMS).
!
! FMS is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! FMS is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
! for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

!> @author Ryan Mulhall
!> @email gfdl.climate.model.info@noaa.gov
!> @description Unit test for mpp_redistribute with both sized integers
!> also tests redistribute with mosaic cubic grid
program test_mpp_redistribute

  use mpp_mod,         only : FATAL, WARNING, NOTE, MPP_INIT_TEST_INIT_TRUE_ONLY
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_error
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,         only : mpp_init, mpp_exit, stdout, stderr
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D, mpp_global_field
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_deallocate_domain
  use mpp_domains_mod, only : mpp_define_mosaic, mpp_nullify_domain_list
  use platform_mod

  implicit none
  integer :: pe, npes, ierr
  integer :: nx=128, ny=128, nz=40, stackmax=4000000

  call mpp_init(MPP_INIT_TEST_INIT_TRUE_ONLY)
  pe = mpp_pe()
  npes = mpp_npes()
  call mpp_domains_init()
  call mpp_domains_set_stack_size(stackmax)

  call mpp_error(NOTE, "----------Starting tests----------")
  call test_redistribute_i4()
  call mpp_error(NOTE, "test_mpp_redistribute: 32-bit integer test passed")
  call test_redistribute_i8()
  call mpp_error(NOTE, "test_mpp_redistribute: 64-bit integer test passed")
  call cubic_grid_redistribute_i4()
  call mpp_error(NOTE, "test_mpp_redistribute: 32-bit integer cubic grid test passed")
  call cubic_grid_redistribute_i8()
  call mpp_error(NOTE, "test_mpp_redistribute: 64-bit integer cubic grid test passed")
  call mpp_error(NOTE, "----------Tests Complete----------")

  call mpp_domains_exit()
  call mpi_finalize(ierr)

contains

  !> redistribute x domain to y with 32-bit integers
  subroutine test_redistribute_i4()
    type(domain2D) :: domainx, domainy
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    integer(i4_kind), allocatable, dimension(:,:,:)       :: gcheck, glbl
    integer(i4_kind), allocatable, dimension(:,:,:) :: x, x1, x2, x3, x4, x5, x6
    integer(i4_kind), allocatable, dimension(:,:,:) :: y, y1, y2, y3, y4, y5, y6
    integer                                      :: k, j, i, layout(2), id
    integer                                      :: is, ie, js, je, isd, ied, jsd, jed
    ! nullify domain list otherwise it retains memory between calls.
    call mpp_nullify_domain_list(domainx)
    call mpp_nullify_domain_list(domainy)

    !fill in glbl array with kiiijjj
    allocate( gcheck(nx,ny,nz), glbl(nx,ny,nz) )
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          glbl(i,j,k) = k*1e6 + i*1e3 + j
        end do
      end do
    end do

    ! set up x arrays
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domainx )
    call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
    allocate( x(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    allocate( x5(isd:ied,jsd:jed,nz) )
    allocate( x6(isd:ied,jsd:jed,nz) )
    x = 0.
    x(is:ie,js:je,:) = glbl(is:ie,js:je,:)
    x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x

    !set up y arrays
    call mpp_define_domains( (/1,nx,1,ny/), (/npes,1/), domainy)
    call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
    call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
    allocate( y(isd:ied,jsd:jed,nz) )
    allocate( y2(isd:ied,jsd:jed,nz) )
    allocate( y3(isd:ied,jsd:jed,nz) )
    allocate( y4(isd:ied,jsd:jed,nz) )
    allocate( y5(isd:ied,jsd:jed,nz) )
    allocate( y6(isd:ied,jsd:jed,nz) )
    y = 0.
    y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.

    !go global and redistribute
    call mpp_broadcast_domain(domainx)
    call mpp_broadcast_domain(domainy)
    call mpp_redistribute( domainx, x, domainy, y )

    !check answers on pelist
    call mpp_global_field( domainy, y(:,:,:), gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL , &
             "test_mpp_redistribute: incorrect results in global array")
    ! redistribute and check x answers
    if(ALLOCATED(y))y=0.
    call mpp_redistribute( domainx, x,  domainy, y,  complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_global_field( domainx, x(:,:,:), gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x")
    call mpp_global_field( domainx, x2, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x2")
    call mpp_global_field( domainx, x3, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x3")
    call mpp_global_field( domainx, x4, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x4")
    call mpp_global_field( domainx, x5, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
             "test_mpp_redistribute: global array differs for x5")
    call mpp_global_field( domainx, x6, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x6")

    ! redistribute and check y answers
    if(ALLOCATED(y))then
      y=0.; y2=0.; y3=0.; y4=0.; y5=0.; y6=0.
    endif
    call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_global_field( domainy, y, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y")
    call mpp_global_field( domainy, y2, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y2")
    call mpp_global_field( domainy, y3, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y3")
    call mpp_global_field( domainy, y4, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y4")
    call mpp_global_field( domainy, y5, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y5")
    call mpp_global_field( domainy, y6, gcheck )
    if(.not. compare_result4( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y6")

    dch =>NULL()

    call mpp_redistribute( domainx, x, domainy, y, free=.true.,list_size=6 )
    deallocate(gcheck, glbl)
    deallocate(x,x2,x3,x4,x5,x6)
    deallocate(y,y2,y3,y4,y5,y6)

  end subroutine test_redistribute_i4

  !> Test redistribute between two domains with 64-bit integers
  subroutine test_redistribute_i8()
    type(domain2D) :: domainx, domainy
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    integer(i8_kind), allocatable, dimension(:,:,:)       :: gcheck, glbl
    integer(i8_kind), allocatable, dimension(:,:,:), save :: x, x1, x2, x3, x4, x5, x6
    integer(i8_kind), allocatable, dimension(:,:,:), save :: y, y1, y2, y3, y4, y5, y6
    integer                                      :: k, j, i, layout(2), id
    integer                                      :: is, ie, js, je, isd, ied, jsd, jed
    ! nullify domain list otherwise it retains memory between calls.
    call mpp_nullify_domain_list(domainx)
    call mpp_nullify_domain_list(domainy)

    !fill in glbl array with kiiijjj
    allocate( gcheck(nx,ny,nz), glbl(nx,ny,nz) )
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          glbl(i,j,k) = k*1e6 + i*1e3 + j
        end do
      end do
    end do

    ! set up x arrays
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domainx )
    call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
    allocate( x(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    allocate( x5(isd:ied,jsd:jed,nz) )
    allocate( x6(isd:ied,jsd:jed,nz) )
    x = 0.
    x(is:ie,js:je,:) = glbl(is:ie,js:je,:)
    x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x

    !set up y arrays
    call mpp_define_domains( (/1,nx,1,ny/), (/npes,1/), domainy)
    call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
    call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
    allocate( y(isd:ied,jsd:jed,nz) )
    allocate( y2(isd:ied,jsd:jed,nz) )
    allocate( y3(isd:ied,jsd:jed,nz) )
    allocate( y4(isd:ied,jsd:jed,nz) )
    allocate( y5(isd:ied,jsd:jed,nz) )
    allocate( y6(isd:ied,jsd:jed,nz) )
    y = 0.
    y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.

    !go global and redistribute
    call mpp_broadcast_domain(domainx)
    call mpp_broadcast_domain(domainy)
    call mpp_redistribute( domainx, x, domainy, y )

    !check answers on pelist
    call mpp_global_field( domainy, y(:,:,:), gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL , &
             "test_mpp_redistribute: incorrect results in global array")
    ! redistribute and check x answers
    if(ALLOCATED(y))y=0.
    call mpp_redistribute( domainx, x,  domainy, y,  complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_global_field( domainx, x(:,:,:), gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x")
    call mpp_global_field( domainx, x2, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x2")
    call mpp_global_field( domainx, x3, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x3")
    call mpp_global_field( domainx, x4, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x4")
    call mpp_global_field( domainx, x5, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
             "test_mpp_redistribute: global array differs for x5")
    call mpp_global_field( domainx, x6, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for x6")

    ! redistribute and check y answers
    if(ALLOCATED(y))then
      y=0.; y2=0.; y3=0.; y4=0.; y5=0.; y6=0.
    endif
    call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_global_field( domainy, y, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y")
    call mpp_global_field( domainy, y2, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y2")
    call mpp_global_field( domainy, y3, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y3")
    call mpp_global_field( domainy, y4, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y4")
    call mpp_global_field( domainy, y5, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y5")
    call mpp_global_field( domainy, y6, gcheck )
    if(.not. compare_result8( glbl(1:nx,1:ny,:), gcheck )) call mpp_error(FATAL,&
            "test_mpp_redistribute: global array differs for y6")

    dch =>NULL()

    call mpp_redistribute( domainx, x, domainy, y, free=.true.,list_size=6 )
    deallocate(gcheck, glbl)
    deallocate(x,x2,x3,x4,x5,x6)
    deallocate(y,y2,y3,y4,y5,y6)

  end subroutine test_redistribute_i8
  !> Tests redistribute with cubic grid and 32-bit ints
  subroutine cubic_grid_redistribute_i4

     integer              :: npes, npes_per_ensemble, npes_per_tile
     integer              :: ensemble_id, tile_id, ensemble_tile_id
     integer              :: i, j, p, n, ntiles, my_root_pe, k
     integer              :: isc_ens, iec_ens, jsc_ens, jec_ens
     integer              :: isd_ens, ied_ens, jsd_ens, jed_ens
     integer              :: isc, iec, jsc, jec
     integer              :: isd, ied, jsd, jed, layout_ensemble(2) = (/0,0/)
     integer, allocatable :: my_ensemble_pelist(:), pe_start(:), pe_end(:)
     integer, allocatable :: global_indices(:,:), layout2D(:,:)
     integer(i4_kind),    allocatable :: x(:,:,:,:), x_ens(:,:,:), y(:,:,:)
     integer              :: layout(2), ensemble_size = 1, layout_cubic(2) = (/0,0/)
     type(domain2D)       :: domain
     type(domain2D), allocatable :: domain_ensemble(:)
     character(len=128)   :: mesg
     integer              :: nx_cubic = 20, ny_cubic = 20

     logical :: check
     ! set up pelist
     npes = mpp_npes()
     if(mod(npes, ensemble_size) .NE. 0) call mpp_error(FATAL, &
         "test_mpp_domains: npes is not divisible by ensemble_size")
     npes_per_ensemble = npes/ensemble_size
     allocate(my_ensemble_pelist(0:npes_per_ensemble-1))
     ensemble_id = mpp_pe()/npes_per_ensemble + 1
     do p = 0, npes_per_ensemble-1
        my_ensemble_pelist(p) = (ensemble_id-1)*npes_per_ensemble + p
     enddo

     call mpp_declare_pelist(my_ensemble_pelist)

     ! set tile count
     ntiles = 6

     if( mod(npes, ntiles) .NE. 0 ) call mpp_error(FATAL, &
          "test_mpp_domains: npes is not divisible by ntiles")

     npes_per_tile = npes/ntiles
     tile_id = mpp_pe()/npes_per_tile + 1
     if( npes_per_tile == layout_cubic(1) * layout_cubic(2) ) then
        layout = layout_cubic
     else
        call mpp_define_layout( (/1,nx_cubic,1,ny_cubic/), npes_per_tile, layout )
     endif
     allocate(global_indices(4, ntiles))
     allocate(layout2D(2, ntiles))
     allocate(pe_start(ntiles), pe_end(ntiles))
     do n = 1, ntiles
       global_indices(:,n) = (/1,nx_cubic,1,ny_cubic/)
       layout2D(:,n)         = layout
     end do

     do n = 1, ntiles
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile-1
     end do

     call define_cubic_mosaic("cubic_grid", domain, (/nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic/), &
                              (/ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic/), &
                                global_indices, layout2D, pe_start, pe_end )

     allocate(domain_ensemble(ensemble_size))
     !-- define domain for each ensemble
     call mpp_set_current_pelist( my_ensemble_pelist )
     if( mod(npes_per_ensemble, ntiles) .NE. 0 ) call mpp_error(FATAL, &
          "test_mpp_domains: npes_per_ensemble is not divisible by ntiles")
     npes_per_tile = npes_per_ensemble/ntiles
     my_root_pe = my_ensemble_pelist(0)
     ensemble_tile_id = (mpp_pe() - my_root_pe)/npes_per_tile + 1

     if( npes_per_tile == layout_ensemble(1) * layout_ensemble(2) ) then
        layout = layout_ensemble
     else
        call mpp_define_layout( (/1,nx_cubic,1,ny_cubic/), npes_per_tile, layout )
     endif
     do n = 1, ntiles
       global_indices(:,n) = (/1,nx_cubic,1,ny_cubic/)
       layout2D(:,n)         = layout
     end do

     do n = 1, ntiles
        pe_start(n) = my_root_pe + (n-1)*npes_per_tile
        pe_end(n)   = my_root_pe + n*npes_per_tile-1
     end do

     call define_cubic_mosaic("cubic_grid",domain_ensemble(ensemble_id),(/nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic/)&
                             ,(/ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic/), &
                                global_indices, layout2D, pe_start, pe_end )

     call mpp_set_current_pelist()
     do n = 1, ensemble_size
        call mpp_broadcast_domain(domain_ensemble(n))
     enddo

     call mpp_get_data_domain( domain_ensemble(ensemble_id), isd_ens, ied_ens, jsd_ens, jed_ens)
     call mpp_get_compute_domain( domain_ensemble(ensemble_id), isc_ens, iec_ens, jsc_ens, jec_ens)
     call mpp_get_data_domain( domain, isd, ied, jsd, jed)
     call mpp_get_compute_domain( domain, isc, iec, jsc, jec)

     allocate(x_ens(isd_ens:ied_ens, jsd_ens:jed_ens, nz))
     allocate(x(isd:ied, jsd:jed, nz, ensemble_size))
     allocate(y(isd:ied, jsd:jed, nz))

     x = 0
     do k = 1, nz
        do j = jsc_ens, jec_ens
           do i = isc_ens, iec_ens
              x_ens(i,j,k) = ensemble_id *1e6 + ensemble_tile_id*1e3 + i + j * 1.e-3 + k * 1.e-6
           enddo
        enddo
     enddo

     do n = 1, ensemble_size
        x = 0
        call mpp_redistribute( domain_ensemble(n), x_ens, domain, x(:,:,:,n) )
        y = 0
        do k = 1, nz
           do j = jsc, jec
              do i = isc, iec
                 y(i,j,k) = n *1e6 + tile_id*1e3 + i + j * 1.e-3 + k * 1.e-6
              enddo
           enddo
        enddo
       write(mesg,'(a,i4)') "cubic_grid redistribute from ensemble", n
       if(.not.compare_result4( x(isc:iec,jsc:jec,:,n), y(isc:iec,jsc:jec,:))) call &
           mpp_error(FATAL, "test_mpp_redistribute: failed cubic grid 32-bit check")
     enddo

     ! redistribute data to each ensemble.
     deallocate(x,y,x_ens)
     allocate(x(isd:ied, jsd:jed, nz, ensemble_size))
     allocate(x_ens(isd_ens:ied_ens, jsd_ens:jed_ens, nz))
     allocate(y(isd_ens:ied_ens, jsd_ens:jed_ens, nz))

     y = 0
     do k = 1, nz
        do j = jsc, jec
           do i = isc, iec
              x(i,j,k,:) = i + j * 1.e-3 + k * 1.e-6
           enddo
        enddo
     enddo

     do n = 1, ensemble_size
        x_ens = 0
        call mpp_redistribute(domain, x(:,:,:,n), domain_ensemble(n), x_ens)
        y = 0
        if( ensemble_id == n ) then
           do k = 1, nz
              do j = jsc_ens, jec_ens
                 do i = isc_ens, iec_ens
                    y(i,j,k) = i + j * 1.e-3 + k * 1.e-6
                 enddo
              enddo
           enddo
        endif
        write(mesg,'(a,i4)') "cubic_grid redistribute to ensemble", n
        if(.not.compare_result4( x_ens(isc_ens:iec_ens,jsc_ens:jec_ens,:), y(isc_ens:iec_ens,jsc_ens:jec_ens,:))) call &
           mpp_error(FATAL, "test_mpp_redistribute: failed cubic grid 32-bit check")
     enddo

     deallocate(x, y, x_ens)
     call mpp_deallocate_domain(domain)
     do n = 1, ensemble_size
        call mpp_deallocate_domain(domain_ensemble(n))
     enddo
     deallocate(domain_ensemble)

  end subroutine cubic_grid_redistribute_i4

  !> Tests redistribute with cubic grid and 64-bit ints
  subroutine cubic_grid_redistribute_i8

     integer              :: npes, npes_per_ensemble, npes_per_tile
     integer              :: ensemble_id, tile_id, ensemble_tile_id
     integer              :: i, j, p, n, ntiles, my_root_pe, k
     integer              :: isc_ens, iec_ens, jsc_ens, jec_ens
     integer              :: isd_ens, ied_ens, jsd_ens, jed_ens
     integer              :: isc, iec, jsc, jec
     integer              :: isd, ied, jsd, jed, layout_ensemble(2) = (/0,0/)
     integer, allocatable :: my_ensemble_pelist(:), pe_start(:), pe_end(:)
     integer, allocatable :: global_indices(:,:), layout2D(:,:)
     integer(i8_kind),    allocatable :: x(:,:,:,:), x_ens(:,:,:), y(:,:,:)
     integer              :: layout(2), ensemble_size = 1, layout_cubic(2) = (/0,0/)
     type(domain2D)       :: domain
     type(domain2D), allocatable :: domain_ensemble(:)
     character(len=128)   :: mesg
     integer              :: nx_cubic = 20, ny_cubic = 20
     logical :: check

     ! --- set up pelist
     npes = mpp_npes()
     if(mod(npes, ensemble_size) .NE. 0) call mpp_error(FATAL, &
         "test_mpp_domains: npes is not divisible by ensemble_size")
     npes_per_ensemble = npes/ensemble_size
     allocate(my_ensemble_pelist(0:npes_per_ensemble-1))
     ensemble_id = mpp_pe()/npes_per_ensemble + 1
     do p = 0, npes_per_ensemble-1
        my_ensemble_pelist(p) = (ensemble_id-1)*npes_per_ensemble + p
     enddo

     call mpp_declare_pelist(my_ensemble_pelist)

     !--- define a mosaic use all the pelist
     ntiles = 6


     if( mod(npes, ntiles) .NE. 0 ) call mpp_error(FATAL, &
          "test_mpp_domains: npes is not divisible by ntiles")

     npes_per_tile = npes/ntiles
     tile_id = mpp_pe()/npes_per_tile + 1
     if( npes_per_tile == layout_cubic(1) * layout_cubic(2) ) then
        layout = layout_cubic
     else
        call mpp_define_layout( (/1,nx_cubic,1,ny_cubic/), npes_per_tile, layout )
     endif
     allocate(global_indices(4, ntiles))
     allocate(layout2D(2, ntiles))
     allocate(pe_start(ntiles), pe_end(ntiles))
     do n = 1, ntiles
       global_indices(:,n) = (/1,nx_cubic,1,ny_cubic/)
       layout2D(:,n)         = layout
     end do

     do n = 1, ntiles
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile-1
     end do

     call define_cubic_mosaic("cubic_grid", domain, (/nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic/), &
                              (/ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic/), &
                                global_indices, layout2D, pe_start, pe_end )

     allocate(domain_ensemble(ensemble_size))
     !-- define domain for each ensemble
     call mpp_set_current_pelist( my_ensemble_pelist )
     if( mod(npes_per_ensemble, ntiles) .NE. 0 ) call mpp_error(FATAL, &
          "test_mpp_domains: npes_per_ensemble is not divisible by ntiles")
     npes_per_tile = npes_per_ensemble/ntiles
     my_root_pe = my_ensemble_pelist(0)
     ensemble_tile_id = (mpp_pe() - my_root_pe)/npes_per_tile + 1

     if( npes_per_tile == layout_ensemble(1) * layout_ensemble(2) ) then
        layout = layout_ensemble
     else
        call mpp_define_layout( (/1,nx_cubic,1,ny_cubic/), npes_per_tile, layout )
     endif
     do n = 1, ntiles
       global_indices(:,n) = (/1,nx_cubic,1,ny_cubic/)
       layout2D(:,n)         = layout
     end do

     do n = 1, ntiles
        pe_start(n) = my_root_pe + (n-1)*npes_per_tile
        pe_end(n)   = my_root_pe + n*npes_per_tile-1
     end do

     call define_cubic_mosaic("cubic_grid",domain_ensemble(ensemble_id),(/nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic/)&
                             ,(/ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic/), &
                                global_indices, layout2D, pe_start, pe_end )

     call mpp_set_current_pelist()
     do n = 1, ensemble_size
        call mpp_broadcast_domain(domain_ensemble(n))
     enddo

     call mpp_get_data_domain( domain_ensemble(ensemble_id), isd_ens, ied_ens, jsd_ens, jed_ens)
     call mpp_get_compute_domain( domain_ensemble(ensemble_id), isc_ens, iec_ens, jsc_ens, jec_ens)
     call mpp_get_data_domain( domain, isd, ied, jsd, jed)
     call mpp_get_compute_domain( domain, isc, iec, jsc, jec)

     allocate(x_ens(isd_ens:ied_ens, jsd_ens:jed_ens, nz))
     allocate(x(isd:ied, jsd:jed, nz, ensemble_size))
     allocate(y(isd:ied, jsd:jed, nz))

     x = 0
     do k = 1, nz
        do j = jsc_ens, jec_ens
           do i = isc_ens, iec_ens
              x_ens(i,j,k) = ensemble_id *1e6 + ensemble_tile_id*1e3 + i + j * 1.e-3 + k * 1.e-6
           enddo
        enddo
     enddo

     do n = 1, ensemble_size
        x = 0
        call mpp_redistribute( domain_ensemble(n), x_ens, domain, x(:,:,:,n) )
        y = 0
        do k = 1, nz
           do j = jsc, jec
              do i = isc, iec
                 y(i,j,k) = n *1e6 + tile_id*1e3 + i + j * 1.e-3 + k * 1.e-6
              enddo
           enddo
        enddo
       write(mesg,'(a,i4)') "cubic_grid redistribute from ensemble", n
       if(.not.compare_result8( x(isc:iec,jsc:jec,:,n), y(isc:iec,jsc:jec,:))) call &
           mpp_error(FATAL, "test_mpp_redistribute: failed cubic grid 32-bit check")
     enddo

     ! redistribute data to each ensemble.
     deallocate(x,y,x_ens)
     allocate(x(isd:ied, jsd:jed, nz, ensemble_size))
     allocate(x_ens(isd_ens:ied_ens, jsd_ens:jed_ens, nz))
     allocate(y(isd_ens:ied_ens, jsd_ens:jed_ens, nz))

     y = 0
     do k = 1, nz
        do j = jsc, jec
           do i = isc, iec
              x(i,j,k,:) = i + j * 1.e-3 + k * 1.e-6
           enddo
        enddo
     enddo

     do n = 1, ensemble_size
        x_ens = 0
        call mpp_redistribute(domain, x(:,:,:,n), domain_ensemble(n), x_ens)
        y = 0
        if( ensemble_id == n ) then
           do k = 1, nz
              do j = jsc_ens, jec_ens
                 do i = isc_ens, iec_ens
                    y(i,j,k) = i + j * 1.e-3 + k * 1.e-6
                 enddo
              enddo
           enddo
        endif
        write(mesg,'(a,i4)') "cubic_grid redistribute to ensemble", n
        if(.not.compare_result8( x_ens(isc_ens:iec_ens,jsc_ens:jec_ens,:), y(isc_ens:iec_ens,jsc_ens:jec_ens,:))) call &
           mpp_error(FATAL, "test_mpp_redistribute: failed cubic grid 32-bit check")
     enddo

     deallocate(x, y, x_ens)
     call mpp_deallocate_domain(domain)
     do n = 1, ensemble_size
        call mpp_deallocate_domain(domain_ensemble(n))
     enddo
     deallocate(domain_ensemble)

  end subroutine cubic_grid_redistribute_i8

  ! define mosaic domain for cubic grid
  subroutine define_cubic_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end)
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)
    integer                       :: nhalo=1, shalo=1, ehalo=1, whalo=1

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

    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2, &
      istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
      pe_start, pe_end, symmetry = .true., whalo=whalo, ehalo=ehalo,   &
      shalo=shalo, nhalo=nhalo, name = trim(type), memory_size = msize  )

  end subroutine define_cubic_mosaic

  !> checks if i4 arrays are equal
  function compare_result4(a, b)
    integer(i4_kind), intent(in), dimension(:,:,:)  :: a, b
    logical                                         :: compare_result4
    if(size(a,1).ne.size(b,1) .or. size(a,2).ne.size(b,2) .or. size(a,3).ne.size(b,3)) call &
                        mpp_error(FATAL, "test_mpp_redistribute: comparing different sized arrays")
    compare_result4 = all(a.eq.b)
  end function compare_result4

  !> checks if i8 arrays are equal
  function compare_result8(a, b)
    integer(i8_kind), intent(in), dimension(:,:,:)  :: a, b
    logical                                         :: compare_result8
    if(size(a,1).ne.size(b,1) .or. size(a,2).ne.size(b,2) .or. size(a,3).ne.size(b,3)) call &
                        mpp_error(FATAL, "test_mpp_redistribute: comparing different sized arrays")
    compare_result8 = all(a.eq.b)
  end function compare_result8
end program test_mpp_redistribute