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
!> @brief Run performance tests using blocking communications with mpp_update_domains,
!! and non-blocking communications with mpp_start_update_domains and mpp_complete_update_domains.
!!
!! The test cases are 'folded_north' and 'cubic_grid'
program test_update_domains_performance
  use compare_data_checksums, only : compare_checksums
  use compare_data_checksums_int, only : compare_checksums_int
  use mpp_mod, only : FATAL, WARNING, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod, only : mpp_init, mpp_pe, mpp_npes, mpp_root_pe, mpp_error
  use mpp_mod, only : mpp_clock_id, mpp_clock_begin, mpp_clock_end, mpp_sync
  use mpp_mod, only : mpp_init_test_requests_allocated
  use mpp_domains_mod, only : BGRID_NE, CGRID_NE, AGRID, SCALAR_PAIR, MPP_DOMAIN_TIME
  use mpp_domains_mod, only : domain2D
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit
  use mpp_domains_mod, only : mpp_update_domains, mpp_get_memory_domain
  use mpp_domains_mod, only : mpp_define_layout
  use mpp_domains_mod, only : mpp_define_mosaic
  use mpp_domains_mod, only : NORTH, SOUTH, WEST, EAST, CENTER
  use mpp_domains_mod, only : mpp_get_global_domain, ZERO
  use mpp_domains_mod, only : mpp_start_update_domains, mpp_complete_update_domains
  use mpp_io_mod, only: mpp_io_init
  use platform_mod

  implicit none

  integer :: ierr, id
  integer :: pe, npes
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
  integer :: num_fields = 4
  integer :: ensemble_size = 1
  integer :: num_iter = 1
  integer :: layout_cubic(2) = (/0,0/)
  integer :: layout_tripolar(2) = (/0,0/)
  integer :: layout_ensemble(2) = (/0,0/)
  logical :: do_sleep = .false.
  logical :: mix_2D_3D = .false.
  !> Initialize mpp and mpp IO modules
  call mpp_init(test_level=mpp_init_test_requests_allocated)
  call mpp_io_init()
  call mpp_domains_init(MPP_DOMAIN_TIME)
  call mpp_domains_set_stack_size(stackmax)
  pe = mpp_pe()
  npes = mpp_npes()
  !> run the tests
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Calling 64-bit real update_domains_performance tests <-------------------'
  call update_domains_performance_r8('Folded-north')
  call update_domains_performance_r8('Cubic-Grid')
  call update_domains_performance_r8('Single-Tile')
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Finished 64-bit real update_domains_performance tests <-------------------'
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Calling 32-bit real update_domains_performance tests <-------------------'
  call update_domains_performance_r4('Folded-north')
  call update_domains_performance_r4('Cubic-Grid')
  call update_domains_performance_r4('Single-Tile')
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Finished 32-bit real update_domains_performance tests'
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Calling 64-bit integer update_domains_performance tests <-------------------'
  call update_domains_performance_i8('Folded-north')
  call update_domains_performance_i8('Cubic-Grid')
  call update_domains_performance_i8('Single-Tile')
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Finished 64-bit integer update_domains_performance tests <-------------------'
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Finished 64-bit integer update_domains_performance tests'
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Calling 32-bit integer update_domains_performance tests <-------------------'
  call update_domains_performance_i4('Folded-north')
  call update_domains_performance_i4('Cubic-Grid')
  call update_domains_performance_i4('Single-Tile')
  if (mpp_pe() == mpp_root_pe()) &
    print *, '--------------------> Finished 32-bit integer update_domains_performance tests <-------------------'
  call mpp_domains_exit()
  !> Finalize mpp
  call MPI_FINALIZE(ierr)
  contains

  !> run performance tests on 64-bit real arrays
  subroutine update_domains_performance_r8( test_type )
    character(len=*), intent(in) :: test_type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile, ntile_per_pe
    integer        :: i, j, k, l, n, shift
    integer        :: ism, iem, jsm, jem
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real(kind=r8_kind),    allocatable, dimension(:,:,:,:) :: x, x1, y, y1, x_save, y_save
    real(kind=r8_kind),    allocatable, dimension(:,:,:,:) :: a, a1, b, b1
    real(kind=r8_kind),     allocatable, dimension(:,:,:  ) :: a1_2D, b1_2D
    integer            :: id_update
    integer            :: id1, id2
    logical            :: folded_north
    logical            :: cubic_grid, single_tile
    character(len=3)   :: text
    integer            :: nx_save, ny_save
    integer            :: id_single, id_update_single

    folded_north       = .false.
    cubic_grid         = .false.
    single_tile        = .false.
    nx_save = nx
    ny_save = ny
    !--- check the test_type
    select case(test_type)
      case ( 'Single-Tile' )   !--- single with cyclic along x- and y-direction
        single_tile = .true.
        ntiles = 1
        num_contact = 2
      case ( 'Folded-north' )
        ntiles = 1
        num_contact = 2
        folded_north = .true.
      case ( 'Cubic-Grid' )
        if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'update_domains_performance_r8: for Cubic_grid mosaic, nx_cubic is zero, '//&
                        'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'update_domains_performance_r8: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                        'No test is done for Cubic-Grid mosaic. ' )
          return
        endif

        nx = nx_cubic
        ny = ny_cubic
        ntiles = 6
        num_contact = 12
        cubic_grid = .true.

      case default
        call mpp_error(FATAL, 'update_domains_performancez_r8: no such test: '//test_type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    if( mod(npes, ntiles) == 0 ) then
      npes_per_tile = npes/ntiles
      text=""
      write(text,"(I3)") npes_per_tile
      call mpp_error(NOTE, 'update_domains_performance_r8: For Mosaic "'//trim(test_type)// &
                    '", each tile will be distributed over '//text//' processors.')
      ntile_per_pe = 1
      allocate(tile(ntile_per_pe))
      tile = pe/npes_per_tile+1
      if(cubic_grid) then
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      else
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      endif
      do n = 1, ntiles
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile-1
      end do
    else if ( mod(ntiles, npes) == 0 ) then
      ntile_per_pe = ntiles/npes
      text=""
      write(text,"(I3)") ntile_per_pe
      call mpp_error(NOTE, 'update_domains_performance_r8: For Mosaic "'//trim(test_type)// &
                      '", there will be '//text//' tiles on each processor.')
      allocate(tile(ntile_per_pe))
      do n = 1, ntile_per_pe
        tile(n) = pe*ntile_per_pe + n
      end do
      do n = 1, ntiles
        pe_start(n) = (n-1)/ntile_per_pe
        pe_end(n)   = pe_start(n)
      end do
      layout = 1
    else
      call mpp_error(NOTE,'update_domains_performance_r8: npes should be multiple of ntiles or ' // &
          'ntiles should be multiple of npes. No test is done for '//trim(test_type) )
      return
    end if

    do n = 1, ntiles
      global_indices(:,n) = (/1,nx,1,ny/)
      layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    if(single_tile) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
      tile1(1) = 1; tile2(1) = 1
      istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
      istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
          istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
          pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
          name = test_type, symmetry = .false. )
    else if(folded_north) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
      tile1(1) = 1; tile2(1) = 1
      istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (NORTH) and tile 1 (NORTH)  --- folded-north-edge
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = ny;  jend1(2) = ny
      istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                            pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                            name = test_type, symmetry = .false.  )
    else if( cubic_grid ) then
      call define_cubic_mosaic(test_type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                              global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x_save (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( a (ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
          do i = isc, iec
              x(i, j, k, l) = tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          enddo
        enddo
      enddo
    enddo

    a  = x
    x_save = x

    if(num_fields<1) then
      call mpp_error(FATAL, "update_domains_performanc_r8: num_fields must be a positive integer")
    endif

    id1 = mpp_clock_id( test_type, flags=MPP_CLOCK_SYNC)
    id_single = mpp_clock_id( test_type//' non-blocking', flags=MPP_CLOCK_SYNC)

    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id1)

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      id_update_single =  mpp_start_update_domains(a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( test_type//' group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( test_type//' group non-blocking', flags=MPP_CLOCK_SYNC )

    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem,jsm:jem,nz, num_fields) )
      allocate( a1(ism:iem,jsm:jem,nz, num_fields) )
      if(mix_2D_3D) allocate( a1_2D(ism:iem,jsm:jem,num_fields) )

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          if(mix_2D_3D) a1_2D(:,:,l) = x_save(:,:,1,1)
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains( x1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id1)

        ! non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
          do l = 1, num_fields
            if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), domain, complete=.false., tile_count=1)
            id_update =  mpp_start_update_domains(a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
          enddo
        else
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), domain, &
              update_id=id_update, complete=.false., tile_count=1)
              id_update = mpp_start_update_domains(a1(:,:,:,l), domain, update_id=id_update, &
                                                   complete=l==num_fields, tile_count=1)
          enddo
        endif
        call mpp_clock_end  (id2)

        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
          if(mix_2D_3D) call mpp_complete_update_domains(id_update, a1_2D(:,:,l), domain, &
                                                         complete=.false., tile_count=1)
          call mpp_complete_update_domains(id_update, a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id2)


        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' X'//text)
        enddo
        if(mix_2D_3D)call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), test_type//' X 2D')
      enddo
      deallocate(x1, a1)
      if(mix_2D_3D) deallocate(a1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type)
    deallocate(x, a, x_save)
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    !--- setup data
    shift = 0
    if(single_tile .or. folded_north) then
      shift = 0
    else
      shift = 1
    endif

    allocate( x (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x_save (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y_save (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( a (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( b (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    x = 0
    y = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec+shift
          do i = isc, iec+shift
            x(i,j,k,l) = 1.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
            y(i,j,k,l) = 2.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
      end do
    enddo
    a  = x; b  = y
    x_save  = x; y_save  = y

    id1 = mpp_clock_id( trim(test_type)//' BGRID', flags=MPP_CLOCK_SYNC )
    id_single = mpp_clock_id( trim(test_type)//' BGRID non-blocking', flags=MPP_CLOCK_SYNC )
    !--- blocking update
    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=BGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id1)

    !--- non-blocking update
    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=BGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if (do_sleep) call sleep(1)

    id1 = mpp_clock_id( trim(test_type)//' BGRID group', flags=MPP_CLOCK_SYNC)
    id2 = mpp_clock_id( trim(test_type)//' BGRID group non-blocking', flags=MPP_CLOCK_SYNC)
    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      allocate( y1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      allocate( a1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      allocate( b1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      if(mix_2D_3D) then
        allocate( a1_2D(ism:iem+shift,jsm:jem+shift,num_fields) )
        allocate( b1_2D(ism:iem+shift,jsm:jem+shift,num_fields) )
      endif

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          y1(:,:,:,l) = y_save(:,:,:,1)
          b1(:,:,:,l) = y_save(:,:,:,1)
          if(mix_2D_3D) then
              a1_2D(:,:,l) = x_save(:,:,1,1)
              b1_2D(:,:,l) = y_save(:,:,1,1)
          endif
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains( x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=BGRID_NE, &
                                  complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id1)

        !--- non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
          do l = 1, num_fields
            if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                          gridtype=BGRID_NE, complete=.false.)
            id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, &
                        gridtype=BGRID_NE, complete=l==num_fields, tile_count=1)
          enddo
        else
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, gridtype=BGRID_NE, &
                              update_id=id_update, complete=.false.)
            id_update =  mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, gridtype=BGRID_NE, &
                         update_id=id_update, complete=l==num_fields, tile_count=1)
          enddo
        endif
        call mpp_clock_end  (id2)

        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
          if(mix_2D_3D)call mpp_complete_update_domains(id_update, a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                                            gridtype=BGRID_NE, complete=.false.)
          call mpp_complete_update_domains(id_update, a1(:,:,:,l), b1(:,:,:,l), domain, &
                                            gridtype=BGRID_NE, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id2)

        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' BGRID X'//text)
          call compare_checksums( y1(:,:,:,l), b1(:,:,:,l), test_type//' BGRID Y'//text)
          if(mix_2D_3D) then
            call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' BGRID X'//text)
            call compare_checksums( y1(:,:,:,1), b1(:,:,:,1), test_type//' BGRID Y'//text)
          endif
        enddo
        if(mix_2D_3D) then
          call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), test_type//' BGRID X 2D')
          call compare_checksums( y1(:,:,1,:), b1_2D(:,:,:), test_type//' BGRID Y 2D')
        endif
      enddo
      deallocate(x1, y1, a1, b1)
      if(mix_2D_3D) deallocate(a1_2D, b1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=BGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end(id_single)

    !--- compare checksums
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' BGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' BGRID Y')

    deallocate(x, y, a, b, x_save, y_save)
    !------------------------------------------------------------------
    !              vector update : CGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    allocate( x (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
    allocate( y (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( a (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
    allocate( b (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x_save (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
    allocate( y_save (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )

    x = 0
    y = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
          do i = isc, iec+shift
            x(i,j,k,l) = 1.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
        do j = jsc, jec+shift
          do i = isc, iec
            y(i,j,k,l) = 2.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
      end do
    enddo

    a  = x; b  = y
    x_save  = x; y_save  = y

    id1 = mpp_clock_id( trim(test_type)//' CGRID', flags=MPP_CLOCK_SYNC )
    id_single = mpp_clock_id( trim(test_type)//' CGRID non-blocking', flags=MPP_CLOCK_SYNC )

    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=CGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id1)

    !--- non-blocking update
    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=CGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( trim(test_type)//' CGRID group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( trim(test_type)//' CGRID group non-blocking', flags=MPP_CLOCK_SYNC )

    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem+shift,jsm:jem      ,nz,num_fields) )
      allocate( y1(ism:iem      ,jsm:jem+shift,nz,num_fields) )
      allocate( a1(ism:iem+shift,jsm:jem      ,nz,num_fields) )
      allocate( b1(ism:iem      ,jsm:jem+shift,nz,num_fields) )
      if(mix_2D_3D) then
        allocate( a1_2D(ism:iem+shift,jsm:jem      ,num_fields) )
        allocate( b1_2D(ism:iem      ,jsm:jem+shift,num_fields) )
      endif

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          y1(:,:,:,l) = y_save(:,:,:,1)
          b1(:,:,:,l) = y_save(:,:,:,1)
          if(mix_2D_3D) then
            a1_2D(:,:,l) = x_save(:,:,1,1)
            b1_2D(:,:,l) = y_save(:,:,1,1)
          endif
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains( x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE, &
                                  complete=l==num_fields, tile_count=1 )
        enddo
        call mpp_clock_end  (id1)

        !--- non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, &
              gridtype=CGRID_NE, complete=.false.)
            id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, &
                        gridtype=CGRID_NE, complete=l==num_fields, tile_count=1)
          enddo
        else
          do l = 1, num_fields
            if(mix_2D_3D)id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, gridtype=CGRID_NE, &
                          update_id=id_update, complete=.false., tile_count=1)
            id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, gridtype=CGRID_NE, &
                          update_id=id_update, complete=l==num_fields, tile_count=1)
          enddo
        endif
        call mpp_clock_end  (id2)

        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
          if(mix_2D_3D) call mpp_complete_update_domains(id_update, a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                gridtype=CGRID_NE, complete=.false., tile_count=1)
          call mpp_complete_update_domains(id_update, a1(:,:,:,l), b1(:,:,:,l), domain, &
                gridtype=CGRID_NE, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id2)

        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' CGRID X'//text)
          call compare_checksums( y1(:,:,:,l), b1(:,:,:,l), test_type//' CGRID Y'//text)
        enddo
        if(mix_2D_3D) then
          call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), test_type//' BGRID X 2D')
          call compare_checksums( y1(:,:,1,:), b1_2D(:,:,:), test_type//' BGRID Y 2D')
        endif
      enddo
      deallocate(x1, y1, a1, b1)
      if(mix_2D_3D) deallocate(a1_2D, b1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=CGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !--- compare checksums

    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' CGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' CGRID Y')

    deallocate(x, y, a, b, x_save, y_save)

    !------------------------------------------------------------------
    !              vector update : AGRID vector and scalar pair
    !------------------------------------------------------------------
    allocate( x (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( y (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( a (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( b (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( x_save (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( y_save (ism:iem,jsm:jem,nz,ntile_per_pe) )

    x = 0
    y = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
          do i = isc, iec+shift
            x(i,j,k,l) = 1.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
        do j = jsc, jec+shift
          do i = isc, iec
            y(i,j,k,l) = 2.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
      end do
    enddo

    a  = x; b  = y
    x_save  = x; y_save  = y
    ! blocking update
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=AGRID, tile_count=l)
    enddo

    id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=AGRID)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=AGRID, tile_count=l)
    enddo

    !--- compare checksum
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' AGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' AGRID Y')

    x = x_save; y = y_save
    a = x_save; b = y_save
    ! blocking update
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=AGRID, flags = SCALAR_PAIR, tile_count=l)

      id_update_single = mpp_start_update_domains(a, b, domain, gridtype=AGRID, &
                                                  flags = SCALAR_PAIR, tile_count=l)
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=AGRID, &
                                       flags = SCALAR_PAIR, tile_count=l)
    enddo
    !--- compare checksums
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' AGRID SCALAR-PAIR X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' AGRID SCALAR-PAIR Y')

    deallocate(x, y, a, b, x_save, y_save)

    nx = nx_save
    ny = ny_save

    deallocate(layout2D, global_indices, pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 )

  end subroutine update_domains_performance_r8

  !> run performance tests on 32-bit real arrays
  subroutine update_domains_performance_r4( test_type )
    character(len=*), intent(in) :: test_type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile, ntile_per_pe
    integer        :: i, j, k, l, n, shift
    integer        :: ism, iem, jsm, jem
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real(kind=r4_kind),    allocatable, dimension(:,:,:,:) :: x, x1, y, y1, x_save, y_save
    real(kind=r4_kind),    allocatable, dimension(:,:,:,:) :: a, a1, b, b1
    real(kind=r4_kind),    allocatable, dimension(:,:,:  ) :: a1_2D, b1_2D
    integer            :: id_update
    integer            :: id1, id2
    logical            :: folded_north
    logical            :: cubic_grid, single_tile
    character(len=3)   :: text
    integer            :: nx_save, ny_save
    integer            :: id_single, id_update_single

    folded_north       = .false.
    cubic_grid         = .false.
    single_tile        = .false.
    nx_save = nx
    ny_save = ny
    !--- check the test_type
    select case(test_type)
      case ( 'Single-Tile' )   !--- single with cyclic along x- and y-direction
        single_tile = .true.
        ntiles = 1
        num_contact = 2
      case ( 'Folded-north' )
        ntiles = 1
        num_contact = 2
        folded_north = .true.
      case ( 'Cubic-Grid' )
        if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'update_domains_performance_r4: for Cubic_grid mosaic, nx_cubic is zero, '//&
                        'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'update_domains_performance_r4: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                         'No test is done for Cubic-Grid mosaic. ' )
          return
        endif

        nx = nx_cubic
        ny = ny_cubic
        ntiles = 6
        num_contact = 12
        cubic_grid = .true.

      case default
        call mpp_error(FATAL, 'update_domains_performancez_r4: no such test: '//test_type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    if( mod(npes, ntiles) == 0 ) then
      npes_per_tile = npes/ntiles
      text=""
      write(text,"(I3)") npes_per_tile
      call mpp_error(NOTE,'update_domains_performance_r4: For Mosaic "'//trim(test_type)// &
                     '", each tile will be distributed over '//text//' processors.')
      ntile_per_pe = 1
      allocate(tile(ntile_per_pe))
      tile = pe/npes_per_tile+1
      if(cubic_grid) then
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      else
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      endif
      do n = 1, ntiles
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile-1
      end do
    else if ( mod(ntiles, npes) == 0 ) then
      ntile_per_pe = ntiles/npes
      text=""
      write(text,"(I3)") ntile_per_pe
      call mpp_error(NOTE,'update_domains_performance_r4: For Mosaic "'//trim(test_type)// &
                      '", there will be '//text//' tiles on each processor.')
      allocate(tile(ntile_per_pe))
      do n = 1, ntile_per_pe
        tile(n) = pe*ntile_per_pe + n
      end do
      do n = 1, ntiles
        pe_start(n) = (n-1)/ntile_per_pe
        pe_end(n)   = pe_start(n)
      end do
      layout = 1
    else
      call mpp_error(NOTE,'update_domains_performance_r4: npes should be multiple of ntiles or ' // &
          'ntiles should be multiple of npes. No test is done for '//trim(test_type) )
      return
    end if

    do n = 1, ntiles
      global_indices(:,n) = (/1,nx,1,ny/)
      layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    if(single_tile) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
      tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
      istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
           istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
           pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
           name = test_type, symmetry = .false. )
    else if(folded_north) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
      tile1(1) = 1; tile2(1) = 1
      istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (NORTH) and tile 1 (NORTH)  --- folded-north-edge
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = ny;  jend1(2) = ny
      istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                            pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                            name = test_type, symmetry = .false.  )
    else if( cubic_grid ) then
      call define_cubic_mosaic(test_type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                              global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x_save (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( a (ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
           do i = isc, iec
              x(i, j, k, l) = tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
           enddo
        enddo
      enddo
    enddo

    a  = x
    x_save = x

    if(num_fields<1) then
      call mpp_error(FATAL, "update_domains_performanc_r4: num_fields must be a positive integer")
    endif

    id1 = mpp_clock_id( test_type, flags=MPP_CLOCK_SYNC)
    id_single = mpp_clock_id( test_type//' non-blocking', flags=MPP_CLOCK_SYNC)

    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, domain, tile_count=l)
    enddo
    call mpp_clock_end(id1)

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      id_update_single =  mpp_start_update_domains(a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( test_type//' group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( test_type//' group non-blocking', flags=MPP_CLOCK_SYNC )

    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem,jsm:jem,nz, num_fields) )
      allocate( a1(ism:iem,jsm:jem,nz, num_fields) )
      if(mix_2D_3D) allocate( a1_2D(ism:iem,jsm:jem,num_fields) )

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          if(mix_2D_3D) a1_2D(:,:,l) = x_save(:,:,1,1)
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains( x1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id1)

        ! non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
           do l = 1, num_fields
             if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), domain, complete=.false., tile_count=1)
             id_update = mpp_start_update_domains(a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
           enddo
        else
           do l = 1, num_fields
             if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), domain, &
                            update_id=id_update, complete=.false., tile_count=1)
             id_update = mpp_start_update_domains(a1(:,:,:,l), domain, update_id=id_update, &
                                                  complete=l==num_fields, tile_count=1)
           enddo
        endif
        call mpp_clock_end  (id2)

        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
          if(mix_2D_3D) call mpp_complete_update_domains(id_update, a1_2D(:,:,l), domain, complete=.false., tile_count=1)
          call mpp_complete_update_domains(id_update, a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id2)


        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' X'//text)
        enddo
        if(mix_2D_3D)call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), test_type//' X 2D')
      enddo
      deallocate(x1, a1)
      if(mix_2D_3D) deallocate(a1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type)
    deallocate(x, a, x_save)
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    !--- set up the data
    shift = 0
    if(single_tile .or. folded_north) then
      shift = 0
    else
      shift = 1
    endif

    allocate( x (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x_save (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y_save (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( a (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( b (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    x = 0
    y = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec+shift
           do i = isc, iec+shift
             x(i,j,k,l) = 1.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             y(i,j,k,l) = 2.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
           end do
        end do
      end do
    enddo
    a  = x; b  = y
    x_save  = x; y_save  = y

    id1 = mpp_clock_id( trim(test_type)//' BGRID', flags=MPP_CLOCK_SYNC )
    id_single = mpp_clock_id( trim(test_type)//' BGRID non-blocking', flags=MPP_CLOCK_SYNC )
    !--- blocking update
    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=BGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id1)

    !--- non-blocking update
    call mpp_clock_begin(id_single)
    id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=BGRID_NE)
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if (do_sleep) call sleep(1)

    id1 = mpp_clock_id( trim(test_type)//' BGRID group', flags=MPP_CLOCK_SYNC)
    id2 = mpp_clock_id( trim(test_type)//' BGRID group non-blocking', flags=MPP_CLOCK_SYNC)
    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      allocate( y1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      allocate( a1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      allocate( b1(ism:iem+shift,jsm:jem+shift,nz,num_fields) )
      if(mix_2D_3D) then
        allocate( a1_2D(ism:iem+shift,jsm:jem+shift,num_fields) )
        allocate( b1_2D(ism:iem+shift,jsm:jem+shift,num_fields) )
      endif

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          y1(:,:,:,l) = y_save(:,:,:,1)
          b1(:,:,:,l) = y_save(:,:,:,1)
          if(mix_2D_3D) then
            a1_2D(:,:,l) = x_save(:,:,1,1)
            b1_2D(:,:,l) = y_save(:,:,1,1)
          endif
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains( x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=BGRID_NE, &
                                  complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id1)

        !--- non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
          do l = 1, num_fields
            if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                           gridtype=BGRID_NE, complete=.false., tile_count=1)
            id_update =  mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, &
                           gridtype=BGRID_NE, complete=l==num_fields, tile_count=1)
          enddo
        else
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, gridtype=BGRID_NE, &
                              update_id=id_update, complete=.false., tile_count=1)
            id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, gridtype=BGRID_NE, &
                              update_id=id_update, complete=l==num_fields, tile_count=1)
          enddo
        endif
        call mpp_clock_end(id2)
        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
          if(mix_2D_3D) call mpp_complete_update_domains(id_update, a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                                            gridtype=BGRID_NE, complete=.false., tile_count=1)
          call mpp_complete_update_domains(id_update, a1(:,:,:,l), b1(:,:,:,l), domain, &
                                            gridtype=BGRID_NE, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end(id2)

        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' BGRID X'//text)
          call compare_checksums( y1(:,:,:,l), b1(:,:,:,l), test_type//' BGRID Y'//text)
          if(mix_2D_3D) then
            call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' BGRID X'//text)
            call compare_checksums( y1(:,:,:,1), b1(:,:,:,1), test_type//' BGRID Y'//text)
          endif
        enddo
        if(mix_2D_3D) then
          call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), test_type//' BGRID X 2D')
          call compare_checksums( y1(:,:,1,:), b1_2D(:,:,:), test_type//' BGRID Y 2D')
        endif
      enddo
      deallocate(x1, y1, a1, b1)
      if(mix_2D_3D) deallocate(a1_2D, b1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=BGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end(id_single)

    !--- compare checksums
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' BGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' BGRID Y')

    deallocate(x, y, a, b, x_save, y_save)
    !------------------------------------------------------------------
    !              vector update : CGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    allocate( x (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
    allocate( y (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( a (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
    allocate( b (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x_save (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
    allocate( y_save (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )

    x = 0
    y = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
          do i = isc, iec+shift
            x(i,j,k,l) = 1.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
        do j = jsc, jec+shift
          do i = isc, iec
            y(i,j,k,l) = 2.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
      end do
    enddo

    a  = x; b  = y
    x_save  = x; y_save  = y

    id1 = mpp_clock_id( trim(test_type)//' CGRID', flags=MPP_CLOCK_SYNC )
    id_single = mpp_clock_id( trim(test_type)//' CGRID non-blocking', flags=MPP_CLOCK_SYNC )

    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=CGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id1)

    !--- non-blocking update
    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=CGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( trim(test_type)//' CGRID group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( trim(test_type)//' CGRID group non-blocking', flags=MPP_CLOCK_SYNC )

    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem+shift,jsm:jem      ,nz,num_fields) )
      allocate( y1(ism:iem      ,jsm:jem+shift,nz,num_fields) )
      allocate( a1(ism:iem+shift,jsm:jem      ,nz,num_fields) )
      allocate( b1(ism:iem      ,jsm:jem+shift,nz,num_fields) )
      if(mix_2D_3D) then
        allocate( a1_2D(ism:iem+shift,jsm:jem      ,num_fields) )
        allocate( b1_2D(ism:iem      ,jsm:jem+shift,num_fields) )
      endif

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          y1(:,:,:,l) = y_save(:,:,:,1)
          b1(:,:,:,l) = y_save(:,:,:,1)
          if(mix_2D_3D) then
            a1_2D(:,:,l) = x_save(:,:,1,1)
            b1_2D(:,:,l) = y_save(:,:,1,1)
          endif
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains(x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE, &
                                  complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end(id1)
        !--- non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                     gridtype=CGRID_NE, complete=.false., tile_count=1)
            id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, &
                     gridtype=CGRID_NE, complete=l==num_fields, tile_count=1)
          enddo
        else
          do l = 1, num_fields
            if(mix_2D_3D)id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, gridtype=CGRID_NE, &
                          update_id=id_update, complete=.false., tile_count=1)
            id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, gridtype=CGRID_NE, &
                          update_id=id_update, complete=l==num_fields, tile_count=1)
          enddo
        endif
        call mpp_clock_end  (id2)

        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
           if(mix_2D_3D)call mpp_complete_update_domains(id_update, a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                 gridtype=CGRID_NE, complete=.false., tile_count=1)
           call mpp_complete_update_domains(id_update, a1(:,:,:,l), b1(:,:,:,l), domain, &
                 gridtype=CGRID_NE, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id2)

        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), test_type//' CGRID X'//text)
          call compare_checksums( y1(:,:,:,l), b1(:,:,:,l), test_type//' CGRID Y'//text)
        enddo
        if(mix_2D_3D) then
          call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), test_type//' BGRID X 2D')
          call compare_checksums( y1(:,:,1,:), b1_2D(:,:,:), test_type//' BGRID Y 2D')
        endif
      enddo
      deallocate(x1, y1, a1, b1)
      if(mix_2D_3D) deallocate(a1_2D, b1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=CGRID_NE, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !--- compare checksum
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' CGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' CGRID Y')

    deallocate(x, y, a, b, x_save, y_save)

    !------------------------------------------------------------------
    !              vector update : AGRID vector and scalar pair
    !------------------------------------------------------------------
    allocate( x (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( y (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( a (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( b (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( x_save (ism:iem,jsm:jem,nz,ntile_per_pe) )
    allocate( y_save (ism:iem,jsm:jem,nz,ntile_per_pe) )


    x = 0
    y = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
          do i = isc, iec+shift
            x(i,j,k,l) = 1.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
        do j = jsc, jec+shift
          do i = isc, iec
            y(i,j,k,l) = 2.0e3 + tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
        end do
      end do
    enddo

    a  = x; b  = y
    x_save  = x; y_save  = y
    ! blocking update
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=AGRID, tile_count=l)

      id_update_single = mpp_start_update_domains(a, b, domain, gridtype=AGRID, tile_count=l)
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=AGRID, tile_count=l)
    enddo
    !--- compare checksum
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' AGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' AGRID Y')

    x = x_save; y = y_save
    a = x_save; b = y_save
    ! blocking update
    do l=1,ntile_per_pe
      call mpp_update_domains( x, y, domain, gridtype=AGRID, flags = SCALAR_PAIR, tile_count=l)

      id_update_single = mpp_start_update_domains(a, b, domain, gridtype=AGRID, &
                                                  flags = SCALAR_PAIR, tile_count=l)
      call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=AGRID, &
                                       flags = SCALAR_PAIR, tile_count=l)
    enddo
    !--- compare checksums
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), test_type//' AGRID SCALAR-PAIR X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), test_type//' AGRID SCALAR-PAIR Y')

    deallocate(x, y, a, b, x_save, y_save)

    nx = nx_save
    ny = ny_save

    deallocate(layout2D, global_indices, pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 )

  end subroutine update_domains_performance_r4

 !> run performance tests on 64-bit integer arrays
  subroutine update_domains_performance_i8( test_type )
    character(len=*), intent(in) :: test_type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile, ntile_per_pe
    integer        :: l, n, shift
    integer(kind=i8_kind) :: i, j, k ! used to define the data arrays
    integer        :: ism, iem, jsm, jem
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    integer(kind=i8_kind),    allocatable, dimension(:,:,:,:) :: x, x1, x_save
    integer(kind=i8_kind),    allocatable, dimension(:,:,:,:) :: a, a1
    integer(kind=i8_kind),     allocatable, dimension(:,:,:  ) :: a1_2D
    integer            :: id_update
    integer            :: id1, id2
    logical            :: folded_north
    logical            :: cubic_grid, single_tile
    character(len=3)   :: text
    integer            :: nx_save, ny_save
    integer            :: id_single, id_update_single

    folded_north       = .false.
    cubic_grid         = .false.
    single_tile        = .false.
    nx_save = nx
    ny_save = ny
    !--- check the test_type
    select case(test_type)
      case ( 'Single-Tile' )   !--- single with cyclic along x- and y-direction
        single_tile = .true.
        ntiles = 1
        num_contact = 2
      case ( 'Folded-north' )
        ntiles = 1
        num_contact = 2
        folded_north = .true.
      case ( 'Cubic-Grid' )
        if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'update_domains_performance_r8: for Cubic_grid mosaic, nx_cubic is zero, '//&
                        'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'update_domains_performance_r8: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                        'No test is done for Cubic-Grid mosaic. ' )
          return
        endif

        nx = nx_cubic
        ny = ny_cubic
        ntiles = 6
        num_contact = 12
        cubic_grid = .true.

      case default
        call mpp_error(FATAL, 'update_domains_performancez_r8: no such test: '//test_type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    if( mod(npes, ntiles) == 0 ) then
      npes_per_tile = npes/ntiles
      text=""
      write(text,"(I3)") npes_per_tile
      call mpp_error(NOTE, 'update_domains_performance_r8: For Mosaic "'//trim(test_type)// &
                    '", each tile will be distributed over '//text//' processors.')
      ntile_per_pe = 1
      allocate(tile(ntile_per_pe))
      tile = pe/npes_per_tile+1
      if(cubic_grid) then
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      else
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      endif
      do n = 1, ntiles
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile-1
      end do
    else if ( mod(ntiles, npes) == 0 ) then
      ntile_per_pe = ntiles/npes
      text=""
      write(text,"(I3)") ntile_per_pe
      call mpp_error(NOTE, 'update_domains_performance_r8: For Mosaic "'//trim(test_type)// &
                      '", there will be '//text//' tiles on each processor.')
      allocate(tile(ntile_per_pe))
      do n = 1, ntile_per_pe
        tile(n) = pe*ntile_per_pe + n
      end do
      do n = 1, ntiles
        pe_start(n) = (n-1)/ntile_per_pe
        pe_end(n)   = pe_start(n)
      end do
      layout = 1
    else
      call mpp_error(NOTE,'update_domains_performance_r8: npes should be multiple of ntiles or ' // &
          'ntiles should be multiple of npes. No test is done for '//trim(test_type) )
      return
    end if

    do n = 1, ntiles
      global_indices(:,n) = (/1,nx,1,ny/)
      layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    if(single_tile) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
      tile1(1) = 1; tile2(1) = 1
      istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
      istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
          istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
          pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
          name = test_type, symmetry = .false. )
    else if(folded_north) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
      tile1(1) = 1; tile2(1) = 1
      istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (NORTH) and tile 1 (NORTH)  --- folded-north-edge
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = ny;  jend1(2) = ny
      istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                            pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                            name = test_type, symmetry = .false.  )
    else if( cubic_grid ) then
      call define_cubic_mosaic(test_type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                              global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x_save (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( a (ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
          do i = isc, iec
              x(i, j, k, l) = tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          enddo
        enddo
      enddo
    enddo

    a  = x
    x_save = x

    if(num_fields<1) then
      call mpp_error(FATAL, "update_domains_performanc_r8: num_fields must be a positive integer")
    endif

    id1 = mpp_clock_id( test_type, flags=MPP_CLOCK_SYNC)
    id_single = mpp_clock_id( test_type//' non-blocking', flags=MPP_CLOCK_SYNC)

    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id1)

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      id_update_single =  mpp_start_update_domains(a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( test_type//' group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( test_type//' group non-blocking', flags=MPP_CLOCK_SYNC )

    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem,jsm:jem,nz, num_fields) )
      allocate( a1(ism:iem,jsm:jem,nz, num_fields) )
      if(mix_2D_3D) allocate( a1_2D(ism:iem,jsm:jem,num_fields) )

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          if(mix_2D_3D) a1_2D(:,:,l) = x_save(:,:,1,1)
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains( x1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id1)

        ! non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
          do l = 1, num_fields
            if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), domain, complete=.false., tile_count=1)
            id_update =  mpp_start_update_domains(a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
          enddo
        else
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), domain, &
              update_id=id_update, complete=.false., tile_count=1)
              id_update = mpp_start_update_domains(a1(:,:,:,l), domain, update_id=id_update, &
                                                   complete=l==num_fields, tile_count=1)
          enddo
        endif
        call mpp_clock_end  (id2)

        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
          if(mix_2D_3D) call mpp_complete_update_domains(id_update, a1_2D(:,:,l), domain, &
                                                         complete=.false., tile_count=1)
          call mpp_complete_update_domains(id_update, a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id2)

        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums_int( x1(:,:,:,l), a1(:,:,:,l), test_type//' X'//text)
        enddo
        if(mix_2D_3D)call compare_checksums_int( x1(:,:,1,:), a1_2D(:,:,:), test_type//' X 2D')
      enddo
      deallocate(x1, a1)
      if(mix_2D_3D) deallocate(a1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)
    call compare_checksums_int( x(:,:,:,1), a(:,:,:,1), test_type)

    deallocate(x, a, x_save)
    deallocate(layout2D, global_indices, pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 )

  end subroutine update_domains_performance_i8

  !> run performance tests on 32-bit integer arrays
  subroutine update_domains_performance_i4( test_type )
    character(len=*), intent(in) :: test_type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile, ntile_per_pe
    integer        :: l, n, shift
    integer(kind=i4_kind) :: i, j, k ! used to define the data arrays
    integer        :: ism, iem, jsm, jem
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    integer(kind=i4_kind),    allocatable, dimension(:,:,:,:) :: x, x1, x_save
    integer(kind=i4_kind),    allocatable, dimension(:,:,:,:) :: a, a1
    integer(kind=i4_kind),    allocatable, dimension(:,:,:  ) :: a1_2D
    integer            :: id_update
    integer            :: id1, id2
    logical            :: folded_north
    logical            :: cubic_grid, single_tile
    character(len=3)   :: text
    integer            :: nx_save, ny_save
    integer            :: id_single, id_update_single

    folded_north       = .false.
    cubic_grid         = .false.
    single_tile        = .false.
    nx_save = nx
    ny_save = ny
    !--- check the test_type
    select case(test_type)
      case ( 'Single-Tile' )   !--- single with cyclic along x- and y-direction
        single_tile = .true.
        ntiles = 1
        num_contact = 2
      case ( 'Folded-north' )
        ntiles = 1
        num_contact = 2
        folded_north = .true.
      case ( 'Cubic-Grid' )
        if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'update_domains_performance_r8: for Cubic_grid mosaic, nx_cubic is zero, '//&
                        'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'update_domains_performance_r8: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                        'No test is done for Cubic-Grid mosaic. ' )
          return
        endif

        nx = nx_cubic
        ny = ny_cubic
        ntiles = 6
        num_contact = 12
        cubic_grid = .true.

      case default
        call mpp_error(FATAL, 'update_domains_performancez_r8: no such test: '//test_type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    if( mod(npes, ntiles) == 0 ) then
      npes_per_tile = npes/ntiles
      text=""
      write(text,"(I3)") npes_per_tile
      call mpp_error(NOTE, 'update_domains_performance_r8: For Mosaic "'//trim(test_type)// &
                    '", each tile will be distributed over '//text//' processors.')
      ntile_per_pe = 1
      allocate(tile(ntile_per_pe))
      tile = pe/npes_per_tile+1
      if(cubic_grid) then
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      else
        call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
      endif
      do n = 1, ntiles
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile-1
      end do
    else if ( mod(ntiles, npes) == 0 ) then
      ntile_per_pe = ntiles/npes
      text=""
      write(text,"(I3)") ntile_per_pe
      call mpp_error(NOTE, 'update_domains_performance_r8: For Mosaic "'//trim(test_type)// &
                      '", there will be '//text//' tiles on each processor.')
      allocate(tile(ntile_per_pe))
      do n = 1, ntile_per_pe
        tile(n) = pe*ntile_per_pe + n
      end do
      do n = 1, ntiles
        pe_start(n) = (n-1)/ntile_per_pe
        pe_end(n)   = pe_start(n)
      end do
      layout = 1
    else
      call mpp_error(NOTE,'update_domains_performance_r8: npes should be multiple of ntiles or ' // &
          'ntiles should be multiple of npes. No test is done for '//trim(test_type) )
      return
    end if

    do n = 1, ntiles
      global_indices(:,n) = (/1,nx,1,ny/)
      layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    if(single_tile) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)
      tile1(1) = 1; tile2(1) = 1
      istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (NORTH)  --- cyclic
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx; jstart1(2) = 1;   jend1(2) = 1
      istart2(2) = 1;  iend2(2) = nx; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
          istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
          pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
          name = test_type, symmetry = .false. )
    else if(folded_north) then
      !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
      tile1(1) = 1; tile2(1) = 1
      istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
      istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
      !--- Contact line 2, between tile 1 (NORTH) and tile 1 (NORTH)  --- folded-north-edge
      tile1(2) = 1; tile2(2) = 1
      istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = ny;  jend1(2) = ny
      istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = ny;  jend2(2) = ny
      call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                            pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                            name = test_type, symmetry = .false.  )
    else if( cubic_grid ) then
      call define_cubic_mosaic(test_type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                              global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x_save (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( a (ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0
    do l = 1, ntile_per_pe
      do k = 1, nz
        do j = jsc, jec
          do i = isc, iec
              x(i, j, k, l) = tile(l) + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          enddo
        enddo
      enddo
    enddo

    a  = x
    x_save = x

    if(num_fields<1) then
      call mpp_error(FATAL, "update_domains_performanc_r8: num_fields must be a positive integer")
    endif

    id1 = mpp_clock_id( test_type, flags=MPP_CLOCK_SYNC)
    id_single = mpp_clock_id( test_type//' non-blocking', flags=MPP_CLOCK_SYNC)

    call mpp_clock_begin(id1)
    do l=1,ntile_per_pe
      call mpp_update_domains( x, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id1)

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      id_update_single =  mpp_start_update_domains(a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( test_type//' group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( test_type//' group non-blocking', flags=MPP_CLOCK_SYNC )

    if(ntile_per_pe == 1) then
      allocate( x1(ism:iem,jsm:jem,nz, num_fields) )
      allocate( a1(ism:iem,jsm:jem,nz, num_fields) )
      if(mix_2D_3D) allocate( a1_2D(ism:iem,jsm:jem,num_fields) )

      do n = 1, num_iter
        do l = 1, num_fields
          x1(:,:,:,l) = x_save(:,:,:,1)
          a1(:,:,:,l) = x_save(:,:,:,1)
          if(mix_2D_3D) a1_2D(:,:,l) = x_save(:,:,1,1)
        enddo

        call mpp_clock_begin(id1)
        do l = 1, num_fields
          call mpp_update_domains( x1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id1)

        ! non-blocking update
        call mpp_clock_begin(id2)
        if( n == 1 ) then
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), domain, complete=.false., tile_count=1)
            id_update =  mpp_start_update_domains(a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
          enddo
        else
          do l = 1, num_fields
            if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), domain, &
              update_id=id_update, complete=.false., tile_count=1)
              id_update = mpp_start_update_domains(a1(:,:,:,l), domain, update_id=id_update, &
                                                   complete=l==num_fields, tile_count=1)
          enddo
        endif
        call mpp_clock_end  (id2)

        !---- sleep some time for non-blocking.
        if(do_sleep) call sleep(1)

        call mpp_clock_begin(id2)
        do l = 1, num_fields
          if(mix_2D_3D) call mpp_complete_update_domains(id_update, a1_2D(:,:,l), domain, &
                                                         complete=.false., tile_count=1)
          call mpp_complete_update_domains(id_update, a1(:,:,:,l), domain, complete=l==num_fields, tile_count=1)
        enddo
        call mpp_clock_end  (id2)


        !--- compare checksum
        do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums_int( x1(:,:,:,l), a1(:,:,:,l), test_type//' X'//text)
        enddo
        if(mix_2D_3D)call compare_checksums_int( x1(:,:,1,:), a1_2D(:,:,:), test_type//' X 2D')
      enddo
      deallocate(x1, a1)
      if(mix_2D_3D) deallocate(a1_2D)
    endif

    call mpp_clock_begin(id_single)
    do l=1,ntile_per_pe
      call mpp_complete_update_domains(id_update_single, a, domain, tile_count=l)
    enddo
    call mpp_clock_end  (id_single)
    call compare_checksums_int( x(:,:,:,1), a(:,:,:,1), test_type)

    deallocate(x, a, x_save)
    deallocate(layout2D, global_indices, pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 )

  end subroutine update_domains_performance_i4

  !> define mosaic domain for cubic grid
  subroutine define_cubic_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end, use_memsize)
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    logical, optional, intent(in) :: use_memsize
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)
    logical                       :: use_memsize_local

    use_memsize_local = .true.
    if(present(use_memsize)) use_memsize_local = use_memsize

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

end program test_update_domains_performance
