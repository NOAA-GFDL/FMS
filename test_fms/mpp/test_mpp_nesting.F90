!***********************************************************************
!*                   Gnu Lesser General Public License
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
!> Tests nested domain operations and routines in mpp_domains
program test_mpp_nesting

  use fms_mod
  use mpp_domains_mod
  use mpp_mod
  use compare_data_checksums
  use test_domains_utility_mod
  use platform_mod

  implicit none

  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, stackmax=4000000
  integer :: stdunit = 6
  logical :: debug=.FALSE., opened

  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  integer :: x_cyclic_offset = 3   ! to be used in test_cyclic_offset
  integer :: y_cyclic_offset = -4  ! to be used in test_cyclic_offset
  character(len=32) :: warn_level = "fatal"
  integer :: wide_halo_x = 0, wide_halo_y = 0
  integer :: nx_cubic = 0, ny_cubic = 0
  ! namelist flags to run each test
  logical :: test_nest = .false.
  logical :: test_nest_regional = .false.
  integer :: ensemble_size = 1
  integer :: layout_cubic(2) = (/0,0/)
  integer :: layout_tripolar(2) = (/0,0/)
  integer :: layout_ensemble(2) = (/0,0/)
  logical :: do_sleep = .false.
  integer :: num_iter = 1
  integer :: num_fields = 4

  logical :: mix_2D_3D = .false.
  logical :: test_subset = .false.
  integer :: nthreads = 1
  logical :: test_adjoint = .false.
  logical :: wide_halo = .false.
  logical :: test_unstruct = .false.

  !--- namelist variable for nest domain
  integer, parameter :: MAX_NNEST=20
  integer, parameter :: MAX_NCONTACT=100
  integer, parameter :: MAX_NTILE=50
  integer :: num_nest = 0
  integer :: nest_level(MAX_NNEST) = 1
  integer :: tile_fine(MAX_NNEST) = 0
  integer :: tile_coarse(MAX_NNEST) = 0
  integer :: ntiles_nest_all = 0
  integer :: npes_nest_tile(MAX_NTILE) = 0
  integer :: refine_ratio(MAX_NNEST) = 1
  integer :: istart_coarse(MAX_NNEST) = 0, icount_coarse(MAX_NNEST) = 0
  integer :: jstart_coarse(MAX_NNEST) = 0, jcount_coarse(MAX_NNEST) = 0
  integer :: extra_halo = 0
  integer :: layout_nest(2,MAX_NNEST) = 0
  character :: cyclic_nest(MAX_NCONTACT) = 'N'
  logical :: test_edge_nonblock = .false. !< enables non-blocking domain updates in edge update test
                                          !! currently fails with gcc

  namelist / test_mpp_domains_nml / nx, ny, nz, stackmax, debug, mpes, &
                               whalo, ehalo, shalo, nhalo, x_cyclic_offset, y_cyclic_offset, &
                               warn_level, wide_halo_x, wide_halo_y, nx_cubic, ny_cubic, &
                               num_fields, do_sleep, num_iter, &
                               test_nest, num_nest, ntiles_nest_all, nest_level, &
                               tile_fine, &
                               tile_coarse, refine_ratio, istart_coarse, icount_coarse, &
                               jstart_coarse, jcount_coarse, extra_halo, npes_nest_tile, &
                               cyclic_nest, mix_2D_3D, &
                               ensemble_size, layout_cubic, &
                               layout_ensemble, nthreads, layout_tripolar, &
                               test_nest_regional, wide_halo
  integer :: i, j, k, n
  integer :: layout(2)
  integer :: id
  integer :: outunit, errunit, io_status
  integer :: omp_get_num_threads, omp_get_thread_num
  integer :: ierr

  call mpp_init()
  if (debug) then
    call mpp_domains_init(MPP_DEBUG)
  else
    call mpp_domains_init()
  endif
  call mpp_domains_set_stack_size(stackmax)

  outunit = stdout()
  errunit = stderr()

  read (input_nml_file, test_mpp_domains_nml, iostat=io_status)
  if (io_status > 0) then
     call mpp_error(FATAL,'=>test_mpp_domains: Error reading input.nml')
  endif

  pe = mpp_pe()
  npes = mpp_npes()


  !! check valid namelist values
  !--- wide_halo_x and wide_halo_y must be either both 0 or both positive.
  if( wide_halo_x < 0 .OR. wide_halo_y < 0) call mpp_error(FATAL, &
     "test_mpp_domain: both wide_halo_x and wide_halo_y should be non-negative")
  if( wide_halo_x == 0 .NEQV. wide_halo_y == 0) call mpp_error(FATAL, &
     "test_mpp_domain: wide_halo_x and wide_halo_y should be both zero or both positive")

  !--- nx_cubic and ny_cubic must be either both 0 or both positive.
  if( nx_cubic < 0 .OR. ny_cubic < 0) call mpp_error(FATAL, &
     "test_mpp_domain: both nx_cubic and ny_cubic should be non-negative")
  if( nx_cubic == 0 .NEQV. ny_cubic == 0) call mpp_error(FATAL, &
     "test_mpp_domain: nx_cubic and ny_cubic should be both zero or both positive")

  ! nest domain update
  if( test_nest .and. (num_nest>0) ) then
     if (mpp_pe() == mpp_root_pe()) &
        print *, '-----------------> Calling test_update_nest_domain Cubic <----------------'
     do n = 1, num_nest
        if( istart_coarse(n) == 0 .OR. jstart_coarse(n) == 0 ) call mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable istart_coarse, jstart_coarse")
        if( icount_coarse(n) == 0 .OR. jcount_coarse(n) == 0 ) call mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable icount_coarse, jcount_coarse")
        if( tile_coarse(n) .LE. 0) call mpp_error(FATAL, &
            "test_mpp_domain: check the setting of namelist variable tile_coarse")
     enddo
     if(ANY(refine_ratio(:).LT.1)) call  mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable refine_ratio")
    if (mpp_pe() == mpp_root_pe()) &
      & print *, '-----------------> Starting test_update_nest_domain_r8 Cubic <----------------'
     call test_update_nest_domain_r8('Cubic-Grid')
    if (mpp_pe() == mpp_root_pe()) &
      & print *, '-----------------> Finished test_update_nest_domain_r8 Cubic <----------------'
    if (mpp_pe() == mpp_root_pe()) &
      & print *, '-----------------> Starting test_update_nest_domain_r4 Cubic <----------------'
     call test_update_nest_domain_r4('Cubic-Grid')
    if (mpp_pe() == mpp_root_pe()) &
      & print *, '-----------------> Finished test_update_nest_domain_r4 Cubic <----------------'
  endif
  ! regional nest domain update
  if( test_nest_regional .and. (num_nest>0) ) then
    if (mpp_pe() == mpp_root_pe()) &
      & print *, '---------------> Calling test_update_nest_domain Regional <--------------'
     do n = 1, num_nest
        if( istart_coarse(n) == 0 .OR. jstart_coarse(n) == 0 ) call mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable istart_coarse, jstart_coarse")
        if( icount_coarse(n) == 0 .OR. jcount_coarse(n) == 0 ) call mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable icount_coarse, jcount_coarse")
        if( tile_coarse(n) .LE. 0) call mpp_error(FATAL, &
            "test_mpp_domain: check the setting of namelist variable tile_coarse")
     enddo
     if(ANY(refine_ratio(:).LT.1)) call  mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable refine_ratio")
    if (mpp_pe() == mpp_root_pe()) &
     & print *, '---------------> Starting test_update_nest_domain_r8 Regional <--------------'
     call test_update_nest_domain_r8('Regional')
    if (mpp_pe() == mpp_root_pe()) &
     & print *, '---------------> Finished test_update_nest_domain_r8 Regional <--------------'
    if (mpp_pe() == mpp_root_pe()) &
     & print *, '---------------> Starting test_update_nest_domain_r4 Regional <--------------'
     call test_update_nest_domain_r4('Regional')
    if (mpp_pe() == mpp_root_pe()) &
     & print *, '---------------> Finished test_update_nest_domain_r4 Regional <--------------'
  endif
  ! single face cubed-sphere
  if(ntiles_nest_all == 4 .and. num_nest > 0) then
    if (mpp_pe() == mpp_root_pe()) &
     & print *, '--------------------> Calling test_update_nest_domain(cubed-sphere, single face) <-------------------'
     do n = 1, num_nest
        if( istart_coarse(n) == 0 .OR. jstart_coarse(n) == 0 ) call mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable istart_coarse, jstart_coarse")
        if( icount_coarse(n) == 0 .OR. jcount_coarse(n) == 0 ) call mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable icount_coarse, jcount_coarse")
        if( tile_coarse(n) .LE. 0) call mpp_error(FATAL, &
            "test_mpp_domain: check the setting of namelist variable tile_coarse")
     enddo
     if(ANY(refine_ratio(:).LT.1)) call  mpp_error(FATAL, &
        "test_mpp_domain: check the setting of namelist variable refine_ratio")
     call test_update_nest_domain_r8('Cubed-sphere, single face')
     call test_update_nest_domain_r4('Cubed-sphere, single face')
    if (mpp_pe() == mpp_root_pe()) &
      & print *, '--------------------> Finished test_update_nest_domain <-------------------'
  endif

  call mpp_exit()

  contains

!###############################################################################
! test halo update for grid nested in global cubic sphere grid. The nested region may cross the edge.
! It is assumed the boundary condition of nested region is solid wall in both direction.
  subroutine test_nest_halo_update_r4( domain )
    type(domain2D), intent(inout) :: domain
    integer        :: i, j, k, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed

    real(r4_kind),    allocatable, dimension(:,:,:) :: x, y
    real(r4_kind),    allocatable, dimension(:,:,:) :: global1, global2, global
    character(len=128) :: type
    integer :: nx, ny

    call mpp_get_global_domain(domain, xsize=nx, ysize=ny)
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !--- setup data
    allocate(global(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )
    global = 0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             global(i,j,k) = i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
       end do
    end do

    allocate( x (isd:ied,jsd:jed,nz) )
    x = 0.
    x(isc:iec,jsc:jec,:) = global(isc:iec,jsc:jec,:)

    type = 'nest grid scalar'
    call mpp_update_domains( x, domain, name=trim(type) )

    call compare_checksums( x(isd:ied,jsd:jed,:), global(isd:ied,jsd:jed,:), trim(type) )

    deallocate(global, x)
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction
    !------------------------------------------------------------------
    !--- setup data
    shift = 1

    allocate(global1(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
    allocate(global2(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
    global1 = 0; global2 = 0
    do k = 1, nz
       do j = 1, ny+shift
          do i = 1, nx+shift
             global1(i,j,k) = 1.0 + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             global2(i,j,k) = 2.0 + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
       end do
    end do

    allocate( x (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y (isd:ied+shift,jsd:jed+shift,nz) )

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec+shift,:) = global1(isc:iec+shift,jsc:jec+shift,:)
    y (isc:iec+shift,jsc:jec+shift,:) = global2(isc:iec+shift,jsc:jec+shift,:)

    type = 'nest grid BGRID_NE'
    call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE, name=trim(type))

    call compare_checksums( x (isd:ied+shift,jsd:jed+shift,:),  global1(isd:ied+shift,jsd:jed+shift,:), &
                          & trim(type)//' X' )
    call compare_checksums( y (isd:ied+shift,jsd:jed+shift,:),  global2(isd:ied+shift,jsd:jed+shift,:), &
                          & trim(type)//' Y' )

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    deallocate(global1, global2, x, y)
    allocate(global1(1-whalo:nx+shift+ehalo,1-shalo:ny+nhalo,nz) )
    allocate(global2(1-whalo:nx+ehalo,1-shalo:ny+shift+nhalo,nz) )
    allocate( x (isd:ied+shift,jsd:jed  ,nz) )
    allocate( y (isd:ied  ,jsd:jed+shift,nz) )
    global1 = 0; global2 = 0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+shift
             global1(i,j,k) = 1.0 + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
       end do
       do j = 1, ny+shift
          do i = 1, nx
             global2(i,j,k) = 2.0 + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
          end do
       end do
    end do

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec,:) = global1(isc:iec+shift,jsc:jec,:)
    y (isc:iec,jsc:jec+shift,:) = global2(isc:iec,jsc:jec+shift,:)

    type = "nest grid CGRID_NE"
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE, name=trim(type))

    call compare_checksums( x(isd:ied+shift,jsd:jed,:), global1(isd:ied+shift,jsd:jed,:), trim(type)//' X' )
    call compare_checksums( y(isd:ied,jsd:jed+shift,:), global2(isd:ied,jsd:jed+shift,:), trim(type)//' Y' )

    deallocate(global1, global2, x, y)

  end subroutine test_nest_halo_update_r4

  subroutine test_nest_halo_update_r8( domain )
    type(domain2D), intent(inout) :: domain
    integer        :: i, j, k, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed

    real(r8_kind),    allocatable, dimension(:,:,:) :: x, y
    real(r8_kind),    allocatable, dimension(:,:,:) :: global1, global2, global
    character(len=128) :: type
    integer :: nx, ny

    call mpp_get_global_domain(domain, xsize=nx, ysize=ny)
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !--- setup data
    allocate(global(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )
    global = 0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             global(i,j,k) = i*1.0e-3_r8_kind + j*1.0e-6_r8_kind + k*1.0e-9_r8_kind
          end do
       end do
    end do

    allocate( x (isd:ied,jsd:jed,nz) )
    x = 0.
    x(isc:iec,jsc:jec,:) = global(isc:iec,jsc:jec,:)

    type = 'nest grid scalar'
    call mpp_update_domains( x, domain, name=trim(type) )

    call compare_checksums( x(isd:ied,jsd:jed,:), global(isd:ied,jsd:jed,:), trim(type) )

    deallocate(global, x)
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction
    !------------------------------------------------------------------
    !--- setup data
    shift = 1

    allocate(global1(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
    allocate(global2(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
    global1 = 0; global2 = 0
    do k = 1, nz
       do j = 1, ny+shift
          do i = 1, nx+shift
             global1(i,j,k) = 1.0_r8_kind + i*1.0e-3_r8_kind + j*1.0e-6_r8_kind + k*1.0e-9_r8_kind
             global2(i,j,k) = 2.0_r8_kind + i*1.0e-3_r8_kind + j*1.0e-6_r8_kind + k*1.0e-9_r8_kind
          end do
       end do
    end do

    allocate( x (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y (isd:ied+shift,jsd:jed+shift,nz) )

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec+shift,:) = global1(isc:iec+shift,jsc:jec+shift,:)
    y (isc:iec+shift,jsc:jec+shift,:) = global2(isc:iec+shift,jsc:jec+shift,:)

    type = 'nest grid BGRID_NE'
    call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE, name=trim(type))

    call compare_checksums( x (isd:ied+shift,jsd:jed+shift,:),  global1(isd:ied+shift,jsd:jed+shift,:), &
                          & trim(type)//' X' )
    call compare_checksums( y (isd:ied+shift,jsd:jed+shift,:),  global2(isd:ied+shift,jsd:jed+shift,:), &
                          & trim(type)//' Y' )

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    deallocate(global1, global2, x, y)
    allocate(global1(1-whalo:nx+shift+ehalo,1-shalo:ny+nhalo,nz) )
    allocate(global2(1-whalo:nx+ehalo,1-shalo:ny+shift+nhalo,nz) )
    allocate( x (isd:ied+shift,jsd:jed  ,nz) )
    allocate( y (isd:ied  ,jsd:jed+shift,nz) )
    global1 = 0; global2 = 0
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx+shift
             global1(i,j,k) = 1.0_r8_kind + i*1.0e-3_r8_kind + j*1.0e-6_r8_kind + k*1.0e-9_r8_kind
          end do
       end do
       do j = 1, ny+shift
          do i = 1, nx
             global2(i,j,k) = 2.0_r8_kind + i*1.0e-3_r8_kind + j*1.0e-6_r8_kind + k*1.0e-9_r8_kind
          end do
       end do
    end do

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec,:) = global1(isc:iec+shift,jsc:jec,:)
    y (isc:iec,jsc:jec+shift,:) = global2(isc:iec,jsc:jec+shift,:)

    type = "nest grid CGRID_NE"
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE, name=trim(type))

    call compare_checksums( x(isd:ied+shift,jsd:jed,:), global1(isd:ied+shift,jsd:jed,:), trim(type)//' X' )
    call compare_checksums( y(isd:ied,jsd:jed+shift,:), global2(isd:ied,jsd:jed+shift,:), trim(type)//' Y' )

    deallocate(global1, global2, x, y)

  end subroutine test_nest_halo_update_r8

  subroutine get_nnest(domain, num_nest, tile_coarse, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, &
                       nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, &
                       je_coarse)
    type(domain2D), intent(inout) :: domain
    integer, intent(in)  :: num_nest, istart_coarse(:), iend_coarse(:), jstart_coarse(:), jend_coarse(:)
    integer, intent(in)  :: tile_coarse(:)
    integer, intent(out) :: nnest, is_coarse(:), ie_coarse(:), js_coarse(:), je_coarse(:)
    integer, intent(out) :: t_coarse(:), iadd_coarse(:), jadd_coarse(:), rotate_coarse(:)
    integer :: is, ie, js, je, tile, isg, ieg, jsg, jeg
    integer :: ncross, rotate, i1, i2, ntiles
    integer :: is2, ie2, js2, je2
    ntiles = 6
    call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
    nnest = 0

    do n = 1, num_nest
       is = istart_coarse(n); ie = iend_coarse(n)
       js = jstart_coarse(n); je = jend_coarse(n)
       tile = tile_coarse(n)
       ncross = 0
       rotate = 0
       do while(is >0 .and. js > 0)
          is2 = max(is,1); ie2 = min(ie,ieg)
          js2 = max(js,1); je2 = min(je,jeg)

          if(ie2 .GE. is2 .and. je2 .GE. js2) then
             nnest = nnest+1
             if(nnest > ntiles) call mpp_error(FATAL, "get_nnest: nnest > ntiles")
             select case(rotate)
             case(0)
                 is_coarse(nnest) = is2; ie_coarse(nnest) = ie2
                 js_coarse(nnest) = js2; je_coarse(nnest) = je2
             case(90)
                 is_coarse(nnest) = js2; ie_coarse(nnest) = je2
                 js_coarse(nnest) = ieg-ie2+1; je_coarse(nnest) = ieg-is2+1
             case(-90)
                is_coarse(nnest) = jeg-je2+1; ie_coarse(nnest) = jeg-js2+1
                js_coarse(nnest) = is2; je_coarse(nnest) = ie2
             end select
             if(jend_coarse(n) > jeg) then
                iadd_coarse(nnest) = 0
                jadd_coarse(nnest) = ncross*jeg
                if( je > jeg ) then
                   js = 1; je = je - jeg
                else
                   is = 0 ; js = 0
                endif
             else if(iend_coarse(n) > ieg) then
                iadd_coarse(nnest) = ncross*ieg
                jadd_coarse(nnest) = 0
                if(ie>ieg) then
                   is = 1; ie = ie - ieg
                else
                   is = 0; js = 0
                endif
             else
                iadd_coarse(nnest) = 0
                jadd_coarse(nnest) = 0
                is = 0; js = 0
             endif
             rotate_coarse(nnest) = rotate
             t_coarse(nnest) = tile
          else
             if(je > jeg) then
                js = js-jeg; je = je-jeg
             else if(ie > ieg) then
                is = is-ieg; ie = ie-ieg
             else
                call mpp_error(FATAL, "get_nnest: je <= jeg and ie <= ieg but (ie2<is2 .or. je2<js2")
             end if
          endif
          !--- figure out tile and rotation
          ncross = ncross+1
          if(jend_coarse(n) > jeg) then
             if(rotate == 0) then ! cross north edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                else  ! rotate 90 degree
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate + 90
                endif
             else if(rotate == 90) then  ! cross east edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate - 90
                else  ! tile 1 3 5
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                endif
             else
                call mpp_error(FATAL, "get_nnest: rotate should be 0 or 90 when je>jeg")
             endif
          else if(iend_coarse(n) > ieg) then
             if(rotate == 0) then ! cross east edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate - 90
                else  ! no rotation
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                endif
             else if(rotate == -90) then  ! cross east edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                else  ! tile 1 3 5
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate + 90
                endif
             else
                call mpp_error(FATAL, "get_nnest: rotate should be 0 or -90 when ie>ieg")
             endif
          endif
       enddo
    enddo

  end subroutine get_nnest


  subroutine get_nnest2(domain, num_nest, tile_coarse, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, &
                       nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, &
                       je_coarse)
    type(domain2D), intent(inout) :: domain
    integer, intent(in)  :: num_nest, istart_coarse(:), iend_coarse(:), jstart_coarse(:), jend_coarse(:)
    integer, intent(in)  :: tile_coarse(:)
    integer, intent(out) :: nnest, is_coarse(:), ie_coarse(:), js_coarse(:), je_coarse(:)
    integer, intent(out) :: t_coarse(:), iadd_coarse(:), jadd_coarse(:), rotate_coarse(:)
    integer :: is, ie, js, je, tile, isg, ieg, jsg, jeg
    integer :: ncross, rotate, i1, i2, ntiles
    integer :: is2, ie2, js2, je2
    ntiles = 6
    call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
    nnest = 0

    do n = 1, num_nest
       is = istart_coarse(n); ie = iend_coarse(n)
       js = jstart_coarse(n); je = jend_coarse(n)
       tile = tile_coarse(n)
       ncross = 0
       rotate = 0
       do while(is >0 .and. js > 0)
          is2 = max(is,1); ie2 = min(ie,ieg)
          js2 = max(js,1); je2 = min(je,jeg)
          if(ie2 .GE. is2 .and. je2 .GE. js2) then
             nnest = nnest+1
             if(nnest > ntiles) call mpp_error(FATAL, "get_nnest2: nnest > ntiles")
             if(jend_coarse(n) > jeg) then
                is_coarse(nnest) = is2; ie_coarse(nnest) = ie2
                js_coarse(nnest) = js2+ncross*jeg; je_coarse(nnest) = je2+ncross*jeg

                iadd_coarse(nnest) = 0
                jadd_coarse(nnest) = ncross*jeg
                if( je > jeg ) then
                   js = 1; je = je - jeg
                else
                   is = 0 ; js = 0
                endif
             else if(iend_coarse(n) > ieg) then
                is_coarse(nnest) = is2+ncross*ieg; ie_coarse(nnest) = ie2+ncross*ieg
                js_coarse(nnest) = js2           ; je_coarse(nnest) = je2
                iadd_coarse(nnest) = ncross*ieg
                jadd_coarse(nnest) = 0
                if(ie>ieg) then
                   is = 1; ie = ie - ieg
                else
                   is = 0; js = 0
                endif
             else
                 is_coarse(nnest) = is2; ie_coarse(nnest) = ie2
                 js_coarse(nnest) = js2; je_coarse(nnest) = je2
                iadd_coarse(nnest) = 0
                jadd_coarse(nnest) = 0
                is = 0; js = 0
             endif
             rotate_coarse(nnest) = rotate
             t_coarse(nnest) = tile
          else
             if(je > jeg) then
                js = js-jeg; je = je-jeg
             else if(ie > ieg) then
                is = is-ieg; ie = ie-ieg
             else
                call mpp_error(FATAL, "get_nnest2: je <= jeg and ie <= ieg but (ie2<is2 .or. je2<js2")
             end if
          endif
          !--- figure out tile and rotation
          ncross = ncross+1
          if(jend_coarse(n) > jeg) then
             if(rotate == 0) then ! cross north edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                else  ! rotate 90 degree
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate + 90
                endif
             else if(rotate == 90) then  ! cross east edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate - 90
                else  ! tile 1 3 5
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                endif
             else
                call mpp_error(FATAL, "get_nnest2: rotate should be 0 or 90 when je>jeg")
             endif
          else if(iend_coarse(n) > ieg) then
             if(rotate == 0) then ! cross east edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate - 90
                else  ! no rotation
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                endif
             else if(rotate == -90) then  ! cross east edge
                if(mod(tile,2) ==0) then ! tile 2 4 6
                   tile = tile + 1
                   if(tile>ntiles) tile=tile-ntiles
                else  ! tile 1 3 5
                   tile = tile + 2
                   if(tile>ntiles) tile=tile-ntiles
                   rotate = rotate + 90
                endif
             else
                call mpp_error(FATAL, "get_nnest2: rotate should be 0 or -90 when ie>ieg")
             endif
          endif
       enddo
    enddo
  end subroutine get_nnest2

!###############################################################################
  subroutine test_update_nest_domain_r8( type )
    character(len=*), intent(in) :: type
    logical                      :: cubic_grid
    logical                      :: is_fine_pe, is_coarse_pe
    integer                      :: n, i, j, k
    integer                      :: ntiles, npes_per_tile
    integer                      :: npes_fine, pos
    integer                      :: isc_coarse, iec_coarse, jsc_coarse, jec_coarse
    integer                      :: isd_coarse, ied_coarse, jsd_coarse, jed_coarse
    integer                      :: isd_fine, ied_fine, jsd_fine, jed_fine
    integer                      :: isc_fine, iec_fine, jsc_fine, jec_fine
    integer                      :: nx_fine, ny_fine, nx_coarse, ny_coarse
    integer                      :: nxc_fine, nyc_fine, nxc_coarse, nyc_coarse
    integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
    integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
    integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
    integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
    integer                      :: isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx
    integer                      :: ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx
    integer                      :: iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx
    integer                      :: isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx
    integer                      :: isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy
    integer                      :: ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy
    integer                      :: iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy
    integer                      :: isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy
    integer                      :: isw_f2, iew_f2, jsw_f2, jew_f2, isw_c2, iew_c2, jsw_c2, jew_c2, tile_w2
    integer                      :: ise_f2, iee_f2, jse_f2, jee_f2, ise_c2, iee_c2, jse_c2, jee_c2, tile_e2
    integer                      :: iss_f2, ies_f2, jss_f2, jes_f2, iss_c2, ies_c2, jss_c2, jes_c2, tile_s2
    integer                      :: isn_f2, ien_f2, jsn_f2, jen_f2, isn_c2, ien_c2, jsn_c2, jen_c2, tile_n2
    integer                      :: isw_fx2, iew_fx2, jsw_fx2, jew_fx2, isw_cx2, iew_cx2, jsw_cx2, jew_cx2, tile_wx2
    integer                      :: ise_fx2, iee_fx2, jse_fx2, jee_fx2, ise_cx2, iee_cx2, jse_cx2, jee_cx2, tile_ex2
    integer                      :: iss_fx2, ies_fx2, jss_fx2, jes_fx2, iss_cx2, ies_cx2, jss_cx2, jes_cx2, tile_sx2
    integer                      :: isn_fx2, ien_fx2, jsn_fx2, jen_fx2, isn_cx2, ien_cx2, jsn_cx2, jen_cx2, tile_nx2
    integer                      :: isw_fy2, iew_fy2, jsw_fy2, jew_fy2, isw_cy2, iew_cy2, jsw_cy2, jew_cy2, tile_wy2
    integer                      :: ise_fy2, iee_fy2, jse_fy2, jee_fy2, ise_cy2, iee_cy2, jse_cy2, jee_cy2, tile_ey2
    integer                      :: iss_fy2, ies_fy2, jss_fy2, jes_fy2, iss_cy2, ies_cy2, jss_cy2, jes_cy2, tile_sy2
    integer                      :: isn_fy2, ien_fy2, jsn_fy2, jen_fy2, isn_cy2, ien_cy2, jsn_cy2, jen_cy2, tile_ny2
    integer                      :: isw_f_T, iew_f_T, jsw_f_T, jew_f_T, isw_c_T, iew_c_T, jsw_c_T, jew_c_T
    integer                      :: ise_f_T, iee_f_T, jse_f_T, jee_f_T, ise_c_T, iee_c_T, jse_c_T, jee_c_T
    integer                      :: iss_f_T, ies_f_T, jss_f_T, jes_f_T, iss_c_T, ies_c_T, jss_c_T, jes_c_T
    integer                      :: isn_f_T, ien_f_T, jsn_f_T, jen_f_T, isn_c_T, ien_c_T, jsn_c_T, jen_c_T
    integer                      :: is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f
    integer                      :: is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx
    integer                      :: is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy
    integer                      :: tile, position, shift
    integer                      :: layout_fine(2), my_fine_id
    integer, allocatable         :: pelist(:), start_pos(:), end_pos(:)
    integer, allocatable         :: my_pelist_fine(:)
    integer, allocatable         :: pe_start(:), pe_end(:)
    integer, allocatable         :: layout2D(:,:), global_indices(:,:)
    real(kind=r8_kind), allocatable :: x(:,:,:), x1(:,:,:), x2(:,:,:)
    real(kind=r8_kind), allocatable :: y(:,:,:), y1(:,:,:), y2(:,:,:)
    real(kind=r8_kind), allocatable :: wbuffer(:,:,:), wbuffer2(:,:,:)
    real(kind=r8_kind), allocatable :: ebuffer(:,:,:), ebuffer2(:,:,:)
    real(kind=r8_kind), allocatable :: sbuffer(:,:,:), sbuffer2(:,:,:)
    real(kind=r8_kind), allocatable :: nbuffer(:,:,:), nbuffer2(:,:,:)
    real(kind=r8_kind), allocatable :: wbufferx(:,:,:), wbufferx2(:,:,:)
    real(kind=r8_kind), allocatable :: ebufferx(:,:,:), ebufferx2(:,:,:)
    real(kind=r8_kind), allocatable :: sbufferx(:,:,:), sbufferx2(:,:,:)
    real(kind=r8_kind), allocatable :: nbufferx(:,:,:), nbufferx2(:,:,:)
    real(kind=r8_kind), allocatable :: wbuffery(:,:,:), wbuffery2(:,:,:)
    real(kind=r8_kind), allocatable :: ebuffery(:,:,:), ebuffery2(:,:,:)
    real(kind=r8_kind), allocatable :: sbuffery(:,:,:), sbuffery2(:,:,:)
    real(kind=r8_kind), allocatable :: nbuffery(:,:,:), nbuffery2(:,:,:)
    integer                      :: x_refine(num_nest), y_refine(num_nest)
    integer                      :: istart_fine(num_nest), iend_fine(num_nest)
    integer                      :: jstart_fine(num_nest), jend_fine(num_nest)
    integer                      :: iend_coarse(num_nest), jend_coarse(num_nest)
    integer                      :: is_fine(6*num_nest), ie_fine(6*num_nest)
    integer                      :: js_fine(6*num_nest), je_fine(6*num_nest)
    integer                      :: is_coarse(6*num_nest), ie_coarse(6*num_nest)
    integer                      :: js_coarse(6*num_nest), je_coarse(6*num_nest)
    integer                      :: t_coarse(6*num_nest), rotate_coarse(6*num_nest)
    integer                      :: iadd_coarse(6*num_nest), jadd_coarse(6*num_nest)
    integer                      :: nnest
    character(len=128)           :: type2
    character(len=32)            :: text, pelist_name
    type(domain2d)               :: domain
    type(domain2d), pointer      :: domain_coarse=>NULL()
    type(domain2d), pointer      :: domain_fine=>NULL()
    type(nest_domain_type)       :: nest_domain
    logical                      :: x_cyclic, y_cyclic
    integer                      :: my_tile_id(1), my_num_nest
    integer, dimension(num_nest) :: my_tile_coarse, my_tile_fine, my_istart_coarse, my_iend_coarse
    integer, dimension(num_nest) :: my_jstart_coarse, my_jend_coarse
    integer                      :: ntiles_nest_top, npes_nest_top, num_nest_level, my_npes, l
    integer                      :: npes_my_fine, npes_my_level
    integer, allocatable         :: my_pelist(:)
    integer, dimension(1)        :: dummy

    x_cyclic = .false.
    y_cyclic = .false.
    if(cyclic_nest(1) == 'X') then
       x_cyclic = .true.
    else if(cyclic_nest(1) == 'Y') then
       y_cyclic = .true.
    endif

    istart_fine = 0; iend_fine = -1
    jstart_fine = 0; jend_fine = -1
    iend_coarse = -1; jend_coarse = -1
    is_fine = 0;  ie_fine = -1
    js_fine = 0;  je_fine = -1
    is_coarse = 0;  ie_coarse = -1
    js_coarse = 0;  je_coarse = -1
    t_coarse = 0; rotate_coarse = -1;
    iadd_coarse = 0; jadd_coarse = 0

    select case(type)
      case ( 'Cubic-Grid' )
        if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic is zero, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        nx = nx_cubic
        ny = ny_cubic
        ntiles_nest_top = 6
        cubic_grid = .true.

      case ( 'Cubed-sphere, single face' )
        if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic is zero, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
        endif
        nx = nx_cubic
        ny = ny_cubic
        ntiles_nest_top = 1
        cubic_grid = .false.

    case ( 'Regional' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Regional grid mosaic, nx_cubic is zero, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Regional grid mosaic, nx_cubic does not equal ny_cubic, '&
                  & //'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       nx = nx_cubic
       ny = ny_cubic
       ntiles_nest_top = 1
       cubic_grid = .false.
    case default
       call mpp_error(FATAL, 'test_update_nest_domain: no such test: '//type)
    end select

    if(ntiles_nest_all > MAX_NTILE) call mpp_error(FATAL, 'test_update_nest_domain: ntiles_nest_all > MAX_NTILE')
    if(ntiles_nest_top .GE. ntiles_nest_all) call mpp_error(FATAL, &
                                            & 'test_update_nest_domain: ntiles_nest_top .GE. ntile_nest_all')
    if(ntiles_nest_all .NE. ntiles_nest_top + num_nest) call mpp_error(FATAL, &
             'test_update_nest_domain: ntiles_nest_all .NE. ntiles_nest_top + num_nest')
    !--- for the ntiles_nest_top, number of processors should be same
    do n = 1, ntiles_nest_all
       if(npes_nest_tile(n) .LE. 0) call mpp_error(FATAL, &
            'test_update_nest_domain: npes_nest_tile is not properly set')
    enddo
    do n = 2, ntiles_nest_top
       if(npes_nest_tile(n) .NE. npes_nest_tile(n-1)) call mpp_error(FATAL, &
            'test_update_nest_domain: each tile of top mosaic grid should use same number of MPI ranks')
    enddo
    npes_nest_top = ntiles_nest_top * npes_nest_tile(1)

    npes = mpp_npes()

    !--- make sure sum(npes_nest_tile) == npes
    if(sum(npes_nest_tile(1:ntiles_nest_all)) .NE. npes ) &
         call mpp_error(FATAL, "test_mpp_domains: sum(npes_nest_tile) .NE. npes")

    !--- make sure tile_fine are monotonically increasing and equal to ntiles_nest_top + nest number
    do n = 1, num_nest
       if(tile_fine(n) .NE. ntiles_nest_top+n) call mpp_error(FATAL, &
           "test_mpp_domains: tile_fine(n) .NE. ntiles_nest_top+n")
    enddo

    !---make sure nest_level is setup properly
    if(nest_level(1) .NE. 1) call mpp_error(FATAL, "test_mpp_domains: nest_level(1) .NE. 1")
    do n = 2, num_nest
       if(nest_level(n) > nest_level(n-1)+1)call mpp_error(FATAL,"test_mpp_domains: nest_level(n) > nest_level(n-1)+1")
       if(nest_level(n) < nest_level(n-1) ) call mpp_error(FATAL, "test_mpp_domains: nest_level(n) < nest_level(n-1)")
    enddo
    num_nest_level = nest_level(num_nest)

    allocate(pelist(npes))
    call mpp_get_current_pelist(pelist)

    !--- compute iend_coarse and jend_coarse
    do n = 1, num_nest
       iend_coarse(n) = istart_coarse(n) + icount_coarse(n) - 1
       jend_coarse(n) = jstart_coarse(n) + jcount_coarse(n) - 1
       istart_fine(n) = 1; iend_fine(n) = icount_coarse(n)*refine_ratio(n)
       jstart_fine(n) = 1; jend_fine(n) = jcount_coarse(n)*refine_ratio(n)
    enddo

    !--- first define the top level grid mosaic domain.

    !--- setup pelist for top level
    allocate(my_pelist(npes_nest_top))
    do n = 1, npes_nest_top
       my_pelist(n) = pelist(n)
    enddo
    call mpp_declare_pelist(my_pelist)
    if(ANY(my_pelist==mpp_pe())) then
       call mpp_set_current_pelist(my_pelist)
       allocate(layout2D(2,ntiles_nest_top), global_indices(4,ntiles_nest_top), pe_start(ntiles_nest_top), &
              & pe_end(ntiles_nest_top) )
       npes_per_tile = npes_nest_tile(1)


       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles_nest_top
          global_indices(:,n) = (/1,nx,1,ny/)
          if (ANY(layout_cubic == 0)) then
             layout2D(:,n)         = layout
          else
             layout2D(:,n)         = layout_cubic(:)
           endif
       end do
       do n = 1, ntiles_nest_top
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do

       if( cubic_grid ) then
         call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                                   global_indices, layout2D, pe_start, pe_end )
       else
          call mpp_define_domains(global_indices(:,ntiles_nest_top), layout2D(:,ntiles_nest_top), domain, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.true., name=trim(type)//' top level grid', tile_id=1  )
       endif
       call mpp_get_compute_domain(domain, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
       call mpp_get_data_domain(domain, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
       deallocate(layout2D, global_indices, pe_start, pe_end )
    endif

    call mpp_set_current_pelist()
    deallocate(my_pelist)
    !--- define domain for all the nest region.
    pos = npes_nest_top
    do n = 1, num_nest
       my_npes = npes_nest_tile(tile_fine(n))
       allocate(my_pelist(my_npes))
       my_pelist(:) = pelist(pos+1:pos+my_npes)
       call mpp_declare_pelist(my_pelist)
       if(ANY(my_pelist==mpp_pe())) then
          call mpp_set_current_pelist(my_pelist)
          nx_fine = iend_fine(n) - istart_fine(n) + 1
          ny_fine = jend_fine(n) - jstart_fine(n) + 1
          if (ANY(layout_nest(:,n) == 0)) then
             call mpp_define_layout( (/1,nx_fine,1,ny_fine/), my_npes, layout )
          else
             layout(:)         = layout_nest(:,n)
          endif
          call mpp_define_domains((/1,nx_fine,1,ny_fine/), layout, domain, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.true., name=trim(type)//' fine grid', tile_id = tile_fine(n) )
          call mpp_get_compute_domain(domain, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain, isd_fine, ied_fine, jsd_fine, jed_fine)
          !--- test halo update for nested region.
          call test_nest_halo_update_r8(domain)
       endif
       pos = pos+my_npes
       deallocate(my_pelist)
       call mpp_set_current_pelist()
    enddo

    !--- reset to the global pelist
    call mpp_set_current_pelist()

    x_refine(:) = refine_ratio(1:num_nest)
    y_refine(:) = refine_ratio(1:num_nest)

    call mpp_define_nest_domains(nest_domain, domain, num_nest, nest_level(1:num_nest), tile_fine(1:num_nest), &
             tile_coarse(1:num_nest), istart_coarse(1:num_nest), icount_coarse(1:num_nest), jstart_coarse(1:num_nest),&
             jcount_coarse(1:num_nest), npes_nest_tile(1:ntiles_nest_all), &
             x_refine(1:num_nest), y_refine(1:num_nest), extra_halo=extra_halo, name="nest_domain")

    !--- loop over nest level
    do l = 1, num_nest_level
       npes_my_level = mpp_get_nest_npes(nest_domain, l)
       npes_my_fine = mpp_get_nest_fine_npes(nest_domain,l)
       allocate(my_pelist(npes_my_level))
       allocate(my_pelist_fine(npes_my_fine))
       call mpp_get_nest_pelist(nest_domain, l, my_pelist)

       call mpp_declare_pelist(my_pelist(:))
       write(type2, '(a,I2)')trim(type)//" nest_level = ",l
       if(ANY(my_pelist(:)==mpp_pe())) then

          call mpp_get_nest_fine_pelist(nest_domain, l, my_pelist_fine)

          call mpp_set_current_pelist(my_pelist)
          my_tile_id = mpp_get_tile_id(domain)
          domain_coarse => mpp_get_nest_coarse_domain(nest_domain, nest_level=l)
          domain_fine => mpp_get_nest_fine_domain(nest_domain, nest_level=l)
          is_fine_pe = mpp_is_nest_fine(nest_domain, l)
          is_coarse_pe = mpp_is_nest_coarse(nest_domain, l)
          if(is_fine_pe .eqv. is_coarse_pe) call mpp_error(FATAL, "test_mpp_domains: is_fine_pe .eqv. is_coarse_pe")
          my_num_nest = 0
          my_fine_id = 0
          do n = 1, num_nest
             if(nest_level(n)==l) then
                my_num_nest = my_num_nest+1
                my_tile_coarse(my_num_nest) = tile_coarse(n)
                my_tile_fine(my_num_nest) = tile_fine(n)
                my_istart_coarse(my_num_nest) = istart_coarse(n)
                my_iend_coarse(my_num_nest) = iend_coarse(n)
                my_jstart_coarse(my_num_nest) = jstart_coarse(n)
                my_jend_coarse(my_num_nest) = jend_coarse(n)
                if(my_tile_id(1) == tile_fine(n)) my_fine_id = n
             endif
          enddo
          !--- each nest region might be over multiple face of cubic sphere grid.
          !---Get the number of nest region with consideration of face.
          call get_nnest(domain_coarse, my_num_nest, my_tile_coarse, my_istart_coarse, my_iend_coarse, &
               my_jstart_coarse, my_jend_coarse, nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, &
               is_coarse, ie_coarse, js_coarse, je_coarse)

          !---------------------------------------------------------------------------
          !
          !                    fine to coarse scalar field, limit to position=CENTER.
          !
          !---------------------------------------------------------------------------
          if(is_fine_pe) then
             call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine)
             call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine)
          endif

          if(is_coarse_pe) then
             call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
             call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          endif

          if(is_fine_pe) then
             call mpp_get_F2C_index(nest_domain, is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f, l, position=CENTER)
             allocate(x(is_c:ie_c, js_c:je_c, nz))
             x = 0
             do k = 1, nz
                do j = js_c, je_c
                   do i = is_c, ie_c
                      x(i,j,k) = i*1.D+6 + j*1.D+3 + k + 0.001_r8_kind
                   enddo
                enddo
             enddo
          else
             allocate(x1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse, nz))
             allocate(x2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse, nz))
             x1 = 0
             x2 = 0
             tile = my_tile_id(1)

             do k = 1, nz
                do j = jsc_coarse, jec_coarse
                   do i = isc_coarse, iec_coarse
                      x1(i,j,k) = i*1.D+6 + j*1.D+3 + k + 0.002_r8_kind
                   enddo
                enddo
             enddo
             x2 = x1
          endif

          if(is_coarse_pe) then
             do n = 1, nnest
                is_c = max(is_coarse(n), isc_coarse)
                ie_c = min(ie_coarse(n),   iec_coarse)
                js_c = max(js_coarse(n), jsc_coarse)
                je_c = min(je_coarse(n),   jec_coarse)
                if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c .GE. js_c ) then
                   call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                        is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, 0, 0.001_r8_kind, &
                        0.001_r8_kind, 1, 1, .false., .false., iend_coarse(1), jend_coarse(1) )
                endif
             enddo
          endif

          call mpp_update_nest_coarse(x, nest_domain, x1, nest_level=l, position=CENTER)

          !--- compare with assumed value
          if( is_coarse_pe) then
             !! initial failure point
             call compare_checksums(x1, x2, trim(type2)//' fine to coarse scalar')
          endif
          if(allocated(x))       deallocate(x)
          if(allocated(x1))      deallocate(x1)
          if(allocated(x2))      deallocate(x2)

       !---------------------------------------------------------------------------
       !
       !                    fine to coarse CGRID scalar pair update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_fine_pe) then
          call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, l, position=EAST)
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l,position=NORTH)
          allocate(x(is_cx:ie_cx, js_cx:je_cx, nz))
          allocate(y(is_cy:ie_cy, js_cy:je_cy, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = js_cx, je_cx
                do i = is_cx, ie_cx
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + dble(k) + 1.0D-6
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = js_cy, je_cy
                do i = is_cy, ie_cy
                   y(i,j,k) = i*1.D+6 + j*1.D+3 + dble(k) + 2.0D-6
                enddo
             enddo
          enddo
          if(x_cyclic) then
             if(ie_cx == iend_coarse(1)+1) then
                i = ie_cx
                do k = 1, nz
                   do j = js_cx, je_cx
                      x(i,j,k) = istart_coarse(1)*1.D+6 + j*1.D+3 + dble(k) + 1.0D-6
                   enddo
                enddo
             endif
          endif
          if(y_cyclic) then
             if(je_cx == jend_coarse(1)+1) then
                j = je_cx
                do k = 1, nz
                   do i = is_cx, ie_cx
                      y(i,j,k) = i*1.D+6 + jstart_coarse(1)*1.D+3 + dble(k) + 1.0D-6
                   enddo
                enddo
             endif
          endif
       else
          allocate(x1(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(x2(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(y2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x1 = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x1(i,j,k) = i*1.D+6 + j*1.D+3 + k + 0.001_r8_kind
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y1(i,j,k) = i*1.D+6 + j*1.D+3 + k + 0.002_r8_kind
                enddo
             enddo
          enddo
          x2 = x1
          y2 = y1
       endif


       if(is_coarse_pe) then
          do n = 1, nnest
             is_c = max(is_coarse(n), isc_coarse)
             ie_c = min(ie_coarse(n),   iec_coarse)
             js_c = max(js_coarse(n), jsc_coarse)
             je_c = min(je_coarse(n),   jec_coarse)
             if( tile == t_coarse(n) .AND. ie_c+shift .GE. is_c .AND. je_c .GE. js_c ) then
                call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0,  &
                     real(1.0D-6,r8_kind), real(2.0D-6,r8_kind), 1, 1, &
                     x_cyclic, .false., iend_coarse(1)+1, jend_coarse(1)+1)
             endif
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift,  &
                     real(2.0D-6, r8_kind),real(1.0D-6,r8_kind), 1, 1, &
                     .false., y_cyclic, iend_coarse(1)+1, jend_coarse(1)+1)
             endif
          enddo
       endif

       call mpp_update_nest_coarse(x, y, nest_domain, x1, y1, nest_level=l, gridtype=CGRID_NE, flags=SCALAR_PAIR)

       !--- compare with assumed value
       if( is_coarse_pe) then
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Scalar_pair X')
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Scalar_pair Y')
       endif
       if(allocated(x))       deallocate(x)
       if(allocated(x1))      deallocate(x1)
       if(allocated(x2))      deallocate(x2)
       if(allocated(y))       deallocate(y)
       if(allocated(y1))      deallocate(y1)
       if(allocated(y2))      deallocate(y2)

       !---------------------------------------------------------------------------
       !
       !                    fine to coarse CGRID vector update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_fine_pe) then
          call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, l, position=EAST)
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l,position=NORTH)
          allocate(x(is_cx:ie_cx, js_cx:je_cx, nz))
          allocate(y(is_cy:ie_cy, js_cy:je_cy, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = js_cx, je_cx
                do i = is_cx, ie_cx
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + k + 1.0D-6
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = js_cy, je_cy
                do i = is_cy, ie_cy
                   y(i,j,k) = i*1.D+6 + j*1.D+3 + k + 2.0D-6
                enddo
             enddo
          enddo
       else
          allocate(x1(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(x2(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(y2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x1 = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x1(i,j,k) = i*1.D+6 + j*1.D+3 + k + 0.001_r8_kind
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y1(i,j,k) = i*1.D+6 + j*1.D+3 + k + 0.001_r8_kind
                enddo
             enddo
          enddo
          x2 = x1
          y2 = y1
       endif


       if(is_coarse_pe) then
          do n = 1, nnest
             is_c = max(is_coarse(n), isc_coarse)
             ie_c = min(ie_coarse(n),   iec_coarse)
             js_c = max(js_coarse(n), jsc_coarse)
             je_c = min(je_coarse(n),   jec_coarse)
             if( tile == t_coarse(n) .AND. ie_c+shift .GE. is_c .AND. je_c .GE. js_c ) then
                call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0, &
                     real(1.0D-6,r8_kind), real(2.0D-6,r8_kind), 1, -1, &
                     x_cyclic, .false., iend_coarse(1)+1, jend_coarse(1)+1)
             endif
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift,  &
                     real(2.0D-6,r8_kind), real(1.0D-6,r8_kind), -1, 1, &
                     .false., y_cyclic, iend_coarse(1)+1, jend_coarse(1)+1)
             endif
          enddo
       endif

       call mpp_update_nest_coarse(x, y, nest_domain, x1, y1, nest_level=l, gridtype=CGRID_NE)

       !--- compare with assumed value
       if( is_coarse_pe) then
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Vector X')
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Vector Y')
       endif
       if(allocated(x))       deallocate(x)
       if(allocated(x1))      deallocate(x1)
       if(allocated(x2))      deallocate(x2)
       if(allocated(y))       deallocate(y)
       if(allocated(y1))      deallocate(y1)
       if(allocated(y2))      deallocate(y2)

       !---------------------------------------------------------------------------
       !
       !                    fine to coarse DGRID vector update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_fine_pe) then
          call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, l,position=NORTH)
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l, position=EAST)
          allocate(x(is_cx:ie_cx, js_cx:je_cx, nz))
          allocate(y(is_cy:ie_cy, js_cy:je_cy, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = js_cx, je_cx
                do i = is_cx, ie_cx
                   x(i,j,k) = i*1.d+6 + j*1.d+3 + k + 1.0d-6
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = js_cy, je_cy
                do i = is_cy, ie_cy
                   y(i,j,k) = i*1.d+6 + j*1.d+3 + k + 2.0d-6
                enddo
             enddo
          enddo
       else
          allocate(x1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(x2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(y1(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y2(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          x1 = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   x1(i,j,k) = i*1.d+6 + j*1.d+3 + dble(k) + 0.001_r8_kind
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   y1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.002
                enddo
             enddo
          enddo
          x2 = x1
          y2 = y1
       endif


       if(is_coarse_pe) then
          do n = 1, nnest
             is_c = max(is_coarse(n), isc_coarse)
             ie_c = min(ie_coarse(n),   iec_coarse)
             js_c = max(js_coarse(n), jsc_coarse)
             je_c = min(je_coarse(n),   jec_coarse)
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift,  &
                     real(1.0D-6,r8_kind), real(2.0D-6,r8_kind), 1, -1, &
                     .false., y_cyclic, iend_coarse(1), jend_coarse(1) )
             endif
             if( tile == t_coarse(n) .AND. ie_c+shift .GE. is_c .AND. je_c .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0,  &
                     real(2.0D-6,r8_kind), real(1.0D-6,r8_kind), -1, 1, &
                     x_cyclic, .false., iend_coarse(1), jend_coarse(1))
             endif
          enddo
       endif

       call mpp_update_nest_coarse(x, y, nest_domain, x1, y1, nest_level=l, gridtype=DGRID_NE)

       !--- compare with assumed value
       if( is_coarse_pe) then
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer DGRID Vector X')
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer DGRID Vector Y')
       endif
       if(allocated(x))       deallocate(x)
       if(allocated(x1))      deallocate(x1)
       if(allocated(x2))      deallocate(x2)
       if(allocated(y))       deallocate(y)
       if(allocated(y1))      deallocate(y1)
       if(allocated(y2))      deallocate(y2)

       !---------------------------------------------------------------------------
       !
       !                 Coarse to Fine scalar field, position = CENTER
       !
       !---------------------------------------------------------------------------
       !--- first check the index is correct or not
       !--- The index from nest domain
       call mpp_get_C2F_index(nest_domain, isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c, WEST, l)
       call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, EAST, l)
       call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, SOUTH, l)
       call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, NORTH, l)

       if(is_fine_pe) then
          call mpp_get_compute_domain(domain, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain, isd_fine, ied_fine, jsd_fine, jed_fine)

          !-- The assumed index
          isw_f2 = 0; iew_f2 = -1; jsw_f2 = 0; jew_f2 = -1
          isw_c2 = 0; iew_c2 = -1; jsw_c2 = 0; jew_c2 = -1
          ise_f2 = 0; iee_f2 = -1; jse_f2 = 0; jee_f2 = -1
          ise_c2 = 0; iee_c2 = -1; jse_c2 = 0; jee_c2 = -1
          iss_f2 = 0; ies_f2 = -1; jss_f2 = 0; jes_f2 = -1
          iss_c2 = 0; ies_c2 = -1; jss_c2 = 0; jes_c2 = -1
          isn_f2 = 0; ien_f2 = -1; jsn_f2 = 0; jen_f2 = -1
          isn_c2 = 0; ien_c2 = -1; jsn_c2 = 0; jen_c2 = -1

          !--- west
          if( isc_fine == 1 ) then
             isw_f2 = isd_fine; iew_f2 = isc_fine - 1
             jsw_f2 = jsd_fine; jew_f2 = jed_fine
             isw_c2 = istart_coarse(my_fine_id)-whalo
             iew_c2 = istart_coarse(my_fine_id)
             jsw_c2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_c2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
          endif
          !--- east
          if( iec_fine == nx_fine ) then
             ise_f2 = iec_fine+1; iee_f2 = ied_fine
             jse_f2 = jsd_fine;   jee_f2 = jed_fine
             ise_c2 = iend_coarse(my_fine_id)
             iee_c2 = iend_coarse(my_fine_id)+ehalo
             jse_c2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_c2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_f2 = isd_fine; ies_f2 = ied_fine
             jss_f2 = jsd_fine; jes_f2 = jsc_fine - 1
             iss_c2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_c2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jss_c2 = jstart_coarse(my_fine_id)-shalo
             jes_c2 = jstart_coarse(my_fine_id)
          endif
          !--- north
          if( jec_fine == ny_fine ) then
             isn_f2 = isd_fine;  ien_f2 = ied_fine
             jsn_f2 = jec_fine+1; jen_f2 = jed_fine
             isn_c2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_c2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jsn_c2 = jend_coarse(my_fine_id)
             jen_c2 = jend_coarse(my_fine_id)+nhalo
          endif

          if( isw_f .NE. isw_f2 .OR. iew_f .NE. iew_f2 .OR. jsw_f .NE. jsw_f2 .OR. jew_f .NE. jew_f2 .OR. &
               isw_c .NE. isw_c2 .OR. iew_c .NE. iew_c2 .OR. jsw_c .NE. jsw_c2 .OR. jew_c .NE. jew_c2 ) then
             write(5000+mpp_pe(),*) "west buffer fine index = ", isw_f, iew_f, jsw_f, jew_f
             write(5000+mpp_pe(),*) "west buffer fine index2 = ", isw_f2, iew_f2, jsw_f2, jew_f2
             write(5000+mpp_pe(),*) "west buffer coarse index = ", isw_c, iew_c, jsw_c, jew_c
             write(5000+mpp_pe(),*) "west buffer coarse index2 = ", isw_c2, iew_c2, jsw_c2, jew_c2
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine scalar")
          endif
          if( ise_f .NE. ise_f2 .OR. iee_f .NE. iee_f2 .OR. jse_f .NE. jse_f2 .OR. jee_f .NE. jee_f2 .OR. &
               ise_c .NE. ise_c2 .OR. iee_c .NE. iee_c2 .OR. jse_c .NE. jse_c2 .OR. jee_c .NE. jee_c2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine scalar")
          endif
          if( iss_f .NE. iss_f2 .OR. ies_f .NE. ies_f2 .OR. jss_f .NE. jss_f2 .OR. jes_f .NE. jes_f2 .OR. &
               iss_c .NE. iss_c2 .OR. ies_c .NE. ies_c2 .OR. jss_c .NE. jss_c2 .OR. jes_c .NE. jes_c2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine scalar")
          endif
          if( isn_f .NE. isn_f2 .OR. ien_f .NE. ien_f2 .OR. jsn_f .NE. jsn_f2 .OR. jen_f .NE. jen_f2 .OR. &
               isn_c .NE. isn_c2 .OR. ien_c .NE. ien_c2 .OR. jsn_c .NE. jsn_c2 .OR. jen_c .NE. jen_c2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine scalar")
          endif
       endif

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse, jsd_coarse:jed_coarse, nz))
          x = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse
                   x(i,j,k) = tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine, jsd_fine:jed_fine, nz))
          x = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + k
                enddo
             enddo
          enddo
       endif

       if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
          allocate(wbuffer(isw_c:iew_c, jsw_c:jew_c,nz))
          allocate(wbuffer2(isw_c:iew_c, jsw_c:jew_c,nz))
       else
          allocate(wbuffer(1,1,1))
          allocate(wbuffer2(1,1,1))
       endif

       if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
          allocate(ebuffer(ise_c:iee_c, jse_c:jee_c,nz))
          allocate(ebuffer2(ise_c:iee_c, jse_c:jee_c,nz))
       else
          allocate(ebuffer(1,1,1))
          allocate(ebuffer2(1,1,1))
       endif

       if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
          allocate(sbuffer(iss_c:ies_c, jss_c:jes_c,nz))
          allocate(sbuffer2(iss_c:ies_c, jss_c:jes_c,nz))
       else
          allocate(sbuffer(1,1,1))
          allocate(sbuffer2(1,1,1))
       endif

       if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
          allocate(nbuffer(isn_c:ien_c, jsn_c:jen_c,nz))
          allocate(nbuffer2(isn_c:ien_c, jsn_c:jen_c,nz))
       else
          allocate(nbuffer(1,1,1))
          allocate(nbuffer2(1,1,1))
       endif
       ebuffer = 0; ebuffer2 = 0
       wbuffer = 0; wbuffer2 = 0
       sbuffer = 0; sbuffer2 = 0
       nbuffer = 0; nbuffer2 = 0

       call mpp_update_nest_fine(x, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=l)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), &
                  (/jsw_c/), (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, &
                  ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbuffer2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r8_kind, 0.0_r8_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(wbuffer, wbuffer2, trim(type2)//' west buffer coarse to fine scalar')

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), &
                  (/jss_c/), (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, &
                  ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbuffer2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r8_kind, 0.0_r8_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(sbuffer, sbuffer2, trim(type2)//' south buffer coarse to fine scalar')

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), &
                  (/jse_c/), (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, &
                   ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(ebuffer2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r8_kind, 0.0_r8_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(ebuffer, ebuffer2, trim(type2)//' east buffer coarse to fine scalar')

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), &
                  (/jsn_c/), (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, &
                  is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbuffer2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r8_kind, 0.0_r8_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(nbuffer, nbuffer2, trim(type2)//' north buffer coarse to fine scalar')
       endif
       if(is_fine_pe) then
          deallocate(wbuffer, ebuffer, sbuffer, nbuffer)
          deallocate(wbuffer2, ebuffer2, sbuffer2, nbuffer2)
       endif
       deallocate(x)

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine BGRID scalar pair update
       !
       !---------------------------------------------------------------------------
       shift = 1
       !--- first check the index is correct or not
       if(is_fine_pe) then
          !--- The index from nest domain
          call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine)
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, &
                               & WEST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, &
                               & EAST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, &
                               & SOUTH, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, &
                               & NORTH, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, &
                               & WEST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, &
                               & EAST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, &
                               & SOUTH, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, &
                               & NORTH, l, position=CORNER)

          !-- The assumed index
          isw_fx2 = 0; iew_fx2 = -1; jsw_fx2 = 0; jew_fx2 = -1
          isw_cx2 = 0; iew_cx2 = -1; jsw_cx2 = 0; jew_cx2 = -1
          ise_fx2 = 0; iee_fx2 = -1; jse_fx2 = 0; jee_fx2 = -1
          ise_cx2 = 0; iee_cx2 = -1; jse_cx2 = 0; jee_cx2 = -1
          iss_fx2 = 0; ies_fx2 = -1; jss_fx2 = 0; jes_fx2 = -1
          iss_cx2 = 0; ies_cx2 = -1; jss_cx2 = 0; jes_cx2 = -1
          isn_fx2 = 0; ien_fx2 = -1; jsn_fx2 = 0; jen_fx2 = -1
          isn_cx2 = 0; ien_cx2 = -1; jsn_cx2 = 0; jen_cx2 = -1
          isw_fy2 = 0; iew_fy2 = -1; jsw_fy2 = 0; jew_fy2 = -1
          isw_cy2 = 0; iew_cy2 = -1; jsw_cy2 = 0; jew_cy2 = -1
          ise_fy2 = 0; iee_fy2 = -1; jse_fy2 = 0; jee_fy2 = -1
          ise_cy2 = 0; iee_cy2 = -1; jse_cy2 = 0; jee_cy2 = -1
          iss_fy2 = 0; ies_fy2 = -1; jss_fy2 = 0; jes_fy2 = -1
          iss_cy2 = 0; ies_cy2 = -1; jss_cy2 = 0; jes_cy2 = -1
          isn_fy2 = 0; ien_fy2 = -1; jsn_fy2 = 0; jen_fy2 = -1
          isn_cy2 = 0; ien_cy2 = -1; jsn_cy2 = 0; jen_cy2 = -1

          !--- west
          if( isc_fine == 1 ) then
             isw_fx2 = isd_fine
             iew_fx2 = isc_fine - 1
             jsw_fx2 = jsd_fine
             jew_fx2 = jed_fine + shift
             isw_cx2 = istart_coarse(my_fine_id)-whalo
             iew_cx2 = istart_coarse(my_fine_id)
             jsw_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
             isw_fy2 = isd_fine
             iew_fy2 = isc_fine - 1
             jsw_fy2 = jsd_fine
             jew_fy2 = jed_fine + shift
             isw_cy2 = istart_coarse(my_fine_id)-whalo
             iew_cy2 = istart_coarse(my_fine_id)
             jsw_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- east
          if( iec_fine == nx_fine ) then
             ise_fx2 = iec_fine+1+shift
             iee_fx2 = ied_fine + shift
             jse_fx2 = jsd_fine
             jee_fx2 = jed_fine + shift
             ise_cx2 = iend_coarse(my_fine_id)+shift
             iee_cx2 = iend_coarse(my_fine_id)+ehalo+shift
             jse_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
             ise_fy2 = iec_fine+1 + shift
             iee_fy2 = ied_fine + shift
             jse_fy2 = jsd_fine
             jee_fy2 = jed_fine + shift
             ise_cy2 = iend_coarse(my_fine_id) + shift
             iee_cy2 = iend_coarse(my_fine_id)+ehalo + shift
             jse_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_fx2 = isd_fine
             ies_fx2 = ied_fine + shift
             jss_fx2 = jsd_fine
             jes_fx2 = jsc_fine - 1
             iss_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jss_cx2 = jstart_coarse(my_fine_id)-shalo
             jes_cx2 = jstart_coarse(my_fine_id)
             iss_fy2 = isd_fine
             ies_fy2 = ied_fine + shift
             jss_fy2 = jsd_fine
             jes_fy2 = jsc_fine - 1
             iss_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jss_cy2 = jstart_coarse(my_fine_id)-shalo
             jes_cy2 = jstart_coarse(my_fine_id)
          endif
          !--- north
          if( jec_fine == ny_fine ) then
             isn_fx2 = isd_fine
             ien_fx2 = ied_fine + shift
             jsn_fx2 = jec_fine+1 + shift
             jen_fx2 = jed_fine + shift
             isn_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jsn_cx2 = jend_coarse(my_fine_id) + shift
             jen_cx2 = jend_coarse(my_fine_id)+nhalo + shift
             isn_fy2 = isd_fine
             ien_fy2 = ied_fine + shift
             jsn_fy2 = jec_fine+1 + shift
             jen_fy2 = jed_fine + shift
             isn_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jsn_cy2 = jend_coarse(my_fine_id) + shift
             jen_cy2 = jend_coarse(my_fine_id)+nhalo + shift
          endif

          if( isw_fx .NE. isw_fx2 .OR. iew_fx .NE. iew_fx2 .OR. jsw_fx .NE. jsw_fx2 .OR. jew_fx .NE. jew_fx2 .OR. &
               isw_cx .NE. isw_cx2 .OR. iew_cx .NE. iew_cx2 .OR. jsw_cx .NE. jsw_cx2 .OR. jew_cx .NE. jew_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine BGRID X")
          endif
          if( ise_fx .NE. ise_fx2 .OR. iee_fx .NE. iee_fx2 .OR. jse_fx .NE. jse_fx2 .OR. jee_fx .NE. jee_fx2 .OR. &
               ise_cx .NE. ise_cx2 .OR. iee_cx .NE. iee_cx2 .OR. jse_cx .NE. jse_cx2 .OR. jee_cx .NE. jee_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine BGRID X")
          endif
          if( iss_fx .NE. iss_fx2 .OR. ies_fx .NE. ies_fx2 .OR. jss_fx .NE. jss_fx2 .OR. jes_fx .NE. jes_fx2 .OR. &
               iss_cx .NE. iss_cx2 .OR. ies_cx .NE. ies_cx2 .OR. jss_cx .NE. jss_cx2 .OR. jes_cx .NE. jes_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine BGRID X")
          endif
          if( isn_fx .NE. isn_fx2 .OR. ien_fx .NE. ien_fx2 .OR. jsn_fx .NE. jsn_fx2 .OR. jen_fx .NE. jen_fx2 .OR. &
               isn_cx .NE. isn_cx2 .OR. ien_cx .NE. ien_cx2 .OR. jsn_cx .NE. jsn_cx2 .OR. jen_cx .NE. jen_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine BGRID X")
          endif

          if( isw_fy .NE. isw_fy2 .OR. iew_fy .NE. iew_fy2 .OR. jsw_fy .NE. jsw_fy2 .OR. jew_fy .NE. jew_fy2 .OR. &
               isw_cy .NE. isw_cy2 .OR. iew_cy .NE. iew_cy2 .OR. jsw_cy .NE. jsw_cy2 .OR. jew_cy .NE. jew_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine BGRID Y")
          endif
          if( ise_fy .NE. ise_fy2 .OR. iee_fy .NE. iee_fy2 .OR. jse_fy .NE. jse_fy2 .OR. jee_fy .NE. jee_fy2 .OR. &
               ise_cy .NE. ise_cy2 .OR. iee_cy .NE. iee_cy2 .OR. jse_cy .NE. jse_cy2 .OR. jee_cy .NE. jee_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine BGRID Y")
          endif
          if( iss_fy .NE. iss_fy2 .OR. ies_fy .NE. ies_fy2 .OR. jss_fy .NE. jss_fy2 .OR. jes_fy .NE. jes_fy2 .OR. &
               iss_cy .NE. iss_cy2 .OR. ies_cy .NE. ies_cy2 .OR. jss_cy .NE. jss_cy2 .OR. jes_cy .NE. jes_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine BGRID Y")
          endif
          if( isn_fy .NE. isn_fy2 .OR. ien_fy .NE. ien_fy2 .OR. jsn_fy .NE. jsn_fy2 .OR. jen_fy .NE. jen_fy2 .OR. &
               isn_cy .NE. isn_cy2 .OR. ien_cy .NE. ien_cy2 .OR. jsn_cy .NE. jsn_cy2 .OR. jen_cy .NE. jen_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine BGRID Y")
          endif
       endif

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse+shift, nz))
          allocate(y(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = 1.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse+shift
                   y(i,j,k) = 2.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine+shift, jsd_fine:jed_fine+shift, nz))
          allocate(y(isd_fine:ied_fine+shift, jsd_fine:jed_fine+shift, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + k + 1.D-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine+shift
                   y(i,j,k) = i*1.D+6 + j*1.D+3 + k + 2.D-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          if( iew_cx .GE. isw_cx .AND. jew_cx .GE. jsw_cx ) then
             allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
             allocate(wbufferx2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery2(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
          else
             allocate(wbufferx(1,1,1))
             allocate(wbuffery(1,1,1))
             allocate(wbufferx2(1,1,1))
             allocate(wbuffery2(1,1,1))
          endif
          if( iee_cx .GE. ise_cx .AND. jee_cx .GE. jse_cx ) then
             allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy,nz))
             allocate(ebufferx2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery2(ise_cy:iee_cy, jse_cy:jee_cy,nz))
          else
             allocate(ebufferx(1,1,1))
             allocate(ebuffery(1,1,1))
             allocate(ebufferx2(1,1,1))
             allocate(ebuffery2(1,1,1))
          endif
          if( ies_cx .GE. iss_cx .AND. jes_cx .GE. jss_cx ) then
             allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy,nz))
             allocate(sbufferx2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery2(iss_cy:ies_cy, jss_cy:jes_cy,nz))
          else
             allocate(sbufferx(1,1,1))
             allocate(sbuffery(1,1,1))
             allocate(sbufferx2(1,1,1))
             allocate(sbuffery2(1,1,1))
          endif
          if( ien_cx .GE. isn_cx .AND. jen_cx .GE. jsn_cx ) then
             allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
             allocate(nbufferx2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery2(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
          else
             allocate(nbufferx(1,1,1))
             allocate(nbuffery(1,1,1))
             allocate(nbufferx2(1,1,1))
             allocate(nbuffery2(1,1,1))
          endif
          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=BGRID_NE, flags=SCALAR_PAIR)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), &
                  (/iew_c/), (/jsw_c/), (/jew_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  real(1.D+3,r8_kind), real(2.D+3,r8_kind), 1, 1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  real(2.D+3,r8_kind), real(1.D+3,r8_kind), 1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), &
                  (/jss_c/), (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, &
                  ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind), &
                  real(2.D+3, r8_kind), 1, 1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind), &
                  real(1.D+3, r8_kind), 1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/),&
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, &
                  real(1.D+3, r8_kind), real(2.D+3, r8_kind), 1, 1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, &
                  real(2.D+3, r8_kind), real(1.D+3, r8_kind), 1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, &
                  real(1.D+3, r8_kind), real(2.D+3, r8_kind), 1, 1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, &
                  real(2.D+3, r8_kind), real(1.D+3, r8_kind), 1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine BGRID scalar pair Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine BGRID scalar pair Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine BGRID scalar pair Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine BGRID scalar pair Y')
       endif
       if(allocated(x)) deallocate(x)
       if(allocated(y)) deallocate(y)
       if(is_fine_pe) then
          deallocate(wbufferx, ebufferx, sbufferx, nbufferx)
          deallocate(wbufferx2, ebufferx2, sbufferx2, nbufferx2)
          deallocate(wbuffery, ebuffery, sbuffery, nbuffery)
          deallocate(wbuffery2, ebuffery2, sbuffery2, nbuffery2)
       endif

       !---------------------------------------------------------------------------
       !
       !                 Coarse to Fine scalar field, position = CORNER
       !
       !---------------------------------------------------------------------------
       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine+shift, jsd_fine:jed_fine+shift, nz))
          x = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + k
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             allocate(wbuffer(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffer2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
          else
             allocate(wbuffer(1,1,1))
             allocate(wbuffer2(1,1,1))
          endif
          wbuffer = 0; wbuffer2 = 0

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             allocate(ebuffer(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffer2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
          else
             allocate(ebuffer(1,1,1))
             allocate(ebuffer2(1,1,1))
          endif
          ebuffer = 0; ebuffer2 = 0

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             allocate(sbuffer(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffer2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
          else
             allocate(sbuffer(1,1,1))
             allocate(sbuffer2(1,1,1))
          endif
          sbuffer = 0; sbuffer2 = 0

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             allocate(nbuffer(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffer2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
          else
             allocate(nbuffer(1,1,1))
             allocate(nbuffer2(1,1,1))
          endif
          nbuffer = 0; nbuffer2 = 0

       endif

       call mpp_update_nest_fine(x, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=l, position=CORNER)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), &
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbuffer2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r8_kind, 0.0_r8_kind, &
                  1, 1, nx, ny)
          endif
          call compare_checksums(wbuffer, wbuffer2, trim(type2)//' west buffer coarse to fine scalar CORNER')

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), &
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbuffer2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r8_kind, 0.0_r8_kind, &
                  1, 1, nx, ny)
          endif
          call compare_checksums(sbuffer, sbuffer2, trim(type2)//' south buffer coarse to fine scalar CORNER')

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/),&
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebuffer2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift,  &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, &
                  0.0_r8_kind, 0.0_r8_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(ebuffer, ebuffer2, trim(type2)//' east buffer coarse to fine scalar CORNER')

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), &
                  (/jsn_c/), (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, &
                  ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbuffer2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse,&
                  0.0_r8_kind, 0.0_r8_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(nbuffer, nbuffer2, trim(type2)//' north buffer coarse to fine scalar CORNER')

       endif
       if(is_fine_pe) then
          deallocate(wbuffer, ebuffer, sbuffer, nbuffer)
          deallocate(wbuffer2, ebuffer2, sbuffer2, nbuffer2)
       endif
       deallocate(x)

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine CGRID scalar pair update
       !
       !---------------------------------------------------------------------------
       shift = 1
       !--- first check the index is correct or not
       if(is_fine_pe) then
          !--- The index from nest domain
          call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine)
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, l,&
                               & position=EAST)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, l,&
                               & position=EAST)
          call mpp_get_C2F_index(nest_domain,iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH, l,&
                               & position=EAST)
          call mpp_get_C2F_index(nest_domain,isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH, l,&
                               & position=EAST)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, l,&
                               & position=NORTH)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, l,&
                               & position=NORTH)
          call mpp_get_C2F_index(nest_domain,iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH, l,&
                               & position=NORTH)
          call mpp_get_C2F_index(nest_domain,isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH, l,&
                               & position=NORTH)

          !-- The assumed index
          isw_fx2 = 0; iew_fx2 = -1; jsw_fx2 = 0; jew_fx2 = -1
          isw_cx2 = 0; iew_cx2 = -1; jsw_cx2 = 0; jew_cx2 = -1
          ise_fx2 = 0; iee_fx2 = -1; jse_fx2 = 0; jee_fx2 = -1
          ise_cx2 = 0; iee_cx2 = -1; jse_cx2 = 0; jee_cx2 = -1
          iss_fx2 = 0; ies_fx2 = -1; jss_fx2 = 0; jes_fx2 = -1
          iss_cx2 = 0; ies_cx2 = -1; jss_cx2 = 0; jes_cx2 = -1
          isn_fx2 = 0; ien_fx2 = -1; jsn_fx2 = 0; jen_fx2 = -1
          isn_cx2 = 0; ien_cx2 = -1; jsn_cx2 = 0; jen_cx2 = -1
          isw_fy2 = 0; iew_fy2 = -1; jsw_fy2 = 0; jew_fy2 = -1
          isw_cy2 = 0; iew_cy2 = -1; jsw_cy2 = 0; jew_cy2 = -1
          ise_fy2 = 0; iee_fy2 = -1; jse_fy2 = 0; jee_fy2 = -1
          ise_cy2 = 0; iee_cy2 = -1; jse_cy2 = 0; jee_cy2 = -1
          iss_fy2 = 0; ies_fy2 = -1; jss_fy2 = 0; jes_fy2 = -1
          iss_cy2 = 0; ies_cy2 = -1; jss_cy2 = 0; jes_cy2 = -1
          isn_fy2 = 0; ien_fy2 = -1; jsn_fy2 = 0; jen_fy2 = -1
          isn_cy2 = 0; ien_cy2 = -1; jsn_cy2 = 0; jen_cy2 = -1

          !--- west
          if( isc_fine == 1 ) then
             isw_fx2 = isd_fine
             iew_fx2 = isc_fine - 1
             jsw_fx2 = jsd_fine
             jew_fx2 = jed_fine
             isw_cx2 = istart_coarse(my_fine_id)-whalo
             iew_cx2 = istart_coarse(my_fine_id)
             jsw_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
             isw_fy2 = isd_fine
             iew_fy2 = isc_fine - 1
             jsw_fy2 = jsd_fine
             jew_fy2 = jed_fine + shift
             isw_cy2 = istart_coarse(my_fine_id)-whalo
             iew_cy2 = istart_coarse(my_fine_id)
             jsw_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- east
          if( iec_fine == nx_fine ) then
             ise_fx2 = iec_fine+1+shift
             iee_fx2 = ied_fine + shift
             jse_fx2 = jsd_fine
             jee_fx2 = jed_fine
             ise_cx2 = iend_coarse(my_fine_id)+shift
             iee_cx2 = iend_coarse(my_fine_id)+ehalo+shift
             jse_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
             ise_fy2 = iec_fine+1
             iee_fy2 = ied_fine
             jse_fy2 = jsd_fine
             jee_fy2 = jed_fine + shift
             ise_cy2 = iend_coarse(my_fine_id)
             iee_cy2 = iend_coarse(my_fine_id)+ehalo
             jse_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_fx2 = isd_fine
             ies_fx2 = ied_fine + shift
             jss_fx2 = jsd_fine
             jes_fx2 = jsc_fine - 1
             iss_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jss_cx2 = jstart_coarse(my_fine_id)-shalo
             jes_cx2 = jstart_coarse(my_fine_id)
             iss_fy2 = isd_fine
             ies_fy2 = ied_fine
             jss_fy2 = jsd_fine
             jes_fy2 = jsc_fine - 1
             iss_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jss_cy2 = jstart_coarse(my_fine_id)-shalo
             jes_cy2 = jstart_coarse(my_fine_id)
          endif
          !--- north
          if( jec_fine == ny_fine ) then
             isn_fx2 = isd_fine
             ien_fx2 = ied_fine + shift
             jsn_fx2 = jec_fine+1
             jen_fx2 = jed_fine
             isn_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jsn_cx2 = jend_coarse(my_fine_id)
             jen_cx2 = jend_coarse(my_fine_id)+nhalo
             isn_fy2 = isd_fine
             ien_fy2 = ied_fine
             jsn_fy2 = jec_fine+1 + shift
             jen_fy2 = jed_fine + shift
             isn_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jsn_cy2 = jend_coarse(my_fine_id) + shift
             jen_cy2 = jend_coarse(my_fine_id)+nhalo + shift
          endif

          if( isw_fx .NE. isw_fx2 .OR. iew_fx .NE. iew_fx2 .OR. jsw_fx .NE. jsw_fx2 .OR. jew_fx .NE. jew_fx2 .OR. &
               isw_cx .NE. isw_cx2 .OR. iew_cx .NE. iew_cx2 .OR. jsw_cx .NE. jsw_cx2 .OR. jew_cx .NE. jew_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine CGRID X")
          endif
          if( ise_fx .NE. ise_fx2 .OR. iee_fx .NE. iee_fx2 .OR. jse_fx .NE. jse_fx2 .OR. jee_fx .NE. jee_fx2 .OR. &
               ise_cx .NE. ise_cx2 .OR. iee_cx .NE. iee_cx2 .OR. jse_cx .NE. jse_cx2 .OR. jee_cx .NE. jee_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine CGRID X")
          endif
          if( iss_fx .NE. iss_fx2 .OR. ies_fx .NE. ies_fx2 .OR. jss_fx .NE. jss_fx2 .OR. jes_fx .NE. jes_fx2 .OR. &
               iss_cx .NE. iss_cx2 .OR. ies_cx .NE. ies_cx2 .OR. jss_cx .NE. jss_cx2 .OR. jes_cx .NE. jes_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine CGRID X")
          endif
          if( isn_fx .NE. isn_fx2 .OR. ien_fx .NE. ien_fx2 .OR. jsn_fx .NE. jsn_fx2 .OR. jen_fx .NE. jen_fx2 .OR. &
               isn_cx .NE. isn_cx2 .OR. ien_cx .NE. ien_cx2 .OR. jsn_cx .NE. jsn_cx2 .OR. jen_cx .NE. jen_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine CGRID X")
          endif

          if( isw_fy .NE. isw_fy2 .OR. iew_fy .NE. iew_fy2 .OR. jsw_fy .NE. jsw_fy2 .OR. jew_fy .NE. jew_fy2 .OR. &
               isw_cy .NE. isw_cy2 .OR. iew_cy .NE. iew_cy2 .OR. jsw_cy .NE. jsw_cy2 .OR. jew_cy .NE. jew_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine CGRID Y")
          endif
          if( ise_fy .NE. ise_fy2 .OR. iee_fy .NE. iee_fy2 .OR. jse_fy .NE. jse_fy2 .OR. jee_fy .NE. jee_fy2 .OR. &
               ise_cy .NE. ise_cy2 .OR. iee_cy .NE. iee_cy2 .OR. jse_cy .NE. jse_cy2 .OR. jee_cy .NE. jee_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine CGRID Y")
          endif
          if( iss_fy .NE. iss_fy2 .OR. ies_fy .NE. ies_fy2 .OR. jss_fy .NE. jss_fy2 .OR. jes_fy .NE. jes_fy2 .OR. &
               iss_cy .NE. iss_cy2 .OR. ies_cy .NE. ies_cy2 .OR. jss_cy .NE. jss_cy2 .OR. jes_cy .NE. jes_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine CGRID Y")
          endif
          if( isn_fy .NE. isn_fy2 .OR. ien_fy .NE. ien_fy2 .OR. jsn_fy .NE. jsn_fy2 .OR. jen_fy .NE. jen_fy2 .OR. &
               isn_cy .NE. isn_cy2 .OR. ien_cy .NE. ien_cy2 .OR. jsn_cy .NE. jsn_cy2 .OR. jen_cy .NE. jen_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine CGRID Y")
          endif
       endif

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = 1.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y(i,j,k) = 2.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine+shift, jsd_fine:jed_fine, nz))
          allocate(y(isd_fine:ied_fine, jsd_fine:jed_fine+shift, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + k + 1.D-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine
                   y(i,j,k) = i*1.D+6 + j*1.D+3 + k + 2.D-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          if( iew_cx .GE. isw_cx .AND. jew_cx .GE. jsw_cx ) then
             allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
             allocate(wbufferx2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery2(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
          else
             allocate(wbufferx(1,1,1))
             allocate(wbuffery(1,1,1))
             allocate(wbufferx2(1,1,1))
             allocate(wbuffery2(1,1,1))
          endif
          if( iee_cx .GE. ise_cx .AND. jee_cx .GE. jse_cx ) then
             allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy,nz))
             allocate(ebufferx2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery2(ise_cy:iee_cy, jse_cy:jee_cy,nz))
          else
             allocate(ebufferx(1,1,1))
             allocate(ebuffery(1,1,1))
             allocate(ebufferx2(1,1,1))
             allocate(ebuffery2(1,1,1))
          endif
          if( ies_cx .GE. iss_cx .AND. jes_cx .GE. jss_cx ) then
             allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy,nz))
             allocate(sbufferx2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery2(iss_cy:ies_cy, jss_cy:jes_cy,nz))
          else
             allocate(sbufferx(1,1,1))
             allocate(sbuffery(1,1,1))
             allocate(sbufferx2(1,1,1))
             allocate(sbuffery2(1,1,1))
          endif
          if( ien_cx .GE. isn_cx .AND. jen_cx .GE. jsn_cx ) then
             allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
             allocate(nbufferx2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery2(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
          else
             allocate(nbufferx(1,1,1))
             allocate(nbuffery(1,1,1))
             allocate(nbufferx2(1,1,1))
             allocate(nbuffery2(1,1,1))
          endif
          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=CGRID_NE, flags=SCALAR_PAIR)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), &
                  (/jsw_c/), (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, &
                  ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind), &
                  real(2.D+3, r8_kind), 1, 1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind), &
                  real(1.D+3, r8_kind), 1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), &
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind), &
                  real(2.D+3, r8_kind), 1, 1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind), &
                  real(1.D+3, r8_kind), 1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/),&
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind),&
                  real(2.D+3, r8_kind), 1, 1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind), &
                  real(1.D+3, r8_kind), 1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind), &
                  real(2.D+3, r8_kind), 1, 1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, real(2.D+3, r8_kind), &
                  real(1.D+3, r8_kind), 1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine CGRID scalar pair Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine CGRID scalar pair Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine CGRID scalar pair Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine CGRID scalar pair Y')
       endif

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine CGRID vector update
       !
       !---------------------------------------------------------------------------
       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = 1.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y(i,j,k) = 2.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
       else
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + k + 1.D-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine
                   y(i,j,k) = i*1.D+6 + j*1.D+3 + k + 2.D-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=CGRID_NE)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/),&
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind), &
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind), &
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), &
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind), &
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind), &
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), &
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind),&
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind),&
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind),&
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, real(2.D+3, r8_kind),&
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine CGRID vector X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine CGRID vector Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine CGRID vector X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine CGRID vector Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine CGRID vector X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine CGRID vector Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine CGRID vector X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine CGRID vector Y')
       endif

       if(allocated(x)) deallocate(x)
       if(allocated(y)) deallocate(y)
       if(is_fine_pe) then
          deallocate(wbufferx, ebufferx, sbufferx, nbufferx)
          deallocate(wbufferx2, ebufferx2, sbufferx2, nbufferx2)
          deallocate(wbuffery, ebuffery, sbuffery, nbuffery)
          deallocate(wbuffery2, ebuffery2, sbuffery2, nbuffery2)
       endif

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine DGRID vector update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(y(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(x(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   x(i,j,k) = 1.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   y(i,j,k) = 2.D+3 + tile + i*1.D-3 + j*1.D-6 + k*1.D-9
                enddo
             enddo
          enddo
       else
          allocate(y(isd_fine:ied_fine+shift, jsd_fine:jed_fine, nz))
          allocate(x(isd_fine:ied_fine, jsd_fine:jed_fine+shift, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine
                   x(i,j,k) = i*1.D+6 + j*1.D+3 + k + 1.D-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine+shift
                   y(i,j,k) = i*1.D+6 + j*1.D+3 + k + 2.D-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, l,&
                               & position=NORTH)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, l,&
                               & position=NORTH)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH,l,&
                               & position=NORTH)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH,l,&
                               & position=NORTH)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, l,&
                               & position=EAST)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, l,&
                               & position=EAST)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH,l,&
                               & position=EAST)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH,l,&
                               & position=EAST)

          if( iew_cx .GE. isw_cx .AND. jew_cx .GE. jsw_cx ) then
             allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
             allocate(wbufferx2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery2(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
          else
             allocate(wbufferx(1,1,1))
             allocate(wbuffery(1,1,1))
             allocate(wbufferx2(1,1,1))
             allocate(wbuffery2(1,1,1))
          endif
          if( iee_cx .GE. ise_cx .AND. jee_cx .GE. jse_cx ) then
             allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy,nz))
             allocate(ebufferx2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery2(ise_cy:iee_cy, jse_cy:jee_cy,nz))
          else
             allocate(ebufferx(1,1,1))
             allocate(ebuffery(1,1,1))
             allocate(ebufferx2(1,1,1))
             allocate(ebuffery2(1,1,1))
          endif
          if( ies_cx .GE. iss_cx .AND. jes_cx .GE. jss_cx ) then
             allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy,nz))
             allocate(sbufferx2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery2(iss_cy:ies_cy, jss_cy:jes_cy,nz))
          else
             allocate(sbufferx(1,1,1))
             allocate(sbuffery(1,1,1))
             allocate(sbufferx2(1,1,1))
             allocate(sbuffery2(1,1,1))
          endif
          if( ien_cx .GE. isn_cx .AND. jen_cx .GE. jsn_cx ) then
             allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
             allocate(nbufferx2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery2(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
          else
             allocate(nbufferx(1,1,1))
             allocate(nbuffery(1,1,1))
             allocate(nbufferx2(1,1,1))
             allocate(nbuffery2(1,1,1))
          endif

          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=DGRID_NE)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/),&
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind),&
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind),&
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/),&
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind),&
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind),&
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/),&
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(1.D+3, r8_kind),&
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind),&
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/),&
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, real(1.D+3, r8_kind),&
                  real(2.D+3, r8_kind), 1, -1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse,&
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, real(2.D+3, r8_kind),&
                  real(1.D+3, r8_kind), -1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine DGRID vector X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine DGRID vector Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine DGRID vector X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine DGRID vector Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine DGRID vector X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine DGRID vector Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine DGRID vector X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine DGRID vector Y')
       endif

       if(allocated(x)) deallocate(x)
       if(allocated(y)) deallocate(y)
       if(is_fine_pe) then
          deallocate(wbufferx, ebufferx, sbufferx, nbufferx)
          deallocate(wbufferx2, ebufferx2, sbufferx2, nbufferx2)
          deallocate(wbuffery, ebuffery, sbuffery, nbuffery)
          deallocate(wbuffery2, ebuffery2, sbuffery2, nbuffery2)
       endif
       endif
       deallocate(my_pelist, my_pelist_fine)
       call mpp_set_current_pelist()

    enddo

    call mpp_set_current_pelist(pelist)
    call mpp_sync(pelist)
    deallocate(pelist)

  end subroutine test_update_nest_domain_r8

  !############################################################################
  !--- this routine will get number of nest.
  subroutine convert_index_up(domain, rotate, ncross, is_coarse, ie_coarse, js_coarse, je_coarse, &
                                is_in, ie_in, js_in, je_in, is_out, ie_out, js_out, je_out)
    type(domain2D), intent(in) :: domain
    integer, intent(in)  :: is_coarse, ie_coarse, js_coarse, je_coarse
    integer, intent(in)  :: is_in, ie_in, js_in, je_in, rotate, ncross
    integer, intent(out) :: is_out, ie_out, js_out, je_out
    integer :: isg, ieg, jsg, jeg

    call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)

    if( je_coarse > jeg .and. ie_coarse > ieg ) then
       call mpp_error(FATAL,"convert_index_up:  je_coarse > jeg .and. ie_convert > ieg")
    else if (je_coarse > jeg) then
       select case(rotate)
       case(0)
          is_out = is_in
          ie_out = ie_in
          js_out = js_in + ncross*jeg
          je_out = je_in + ncross*jeg
       case(90)
          is_out = js_in + ncross*jeg
          ie_out = je_in + ncross*jeg
          js_out = jeg - ie_in
          je_out = jeg - is_in
       case default
          call mpp_error(FATAL, "convert_index_back: rotate should be 0 or 90 when je_in>jeg")
       end select
    else if (ie_coarse > ieg) then
       select case(rotate)
       case(0)
          is_out = is_in + ncross*ieg
          ie_out = ie_in + ncross*ieg
          js_out = js_in
          je_out = je_in
       case(-90)
          js_out = is_in + ncross*ieg
          je_out = ie_in + ncross*ieg
          is_out = ieg - je_in
          ie_out = ieg - js_in
       case default
          call mpp_error(FATAL, "convert_index_back: rotate should be 0 or -90 when ie_in>ieg")
       end select
    else
       is_out = is_in
       ie_out = ie_in
       js_out = js_in
       je_out = je_in
    endif

  end subroutine convert_index_up

!###############################################################################
  subroutine test_update_nest_domain_r4( type )
    character(len=*), intent(in) :: type
    logical                      :: cubic_grid
    logical                      :: is_fine_pe, is_coarse_pe
    integer                      :: n, i, j, k
    integer                      :: ntiles, npes_per_tile
    integer                      :: npes_fine, pos
    integer                      :: isc_coarse, iec_coarse, jsc_coarse, jec_coarse
    integer                      :: isd_coarse, ied_coarse, jsd_coarse, jed_coarse
    integer                      :: isd_fine, ied_fine, jsd_fine, jed_fine
    integer                      :: isc_fine, iec_fine, jsc_fine, jec_fine
    integer                      :: nx_fine, ny_fine, nx_coarse, ny_coarse
    integer                      :: nxc_fine, nyc_fine, nxc_coarse, nyc_coarse
    integer                      :: isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c
    integer                      :: ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c
    integer                      :: iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c
    integer                      :: isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c
    integer                      :: isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx
    integer                      :: ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx
    integer                      :: iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx
    integer                      :: isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx
    integer                      :: isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy
    integer                      :: ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy
    integer                      :: iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy
    integer                      :: isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy
    integer                      :: isw_f2, iew_f2, jsw_f2, jew_f2, isw_c2, iew_c2, jsw_c2, jew_c2, tile_w2
    integer                      :: ise_f2, iee_f2, jse_f2, jee_f2, ise_c2, iee_c2, jse_c2, jee_c2, tile_e2
    integer                      :: iss_f2, ies_f2, jss_f2, jes_f2, iss_c2, ies_c2, jss_c2, jes_c2, tile_s2
    integer                      :: isn_f2, ien_f2, jsn_f2, jen_f2, isn_c2, ien_c2, jsn_c2, jen_c2, tile_n2
    integer                      :: isw_fx2, iew_fx2, jsw_fx2, jew_fx2, isw_cx2, iew_cx2, jsw_cx2, jew_cx2, tile_wx2
    integer                      :: ise_fx2, iee_fx2, jse_fx2, jee_fx2, ise_cx2, iee_cx2, jse_cx2, jee_cx2, tile_ex2
    integer                      :: iss_fx2, ies_fx2, jss_fx2, jes_fx2, iss_cx2, ies_cx2, jss_cx2, jes_cx2, tile_sx2
    integer                      :: isn_fx2, ien_fx2, jsn_fx2, jen_fx2, isn_cx2, ien_cx2, jsn_cx2, jen_cx2, tile_nx2
    integer                      :: isw_fy2, iew_fy2, jsw_fy2, jew_fy2, isw_cy2, iew_cy2, jsw_cy2, jew_cy2, tile_wy2
    integer                      :: ise_fy2, iee_fy2, jse_fy2, jee_fy2, ise_cy2, iee_cy2, jse_cy2, jee_cy2, tile_ey2
    integer                      :: iss_fy2, ies_fy2, jss_fy2, jes_fy2, iss_cy2, ies_cy2, jss_cy2, jes_cy2, tile_sy2
    integer                      :: isn_fy2, ien_fy2, jsn_fy2, jen_fy2, isn_cy2, ien_cy2, jsn_cy2, jen_cy2, tile_ny2
    integer                      :: isw_f_T, iew_f_T, jsw_f_T, jew_f_T, isw_c_T, iew_c_T, jsw_c_T, jew_c_T
    integer                      :: ise_f_T, iee_f_T, jse_f_T, jee_f_T, ise_c_T, iee_c_T, jse_c_T, jee_c_T
    integer                      :: iss_f_T, ies_f_T, jss_f_T, jes_f_T, iss_c_T, ies_c_T, jss_c_T, jes_c_T
    integer                      :: isn_f_T, ien_f_T, jsn_f_T, jen_f_T, isn_c_T, ien_c_T, jsn_c_T, jen_c_T
    integer                      :: is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f
    integer                      :: is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx
    integer                      :: is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy
    integer                      :: tile, position, shift
    integer                      :: layout_fine(2), my_fine_id
    integer, allocatable         :: pelist(:), start_pos(:), end_pos(:)
    integer, allocatable         :: my_pelist_fine(:)
    integer, allocatable         :: pe_start(:), pe_end(:)
    integer, allocatable         :: layout2D(:,:), global_indices(:,:)
    real(kind=r4_kind), allocatable :: x(:,:,:), x1(:,:,:), x2(:,:,:)
    real(kind=r4_kind), allocatable :: y(:,:,:), y1(:,:,:), y2(:,:,:)
    real(kind=r4_kind), allocatable :: wbuffer(:,:,:), wbuffer2(:,:,:)
    real(kind=r4_kind), allocatable :: ebuffer(:,:,:), ebuffer2(:,:,:)
    real(kind=r4_kind), allocatable :: sbuffer(:,:,:), sbuffer2(:,:,:)
    real(kind=r4_kind), allocatable :: nbuffer(:,:,:), nbuffer2(:,:,:)
    real(kind=r4_kind), allocatable :: wbufferx(:,:,:), wbufferx2(:,:,:)
    real(kind=r4_kind), allocatable :: ebufferx(:,:,:), ebufferx2(:,:,:)
    real(kind=r4_kind), allocatable :: sbufferx(:,:,:), sbufferx2(:,:,:)
    real(kind=r4_kind), allocatable :: nbufferx(:,:,:), nbufferx2(:,:,:)
    real(kind=r4_kind), allocatable :: wbuffery(:,:,:), wbuffery2(:,:,:)
    real(kind=r4_kind), allocatable :: ebuffery(:,:,:), ebuffery2(:,:,:)
    real(kind=r4_kind), allocatable :: sbuffery(:,:,:), sbuffery2(:,:,:)
    real(kind=r4_kind), allocatable :: nbuffery(:,:,:), nbuffery2(:,:,:)
    integer                      :: x_refine(num_nest), y_refine(num_nest)
    integer                      :: istart_fine(num_nest), iend_fine(num_nest)
    integer                      :: jstart_fine(num_nest), jend_fine(num_nest)
    integer                      :: iend_coarse(num_nest), jend_coarse(num_nest)
    integer                      :: is_fine(6*num_nest), ie_fine(6*num_nest)
    integer                      :: js_fine(6*num_nest), je_fine(6*num_nest)
    integer                      :: is_coarse(6*num_nest), ie_coarse(6*num_nest)
    integer                      :: js_coarse(6*num_nest), je_coarse(6*num_nest)
    integer                      :: t_coarse(6*num_nest), rotate_coarse(6*num_nest)
    integer                      :: iadd_coarse(6*num_nest), jadd_coarse(6*num_nest)
    integer                      :: nnest
    character(len=128)           :: type2
    character(len=32)            :: text, pelist_name
    type(domain2d)               :: domain
    type(domain2d), pointer      :: domain_coarse=>NULL()
    type(domain2d), pointer      :: domain_fine=>NULL()
    type(nest_domain_type)       :: nest_domain
    logical                      :: x_cyclic, y_cyclic
    integer                      :: my_tile_id(1), my_num_nest
    integer, dimension(num_nest) :: my_tile_coarse, my_tile_fine, my_istart_coarse, my_iend_coarse
    integer, dimension(num_nest) :: my_jstart_coarse, my_jend_coarse
    integer                      :: ntiles_nest_top, npes_nest_top, num_nest_level, my_npes, l
    integer                      :: npes_my_fine, npes_my_level
    integer, allocatable         :: my_pelist(:)
    integer, dimension(1)        :: dummy

    x_cyclic = .false.
    y_cyclic = .false.
    if(cyclic_nest(1) == 'X') then
       x_cyclic = .true.
    else if(cyclic_nest(1) == 'Y') then
       y_cyclic = .true.
    endif

    istart_fine = 0; iend_fine = -1
    jstart_fine = 0; jend_fine = -1
    iend_coarse = -1; jend_coarse = -1
    is_fine = 0;  ie_fine = -1
    js_fine = 0;  je_fine = -1
    is_coarse = 0;  ie_coarse = -1
    js_coarse = 0;  je_coarse = -1
    t_coarse = 0; rotate_coarse = -1;
    iadd_coarse = 0; jadd_coarse = 0

    select case(type)
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic is zero, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       nx = nx_cubic
       ny = ny_cubic
       ntiles_nest_top = 6
       cubic_grid = .true.

    case ( 'Cubed-sphere, single face' )
      if( nx_cubic == 0 ) then
        call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic is zero, '//&
                'No test is done for Cubic-Grid mosaic. ' )
        return
      endif
      if( nx_cubic .NE. ny_cubic ) then
        call mpp_error(NOTE,'test_update_nest_domain: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                'No test is done for Cubic-Grid mosaic. ' )
        return
      endif
      nx = nx_cubic
      ny = ny_cubic
      ntiles_nest_top = 1
      cubic_grid = .false.

    case ( 'Regional' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Regional grid mosaic, nx_cubic is zero, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_update_nest_domain: for Regional grid mosaic, nx_cubic does not equal ny_cubic,'//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       nx = nx_cubic
       ny = ny_cubic
       ntiles_nest_top = 1
       cubic_grid = .false.
    case default
       call mpp_error(FATAL, 'test_update_nest_domain: no such test: '//type)
    end select

    if(ntiles_nest_all > MAX_NTILE) call mpp_error(FATAL, 'test_update_nest_domain: ntiles_nest_all > MAX_NTILE')
    if(ntiles_nest_top .GE. ntiles_nest_all) call mpp_error(FATAL, &
                            'test_update_nest_domain: ntiles_nest_top .GE. ntile_nest_all')
    if(ntiles_nest_all .NE. ntiles_nest_top + num_nest) call mpp_error(FATAL, &
             'test_update_nest_domain: ntiles_nest_all .NE. ntiles_nest_top + num_nest')
    !--- for the ntiles_nest_top, number of processors should be same
    do n = 1, ntiles_nest_all
       if(npes_nest_tile(n) .LE. 0) call mpp_error(FATAL, &
            'test_update_nest_domain: npes_nest_tile is not properly set')
    enddo
    do n = 2, ntiles_nest_top
       if(npes_nest_tile(n) .NE. npes_nest_tile(n-1)) call mpp_error(FATAL, &
            'test_update_nest_domain: each tile of top mosaic grid should use same number of MPI ranks')
    enddo
    npes_nest_top = ntiles_nest_top * npes_nest_tile(1)

    npes = mpp_npes()

    !--- make sure sum(npes_nest_tile) == npes
    if(sum(npes_nest_tile(1:ntiles_nest_all)) .NE. npes ) &
         call mpp_error(FATAL, "test_mpp_domains: sum(npes_nest_tile) .NE. npes")

    !--- make sure tile_fine are monotonically increasing and equal to ntiles_nest_top + nest number
    do n = 1, num_nest
       if(tile_fine(n) .NE. ntiles_nest_top+n) call mpp_error(FATAL, &
           "test_mpp_domains: tile_fine(n) .NE. ntiles_nest_top+n")
    enddo

    !---make sure nest_level is setup properly
    if(nest_level(1) .NE. 1) call mpp_error(FATAL, "test_mpp_domains: nest_level(1) .NE. 1")
    do n = 2, num_nest
       if(nest_level(n) > nest_level(n-1)+1)call mpp_error(FATAL,"test_mpp_domains: nest_level(n) > nest_level(n-1)+1")
       if(nest_level(n) < nest_level(n-1) ) call mpp_error(FATAL, "test_mpp_domains: nest_level(n) < nest_level(n-1)")
    enddo
    num_nest_level = nest_level(num_nest)

    allocate(pelist(npes))
    call mpp_get_current_pelist(pelist)

    !--- compute iend_coarse and jend_coarse
    do n = 1, num_nest
       iend_coarse(n) = istart_coarse(n) + icount_coarse(n) - 1
       jend_coarse(n) = jstart_coarse(n) + jcount_coarse(n) - 1
       istart_fine(n) = 1; iend_fine(n) = icount_coarse(n)*refine_ratio(n)
       jstart_fine(n) = 1; jend_fine(n) = jcount_coarse(n)*refine_ratio(n)
    enddo

    !--- first define the top level grid mosaic domain.

    !--- setup pelist for top level
    allocate(my_pelist(npes_nest_top))
    do n = 1, npes_nest_top
       my_pelist(n) = pelist(n)
    enddo
    call mpp_declare_pelist(my_pelist)
    if(ANY(my_pelist==mpp_pe())) then
       call mpp_set_current_pelist(my_pelist)

       allocate(layout2D(2,ntiles_nest_top), global_indices(4,ntiles_nest_top), pe_start(ntiles_nest_top), &
                pe_end(ntiles_nest_top) )
       npes_per_tile = npes_nest_tile(1)

       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles_nest_top
          global_indices(:,n) = (/1,nx,1,ny/)
          if (ANY(layout_cubic == 0)) then
             layout2D(:,n)         = layout
          else
             layout2D(:,n)         = layout_cubic(:)
           endif
       end do
       do n = 1, ntiles_nest_top
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do

       if( cubic_grid ) then
         call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                                   global_indices, layout2D, pe_start, pe_end )
       else
          call mpp_define_domains(global_indices(:,ntiles_nest_top), layout2D(:,ntiles_nest_top), domain, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.true., name=trim(type)//' top level grid', tile_id=1  )
       endif
       call mpp_get_compute_domain(domain, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
       call mpp_get_data_domain(domain, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
       deallocate(layout2D, global_indices, pe_start, pe_end )
    endif

    call mpp_set_current_pelist()
    deallocate(my_pelist)
    !--- define domain for all the nest regoin.
    pos = npes_nest_top
    do n = 1, num_nest
       my_npes = npes_nest_tile(tile_fine(n))
       allocate(my_pelist(my_npes))
       my_pelist(:) = pelist(pos+1:pos+my_npes)
       call mpp_declare_pelist(my_pelist)
       if(ANY(my_pelist==mpp_pe())) then
          call mpp_set_current_pelist(my_pelist)
          nx_fine = iend_fine(n) - istart_fine(n) + 1
          ny_fine = jend_fine(n) - jstart_fine(n) + 1
          if (ANY(layout_nest(:,n) == 0)) then
             call mpp_define_layout( (/1,nx_fine,1,ny_fine/), my_npes, layout )
          else
             layout(:)         = layout_nest(:,n)
          endif
          call mpp_define_domains((/1,nx_fine,1,ny_fine/), layout, domain, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.true., name=trim(type)//' fine grid', tile_id = tile_fine(n) )
          call mpp_get_compute_domain(domain, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain, isd_fine, ied_fine, jsd_fine, jed_fine)
          !--- test halo update for nested region.
          call test_nest_halo_update_r4(domain)
       endif
       pos = pos+my_npes
       deallocate(my_pelist)
       call mpp_set_current_pelist()
    enddo

    !--- reset to the global pelist
    call mpp_set_current_pelist()

    x_refine(:) = refine_ratio(1:num_nest)
    y_refine(:) = refine_ratio(1:num_nest)

    call mpp_define_nest_domains(nest_domain, domain, num_nest, nest_level(1:num_nest), tile_fine(1:num_nest), &
             tile_coarse(1:num_nest), istart_coarse(1:num_nest), icount_coarse(1:num_nest), jstart_coarse(1:num_nest),&
             jcount_coarse(1:num_nest), npes_nest_tile(1:ntiles_nest_all), &
             x_refine(1:num_nest), y_refine(1:num_nest), extra_halo=extra_halo, name="nest_domain")

    !--- loop over nest level
    do l = 1, num_nest_level
       npes_my_level = mpp_get_nest_npes(nest_domain, l)
       npes_my_fine = mpp_get_nest_fine_npes(nest_domain,l)
       allocate(my_pelist(npes_my_level))
       allocate(my_pelist_fine(npes_my_fine))
       call mpp_get_nest_pelist(nest_domain, l, my_pelist)

       call mpp_declare_pelist(my_pelist(:))
       write(type2, '(a,I2)')trim(type)//" nest_level = ",l
       if(ANY(my_pelist(:) == mpp_pe())) then

          call mpp_get_nest_fine_pelist(nest_domain, l, my_pelist_fine)

          call mpp_set_current_pelist(my_pelist)
          my_tile_id = mpp_get_tile_id(domain)
          domain_coarse => mpp_get_nest_coarse_domain(nest_domain, nest_level=l)
          domain_fine => mpp_get_nest_fine_domain(nest_domain, nest_level=l)
          is_fine_pe = mpp_is_nest_fine(nest_domain, l)
          is_coarse_pe = mpp_is_nest_coarse(nest_domain, l)
          if(is_fine_pe .eqv. is_coarse_pe) call mpp_error(FATAL, "test_mpp_domains: is_fine_pe .eqv. is_coarse_pe")
          my_num_nest = 0
          my_fine_id = 0
          do n = 1, num_nest
             if(nest_level(n)==l) then
                my_num_nest = my_num_nest+1
                my_tile_coarse(my_num_nest) = tile_coarse(n)
                my_tile_fine(my_num_nest) = tile_fine(n)
                my_istart_coarse(my_num_nest) = istart_coarse(n)
                my_iend_coarse(my_num_nest) = iend_coarse(n)
                my_jstart_coarse(my_num_nest) = jstart_coarse(n)
                my_jend_coarse(my_num_nest) = jend_coarse(n)
                if(my_tile_id(1) == tile_fine(n)) my_fine_id = n
             endif
          enddo
          !--- each nest region might be over multiple face of cubic sphere grid.
          !---Get the number of nest region with consideration of face.
          call get_nnest(domain_coarse, my_num_nest, my_tile_coarse, my_istart_coarse, my_iend_coarse, &
               my_jstart_coarse, my_jend_coarse, nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, &
               is_coarse, ie_coarse, js_coarse, je_coarse)

          !---------------------------------------------------------------------------
          !
          !                    fine to coarse scalar field, limit to position=CENTER.
          !
          !---------------------------------------------------------------------------
          if(is_fine_pe) then
             call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine)
             call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine)
          endif

          if(is_coarse_pe) then
             call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
             call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          endif

          if(is_fine_pe) then
             call mpp_get_F2C_index(nest_domain, is_c, ie_c, js_c, je_c, is_f, ie_f, js_f, je_f, l, position=CENTER)
             allocate(x(is_c:ie_c, js_c:je_c, nz))
             x = 0
             do k = 1, nz
                do j = js_c, je_c
                   do i = is_c, ie_c
                      x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.001
                   enddo
                enddo
             enddo
          else
             allocate(x1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse, nz))
             allocate(x2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse, nz))
             x1 = 0
             tile = my_tile_id(1)

             do k = 1, nz
                do j = jsc_coarse, jec_coarse
                   do i = isc_coarse, iec_coarse
                      x1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.002
                   enddo
                enddo
             enddo
             x2 = x1
          endif


          if(is_coarse_pe) then
             do n = 1, nnest
                is_c = max(is_coarse(n), isc_coarse)
                ie_c = min(ie_coarse(n),   iec_coarse)
                js_c = max(js_coarse(n), jsc_coarse)
                je_c = min(je_coarse(n),   jec_coarse)
                if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c .GE. js_c ) then
                   call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), is_c, ie_c, js_c, je_c,&
                         nz, isd_coarse, jsd_coarse, nx, ny, 0, 0, 0.001_r4_kind, 0.001_r4_kind, 1, 1, &
                        .false., .false., iend_coarse(1), jend_coarse(1) )
                endif
             enddo
          endif

          call mpp_update_nest_coarse(x, nest_domain, x1, nest_level=l, position=CENTER)

          !--- compare with assumed value
          if( is_coarse_pe) then
             call compare_checksums(x1, x2, trim(type2)//' fine to coarse scalar')
          endif
          if(allocated(x))       deallocate(x)
          if(allocated(x1))      deallocate(x1)
          if(allocated(x2))      deallocate(x2)

       !---------------------------------------------------------------------------
       !
       !                    fine to coarse CGRID scalar pair update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_fine_pe) then
          call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, l, position=EAST)
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l,position=NORTH)
          allocate(x(is_cx:ie_cx, js_cx:je_cx, nz))
          allocate(y(is_cy:ie_cy, js_cy:je_cy, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = js_cx, je_cx
                do i = is_cx, ie_cx
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 1.0E-6
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = js_cy, je_cy
                do i = is_cy, ie_cy
                   y(i,j,k) = i*1.e+6 + j*1.e+3 + k + 2.0E-6
                enddo
             enddo
          enddo
          if(x_cyclic) then
             if(ie_cx == iend_coarse(1)+1) then
                i = ie_cx
                do k = 1, nz
                   do j = js_cx, je_cx
                      x(i,j,k) = istart_coarse(1)*1.e+6 + j*1.e+3 + k + 1.0E-6
                   enddo
                enddo
             endif
          endif
          if(y_cyclic) then
             if(je_cx == jend_coarse(1)+1) then
                j = je_cx
                do k = 1, nz
                   do i = is_cx, ie_cx
                      y(i,j,k) = i*1.e+6 + jstart_coarse(1)*1.e+3 + k + 1.0E-6
                   enddo
                enddo
             endif
          endif
       else
          allocate(x1(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(x2(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(y2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x1 = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.001
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.002
                enddo
             enddo
          enddo
          x2 = x1
          y2 = y1
       endif


       if(is_coarse_pe) then
          do n = 1, nnest
             is_c = max(is_coarse(n), isc_coarse)
             ie_c = min(ie_coarse(n),   iec_coarse)
             js_c = max(js_coarse(n), jsc_coarse)
             je_c = min(je_coarse(n),   jec_coarse)
             if( tile == t_coarse(n) .AND. ie_c+shift .GE. is_c .AND. je_c .GE. js_c ) then
                call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0, 1.0E-6_r4_kind, &
                     2.0E-6_r4_kind, 1, 1, x_cyclic, .false., iend_coarse(1)+1, jend_coarse(1)+1)
             endif
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift, 2.0E-6_r4_kind, &
                     1.0E-6_r4_kind, 1, 1, .false., y_cyclic, iend_coarse(1)+1, jend_coarse(1)+1)
             endif
          enddo
       endif

       call mpp_update_nest_coarse(x, y, nest_domain, x1, y1, nest_level=l, gridtype=CGRID_NE, flags=SCALAR_PAIR)

       !--- compare with assumed value
       if( is_coarse_pe) then
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Scalar_pair X')
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Scalar_pair Y')
       endif
       if(allocated(x))       deallocate(x)
       if(allocated(x1))      deallocate(x1)
       if(allocated(x2))      deallocate(x2)
       if(allocated(y))       deallocate(y)
       if(allocated(y1))      deallocate(y1)
       if(allocated(y2))      deallocate(y2)

       !---------------------------------------------------------------------------
       !
       !                    fine to coarse CGRID vector update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_fine_pe) then
          call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, l, position=EAST)
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l,position=NORTH)
          allocate(x(is_cx:ie_cx, js_cx:je_cx, nz))
          allocate(y(is_cy:ie_cy, js_cy:je_cy, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = js_cx, je_cx
                do i = is_cx, ie_cx
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 1.0E-6
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = js_cy, je_cy
                do i = is_cy, ie_cy
                   y(i,j,k) = i*1.e+6 + j*1.e+3 + k + 2.0E-6
                enddo
             enddo
          enddo
       else
          allocate(x1(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(x2(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(y2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x1 = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.001
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.002
                enddo
             enddo
          enddo
          x2 = x1
          y2 = y1
       endif


       if(is_coarse_pe) then
          do n = 1, nnest
             is_c = max(is_coarse(n), isc_coarse)
             ie_c = min(ie_coarse(n),   iec_coarse)
             js_c = max(js_coarse(n), jsc_coarse)
             je_c = min(je_coarse(n),   jec_coarse)
             if( tile == t_coarse(n) .AND. ie_c+shift .GE. is_c .AND. je_c .GE. js_c ) then
                call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0, 1.0E-6_r4_kind, &
                     2.0E-6_r4_kind, 1, -1, x_cyclic, .false., iend_coarse(1)+1, jend_coarse(1)+1)
             endif
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift, 2.0E-6_r4_kind, &
                     1.0E-6_r4_kind, -1, 1, .false., y_cyclic, iend_coarse(1)+1, jend_coarse(1)+1)
             endif
          enddo
       endif

       call mpp_update_nest_coarse(x, y, nest_domain, x1, y1, nest_level=l, gridtype=CGRID_NE)

       !--- compare with assumed value
       if( is_coarse_pe) then
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Vector X')
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer CGRID Vector Y')
       endif
       if(allocated(x))       deallocate(x)
       if(allocated(x1))      deallocate(x1)
       if(allocated(x2))      deallocate(x2)
       if(allocated(y))       deallocate(y)
       if(allocated(y1))      deallocate(y1)
       if(allocated(y2))      deallocate(y2)

       !---------------------------------------------------------------------------
       !
       !                    fine to coarse DGRID vector update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_fine_pe) then
          call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, l,position=NORTH)
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l, position=EAST)
          allocate(x(is_cx:ie_cx, js_cx:je_cx, nz))
          allocate(y(is_cy:ie_cy, js_cy:je_cy, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = js_cx, je_cx
                do i = is_cx, ie_cx
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 1.0E-6
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = js_cy, je_cy
                do i = is_cy, ie_cy
                   y(i,j,k) = i*1.e+6 + j*1.e+3 + k + 2.0E-6
                enddo
             enddo
          enddo
       else
          allocate(x1(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(x2(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          allocate(y1(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y2(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          x1 = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   x1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.001
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   y1(i,j,k) = i*1.e+6 + j*1.e+3 + k + 0.002
                enddo
             enddo
          enddo
          x2 = x1
          y2 = y1
       endif


       if(is_coarse_pe) then
          do n = 1, nnest
             is_c = max(is_coarse(n), isc_coarse)
             ie_c = min(ie_coarse(n),   iec_coarse)
             js_c = max(js_coarse(n), jsc_coarse)
             je_c = min(je_coarse(n),   jec_coarse)
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift, 1.0E-6_r4_kind, &
                     2.0E-6_r4_kind, 1, -1, .false., y_cyclic, iend_coarse(1), jend_coarse(1) )
             endif
             if( tile == t_coarse(n) .AND. ie_c+shift .GE. is_c .AND. je_c .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0, 2.0E-6_r4_kind, &
                     1.0E-6_r4_kind, -1, 1, x_cyclic, .false., iend_coarse(1), jend_coarse(1))
             endif
          enddo
       endif

       call mpp_update_nest_coarse(x, y, nest_domain, x1, y1, nest_level=l, gridtype=DGRID_NE)

       !--- compare with assumed value
       if( is_coarse_pe) then
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer DGRID Vector X')
          call compare_checksums(x1, x2, trim(type2)//' fine to coarse buffer DGRID Vector Y')
       endif
       if(allocated(x))       deallocate(x)
       if(allocated(x1))      deallocate(x1)
       if(allocated(x2))      deallocate(x2)
       if(allocated(y))       deallocate(y)
       if(allocated(y1))      deallocate(y1)
       if(allocated(y2))      deallocate(y2)

       !---------------------------------------------------------------------------
       !
       !                 Coarse to Fine scalar field, position = CENTER
       !
       !---------------------------------------------------------------------------
       !--- first check the index is correct or not
       !--- The index from nest domain
       call mpp_get_C2F_index(nest_domain, isw_f, iew_f, jsw_f, jew_f, isw_c, iew_c, jsw_c, jew_c, WEST, l)
       call mpp_get_C2F_index(nest_domain, ise_f, iee_f, jse_f, jee_f, ise_c, iee_c, jse_c, jee_c, EAST, l)
       call mpp_get_C2F_index(nest_domain, iss_f, ies_f, jss_f, jes_f, iss_c, ies_c, jss_c, jes_c, SOUTH, l)
       call mpp_get_C2F_index(nest_domain, isn_f, ien_f, jsn_f, jen_f, isn_c, ien_c, jsn_c, jen_c, NORTH, l)

       if(is_fine_pe) then
          call mpp_get_compute_domain(domain, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain, isd_fine, ied_fine, jsd_fine, jed_fine)

          !-- The assumed index
          isw_f2 = 0; iew_f2 = -1; jsw_f2 = 0; jew_f2 = -1
          isw_c2 = 0; iew_c2 = -1; jsw_c2 = 0; jew_c2 = -1
          ise_f2 = 0; iee_f2 = -1; jse_f2 = 0; jee_f2 = -1
          ise_c2 = 0; iee_c2 = -1; jse_c2 = 0; jee_c2 = -1
          iss_f2 = 0; ies_f2 = -1; jss_f2 = 0; jes_f2 = -1
          iss_c2 = 0; ies_c2 = -1; jss_c2 = 0; jes_c2 = -1
          isn_f2 = 0; ien_f2 = -1; jsn_f2 = 0; jen_f2 = -1
          isn_c2 = 0; ien_c2 = -1; jsn_c2 = 0; jen_c2 = -1

          !--- west
          if( isc_fine == 1 ) then
             isw_f2 = isd_fine; iew_f2 = isc_fine - 1
             jsw_f2 = jsd_fine; jew_f2 = jed_fine
             isw_c2 = istart_coarse(my_fine_id)-whalo
             iew_c2 = istart_coarse(my_fine_id)
             jsw_c2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_c2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
          endif
          !--- east
          if( iec_fine == nx_fine ) then
             ise_f2 = iec_fine+1; iee_f2 = ied_fine
             jse_f2 = jsd_fine;   jee_f2 = jed_fine
             ise_c2 = iend_coarse(my_fine_id)
             iee_c2 = iend_coarse(my_fine_id)+ehalo
             jse_c2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_c2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_f2 = isd_fine; ies_f2 = ied_fine
             jss_f2 = jsd_fine; jes_f2 = jsc_fine - 1
             iss_c2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_c2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jss_c2 = jstart_coarse(my_fine_id)-shalo
             jes_c2 = jstart_coarse(my_fine_id)
          endif
          !--- north
          if( jec_fine == ny_fine ) then
             isn_f2 = isd_fine;  ien_f2 = ied_fine
             jsn_f2 = jec_fine+1; jen_f2 = jed_fine
             isn_c2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_c2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jsn_c2 = jend_coarse(my_fine_id)
             jen_c2 = jend_coarse(my_fine_id)+nhalo
          endif

          if( isw_f .NE. isw_f2 .OR. iew_f .NE. iew_f2 .OR. jsw_f .NE. jsw_f2 .OR. jew_f .NE. jew_f2 .OR. &
               isw_c .NE. isw_c2 .OR. iew_c .NE. iew_c2 .OR. jsw_c .NE. jsw_c2 .OR. jew_c .NE. jew_c2 ) then
             write(5000+mpp_pe(),*) "west buffer fine index = ", isw_f, iew_f, jsw_f, jew_f
             write(5000+mpp_pe(),*) "west buffer fine index2 = ", isw_f2, iew_f2, jsw_f2, jew_f2
             write(5000+mpp_pe(),*) "west buffer coarse index = ", isw_c, iew_c, jsw_c, jew_c
             write(5000+mpp_pe(),*) "west buffer coarse index2 = ", isw_c2, iew_c2, jsw_c2, jew_c2
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine scalar")
          endif
          if( ise_f .NE. ise_f2 .OR. iee_f .NE. iee_f2 .OR. jse_f .NE. jse_f2 .OR. jee_f .NE. jee_f2 .OR. &
               ise_c .NE. ise_c2 .OR. iee_c .NE. iee_c2 .OR. jse_c .NE. jse_c2 .OR. jee_c .NE. jee_c2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine scalar")
          endif
          if( iss_f .NE. iss_f2 .OR. ies_f .NE. ies_f2 .OR. jss_f .NE. jss_f2 .OR. jes_f .NE. jes_f2 .OR. &
               iss_c .NE. iss_c2 .OR. ies_c .NE. ies_c2 .OR. jss_c .NE. jss_c2 .OR. jes_c .NE. jes_c2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine scalar")
          endif
          if( isn_f .NE. isn_f2 .OR. ien_f .NE. ien_f2 .OR. jsn_f .NE. jsn_f2 .OR. jen_f .NE. jen_f2 .OR. &
               isn_c .NE. isn_c2 .OR. ien_c .NE. ien_c2 .OR. jsn_c .NE. jsn_c2 .OR. jen_c .NE. jen_c2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine scalar")
          endif
       endif

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse, jsd_coarse:jed_coarse, nz))
          x = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse
                   x(i,j,k) = tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine, jsd_fine:jed_fine, nz))
          x = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k
                enddo
             enddo
          enddo
       endif

       if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
          allocate(wbuffer(isw_c:iew_c, jsw_c:jew_c,nz))
          allocate(wbuffer2(isw_c:iew_c, jsw_c:jew_c,nz))
       else
          allocate(wbuffer(1,1,1))
          allocate(wbuffer2(1,1,1))
       endif

       if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
          allocate(ebuffer(ise_c:iee_c, jse_c:jee_c,nz))
          allocate(ebuffer2(ise_c:iee_c, jse_c:jee_c,nz))
       else
          allocate(ebuffer(1,1,1))
          allocate(ebuffer2(1,1,1))
       endif

       if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
          allocate(sbuffer(iss_c:ies_c, jss_c:jes_c,nz))
          allocate(sbuffer2(iss_c:ies_c, jss_c:jes_c,nz))
       else
          allocate(sbuffer(1,1,1))
          allocate(sbuffer2(1,1,1))
       endif

       if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
          allocate(nbuffer(isn_c:ien_c, jsn_c:jen_c,nz))
          allocate(nbuffer2(isn_c:ien_c, jsn_c:jen_c,nz))
       else
          allocate(nbuffer(1,1,1))
          allocate(nbuffer2(1,1,1))
       endif
       ebuffer = 0; ebuffer2 = 0
       wbuffer = 0; wbuffer2 = 0
       sbuffer = 0; sbuffer2 = 0
       nbuffer = 0; nbuffer2 = 0

       call mpp_update_nest_fine(x, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=l)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/),&
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbuffer2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r4_kind, 0.0_r4_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(wbuffer, wbuffer2, trim(type2)//' west buffer coarse to fine scalar')

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/),&
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbuffer2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r4_kind, 0.0_r4_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(sbuffer, sbuffer2, trim(type2)//' south buffer coarse to fine scalar')

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/),&
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebuffer2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r4_kind, 0.0_r4_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(ebuffer, ebuffer2, trim(type2)//' east buffer coarse to fine scalar')

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/),&
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbuffer2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, 0, 0, iadd_coarse,jadd_coarse,&
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r4_kind, 0.0_r4_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(nbuffer, nbuffer2, trim(type2)//' north buffer coarse to fine scalar')
       endif
       if(is_fine_pe) then
          deallocate(wbuffer, ebuffer, sbuffer, nbuffer)
          deallocate(wbuffer2, ebuffer2, sbuffer2, nbuffer2)
       endif
       deallocate(x)

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine BGRID scalar pair update
       !
       !---------------------------------------------------------------------------
       shift = 1
       !--- first check the index is correct or not
       if(is_fine_pe) then
          !--- The index from nest domain
          call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine)
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, l,&
                                 position=CORNER)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, l,&
                                 position=CORNER)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH,l,&
                                 position=CORNER)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH,l,&
                                 position=CORNER)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, l,&
                                 position=CORNER)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, l,&
                                 position=CORNER)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH,l,&
                                 position=CORNER)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH,l,&
                                 position=CORNER)

          !-- The assumed index
          isw_fx2 = 0; iew_fx2 = -1; jsw_fx2 = 0; jew_fx2 = -1
          isw_cx2 = 0; iew_cx2 = -1; jsw_cx2 = 0; jew_cx2 = -1
          ise_fx2 = 0; iee_fx2 = -1; jse_fx2 = 0; jee_fx2 = -1
          ise_cx2 = 0; iee_cx2 = -1; jse_cx2 = 0; jee_cx2 = -1
          iss_fx2 = 0; ies_fx2 = -1; jss_fx2 = 0; jes_fx2 = -1
          iss_cx2 = 0; ies_cx2 = -1; jss_cx2 = 0; jes_cx2 = -1
          isn_fx2 = 0; ien_fx2 = -1; jsn_fx2 = 0; jen_fx2 = -1
          isn_cx2 = 0; ien_cx2 = -1; jsn_cx2 = 0; jen_cx2 = -1
          isw_fy2 = 0; iew_fy2 = -1; jsw_fy2 = 0; jew_fy2 = -1
          isw_cy2 = 0; iew_cy2 = -1; jsw_cy2 = 0; jew_cy2 = -1
          ise_fy2 = 0; iee_fy2 = -1; jse_fy2 = 0; jee_fy2 = -1
          ise_cy2 = 0; iee_cy2 = -1; jse_cy2 = 0; jee_cy2 = -1
          iss_fy2 = 0; ies_fy2 = -1; jss_fy2 = 0; jes_fy2 = -1
          iss_cy2 = 0; ies_cy2 = -1; jss_cy2 = 0; jes_cy2 = -1
          isn_fy2 = 0; ien_fy2 = -1; jsn_fy2 = 0; jen_fy2 = -1
          isn_cy2 = 0; ien_cy2 = -1; jsn_cy2 = 0; jen_cy2 = -1

          !--- west
          if( isc_fine == 1 ) then
             isw_fx2 = isd_fine
             iew_fx2 = isc_fine - 1
             jsw_fx2 = jsd_fine
             jew_fx2 = jed_fine + shift
             isw_cx2 = istart_coarse(my_fine_id)-whalo
             iew_cx2 = istart_coarse(my_fine_id)
             jsw_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
             isw_fy2 = isd_fine
             iew_fy2 = isc_fine - 1
             jsw_fy2 = jsd_fine
             jew_fy2 = jed_fine + shift
             isw_cy2 = istart_coarse(my_fine_id)-whalo
             iew_cy2 = istart_coarse(my_fine_id)
             jsw_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- east
          if( iec_fine == nx_fine ) then
             ise_fx2 = iec_fine+1+shift
             iee_fx2 = ied_fine + shift
             jse_fx2 = jsd_fine
             jee_fx2 = jed_fine + shift
             ise_cx2 = iend_coarse(my_fine_id)+shift
             iee_cx2 = iend_coarse(my_fine_id)+ehalo+shift
             jse_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
             ise_fy2 = iec_fine+1 + shift
             iee_fy2 = ied_fine + shift
             jse_fy2 = jsd_fine
             jee_fy2 = jed_fine + shift
             ise_cy2 = iend_coarse(my_fine_id) + shift
             iee_cy2 = iend_coarse(my_fine_id)+ehalo + shift
             jse_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_fx2 = isd_fine
             ies_fx2 = ied_fine + shift
             jss_fx2 = jsd_fine
             jes_fx2 = jsc_fine - 1
             iss_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jss_cx2 = jstart_coarse(my_fine_id)-shalo
             jes_cx2 = jstart_coarse(my_fine_id)
             iss_fy2 = isd_fine
             ies_fy2 = ied_fine + shift
             jss_fy2 = jsd_fine
             jes_fy2 = jsc_fine - 1
             iss_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jss_cy2 = jstart_coarse(my_fine_id)-shalo
             jes_cy2 = jstart_coarse(my_fine_id)
          endif
          !--- north
          if( jec_fine == ny_fine ) then
             isn_fx2 = isd_fine
             ien_fx2 = ied_fine + shift
             jsn_fx2 = jec_fine+1 + shift
             jen_fx2 = jed_fine + shift
             isn_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jsn_cx2 = jend_coarse(my_fine_id) + shift
             jen_cx2 = jend_coarse(my_fine_id)+nhalo + shift
             isn_fy2 = isd_fine
             ien_fy2 = ied_fine + shift
             jsn_fy2 = jec_fine+1 + shift
             jen_fy2 = jed_fine + shift
             isn_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jsn_cy2 = jend_coarse(my_fine_id) + shift
             jen_cy2 = jend_coarse(my_fine_id)+nhalo + shift
          endif

          if( isw_fx .NE. isw_fx2 .OR. iew_fx .NE. iew_fx2 .OR. jsw_fx .NE. jsw_fx2 .OR. jew_fx .NE. jew_fx2 .OR. &
               isw_cx .NE. isw_cx2 .OR. iew_cx .NE. iew_cx2 .OR. jsw_cx .NE. jsw_cx2 .OR. jew_cx .NE. jew_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine BGRID X")
          endif
          if( ise_fx .NE. ise_fx2 .OR. iee_fx .NE. iee_fx2 .OR. jse_fx .NE. jse_fx2 .OR. jee_fx .NE. jee_fx2 .OR. &
               ise_cx .NE. ise_cx2 .OR. iee_cx .NE. iee_cx2 .OR. jse_cx .NE. jse_cx2 .OR. jee_cx .NE. jee_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine BGRID X")
          endif
          if( iss_fx .NE. iss_fx2 .OR. ies_fx .NE. ies_fx2 .OR. jss_fx .NE. jss_fx2 .OR. jes_fx .NE. jes_fx2 .OR. &
               iss_cx .NE. iss_cx2 .OR. ies_cx .NE. ies_cx2 .OR. jss_cx .NE. jss_cx2 .OR. jes_cx .NE. jes_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine BGRID X")
          endif
          if( isn_fx .NE. isn_fx2 .OR. ien_fx .NE. ien_fx2 .OR. jsn_fx .NE. jsn_fx2 .OR. jen_fx .NE. jen_fx2 .OR. &
               isn_cx .NE. isn_cx2 .OR. ien_cx .NE. ien_cx2 .OR. jsn_cx .NE. jsn_cx2 .OR. jen_cx .NE. jen_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine BGRID X")
          endif

          if( isw_fy .NE. isw_fy2 .OR. iew_fy .NE. iew_fy2 .OR. jsw_fy .NE. jsw_fy2 .OR. jew_fy .NE. jew_fy2 .OR. &
               isw_cy .NE. isw_cy2 .OR. iew_cy .NE. iew_cy2 .OR. jsw_cy .NE. jsw_cy2 .OR. jew_cy .NE. jew_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine BGRID Y")
          endif
          if( ise_fy .NE. ise_fy2 .OR. iee_fy .NE. iee_fy2 .OR. jse_fy .NE. jse_fy2 .OR. jee_fy .NE. jee_fy2 .OR. &
               ise_cy .NE. ise_cy2 .OR. iee_cy .NE. iee_cy2 .OR. jse_cy .NE. jse_cy2 .OR. jee_cy .NE. jee_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine BGRID Y")
          endif
          if( iss_fy .NE. iss_fy2 .OR. ies_fy .NE. ies_fy2 .OR. jss_fy .NE. jss_fy2 .OR. jes_fy .NE. jes_fy2 .OR. &
               iss_cy .NE. iss_cy2 .OR. ies_cy .NE. ies_cy2 .OR. jss_cy .NE. jss_cy2 .OR. jes_cy .NE. jes_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine BGRID Y")
          endif
          if( isn_fy .NE. isn_fy2 .OR. ien_fy .NE. ien_fy2 .OR. jsn_fy .NE. jsn_fy2 .OR. jen_fy .NE. jen_fy2 .OR. &
               isn_cy .NE. isn_cy2 .OR. ien_cy .NE. ien_cy2 .OR. jsn_cy .NE. jsn_cy2 .OR. jen_cy .NE. jen_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine BGRID Y")
          endif
       endif

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse+shift, nz))
          allocate(y(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = 1e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse+shift
                   y(i,j,k) = 2e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine+shift, jsd_fine:jed_fine+shift, nz))
          allocate(y(isd_fine:ied_fine+shift, jsd_fine:jed_fine+shift, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 1e-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine+shift
                   y(i,j,k) = i*1.e+6 + j*1.e+3 + k + 2e-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          if( iew_cx .GE. isw_cx .AND. jew_cx .GE. jsw_cx ) then
             allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
             allocate(wbufferx2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery2(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
          else
             allocate(wbufferx(1,1,1))
             allocate(wbuffery(1,1,1))
             allocate(wbufferx2(1,1,1))
             allocate(wbuffery2(1,1,1))
          endif
          if( iee_cx .GE. ise_cx .AND. jee_cx .GE. jse_cx ) then
             allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy,nz))
             allocate(ebufferx2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery2(ise_cy:iee_cy, jse_cy:jee_cy,nz))
          else
             allocate(ebufferx(1,1,1))
             allocate(ebuffery(1,1,1))
             allocate(ebufferx2(1,1,1))
             allocate(ebuffery2(1,1,1))
          endif
          if( ies_cx .GE. iss_cx .AND. jes_cx .GE. jss_cx ) then
             allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy,nz))
             allocate(sbufferx2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery2(iss_cy:ies_cy, jss_cy:jes_cy,nz))
          else
             allocate(sbufferx(1,1,1))
             allocate(sbuffery(1,1,1))
             allocate(sbufferx2(1,1,1))
             allocate(sbuffery2(1,1,1))
          endif
          if( ien_cx .GE. isn_cx .AND. jen_cx .GE. jsn_cx ) then
             allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
             allocate(nbufferx2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery2(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
          else
             allocate(nbufferx(1,1,1))
             allocate(nbuffery(1,1,1))
             allocate(nbufferx2(1,1,1))
             allocate(nbuffery2(1,1,1))
          endif
          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=BGRID_NE, flags=SCALAR_PAIR)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/),&
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, 1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/),&
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, 1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/),&
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, &
                  1E3_r4_kind, 2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, &
                  2E3_r4_kind, 1E3_r4_kind, 1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, &
                  1E3_r4_kind, 2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, &
                  2E3_r4_kind, 1E3_r4_kind, 1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine BGRID scalar pair Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine BGRID scalar pair Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine BGRID scalar pair Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine BGRID scalar pair X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine BGRID scalar pair Y')
       endif
       if(allocated(x)) deallocate(x)
       if(allocated(y)) deallocate(y)
       if(is_fine_pe) then
          deallocate(wbufferx, ebufferx, sbufferx, nbufferx)
          deallocate(wbufferx2, ebufferx2, sbufferx2, nbufferx2)
          deallocate(wbuffery, ebuffery, sbuffery, nbuffery)
          deallocate(wbuffery2, ebuffery2, sbuffery2, nbuffery2)
       endif

       !---------------------------------------------------------------------------
       !
       !                 Coarse to Fine scalar field, position = CORNER
       !
       !---------------------------------------------------------------------------
       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine+shift, jsd_fine:jed_fine+shift, nz))
          x = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             allocate(wbuffer(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffer2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
          else
             allocate(wbuffer(1,1,1))
             allocate(wbuffer2(1,1,1))
          endif
          wbuffer = 0; wbuffer2 = 0

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             allocate(ebuffer(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffer2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
          else
             allocate(ebuffer(1,1,1))
             allocate(ebuffer2(1,1,1))
          endif
          ebuffer = 0; ebuffer2 = 0

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             allocate(sbuffer(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffer2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
          else
             allocate(sbuffer(1,1,1))
             allocate(sbuffer2(1,1,1))
          endif
          sbuffer = 0; sbuffer2 = 0

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             allocate(nbuffer(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffer2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
          else
             allocate(nbuffer(1,1,1))
             allocate(nbuffer2(1,1,1))
          endif
          nbuffer = 0; nbuffer2 = 0

       endif

       call mpp_update_nest_fine(x, nest_domain, wbuffer, sbuffer, ebuffer, nbuffer, nest_level=l, position=CORNER)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), &
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbuffer2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r4_kind, &
                  0.0_r4_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(wbuffer, wbuffer2, trim(type2)//' west buffer coarse to fine scalar CORNER')

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), &
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbuffer2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0_r4_kind, 0.0_r4_kind, &
                  1, 1, nx, ny)
          endif
          call compare_checksums(sbuffer, sbuffer2, trim(type2)//' south buffer coarse to fine scalar CORNER')

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), &
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebuffer2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, &
                  0.0_r4_kind, 0.0_r4_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(ebuffer, ebuffer2, trim(type2)//' east buffer coarse to fine scalar CORNER')

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbuffer2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse,&
                  0.0_r4_kind, 0.0_r4_kind, 1, 1, nx, ny)
          endif
          call compare_checksums(nbuffer, nbuffer2, trim(type2)//' north buffer coarse to fine scalar CORNER')

       endif
       if(is_fine_pe) then
          deallocate(wbuffer, ebuffer, sbuffer, nbuffer)
          deallocate(wbuffer2, ebuffer2, sbuffer2, nbuffer2)
       endif
       deallocate(x)

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine CGRID scalar pair update
       !
       !---------------------------------------------------------------------------
       shift = 1
       !--- first check the index is correct or not
       if(is_fine_pe) then
          !--- The index from nest domain
          call mpp_get_compute_domain(domain_fine, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain_fine, isd_fine, ied_fine, jsd_fine, jed_fine)
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, l,&
                                 position=EAST)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, l,&
                                 position=EAST)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH,l,&
                                 position=EAST)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH,l,&
                                 position=EAST)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, l,&
                                 position=NORTH)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, l,&
                                 position=NORTH)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH,l,&
                                 position=NORTH)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH,l,&
                                 position=NORTH)

          !-- The assumed index
          isw_fx2 = 0; iew_fx2 = -1; jsw_fx2 = 0; jew_fx2 = -1
          isw_cx2 = 0; iew_cx2 = -1; jsw_cx2 = 0; jew_cx2 = -1
          ise_fx2 = 0; iee_fx2 = -1; jse_fx2 = 0; jee_fx2 = -1
          ise_cx2 = 0; iee_cx2 = -1; jse_cx2 = 0; jee_cx2 = -1
          iss_fx2 = 0; ies_fx2 = -1; jss_fx2 = 0; jes_fx2 = -1
          iss_cx2 = 0; ies_cx2 = -1; jss_cx2 = 0; jes_cx2 = -1
          isn_fx2 = 0; ien_fx2 = -1; jsn_fx2 = 0; jen_fx2 = -1
          isn_cx2 = 0; ien_cx2 = -1; jsn_cx2 = 0; jen_cx2 = -1
          isw_fy2 = 0; iew_fy2 = -1; jsw_fy2 = 0; jew_fy2 = -1
          isw_cy2 = 0; iew_cy2 = -1; jsw_cy2 = 0; jew_cy2 = -1
          ise_fy2 = 0; iee_fy2 = -1; jse_fy2 = 0; jee_fy2 = -1
          ise_cy2 = 0; iee_cy2 = -1; jse_cy2 = 0; jee_cy2 = -1
          iss_fy2 = 0; ies_fy2 = -1; jss_fy2 = 0; jes_fy2 = -1
          iss_cy2 = 0; ies_cy2 = -1; jss_cy2 = 0; jes_cy2 = -1
          isn_fy2 = 0; ien_fy2 = -1; jsn_fy2 = 0; jen_fy2 = -1
          isn_cy2 = 0; ien_cy2 = -1; jsn_cy2 = 0; jen_cy2 = -1

          !--- west
          if( isc_fine == 1 ) then
             isw_fx2 = isd_fine
             iew_fx2 = isc_fine - 1
             jsw_fx2 = jsd_fine
             jew_fx2 = jed_fine
             isw_cx2 = istart_coarse(my_fine_id)-whalo
             iew_cx2 = istart_coarse(my_fine_id)
             jsw_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
             isw_fy2 = isd_fine
             iew_fy2 = isc_fine - 1
             jsw_fy2 = jsd_fine
             jew_fy2 = jed_fine + shift
             isw_cy2 = istart_coarse(my_fine_id)-whalo
             iew_cy2 = istart_coarse(my_fine_id)
             jsw_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- east
          if( iec_fine == nx_fine ) then
             ise_fx2 = iec_fine+1+shift
             iee_fx2 = ied_fine + shift
             jse_fx2 = jsd_fine
             jee_fx2 = jed_fine
             ise_cx2 = iend_coarse(my_fine_id)+shift
             iee_cx2 = iend_coarse(my_fine_id)+ehalo+shift
             jse_cx2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
             ise_fy2 = iec_fine+1
             iee_fy2 = ied_fine
             jse_fy2 = jsd_fine
             jee_fy2 = jed_fine + shift
             ise_cy2 = iend_coarse(my_fine_id)
             iee_cy2 = iend_coarse(my_fine_id)+ehalo
             jse_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) &
                     & + nhalo + shift
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_fx2 = isd_fine
             ies_fx2 = ied_fine + shift
             jss_fx2 = jsd_fine
             jes_fx2 = jsc_fine - 1
             iss_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jss_cx2 = jstart_coarse(my_fine_id)-shalo
             jes_cx2 = jstart_coarse(my_fine_id)
             iss_fy2 = isd_fine
             ies_fy2 = ied_fine
             jss_fy2 = jsd_fine
             jes_fy2 = jsc_fine - 1
             iss_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jss_cy2 = jstart_coarse(my_fine_id)-shalo
             jes_cy2 = jstart_coarse(my_fine_id)
          endif
          !--- north
          if( jec_fine == ny_fine ) then
             isn_fx2 = isd_fine
             ien_fx2 = ied_fine + shift
             jsn_fx2 = jec_fine+1
             jen_fx2 = jed_fine
             isn_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) &
                     & + ehalo + shift
             jsn_cx2 = jend_coarse(my_fine_id)
             jen_cx2 = jend_coarse(my_fine_id)+nhalo
             isn_fy2 = isd_fine
             ien_fy2 = ied_fine
             jsn_fy2 = jec_fine+1 + shift
             jen_fy2 = jed_fine + shift
             isn_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo
             jsn_cy2 = jend_coarse(my_fine_id) + shift
             jen_cy2 = jend_coarse(my_fine_id)+nhalo + shift
          endif

          if( isw_fx .NE. isw_fx2 .OR. iew_fx .NE. iew_fx2 .OR. jsw_fx .NE. jsw_fx2 .OR. jew_fx .NE. jew_fx2 .OR. &
               isw_cx .NE. isw_cx2 .OR. iew_cx .NE. iew_cx2 .OR. jsw_cx .NE. jsw_cx2 .OR. jew_cx .NE. jew_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine CGRID X")
          endif
          if( ise_fx .NE. ise_fx2 .OR. iee_fx .NE. iee_fx2 .OR. jse_fx .NE. jse_fx2 .OR. jee_fx .NE. jee_fx2 .OR. &
               ise_cx .NE. ise_cx2 .OR. iee_cx .NE. iee_cx2 .OR. jse_cx .NE. jse_cx2 .OR. jee_cx .NE. jee_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine CGRID X")
          endif
          if( iss_fx .NE. iss_fx2 .OR. ies_fx .NE. ies_fx2 .OR. jss_fx .NE. jss_fx2 .OR. jes_fx .NE. jes_fx2 .OR. &
               iss_cx .NE. iss_cx2 .OR. ies_cx .NE. ies_cx2 .OR. jss_cx .NE. jss_cx2 .OR. jes_cx .NE. jes_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine CGRID X")
          endif
          if( isn_fx .NE. isn_fx2 .OR. ien_fx .NE. ien_fx2 .OR. jsn_fx .NE. jsn_fx2 .OR. jen_fx .NE. jen_fx2 .OR. &
               isn_cx .NE. isn_cx2 .OR. ien_cx .NE. ien_cx2 .OR. jsn_cx .NE. jsn_cx2 .OR. jen_cx .NE. jen_cx2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine CGRID X")
          endif

          if( isw_fy .NE. isw_fy2 .OR. iew_fy .NE. iew_fy2 .OR. jsw_fy .NE. jsw_fy2 .OR. jew_fy .NE. jew_fy2 .OR. &
               isw_cy .NE. isw_cy2 .OR. iew_cy .NE. iew_cy2 .OR. jsw_cy .NE. jsw_cy2 .OR. jew_cy .NE. jew_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: west buffer index mismatch for coarse to fine CGRID Y")
          endif
          if( ise_fy .NE. ise_fy2 .OR. iee_fy .NE. iee_fy2 .OR. jse_fy .NE. jse_fy2 .OR. jee_fy .NE. jee_fy2 .OR. &
               ise_cy .NE. ise_cy2 .OR. iee_cy .NE. iee_cy2 .OR. jse_cy .NE. jse_cy2 .OR. jee_cy .NE. jee_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: east buffer index mismatch for coarse to fine CGRID Y")
          endif
          if( iss_fy .NE. iss_fy2 .OR. ies_fy .NE. ies_fy2 .OR. jss_fy .NE. jss_fy2 .OR. jes_fy .NE. jes_fy2 .OR. &
               iss_cy .NE. iss_cy2 .OR. ies_cy .NE. ies_cy2 .OR. jss_cy .NE. jss_cy2 .OR. jes_cy .NE. jes_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: south buffer index mismatch for coarse to fine CGRID Y")
          endif
          if( isn_fy .NE. isn_fy2 .OR. ien_fy .NE. ien_fy2 .OR. jsn_fy .NE. jsn_fy2 .OR. jen_fy .NE. jen_fy2 .OR. &
               isn_cy .NE. isn_cy2 .OR. ien_cy .NE. ien_cy2 .OR. jsn_cy .NE. jsn_cy2 .OR. jen_cy .NE. jen_cy2 ) then
             call mpp_error(FATAL, "test_mpp_domains: north buffer index mismatch for coarse to fine CGRID Y")
          endif
       endif

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(x(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(y(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = 1e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y(i,j,k) = 2e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
       else
          allocate(x(isd_fine:ied_fine+shift, jsd_fine:jed_fine, nz))
          allocate(y(isd_fine:ied_fine, jsd_fine:jed_fine+shift, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 1e-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine
                   y(i,j,k) = i*1.e+6 + j*1.e+3 + k + 2e-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          if( iew_cx .GE. isw_cx .AND. jew_cx .GE. jsw_cx ) then
             allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
             allocate(wbufferx2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery2(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
          else
             allocate(wbufferx(1,1,1))
             allocate(wbuffery(1,1,1))
             allocate(wbufferx2(1,1,1))
             allocate(wbuffery2(1,1,1))
          endif
          if( iee_cx .GE. ise_cx .AND. jee_cx .GE. jse_cx ) then
             allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy,nz))
             allocate(ebufferx2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery2(ise_cy:iee_cy, jse_cy:jee_cy,nz))
          else
             allocate(ebufferx(1,1,1))
             allocate(ebuffery(1,1,1))
             allocate(ebufferx2(1,1,1))
             allocate(ebuffery2(1,1,1))
          endif
          if( ies_cx .GE. iss_cx .AND. jes_cx .GE. jss_cx ) then
             allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy,nz))
             allocate(sbufferx2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery2(iss_cy:ies_cy, jss_cy:jes_cy,nz))
          else
             allocate(sbufferx(1,1,1))
             allocate(sbuffery(1,1,1))
             allocate(sbufferx2(1,1,1))
             allocate(sbuffery2(1,1,1))
          endif
          if( ien_cx .GE. isn_cx .AND. jen_cx .GE. jsn_cx ) then
             allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
             allocate(nbufferx2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery2(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
          else
             allocate(nbufferx(1,1,1))
             allocate(nbuffery(1,1,1))
             allocate(nbufferx2(1,1,1))
             allocate(nbuffery2(1,1,1))
          endif
          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=CGRID_NE, flags=SCALAR_PAIR)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), &
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  1E3_r4_kind, 2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  2E3_r4_kind, 1E3_r4_kind, 1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), &
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  1E3_r4_kind, 2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  2E3_r4_kind, 1E3_r4_kind, 1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), &
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, &
                  1E3_r4_kind, 2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  2E3_r4_kind, 1E3_r4_kind, 1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, &
                  1E3_r4_kind, 2E3_r4_kind, 1, 1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, &
                  iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, &
                  2E3_r4_kind, 1E3_r4_kind, 1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine CGRID scalar pair Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine CGRID scalar pair Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine CGRID scalar pair Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine CGRID scalar pair X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine CGRID scalar pair Y')
       endif

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine CGRID vector update
       !
       !---------------------------------------------------------------------------
       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   x(i,j,k) = 1e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   y(i,j,k) = 2e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
       else
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine+shift
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 1e-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine
                   y(i,j,k) = i*1.e+6 + j*1.e+3 + k + 2e-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=CGRID_NE)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), &
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, 2E3_r4_kind, &
                  1, -1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, 1E3_r4_kind, &
                  -1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), &
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, -1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, -1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), &
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, -1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, -1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, -1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, -1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine CGRID vector X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine CGRID vector Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine CGRID vector X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine CGRID vector Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine CGRID vector X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine CGRID vector Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine CGRID vector X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine CGRID vector Y')
       endif

       if(allocated(x)) deallocate(x)
       if(allocated(y)) deallocate(y)
       if(is_fine_pe) then
          deallocate(wbufferx, ebufferx, sbufferx, nbufferx)
          deallocate(wbufferx2, ebufferx2, sbufferx2, nbufferx2)
          deallocate(wbuffery, ebuffery, sbuffery, nbuffery)
          deallocate(wbuffery2, ebuffery2, sbuffery2, nbuffery2)
       endif

       !---------------------------------------------------------------------------
       !
       !                    coarse to fine DGRID vector update
       !
       !---------------------------------------------------------------------------
       shift = 1

       if(is_coarse_pe) then
          call mpp_get_compute_domain(domain_coarse, isc_coarse, iec_coarse, jsc_coarse, jec_coarse)
          call mpp_get_data_domain(domain_coarse, isd_coarse, ied_coarse, jsd_coarse, jed_coarse)
          allocate(y(isd_coarse:ied_coarse+shift, jsd_coarse:jed_coarse, nz))
          allocate(x(isd_coarse:ied_coarse, jsd_coarse:jed_coarse+shift, nz))
          x = 0
          y = 0
          tile = my_tile_id(1)
          do k = 1, nz
             do j = jsc_coarse, jec_coarse+shift
                do i = isc_coarse, iec_coarse
                   x(i,j,k) = 1e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_coarse, jec_coarse
                do i = isc_coarse, iec_coarse+shift
                   y(i,j,k) = 2e3 + tile + i*1.e-3 + j*1.e-6 + k*1.e-9
                enddo
             enddo
          enddo
       else
          allocate(y(isd_fine:ied_fine+shift, jsd_fine:jed_fine, nz))
          allocate(x(isd_fine:ied_fine, jsd_fine:jed_fine+shift, nz))
          x = 0
          y = 0
          do k = 1, nz
             do j = jsc_fine, jec_fine+shift
                do i = isc_fine, iec_fine
                   x(i,j,k) = i*1.e+6 + j*1.e+3 + k + 1e-3
                enddo
             enddo
          enddo
          do k = 1, nz
             do j = jsc_fine, jec_fine
                do i = isc_fine, iec_fine+shift
                   y(i,j,k) = i*1.e+6 + j*1.e+3 + k + 2e-3
                enddo
             enddo
          enddo
       endif

       if(is_fine_pe) then
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST,l,&
                                 position=NORTH)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST,l,&
                                 position=NORTH)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH,l,&
                                 position=NORTH)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH,l,&
                                 position=NORTH)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST,l,&
                                 position=EAST)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST,l,&
                                 position=EAST)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH,l,&
                                 position=EAST)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH,l,&
                                 position=EAST)

          if( iew_cx .GE. isw_cx .AND. jew_cx .GE. jsw_cx ) then
             allocate(wbufferx(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
             allocate(wbufferx2(isw_cx:iew_cx, jsw_cx:jew_cx,nz))
             allocate(wbuffery2(isw_cy:iew_cy, jsw_cy:jew_cy,nz))
          else
             allocate(wbufferx(1,1,1))
             allocate(wbuffery(1,1,1))
             allocate(wbufferx2(1,1,1))
             allocate(wbuffery2(1,1,1))
          endif
          if( iee_cx .GE. ise_cx .AND. jee_cx .GE. jse_cx ) then
             allocate(ebufferx(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery(ise_cy:iee_cy, jse_cy:jee_cy,nz))
             allocate(ebufferx2(ise_cx:iee_cx, jse_cx:jee_cx,nz))
             allocate(ebuffery2(ise_cy:iee_cy, jse_cy:jee_cy,nz))
          else
             allocate(ebufferx(1,1,1))
             allocate(ebuffery(1,1,1))
             allocate(ebufferx2(1,1,1))
             allocate(ebuffery2(1,1,1))
          endif
          if( ies_cx .GE. iss_cx .AND. jes_cx .GE. jss_cx ) then
             allocate(sbufferx(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery(iss_cy:ies_cy, jss_cy:jes_cy,nz))
             allocate(sbufferx2(iss_cx:ies_cx, jss_cx:jes_cx,nz))
             allocate(sbuffery2(iss_cy:ies_cy, jss_cy:jes_cy,nz))
          else
             allocate(sbufferx(1,1,1))
             allocate(sbuffery(1,1,1))
             allocate(sbufferx2(1,1,1))
             allocate(sbuffery2(1,1,1))
          endif
          if( ien_cx .GE. isn_cx .AND. jen_cx .GE. jsn_cx ) then
             allocate(nbufferx(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
             allocate(nbufferx2(isn_cx:ien_cx, jsn_cx:jen_cx,nz))
             allocate(nbuffery2(isn_cy:ien_cy, jsn_cy:jen_cy,nz))
          else
             allocate(nbufferx(1,1,1))
             allocate(nbuffery(1,1,1))
             allocate(nbufferx2(1,1,1))
             allocate(nbuffery2(1,1,1))
          endif

          wbufferx = 0; wbufferx2 = 0
          wbuffery = 0; wbuffery2 = 0
          sbufferx = 0; sbufferx2 = 0
          sbuffery = 0; sbuffery2 = 0
          ebufferx = 0; ebufferx2 = 0
          ebuffery = 0; ebuffery2 = 0
          nbufferx = 0; nbufferx2 = 0
          nbuffery = 0; nbuffery2 = 0
       endif
       call mpp_update_nest_fine(x, y, nest_domain, wbufferx, wbuffery, sbufferx, sbuffery, ebufferx, ebuffery, &
            nbufferx, nbuffery, nest_level=l, gridtype=DGRID_NE)

       !--- compare with the assumed value.
       if( is_fine_pe ) then
          call mpp_set_current_pelist(my_pelist_fine)
          if( iew_c .GE. isw_c .AND. jew_c .GE. jsw_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), &
                  (/jew_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, -1, nx, ny)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, -1, 1, nx, ny)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), &
                  (/jes_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, -1, nx, ny)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, -1, 1, nx, ny)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), &
                  (/jee_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, -1, nx, ny)
             call fill_nest_data(ebuffery2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, -1, 1, nx, ny)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), &
                  (/jen_c/), nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, &
                  js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 1E3_r4_kind, &
                  2E3_r4_kind, 1, -1, nx, ny)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse, &
                  jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2E3_r4_kind, &
                  1E3_r4_kind, -1, 1, nx, ny)
          endif

          call compare_checksums(wbufferx, wbufferx2, trim(type2)//' west buffer coarse to fine DGRID vector X')
          call compare_checksums(wbuffery, wbuffery2, trim(type2)//' west buffer coarse to fine DGRID vector Y')
          call compare_checksums(sbufferx, sbufferx2, trim(type2)//' south buffer coarse to fine DGRID vector X')
          call compare_checksums(sbuffery, sbuffery2, trim(type2)//' south buffer coarse to fine DGRID vector Y')
          call compare_checksums(ebufferx, ebufferx2, trim(type2)//' east buffer coarse to fine DGRID vector X')
          call compare_checksums(ebuffery, ebuffery2, trim(type2)//' east buffer coarse to fine DGRID vector Y')
          call compare_checksums(nbufferx, nbufferx2, trim(type2)//' north buffer coarse to fine DGRID vector X')
          call compare_checksums(nbuffery, nbuffery2, trim(type2)//' north buffer coarse to fine DGRID vector Y')
       endif



       if(allocated(x)) deallocate(x)
       if(allocated(y)) deallocate(y)
       if(is_fine_pe) then
          deallocate(wbufferx, ebufferx, sbufferx, nbufferx)
          deallocate(wbufferx2, ebufferx2, sbufferx2, nbufferx2)
          deallocate(wbuffery, ebuffery, sbuffery, nbuffery)
          deallocate(wbuffery2, ebuffery2, sbuffery2, nbuffery2)
       endif
       endif

       deallocate(my_pelist, my_pelist_fine)
       call mpp_set_current_pelist()

    enddo

    call mpp_set_current_pelist(pelist)
    call mpp_sync(pelist)
    deallocate(pelist)

  end subroutine test_update_nest_domain_r4

end program
