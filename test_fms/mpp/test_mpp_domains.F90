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
program test_mpp_domains
  use mpp_mod,         only : FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_node, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_mod,         only : input_nml_file
  use mpp_mod,         only : mpp_get_current_pelist, mpp_broadcast
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE, DGRID_NE, AGRID
  use mpp_domains_mod, only : FOLD_SOUTH_EDGE, FOLD_NORTH_EDGE, FOLD_WEST_EDGE, FOLD_EAST_EDGE
  use mpp_domains_mod, only : MPP_DOMAIN_TIME, CYCLIC_GLOBAL_DOMAIN, NUPDATE,EUPDATE, XUPDATE, YUPDATE, SCALAR_PAIR
  use mpp_domains_mod, only : domain1D, domain2D, DomainCommunicator2D, BITWISE_EFP_SUM
  use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_data_domain, mpp_domains_set_stack_size
  use mpp_domains_mod, only : mpp_global_field, mpp_global_sum, mpp_global_max, mpp_global_min
  use mpp_domains_mod, only : mpp_domains_init, mpp_domains_exit, mpp_broadcast_domain
  use mpp_domains_mod, only : mpp_update_domains, mpp_check_field, mpp_redistribute, mpp_get_memory_domain
  use mpp_domains_mod, only : mpp_define_layout, mpp_define_domains, mpp_modify_domain
  use mpp_domains_mod, only : mpp_get_neighbor_pe, mpp_define_mosaic, mpp_nullify_domain_list
  use mpp_domains_mod, only : NORTH, NORTH_EAST, EAST, SOUTH_EAST, CORNER, CENTER
  use mpp_domains_mod, only : SOUTH, SOUTH_WEST, WEST, NORTH_WEST, mpp_define_mosaic_pelist
  use mpp_domains_mod, only : mpp_get_global_domain, ZERO, NINETY, MINUS_NINETY
  use mpp_domains_mod, only : mpp_get_boundary, mpp_start_update_domains, mpp_complete_update_domains
  use mpp_domains_mod, only : mpp_define_nest_domains, nest_domain_type
  use mpp_domains_mod, only : mpp_get_C2F_index, mpp_update_nest_fine
  use mpp_domains_mod, only : mpp_get_F2C_index, mpp_update_nest_coarse
  use mpp_domains_mod, only : mpp_get_nest_coarse_domain, mpp_get_nest_fine_domain
  use mpp_domains_mod, only : mpp_is_nest_fine, mpp_is_nest_coarse
  use mpp_domains_mod, only : mpp_get_nest_pelist, mpp_get_nest_npes
  use mpp_domains_mod, only : mpp_get_nest_fine_pelist, mpp_get_nest_fine_npes
  use mpp_domains_mod, only : mpp_get_domain_shift, EDGEUPDATE, mpp_deallocate_domain
  use mpp_domains_mod, only : mpp_group_update_type, mpp_create_group_update
  use mpp_domains_mod, only : mpp_do_group_update, mpp_clear_group_update
  use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
  use mpp_domains_mod, only : WUPDATE, SUPDATE, mpp_get_compute_domains, NONSYMEDGEUPDATE
  use mpp_domains_mod, only : domainUG, mpp_define_unstruct_domain, mpp_get_UG_domain_tile_id
  use mpp_domains_mod, only : mpp_get_UG_compute_domain, mpp_pass_SG_to_UG, mpp_pass_UG_to_SG
  use mpp_domains_mod, only : mpp_get_ug_global_domain, mpp_global_field_ug, mpp_get_tile_id
  use mpp_memutils_mod, only : mpp_memuse_begin, mpp_memuse_end
  use fms_affinity_mod, only : fms_affinity_set


  implicit none
#include "../../include/fms_platform.h"
  integer :: pe, npes
  integer :: nx=128, ny=128, nz=40, stackmax=4000000
  integer :: unit=7
  integer :: stdunit = 6
  logical :: debug=.FALSE., opened

  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  integer :: x_cyclic_offset = 3   ! to be used in test_cyclic_offset
  integer :: y_cyclic_offset = -4  ! to be used in test_cyclic_offset
  character(len=32) :: warn_level = "fatal"
  integer :: wide_halo_x = 0, wide_halo_y = 0
  integer :: nx_cubic = 0, ny_cubic = 0
  logical :: test_nest = .false.
  logical :: test_performance = .false.
  logical :: test_interface = .true.
  logical :: test_edge_update = .false.
  logical :: test_nonsym_edge = .false.
  logical :: test_group = .false.
  logical :: test_cubic_grid_redistribute = .false.
  logical :: check_parallel = .FALSE.  ! when check_parallel set to false,
  logical :: test_get_nbr = .FALSE.
  logical :: test_boundary = .false.
  logical :: test_global_sum = .false.
  logical :: test_halosize_performance = .false.
  integer :: ensemble_size = 1
  integer :: layout_cubic(2) = (/0,0/)
  integer :: layout_tripolar(2) = (/0,0/)
  integer :: layout_ensemble(2) = (/0,0/)
  logical :: do_sleep = .false.
  integer :: num_iter = 1
  integer :: num_fields = 4

  logical :: mix_2D_3D = .false.
  logical :: test_subset = .false.
  logical :: test_unstruct = .false.
  integer :: nthreads = 1
  logical :: test_adjoint = .false.
  logical :: wide_halo = .false.

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
  character :: cyclic_nest(MAX_NCONTACT) = 'N'

  namelist / test_mpp_domains_nml / nx, ny, nz, stackmax, debug, mpes, check_parallel, &
                               whalo, ehalo, shalo, nhalo, x_cyclic_offset, y_cyclic_offset, &
                               warn_level, wide_halo_x, wide_halo_y, nx_cubic, ny_cubic, &
                               test_performance, test_interface, num_fields, do_sleep, num_iter, &
                               test_nest, num_nest, ntiles_nest_all, nest_level, tile_fine, tile_coarse, &
                               refine_ratio, istart_coarse, icount_coarse, jstart_coarse, jcount_coarse, &
                               extra_halo, npes_nest_tile, cyclic_nest, mix_2D_3D, test_get_nbr, &
                               test_edge_update, test_cubic_grid_redistribute, ensemble_size, &
                               layout_cubic, layout_ensemble, nthreads, test_boundary, &
                               layout_tripolar, test_group, test_global_sum, test_subset, test_unstruct, &
                               test_nonsym_edge, test_halosize_performance, test_adjoint, wide_halo
  integer :: i, j, k, n
  integer :: layout(2)
  integer :: id
  integer :: outunit, errunit, io_status
  integer :: omp_get_num_threads, omp_get_thread_num

  call mpp_memuse_begin()
  call mpp_init()

  outunit = stdout()
  errunit = stderr()
#ifdef INTERNAL_FILE_NML
  read (input_nml_file, test_mpp_domains_nml, iostat=io_status)
#else
  do
     inquire( unit=unit, opened=opened )
     if( .NOT.opened )exit
     unit = unit + 1
     if( unit.EQ.100 )call mpp_error( FATAL, 'Unable to locate unit number.' )
  end do
  open( unit=unit, file='input.nml', iostat=io_status )
  read( unit,test_mpp_domains_nml, iostat=io_status )
  close(unit)
#endif

  if (io_status > 0) then
     call mpp_error(FATAL,'=>test_mpp_domains: Error reading input.nml')
  endif

  select case(trim(warn_level))
  case("fatal")
     call mpp_set_warn_level(FATAL)
  case("warning")
     call mpp_set_warn_level(WARNING)
  case default
     call mpp_error(FATAL, "test_mpp_domains: warn_level should be fatal or warning")
  end select

  pe = mpp_pe()
  npes = mpp_npes()

  !--- initialize mpp domains
  if( (.not.debug) .and. test_nest ) then
      call mpp_domains_init()
  elseif( debug )then
      call mpp_domains_init(MPP_DEBUG)
  else
      call mpp_domains_init(MPP_DOMAIN_TIME)
  end if
  call mpp_domains_set_stack_size(stackmax)

!$  call omp_set_num_threads(nthreads)
!$OMP PARALLEL
!$  call fms_affinity_set("test_mpp_domains", .FALSE., omp_get_num_threads())
!$OMP END PARALLEL

  if( pe.EQ.mpp_root_pe() )print '(a,9i6)', 'npes, mpes, nx, ny, nz, whalo, ehalo, shalo, nhalo =', &
                           npes, mpes, nx, ny, nz, whalo, ehalo, shalo, nhalo
  call mpp_memuse_end("in the begining", outunit)

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

  if( test_nest .and. (num_nest>0) ) then
    if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_update_nest_domain <-------------------'
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
     call test_update_nest_domain('Cubic-Grid')
    if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_update_nest_domain <-------------------'
  endif

  if(test_subset) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_subset_update <-------------------'
      call test_subset_update()
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_subset_update <-------------------'
  endif

  if( test_halosize_performance ) then
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_halosize_performance <-------------------'
     call test_halosize_update( 'Folded-north' )
     call test_halosize_update( 'Folded-north symmetry' )
     call test_halosize_update( 'Cubic-Grid' )
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_halosize_performance <-------------------'
  endif

  if( test_edge_update ) then
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_edge_update <-------------------'
      call test_update_edge( 'Cyclic' )
      call test_update_edge( 'Folded-north' ) !includes vector field test
      call test_update_edge( 'Folded-north symmetry' )
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_edge_update <-------------------'
  endif

  if( test_nonsym_edge ) then
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_nonsym_edge <-------------------'
      call test_update_nonsym_edge( 'Folded-north' ) !includes vector field test
      call test_update_nonsym_edge( 'Folded-north symmetry' )
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_nonsym_edge <-------------------'
  endif

  if( test_performance) then
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_performance <-------------------'
      call update_domains_performance('Folded-north')
      call update_domains_performance('Cubic-Grid')
     if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_performance <-------------------'
  endif

  if( test_global_sum ) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_mpp_global_sum <-------------------'
      call test_mpp_global_sum('Folded-north')
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_mpp_global_sum <-------------------'
  endif

  if( test_cubic_grid_redistribute ) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling cubic_grid_redistribute <-------------------'
     call cubic_grid_redistribute()
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished cubic_grid_redistribute <-------------------'
  endif

  if(test_boundary) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_boundary <-------------------'
      call test_get_boundary('torus')
      call test_get_boundary('Four-Tile')
      call test_get_boundary('Cubic-Grid')
      call test_get_boundary('Folded-north')
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_boundary <-------------------'
  endif

! Adjoint Dot Test ------------------------------------------
  if (test_adjoint) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_adjoint <-------------------'
       call test_get_boundary_ad('Four-Tile')
       call test_halo_update_ad( 'Simple' )
       call test_global_reduce_ad( 'Simple')
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finished test_adjoint <-------------------'
  endif

  if( test_unstruct) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_unstruct <-------------------'
     call test_unstruct_update( 'Cubic-Grid' )
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_unstruct <-------------------'
  endif

  if( test_group) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_group <-------------------'
     call test_group_update( 'Folded-north' )
     call test_group_update( 'Cubic-Grid' )
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_group <-------------------'
  endif

  if( test_interface ) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_interface  <-------------------'
      call test_modify_domain()
!!$      call test_cyclic_offset('x_cyclic_offset')
!!$      call test_cyclic_offset('y_cyclic_offset')
!!$      call test_cyclic_offset('torus_x_offset')
!!$      call test_cyclic_offset('torus_y_offset')
      if(.not. wide_halo) call test_uniform_mosaic('Single-Tile')
      call test_uniform_mosaic('Folded-north mosaic') ! one-tile tripolar grid
      call test_uniform_mosaic('Folded-north symmetry mosaic') ! one-tile tripolar grid
      if(.not. wide_halo) then
         call test_uniform_mosaic('Folded-south symmetry mosaic') ! one-tile tripolar grid
         call test_uniform_mosaic('Folded-west symmetry mosaic') ! one-tile tripolar grid
         call test_uniform_mosaic('Folded-east symmetry mosaic') ! one-tile tripolar grid
         call test_uniform_mosaic('Four-Tile')
      endif
      call test_uniform_mosaic('Cubic-Grid') ! 6 tiles.
      call test_nonuniform_mosaic('Five-Tile')

      call test_halo_update( 'Simple' ) !includes global field, global sum tests
      call test_halo_update( 'Cyclic' )
      call test_halo_update( 'Folded-north' ) !includes vector field test
!      call test_halo_update( 'Masked' ) !includes vector field test
      call test_halo_update( 'Folded xy_halo' ) !
      if(.not. wide_halo) then
         call test_halo_update( 'Simple symmetry' ) !includes global field, global sum tests
         call test_halo_update( 'Cyclic symmetry' )
      endif
      call test_halo_update( 'Folded-north symmetry' ) !includes vector field test
      if(.not. wide_halo) then
         call test_halo_update( 'Folded-south symmetry' ) !includes vector field test
         call test_halo_update( 'Folded-west symmetry' ) !includes vector field test
         call test_halo_update( 'Folded-east symmetry' ) !includes vector field test
      endif

      !--- z1l: The following will not work due to symmetry and domain%x is cyclic.
      !--- Will solve this problem in the future if needed.
      ! call test_halo_update( 'Masked symmetry' ) !includes vector field test

      call test_global_field( 'Non-symmetry' )
      call test_global_field( 'Symmetry center' )
      call test_global_field( 'Symmetry corner' )
      call test_global_field( 'Symmetry east' )
      call test_global_field( 'Symmetry north' )

      if(.not. wide_halo) then
         call test_global_reduce( 'Simple')
         call test_global_reduce( 'Simple symmetry center')
         call test_global_reduce( 'Simple symmetry corner')
         call test_global_reduce( 'Simple symmetry east')
         call test_global_reduce( 'Simple symmetry north')
         call test_global_reduce( 'Cyclic symmetry center')
         call test_global_reduce( 'Cyclic symmetry corner')
         call test_global_reduce( 'Cyclic symmetry east')
         call test_global_reduce( 'Cyclic symmetry north')
      endif

      call test_redistribute( 'Complete pelist' )
!      call test_redistribute( 'Overlap  pelist' )
!      call test_redistribute( 'Disjoint pelist' )
      if(.not. wide_halo) then
         call test_define_mosaic_pelist('One tile', 1)
         call test_define_mosaic_pelist('Two uniform tile', 2)
         call test_define_mosaic_pelist('Two nonuniform tile', 2)
         call test_define_mosaic_pelist('Ten tile', 10)
         call test_define_mosaic_pelist('Ten tile with nonuniform cost', 10)
      endif
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finish test_interface  <-------------------'
  endif

  if( check_parallel) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_check_parallel <-------------------'
     call test_parallel( )
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finish test_check_parallel <-------------------'
  endif

!!$!Balaji adding openMP tests
!!$  call test_openmp()
!!$! Alewxander.Pletzer get_neighbor tests
  if( test_get_nbr ) then
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Calling test_get_nbr <-------------------'
     call test_get_neighbor_1d
     call test_get_neighbor_non_cyclic
     call test_get_neighbor_cyclic
     call test_get_neighbor_folded_north
     call test_get_neighbor_mask
     call mpp_sync()
      if (mpp_pe() == mpp_root_pe())  print *, '--------------------> Finish test_get_nbr <-------------------'
  endif

  call mpp_domains_exit()
  call mpp_exit()

contains
  subroutine test_openmp()
#ifdef _OPENMP_TEST
    integer :: omp_get_num_thread, omp_get_max_threads, omp_get_thread_num
    real, allocatable :: a(:,:,:)
    type(domain2D) :: domain
    integer :: layout(2)
    integer :: i,j,k, jthr
    integer :: thrnum, maxthr
    integer(LONG_KIND) :: sum1, sum2

    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain )
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( a(isd:ied,jsd:jed,nz) )
    maxthr = omp_get_max_threads()
    write( outunit,'(a,4i6)' )'pe,js,je,maxthr=', pe, js, je, maxthr
    if( mod(je-js+1,maxthr).NE.0 ) &
         call mpp_error( FATAL, 'maxthr must divide domain (TEMPORARY).' )
    jthr = (je-js+1)/maxthr
!$OMP PARALLEL PRIVATE(i,j,k,thrnum)
    thrnum = omp_get_thread_num()
    write( outunit,'(a,4i6)' )'pe,thrnum,js,je=', &
         pe, thrnum, js+thrnum*jthr,js+(thrnum+1)*jthr-1
    write( outunit,'(a,3i6)' )'pe,thrnum,node=', pe, thrnum, mpp_node()
!!$OMP DO
    do k = 1,nz
!when omp DO is commented out, user must compute j loop limits
!with omp DO, let OMP figure it out
       do j = js+thrnum*jthr,js+(thrnum+1)*jthr-1
!       do j = js,je
          do i = is,ie
             a(i,j,k) = global(i,j,k)
          end do
       end do
    end do
!!$OMP END DO
!$OMP END PARALLEL
    sum1 = mpp_chksum( a(is:ie,js:je,:) )
    sum2 = mpp_chksum( global(is:ie,js:je,:) )
    if( sum1.EQ.sum2 )then
        call mpp_error( NOTE, 'OMP parallel test OK.' )
    else
        if( mpp_pe().EQ.mpp_root_pe() )write( errunit,'(a,2z18)' )'OMP checksums: ', sum1, sum2
        call mpp_error( FATAL, 'OMP parallel test failed.' )
    end if
#endif
    return
  end subroutine test_openmp

  subroutine test_redistribute( type )
!test redistribute between two domains
    character(len=*), intent(in) :: type
    type(domain2D) :: domainx, domainy
    type(DomainCommunicator2D), pointer, save :: dch =>NULL()
    real, allocatable, dimension(:,:,:)       :: gcheck, global
    real, allocatable, dimension(:,:,:), save :: x, y
    real, allocatable, dimension(:,:,:), save :: x2, y2
    real, allocatable, dimension(:,:,:), save :: x3, y3
    real, allocatable, dimension(:,:,:), save :: x4, y4
    real, allocatable, dimension(:,:,:), save :: x5, y5
    real, allocatable, dimension(:,:,:), save :: x6, y6
    integer, allocatable :: pelist(:)
    integer :: pemax
    integer :: is, ie, js, je, isd, ied, jsd, jed

    pemax = npes/2              !the partial pelist will run from 0...pemax
    !--- nullify domain list otherwise it retains memory between calls.
    call mpp_nullify_domain_list(domainx)
    call mpp_nullify_domain_list(domainy)

    allocate( gcheck(nx,ny,nz), global(nx,ny,nz) )
    !fill in global array: with k.iiijjj
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             global(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

!select pelists
    select case(type)
    case( 'Complete pelist' )
!both pelists run from 0...npes-1
        if(nx < npes) then
           call mpp_error(NOTE, &
              "test_mpp_domains(test_redistribute): nx is less than npes, no test will be done for complete pelist")
           return
        endif
        allocate( pelist(0:npes-1) )
        pelist = (/ (i,i=0,npes-1) /)
        call mpp_declare_pelist( pelist )
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)
        call mpp_declare_pelist( pelist )
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1
        if( pemax+1.GE.npes )return
        allocate( pelist(0:pemax) )
        pelist = (/ (i,i=0,pemax) /)

        call mpp_declare_pelist( pelist )
        ! z1l: the follwing will cause deadlock will happen
        ! for npes = 6, x- mpp_global_field will call mpp_sync
        call mpp_declare_pelist( (/ (i,i=pemax+1,npes-1) /))
    case default
        call mpp_error( FATAL, 'TEST_REDISTRIBUTE: no such test: '//type )
    end select

!set up x and y arrays
    select case(type)
    case( 'Complete pelist' )
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        call mpp_define_domains( (/1,nx,1,ny/), (/npes,1/), domainy, name=type )
        call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
        allocate( y(isd:ied,jsd:jed,nz) )
        allocate( y2(isd:ied,jsd:jed,nz) )
        allocate( y3(isd:ied,jsd:jed,nz) )
        allocate( y4(isd:ied,jsd:jed,nz) )
        allocate( y5(isd:ied,jsd:jed,nz) )
        allocate( y6(isd:ied,jsd:jed,nz) )
        y = 0.
        y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
    case( 'Overlap  pelist' )
!one pelist from 0...pemax, other from 0...npes-1
!set up x array
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
        call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
        call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
        allocate( x(isd:ied,jsd:jed,nz) )
        allocate( x2(isd:ied,jsd:jed,nz) )
        allocate( x3(isd:ied,jsd:jed,nz) )
        allocate( x4(isd:ied,jsd:jed,nz) )
        allocate( x5(isd:ied,jsd:jed,nz) )
        allocate( x6(isd:ied,jsd:jed,nz) )
        x = 0.
        x(is:ie,js:je,:) = global(is:ie,js:je,:)
        x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        end if
    case( 'Disjoint pelist' )
!one pelist from 0...pemax, other from pemax+1...npes-1

!set up y array
        if( ANY(pelist.EQ.pe) )then
            call mpp_set_current_pelist(pelist)
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainy, name=type )
            call mpp_get_compute_domain( domainy, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainy, isd, ied, jsd, jed )
            allocate( y(isd:ied,jsd:jed,nz) )
            allocate( y2(isd:ied,jsd:jed,nz) )
            allocate( y3(isd:ied,jsd:jed,nz) )
            allocate( y4(isd:ied,jsd:jed,nz) )
            allocate( y5(isd:ied,jsd:jed,nz) )
            allocate( y6(isd:ied,jsd:jed,nz) )
            y = 0.
            y2 = 0.;y3 = 0.;y4 = 0.;y5 = 0.;y6 = 0.
        else
!set up x array
            call mpp_set_current_pelist( (/ (i,i=pemax+1,npes-1) /) )
            call mpp_define_layout( (/1,nx,1,ny/), mpp_npes(), layout )
            call mpp_define_domains( (/1,nx,1,ny/), layout, domainx, name=type )
            call mpp_get_compute_domain( domainx, is,  ie,  js,  je  )
            call mpp_get_data_domain   ( domainx, isd, ied, jsd, jed )
            allocate( x(isd:ied,jsd:jed,nz) )
            allocate( x2(isd:ied,jsd:jed,nz) )
            allocate( x3(isd:ied,jsd:jed,nz) )
            allocate( x4(isd:ied,jsd:jed,nz) )
            allocate( x5(isd:ied,jsd:jed,nz) )
            allocate( x6(isd:ied,jsd:jed,nz) )
            x = 0.
            x(is:ie,js:je,:) = global(is:ie,js:je,:)
            x2 = x;  x3 = x; x4 = x; x5 = x; x6 = x
         end if
    end select

!go global and redistribute
    call mpp_set_current_pelist()
    call mpp_broadcast_domain(domainx)
    call mpp_broadcast_domain(domainy)

    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_redistribute( domainx, x, domainy, y )
    call mpp_clock_end  (id)

!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if

    call mpp_set_current_pelist()

    call mpp_clock_begin(id)
    if(ALLOCATED(y))y=0.
    call mpp_redistribute( domainx, x,  domainy, y,  complete=.false. )
    call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
    call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
    call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
    call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
    call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
    call mpp_clock_end  (id)

!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if

    call mpp_set_current_pelist()

    if(type == 'Complete pelist')then
      write(outunit,*) 'Use domain communicator handle'
      call mpp_clock_begin(id)
      if(ALLOCATED(y))then
         y=0.; y2=0.; y3=0.; y4=0.; y5=0.; y6=0.
      endif
      call mpp_redistribute( domainx, x, domainy, y, complete=.false. )
      call mpp_redistribute( domainx, x2, domainy, y2, complete=.false. )
      call mpp_redistribute( domainx, x3, domainy, y3, complete=.false. )
      call mpp_redistribute( domainx, x4, domainy, y4, complete=.false. )
      call mpp_redistribute( domainx, x5, domainy, y5, complete=.false. )
      call mpp_redistribute( domainx, x6, domainy, y6, complete=.true., dc_handle=dch )
      call mpp_clock_end  (id)

!check answers on pelist
    if( ANY(pelist.EQ.pe) )then
        call mpp_set_current_pelist(pelist)
        call mpp_global_field( domainy, y, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y2, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y3, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y4, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y5, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
        call mpp_global_field( domainy, y6, gcheck )
        call compare_checksums( global(1:nx,1:ny,:), gcheck, type )
    end if
    endif
    dch =>NULL()

    call mpp_set_current_pelist()

    deallocate(gcheck, global)
    if(ALLOCATED(pelist)) deallocate(pelist)

    if(ALLOCATED(x))then
      call mpp_redistribute( domainx, x, domainy, y, free=.true.,list_size=6 )
      deallocate(x,x2,x3,x4,x5,x6)
    endif
    if(ALLOCATED(y))deallocate(y,y2,y3,y4,y5,y6)
  end subroutine test_redistribute

  subroutine cubic_grid_redistribute

     integer              :: npes, npes_per_ensemble, npes_per_tile
     integer              :: ensemble_id, tile_id, ensemble_tile_id
     integer              :: i, j, p, n, ntiles, my_root_pe
     integer              :: isc_ens, iec_ens, jsc_ens, jec_ens
     integer              :: isd_ens, ied_ens, jsd_ens, jed_ens
     integer              :: isc, iec, jsc, jec
     integer              :: isd, ied, jsd, jed
     integer, allocatable :: my_ensemble_pelist(:), pe_start(:), pe_end(:)
     integer, allocatable :: global_indices(:,:), layout2D(:,:)
     real,    allocatable :: x(:,:,:,:), x_ens(:,:,:), y(:,:,:)
     integer              :: layout(2)
     type(domain2D)       :: domain
     type(domain2D), allocatable :: domain_ensemble(:)
     character(len=128)   :: mesg

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

     call define_cubic_mosaic("cubic_grid", domain_ensemble(ensemble_id), (/nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic,nx_cubic/), &
                              (/ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic,ny_cubic/), &
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
        call compare_checksums( x(isc:iec,jsc:jec,:,n), y(isc:iec,jsc:jec,:), trim(mesg) )
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
        call compare_checksums( x_ens(isc_ens:iec_ens,jsc_ens:jec_ens,:), y(isc_ens:iec_ens,jsc_ens:jec_ens,:), trim(mesg) )
     enddo

     deallocate(x, y, x_ens)
     call mpp_deallocate_domain(domain)
     do n = 1, ensemble_size
        call mpp_deallocate_domain(domain_ensemble(n))
     enddo
     deallocate(domain_ensemble)

  end subroutine cubic_grid_redistribute


  subroutine test_uniform_mosaic( type )
    character(len=*), intent(in) :: type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile, ntile_per_pe, update_flags
    integer        :: i, j, k, l, n, shift, tw, te, ts, tn, tsw, tnw, tse, tne
    integer        :: ism, iem, jsm, jem, wh, eh, sh, nh
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    real           :: gsum, lsum

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:)     :: global2D
    real,    allocatable, dimension(:,:,:)   :: local1, local2
    real,    allocatable, dimension(:,:,:,:) :: x, y, x1, x2, x3, x4, y1, y2, y3, y4
    real,    allocatable, dimension(:,:,:,:) :: global1, global2, gcheck
    real,    allocatable, dimension(:,:,:,:) :: global1_all, global2_all, global_all
    character(len=256) :: type2, type3
    logical            :: folded_north, folded_north_sym, folded_north_nonsym
    logical            :: folded_south_sym, folded_west_sym, folded_east_sym
    logical            :: cubic_grid, single_tile, four_tile
    integer            :: whalo_save, ehalo_save, nhalo_save, shalo_save
    integer            :: nx_save, ny_save
    logical            :: same_layout = .false.


    nx_save = nx
    ny_save = ny
    if(type == 'Cubic-Grid' .and. nx_cubic >0) then
       nx = nx_cubic
       ny = ny_cubic
    endif

    if(wide_halo_x > 0) then
       whalo_save = whalo
       ehalo_save = ehalo
       shalo_save = shalo
       nhalo_save = nhalo
       if(type == 'Single-Tile' .OR. type == 'Folded-north mosaic' .OR. type == 'Cubic-Grid') then
          whalo = wide_halo_x
          ehalo = wide_halo_x
          shalo = wide_halo_y
          nhalo = wide_halo_y
       endif
    endif

    folded_north_nonsym = .false.
    folded_north_sym    = .false.
    folded_north        = .false.
    folded_south_sym    = .false.
    folded_west_sym     = .false.
    folded_east_sym     = .false.
    cubic_grid        = .false.
    single_tile        = .false.
    four_tile          = .false.
    !--- check the type
    select case(type)
    case ( 'Single-Tile' )   !--- single with cyclic along x- and y-direction
       single_tile = .true.
       ntiles = 1
       num_contact = 2
    case ( 'Folded-north mosaic' )
       ntiles = 1
       num_contact = 2
       folded_north_nonsym = .true.
    case ( 'Folded-north symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_north_sym = .true.
    case ( 'Folded-south symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_south_sym = .true.
    case ( 'Folded-west symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_west_sym = .true.
    case ( 'Folded-east symmetry mosaic' )
       ntiles = 1
       num_contact = 2
       folded_east_sym = .true.
    case ( 'Four-Tile' ) !--- cyclic along both x- and y-direction.
       ntiles = 4
       num_contact = 8
       four_tile = .true.
    case ( 'Cubic-Grid' )
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( nx .NE. ny) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: for Cubic_grid mosaic, nx should equal ny, '//&
                   'No test is done for Cubic-Grid mosaic. ' )
          if(wide_halo_x > 0) then
             whalo = whalo_save
             ehalo = ehalo_save
             shalo = shalo_save
             nhalo = nhalo_save
             nx    = nx_save
             ny    = ny_save
          endif
          return
       end if
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS: no such test: '//type)
    end select

    folded_north = folded_north_nonsym .OR. folded_north_sym

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    if( mod(npes, ntiles) == 0 ) then
       npes_per_tile = npes/ntiles
       write(outunit,*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                       '", each tile will be distributed over ', npes_per_tile, ' processors.'
       ntile_per_pe = 1
       allocate(tile(ntile_per_pe))
       tile = pe/npes_per_tile+1
       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do
    else if ( mod(ntiles, npes) == 0 ) then
       ntile_per_pe = ntiles/npes
       write(outunit,*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                        '", there will be ', ntile_per_pe, ' tiles on each processor.'
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
       call mpp_error(NOTE,'TEST_MPP_DOMAINS: npes should be multiple of ntiles or ' // &
            'ntiles should be multiple of npes. No test is done for '//trim(type) )
       nx = nx_save
       ny = ny_save
       if(wide_halo_x > 0) then
          whalo = whalo_save
          ehalo = ehalo_save
          shalo = shalo_save
          nhalo = nhalo_save
       endif
       return
    end if

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do
    same_layout = .false.
    if(layout(1) == layout(2)) same_layout = .true.

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    call mpp_memuse_begin()
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
            name = type, symmetry = .false. )
    else if(folded_north) then
       !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (NORTH) and tile 1 (NORTH)  --- folded-north-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = ny;  jend1(2) = ny
       istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = ny;  jend2(2) = ny
       if(folded_north_nonsym) then
          call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                                 istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                                 pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                                 name = type, symmetry = .false.  )
       else
          call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                                 istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                                 pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                                 name = type, symmetry = .true.  )
       endif
    else if(folded_south_sym) then
       !--- Contact line 1, between tile 1 (EAST) and tile 1 (WEST)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = nx; iend1(1) = nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;  iend2(1) = 1;  jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (SOUTH) and tile 1 (SOUTH)  --- folded-south-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = nx/2;   jstart1(2) = 1;  jend1(2) = 1
       istart2(2) = nx; iend2(2) = nx/2+1; jstart2(2) = 1;  jend2(2) = 1
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                              name = type, symmetry = .true.  )
    else if(folded_west_sym) then
       !--- Contact line 1, between tile 1 (NORTH) and tile 1 (SOUTH)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = 1; iend1(1) = nx; jstart1(1) = ny;  jend1(1) = ny
       istart2(1) = 1; iend2(1) = nx; jstart2(1) = 1;   jend2(1) = 1
       !--- Contact line 2, between tile 1 (WEST) and tile 1 (WEST)  --- folded-west-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = 1;  iend1(2) = 1; jstart1(2) = 1;  jend1(2) = ny/2
       istart2(2) = 1;  iend2(2) = 1; jstart2(2) = ny; jend2(2) = ny/2+1
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                              name = type, symmetry = .true.  )
    else if(folded_east_sym) then
       !--- Contact line 1, between tile 1 (NORTH) and tile 1 (SOUTH)  --- cyclic
       tile1(1) = 1; tile2(1) = 1
       istart1(1) = 1; iend1(1) = nx; jstart1(1) = ny;  jend1(1) = ny
       istart2(1) = 1; iend2(1) = nx; jstart2(1) = 1;   jend2(1) = 1
       !--- Contact line 2, between tile 1 (EAST) and tile 1 (EAST)  --- folded-west-edge
       tile1(2) = 1; tile2(2) = 1
       istart1(2) = nx;  iend1(2) = nx; jstart1(2) = 1;  jend1(2) = ny/2
       istart2(2) = nx;  iend2(2) = nx; jstart2(2) = ny; jend2(2) = ny/2+1
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
                              istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
                              pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                              name = type, symmetry = .true.  )
    else if( four_tile ) then
       call define_fourtile_mosaic(type, domain, (/nx,nx,nx,nx/), (/ny,ny,ny,ny/), global_indices, &
                                   layout2D, pe_start, pe_end, symmetry = .false. )
    else if( cubic_grid ) then
       call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                                global_indices, layout2D, pe_start, pe_end )
    endif
    call mpp_memuse_end(trim(type)//" mpp_define_mosaic", outunit )

    !--- setup data
    allocate(global2(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz, ntile_per_pe) )
    allocate(global_all(1:nx,1:ny,nz, ntiles) )
    global2 = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                global_all(i,j,k,l) = l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global2(1:nx,1:ny,:,n) = global_all(:,:,:,tile(n))
    end do

    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( gcheck(nx, ny, nz, ntile_per_pe) )
    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x1(ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x2(ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x3(ism:iem,jsm:jem,nz, ntile_per_pe) )
    allocate( x4(ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0.
    x(isc:iec,jsc:jec,:,:) = global2(isc:iec,jsc:jec,:,:)
    x1 = x; x2 = x; x3 = x; x4 = x;

    !--- test mpp_global_sum
    gsum = 0
    allocate(global2D(nx,ny))
    do n = 1, ntiles
       do j = 1, ny
          do i = 1, nx
             global2D(i,j) = sum(global_all(i,j,:,n))
          end do
       end do
       gsum = gsum + sum(global2D)
    end do

    do n = 1, ntile_per_pe
       lsum = mpp_global_sum( domain, x(:,:,:,n), tile_count=n )
    end do
    if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', type//' Fast sum=', lsum, gsum

    !test exact mpp_global_sum
    do n = 1, ntile_per_pe
       lsum = mpp_global_sum( domain, x(:,:,:,n), BITWISE_EXACT_SUM, tile_count=n)
    end do
    call compare_data_scalar(lsum, gsum, FATAL, type//' mpp_global_exact_sum')

    !--- test mpp_global_field
    gcheck = 0.
    id = mpp_clock_id( type//' global field ', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_global_field( domain, x(:,:,:,n), gcheck(:,:,:,n), tile_count=n)
    end do
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    do n = 1, ntile_per_pe
       call compare_checksums( global2(1:nx,1:ny,:,n), gcheck(:,:,:,n), type//' mpp_global_field ' )
    end do

    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    do n = 1, ntile_per_pe
       !--- fill up the value at halo points.
       if(single_tile) then
          call fill_regular_mosaic_halo(global2(:,:,:,n), global_all, 1, 1, 1, 1, 1, 1, 1, 1)
       else if(folded_north) then
          call fill_folded_north_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(folded_south_sym) then
          call fill_folded_south_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(folded_west_sym) then
          call fill_folded_west_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(folded_east_sym) then
          call fill_folded_east_halo(global2(:,:,:,n), 0, 0, 0, 0, 1)
       else if(four_tile) then
          select case ( tile(n) )
          case (1)
             tw = 2; ts = 3; tsw = 4
          case (2)
             tw = 1; ts = 4; tsw = 3
          case (3)
             tw = 4; ts = 1; tsw = 2
          case (4)
             tw = 3; ts = 2; tsw = 1
          end select
          te = tw; tn = ts; tse = tsw; tnw = tsw; tne = tsw
          call fill_regular_mosaic_halo(global2(:,:,:,n), global_all, te, tse, ts, tsw, tw, tnw, tn, tne )
       else if(cubic_grid) then
          call fill_cubic_grid_halo(global2(:,:,:,n), global_all, global_all, tile(n), 0, 0, 1, 1 )
       endif

       !full update
       call mpp_clock_begin(id)
       if(ntile_per_pe == 1) then
          call mpp_update_domains( x(:,:,:,n), domain )
       else
          call mpp_update_domains( x(:,:,:,n), domain, tile_count = n )
       end if
       call mpp_clock_end  (id)
    end do
    type2 = type
    do n = 1, ntile_per_pe
       if(ntile_per_pe>1)   write(type2, *)type, " at tile_count = ",n
       call compare_checksums( x(ism:ism+ied-isd,jsm:jsm+jed-jsd,:,n), global2(isd:ied,jsd:jed,:,n), trim(type2) )
    end do

    !partial update only be done when there is at most one tile on each pe
    if(ntile_per_pe == 1 ) then
       id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
       call mpp_clock_begin(id)
       call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
       call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
       call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
       call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
       call mpp_clock_end  (id)
       call compare_checksums( x1(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x1' )
       call compare_checksums( x2(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x2' )
       call compare_checksums( x3(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x3' )
       call compare_checksums( x4(isc:ied,jsc:jed,:,1), global2(isc:ied,jsc:jed,:,1), type//' partial x4' )

       !arbitrary halo update. not for tripolar grid
       if(wide_halo_x == 0) then
          if(single_tile .or. four_tile .or. (cubic_grid .and. same_layout) .or. folded_north ) then
             allocate(local2(isd:ied,jsd:jed,nz) )
             do wh = 1, whalo
                do eh = 1, ehalo
                   if(wh .NE. eh) cycle
                   do sh = 1, shalo
                      do nh = 1, nhalo
                         if(sh .NE. nh) cycle
                         local2(isd:ied,jsd:jed,:) = global2(isd:ied,jsd:jed,:,1)
                         x = 0.
                         x(isc:iec,jsc:jec,:,1) = local2(isc:iec,jsc:jec,:)
                         call fill_halo_zero(local2, wh, eh, sh, nh, 0, 0, isc, iec, jsc, jec, isd, ied, jsd, jed)

                         write(type2,'(a,a,i2,a,i2,a,i2,a,i2)') trim(type), ' with whalo = ', wh, &
                              ', ehalo = ',eh, ', shalo = ', sh, ', nhalo = ', nh
                         call mpp_update_domains( x, domain, whalo=wh, ehalo=eh, shalo=sh, nhalo=nh, name = type2  )
                         call compare_checksums( x(isd:ied,jsd:jed,:,1), local2, trim(type2) )
                      end do
                   end do
                end do
             end do
             deallocate(local2)
          end if
       endif
    end if

    deallocate(global2, global_all, x, x1, x2, x3, x4)
    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    !--- setup data
    shift = 0
    if(single_tile .or. four_tile .or. folded_north_nonsym) then
       shift = 0
    else
       shift = 1
    endif

    allocate(global1(1-whalo:nx+shift+ehalo,1-shalo:ny+shift+nhalo,nz,ntile_per_pe) )
    allocate(global2(1-whalo:nx+shift+ehalo,1-shalo:ny+shift+nhalo,nz,ntile_per_pe) )
    allocate(global1_all(nx+shift,ny+shift,nz, ntiles),  global2_all(nx+shift,ny+shift,nz, ntiles))
    global1 = 0; global2 = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+shift
             do i = 1, nx+shift
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !-----------------------------------------------------------------------
    !--- make sure consistency on the boundary for cubic grid
    !--- east boundary will take the value of neighbor tile ( west/south),
    !--- north boundary will take the value of neighbor tile ( south/west).
    !--- for the point on the corner, the 12 corner take the following value
    !--- corner between 1, 2, 3 takes the value at 3,
    !--- corner between 1, 3, 5 takes the value at 3
    !-----------------------------------------------------------------------
    if( cubic_grid ) then
       do l = 1, ntiles
          if(mod(l,2) == 0) then ! tile 2, 4, 6
             te = l + 2
             tn = l + 1
             if(te>6) te = te - 6
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,1:ny+1,:,l) = global2_all(nx+shift:1:-1,1,:,te)  ! east
             global2_all(nx+shift,1:ny+1,:,l) = global1_all(nx+shift:1:-1,1,:,te)  ! east
             global1_all(1:nx,ny+shift,:,l)    = global1_all(1:nx,1,:,tn) ! north
             global2_all(1:nx,ny+shift,:,l)    = global2_all(1:nx,1,:,tn) ! north
          else                   ! tile 1, 3, 5
             te = l + 1
             tn = l + 2
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,:,:,l)    = global1_all(1,:,:,te)  ! east
             global2_all(nx+shift,:,:,l)    = global2_all(1,:,:,te)  ! east
             global1_all(1:nx+1,ny+shift,:,l) = global2_all(1,ny+shift:1:-1,:,tn) ! north
             global2_all(1:nx+1,ny+shift,:,l) = global1_all(1,ny+shift:1:-1,:,tn) ! north
          end if
       end do
       ! set the corner value to 0
       global1_all(1,ny+1,:,:) = 0; global1_all(nx+1,1,:,:) = 0; global1_all(1,1,:,:) = 0; global1_all(nx+1,ny+1,:,:) = 0
       global2_all(1,ny+1,:,:) = 0; global2_all(nx+1,1,:,:) = 0; global2_all(1,1,:,:) = 0; global2_all(nx+1,ny+1,:,:) = 0
    end if

    do n = 1, ntile_per_pe
       global1(1:nx+shift,1:ny+shift,:,n) = global1_all(:,:,:,tile(n))
       global2(1:nx+shift,1:ny+shift,:,n) = global2_all(:,:,:,tile(n))
    end do

    if(folded_north) then
       call fill_folded_north_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)
       call fill_folded_north_halo(global2(:,:,:,1), 1, 1, shift, shift, -1)
    else if(folded_south_sym) then
       call fill_folded_south_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)
       call fill_folded_south_halo(global2(:,:,:,1), 1, 1, shift, shift, -1)
    else if(folded_west_sym) then
       call fill_folded_west_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)
       call fill_folded_west_halo(global2(:,:,:,1), 1, 1, shift, shift, -1)
    else if(folded_east_sym) then
       call fill_folded_east_halo(global1(:,:,:,1), 1, 1, shift, shift, -1)
       call fill_folded_east_halo(global2(:,:,:,1), 1, 1, shift, shift, -1)
    endif

    allocate( x (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x1(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x2(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x3(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( x4(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y1(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y2(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y3(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y4(ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec+shift,:,:) = global1(isc:iec+shift,jsc:jec+shift,:,:)
    y (isc:iec+shift,jsc:jec+shift,:,:) = global2(isc:iec+shift,jsc:jec+shift,:,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points.
    !-----------------------------------------------------------------------
    if(cubic_grid) then
       type2 = type//' paired-scalar BGRID_NE'
       update_flags = SCALAR_PAIR
    else
       type2 = type//' vector BGRID_NE'
       update_flags = XUPDATE + YUPDATE
    endif

    id = mpp_clock_id( trim(type2), flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    type3 = type2

    do n = 1, ntile_per_pe
       if(single_tile) then
          call fill_regular_mosaic_halo(global1(:,:,:,n), global1_all, 1, 1, 1, 1, 1, 1, 1, 1)
          call fill_regular_mosaic_halo(global2(:,:,:,n), global2_all, 1, 1, 1, 1, 1, 1, 1, 1)
       else if(folded_north) then
          !redundant points must be equal and opposite for tripolar grid
          global1(nx/2+shift,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global1(nx+shift  ,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global1(nx/2+1+shift:nx-1+shift,   ny+shift,:,:) = -global1(nx/2-1+shift:1+shift:-1, ny+shift,:,:)
          global1(1-whalo:shift,             ny+shift,:,:) = -global1(nx-whalo+1:nx+shift,     ny+shift,:,:)
          global1(nx+1+shift:nx+ehalo+shift, ny+shift,:,:) = -global1(1+shift:ehalo+shift,     ny+shift,:,:)
          global2(nx/2+shift,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global2(nx+shift  ,                ny+shift,:,:) = 0.  !pole points must have 0 velocity
          global2(nx/2+1+shift:nx-1+shift,   ny+shift,:,:) = -global2(nx/2-1+shift:1+shift:-1, ny+shift,:,:)
          global2(1-whalo:shift,             ny+shift,:,:) = -global2(nx-whalo+1:nx+shift,     ny+shift,:,:)
          global2(nx+1+shift:nx+ehalo+shift, ny+shift,:,:) = -global2(1+shift:ehalo+shift,     ny+shift,:,:)
          !--- the following will fix the +0/-0 problem on altix
          if(nhalo >0) then
             global1(shift,ny+shift,:,:) = 0.  !pole points must have 0 velocity
             global2(shift,ny+shift,:,:) = 0.  !pole points must have 0 velocity
          end if
       else if(folded_south_sym) then
          global1(nx/2+shift,                1,:,:) = 0.  !pole points must have 0 velocity
          global1(nx+shift  ,                1,:,:) = 0.  !pole points must have 0 velocity
          global1(nx/2+1+shift:nx-1+shift,   1,:,:) = -global1(nx/2-1+shift:1+shift:-1, 1,:,:)
          global1(1-whalo:shift,             1,:,:) = -global1(nx-whalo+1:nx+shift,     1,:,:)
          global1(nx+1+shift:nx+ehalo+shift, 1,:,:) = -global1(1+shift:ehalo+shift,     1,:,:)
          global2(nx/2+shift,                1,:,:) = 0.  !pole points must have 0 velocity
          global2(nx+shift  ,                1,:,:) = 0.  !pole points must have 0 velocity
          global2(nx/2+1+shift:nx-1+shift,   1,:,:) = -global2(nx/2-1+shift:1+shift:-1, 1,:,:)
          global2(1-whalo:shift,             1,:,:) = -global2(nx-whalo+1:nx+shift,     1,:,:)
          global2(nx+1+shift:nx+ehalo+shift, 1,:,:) = -global2(1+shift:ehalo+shift,     1,:,:)
          !--- the following will fix the +0/-0 problem on altix
          if(shalo >0) then
             global1(shift,1,:,:) = 0.  !pole points must have 0 velocity
             global2(shift,1,:,:) = 0.  !pole points must have 0 velocity
          endif
       else if(folded_west_sym) then
          global1(1, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global1(1, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global1(1, ny/2+1+shift:ny-1+shift,   :,:) = -global1(1, ny/2-1+shift:1+shift:-1, :,:)
          global1(1, 1-shalo:shift,             :,:) = -global1(1, ny-shalo+1:ny+shift,     :,:)
          global1(1, ny+1+shift:ny+nhalo+shift, :,:) = -global1(1, 1+shift:nhalo+shift,     :,:)
          global2(1, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global2(1, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global2(1, ny/2+1+shift:ny-1+shift,   :,:) = -global2(1, ny/2-1+shift:1+shift:-1, :,:)
          global2(1, 1-shalo:shift,             :,:) = -global2(1, ny-shalo+1:ny+shift,     :,:)
          global2(1, ny+1+shift:ny+nhalo+shift, :,:) = -global2(1, 1+shift:nhalo+shift,     :,:)
          !--- the following will fix the +0/-0 problem on altix
          if(whalo>0) then
             global1(1, shift, :, :) = 0.  !pole points must have 0 velocity
             global2(1, shift, :, :) = 0.  !pole points must have 0 velocity
          endif
       else if(folded_east_sym) then
          global1(nx+shift, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global1(nx+shift, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global1(nx+shift, ny/2+1+shift:ny-1+shift,   :,:) = -global1(nx+shift, ny/2-1+shift:1+shift:-1, :,:)
          global1(nx+shift, 1-shalo:shift,             :,:) = -global1(nx+shift, ny-shalo+1:ny+shift,     :,:)
          global1(nx+shift, ny+1+shift:ny+nhalo+shift, :,:) = -global1(nx+shift, 1+shift:nhalo+shift,     :,:)
          global2(nx+shift, ny/2+shift, :,:) = 0. !pole points must have 0 velocity
          global2(nx+shift, ny+shift,   :,:) = 0. !pole points must have 0 velocity
          global2(nx+shift, ny/2+1+shift:ny-1+shift,   :,:) = -global2(nx+shift, ny/2-1+shift:1+shift:-1, :,:)
          global2(nx+shift, 1-shalo:shift,             :,:) = -global2(nx+shift, ny-shalo+1:ny+shift,     :,:)
          global2(nx+shift, ny+1+shift:ny+nhalo+shift, :,:) = -global2(nx+shift, 1+shift:nhalo+shift,     :,:)
          !--- the following will fix the +0/-0 problem on altix
          if(ehalo >0) then
             global1(nx+shift, shift, :,:) = 0.  !pole points must have 0 velocity
             global2(nx+shift, shift, :,:) = 0.  !pole points must have 0 velocity
          end if
       else if(four_tile) then
          select case ( tile(n) )
          case (1)
             tw = 2; ts = 3; tsw = 4
          case (2)
             tw = 1; ts = 4; tsw = 3
          case (3)
             tw = 4; ts = 1; tsw = 2
          case (4)
             tw = 3; ts = 2; tsw = 1
          end select
          te = tw; tn = ts; tse = tsw; tnw = tsw; tne = tsw
          call fill_regular_mosaic_halo(global1(:,:,:,n), global1_all, te, tse, ts, tsw, tw, tnw, tn, tne )
          call fill_regular_mosaic_halo(global2(:,:,:,n), global2_all, te, tse, ts, tsw, tw, tnw, tn, tne )
       else if(cubic_grid) then
          call fill_cubic_grid_halo(global1(:,:,:,n), global1_all, global2_all, tile(n), 1, 1, 1, 1 )
          call fill_cubic_grid_halo(global2(:,:,:,n), global2_all, global1_all, tile(n), 1, 1, 1, 1 )
       endif

       if(ntile_per_pe > 1) write(type3, *)trim(type2), " at tile_count = ",n
       call mpp_clock_begin(id)
       if(ntile_per_pe == 1) then
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, flags=update_flags, gridtype=BGRID_NE, name=type3)
       else
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, flags=update_flags, gridtype=BGRID_NE, &
               name=type3, tile_count = n)
       end if
       call mpp_clock_end  (id)
    end do

    do n = 1, ntile_per_pe
       if(ntile_per_pe > 1) write(type3, *)trim(type2), " at tile_count = ", n
       call compare_checksums( x (isd:ied+shift,jsd:jed+shift,:,n),  global1(isd:ied+shift,jsd:jed+shift,:,n), trim(type3)//' X' )
       call compare_checksums( y (isd:ied+shift,jsd:jed+shift,:,n),  global2(isd:ied+shift,jsd:jed+shift,:,n), trim(type3)//' Y' )
    end do

    if(ntile_per_pe == 1) then
       call mpp_clock_begin(id)
       call mpp_update_domains( x1, y1, domain, flags=update_flags, gridtype=BGRID_NE, complete=.false., name=type2)
       call mpp_update_domains( x2, y2, domain, flags=update_flags, gridtype=BGRID_NE, complete=.false.,  name=type2)
       call mpp_update_domains( x3, y3, domain, flags=update_flags, gridtype=BGRID_NE, complete=.false., name=type2)
       call mpp_update_domains( x4, y4, domain, flags=update_flags, gridtype=BGRID_NE, complete=.true.,  name=type2)
       call mpp_clock_end  (id)

       call compare_checksums( x1(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X1')
       call compare_checksums( x2(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X2')
       call compare_checksums( x3(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X3')
       call compare_checksums( x4(isd:ied+shift,jsd:jed+shift,:,1), global1(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' X4')
       call compare_checksums( y1(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y1')
       call compare_checksums( y2(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y2')
       call compare_checksums( y3(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y3')
       call compare_checksums( y4(isd:ied+shift,jsd:jed+shift,:,1), global2(isd:ied+shift,jsd:jed+shift,:,1), trim(type2)//' Y4')

       !--- arbitrary halo updates ---------------------------------------
       if(wide_halo_x == 0) then
          if(single_tile .or. four_tile .or. (cubic_grid .and. same_layout) .or. folded_north) then
             allocate(local1(isd:ied+shift,jsd:jed+shift,nz) )
             allocate(local2(isd:ied+shift,jsd:jed+shift,nz) )
             do wh = 1, whalo
                do eh = 1, ehalo
                   if(wh .NE. eh) cycle
                   do sh = 1, shalo
                      do nh = 1, nhalo
                         if(nh .NE. sh) cycle

                         local1(isd:ied+shift,jsd:jed+shift,:) = global1(isd:ied+shift,jsd:jed+shift,:,1)
                         local2(isd:ied+shift,jsd:jed+shift,:) = global2(isd:ied+shift,jsd:jed+shift,:,1)
                         x = 0.; y = 0.
                         x(isc:iec+shift,jsc:jec+shift,:,1) = global1(isc:iec+shift,jsc:jec+shift,:,1)
                         y(isc:iec+shift,jsc:jec+shift,:,1) = global2(isc:iec+shift,jsc:jec+shift,:,1)

                         call fill_halo_zero(local1, wh, eh, sh, nh, shift, shift, isc, iec, jsc, jec, isd, ied, jsd, jed)
                         call fill_halo_zero(local2, wh, eh, sh, nh, shift, shift, isc, iec, jsc, jec, isd, ied, jsd, jed)

                         write(type3,'(a,a,i2,a,i2,a,i2,a,i2)') trim(type2), ' with whalo = ', wh, &
                              ', ehalo = ',eh, ', shalo = ', sh, ', nhalo = ', nh
                         call mpp_update_domains( x,  y,  domain, flags=update_flags, gridtype=BGRID_NE, &
                              whalo=wh, ehalo=eh, shalo=sh, nhalo=nh, name=type3)
                         call compare_checksums( x(isd:ied+shift,jsd:jed+shift,:,1),  local1, trim(type3)//' X' )
                         call compare_checksums( y(isd:ied+shift,jsd:jed+shift,:,1),  local2, trim(type3)//' Y' )
                      end do
                   end do
                end do
             end do
             deallocate(local1, local2)
          end if
       endif
    end if
    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- setup data
    if(cubic_grid .or. folded_north .or. folded_south_sym .or. folded_west_sym .or. folded_east_sym ) then
       deallocate(global1_all, global2_all)
       allocate(global1_all(nx+shift,ny,nz, ntiles),  global2_all(nx,ny+shift,nz, ntiles))
       deallocate(global1, global2, x, y, x1, x2, x3, x4, y1, y2, y3, y4)
       allocate(global1(1-whalo:nx+shift+ehalo,1-shalo:ny  +nhalo,nz,ntile_per_pe) )
       allocate( x (ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( y (ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( x1(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( x2(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( x3(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( x4(ism:iem+shift,jsm:jem  ,nz,ntile_per_pe) )
       allocate( y1(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( y2(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( y3(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate( y4(ism:iem  ,jsm:jem+shift,nz,ntile_per_pe) )
       allocate(global2(1-whalo:nx  +ehalo,1-shalo:ny+shift+nhalo,nz,ntile_per_pe) )
       global1 = 0; global2 = 0
       do l = 1, ntiles
          do k = 1, nz
             do j = 1, ny
                do i = 1, nx+shift
                   global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                end do
             end do
             do j = 1, ny+shift
                do i = 1, nx
                   global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                end do
             end do
          end do
       end do
    endif
    if( folded_north .or. folded_south_sym .or. folded_west_sym .or. folded_east_sym ) then
       do n = 1, ntile_per_pe
          global1(1:nx+shift,1:ny  ,:,n) = global1_all(1:nx+shift,1:ny,  :,tile(n))
          global2(1:nx  ,1:ny+shift,:,n) = global2_all(1:nx  ,1:ny+shift,:,tile(n))
       end do
    endif

    if( cubic_grid ) then
       !-----------------------------------------------------------------------
       !--- make sure consistency on the boundary for cubic grid
       !--- east boundary will take the value of neighbor tile ( west/south),
       !--- north boundary will take the value of neighbor tile ( south/west).
       !-----------------------------------------------------------------------
       do l = 1, ntiles
          if(mod(l,2) == 0) then ! tile 2, 4, 6
             te = l + 2
             tn = l + 1
             if(te>6) te = te - 6
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,1:ny,:,l) = global2_all(nx:1:-1,1,:,te)  ! east
             global2_all(1:nx,ny+shift,:,l) = global2_all(1:nx,1,:,tn) ! north
          else                   ! tile 1, 3, 5
             te = l + 1
             tn = l + 2
             if(tn > 6) tn = tn - 6
             global1_all(nx+shift,:,:,l)    = global1_all(1,:,:,te)  ! east
             global2_all(1:nx,ny+shift,:,l) = global1_all(1,ny:1:-1,:,tn) ! north
          end if
       end do
       do n = 1, ntile_per_pe
          global1(1:nx+shift,1:ny  ,:,n) = global1_all(1:nx+shift,1:ny,  :,tile(n))
          global2(1:nx  ,1:ny+shift,:,n) = global2_all(1:nx  ,1:ny+shift,:,tile(n))
       end do
    else if( folded_north ) then
       call fill_folded_north_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_north_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    else if(folded_south_sym ) then
       call fill_folded_south_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_south_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    else if(folded_west_sym ) then
       call fill_folded_west_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_west_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    else if(folded_east_sym ) then
       call fill_folded_east_halo(global1(:,:,:,1), 1, 0, shift, 0, -1)
       call fill_folded_east_halo(global2(:,:,:,1), 0, 1, 0, shift, -1)
    endif
    x = 0.; y = 0.
    x (isc:iec+shift,jsc:jec  ,:,:) = global1(isc:iec+shift,jsc:jec  ,:,:)
    y (isc:iec  ,jsc:jec+shift,:,:) = global2(isc:iec  ,jsc:jec+shift,:,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points for cubic-grid.
    !   On the contact line, the following relation will be used to
    !   --- fill the value on contact line ( balance send and recv).
    !       2W --> 1E, 1S --> 6N, 3W --> 1N, 4S --> 2E
    !       4W --> 3E, 3S --> 2N, 1W --> 5N, 2S --> 6E
    !       6W --> 5E, 5S --> 4N, 5W --> 3N, 6S --> 4E
    !---------------------------------------------------------------------------
    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    type2 = type
    do n = 1, ntile_per_pe
       if( cubic_grid ) then
          call fill_cubic_grid_halo(global1(:,:,:,n), global1_all, global2_all, tile(n), 1, 0, 1, -1 )
          call fill_cubic_grid_halo(global2(:,:,:,n), global2_all, global1_all, tile(n), 0, 1, -1, 1 )
       else if( folded_north ) then
          !redundant points must be equal and opposite
          global2(nx/2+1:nx,     ny+shift,:,:) = -global2(nx/2:1:-1, ny+shift,:,:)
          global2(1-whalo:0,     ny+shift,:,:) = -global2(nx-whalo+1:nx, ny+shift,:,:)
          global2(nx+1:nx+ehalo, ny+shift,:,:) = -global2(1:ehalo,       ny+shift,:,:)
       else if( folded_south_sym ) then
          global2(nx/2+1:nx,     1,:,:) = -global2(nx/2:1:-1, 1,:,:)
          global2(1-whalo:0,     1,:,:) = -global2(nx-whalo+1:nx, 1, :,:)
          global2(nx+1:nx+ehalo, 1,:,:) = -global2(1:ehalo,       1, :,:)
       else if( folded_west_sym ) then
          global1(1, ny/2+1:ny,     :,:) = -global1(1, ny/2:1:-1,     :,:)
          global1(1, 1-shalo:0,     :,:) = -global1(1, ny-shalo+1:ny, :,:)
          global1(1, ny+1:ny+nhalo, :,:) = -global1(1, 1:nhalo,       :,:)
       else if( folded_east_sym ) then
          global1(nx+shift, ny/2+1:ny,     :,:) = -global1(nx+shift, ny/2:1:-1,     :,:)
          global1(nx+shift, 1-shalo:0,     :,:) = -global1(nx+shift, ny-shalo+1:ny, :,:)
          global1(nx+shift, ny+1:ny+nhalo, :,:) = -global1(nx+shift, 1:nhalo,       :,:)
       end if

       if(ntile_per_pe > 1) write(type2, *)type, " at tile_count = ",n
       call mpp_clock_begin(id)
       if(ntile_per_pe == 1) then
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, gridtype=CGRID_NE, name=type2//' vector CGRID_NE')
       else
          call mpp_update_domains( x(:,:,:,n),  y(:,:,:,n),  domain, gridtype=CGRID_NE, &
               name=type2//' vector CGRID_NE', tile_count = n)
       end if
       call mpp_clock_end  (id)
    end do



    do n = 1, ntile_per_pe
       if(ntile_per_pe > 1) write(type2, *)type, " at tile_count = ",n
       call compare_checksums( x(isd:ied+shift,jsd:jed,:,n), global1(isd:ied+shift,jsd:jed,  :,n), &
                               trim(type2)//' CGRID_NE X')
       call compare_checksums( y(isd:ied,jsd:jed+shift,:,n), global2(isd:ied,  jsd:jed+shift,:,n), &
                               trim(type2)//' CGRID_NE Y')
    end do

    if(ntile_per_pe == 1) then
       call mpp_clock_begin(id)
       call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false., name=type//' vector CGRID_NE' )
       call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.false., name=type//' vector CGRID_NE')
       call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false., name=type//' vector CGRID_NE' )
       call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true. , name=type//' vector CGRID_NE')
       call mpp_clock_end  (id)

       call compare_checksums( x1(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X1')
       call compare_checksums( x2(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X2')
       call compare_checksums( x3(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X3')
       call compare_checksums( x4(isd:ied+shift,jsd:jed,:,1), global1(isd:ied+shift,jsd:jed,:,1), type//' CGRID_NE X4')
       call compare_checksums( y1(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y1')
       call compare_checksums( y2(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y2')
       call compare_checksums( y3(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y3')
       call compare_checksums( y4(isd:ied,jsd:jed+shift,:,1), global2(isd:ied,jsd:jed+shift,:,1), type//' CGRID_NE Y4')

       !--- arbitrary halo updates ---------------------------------------
       if(wide_halo_x ==0) then
          if(single_tile .or. four_tile .or.  (cubic_grid .and. same_layout) .or. folded_north ) then
             allocate(local1(isd:ied+shift,jsd:jed,      nz) )
             allocate(local2(isd:ied,      jsd:jed+shift,nz) )

             do wh = 1, whalo
                do eh = 1, ehalo
                   if(wh .NE. eh) cycle
                   do sh = 1, shalo
                      do nh = 1, nhalo
                         if(sh .NE. nh) cycle
                         local1(isd:ied+shift,jsd:jed,      :) = global1(isd:ied+shift,jsd:jed,      :,1)
                         local2(isd:ied,      jsd:jed+shift,:) = global2(isd:ied,      jsd:jed+shift,:,1)
                         x = 0.; y = 0.
                         x(isc:iec+shift,jsc:jec,      :,1) = global1(isc:iec+shift,jsc:jec,      :,1)
                         y(isc:iec,      jsc:jec+shift,:,1) = global2(isc:iec,      jsc:jec+shift,:,1)
                         call fill_halo_zero(local1, wh, eh, sh, nh, shift, 0, isc, iec, jsc, jec, isd, ied, jsd, jed)
                         call fill_halo_zero(local2, wh, eh, sh, nh, 0, shift, isc, iec, jsc, jec, isd, ied, jsd, jed)

                         write(type3,'(a,a,i2,a,i2,a,i2,a,i2)') trim(type), ' vector CGRID_NE with whalo = ', &
                              wh, ', ehalo = ',eh, ', shalo = ', sh, ', nhalo = ', nh
                         call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE, whalo=wh, ehalo=eh, &
                              shalo=sh, nhalo=nh, name=type3)
                         call compare_checksums( x(isd:ied+shift,jsd:jed, :,1),  local1, trim(type3)//' X' )
                         call compare_checksums( y(isd:ied,jsd:jed+shift, :,1),  local2, trim(type3)//' Y' )
                      end do
                   end do
                end do
             end do
             deallocate(local1, local2)
          end if
       endif
    end if

    deallocate(global1, global2, x, y, x1, x2, x3, x4, y1, y2, y3, y4, global1_all, global2_all)
    deallocate(layout2D, global_indices, pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 )

    if(wide_halo_x > 0) then
       whalo = whalo_save
       ehalo = ehalo_save
       shalo = shalo_save
       nhalo = nhalo_save
    endif
    nx = nx_save
    ny = ny_save

  end subroutine test_uniform_mosaic

  !#################################################################################
  subroutine update_domains_performance( type )
    character(len=*), intent(in) :: type

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
    real,    allocatable, dimension(:,:,:,:) :: x, x1, y, y1, x_save, y_save
    real,    allocatable, dimension(:,:,:,:) :: a, a1, b, b1
    real,    allocatable, dimension(:,:,:  ) :: a1_2D, b1_2D
    integer            :: id_update
    integer            :: id1, id2
    logical            :: folded_north
    logical            :: cubic_grid, single_tile, four_tile
    character(len=3)   :: text
    integer            :: nx_save, ny_save
    integer            :: id_single, id_update_single

    folded_north       = .false.
    cubic_grid         = .false.
    single_tile        = .false.
    four_tile          = .false.
    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Single-Tile' )   !--- single with cyclic along x- and y-direction
       single_tile = .true.
       ntiles = 1
       num_contact = 2
    case ( 'Folded-north' )
       ntiles = 1
       num_contact = 2
       folded_north = .true.
    case ( 'Four-Tile' ) !--- cyclic along both x- and y-direction.
       ntiles = 4
       num_contact = 8
       four_tile = .true.
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'update_domains_performance: for Cubic_grid mosaic, nx_cubic is zero, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'update_domains_performance: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
                  'No test is done for Cubic-Grid mosaic. ' )
          return
       endif

       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.

    case default
       call mpp_error(FATAL, 'update_domains_performance: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    if( mod(npes, ntiles) == 0 ) then
       npes_per_tile = npes/ntiles
       write(outunit,*)'NOTE from update_domains_performance ==> For Mosaic "', trim(type), &
                       '", each tile will be distributed over ', npes_per_tile, ' processors.'
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
       write(outunit,*)'NOTE from update_domains_performance ==> For Mosaic "', trim(type), &
                        '", there will be ', ntile_per_pe, ' tiles on each processor.'
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
       call mpp_error(NOTE,'update_domains_performance: npes should be multiple of ntiles or ' // &
            'ntiles should be multiple of npes. No test is done for '//trim(type) )
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
            name = type, symmetry = .false. )
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
                              name = type, symmetry = .false.  )
    else if( four_tile ) then
       call define_fourtile_mosaic(type, domain, (/nx,nx,nx,nx/), (/ny,ny,ny,ny/), global_indices, &
                                   layout2D, pe_start, pe_end, symmetry = .false. )
    else if( cubic_grid ) then
       call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
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
       call mpp_error(FATAL, "test_mpp_domains: num_fields must be a positive integer")
    endif

    id1 = mpp_clock_id( type, flags=MPP_CLOCK_SYNC)
    id_single = mpp_clock_id( type//' non-blocking', flags=MPP_CLOCK_SYNC)


    call mpp_clock_begin(id1)
    call mpp_update_domains( x, domain)
    call mpp_clock_end  (id1)

    call mpp_clock_begin(id_single)
    id_update_single =  mpp_start_update_domains(a, domain)
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( type//' group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( type//' group non-blocking', flags=MPP_CLOCK_SYNC )


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
             call mpp_update_domains( x1(:,:,:,l), domain, complete=l==num_fields )
          enddo
          call mpp_clock_end  (id1)

          ! non-blocking update
          call mpp_clock_begin(id2)
          if( n == 1 ) then
             do l = 1, num_fields
                if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), domain, complete=.false.)
                id_update =  mpp_start_update_domains(a1(:,:,:,l), domain, complete=l==num_fields)
             enddo
          else
             do l = 1, num_fields
                if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), domain, update_id=id_update, complete=.false.)
                id_update =  mpp_start_update_domains(a1(:,:,:,l), domain, update_id=id_update, complete=l==num_fields)
             enddo
          endif
          call mpp_clock_end  (id2)

          !---- sleep some time for non-blocking.
          if(do_sleep) call sleep(1)

          call mpp_clock_begin(id2)
          do l = 1, num_fields
             if(mix_2D_3D) call mpp_complete_update_domains(id_update, a1_2D(:,:,l), domain, complete=.false.)
             call mpp_complete_update_domains(id_update, a1(:,:,:,l), domain, complete=l==num_fields)
          enddo
          call mpp_clock_end  (id2)


          !--- compare checksum
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), type//' X'//text)
          enddo
          if(mix_2D_3D)call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), type//' X 2D')
       enddo
       deallocate(x1, a1)
       if(mix_2D_3D) deallocate(a1_2D)
    endif

    call mpp_clock_begin(id_single)
    call mpp_complete_update_domains(id_update_single, a, domain)
    call mpp_clock_end  (id_single)
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), type)
    deallocate(x, a, x_save)


    !------------------------------------------------------------------
    !              vector update : BGRID_NE, one extra point in each direction for cubic-grid
    !------------------------------------------------------------------
    !--- setup data
    shift = 0
    if(single_tile .or. four_tile .or. folded_north) then
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

    id1 = mpp_clock_id( trim(type)//' BGRID', flags=MPP_CLOCK_SYNC )
    id_single = mpp_clock_id( trim(type)//' BGRID non-blocking', flags=MPP_CLOCK_SYNC )

    call mpp_clock_begin(id1)
    call mpp_update_domains( x, y, domain, gridtype=BGRID_NE)
    call mpp_clock_end  (id1)

    !--- non-blocking update
    call mpp_clock_begin(id_single)
    id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=BGRID_NE)
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( trim(type)//' BGRID group', flags=MPP_CLOCK_SYNC)
    id2 = mpp_clock_id( trim(type)//' BGRID group non-blocking', flags=MPP_CLOCK_SYNC)
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
             call mpp_update_domains( x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=BGRID_NE, complete=l==num_fields )
          enddo
          call mpp_clock_end  (id1)

          !--- non-blocking update
          call mpp_clock_begin(id2)
          if( n == 1 ) then
             do l = 1, num_fields
                if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                             gridtype=BGRID_NE, complete=.false.)
                id_update =  mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, &
                             gridtype=BGRID_NE, complete=l==num_fields)
             enddo
          else
             do l = 1, num_fields
                if(mix_2D_3D) id_update =  mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, gridtype=BGRID_NE, &
                                update_id=id_update, complete=.false.)
                id_update =  mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, gridtype=BGRID_NE, &
                                update_id=id_update, complete=l==num_fields)
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
                                              gridtype=BGRID_NE, complete=l==num_fields)
          enddo
          call mpp_clock_end  (id2)

          !--- compare checksum
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), type//' BGRID X'//text)
             call compare_checksums( y1(:,:,:,l), b1(:,:,:,l), type//' BGRID Y'//text)
             if(mix_2D_3D) then
                call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), type//' BGRID X'//text)
                call compare_checksums( y1(:,:,:,1), b1(:,:,:,1), type//' BGRID Y'//text)
             endif
          enddo
          if(mix_2D_3D) then
             call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), type//' BGRID X 2D')
             call compare_checksums( y1(:,:,1,:), b1_2D(:,:,:), type//' BGRID Y 2D')
          endif
       enddo
       deallocate(x1, y1, a1, b1)
       if(mix_2D_3D) deallocate(a1_2D, b1_2D)
    endif

    call mpp_clock_begin(id_single)
    call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=BGRID_NE)
    call mpp_clock_end  (id_single)


    !--- compare checksum

    call compare_checksums( x(:,:,:,1), a(:,:,:,1), type//' BGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), type//' BGRID Y')


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

    id1 = mpp_clock_id( trim(type)//' CGRID', flags=MPP_CLOCK_SYNC )
    id_single = mpp_clock_id( trim(type)//' CGRID non-blocking', flags=MPP_CLOCK_SYNC )

    call mpp_clock_begin(id1)
    call mpp_update_domains( x, y, domain, gridtype=CGRID_NE)
    call mpp_clock_end  (id1)

    !--- non-blocking update
    call mpp_clock_begin(id_single)
    id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=CGRID_NE)
    call mpp_clock_end  (id_single)

    !---- sleep some time for non-blocking.
    if(do_sleep) call sleep(1)

    id1 = mpp_clock_id( trim(type)//' CGRID group', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( trim(type)//' CGRID group non-blocking', flags=MPP_CLOCK_SYNC )

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
             call mpp_update_domains( x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE, complete=l==num_fields )
          enddo
          call mpp_clock_end  (id1)

          !--- non-blocking update
          call mpp_clock_begin(id2)
          if( n == 1 ) then
             do l = 1, num_fields
                if(mix_2D_3D) id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                       gridtype=CGRID_NE, complete=.false.)
                id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, &
                       gridtype=CGRID_NE, complete=l==num_fields)
             enddo
          else
             do l = 1, num_fields
             if(mix_2D_3D)id_update = mpp_start_update_domains(a1_2D(:,:,l), b1_2D(:,:,l), domain, gridtype=CGRID_NE, &
                            update_id=id_update, complete=.false.)
             id_update = mpp_start_update_domains(a1(:,:,:,l), b1(:,:,:,l), domain, gridtype=CGRID_NE, &
                            update_id=id_update, complete=l==num_fields)
             enddo
          endif
          call mpp_clock_end  (id2)

          !---- sleep some time for non-blocking.
          if(do_sleep) call sleep(1)

          call mpp_clock_begin(id2)
          do l = 1, num_fields
             if(mix_2D_3D)call mpp_complete_update_domains(id_update, a1_2D(:,:,l), b1_2D(:,:,l), domain, &
                   gridtype=CGRID_NE, complete=.false.)
             call mpp_complete_update_domains(id_update, a1(:,:,:,l), b1(:,:,:,l), domain, &
                   gridtype=CGRID_NE, complete=l==num_fields)
          enddo
          call mpp_clock_end  (id2)

          !--- compare checksum
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums( x1(:,:,:,l), a1(:,:,:,l), type//' CGRID X'//text)
             call compare_checksums( y1(:,:,:,l), b1(:,:,:,l), type//' CGRID Y'//text)
          enddo
          if(mix_2D_3D) then
             call compare_checksums( x1(:,:,1,:), a1_2D(:,:,:), type//' BGRID X 2D')
             call compare_checksums( y1(:,:,1,:), b1_2D(:,:,:), type//' BGRID Y 2D')
          endif
       enddo
       deallocate(x1, y1, a1, b1)
       if(mix_2D_3D) deallocate(a1_2D, b1_2D)
    endif

    call mpp_clock_begin(id_single)
    call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=CGRID_NE)
    call mpp_clock_end  (id_single)

    !--- compare checksum

    call compare_checksums( x(:,:,:,1), a(:,:,:,1), type//' CGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), type//' CGRID Y')

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

    call mpp_update_domains( x, y, domain, gridtype=AGRID)

    id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=AGRID)
    call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=AGRID)

    !--- compare checksum
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), type//' AGRID X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), type//' AGRID Y')

    x = x_save; y = y_save
    a = x_save; b = y_save

    call mpp_update_domains( x, y, domain, gridtype=AGRID, flags = SCALAR_PAIR)

    id_update_single =  mpp_start_update_domains(a, b, domain, gridtype=AGRID, flags = SCALAR_PAIR)
    call mpp_complete_update_domains(id_update_single, a, b, domain, gridtype=AGRID, flags = SCALAR_PAIR)

    !--- compare checksum
    call compare_checksums( x(:,:,:,1), a(:,:,:,1), type//' AGRID SCALAR-PAIR X')
    call compare_checksums( y(:,:,:,1), b(:,:,:,1), type//' AGRID SCALAR-PAIR Y')

    deallocate(x, y, a, b, x_save, y_save)

    nx = nx_save
    ny = ny_save

    deallocate(layout2D, global_indices, pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2 )


  end subroutine update_domains_performance


  !###############################################################
  subroutine test_mpp_global_sum( type )
    character(len=*), intent(in) :: type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile
    integer        :: i, j, k, l, n, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem

    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:,:)   :: data_3D
    real,    allocatable, dimension(:,:)     :: data_2D

    integer(kind=8)    :: mold
    logical            :: folded_north, cubic_grid
    character(len=3)   :: text
    integer            :: nx_save, ny_save
    integer            :: id1, id2, id3, id4
    real               :: gsum1, gsum2, gsum3, gsum4

    folded_north       = .false.
    cubic_grid         = .false.

    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Folded-north' )
       ntiles = 1
       shift = 0
       num_contact = 2
       folded_north = .true.
       npes_per_tile = npes
       if(layout_tripolar(1)*layout_tripolar(2) == npes ) then
          layout = layout_tripolar
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_group_update: for Cubic_grid mosaic, nx_cubic is zero, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_group_update: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       shift = 1
       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( mod(npes, ntiles) == 0 ) then
          npes_per_tile = npes/ntiles
          write(outunit,*)'NOTE from test_mpp_global_sum ==> For Mosaic "', trim(type), &
               '", each tile will be distributed over ', npes_per_tile, ' processors.'
       else
          call mpp_error(NOTE,'test_group_update: npes should be multiple of ntiles No test is done for '//trim(type))
          return
       endif
       if(layout_cubic(1)*layout_cubic(2) == npes_per_tile) then
          layout = layout_cubic
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
    case default
       call mpp_error(FATAL, 'test_mpp_global_sum: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       pe_start(n) = (n-1)*npes_per_tile
       pe_end(n)   = n*npes_per_tile-1
    end do

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    if(folded_north) then
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
            name = type, symmetry = .false.  )
    else if( cubic_grid ) then
       call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    allocate(data_2d(isd:ied,jsd:jed))
    allocate(data_3d(isd:ied,jsd:jed,nz))

    do k = 1, nz
       do j = jsd, jed
          do i = isd, ied
             data_3d(i,j,k) =  k*1e3 + i + j*1e-3
          enddo
       enddo
    enddo

    do j = jsd, jed
       do i = isd, ied
          data_2d(i,j) =  i*1e3 + j*1e-3
       enddo
    enddo

    id1 = mpp_clock_id( type//' bitwise sum 3D', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( type//' EFP sum 3D', flags=MPP_CLOCK_SYNC )
    id3 = mpp_clock_id( type//' EFP sum 3D check', flags=MPP_CLOCK_SYNC )
    id4 = mpp_clock_id( type//' non-bitwise sum 3D', flags=MPP_CLOCK_SYNC )

    call mpp_clock_begin(id1)
    do n = 1, num_iter
       gsum1 = mpp_global_sum(domain, data_3d, flags=BITWISE_EXACT_SUM)
    enddo
    call mpp_clock_end(id1)

    call mpp_clock_begin(id2)
    do n = 1, num_iter
       gsum2 = mpp_global_sum(domain, data_3d, flags=BITWISE_EFP_SUM)
    enddo
    call mpp_clock_end(id2)

    call mpp_clock_begin(id3)
    do n = 1, num_iter
       gsum3 = mpp_global_sum(domain, data_3d, flags=BITWISE_EFP_SUM, overflow_check=.true. )
    enddo
    call mpp_clock_end(id3)

    call mpp_clock_begin(id4)
    do n = 1, num_iter
       gsum4= mpp_global_sum(domain, data_3d)
    enddo
    call mpp_clock_end(id4)

    write(outunit, *) " ********************************************************************************"
    write(outunit, *) " global sum for "//type//' bitwise exact sum 3D = ', gsum1
    write(outunit, *) " global sum for "//type//' bitwise EFP sum 3D = ', gsum2
    write(outunit, *) " global sum for "//type//' bitwise EFP sum 3D with overflow_check = ', gsum3
    write(outunit, *) " global sum for "//type//' non-bitwise sum 3D = ', gsum4
    write(outunit, *) " "
    write(outunit, *) " chksum for "//type//' bitwise exact sum 3D = ', transfer(gsum1, mold)
    write(outunit, *) " chksum for "//type//' bitwise EFP sum 3D = ', transfer(gsum2, mold)
    write(outunit, *) " chksum for "//type//' bitwise EFP sum 3D with overflow_check = ', transfer(gsum3, mold)
    write(outunit, *) " chksum for "//type//' non-bitwise sum 3D = ', transfer(gsum4, mold)
    write(outunit, *) " ********************************************************************************"

    id1 = mpp_clock_id( type//' bitwise sum 2D', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( type//' EFP sum 2D', flags=MPP_CLOCK_SYNC )
    id3 = mpp_clock_id( type//' EFP sum 2D check', flags=MPP_CLOCK_SYNC )
    id4 = mpp_clock_id( type//' non-bitwise sum 2D', flags=MPP_CLOCK_SYNC )

    call mpp_clock_begin(id1)
    do n = 1, num_iter
       gsum1 = mpp_global_sum(domain, data_2d, flags=BITWISE_EXACT_SUM)
    enddo
    call mpp_clock_end(id1)

    call mpp_clock_begin(id2)
    do n = 1, num_iter
       gsum2 = mpp_global_sum(domain, data_2d, flags=BITWISE_EFP_SUM)
    enddo
    call mpp_clock_end(id2)

    call mpp_clock_begin(id3)
    do n = 1, num_iter
       gsum3 = mpp_global_sum(domain, data_2d, flags=BITWISE_EFP_SUM, overflow_check=.true. )
    enddo
    call mpp_clock_end(id3)

    call mpp_clock_begin(id4)
    do n = 1, num_iter
       gsum4= mpp_global_sum(domain, data_2d)
    enddo
    call mpp_clock_end(id4)

    write(outunit, *) " ********************************************************************************"
    write(outunit, *) " global sum for "//type//' bitwise exact sum 2D = ', gsum1
    write(outunit, *) " global sum for "//type//' bitwise EFP sum 2D = ', gsum2
    write(outunit, *) " global sum for "//type//' bitwise EFP sum 2D with overflow_check = ', gsum3
    write(outunit, *) " global sum for "//type//' non-bitwise sum 2D = ', gsum4
    write(outunit, *) " "
    write(outunit, *) " chksum for "//type//' bitwise exact sum 2D = ', transfer(gsum1, mold)
    write(outunit, *) " chksum for "//type//' bitwise EFP sum 2D = ', transfer(gsum2, mold)
    write(outunit, *) " chksum for "//type//' bitwise EFP sum 2D with overflow_check = ', transfer(gsum3, mold)
    write(outunit, *) " chksum for "//type//' non-bitwise sum 2D = ', transfer(gsum4, mold)
    write(outunit, *) " ********************************************************************************"



    nx = nx_save
    ny = ny_save

  end subroutine test_mpp_global_sum

  !###############################################################
  subroutine test_group_update( type )
    character(len=*), intent(in) :: type

    type(domain2D) :: domain
    integer        :: num_contact, ntiles, npes_per_tile
    integer        :: i, j, k, l, n, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem

    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:,:,:) :: x1, y1, x2, y2
    real,    allocatable, dimension(:,:,:,:) :: a1, a2
    real,    allocatable, dimension(:,:,:)   :: base
    integer            :: id1, id2, id3
    logical            :: folded_north
    logical            :: cubic_grid
    character(len=3)   :: text
    integer            :: nx_save, ny_save
    type(mpp_group_update_type) :: group_update
    type(mpp_group_update_type), allocatable :: update_list(:)

    folded_north       = .false.
    cubic_grid         = .false.

    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Folded-north' )
       ntiles = 1
       shift = 0
       num_contact = 2
       folded_north = .true.
       npes_per_tile = npes
       if(layout_tripolar(1)*layout_tripolar(2) == npes ) then
          layout = layout_tripolar
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_group_update: for Cubic_grid mosaic, nx_cubic is zero, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_group_update: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       shift = 1
       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( mod(npes, ntiles) == 0 ) then
          npes_per_tile = npes/ntiles
          write(outunit,*)'NOTE from update_domains_performance ==> For Mosaic "', trim(type), &
               '", each tile will be distributed over ', npes_per_tile, ' processors.'
       else
          call mpp_error(NOTE,'test_group_update: npes should be multiple of ntiles No test is done for '//trim(type))
          return
       endif
       if(layout_cubic(1)*layout_cubic(2) == npes_per_tile) then
          layout = layout_cubic
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
    case default
       call mpp_error(FATAL, 'test_group_update: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       pe_start(n) = (n-1)*npes_per_tile
       pe_end(n)   = n*npes_per_tile-1
    end do

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    if(folded_north) then
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
            name = type, symmetry = .false.  )
    else if( cubic_grid ) then
       call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )

    if(num_fields<1) then
       call mpp_error(FATAL, "test_mpp_domains: num_fields must be a positive integer")
    endif

    allocate(update_list(num_fields))

    id1 = mpp_clock_id( type//' group 2D', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( type//' non-group 2D', flags=MPP_CLOCK_SYNC )
    id3 = mpp_clock_id( type//' non-block group 2D', flags=MPP_CLOCK_SYNC )

    allocate( a1(ism:iem,      jsm:jem,       nz, num_fields) )
    allocate( x1(ism:iem+shift,jsm:jem,       nz, num_fields) )
    allocate( y1(ism:iem,      jsm:jem+shift, nz, num_fields) )
    allocate( a2(ism:iem,      jsm:jem,       nz, num_fields) )
    allocate( x2(ism:iem+shift,jsm:jem,       nz, num_fields) )
    allocate( y2(ism:iem,      jsm:jem+shift, nz, num_fields) )
    allocate( base(isc:iec+shift,jsc:jec+shift,nz) )
    a1 = 0; x1 = 0; y1 = 0

    base = 0
    do k = 1,nz
       do j = jsc, jec+shift
          do i = isc, iec+shift
             base(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    !--- Test for partial direction update
    do l =1, num_fields
       call mpp_create_group_update(group_update, a1(:,:,:,l), domain, flags=WUPDATE+SUPDATE)
    end do

    do l = 1, num_fields
       a1(isc:iec,jsc:jec,:,l) = base(isc:iec,jsc:jec,:) + l*1e3
       do k=1,nz
          do i=isc-1,iec+1
             a1(i,jsc-1,k,l) = 999;
             a1(i,jec+1,k,l) = 999;
          enddo
          do j=jsc,jec
             a1(isc-1,j,k,l) = 999
             a1(iec+1,j,k,l) = 999
          enddo
       enddo
    enddo

    a2 = a1
    call mpp_do_group_update(group_update, domain, a1(isc,jsc,1,1))

    do l = 1, num_fields
       call mpp_update_domains( a2(:,:,:,l), domain, flags=WUPDATE+SUPDATE, complete=l==num_fields )
    enddo

    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(isd:ied,jsd:jed,:,l),a2(isd:ied,jsd:jed,:,l),type//' CENTER South West '//text)
    enddo

    call mpp_clear_group_update(group_update)

    !--- Test for DGRID update
    if(type == 'Cubic-Grid' ) then
       x1 = 0; y1 = 0
       do l =1, num_fields
          call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=DGRID_NE)
       end do

       do l = 1, num_fields
          y1(isc:iec+shift,jsc:jec,      :,l) = base(isc:iec+shift,jsc:jec,      :) + l*1e3 + 1e6
          x1(isc:iec,      jsc:jec+shift,:,l) = base(isc:iec,      jsc:jec+shift,:) + l*1e3 + 2*1e6
       enddo
       x2 = x1; y2 = y1
       call mpp_start_group_update(group_update, domain, x1(isc,jsc,1,1))
       call mpp_complete_group_update(group_update, domain, x1(isc,jsc,1,1))

       do l = 1, num_fields
          call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=DGRID_NE, complete=l==num_fields )
       enddo

    !--- compare checksum
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(x1(isd:ied+shift,jsd:jed,      :,l),x2(isd:ied+shift,jsd:jed,      :,l),type//' DGRID X'//text)
          call compare_checksums(y1(isd:ied,      jsd:jed+shift,:,l),y2(isd:ied,      jsd:jed+shift,:,l),type//' DGRID Y'//text)
       enddo

       call mpp_clear_group_update(group_update)
    endif
    !--- Test for CGRID
    a1 = 0; x1 = 0; y1 = 0
    do l =1, num_fields
       call mpp_create_group_update(group_update, a1(:,:,:,l), domain)
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE)
    end do

    do n = 1, num_iter
       a1 = 0; x1 = 0; y1 = 0
       do l = 1, num_fields
          a1(isc:iec,      jsc:jec,      :,l) = base(isc:iec,      jsc:jec,      :) + l*1e3
          x1(isc:iec+shift,jsc:jec,      :,l) = base(isc:iec+shift,jsc:jec,      :) + l*1e3 + 1e6
          y1(isc:iec,      jsc:jec+shift,:,l) = base(isc:iec,      jsc:jec+shift,:) + l*1e3 + 2*1e6
       enddo
       a2 = a1; x2 = x1; y2 = y1
       call mpp_clock_begin(id1)
       call mpp_do_group_update(group_update, domain, a1(isc,jsc,1,1))
       call mpp_clock_end  (id1)

       call mpp_clock_begin(id2)
       do l = 1, num_fields
          call mpp_update_domains( a2(:,:,:,l), domain, complete=l==num_fields )
       enddo
       do l = 1, num_fields
          call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=CGRID_NE, complete=l==num_fields )
       enddo
       call mpp_clock_end(id2)

       !--- compare checksum
       if( n == num_iter ) then
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(a1(isd:ied,      jsd:jed,      :,l),a2(isd:ied,      jsd:jed,      :,l),type//' CENTER '//text)
          call compare_checksums(x1(isd:ied+shift,jsd:jed,      :,l),x2(isd:ied+shift,jsd:jed,      :,l),type//' CGRID X'//text)
          call compare_checksums(y1(isd:ied,      jsd:jed+shift,:,l),y2(isd:ied,      jsd:jed+shift,:,l),type//' CGRID Y'//text)
       enddo
       endif
       a1 = 0; x1 = 0; y1 = 0
       do l = 1, num_fields
          a1(isc:iec,      jsc:jec,      :,l) = base(isc:iec,      jsc:jec,      :) + l*1e3
          x1(isc:iec+shift,jsc:jec,      :,l) = base(isc:iec+shift,jsc:jec,      :) + l*1e3 + 1e6
          y1(isc:iec,      jsc:jec+shift,:,l) = base(isc:iec,      jsc:jec+shift,:) + l*1e3 + 2*1e6
       enddo
       call mpp_clock_begin(id3)
       call mpp_start_group_update(group_update, domain, a1(isc,jsc,1,1))
       call mpp_complete_group_update(group_update, domain, a1(isc,jsc,1,1))
       call mpp_clock_end  (id3)
       !--- compare checksum
       if( n == num_iter ) then
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(a1(isd:ied,      jsd:jed,      :,l),a2(isd:ied,      jsd:jed,      :,l), &
                                 type//' nonblock CENTER '//text)
          call compare_checksums(x1(isd:ied+shift,jsd:jed,      :,l),x2(isd:ied+shift,jsd:jed,      :,l), &
                                 type//' nonblock CGRID X'//text)
          call compare_checksums(y1(isd:ied,      jsd:jed+shift,:,l),y2(isd:ied,      jsd:jed+shift,:,l), &
                                 type//' nonblock CGRID Y'//text)
       enddo
       endif
    enddo

    call mpp_clear_group_update(group_update)

    !--- The following is to test overlapping start and complete
    if( num_fields > 1 ) then
       do l =1, num_fields
          call mpp_create_group_update(update_list(l), a1(:,:,:,l), domain)
          call mpp_create_group_update(update_list(l), x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE)
       end do

       do n = 1, num_iter

          a1 = 0; x1 = 0; y1 = 0
          do l = 1, num_fields
             a1(isc:iec,      jsc:jec,      :,l) = base(isc:iec,      jsc:jec,      :) + l*1e3
             x1(isc:iec+shift,jsc:jec,      :,l) = base(isc:iec+shift,jsc:jec,      :) + l*1e3 + 1e6
             y1(isc:iec,      jsc:jec+shift,:,l) = base(isc:iec,      jsc:jec+shift,:) + l*1e3 + 2*1e6
          enddo
          do l = 1, num_fields-1
             call mpp_start_group_update(update_list(l), domain, a1(isc,jsc,1,1))
          enddo

          call mpp_complete_group_update(update_list(1), domain, a1(isc,jsc,1,1))
          call mpp_start_group_update(update_list(num_fields), domain, a1(isc,jsc,1,1))
          do l = 2, num_fields
             call mpp_complete_group_update(update_list(l), domain, a1(isc,jsc,1,1))
          enddo
          !--- compare checksum
          if( n == num_iter ) then
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums(a1(isd:ied,      jsd:jed,      :,l),a2(isd:ied,      jsd:jed,      :,l), &
                                    type//' multiple nonblock CENTER '//text)
             call compare_checksums(x1(isd:ied+shift,jsd:jed,      :,l),x2(isd:ied+shift,jsd:jed,      :,l), &
                                    type//' multiple nonblock CGRID X'//text)
             call compare_checksums(y1(isd:ied,      jsd:jed+shift,:,l),y2(isd:ied,      jsd:jed+shift,:,l), &
                                    type//' multiple nonblock CGRID Y'//text)
          enddo
          endif
       enddo
    endif

    do l =1, num_fields
      call mpp_clear_group_update(update_list(l))
    enddo
    deallocate(update_list)

    !--- test scalar 4-D variable
    call mpp_create_group_update(group_update, a1(:,:,:,:), domain)

    a1 = 0; x1 = 0; y1 = 0
    do l = 1, num_fields
       a1(isc:iec,      jsc:jec,      :,l) = base(isc:iec,      jsc:jec,      :) + l*1e3
    enddo
    a2 = a1; x2 = x1; y2 = y1
    call mpp_clock_begin(id1)
    call mpp_do_group_update(group_update, domain, a1(isc,jsc,1,1))
    call mpp_clock_end  (id1)

    call mpp_clock_begin(id2)
    call mpp_update_domains( a2(:,:,:,:), domain )
    call mpp_clock_end(id2)

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(isd:ied, jsd:jed, :,l),a2(isd:ied, jsd:jed, :,l),type//' 4D CENTER '//text)
    enddo

    a1 = 0
    do l = 1, num_fields
       a1(isc:iec,      jsc:jec,      :,l) = base(isc:iec,      jsc:jec,      :) + l*1e3
    enddo
    call mpp_clock_begin(id3)
    call mpp_start_group_update(group_update, domain, a1(isc,jsc,1,1))
    call mpp_complete_group_update(group_update, domain, a1(isc,jsc,1,1))
    call mpp_clock_end  (id3)

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(isd:ied, jsd:jed, :,l),a2(isd:ied, jsd:jed, :,l), &
                              type//' nonblock 4D CENTER '//text)
    enddo



    !--- test for BGRID.
    deallocate(a1, x1, y1)
    deallocate(a2, x2, y2)
    call mpp_clear_group_update(group_update)

    allocate( a1(ism:iem+shift,jsm:jem+shift, nz, num_fields) )
    allocate( x1(ism:iem+shift,jsm:jem+shift, nz, num_fields) )
    allocate( y1(ism:iem+shift,jsm:jem+shift, nz, num_fields) )
    allocate( a2(ism:iem+shift,jsm:jem+shift, nz, num_fields) )
    allocate( x2(ism:iem+shift,jsm:jem+shift, nz, num_fields) )
    allocate( y2(ism:iem+shift,jsm:jem+shift, nz, num_fields) )

    do l =1, num_fields
       call mpp_create_group_update(group_update, a1(:,:,:,l), domain, position=CORNER)
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=BGRID_NE)
    end do

    do n = 1, num_iter
       a1 = 0; x1 = 0; y1 = 0
       do l = 1, num_fields
          a1(isc:iec+shift,jsc:jec+shift,:,l) = base(isc:iec+shift,jsc:jec+shift,:) + l*1e3
          x1(isc:iec+shift,jsc:jec+shift,:,l) = base(isc:iec+shift,jsc:jec+shift,:) + l*1e3 + 1e6
          y1(isc:iec+shift,jsc:jec+shift,:,l) = base(isc:iec+shift,jsc:jec+shift,:) + l*1e3 + 2*1e6
       enddo
       a2 = a1; x2 = x1; y2 = y1
       call mpp_clock_begin(id1)
       call mpp_do_group_update(group_update, domain, a1(isc,jsc,1,1))
       call mpp_clock_end  (id1)

       call mpp_clock_begin(id2)
       do l = 1, num_fields
          call mpp_update_domains( a2(:,:,:,l), domain, position=CORNER, complete=l==num_fields )
       enddo
       do l = 1, num_fields
          call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=BGRID_NE, complete=l==num_fields )
       enddo
       call mpp_clock_end(id2)

       !--- compare checksum
       if( n == num_iter ) then
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(a1(isd:ied+shift,jsd:jed+shift,:,l),a2(isd:ied+shift,jsd:jed+shift,:,l),type//' CORNER '//text)
          call compare_checksums(x1(isd:ied+shift,jsd:jed+shift,:,l),x2(isd:ied+shift,jsd:jed+shift,:,l),type//' BGRID X'//text)
          call compare_checksums(y1(isd:ied+shift,jsd:jed+shift,:,l),y2(isd:ied+shift,jsd:jed+shift,:,l),type//' BGRID Y'//text)
       enddo
       endif

       a1 = 0; x1 = 0; y1 = 0
       do l = 1, num_fields
          a1(isc:iec+shift,jsc:jec+shift,:,l) = base(isc:iec+shift,jsc:jec+shift,:) + l*1e3
          x1(isc:iec+shift,jsc:jec+shift,:,l) = base(isc:iec+shift,jsc:jec+shift,:) + l*1e3 + 1e6
          y1(isc:iec+shift,jsc:jec+shift,:,l) = base(isc:iec+shift,jsc:jec+shift,:) + l*1e3 + 2*1e6
       enddo
       call mpp_clock_begin(id3)
       call mpp_start_group_update(group_update, domain, a1(isc,jsc,1,1))
       call mpp_complete_group_update(group_update, domain, a1(isc,jsc,1,1))
       call mpp_clock_end  (id3)
       !--- compare checksum
       if( n == num_iter ) then
       do l = 1, num_fields
          write(text, '(i3.3)') l
          call compare_checksums(a1(isd:ied+shift,jsd:jed+shift,:,l),a2(isd:ied+shift,jsd:jed+shift,:,l), &
                                 type//' nonblockCORNER '//text)
          call compare_checksums(x1(isd:ied+shift,jsd:jed+shift,:,l),x2(isd:ied+shift,jsd:jed+shift,:,l), &
                                 type//' nonblock BGRID X'//text)
          call compare_checksums(y1(isd:ied+shift,jsd:jed+shift,:,l),y2(isd:ied+shift,jsd:jed+shift,:,l), &
                                 type//' nonblock BGRID Y'//text)
       enddo
       endif
    enddo

    call mpp_clear_group_update(group_update)

    !-----------------------------------------------------------------------------
    !                   test for AGRID vector and scalar pair
    !-----------------------------------------------------------------------------
    deallocate(x1, y1)
    deallocate(x2, y2)

    allocate( x1(ism:iem,jsm:jem, nz, num_fields) )
    allocate( y1(ism:iem,jsm:jem, nz, num_fields) )
    allocate( x2(ism:iem,jsm:jem, nz, num_fields) )
    allocate( y2(ism:iem,jsm:jem, nz, num_fields) )

    x1 = 0; y1 = 0
    do l = 1, num_fields
       x1(isc:iec,jsc:jec,:,l) = base(isc:iec,jsc:jec,:) + l*1e3 + 1e6
       y1(isc:iec,jsc:jec,:,l) = base(isc:iec,jsc:jec,:) + l*1e3 + 2*1e6
    enddo
    x2 = x1; y2 = y1

    do l =1, num_fields
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=AGRID)
    end do

    do l = 1, num_fields
       call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=AGRID, complete=l==num_fields )
    enddo

    call mpp_start_group_update(group_update, domain, a1(isc,jsc,1,1))
    call mpp_complete_group_update(group_update, domain, a1(isc,jsc,1,1))

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(x1(isd:ied,jsd:jed,:,l),x2(isd:ied,jsd:jed,:,l),type//' AGRID X'//text)
       call compare_checksums(y1(isd:ied,jsd:jed,:,l),y2(isd:ied,jsd:jed,:,l),type//' AGRID Y'//text)
    enddo

    call mpp_clear_group_update(group_update)

    x1 = 0; y1 = 0
    do l = 1, num_fields
       x1(isc:iec,jsc:jec,:,l) = base(isc:iec,jsc:jec,:) + l*1e3 + 1e6
       y1(isc:iec,jsc:jec,:,l) = base(isc:iec,jsc:jec,:) + l*1e3 + 2*1e6
    enddo
    x2 = x1; y2 = y1

    do l =1, num_fields
       call mpp_create_group_update(group_update, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=AGRID, flags=SCALAR_PAIR)
    end do

    do l = 1, num_fields
       call mpp_update_domains( x2(:,:,:,l), y2(:,:,:,l), domain, gridtype=AGRID, flags=SCALAR_PAIR, complete=l==num_fields)
    enddo

    call mpp_start_group_update(group_update, domain, x1(isc,jsc,1,1))
    call mpp_complete_group_update(group_update, domain, x1(isc,jsc,1,1))

    !--- compare checksum
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(x1(isd:ied,jsd:jed,:,l),x2(isd:ied,jsd:jed,:,l),type//' AGRID SCALAR_PAIR X'//text)
       call compare_checksums(y1(isd:ied,jsd:jed,:,l),y2(isd:ied,jsd:jed,:,l),type//' AGRID SCALAR_PAIR Y'//text)
    enddo

    call mpp_clear_group_update(group_update)

    deallocate(pe_start, pe_end, tile1, tile2)
    deallocate(istart1, iend1, jstart1, jend1)
    deallocate(istart2, iend2, jstart2, jend2)
    deallocate(layout2D, global_indices)

    deallocate(a1, x1, y1)
    deallocate(a2, x2, y2)
    deallocate(base)
    call mpp_deallocate_domain(domain)

end subroutine test_group_update


  !###############################################################
  !--- This will test scalar and CGRID performance between halo=1 and halo=3
  subroutine test_halosize_update( type )
    character(len=*), intent(in) :: type

    type(domain2D) :: domain
    integer        :: ntiles, npes_per_tile
    integer        :: i, j, k, l, n, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem

    integer, allocatable, dimension(:)       :: pe_start, pe_end
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:,:,:) :: x1, y1, x2, y2
    real,    allocatable, dimension(:,:,:,:) :: a1, a2
    real,    allocatable, dimension(:,:,:,:) :: base, global_all, global
    real,    allocatable, dimension(:,:,:,:) :: base1, global1_all, global1
    real,    allocatable, dimension(:,:,:,:) :: base2, global2_all, global2

    integer            :: id1, id2, id3
    logical            :: folded_north
    logical            :: cubic_grid, is_symmetry
    character(len=3)   :: text
    character(len=1)   :: halostr
    integer            :: nx_save, ny_save
    integer            :: mytile
    type(mpp_group_update_type) :: group_update1, group_update2
    type(mpp_group_update_type), allocatable :: update_list(:)

    if(whalo .ne. ehalo .or. whalo .ne. shalo .or. whalo .ne. nhalo) then
       call mpp_error(FATAL,"test_mpp_domains: whalo, ehalo, shalo, nhalo must be the same when test_halosize_performance=true")
    endif

    folded_north       = .false.
    cubic_grid         = .false.

    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Folded-north', 'Folded-north symmetry' )
       ntiles = 1
       mytile = 1
       folded_north = .true.
       npes_per_tile = npes
       if(layout_tripolar(1)*layout_tripolar(2) == npes ) then
          layout = layout_tripolar
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
       if(index(type, 'symmetry') == 0) then
          is_symmetry = .false.
       else
          is_symmetry = .true.
       end if
    case ( 'Cubic-Grid' )
       is_symmetry = .true.
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_halosize_update: for Cubic_grid mosaic, nx_cubic is zero, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_halosize_update: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       if( mod(npes, ntiles) .ne. 0 ) then
          call mpp_error(NOTE,'test_halosize_update: npes is not divisible by ntiles, no test is done for '//trim(type) )
          return
       endif
       npes_per_tile = npes/ntiles
       mytile = mpp_pe()/npes_per_tile + 1
       cubic_grid = .true.
       if(layout_cubic(1)*layout_cubic(2) == npes_per_tile) then
          layout = layout_cubic
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
    case default
       call mpp_error(FATAL, 'test_group_update: no such test: '//type)
    end select

    shift = 0
    if(is_symmetry) shift = 1

    !--- define domain
    if(folded_north) then
        call mpp_define_domains((/1,nx,1,ny/), layout, domain, &
                          xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=is_symmetry, name=type )
    else if( cubic_grid ) then
       allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
       do n = 1, ntiles
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do

       do n = 1, ntiles
          global_indices(:,n) = (/1,nx,1,ny/)
          layout2D(:,n)         = layout
       end do

       call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end, use_memsize=.false.)
       deallocate(pe_start, pe_end)
       deallocate(layout2D, global_indices)
    endif

    !--- setup data
    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    if(num_fields<1) then
       call mpp_error(FATAL, "test_mpp_domains: num_fields must be a positive integer")
    endif

    !--- scalar update
    write(halostr,'(I1)') whalo
    id1 = mpp_clock_id( type//' halo='//halostr//' scalar', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( type//' halo=1 scalar', flags=MPP_CLOCK_SYNC )

    allocate( a1(isd:ied,      jsd:jed,       nz, num_fields) )
    allocate( a2(isd:ied,      jsd:jed,       nz, num_fields) )
    allocate(base(isc:iec, jsc:jec, nz, num_fields))
    allocate(global_all(1:nx,1:ny,nz,ntiles) )
    allocate(global(1-whalo:nx+ehalo, 1-shalo:ny+nhalo, nz, num_fields))

    do n = 1, ntiles
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                global_all(i,j,k,n) = n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do
    global = 0.0
    do l = 1, num_fields
       global(1:nx,1:ny,:,l) = global_all(:,:,:,mytile)
    enddo

    base(isc:iec,jsc:jec,:,:) = global(isc:iec,jsc:jec,:,:)

    !--- fill up the value at halo points
    do l = 1, num_fields
       if(folded_north) then
          call fill_folded_north_halo(global(:,:,:,l), 0, 0, 0, 0, 1)
       else if(cubic_grid) then
          call fill_cubic_grid_halo(global(:,:,:,l), global_all, global_all, mytile, 0, 0, 1, 1 )
       endif
    enddo
    a1 = 0.0
    a2(isd:ied,jsd:jed,:,:) = global(isd:ied,jsd:jed,:,:)

    do l =1, num_fields
       call mpp_create_group_update(group_update1, a1(:,:,:,l), domain)
    end do
    do l =1, num_fields
       call mpp_create_group_update(group_update2, a1(:,:,:,l), domain, whalo=1, ehalo=1, shalo=1, nhalo=1)
    end do

    do n = 1, num_iter
       a1 = 0.0
       a1(isc:iec,jsc:jec,:,:) = base(isc:iec,jsc:jec,:,:)

       call mpp_clock_begin(id1)
       call mpp_do_group_update(group_update1, domain, a1(isc,jsc,1,1))
       call mpp_clock_end(id1)
       if(n==num_iter) then
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums(a1(:,:,:,l),a2(:,:,:,l),type//' halo='//halostr//' scalar'//text)
          enddo
       endif
    enddo

    !--- make sure mpp_start_group_update/mpp_complete_group_update is OK
    a1 = 0.0
    a1(isc:iec,jsc:jec,:,:) = base(isc:iec,jsc:jec,:,:)

    call mpp_start_group_update(group_update1, domain, a1(isc,jsc,1,1))
    call mpp_complete_group_update(group_update1, domain, a1(isc,jsc,1,1))
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(:,:,:,l),a2(:,:,:,l),type//'nonblock halo='//halostr//' scalar'//text)
    enddo


    a2 = 0
    a2(isc-1:iec+1,jsc-1:jec+1,:,:) = global(isc-1:iec+1,jsc-1:jec+1,:,:)

    do n = 1, num_iter
       a1 = 0.0
       a1(isc:iec,jsc:jec,:,:) = base(isc:iec,jsc:jec,:,:)
       call mpp_clock_begin(id2)
       call mpp_do_group_update(group_update2, domain, a1(isc,jsc,1,1))
       call mpp_clock_end(id2)
       if(n==num_iter) then
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums(a1(:,:,:,l),a2(:,:,:,l),type//' halo=1 scalar'//text)
          enddo
       endif
    enddo

    a1 = 0.0
    a1(isc:iec,jsc:jec,:,:) = base(isc:iec,jsc:jec,:,:)
    call mpp_start_group_update(group_update2, domain, a1(isc,jsc,1,1))
    call mpp_complete_group_update(group_update2, domain, a1(isc,jsc,1,1))
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(a1(:,:,:,l),a2(:,:,:,l),type//' nonblock halo=1 scalar'//text)
    enddo

    call mpp_clear_group_update(group_update1)
    call mpp_clear_group_update(group_update2)
    deallocate(a1,a2,global,global_all,base)

    !--- CGRID vector update -------------------------
    id1 = mpp_clock_id( type//' halo='//halostr//' CGRID', flags=MPP_CLOCK_SYNC )
    id2 = mpp_clock_id( type//' halo=1 CGRID', flags=MPP_CLOCK_SYNC )

    allocate( x1(isd:ied+shift,jsd:jed,       nz, num_fields) )
    allocate( y1(isd:ied,      jsd:jed+shift, nz, num_fields) )
    allocate( x2(isd:ied+shift,jsd:jed,       nz, num_fields) )
    allocate( y2(isd:ied,      jsd:jed+shift, nz, num_fields) )
    allocate(base1(isc:iec+shift, jsc:jec, nz, num_fields))
    allocate(base2(isc:iec, jsc:jec+shift, nz, num_fields))
    allocate(global1_all(1:nx+shift,1:ny,nz,ntiles) )
    allocate(global2_all(1:nx,1:ny+shift,nz,ntiles) )
    allocate(global1(1-whalo:nx+ehalo+shift, 1-shalo:ny+nhalo, nz, num_fields))
    allocate(global2(1-whalo:nx+ehalo, 1-shalo:ny+nhalo+shift, nz, num_fields))
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx+shift
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
          do j = 1, ny+shift
             do i = 1, nx
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    global1 = 0.0; global2 = 0.0
    do l = 1, num_fields
       global1(1:nx+shift,1:ny,:,l) = global1_all(1:nx+shift,1:ny,:,mytile)
       global2(1:nx,1:ny+shift,:,l) = global2_all(1:nx,1:ny+shift,:,mytile)
    end do

   if(folded_north) then
      do l = 1, num_fields
         call fill_folded_north_halo(global1(:,:,:,l), 1, 0, shift, 0, -1)
         call fill_folded_north_halo(global2(:,:,:,l), 0, 1, 0, shift, -1)
      enddo
   endif

    base1(isc:iec+shift,jsc:jec,:,:) = global1(isc:iec+shift,jsc:jec,:,:)
    base2(isc:iec,jsc:jec+shift,:,:) = global2(isc:iec,jsc:jec+shift,:,:)

    if(folded_north) then
       !redundant points must be equal and opposite for tripolar grid
        global2(nx/2+1:nx,     ny+shift,:,:) = -global2(nx/2:1:-1, ny+shift,:,:)
        global2(1-whalo:0,     ny+shift,:,:) = -global2(nx-whalo+1:nx, ny+shift,:,:)
        global2(nx+1:nx+ehalo, ny+shift,:,:) = -global2(1:ehalo,       ny+shift,:,:)
    else if(cubic_grid) then
       do l = 1, num_fields
          call fill_cubic_grid_halo(global1(:,:,:,l), global1_all, global2_all, mytile, 1, 0, 1, -1 )
          call fill_cubic_grid_halo(global2(:,:,:,l), global2_all, global1_all, mytile, 0, 1, -1, 1 )
       enddo
    endif

    x1 = 0; y1 = 0
    do l =1, num_fields
       call mpp_create_group_update(group_update1, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE)
    end do
    do l =1, num_fields
       call mpp_create_group_update(group_update2, x1(:,:,:,l), y1(:,:,:,l), domain, gridtype=CGRID_NE, &
            whalo=1, ehalo=1, shalo=1, nhalo=1 )
    end do

    x2(:,:,:,:) = global1(isd:ied+shift,jsd:jed,:,:)
    y2(:,:,:,:) = global2(isd:ied,jsd:jed+shift,:,:)

    do n = 1, num_iter
       x1 = 0.0; y1 = 0.0
       x1(isc:iec+shift,jsc:jec,      :,:) = base1(isc:iec+shift,jsc:jec,      :,:)
       y1(isc:iec,      jsc:jec+shift,:,:) = base2(isc:iec,      jsc:jec+shift,:,:)
       call mpp_clock_begin(id1)
       call mpp_do_group_update(group_update1, domain, x1(isc,jsc,1,1))
       call mpp_clock_end(id1)
       if(n==num_iter) then
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums(x1(:,:,:,l),x2(:,:,:,l),type//' halo='//halostr//' CGRID X'//text)
             call compare_checksums(y1(:,:,:,l),y2(:,:,:,l),type//' halo='//halostr//' CGRID Y'//text)
          enddo
       endif
    enddo

    !--- make sure non-blocking call is OK
    x1 = 0.0; y1 = 0.0
    x1(isc:iec+shift,jsc:jec,      :,:) = base1(isc:iec+shift,jsc:jec,      :,:)
    y1(isc:iec,      jsc:jec+shift,:,:) = base2(isc:iec,      jsc:jec+shift,:,:)
    call mpp_start_group_update(group_update1, domain, x1(isc,jsc,1,1))
    call mpp_complete_group_update(group_update1, domain, x1(isc,jsc,1,1))
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(x1(:,:,:,l),x2(:,:,:,l),type//' nonblock halo='//halostr//' CGRID X'//text)
       call compare_checksums(y1(:,:,:,l),y2(:,:,:,l),type//' nonblock halo='//halostr//' CGRID Y'//text)
    enddo

    x2 = 0; y2 = 0
    x2(isc-1:iec+1+shift,jsc-1:jec+1,:,:) = global1(isc-1:iec+1+shift,jsc-1:jec+1,:,:)
    y2(isc-1:iec+1,jsc-1:jec+1+shift,:,:) = global2(isc-1:iec+1,jsc-1:jec+1+shift,:,:)

    do n = 1, num_iter
       x1 = 0.0; y1 = 0.0
       x1(isc:iec+shift,jsc:jec,      :,:) = base1(isc:iec+shift,jsc:jec,      :,:)
       y1(isc:iec,      jsc:jec+shift,:,:) = base2(isc:iec,      jsc:jec+shift,:,:)
       call mpp_clock_begin(id2)
       call mpp_do_group_update(group_update2, domain, x1(isc,jsc,1,1))
       call mpp_clock_end(id2)
       if(n==num_iter) then
          do l = 1, num_fields
             write(text, '(i3.3)') l
             call compare_checksums(x1(:,:,:,l),x2(:,:,:,l),type//' halo=1 CGRID X'//text)
             call compare_checksums(y1(:,:,:,l),y2(:,:,:,l),type//' halo=1 CGRID Y'//text)
          enddo
       endif
    enddo

    x1 = 0.0; y1 = 0.0
    x1(isc:iec+shift,jsc:jec,      :,:) = base1(isc:iec+shift,jsc:jec,      :,:)
    y1(isc:iec,      jsc:jec+shift,:,:) = base2(isc:iec,      jsc:jec+shift,:,:)
    call mpp_start_group_update(group_update2, domain, x1(isc,jsc,1,1))
    call mpp_complete_group_update(group_update2, domain, x1(isc,jsc,1,1))
    do l = 1, num_fields
       write(text, '(i3.3)') l
       call compare_checksums(x1(:,:,:,l),x2(:,:,:,l),type//' nonblock halo=1 CGRID X'//text)
       call compare_checksums(y1(:,:,:,l),y2(:,:,:,l),type//' nonblock halo=1 CGRID Y'//text)
    enddo

    call mpp_clear_group_update(group_update1)
    call mpp_clear_group_update(group_update2)

    deallocate(x1, y1, global1, global2)
    deallocate(x2, y2, global1_all, global2_all)
    deallocate(base1, base2)
    call mpp_deallocate_domain(domain)

end subroutine test_halosize_update

  !###############################################################
  subroutine test_unstruct_update( type )
    character(len=*), intent(in) :: type

    type(domain2D) :: SG_domain
    type(domainUG) :: UG_domain
    integer        :: num_contact, ntiles, npes_per_tile
    integer        :: i, j, k, l, n, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem, lsg, leg

    integer, allocatable, dimension(:)       :: pe_start, pe_end, npts_tile, grid_index, ntiles_grid
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:)     :: x1, x2, g1, g2
    real,    allocatable, dimension(:,:,:)   :: a1, a2, gdata
    real,    allocatable, dimension(:,:)     :: rmask
    real,    allocatable, dimension(:)       :: frac_crit
    logical, allocatable, dimension(:,:,:)   :: lmask
    integer, allocatable, dimension(:)       :: isl, iel, jsl, jel
    logical            :: cubic_grid
    character(len=3)   :: text
    integer            :: nx_save, ny_save, tile
    integer            :: ntotal_land, istart, iend, pos

    cubic_grid         = .false.

    nx_save = nx
    ny_save = ny
    !--- check the type
    select case(type)
    case ( 'Cubic-Grid' )
       if( nx_cubic == 0 ) then
          call mpp_error(NOTE,'test_unstruct_update: for Cubic_grid mosaic, nx_cubic is zero, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       if( nx_cubic .NE. ny_cubic ) then
          call mpp_error(NOTE,'test_unstruct_update: for Cubic_grid mosaic, nx_cubic does not equal ny_cubic, '//&
               'No test is done for Cubic-Grid mosaic. ' )
          return
       endif
       nx = nx_cubic
       ny = ny_cubic
       ntiles = 6
       num_contact = 12
       cubic_grid = .true.
       if( mod(npes, ntiles) == 0 ) then
          npes_per_tile = npes/ntiles
          write(outunit,*)'NOTE from test_unstruct_update ==> For Mosaic "', trim(type), &
               '", each tile will be distributed over ', npes_per_tile, ' processors.'
       else
          call mpp_error(NOTE,'test_unstruct_update: npes should be multiple of ntiles No test is done for '//trim(type))
          return
       endif
       if(layout_cubic(1)*layout_cubic(2) == npes_per_tile) then
          layout = layout_cubic
       else
          call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       endif
       allocate(frac_crit(ntiles))
       frac_crit(1) = 0.3; frac_crit(2) = 0.1; frac_crit(3) = 0.6
       frac_crit(4) = 0.2; frac_crit(5) = 0.4; frac_crit(6) = 0.5

    case default
       call mpp_error(FATAL, 'test_group_update: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    do n = 1, ntiles
       pe_start(n) = (n-1)*npes_per_tile
       pe_end(n)   = n*npes_per_tile-1
    end do

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    !--- define domain
    if( cubic_grid ) then
       call define_cubic_mosaic(type, SG_domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end )
    endif

    !--- setup data
    call mpp_get_compute_domain( SG_domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( SG_domain, isd, ied, jsd, jed )

    allocate(lmask(nx,ny,ntiles))
    allocate(npts_tile(ntiles))
    lmask = .false.
    if(mpp_pe() == mpp_root_pe() ) then
       allocate(rmask(nx,ny))
       !--- construct gmask.
       do n = 1, ntiles
          call random_number(rmask)
          do j = 1, ny
             do i = 1, nx
                if(rmask(i,j) > frac_crit(n)) then
                   lmask(i,j,n) = .true.
                endif
             enddo
          enddo
          npts_tile(n) = count(lmask(:,:,n))
       enddo
       ntotal_land = sum(npts_tile)
       allocate(grid_index(ntotal_land))
       l = 0
       allocate(isl(0:mpp_npes()-1), iel(0:mpp_npes()-1))
       allocate(jsl(0:mpp_npes()-1), jel(0:mpp_npes()-1))
       call mpp_get_compute_domains(SG_domain,xbegin=isl,xend=iel,ybegin=jsl,yend=jel)

       do n = 1, ntiles
          do j = 1, ny
             do i = 1, nx
                if(lmask(i,j,n)) then
                   l = l + 1
                   grid_index(l) = (j-1)*nx+i
                endif
             enddo
          enddo
       enddo
       deallocate(rmask, isl, iel, jsl, jel)
    endif
    call mpp_broadcast(npts_tile, ntiles, mpp_root_pe())
    if(mpp_pe() .NE. mpp_root_pe()) then
       ntotal_land = sum(npts_tile)
       allocate(grid_index(ntotal_land))
    endif
    call mpp_broadcast(grid_index, ntotal_land, mpp_root_pe())

    allocate(ntiles_grid(ntotal_land))
    ntiles_grid = 1
   !--- define the unstructured grid domain
    call mpp_define_unstruct_domain(UG_domain, SG_domain, npts_tile, ntiles_grid, mpp_npes(), 1, grid_index, name="LAND unstruct")
    call mpp_get_UG_compute_domain(UG_domain, istart, iend)

    !--- figure out lmask according to grid_index
    pos = 0
    do n = 1, ntiles
       do l = 1, npts_tile(n)
          pos = pos + 1
          j = (grid_index(pos)-1)/nx + 1
          i = mod((grid_index(pos)-1),nx) + 1
          lmask(i,j,n) = .true.
       enddo
    enddo

    !--- set up data
    allocate(gdata(nx,ny,ntiles))
    gdata = -999
    do n = 1, ntiles
       do j = 1, ny
          do i = 1, nx
             if(lmask(i,j,n)) then
                gdata(i,j,n) = n*1.e+3 + i + j*1.e-3
             endif
          end do
       end do
    end do

    !--- test the 2-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,1), a2(isc:iec,jsc:jec,1 ) )

    tile = mpp_pe()/npes_per_tile + 1
    do j = jsc, jec
       do i = isc, iec
          a1(i,j,1) = gdata(i,j,tile)
       enddo
    enddo
    a2 = -999

    allocate(x1(istart:iend,1), x2(istart:iend,1))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       x2(l,1) = gdata(i,j,tile)
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1(:,:,1), x1(:,1))
    call compare_checksums_2D(x1, x2, type//' SG2UG 2-D compute domain')
    call mpp_pass_UG_to_SG(UG_domain, x1(:,1), a2(:,:,1))

    call compare_checksums(a1(:,:,1:1),a2(:,:,1:1),type//' UG2SG 2-D compute domain')
    deallocate(a1,a2,x1,x2)

    !--- test the 3-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,nz), a2(isc:iec,jsc:jec,nz ) )

    tile = mpp_pe()/npes_per_tile + 1
    do k = 1, nz
       do j = jsc, jec
          do i = isc, iec
             a1(i,j,k) = gdata(i,j,tile)
             if(a1(i,j,k) .NE. -999) a1(i,j,k) = a1(i,j,k) + k*1.e-6
          enddo
       enddo
    enddo
    a2 = -999

    allocate(x1(istart:iend,nz), x2(istart:iend,nz))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       do k = 1, nz
          x2(l,k) = gdata(i,j,tile) + k*1.e-6
       enddo
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1, x1)
    call compare_checksums_2D(x1, x2, type//' SG2UG 3-D compute domain')
    call mpp_pass_UG_to_SG(UG_domain, x1, a2)

    call compare_checksums(a1,a2,type//' UG2SG 3-D compute domain')
    deallocate(a1,a2,x1,x2)

    !--- test the 2-D data is on data domain
    allocate( a1(isd:ied, jsd:jed,1), a2(isd:ied,jsd:jed,1 ) )
    a1 = -999; a2 = -999

    tile = mpp_pe()/npes_per_tile + 1
    do j = jsc, jec
       do i = isc, iec
          a1(i,j,1) = gdata(i,j,tile)
       enddo
    enddo
    a2 = -999

    allocate(x1(istart:iend,1), x2(istart:iend,1))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       x2(l,1) = gdata(i,j,tile)
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1(:,:,1), x1(:,1))
    call compare_checksums_2D(x1, x2, type//' SG2UG 2-D data domain')
    call mpp_pass_UG_to_SG(UG_domain, x1(:,1), a2(:,:,1))

    call compare_checksums(a1(:,:,1:1),a2(:,:,1:1),type//' UG2SG 2-D data domain')
    deallocate(a1,a2,x1,x2)

    !--- test the 3-D data is on computing domain
    allocate( a1(isd:ied, jsd:jed,nz), a2(isd:ied,jsd:jed,nz ) )
    a1 = -999; a2 = -999

    tile = mpp_pe()/npes_per_tile + 1
    do k = 1, nz
       do j = jsc, jec
          do i = isc, iec
             a1(i,j,k) = gdata(i,j,tile)
             if(a1(i,j,k) .NE. -999) a1(i,j,k) = a1(i,j,k) + k*1.e-6
          enddo
       enddo
    enddo
    a2 = -999

    allocate(x1(istart:iend,nz), x2(istart:iend,nz))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    pos = 0
    do n = 1, tile-1
       pos = pos + npts_tile(n)
    enddo
    do l = istart, iend
       i = mod((grid_index(pos+l)-1), nx) + 1
       j = (grid_index(pos+l)-1)/nx + 1
       do k = 1, nz
          x2(l,k) = gdata(i,j,tile) + k*1.e-6
       enddo
    enddo

    call mpp_pass_SG_to_UG(UG_domain, a1, x1)
    call compare_checksums_2D(x1, x2, type//' SG2UG 3-D data domain')
    call mpp_pass_UG_to_SG(UG_domain, x1, a2)

    call compare_checksums(a1,a2,type//' UG2SG 3-D data domain')
    deallocate(a1,a2,x1,x2)

    !----------------------------------------------------------------
    !    test mpp_global_field_ug
    !----------------------------------------------------------------
    call mpp_get_UG_global_domain(UG_domain, lsg, leg)
    tile = mpp_get_UG_domain_tile_id(UG_domain)
    allocate(g1(lsg:leg,nz), g2(lsg:leg,nz), x1(istart:iend,nz))
    g1 = 0
    g2 = 0
    x1 = 0
    do k = 1, nz
       do l = lsg, leg
          g1(l,k) = tile*1e6 + l + k*1.e-3
       enddo
       do l = istart, iend
          x1(l,k) = g1(l,k)
       enddo
    enddo

    call mpp_global_field_ug(UG_domain, x1, g2)
    call compare_checksums_2D(g1,g2,type//' global_field_ug 3-D')

    g2 = 0.0
    call mpp_global_field_ug(UG_domain, x1(:,1), g2(:,1))
    call compare_checksums_2D(g1(:,1:1),g2(:,1:1),type//' global_field_ug 2-D')

    deallocate(g1,g2,x1)

  end subroutine test_unstruct_update



  !#################################################################################

  subroutine fill_halo_zero(data, whalo, ehalo, shalo, nhalo, xshift, yshift, isc, iec, jsc, jec, isd, ied, jsd, jed)
    integer,                         intent(in) :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer,                         intent(in) :: whalo, ehalo, shalo, nhalo, xshift, yshift
    real, dimension(isd:,jsd:,:), intent(inout) :: data

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

  end subroutine fill_halo_zero

  !##############################################################################
  ! this routine fill the halo points for the regular mosaic.
  subroutine fill_regular_mosaic_halo(data, data_all, te, tse, ts, tsw, tw, tnw, tn, tne)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data_all
    integer,                              intent(in)    :: te, tse, ts, tsw, tw, tnw, tn, tne

       data(nx+1:nx+ehalo, 1:ny,          :) = data_all(1:ehalo,       1:ny,          :, te) ! east
       data(1:nx,          1-shalo:0,     :) = data_all(1:nx,          ny-shalo+1:ny, :, ts) ! south
       data(1-whalo:0,     1:ny,          :) = data_all(nx-whalo+1:nx, 1:ny,          :, tw) ! west
       data(1:nx,          ny+1:ny+nhalo, :) = data_all(1:nx,          1:nhalo,       :, tn) ! north
       data(nx+1:nx+ehalo, 1-shalo:0,     :) = data_all(1:ehalo,       ny-shalo+1:ny, :,tse) ! southeast
       data(1-whalo:0,     1-shalo:0,     :) = data_all(nx-whalo+1:nx, ny-shalo+1:ny, :,tsw) ! southwest
       data(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = data_all(1:ehalo,       1:nhalo,       :,tnw) ! northeast
       data(1-whalo:0,     ny+1:ny+nhalo, :) = data_all(nx-whalo+1:nx, 1:nhalo,       :,tne) ! northwest



  end subroutine fill_regular_mosaic_halo

  !################################################################################
  subroutine fill_folded_north_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff

    data(1-whalo:0,                  1:nyp,:) =      data(nx-whalo+1:nx,        1:ny+jshift,:) ! west
    data(nx+1:nx+ehalo+ishift,       1:nyp,:) =      data(1:ehalo+ishift,       1:ny+jshift,:) ! east
    if(m1 .GE. 1-whalo) data(1-whalo:m1,  nyp+1:nyp+nhalo,:) = sign*data(whalo+m2:1+ishift:-1, nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(m1+1:nx+m2,       nyp+1:nyp+nhalo,:) = sign*data(nx+ishift:1:-1,       nyp-joff:nyp-nhalo-joff+1:-1,:)
    data(nx+m2+1:nxp+ehalo,nyp+1:nyp+nhalo,:) = sign*data(nx:nx-ehalo+m1+1:-1,  nyp-joff:nyp-nhalo-joff+1:-1,:)

  end subroutine fill_folded_north_halo

  !################################################################################
  subroutine fill_folded_south_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = ishift - ioff
    m2 = 2*ishift - ioff


    data(1-whalo:0,                  1:nyp,:) =      data(nx-whalo+1:nx,        1:nyp,:) ! west
    data(nx+1:nx+ehalo+ishift,       1:nyp,:) =      data(1:ehalo+ishift,       1:nyp,:) ! east
    if(m1 .GE. 1-whalo)data(1-whalo:m1, 1-shalo:0,:) = sign*data(whalo+m2:1+ishift:-1, shalo+jshift:1+jshift:-1,:)
    data(m1+1:nx+m2,       1-shalo:0,:) = sign*data(nxp:1:-1,             shalo+jshift:1+jshift:-1,:)
    data(nx+m2+1:nxp+ehalo,1-shalo:0,:) = sign*data(nx:nx-ehalo+m1+1:-1,  shalo+jshift:1+jshift:-1,:)

  end subroutine fill_folded_south_halo

  !################################################################################
  subroutine fill_folded_west_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)      = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) data(1-whalo:0, 1-shalo:m1, :) = sign*data(whalo+ishift:1+ishift:-1, shalo+m2:1+jshift:-1,:)
    data(1-whalo:0, m1+1:ny+m2, :) = sign*data(whalo+ishift:1+ishift:-1, nyp:1:-1, :)
    data(1-whalo:0, ny+m2+1:nyp+nhalo,:) = sign*data(whalo+ishift:1+ishift:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_west_halo

  !################################################################################
  subroutine fill_folded_east_halo(data, ioff, joff, ishift, jshift, sign)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    integer,                              intent(in   ) :: ioff, joff, ishift, jshift, sign
    integer  :: nxp, nyp, m1, m2

    nxp = nx+ishift
    nyp = ny+jshift
    m1 = jshift - joff
    m2 = 2*jshift - joff

    data(1:nxp, 1-shalo:0, :)      = data(1:nxp, ny-shalo+1:ny, :) ! south
    data(1:nxp, ny+1:nyp+nhalo, :) = data(1:nxp, 1:nhalo+jshift,:) ! north
    if(m1 .GE. 1-shalo) data(nxp+1:nxp+ehalo, 1-shalo:m1, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, shalo+m2:1+jshift:-1,:)
    data(nxp+1:nxp+ehalo, m1+1:ny+m2, :) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, nyp:1:-1, :)
    data(nxp+1:nxp+ehalo, ny+m2+1:nyp+nhalo,:) = sign*data(nxp-ioff:nxp-ehalo-ioff+1:-1, ny:ny-nhalo+m1+1:-1,:)

  end subroutine fill_folded_east_halo

  !################################################################################
  subroutine fill_four_tile_bound(data_all, is, ie, js, je, ioff, joff, tile, &
                                   ebound, sbound, wbound, nbound )
    real, dimension(:,:,:,:),       intent(in)    :: data_all
    integer,                        intent(in)    :: is, ie, js, je
    integer,                        intent(in)    :: tile, ioff, joff
    real, dimension(:,:), optional, intent(inout) :: ebound, sbound, wbound, nbound
    integer                                       :: tw, te, ts, tn

    if(tile == 1 .OR. tile == 3) te = tile + 1
    if(tile == 2 .OR. tile == 4) te = tile - 1
    if(tile == 1 .OR. tile == 2) ts = tile + 2
    if(tile == 3 .OR. tile == 4) ts = tile - 2
    tw = te;   tn = ts
    if(present(ebound)) then
       if( ie == nx ) then
          ebound(:,:) = data_all(1, js:je+joff, :, te)
       else
          ebound(:,:) = data_all(ie+ioff, js:je+joff, :, tile)
       end if
    end if

    if(present(wbound)) then
       if( is == 1 ) then
          wbound(:,:) = data_all(nx+ioff, js:je+joff, :, tw)
       else
          wbound(:,:) = data_all(is, js:je+joff, :, tile)
       end if
    end if

    if(present(sbound)) then
       if( js == 1 ) then
          sbound(:,:) = data_all(is:ie+ioff, ny+joff, :, ts)
       else
          sbound(:,:) = data_all(is:ie+ioff, js, :, tile)
       end if
    end if

    if(present(nbound)) then
       if( je == ny ) then
          nbound(:,:) = data_all(is:ie+ioff, 1, :, tn)
       else
          nbound(:,:) = data_all(is:ie+ioff, je+joff, :, tile)
       end if
    end if

    return

  end subroutine fill_four_tile_bound


  !################################################################################
  subroutine fill_torus_bound(data_all, is, ie, js, je, ioff, joff, tile, &
                                     sbound, wbound)
    real, dimension(:,:,:),       intent(in)    :: data_all
    integer,                        intent(in)    :: is, ie, js, je
    integer,                        intent(in)    :: tile, ioff, joff
    real, dimension(:,:), optional, intent(inout) :: sbound, wbound
    integer                                       :: tw, te, ts, tn
    integer                                       :: js1, js2, is1, is2

    if(tile .NE. 1) call mpp_error(FATAL, "fill_torus_bound: tile must be 1")

    js2 = js
    js1 = 1
    if( js == 1 .AND. joff==1 )  then
      js1 = 2
      js2 = js+1
    endif
    is2 = is
    is1 = 1
    if( is == 1 .AND. ioff==1 )  then
      is1 = 2
      is2 = is+1
    endif

    if(present(wbound)) then
       if(ioff .NE. 1) call mpp_error(FATAL, "fill_torus_bound: ioff must be 1 when wbound present")
       if( is == 1 ) then
          wbound(js1:,:) = data_all(nx+ioff, js2:je+joff, :)
       else
          wbound(js1:,:) = data_all(is, js2:je+joff, :)
       end if
       if(js1 == 2) then
          if( is == 1 ) then
             wbound(1,:) = data_all(nx+1, ny+1, :)
          else
             wbound(1,:) = data_all(is, ny+1, :)
          endif
       endif
    end if

    if(present(sbound)) then
       if(joff .NE. 1) call mpp_error(FATAL, "fill_torus_bound: joff must be 1 when sbound present")
       if( js == 1 ) then
          sbound(is1:,:) = data_all(is2:ie+ioff, ny+joff, :)
       else
          sbound(is1:,:) = data_all(is2:ie+ioff, js, :)
       end if
       if(is1 == 2) then
          if( js == 1 ) then
             sbound(1,:) = data_all(nx+1, ny+1, :)
          else
             sbound(1,:) = data_all(nx+1, js, :)
          endif
       endif
    end if

    return

  end subroutine fill_torus_bound

  !################################################################################
  subroutine fill_folded_north_bound(data_all, is, ie, js, je, ioff, joff, tile, &
                                     sbound, wbound)
    real, dimension(:,:,:),       intent(in)    :: data_all
    integer,                        intent(in)    :: is, ie, js, je
    integer,                        intent(in)    :: tile, ioff, joff
    real, dimension(:,:), optional, intent(inout) :: sbound, wbound
    integer                                       :: tw, te, ts, tn
    integer                                       :: js1, js2

    if(tile .NE. 1) call mpp_error(FATAL, "fill_folded_north_bound: tile must be 1")

    js2 = js
    js1 = 1
    if( js == 1 .AND. joff==1 )  then
      js1 = 2
      js2 = js+1
    endif

    if(present(wbound)) then
       if( is == 1 ) then
          wbound(js1:,:) = data_all(nx+ioff, js2:je+joff, :)
       else
          wbound(js1:,:) = data_all(is, js2:je+joff, :)
       end if
    end if

    if(present(sbound)) then
       if( js == 1 ) then
          sbound(:,:) = 0
       else
          if( is == 1 .AND. ioff == 1 ) then
             sbound(1,:)  = data_all(nx+1, js, :)
             sbound(2:,:) = data_all(is+1:ie+ioff, js, :)
          else
             sbound(:,:) = data_all(is:ie+ioff, js, :)
          endif
       end if
    end if

    return

  end subroutine fill_folded_north_bound

  !################################################################################
  subroutine fill_cubic_grid_bound(data1_all, data2_all, is, ie, js, je, ioff, joff, tile, sign1, sign2, &
                                   ebound, sbound, wbound, nbound )
    real, dimension(:,:,:,:),       intent(in)    :: data1_all, data2_all
    integer,                        intent(in)    :: is, ie, js, je
    integer,                        intent(in)    :: tile, ioff, joff, sign1, sign2
    real, dimension(:,:), optional, intent(inout) :: ebound, sbound, wbound, nbound
    integer                                       :: tw, te, ts, tn

    if(mod(tile,2) == 0) then ! tile 2, 4, 6
       tw = tile - 1; te = tile + 2; ts = tile - 2; tn = tile + 1
       if(te > 6 ) te = te - 6
       if(ts < 1 ) ts = ts + 6
       if(tn > 6 ) tn = tn - 6
       !--- East bound
       if(present(ebound)) then
          if(ie == nx) then
             ebound(:,:) = sign1*data2_all(nx+joff-js+1:nx-je+1:-1,1,:,te)
          else
             ebound(:,:) = data1_all(ie+ioff, js:je+joff, :,tile)
          end if
       end if
       !--- South bound
       if(present(sbound)) then
          if(js == 1) then
             sbound(:,:) = sign2*data2_all(nx+joff, ny+ioff-is+1:ny-ie+1:-1,:,ts)
          else
             sbound(:,:) = data1_all(is:ie+ioff, js, :,tile)
          end if
       end if

       !--- West bound
       if(present(wbound)) then
          if(is == 1) then
             wbound(:,:) = data1_all(nx+ioff, js:je+joff,:,tw)
          else
             wbound(:,:) = data1_all(is, js:je+joff,:,tile)
          end if
       end if

       !--- north bound
       if(present(nbound)) then
          if(je == ny) then
             nbound(:,:) = data1_all(is:ie+ioff, 1,:,tn)
          else
             nbound(:,:) = data1_all(is:ie+ioff, je+joff, :,tile)
          end if
       end if
    else ! tile 1, 3, 5
       tw = tile - 2; te = tile + 1; ts = tile - 1; tn = tile + 2
       if(tw < 1 ) tw = tw + 6
       if(ts < 1 ) ts = ts + 6
       if(tn > 6 ) tn = tn - 6
       !--- East bound
       if(present(ebound)) then
          if(ie == nx) then
             ebound(:,:) = data1_all(1, js:je+joff, :,te)
          else
             ebound(:,:) = data1_all(ie+ioff, js:je+joff, :,tile)
          end if
       end if
       !--- South bound
       if(present(sbound)) then
          if(js == 1) then
             sbound(:,:) = data1_all(is:ie+ioff,ny+joff,:,ts)
          else
             sbound(:,:) = data1_all(is:ie+ioff, js, :,tile)
          end if
       end if

       !--- West bound
       if(present(wbound)) then
          if(is == 1) then
             wbound(:,:) = sign1*data2_all(nx+joff-js+1:nx-je+1:-1,ny+ioff,:,tw)
          else
             wbound(:,:) = data1_all(is, js:je+joff,:,tile)
          end if
       end if

       !--- north bound
       if(present(nbound)) then
          if(je == ny) then
             nbound(:,:) = sign2*data2_all(1, ny+ioff-is+1:ny-ie+1:-1,:,tn)
          else
             nbound(:,:) = data1_all(is:ie+ioff, je+joff, :,tile)
          end if
       end if

    end if

  end subroutine fill_cubic_grid_bound

  !##############################################################################
  ! this routine fill the halo points for the cubic grid. ioff and joff is used to distinguish
  ! T, C, E, or N-cell
  subroutine fill_cubic_grid_halo(data, data1_all, data2_all, tile, ioff, joff, sign1, sign2)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data1_all, data2_all
    integer,                              intent(in)    :: tile, ioff, joff, sign1, sign2
    integer                                             :: lw, le, ls, ln

    if(mod(tile,2) == 0) then ! tile 2, 4, 6
       lw = tile - 1; le = tile + 2; ls = tile - 2; ln = tile + 1
       if(le > 6 ) le = le - 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       data(1-whalo:0, 1:ny+joff, :) = data1_all(nx-whalo+1:nx, 1:ny+joff, :, lw) ! west
       do i = 1, ehalo
          data(nx+i+ioff, 1:ny+joff, :)    = sign1*data2_all(nx+joff:1:-1, i+ioff, :, le) ! east
       end do
       do i = 1, shalo
          data(1:nx+ioff, 1-i, :)     = sign2*data2_all(nx-i+1, ny+ioff:1:-1, :, ls) ! south
       end do
       data(1:nx+ioff, ny+1+joff:ny+nhalo+joff, :) = data1_all(1:nx+ioff, 1+joff:nhalo+joff, :, ln) ! north
    else ! tile 1, 3, 5
       lw = tile - 2; le = tile + 1; ls = tile - 1; ln = tile + 2
       if(lw < 1 ) lw = lw + 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       do i = 1, whalo
          data(1-i, 1:ny+joff, :)     = sign1*data2_all(nx+joff:1:-1, ny-i+1, :, lw) ! west
       end do
       data(nx+1+ioff:nx+ehalo+ioff, 1:ny+joff, :) = data1_all(1+ioff:ehalo+ioff, 1:ny+joff, :, le) ! east
       data(1:nx+ioff, 1-shalo:0, :)     = data1_all(1:nx+ioff, ny-shalo+1:ny, :, ls) ! south
       do i = 1, nhalo
          data(1:nx+ioff, ny+i+joff, :)    = sign2*data2_all(i+joff, ny+ioff:1:-1, :, ln) ! north
       end do
    end if

  end subroutine fill_cubic_grid_halo

   !#####################################################################
  subroutine test_nonuniform_mosaic( type )
    character(len=*), intent(in) :: type

    type(domain2D)               :: domain
    integer                      :: num_contact, ntiles, ntile_per_pe
    integer                      :: i, j, k, n, nxm, nym, ni, nj, shift
    integer                      :: ism, iem, jsm, jem, isc, iec, jsc, jec
    integer                      :: isd, ied, jsd, jed
    integer                      :: indices(4), msize(2)
    character(len=128)           :: type2

    integer, allocatable, dimension(:)       :: tile
    integer, allocatable, dimension(:)       :: pe_start, pe_end, tile1, tile2
    integer, allocatable, dimension(:)       :: istart1, iend1, jstart1, jend1
    integer, allocatable, dimension(:)       :: istart2, iend2, jstart2, jend2
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:,:,:) :: global1_all, global2_all
    real,    allocatable, dimension(:,:,:,:) :: global1, global2, x, y

    shift = 0
    select case(type)
    case('Five-Tile') ! one tile will run on pe 0 and other four tiles will run on pe 1
       shift = 1      ! one extra point for symmetry domain
       ntiles = 5     ! tile 1 with resolution 2*nx and 2*ny and the tiles are nx and ny.
       num_contact = 11
       if(npes .NE. 2) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: Five-Tile mosaic will not be tested because npes is not 2')
          return
       end if
       nxm = 2*nx; nym = 2*ny
       layout = 1
       if( pe == 0) then
          ntile_per_pe = 1
          allocate(tile(ntile_per_pe))
          tile = 1
          indices = (/1,2*nx,1,2*ny/)
          ni = 2*nx; nj = 2*ny
       else
          ntile_per_pe = 4
          allocate(tile(ntile_per_pe))
          do n = 1, ntile_per_pe
             tile(n) = n + 1
          end do
          indices = (/1,nx,1,ny/)
          ni = nx; nj = ny
       end if
       allocate(pe_start(ntiles), pe_end(ntiles) )
       pe_start(1) = 0; pe_start(2:) = 1
       pe_end = pe_start
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles) )

    do n = 1, ntiles
       if(n==1) then
          global_indices(:,n) = (/1,2*nx,1,2*ny/)
       else
          global_indices(:,n) = (/1,nx,1,ny/)
       endif
!       global_indices(:,n) = indices
       layout2D(:,n)       = layout
    end do

    allocate(tile1(num_contact), tile2(num_contact) )
    allocate(istart1(num_contact), iend1(num_contact), jstart1(num_contact), jend1(num_contact) )
    allocate(istart2(num_contact), iend2(num_contact), jstart2(num_contact), jend2(num_contact) )

    !--- define domain
    select case(type)
    case( 'Five-Tile' )
       !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
       tile1(1) = 1; tile2(1) = 2
       istart1(1) = 2*nx; iend1(1) = 2*nx; jstart1(1) = 1;  jend1(1) = ny
       istart2(1) = 1;    iend2(1) = 1;    jstart2(1) = 1;  jend2(1) = ny
       !--- Contact line 2, between tile 1 (EAST) and tile 4 (WEST)
       tile1(2) = 1; tile2(2) = 4
       istart1(2) = 2*nx; iend1(2) = 2*nx; jstart1(2) = ny+1; jend1(2) = 2*ny
       istart2(2) = 1;    iend2(2) = 1;    jstart2(2) = 1;    jend2(2) = ny
       !--- Contact line 3, between tile 1 (SOUTH) and tile 1 (NORTH)
       tile1(3) = 1; tile2(3) = 1
       istart1(3) = 1; iend1(3) = 2*nx; jstart1(3) = 1;    jend1(3) = 1
       istart2(3) = 1; iend2(3) = 2*nx; jstart2(3) = 2*ny; jend2(3) = 2*ny
       !--- Contact line 4, between tile 1 (WEST) and tile 3 (EAST)
       tile1(4) = 1; tile2(4) = 3
       istart1(4) = 1;  iend1(4) = 1;  jstart1(4) = 1;  jend1(4) = ny
       istart2(4) = nx; iend2(4) = nx; jstart2(4) = 1;  jend2(4) = ny
       !--- Contact line 5, between tile 1 (WEST) and tile 5 (EAST)
       tile1(5) = 1; tile2(5) = 5
       istart1(5) = 1;  iend1(5) = 1;  jstart1(5) = ny+1;  jend1(5) = 2*ny
       istart2(5) = nx; iend2(5) = nx; jstart2(5) = 1;     jend2(5) = ny
       !--- Contact line 6, between tile 2 (EAST) and tile 3 (WEST)
       tile1(6) = 2; tile2(6) = 3
       istart1(6) = nx; iend1(6) = nx; jstart1(6) = 1;  jend1(6) = ny
       istart2(6) = 1;  iend2(6) = 1;  jstart2(6) = 1;  jend2(6) = ny
       !--- Contact line 7, between tile 2 (SOUTH) and tile 4 (NORTH)  --- cyclic
       tile1(7) = 2; tile2(7) = 4
       istart1(7) = 1;  iend1(7) = nx; jstart1(7) = 1;   jend1(7) = 1
       istart2(7) = 1;  iend2(7) = nx; jstart2(7) = ny;  jend2(7) = ny
       !--- Contact line 8, between tile 2 (NORTH) and tile 4 (SOUTH)
       tile1(8) = 2; tile2(8) = 4
       istart1(8) = 1;  iend1(8) = nx; jstart1(8) = ny;  jend1(8) = ny
       istart2(8) = 1;  iend2(8) = nx; jstart2(8) = 1;   jend2(8) = 1
       !--- Contact line 9, between tile 3 (SOUTH) and tile 5 (NORTH)  --- cyclic
       tile1(9) = 3; tile2(9) = 5
       istart1(9) = 1;  iend1(9) = nx; jstart1(9) = 1;   jend1(9) = 1
       istart2(9) = 1;  iend2(9) = nx; jstart2(9) = ny;  jend2(9) = ny
       !--- Contact line 10, between tile 3 (NORTH) and tile 5 (SOUTH)
       tile1(10) = 3; tile2(10) = 5
       istart1(10) = 1;  iend1(10) = nx; jstart1(10) = ny;  jend1(10) = ny
       istart2(10) = 1;  iend2(10) = nx; jstart2(10) = 1;   jend2(10) = 1
       !--- Contact line 11, between tile 4 (EAST) and tile 5 (WEST)
       tile1(11) = 4; tile2(11) = 5
       istart1(11) = nx; iend1(11) = nx; jstart1(11) = 1;  jend1(11) = ny
       istart2(11) = 1;  iend2(11) = 1;  jstart2(11) = 1;  jend2(11) = ny
       msize(1) = 2*nx + whalo + ehalo
       msize(2) = 2*ny + shalo + nhalo
       call mpp_define_mosaic(global_indices, layout2D, domain, ntiles, num_contact, tile1, tile2, &
            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,      &
            pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,       &
            name = type, memory_size = msize, symmetry = .true.  )
    end select

    !--- setup data
    allocate(global1_all(1:nxm,1:nym,nz, ntiles) )
    allocate(global1(1-whalo:ni+ehalo,1-shalo:nj+nhalo,nz, ntile_per_pe) )
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, nym
             do i = 1, nxm
                global1_all(i,j,k,n) = n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global1(1:ni,1:nj,:,n) = global1_all(1:ni,1:nj,:,tile(n))
    end do

    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )

    allocate( x (ism:iem,jsm:jem,nz, ntile_per_pe) )
    x = 0.
    x(isc:iec,jsc:jec,:,:) = global1(isc:iec,jsc:jec,:,:)

    !--- fill up the value at halo points
    do n = 1, ntile_per_pe
       call fill_five_tile_halo(global1(:,:,:,n), global1_all, tile(n), 0, 0 )
    end do

    ! full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_update_domains( x(:,:,:,n), domain, tile_count = n )
    end do
    call mpp_clock_end(id)

   do n = 1, ntile_per_pe
      write(type2, *)type, " at tile_count = ",n
      call compare_checksums( x(isd:ied,jsd:jed,:,n), global1(isd:ied,jsd:jed,:,n), trim(type2) )
   end do

   deallocate(global1_all, global1, x)

    !------------------------------------------------------------------
    !  vector update : BGRID_NE, one extra point in each direction for Five-Tile
    !------------------------------------------------------------------
    !--- setup data
    allocate(global1_all(nxm+shift,nym+shift,nz, ntiles), global2_all(nxm+shift,nym+shift,nz, ntiles) )
    allocate(global1(1-whalo:ni+ehalo+shift,1-shalo:nj+nhalo+shift,nz, ntile_per_pe) )
    allocate(global2(1-whalo:ni+ehalo+shift,1-shalo:nj+nhalo+shift,nz, ntile_per_pe) )
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, nym+shift
             do i = 1, nxm+shift
                global1_all(i,j,k,n) = 1.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,n) = 2.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !------------------------------------------------------------------------
    ! --- make sure consisency on the boundary for Five-Tile mosaic
    ! --- east boundary will take the value of neighbor tile west,
    ! --- north boundary will take the value of neighbor tile south.
    !------------------------------------------------------------------------
    if(type == 'Five-Tile') then
       global1_all(nxm+1,    1:ny,:,1) = global1_all(1,    1:ny,:,2)  ! east
       global1_all(nxm+1,ny+1:nym,:,1) = global1_all(1,    1:ny,:,4)  ! east
       global1_all(1:nxm+1, nym+1,:,1) = global1_all(1:nxm+1, 1,:,1)  ! north
       global1_all(nx+1,     1:ny,:,2) = global1_all(1,    1:ny,:,3)  ! east
       global1_all(1:nx+1,   ny+1,:,2) = global1_all(1:nx+1,  1,:,4)  ! north
       global1_all(nx+1,     1:ny,:,3) = global1_all(1,    1:ny,:,1)  ! east
       global1_all(1:nx+1,   ny+1,:,3) = global1_all(1:nx+1,  1,:,5)  ! north
       global1_all(nx+1,     1:ny,:,4) = global1_all(1,    1:ny,:,5)  ! east
       global1_all(1:nx+1,   ny+1,:,4) = global1_all(1:nx+1,  1,:,2)  ! north
       global1_all(nx+1,     1:ny,:,5) = global1_all(1,ny+1:nym,:,1)  ! east
       global1_all(1:nx+1,   ny+1,:,5) = global1_all(1:nx+1,  1,:,3)  ! north
       global1_all(nx+1,     ny+1,:,2) = global1_all(1,       1,:,5)  ! northeast
       global1_all(nx+1,     ny+1,:,3) = global1_all(1,    ny+1,:,1)  ! northeast
       global2_all(nxm+1,    1:ny,:,1) = global2_all(1,    1:ny,:,2)  ! east
       global2_all(nxm+1,ny+1:nym,:,1) = global2_all(1,    1:ny,:,4)  ! east
       global2_all(1:nxm+1, nym+1,:,1) = global2_all(1:nxm+1, 1,:,1)  ! north
       global2_all(nx+1,     1:ny,:,2) = global2_all(1,    1:ny,:,3)  ! east
       global2_all(1:nx+1,   ny+1,:,2) = global2_all(1:nx+1,  1,:,4)  ! north
       global2_all(nx+1,     1:ny,:,3) = global2_all(1,    1:ny,:,1)  ! east
       global2_all(1:nx+1,   ny+1,:,3) = global2_all(1:nx+1,  1,:,5)  ! north
       global2_all(nx+1,     1:ny,:,4) = global2_all(1,    1:ny,:,5)  ! east
       global2_all(1:nx+1,   ny+1,:,4) = global2_all(1:nx+1,  1,:,2)  ! north
       global2_all(nx+1,     1:ny,:,5) = global2_all(1,ny+1:nym,:,1)  ! east
       global2_all(1:nx+1,   ny+1,:,5) = global2_all(1:nx+1,  1,:,3)  ! north
       global2_all(nx+1,     ny+1,:,2) = global2_all(1,       1,:,5)  ! northeast
       global2_all(nx+1,     ny+1,:,3) = global2_all(1,    ny+1,:,1)  ! northeast
    end if

    do n = 1, ntile_per_pe
       global1(1:ni+shift,1:nj+shift,:,n) = global1_all(1:ni+shift,1:nj+shift,:,tile(n))
       global2(1:ni+shift,1:nj+shift,:,n) = global2_all(1:ni+shift,1:nj+shift,:,tile(n))
    end do

    allocate( x (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )
    allocate( y (ism:iem+shift,jsm:jem+shift,nz,ntile_per_pe) )

    x = 0.; y = 0
    x (isc:iec+shift,jsc:jec+shift,:,:) = global1(isc:iec+shift,jsc:jec+shift,:,:)
    y (isc:iec+shift,jsc:jec+shift,:,:) = global2(isc:iec+shift,jsc:jec+shift,:,:)

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points.
    !-----------------------------------------------------------------------
    do n = 1, ntile_per_pe
       call fill_five_tile_halo(global1(:,:,:,n), global1_all, tile(n), shift, shift)
       call fill_five_tile_halo(global2(:,:,:,n), global2_all, tile(n), shift, shift)
    end do

    id = mpp_clock_id( type//' BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_update_domains( x(:,:,:,n), y(:,:,:,n), domain, gridtype=BGRID_NE, tile_count = n )
    end do
    call mpp_clock_end(id)

   do n = 1, ntile_per_pe
      write(type2, *)type, " at tile_count = ",n
      call compare_checksums( x(isd:ied+shift,jsd:jed+shift,:,n), global1(isd:ied+shift,jsd:jed+shift,:,n), &
                              trim(type2)//' BGRID_NE X')
      call compare_checksums( y(isd:ied+shift,jsd:jed+shift,:,n), global2(isd:ied+shift,jsd:jed+shift,:,n), &
                              trim(type2)//' BGRID_NE Y')
   end do

   deallocate(global1_all, global2_all, global1, global2, x, y)

    !------------------------------------------------------------------
    !  vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- setup data
    allocate(global1_all(nxm+shift,nym,nz, ntiles), global2_all(nxm,nym+shift,nz, ntiles) )
    allocate(global1(1-whalo:ni+ehalo+shift, 1-shalo:nj+nhalo,       nz, ntile_per_pe) )
    allocate(global2(1-whalo:ni+ehalo,       1-shalo:nj+nhalo+shift, nz, ntile_per_pe) )
    do n = 1, ntiles
       do k = 1, nz
          do j = 1, nym
             do i = 1, nxm+shift
                global1_all(i,j,k,n) = 1.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
          do j = 1, nym+shift
             do i = 1, nxm
                global2_all(i,j,k,n) = 2.0e3 + n + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    !------------------------------------------------------------------------
    ! --- make sure consisency on the boundary for Five-Tile mosaic
    ! --- east boundary will take the value of neighbor tile west,
    ! --- north boundary will take the value of neighbor tile south.
    !------------------------------------------------------------------------
    if(type == 'Five-Tile') then
       global1_all(nxm+1,    1:ny,:,1) = global1_all(1,    1:ny,:,2)  ! east
       global1_all(nxm+1,ny+1:nym,:,1) = global1_all(1,    1:ny,:,4)  ! east
       global1_all(nx+1,     1:ny,:,2) = global1_all(1,    1:ny,:,3)  ! east
       global1_all(nx+1,     1:ny,:,3) = global1_all(1,    1:ny,:,1)  ! east
       global1_all(nx+1,     1:ny,:,4) = global1_all(1,    1:ny,:,5)  ! east
       global1_all(nx+1,     1:ny,:,5) = global1_all(1,ny+1:nym,:,1)  ! east
       global2_all(1:nxm,   nym+1,:,1) = global2_all(1:nxm,   1,:,1)  ! north
       global2_all(1:nx,     ny+1,:,2) = global2_all(1:nx,    1,:,4)  ! north
       global2_all(1:nx,     ny+1,:,3) = global2_all(1:nx,    1,:,5)  ! north
       global2_all(1:nx,     ny+1,:,4) = global2_all(1:nx,    1,:,2)  ! north
       global2_all(1:nx,     ny+1,:,5) = global2_all(1:nx,    1,:,3)  ! north
    end if

    do n = 1, ntile_per_pe
       global1(1:ni+shift,      1:nj,:,n) = global1_all(1:ni+shift,      1:nj,:,tile(n))
       global2(1:ni,      1:nj+shift,:,n) = global2_all(1:ni,      1:nj+shift,:,tile(n))
    end do

    allocate( x (ism:iem+shift,      jsm:jem,nz,ntile_per_pe) )
    allocate( y (ism:iem,      jsm:jem+shift,nz,ntile_per_pe) )

    x = 0.; y = 0
    x (isc:iec+shift,      jsc:jec,:,:) = global1(isc:iec+shift,      jsc:jec,:,:)
    y (isc:iec,      jsc:jec+shift,:,:) = global2(isc:iec,      jsc:jec+shift,:,:)

    !-----------------------------------------------------------------------
    !                   fill up the value at halo points.
    !-----------------------------------------------------------------------
    do n = 1, ntile_per_pe
       call fill_five_tile_halo(global1(:,:,:,n), global1_all, tile(n), shift, 0)
       call fill_five_tile_halo(global2(:,:,:,n), global2_all, tile(n), 0, shift)
    end do

    id = mpp_clock_id( type//' CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    do n = 1, ntile_per_pe
       call mpp_update_domains( x(:,:,:,n), y(:,:,:,n), domain, gridtype=CGRID_NE, tile_count = n )
    end do
    call mpp_clock_end(id)

   do n = 1, ntile_per_pe
      write(type2, *)type, " at tile_count = ",n
      call compare_checksums( x(isd:ied+shift,jsd:jed,:,n), global1(isd:ied+shift,jsd:jed,:,n), &
                              trim(type2)//' CGRID_NE X')
      call compare_checksums( y(isd:ied,jsd:jed+shift,:,n), global2(isd:ied,jsd:jed+shift,:,n), &
                              trim(type2)//' CGRID_NE Y')
   end do

   deallocate(global1_all, global2_all, global1, global2, x, y)

  end subroutine test_nonuniform_mosaic

  subroutine fill_five_tile_halo(data, data_all, tile, ioff, joff)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data_all
    integer,                              intent(in)    :: tile, ioff, joff
    integer                                             :: nxm, nym

    nxm = 2*nx; nym = 2*ny

    select case(tile)
    case(1)
       data(nxm+1+ioff:nxm+ehalo+ioff,                     1:ny,:) = data_all(1+ioff:ehalo+ioff,              1:ny,:,2) ! east
       data(nxm+1+ioff:nxm+ehalo+ioff,            ny+1:nym+joff,:) = data_all(1+ioff:ehalo+ioff,         1:ny+joff,:,4) ! east
       data(1-whalo:0,                                     1:ny,:) = data_all(nx-whalo+1:nx,                  1:ny,:,3) ! west
       data(1-whalo:0,                            ny+1:nym+joff,:) = data_all(nx-whalo+1:nx,             1:ny+joff,:,5) ! west
       data(1:nxm+ioff,                               1-shalo:0,:) = data_all(1:nxm+ioff,          nym-shalo+1:nym,:,1) ! south
       data(1:nxm+ioff,               nym+1+joff:nym+nhalo+joff,:) = data_all(1:nxm+ioff,        1+joff:nhalo+joff,:,1) ! north
       data(nxm+1+ioff:nxm+ehalo+ioff,                1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,     ny-shalo+1:ny,:,4) ! southeast
       data(1-whalo:0,                                1-shalo:0,:) = data_all(nx-whalo+1:nx,         ny-shalo+1:ny,:,5) ! southwest
       data(nxm+1+ioff:nxm+ehalo+ioff,nym+1+joff:nym+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff, 1+joff:nhalo+joff,:,2) ! northeast
       data(1-whalo:0,                nym+1+joff:nym+nhalo+joff,:) = data_all(nx-whalo+1:nx,     1+joff:nhalo+joff,:,3) ! northwest
    case(2)
       data(nx+1+ioff:nx+ehalo+ioff,              1:ny+joff,:) = data_all(1+ioff:ehalo+ioff,              1:ny+joff,:,3) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nxm-whalo+1:nxm,                1:ny+joff,:,1) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,                  ny-shalo+1:ny,:,4) ! south
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,              1+joff:nhalo+joff,:,4) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,          ny-shalo+1:ny,:,5) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nxm-whalo+1:nxm,          nym-shalo+1:nym,:,1) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,      1+joff:nhalo+joff,:,5) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nxm-whalo+1:nxm,  ny+1+joff:ny+nhalo+joff,:,1) ! northwest
    case(3)
       data(nx+1+ioff:nx+ehalo+ioff,              1:ny+joff,:) = data_all(1+ioff:ehalo+ioff,              1:ny+joff,:,1) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nx-whalo+1:nx,                  1:ny+joff,:,2) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,                  ny-shalo+1:ny,:,5) ! south
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,              1+joff:nhalo+joff,:,5) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,        nym-shalo+1:nym,:,1) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nx-whalo+1:nx,              ny-shalo+1:ny,:,4) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,ny+1+joff:ny+nhalo+joff,:,1) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nx-whalo+1:nx,          1+joff:nhalo+joff,:,4) ! northwest
    case(4)
       data(nx+1+ioff:nx+ehalo+ioff,              1:ny+joff,:) = data_all(1+ioff:ehalo+ioff,        1:ny+joff,:,5) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nxm-whalo+1:nxm,     ny+1:2*ny+joff,:,1) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,            ny-shalo+1:ny,:,2) ! south
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,        1+joff:nhalo+joff,:,2) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,    ny-shalo+1:ny,:,3) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nxm-whalo+1:nxm,      ny-shalo+1:ny,:,1) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,1+joff:nhalo+joff,:,3) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nxm-whalo+1:nxm,  1+joff:nhalo+joff,:,1) ! northwest
    case(5)
       data(nx+1+ioff:nx+ehalo+ioff,            1:  ny+joff,:) = data_all(1+ioff:ehalo+ioff,   ny+1:2*ny+joff,:,1) ! east
       data(1-whalo:0,                            1:ny+joff,:) = data_all(nx-whalo+1:nx,            1:ny+joff,:,4) ! west
       data(1:nx+ioff,                            1-shalo:0,:) = data_all(1:nx+ioff,            ny-shalo+1:ny,:,3) ! south
       data(1:nx+ioff,              ny+1+joff:ny+nhalo+joff,:) = data_all(1:nx+ioff,        1+joff:nhalo+joff,:,3) ! north
       data(nx+1+ioff:nx+ehalo+ioff,              1-shalo:0,:) = data_all(1+ioff:ehalo+ioff,    ny-shalo+1:ny,:,1) ! southeast
       data(1-whalo:0,                            1-shalo:0,:) = data_all(nx-whalo+1:nx,        ny-shalo+1:ny,:,2) ! southwest
       data(nx+1+ioff:nx+ehalo+ioff,ny+1+joff:ny+nhalo+joff,:) = data_all(1+ioff:ehalo+ioff,1+joff:nhalo+joff,:,1) ! northeast
       data(1-whalo:0,              ny+1+joff:ny+nhalo+joff,:) = data_all(nx-whalo+1:nx,    1+joff:nhalo+joff,:,2) ! northwest
    end select

  end subroutine fill_five_tile_halo

  !#######################################################################################
  subroutine test_get_boundary(type)
     character(len=*), intent(in)  :: type

     type(domain2D)       :: domain, domain_nonsym
     integer              :: ntiles, num_contact, npes_per_tile, ntile_per_pe, layout(2)
     integer              :: n, l, isc, iec, jsc, jec, ism, iem, jsm, jem
     integer, allocatable, dimension(:)       :: tile, ni, nj, pe_start, pe_end
     integer, allocatable, dimension(:,:)     :: layout2D, global_indices
     real,    allocatable, dimension(:,:,:)   :: ebuffer,   sbuffer,   wbuffer,   nbuffer
     real,    allocatable, dimension(:,:,:)   :: ebuffer1,  sbuffer1,  wbuffer1,  nbuffer1
     real,    allocatable, dimension(:,:,:)   :: ebuffer2,  sbuffer2,  wbuffer2,  nbuffer2
     real,    allocatable, dimension(:,:,:)   :: ebound,    sbound,    wbound,    nbound
     real,    allocatable, dimension(:,:,:)   :: ebufferx,  sbufferx,  wbufferx,  nbufferx
     real,    allocatable, dimension(:,:,:)   :: ebufferx1, sbufferx1, wbufferx1, nbufferx1
     real,    allocatable, dimension(:,:,:)   :: ebufferx2, sbufferx2, wbufferx2, nbufferx2
     real,    allocatable, dimension(:,:,:)   :: eboundx,   sboundx,   wboundx,   nboundx
     real,    allocatable, dimension(:,:,:)   :: ebuffery,  sbuffery,  wbuffery,  nbuffery
     real,    allocatable, dimension(:,:,:)   :: ebuffery1, sbuffery1, wbuffery1, nbuffery1
     real,    allocatable, dimension(:,:,:)   :: ebuffery2, sbuffery2, wbuffery2, nbuffery2
     real,    allocatable, dimension(:,:,:)   :: eboundy,   sboundy,   wboundy,   nboundy
     real,    allocatable, dimension(:,:,:,:) :: global_all, global1_all, global2_all
     real,    allocatable, dimension(:,:,:,:) :: global, global1, global2
     real,    allocatable, dimension(:,:,:,:) :: x, x1, x2, y, y1, y2
     real,    allocatable, dimension(:,:)     :: u_nonsym, v_nonsym
     logical    :: folded_north = .false.
     logical    :: is_torus = .false.
     integer    :: nx_save, ny_save

     nx_save    = nx
     ny_save    = ny

     !--- check the type
    select case(type)
    case ( 'Four-Tile' ) !--- cyclic along both x- and y-direction.
       ntiles = 4
       num_contact = 8
    case ( 'Cubic-Grid' )
       ntiles = 6
       num_contact = 12
       nx = nx_cubic
       ny = nx
    case ( 'Folded-north' )
       folded_north = .true.
       ntiles = 1
    case ( 'torus' )
       is_torus = .true.
       ntiles = 1
    case default
       call mpp_error(FATAL, 'TEST_MPP_DOMAINS: no such test: '//type)
    end select

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    allocate(ni(ntiles), nj(ntiles))
    ni(:) = nx; nj(:) = ny
    if( mod(npes, ntiles) == 0 ) then
       npes_per_tile = npes/ntiles
       write(outunit,*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                       '", each tile will be distributed over ', npes_per_tile, ' processors.'
       ntile_per_pe = 1
       allocate(tile(ntile_per_pe))
       tile = pe/npes_per_tile+1
       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do
    else if ( mod(ntiles, npes) == 0 ) then
       ntile_per_pe = ntiles/npes
       write(outunit,*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                        '", there will be ', ntile_per_pe, ' tiles on each processor.'
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
       call mpp_error(NOTE,'TEST_MPP_DOMAINS: npes should be multiple of ntiles or ' // &
            'ntiles should be multiple of npes. No test is done for '//trim(type) )
       return
    end if

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

     select case(type)
     case("Four-Tile")
        call define_fourtile_mosaic(type, domain, (/nx,nx,nx,nx/), (/ny,ny,ny,ny/), global_indices, &
                                    layout2D, pe_start, pe_end, .true. )
     case("Cubic-Grid")
        call define_cubic_mosaic(type, domain, ni, nj, global_indices, layout2D, pe_start, pe_end )
     case("Folded-north")
        call mpp_define_domains((/1,nx,1,ny/), layout, domain, &
                          xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.true., name='tripolar' )
        call mpp_define_domains((/1,nx,1,ny/), layout, domain_nonsym, &
                          xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.false., name='tripolar' )
     case("torus")
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 yflags=CYCLIC_GLOBAL_DOMAIN, symmetry=.true., name=type)
     end select

    !--- Test the get_boundary of the data at C-cell center.
    allocate(global_all(1:nx+1,1:ny+1,nz, ntiles) )
    allocate(global(1:nx+1,1:ny+1,nz, ntile_per_pe) )
    global = 0
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+1
             do i = 1, nx+1
                global_all(i,j,k,l) = l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global(:,:,:,n) = global_all(:,:,:,tile(n))
    end do

    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_memory_domain   ( domain, ism, iem, jsm, jem )
    allocate( x (ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x1(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x2(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    x = 0.
    x(isc:iec+1,jsc:jec+1,:,:) = global(isc:iec+1,jsc:jec+1,:,:)
    x1 = x; x2 = x*10

    !--- buffer allocation
    allocate(ebuffer(jsc:jec+1, nz, ntile_per_pe), wbuffer(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbuffer(isc:iec+1, nz, ntile_per_pe), nbuffer(isc:iec+1, nz, ntile_per_pe))
    allocate(ebuffer1(jsc:jec+1, nz, ntile_per_pe), wbuffer1(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbuffer1(isc:iec+1, nz, ntile_per_pe), nbuffer1(isc:iec+1, nz, ntile_per_pe))
    allocate(ebuffer2(jsc:jec+1, nz, ntile_per_pe), wbuffer2(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbuffer2(isc:iec+1, nz, ntile_per_pe), nbuffer2(isc:iec+1, nz, ntile_per_pe))
    allocate(ebound(jsc:jec+1, nz, ntile_per_pe), wbound(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbound(isc:iec+1, nz, ntile_per_pe), nbound(isc:iec+1, nz, ntile_per_pe))
    ebound  = 0; ebuffer = 0; ebuffer1 = 0; ebuffer2 = 0
    sbound  = 0; sbuffer = 0; sbuffer1 = 0; sbuffer2 = 0
    wbound  = 0; wbuffer = 0; wbuffer1 = 0; wbuffer2 = 0
    nbound  = 0; nbuffer = 0; nbuffer1 = 0; nbuffer2 = 0

    do n = 1, ntile_per_pe
       if(folded_north .or. is_torus ) then
          call mpp_get_boundary(x(:,:,:,n), domain, sbuffer=sbuffer(:,:,n), wbuffer=wbuffer(:,:,n), &
                                position=CORNER, tile_count=n  )
       else
          call mpp_get_boundary(x(:,:,:,n), domain, ebuffer=ebuffer(:,:,n), sbuffer=sbuffer(:,:,n), wbuffer=wbuffer(:,:,n), &
                                nbuffer=nbuffer(:,:,n), position=CORNER, tile_count=n  )
       endif
    end do

    !--- multiple variable
    do n = 1, ntile_per_pe
       if(folded_north .or. is_torus) then
          call mpp_get_boundary(x1(:,:,:,n), domain, sbuffer=sbuffer1(:,:,n), wbuffer=wbuffer1(:,:,n), &
               position=CORNER, tile_count=n, complete = .false.  )
          call mpp_get_boundary(x2(:,:,:,n), domain, sbuffer=sbuffer2(:,:,n), wbuffer=wbuffer2(:,:,n), &
               position=CORNER, tile_count=n, complete = .true.  )
       else
          call mpp_get_boundary(x1(:,:,:,n), domain, ebuffer=ebuffer1(:,:,n), sbuffer=sbuffer1(:,:,n), wbuffer=wbuffer1(:,:,n), &
               nbuffer=nbuffer1(:,:,n), position=CORNER, tile_count=n, complete = .false.  )
          call mpp_get_boundary(x2(:,:,:,n), domain, ebuffer=ebuffer2(:,:,n), sbuffer=sbuffer2(:,:,n), wbuffer=wbuffer2(:,:,n), &
               nbuffer=nbuffer2(:,:,n), position=CORNER, tile_count=n, complete = .true.  )
       endif
    end do

    !--- compare the buffer.
    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), ebound(:,:,n), sbound(:,:,n), wbound(:,:,n), nbound(:,:,n) )
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global_all, global_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, ebound(:,:,n), sbound(:,:,n), wbound(:,:,n), nbound(:,:,n) )
       end do
    case("Folded-north")
       !---- folded line update
       global_all(nx/2+2:nx, ny+1,:,1) = global_all(nx/2:2:-1, ny+1,:,1)
       do n = 1, ntile_per_pe
          call fill_folded_north_bound(global_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sbound(:,:,n), wbound(:,:,n) )
       end do
    case("torus")
       do n = 1, ntile_per_pe
          call fill_torus_bound(global_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sbound(:,:,n), wbound(:,:,n) )
       end do
    end select

    if(.not. folded_north .AND. .not. is_torus) then
       call compare_checksums( ebound, ebuffer(:,:,:),  "east bound of "//trim(type) )
       call compare_checksums( nbound, nbuffer(:,:,:),  "north bound of "//trim(type) )
       call compare_checksums( ebound, ebuffer1(:,:,:),  "east bound of "//trim(type)//" X1" )
       call compare_checksums( nbound, nbuffer1(:,:,:),  "north bound of "//trim(type)//" X1" )
       call compare_checksums( ebound*10, ebuffer2(:,:,:),  "east bound of "//trim(type)//" X2" )
       call compare_checksums( nbound*10, nbuffer2(:,:,:),  "north bound of "//trim(type)//" X2" )
    endif
    call compare_checksums( sbound, sbuffer(:,:,:),  "south bound of "//trim(type) )
    call compare_checksums( wbound, wbuffer(:,:,:),  "west bound of "//trim(type) )
    call compare_checksums( sbound, sbuffer1(:,:,:),  "south bound of "//trim(type)//" X1" )
    call compare_checksums( wbound, wbuffer1(:,:,:),  "west bound of "//trim(type)//" X1" )
    call compare_checksums( sbound*10, sbuffer2(:,:,:),  "south bound of "//trim(type)//" X2" )
    call compare_checksums( wbound*10, wbuffer2(:,:,:),  "west bound of "//trim(type)//" X2" )

    !--- release memory
    deallocate(global, global_all, x, x1, x2)
    deallocate(ebuffer, sbuffer, wbuffer, nbuffer)
    deallocate(ebuffer1, sbuffer1, wbuffer1, nbuffer1)
    deallocate(ebuffer2, sbuffer2, wbuffer2, nbuffer2)
    deallocate(ebound, sbound, wbound, nbound )

    !-------------------------------------------------------------------------------------------
    !
    !             Test SCALAR_PAIR BGRID
    !
    !-------------------------------------------------------------------------------------------
    allocate(global1_all(1:nx+1,1:ny+1,nz, ntiles) )
    allocate(global2_all(1:nx+1,1:ny+1,nz, ntiles) )
    allocate(global1(1:nx+1,1:ny+1,nz, ntile_per_pe) )
    allocate(global2(1:nx+1,1:ny+1,nz, ntile_per_pe) )
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+1
             do i = 1, nx+1
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global1(:,:,:,n) = global1_all(:,:,:,tile(n))
       global2(:,:,:,n) = global2_all(:,:,:,tile(n))
    end do
    allocate( x (ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x1(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( x2(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( y (ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( y1(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    allocate( y2(ism:iem+1,jsm:jem+1,nz, ntile_per_pe) )
    x = 0.; y = 0
    if( trim(type) == "Folded-north" ) then
       x(isc+1:iec+1,jsc+1:jec+1,:,:) = global1(isc+1:iec+1,jsc+1:jec+1,:,:)
       y(isc+1:iec+1,jsc+1:jec+1,:,:) = global2(isc+1:iec+1,jsc+1:jec+1,:,:)
    else
       x(isc:iec+1,jsc:jec+1,:,:) = global1(isc:iec+1,jsc:jec+1,:,:)
       y(isc:iec+1,jsc:jec+1,:,:) = global2(isc:iec+1,jsc:jec+1,:,:)
    endif
    x1 = x; x2 = x*10
    y1 = y; y2 = y*10

    !--- buffer allocation
    allocate(ebufferx(jsc:jec+1, nz, ntile_per_pe), wbufferx(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbufferx(isc:iec+1, nz, ntile_per_pe), nbufferx(isc:iec+1, nz, ntile_per_pe))
    allocate(ebufferx1(jsc:jec+1, nz, ntile_per_pe), wbufferx1(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbufferx1(isc:iec+1, nz, ntile_per_pe), nbufferx1(isc:iec+1, nz, ntile_per_pe))
    allocate(ebufferx2(jsc:jec+1, nz, ntile_per_pe), wbufferx2(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbufferx2(isc:iec+1, nz, ntile_per_pe), nbufferx2(isc:iec+1, nz, ntile_per_pe))
    allocate(eboundx(jsc:jec+1, nz, ntile_per_pe), wboundx(jsc:jec+1, nz, ntile_per_pe))
    allocate(sboundx(isc:iec+1, nz, ntile_per_pe), nboundx(isc:iec+1, nz, ntile_per_pe))
    allocate(ebuffery(jsc:jec+1, nz, ntile_per_pe), wbuffery(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbuffery(isc:iec+1, nz, ntile_per_pe), nbuffery(isc:iec+1, nz, ntile_per_pe))
    allocate(ebuffery1(jsc:jec+1, nz, ntile_per_pe), wbuffery1(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbuffery1(isc:iec+1, nz, ntile_per_pe), nbuffery1(isc:iec+1, nz, ntile_per_pe))
    allocate(ebuffery2(jsc:jec+1, nz, ntile_per_pe), wbuffery2(jsc:jec+1, nz, ntile_per_pe))
    allocate(sbuffery2(isc:iec+1, nz, ntile_per_pe), nbuffery2(isc:iec+1, nz, ntile_per_pe))
    allocate(eboundy(jsc:jec+1, nz, ntile_per_pe), wboundy(jsc:jec+1, nz, ntile_per_pe))
    allocate(sboundy(isc:iec+1, nz, ntile_per_pe), nboundy(isc:iec+1, nz, ntile_per_pe))
    eboundx  = 0; ebufferx = 0; ebufferx1 = 0; ebufferx2 = 0
    sboundx  = 0; sbufferx = 0; sbufferx1 = 0; sbufferx2 = 0
    wboundx  = 0; wbufferx = 0; wbufferx1 = 0; wbufferx2 = 0
    nboundx  = 0; nbufferx = 0; nbufferx1 = 0; nbufferx2 = 0
    eboundy  = 0; ebuffery = 0; ebuffery1 = 0; ebuffery2 = 0
    sboundy  = 0; sbuffery = 0; sbuffery1 = 0; sbuffery2 = 0
    wboundy  = 0; wbuffery = 0; wbuffery1 = 0; wbuffery2 = 0
    nboundy  = 0; nbuffery = 0; nbuffery1 = 0; nbuffery2 = 0


    do n = 1, ntile_per_pe
       if(folded_north .or. is_torus) then
          call mpp_get_boundary(x(:,:,:,n), y(:,:,:,n), domain, sbufferx=sbufferx(:,:,n), wbufferx=wbufferx(:,:,n), &
               sbuffery=sbuffery(:,:,n), wbuffery=wbuffery(:,:,n), gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR  )
       else
          call mpp_get_boundary(x(:,:,:,n), y(:,:,:,n), domain, ebufferx=ebufferx(:,:,n), sbufferx=sbufferx(:,:,n), &
               wbufferx=wbufferx(:,:,n), nbufferx=nbufferx(:,:,n), ebuffery=ebuffery(:,:,n),       &
               sbuffery=sbuffery(:,:,n), wbuffery=wbuffery(:,:,n), nbuffery=nbuffery(:,:,n),       &
               gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR  )
       endif
    end do

    do n = 1, ntile_per_pe
       if(folded_north .or. is_torus) then
          call mpp_get_boundary(x1(:,:,:,n), y1(:,:,:,n), domain, sbufferx=sbufferx1(:,:,n), wbufferx=wbufferx1(:,:,n), &
               sbuffery=sbuffery1(:,:,n), wbuffery=wbuffery1(:,:,n),                           &
               gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR, complete = .false.  )
          call mpp_get_boundary(x2(:,:,:,n), y2(:,:,:,n), domain, sbufferx=sbufferx2(:,:,n), wbufferx=wbufferx2(:,:,n), &
               sbuffery=sbuffery2(:,:,n), wbuffery=wbuffery2(:,:,n),       &
               gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR, complete = .true.  )
       else
          call mpp_get_boundary(x1(:,:,:,n), y1(:,:,:,n), domain, ebufferx=ebufferx1(:,:,n), sbufferx=sbufferx1(:,:,n), &
               wbufferx=wbufferx1(:,:,n), nbufferx=nbufferx1(:,:,n), ebuffery=ebuffery1(:,:,n),       &
               sbuffery=sbuffery1(:,:,n), wbuffery=wbuffery1(:,:,n), nbuffery=nbuffery1(:,:,n),       &
               gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR, complete = .false.  )
          call mpp_get_boundary(x2(:,:,:,n), y2(:,:,:,n), domain, ebufferx=ebufferx2(:,:,n), sbufferx=sbufferx2(:,:,n), &
               wbufferx=wbufferx2(:,:,n), nbufferx=nbufferx2(:,:,n), ebuffery=ebuffery2(:,:,n),       &
               sbuffery=sbuffery2(:,:,n), wbuffery=wbuffery2(:,:,n), nbuffery=nbuffery2(:,:,n),       &
               gridtype=BGRID_NE, tile_count=n, flags = SCALAR_PAIR, complete = .true.  )
       endif
    end do

    !--- compare the buffer.
    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )
          call fill_four_tile_bound(global2_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all, global2_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )
          call fill_cubic_grid_bound(global2_all, global1_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )
       end do
    case("Folded-north")
       global1_all(nx/2+2:nx, ny+1,:,1) = global1_all(nx/2:2:-1, ny+1,:,1)
       global2_all(nx/2+2:nx, ny+1,:,1) = global2_all(nx/2:2:-1, ny+1,:,1)
       do n = 1, ntile_per_pe
          call fill_folded_north_bound(global1_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundx(:,:,n), wboundx(:,:,n) )
          call fill_folded_north_bound(global2_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundy(:,:,n), wboundy(:,:,n) )
       end do
    case("torus")
       do n = 1, ntile_per_pe
          call fill_torus_bound(global1_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundx(:,:,n), wboundx(:,:,n) )
          call fill_torus_bound(global2_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundy(:,:,n), wboundy(:,:,n) )
       end do
    end select

    if(.not. folded_north .AND. .not. is_torus ) then
       call compare_checksums( eboundx, ebufferx(:,:,:),   "east bound of SCALAR_PAIR BGRID " //trim(type)//" X" )
       call compare_checksums( nboundx, nbufferx(:,:,:),   "north bound of SCALAR_PAIR BGRID "//trim(type)//" X" )
       call compare_checksums( eboundy, ebuffery(:,:,:),   "east bound of SCALAR_PAIR BGRID " //trim(type)//" Y" )
       call compare_checksums( nboundy, nbuffery(:,:,:),   "north bound of SCALAR_PAIR BGRID "//trim(type)//" Y" )
       call compare_checksums( eboundx, ebufferx1(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" X1" )
       call compare_checksums( nboundx, nbufferx1(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" X1" )
       call compare_checksums( eboundy, ebuffery1(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" Y1" )
       call compare_checksums( nboundy, nbuffery1(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" Y1" )
    endif

    call compare_checksums( sboundx, sbufferx(:,:,:),   "south bound of SCALAR_PAIR BGRID "//trim(type)//" X" )
    call compare_checksums( wboundx, wbufferx(:,:,:),   "west bound of SCALAR_PAIR BGRID " //trim(type)//" X" )
    call compare_checksums( sboundy, sbuffery(:,:,:),   "south bound of SCALAR_PAIR BGRID "//trim(type)//" Y" )
    call compare_checksums( wboundy, wbuffery(:,:,:),   "west bound of SCALAR_PAIR BGRID " //trim(type)//" Y" )
    call compare_checksums( sboundx, sbufferx1(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" X1" )
    call compare_checksums( wboundx, wbufferx1(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" X1" )
    call compare_checksums( sboundy, sbuffery1(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" Y1" )
    call compare_checksums( wboundy, wbuffery1(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" Y1" )

    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )
          call fill_four_tile_bound(global2_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all*10, global2_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )
          call fill_cubic_grid_bound(global2_all*10, global1_all*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, 1, eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )
       end do
    case("Folded-north")
       do n = 1, ntile_per_pe
          call fill_folded_north_bound(global1_all(:,:,:,1)*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundx(:,:,n), wboundx(:,:,n) )
          call fill_folded_north_bound(global2_all(:,:,:,1)*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundy(:,:,n), wboundy(:,:,n) )
       end do
    case("torus")
       do n = 1, ntile_per_pe
          call fill_torus_bound(global1_all(:,:,:,1)*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundx(:,:,n), wboundx(:,:,n) )
          call fill_torus_bound(global2_all(:,:,:,1)*10, isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundy(:,:,n), wboundy(:,:,n) )
       end do
    end select

    if(.not. folded_north .AND. .not. is_torus ) then
       call compare_checksums( eboundx, ebufferx2(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" X2" )
       call compare_checksums( nboundx, nbufferx2(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" X2" )
       call compare_checksums( eboundy, ebuffery2(:,:,:),  "east bound of SCALAR_PAIR BGRID " //trim(type)//" Y2" )
       call compare_checksums( nboundy, nbuffery2(:,:,:),  "north bound of SCALAR_PAIR BGRID "//trim(type)//" Y2" )
    endif
    call compare_checksums( sboundx, sbufferx2(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" X2" )
    call compare_checksums( wboundx, wbufferx2(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" X2" )
    call compare_checksums( sboundy, sbuffery2(:,:,:),  "south bound of SCALAR_PAIR BGRID "//trim(type)//" Y2" )
    call compare_checksums( wboundy, wbuffery2(:,:,:),  "west bound of SCALAR_PAIR BGRID " //trim(type)//" Y2" )

    !-------------------------------------------------------------------------------------------
    !
    !             Test 2-D Vector BGRID
    !
    !-------------------------------------------------------------------------------------------
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny+1
             do i = 1, nx+1
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    x = 0.; y = 0
    eboundx  = 0; ebufferx = 0; ebufferx1 = 0; ebufferx2 = 0
    sboundx  = 0; sbufferx = 0; sbufferx1 = 0; sbufferx2 = 0
    wboundx  = 0; wbufferx = 0; wbufferx1 = 0; wbufferx2 = 0
    nboundx  = 0; nbufferx = 0; nbufferx1 = 0; nbufferx2 = 0
    eboundy  = 0; ebuffery = 0; ebuffery1 = 0; ebuffery2 = 0
    sboundy  = 0; sbuffery = 0; sbuffery1 = 0; sbuffery2 = 0
    wboundy  = 0; wbuffery = 0; wbuffery1 = 0; wbuffery2 = 0
    nboundy  = 0; nbuffery = 0; nbuffery1 = 0; nbuffery2 = 0

    x(isc:iec+1,jsc:jec+1,1,:) = global1(isc:iec+1,jsc:jec+1,1,:)
    y(isc:iec+1,jsc:jec+1,1,:) = global2(isc:iec+1,jsc:jec+1,1,:)

    do n = 1, ntile_per_pe
       if(folded_north .or. is_torus ) then
          call mpp_get_boundary(x(:,:,1,n), y(:,:,1,n), domain, sbufferx=sbufferx(:,1,n), wbufferx=wbufferx(:,1,n), &
               sbuffery=sbuffery(:,1,n), wbuffery=wbuffery(:,1,n), gridtype=BGRID_NE, tile_count=n)

       else
          call mpp_get_boundary(x(:,:,1,n), y(:,:,1,n), domain, ebufferx=ebufferx(:,1,n), sbufferx=sbufferx(:,1,n), &
               wbufferx=wbufferx(:,1,n), nbufferx=nbufferx(:,1,n), ebuffery=ebuffery(:,1,n),       &
               sbuffery=sbuffery(:,1,n), wbuffery=wbuffery(:,1,n), nbuffery=nbuffery(:,1,n),       &
               gridtype=BGRID_NE, tile_count=n)
       endif
    end do

    if(folded_north) then
       allocate(u_nonsym(ism:iem,jsm:jem), v_nonsym(ism:iem,jsm:jem))
       u_nonsym = 0.0; v_nonsym = 0.0
       u_nonsym(isc:iec,jsc:jec) = global1(isc+1:iec+1,jsc+1:jec+1,1,1)
       v_nonsym(isc:iec,jsc:jec) = global2(isc+1:iec+1,jsc+1:jec+1,1,1)
       call mpp_update_domains(u_nonsym, v_nonsym, domain_nonsym, gridtype=BGRID_NE)
       !--- comparing boundary data
       do i = isc,iec+1
          if(i==1) cycle
          if(sbufferx(i,1,1) .NE. u_nonsym(i-1,jsc-1)) then
             print*,"pe ", mpp_pe(), i, jsc-1, sbufferx(i,1,1), u_nonsym(i-1,jsc-1)
             call mpp_error(FATAL, "test_get_boundary: mismatch of sbufferx")
          endif
       enddo
       call mpp_error(NOTE,"test_get_boundary: reproduce non-symmetric halo update for sbufferx")

       do i = isc,iec+1
          if(i==1) cycle
          if(sbuffery(i,1,1) .NE. v_nonsym(i-1,jsc-1)) then
             print*,"pe ", mpp_pe(), i, jsc-1, sbufferx(i,1,1), v_nonsym(i-1,jsc-1)
             call mpp_error(FATAL, "test_get_boundary: mismatch of sbuffery")
          endif
       enddo
      call mpp_error(NOTE,"test_get_boundary: reproduce non-symmetric halo update for sbuffery")

       do j = jsc,jec+1
          if(j == 1) cycle
          if(wbufferx(j,1,1) .NE. u_nonsym(isc-1,j-1)) then
             print*,"pe ", mpp_pe(), isc-1, j, wbufferx(j,1,1), u_nonsym(isc-1,j-1)
             call mpp_error(FATAL, "test_get_boundary: mismatch of wbufferx")
          endif
       enddo
      call mpp_error(NOTE,"test_get_boundary: reproduce non-symmetric halo update for wbufferx")

       do j = jsc,jec+1
          if(j==1) cycle
          if(wbuffery(j,1,1) .NE. v_nonsym(isc-1,j-1)) then
             print*,"pe ", mpp_pe(), isc-1, j, wbuffery(j,1,1), v_nonsym(isc-1,j-1)
             call mpp_error(FATAL, "test_get_boundary: mismatch of wbuffery")
          endif
       enddo
       call mpp_error(NOTE,"test_get_boundary: reproduce non-symmetric halo update for wbuffery")

       deallocate(u_nonsym, v_nonsym)

    endif

    !--- compare the buffer.
    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )
          call fill_four_tile_bound(global2_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all, global2_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), 1, -1, eboundx(:,:,n), sboundx(:,:,n), wboundx(:,:,n), nboundx(:,:,n) )
          call fill_cubic_grid_bound(global2_all, global1_all, isc, iec, jsc, jec, 1, 1, &
               tile(n), -1, 1, eboundy(:,:,n), sboundy(:,:,n), wboundy(:,:,n), nboundy(:,:,n) )
       end do
    case("Folded-north")
       global1_all(nx/2+2:nx, ny+1,:,1) = -global1_all(nx/2:2:-1, ny+1,:,1)
       global2_all(nx/2+2:nx, ny+1,:,1) = -global2_all(nx/2:2:-1, ny+1,:,1)
       global1_all(1, ny+1,:,1) = 0
       global2_all(1, ny+1,:,1) = 0
       global1_all(nx/2+1, ny+1,:,1) = 0
       global2_all(nx/2+1, ny+1,:,1) = 0
       global1_all(nx+1, ny+1,:,1) = 0
       global2_all(nx+1, ny+1,:,1) = 0


       do n = 1, ntile_per_pe
          call fill_folded_north_bound(global1_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundx(:,:,n), wboundx(:,:,n) )
          call fill_folded_north_bound(global2_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundy(:,:,n), wboundy(:,:,n) )
          ! set wboundx and wbouny to zero at pole (i=1, nx/2+1, nx+1)
!          if( jec == ny ) then
!             if( isc == 1 .OR. isc == nx/2+1 .OR. isc == nx+1 ) then
!                wboundx(jec+1,:,n) = 0
!                wboundy(jec+1,:,n) = 0
!             endif
!         endif
       end do
    case("torus")
       do n = 1, ntile_per_pe
          call fill_torus_bound(global1_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundx(:,:,n), wboundx(:,:,n) )
          call fill_torus_bound(global2_all(:,:,:,1), isc, iec, jsc, jec, 1, 1, &
               tile(n), sboundy(:,:,n), wboundy(:,:,n) )
       enddo
    end select

    if(.not. folded_north .AND. .not. is_torus ) then
       call compare_checksums( eboundx(:,1:1,:), ebufferx(:,1:1,:),   "east bound of 2-D BGRID " //trim(type)//" X" )
       call compare_checksums( nboundx(:,1:1,:), nbufferx(:,1:1,:),   "north bound of 2-D BGRID "//trim(type)//" X" )
       call compare_checksums( eboundy(:,1:1,:), ebuffery(:,1:1,:),   "east bound of 2-D BGRID " //trim(type)//" Y" )
       call compare_checksums( nboundy(:,1:1,:), nbuffery(:,1:1,:),   "north bound of 2-D BGRID "//trim(type)//" Y" )
    endif

    call compare_checksums( sboundx(:,1:1,:), sbufferx(:,1:1,:),   "south bound of 2-D BGRID "//trim(type)//" X" )
    call compare_checksums( wboundx(:,1:1,:), wbufferx(:,1:1,:),   "west bound of 2-D BGRID " //trim(type)//" X" )
    call compare_checksums( sboundy(:,1:1,:), sbuffery(:,1:1,:),   "south bound of 2-D BGRID "//trim(type)//" Y" )
    call compare_checksums( wboundy(:,1:1,:), wbuffery(:,1:1,:),   "west bound of 2-D BGRID " //trim(type)//" Y" )


    !--- release memory
    deallocate(global1, global1_all, global2, global2_all)
    deallocate(x, y, x1, y1, x2, y2)
    deallocate(ebufferx, sbufferx, wbufferx, nbufferx)
    deallocate(ebufferx1, sbufferx1, wbufferx1, nbufferx1)
    deallocate(ebufferx2, sbufferx2, wbufferx2, nbufferx2)
    deallocate(ebuffery, sbuffery, wbuffery, nbuffery)
    deallocate(ebuffery1, sbuffery1, wbuffery1, nbuffery1)
    deallocate(ebuffery2, sbuffery2, wbuffery2, nbuffery2)
    deallocate(eboundx, sboundx, wboundx, nboundx )
    deallocate(eboundy, sboundy, wboundy, nboundy )

    !-------------------------------------------------------------------------------------------
    !
    !             Test VECTOR CGRID
    !
    !-------------------------------------------------------------------------------------------
    allocate(global1_all(1:nx+1,1:ny,  nz, ntiles) )
    allocate(global2_all(1:nx,  1:ny+1,nz, ntiles) )
    allocate(global1(1:nx+1,1:ny,  nz, ntile_per_pe) )
    allocate(global2(1:nx,  1:ny+1,nz, ntile_per_pe) )
    do l = 1, ntiles
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx+1
                global1_all(i,j,k,l) = 1.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
          do j = 1, ny+1
             do i = 1, nx
                global2_all(i,j,k,l) = 2.0e3 + l + i*1.0e-3 + j*1.0e-6 + k*1.0e-9
             end do
          end do
       end do
    end do

    do n = 1, ntile_per_pe
       global1(:,:,:,n) = global1_all(:,:,:,tile(n))
       global2(:,:,:,n) = global2_all(:,:,:,tile(n))
    end do
    allocate( x (ism:iem+1,jsm:jem,  nz, ntile_per_pe) )
    allocate( x1(ism:iem+1,jsm:jem,  nz, ntile_per_pe) )
    allocate( x2(ism:iem+1,jsm:jem,  nz, ntile_per_pe) )
    allocate( y (ism:iem,  jsm:jem+1,nz, ntile_per_pe) )
    allocate( y1(ism:iem,  jsm:jem+1,nz, ntile_per_pe) )
    allocate( y2(ism:iem,  jsm:jem+1,nz, ntile_per_pe) )
    x = 0.; y = 0
    x(isc:iec+1,jsc:jec,  :,:) = global1(isc:iec+1,jsc:jec,  :,:)
    y(isc:iec,  jsc:jec+1,:,:) = global2(isc:iec,  jsc:jec+1,:,:)
    x1 = x; x2 = x*10
    y1 = y; y2 = y*10

    !--- buffer allocation
    allocate(ebufferx(jec-jsc+1, nz, ntile_per_pe), wbufferx(jec-jsc+1, nz, ntile_per_pe))
    allocate(sbufferx(iec-isc+2, nz, ntile_per_pe), nbufferx(iec-isc+2, nz, ntile_per_pe))
    allocate(ebufferx1(jec-jsc+1, nz, ntile_per_pe), wbufferx1(jec-jsc+1, nz, ntile_per_pe))
    allocate(sbufferx1(iec-isc+2, nz, ntile_per_pe), nbufferx1(iec-isc+2, nz, ntile_per_pe))
    allocate(ebufferx2(jec-jsc+1, nz, ntile_per_pe), wbufferx2(jec-jsc+1, nz, ntile_per_pe))
    allocate(sbufferx2(iec-isc+2, nz, ntile_per_pe), nbufferx2(iec-isc+2, nz, ntile_per_pe))
    allocate(ebuffery(jec-jsc+2, nz, ntile_per_pe), wbuffery(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery(iec-isc+1, nz, ntile_per_pe), nbuffery(iec-isc+1, nz, ntile_per_pe))
    allocate(ebuffery1(jec-jsc+2, nz, ntile_per_pe), wbuffery1(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery1(iec-isc+1, nz, ntile_per_pe), nbuffery1(iec-isc+1, nz, ntile_per_pe))
    allocate(ebuffery2(jec-jsc+2, nz, ntile_per_pe), wbuffery2(jec-jsc+2, nz, ntile_per_pe))
    allocate(sbuffery2(iec-isc+1, nz, ntile_per_pe), nbuffery2(iec-isc+1, nz, ntile_per_pe))
    allocate(eboundx(jec-jsc+1, nz, ntile_per_pe), wboundx(jec-jsc+1, nz, ntile_per_pe))
    allocate(sboundy(iec-isc+1, nz, ntile_per_pe), nboundy(iec-isc+1, nz, ntile_per_pe))
    eboundx  = 0; ebufferx = 0; ebufferx1 = 0; ebufferx2 = 0
    wboundx  = 0; wbufferx = 0; wbufferx1 = 0; wbufferx2 = 0
    sboundy  = 0; sbuffery = 0; sbuffery1 = 0; sbuffery2 = 0
    nboundy  = 0; nbuffery = 0; nbuffery1 = 0; nbuffery2 = 0


    do n = 1, ntile_per_pe
       if(folded_north .or. is_torus) then
          call mpp_get_boundary(x(:,:,:,n), y(:,:,:,n), domain, wbufferx=wbufferx(:,:,n), &
               sbuffery=sbuffery(:,:,n), gridtype=CGRID_NE, tile_count=n  )
       else
          call mpp_get_boundary(x(:,:,:,n), y(:,:,:,n), domain, ebufferx=ebufferx(:,:,n), wbufferx=wbufferx(:,:,n), &
               sbuffery=sbuffery(:,:,n), nbuffery=nbuffery(:,:,n), gridtype=CGRID_NE, tile_count=n  )
       endif
    end do

    do n = 1, ntile_per_pe
       if( folded_north .or. is_torus ) then
          call mpp_get_boundary(x1(:,:,:,n), y1(:,:,:,n), domain, wbufferx=wbufferx1(:,:,n), &
               sbuffery=sbuffery1(:,:,n), gridtype=CGRID_NE, tile_count=n,  &
               complete = .false.  )
          call mpp_get_boundary(x2(:,:,:,n), y2(:,:,:,n), domain, wbufferx=wbufferx2(:,:,n), &
               sbuffery=sbuffery2(:,:,n), gridtype=CGRID_NE, tile_count=n,  &
               complete = .true.  )
       else
          call mpp_get_boundary(x1(:,:,:,n), y1(:,:,:,n), domain, ebufferx=ebufferx1(:,:,n), wbufferx=wbufferx1(:,:,n), &
               sbuffery=sbuffery1(:,:,n), nbuffery=nbuffery1(:,:,n), gridtype=CGRID_NE, tile_count=n,  &
               complete = .false.  )
          call mpp_get_boundary(x2(:,:,:,n), y2(:,:,:,n), domain, ebufferx=ebufferx2(:,:,n), wbufferx=wbufferx2(:,:,n), &
               sbuffery=sbuffery2(:,:,n), nbuffery=nbuffery2(:,:,n), gridtype=CGRID_NE, tile_count=n,  &
               complete = .true.  )
       endif
    end do

    !--- compare the buffer.
    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all, isc, iec, jsc, jec, 1, 0, &
               tile(n), ebound=eboundx(:,:,n), wbound=wboundx(:,:,n) )
          call fill_four_tile_bound(global2_all, isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all, global2_all, isc, iec, jsc, jec, 1, 0, &
               tile(n), 1, -1, ebound=eboundx(:,:,n), wbound=wboundx(:,:,n)  )
          call fill_cubic_grid_bound(global2_all, global1_all, isc, iec, jsc, jec, 0, 1, &
               tile(n), -1, 1, sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )
       end do
    case("Folded-north")
       do n = 1, ntile_per_pe
          call fill_folded_north_bound(global1_all(:,:,:,1), isc, iec, jsc, jec, 1, 0, &
               tile(n), wbound=wboundx(:,:,n) )
          call fill_folded_north_bound(global2_all(:,:,:,1), isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n) )
       end do
    case("torus")
       do n = 1, ntile_per_pe
          call fill_torus_bound(global1_all(:,:,:,1), isc, iec, jsc, jec, 1, 0, &
               tile(n), wbound=wboundx(:,:,n) )
          call fill_torus_bound(global2_all(:,:,:,1), isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n) )
       end do
    end select

    if(.not. folded_north .and. .not. is_torus ) then
       call compare_checksums( eboundx, ebufferx(:,:,:),   "east bound of CGRID " //trim(type)//" X" )
       call compare_checksums( nboundy, nbuffery(:,:,:),   "north bound of CGRID "//trim(type)//" Y" )
       call compare_checksums( eboundx, ebufferx1(:,:,:),  "east bound of CGRID " //trim(type)//" X1" )
       call compare_checksums( nboundy, nbuffery1(:,:,:),  "north bound of CGRID "//trim(type)//" Y1" )
    endif
    call compare_checksums( wboundx, wbufferx(:,:,:),   "west bound of CGRID " //trim(type)//" X" )
    call compare_checksums( sboundy, sbuffery(:,:,:),   "south bound of CGRID "//trim(type)//" Y" )
    call compare_checksums( wboundx, wbufferx1(:,:,:),  "west bound of CGRID " //trim(type)//" X1" )
    call compare_checksums( sboundy, sbuffery1(:,:,:),  "south bound of CGRID "//trim(type)//" Y1" )

    select case(type)
    case("Four-Tile")
       do n = 1, ntile_per_pe
          call fill_four_tile_bound(global1_all*10, isc, iec, jsc, jec, 1, 0, &
               tile(n), ebound=eboundx(:,:,n), wbound=wboundx(:,:,n) )
          call fill_four_tile_bound(global2_all*10, isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )
       end do
    case("Cubic-Grid")
       do n = 1, ntile_per_pe
          call fill_cubic_grid_bound(global1_all*10, global2_all*10, isc, iec, jsc, jec, 1, 0, &
               tile(n), 1, -1, ebound=eboundx(:,:,n), wbound=wboundx(:,:,n) )
          call fill_cubic_grid_bound(global2_all*10, global1_all*10, isc, iec, jsc, jec, 0, 1, &
               tile(n), -1, 1, sbound=sboundy(:,:,n), nbound=nboundy(:,:,n) )
       end do
    case("Folded-north")
       do n = 1, ntile_per_pe
          call fill_folded_north_bound(global1_all(:,:,:,1)*10, isc, iec, jsc, jec, 1, 0, &
               tile(n), wbound=wboundx(:,:,n) )
          call fill_folded_north_bound(global2_all(:,:,:,1)*10, isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n) )
       end do
    case("torus")
       do n = 1, ntile_per_pe
          call fill_torus_bound(global1_all(:,:,:,1)*10, isc, iec, jsc, jec, 1, 0, &
               tile(n), wbound=wboundx(:,:,n) )
          call fill_torus_bound(global2_all(:,:,:,1)*10, isc, iec, jsc, jec, 0, 1, &
               tile(n), sbound=sboundy(:,:,n) )
       end do
    end select

    if(.not. folded_north .and. .not. is_torus ) then
       call compare_checksums( eboundx, ebufferx2(:,:,:),  "east bound of CGRID " //trim(type)//" X2" )
       call compare_checksums( nboundy, nbuffery2(:,:,:),  "north bound of CGRID "//trim(type)//" Y2" )
    endif
    call compare_checksums( wboundx, wbufferx2(:,:,:),  "west bound of CGRID " //trim(type)//" X2" )
    call compare_checksums( sboundy, sbuffery2(:,:,:),  "south bound of CGRID "//trim(type)//" Y2" )

    !--- release memory
    deallocate(global1, global1_all, global2, global2_all)
    deallocate(x, y, x1, y1, x2, y2)
    deallocate(ebufferx, sbufferx, wbufferx, nbufferx)
    deallocate(ebufferx1, sbufferx1, wbufferx1, nbufferx1)
    deallocate(ebufferx2, sbufferx2, wbufferx2, nbufferx2)
    deallocate(ebuffery, sbuffery, wbuffery, nbuffery)
    deallocate(ebuffery1, sbuffery1, wbuffery1, nbuffery1)
    deallocate(ebuffery2, sbuffery2, wbuffery2, nbuffery2)
    deallocate(eboundx, sboundy, wboundx, nboundy )

    nx = nx_save
    ny = ny_save

  end subroutine test_get_boundary

  !######################################################################################
  subroutine define_fourtile_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end, symmetry )
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    logical,        intent(in)    :: symmetry
    integer, dimension(8)         :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(8)         :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)

    ntiles = 4
    num_contact = 8
    if(size(pe_start(:)) .NE. 4 .OR. size(pe_end(:)) .NE. 4 ) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of pe_start and pe_end should be 4")
    if(size(global_indices,1) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of first dimension of global_indices should be 4")
    if(size(global_indices,2) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of second dimension of global_indices should be 4")
    if(size(layout,1) .NE. 2) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of first dimension of layout should be 2")
    if(size(layout,2) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of second dimension of layout should be 4")
    if(size(ni(:)) .NE. 4 .OR. size(nj(:)) .NE. 4) call mpp_error(FATAL, &
         "define_fourtile_mosaic: size of ni and nj should be 4")

    !--- Contact line 1, between tile 1 (EAST) and tile 2 (WEST)
    tile1(1) = 1; tile2(1) = 2
    istart1(1) = ni(1); iend1(1) = ni(1); jstart1(1) = 1;     jend1(1) = nj(1)
    istart2(1) = 1;     iend2(1) = 1;     jstart2(1) = 1;     jend2(1) = nj(2)
    !--- Contact line 2, between tile 1 (SOUTH) and tile 3 (NORTH)  --- cyclic
    tile1(2) = 1; tile2(2) = 3
    istart1(2) = 1;     iend1(2) = ni(1); jstart1(2) = 1;     jend1(2) = 1
    istart2(2) = 1;     iend2(2) = ni(3); jstart2(2) = nj(3); jend2(2) = nj(3)
    !--- Contact line 3, between tile 1 (WEST) and tile 2 (EAST) --- cyclic
    tile1(3) = 1; tile2(3) = 2
    istart1(3) = 1;     iend1(3) = 1;     jstart1(3) = 1;     jend1(3) = nj(1)
    istart2(3) = ni(2); iend2(3) = ni(2); jstart2(3) = 1;     jend2(3) = nj(2)
    !--- Contact line 4, between tile 1 (NORTH) and tile 3 (SOUTH)
    tile1(4) = 1; tile2(4) = 3
    istart1(4) = 1;     iend1(4) = ni(1); jstart1(4) = nj(1); jend1(4) = nj(1)
    istart2(4) = 1;     iend2(4) = ni(3); jstart2(4) = 1;     jend2(4) = 1
    !--- Contact line 5, between tile 2 (SOUTH) and tile 4 (NORTH) --- cyclic
    tile1(5) = 2; tile2(5) = 4
    istart1(5) = 1;     iend1(5) = ni(2); jstart1(5) = 1;     jend1(5) = 1
    istart2(5) = 1;     iend2(5) = ni(4); jstart2(5) = nj(4); jend2(5) = nj(4)
    !--- Contact line 6, between tile 2 (NORTH) and tile 4 (SOUTH)
    tile1(6) = 2; tile2(6) = 4
    istart1(6) = 1;     iend1(6) = ni(2); jstart1(6) = nj(2); jend1(6) = nj(2)
    istart2(6) = 1;     iend2(6) = ni(4); jstart2(6) = 1;     jend2(6) = 1
    !--- Contact line 7, between tile 3 (EAST) and tile 4 (WEST)
    tile1(7) = 3; tile2(7) = 4
    istart1(7) = ni(3); iend1(7) = ni(3); jstart1(7) = 1;     jend1(7) = nj(3)
    istart2(7) = 1;     iend2(7) = 1;     jstart2(7) = 1;     jend2(7) = nj(4)
    !--- Contact line 8, between tile 3 (WEST) and tile 4 (EAST) --- cyclic
    tile1(8) = 3; tile2(8) = 4
    istart1(8) = 1;     iend1(8) = 1;     jstart1(8) = 1;     jend1(8) = nj(3)
    istart2(8) = ni(4); iend2(8) = ni(4); jstart2(8) = 1;     jend2(8) = nj(4)
    msize(1) = maxval(ni(:)/layout(1,:)) + whalo + ehalo + 1 ! make sure memory domain size is no smaller than
    msize(2) = maxval(nj(:)/layout(2,:)) + shalo + nhalo + 1 ! data domain size
    call mpp_define_mosaic(global_indices, layout, domain, ntiles, num_contact, tile1, tile2,       &
         istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,          &
         pe_start, pe_end, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo,    &
         name = type, memory_size = msize, symmetry = symmetry )

    return

  end subroutine define_fourtile_mosaic

  !#######################################################################################
  !--- define mosaic domain for cubic grid
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

  !#######################################################################################
  subroutine fill_regular_refinement_halo( data, data_all, ni, nj, tm, te, tse, ts, tsw, tw, tnw, tn, tne, ioff, joff )
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data_all
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

  end subroutine fill_regular_refinement_halo

  !##############################################################################
  ! this routine fill the halo points for the refined cubic grid. ioff and joff is used to distinguish
  ! T, C, E, or N-cell
  subroutine fill_cubicgrid_refined_halo(data, data1_all, data2_all, ni, nj, tile, ioff, joff, sign1, sign2)
    real, dimension(1-whalo:,1-shalo:,:), intent(inout) :: data
    real, dimension(:,:,:,:),             intent(in)    :: data1_all, data2_all
    integer, dimension(:),                intent(in)    :: ni, nj
    integer,                              intent(in)    :: tile, ioff, joff, sign1, sign2
    integer                                             :: lw, le, ls, ln

    if(mod(tile,2) == 0) then ! tile 2, 4, 6
       lw = tile - 1; le = tile + 2; ls = tile - 2; ln = tile + 1
       if(le > 6 ) le = le - 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       if( nj(tile) == nj(lw) ) then
          data(1-whalo:0, 1:nj(tile)+joff, :) = data1_all(ni(lw)-whalo+1:ni(lw), 1:nj(lw)+joff, :, lw) ! west
       end if
       if( nj(tile) == ni(le) ) then
          do i = 1, ehalo
             data(ni(tile)+i+ioff, 1:nj(tile)+joff, :)    = sign1*data2_all(ni(le)+joff:1:-1, i+ioff, :, le) ! east
          end do
       end if
       if(ni(tile) == nj(ls) ) then
          do i = 1, shalo
             data(1:ni(tile)+ioff, 1-i, :)     = sign2*data2_all(ni(ls)-i+1, nj(ls)+ioff:1:-1, :, ls) ! south
          end do
       end if
       if(ni(tile) == ni(ln) ) then
          data(1:ni(tile)+ioff, nj(tile)+1+joff:nj(tile)+nhalo+joff, :) = data1_all(1:ni(ln)+ioff, 1+joff:nhalo+joff, :, ln) ! north
       end if
    else ! tile 1, 3, 5
       lw = tile - 2; le = tile + 1; ls = tile - 1; ln = tile + 2
       if(lw < 1 ) lw = lw + 6
       if(ls < 1 ) ls = ls + 6
       if(ln > 6 ) ln = ln - 6
       if(nj(tile) == ni(lw) ) then
          do i = 1, whalo
             data(1-i, 1:nj(tile)+joff, :)     = sign1*data2_all(ni(lw)+joff:1:-1, nj(lw)-i+1, :, lw) ! west
          end do
       end if
       if(nj(tile) == nj(le) ) then
          data(ni(tile)+1+ioff:ni(tile)+ehalo+ioff, 1:nj(tile)+joff, :) = data1_all(1+ioff:ehalo+ioff, 1:nj(le)+joff, :, le) ! east
       end if
       if(ni(tile) == ni(ls) ) then
          data(1:ni(tile)+ioff, 1-shalo:0, :)     = data1_all(1:ni(ls)+ioff, nj(ls)-shalo+1:nj(ls), :, ls) ! south
       end if
       if(ni(tile) == nj(ln) ) then
          do i = 1, nhalo
             data(1:ni(tile)+ioff, nj(tile)+i+joff, :)    = sign2*data2_all(i+joff, nj(ln)+ioff:1:-1, :, ln) ! north
          end do
       end if
    end if

  end subroutine fill_cubicgrid_refined_halo

 !##################################################################################
  subroutine test_subset_update( )
    real, allocatable, dimension(:,:,:) :: x
    type(domain2D) :: domain
    real,    allocatable :: global(:,:,:)
    integer              :: i, xhalo, yhalo
    integer              :: is, ie, js, je, isd, ied, jsd, jed
!   integer :: pes9(9)=(/1,2,3,4,5,6,7,8,9/)
    integer :: pes9(9)=(/0,2,4,10,12,14,20,22,24/)
    integer :: ni, nj

    if(mpp_npes() < 25) then
       call mpp_error(FATAL,"test_mpp_domains: test_subset_update will&
            & not be done when npes < 25")
       return
    endif

    call mpp_declare_pelist(pes9)
    if(any(mpp_pe()==pes9)) then
       call mpp_set_current_pelist(pes9)
       layout = (/3,3/)
       ni = 3; nj =3
       call mpp_define_domains((/1,ni,1,nj/), layout, domain, xhalo=1&
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

    !set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x (isd:ied,jsd:jed,nz) )

    x = 0.
    x (is:ie,js:je,:) = global(is:ie,js:je,:)

!full update
    call mpp_update_domains( x, domain )
    call compare_checksums( x, global(isd:ied,jsd:jed,:), '9pe subset' )

    deallocate(x, global)
    call mpp_deallocate_domain(domain)
  endif

   call mpp_set_current_pelist()

  end subroutine test_subset_update

  !##################################################################################
  subroutine test_halo_update( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x1, x2, x3, x4
    real, allocatable, dimension(:,:,:) :: y, y1, y2, y3, y4
    type(domain2D) :: domain
    real,    allocatable :: global1(:,:,:), global2(:,:,:), global(:,:,:)
    logical, allocatable :: maskmap(:,:)
    integer              :: shift, i, xhalo, yhalo
    logical              :: is_symmetry, folded_south, folded_west, folded_east
    integer              :: is, ie, js, je, isd, ied, jsd, jed

    ! when testing maskmap option, nx*ny should be able to be divided by both npes and npes+1
    if(type == 'Masked' .or. type == 'Masked symmetry') then
       if(mod(nx*ny, npes) .NE. 0 .OR. mod(nx*ny, npes+1) .NE. 0 ) then
          call mpp_error(NOTE,'TEST_MPP_DOMAINS: nx*ny can not be divided by both npes and npes+1, '//&
               'Masked test_halo_update will not be tested')
          return
       end if
    end if

    if(type == 'Folded xy_halo' ) then
       xhalo = max(whalo, ehalo); yhalo = max(shalo, nhalo)
       allocate(global(1-xhalo:nx+xhalo,1-yhalo:ny+yhalo,nz) )
    else
       allocate(global(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )
    end if

    global = 0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             global(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    if(index(type, 'symmetry') == 0) then
       is_symmetry = .false.
    else
       is_symmetry = .true.
    end if
    select case(type)
    case( 'Simple', 'Simple symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, name=type, symmetry = is_symmetry )
    case( 'Cyclic', 'Cyclic symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,        &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry )
        global(1-whalo:0,                 1:ny,:) = global(nx-whalo+1:nx,             1:ny,:)
        global(nx+1:nx+ehalo,             1:ny,:) = global(1:ehalo,                   1:ny,:)
        global(1-whalo:nx+ehalo,     1-shalo:0,:) = global(1-whalo:nx+ehalo, ny-shalo+1:ny,:)
        global(1-whalo:nx+ehalo, ny+1:ny+nhalo,:) = global(1-whalo:nx+ehalo,       1:nhalo,:)
    case( 'Folded-north', 'Folded-north symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_north_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded-south symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_SOUTH_EDGE, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_south_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded-west symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=FOLD_WEST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_west_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded-east symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=FOLD_EAST_EDGE, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_east_halo(global, 0, 0, 0, 0, 1)
    case( 'Folded xy_halo' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, xhalo=xhalo, yhalo=yhalo,   &
             xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, name=type, symmetry = is_symmetry  )
        global(1-xhalo:0,                1:ny,:) = global(nx-xhalo+1:nx,                   1:ny,:)
        global(nx+1:nx+xhalo,            1:ny,:) = global(1:xhalo,                         1:ny,:)
        global(1-xhalo:nx+xhalo,ny+1:ny+yhalo,:) = global(nx+xhalo:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
    case( 'Masked', 'Masked symmetry' )
!with fold and cyclic, assign to npes+1 and mask out the top-rightdomain
        call mpp_define_layout( (/1,nx,1,ny/), npes+1, layout )
        allocate( maskmap(layout(1),layout(2)) )
        maskmap(:,:) = .TRUE.; maskmap(layout(1),layout(2)) = .FALSE.
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
             maskmap=maskmap, name=type, symmetry = is_symmetry  )
        deallocate(maskmap)
       !we need to zero out the global data on the missing domain.
       !this logic assumes top-right, in an even division
        if( mod(nx,layout(1)).NE.0 .OR. mod(ny,layout(2)).NE.0 )call mpp_error( FATAL, &
             'TEST_MPP_DOMAINS: test for masked domains needs (nx,ny) to divide evenly on npes+1 PEs.' )
        global(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny,:) = 0
        call fill_folded_north_halo(global, 0, 0, 0, 0, 1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( x1(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    x = 0.
    x (is:ie,js:je,:) = global(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x

!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)
    call compare_checksums( x, global(isd:ied,jsd:jed,:), type )

!partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x1(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x1' )
    call compare_checksums( x2(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x2' )
    call compare_checksums( x3(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x3' )
    call compare_checksums( x4(is:ied,js:jed,:), global(is:ied,js:jed,:), type//' partial x4' )

    !--- test vector update for FOLDED and MASKED case.
    if(type == 'Simple' .or. type == 'Simple symmetry' .or. type == 'Cyclic' .or. type == 'Cyclic symmetry') then
       deallocate(x,x1,x2,x3,x4)
       return
    end if

    !------------------------------------------------------------------
    !              vector update : BGRID_NE
    !------------------------------------------------------------------
    shift = 0
    if(is_symmetry) then
       shift = 1
       deallocate(global)
       allocate(global(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
       global = 0.0
       do k = 1,nz
          do j = 1,ny+1
             do i = 1,nx+1
                global(i,j,k) = k + i*1e-3 + j*1e-6
             end do
          end do
       end do
       if(type == 'Masked symmetry') then
           global(nx-nx/layout(1)+1:nx+1,ny-ny/layout(2)+1:ny+1,:) = 0
       endif
       deallocate(x, x1, x2, x3, x4)
       allocate( x (isd:ied+1,jsd:jed+1,nz) )
       allocate( x1(isd:ied+1,jsd:jed+1,nz) )
       allocate( x2(isd:ied+1,jsd:jed+1,nz) )
       allocate( x3(isd:ied+1,jsd:jed+1,nz) )
       allocate( x4(isd:ied+1,jsd:jed+1,nz) )
    endif

    folded_south = .false.
    folded_west  = .false.
    folded_east  = .false.
    select case (type)
    case ('Folded-north', 'Masked')
       !fill in folded north edge, cyclic east and west edge
       call fill_folded_north_halo(global, 1, 1, 0, 0, -1)
    case ('Folded xy_halo')
       !fill in folded north edge, cyclic east and west edge
       global(1-xhalo:0,                  1:ny,:) =  global(nx-xhalo+1:nx,                     1:ny,:)
       global(nx+1:nx+xhalo,              1:ny,:) =  global(1:xhalo,                           1:ny,:)
       global(1-xhalo:nx+xhalo-1,ny+1:ny+yhalo,:) = -global(nx+xhalo-1:1-xhalo:-1,ny-1:ny-yhalo:-1,:)
       global(nx+xhalo,          ny+1:ny+yhalo,:) = -global(nx-xhalo,             ny-1:ny-yhalo:-1,:)
    case ('Folded-north symmetry', 'Masked symmetry' )
       call fill_folded_north_halo(global, 1, 1, 1, 1, -1)
    case ('Folded-south symmetry' )
       folded_south = .true.
       call fill_folded_south_halo(global, 1, 1, 1, 1, -1)
    case ('Folded-west symmetry' )
       folded_west = .true.
       call fill_folded_west_halo(global, 1, 1, 1, 1, -1)
    case ('Folded-east symmetry' )
       folded_east = .true.
       call fill_folded_east_halo(global, 1, 1, 1, 1, -1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    x = 0.
    x(is:ie+shift,js:je+shift,:) = global(is:ie+shift,js:je+shift,:)
    !set up y array
    allocate( y (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y1(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y2(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y3(isd:ied+shift,jsd:jed+shift,nz) )
    allocate( y4(isd:ied+shift,jsd:jed+shift,nz) )
    y = x; x1 = x; x2 = x; x3 = x; x4 = x
    y = x; y1 = x; y2 = x; y3 = x; y4 = x

    id = mpp_clock_id( type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=BGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    !redundant points must be equal and opposite

    if(folded_south) then
       global(nx/2+shift,                1,:) = 0.  !pole points must have 0 velocity
       global(nx+shift  ,                1,:) = 0.  !pole points must have 0 velocity
       global(nx/2+1+shift:nx-1+shift,   1,:) = -global(nx/2-1+shift:1+shift:-1, 1,:)
       global(1-whalo:shift,             1,:) = -global(nx-whalo+1:nx+shift,     1,:)
       global(nx+1+shift:nx+ehalo+shift, 1,:) = -global(1+shift:ehalo+shift,     1,:)
       !--- the following will fix the +0/-0 problem on altix
       if(shalo >0) global(shift,1,:) = 0.  !pole points must have 0 velocity
    else if(folded_west) then
       global(1, ny/2+shift, :) = 0. !pole points must have 0 velocity
       global(1, ny+shift,   :) = 0. !pole points must have 0 velocity
       global(1, ny/2+1+shift:ny-1+shift,   :) = -global(1, ny/2-1+shift:1+shift:-1, :)
       global(1, 1-shalo:shift,             :) = -global(1, ny-shalo+1:ny+shift,     :)
       global(1, ny+1+shift:ny+nhalo+shift, :) = -global(1, 1+shift:nhalo+shift,     :)
       !--- the following will fix the +0/-0 problem on altix
       if(whalo>0) global(1, shift, :) = 0.  !pole points must have 0 velocity
    else if(folded_east) then
       global(nx+shift, ny/2+shift, :) = 0. !pole points must have 0 velocity
       global(nx+shift, ny+shift,   :) = 0. !pole points must have 0 velocity
       global(nx+shift, ny/2+1+shift:ny-1+shift,   :) = -global(nx+shift, ny/2-1+shift:1+shift:-1, :)
       global(nx+shift, 1-shalo:shift,             :) = -global(nx+shift, ny-shalo+1:ny+shift,     :)
       global(nx+shift, ny+1+shift:ny+nhalo+shift, :) = -global(nx+shift, 1+shift:nhalo+shift,     :)
       if(ehalo >0) global(nx+shift, shift, :) = 0.  !pole points must have 0 velocity
    else
       global(nx/2+shift,                ny+shift,:) = 0.  !pole points must have 0 velocity
       global(nx+shift  ,                ny+shift,:) = 0.  !pole points must have 0 velocity
       global(nx/2+1+shift:nx-1+shift,   ny+shift,:) = -global(nx/2-1+shift:1+shift:-1, ny+shift,:)
       if(type == 'Folded xy_halo') then
          global(1-xhalo:shift,             ny+shift,:) = -global(nx-xhalo+1:nx+shift,     ny+shift,:)
          global(nx+1+shift:nx+xhalo+shift, ny+shift,:) = -global(1+shift:xhalo+shift,     ny+shift,:)
       else
          global(1-whalo:shift,             ny+shift,:) = -global(nx-whalo+1:nx+shift,     ny+shift,:)
          global(nx+1+shift:nx+ehalo+shift, ny+shift,:) = -global(1+shift:ehalo+shift,     ny+shift,:)
       end if
       !--- the following will fix the +0/-0 problem on altix
       if(nhalo >0) global(shift,ny+shift,:) = 0.  !pole points must have 0 velocity
    endif

    call compare_checksums( x,  global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X' )
    call compare_checksums( y,  global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y' )
    call compare_checksums( x1, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X1' )
    call compare_checksums( x2, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X2' )
    call compare_checksums( x3, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X3' )
    call compare_checksums( x4, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE X4' )
    call compare_checksums( y1, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y1' )
    call compare_checksums( y2, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y2' )
    call compare_checksums( y3, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y3' )
    call compare_checksums( y4, global(isd:ied+shift,jsd:jed+shift,:), type//' BGRID_NE Y4' )

    deallocate(global, x, x1, x2, x3, x4, y, y1, y2, y3, y4)

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- global1 is x-component and global2 is y-component
    if(type == 'Folded xy_halo') then
       allocate(global1(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
       allocate(global2(1-xhalo:nx+xhalo, 1-yhalo:ny+yhalo, nz))
    else
       allocate(global1(1-whalo:nx+ehalo+shift, 1-shalo:ny+nhalo, nz))
       allocate(global2(1-whalo:nx+ehalo, 1-shalo:ny+nhalo+shift, nz))
    end if
    allocate(x (isd:ied+shift,jsd:jed,nz), y (isd:ied,jsd:jed+shift,nz) )
    allocate(x1(isd:ied+shift,jsd:jed,nz), y1(isd:ied,jsd:jed+shift,nz) )
    allocate(x2(isd:ied+shift,jsd:jed,nz), y2(isd:ied,jsd:jed+shift,nz) )
    allocate(x3(isd:ied+shift,jsd:jed,nz), y3(isd:ied,jsd:jed+shift,nz) )
    allocate(x4(isd:ied+shift,jsd:jed,nz), y4(isd:ied,jsd:jed+shift,nz) )

    global1 = 0.0
    global2 = 0.0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx+shift
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
       do j = 1,ny+shift
          do i = 1,nx
             global2(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    if(type == 'Masked' .or. type == 'Masked symmetry') then
       global1(nx-nx/layout(1)+1:nx+shift,ny-ny/layout(2)+1:ny,:) = 0
       global2(nx-nx/layout(1)+1:nx,ny-ny/layout(2)+1:ny+shift,:) = 0
    end if

    select case (type)
    case ('Folded-north', 'Masked')
       !fill in folded north edge, cyclic east and west edge
       call fill_folded_north_halo(global1, 1, 0, 0, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 0, -1)
    case ('Folded xy_halo')
       global1(1-xhalo:0,                   1:ny,:) =  global1(nx-xhalo+1:nx,                     1:ny,:)
       global1(nx+1:nx+xhalo,               1:ny,:) =  global1(1:xhalo,                           1:ny,:)
       global2(1-xhalo:0,                   1:ny,:) =  global2(nx-xhalo+1:nx,                     1:ny,:)
       global2(nx+1:nx+xhalo,               1:ny,:) =  global2(1:xhalo,                           1:ny,:)
       global1(1-xhalo:nx+xhalo-1, ny+1:ny+yhalo,:) = -global1(nx+xhalo-1:1-xhalo:-1, ny:ny-yhalo+1:-1,:)
       global1(nx+xhalo,           ny+1:ny+yhalo,:) = -global1(nx-xhalo,              ny:ny-yhalo+1:-1,:)
       global2(1-xhalo:nx+xhalo,   ny+1:ny+yhalo,:) = -global2(nx+xhalo:1-xhalo:-1,   ny-1:ny-yhalo:-1,:)
    case ('Folded-north symmetry')
       call fill_folded_north_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 1, -1)
    case ('Folded-south symmetry')
       call fill_folded_south_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_south_halo(global2, 0, 1, 0, 1, -1)
    case ('Folded-west symmetry')
       call fill_folded_west_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_west_halo(global2, 0, 1, 0, 1, -1)
    case ('Folded-east symmetry')
       call fill_folded_east_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_east_halo(global2, 0, 1, 0, 1, -1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    x = 0.; y = 0.
    x(is:ie+shift,js:je,      :) = global1(is:ie+shift,js:je,      :)
    y(is:ie      ,js:je+shift,:) = global2(is:ie,      js:je+shift,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    !redundant points must be equal and opposite
    if(folded_south) then
       global2(nx/2+1:nx,     1,:) = -global2(nx/2:1:-1, 1,:)
       global2(1-whalo:0,     1,:) = -global2(nx-whalo+1:nx, 1, :)
       global2(nx+1:nx+ehalo, 1,:) = -global2(1:ehalo,       1, :)
    else if(folded_west) then
       global1(1, ny/2+1:ny,     :) = -global1(1, ny/2:1:-1,     :)
       global1(1, 1-shalo:0,     :) = -global1(1, ny-shalo+1:ny, :)
       global1(1, ny+1:ny+nhalo, :) = -global1(1, 1:nhalo,       :)
    else if(folded_east) then
       global1(nx+shift, ny/2+1:ny,     :) = -global1(nx+shift, ny/2:1:-1,     :)
       global1(nx+shift, 1-shalo:0,     :) = -global1(nx+shift, ny-shalo+1:ny, :)
       global1(nx+shift, ny+1:ny+nhalo, :) = -global1(nx+shift, 1:nhalo,       :)
    else
       global2(nx/2+1:nx,     ny+shift,:) = -global2(nx/2:1:-1, ny+shift,:)
       if(type == 'Folded xy_halo') then
          global2(1-xhalo:0,     ny+shift,:) = -global2(nx-xhalo+1:nx, ny+shift,:)
          global2(nx+1:nx+xhalo, ny+shift,:) = -global2(1:xhalo,       ny+shift,:)
       else
          global2(1-whalo:0,     ny+shift,:) = -global2(nx-whalo+1:nx, ny+shift,:)
          global2(nx+1:nx+ehalo, ny+shift,:) = -global2(1:ehalo,       ny+shift,:)
       end if
    endif

    call compare_checksums( x,  global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied+shift,jsd:jed,      :), type//' CGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,      jsd:jed+shift,:), type//' CGRID_NE Y4' )

    deallocate(global1, global2, x, x1, x2, x3, x4, y, y1, y2, y3, y4)


  end subroutine test_halo_update

  subroutine set_corner_zero( data, isd, ied, jsd, jed, isc, iec, jsc, jec )
     integer,                               intent(in) :: isd, ied, jsd, jed
     integer,                               intent(in) :: isc, iec, jsc, jec
     real, dimension(isd:,jsd:,:), intent(inout) :: data

    data (isd  :isc-1, jsd  :jsc-1,:) = 0
    data (isd  :isc-1, jec+1:jed,  :) = 0
    data (iec+1:ied  , jsd  :jsc-1,:) = 0
    data (iec+1:ied  , jec+1:jed,  :) = 0


  end subroutine set_corner_zero

  !##################################################################################
  subroutine test_update_edge( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x2, a
    real, allocatable, dimension(:,:,:) :: y, y2, b
    type(domain2D) :: domain
    real,    allocatable :: global1(:,:,:), global2(:,:,:), global(:,:,:)
    logical, allocatable :: maskmap(:,:)
    integer              :: shift, i, xhalo, yhalo
    logical              :: is_symmetry, folded_south, folded_west, folded_east
    integer              :: is, ie, js, je, isd, ied, jsd, jed
    integer              :: id_update

    allocate(global(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz) )

    global = 0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             global(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    if(index(type, 'symmetry') == 0) then
       is_symmetry = .false.
    else
       is_symmetry = .true.
    end if
    select case(type)
    case( 'Cyclic' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,        &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry )
        global(1-whalo:0,          1:ny,:) = global(nx-whalo+1:nx,             1:ny,:)
        global(nx+1:nx+ehalo,      1:ny,:) = global(1:ehalo,                   1:ny,:)
        global(1:nx,          1-shalo:0,:) = global(1:nx,             ny-shalo+1:ny,:)
        global(1:nx,      ny+1:ny+nhalo,:) = global(1:nx,             1:nhalo,      :)
    case( 'Folded-north', 'Folded-north symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
             name=type, symmetry = is_symmetry  )
        call fill_folded_north_halo(global, 0, 0, 0, 0, 1)
        !--- set the corner to 0
        call set_corner_zero(global, 1-whalo, nx+ehalo, 1-shalo, ny+ehalo, 1, nx, 1, ny)
    case default
        call mpp_error( FATAL, 'test_update_edge: no such test: '//type )
    end select

!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( a (isd:ied,jsd:jed,nz) )
    allocate( x2 (isd:ied,jsd:jed,nz) )
    x2 (isd:ied,jsd:jed,:) = global(isd:ied,jsd:jed,:)
    call set_corner_zero(x2, isd, ied, jsd, jed, is, ie, js, je)

    x = 0.
    x (is:ie,js:je,:) = global(is:ie,js:je,:)

!full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain, flags=EDGEUPDATE)
    call mpp_clock_end  (id)
    call compare_checksums( x, x2, type )
    deallocate(x2)

    a = 0
    a(is:ie,js:je,:) = global(is:ie,js:je,:)
    id_update = mpp_start_update_domains( a, domain, flags=EDGEUPDATE)
    call mpp_complete_update_domains(id_update, a, domain, flags=EDGEUPDATE)
    call compare_checksums( x, a, type//" nonblock")

        !--- test vector update for FOLDED and MASKED case.
    if( type == 'Cyclic' ) then
       deallocate(global, x, a)
       return
    end if

    !------------------------------------------------------------------
    !              vector update : BGRID_NE
    !------------------------------------------------------------------
    shift = 0
    if(is_symmetry) then
       shift = 1
       deallocate(global)
       allocate(global(1-whalo:nx+ehalo+shift,1-shalo:ny+nhalo+shift,nz) )
       global = 0.0
       do k = 1,nz
          do j = 1,ny+1
             do i = 1,nx+1
                global(i,j,k) = k + i*1e-3 + j*1e-6
             end do
          end do
       end do
       deallocate(x,a)
       allocate( x (isd:ied+1,jsd:jed+1,nz) )
       allocate( a (isd:ied+1,jsd:jed+1,nz) )
    endif

    select case (type)
    case ('Folded-north')
       !fill in folded north edge, cyclic east and west edge
       call fill_folded_north_halo(global, 1, 1, 0, 0, -1)
    case ('Folded-north symmetry')
       call fill_folded_north_halo(global, 1, 1, 1, 1, -1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    x = 0.
    a = 0.
    x(is:ie+shift,js:je+shift,:) = global(is:ie+shift,js:je+shift,:)
    a(is:ie+shift,js:je+shift,:) = global(is:ie+shift,js:je+shift,:)
    !set up y array
    allocate( y (isd:ied+shift,jsd:jed+shift,nz) )
    allocate( b (isd:ied+shift,jsd:jed+shift,nz) )
    b = x
    y = x
    id = mpp_clock_id( type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y, domain, flags=EDGEUPDATE, gridtype=BGRID_NE)
    call mpp_clock_end  (id)

    !--nonblocking update
    id_update = mpp_start_update_domains(a,b, domain, flags=EDGEUPDATE, gridtype=BGRID_NE)
    call mpp_complete_update_domains(id_update, a,b, domain, flags=EDGEUPDATE, gridtype=BGRID_NE)

    !redundant points must be equal and opposite


    global(nx/2+shift,                ny+shift,:) = 0.  !pole points must have 0 velocity
    global(nx+shift  ,                ny+shift,:) = 0.  !pole points must have 0 velocity
    global(nx/2+1+shift:nx-1+shift,   ny+shift,:) = -global(nx/2-1+shift:1+shift:-1, ny+shift,:)

    global(1-whalo:shift,             ny+shift,:) = -global(nx-whalo+1:nx+shift,     ny+shift,:)
    global(nx+1+shift:nx+ehalo+shift, ny+shift,:) = -global(1+shift:ehalo+shift,     ny+shift,:)
    !--- the following will fix the +0/-0 problem on altix
    if(nhalo >0) global(shift,ny+shift,:) = 0.  !pole points must have 0 velocity

    allocate( x2 (isd:ied+shift,jsd:jed+shift,nz) )
    x2 (isd:ied+shift,jsd:jed+shift,:) = global(isd:ied+shift,jsd:jed+shift,:)
    call set_corner_zero(x2, isd, ied+shift, jsd, jed+shift, is, ie+shift, js, je+shift)

    call compare_checksums( x,  x2, type//' BGRID_NE X' )
    call compare_checksums( y,  x2, type//' BGRID_NE Y' )
    call compare_checksums( a,  x2, type//' BGRID_NE X nonblock' )
    call compare_checksums( b,  x2, type//' BGRID_NE Y nonblock' )

    deallocate(global, x, y, x2, a, b)

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- global1 is x-component and global2 is y-component
    allocate(global1(1-whalo:nx+ehalo+shift, 1-shalo:ny+nhalo, nz))
    allocate(global2(1-whalo:nx+ehalo, 1-shalo:ny+nhalo+shift, nz))
    allocate(x  (isd:ied+shift,jsd:jed,nz), y (isd:ied,jsd:jed+shift,nz) )
    allocate(x2 (isd:ied+shift,jsd:jed,nz), y2 (isd:ied,jsd:jed+shift,nz) )
    allocate(a  (isd:ied+shift,jsd:jed,nz), b (isd:ied,jsd:jed+shift,nz) )

    global1 = 0.0
    global2 = 0.0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx+shift
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
       do j = 1,ny+shift
          do i = 1,nx
             global2(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    select case (type)
    case ('Folded-north')
       !fill in folded north edge, cyclic east and west edge
       call fill_folded_north_halo(global1, 1, 0, 0, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 0, -1)
       !--- set the corner to 0
       global1(1-whalo:0,     1-shalo:0,     :) = 0
       global1(1-whalo:0,     ny+1:ny+nhalo, :) = 0
       global1(nx+1:nx+ehalo, 1-shalo:0,     :) = 0
       global1(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = 0
       global2(1-whalo:0,     1-shalo:0,     :) = 0
       global2(1-whalo:0,     ny+1:ny+nhalo, :) = 0
       global2(nx+1:nx+ehalo, 1-shalo:0,     :) = 0
       global2(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = 0
    case ('Folded-north symmetry')
       call fill_folded_north_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 1, -1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    x = 0.; y = 0.
    x(is:ie+shift,js:je,      :) = global1(is:ie+shift,js:je,      :)
    y(is:ie      ,js:je+shift,:) = global2(is:ie,      js:je+shift,:)
    a = 0.; b = 0.
    a(is:ie+shift,js:je,      :) = global1(is:ie+shift,js:je,      :)
    b(is:ie      ,js:je+shift,:) = global2(is:ie,      js:je+shift,:)

    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y, domain, flags=EDGEUPDATE, gridtype=CGRID_NE)
    call mpp_clock_end  (id)

    !--nonblocking
    id_update = mpp_start_update_domains( a,  b, domain, flags=EDGEUPDATE, gridtype=CGRID_NE)
    call mpp_complete_update_domains(id_update, a,  b, domain, flags=EDGEUPDATE, gridtype=CGRID_NE)

    !redundant points must be equal and opposite
    global2(nx/2+1:nx,     ny+shift,:) = -global2(nx/2:1:-1, ny+shift,:)
    global2(1-whalo:0,     ny+shift,:) = -global2(nx-whalo+1:nx, ny+shift,:)
    global2(nx+1:nx+ehalo, ny+shift,:) = -global2(1:ehalo,       ny+shift,:)

    x2(isd:ied+shift,jsd:jed,:) = global1(isd:ied+shift,jsd:jed,:)
    y2(isd:ied,jsd:jed+shift,:) = global2(isd:ied,jsd:jed+shift,:)
    call set_corner_zero(x2, isd, ied+shift, jsd, jed, is, ie+shift, js, je)
    call set_corner_zero(y2, isd, ied, jsd, jed+shift, is, ie, js, je+shift)

    call compare_checksums( x,  x2, type//' CGRID_NE X' )
    call compare_checksums( y,  y2, type//' CGRID_NE Y' )
    call compare_checksums( a,  x2, type//' CGRID_NE X nonblock' )
    call compare_checksums( b,  y2, type//' CGRID_NE Y nonblock' )

    deallocate(global1, global2, x, y, x2, y2, a, b)


  end subroutine test_update_edge


  !##################################################################################
  subroutine test_update_nonsym_edge( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x2
    real, allocatable, dimension(:,:,:) :: y, y2
    type(domain2D) :: domain
    real,    allocatable :: global1(:,:,:), global2(:,:,:)
    integer              :: shift, i, xhalo, yhalo
    logical              :: is_symmetry
    integer              :: is, ie, js, je, isd, ied, jsd, jed
    type(mpp_group_update_type) :: group_update

    if(index(type, 'symmetry') == 0) then
       shift = 0
       is_symmetry = .false.
    else
       shift = 1
       is_symmetry = .true.
    end if
    select case(type)
    case( 'Folded-north', 'Folded-north symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,   &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE, &
             name=type, symmetry = is_symmetry  )
    case default
        call mpp_error( FATAL, 'test_update_edge: no such test: '//type )
    end select

    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------
    !--- global1 is x-component and global2 is y-component
    allocate(global1(1-whalo:nx+ehalo+shift, 1-shalo:ny+nhalo, nz))
    allocate(global2(1-whalo:nx+ehalo, 1-shalo:ny+nhalo+shift, nz))
    allocate(x  (isd:ied+shift,jsd:jed,nz), y (isd:ied,jsd:jed+shift,nz) )
    allocate(x2 (isd:ied+shift,jsd:jed,nz), y2 (isd:ied,jsd:jed+shift,nz) )

    global1 = 0.0
    global2 = 0.0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx+shift
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
       do j = 1,ny+shift
          do i = 1,nx
             global2(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    select case (type)
    case ('Folded-north')
       !fill in folded north edge, cyclic east and west edge
       call fill_folded_north_halo(global1, 1, 0, 0, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 0, -1)
       !--- set the corner to 0
       global1(1-whalo:0,     1-shalo:0,     :) = 0
       global1(1-whalo:0,     ny+1:ny+nhalo, :) = 0
       global1(nx+1:nx+ehalo, 1-shalo:0,     :) = 0
       global1(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = 0
       global2(1-whalo:0,     1-shalo:0,     :) = 0
       global2(1-whalo:0,     ny+1:ny+nhalo, :) = 0
       global2(nx+1:nx+ehalo, 1-shalo:0,     :) = 0
       global2(nx+1:nx+ehalo, ny+1:ny+nhalo, :) = 0
    case ('Folded-north symmetry')
       call fill_folded_north_halo(global1, 1, 0, 1, 0, -1)
       call fill_folded_north_halo(global2, 0, 1, 0, 1, -1)
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    !redundant points must be equal and opposite
    global2(nx/2+1:nx,     ny+shift,:) = -global2(nx/2:1:-1, ny+shift,:)
    global2(1-whalo:0,     ny+shift,:) = -global2(nx-whalo+1:nx, ny+shift,:)
!    global2(nx+1:nx+ehalo, ny+shift,:) = -global2(1:ehalo,       ny+shift,:)

    x2 = 0.0; y2 = 0.0
    if(is_symmetry) then
       x2(isd:ie+shift,jsd:je,:) = global1(isd:ie+shift,jsd:je,:)
       y2(isd:ie,jsd:je+shift,:) = global2(isd:ie,jsd:je+shift,:)
    else
       x2(isd:ie+shift,js:je,:) = global1(isd:ie+shift,js:je,:)
       y2(is:ie,jsd:je+shift,:) = global2(is:ie,jsd:je+shift,:)
    endif

    x = 0.; y = 0.
    x(is:ie+shift,js:je,      :) = global1(is:ie+shift,js:je,      :)
    y(is:ie      ,js:je+shift,:) = global2(is:ie,      js:je+shift,:)

    call mpp_create_group_update(group_update, x, y, domain, gridtype=CGRID_NE, &
                                 flags=WUPDATE+SUPDATE+NONSYMEDGEUPDATE, whalo=1, ehalo=1, shalo=1, nhalo=1)
    call mpp_do_group_update(group_update, domain, x(is,js,1))

    call compare_checksums( x,  x2, type//' CGRID_NE X' )
    call compare_checksums( y,  y2, type//' CGRID_NE Y' )

    call mpp_sync()

    x = 0.; y = 0.
    x(is:ie+shift,js:je,      :) = global1(is:ie+shift,js:je,      :)
    y(is:ie      ,js:je+shift,:) = global2(is:ie,      js:je+shift,:)
    call mpp_start_group_update(group_update, domain, x(is,js,1))
    call mpp_complete_group_update(group_update, domain, x(is,js,1))

    call compare_checksums( x,  x2, type//' CGRID_NE X nonblock' )
    call compare_checksums( y,  y2, type//' CGRID_NE Y nonblock' )

    deallocate(global1, global2, x, y, x2, y2)
    call mpp_clear_group_update(group_update)

  end subroutine test_update_nonsym_edge


  !##################################################################################
  subroutine test_cyclic_offset( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, x1, x2, x3, x4
    real, allocatable, dimension(:,:,:) :: y, y1, y2, y3, y4
    type(domain2D) :: domain
    real,    allocatable :: global1(:,:,:), global2(:,:,:), global(:,:,:)
    integer              :: i, j, k, jj, ii
    integer              :: is, ie, js, je, isd, ied, jsd, jed
    character(len=128)   :: type2

    allocate(global(1-whalo:nx+ehalo,1-shalo:ny+nhalo,nz))

    global = 0
    do k = 1,nz
       do j = 1,ny
          do i = 1,nx
             global(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    end do

    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'x_cyclic_offset' )
        write(type2, *)type, ' x_cyclic=', x_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 name=type, x_cyclic_offset = x_cyclic_offset)
        do j = 1, ny
           jj = mod(j + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0,j,:) = global(nx-whalo+1:nx, jj,:) ! West
           jj = mod(j - x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo,j,:) = global(1:ehalo,jj,:)    ! East
        end do
    case( 'y_cyclic_offset' )
        write(type2, *)type, ' y_cyclic = ', y_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, yflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 name=type, y_cyclic_offset = y_cyclic_offset)
        do i = 1, nx
           ii = mod(i + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! South
           ii = mod(i - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! NORTH
        end do
    case( 'torus_x_offset' )
        write(type2, *)type, ' x_cyclic = ', x_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 yflags=CYCLIC_GLOBAL_DOMAIN, name=type,                  &
                                 x_cyclic_offset = x_cyclic_offset)
        do j = 1, ny
           jj = mod(j + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0,j,:) = global(nx-whalo+1:nx, jj,:) ! West
           jj = mod(j - x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo,j,:) = global(1:ehalo,jj,:)    ! East
        end do
        global(1:nx,1-shalo:0,:)     = global(1:nx, ny-shalo+1:ny,:) ! South
        global(1:nx,ny+1:ny+nhalo,:) = global(1:nx, 1:nhalo, :)    ! NORTH

        do j = 1, shalo
           jj = mod(ny-j+1 + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0, 1-j,:) = global(nx-whalo+1:nx, jj, :)  ! Southwest
           jj = mod(ny-j+1-x_cyclic_offset+ny,ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo, 1-j,:) = global(1:ehalo, jj, :)    ! Southeast
        end do
        do j = 1, nhalo
           jj = mod(j + x_cyclic_offset + ny, ny)
           if(jj==0) jj = ny
           global(1-whalo:0, ny+j,:) = global(nx-whalo+1:nx, jj, :)  ! northwest
           jj = mod(j - x_cyclic_offset+ny,ny)
           if(jj==0) jj = ny
           global(nx+1:nx+ehalo, ny+j,:) = global(1:ehalo, jj, :)    ! northeast
        end do

    case( 'torus_y_offset' )
        write(type2, *)type, ' y_cyclic = ', y_cyclic_offset
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN,   &
                                 yflags=CYCLIC_GLOBAL_DOMAIN, name=type,                  &
                                 y_cyclic_offset = y_cyclic_offset)
        do i = 1, nx
           ii = mod(i + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! South
           ii = mod(i - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! NORTH
        end do
        global(1-whalo:0,1:ny,:)     = global(nx-whalo+1:nx, 1:ny,:) ! West
        global(nx+1:nx+ehalo,1:ny,:) = global(1:ehalo, 1:ny, :)      ! East
        do i = 1, whalo
           ii = mod(nx-i+1 + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(1-i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! southwest
           ii = mod(nx-i+1 - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(1-i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! northwest
        end do
        do i = 1, ehalo
           ii = mod(i + y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(nx+i, 1-shalo:0,:) = global(ii, ny-shalo+1:ny,:) ! southeast
           ii = mod(i - y_cyclic_offset + nx, nx)
           if(ii==0) ii = nx
           global(nx+i,ny+1:ny+nhalo,:) = global(ii,1:nhalo,:)    ! northeast
        end do
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

    !set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )
    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( x1(isd:ied,jsd:jed,nz) )
    allocate( x2(isd:ied,jsd:jed,nz) )
    allocate( x3(isd:ied,jsd:jed,nz) )
    allocate( x4(isd:ied,jsd:jed,nz) )
    x = 0.
    x (is:ie,js:je,:) = global(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x

    !full update
    id = mpp_clock_id( type, flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x, domain )
    call mpp_clock_end  (id)
    call compare_checksums( x, global(isd:ied,jsd:jed,:), trim(type2) )

    !partial update
    id = mpp_clock_id( type//' partial', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x1, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x2, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x3, domain, NUPDATE+EUPDATE, complete=.false. )
    call mpp_update_domains( x4, domain, NUPDATE+EUPDATE, complete=.true. )
    call mpp_clock_end  (id)
    call compare_checksums( x1(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x1' )
    call compare_checksums( x2(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x2' )
    call compare_checksums( x3(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x3' )
    call compare_checksums( x4(is:ied,js:jed,:), global(is:ied,js:jed,:), trim(type2)//' partial x4' )

    !--- test vector update for FOLDED and MASKED case.
    deallocate(x,x1,x2,x3,x4)


    !------------------------------------------------------------------
    !              vector update : BGRID_NE
    !------------------------------------------------------------------
    !--- global1 is x-component and global2 is y-component
    allocate(global1(1-whalo:nx+ehalo, 1-shalo:ny+nhalo, nz))
    allocate(global2(1-whalo:nx+ehalo, 1-shalo:ny+nhalo, nz))
    allocate(x (isd:ied,jsd:jed,nz), y (isd:ied,jsd:jed,nz) )
    allocate(x1(isd:ied,jsd:jed,nz), y1(isd:ied,jsd:jed,nz) )
    allocate(x2(isd:ied,jsd:jed,nz), y2(isd:ied,jsd:jed,nz) )
    allocate(x3(isd:ied,jsd:jed,nz), y3(isd:ied,jsd:jed,nz) )
    allocate(x4(isd:ied,jsd:jed,nz), y4(isd:ied,jsd:jed,nz) )
    where (global >0)
       global1 = 1000 + global
       global2 = 2000 + global
    elsewhere
       global1 = 0
       global2 = 0
    end where
    x = 0.; y = 0
    x(is:ie,js:je,:) = global1(is:ie,js:je,:)
    y(is:ie,js:je,:) = global2(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector BGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=BGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=BGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=BGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    !redundant points must be equal and opposite

    call compare_checksums( x,  global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,jsd:jed,:), trim(type2)//' BGRID_NE Y4' )

    !------------------------------------------------------------------
    !              vector update : CGRID_NE
    !------------------------------------------------------------------

    x = 0.; y = 0.
    x(is:ie,js:je,:) = global1(is:ie,js:je,:)
    y(is:ie,js:je,:) = global2(is:ie,js:je,:)
    x1 = x; x2 = x; x3 = x; x4 = x
    y1 = y; y2 = y; y3 = y; y4 = y

    id = mpp_clock_id( type//' vector CGRID_NE', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_update_domains( x,  y,  domain, gridtype=CGRID_NE)
    call mpp_update_domains( x1, y1, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x2, y2, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x3, y3, domain, gridtype=CGRID_NE, complete=.false. )
    call mpp_update_domains( x4, y4, domain, gridtype=CGRID_NE, complete=.true.  )
    call mpp_clock_end  (id)

    call compare_checksums( x,  global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X' )
    call compare_checksums( y,  global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y' )
    call compare_checksums( x1, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X1' )
    call compare_checksums( x2, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X2' )
    call compare_checksums( x3, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X3' )
    call compare_checksums( x4, global1(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE X4' )
    call compare_checksums( y1, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y1' )
    call compare_checksums( y2, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y2' )
    call compare_checksums( y3, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y3' )
    call compare_checksums( y4, global2(isd:ied,jsd:jed,:), trim(type2)//' CGRID_NE Y4' )

    deallocate(global1, global2, x, x1, x2, x3, x4, y, y1, y2, y3, y4)


  end subroutine test_cyclic_offset


  subroutine test_global_field( type )
    character(len=*), intent(in) :: type
    real, allocatable, dimension(:,:,:) :: x, gcheck
    type(domain2D) :: domain
    real, allocatable    :: global1(:,:,:)
    integer              :: ishift, jshift, ni, nj, i, j, position
    integer, allocatable :: pelist(:)
    integer              :: is, ie, js, je, isd, ied, jsd, jed

    !--- set up domain
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
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !--- determine if an extra point is needed
    ishift = 0; jshift = 0
    position = CENTER
    select case(type)
    case ('Symmetry corner')
       ishift = 1; jshift = 1; position=CORNER
    case ('Symmetry east')
       ishift = 1; jshift = 0; position=EAST
    case ('Symmetry north')
       ishift = 0; jshift = 1; position=NORTH
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-whalo:ni+ehalo, 1-shalo:nj+nhalo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    allocate( gcheck(ni, nj, nz) )
    allocate( x (isd:ied,jsd:jed,nz) )

    x(:,:,:) = global1(isd:ied,jsd:jed,:)

    !--- test the data on data domain
    gcheck = 0.
    id = mpp_clock_id( type//' global field on data domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on data domain' )

    !--- Since in the disjoint redistribute mpp test, pelist1 = (npes/2+1 .. npes-1)
    !--- will be declared. But for the x-direction global field, mpp_sync_self will
    !--- be called. For some pe count, pelist1 will be set ( only on pe of pelist1 )
    !--- in the mpp_sync_self call, later when calling mpp_declare_pelist(pelist1),
    !--- deadlock will happen. For example npes = 6 and layout = (2,3), pelist = (4,5)
    !--- will be set in mpp_sync_self. To solve the problem, some explicit mpp_declare_pelist
    !--- on all pe is needed for those partial pelist. But for y-update, it is ok.
    !--- because the pelist in y-update is not continous.
    allocate(pelist(0:layout(1)-1))
    do j = 0, layout(2)-1
       do i = 0, layout(1)-1
          pelist(i) = j*layout(1) + i
       end do
       call mpp_declare_pelist(pelist)
    end do
    deallocate(pelist)

    !xupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = XUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on data domain' )

    !yupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, flags = YUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on data domain' )

    call mpp_clock_begin(id)
    call mpp_global_field( domain, x, gcheck, position=position )

    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, &
                            type//' mpp_global_field on data domain' )

    !--- test the data on compute domain
    gcheck = 0.
    id = mpp_clock_id( type//' global field on compute domain', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je, :), gcheck, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,1:nj,:), gcheck, type//' mpp_global_field on compute domain' )

    !xupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = XUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(1:ni,js:je,:), gcheck(1:ni,js:je,:), &
                            type//' mpp_global_field xupdate only on compute domain' )

    !yupdate
    gcheck = 0.
    call mpp_clock_begin(id)
    call mpp_global_field( domain, x(is:ie, js:je,:), gcheck, flags = YUPDATE, position=position )
    call mpp_clock_end  (id)
    !compare checksums between global and x arrays
    call compare_checksums( global1(is:ie,1:nj,:), gcheck(is:ie,1:nj,:), &
                            type//' mpp_global_field yupdate only on compute domain' )


    deallocate(global1, gcheck, x)

  end subroutine test_global_field

    !--- test mpp_global_sum, mpp_global_min and mpp_global_max
  subroutine test_global_reduce (type)
    character(len=*), intent(in) :: type
    real    :: lsum, gsum, lmax, gmax, lmin, gmin
    integer :: ni, nj, ishift, jshift, position
    integer              :: is, ie, js, je, isd, ied, jsd, jed

    type(domain2D) :: domain
    real, allocatable, dimension(:,:,:) :: global1, x
    real, allocatable, dimension(:,:)   :: global2D
    !--- set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Simple' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type )
    case( 'Simple symmetry center', 'Simple symmetry corner', 'Simple symmetry east', 'Simple symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case( 'Cyclic symmetry center', 'Cyclic symmetry corner', 'Cyclic symmetry east', 'Cyclic symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                                    name=type, symmetry = .true., xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !--- determine if an extra point is needed
    ishift = 0; jshift = 0; position = CENTER
    select case(type)
    case ('Simple symmetry corner', 'Cyclic symmetry corner')
       ishift = 1; jshift = 1; position = CORNER
    case ('Simple symmetry east', 'Cyclic symmetry east' )
       ishift = 1; jshift = 0; position = EAST
    case ('Simple symmetry north', 'Cyclic symmetry north')
       ishift = 0; jshift = 1; position = NORTH
    end select

    ie  = ie+ishift;  je  = je+jshift
    ied = ied+ishift; jed = jed+jshift
    ni  = nx+ishift;  nj  = ny+jshift
    allocate(global1(1-whalo:ni+ehalo, 1-shalo:nj+nhalo, nz))
    global1 = 0.0
    do k = 1,nz
       do j = 1,nj
          do i = 1,ni
             global1(i,j,k) = k + i*1e-3 + j*1e-6
          end do
       end do
    enddo

    !--- NOTE: even though the domain is cyclic, no need to apply cyclic condition on the global data

    allocate( x (isd:ied,jsd:jed,nz) )
    allocate( global2D(ni,nj))

    x(:,:,:) = global1(isd:ied,jsd:jed,:)
    do j = 1, nj
       do i = 1, ni
          global2D(i,j) = sum(global1(i,j,:))
       enddo
    enddo
    !test mpp_global_sum

    if(type(1:6) == 'Simple') then
       gsum = sum( global2D(1:ni,1:nj) )
    else
       gsum = sum( global2D(1:nx, 1:ny) )
    endif
    id = mpp_clock_id( type//' sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x, position = position  )
    call mpp_clock_end  (id)
    if( pe.EQ.mpp_root_pe() )print '(a,2es15.8,a,es12.4)', type//' Fast sum=', lsum, gsum

    !test exact mpp_global_sum
    id = mpp_clock_id( type//' exact sum', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lsum = mpp_global_sum( domain, x, BITWISE_EXACT_SUM, position = position )
    call mpp_clock_end  (id)
    !--- The following check will fail on altix in normal mode, but it is ok
    !--- in debugging mode. It is ok on irix.
    call compare_data_scalar(lsum, gsum, FATAL, type//' mpp_global_exact_sum')

    !test mpp_global_min
    gmin = minval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' min', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmin = mpp_global_min( domain, x, position = position )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmin, gmin, FATAL, type//' mpp_global_min')

    !test mpp_global_max
    gmax = maxval(global1(1:ni, 1:nj, :))
    id = mpp_clock_id( type//' max', flags=MPP_CLOCK_SYNC+MPP_CLOCK_DETAILED )
    call mpp_clock_begin(id)
    lmax = mpp_global_max( domain, x, position = position )
    call mpp_clock_end  (id)
    call compare_data_scalar(lmax, gmax, FATAL, type//' mpp_global_max' )

    deallocate(global1, x)

  end subroutine test_global_reduce


  subroutine test_parallel ( )

    integer :: npes, layout(2), i, j, k,is, ie, js, je, isd, ied, jsd, jed
    real, dimension(:,:), allocatable :: field, lfield
    real, dimension(:,:,:), allocatable :: field3d, lfield3d
    type(domain2d) :: domain
    integer, dimension(:), allocatable :: pelist1 , pelist2
    logical :: group1, group2
    character(len=128)  :: mesg

    npes = mpp_npes()
    allocate(pelist1(npes-mpes), pelist2(mpes))
    pelist1 = (/(i, i = 0, npes-mpes -1)/)
    pelist2 = (/(i, i = npes-mpes, npes - 1)/)
    call mpp_declare_pelist(pelist1)
    call mpp_declare_pelist(pelist2)
    group1 = .FALSE. ; group2 = .FALSE.
    if(any(pelist1==pe)) group1 = .TRUE.
    if(any(pelist2==pe)) group2 = .TRUE.
    mesg = 'parallel checking'

    if(group1) then
       call mpp_set_current_pelist(pelist1)
       call mpp_define_layout( (/1,nx,1,ny/), npes-mpes, layout )
    else if(group2) then
       call mpp_set_current_pelist(pelist2)
       call mpp_define_layout( (/1,nx,1,ny/), mpes, layout )
    endif
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo)

    call mpp_set_current_pelist()

     call mpp_get_compute_domain(domain, is, ie, js, je)
     call mpp_get_data_domain(domain, isd, ied, jsd, jed)
     allocate(lfield(is:ie,js:je),field(isd:ied,jsd:jed))
     allocate(lfield3d(is:ie,js:je,nz),field3d(isd:ied,jsd:jed,nz))

     do i = is, ie
     do j = js, je
        lfield(i,j) = real(i)+real(j)*0.001
     enddo
     enddo
     do i = is, ie
     do j = js, je
     do k = 1, nz
        lfield3d(i,j,k) = real(i)+real(j)*0.001+real(k)*0.00001
     enddo
     enddo
     enddo
     field = 0.0
     field3d = 0.0
     field(is:ie,js:je)= lfield(is:ie,js:je)
     field3d(is:ie,js:je,:) = lfield3d(is:ie,js:je,:)
     call mpp_update_domains(field,domain)
     call mpp_update_domains(field3d,domain)

    call mpp_check_field(field, pelist1, pelist2,domain, '2D '//mesg, w_halo = whalo, &
                            s_halo = shalo, e_halo = ehalo, n_halo = nhalo)
    call mpp_check_field(field3d, pelist1, pelist2,domain, '3D '//mesg, w_halo = whalo, &
                            s_halo = shalo, e_halo = ehalo, n_halo = nhalo)

  end subroutine test_parallel

  subroutine test_modify_domain( )

    type(domain2D) :: domain2d_no_halo, domain2d_with_halo
    integer :: is1, ie1, js1, je1, isd1, ied1, jsd1, jed1
    integer :: is2, ie2, js2, je2, isd2, ied2, jsd2, jed2

    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    call mpp_define_domains( (/1,nx,1,ny/), layout, domain2d_no_halo,   &
                            yflags=CYCLIC_GLOBAL_DOMAIN, xhalo=0, yhalo=0)

    call mpp_get_compute_domain(domain2d_no_halo, is1, ie1, js1, je1)
    call mpp_get_data_domain(domain2d_no_halo, isd1, ied1, jsd1, jed1)
    call mpp_modify_domain(domain2d_no_halo, domain2d_with_halo, whalo=whalo,ehalo=ehalo,shalo=shalo,nhalo=nhalo)
    call mpp_get_compute_domain(domain2d_with_halo, is2, ie2, js2, je2)
    call mpp_get_data_domain(domain2d_with_halo, isd2, ied2, jsd2, jed2)
    if( is1 .NE. is2 .OR. ie1 .NE. ie2 .OR. js1 .NE. js2 .OR. je1 .NE. je2 ) then
        print*, "at pe ", pe, " compute domain without halo: ", is1, ie1, js1, je1, &
                " is not equal to the domain with halo ", is2, ie2, js2, je2
        call mpp_error(FATAL, "compute domain mismatch between domain without halo and domain with halo")
    end if

    if( isd1-whalo .NE. isd2 .OR. ied1+ehalo .NE. ied2 .OR. jsd1-shalo .NE. jsd2 .OR. jed1+nhalo .NE. jed2 ) then
        print*, "at pe ", pe, "halo is w=",whalo,",e=",ehalo,",s=",shalo,"n=",nhalo, &
               ",data domain without halo is ",isd1, ied1, jsd1, jed1,                     &
               ", data domain with halo is ", isd2, ied2, jsd2, jed2
        call mpp_error(FATAL, "compute domain mismatch between data domain without halo and data domain with halo")
    else
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, 'test_modify_domain: OK.' )
    end if

    return

  end subroutine test_modify_domain

  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(LONG_KIND) :: sum1, sum2
    integer :: i, j, k

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    call mpp_sync_self()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_chksum: size of a and b does not match')

    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
                write(*,'(a,i3,a,i3,a,i3,a,i3,a,f20.9,a,f20.9)') trim(string)//" at pe ", mpp_pe(), &
                     ", at point (",i,", ", j, ", ", k, "), a = ", a(i,j,k), ", b = ", b(i,j,k)
                call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
             endif
          enddo
       enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
    else
        write(stdunit, *)"sum1 =", sum1, mpp_pe()
        write(stdunit, *)"sum2 =", sum2, mpp_pe()
        write(stdunit,'(a,i3,a,i20,a,i20)')" at pe ", mpp_pe(), " sum(a)=", sum1, " sum(b)=", sum2
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums

  !###########################################################################
  subroutine compare_checksums_2D( a, b, string )
    real, intent(in), dimension(:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(LONG_KIND) :: sum1, sum2
    integer :: i, j

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    ! mpp_sync()
    call mpp_sync_self()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) ) &
         call mpp_error(FATAL,'compare_chksum_2D: size of a and b does not match')

    do j = 1, size(a,2)
       do i = 1, size(a,1)
          if(a(i,j) .ne. b(i,j)) then
            print*, "a =", a(i,j)
            print*, "b =", b(i,j)
             write(*,'(a,i3,a,i3,a,i3,a,f20.9,a,f20.9)')"at the pe ", mpp_pe(), &
                  ", at point (",i,", ", j, "),a=", a(i,j), ",b=", b(i,j)
             call mpp_error(FATAL, trim(string)//': point by point comparison are not OK.')
          endif
       enddo
    enddo

    sum1 = mpp_chksum( a, (/pe/) )
    sum2 = mpp_chksum( b, (/pe/) )

    if( sum1.EQ.sum2 )then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': OK.' )
        !--- in some case, even though checksum agree, the two arrays
        !    actually are different, like comparing (1.1,-1.2) with (-1.1,1.2)
        !--- hence we need to check the value point by point.
    else
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums_2D


  !###########################################################################

  subroutine compare_data_scalar( a, b, action, string )
    real,             intent(in) :: a, b
    integer,          intent(in) :: action
    character(len=*), intent(in) :: string
    if( a .EQ. b)then
        if( pe.EQ.mpp_root_pe() )call mpp_error( NOTE, trim(string)//': data comparison are OK.' )
    else
        write(stdunit,'(a,i3,a,es12.4,a,es12.4,a,es12.4)')' on pe ', mpp_pe(),' a = ', a, ', b = ', b, ', a - b =', a-b
        call mpp_error( action, trim(string)//': data comparison are not OK.' )
    end if

  end subroutine compare_data_scalar

  subroutine test_get_neighbor_1d
    type(domain1d) :: dmn1d
    integer npes, peN, peS
    npes = mpp_npes()
    call mpp_define_domains((/1,npes/), npes, dmn1d)
    call mpp_get_neighbor_pe(dmn1d, direction=+1, pe=peN)
    call mpp_get_neighbor_pe(dmn1d, direction=-1, pe=peS)
    print '(a,i2,a,2i3)', 'PE: ', mpp_pe(), ' R/L pes: ', peN, peS
  end subroutine test_get_neighbor_1d

  subroutine test_get_neighbor_non_cyclic
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_non_cyclic '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (non-cyclic): ', layout,  &
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_non_cyclic

  subroutine test_get_neighbor_cyclic
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_cyclic '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (cyclic)    : ', layout, &
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_cyclic

  subroutine test_get_neighbor_folded_north
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()
    if( npes .NE. 8 ) then
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_folded_north '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes, layout )
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         xflags=CYCLIC_GLOBAL_DOMAIN, yflags=FOLD_NORTH_EDGE)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i2,a,2i2,a,8i3)','PE: ', mpp_pe(), ' layout (folded N)  : ', layout, &
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_folded_north

  subroutine test_get_neighbor_mask
    logical, allocatable ::  mask(:,:)
    integer :: im, jm, n_remove
    type(domain2d) :: domain
    integer nx, ny,layout(2), halo, peN, peS, peE, peW, peNE, peNW, peSE, peSW, npes
    nx = 10
    ny = 20
    halo = 2
    npes = mpp_npes()

    n_remove = 2
    if( npes .NE. 8 ) then
       call mpp_error(NOTE, 'test_mpp_domains: test_get_neighbor_mask '// &
                            ' will be performed only when npes = 8')
      return
    end if
    call mpp_define_layout( (/1,nx, 1,ny/), npes+n_remove, layout )
    allocate(mask(layout(1), layout(2)))
    mask = .TRUE.  ! activate domains
    im = min(layout(1), ceiling(layout(1)/2.0))
    jm = min(layout(2), ceiling(layout(2)/2.0))
    mask(im  ,jm  ) = .FALSE. ! deactivate domain
    mask(im  ,jm-1) = .FALSE. ! deactivate domain
    print '(a,2i3,a,2i3)', 'Masked out domains ', im, jm, ' and ', im,jm-1
    call mpp_define_domains((/1,nx, 1,ny/), layout, domain, xhalo=halo, yhalo=halo, &
         maskmap=mask)
    call mpp_get_neighbor_pe(domain, direction=NORTH, pe=peN)
    call mpp_get_neighbor_pe(domain, direction=SOUTH, pe=peS)
    call mpp_get_neighbor_pe(domain, direction=EAST, pe=peE)
    call mpp_get_neighbor_pe(domain, direction=WEST, pe=peW)
    call mpp_get_neighbor_pe(domain, direction=NORTH_EAST, pe=peNE)
    call mpp_get_neighbor_pe(domain, direction=NORTH_WEST, pe=peNW)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_EAST, pe=peSE)
    call mpp_get_neighbor_pe(domain, direction=SOUTH_WEST, pe=peSW)
    print '(a,i3,a,2i3,a,8i3)','PE: ', mpp_pe(), ' layout (mask   )  : ', layout, &
         & ' N/S/E/W/NE/SE/SW/NW pes: ', peN, peS, peE, peW, peNE, peSE, peSW, peNW
  end subroutine test_get_neighbor_mask

  subroutine test_define_mosaic_pelist(type, ntile)
    character(len=*),       intent(in) :: type
    integer,                intent(in) :: ntile
    integer                            :: npes, root_pe, start_pe, n, ntile_per_pe
    integer, dimension(:), allocatable :: pe1_start, pe1_end, pe2_start, pe2_end
    integer, dimension(:), allocatable :: sizes, costpertile

    root_pe = mpp_root_pe()
    npes = mpp_npes()

    allocate(sizes(ntile), pe1_start(ntile), pe1_end(ntile), pe2_start(ntile), pe2_end(ntile),costpertile(ntile) )
    costpertile = 1
    sizes = nx*ny
    if(npes ==1) then
       pe1_start = root_pe; pe1_end = root_pe
    end if
    select case(type)
    case('One tile')
       pe1_start = root_pe; pe1_end = npes+root_pe-1
    case('Two uniform tile')
       if(mod(npes,2) .NE. 0 .AND. npes .NE. 1) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 2, no test for '//type )
          return
       end if
       if(npes .NE. 1) then
          pe1_start(1) = root_pe;        pe1_end(1) = npes/2+root_pe-1
          pe1_start(2) = npes/2+root_pe; pe1_end(2) = npes+root_pe-1
       end if
    case('Two nonuniform tile')
       if(mod(npes,3) .NE. 0 .AND. npes .NE. 1) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 3, no test for '//type )
          return
       end if
       sizes(1) = 2*nx*ny
       if(npes .NE. 1) then
          pe1_start(1) = root_pe;          pe1_end(1) = npes/3*2+root_pe-1
          pe1_start(2) = npes/3*2+root_pe; pe1_end(2) = npes+root_pe-1
       end if
    case('Ten tile')
       if(mod(npes,10) .NE. 0 .AND. npes .NE. 1 .AND. mod(10,npes) .NE. 0) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 10(or reverse), no test for '//type )
          return
       end if
       if(mod(10, npes)==0) then
          ntile_per_pe = ntile/npes
          do n = 1, ntile
             pe1_start(n) = root_pe+(n-1)/ntile_per_pe; pe1_end(n) = pe1_start(n)
          end do
       else if(mod(npes,10) == 0) then
          do n = 1, ntile
             pe1_start(n) = npes/10*(n-1)+root_pe; pe1_end(n) = npes/10*n+root_pe-1
          end do
       end if
    case('Ten tile with nonuniform cost')
       if(mod(npes,15) .NE. 0 .AND. npes .NE. 1) then
          call mpp_error(NOTE, 'test_define_mosaic_pelist: npes can not be divided by 15, no test for '//type )
          return
       end if
       costpertile(1:5) = 2; costpertile(6:ntile) = 1
       if(npes .NE. 1) then
          start_pe = root_pe
          do n = 1, ntile
             pe1_start(n) = start_pe
             pe1_end(n)   = start_pe + npes/15*costpertile(n)-1
             start_pe = pe1_end(n) + 1
          end do
       end if
    case default
       call mpp_error(FATAL,"test_define_mosaic_pelist: "//type//" is an invalid type")
    end select

    call mpp_define_mosaic_pelist( sizes, pe2_start, pe2_end, costpertile=costpertile)
    if( ANY(pe1_start .NE. pe2_start) .OR. ANY(pe1_end .NE. pe2_end) ) then
       call mpp_error(FATAL,"test_define_mosaic_pelist: test failed for "//trim(type) )
    else
       call mpp_error(NOTE,"test_define_mosaic_pelist: test successful for "//trim(type) )
    end if

  end subroutine test_define_mosaic_pelist

!###############################################################################
! test halo update for grid nested in global cubic sphere grid. The nested region may cross the edge.
! It is assumed the boundary condition of nested region is solid wall in both direction.
  subroutine test_nest_halo_update( domain )
    type(domain2D), intent(inout) :: domain
    integer        :: i, j, k, shift
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed

    real,    allocatable, dimension(:,:,:) :: x, y
    real,    allocatable, dimension(:,:,:) :: global1, global2, global
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

    call compare_checksums( x (isd:ied+shift,jsd:jed+shift,:),  global1(isd:ied+shift,jsd:jed+shift,:), trim(type)//' X' )
    call compare_checksums( y (isd:ied+shift,jsd:jed+shift,:),  global2(isd:ied+shift,jsd:jed+shift,:), trim(type)//' Y' )

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

  end subroutine test_nest_halo_update

  subroutine get_nnest(domain, num_nest, tile_coarse, istart_coarse, iend_coarse, jstart_coarse, jend_coarse, &
                       nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
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
                       nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
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
  subroutine fill_nest_data(buffer, is, ie, js, je, nnest, tile, ishift, jshift, iadd, jadd, rotate, &
                            isl, iel, jsl, jel, xadd, yadd, sign1, sign2)
     real, dimension(is:,js:,:), intent(inout) :: buffer
     integer,                       intent(in) :: is, ie, js, je, nnest
     integer,                       intent(in) :: ishift, jshift
     integer, dimension(:),         intent(in) :: tile, iadd, jadd, rotate, isl, iel, jsl, jel
     real,                          intent(in) :: xadd, yadd
     integer,                       intent(in) :: sign1, sign2
     integer :: i, j, k, n, nk
     integer :: ioff, joff

     ioff = 0
     joff = 0
     nk = size(buffer,3)
     do k = 1, nk
        do n = 1, nnest
           if(iel(n) == ie) ioff = ishift
           if(jel(n) == je) joff = jshift

           select case (rotate(n))
           case(ZERO)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = xadd + tile(n) + (i-iadd(n))*1.e-3 + (j-jadd(n))*1.e-6 + k*1.e-9
                 enddo
              enddo
           case (NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign2*(yadd + tile(n) + (j-jadd(n))*1.e-3 + (nx-i+iadd(n)+1+ioff)*1.e-6 + k*1.e-9)
                 enddo
              enddo
           case (MINUS_NINETY)
              do j = jsl(n), jel(n)+joff
                 do i = isl(n), iel(n)+ioff
                    buffer(i,j,k) = sign1*(yadd + tile(n) + (ny-j+jadd(n)+1+joff)*1.e-3 + (i-iadd(n))*1.e-6 + k*1.e-9)
                 enddo
              enddo
           case default
              call mpp_error(FATAL,"fill_nest_data: rotate must be ZERO, NINETY, MINUS_NINETY")
           end select
        enddo
     enddo

  end subroutine fill_nest_data

!###############################################################################


  subroutine test_update_nest_domain( type )
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
    real,    allocatable         :: x(:,:,:), x1(:,:,:), x2(:,:,:)
    real,    allocatable         :: y(:,:,:), y1(:,:,:), y2(:,:,:)
    real,    allocatable         :: wbuffer(:,:,:), wbuffer2(:,:,:)
    real,    allocatable         :: ebuffer(:,:,:), ebuffer2(:,:,:)
    real,    allocatable         :: sbuffer(:,:,:), sbuffer2(:,:,:)
    real,    allocatable         :: nbuffer(:,:,:), nbuffer2(:,:,:)
    real,    allocatable         :: wbufferx(:,:,:), wbufferx2(:,:,:)
    real,    allocatable         :: ebufferx(:,:,:), ebufferx2(:,:,:)
    real,    allocatable         :: sbufferx(:,:,:), sbufferx2(:,:,:)
    real,    allocatable         :: nbufferx(:,:,:), nbufferx2(:,:,:)
    real,    allocatable         :: wbuffery(:,:,:), wbuffery2(:,:,:)
    real,    allocatable         :: ebuffery(:,:,:), ebuffery2(:,:,:)
    real,    allocatable         :: sbuffery(:,:,:), sbuffery2(:,:,:)
    real,    allocatable         :: nbuffery(:,:,:), nbuffery2(:,:,:)
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
    case default
       call mpp_error(FATAL, 'test_update_nest_domain: no such test: '//type)
    end select

    if(ntiles_nest_all > MAX_NTILE) call mpp_error(FATAL, 'test_update_nest_domain: ntiles_nest_all > MAX_NTILE')
    if(ntiles_nest_top .GE. ntiles_nest_all) call mpp_error(FATAL, 'test_update_nest_domain: ntiles_nest_top .GE. ntile_nest_all')
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
       if(nest_level(n) > nest_level(n-1)+1) call mpp_error(FATAL, "test_mpp_domains: nest_level(n) > nest_level(n-1)+1")
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

       allocate(layout2D(2,ntiles_nest_top), global_indices(4,ntiles_nest_top), pe_start(ntiles_nest_top), pe_end(ntiles_nest_top) )
       npes_per_tile = npes_nest_tile(1)

       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles_nest_top
          global_indices(:,n) = (/1,nx,1,ny/)
          layout2D(:,n)         = layout
       end do
       do n = 1, ntiles_nest_top
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do

       if( cubic_grid ) then
          call define_cubic_mosaic(type, domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
                                   global_indices, layout2D, pe_start, pe_end )
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
          call mpp_define_layout( (/1,nx_fine,1,ny_fine/), my_npes, layout )
          call mpp_define_domains((/1,nx_fine,1,ny_fine/), layout, domain, &
                          whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                          symmetry=.true., name=trim(type)//' fine grid', tile_id = tile_fine(n) )
          call mpp_get_compute_domain(domain, isc_fine, iec_fine, jsc_fine, jec_fine)
          call mpp_get_data_domain(domain, isd_fine, ied_fine, jsd_fine, jed_fine)
          !--- test halo update for nested region.
          call test_nest_halo_update(domain)
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
             tile_coarse(1:num_nest), istart_coarse(1:num_nest), icount_coarse(1:num_nest), jstart_coarse(1:num_nest), &
             jcount_coarse(1:num_nest), npes_nest_tile(1:ntiles_nest_all), &
             x_refine(1:num_nest), y_refine(1:num_nest), extra_halo=extra_halo, name="nest_domain")

    !--- loop over nest level
    do l = 1, num_nest_level
       npes_my_level = mpp_get_nest_npes(nest_domain, l)
       npes_my_fine = mpp_get_nest_fine_npes(nest_domain,l)
       allocate(my_pelist(npes_my_level))
       allocate(my_pelist_fine(npes_my_fine))
       call mpp_get_nest_pelist(nest_domain, l, my_pelist)
       call mpp_get_nest_fine_pelist(nest_domain, l, my_pelist_fine)

       call mpp_declare_pelist(my_pelist(:))
       write(type2, '(a,I2)')trim(type)//" nest_level = ",l
       if(ANY(my_pelist(:)==mpp_pe())) then
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
                   call fill_coarse_data(x2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                        is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, 0, 0.001, 0.001, 1, 1, &
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
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l, position=NORTH)
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
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0, 1.0E-6, 2.0E-6, 1, 1, &
                     x_cyclic, .false., iend_coarse(1)+1, jend_coarse(1)+1)
             endif
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift, 2.0E-6, 1.0E-6, 1, 1, &
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
          call mpp_get_F2C_index(nest_domain, is_cy, ie_cy, js_cy, je_cy, is_fy, ie_fy, js_fy, je_fy, l, position=NORTH)
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
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0, 1.0E-6, 2.0E-6, 1, -1, &
                     x_cyclic, .false., iend_coarse(1)+1, jend_coarse(1)+1)
             endif
             if( tile == t_coarse(n) .AND. ie_c .GE. is_c .AND. je_c+shift .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift, 2.0E-6, 1.0E-6, -1, 1, &
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
          call mpp_get_F2C_index(nest_domain, is_cx, ie_cx, js_cx, je_cx, is_fx, ie_fx, js_fx, je_fx, l, position=NORTH)
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
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, 0, shift, 1.0E-6, 2.0E-6, 1, -1, &
                     .false., y_cyclic, iend_coarse(1), jend_coarse(1) )
             endif
             if( tile == t_coarse(n) .AND. ie_c+shift .GE. is_c .AND. je_c .GE. js_c ) then
                call fill_coarse_data(y2, rotate_coarse(n), iadd_coarse(n), jadd_coarse(n), &
                     is_c, ie_c, js_c, je_c, nz, isd_coarse, jsd_coarse, nx, ny, shift, 0, 2.0E-6, 1.0E-6, -1, 1, &
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
             write(5000+mpp_pe(),*), "west buffer fine index = ", isw_f, iew_f, jsw_f, jew_f
             write(5000+mpp_pe(),*), "west buffer fine index2 = ", isw_f2, iew_f2, jsw_f2, jew_f2
             write(5000+mpp_pe(),*), "west buffer coarse index = ", isw_c, iew_c, jsw_c, jew_c
             write(5000+mpp_pe(),*), "west buffer coarse index2 = ", isw_c2, iew_c2, jsw_c2, jew_c2
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
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), (/jew_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbuffer2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0, 0.0, 1, 1)
          endif
          call compare_checksums(wbuffer, wbuffer2, trim(type2)//' west buffer coarse to fine scalar')

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), (/jes_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbuffer2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0, 0.0, 1, 1)
          endif
          call compare_checksums(sbuffer, sbuffer2, trim(type2)//' south buffer coarse to fine scalar')

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), (/jee_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(ebuffer2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0, 0.0, 1, 1)
          endif
          call compare_checksums(ebuffer, ebuffer2, trim(type2)//' east buffer coarse to fine scalar')

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), (/jen_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbuffer2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0, 0.0, 1, 1)
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
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH, l, position=CORNER)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH, l, position=CORNER)

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
             jew_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo + shift
             isw_fy2 = isd_fine
             iew_fy2 = isc_fine - 1
             jsw_fy2 = jsd_fine
             jew_fy2 = jed_fine + shift
             isw_cy2 = istart_coarse(my_fine_id)-whalo
             iew_cy2 = istart_coarse(my_fine_id)
             jsw_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jew_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo + shift
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
             jee_cx2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo + shift
             ise_fy2 = iec_fine+1 + shift
             iee_fy2 = ied_fine + shift
             jse_fy2 = jsd_fine
             jee_fy2 = jed_fine + shift
             ise_cy2 = iend_coarse(my_fine_id) + shift
             iee_cy2 = iend_coarse(my_fine_id)+ehalo + shift
             jse_cy2 = jstart_coarse(my_fine_id) + (jsc_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) - shalo
             jee_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo + shift
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_fx2 = isd_fine
             ies_fx2 = ied_fine + shift
             jss_fx2 = jsd_fine
             jes_fx2 = jsc_fine - 1
             iss_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo + shift
             jss_cx2 = jstart_coarse(my_fine_id)-shalo
             jes_cx2 = jstart_coarse(my_fine_id)
             iss_fy2 = isd_fine
             ies_fy2 = ied_fine + shift
             jss_fy2 = jsd_fine
             jes_fy2 = jsc_fine - 1
             iss_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo + shift
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
             ien_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo + shift
             jsn_cx2 = jend_coarse(my_fine_id) + shift
             jen_cx2 = jend_coarse(my_fine_id)+nhalo + shift
             isn_fy2 = isd_fine
             ien_fy2 = ied_fine + shift
             jsn_fy2 = jec_fine+1 + shift
             jen_fy2 = jed_fine + shift
             isn_cy2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ien_cy2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo + shift
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
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), (/jew_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, 1, 1)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), (/jes_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, 1, 1)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), (/jee_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(ebuffery2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, 1, 1)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), (/jen_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 2e3, 1e3, 1, 1)
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
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), (/jew_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbuffer2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0, 0.0, 1, 1)
          endif
          call compare_checksums(wbuffer, wbuffer2, trim(type2)//' west buffer coarse to fine scalar CORNER')

          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), (/jes_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbuffer2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 0.0, 0.0, 1, 1)
          endif
          call compare_checksums(sbuffer, sbuffer2, trim(type2)//' south buffer coarse to fine scalar CORNER')

          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), (/jee_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(ebuffer2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 0.0, 0.0, 1, 1)
          endif
          call compare_checksums(ebuffer, ebuffer2, trim(type2)//' east buffer coarse to fine scalar CORNER')

          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), (/jen_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbuffer2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, shift, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 0.0, 0.0, 1, 1)
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
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, l, position=EAST)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, l, position=EAST)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH, l, position=EAST)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH, l, position=EAST)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, l, position=NORTH)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, l, position=NORTH)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH, l, position=NORTH)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH, l, position=NORTH)

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
             jew_cy2 = jstart_coarse(my_fine_id) + (jec_fine + shift - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo
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
             jee_cy2 = jstart_coarse(my_fine_id) + (jec_fine - jstart_fine(my_fine_id))/y_refine(my_fine_id) + nhalo + shift
          endif
          !--- south
          if( jsc_fine == 1 ) then
             iss_fx2 = isd_fine
             ies_fx2 = ied_fine + shift
             jss_fx2 = jsd_fine
             jes_fx2 = jsc_fine - 1
             iss_cx2 = istart_coarse(my_fine_id) + (isc_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) - whalo
             ies_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo + shift
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
             ien_cx2 = istart_coarse(my_fine_id) + (iec_fine - istart_fine(my_fine_id))/x_refine(my_fine_id) + ehalo + shift
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
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), (/jew_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, 1, 1)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), (/jes_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, 1, 1)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), (/jee_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(ebuffery2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, 1, 1)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), (/jen_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, 1)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 2e3, 1e3, 1, 1)
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
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), (/jew_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, -1, 1)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), (/jes_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, -1, 1)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), (/jee_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(ebuffery2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, -1, 1)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), (/jen_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 2e3, 1e3, -1, 1)
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
          call mpp_get_C2F_index(nest_domain, isw_fx, iew_fx, jsw_fx, jew_fx, isw_cx, iew_cx, jsw_cx, jew_cx, WEST, l, position=NORTH)
          call mpp_get_C2F_index(nest_domain, ise_fx, iee_fx, jse_fx, jee_fx, ise_cx, iee_cx, jse_cx, jee_cx, EAST, l, position=NORTH)
          call mpp_get_C2F_index(nest_domain, iss_fx, ies_fx, jss_fx, jes_fx, iss_cx, ies_cx, jss_cx, jes_cx, SOUTH, l, position=NORTH)
          call mpp_get_C2F_index(nest_domain, isn_fx, ien_fx, jsn_fx, jen_fx, isn_cx, ien_cx, jsn_cx, jen_cx, NORTH, l, position=NORTH)
          call mpp_get_C2F_index(nest_domain, isw_fy, iew_fy, jsw_fy, jew_fy, isw_cy, iew_cy, jsw_cy, jew_cy, WEST, l, position=EAST)
          call mpp_get_C2F_index(nest_domain, ise_fy, iee_fy, jse_fy, jee_fy, ise_cy, iee_cy, jse_cy, jee_cy, EAST, l, position=EAST)
          call mpp_get_C2F_index(nest_domain, iss_fy, ies_fy, jss_fy, jes_fy, iss_cy, ies_cy, jss_cy, jes_cy, SOUTH, l, position=EAST)
          call mpp_get_C2F_index(nest_domain, isn_fy, ien_fy, jsn_fy, jen_fy, isn_cy, ien_cy, jsn_cy, jen_cy, NORTH, l, position=EAST)

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
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isw_c/), (/iew_c/), (/jsw_c/), (/jew_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(wbufferx2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(wbuffery2, isw_c, iew_c, jsw_c, jew_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, -1, 1)
          endif
          if( ies_c .GE. iss_c .AND. jes_c .GE. jss_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/iss_c/), (/ies_c/), (/jss_c/), (/jes_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(sbufferx2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, 0, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(sbuffery2, iss_c, ies_c, jss_c, jes_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, -1, 1)
          endif
          if( iee_c .GE. ise_c .AND. jee_c .GE. jse_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/ise_c/), (/iee_c/), (/jse_c/), (/jee_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(ebufferx2, ise_c, iee_c, jse_c, jee_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(ebuffery2, ise_c+shift, iee_c, jse_c, jee_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse+shift, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, -1, 1)
          endif
          if( ien_c .GE. isn_c .AND. jen_c .GE. jsn_c ) then
             call get_nnest2(domain_coarse, 1, tile_coarse(my_fine_id:my_fine_id), (/isn_c/), (/ien_c/), (/jsn_c/), (/jen_c/), &
                  nnest, t_coarse, iadd_coarse, jadd_coarse, rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse)
             call fill_nest_data(nbufferx2, isn_c, ien_c, jsn_c+shift, jen_c, nnest, t_coarse, 0, shift, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse+shift, je_coarse, 1e3, 2e3, 1, -1)
             call fill_nest_data(nbuffery2, isn_c, ien_c, jsn_c, jen_c, nnest, t_coarse, shift, 0, iadd_coarse, jadd_coarse, &
                  rotate_coarse, is_coarse, ie_coarse, js_coarse, je_coarse, 2e3, 1e3, -1, 1)
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

    call mpp_set_current_pelist()
    deallocate(pelist)

  end subroutine test_update_nest_domain

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

  !############################################################################
  subroutine fill_coarse_data(data, rotate, iadd, jadd, is_c, ie_c, js_c, je_c, nz, isd, jsd, nx, ny, &
                              ishift, jshift, x_add, y_add, sign1, sign2, x_cyclic, y_cyclic, ieg, jeg)
    integer, intent(in)    :: rotate, is_c, ie_c, js_c, je_c, nz, isd, jsd, iadd, jadd, nx, ny, ishift, jshift
    integer, intent(in)    :: sign1, sign2
    real,    intent(inout) :: data(isd:, jsd:, :)
    real,    intent(in)    :: x_add, y_add
    logical, intent(in)    :: x_cyclic, y_cyclic
    integer, intent(in)    :: ieg, jeg
    integer :: i, j, k

    select case (rotate)
    case (ZERO)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = (i+iadd)*1.e+6 + (j+jadd)*1.e+3 + k + x_add
             enddo
          enddo
       enddo
    case (NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign1*((nx-j+1+iadd+jshift)*1.e+6 + (i+jadd)*1.e+3 + k + y_add)
             enddo
          enddo
       enddo
    case (MINUS_NINETY)
       ! convert the index to be consistent with the fine grid.
       do k = 1, nz
          do j = js_c, je_c+jshift
             do i = is_c, ie_c+ishift
                data(i,j,k) = sign2*((j+iadd)*1.e+6 + (ny-i+1+jadd+ishift)*1.e+3 + k + y_add)
             enddo
          enddo
       enddo
    case default
       call mpp_error(FATAL,"fill_coarse_data: rotate_coarse must be ZERO, NINETY, MINUS_NINETY")
    end select

    !---handle cyclic condition
    if(x_cyclic) then
       if(ie_c+ishift+iadd == ieg) then
          i = ie_c+ishift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = i*1.e+6 + (j+jadd)*1.e+3 + k + x_add
             enddo
          enddo
       endif
    endif


    if(y_cyclic) then
       if(je_c+jshift+jadd == jeg) then
          j = je_c+jshift
          do k = 1, nz
             do j = js_c, je_c+jshift
                data(i,j,k) = (i+iadd)*1.e+6 + j*1.e+3 + k + x_add
             enddo
          enddo
       endif
    endif

  end subroutine fill_coarse_data

  subroutine test_get_boundary_ad(type)
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sum
  use mpp_domains_mod, only : CGRID_NE
  use mpp_domains_mod, only : mpp_get_boundary
  use mpp_domains_mod, only : mpp_get_boundary_ad

     character(len=*), intent(in)  :: type

     type(domain2D)       :: domain
     integer              :: ntiles, num_contact, npes_per_tile, ntile_per_pe, layout(2)
     integer              :: n, l, isc, iec, jsc, jec, ism, iem, jsm, jem
     integer, allocatable, dimension(:)       :: tile, ni, nj, pe_start, pe_end
     integer, allocatable, dimension(:,:)     :: layout2D, global_indices

     real*8,  allocatable, dimension(:,:,:) :: x_ad, y_ad, x_fd, y_fd, x_save, y_save
     real*8,  allocatable, dimension(:,:) :: ebufferx2_ad, wbufferx2_ad
     real*8,  allocatable, dimension(:,:) :: sbuffery2_ad, nbuffery2_ad
     real*8 :: ad_sum, fd_sum
     integer :: shift,i,j,k,pe

    !--- check the type
    ntiles = 4
    num_contact = 8

    allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
    allocate(ni(ntiles), nj(ntiles))
    ni(:) = nx; nj(:) = ny
    if( mod(npes, ntiles) == 0 ) then
       npes_per_tile = npes/ntiles
       write(outunit,*)'NOTE from test_uniform_mosaic ==> For Mosaic "', trim(type), &
                       '", each tile will be distributed over ', npes_per_tile, ' processors.'
       ntile_per_pe = 1
       allocate(tile(ntile_per_pe))
       tile = pe/npes_per_tile+1
       call mpp_define_layout( (/1,nx,1,ny/), npes_per_tile, layout )
       do n = 1, ntiles
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do
    else
       call mpp_error(NOTE,'TEST_MPP_DOMAINS: npes should be multiple of ntiles or ' // &
            'ntiles should be multiple of npes. No test is done for '//trim(type) )
       return
    end if

    do n = 1, ntiles
       global_indices(:,n) = (/1,nx,1,ny/)
       layout2D(:,n)         = layout
    end do

    call define_fourtile_mosaic(type, domain, (/nx,nx,nx,nx/), (/ny,ny,ny,ny/), global_indices, &
                                layout2D, pe_start, pe_end, .true. )

    call mpp_get_compute_domain( domain, isc, iec, jsc, jec )
    call mpp_get_memory_domain( domain, ism, iem, jsm, jem )

    deallocate(layout2D, global_indices, pe_start, pe_end )
    deallocate(ni, nj)

    shift = 1
    allocate( x_ad  (ism:iem+shift,jsm:jem  ,nz) )
    allocate( x_fd  (ism:iem+shift,jsm:jem  ,nz) )
    allocate( x_save(ism:iem+shift,jsm:jem  ,nz) )
    allocate( y_ad  (ism:iem  ,jsm:jem+shift,nz) )
    allocate( y_fd  (ism:iem  ,jsm:jem+shift,nz) )
    allocate( y_save(ism:iem  ,jsm:jem+shift,nz) )
    allocate(ebufferx2_ad(jec-jsc+1, nz), wbufferx2_ad(jec-jsc+1, nz))
    allocate(sbuffery2_ad(iec-isc+1, nz), nbuffery2_ad(iec-isc+1, nz))

    pe = mpp_pe()

    x_fd=0; y_fd=0
    do k = 1,nz
      do j = jsc,jec
        do i = isc,iec
            x_fd(i,j,k)= i*j
            y_fd(i,j,k)= i*j
        end do
      end do
    end do

    x_save=x_fd
    y_save=y_fd

    ebufferx2_ad = 0
    wbufferx2_ad = 0
    sbuffery2_ad = 0
    nbuffery2_ad = 0

    call mpp_get_boundary(x_fd, y_fd, domain, ebufferx=ebufferx2_ad(:,:), wbufferx=wbufferx2_ad(:,:), &
                             sbuffery=sbuffery2_ad(:,:), nbuffery=nbuffery2_ad(:,:), gridtype=CGRID_NE,  &
                             complete = .true.  )
    fd_sum = 0.
    do k = 1,nz
      do j = jsc,jec
        do i = isc,iec
           fd_sum = fd_sum + x_fd(i,j,k)*x_fd(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsc,jec
        do i = isc,iec
           fd_sum = fd_sum + y_fd(i,j,k)*y_fd(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
        do i = 1,jec-jsc+1
           fd_sum = fd_sum + ebufferx2_ad(i,k)*ebufferx2_ad(i,k)
        end do
    end do
    do k = 1,nz
        do i = 1,jec-jsc+1
           fd_sum = fd_sum + wbufferx2_ad(i,k)*wbufferx2_ad(i,k)
        end do
    end do
    do k = 1,nz
        do i = 1,iec-isc+1
           fd_sum = fd_sum + sbuffery2_ad(i,k)*sbuffery2_ad(i,k)
        end do
    end do
    do k = 1,nz
        do i = 1,iec-isc+1
           fd_sum = fd_sum + nbuffery2_ad(i,k)*nbuffery2_ad(i,k)
        end do
    end do
    call mpp_sum( fd_sum )

    x_ad = x_fd
    y_ad = y_fd

    call mpp_get_boundary_ad(x_ad, y_ad, domain, ebufferx=ebufferx2_ad(:,:), wbufferx=wbufferx2_ad(:,:), &
                             sbuffery=sbuffery2_ad(:,:), nbuffery=nbuffery2_ad(:,:), gridtype=CGRID_NE,  &
                             complete = .true.  )

    ad_sum = 0.
    do k = 1,nz
      do j = jsc,jec
        do i = isc,iec
           ad_sum = ad_sum + x_ad(i,j,k)*x_save(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsc,jec
        do i = isc,iec
           ad_sum = ad_sum + y_ad(i,j,k)*y_save(i,j,k)
        end do
      end do
    end do
    call mpp_sum( ad_sum )

    if( pe.EQ.mpp_root_pe() ) then
       if (abs(ad_sum-fd_sum)/fd_sum.lt.1E-7) then
           print*, "Passed Adjoint Dot Test: mpp_get_boundary_ad"
       endif
    endif

    deallocate (x_ad, y_ad, x_fd, y_fd, x_save, y_save)
    deallocate (ebufferx2_ad, wbufferx2_ad)
    deallocate (sbuffery2_ad, nbuffery2_ad)

  end subroutine test_get_boundary_ad

  subroutine test_halo_update_ad( type )
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sum
  use mpp_domains_mod, only : CGRID_NE
  use mpp_domains_mod, only : mpp_update_domains, mpp_update_domains_ad

    character(len=*), intent(in) :: type
    type(domain2D) :: domain

    integer              :: shift, i, j, k
    logical              :: is_symmetry
    integer              :: is, ie, js, je, isd, ied, jsd, jed, pe

    real*8,  allocatable, dimension(:,:,:) :: x_ad, y_ad, x_fd, y_fd, x_save, y_save
    real*8 :: ad_sum, fd_sum

    if(index(type, 'symmetry') == 0) then
       is_symmetry = .false.
    else
       is_symmetry = .true.
    end if
    select case(type)
    case( 'Simple', 'Simple symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                 shalo=shalo, nhalo=nhalo, name=type, symmetry = is_symmetry )
    case( 'Cyclic', 'Cyclic symmetry' )
        call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
        call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo,        &
             shalo=shalo, nhalo=nhalo, xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN, &
             name=type, symmetry = is_symmetry )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type )
    end select

!set up x array
    call mpp_get_compute_domain( domain, is,  ie,  js,  je  )
    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    shift=1
!---test 3d single fields----------------------------------------------------------
    allocate( x_fd(isd:ied,jsd:jed,nz) )
    allocate( x_ad(isd:ied,jsd:jed,nz) )
    allocate( x_save(isd:ied,jsd:jed,nz) )
    x_fd = 0.; x_ad = 0.; x_save = 0.

    do k = 1,nz
       do j = js,je
          do i = is,ie
             x_fd(i,j,k) = i*j
          end do
       end do
    end do
    x_save = x_fd

!full update
    call mpp_update_domains( x_fd, domain )

    fd_sum = 0.
    do k = 1,nz
       do j = jsd,jed
          do i = isd,ied
             fd_sum = fd_sum + x_fd(i,j,k)*x_fd(i,j,k)
          end do
       end do
    end do
    call mpp_sum( fd_sum )

    x_ad = x_fd
    call mpp_update_domains_ad( x_ad, domain )

    ad_sum = 0.
    do k = 1,nz
       do j = jsd,jed
          do i = isd,ied
             ad_sum = ad_sum + x_ad(i,j,k)*x_save(i,j,k)
          end do
       end do
    end do
    call mpp_sum( ad_sum )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
       if (abs(ad_sum-fd_sum)/fd_sum.lt.1E-7) then
           print*, "Passed Adjoint Dot Test: mpp_update_domains_ad(single 3D field)"
       endif
    endif

    deallocate (x_ad, x_fd, x_save)

!---test 3d vector fields----------------------------------------------------------
    allocate( x_ad  (isd:ied+shift,jsd:jed  ,nz) )
    allocate( x_fd  (isd:ied+shift,jsd:jed  ,nz) )
    allocate( x_save(isd:ied+shift,jsd:jed  ,nz) )
    allocate( y_ad  (isd:ied  ,jsd:jed+shift,nz) )
    allocate( y_fd  (isd:ied  ,jsd:jed+shift,nz) )
    allocate( y_save(isd:ied  ,jsd:jed+shift,nz) )

    x_fd=0; y_fd=0
    do k = 1,nz
      do j = js,je
        do i = is,ie
           x_fd(i,j,k)=i*j
           y_fd(i,j,k)=i*j
        end do
      end do
    end do

    call mpp_update_domains( x_fd, y_fd, domain, gridtype=CGRID_NE)
    x_save=x_fd
    y_save=y_fd

    fd_sum = 0.
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied+shift
           fd_sum = fd_sum + x_fd(i,j,k)*x_fd(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsd,jed+shift
        do i = isd,ied
           fd_sum = fd_sum + y_fd(i,j,k)*y_fd(i,j,k)
        end do
      end do
    end do
    call mpp_sum( fd_sum )

    x_ad = x_fd
    y_ad = y_fd
    call mpp_update_domains_ad( x_ad, y_ad, domain, gridtype=CGRID_NE)

    ad_sum = 0.
    do k = 1,nz
      do j = jsd,jed
        do i = isd,ied+shift
           ad_sum = ad_sum + x_ad(i,j,k)*x_save(i,j,k)
        end do
      end do
    end do
    do k = 1,nz
      do j = jsd,jed+shift
        do i = isd,ied
           ad_sum = ad_sum + y_ad(i,j,k)*y_save(i,j,k)
        end do
      end do
    end do
    call mpp_sum( ad_sum )

    if( pe.EQ.mpp_root_pe() ) then
       if (abs(ad_sum-fd_sum)/fd_sum.lt.1E-7) then
           print*, "Passed Adjoint Dot Test: mpp_update_domains_ad(vector 3D fields)"
       endif
    endif
    deallocate (x_ad, y_ad, x_fd, y_fd, x_save, y_save)

  end subroutine test_halo_update_ad

  subroutine test_global_reduce_ad (type)
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe, mpp_sum
  use mpp_domains_mod, only : mpp_global_sum_tl, mpp_global_sum_ad
    character(len=*), intent(in) :: type
    real    :: gsum_tl, gsum_ad
    real*8  :: gsum_tl_save, gsum_ad_save
    real    :: gsum_tl_bit, gsum_ad_bit
    real*8  :: gsum_tl_save_bit, gsum_ad_save_bit
    integer :: i,j,k, ishift, jshift, position
    integer :: isd, ied, jsd, jed

    type(domain2D) :: domain
    real, allocatable, dimension(:,:,:) :: x, x_ad, x_ad_bit

    !--- set up domain
    call mpp_define_layout( (/1,nx,1,ny/), npes, layout )
    select case(type)
    case( 'Simple' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type )
    case( 'Simple symmetry center', 'Simple symmetry corner', 'Simple symmetry east', 'Simple symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, &
                                    shalo=shalo, nhalo=nhalo, name=type, symmetry = .true. )
    case( 'Cyclic symmetry center', 'Cyclic symmetry corner', 'Cyclic symmetry east', 'Cyclic symmetry north' )
           call mpp_define_domains( (/1,nx,1,ny/), layout, domain, whalo=whalo, ehalo=ehalo, shalo=shalo, nhalo=nhalo, &
                                    name=type, symmetry = .true., xflags=CYCLIC_GLOBAL_DOMAIN, yflags=CYCLIC_GLOBAL_DOMAIN )
    case default
        call mpp_error( FATAL, 'TEST_MPP_DOMAINS: no such test: '//type//' in test_global_field' )
    end select

    call mpp_get_data_domain   ( domain, isd, ied, jsd, jed )

    !--- determine if an extra point is needed
    ishift = 0; jshift = 0; position = CENTER
    select case(type)
    case ('Simple symmetry corner', 'Cyclic symmetry corner')
       ishift = 1; jshift = 1; position = CORNER
    case ('Simple symmetry east', 'Cyclic symmetry east' )
       ishift = 1; jshift = 0; position = EAST
    case ('Simple symmetry north', 'Cyclic symmetry north')
       ishift = 0; jshift = 1; position = NORTH
    end select

    ied = ied+ishift; jed = jed+jshift

    allocate( x(isd:ied,jsd:jed,nz), x_ad(isd:ied,jsd:jed,nz), x_ad_bit(isd:ied,jsd:jed,nz) )

    x=0.
    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           x(i,j,k) = i+j+k
         enddo
       enddo
    enddo

    gsum_tl      = mpp_global_sum( domain, x, position = position  )
    gsum_tl_bit  = mpp_global_sum( domain, x, flags=BITWISE_EXACT_SUM  )
    gsum_tl_save = gsum_tl*gsum_tl
    gsum_tl_save_bit = gsum_tl_bit*gsum_tl_bit

    gsum_ad      = gsum_tl
    gsum_ad_bit  = gsum_tl_bit

    x_ad     = 0.
    x_ad_bit = 0.
    call mpp_global_sum_ad( domain, x_ad, gsum_ad, position = position )
    call mpp_global_sum_ad( domain, x_ad_bit, gsum_ad_bit, flags = BITWISE_EXACT_SUM )

    gsum_ad_save     = 0.
    gsum_ad_save_bit = 0.

    do k = 1,nz
       do j = jsd, jed
         do i = isd, ied
           gsum_ad_save     = gsum_ad_save + x_ad(i,j,k)*x(i,j,k)
           gsum_ad_save_bit = gsum_ad_save_bit + x_ad_bit(i,j,k)*x(i,j,k)
         enddo
       enddo
    enddo

    call mpp_sum( gsum_ad_save )
    call mpp_sum( gsum_ad_save_bit )

    pe = mpp_pe()
    if( pe.EQ.mpp_root_pe() ) then
       if (abs(gsum_ad_save-gsum_tl_save)/gsum_tl_save.lt.1E-7) then
           print*, "Passed Adjoint Dot Test: mpp_global_sum_ad"
       endif
       if (abs(gsum_ad_save_bit-gsum_tl_save_bit)/gsum_tl_save_bit.lt.1E-7) then
           print*, "Passed Adjoint Dot Test: mpp_global_sum_ad, flags=BITWISE_EXACT_SUM"
       endif
    endif

    deallocate(x, x_ad, x_ad_bit)

  end subroutine test_global_reduce_ad

end program test_mpp_domains
