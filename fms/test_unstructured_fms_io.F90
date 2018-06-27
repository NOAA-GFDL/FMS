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
#ifdef TEST_UNSTRUCTURED_FMS_IO
program test_unstructured_fms_io
#include <fms_platform.h>
    use,intrinsic :: iso_fortran_env, only: output_unit
    use mpp_parameter_mod,            only: FATAL, &
                                            NOTE, &
                                            MPP_DEBUG, &
                                            MPP_CLOCK_SYNC
    use mpp_mod,                      only: mpp_init, &
                                            mpp_pe, &
                                            mpp_npes, &
                                            mpp_root_pe, &
                                            mpp_error, &
                                            mpp_set_stack_size, &
                                            mpp_exit, &
                                            mpp_clock_begin, &
                                            mpp_clock_end, &
                                            mpp_clock_id
#ifdef INTERNAL_FILE_NML
    use mpp_mod,                      only: input_nml_file
#endif
    use mpp_domains_mod,              only: mpp_domains_init, &
                                            mpp_domains_set_stack_size, &
                                            mpp_domains_exit, &
                                            domain2D, &
                                            domainUG
    use mpp_io_mod,                   only: mpp_io_init, &
                                            mpp_io_exit
    use fms_io_mod,                   only: fms_io_init, &
                                            fms_io_exit
    implicit none

#ifdef use_netCDF
#include <netcdf.inc>
#endif

   !Local variables
    integer(INT_KIND)              :: nx = 8                               !<Total number of grid points in the x-dimension (longitude?)
    integer(INT_KIND)              :: ny = 8                               !<Total number of grid points in the y-dimension (latitude?)
    integer(INT_KIND)              :: nz = 2                               !<Total number of grid points in the z-dimension (height)
!   integer(INT_KIND)              :: nx = 128                             !<Total number of grid points in the x-dimension (longitude?)
!   integer(INT_KIND)              :: ny = 128                             !<Total number of grid points in the y-dimension (latitude?)
!   integer(INT_KIND)              :: nz = 40                              !<Total number of grid points in the z-dimension (height)
    integer(INT_KIND)              :: nt = 2                               !<Total number of time grid points.
    integer(INT_KIND)              :: halo = 2                             !<Number of grid points in the halo???
    integer(INT_KIND)              :: ntiles_x = 1                         !<Number of tiles in the x-direction (A 2D grid of tiles is used in this test.)
    integer(INT_KIND)              :: ntiles_y = 2                         !<Number of tiles in the y-direction (A 2D grid of tiles is used in this test.)
    integer(INT_KIND)              :: total_num_tiles                      !<The total number of tiles for the run (= ntiles_x*ntiles_y)
    integer(INT_KIND),dimension(2) :: layout(2) = (/1,1/)                  !<Rank layout (x,y of a 2D grid) for each tile.  Total number of ranks = total_num_tiles*layout(1)*layout(2)
    integer(INT_KIND),dimension(2) :: io_layout(2) = (/1,1/)               !<Layout (x,y of a 2D grid) for ranks that will perform I/O for each tile.  These dimensions must divide equally
                                                                           !!into the dimensions of the layout array.
    integer(INT_KIND)              :: stackmax = 1500000                   !<Default size to which the mpp stack will be set.
    integer(INT_KIND)              :: stackmaxd = 500000                   !<Default size to which the mpp_domains stack will be set.
    logical(INT_KIND)              :: debug = .false.                      !<Flag to print debugging information.
    character(len=64)              :: test_file = "test_unstructured_grid" !<Base filename for the unit tests.
    character(len=64)              :: iospec = '-F cachea'                 !<Something cray related ???
    integer(INT_KIND)              :: pack_size = 1                        !<(Number of bits in real(DOUBLE_KIND))/(Number of bits in real)
    integer(INT_KIND)              :: npes                                 !<Total number of ranks in the current pelist.
    integer(INT_KIND)              :: io_status                            !<Namelist read error code.
    real(DOUBLE_KIND)              :: doubledata = 0.0                     !<Used to determine pack_size.  This must be kind=DOUBLE_KIND.
    real                           :: realdata = 0.0                       !<Used to determine pack_size.  Do not specify a kind parameter.
    integer(INT_KIND)              :: funit = 7                            !<File unit.
    logical(INT_KIND)              :: fopened                              !<Flag telling if a file is already open.
    type(domain2D)                 :: structured_domain                    !<A structured 2D domain.
    type(domainUG)                 :: unstructured_domain                  !<An unstructured mpp domain.
    integer(INT_KIND)              :: test_num                             !<Which test to perform.
    integer(INT_KIND),parameter    :: test_1_id = 1                        !<Test 1 id.
!   integer(INT_KIND)              :: id_single_tile_mult_file             !<Mpp timer id.
!   integer(INT_KIND)              :: id_mult_tile                         !<Mpp timer id.
!   integer(INT_KIND)              :: id_single_tile_with_group            !<Mpp timer id.
!   integer(INT_KIND)              :: id_mult_tile_with_group              !<Mpp timer id.

    namelist /test_unstructured_io_nml/ nx, &
                                        ny, &
                                        nz, &
                                        nt, &
                                        halo, &
                                        stackmax, &
                                        stackmaxd, &
                                        debug, &
                                        test_file, &
                                        iospec, &
                                        ntiles_x, &
                                        ntiles_y, &
                                        layout, &
                                        io_layout

   !Initialize mpp, get the rank of the current process and the number
   !of ranks in the current pelist.
    call mpp_init() 
    npes = mpp_npes()

   !Get the values from the namelist file.
#ifdef INTERNAL_FILE_NML
    read (input_nml_file,test_unstructured_io_nml,iostat=io_status)
#else
    do
        inquire(unit=funit,opened=fopened)
        if (.not. fopened) then
            exit
        endif
        funit = funit + 1
        if (funit .eq. 100) then
            call mpp_error(FATAL, &
                           "test_unstructured_io: Unable to locate unit" &
                           //" number for the input.nml file.")
        endif
    enddo
    open(unit=funit,file='input.nml',iostat=io_status)
    read(funit,test_unstructured_io_nml,iostat=io_status)
    close(funit)
#endif

   !Check the namelist read error code.
    if (io_status > 0) then
        call mpp_error(FATAL, &
                       "test_unstructured_io: Error reading input.nml")
    endif

   !Initialize the mpp_domains module
    if (debug) then
        call mpp_domains_init(MPP_DEBUG)
    else
        call mpp_domains_init()
    endif

   !Initialize the mpp_io module.
    if (debug) then
        call mpp_io_init(MPP_DEBUG)
    else
        call mpp_io_init()
    endif

   !Initialize the fms_io module.
    call fms_io_init()

   !Set the mpp and mpp_domains stack sizes.
    call mpp_set_stack_size(stackmax)
    call mpp_domains_set_stack_size(stackmaxd)

   !Write out test configuration parameters.
    if (mpp_pe() .eq. mpp_root_pe()) then
        write(output_unit,*)
        write(output_unit,*) "Performing unstructured_io unit test with:"
        write(output_unit,*) "Total number of ranks:                          ", &
                             npes
        write(output_unit,*) "Total number of grid points in the x-dimension: ", &
                             nx
        write(output_unit,*) "Total number of grid points in the y-dimension: ", &
                             ny
        write(output_unit,*) "Total number of grid points in the z-dimension: ", &
                             nz
        write(output_unit,*) "Total number of grid points in the t-dimension: ", &
                             nt
        write(output_unit,*) "Halo width (# of grid points):                  ", &
                             halo
        write(output_unit,*) "Using NEW domaintypes and calls..."
    endif

   !Add a suffix to the test file.
    write(test_file,'(a,i3.3)') trim(test_file),npes

   !Perform the chosen test.
    test_num = test_1_id
    select case (test_num)
        case (test_1_id)
            !Test 1: This test simulates the model running, saving the restart
            !        files, and then restarting a user specified number of
            !        times.  Each time the model is "restarted", the data is
            !        checked to make sure that it matches the data that saved
            !        in the restart.
             if (mpp_pe() .eq. mpp_root_pe()) then
                 write(output_unit,*)
                 write(output_unit,*) "///////////////////////////////////////"
                 write(output_unit,*) "Performing test 1: ..."
             endif

            !Create a structured 2D mpp domain and an unstructured mpp domain
            !for the test.
             call create_mpp_domains(nx, &
                                     ny, &
                                     nz, &
                                     npes, &
                                     ntiles_x, &
                                     ntiles_y, &
                                     structured_domain, &
                                     unstructured_domain)
            !Perform test 1.
             call test_1(unstructured_domain, &
                         10, &
                         nx, &
                         ny, &
                         nz, &
                         npes)

            !Free memory allocated to the domains for this test.
             call destroy_mpp_domains(structured_domain, &
                                      unstructured_domain)

            !Test 1 complete.
             if (mpp_pe() .eq. mpp_root_pe()) then
                 write(output_unit,*) "Test 1 complete."
                 write(output_unit,*) "///////////////////////////////////////"
                 write(output_unit,*)
             endif

        case default
            !No test was selected, so throw an error
             call mpp_error(FATAL, &
                            "test_unstructured_io: invalid test specified.")
    end select

   !Finalize the fms_io, mpp_io, mpp_domains, and mpp modules.
    call fms_io_exit()
    call mpp_io_exit()
    call mpp_domains_exit()
    call mpp_exit()

contains

   !---------------------------------------------------------------------------
   !Create a structured 2D mpp domain and an unstructured mpp domain.
    subroutine create_mpp_domains(nx, &
                                  ny, &
                                  nz, &
                                  npes, &
                                  num_domain_tiles_x, &
                                  num_domain_tiles_y, &
                                  structured_domain, &
                                  unstructured_domain)
        use, intrinsic :: iso_fortran_env, only: output_unit
        use mpp_parameter_mod,             only: FATAL
        use mpp_mod,                       only: mpp_error, &
                                                 mpp_pe, &
                                                 mpp_root_pe, &
                                                 mpp_sync
        use mpp_domains_mod,               only: domain2D, &
                                                 mpp_define_mosaic, &
                                                 domainUG, &
                                                 mpp_define_unstruct_domain

       !Inputs/Ouputs
        integer(INT_KIND),intent(in) :: nx                  !<The number of grid points in the x-direction.
        integer(INT_KIND),intent(in) :: ny                  !<The number of grid points in the y-direction.
        integer(INT_KIND),intent(in) :: nz                  !<The number of grid points in the z-direction.
        integer(INT_KIND),intent(in) :: npes                !<The total number of ranks used in this test.
        integer(INT_KIND),intent(in) :: num_domain_tiles_x  !<The total number of domain tiles in the x-dimension for the 2D structured domain in this test.
        integer(INT_KIND),intent(in) :: num_domain_tiles_y  !<The total number of domain tiles in the y-dimension for the 2D structured domain in this test.
        type(domain2D),intent(inout) :: structured_domain   !<A structured 2D domain.
        type(domainUG),intent(inout) :: unstructured_domain !<An unstructured mpp domain.

       !Local variables
        integer(INT_KIND)                              :: num_domain_tiles                           !<The total number of domain tiles for the 2D structured domain in this test.
        integer(INT_KIND)                              :: npes_per_domain_tile                       !<The number of ranks per domain tile for the 2D structured domain.
        integer(INT_KIND)                              :: my_domain_tile_id                          !<The 2D structured domain tile id for the current rank.
        logical(INT_KIND)                              :: is_domain_tile_root                        !<Flag telling if the current rank is the root rank of its associated
                                                                                                     !!2D structured domain tile.
        integer(INT_KIND),dimension(2)                 :: layout_for_full_domain                     !<Rank layout (2D grid) for the full 2D structured domain.
                                                                                                     !!Example: 16 ranks -> (16,1) or (8,2) or (4,4) or (2,8) or (1,16)
        integer(INT_KIND),dimension(:),allocatable     :: pe_start                                   !<Array holding the smallest rank id assigned to each 2D structured domain tile.
        integer(INT_KIND),dimension(:),allocatable     :: pe_end                                     !<Array holding the largest rank id assigned to each 2D structured domain tile.
        integer(INT_KIND)                              :: x_grid_points_per_domain_tile              !<The number of grid points in the x-dimension on each 2D structured domain tile.
        integer(INT_KIND)                              :: y_grid_points_per_domain_tile              !<The number of grid points in the y-dimension on each 2D structured domain tile.
        integer(INT_KIND),dimension(:,:),allocatable   :: global_indices                             !<Required to define the 2D structured domain.
        integer(INT_KIND),dimension(:,:),allocatable   :: layout2D                                   !<Required to define the 2D structured domain.
        logical(INT_KIND),dimension(:,:,:),allocatable :: land_mask                                  !<A toy mask.
        integer(INT_KIND),dimension(:),allocatable     :: num_non_masked_grid_points_per_domain_tile !<Total number of non-masked grid points on each 2D structured domain tile.
        integer(INT_KIND)                              :: mask_counter                               !<Counting variable.
        integer(INT_KIND)                              :: num_non_masked_grid_points                 !<Total number of non-masked grid points for the 2D structured domain.
        integer(INT_KIND),dimension(:),allocatable     :: num_land_tiles_per_non_masked_grid_point   !<Number of land tiles per non-masked grid point for the 2D structured domain.
        integer(INT_KIND)                              :: num_ranks_using_unstructured_grid          !<Number of ranks using the unstructured domain.
        integer(INT_KIND)                              :: io_tile_factor                             !<I/O tile factor.  See below.
        integer(INT_KIND),dimension(:),allocatable     :: unstructured_grid_point_index_map          !<Array that maps indices between the 2D structured and unstructured domains.
        integer(INT_KIND)                              :: i                                          !<Loop variable.
        integer(INT_KIND)                              :: j                                          !<Loop variable.
        integer(INT_KIND)                              :: k                                          !<Loop variable.
        integer(INT_KIND)                              :: p                                          !<Counting variable.

       !Needed to define the 2D structured domain but not otherwised used.
        integer(INT_KIND)              :: ncontacts
        integer(INT_KIND),dimension(1) :: tile1
        integer(INT_KIND),dimension(1) :: tile2
        integer(INT_KIND),dimension(1) :: istart1
        integer(INT_KIND),dimension(1) :: iend1
        integer(INT_KIND),dimension(1) :: jstart1
        integer(INT_KIND),dimension(1) :: jend1
        integer(INT_KIND),dimension(1) :: istart2
        integer(INT_KIND),dimension(1) :: iend2
        integer(INT_KIND),dimension(1) :: jstart2
        integer(INT_KIND),dimension(1) :: jend2

       !Print out a message that the routine is starting.
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*)
            write(output_unit,*) "Creating a structured and unstructured" &
                                 //" domain ..."
        endif

       !Synchronize all ranks.
        call mpp_sync()

       !Make sure that valid inputs were passed in.
        if (nx .lt. 1 .or. ny .lt. 1) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" there must be at least on grid point in the" &
                           //" x- and y- dimensions.")
        endif
        if (npes .gt. nx*ny) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" the total number of ranks cannot be greater" &
                           //" than the total number of grid points in the" &
                           //" x-y plane.")
        endif
        if (num_domain_tiles_x .lt. 1 .or. num_domain_tiles_y .lt. 1) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" there must be at least on domain tile in the" &
                           //" x- and y- dimensions.")
        endif
        if (mod(nx,num_domain_tiles_x) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" the total number of grid points in the" &
                           //" x-dimension must be evenly divisible by the" &
                           //" total number of domain tiles in the" &
                           //" x-dimension.")
        endif
        if (mod(ny,num_domain_tiles_y) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" the total number of grid points in the" &
                           //" y-dimension must be evenly divisible by the" &
                           //" total number of domain tiles in the" &
                           //" y-dimension.")
        endif
        if (num_domain_tiles_x*num_domain_tiles_y .gt. npes) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" the total number of domain tiles cannot be" &
                           //" greater than the total number of ranks.")
        endif
        if (mod(npes,num_domain_tiles_x) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" the total number of ranks must be evenly" &
                           //" divisible by the total number of domain" &
                           //" tiles in the x-dimension.")
        endif
        if (mod(npes,num_domain_tiles_y) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" the total number of ranks must be evenly" &
                           //" divisible by the total number of domain" &
                           //" tiles in the y-dimension.")
        endif

       !Set domain tile values for the 2D structured domain.
        num_domain_tiles = num_domain_tiles_x*num_domain_tiles_y
        npes_per_domain_tile = npes/num_domain_tiles
        my_domain_tile_id = (mpp_pe())/npes_per_domain_tile + 1
        if (mpp_pe() .eq. (my_domain_tile_id-1)*npes_per_domain_tile) then
            is_domain_tile_root = .true.
        else
            is_domain_tile_root = .false.
        endif
        layout_for_full_domain(1) = num_domain_tiles_x
        layout_for_full_domain(2) = npes/layout_for_full_domain(1)

       !For each 2D structured domain tile, store the beginning and ending
       !rank ids assigned to it.  For example, if there are 8 ranks and 2
       !domain tiles, then tile 1 will be assigned ranks 0 - 3 and tile 2
       !will be assigned ranks 4 - 7.
        allocate(pe_start(num_domain_tiles))
        allocate(pe_end(num_domain_tiles))
        do i = 1,num_domain_tiles
            pe_start(i) = (i-1)*npes_per_domain_tile
            pe_end(i) = i*npes_per_domain_tile - 1
        enddo

       !Calculate parameters needed to construct the 2D structured domain.
       !All domain tiles are assumed to be the same size.
        x_grid_points_per_domain_tile = nx/num_domain_tiles_x
        y_grid_points_per_domain_tile = ny/num_domain_tiles_y
        allocate(global_indices(4,num_domain_tiles))
        do i = 1,num_domain_tiles
            global_indices(:,i) = (/1,x_grid_points_per_domain_tile, &
                                    1,y_grid_points_per_domain_tile/)
        enddo
        allocate(layout2D(2,num_domain_tiles))
        do i = 1,num_domain_tiles
            layout2D(1,i) = layout_for_full_domain(1)/num_domain_tiles_x
            layout2D(2,i) = layout_for_full_domain(2)/num_domain_tiles_y
        enddo

       !This test does not use the "contact" region between tiles, but
       !the 2D structured domain requires these inputs, so just set them
       !all equal to 1.
        ncontacts = 0
        tile1 = 1
        tile2 = 1
        istart1 = 1
        iend1 = 1
        jstart1 = 1
        jend1 = 1
        istart2 = 1
        iend2 = 1
        jstart2 = 1
        jend2 = 1

       !Define the 2D structured domain.
        call mpp_define_mosaic(global_indices, &
                               layout2D, &
                               structured_domain, &
                               num_domain_tiles, &
                               ncontacts, &
                               tile1, &
                               tile2, &
                               istart1, &
                               iend1, &
                               jstart1, &
                               jend1, &
                               istart2, &
                               iend2, &
                               jstart2, &
                               jend2, &
                               pe_start, &
                               pe_end)

       !Define a toy mask to mimic what happens in the land model.
        allocate(land_mask(x_grid_points_per_domain_tile, &
                           y_grid_points_per_domain_tile, &
                           num_domain_tiles))
        allocate(num_non_masked_grid_points_per_domain_tile(num_domain_tiles))
        land_mask = .false.
        do k = 1,num_domain_tiles
            mask_counter = 0
            do j = 1,y_grid_points_per_domain_tile
                do i = 1,x_grid_points_per_domain_tile
                    if (mod((k-1)*y_grid_points_per_domain_tile*x_grid_points_per_domain_tile + &
                            (j-1)*x_grid_points_per_domain_tile + &
                            (i-1),2) .eq. 0) then
                        land_mask(i,j,k) = .true.
                        mask_counter = mask_counter + 1
                    endif
                enddo
            enddo
            num_non_masked_grid_points_per_domain_tile(k) = mask_counter
        enddo

       !Set the number of land tiles allowed per non-masked grid point.
        num_non_masked_grid_points = sum(num_non_masked_grid_points_per_domain_tile)
        allocate(num_land_tiles_per_non_masked_grid_point(num_non_masked_grid_points))
        num_land_tiles_per_non_masked_grid_point = 1

       !Set the number of ranks to use with the unstructured domain.  There
       !must be at least one grid point per rank.
        num_ranks_using_unstructured_grid = npes
        if (num_ranks_using_unstructured_grid .gt. num_non_masked_grid_points) then
            call mpp_error(FATAL, &
                           "create_mpp_domains:" &
                           //" the number of ranks exceeds the number of" &
                           //" non-masked grid points for the unstructured" &
                           //" domain.")
        endif

       !Define the number of "I/O tile factor".  The number of ranks that
       !participate in I/O for a tile is equal to: 
       !
       ! num_io_ranks_on_a_tile = num_ranks_on_the_tile / "I/O tile factor".
       !
       !so for:
       !
       ! io_tile_factor = 1, all of the ranks on a tile participate in the I/O
       ! io_tile_factor = 2, 1/2 of the ranks on a tile participate in the I/O
       ! io_tile_factor = 3, 1/3 of the ranks on a tile participate in the I/O
       ! ...
       ! io_tile_factor = 0 is a special case where only one rank participates
       !                  in the I/O for a tile.
        io_tile_factor = 0

       !Define an array used to map grid points from the "structured" 2D grid
       !to the "unstructured" 1D grid.  The mapping goes as follows (fortran
       !ording so first index is fastest):
       !
       ! 2D "structured" grid (lon,lat,tile) => 1D "unstructured" grid (p)
       !
       !where masked points are skipped.
        allocate(unstructured_grid_point_index_map(num_non_masked_grid_points))
        p = 0
        do k = 1,num_domain_tiles
            do j = 1,y_grid_points_per_domain_tile
                do i = 1,x_grid_points_per_domain_tile
                    if (land_mask(i,j,k)) then
                        p = p + 1
                        unstructured_grid_point_index_map(p) = (j-1)*x_grid_points_per_domain_tile + i
                    endif
                enddo
            enddo
        enddo

       !Define the "unstructured" domain decomposition.
        call mpp_define_unstruct_domain(unstructured_domain, &
                                        structured_domain, &
                                        num_non_masked_grid_points_per_domain_tile, &
                                        num_land_tiles_per_non_masked_grid_point, &
                                        num_ranks_using_unstructured_grid, &
                                        io_tile_factor, &
                                        unstructured_grid_point_index_map)

       !Print out information about the unstructured domain.
!       if (mpp_pe() .eq. mpp_root_pe()) then
!           write(*,*) 'num tiles=',mpp_get_UG_domain_ntiles(domain_ug)
!       endif

       !Deallocate local allocatables.
        deallocate(pe_start)
        deallocate(pe_end)
        deallocate(global_indices)
        deallocate(layout2D)
        deallocate(land_mask)
        deallocate(num_non_masked_grid_points_per_domain_tile)
        deallocate(num_land_tiles_per_non_masked_grid_point)
        deallocate(unstructured_grid_point_index_map)

       !Print out a message that the routine is done.
        call mpp_sync()
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*) "Domains created."
            write(output_unit,*)
        endif

        return
    end subroutine create_mpp_domains

   !---------------------------------------------------------------------------
   !Destroy a 2D mpp domain and an unstructured mpp domain.
    subroutine destroy_mpp_domains(structured_domain, &
                                   unstructured_domain)
        use, intrinsic :: iso_fortran_env, only: output_unit
        use mpp_mod,                       only: mpp_pe, &
                                                 mpp_root_pe, &
                                                 mpp_sync
        use mpp_domains_mod,               only: domain2D, &
                                                 mpp_deallocate_domain, &
                                                 domainUG, &
                                                 mpp_deallocate_domainUG

       !Inputs/Ouputs
        type(domain2D),intent(inout) :: structured_domain   !<A structured 2D domain.
        type(domainUG),intent(inout) :: unstructured_domain !<An unstructured mpp domain.

       !Print out a message that the routine is starting.
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*)
            write(output_unit,*) "Creating a structured and unstructured" &
                                 //" domain ..."
        endif

       !Synchronize all ranks.
        call mpp_sync()

       !Deallocate the unstructured domain.
        call mpp_deallocate_domainUG(unstructured_domain)

       !Deallocate the 2D structured domain.
        call mpp_deallocate_domain(structured_domain)

       !Print out a message that the routine is done.
        call mpp_sync()
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*) "Domains destroyed."
            write(output_unit,*)
        endif

        return
    end subroutine destroy_mpp_domains

   !---------------------------------------------------------------------------
   !Test 1:
   !Register axes and fields to a restart file, and then save the restart
   !file.  After the restart file is saved, reset the data and then restore
   !the state (i.e., read in the data from the written restart file).  Make
   !sure that the restored data matches the written data.
    subroutine test_1(unstructured_domain, &
                      num_restarts, &
                      nx, &
                      ny, &
                      nz, &
                      npes)
        use, intrinsic :: iso_fortran_env, only: output_unit
        use mpp_parameter_mod,             only: FATAL, &
                                                 COMM_TAG_1, &
                                                 EVENT_RECV
        use mpp_mod,                       only: mpp_error, &
                                                 mpp_pe, &
                                                 mpp_root_pe, &
                                                 mpp_sync, &
                                                 mpp_chksum, &
                                                 mpp_send, &
                                                 mpp_recv, &
                                                 mpp_sync_self
        use mpp_domains_mod,               only: domainUG, &
                                                 mpp_get_UG_io_domain, &
                                                 mpp_get_UG_domain_npes, &
                                                 mpp_get_UG_domain_pelist
        use mpp_io_mod,                    only: mpp_close
        use fms_io_mod,                    only: restart_file_type, &
                                                 fms_io_unstructured_register_restart_axis, &
                                                 fms_io_unstructured_register_restart_field, &
                                                 CIDX, &
                                                 HIDX, &
                                                 ZIDX, &
                                                 CCIDX, &
                                                 fms_io_unstructured_save_restart, &
                                                 fms_io_unstructured_get_field_size, &
                                                 fms_io_unstructured_read, &
                                                 fms_io_unstructured_file_unit

       !Inputs/Ouputs
        type(domainUG),intent(in)    :: unstructured_domain !<An unstructured mpp domain.
        integer(INT_KIND),intent(in) :: num_restarts        !<Number of times to "restart" the "model run".
        integer(INT_KIND),intent(in) :: nx                  !<The number of grid points in the x-direction.
        integer(INT_KIND),intent(in) :: ny                  !<The number of grid points in the y-direction.
        integer(INT_KIND),intent(in) :: nz                  !<The number of grid points in the z-direction.
        integer(INT_KIND),intent(in) :: npes                !<The total number of ranks used in this test.

       !Local variables
        type(domainUG),pointer                     :: io_domain                                !<Pointer to unstructured domain's I/O domain.
        integer(INT_KIND)                          :: io_domain_npes                           !<The total number of ranks in the unstructured I/O domain pelist.
        integer(INT_KIND),dimension(:),allocatable :: pelist                                   !<A pelist.
        character(len=256)                         :: restart_file_name                        !<Name for the restart file.
        real,dimension(:),allocatable              :: x_axis_data                              !<Data for the x-axis that is registered to the restart file.
        real,dimension(:),allocatable              :: y_axis_data                              !<Data for the y-axis that is registered to the restart file.
        real,dimension(:),allocatable              :: z_axis_data                              !<Data for the z-axis that is registered to the restart file.
        integer(INT_KIND)                          :: cc_axis_size                             !<Size of the cc-axis (???).
        real,dimension(:),allocatable              :: cc_axis_data                             !<Data for the cc-axis (???) that is registered to the restart file.
        integer(INT_KIND)                          :: compressed_c_axis_size                   !<Size of the compressed c (???) axis.:
        integer(INT_KIND),dimension(:),allocatable :: compressed_c_axis_data                   !<Data that is registered to the restart file for the compressed c (???) axis.
        integer(INT_KIND)                          :: compressed_h_axis_size                   !<Size of the compressed c (???) axis.:
        integer(INT_KIND),dimension(:),allocatable :: compressed_h_axis_data                   !<Data that is registered to the restart file for the compressed c (???) axis.
        integer(INT_KIND),dimension(:),allocatable :: compressed_c_axis_size_per_rank          !<Array of "compressed c" axis sizes for each rank.
        integer(INT_KIND),dimension(:),allocatable :: compressed_h_axis_size_per_rank          !<Array of "compressed h" axis sizes for each rank.
        type(restart_file_type)                    :: restart_file                             !<A restart file.
        integer(INT_KIND)                          :: register_id                              !<Id returned from a registered field.
        character(len=256)                         :: real_scalar_field_name                   !<Name for a real scalar field.
        real                                       :: real_scalar_field_data                   !<Data for a real scalar field.
        character(len=256)                         :: compressed_c_real_1D_field_name          !<Name for a "compressed c" real 1D field.
        real,dimension(:),allocatable              :: compressed_c_real_1D_field_data          !<Data for a "compressed c" real 1D field.
        character(len=256)                         :: compressed_h_real_1D_field_name          !<Name for a "compressed h" real 1D field.
        real,dimension(:),allocatable              :: compressed_h_real_1D_field_data          !<Data for a "compressed h" real 1D field.
        character(len=256)                         :: compressed_c_z_real_2D_field_name        !<Name for a "compressed c, z" real 2D field.
        real,dimension(:,:),allocatable            :: compressed_c_z_real_2D_field_data        !<Data for a "compressed c, z" real 2D field.
        character(len=256)                         :: compressed_h_z_real_2D_field_name        !<Name for a "compressed h, z" real 2D field.
        real,dimension(:,:),allocatable            :: compressed_h_z_real_2D_field_data        !<Data for a "compressed h, z" real 2D field.
        character(len=256)                         :: compressed_c_cc_real_2D_field_name       !<Name for a "compressed c, cc" real 2D field.
        real,dimension(:,:),allocatable            :: compressed_c_cc_real_2D_field_data       !<Data for a "compressed c, cc" real 2D field.
        character(len=256)                         :: compressed_h_cc_real_2D_field_name       !<Name for a "compressed h, cc" real 2D field.
        real,dimension(:,:),allocatable            :: compressed_h_cc_real_2D_field_data       !<Data for a "compressed h, cc" real 2D field.
        character(len=256)                         :: compressed_c_z_cc_real_3D_field_name     !<Name for a "compressed c, z, cc" real 3D field.
        real,dimension(:,:,:),allocatable          :: compressed_c_z_cc_real_3D_field_data     !<Data for a "compressed c, z, cc" real 3D field.
        character(len=256)                         :: compressed_h_z_cc_real_3D_field_name     !<Name for a "compressed h, z, cc" real 3D field.
        real,dimension(:,:,:),allocatable          :: compressed_h_z_cc_real_3D_field_data     !<Data for a "compressed h, z, cc" real 3D field.
        character(len=256)                         :: compressed_c_cc_z_real_3D_field_name     !<Name for a "compressed c, cc, z" real 3D field.
        real,dimension(:,:,:),allocatable          :: compressed_c_cc_z_real_3D_field_data     !<Data for a "compressed c, cc, z" real 3D field.
        character(len=256)                         :: compressed_h_cc_z_real_3D_field_name     !<Name for a "compressed h, cc, z" real 3D field.
        real,dimension(:,:,:),allocatable          :: compressed_h_cc_z_real_3D_field_data     !<Data for a "compressed h, cc, z" real 3D field.
        character(len=256)                         :: int_scalar_field_name                    !<Name for an integer scalar field.
        integer                                    :: int_scalar_field_data                    !<Data for an integer scalar field.
        character(len=256)                         :: compressed_c_int_1D_field_name           !<Name for a "compressed c" integer 1D field.
        integer,dimension(:),allocatable           :: compressed_c_int_1D_field_data           !<Data for a "compressed c" integer 1D field.
        character(len=256)                         :: compressed_h_int_1D_field_name           !<Name for a "compressed h" integer 1D field.
        integer,dimension(:),allocatable           :: compressed_h_int_1D_field_data           !<Data for a "compressed h" integer 1D field.
        character(len=256)                         :: compressed_c_z_int_2D_field_name         !<Name for a "compressed c, z" integer 2D field.
        integer,dimension(:,:),allocatable         :: compressed_c_z_int_2D_field_data         !<Data for a "compressed c, z" integer 2D field.
        character(len=256)                         :: compressed_h_z_int_2D_field_name         !<Name for a "compressed h, z" integer 2D field.
        integer,dimension(:,:),allocatable         :: compressed_h_z_int_2D_field_data         !<Data for a "compressed h, z" integer 2D field.
        character(len=256)                         :: compressed_c_cc_int_2D_field_name        !<Name for a "compressed c, cc" integer 2D field.
        integer,dimension(:,:),allocatable         :: compressed_c_cc_int_2D_field_data        !<Data for a "compressed c, cc" integer 2D field.
        character(len=256)                         :: compressed_h_cc_int_2D_field_name        !<Name for a "compressed h, cc" integer 2D field.
        integer,dimension(:,:),allocatable         :: compressed_h_cc_int_2D_field_data        !<Data for a "compressed h, cc" integer 2D field.
        real,dimension(:),allocatable              :: real_scalar_field_data_ref               !<Reference test data for a real scalar field.
        real,dimension(:,:),allocatable            :: compressed_c_real_1D_field_data_ref      !<Reference test data for a "compressed c" real 1D field.
        real,dimension(:,:),allocatable            :: compressed_h_real_1D_field_data_ref      !<Reference test data for a "compressed c" real 1D field.
        real,dimension(:,:,:),allocatable          :: compressed_c_z_real_2D_field_data_ref    !<Reference test data for a "compressed c, z" real 2D field.
        real,dimension(:,:,:),allocatable          :: compressed_h_z_real_2D_field_data_ref    !<Reference test data for a "compressed h, z" real 2D field.
        real,dimension(:,:,:),allocatable          :: compressed_c_cc_real_2D_field_data_ref   !<Reference test data for a "compressed c, cc" real 2D field.
        real,dimension(:,:,:),allocatable          :: compressed_h_cc_real_2D_field_data_ref   !<Reference test data for a "compressed h, cc" real 2D field.
        real,dimension(:,:,:,:),allocatable        :: compressed_c_z_cc_real_3D_field_data_ref !<Reference test data for a compressed real 3D field.
        real,dimension(:,:,:,:),allocatable        :: compressed_h_z_cc_real_3D_field_data_ref !<Reference test data for a compressed real 3D field.
        real,dimension(:,:,:,:),allocatable        :: compressed_c_cc_z_real_3D_field_data_ref !<Reference test data for a compressed real 3D field.
        real,dimension(:,:,:,:),allocatable        :: compressed_h_cc_z_real_3D_field_data_ref !<Reference test data for a compressed real 3D field.
        integer,dimension(:),allocatable           :: int_scalar_field_data_ref                !<Reference test data for an integer scalar field.
        integer,dimension(:,:),allocatable         :: compressed_c_int_1D_field_data_ref       !<Reference test data for a "compressed c" integer 1D field.
        integer,dimension(:,:),allocatable         :: compressed_h_int_1D_field_data_ref       !<Reference test data for a "compressed c" integer 1D field.
        integer,dimension(:,:,:),allocatable       :: compressed_c_z_int_2D_field_data_ref     !<Reference test data for a "compressed c, z" integer 2D field.
        integer,dimension(:,:,:),allocatable       :: compressed_h_z_int_2D_field_data_ref     !<Reference test data for a "compressed h, z" integer 2D field.
        integer,dimension(:,:,:),allocatable       :: compressed_c_cc_int_2D_field_data_ref    !<Reference test data for a "compressed c, cc" integer 2D field.
        integer,dimension(:,:,:),allocatable       :: compressed_h_cc_int_2D_field_data_ref    !<Reference test data for a "compressed h, cc" integer 2D field.
        integer(INT_KIND),dimension(5)             :: field_dimension_sizes                    !<Array to hold the dimensions of fields when they are read back in.
        logical(INT_KIND)                          :: field_found_in_file                      !<Flag telling if a field was found in a file.
        real,dimension(:),allocatable              :: real_buffer_1D                           !<Buffer used to read back in the data.
        real,dimension(:,:),allocatable            :: real_buffer_2D                           !<Buffer used to read back in the data.
        real,dimension(:,:,:),allocatable          :: real_buffer_3D                           !<Buffer used to read back in the data.
        integer,dimension(:),allocatable           :: int_buffer_1D                            !<Buffer used to read back in the data.
        integer,dimension(:,:),allocatable         :: int_buffer_2D                            !<Buffer used to read back in the data.
        integer(INT_KIND)                          :: offset                                   !<Offset used to check the read in data.
        real                                       :: rmax_error                               !<Error for real data.
        integer                                    :: imax_error                               !<Error for integer data.
        integer(LONG_KIND)                         :: read_in_chksum                           !<Check-sum for read in data.
        integer(LONG_KIND)                         :: ref_chksum                               !<Check-sum for reference test data.
        integer(INT_KIND)                          :: funit                                    !<File unit used to close a file.
        integer(INT_KIND)                          :: i                                        !<Loop variable.
        integer(INT_KIND)                          :: j                                        !<Loop variable.
        integer(INT_KIND)                          :: k                                        !<Loop variable.
        integer(INT_KIND)                          :: q                                        !<Loop variable.

       !Print out a message that the test is starting.
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*)
            write(output_unit,*) "Executing test 1 body ..."
        endif

       !Make sure that the inputted number of restarts is less than or
       !equal to 10, to keep the size of the output files reasonable.
        if (num_restarts .gt. 10) then
            call mpp_error(FATAL, &
                           "test 1: the inputted number of restarts should" &
                           //" be less than or equal to 10.")
        endif

       !Synchronize all ranks.
        call mpp_sync()

       !Get the pelist associated with the unstructured I/O domain.
        io_domain => null()
        io_domain => mpp_get_UG_io_domain(unstructured_domain)
        io_domain_npes = mpp_get_UG_domain_npes(io_domain)
        allocate(pelist(io_domain_npes))
        call mpp_get_UG_domain_pelist(io_domain, &
                                      pelist)
        io_domain => null()

       !Set the name of the restart file.
        restart_file_name = "test_1_unstructured_restart_file.nc"

       !Add a x-axis to the restart file.  This must be done before any fields
       !are registered.
        allocate(x_axis_data(nx))
        do i = 1,nx
            x_axis_data(i) = (real(i-1))*(360.0/real(nx))
        enddo
        call fms_io_unstructured_register_restart_axis(restart_file, &
                                                       restart_file_name, &
                                                       "lon", &
                                                       x_axis_data, &
                                                       "X", &
                                                       unstructured_domain)

       !Add a y-axis to the restart file.  This must be done before any fields
       !are registered.
        allocate(y_axis_data(ny))
        do i = 1,ny
            y_axis_data(i) = (real(i-1))*(180.0/real(ny))
        enddo
        call fms_io_unstructured_register_restart_axis(restart_file, &
                                                       restart_file_name, &
                                                       "lat", &
                                                       y_axis_data, &
                                                       "Y", &
                                                       unstructured_domain)

       !Add a z-axis to the restart file.  This must be done before any fields
       !are registered.
        allocate(z_axis_data(nz))
        do i = 1,nz
            z_axis_data(i) = real(i*5.0)
        enddo
        call fms_io_unstructured_register_restart_axis(restart_file, &
                                                       restart_file_name, &
                                                       "height", &
                                                       z_axis_data, &
                                                       "Z", &
                                                       unstructured_domain)

       !Add a cc-axis (???)to the restart file.  This must be done before any
       !fields are registered.
        cc_axis_size = 9
        allocate(cc_axis_data(cc_axis_size))
        do i = 1,cc_axis_size
            cc_axis_data(i) = real(i)*12.0 - real(i)
        enddo
        call fms_io_unstructured_register_restart_axis(restart_file, &
                                                       restart_file_name, &
                                                       "cc_axis", &
                                                       cc_axis_data, &
                                                       "CC", &
                                                       unstructured_domain)

       !Add a "C" type "compressed" axis to the restart file. This must be
       !done before any fields are registered.  Get the "compressed c" axis
       !size from all other ranks on the same I/O domain pelist.  This is
       !needed to check to data after it is read back in.
        compressed_c_axis_size = mpp_pe() + 1
        allocate(compressed_c_axis_data(compressed_c_axis_size))
        do i = 1,compressed_c_axis_size
            compressed_c_axis_data(i) = real(mpp_pe()*mpp_pe() + i)
        enddo
        call fms_io_unstructured_register_restart_axis(restart_file, &
                                                       restart_file_name, &
                                                       "compressed_c_axis", &
                                                       compressed_c_axis_data, &
                                                       "compressed c", &
                                                       "C", &
                                                       compressed_c_axis_size, &
                                                       unstructured_domain)
        allocate(compressed_c_axis_size_per_rank(io_domain_npes))
        compressed_c_axis_size_per_rank = 0
        do i = 1,io_domain_npes
            if (pelist(i) .ne. mpp_pe()) then
                call mpp_recv(compressed_c_axis_size_per_rank(i), &
                              pelist(i), &
                              block = .false., &
                              tag=COMM_TAG_1)
                call mpp_send(compressed_c_axis_size, &
                              pelist(i), &
                              tag=COMM_TAG_1)
            else
                compressed_c_axis_size_per_rank(i) = compressed_c_axis_size
            endif
        enddo
        call mpp_sync_self(check=EVENT_RECV)
        call mpp_sync_self()
        call mpp_sync()

       !Add a "H" type "compressed" axis to the restart file. This must be
       !done before any fields are registered.
        compressed_h_axis_size = npes
        allocate(compressed_h_axis_data(compressed_h_axis_size))
        do i = 1,compressed_h_axis_size
            compressed_h_axis_data(i) = real((mpp_pe())*npes + i)
        enddo
        call fms_io_unstructured_register_restart_axis(restart_file, &
                                                       restart_file_name, &
                                                       "compressed_h_axis", &
                                                       compressed_h_axis_data, &
                                                       "compressed h", &
                                                       "H", &
                                                       compressed_h_axis_size, &
                                                       unstructured_domain)
        allocate(compressed_h_axis_size_per_rank(io_domain_npes))
        compressed_h_axis_size_per_rank = 0
        do i = 1,io_domain_npes
            if (pelist(i) .ne. mpp_pe()) then
                call mpp_recv(compressed_h_axis_size_per_rank(i), &
                              pelist(i), &
                              block = .false., &
                              tag=COMM_TAG_1)
                call mpp_send(compressed_h_axis_size, &
                              pelist(i), &
                              tag=COMM_TAG_1)
            else
                compressed_h_axis_size_per_rank(i) = compressed_h_axis_size
            endif
        enddo
        call mpp_sync_self(check=EVENT_RECV)
        call mpp_sync_self()
        call mpp_sync()

       !Create a real scalar field and register it to the restart file.  If
       !a scalar field value is different across ranks, only the value on the
       !root of the I/O domain pelist gets written to the file.
       !Should we check for this?  Should it be a fatal error?
        real_scalar_field_name = "real_scalar_field_1"
        real_scalar_field_data = 1234.5678
!       real_scalar_field_data = 1234.5678 + real(mpp_pe()*1000)
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 real_scalar_field_name, &
                                                                 real_scalar_field_data, &
                                                                 unstructured_domain, &
                                                                 longname="rsf1", &
                                                                 units="ergs")

       !Create a "compressed c" real 1D field and register it to the restart
       !file.  This field is of the form:
       !field = field(compressed c).
        compressed_c_real_1D_field_name = "compressed_c_real_1D_field_1"
        allocate(compressed_c_real_1D_field_data(compressed_c_axis_size))
        do i = 1,compressed_c_axis_size
            compressed_c_real_1D_field_data(i) = real(mpp_pe())
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_real_1D_field_name, &
                                                                 compressed_c_real_1D_field_data, &
                                                                 (/CIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r1Dcompcf1", &
                                                                 units="cm")

       !Create a "compressed h" real 1D field and register it to the restart
       !file.  This field is of the form:
       !field = field(compressed h).
        compressed_h_real_1D_field_name = "compressed_h_real_1D_field_1"
        allocate(compressed_h_real_1D_field_data(compressed_h_axis_size))
        do i = 1,compressed_h_axis_size
            compressed_h_real_1D_field_data(i) = real(i) + real(mpp_pe()) &
                                                 + 1111.1111
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_real_1D_field_name, &
                                                                 compressed_h_real_1D_field_data, &
                                                                 (/HIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r1Dcomphf1", &
                                                                 units="km")

       !Create a "compressed c, z" real 2D field and register it to the restart
       !file.  The field is of the form:
       !field = field(compressed c,z).
        compressed_c_z_real_2D_field_name = "compressed_c_z_real_2D_field_1"
        allocate(compressed_c_z_real_2D_field_data(compressed_c_axis_size,nz))
        do j = 1,nz
            do i = 1,compressed_c_axis_size
                compressed_c_z_real_2D_field_data(i,j) = real(mpp_pe()*1000) &
                                                         + real(100*j) &
                                                         + real(10*i)
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_z_real_2D_field_name, &
                                                                 compressed_c_z_real_2D_field_data, &
                                                                 (/CIDX,ZIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r2Dcompczf1", &
                                                                 units="N")

       !Create a "compressed h, z" real 2D field and register it to the restart
       !file.  The field is of the form:
       !field = field(compressed h,z).
        compressed_h_z_real_2D_field_name = "compressed_h_z_real_2D_field_1"
        allocate(compressed_h_z_real_2D_field_data(compressed_h_axis_size,nz))
        do j = 1,nz
            do i = 1,compressed_h_axis_size
                compressed_h_z_real_2D_field_data(i,j) = real(mpp_pe()*1000) &
                                                         - 1.0*real((j-1)* &
                                                         compressed_h_axis_size+i) &
                                                         + 2.2222222
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_z_real_2D_field_name, &
                                                                 compressed_h_z_real_2D_field_data, &
                                                                 (/HIDX,ZIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r2Dcomphzf1", &
                                                                 units="kN")

       !Create a "compressed c, cc" real 2D field and register it to the restart
       !file.  The field is of the form:
       !field = field(compressed c,cc).
        compressed_c_cc_real_2D_field_name = "compressed_c_cc_real_2D_field_1"
        allocate(compressed_c_cc_real_2D_field_data(compressed_c_axis_size, &
                                                    cc_axis_size))
        do j = 1,cc_axis_size
            do i = 1,compressed_c_axis_size
                compressed_c_cc_real_2D_field_data(i,j) = real(mpp_pe()*1111) &
                                                         + real(111*j) &
                                                         + real(11*i)
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_cc_real_2D_field_name, &
                                                                 compressed_c_cc_real_2D_field_data, &
                                                                 (/CIDX,CCIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r2Dcompcccf1", &
                                                                 units="T")

       !Create a "compressed h, z" real 2D field and register it to the restart
       !file.  The field is of the form:
       !field = field(compressed h,cc).
        compressed_h_cc_real_2D_field_name = "compressed_h_cc_real_2D_field_1"
        allocate(compressed_h_cc_real_2D_field_data(compressed_h_axis_size, &
                                                    cc_axis_size))
        do j = 1,cc_axis_size
            do i = 1,compressed_h_axis_size
                compressed_h_cc_real_2D_field_data(i,j) = real(mpp_pe()*1111) &
                                                          - 5.0*real((j-1)* &
                                                          compressed_h_axis_size+i) &
                                                          + 2.2222222
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_cc_real_2D_field_name, &
                                                                 compressed_h_cc_real_2D_field_data, &
                                                                 (/HIDX,CCIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r2Dcomphccf1", &
                                                                 units="kT")

       !Create a "compressed c, z, cc" real 3D field and register it to the
       !restart file.  This field is of the form:
       !field = field(compressed c,z,cc).
        compressed_c_z_cc_real_3D_field_name = "compressed_c_z_cc_real_3D_field_1"
        allocate(compressed_c_z_cc_real_3D_field_data(compressed_c_axis_size, &
                                                      nz, &
                                                      cc_axis_size))
        do k = 1,cc_axis_size
            do j = 1,nz
                do i = 1,compressed_c_axis_size
                    compressed_c_z_cc_real_3D_field_data(i,j,k) = real(mpp_pe()*10000) &
                                                                  + real(k*1000) &
                                                                  + real(j*100) &
                                                                  + real(i*10)
                enddo
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_z_cc_real_3D_field_name, &
                                                                 compressed_c_z_cc_real_3D_field_data, &
                                                                 (/CIDX,ZIDX,CCIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r3Dcompczccf1", &
                                                                 units="nm")

       !Create a "compressed h, z, cc" real 3D field and register it to the
       !restart file.  This field is of the form:
       !field = field(compressed h,z,cc).
        compressed_h_z_cc_real_3D_field_name = "compressed_h_z_cc_real_3D_field_1"
        allocate(compressed_h_z_cc_real_3D_field_data(compressed_h_axis_size, &
                                                      nz, &
                                                      cc_axis_size))
        do k = 1,cc_axis_size
            do j = 1,nz
                do i = 1,compressed_h_axis_size
                    compressed_h_z_cc_real_3D_field_data(i,j,k) = real(mpp_pe()*1000) &
                                                                  - 1.0*real((j-1)* &
                                                                  compressed_h_axis_size+i) &
                                                                  + 2.2222222*(real(k-mpp_pe()))
                enddo
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_z_cc_real_3D_field_name, &
                                                                 compressed_h_z_cc_real_3D_field_data, &
                                                                 (/HIDX,ZIDX,CCIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r3Dcomphzccf1", &
                                                                 units="mB")

       !Create a "compressed c, cc, z" real 3D field and register it to the
       !restart file.  This field is of the form:
       !field = field(compressed c,cc,z).
        compressed_c_cc_z_real_3D_field_name = "compressed_c_cc_z_real_3D_field_1"
        allocate(compressed_c_cc_z_real_3D_field_data(compressed_c_axis_size, &
                                                      cc_axis_size, &
                                                      nz))
        do k = 1,nz
            do j = 1,cc_axis_size
                do i = 1,compressed_c_axis_size
                    compressed_c_cc_z_real_3D_field_data(i,j,k) = real(mpp_pe()*11111) &
                                                                  + real(k*1000) &
                                                                  + real(j*100) &
                                                                  + real(i*10)
                enddo
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_cc_z_real_3D_field_name, &
                                                                 compressed_c_cc_z_real_3D_field_data, &
                                                                 (/CIDX,CCIDX,ZIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r3Dcompccczf1", &
                                                                 units="yr")

       !Create a "compressed h, cc, z" real 3D field and register it to the
       !restart file.  This field is of the form:
       !field = field(compressed h,cc,z).
        compressed_h_cc_z_real_3D_field_name = "compressed_h_cc_z_real_3D_field_1"
        allocate(compressed_h_cc_z_real_3D_field_data(compressed_h_axis_size, &
                                                      cc_axis_size, &
                                                      nz))
        do k = 1,nz
            do j = 1,cc_axis_size
                do i = 1,compressed_h_axis_size
                    compressed_h_cc_z_real_3D_field_data(i,j,k) = real(mpp_pe()*11111) &
                                                                  - 1.0*real((j-1)* &
                                                                  compressed_h_axis_size+i) &
                                                                  + 2.2222222*(real(k-mpp_pe()))
                enddo
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_cc_z_real_3D_field_name, &
                                                                 compressed_h_cc_z_real_3D_field_data, &
                                                                 (/HIDX,CCIDX,ZIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="r3Dcomphcczf1", &
                                                                 units="hour")

       !Create an integer scalar field and register it to the restart file.  If
       !a scalar field value is different across ranks, only the value on the
       !root of the I/O domain pelist gets written to the file.
       !Should we check for this?  Should it be a fatal error?
        int_scalar_field_name = "int_scalar_field_1"
        int_scalar_field_data = 4321
!       int_scalar_field_data = 4321 + (mpp_pe()*1000)
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 int_scalar_field_name, &
                                                                 int_scalar_field_data, &
                                                                 unstructured_domain, &
                                                                 longname="isf1", &
                                                                 units="Watt")

       !Create a "compressed c" integer 1D field and register it to the restart
       !file.  This field is of the form:
       !field = field(compressed c).
        compressed_c_int_1D_field_name = "compressed_c_int_1D_field_1"
        allocate(compressed_c_int_1D_field_data(compressed_c_axis_size))
        do i = 1,compressed_c_axis_size
            compressed_c_int_1D_field_data(i) = mpp_pe()
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_int_1D_field_name, &
                                                                 compressed_c_int_1D_field_data, &
                                                                 (/CIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="i1Dcompcf1", &
                                                                 units="pc")

       !Create a "compressed h" integer 1D field and register it to the restart
       !file.  This field is of the form:
       !field = field(compressed h).
        compressed_h_int_1D_field_name = "compressed_h_int_1D_field_1"
        allocate(compressed_h_int_1D_field_data(compressed_h_axis_size))
        do i = 1,compressed_h_axis_size
            compressed_h_int_1D_field_data(i) = i + mpp_pe() + 1111
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_int_1D_field_name, &
                                                                 compressed_h_int_1D_field_data, &
                                                                 (/HIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="i1Dcomphf1", &
                                                                 units="kpc")

       !Create a "compressed c, z" integer 2D field and register it to the
       !restart file.  The field is of the form:
       !field = field(compressed c,z).
        compressed_c_z_int_2D_field_name = "compressed_c_z_int_2D_field_1"
        allocate(compressed_c_z_int_2D_field_data(compressed_c_axis_size,nz))
        do j = 1,nz
            do i = 1,compressed_c_axis_size
                compressed_c_z_int_2D_field_data(i,j) = mpp_pe()*1000 &
                                                        + 100*j + 10
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_z_int_2D_field_name, &
                                                                 compressed_c_z_int_2D_field_data, &
                                                                 (/CIDX,ZIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="i2Dcompczf1", &
                                                                 units="au")

       !Create a "compressed h, z" integer 2D field and register it to the
       !restart file.  The field is of the form:
       !field = field(compressed h,z).
        compressed_h_z_int_2D_field_name = "compressed_h_z_int_2D_field_1"
        allocate(compressed_h_z_int_2D_field_data(compressed_h_axis_size,nz))
        do j = 1,nz
            do i = 1,compressed_h_axis_size
                compressed_h_z_int_2D_field_data(i,j) = mpp_pe()*1000 &
                                                        - 1.0*(j-1)* &
                                                         compressed_h_axis_size+i &
                                                         + 2
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_z_int_2D_field_name, &
                                                                 compressed_h_z_int_2D_field_data, &
                                                                 (/HIDX,ZIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="i2Dcomphzf1", &
                                                                 units="kau")

       !Create a "compressed c, cc" integer 2D field and register it to the
       !restart file.  The field is of the form:
       !field = field(compressed c,cc).
        compressed_c_cc_int_2D_field_name = "compressed_c_cc_int_2D_field_1"
        allocate(compressed_c_cc_int_2D_field_data(compressed_c_axis_size, &
                                                   cc_axis_size))
        do j = 1,cc_axis_size
            do i = 1,compressed_c_axis_size
                compressed_c_cc_int_2D_field_data(i,j) = mpp_pe()*1111 &
                                                         + 111*j + 11*i
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_c_cc_int_2D_field_name, &
                                                                 compressed_c_cc_int_2D_field_data, &
                                                                 (/CIDX,CCIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="i2Dcompcccf1", &
                                                                 units="mol")

       !Create a "compressed h, z" integer 2D field and register it to the
       !restart file.  The field is of the form:
       !field = field(compressed h,cc).
        compressed_h_cc_int_2D_field_name = "compressed_h_cc_int_2D_field_1"
        allocate(compressed_h_cc_int_2D_field_data(compressed_h_axis_size, &
                                                   cc_axis_size))
        do j = 1,cc_axis_size
            do i = 1,compressed_h_axis_size
                compressed_h_cc_int_2D_field_data(i,j) = mpp_pe()*1111 &
                                                         - 5.0*(j-1)* &
                                                         compressed_h_axis_size+i &
                                                         + 2
            enddo
        enddo
        register_id = fms_io_unstructured_register_restart_field(restart_file, &
                                                                 restart_file_name, &
                                                                 compressed_h_cc_int_2D_field_name, &
                                                                 compressed_h_cc_int_2D_field_data, &
                                                                 (/HIDX,CCIDX/), &
                                                                 unstructured_domain, &
                                                                 longname="i2Dcomphccf1", &
                                                                 units="kmol")

       !Allocate arrays used to check the results on a "restart".
        allocate(real_scalar_field_data_ref(num_restarts))
        allocate(compressed_c_real_1D_field_data_ref(compressed_c_axis_size,num_restarts))
        allocate(compressed_h_real_1D_field_data_ref(compressed_h_axis_size,num_restarts))
        allocate(compressed_c_z_real_2D_field_data_ref(compressed_c_axis_size,nz,num_restarts))
        allocate(compressed_h_z_real_2D_field_data_ref(compressed_h_axis_size,nz,num_restarts))
        allocate(compressed_c_cc_real_2D_field_data_ref(compressed_c_axis_size,cc_axis_size,num_restarts))
        allocate(compressed_h_cc_real_2D_field_data_ref(compressed_h_axis_size,cc_axis_size,num_restarts))
        allocate(compressed_c_z_cc_real_3D_field_data_ref(compressed_c_axis_size,nz,cc_axis_size,num_restarts))
        allocate(compressed_h_z_cc_real_3D_field_data_ref(compressed_h_axis_size,nz,cc_axis_size,num_restarts))
        allocate(compressed_c_cc_z_real_3D_field_data_ref(compressed_c_axis_size,cc_axis_size,nz,num_restarts))
        allocate(compressed_h_cc_z_real_3D_field_data_ref(compressed_h_axis_size,cc_axis_size,nz,num_restarts))
        allocate(int_scalar_field_data_ref(num_restarts))
        allocate(compressed_c_int_1D_field_data_ref(compressed_c_axis_size,num_restarts))
        allocate(compressed_h_int_1D_field_data_ref(compressed_h_axis_size,num_restarts))
        allocate(compressed_c_z_int_2D_field_data_ref(compressed_c_axis_size,nz,num_restarts))
        allocate(compressed_h_z_int_2D_field_data_ref(compressed_h_axis_size,nz,num_restarts))
        allocate(compressed_c_cc_int_2D_field_data_ref(compressed_c_axis_size,cc_axis_size,num_restarts))
        allocate(compressed_h_cc_int_2D_field_data_ref(compressed_h_axis_size,cc_axis_size,num_restarts))
        call mpp_sync()

       !Write out the restart data at an inputted number of time levels.
        do q = 1,num_restarts

           !Simulate the data evolving in time.
            if (q .gt. 1) then
                real_scalar_field_data = real_scalar_field_data + 8.0
                compressed_c_real_1D_field_data = compressed_c_real_1D_field_data + 9.0
                compressed_h_real_1D_field_data = compressed_h_real_1D_field_data + 10.0
                compressed_c_z_real_2D_field_data = compressed_c_z_real_2D_field_data + 11.0
                compressed_h_z_real_2D_field_data = compressed_h_z_real_2D_field_data + 12.0
                compressed_c_cc_real_2D_field_data = compressed_c_cc_real_2D_field_data + 13.0
                compressed_h_cc_real_2D_field_data = compressed_h_cc_real_2D_field_data + 14.0
                compressed_c_z_cc_real_3D_field_data = compressed_c_z_cc_real_3D_field_data + 15.0
                compressed_h_z_cc_real_3D_field_data = compressed_h_z_cc_real_3D_field_data + 16.0
                compressed_c_cc_z_real_3D_field_data = compressed_c_cc_z_real_3D_field_data + 17.0
                compressed_h_cc_z_real_3D_field_data = compressed_h_cc_z_real_3D_field_data + 18.0
                int_scalar_field_data = int_scalar_field_data + 19
                compressed_c_int_1D_field_data = compressed_c_int_1D_field_data + 20
                compressed_h_int_1D_field_data = compressed_h_int_1D_field_data + 21
                compressed_c_z_int_2D_field_data = compressed_c_z_int_2D_field_data + 22
                compressed_h_z_int_2D_field_data = compressed_h_z_int_2D_field_data + 23
                compressed_c_cc_int_2D_field_data = compressed_c_cc_int_2D_field_data + 24
                compressed_h_cc_int_2D_field_data = compressed_h_cc_int_2D_field_data + 25
            endif

           !Save the state that will be written to the restarts.
            real_scalar_field_data_ref(q) = real_scalar_field_data
            compressed_c_real_1D_field_data_ref(:,q) = compressed_c_real_1D_field_data
            compressed_h_real_1D_field_data_ref(:,q) = compressed_h_real_1D_field_data
            compressed_c_z_real_2D_field_data_ref(:,:,q) = compressed_c_z_real_2D_field_data
            compressed_h_z_real_2D_field_data_ref(:,:,q) = compressed_h_z_real_2D_field_data
            compressed_c_cc_real_2D_field_data_ref(:,:,q) = compressed_c_cc_real_2D_field_data
            compressed_h_cc_real_2D_field_data_ref(:,:,q) = compressed_h_cc_real_2D_field_data
            compressed_c_z_cc_real_3D_field_data_ref(:,:,:,q) = compressed_c_z_cc_real_3D_field_data
            compressed_h_z_cc_real_3D_field_data_ref(:,:,:,q) = compressed_h_z_cc_real_3D_field_data
            compressed_c_cc_z_real_3D_field_data_ref(:,:,:,q) = compressed_c_cc_z_real_3D_field_data
            compressed_h_cc_z_real_3D_field_data_ref(:,:,:,q) = compressed_h_cc_z_real_3D_field_data
            int_scalar_field_data_ref(q) = int_scalar_field_data
            compressed_c_int_1D_field_data_ref(:,q) = compressed_c_int_1D_field_data
            compressed_h_int_1D_field_data_ref(:,q) = compressed_h_int_1D_field_data
            compressed_c_z_int_2D_field_data_ref(:,:,q) = compressed_c_z_int_2D_field_data
            compressed_h_z_int_2D_field_data_ref(:,:,q) = compressed_h_z_int_2D_field_data
            compressed_c_cc_int_2D_field_data_ref(:,:,q) = compressed_c_cc_int_2D_field_data
            compressed_h_cc_int_2D_field_data_ref(:,:,q) = compressed_h_cc_int_2D_field_data

           !Write out the restart file.
            call mpp_sync()
            if (q .gt. 1) then
                call fms_io_unstructured_save_restart(restart_file, &
                                                      directory="RESTART", &
                                                      append=.true., &
                                                      time_level=real(q))
            else
                call fms_io_unstructured_save_restart(restart_file, &
                                                      directory="RESTART")
            endif
            call mpp_sync()
        enddo

       !For each time level, read the data back in and compare it to the
       !stored reference data.
        do q = 1,num_restarts

           !Write a message specifiying the time level at which the data
           !is being compared.
            if (mpp_pe() .eq. mpp_root_pe()) then
                write(output_unit,*)
                write(output_unit,*) "Checking restart data at timelevel:",q
            endif

           !Zero out all field data to simulate a restart run starting.
            call mpp_sync()
            real_scalar_field_data = 0.0
            compressed_c_real_1D_field_data = 0.0
            compressed_h_real_1D_field_data = 0.0
            compressed_c_z_real_2D_field_data = 0.0
            compressed_h_z_real_2D_field_data = 0.0
            compressed_c_cc_real_2D_field_data = 0.0
            compressed_h_cc_real_2D_field_data = 0.0
            compressed_c_z_cc_real_3D_field_data = 0.0
            compressed_h_z_cc_real_3D_field_data = 0.0
            compressed_c_cc_z_real_3D_field_data = 0.0
            compressed_h_cc_z_real_3D_field_data = 0.0
            int_scalar_field_data = 0
            compressed_c_int_1D_field_data = 0
            compressed_h_int_1D_field_data = 0
            compressed_c_z_int_2D_field_data = 0
            compressed_h_z_int_2D_field_data = 0
            compressed_c_cc_int_2D_field_data = 0
            compressed_h_cc_int_2D_field_data = 0

           !Read the data back in.  The read will read-in an entire I/O domain
           !tile's worth of data.  Each rank is then required to copy its own
           !part of that data back into the appropriate arrays so that the
           !data can be checked.

           !-------------------------------------------------------------------
           !Real scalar field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    real_scalar_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(real_scalar_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_1D(1))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          real_scalar_field_name, &
                                          real_buffer_1D, &
                                          unstructured_domain, &
                                          timelevel=q)
            real_scalar_field_data = real_buffer_1D(1)
            deallocate(real_buffer_1D)
            rmax_error = abs(real_scalar_field_data - &
                             real_scalar_field_data_ref(q))
            if (rmax_error .ne. 0.0) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                     ",mpp_pe()
                write(output_unit,*) "Restart iteration:           ",q
                write(output_unit,*) "Max real scalar field error: ",rmax_error
                call mpp_error(FATAL, &
                               "test 1: real scalar field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,*) "Real scalar field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c real 1D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_real_1D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_real_1D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_1D(field_dimension_sizes(1)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_real_1D_field_name, &
                                          real_buffer_1D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do i = 1,compressed_c_axis_size
                compressed_c_real_1D_field_data(i) = real_buffer_1D(i+offset)
            enddo
            deallocate(real_buffer_1D)
            rmax_error = maxval(abs(compressed_c_real_1D_field_data - &
                                    compressed_c_real_1D_field_data_ref(:,q)))
            read_in_chksum = mpp_chksum(compressed_c_real_1D_field_data)
            ref_chksum = mpp_chksum(compressed_c_real_1D_field_data_ref(:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                              ",mpp_pe()
                write(output_unit,*) "Restart iteration:                    ",q
                write(output_unit,*) "Max compressed c real 1D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:             ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:           ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c real 1D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c real 1D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h real 1D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_real_1D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_real_1D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_1D(field_dimension_sizes(1)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_real_1D_field_name, &
                                          real_buffer_1D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do i = 1,compressed_h_axis_size
                compressed_h_real_1D_field_data(i) = real_buffer_1D(i+offset)
            enddo
            deallocate(real_buffer_1D)
            rmax_error = maxval(abs(compressed_h_real_1D_field_data - &
                                    compressed_h_real_1D_field_data_ref(:,q)))
            read_in_chksum = mpp_chksum(compressed_h_real_1D_field_data)
            ref_chksum = mpp_chksum(compressed_h_real_1D_field_data_ref(:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                              ",mpp_pe()
                write(output_unit,*) "Restart iteration:                    ",q
                write(output_unit,*) "Max compressed h real 1D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:             ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:           ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h real 1D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h real 1D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c, z real 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_z_real_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_z_real_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_2D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_z_real_2D_field_name, &
                                          real_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_c_z_real_2D_field_data,2)
                do i = 1,size(compressed_c_z_real_2D_field_data,1)
                    compressed_c_z_real_2D_field_data(i,j) = real_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(real_buffer_2D)
            rmax_error = maxval(abs(compressed_c_z_real_2D_field_data - &
                                    compressed_c_z_real_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_c_z_real_2D_field_data)
            ref_chksum = mpp_chksum(compressed_c_z_real_2D_field_data_ref(:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                 ",mpp_pe()
                write(output_unit,*) "Restart iteration:                       ",q
                write(output_unit,*) "Max compressed c, z real 2D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:              ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c, z real 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c, z real 2D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h, z real 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_z_real_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_z_real_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_2D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_z_real_2D_field_name, &
                                          real_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_h_z_real_2D_field_data,2)
                do i = 1,size(compressed_h_z_real_2D_field_data,1)
                    compressed_h_z_real_2D_field_data(i,j) = real_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(real_buffer_2D)
            rmax_error = maxval(abs(compressed_h_z_real_2D_field_data - &
                                    compressed_h_z_real_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_h_z_real_2D_field_data)
            ref_chksum = mpp_chksum(compressed_h_z_real_2D_field_data_ref(:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                 ",mpp_pe()
                write(output_unit,*) "Restart iteration:                       ",q
                write(output_unit,*) "Max compressed h, z real 2D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:              ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h, z real 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h, z real 2D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c, cc real 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_cc_real_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_cc_real_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_2D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_cc_real_2D_field_name, &
                                          real_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_c_cc_real_2D_field_data,2)
                do i = 1,size(compressed_c_cc_real_2D_field_data,1)
                    compressed_c_cc_real_2D_field_data(i,j) = real_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(real_buffer_2D)
            rmax_error = maxval(abs(compressed_c_cc_real_2D_field_data - &
                                    compressed_c_cc_real_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_c_cc_real_2D_field_data)
            ref_chksum = mpp_chksum(compressed_c_cc_real_2D_field_data_ref(:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                  ",mpp_pe()
                write(output_unit,*) "Restart iteration:                        ",q
                write(output_unit,*) "Max compressed c, cc real 2D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                 ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:               ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c, z real 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c, cc real 2D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h, cc real 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_cc_real_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_cc_real_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_2D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_cc_real_2D_field_name, &
                                          real_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_h_cc_real_2D_field_data,2)
                do i = 1,size(compressed_h_cc_real_2D_field_data,1)
                    compressed_h_cc_real_2D_field_data(i,j) = real_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(real_buffer_2D)
            rmax_error = maxval(abs(compressed_h_cc_real_2D_field_data - &
                                    compressed_h_cc_real_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_h_cc_real_2D_field_data)
            ref_chksum = mpp_chksum(compressed_h_cc_real_2D_field_data_ref(:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                  ",mpp_pe()
                write(output_unit,*) "Restart iteration:                        ",q
                write(output_unit,*) "Max compressed h, cc real 2D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                 ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:               ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h, z real 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h, cc real 2D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c, z, cc real 3D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_z_cc_real_3D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_z_cc_real_3D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_3D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2), &
                                    field_dimension_sizes(3)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_z_cc_real_3D_field_name, &
                                          real_buffer_3D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do k = 1,size(compressed_c_z_cc_real_3D_field_data,3)
                do j = 1,size(compressed_c_z_cc_real_3D_field_data,2)
                    do i = 1,size(compressed_c_z_cc_real_3D_field_data,1)
                        compressed_c_z_cc_real_3D_field_data(i,j,k) = real_buffer_3D(i+offset,j,k)
                    enddo
                enddo
            enddo
            deallocate(real_buffer_3D)
            rmax_error = maxval(abs(compressed_c_z_cc_real_3D_field_data - &
                                    compressed_c_z_cc_real_3D_field_data_ref(:,:,:,q)))
            read_in_chksum = mpp_chksum(compressed_c_z_cc_real_3D_field_data)
            ref_chksum = mpp_chksum(compressed_c_z_cc_real_3D_field_data_ref(:,:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                     ",mpp_pe()
                write(output_unit,*) "Restart iteration:                           ",q
                write(output_unit,*) "Max compressed c, z, cc real 3D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                    ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                  ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c, z, cc real 3D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c, z, cc real 3D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h, z, cc real 3D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_z_cc_real_3D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_z_cc_real_3D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_3D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2), &
                                    field_dimension_sizes(3)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_z_cc_real_3D_field_name, &
                                          real_buffer_3D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do k = 1,size(compressed_h_z_cc_real_3D_field_data,3)
                do j = 1,size(compressed_h_z_cc_real_3D_field_data,2)
                    do i = 1,size(compressed_h_z_cc_real_3D_field_data,1)
                        compressed_h_z_cc_real_3D_field_data(i,j,k) = real_buffer_3D(i+offset,j,k)
                    enddo
                enddo
            enddo
            deallocate(real_buffer_3D)
            rmax_error = maxval(abs(compressed_h_z_cc_real_3D_field_data - &
                                    compressed_h_z_cc_real_3D_field_data_ref(:,:,:,q)))
            read_in_chksum = mpp_chksum(compressed_h_z_cc_real_3D_field_data)
            ref_chksum = mpp_chksum(compressed_h_z_cc_real_3D_field_data_ref(:,:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                     ",mpp_pe()
                write(output_unit,*) "Restart iteration:                           ",q
                write(output_unit,*) "Max compressed h, z, cc real 3D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                    ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                  ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h, z, cc real 3D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h, z, cc real 3D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c, cc, z real 3D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_cc_z_real_3D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_cc_z_real_3D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_3D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2), &
                                    field_dimension_sizes(3)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_cc_z_real_3D_field_name, &
                                          real_buffer_3D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do k = 1,size(compressed_c_cc_z_real_3D_field_data,3)
                do j = 1,size(compressed_c_cc_z_real_3D_field_data,2)
                    do i = 1,size(compressed_c_cc_z_real_3D_field_data,1)
                        compressed_c_cc_z_real_3D_field_data(i,j,k) = real_buffer_3D(i+offset,j,k)
                    enddo
                enddo
            enddo
            deallocate(real_buffer_3D)
            rmax_error = maxval(abs(compressed_c_cc_z_real_3D_field_data - &
                                    compressed_c_cc_z_real_3D_field_data_ref(:,:,:,q)))
            read_in_chksum = mpp_chksum(compressed_c_cc_z_real_3D_field_data)
            ref_chksum = mpp_chksum(compressed_c_cc_z_real_3D_field_data_ref(:,:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                     ",mpp_pe()
                write(output_unit,*) "Restart iteration:                           ",q
                write(output_unit,*) "Max compressed c, cc, z real 3D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                    ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                  ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c, cc, z real 3D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c, cc, z real 3D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h, cc, z real 3D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_cc_z_real_3D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_cc_z_real_3D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(real_buffer_3D(field_dimension_sizes(1), &
                                    field_dimension_sizes(2), &
                                    field_dimension_sizes(3)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_cc_z_real_3D_field_name, &
                                          real_buffer_3D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do k = 1,size(compressed_h_cc_z_real_3D_field_data,3)
                do j = 1,size(compressed_h_cc_z_real_3D_field_data,2)
                    do i = 1,size(compressed_h_cc_z_real_3D_field_data,1)
                        compressed_h_cc_z_real_3D_field_data(i,j,k) = real_buffer_3D(i+offset,j,k)
                    enddo
                enddo
            enddo
            deallocate(real_buffer_3D)
            rmax_error = maxval(abs(compressed_h_cc_z_real_3D_field_data - &
                                    compressed_h_cc_z_real_3D_field_data_ref(:,:,:,q)))
            read_in_chksum = mpp_chksum(compressed_h_cc_z_real_3D_field_data)
            ref_chksum = mpp_chksum(compressed_h_cc_z_real_3D_field_data_ref(:,:,:,q))
            if (rmax_error .ne. 0.0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                     ",mpp_pe()
                write(output_unit,*) "Restart iteration:                           ",q
                write(output_unit,*) "Max compressed h, cc, z real 3D field error: ",rmax_error
                write(output_unit,'(a,z16)') "Read in checksum:                    ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                  ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h, cc, z real 3D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h, cc, z real 3D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Integer scalar field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    int_scalar_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(int_scalar_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(int_buffer_1D(1))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          int_scalar_field_name, &
                                          int_buffer_1D, &
                                          unstructured_domain, &
                                          timelevel=q)
            int_scalar_field_data = int_buffer_1D(1)
            deallocate(int_buffer_1D)
            imax_error = abs(int_scalar_field_data - &
                             int_scalar_field_data_ref(q))
            if (imax_error .ne. 0) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                        ",mpp_pe()
                write(output_unit,*) "Restart iteration:              ",q
                write(output_unit,*) "Max integer scalar field error: ",imax_error
                call mpp_error(FATAL, &
                               "test 1: integer scalar field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,*) "Integer scalar field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c integer 1D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_int_1D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_int_1D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(int_buffer_1D(field_dimension_sizes(1)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_int_1D_field_name, &
                                          int_buffer_1D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do i = 1,compressed_c_axis_size
                compressed_c_int_1D_field_data(i) = int_buffer_1D(i+offset)
            enddo
            deallocate(int_buffer_1D)
            imax_error = maxval(abs(compressed_c_int_1D_field_data - &
                                    compressed_c_int_1D_field_data_ref(:,q)))
            read_in_chksum = mpp_chksum(compressed_c_int_1D_field_data)
            ref_chksum = mpp_chksum(compressed_c_int_1D_field_data_ref(:,q))
            if (imax_error .ne. 0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                 ",mpp_pe()
                write(output_unit,*) "Restart iteration:                       ",q
                write(output_unit,*) "Max compressed c integer 1D field error: ",imax_error
                write(output_unit,'(a,z16)') "Read in checksum:                ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:              ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c integer 1D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c integer 1D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h integer 1D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_int_1D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_int_1D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(int_buffer_1D(field_dimension_sizes(1)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_int_1D_field_name, &
                                          int_buffer_1D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do i = 1,compressed_h_axis_size
                compressed_h_int_1D_field_data(i) = int_buffer_1D(i+offset)
            enddo
            deallocate(int_buffer_1D)
            imax_error = maxval(abs(compressed_h_int_1D_field_data - &
                                    compressed_h_int_1D_field_data_ref(:,q)))
            read_in_chksum = mpp_chksum(compressed_h_int_1D_field_data)
            ref_chksum = mpp_chksum(compressed_h_int_1D_field_data_ref(:,q))
            if (imax_error .ne. 0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                 ",mpp_pe()
                write(output_unit,*) "Restart iteration:                       ",q
                write(output_unit,*) "Max compressed h integer 1D field error: ",imax_error
                write(output_unit,'(a,z16)') "Read in checksum:                ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:              ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h integer 1D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h integer 1D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c, z integer 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_z_int_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_z_int_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(int_buffer_2D(field_dimension_sizes(1), &
                                   field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_z_int_2D_field_name, &
                                          int_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_c_z_int_2D_field_data,2)
                do i = 1,size(compressed_c_z_int_2D_field_data,1)
                    compressed_c_z_int_2D_field_data(i,j) = int_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(int_buffer_2D)
            imax_error = maxval(abs(compressed_c_z_int_2D_field_data - &
                                    compressed_c_z_int_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_c_z_int_2D_field_data)
            ref_chksum = mpp_chksum(compressed_c_z_int_2D_field_data_ref(:,:,q))
            if (imax_error .ne. 0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                    ",mpp_pe()
                write(output_unit,*) "Restart iteration:                          ",q
                write(output_unit,*) "Max compressed c, z integer 2D field error: ",imax_error
                write(output_unit,'(a,z16)') "Read in checksum:                   ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                 ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c, z integer 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c, z integer 2D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h, z integer 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_z_int_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_z_int_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(int_buffer_2D(field_dimension_sizes(1), &
                                   field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_z_int_2D_field_name, &
                                          int_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_h_z_int_2D_field_data,2)
                do i = 1,size(compressed_h_z_int_2D_field_data,1)
                    compressed_h_z_int_2D_field_data(i,j) = int_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(int_buffer_2D)
            imax_error = maxval(abs(compressed_h_z_int_2D_field_data - &
                                    compressed_h_z_int_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_h_z_int_2D_field_data)
            ref_chksum = mpp_chksum(compressed_h_z_int_2D_field_data_ref(:,:,q))
            if (imax_error .ne. 0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                    ",mpp_pe()
                write(output_unit,*) "Restart iteration:                          ",q
                write(output_unit,*) "Max compressed h, z integer 2D field error: ",imax_error
                write(output_unit,'(a,z16)') "Read in checksum:                   ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                 ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h, z integer 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h, z integer 2D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed c, cc integer 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_c_cc_int_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_c_cc_int_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(int_buffer_2D(field_dimension_sizes(1), &
                                   field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_c_cc_int_2D_field_name, &
                                          int_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_c_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_c_cc_int_2D_field_data,2)
                do i = 1,size(compressed_c_cc_int_2D_field_data,1)
                    compressed_c_cc_int_2D_field_data(i,j) = int_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(int_buffer_2D)
            imax_error = maxval(abs(compressed_c_cc_int_2D_field_data - &
                                    compressed_c_cc_int_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_c_cc_int_2D_field_data)
            ref_chksum = mpp_chksum(compressed_c_cc_int_2D_field_data_ref(:,:,q))
            if (imax_error .ne. 0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                     ",mpp_pe()
                write(output_unit,*) "Restart iteration:                           ",q
                write(output_unit,*) "Max compressed c, cc integer 2D field error: ",imax_error
                write(output_unit,'(a,z16)') "Read in checksum:                    ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                  ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed c, z integer 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed c, cc integer 2D field data correct."
                endif
            endif
            call mpp_sync()

           !-------------------------------------------------------------------
           !Compressed h, cc integer 2D field.
            call mpp_sync()
            call fms_io_unstructured_get_field_size("RESTART/"//trim(restart_file_name), &
                                                    compressed_h_cc_int_2D_field_name, &
                                                    field_dimension_sizes, &
                                                    unstructured_domain, &
                                                    field_found_in_file)
            if (.not. field_found_in_file) then
                call mpp_error(FATAL, &
                               "test 1: field "//trim(compressed_h_cc_int_2D_field_name) &
                               //" was not found in file " &
                               //trim(restart_file_name))
            endif
            allocate(int_buffer_2D(field_dimension_sizes(1), &
                                   field_dimension_sizes(2)))
            call fms_io_unstructured_read("RESTART/"//trim(restart_file_name), &
                                          compressed_h_cc_int_2D_field_name, &
                                          int_buffer_2D, &
                                          unstructured_domain, &
                                          timelevel=q)
            if (mpp_pe() .eq. pelist(1)) then
                offset = 0
            else
                offset = 0
                do i = 1,io_domain_npes
                    if (mpp_pe() .eq. pelist(i)) then
                        exit
                    else
                        offset = offset + compressed_h_axis_size_per_rank(i)
                    endif
                enddo
            endif
            do j = 1,size(compressed_h_cc_int_2D_field_data,2)
                do i = 1,size(compressed_h_cc_int_2D_field_data,1)
                    compressed_h_cc_int_2D_field_data(i,j) = int_buffer_2D(i+offset,j)
                enddo
            enddo
            deallocate(int_buffer_2D)
            imax_error = maxval(abs(compressed_h_cc_int_2D_field_data - &
                                    compressed_h_cc_int_2D_field_data_ref(:,:,q)))
            read_in_chksum = mpp_chksum(compressed_h_cc_int_2D_field_data)
            ref_chksum = mpp_chksum(compressed_h_cc_int_2D_field_data_ref(:,:,q))
            if (imax_error .ne. 0 .or. read_in_chksum .ne. ref_chksum) then
                write(output_unit,*)
                write(output_unit,*) "My rank:                                    ",mpp_pe()
                write(output_unit,*) "Restart iteration:                           ",q
                write(output_unit,*) "Max compressed h, cc integer 2D field error: ",imax_error
                write(output_unit,'(a,z16)') "Read in checksum:                    ",read_in_chksum
                write(output_unit,'(a,z16)') "Reference checksum:                  ",ref_chksum
                call mpp_error(FATAL, &
                               "test 1: compressed h, z integer 2D field data incorrect.")
            else
                if (mpp_pe() .eq. mpp_root_pe()) then
                    write(output_unit,*)
                    write(output_unit,'(a,z16)') "Read-in data check-sum:  ",read_in_chksum
                    write(output_unit,'(a,z16)') "Reference data check-sum:",ref_chksum
                    write(output_unit,*) "Compressed h, cc integer 2D field data correct."
                endif
            endif
            call mpp_sync()
        enddo

       !Close the file.
        call fms_io_unstructured_file_unit("RESTART/"//trim(restart_file_name), &
                                           funit, &
                                           unstructured_domain)
        call mpp_close(funit)

       !Deallocate local allocatables.
        if (allocated(pelist)) then
            deallocate(pelist)
        endif
        if (allocated(compressed_c_axis_size_per_rank)) then
            deallocate(compressed_c_axis_size_per_rank)
        endif
        if (allocated(compressed_h_axis_size_per_rank)) then
            deallocate(compressed_h_axis_size_per_rank)
        endif
        if (allocated(x_axis_data)) then
            deallocate(x_axis_data)
        endif
        if (allocated(y_axis_data)) then
            deallocate(y_axis_data)
        endif
        if (allocated(z_axis_data)) then
            deallocate(z_axis_data)
        endif
        if (allocated(cc_axis_data)) then
            deallocate(cc_axis_data)
        endif
        if (allocated(compressed_c_axis_data)) then
            deallocate(compressed_c_axis_data)
        endif
        if (allocated(compressed_h_axis_data)) then
            deallocate(compressed_h_axis_data)
        endif
        if (allocated(compressed_c_real_1D_field_data)) then
            deallocate(compressed_c_real_1D_field_data)
        endif
        if (allocated(compressed_h_real_1D_field_data)) then
            deallocate(compressed_h_real_1D_field_data)
        endif
        if (allocated(compressed_c_z_real_2D_field_data)) then
            deallocate(compressed_c_z_real_2D_field_data)
        endif
        if (allocated(compressed_h_z_real_2D_field_data)) then
            deallocate(compressed_h_z_real_2D_field_data)
        endif
        if (allocated(compressed_c_cc_real_2D_field_data)) then
            deallocate(compressed_c_cc_real_2D_field_data)
        endif
        if (allocated(compressed_h_cc_real_2D_field_data)) then
            deallocate(compressed_h_cc_real_2D_field_data)
        endif
        if (allocated(compressed_c_z_cc_real_3D_field_data)) then
            deallocate(compressed_c_z_cc_real_3D_field_data)
        endif
        if (allocated(compressed_h_z_cc_real_3D_field_data)) then
            deallocate(compressed_h_z_cc_real_3D_field_data)
        endif
        if (allocated(compressed_c_cc_z_real_3D_field_data)) then
            deallocate(compressed_c_cc_z_real_3D_field_data)
        endif
        if (allocated(compressed_h_cc_z_real_3D_field_data)) then
            deallocate(compressed_h_cc_z_real_3D_field_data)
        endif
        if (allocated(compressed_c_int_1D_field_data)) then
            deallocate(compressed_c_int_1D_field_data)
        endif
        if (allocated(compressed_h_int_1D_field_data)) then
            deallocate(compressed_h_int_1D_field_data)
        endif
        if (allocated(compressed_c_z_int_2D_field_data)) then
            deallocate(compressed_c_z_int_2D_field_data)
        endif
        if (allocated(compressed_h_z_int_2D_field_data)) then
            deallocate(compressed_h_z_int_2D_field_data)
        endif
        if (allocated(compressed_c_cc_int_2D_field_data)) then
            deallocate(compressed_c_cc_int_2D_field_data)
        endif
        if (allocated(compressed_h_cc_int_2D_field_data)) then
            deallocate(compressed_h_cc_int_2D_field_data)
        endif
        if (allocated(real_scalar_field_data_ref)) then
            deallocate(real_scalar_field_data_ref)
        endif
        if (allocated(compressed_c_real_1D_field_data_ref)) then
            deallocate(compressed_c_real_1D_field_data_ref)
        endif
        if (allocated(compressed_h_real_1D_field_data_ref)) then
            deallocate(compressed_h_real_1D_field_data_ref)
        endif
        if (allocated(compressed_c_z_real_2D_field_data_ref)) then
            deallocate(compressed_c_z_real_2D_field_data_ref)
        endif
        if (allocated(compressed_h_z_real_2D_field_data_ref)) then
            deallocate(compressed_h_z_real_2D_field_data_ref)
        endif
        if (allocated(compressed_c_cc_real_2D_field_data_ref)) then
            deallocate(compressed_c_cc_real_2D_field_data_ref)
        endif
        if (allocated(compressed_h_cc_real_2D_field_data_ref)) then
            deallocate(compressed_h_cc_real_2D_field_data_ref)
        endif
        if (allocated(compressed_c_z_cc_real_3D_field_data_ref)) then
            deallocate(compressed_c_z_cc_real_3D_field_data_ref)
        endif
        if (allocated(compressed_h_z_cc_real_3D_field_data_ref)) then
            deallocate(compressed_h_z_cc_real_3D_field_data_ref)
        endif
        if (allocated(compressed_c_cc_z_real_3D_field_data_ref)) then
            deallocate(compressed_c_cc_z_real_3D_field_data_ref)
        endif
        if (allocated(compressed_h_cc_z_real_3D_field_data_ref)) then
            deallocate(compressed_h_cc_z_real_3D_field_data_ref)
        endif
        if (allocated(int_scalar_field_data_ref)) then
            deallocate(int_scalar_field_data_ref)
        endif
        if (allocated(compressed_c_int_1D_field_data_ref)) then
            deallocate(compressed_c_int_1D_field_data_ref)
        endif
        if (allocated(compressed_h_int_1D_field_data_ref)) then
            deallocate(compressed_h_int_1D_field_data_ref)
        endif
        if (allocated(compressed_c_z_int_2D_field_data_ref)) then
            deallocate(compressed_c_z_int_2D_field_data_ref)
        endif
        if (allocated(compressed_h_z_int_2D_field_data_ref)) then
            deallocate(compressed_h_z_int_2D_field_data_ref)
        endif
        if (allocated(compressed_c_cc_int_2D_field_data_ref)) then
            deallocate(compressed_c_cc_int_2D_field_data_ref)
        endif
        if (allocated(compressed_h_cc_int_2D_field_data_ref)) then
            deallocate(compressed_h_cc_int_2D_field_data_ref)
        endif

       !Print out a message that the test is done.
        call mpp_sync()
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*) "Test 1 body complete."
            write(output_unit,*)
        endif

        return
    end subroutine test_1

   !--------------------------------------------------------------------------

end program test_unstructured_fms_io
#endif
