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

! This program runs only one of many possible tests with each execution.
! Each test ends with an intentional fatal error.
! diag_manager_mod is not a stateless module, and there are situations
! where a fatal error leaves the module in a state that does not allow
! it to function properly if used again. Therefore, the program must
! be terminated after each intentional fatal error.

! Each test is dependent on the diag_table, and different diag_tables
! exist for each test. Depending on the test, an intentional fatal error
! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
! Because of this, the calls to all of those routines differ depending on the test.

! The diag_table for each test is included below.

!--------------------------------------------------------------------------------------------------
! diag_table for test 1

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 2

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 3

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 4

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 5

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 6

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 7

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat1", "dat1", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 8

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 9

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 10

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "bk",   "bk",   "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 11

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 12

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
! # Test of the error check that duplicate field names do not appear in same file,
!  "test_mod",              "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 13

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test",  1, "days", 1, "days", "time",
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test",  "all", .false., "none", 2,
! # Test of WARNING message that no data is written when run length is less than output interval
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 14

! test_diag_manager
! 1990 1 29 0 0 0
! #output files
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
! # Test of check for invalid date. (Jan 29 1990 + one month = Feb 29 1990)
!  "test_mod",              "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 16

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test2", 1, "months", 1, "days", "time",
! #output variables
! # Test for output file name to be modified with appended string
!  "test_diag_manager_mod", "dat2", "dat2", "diag_test2", "all", .false., "none", 2,
!--------------------------------------------------------------------------------------------------
! diag_table for test 17

! test_diag_manager
! 1 3 1 0 0 0
! #output files
!  "diag_test2", 1, "days", 1, "days", "time",
! #output variables
!  "test_diag_manager_mod", "dat2", "dat2_rms", "diag_test2", "all", "rms",  "none", 2,
!  "test_diag_manager_mod", "dat2", "dat2",     "diag_test2", "all", .true., "none", 2,
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
!> diag_table for test 23 (unstructured grid)
!!
!test_diag_manager_23
!1990 1 1 0 0 0
!#output files
!"unstructured_diag_test", 2, "days", 2, "days", "time",
!#output variables
!"UG_unit_test", "unstructured_real_scalar_field_data", "rsf_diag_1", "unstructured_diag_test", "all", .TRUE., "none", 1,
!"UG_unit_test", "unstructured_real_1D_field_data", "unstructured_real_1D_field_data", "unstructured_diag_test", "all", .TRUE., "none", 1,
!"UG_unit_test", "unstructured_real_2D_field_data", "unstructured_real_2D_field_data", "unstructured_diag_test", "all", .TRUE., "none", 1,
!"UG_unit_test", "lon", "grid_xt", "unstructured_diag_test", "all", .TRUE., "none", 1,
!"UG_unit_test", "lat", "grid_yt", "unstructured_diag_test", "all", .TRUE., "none", 1,
!--------------------------------------------------------------------------------------------------
PROGRAM test
  ! This program runs only one of many possible tests with each execution.
  ! Each test ends with an intentional fatal error.
  ! diag_manager_mod is not a stateless module, and there are situations
  ! where a fatal error leaves the module in a state that does not allow
  ! it to function properly if used again. Therefore, the program must
  ! be terminated after each intentional fatal error.

  ! Each test is dependent on the diag_table, and different diag_tables
  ! exist for each test. Depending on the test, an intentional fatal error
  ! may be triggered upon the call to diag_manager_init, register_diag_field or send_data.
  ! Because of this, the calls to all of those routines differ depending on the test.

  USE mpp_mod, ONLY: mpp_pe, mpp_root_pe, mpp_debug, mpp_set_stack_size
  USE mpp_io_mod, ONLY: mpp_io_init
  USE mpp_domains_mod, ONLY: domain2d, mpp_define_domains, mpp_get_compute_domain
  USE mpp_domains_mod, ONLY: mpp_define_io_domain, mpp_define_layout
  USE mpp_domains_mod, ONLY: mpp_domains_init, mpp_domains_set_stack_size
  USE fms_mod, ONLY: fms_init, fms_end, mpp_npes, file_exist, check_nml_error, open_file
  USE fms_mod, ONLY: error_mesg, FATAL, WARNING, NOTE, stdlog, stdout
#ifdef INTERNAL_FILE_NML
  USE mpp_mod, ONLY: input_nml_file
#else
  USE fms_mod, ONLY:  open_namelist_file, close_file
#endif
  USE fms_io_mod, ONLY: fms_io_init
  USE fms_io_mod, ONLY: fms_io_exit, set_filename_appendix
  USE constants_mod, ONLY: constants_init, PI, RAD_TO_DEG

  USE time_manager_mod, ONLY: time_type, set_calendar_type, set_date, decrement_date, OPERATOR(+), set_time
  USE time_manager_mod, ONLY: NOLEAP, JULIAN, GREGORIAN, THIRTY_DAY_MONTHS, OPERATOR(*), assignment(=)
  use time_manager_mod, ONLY: OPERATOR(+), OPERATOR(-), OPERATOR(/), days_in_month

  USE diag_manager_mod, ONLY: diag_manager_init, send_data, diag_axis_init, diag_manager_end
  USE diag_manager_mod, ONLY: register_static_field, register_diag_field, diag_send_complete
  USE diag_manager_mod, ONLY: diag_manager_set_time_end, diag_field_add_attribute, diag_axis_add_attribute
  USE diag_manager_mod, ONLY: diag_field_add_cell_measures
  USE diag_manager_mod, ONLY: get_diag_field_id, DIAG_FIELD_NOT_FOUND
  USE diag_axis_mod, ONLY: get_axis_num
#include "fms_platform.h"
  IMPLICIT NONE

  TYPE(domain2d) :: Domain1
  TYPE(domain2d) :: Domain2

  REAL, ALLOCATABLE, DIMENSION(:) :: lon_global1, lonb_global1
  REAL, ALLOCATABLE, DIMENSION(:) :: lat_global1, latb_global1
  REAL, ALLOCATABLE, DIMENSION(:) :: lon_global2, lonb_global2
  REAL, ALLOCATABLE, DIMENSION(:) :: lat_global2, latb_global2
  REAL, ALLOCATABLE, DIMENSION(:) :: pfull, bk, phalf
  REAL, ALLOCATABLE, DIMENSION(:) :: lon1, lat1, lonb1, latb1
  REAL, ALLOCATABLE, DIMENSION(:) :: lon2, lat2, lonb2, latb2
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dat1, dat1h
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: dat2, dat2h
  REAL, ALLOCATABLE, DIMENSION(:,:) :: dat2_2d
  REAL :: solar_constant = 1600
  REAL :: surf_press = 1.e5
  REAL :: dp
  INTEGER :: id_phalf, id_pfull, id_bk
  INTEGER :: id_lon1, id_lonb1, id_latb1, id_lat1, id_dat1
  INTEGER :: id_lon2, id_lat2, id_dat2, id_dat2_2d, id_sol_con, id_dat2h, id_dat2h_2
  INTEGER :: id_dat2_got, id_none_got
  INTEGER :: i, j, k, is1, ie1, js1, je1, nml_unit, ierr, log_unit, out_unit, m
  INTEGER :: is_in, ie_in, js_in, je_in
  INTEGER :: is2, ie2, js2, je2, hi=1, hj=1
  INTEGER :: nlon1, nlat1, nlon2, nlat2
  INTEGER, DIMENSION(2) :: layout = (/1,1/)
  INTEGER :: test_number=1
  INTEGER :: nlon=10, nlat=10, nlev=10
  INTEGER :: io_layout(2) = (/1,1/)
  INTEGER :: nstep = 2
  TYPE(time_type) :: Time, Time_step, Time_end, Time_start, Run_length
  LOGICAL :: used, test_successful
  CHARACTER(len=256) :: err_msg
  INTEGER :: omp_get_num_threads

  INTEGER :: nyc1, n, jsw, jew, isw, iew
  INTEGER :: numthreads=1, ny_per_thread, idthread
  INTEGER :: months=0, days=1, dt_step=1

  ! Variables needed for test 22
  INTEGER :: id_nv, id_nv_init

!!!!!! Stuff for unstrctured grid
    integer(INT_KIND)              :: nx = 8                               !<Total number of grid points in the x-dimension (longitude?)
    integer(INT_KIND)              :: ny = 8                               !<Total number of grid points in the y-dimension (latitude?)
    integer(INT_KIND)              :: nz = 2                               !<Total number of grid points in the z-dimension (height)
    integer(INT_KIND)              :: nt = 2                               !<Total number of time grid points.
    integer(INT_KIND)              :: io_tile_factor = 1                   !< The IO tile factor
    integer(INT_KIND)              :: halo = 2                             !<Number of grid points in the halo???
    integer(INT_KIND)              :: ntiles_x = 1                         !<Number of tiles in the x-direction (A 2D grid of tiles is used in this test.)
    integer(INT_KIND)              :: ntiles_y = 2                         !<Number of tiles in the y-direction (A 2D grid of tiles is used in this test.)
    integer(INT_KIND)              :: total_num_tiles                      !<The total number of tiles for the run (= ntiles_x*ntiles_y)
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
    type(time_type)                :: diag_time                            !<

    integer(INT_KIND)              :: output_unit=6
!!!!!!



  NAMELIST /test_diag_manager_nml/ layout, test_number, nlon, nlat, nlev, io_layout, numthreads, &
                                   dt_step, months, days
  NAMELIST /utest_nml/nx,ny,nz,nt,ntiles_x,ntiles_y,io_tile_factor
  ! Initialize all id* vars to be -1
  id_nv = -1
  id_nv_init = -1
  id_phalf = -1
  id_pfull = -1
  id_bk = -1
  id_lon1 = -1
  id_lonb1 = -1
  id_latb1 = -1
  id_lat1 = -1
  id_dat1 = -1
  id_lon2 = -1
  id_lat2 = -1
  id_dat2 = -1
  id_dat2_2d = -1
  id_sol_con = -1
  id_dat2h = -1
  id_dat2h_2 = -1
  id_dat2_got = -1
  id_none_got = -1

  CALL fms_init
  log_unit = stdlog()
  out_unit = stdout()
  CALL constants_init
  CALL set_calendar_type(JULIAN)
  npes = mpp_npes()
#ifdef INTERNAL_FILE_NML
  READ (input_nml_file, NML=test_diag_manager_nml, IOSTAT=ierr)
  READ (input_nml_file, NML=utest_nml, IOSTAT=i)
#else
  IF ( file_exist('input.nml') ) THEN
     nml_unit = open_namelist_file()
     READ(nml_unit, nml=test_diag_manager_nml, iostat=ierr)
     READ(nml_unit, nml=utest_nml, iostat=i)
     CALL close_file(nml_unit)
  ELSE
     ! Set ierr to an arbitrary positive number if input.nml does not exist.
     ierr = 100
  END IF
#endif
  ! Check the status of reading the diag_manager_nml
  IF ( check_nml_error(IOSTAT=ierr, NML_NAME='DIAG_MANAGER_NML') < 0 ) THEN
     IF ( mpp_pe() == mpp_root_pe() ) THEN
        CALL error_mesg('diag_manager_mod::diag_manager_init', 'TEST_DIAG_MANAGER_NML not found in input.nml.  Using defaults.',&
             & WARNING)
     END IF
  END IF
  WRITE (log_unit,test_diag_manager_nml)

SELECT CASE ( test_number ) ! Closes just before the CONTAINS block.
  ! If the test_number == 23, then call the unstrcutured grid unit test and skip everything else.
  CASE ( 23 )
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
        write(output_unit,*) "Using Unstructured domaintypes and calls..."
    endif

   !Add a suffix to the test file.
    write(test_file,'(a,i3.3)') trim(test_file),npes

   !Initialize the diag manager module.
    call diag_manager_init()

   !Set the diag_time variable to be 01/01/1990 at 00:00:00 (midnight).
    call set_calendar_type(JULIAN)
    time = set_date(1990,1,1,0,0,0)
   CALL unstruct_test (nx, ny, nz, npes, ntiles_x, 1, time,io_tile_factor)

   ! If the test_number == 12, check for the correct error and skip everything else.
   CASE ( 12 ) 
     CALL diag_manager_init(err_msg=err_msg)
     IF ( err_msg /= '' ) THEN
        WRITE (out_unit,'(a)') 'test12 successful: err_msg='//TRIM(err_msg)
        CALL error_mesg('test_diag_manager','test12 successful.',NOTE)
     ELSE
        WRITE (out_unit,'(a)') 'test12 fails'
        CALL error_mesg('test_diag_manager','test12 fails',FATAL)
     END IF

  ! If the test number is not 12 or 23, run all other tests.
  CASE DEFAULT ! Contains all remaining code up to CONTAINS block.
     CALL diag_manager_init

  IF ( layout(1)*layout(2) .NE. mpp_npes() ) THEN
     CALL mpp_define_layout((/1,nlon,1,nlat/), mpp_npes(), layout )
  END IF

  nlon1 = nlon
  nlat1 = nlat
  nlon2 = nlon * 2
  nlat2 = nlat * 2

  CALL mpp_define_domains((/1,nlon1,1,nlat1/), layout, Domain1, name='test_diag_manager')
  CALL mpp_get_compute_domain(Domain1, is1, ie1, js1, je1)
  ALLOCATE(lon_global1(nlon1), lonb_global1(nlon1+1))
  ALLOCATE(lat_global1(nlat1), latb_global1(nlat1+1))
  ALLOCATE(lon_global2(nlon2), lonb_global2(nlon2+1))
  ALLOCATE(lat_global2(nlat2), latb_global2(nlat2+1))
  ALLOCATE(pfull(nlev), bk(nlev), phalf(nlev+1))

  ALLOCATE(lon1(is1:ie1), lat1(js1:je1), lonb1(is1:ie1+1), latb1(js1:je1+1))
  CALL compute_grid(nlon1, nlat1, is1, ie1, js1, je1, lon_global1, lat_global1, lonb_global1, latb_global1, lon1, lat1, lonb1, latb1)
  CALL mpp_define_domains((/1,nlon2,1,nlat2/), layout, Domain2, name='test_diag_manager')
  CALL mpp_get_compute_domain(Domain2, is2, ie2, js2, je2)
  CALL mpp_define_io_domain(Domain1, io_layout)
  CALL mpp_define_io_domain(Domain2, io_layout)

  ALLOCATE(lon2(is2:ie2), lat2(js2:je2), lonb2(is2:ie2+1), latb2(js2:je2+1))
  CALL compute_grid(nlon2, nlat2, is2, ie2, js2, je2, lon_global2, lat_global2, lonb_global2, latb_global2, lon2, lat2, lonb2, latb2)
  dp = surf_press/nlev
  DO k=1, nlev+1
     phalf(k) = dp*(k-1)
  END DO
  DO k=1, nlev
     pfull(k) = .5*(phalf(k) + phalf(k+1))
     bk(k) = pfull(k)/surf_press
  END DO

  ALLOCATE(dat1(is1:ie1,js1:je1,nlev))
  ALLOCATE(dat1h(is1-hi:ie1+hi,js1-hj:je1+hj,nlev))
  dat1h = 0.
  DO j=js1, je1
     DO i=is1, ie1
        dat1(i,j,1) = SIN(lon1(i))*COS(lat1(j))
     END DO
  END DO
  dat1h(is1:ie1,js1:je1,1) = dat1(:,:,1)
  dat1(:,:,2) = -dat1(:,:,1)
  dat1h(:,:,2) = -dat1h(:,:,1)

  ALLOCATE(dat2(is2:ie2,js2:je2,nlev))
  ALLOCATE(dat2_2d(is2:ie2,js2:je2))
  ALLOCATE(dat2h(is2-hi:ie2+hi,js2-hj:je2+hj,nlev))
  dat2h = 0.
  dat2 = 0.
  DO j=js2, je2
     DO i=is2, ie2
        dat2(i,j,1) = SIN(lon2(i))*COS(lat2(j))
     END DO
  END DO
  dat2h(is2:ie2,js2:je2,1) = dat2(:,:,1)
  dat2(:,:,2) = -dat2(:,:,1)
  dat2h(:,:,2) = -dat2h(:,:,1)
  dat2_2d = dat2(:,:,1)

  id_lonb1 = diag_axis_init('lonb1', RAD_TO_DEG*lonb_global1, 'degrees_E', 'x', long_name='longitude edges', Domain2=Domain1)
  id_latb1 = diag_axis_init('latb1', RAD_TO_DEG*latb_global1, 'degrees_N', 'y', long_name='latitude edges',  Domain2=Domain1)

  id_lon1  = diag_axis_init('lon1',  RAD_TO_DEG*lon_global1, 'degrees_E','x',long_name='longitude',Domain2=Domain1,edges=id_lonb1)
  id_lat1  = diag_axis_init('lat1',  RAD_TO_DEG*lat_global1, 'degrees_N','y',long_name='latitude', Domain2=Domain1,edges=id_latb1)

  id_phalf= diag_axis_init('phalf', phalf, 'Pa', 'z', long_name='half pressure level', direction=-1)
  id_pfull= diag_axis_init('pfull', pfull, 'Pa', 'z', long_name='full pressure level', direction=-1, edges=id_phalf)

  id_lon2 = diag_axis_init('lon2',  RAD_TO_DEG*lon_global2,  'degrees_E', 'x', long_name='longitude', Domain2=Domain2)
  id_lat2 = diag_axis_init('lat2',  RAD_TO_DEG*lat_global2,  'degrees_N', 'y', long_name='latitude',  Domain2=Domain2)

  IF ( test_number == 22 ) THEN
     ! Can we get the 'nv' axis ID?
     id_nv = get_axis_num('nv', 'nv')
     IF ( id_nv .GT. 0 ) THEN
        write (out_unit,'(a)') 'test22.1 Passes: id_nv has a positive value'
     ELSE
        write (out_unit,'(a)') 'test22.1 Failed: id_nv does not have a positive value'
     END IF

     ! Can I call diag_axis_init on 'nv' again, and get the same ID back?
     id_nv_init = diag_axis_init( 'nv',(/1.,2./),'none','N','vertex number', set_name='nv')
     IF ( id_nv_init .EQ. id_nv ) THEN
        write (out_unit,'(a)') 'test22.2 Passes: Can call diag_axis_init on "nv" and get same ID'
     ELSE
        write (out_unit,'(a)') 'test22.2 Failed: Cannot call diag_axis_init on "nv" and get same ID'
     END IF
  END IF

  IF ( test_number == 21 ) THEN
     ! Testing addition of axis attributes
     CALL diag_axis_add_attribute(id_lon1, 'real_att', 2.3)
     CALL diag_axis_add_attribute(id_lat1, 'int_att', (/ 2, 3 /))
     CALL diag_axis_add_attribute(id_pfull, 'char_att', 'Some string')
  END IF

  IF ( test_number == 14 ) THEN
     Time = set_date(1990,1,29,0,0,0)
  ELSE
     Time = set_date(1990,1,1,0,0,0)
  END IF

  IF ( test_number == 16 ) THEN
     ! Test 16 tests the filename appendix
     CALL set_filename_appendix('g01')
  END IF
  id_dat1 = register_diag_field('test_diag_manager_mod', 'dat1', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K')
  IF ( test_number == 18 ) THEN
     CALL diag_field_add_attribute(id_dat1, 'real_att', 2.3)
     CALL diag_field_add_attribute(id_dat1, 'cell_methods', 'area: mean')
     CALL diag_field_add_attribute(id_dat1, 'cell_methods', 'lon: mean')
  END IF
  IF ( test_number == 18 .OR. test_number == 19 ) THEN
     id_dat2 = register_diag_field('test_diag_manager_mod', 'dat2', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K')
     CALL diag_field_add_attribute(id_dat2, 'interp_method', 'none')
     CALL diag_field_add_attribute(id_dat2, 'int_att', (/ 1, 2 /) )
  ELSE
     id_dat2 = register_diag_field('test_diag_manager_mod', 'dat2', (/id_lon2,id_lat2,id_pfull/), Time, 'sample data', 'K')
  END IF
  id_sol_con = register_diag_field ('test_diag_manager_mod', 'solar_constant', Time, &
                  'solar constant', 'watts/m2')

  IF ( test_number == 20 ) THEN
     id_dat2_got = get_diag_field_id('test_diag_manager_mod', 'dat2')
     IF ( id_dat2_got == id_dat2 ) THEN
        WRITE (out_unit,'(a)') 'test20.1 Passes, id_dat2.EQ.id_dat2_got'
     ELSE
        WRITE (out_unit,'(a)') 'test20.1 Failed, id_dat2.NE.id_dat2_got'
     END IF

     id_none_got = get_diag_field_id('no_mod', 'no_var')
     IF ( id_none_got == DIAG_FIELD_NOT_FOUND ) THEN
        write (out_unit,'(a)') 'test20.2 Passes, id_none_got.EQ.DIAG_FIELD_NOT_FOUND'
     ELSE
        write (out_unit,'(a)') 'test20.2 Failed, id_none_got.NE.DIAG_FIELD_NOT_FOUND'
     END IF
  END IF

  IF ( dt_step == 0 ) CALL error_mesg ('test_diag_manager',&
       & 'dt_step is not set', FATAL)

  Time_step = set_time(dt_step,0)
  Time_start = Time
  Time_end = Time
  DO m = 1,months
     Time_end = Time_end + set_time(0,days_in_month(Time_end))
  END DO
  Time_end   = Time_end + set_time(0, days)
  Run_length = Time_end - Time_start
  nstep = Run_length / Time_step

  IF ( test_number == 18 ) THEN
     id_dat2h = register_diag_field('test_mod', 'dat2h', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K',&
          & volume=id_dat1, area=id_dat2, realm='myRealm', err_msg=err_msg)
     IF ( err_msg /= '' .OR. id_dat2h <= 0 ) THEN
        CALL error_mesg ('test_diag_manager',&
             & 'Unexpected error registering dat2h '//err_msg, FATAL)
     END IF
     id_dat2h_2 = register_diag_field('test_mod', 'dat2h_2', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K',&
          & err_msg=err_msg)
     CALL diag_field_add_cell_measures(id_dat2h_2, area=id_dat2, volume=id_dat1)
  ELSE IF ( test_number == 19 ) THEN
     id_dat2h = register_diag_field('test_mod', 'dat2h', (/id_lon1,id_lat1,id_pfull/), Time, 'sample data', 'K',&
          & volume=id_dat1, area=id_dat1, err_msg=err_msg)
     IF ( err_msg /= '' .OR. id_dat2h <= 0 ) THEN
        CALL error_mesg ('test_diag_manager',&
             & 'Expected error registering dat2h '//err_msg, NOTE)
     END IF
  END IF

  IF ( test_number == 16 .OR. test_number == 17 .OR. test_number == 18 .OR. test_number == 21 .OR. test_number == 22 ) THEN
     is_in = 1
     js_in = 1
     ie_in = nlon
     je_in = nlat

     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     IF ( id_dat2h > 0 ) used = send_data(id_dat2h, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
     IF ( id_dat2h_2 > 0 ) used = send_data(id_dat2h_2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     IF ( id_dat2h > 0 ) used = send_data(id_dat2h, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
     IF ( id_dat2h_2 > 0 ) used = send_data(id_dat2h_2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
  END IF

  !-- The following is used to test openMP
  IF ( test_number == 15 ) THEN
!$      call omp_set_num_threads(numthreads)
      nyc1 = je1 - js1 + 1
      IF (MOD(nyc1, numthreads ) /= 0) THEN
         CALL error_mesg ('test_diag_manager',&
              & 'The number of OpenMP threads must be an integral multiple &
              &of the number of rows in the compute domain', FATAL)
     END IF
     ny_per_thread = nyc1/numthreads

     dat1 = 1
     CALL diag_manager_set_time_end(Time_end)
     DO n = 1, nstep

        Time = Time + Time_step
        !$OMP parallel do default(shared) private(isw, iew, jsw, jew )

        DO jsw = js1, je1, ny_per_thread
           jew = jsw + ny_per_thread -1
           isw = is1
           iew = ie1
           if(id_dat1>0) used = send_data(id_dat1, dat1(isw:iew, jsw:jew,:), Time, &
                                is_in=isw-is1+1, js_in=jsw-js1+1,err_msg=err_msg)
           if(id_sol_con>0) used = send_data(id_sol_con, solar_constant, Time, err_msg=err_msg )
        END DO
        !$OMP END parallel do
        !CALL diag_send_complete(Time_step)
     END DO
 END IF


  IF ( test_number == 14 ) THEN
     id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data', 'K', err_msg=err_msg)
     IF ( err_msg /= '' ) THEN
        WRITE (out_unit,'(a)') 'test14 successful. err_msg='//TRIM(err_msg)
     ELSE
        WRITE (out_unit,'(a)') 'test14 fails.'
     END IF
  ELSE
     id_dat2_2d = register_diag_field('test_mod', 'dat2', (/id_lon2,id_lat2/), Time, 'sample data', 'K')
  END IF

  id_bk = register_static_field('test_diag_manager_mod', 'bk', (/id_pfull/), 'half level sigma', 'none')

  IF ( test_number == 13 ) THEN
     IF ( id_dat2_2d > 0 ) used=send_data(id_dat2_2d, dat2(:,:,1), Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test13: successful if a WARNING message appears that refers to output interval greater than runlength'
     ELSE
        WRITE (out_unit,'(a)') 'test13 fails: err_msg='//TRIM(err_msg)
     END IF
  END IF

  ! Note: test12 involves diag_manager_init, it does not require a call to send_data.
  !       See call to diag_manager_init above.

  IF ( test_number == 11 ) THEN
     is_in = 1+hi
     js_in = 1+hj
     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj

     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, je_in=je_in, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test11.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test11.1 fails. err_msg='//TRIM(err_msg)
     END IF

     ! intentional_error: je_in is missing
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in, ie_in=ie_in, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test11.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test11.2 successful. err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 10 ) THEN
     !  1 window, no halos, static, 1 dimension, global data.

     IF ( id_bk > 0 ) used = send_data(id_bk, bk, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test10.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test10.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too large.
     IF ( id_bk > 0 ) used = send_data(id_bk, phalf, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE(out_unit,'(a)') 'test10.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test10.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 9 ) THEN
     !  1 window, no halos, static, 1 dimension, global data
     IF ( id_bk > 0 ) used = send_data(id_bk, bk, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test9.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test9.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small
     IF ( id_bk > 0 ) used = send_data(id_bk, bk(1:nlev-1), err_msg=err_msg) ! intentional_error
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test9.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test9.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 8 ) THEN
     !  1 window with halos
     is_in = 1+hi
     js_in = 1+hj

     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat2h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test8.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test8.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small in both x and y directions
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     ie_in = ie1-is1+1+hi
     je_in = je1-js1+1+hj
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used=send_data(id_dat2, dat1h, Time, is_in=is_in, js_in=js_in, &
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test8.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test8.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 7 ) THEN
     !  1 window with halos
     is_in = 1+hi
     js_in = 1+hj

     ie_in = ie1-is1+1+hi
     je_in = je1-js1+1+hj
     IF ( id_dat1 > 0 ) used=send_data(id_dat1, dat1h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test7.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test7.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too large in both x and y directions
     ie_in = ie2-is2+1+hi
     je_in = je2-js2+1+hj
     IF ( id_dat1 > 0 ) used=send_data(id_dat1, dat2h, Time, is_in=is_in, js_in=js_in,&
          & ks_in=1, ie_in=ie_in, je_in=je_in, ke_in=nlev, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test7.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test7.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 6 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within either do loop for test6.1
     test_successful = .TRUE.
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test6.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,:,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test6.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test6.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test6.1 fails.'
     END IF

     !  intentional_error: data array too small in y direction
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1)
     END DO
     Time = Time + set_time(0,1)
     DO i=is2, ie2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(i:i,js2:je2-1,:), Time, i-is2+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test6.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test6.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 5 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within either do loop for test5.1
     test_successful = .TRUE.
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test5.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(:,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test5.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test5.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test5.1 fails.'
     END IF

     !  intentional_error: data array too small in x direction.
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1)
     END DO
     Time = Time + set_time(0,1)
     DO j=js2, je2
        IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2(is2:ie2-1,j:j,:), Time, 1, j-js2+1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test5.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test5.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 4 ) THEN
     !  1 window, no halos
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat2, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test4.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test4.1 fails: err_msg='//TRIM(err_msg)
     END IF

     !  intentional_error: data array too small in both x and y directions
     !  Error check is done on second call to send_data. Change in value of Time triggers the check.
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     Time = Time + set_time(0,1)
     IF ( id_dat2 > 0 ) used = send_data(id_dat2, dat1, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test4.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test4.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 3 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within do loop for test3.1
     test_successful = .TRUE.
     DO i=is1, ie1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test3.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test3.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test3.1 fails.'
     END IF

     !  intentional_error: data array too large in y direction
     DO i=is1, ie1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2(i:i,:,:), Time, i-is1+1, 1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test3.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test3.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 2 ) THEN
     !  multiple windows, no halos
     !  No error messages should appear at any point within do loop for test2.1
     test_successful = .TRUE.
     DO j=js1, je1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
        IF ( err_msg /= '' ) THEN
           WRITE (out_unit,'(a)') 'test2.1 fails: err_msg='//TRIM(err_msg)
           test_successful = .FALSE.
        END IF
     END DO
     IF ( test_successful ) THEN
        WRITE (out_unit,'(a)') 'test2.1 successful.'
     ELSE
        WRITE (out_unit,'(a)') 'test2.1 fails.'
     END IF

     !  intentional_error: data array too large in x direction
     DO j=js1, je1
        IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2(:,j:j,:), Time, 1, j-js1+1, err_msg=err_msg)
        IF ( err_msg /= '' ) EXIT ! exit immediately after error is detected. No need to continue.
     END DO
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test2.2 fails.'
     ELSE
        WRITE (out_unit,'(a)') 'test2.2 successful: err_msg='//TRIM(err_msg)
     END IF
  END IF

  IF ( test_number == 1 ) THEN
     !  1 window, no halos
     ! Here dat2 is too large in both x and y direction so you should get an error.
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat2, Time, err_msg=err_msg)
     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test1.1 fails: Intentional error not detected'
     ELSE
        WRITE (out_unit,'(a)') 'test1.1 successful: '//TRIM(err_msg)
     END IF

     ! Here dat1 has the correct shape, so no error
     IF ( id_dat1 > 0 ) used = send_data(id_dat1, dat1, Time, err_msg=err_msg)

     IF ( err_msg == '' ) THEN
        WRITE (out_unit,'(a)') 'test1.2 successful'
     ELSE
        WRITE (out_unit,'(a)') 'test1.2 fails: '//TRIM(err_msg)
     END IF
  END IF
  CALL diag_manager_end(Time)
END SELECT ! End of case handling opened for test 12. 

  CALL fms_io_exit
  CALL fms_end

CONTAINS

  SUBROUTINE compute_grid(nlon, nlat, is, ie, js, je, lon_global, lat_global, lonb_global, latb_global, lon, lat, lonb, latb)
    INTEGER, INTENT(in) :: nlon, nlat, is, ie, js, je
    REAL, INTENT(out), DIMENSION(:) :: lon_global, lat_global, lonb_global, latb_global, lon, lat, lonb, latb

    REAL :: dlon, dlat
    INTEGER :: i, j

    dlon = 2*PI/nlon
    dlat = PI/nlat

    DO i=1, nlon+1
       lonb_global(i) = dlon*(i-1)
    END DO
    DO j=1,nlat+1
       latb_global(j) = dlat*(j-1) - .5*PI
    END DO
    DO i=1,nlon
       lon_global(i) = .5*(lonb_global(i) + lonb_global(i+1))
    END DO
    DO j=1,nlat
       lat_global(j) = .5*(latb_global(j) + latb_global(j+1))
    END DO
    lon  = lon_global(is:ie)
    lat  = lat_global(js:je)
    lonb = lonb_global(is:ie+1)
    latb = latb_global(js:je+1)
  END SUBROUTINE compute_grid

  SUBROUTINE unstruct_test(nx, ny, nz, npes, num_domain_tiles_x, num_domain_tiles_y, diag_time,io_tile_factor)
        use, intrinsic :: iso_fortran_env, only: output_unit
        use mpp_parameter_mod,             only: FATAL
        use mpp_mod,                       only: mpp_error, &
                                                 mpp_pe, &
                                                 mpp_root_pe, &
                                                 mpp_sync, &
                                                 mpp_chksum
        use mpp_domains_mod,               only: domain2D, &
                                                 mpp_define_mosaic, &
                                                 mpp_deallocate_domain, &
                                                 domainUG, &
                                                 mpp_define_unstruct_domain, &
!                                                mpp_deallocate_domainUG, &
                                                 mpp_get_UG_compute_domain, &
                                                 mpp_get_UG_domain_grid_index, &
                                                 mpp_get_UG_domain_ntiles
        use diag_axis_mod,                 only: diag_axis_init, diag_axis_add_attribute
        use diag_manager_mod,              only: register_diag_field, &
                                                 send_data
        use time_manager_mod,              only: time_type, &
                                                 set_time, &
                                                 operator(+), &
                                                 assignment(=)
        implicit none

       !Inputs/Ouputs
        integer(INT_KIND),intent(in)  :: nx                 !<The number of grid points in the x-direction.
        integer(INT_KIND),intent(in)  :: ny                 !<The number of grid points in the y-direction.
        integer(INT_KIND),intent(in)  :: nz                 !<The number of grid points in the z-direction.
        integer(INT_KIND),intent(in)  :: npes               !<The total number of ranks used in this test.
        integer(INT_KIND),intent(in)  :: num_domain_tiles_x !<The total number of domain tiles in the x-dimension for the 2D structured domain in this test.
        integer(INT_KIND),intent(in)  :: num_domain_tiles_y !<The total number of domain tiles in the y-dimension for the 2D structured domain in this test.
        type(time_type),intent(inout) :: diag_time          !<Time for diag_manager.
        integer(INT_KIND),intent(in)  :: io_tile_factor     !<I/O tile factor.  See below.

       !Local variables
        integer(INT_KIND)                              :: num_domain_tiles                           !<The total number of domain tiles for the 2D structured domain in this test.
        integer(INT_KIND)                              :: npes_per_domain_tile                       !<The number of ranks per domain tile for the 2D structured domain.
        integer(INT_KIND)                              :: my_domain_tile_id                          !<The 2D structured domain tile id for the current rank.
        logical(INT_KIND)                              :: is_domain_tile_root                        !<Flag telling if the current rank is the root rank of its associated 2D structured domain tile.
        integer(INT_KIND),dimension(2)                 :: layout_for_full_domain                     !<Rank layout (2D grid) for the full 2D structured domain. Example: 16 ranks -> (16,1) or (8,2) or (4,4) or (2,8) or (1,16)
        integer(INT_KIND),dimension(:),allocatable     :: pe_start                                   !<Array holding the smallest rank id assigned to each 2D structured domain tile.
        integer(INT_KIND),dimension(:),allocatable     :: pe_end                                     !<Array holding the largest rank id assigned to each 2D structured domain tile.
        integer(INT_KIND)                              :: x_grid_points_per_domain_tile              !<The number of grid points in the x-dimension on each 2D structured domain tile.
        integer(INT_KIND)                              :: y_grid_points_per_domain_tile              !<The number of grid points in the y-dimension on each 2D structured domain tile.
        integer(INT_KIND),dimension(:,:),allocatable   :: global_indices                             !<Required to define the 2D structured domain.
        integer(INT_KIND),dimension(:,:),allocatable   :: layout2D                                   !<Required to define the 2D structured domain.
        type(domain2D)                                 :: domain_2D                                  !<A structured 2D domain.
        logical(INT_KIND),dimension(:,:,:),allocatable :: land_mask                                  !<A toy mask.
        integer(INT_KIND),dimension(:),allocatable     :: num_non_masked_grid_points_per_domain_tile !<Total number of non-masked grid points on each 2D structured domain tile.
        integer(INT_KIND)                              :: mask_counter                               !<Counting variable.
        integer(INT_KIND)                              :: num_non_masked_grid_points                 !<Total number of non-masked grid points for the 2D structured domain.
        integer(INT_KIND),dimension(:),allocatable     :: num_land_tiles_per_non_masked_grid_point   !<Number of land tiles per non-masked grid point for the 2D structured domain.
        integer(INT_KIND)                              :: num_ranks_using_unstructured_grid          !<Number of ranks using the unstructured domain.
        integer(INT_KIND),dimension(:),allocatable     :: unstructured_grid_point_index_map          !<Array that maps indices between the 2D structured and unstructured domains.
        type(domainUG)                                 :: domain_ug                                  !<An unstructured mpp domain.
        integer(INT_KIND),dimension(:),allocatable     :: unstructured_axis_data                     !<Data that is registered to the restart file for the unstructured axis.
        integer(INT_KIND)                              :: unstructured_axis_data_size                !<Size of the unstructured axis data array.
        character(len=256)                             :: unstructured_axis_name                     !<Name for the unstructured axis.
        real,dimension(:),allocatable                  :: x_axis_data                                !<Data for the x-axis that is registered to the restart file.
        real,dimension(:),allocatable                  :: y_axis_data                                !<Data for the y-axis that is registered to the restart file.
        real,dimension(:),allocatable                  :: z_axis_data                                !<Data for the z-axis that is registered to the restart file.
        real                                           :: unstructured_real_scalar_field_data_ref    !<Reference test data for an unstructured real scalar field.
        real,dimension(:),allocatable                  :: unstructured_real_1D_field_data_ref        !<Reference test data for an unstructured real 1D field.
        real,dimension(:,:),allocatable                :: unstructured_real_2D_field_data_ref        !<Reference test data for an unstructured real 2D field.
        real,dimension(:,:,:),allocatable              :: unstructured_real_3D_field_data_ref        !<Reference test data for an unstructured real 3D field.
        integer                                        :: unstructured_int_scalar_field_data_ref     !<Reference test data for an unstructured integer scalar field.
        integer,dimension(:),allocatable               :: unstructured_int_1D_field_data_ref         !<Reference test data for an unstructured integer 1D field.
        integer,dimension(:,:),allocatable             :: unstructured_int_2D_field_data_ref         !<Reference test data for an unstructured integer 2D field.
        character(len=256)                             :: unstructured_real_scalar_field_name        !<Name for an unstructured real scalar field.
        real                                           :: unstructured_real_scalar_field_data        !<Data for an unstructured real scalar field.
        character(len=256)                             :: unstructured_real_1D_field_name            !<Name for an unstructured real 1D field.
        real,dimension(:),allocatable                  :: unstructured_real_1D_field_data            !<Data for an unstructured real 1D field.
        character(len=256)                             :: unstructured_real_2D_field_name            !<Name for an unstructured real 2D field.
        real,dimension(:,:),allocatable                :: unstructured_real_2D_field_data            !<Data for an unstructured real 2D field.
        character(len=256)                             :: unstructured_real_3D_field_name            !<Name for an unstructured real 3D field.
        real,dimension(:,:,:),allocatable              :: unstructured_real_3D_field_data            !<Data for an unstructured real 3D field.
        character(len=256)                             :: unstructured_int_scalar_field_name         !<Name for an unstructured integer scalar field.
        integer                                        :: unstructured_int_scalar_field_data         !<Data for an unstructured integer scalar field.
        character(len=256)                             :: unstructured_int_1D_field_name             !<Name for an unstructured integer 1D field.
        integer,dimension(:),allocatable               :: unstructured_int_1D_field_data             !<Data for an unstructured integer 1D field.
        character(len=256)                             :: unstructured_int_2D_field_name             !<Name for an unstructured integer 2D field.
        character(len=100)                             :: unstructured_1d_alt                       !<Name of the unstrucutred 1D field if L>1
        integer,dimension(:,:),allocatable             :: unstructured_int_2D_field_data             !<Data for an unstructured integer 2D field.
       integer(INT_KIND),allocatable,dimension(:)      :: unstructured_axis_diag_id                  !<Id returned for the unstructured axis by diag_axis_init.
       integer(INT_KIND)                               :: x_axis_diag_id                             !<Id returned for the x-axis by diag_axis_init.
       integer(INT_KIND)                              :: y_axis_diag_id                             !<Id returned for the y-axis by diag_axis_init.
       integer(INT_KIND)                              :: z_axis_diag_id                             !<Id returned for the z-axis by diag_axis_init.
       real,allocatable,dimension(:) :: lat, lon
       integer(INT_KIND)             :: idlat
       integer(INT_KIND)                              :: idlon
       integer(INT_KIND)                              :: rsf_diag_id                                !<Id returned for a real scalar field associated with the unstructured grid by
                                !!register_diag_field.
       integer(INT_KIND),allocatable,dimension(:)     :: rsf_diag_1d_id                             !<Id returned for a real 1D array  field associated with the unstructured grid by                                                                                                     !!register_diag_field.
       integer(INT_KIND)                              :: rsf_diag_2d_id                             !<Id returned for a real 2D array  field associated with the unstructured grid by                                                                                                     !!register_diag_field.
        integer(INT_KIND)                              :: num_diag_time_steps                        !<Number of timesteps (to simulate the model running).
        type(time_type)                                :: diag_time_start                            !<Starting time for the test.
        type(time_type)                                :: diag_time_step                             !<Time step for the test.
        logical(INT_KIND)                              :: used                                       !<Return value from send data.

        integer(INT_KIND)                              :: i                                          !<Loop variable.
        integer(INT_KIND)                              :: j                                          !<Loop variable.
        integer(INT_KIND)                              :: k,l=1                                          !<Loop variable.
        integer(INT_KIND)                              :: p                                          !<Counting variable.

       !Needed to define the 2D structured domain but never used.
        integer(INT_KIND)              :: ncontacts
        integer(INT_KIND),dimension(20) :: tile1
        integer(INT_KIND),dimension(20) :: tile2
        integer(INT_KIND),dimension(20) :: istart1
        integer(INT_KIND),dimension(20) :: iend1
        integer(INT_KIND),dimension(20) :: jstart1
        integer(INT_KIND),dimension(20) :: jend1
        integer(INT_KIND),dimension(20) :: istart2
        integer(INT_KIND),dimension(20) :: iend2
        integer(INT_KIND),dimension(20) :: jstart2
        integer(INT_KIND),dimension(20) :: jend2

        integer(INT_KIND),dimension(3)  :: npes_io_group

       !Print out a message that the test is starting.
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*)
            write(output_unit,*) "</----------------------------------------"
            write(output_unit,*) "Test create_unstructured_test_restart_file" &
                                 //" starting ..."
            write(output_unit,*)
        endif

       !Synchronize all ranks.
        call mpp_sync()

       !Make sure that valid inputs were passed in.
        if (nx .lt. 1 .or. ny .lt. 1) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
                           //" there must be at least on grid point in the" &
                           //" x- and y- dimensions.")
        endif
        if (npes .gt. nx*ny) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
                           //" the total number of ranks cannot be greater" &
                           //" than the total number of grid points in the" &
                           //" x-y plane.")
        endif
        if (num_domain_tiles_x .lt. 1 .or. num_domain_tiles_y .lt. 1) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
                           //" there must be at least on domain tile in the" &
                           //" x- and y- dimensions.")
        endif
        if (mod(nx,num_domain_tiles_x) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
                           //" the total number of grid points in the" &
                           //" x-dimension must be evenly divisible by the" &
                           //" total number of domain tiles in the" &
                           //" x-dimension.")
        endif
        if (mod(ny,num_domain_tiles_y) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
                           //" the total number of grid points in the" &
                           //" y-dimension must be evenly divisible by the" &
                           //" total number of domain tiles in the" &
                           //" y-dimension.")
        endif
        if (num_domain_tiles_x*num_domain_tiles_y .gt. npes) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
                           //" the total number of domain tiles cannot be" &
                           //" greater than the total number of ranks.")
        endif
        if (mod(npes,num_domain_tiles_x) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
                           //" the total number of ranks must be evenly" &
                           //" divisible by the total number of domain" &
                           //" tiles in the x-dimension.")
        endif
        if (mod(npes,num_domain_tiles_y) .ne. 0) then
            call mpp_error(FATAL, &
                           "create_unstructured_test_restart_file:" &
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
        ncontacts = 1
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
!write (6,*)size(tile1)
       !Define the 2D structured domain.
        call mpp_define_mosaic(global_indices, &
                               layout2D, &
                               domain_2D, &
                               num_domain_tiles, &
                               0, &
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
                           "create_unstructured_test_restart_file:" &
                           //" the number of ranks exceeds the number of" &
                           //" non-masked grid points for the unstructured" &
                           //" domain.")
        endif

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
       !> Set in namelist is "I/O tile factor".  The number of ranks that
       !! participate in I/O for a tile is equal to:
       !!
       !! num_io_ranks_on_a_tile = num_ranks_on_the_tile / "I/O tile factor".
       !!
       !!so for:
       !!
       !! io_tile_factor = 1, all of the ranks on a tile participate in the I/O
       !! io_tile_factor = 2, 1/2 of the ranks on a tile participate in the I/O
       !! io_tile_factor = 3, 1/3 of the ranks on a tile participate in the I/O
       !! ...
       !! io_tile_factor = 0 is a special case where only one rank participates
       !!                  in the I/O for a tile.
       !! io_tile_factor = 1
if (mpp_pe() == mpp_root_pe()) write(6,*) "IO_TILE_FACTOR is ",io_tile_factor
allocate(unstructured_axis_diag_id(1))
allocate(rsf_diag_1d_id(1))

       !Define the "unstructured" domain decomposition.
        call mpp_define_unstruct_domain(domain_ug, &
                                        domain_2D, &
                                        num_non_masked_grid_points_per_domain_tile, &
                                        num_land_tiles_per_non_masked_grid_point, &
                                        num_ranks_using_unstructured_grid, &
                                        io_tile_factor, &
                                        unstructured_grid_point_index_map)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      !Don't need to modify above here!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

       !Get the that will be registered for the unstructured axis. This should
       !be each rank's unstructured compute domain (I think, because a gather
       !is performed by the root of each I/O domain pelist.
        call mpp_get_UG_compute_domain(domain_ug,size=unstructured_axis_data_size)
        if(.not.allocated(unstructured_axis_data))allocate(unstructured_axis_data(unstructured_axis_data_size))
!! THIS IS A PROBLEM !!
        call mpp_get_UG_domain_grid_index(domain_ug,unstructured_axis_data)
!write(6,*)"ID:",mpp_pe()," DATA: ",unstructured_axis_data
       !Initialize the "unstructured" axis for the diagnostics.
        unstructured_axis_name = "ug_axis"

        unstructured_axis_diag_id(l) = diag_axis_init(trim(unstructured_axis_name), &
                                                   real(unstructured_axis_data), &
                                                   "none", &
                                                   "U", &
                                                   long_name="mapping indices", &
                                                   domainU=domain_ug)
   call diag_axis_add_attribute(unstructured_axis_diag_id(l),'compress','grid_xt grid_yt')

!write(6,*) "ID U",unstructured_axis_diag_id
       !Add the x-, y-, and z-axes to the restart file.  Until a bug in
       !the code is resolved, I must register the unstructured axis first.
       !Also initialize the axes for the diagnostics.
        if (.not.allocated(x_axis_data)) allocate(x_axis_data(nx))
!        if (.not.allocated(y_axis_data))allocate(y_axis_data(ny))
!! ASSUMES 4 PEs!!!
! if (mpp_pe() > 4) call error_mesg("Diag_test_unstruct","Only 4 PEs please",fatal)
     do i=1,nx
          x_axis_data(i) = real(i)
     enddo
!     if (mod(mpp_pe(),2).eq.0) then
!        do j = 1,ny/4
!            y_axis_data(j) = real(j)
!        enddo
!
!     else
!        do j = 1,ny/4
!            y_axis_data(j) = real(j+ny/4)
!        enddo
!     endif

       x_axis_diag_id = diag_axis_init("grid_xt", &
                                       x_axis_data, &
                                       "degrees", &
                                       "X", &
                                       long_name="longitude")

        if (.not.allocated(y_axis_data))allocate(y_axis_data(ny/num_domain_tiles_y))
        do i = 1,ny/num_domain_tiles_y
            y_axis_data(i) = real(i)
        enddo
       y_axis_diag_id = diag_axis_init("grid_yt", &
                                       y_axis_data, &
                                       "degrees", &
                                       "Y", &
                                       long_name="latitude")

        if (.not.allocated(z_axis_data))allocate(z_axis_data(nz))
        do i = 1,nz
            z_axis_data(i) = real(i*5.0)
        enddo
       z_axis_diag_id = diag_axis_init("zfull", &
                                       z_axis_data, &
                                       "km", &
                                       "Z", &
                                       long_name="dont look down")
!write (6,*) z_axis_diag_id

       !Define some reference test data.

       !real scalar field.
        unstructured_real_scalar_field_data_ref = 1234.5678*real(l)

       !real 1D field.
        if (.not.allocated(unstructured_real_1D_field_data_ref)) allocate(unstructured_real_1D_field_data_ref(unstructured_axis_data_size))
        do i = 1,unstructured_axis_data_size
            unstructured_real_1D_field_data_ref(i) = real(i) *real(i)+0.1*(mpp_pe()+1)
        enddo

       !real 2D field.
        if (.not.allocated(unstructured_real_2D_field_data_ref)) allocate(unstructured_real_2D_field_data_ref(unstructured_axis_data_size,nz))
        do j = 1,nz
            do i = 1,unstructured_axis_data_size
                unstructured_real_2D_field_data_ref(i,j) = real(j)+0.1*(mpp_pe()+1.0)
                                                           !-1.0*real((j-1)* &
                                                           !unstructured_axis_data_size+i) &
                                                           !+ 1.1111111*real(l)
            enddo
        enddo

       !real 3D field.
!       if(.not.allocated(unstructured_real_3D_field_data_ref) allocate(unstructured_real_3D_field_data_ref(unstructured_axis_data_size,nz,cc_axis_size))
!       do k = 1,cc_axis_size
!           do j = 1,nz
!               do i = 1,unstructured_axis_data_size
!                   unstructured_real_3D_field_data_ref(i,j,k) = -1.0*real((k-1)*nz* &
!                                                                unstructured_axis_data_size+(j-1)* &
!                                                                unstructured_axis_data_size+i) &
!                                                                + 2.2222222
!               enddo
!           enddo
!       enddo

       !integer scalar field.
        unstructured_int_scalar_field_data_ref = 7654321*L

       !integer 1D field.
        if (.not.allocated(unstructured_int_1D_field_data_ref)) allocate(unstructured_int_1D_field_data_ref(unstructured_axis_data_size))
        do i = 1,unstructured_axis_data_size
            unstructured_int_1D_field_data_ref(i) = i - 8*l
        enddo

       !integer 2D field.
        if (.not.allocated(unstructured_int_2D_field_data_ref)) allocate(unstructured_int_2D_field_data_ref(unstructured_axis_data_size,nz))
        do j = 1,nz
            do i = 1,unstructured_axis_data_size
                unstructured_int_2D_field_data_ref(i,j) = -1*((j-1)*unstructured_axis_data_size+i) + 2*L
            enddo
        enddo

     !> Latitude and Longitude
     allocate(lat(ny/num_domain_tiles_y),lon(nx))
     do i=1,nx
          lon(i) = real(i)*360.0/real(nx)
     enddo
     do j=1,ny/num_domain_tiles_y
          lat(j) = real(j)*180.8/real(ny)
     enddo

       !Add a real scalar field to the restart file.  Initialize it as a
       !diagnostic.
        unstructured_real_scalar_field_name = "unstructured_real_scalar_field_1"
        unstructured_real_scalar_field_data = unstructured_real_scalar_field_data_ref

       idlon = register_diag_field("UG_unit_test", &
                                         "lon", &
                                         (/x_axis_diag_id/),&
                                         init_time=diag_time, &
                                         long_name="E-W longitude", &
                                         units="degrees")
l=SIZE(unstructured_axis_diag_id)

       rsf_diag_id = register_diag_field("UG_unit_test", &
                                         "unstructured_real_scalar_field_data", &
                                         init_time=diag_time, &
                                         long_name="rsf_diag_1", &
                                         units="ergs")
       rsf_diag_1d_id(1) = register_diag_field("UG_unit_test", &
                                         "unstructured_real_1D_field_data", &
                                         (/unstructured_axis_diag_id(1)/),&
                                         init_time=diag_time, &
                                         long_name="ONE_D_ARRAY", &
                                         units="ergs")

       rsf_diag_2d_id = register_diag_field("UG_unit_test", &
                                         "unstructured_real_2D_field_data", &
                                         (/unstructured_axis_diag_id(1), z_axis_diag_id/),&
                                         init_time=diag_time, &
                                         long_name="TWO_D_ARRAY", &
                                         units="ergs")

       idlat = register_diag_field("UG_unit_test", &
                                         "lat", &
                                         (/y_axis_diag_id/),&
                                         init_time=diag_time, &
                                         long_name="S-N latitude", &
                                         units="degrees")


IF (l .NE. 1) THEN
  do l=2,3
   write(unstructured_1d_alt,'(a,I0)') "unstructured_real_1D",L
   rsf_diag_1d_id(L) = register_diag_field ("UG_unit_test", trim(unstructured_1d_alt),&
                                          (/unstructured_axis_diag_id(L)/),&
                                           init_time=diag_time, &
                                           long_name="OTHER"//trim(unstructured_1d_alt), &
                                           units="kg")
  enddo
ENDIF !L.ne.1
       !Add a real 1D field to the restart file.  This field is of the form:
       !field = field(unstructured).
        unstructured_real_1D_field_name = "unstructured_real_1D_field_1"
        if (.not.allocated(unstructured_real_1D_field_data)) allocate(unstructured_real_1D_field_data(unstructured_axis_data_size))
        unstructured_real_1D_field_data = unstructured_real_1D_field_data_ref

       !Add a real 2D field to the restart file.  This field is of the form:
       !field = field(unstructured,z).
        unstructured_real_2D_field_name = "unstructured_real_2D_field_1"
       if (.not.allocated(unstructured_real_2D_field_data)) allocate(unstructured_real_2D_field_data(unstructured_axis_data_size,nz))
       unstructured_real_2D_field_data = unstructured_real_2D_field_data_ref
!       allocate(unstructured_real_2D_field_data(unstructured_axis_data_size,nx))
!       unstructured_real_2D_field_data = 1

       !Add a real 3D field to the restart file.  This field is of the form:
       !field = field(unstructured,z,cc).
!       unstructured_real_3D_field_name = "unstructured_real_3D_field_1"
!       if (.not.allocated(unstructured_real_3D_field_data)) allocate(unstructured_real_3D_field_data(unstructured_axis_data_size,nz,cc_axis_size))
!       unstructured_real_3D_field_data = unstructured_real_3D_field_data_ref

       !Add an integer scalar field to the restart file.
        unstructured_int_scalar_field_name = "unstructured_int_scalar_field_1"
        unstructured_int_scalar_field_data = unstructured_int_scalar_field_data_ref

       !Add an integer 1D field to the restart file.  This field is of the
       !from: field = field(unstructured).
        unstructured_int_1D_field_name = "unstructured_int_1D_field_1"
        if (.not.allocated(unstructured_int_1D_field_data)) allocate(unstructured_int_1D_field_data(unstructured_axis_data_size))
        unstructured_int_1D_field_data = unstructured_int_1D_field_data_ref

       !Add an integer 2D field to the restart file.  This field is of the
       !form: field = field(unstructured,z).
        unstructured_int_2D_field_name = "unstructured_int_2D_field_1"
        if (.not.allocated(unstructured_int_2D_field_data)) allocate(unstructured_int_2D_field_data(unstructured_axis_data_size,nz))
        unstructured_int_2D_field_data = unstructured_int_2D_field_data_ref

       !Simulate the model timesteps, so that diagnostics may be written
       !out.
        num_diag_time_steps = 4
        diag_time_step = set_time(12*3600, 0)
        diag_time_start = diag_time
! used = send_data(idlat,lat,diag_time)
! used = send_data(idlon,lon,diag_time)
        do i = 1,num_diag_time_steps

           !Update the current time.
            diag_time = diag_time  + diag_time_step

           !"Evolve" the test data.
            unstructured_real_scalar_field_data_ref = unstructured_real_scalar_field_data_ref + &
                                                      real(1)
            unstructured_real_scalar_field_data = unstructured_real_scalar_field_data_ref

           !Update the data.
           if (rsf_diag_id .gt. 0) then
               used = send_data(rsf_diag_id, &
                                unstructured_real_scalar_field_data, &
                                diag_time)
           endif

        IF (SIZE(rsf_diag_1d_id) == 1) THEN
          used = send_data(rsf_diag_1d_id(1), &
                                unstructured_real_1D_field_data, &
                                diag_time)
         ELSE
          DO L=1,3
           used = send_data(rsf_diag_1d_id(L), &
                                unstructured_real_1D_field_data, &
                                diag_time)
          ENDDO
         ENDIF
          used = send_data(rsf_diag_2d_id, &
                                unstructured_real_2D_field_data, &
                                diag_time)
 used = send_data(idlat,lat,diag_time)
 used = send_data(idlon,lon,diag_time)

        enddo
       !Deallocate the unstructured domain.
        call mpp_sync()
!       call mpp_deallocate_domainUG(domain_ug)

       !Deallocate the 2D structured domain.
        call mpp_deallocate_domain(domain_2D)

       !Deallocate local allocatables.
        deallocate(pe_start)
        deallocate(pe_end)
        deallocate(global_indices)
        deallocate(layout2D)
        deallocate(land_mask)
        deallocate(num_non_masked_grid_points_per_domain_tile)
        deallocate(num_land_tiles_per_non_masked_grid_point)
        deallocate(unstructured_grid_point_index_map)
        deallocate(x_axis_data)
        deallocate(y_axis_data)
        deallocate(z_axis_data)
        deallocate(unstructured_axis_data)
        deallocate(unstructured_real_1D_field_data_ref)
        deallocate(unstructured_real_2D_field_data_ref)
!       deallocate(unstructured_real_3D_field_data_ref)
        deallocate(unstructured_int_1D_field_data_ref)
        deallocate(unstructured_int_2D_field_data_ref)
        deallocate(unstructured_real_1D_field_data)
        deallocate(unstructured_real_2D_field_data)
!       deallocate(unstructured_real_3D_field_data)
        deallocate(unstructured_int_1D_field_data)
        deallocate(unstructured_int_2D_field_data)



       !Print out a message that the test is done.
        call mpp_sync()
        if (mpp_pe() .eq. mpp_root_pe()) then
            write(output_unit,*)
            write(output_unit,*) "Test create_unstructured_test_restart_file" &
                                 //" complete."
            write(output_unit,*) "----------------------------------------/>"
            write(output_unit,*)
        endif


        return
  END SUBROUTINE unstruct_test

END PROGRAM test
