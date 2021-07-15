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

! The diag_table needed for the test is included below.

!--------------------------------------------------------------------------------------------------
! test_data_override
! 1 3 1 0 0 0
!
! #output files
! "test_data_override",  -1, "days", 1, "days", "time"
!
! #output variables
! "test_data_override_mod", "sst", "sst", "test_data_override",  "all", .false., "none", 2
! "test_data_override_mod", "ice", "ice", "test_data_override",  "all", .false., "none", 2
!--------------------------------------------------------------------------------------------------

! The data table needed for the test is included below.
!--------------------------------------------------------------------------------------------------
! "ICE", "sst_obs",  "SST", "INPUT/sst_ice_clim.nc", .false., 300.0
! "ICE", "sic_obs",  "SIC", "INPUT/sst_ice_clim.nc", .false., 300.0
! "OCN", "sst_obs",  "SST", "INPUT/sst_ice_clim.nc", .false., 300.0
! "LND", "sst_obs",  "SST", "INPUT/sst_ice_clim.nc", .false., 300.0
!--------------------------------------------------------------------------------------------------

program test

  ! Input data and path_names file for this program is in:
  ! /archive/pjp/unit_tests/test_data_override/lima/exp1
 use           mpp_mod, only: input_nml_file, stdout, mpp_chksum
 use   mpp_domains_mod, only: domain2d, mpp_define_domains, mpp_define_io_domain, mpp_get_compute_domain, mpp_define_layout
 use           fms_mod, only: fms_init, fms_end, mpp_npes, file_exist, check_nml_error
 use           fms_mod, only: error_mesg, FATAL, file_exist, field_exist, field_size
 use  fms_affinity_mod, only: fms_affinity_set
 use        fms_io_mod, only: read_data, fms_io_exit
 use     constants_mod, only: constants_init, pi
 use  time_manager_mod, only: time_type, set_calendar_type, set_date, NOLEAP, JULIAN, operator(+), set_time, print_time
 use  diag_manager_mod, only: diag_manager_init, diag_manager_end, register_static_field, register_diag_field
 use  diag_manager_mod, only: send_data, diag_axis_init
 use data_override_mod, only: data_override_init, data_override, data_override_UG
  use mpp_mod,         only : FATAL, WARNING, MPP_DEBUG, NOTE, MPP_CLOCK_SYNC,MPP_CLOCK_DETAILED
  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_root_pe, mpp_error, mpp_set_warn_level
  use mpp_mod,         only : mpp_declare_pelist, mpp_set_current_pelist, mpp_sync, mpp_sync_self
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id
  use mpp_mod,         only : mpp_init, mpp_exit, mpp_chksum, stdout, stderr
  use mpp_mod,         only : input_nml_file
  use mpp_mod,         only : mpp_get_current_pelist, mpp_broadcast
  use mpp_domains_mod, only : GLOBAL_DATA_DOMAIN, BITWISE_EXACT_SUM, BGRID_NE, CGRID_NE, DGRID_NE
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
  use mpp_domains_mod, only : mpp_get_domain_shift, EDGEUPDATE, mpp_deallocate_domain
  use mpp_domains_mod, only : mpp_group_update_type, mpp_create_group_update
  use mpp_domains_mod, only : mpp_do_group_update, mpp_clear_group_update
  use mpp_domains_mod, only : mpp_start_group_update, mpp_complete_group_update
  use mpp_domains_mod, only : WUPDATE, SUPDATE, mpp_get_compute_domains
  use mpp_domains_mod, only : domainUG, mpp_define_unstruct_domain, mpp_get_UG_domain_tile_id
  use mpp_domains_mod, only : mpp_get_UG_compute_domain, mpp_pass_SG_to_UG, mpp_pass_UG_to_SG
  use mpp_domains_mod, only : mpp_get_ug_global_domain, mpp_global_field_ug
  use mpp_memutils_mod, only : mpp_memuse_begin, mpp_memuse_end
  use platform_mod

 implicit none

 integer                           :: stdoutunit
 integer                           :: num_threads = 1
 integer                           :: omp_get_num_threads
 integer                           :: isw, iew, jsw, jew
 integer, allocatable              :: is_win(:), js_win(:)
 integer                           :: nx_dom, ny_dom, nx_win, ny_win
 type(domain2d)                    :: Domain
 integer                           :: nlon, nlat, siz(4)
 real, allocatable, dimension(:)   :: x, y
 real, allocatable, dimension(:,:) :: lon, lat
 real, allocatable, dimension(:,:) :: sst, ice
 integer                           :: id_x, id_y, id_lon, id_lat, id_sst, id_ice
 integer                           :: i, j, is, ie, js, je, io, ierr, n
 real                              :: rad_to_deg
 character(len=36)                 :: message
 type(time_type)                   :: Time
 logical                           :: used
 logical, allocatable              :: ov_sst(:), ov_ice(:)
 integer, dimension(2)             :: layout = (/0,0/)
 character(len=256)                :: solo_mosaic_file, tile_file
 character(len=128)                :: grid_file   = "INPUT/grid_spec.nc"
 integer                           :: window(2) = (/1,1/)
 integer                           :: nthreads=1
 integer                           :: nwindows
 integer                           :: nx_cubic=90, ny_cubic=90, nx_latlon=90, ny_latlon=90
 integer                           :: test_num=1 !* 1 for unstruct cubic grid, 2 for unstruct latlon-grid
 namelist / test_data_override_nml / layout, window, nthreads, nx_cubic, ny_cubic, nx_latlon, ny_latlon, test_num

 call fms_init
 call constants_init
 call set_calendar_type(NOLEAP)
 call diag_manager_init

 rad_to_deg = 180./pi

 read (input_nml_file, test_data_override_nml, iostat=io)
 ierr = check_nml_error(io, 'test_data_override_nml')

 if(field_exist(grid_file, "x_T" ) ) then
    call field_size(grid_file, 'x_T', siz)
    nlon = siz(1)
    nlat = siz(2)
 else if(field_exist(grid_file, "geolon_t" ) ) then
    call field_size(grid_file, 'geolon_t', siz)
    nlon = siz(1)
    nlat = siz(2)
 else if (field_exist(grid_file, "ocn_mosaic_file" )) then
    call read_data(grid_file, 'ocn_mosaic_file', solo_mosaic_file)
    solo_mosaic_file = 'INPUT/'//trim(solo_mosaic_file)
    call field_size(solo_mosaic_file, 'gridfiles', siz)
    if( siz(2) .NE. 1) call error_mesg('test_data_override', 'only support single tile mosaic, contact developer', FATAL)
    call read_data(solo_mosaic_file, 'gridfiles', tile_file)
    tile_file = 'INPUT/'//trim(tile_file)
    call field_size(tile_file, 'area', siz)
    if(mod(siz(1),2) .NE. 0 .OR. mod(siz(2),2) .NE. 0 ) call error_mesg('test_data_override', &
        "test_data_override: supergrid size can not be divided by 2", FATAL)
    nlon = siz(1)/2
    nlat = siz(2)/2
 else
    call error_mesg('test_data_override', 'x_T, geolon_t and ocn_mosaic_file does not exist', FATAL)
 end if

 if(layout(1)*layout(2) .NE. mpp_npes() ) then
    call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes(), layout )
 end if


 call mpp_define_domains( (/1,nlon,1,nlat/), layout, Domain, name='test_data_override')
 call mpp_define_io_domain(Domain, (/1,1/))
 call data_override_init(Ice_domain_in=Domain, Ocean_domain_in=Domain)
 call data_override_init(Ice_domain_in=Domain, Ocean_domain_in=Domain)
 call mpp_get_compute_domain(Domain, is, ie, js, je)
 call get_grid

 allocate(x(nlon), y(nlat))

 do i=1,nlon
   x(i) = i
 enddo
 do j=1,nlat
   y(j) = j
 enddo

 Time = set_date(2,1,1,0,0,0)

 allocate(sst(is:ie,js:je), ice(is:ie,js:je))
 sst = 0
 ice = 0

 id_x  = diag_axis_init('x',  x,  'point_E', 'x', long_name='point_E', Domain2=Domain)
 id_y  = diag_axis_init('y',  y,  'point_N', 'y', long_name='point_N', Domain2=Domain)

 Time = Time + set_time(0,1)

 id_lon = register_static_field('test_data_override_mod', 'lon', (/id_x,id_y/), 'longitude', 'Degrees')
 id_lat = register_static_field('test_data_override_mod', 'lat', (/id_x,id_y/), 'longitude', 'Degrees')
 id_sst = register_diag_field  ('test_data_override_mod', 'sst', (/id_x,id_y/), Time, 'SST', 'K')
 id_ice = register_diag_field  ('test_data_override_mod', 'ice', (/id_x,id_y/), Time, 'ICE', ' ')
 used = send_data(id_lon, lon, Time)
 used = send_data(id_lat, lat, Time)

 sst = 0.
 ice = 0.

! get number of threads

nx_dom = ie - is + 1
ny_dom = je - js + 1
if( mod( nx_dom, window(1) ) .NE. 0 ) call error_mesg('test_data_override', &
        "nx_dom is not divisible by window(1)", FATAL)
if( mod( ny_dom, window(2) ) .NE. 0 ) call error_mesg('test_data_override', &
        "ny_dom is not divisible by window(2)", FATAL)

nwindows = window(1)*window(2)
!$ call omp_set_num_threads(nthreads)
!$OMP PARALLEL
!$ call fms_affinity_set("test_data_override", .FALSE., omp_get_num_threads() )
!$OMP END PARALLEL

nx_win = nx_dom/window(1)
ny_win = ny_dom/window(2)
allocate(is_win(nwindows), js_win(nwindows))

i = 1
do jsw = js,je,ny_win
   do isw = is,ie,nx_win
      is_win(i) = isw
      js_win(i) = jsw
      i = i + 1
   enddo
enddo

allocate(ov_sst(nwindows), ov_ice(nwindows))
!$OMP parallel  do schedule(static) default(shared) private(isw, iew, jsw, jew)
do n = 1, nwindows
   isw = is_win(n)
   iew = isw + nx_win - 1
   jsw = js_win(n)
   jew = jsw + ny_win - 1
   call data_override('OCN','sst_obs',sst(isw:iew,jsw:jew),Time,override=ov_sst(n), &
                      is_in=isw-is+1, ie_in=iew-is+1, js_in=jsw-js+1, je_in=jew-js+1)
   call data_override('ICE', 'sic_obs', ice(isw:iew,jsw:jew), Time, override=ov_ice(n), &
                      is_in=isw-is+1, ie_in=iew-is+1, js_in=jsw-js+1, je_in=jew-js+1)
enddo

 if(ANY(.NOT. ov_sst) .or. ANY(.not.ov_ice)) then
   if(ANY(.NOT. ov_sst)) then
     message = 'override failed for ice'
   else if(ANY(.NOT. ov_ice)) then
     message = 'override failed for sst'
   else
     message = 'override failed for both sst and ice'
   endif
   call error_mesg('test_data_override', trim(message), FATAL)
 endif

 stdoutunit = stdout()
 write(stdoutunit,*)"===>NOTE from test_data_override: sst chksum = ", mpp_chksum(sst)
 write(stdoutunit,*)"===>NOTE from test_data_override: ice chksum = ", mpp_chksum(ice)


 if(id_sst > 0) used = send_data(id_sst, sst, Time)
 if(id_ice > 0) used = send_data(id_ice, ice, Time)

 if(test_num == 1 .and. nx_cubic > 0 .and. ny_cubic > 0) then
    call mpp_memuse_begin()
    if(mpp_pe() == mpp_root_pe() ) print *, '---------> Testing cubic grid'
    call test_unstruct_grid( 'Cubic-Grid', Time )
    call mpp_memuse_end('Cubic-grid')
    if(mpp_pe() == mpp_root_pe() ) print *, '---------> Done testing cubic grid'
 endif
 if(test_num ==2 .and. nx_latlon > 0 .and. ny_latlon > 0) then
    call mpp_memuse_begin
    if(mpp_pe() == mpp_root_pe() ) print *, '---------> Testing latlon-grid'
    call test_unstruct_grid( 'Latlon-Grid', Time )
    call mpp_memuse_end('Latlon-Grid')
    if(mpp_pe() == mpp_root_pe() ) print *, '---------> Finish testing latlon-grid'
 endif




!-------------------------------------------------------------------------------------------------------
! What follows is a test of calendar conversion

!Time = set_date(1980,2,27,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)

!Time = set_date(1980,2,28,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)

!Time = set_date(1980,2,29,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)

!Time = set_date(1980,3,1,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)

!Time = set_date(1980,3,2,0,0,0)
!call print_time(Time)
!call data_override('OCN','sst_obs',sst,Time)
!if(id_sst > 0) used = send_data(id_sst, sst, Time)
!-------------------------------------------------------------------------------------------------------

 call diag_manager_end(Time)
 call fms_io_exit
 call fms_end

contains

!=================================================================================================================================
 subroutine get_grid
   real, allocatable, dimension(:,:,:) :: lon_vert_glo, lat_vert_glo
   real, allocatable, dimension(:,:)   :: lon_global, lat_global
   integer, dimension(4)  :: siz
   character(len=128) :: message


   if(field_exist(grid_file, 'x_T')) then
      call field_size(grid_file, 'x_T', siz)
      if(siz(1) /= nlon .or. siz(2) /= nlat) then
         write(message,'(a,2i4)') 'x_T is wrong shape. shape(x_T)=',siz(1:2)
         call error_mesg('test_data_override', trim(message), FATAL)
      endif
      allocate(lon_vert_glo(nlon,nlat,4), lat_vert_glo(nlon,nlat,4) )
      allocate(lon_global  (nlon,nlat  ), lat_global  (nlon,nlat  ) )
      call read_data(trim(grid_file), 'x_vert_T', lon_vert_glo, no_domain=.true.)
      call read_data(trim(grid_file), 'y_vert_T', lat_vert_glo, no_domain=.true.)
      lon_global(:,:)  = (lon_vert_glo(:,:,1) + lon_vert_glo(:,:,2) + lon_vert_glo(:,:,3) + lon_vert_glo(:,:,4))*0.25
      lat_global(:,:) =  (lat_vert_glo(:,:,1) + lat_vert_glo(:,:,2) + lat_vert_glo(:,:,3) + lat_vert_glo(:,:,4))*0.25
   else  if(field_exist(grid_file, "geolon_t" ) ) then
      call field_size(grid_file, 'geolon_vert_t', siz)
      if(siz(1) /= nlon+1 .or. siz(2) /= nlat+1) then
         write(message,'(a,2i4)') 'geolon_vert_t is wrong shape. shape(geolon_vert_t)=',siz(1:2)
         call error_mesg('test_data_override', trim(message), FATAL)
      endif
      allocate(lon_vert_glo(nlon+1,nlat+1,1), lat_vert_glo(nlon+1,nlat+1,1))
      allocate(lon_global  (nlon,  nlat    ), lat_global  (nlon,  nlat    ))
      call read_data(trim(grid_file), 'geolon_vert_t', lon_vert_glo, no_domain=.true.)
      call read_data(trim(grid_file), 'geolat_vert_t', lat_vert_glo, no_domain=.true.)

      do i = 1, nlon
         do j = 1, nlat
            lon_global(i,j) = (lon_vert_glo(i,j,1) + lon_vert_glo(i+1,j,1) + &
                 lon_vert_glo(i+1,j+1,1) + lon_vert_glo(i,j+1,1))*0.25
            lat_global(i,j) = (lat_vert_glo(i,j,1) + lat_vert_glo(i+1,j,1) + &
                 lat_vert_glo(i+1,j+1,1) + lat_vert_glo(i,j+1,1))*0.25
         enddo
      enddo
   else if( field_exist(grid_file, "ocn_mosaic_file") ) then ! reading from mosaic file
      call field_size(tile_file, 'area', siz)
      if(siz(1) /= nlon*2 .or. siz(2) /= nlat*2) then
         write(message,'(a,2i4)') 'area is wrong shape. shape(area)=',siz(1:2)
         call error_mesg('test_data_override', trim(message), FATAL)
      endif
      allocate(lon_vert_glo(siz(1)+1,siz(2)+1,1), lat_vert_glo(siz(1)+1,siz(2)+1,1))
      allocate(lon_global  (nlon,  nlat    ), lat_global  (nlon,  nlat    ))
      call read_data( tile_file, 'x', lon_vert_glo, no_domain=.true.)
      call read_data( tile_file, 'y', lat_vert_glo, no_domain=.true.)
      do j = 1, nlat
         do i = 1, nlon
            lon_global(i,j) = lon_vert_glo(i*2,j*2,1)
            lat_global(i,j) = lat_vert_glo(i*2,j*2,1)
         end do
      end do
   end if

  allocate(lon(is:ie,js:je), lat(is:ie,js:je))
  lon = lon_global(is:ie,js:je)
  lat = lat_global(is:ie,js:je)

  deallocate(lon_vert_glo)
  deallocate(lat_vert_glo)
  deallocate(lon_global)
  deallocate(lat_global)

 end subroutine get_grid

  subroutine test_unstruct_grid( type, Time )

   character(len=*), intent(in) :: type
  type(time_type),              intent(in) :: Time !(target) model time

  integer :: pe, npes
  integer :: nx, ny, nz=40, stackmax=4000000
  integer :: stdunit = 6
  logical :: debug=.FALSE., opened

  integer :: mpes = 0
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2
  character(len=32) :: warn_level = "fatal"
  integer :: layout_cubic(2) = (/0,0/)
  integer :: layout_tripolar(2) = (/0,0/)
  integer :: layout_ensemble(2) = (/0,0/)
  logical :: do_sleep = .false.
  integer :: num_iter = 1
  integer :: num_fields = 4

    type(domain2D) :: SG_domain
    type(domainUG) :: UG_domain
    integer        :: num_contact, ntiles, npes_per_tile
    integer        :: isc, iec, jsc, jec, isd, ied, jsd, jed
    integer        :: ism, iem, jsm, jem, lsg, leg

    integer, allocatable, dimension(:)       :: pe_start, pe_end, npts_tile, grid_index, ntiles_grid
    integer, allocatable, dimension(:,:)     :: layout2D, global_indices
    real,    allocatable, dimension(:,:)     :: x1, x2, g1, g2
    real,    allocatable, dimension(:,:,:)   :: a1, a2, gdata
    real,    allocatable, dimension(:,:)     :: rmask
    real,    allocatable, dimension(:)       :: frac_crit
    logical, allocatable, dimension(:,:,:)   :: lmask,msk
    integer, allocatable, dimension(:)       :: isl, iel, jsl, jel
    character(len=3)   :: text
    integer            :: tile
    integer            :: ntotal_land, istart, iend, pos
    integer            :: outunit, errunit, k, l

  npes = mpp_npes()

  outunit = stdout()
  errunit = stderr()
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
       allocate(layout2D(2,ntiles), global_indices(4,ntiles), pe_start(ntiles), pe_end(ntiles) )
       do n = 1, ntiles
          pe_start(n) = (n-1)*npes_per_tile
          pe_end(n)   = n*npes_per_tile-1
       end do

       do n = 1, ntiles
          global_indices(:,n) = (/1,nx,1,ny/)
          layout2D(:,n)         = layout
       end do

       call define_cubic_mosaic(type, SG_domain, (/nx,nx,nx,nx,nx,nx/), (/ny,ny,ny,ny,ny,ny/), &
            global_indices, layout2D, pe_start, pe_end )
    case ( 'Latlon-Grid' )
       if(nx_latlon == 0 .OR. ny_latlon == 0 ) then
          call mpp_error(NOTE,'test_unstruct_update: for latlon mosaic, nx_latlon and ny_latlon are zero, '//&
               'No test is done for Lalton-Grid mosaic. ' )
          return
       endif
       nx = nx_latlon
       ny = ny_latlon
       ntiles = 1
       npes_per_tile = npes
       allocate(frac_crit(ntiles))
       frac_crit(1) = 0.3
       call mpp_define_layout((/1,nx,1,ny/), npes, layout)
       call mpp_define_domains((/1,nx,1,ny/), layout, SG_domain, xflags = cyclic_global_domain)
    case default
       call mpp_error(FATAL, 'test_group_update: no such test: '//type)
    end select

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

    !--- test the 2-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,1), a2(isc:iec,jsc:jec,1 ) )
    allocate(msk(isc:iec, jsc:jec,1)); msk = .false.

    tile = mpp_pe()/npes_per_tile + 1
    do j = jsc, jec
       do i = isc, iec
          msk(i,j,1) = lmask(i,j,tile)
       enddo
    enddo
    !First override the test SG data from file/field
    call data_override_init(Land_domain_in=SG_domain)
    call data_override('LND','sst_obs',a1(:,:,1),Time)

    !Create the test UG data
    a2 = -9999
    do j = jsc, jec
       do i = isc, iec
          if(.NOT. msk(i,j,1)) a2(i,j,1)=a1(i,j,1)
       enddo
    enddo

    allocate(x1(istart:iend,1), x2(istart:iend,1))
    x1 = -99999
    x2 = -999999
    !--- fill the value of x2

    !Now override the test UG data from the same file/field
    call data_override_init(Land_domainUG_in=UG_domain)
    call data_override_UG('LND','sst_obs',x2(:,1),Time)

    !Ensure you get the same UG data from the SG data
    call mpp_pass_SG_to_UG(UG_domain, a1(:,:,1), x1(:,1))
    call compare_checksums_2D(x1, x2, type//' SG2UG 2-D compute domain')

    !Ensure you get the same SG data from the UG data if you transform back
    call mpp_pass_UG_to_SG(UG_domain, x1(:,1), a2(:,:,1))
    call compare_checksums(a1(:,:,1:1),a2(:,:,1:1),type//' UG2SG 2-D compute domain')

    deallocate(a1,a2,x1,x2)

    !--- test the 3-D data is on computing domain
    allocate( a1(isc:iec, jsc:jec,nz), a2(isc:iec,jsc:jec,nz ) )

!    tile = mpp_pe()/npes_per_tile + 1
!    do k = 1, nz
!       do j = jsc, jec
!          do i = isc, iec
!             a1(i,j,k) = gdata(i,j,tile)
!             if(a1(i,j,k) .NE. -999) a1(i,j,k) = a1(i,j,k) + k*1.e-6
!          enddo
!       enddo
!    enddo

    !First override the test SG data from file/field
    call data_override_init(Land_domain_in=SG_domain)
    call data_override('LND','sst_obs',a1,Time)

    a2 = -999
    !For this test on non-Land points a2 must match a1
    do k = 1, nz
    do j = jsc, jec
       do i = isc, iec
          if(.NOT. msk(i,j,1)) a2(i,j,k)=a1(i,j,k)
       enddo
    enddo
    enddo

    allocate(x1(istart:iend,nz), x2(istart:iend,nz))
    x1 = -999
    x2 = -999
    !--- fill the value of x2
 !   tile = mpp_get_UG_domain_tile_id(UG_domain)
 !   pos = 0
 !   do n = 1, tile-1
 !      pos = pos + npts_tile(n)
 !   enddo
 !   do l = istart, iend
 !      i = mod((grid_index(pos+l)-1), nx) + 1
 !      j = (grid_index(pos+l)-1)/nx + 1
 !      do k = 1, nz
 !         x2(l,k) = gdata(i,j,tile) + k*1.e-6
 !      enddo
 !   enddo

    !Now override the test UG data from the same file/field
    call data_override_init(Land_domainUG_in=UG_domain)
    call data_override_UG('LND','sst_obs',x2,Time)

    !Ensure you get the same UG data from the SG data
    call mpp_pass_SG_to_UG(UG_domain, a1, x1)
    call compare_checksums_2D(x1, x2, type//' SG2UG 3-D compute domain')
    !Ensure you get the same SG data from the UG data if you transform back
    call mpp_pass_UG_to_SG(UG_domain, x1, a2)
    call compare_checksums(a1,a2,type//' UG2SG 3-D compute domain')
    deallocate(a1,a2,x1,x2)


  end subroutine test_unstruct_grid

  subroutine compare_checksums( a, b, string )
    real, intent(in), dimension(:,:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(i8_kind) :: sum1, sum2
    integer :: i, j, k,pe

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    ! mpp_sync()
    call mpp_sync_self()
  pe = mpp_pe()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) .or. size(a,3) .ne. size(b,3) ) &
         call mpp_error(FATAL,'compare_chksum: size of a and b does not match')

!print*, "pe,a(1,1,1),b(1,1,1)=", pe,a(1,1,1),b(1,1,1)

    do k = 1, size(a,3)
       do j = 1, size(a,2)
          do i = 1, size(a,1)
             if(a(i,j,k) .ne. b(i,j,k)) then
            print*, "pe,i,j,k", pe,i,j,k
            print*, "a =", a(i,j,k)
            print*, "b =", b(i,j,k)
!                write(stdunit,'(a,i3,a,i3,a,i3,a,i3,a,f20.9,a,f20.9)')" at pe ", mpp_pe(), &
!                     ", at point (",i,", ", j, ", ", k, "), a = ", a(i,j,k), ", b = ", b(i,j,k)
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
        call mpp_error( FATAL, trim(string)//': chksums are not OK.' )
    end if
  end subroutine compare_checksums

  !###########################################################################
  subroutine compare_checksums_2D( a, b, string )
    real, intent(in), dimension(:,:) :: a, b
    character(len=*), intent(in) :: string
    integer(i8_kind) :: sum1, sum2
    integer :: i, j,pe

    ! z1l can not call mpp_sync here since there might be different number of tiles on each pe.
    ! mpp_sync()
    call mpp_sync_self()
  pe = mpp_pe()

    if(size(a,1) .ne. size(b,1) .or. size(a,2) .ne. size(b,2) ) &
         call mpp_error(FATAL,'compare_chksum_2D: size of a and b does not match')

!print*, "a(1,1),b(1,1)=", a(1,1),b(1,1)

    do j = 1, size(a,2)
       do i = 1, size(a,1)
          if(a(i,j) .ne. b(i,j)) then
            print*, "i,j= ", i,j
            print*, "a =", a(i,j)
            print*, "b =", b(i,j)
!             write(stdunit,'(a,i3,a,i3,a,i3,a,f20.9,a,f20.9)')"at pe ", mpp_pe(), &
!                  ", at point (",i,", ", j, "),a=", a(i,j), ",b=", b(i,j)
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
  subroutine define_cubic_mosaic(type, domain, ni, nj, global_indices, layout, pe_start, pe_end)
    character(len=*), intent(in)  :: type
    type(domain2d), intent(inout) :: domain
    integer,        intent(in)    :: global_indices(:,:), layout(:,:)
    integer,        intent(in)    :: ni(:), nj(:)
    integer,        intent(in)    :: pe_start(:), pe_end(:)
    integer, dimension(12)        :: istart1, iend1, jstart1, jend1, tile1
    integer, dimension(12)        :: istart2, iend2, jstart2, jend2, tile2
    integer                       :: ntiles, num_contact, msize(2)
  integer :: whalo = 2, ehalo = 2, shalo = 2, nhalo = 2


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

    return

  end subroutine define_cubic_mosaic

!=================================================================================================================================
 end program test
