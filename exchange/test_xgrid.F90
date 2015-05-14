#ifdef TEST_XGRID
! Now only test some simple test, will test cubic grid mosaic in the future.

program xgrid_test

  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_error, FATAL, mpp_chksum, mpp_min, mpp_max
  use mpp_mod,         only : mpp_set_current_pelist, mpp_declare_pelist
  use mpp_domains_mod, only : mpp_define_domains, mpp_define_layout, mpp_domains_exit
  use mpp_domains_mod, only : mpp_get_compute_domain, domain2d, mpp_domains_init
  use mpp_domains_mod, only : mpp_define_mosaic_pelist, mpp_define_mosaic, mpp_global_sum
  use mpp_domains_mod, only : mpp_get_data_domain, mpp_get_global_domain, mpp_update_domains
  use mpp_io_mod,      only : mpp_open, MPP_RDONLY,MPP_NETCDF, MPP_MULTI, MPP_SINGLE, mpp_close
  use mpp_io_mod,      only : mpp_get_att_value
  use fms_mod,         only : fms_init, file_exist, field_exist, field_size, open_namelist_file
  use fms_mod,         only : check_nml_error, close_file, read_data, stdout, fms_end
  use fms_mod,         only : get_mosaic_tile_grid, write_data, set_domain
  use fms_io_mod,      only : fms_io_exit
  use constants_mod,   only : DEG_TO_RAD
  use xgrid_mod,       only : xgrid_init, setup_xmap, put_to_xgrid, get_from_xgrid
  use xgrid_mod,       only : xmap_type, xgrid_count, grid_box_type, SECOND_ORDER
  use xgrid_mod,       only : get_xmap_grid_area, set_frac_area
  use mosaic_mod,      only : get_mosaic_ntiles, get_mosaic_grid_sizes
  use mosaic_mod,      only : get_mosaic_ncontacts, get_mosaic_contact
  use gradient_mod,    only : calc_cubic_grid_info
  use ensemble_manager_mod, only : ensemble_manager_init, ensemble_pelist_setup
  use ensemble_manager_mod, only : get_ensemble_size

implicit none

  real, parameter :: EPSLN = 1.0e-10
  character(len=256) :: atm_input_file  = "INPUT/atmos_input.nc"
  character(len=256) :: atm_output_file = "atmos_output.nc"
  character(len=256) :: lnd_output_file = "land_output.nc"
  character(len=256) :: ice_output_file = "ocean_output.nc"
  character(len=256) :: atm_field_name  = "none"

  character(len=256) :: runoff_input_file  = "INPUT/land_runoff.nc"
  character(len=256) :: runoff_output_file  = "land_runoff.nc"
  character(len=256) :: runoff_field_name  = "none"
  integer            :: num_iter           = 0 
  integer            :: nk_lnd = 1, nk_ice = 1
  integer            :: atm_layout(2) = (/0,0/)
  integer            :: lnd_layout(2) = (/0,0/)
  integer            :: ice_layout(2) = (/0,0/)
  integer            :: atm_nest_layout(2) = (/0,0/)
  integer            :: atm_npes = 0
  integer            :: lnd_npes = 0
  integer            :: ice_npes = 0
  integer            :: ocn_npes = 0
  integer            :: atm_nest_npes = 0
  logical            :: concurrent = .false.

  namelist /xgrid_test_nml/ atm_input_file, atm_field_name, runoff_input_file, runoff_field_name, num_iter, &
                            nk_lnd, nk_ice, atm_layout, ice_layout, lnd_layout, atm_nest_layout, &
                            atm_nest_npes, atm_npes, lnd_npes, ice_npes

  integer              :: remap_method
  integer              :: pe, npes, ierr, nml_unit, io, n
  integer              :: siz(4), ntile_lnd, ntile_atm, ntile_ice, ncontact
  integer, allocatable :: layout(:,:), global_indices(:,:)
  integer, allocatable :: atm_nx(:), atm_ny(:), ice_nx(:), ice_ny(:), lnd_nx(:), lnd_ny(:)
  integer, allocatable :: pe_start(:), pe_end(:), dummy(:)
  integer, allocatable :: istart1(:), iend1(:), jstart1(:), jend1(:), tile1(:)
  integer, allocatable :: istart2(:), iend2(:), jstart2(:), jend2(:), tile2(:)
  character(len=256)   :: grid_file = "INPUT/grid_spec.nc"
  character(len=256)   :: atm_mosaic, ocn_mosaic, lnd_mosaic
  character(len=256)   :: atm_mosaic_file, ocn_mosaic_file, lnd_mosaic_file
  character(len=256)   :: grid_descriptor, tile_file
  integer              :: isc_atm, iec_atm, jsc_atm, jec_atm, nxc_atm, nyc_atm
  integer              :: isc_lnd, iec_lnd, jsc_lnd, jec_lnd
  integer              :: isc_ice, iec_ice, jsc_ice, jec_ice
  integer              :: isd_atm, ied_atm, jsd_atm, jed_atm
  integer              :: unit, i, j, nxa, nya, nxgrid, nxl, nyl, out_unit, k
  type(domain2d)       :: Atm_domain, Ice_domain, Lnd_domain
  type(xmap_type)      :: Xmap, Xmap_runoff
  type(grid_box_type)  :: atm_grid
  real, allocatable    :: xt(:,:), yt(:,:)  ! on T-cell data domain
  real, allocatable    :: xc(:,:), yc(:,:)  ! on C-cell compute domain
  real, allocatable    :: tmpx(:,:), tmpy(:,:)
  real, allocatable    :: atm_data_in(:,:), atm_data_out(:,:)
  real, allocatable    :: atm_data_out_1(:,:), atm_data_out_2(:,:), atm_data_out_3(:,:)
  real, allocatable    :: lnd_data_out(:,:,:), ice_data_out(:,:,:)
  real, allocatable    :: runoff_data_in(:,:), runoff_data_out(:,:,:)
  real, allocatable    :: atm_area(:,:), lnd_area(:,:), ice_area(:,:)
  real, allocatable    :: lnd_frac(:,:,:), ice_frac(:,:,:)
  real, allocatable    :: x_1(:), x_2(:), x_3(:), x_4(:)
  real                 :: sum_atm_in, sum_ice_out, sum_lnd_out, sum_atm_out
  real                 :: sum_runoff_in, sum_runoff_out, tot
  real                 :: min_atm_in, max_atm_in, min_atm_out, max_atm_out
  real                 :: min_x, max_x
  logical              :: atm_input_file_exist, runoff_input_file_exist
  integer              :: npes_per_tile
  integer              :: id_put_side1_to_xgrid, id_get_side1_from_xgrid
  integer              :: id_put_side2_to_xgrid, id_get_side2_from_xgrid
  integer              :: ens_siz(6), ensemble_size
  integer              :: atm_root_pe, lnd_root_pe, ocn_root_pe, ice_root_pe, atm_nest_root_pe
  integer              :: atm_global_npes, ntile_atm_global, ncontact_global
  integer              :: atm_nest_nx, atm_nest_ny, ntile_atm_nest
  logical              :: atm_pe, lnd_pe, ice_pe, ocn_pe, atm_global_pe, atm_nest_pe
  integer, allocatable :: atm_pelist(:), ocn_pelist(:), ice_pelist(:), lnd_pelist(:)
  integer, allocatable :: atm_global_pelist(:), atm_nest_pelist(:)
  integer, allocatable :: tile_id(:)

  call fms_init

  call mpp_domains_init

  call xgrid_init(remap_method)
  call ensemble_manager_init() 

  npes     = mpp_npes()
  pe       = mpp_pe()
  out_unit = stdout()

#ifdef INTERNAL_FILE_NML
      read (input_nml_file, xgrid_test_nml, iostat=io)
      ierr = check_nml_error(io, 'xgrid_test_nml')
#else
  if (file_exist('input.nml')) then
     ierr=1
     nml_unit = open_namelist_file()
     do while (ierr /= 0)
        read(nml_unit, nml=xgrid_test_nml, iostat=io, end=10)
        ierr = check_nml_error(io, 'xgrid_test_nml')
     enddo
10   call close_file(nml_unit)
  endif
#endif

  !--- get ensemble size
  ens_siz = get_ensemble_size()   
  ensemble_size = ens_siz(1)      
  npes = ens_siz(2)              

  if( atm_npes == 0 ) atm_npes = mpp_npes()
  if( lnd_npes == 0 ) lnd_npes = atm_npes
  if( ocn_npes == 0 ) ocn_npes = atm_npes
  if( ice_npes == 0 ) ice_npes = atm_npes
  if(lnd_npes > atm_npes) call mpp_error(FATAL, 'xgrid_test: lnd_npes > atm_npes')
  if(ocn_npes > atm_npes) call mpp_error(FATAL, 'xgrid_test: ocn_npes > atm_npes')
  
  write(out_unit, *) "NOTE from test_xgrid: atm_npes = ", atm_npes
  write(out_unit, *) "NOTE from test_xgrid: atm_nest_npes = ", atm_nest_npes
  write(out_unit, *) "NOTE from test_xgrid: lnd_npes = ", lnd_npes
  write(out_unit, *) "NOTE from test_xgrid: ice_npes = ", ice_npes
  write(out_unit, *) "NOTE from test_xgrid: ocn_npes = ", ocn_npes
  
  allocate(atm_pelist(atm_npes))
  allocate(ice_pelist(ice_npes))
  allocate(lnd_pelist(lnd_npes))
  allocate(ocn_pelist(ocn_npes))
  concurrent = .false.

  call ensemble_pelist_setup(concurrent, atm_npes, ocn_npes, lnd_npes, ice_npes, &
                             atm_pelist, ocn_pelist, lnd_pelist, ice_pelist)
 
  if(atm_nest_npes > 0) then
    atm_global_npes = atm_npes - atm_nest_npes
    allocate(atm_global_pelist(atm_global_npes))
    allocate(atm_nest_pelist(atm_nest_npes))
    atm_global_pelist(1:atm_global_npes) = atm_pelist(1:atm_global_npes)
    atm_nest_pelist(1:atm_nest_npes)     = atm_pelist(atm_global_npes+1:atm_npes)
    atm_global_pe = ANY(atm_global_pelist == mpp_pe())
    atm_nest_pe = ANY(atm_nest_pelist == mpp_pe())
    atm_nest_root_pe = atm_nest_pelist(1)
    call mpp_declare_pelist(atm_global_pelist, "atm global pelist")
    call mpp_declare_pelist(atm_nest_pelist, "atm nest pelist")
  else
    atm_global_npes = atm_npes
    allocate(atm_global_pelist(atm_global_npes))
    atm_global_pelist(1:atm_global_npes) = atm_pelist(1:atm_global_npes)
    atm_global_pe = ANY(atm_global_pelist == mpp_pe())
    atm_nest_pe = .FALSE.
  endif

  atm_pe = ANY(atm_pelist   .EQ. mpp_pe()) 
  lnd_pe = ANY(lnd_pelist   .EQ. mpp_pe()) 
  ocn_pe = ANY(ocn_pelist   .EQ. mpp_pe()) 
  ice_pe = ANY(ice_pelist   .EQ. mpp_pe())
  atm_root_pe = atm_pelist(1)
  lnd_root_pe = lnd_pelist(1)
  ice_root_pe = ice_pelist(1)
  ocn_root_pe = ocn_pelist(1)

  ntile_atm = 1
  ntile_ice = 1
  ntile_lnd = 1

  if(field_exist(grid_file, "AREA_ATM" ) ) then
     if( atm_nest_npes > 0 ) call mpp_error(FATAL,  &
           'xgrid_test: nested atmosphere model is only supported for mosaic grid')
     allocate(atm_nx(1), atm_ny(1))
     allocate(lnd_nx(1), lnd_ny(1))
     allocate(ice_nx(1), ice_ny(1))
     call field_size(grid_file, "AREA_ATM", siz )
     atm_nx = siz(1); atm_ny = siz(2)
     call field_size(grid_file, "AREA_OCN", siz )
     ice_nx = siz(1); ice_ny = siz(2)
     call field_size(grid_file, "AREA_LND", siz )
     lnd_nx = siz(1); lnd_ny = siz(2)
     if( atm_layout(1)*atm_layout(2) .NE. npes ) then
        call mpp_define_layout( (/1,atm_nx,1,atm_ny/), npes, atm_layout) 
     endif
     call mpp_define_domains( (/1,atm_nx,1,atm_ny/), atm_layout, Atm_domain, name="atmosphere")
     if( lnd_layout(1)*lnd_layout(2) .NE. npes ) then
        call mpp_define_layout( (/1,lnd_nx,1,lnd_ny/), npes, lnd_layout) 
     endif
     call mpp_define_domains( (/1,lnd_nx,1,lnd_ny/), lnd_layout, Lnd_domain, name="land") 
     if( ice_layout(1)*ice_layout(2) .NE. npes ) then
        call mpp_define_layout( (/1,ice_nx,1,ice_ny/), npes, ice_layout)
     endif
     call mpp_define_domains( (/1,ice_nx,1,ice_ny/), ice_layout, Ice_domain, name="Ice")
  else if (field_exist(grid_file, "atm_mosaic" ) ) then
     !--- Get the mosaic data of each component model 
     call read_data(grid_file, 'atm_mosaic', atm_mosaic)
     call read_data(grid_file, 'lnd_mosaic', lnd_mosaic)
     call read_data(grid_file, 'ocn_mosaic', ocn_mosaic)
     atm_mosaic_file = 'INPUT/'//trim(atm_mosaic)//'.nc'
     lnd_mosaic_file = 'INPUT/'//trim(lnd_mosaic)//'.nc'
     ocn_mosaic_file = 'INPUT/'//trim(ocn_mosaic)//'.nc'

     ntile_lnd = get_mosaic_ntiles(lnd_mosaic_file);
     ntile_ice = get_mosaic_ntiles(ocn_mosaic_file);
     ntile_atm = get_mosaic_ntiles(atm_mosaic_file);

     if(ntile_ice > 1) call mpp_error(FATAL,  &
           'xgrid_test: there is more than one tile in ocn_mosaic, which is not implemented yet')

     write(out_unit,*)" There is ", ntile_atm, " tiles in atmos mosaic"
     write(out_unit,*)" There is ", ntile_lnd, " tiles in land  mosaic"
     write(out_unit,*)" There is ", ntile_ice, " tiles in ocean mosaic"
     allocate(atm_nx(ntile_atm), atm_ny(ntile_atm))
     allocate(lnd_nx(ntile_lnd), lnd_ny(ntile_lnd))
     allocate(ice_nx(ntile_ice), ice_ny(ntile_ice))

     call get_mosaic_grid_sizes(atm_mosaic_file, atm_nx, atm_ny)
     call get_mosaic_grid_sizes(lnd_mosaic_file, lnd_nx, lnd_ny)
     call get_mosaic_grid_sizes(ocn_mosaic_file, ice_nx, ice_ny)

     ncontact = get_mosaic_ncontacts(atm_mosaic_file)

     if(ncontact > 0) then
        allocate(tile1(ncontact),   tile2(ncontact) )
        allocate(istart1(ncontact), iend1(ncontact) )
        allocate(jstart1(ncontact), jend1(ncontact) )
        allocate(istart2(ncontact), iend2(ncontact) )
        allocate(jstart2(ncontact), jend2(ncontact) )
        call get_mosaic_contact( atm_mosaic_file, tile1, tile2, istart1, iend1, jstart1, jend1, &
                                 istart2, iend2, jstart2, jend2)
     endif

     ntile_atm_global = ntile_atm
     ncontact_global = ncontact
     if( atm_nest_npes > 0 ) then
        if(ntile_atm .NE. 7) call mpp_error(FATAL,  &
           'xgrid_test: ntile_atm should be 7 when atmos_nest_npes > 0')
        if(ncontact .NE. 13 ) call mpp_error(FATAL,  &
           'xgrid_test: ncontact_atm should be 13 when atmos_nest_npes > 0')
        ntile_atm_global = ntile_atm_global - 1
        ncontact_global  = ncontact_global  - 1
     endif

     if(atm_global_pe) then         
        call mpp_set_current_pelist(atm_global_pelist)
        if(mod(atm_global_npes, ntile_atm_global) .NE. 0 ) call mpp_error(FATAL, &
             "atm_npes_global should be divided by ntile_atm_global")

        allocate(pe_start(ntile_atm_global), pe_end(ntile_atm_global) )
        allocate(global_indices(4, ntile_atm_global), layout(2,ntile_atm_global))
        npes_per_tile = atm_global_npes/ntile_atm_global
        do n = 1, ntile_atm_global
           pe_start(n) = atm_root_pe + (n-1)*npes_per_tile
           pe_end(n)   = atm_root_pe + n*npes_per_tile - 1
           global_indices(:,n) = (/1, atm_nx(n), 1, atm_ny(n)/)
           if(atm_layout(1)*atm_layout(2) == npes_per_tile ) then
              layout(:,n) = atm_layout(:)
           else
              call mpp_define_layout( global_indices(:,n), pe_end(n)-pe_start(n)+1, layout(:,n))
           endif
        end do
 
        allocate(tile_id(ntile_atm_global))
        do n = 1, ntile_atm_global
           tile_id(n) = n
        enddo

        call mpp_define_mosaic(global_indices, layout, Atm_domain, ntile_atm_global, ncontact_global, &
                               tile1(1:ncontact_global), tile2(1:ncontact_global),                    &
                               istart1(1:ncontact_global), iend1(1:ncontact_global),                  &
                               jstart1(1:ncontact_global), jend1(1:ncontact_global),                  &
                               istart2(1:ncontact_global), iend2(1:ncontact_global),                  &
                               jstart2(1:ncontact_global), jend2(1:ncontact_global),                  &
                               pe_start, pe_end, whalo=1, ehalo=1, shalo=1, nhalo=1,                  &
                               tile_id=tile_id, name="atmosphere")
        deallocate( pe_start, pe_end, global_indices, layout, tile_id )
     endif
     if( atm_nest_pe ) then
        atm_nest_nx = atm_nx(ntile_atm)
        atm_nest_ny = atm_ny(ntile_atm)
        call mpp_set_current_pelist(atm_nest_pelist)
        ntile_atm_nest = 1
        ncontact = 0
        allocate(pe_start(1), pe_end(1) )
        allocate(global_indices(4, 1), layout(2,1))
        pe_start(1) = atm_nest_root_pe 
        pe_end(1)   = atm_nest_root_pe + atm_nest_npes - 1
        global_indices(:,1) = (/1, atm_nest_nx, 1, atm_nest_ny/)  !-- the last tile is the nested tile
        if(atm_nest_layout(1)*atm_nest_layout(2) == atm_nest_npes ) then
           layout(:,1) = atm_nest_layout(:)
        else
           call mpp_define_layout( global_indices(:,1), atm_nest_npes, layout(:,1))
        endif
        allocate(tile_id(1))
        tile_id(1) = ntile_atm
        call mpp_define_mosaic(global_indices, layout, Atm_domain, ntile_atm_nest, ncontact, dummy, dummy, &
                               dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end, &
                               whalo=1, ehalo=1, shalo=1, nhalo=1, tile_id=tile_id, name="atmos nest")
        deallocate( pe_start, pe_end, global_indices, layout, tile_id )
     endif

     if( lnd_pe ) then
        call mpp_set_current_pelist(lnd_pelist)
        ncontact = 0 ! no update is needed for land model.
        allocate(pe_start(ntile_lnd), pe_end(ntile_lnd) )
        allocate(global_indices(4,ntile_lnd), layout(2,ntile_lnd))
        if(mod(lnd_npes, ntile_lnd) .NE. 0 ) call mpp_error(FATAL,"lnd_npes should be divided by ntile_lnd")
        npes_per_tile = lnd_npes/ntile_lnd  
        do n = 1, ntile_lnd
           pe_start(n) = lnd_root_pe + (n-1)*npes_per_tile
           pe_end(n)   = lnd_root_pe + n*npes_per_tile - 1
           global_indices(:,n) = (/1, lnd_nx(n), 1, lnd_ny(n)/)
           if(lnd_layout(1)*lnd_layout(2) == npes_per_tile ) then
              layout(:,n) = lnd_layout(:)
           else
              call mpp_define_layout( global_indices(:,n), npes_per_tile, layout(:,n))
           endif
        end do

        ncontact = 0 ! no update is needed for land and ocean model.

        allocate(tile_id(ntile_lnd))
        do n = 1, ntile_lnd
           tile_id(n) = n      
        enddo

        call mpp_define_mosaic(global_indices, layout, Lnd_domain, ntile_lnd, ncontact, dummy, dummy, &
             dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end, tile_id=tile_id, name="land")
        deallocate( pe_start, pe_end, global_indices, layout, tile_id )
     endif

     if( ice_pe ) then
        call mpp_set_current_pelist(ice_pelist)
        ncontact = 0 ! no update is needed for ocn model.
        allocate(pe_start(ntile_ice), pe_end(ntile_ice) )
        allocate(global_indices(4, ntile_ice), layout(2, ntile_ice))
        npes_per_tile = ice_npes/ntile_ice
        do n = 1, ntile_ice
           pe_start(n) = ice_root_pe + (n-1)*npes_per_tile
           pe_end(n)   = ice_root_pe + n*npes_per_tile - 1
           global_indices(:,n) = (/1, ice_nx(n), 1, ice_ny(n)/)
           if(ice_layout(1)*ice_layout(2) == npes_per_tile ) then
              layout(:,n) = ice_layout(:)
           else
              call mpp_define_layout( global_indices(:,n), ice_npes, layout(:,n))
           endif
        end do
        allocate(tile_id(ntile_ice))
        do n = 1, ntile_ice
           tile_id(n) = n
        enddo

        call mpp_define_mosaic(global_indices, layout, Ice_domain, ntile_ice, ncontact, dummy, dummy, &
             dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end, tile_id=tile_id, name="Ice")
        deallocate( pe_start, pe_end, global_indices, layout, tile_id )
     endif
  else
     call mpp_error(FATAL, 'test_xgrid:both AREA_ATM and atm_mosaic does not exist in '//trim(grid_file))
  end if

  deallocate(atm_nx, atm_ny, lnd_nx, lnd_ny, ice_nx, ice_ny)

  if( atm_pe ) then
     call mpp_get_compute_domain(atm_domain, isc_atm, iec_atm, jsc_atm, jec_atm)
     call mpp_get_data_domain(atm_domain, isd_atm, ied_atm, jsd_atm, jed_atm)
     call mpp_get_global_domain(atm_domain, xsize = nxa, ysize = nya)
     nxc_atm = iec_atm - isc_atm + 1
     nyc_atm = jec_atm - jsc_atm + 1
  endif
  if( lnd_pe ) then
     call mpp_get_compute_domain(lnd_domain, isc_lnd, iec_lnd, jsc_lnd, jec_lnd)
     call mpp_get_global_domain(lnd_domain, xsize = nxl, ysize = nyl)
  endif
  if( ice_pe ) then
     call mpp_get_compute_domain(Ice_domain, isc_ice, iec_ice, jsc_ice, jec_ice)
  endif

  ! set up atm_grid for second order conservative interpolation and atm grid is cubic grid.
  if(remap_method == SECOND_ORDER ) then  
      if( atm_nest_npes > 0 ) call mpp_error(FATAL, &
            'test_xgrid: remap_method could not be SECOND_ORDER when atmos_nest_npes > 0')
      if(ntile_atm == 6) then ! 6 tile for cubic grid
        allocate(xt         (isd_atm:ied_atm,jsd_atm:jed_atm),     yt         (isd_atm:ied_atm,jsd_atm:jed_atm)   )
        allocate(xc         (isc_atm:ied_atm,jsc_atm:jed_atm),     yc         (isc_atm:ied_atm,jsc_atm:jed_atm)   )
        allocate(atm_grid%dx(isc_atm:iec_atm,jsc_atm:jed_atm),     atm_grid%dy(isc_atm:iec_atm+1,jsc_atm:jec_atm) )
        allocate(atm_grid%edge_w(jsc_atm:jed_atm),               atm_grid%edge_e(jsc_atm:jed_atm))
        allocate(atm_grid%edge_s(isc_atm:ied_atm),               atm_grid%edge_n(isc_atm:ied_atm))
        allocate(atm_grid%en1 (3,isc_atm:iec_atm,jsc_atm:jed_atm), atm_grid%en2 (3,isc_atm:ied_atm,jsc_atm:jec_atm) )
        allocate(atm_grid%vlon(3,isc_atm:iec_atm,jsc_atm:jec_atm), atm_grid%vlat(3,isc_atm:iec_atm,jsc_atm:jec_atm) )
        allocate(atm_grid%area(isc_atm:iec_atm,jsc_atm:jec_atm) )

        ! first get grid from grid file 
        call get_mosaic_tile_grid(tile_file, atm_mosaic_file, atm_domain)
        allocate(tmpx(nxa*2+1, nya*2+1), tmpy(nxa*2+1, nya*2+1))
        call read_data( tile_file, 'x', tmpx, no_domain=.true.)
        call read_data( tile_file, 'y', tmpy, no_domain=.true.) 
        xt = 0; yt = 0;
        do j = jsc_atm, jec_atm
           do i = isc_atm, iec_atm
              xt(i,j) = tmpx(2*i, 2*j)*DEG_TO_RAD
              yt(i,j) = tmpy(2*i, 2*j)*DEG_TO_RAD
           end do
        end do
        do j = jsc_atm, jed_atm
           do i = isc_atm, ied_atm
              xc(i,j) = tmpx(2*i-1, 2*j-1)*DEG_TO_RAD
              yc(i,j) = tmpy(2*i-1, 2*j-1)*DEG_TO_RAD
           end do
        end do        
        call mpp_update_domains(xt, atm_domain)
        call mpp_update_domains(yt, atm_domain)

        call calc_cubic_grid_info(xt, yt, xc, yc, atm_grid%dx, atm_grid%dy, atm_grid%area, atm_grid%edge_w, &
                                  atm_grid%edge_e, atm_grid%edge_s, atm_grid%edge_n, atm_grid%en1,          &
                                  atm_grid%en2, atm_grid%vlon, atm_grid%vlat, isc_atm==1, iec_atm==nxa,     &
                                  jsc_atm==1, jec_atm==nya                               )
     end if
  end if

  call mpp_set_current_pelist()

  !--- conservation check is done in setup_xmap. 
  call setup_xmap(Xmap, (/ 'ATM', 'OCN', 'LND' /), (/ Atm_domain, Ice_domain, Lnd_domain /), grid_file, atm_grid)
  call setup_xmap(Xmap_runoff, (/ 'LND', 'OCN'/), (/ Lnd_domain, Ice_domain/), grid_file )
  !--- set frac area if nk_lnd or nk_ocn is greater than 1.
  if(nk_lnd > 0 .AND. lnd_pe) then
    allocate(lnd_frac(isc_lnd:iec_lnd, jsc_lnd:jec_lnd, nk_lnd))
    call random_number(lnd_frac)
    lnd_frac = lnd_frac + 0.5
    do j = jsc_lnd, jec_lnd
       do i = isc_lnd, iec_lnd
          tot = sum(lnd_frac(i,j,:))
          do k = 1, nk_lnd
             lnd_frac(i,j,k)=lnd_frac(i,j,k)/tot
          enddo
       enddo
    enddo
    call set_frac_area(lnd_frac, 'LND', xmap)
  endif
  if(nk_ice > 0 ) then
    if( ice_pe ) then
       allocate(ice_frac(isc_ice:iec_ice, jsc_ice:jec_ice, nk_ice))
       call random_number(ice_frac)
       ice_frac = ice_frac + 0.5
       do j = jsc_ice, jec_ice
          do i = isc_ice, iec_ice
             tot = sum(ice_frac(i,j,:))
             do k = 1, nk_ice
                ice_frac(i,j,k)=ice_frac(i,j,k)/tot
             enddo
          enddo
       enddo
    endif
    call set_frac_area(ice_frac, 'OCN', xmap)
  endif

  !--- remap realistic data and write the output file when atmos_input_file does exist
  atm_input_file_exist = file_exist(atm_input_file, domain=atm_domain)
  if( atm_input_file_exist ) then
     if(trim(atm_input_file) == trim(atm_output_file) ) call mpp_error(FATAL, &
          "test_xgrid: atm_input_file should have a different name from atm_output_file")
     call field_size(atm_input_file, atm_field_name, siz, domain=Atm_domain )
     if(siz(1) .NE. nxa .OR. siz(2) .NE. nya ) call mpp_error(FATAL,"test_xgrid: x- and y-size of field "//trim(atm_field_name) &
            //" in file "//trim(atm_input_file) //" does not compabile with the grid size" )
     if(siz(3) > 1) call mpp_error(FATAL,"test_xgrid: number of vertical level of field "//trim(atm_field_name) &
            //" in file "//trim(atm_input_file) //" should be no larger than 1")

     allocate(atm_data_in (isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(atm_data_out(isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(lnd_data_out(isc_lnd:iec_lnd, jsc_lnd:jec_lnd, nk_lnd) )
     allocate(ice_data_out(isc_ice:iec_ice, jsc_ice:jec_ice, nk_ice) )
     allocate(atm_data_out_1(isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(atm_data_out_2(isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(atm_data_out_3(isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     nxgrid = max(xgrid_count(Xmap), 1)
     allocate(x_1(nxgrid), x_2(nxgrid))
     allocate(x_3(nxgrid), x_4(nxgrid))
     x_1 = 0
     x_2 = 0
     x_3 = 0
     x_4 = 0

     atm_data_in  = 0
     atm_data_out = 0
     lnd_data_out = 0
     ice_data_out = 0
     atm_data_out_1 = 0
     atm_data_out_2 = 0
     atm_data_out_3 = 0
     ! test one time level should be sufficient
     call read_data(atm_input_file, atm_field_name, atm_data_in, atm_domain)
     call put_to_xgrid(atm_data_in, 'ATM', x_1, Xmap, remap_method=remap_method)
     call put_to_xgrid(atm_data_in, 'ATM', x_2, Xmap, remap_method=remap_method, complete=.false.)
     call put_to_xgrid(atm_data_in, 'ATM', x_3, Xmap, remap_method=remap_method, complete=.false.)
     call put_to_xgrid(atm_data_in, 'ATM', x_4, Xmap, remap_method=remap_method, complete=.true.)
     min_x = minval(x_1)
     max_x = maxval(x_1)
     !--- check make sure x_2, x_3 and x_4 are the same as x_1
     if(ANY(x_1 .NE. x_2)) call mpp_error(FATAL,"test_xgrid: x_1 and x_2 are not equal")
     if(ANY(x_1 .NE. x_3)) call mpp_error(FATAL,"test_xgrid: x_1 and x_3 are not equal")
     if(ANY(x_1 .NE. x_4)) call mpp_error(FATAL,"test_xgrid: x_1 and x_4 are not equal")
 
     deallocate(x_3,x_4)
     x_2 = 0
     call get_from_xgrid(lnd_data_out, 'LND', x_1, xmap)
     call get_from_xgrid(ice_data_out, 'OCN', x_1, xmap)
     call put_to_xgrid(lnd_data_out, 'LND', x_2, xmap)
     call put_to_xgrid(ice_data_out, 'OCN', x_2, xmap)
     call get_from_xgrid(atm_data_out, 'ATM', x_2, xmap)
     call get_from_xgrid(atm_data_out_1, 'ATM', x_2, xmap, complete=.false.)
     call get_from_xgrid(atm_data_out_2, 'ATM', x_2, xmap, complete=.false.)
     call get_from_xgrid(atm_data_out_3, 'ATM', x_2, xmap, complete=.true.)
     if(ANY(atm_data_out .NE. atm_data_out_1)) &
        call mpp_error(FATAL,"test_xgrid: atm_data_out and atm_data_out_1 are not equal")
     if(ANY(atm_data_out .NE. atm_data_out_2)) &
        call mpp_error(FATAL,"test_xgrid: atm_data_out and atm_data_out_2 are not equal")
     if(ANY(atm_data_out .NE. atm_data_out_3)) &
        call mpp_error(FATAL,"test_xgrid: atm_data_out and atm_data_out_3 are not equal")

     call write_data( atm_output_file, atm_field_name, atm_data_out, atm_domain)
     call write_data( lnd_output_file, atm_field_name, lnd_data_out, lnd_domain)
     call write_data( ice_output_file, atm_field_name, ice_data_out, Ice_domain)
     !--- print out checksum
     write(out_unit,*) "chksum for atm_data_in",  mpp_chksum(atm_data_in)
     write(out_unit,*) "chksum for lnd_data_out", mpp_chksum(lnd_data_out)
     write(out_unit,*) "chksum for ice_data_out", mpp_chksum(ice_data_out)
     write(out_unit,*) "chksum for atm_data_out", mpp_chksum(atm_data_out)

     ! conservation check 
     allocate(atm_area(isc_atm:iec_atm, jsc_atm:jec_atm ) )
     allocate(lnd_area(isc_lnd:iec_lnd, jsc_lnd:jec_lnd ) )
     allocate(ice_area(isc_ice:iec_ice, jsc_ice:jec_ice ) )
     call get_xmap_grid_area("ATM", Xmap, atm_area)
     call get_xmap_grid_area("LND", Xmap, lnd_area)
     call get_xmap_grid_area("OCN", Xmap, ice_area)

     min_atm_in   = minval(atm_data_in)
     max_atm_in   = maxval(atm_data_in)
     min_atm_out  = minval(atm_data_out)
     max_atm_out  = maxval(atm_data_out)
     call mpp_min(min_atm_in)
     call mpp_max(max_atm_in)
     call mpp_min(min_atm_out)
     call mpp_max(max_atm_out)
     call mpp_min(min_x)
     call mpp_max(max_x)

     sum_atm_in  = mpp_global_sum(atm_domain, atm_area * atm_data_in)
     sum_lnd_out = 0
     do k = 1, nk_lnd
        sum_lnd_out = sum_lnd_out + mpp_global_sum(lnd_domain, lnd_area * lnd_data_out(:,:,k))
     enddo
     sum_ice_out = 0
     do k = 1, nk_ice
        sum_ice_out = sum_ice_out + mpp_global_sum(ice_domain, ice_area * ice_data_out(:,:,k))
     enddo
     sum_atm_out = mpp_global_sum(atm_domain, atm_area * atm_data_out)
     write(out_unit,*) "********************** check conservation *********************** "
     write(out_unit,*) "the global area sum of atmos input data is                    : ", sum_atm_in 
     write(out_unit,*) "the global area sum of atmos output data is                   : ", sum_atm_out
     write(out_unit,*) "the global area sum of land output data + ocean output data is: ", sum_lnd_out+sum_ice_out
     write(out_unit,*) "The min of atmos input   data is ", min_atm_in
     write(out_unit,*) "The min of xgrid         data is ", min_x
     write(out_unit,*) "The min of atmos output  data is ", min_atm_out
     write(out_unit,*) "The max of atmos input   data is ", max_atm_in
     write(out_unit,*) "The max of xgrid         data is ", max_x
     write(out_unit,*) "The max of atmos output  data is ", max_atm_out


     deallocate(atm_area, lnd_area, ice_area, atm_data_in, atm_data_out, lnd_data_out, ice_data_out)
     deallocate(atm_data_out_1, atm_data_out_2, atm_data_out_3)
     deallocate(x_1, x_2)
  else
     write(out_unit,*) "NOTE from test_xgrid ==> file "//trim(atm_input_file)//" does not exist, no check is done for real data sets."
  end if           

  runoff_input_file_exist = file_exist(runoff_input_file, domain=atm_domain)     
  if( runoff_input_file_exist ) then
     if( atm_nest_npes > 0 ) call mpp_error(FATAL, &
          "test_xgrid: runoff_input_file_exist should be false when atmos_nest_npes > 0")
     if(trim(runoff_input_file) == trim(runoff_output_file) ) call mpp_error(FATAL, &
          "test_xgrid: runoff_input_file should have a different name from runoff_output_file")
     call field_size(runoff_input_file, runoff_field_name, siz )
     if(siz(1) .NE. nxl .OR. siz(2) .NE. nyl ) call mpp_error(FATAL,"test_xgrid: x- and y-size of field "//trim(runoff_field_name) &
            //" in file "//trim(runoff_input_file) //" does not compabile with the grid size" )
     if(siz(3) > 1) call mpp_error(FATAL,"test_xgrid: number of vertical level of field "//trim(runoff_field_name) &
            //" in file "//trim(runoff_input_file) //" should be no larger than 1")

     allocate(runoff_data_in (isc_lnd:iec_lnd, jsc_lnd:jec_lnd   ) ) 
     allocate(runoff_data_out(isc_ice:iec_ice, jsc_ice:jec_ice, 1) )
     nxgrid = max(xgrid_count(Xmap_runoff), 1)
     allocate(x_1(nxgrid), x_2(nxgrid))

     runoff_data_in  = 0
     runoff_data_out = 0
     ! test one time level should be sufficient
     call read_data(runoff_input_file, runoff_field_name, runoff_data_in, lnd_domain)
     call put_to_xgrid(runoff_data_in, 'LND', x_1, Xmap_runoff)
     call get_from_xgrid(runoff_data_out, 'OCN', x_1, xmap_runoff)
     call write_data( runoff_output_file, runoff_field_name, runoff_data_out, ice_domain)
     ! conservation check 
     allocate(lnd_area(isc_lnd:iec_lnd, jsc_lnd:jec_lnd ) )
     allocate(ice_area(isc_ice:iec_ice, jsc_ice:jec_ice ) )
     call get_xmap_grid_area("LND", Xmap_runoff, lnd_area)
     call get_xmap_grid_area("OCN", Xmap_runoff, ice_area)

     sum_runoff_in  = mpp_global_sum(lnd_domain, lnd_area * runoff_data_in)
     sum_runoff_out = mpp_global_sum(Ice_domain, ice_area * runoff_data_out(:,:,1))
     write(out_unit,*) "********************** check conservation *********************** "
     write(out_unit,*) "the global area sum of runoff input data is                    : ", sum_runoff_in 
     write(out_unit,*) "the global area sum of runoff output data is                   : ", sum_runoff_out
  else
     write(out_unit,*) "NOTE from test_xgrid ==> file "//trim(runoff_input_file)//" does not exist, no check is done for real data sets."
  end if           

  ! when num_iter is greater than 0, create random number as input to test the performance of xgrid_mod.
  if(num_iter > 0) then
    if( atm_nest_npes > 0 ) call mpp_error(FATAL, &
          "test_xgrid: num_iter > 0 when atm_nest_npes > 0")

     allocate(atm_data_in (isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(atm_data_out(isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(lnd_data_out(isc_lnd:iec_lnd, jsc_lnd:jec_lnd, nk_lnd) )
     allocate(ice_data_out(isc_ice:iec_ice, jsc_ice:jec_ice, nk_ice) )
     nxgrid = max(xgrid_count(Xmap), 1)
     allocate(x_1(nxgrid), x_2(nxgrid))  
     atm_data_in  = 0
     atm_data_out = 0
     lnd_data_out = 0
     ice_data_out = 0   
     allocate(atm_area(isc_atm:iec_atm, jsc_atm:jec_atm ) )
     allocate(lnd_area(isc_lnd:iec_lnd, jsc_lnd:jec_lnd ) )
     allocate(ice_area(isc_ice:iec_ice, jsc_ice:jec_ice ) )
     call get_xmap_grid_area("ATM", Xmap, atm_area)
     call get_xmap_grid_area("LND", Xmap, lnd_area)
     call get_xmap_grid_area("OCN", Xmap, ice_area)
     do n = 1, num_iter
        call random_number(atm_data_in)
        call put_to_xgrid(atm_data_in, 'ATM', x_1, Xmap, remap_method=remap_method)

        call get_from_xgrid(lnd_data_out, 'LND', x_1, xmap)
        call get_from_xgrid(ice_data_out, 'OCN', x_1, xmap)

        call put_to_xgrid(lnd_data_out, 'LND', x_2, xmap)
        call put_to_xgrid(ice_data_out, 'OCN', x_2, xmap)

        call get_from_xgrid(atm_data_out, 'ATM', x_2, xmap)
        sum_atm_in  = mpp_global_sum(atm_domain, atm_area * atm_data_in)
        sum_lnd_out = 0
        do k = 1, nk_lnd
           sum_lnd_out = sum_lnd_out + mpp_global_sum(lnd_domain, lnd_area * lnd_data_out(:,:,k))
        enddo
        sum_ice_out = 0
        do k = 1, nk_ice
           sum_ice_out = sum_ice_out + mpp_global_sum(Ice_domain, ice_area * ice_data_out(:,:,k))
        enddo
        sum_atm_out = mpp_global_sum(atm_domain, atm_area * atm_data_out)
        write(out_unit,*) "********************** check conservation *********************** "
        write(out_unit,*) "the global area sum of atmos input data is                    : ", sum_atm_in 
        write(out_unit,*) "the global area sum of atmos output data is                   : ", sum_atm_out
        write(out_unit,*) "the global area sum of land output data + ocean output data is: ", sum_lnd_out+sum_ice_out
     enddo
     deallocate(atm_area, lnd_area, ice_area, atm_data_in, atm_data_out, lnd_data_out, ice_data_out)
     deallocate(x_1, x_2)
  endif  

  write(out_unit,*) "************************************************************************"
  write(out_unit,*) "***********      Finish running program test_xgrid         *************"
  write(out_unit,*) "************************************************************************"

  call fms_io_exit
  call fms_end

end program xgrid_test

#else
module null_test_xgrid
end module  

#endif /* test_mpp */
