#ifdef TEST_XGRID
! Now only test some simple test, will test cubic grid mosaic in the future.

program xgrid_test

  use mpp_mod,         only : mpp_pe, mpp_npes, mpp_error, FATAL, mpp_chksum, mpp_min, mpp_max
  use mpp_mod,         only : mpp_clock_begin, mpp_clock_end, mpp_clock_id, MPP_CLOCK_SYNC
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

implicit none

  real, parameter :: EPSLN = 1.0e-10
  character(len=256) :: atm_input_file  = "INPUT/atmos_input.nc"
  character(len=256) :: atm_output_file = "atmos_output.nc"
  character(len=256) :: lnd_output_file = "land_output.nc"
  character(len=256) :: ocn_output_file = "ocean_output.nc"
  character(len=256) :: atm_field_name  = "none"

  character(len=256) :: runoff_input_file  = "INPUT/land_runoff.nc"
  character(len=256) :: runoff_output_file  = "land_runoff.nc"
  character(len=256) :: runoff_field_name  = "none"
  integer            :: num_iter           = 0 
  integer            :: nk_lnd = 1, nk_ocn = 1
  integer            :: atm_layout(2) = (/0,0/)
  integer            :: lnd_layout(2) = (/0,0/)
  integer            :: ocn_layout(2) = (/0,0/)


  namelist /xgrid_test_nml/ atm_input_file, atm_field_name, runoff_input_file, runoff_field_name, num_iter, &
                            nk_lnd, nk_ocn, atm_layout, ocn_layout, lnd_layout

  integer              :: remap_method
  integer              :: pe, npes, ierr, nml_unit, io, n
  integer              :: siz(4), ntile_lnd, ntile_atm, ntile_ocn, ncontact
  integer, allocatable :: layout(:,:), global_indices(:,:)
  integer, allocatable :: atm_nx(:), atm_ny(:), ocn_nx(:), ocn_ny(:), lnd_nx(:), lnd_ny(:)
  integer, allocatable :: pe_start(:), pe_end(:), dummy(:)
  integer, allocatable :: istart1(:), iend1(:), jstart1(:), jend1(:), tile1(:)
  integer, allocatable :: istart2(:), iend2(:), jstart2(:), jend2(:), tile2(:)
  character(len=256)   :: grid_file = "INPUT/grid_spec.nc"
  character(len=256)   :: atm_mosaic, ocn_mosaic, lnd_mosaic
  character(len=256)   :: atm_mosaic_file, ocn_mosaic_file, lnd_mosaic_file
  character(len=256)   :: grid_descriptor, tile_file
  integer              :: isc_atm, iec_atm, jsc_atm, jec_atm, nxc_atm, nyc_atm
  integer              :: isc_lnd, iec_lnd, jsc_lnd, jec_lnd
  integer              :: isc_ocn, iec_ocn, jsc_ocn, jec_ocn
  integer              :: isd_atm, ied_atm, jsd_atm, jed_atm
  integer              :: unit, i, j, nxa, nya, nxgrid, nxl, nyl, out_unit, k
  type(domain2d)       :: Atm_domain, Ocn_domain, Lnd_domain
  type(xmap_type)      :: Xmap, Xmap_runoff
  type(grid_box_type)  :: atm_grid
  real, allocatable    :: xt(:,:), yt(:,:)  ! on T-cell data domain
  real, allocatable    :: xc(:,:), yc(:,:)  ! on C-cell compute domain
  real, allocatable    :: tmpx(:,:), tmpy(:,:)
  real, allocatable    :: atm_data_in(:,:), atm_data_out(:,:)
  real, allocatable    :: atm_data_out_1(:,:), atm_data_out_2(:,:), atm_data_out_3(:,:)
  real, allocatable    :: lnd_data_out(:,:,:), ocn_data_out(:,:,:)
  real, allocatable    :: runoff_data_in(:,:), runoff_data_out(:,:,:)
  real, allocatable    :: atm_area(:,:), lnd_area(:,:), ocn_area(:,:)
  real, allocatable    :: lnd_frac(:,:,:), ocn_frac(:,:,:)
  real, allocatable    :: x_1(:), x_2(:), x_3(:), x_4(:)
  real                 :: sum_atm_in, sum_ocn_out, sum_lnd_out, sum_atm_out
  real                 :: sum_runoff_in, sum_runoff_out, tot
  real                 :: min_atm_in, max_atm_in, min_atm_out, max_atm_out
  real                 :: min_x, max_x
  logical              :: atm_input_file_exist, runoff_input_file_exist
  integer              :: npes_per_tile
  integer              :: id_put_side1_to_xgrid, id_get_side1_from_xgrid
  integer              :: id_put_side2_to_xgrid, id_get_side2_from_xgrid
  integer              :: id_setup_xmap, id_other, id_frac_area

  call fms_init
  id_setup_xmap = mpp_clock_id("setup_xmap", flags=MPP_CLOCK_SYNC)
  id_put_side1_to_xgrid = mpp_clock_id("put_side1_to_xgrid", flags=MPP_CLOCK_SYNC)
  id_get_side1_from_xgrid = mpp_clock_id("get_side1_from_xgrid", flags=MPP_CLOCK_SYNC)
  id_put_side2_to_xgrid = mpp_clock_id("put_side2_to_xgrid", flags=MPP_CLOCK_SYNC)
  id_get_side2_from_xgrid = mpp_clock_id("get_side2_from_xgrid", flags=MPP_CLOCK_SYNC)
  id_other = mpp_clock_id("other", flags=MPP_CLOCK_SYNC)
  id_frac_area = mpp_clock_id("set_frac_area", flags=MPP_CLOCK_SYNC)

  call mpp_clock_begin(id_other)
  call mpp_domains_init
  call xgrid_init(remap_method)

  npes     = mpp_npes()
  pe       = mpp_pe()
  out_unit = stdout()

 if (file_exist('input.nml')) then
   ierr=1
   nml_unit = open_namelist_file()
   do while (ierr /= 0)
     read(nml_unit, nml=xgrid_test_nml, iostat=io, end=10)
          ierr = check_nml_error(io, 'xgrid_test_nml')
   enddo
10 call close_file(nml_unit)
 endif

  ntile_atm = 1
  ntile_ocn = 1
  ntile_lnd = 1

  if(field_exist(grid_file, "AREA_ATM" ) ) then
     allocate(atm_nx(1), atm_ny(1))
     allocate(lnd_nx(1), lnd_ny(1))
     allocate(ocn_nx(1), ocn_ny(1))
     call field_size(grid_file, "AREA_ATM", siz )
     atm_nx = siz(1); atm_ny = siz(2)
     call field_size(grid_file, "AREA_OCN", siz )
     ocn_nx = siz(1); ocn_ny = siz(2)
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
     if( ocn_layout(1)*ocn_layout(2) .NE. npes ) then
        call mpp_define_layout( (/1,ocn_nx,1,ocn_ny/), npes, ocn_layout)
     endif
     call mpp_define_domains( (/1,ocn_nx,1,ocn_ny/), ocn_layout, Ocn_domain, name="ocean")
  else if (field_exist(grid_file, "atm_mosaic" ) ) then
     !--- Get the mosaic data of each component model 
     call read_data(grid_file, 'atm_mosaic', atm_mosaic)
     call read_data(grid_file, 'lnd_mosaic', lnd_mosaic)
     call read_data(grid_file, 'ocn_mosaic', ocn_mosaic)
     atm_mosaic_file = 'INPUT/'//trim(atm_mosaic)//'.nc'
     lnd_mosaic_file = 'INPUT/'//trim(lnd_mosaic)//'.nc'
     ocn_mosaic_file = 'INPUT/'//trim(ocn_mosaic)//'.nc'

     ntile_lnd = get_mosaic_ntiles(lnd_mosaic_file);
     ntile_ocn = get_mosaic_ntiles(ocn_mosaic_file);
     ntile_atm = get_mosaic_ntiles(atm_mosaic_file);
     if(ntile_ocn > 1) call mpp_error(FATAL,  &
           'xgrid_test: there is more than one tile in ocn_mosaic, which is not implemented yet')

     write(out_unit,*)" There is ", ntile_atm, " tiles in atmos mosaic"
     write(out_unit,*)" There is ", ntile_lnd, " tiles in land  mosaic"
     write(out_unit,*)" There is ", ntile_ocn, " tiles in ocean mosaic"
     allocate(atm_nx(ntile_atm), atm_ny(ntile_atm))
     allocate(lnd_nx(ntile_lnd), lnd_ny(ntile_lnd))
     allocate(ocn_nx(ntile_ocn), ocn_ny(ntile_ocn))

     call get_mosaic_grid_sizes(atm_mosaic_file, atm_nx, atm_ny)
     call get_mosaic_grid_sizes(lnd_mosaic_file, lnd_nx, lnd_ny)
     call get_mosaic_grid_sizes(ocn_mosaic_file, ocn_nx, ocn_ny)

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

     if(mod(npes, ntile_atm) .NE. 0 ) call mpp_error(FATAL,"npes should be divided by ntile_atm")

     allocate(pe_start(ntile_atm), pe_end(ntile_atm) )
     allocate(global_indices(4, ntile_atm), layout(2,ntile_atm))
     npes_per_tile = npes/ntile_atm
     do n = 1, ntile_atm
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile - 1
        global_indices(:,n) = (/1, atm_nx(n), 1, atm_ny(n)/)
        if(atm_layout(1)*atm_layout(2) == npes_per_tile ) then
           layout(:,n) = atm_layout(:)
        else
           call mpp_define_layout( global_indices(:,n), pe_end(n)-pe_start(n)+1, layout(:,n))
        endif
     end do
 
     call mpp_define_mosaic(global_indices, layout, Atm_domain, ntile_atm, ncontact, tile1, tile2, &
                            istart1, iend1, jstart1, jend1, istart2, iend2, jstart2, jend2,        &
                            pe_start, pe_end, whalo=1, ehalo=1, shalo=1, nhalo=1, name="atmosphere")
     deallocate( pe_start, pe_end, global_indices, layout )

     allocate(pe_start(ntile_lnd), pe_end(ntile_lnd) )
     allocate(global_indices(4,ntile_lnd), layout(2,ntile_lnd))
     npes_per_tile = npes/ntile_lnd  
     do n = 1, ntile_lnd
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile - 1
        global_indices(:,n) = (/1, lnd_nx(n), 1, lnd_ny(n)/)
        if(lnd_layout(1)*lnd_layout(2) == npes_per_tile ) then
           layout(:,n) = lnd_layout(:)
        else
           call mpp_define_layout( global_indices(:,n), pe_end(n)-pe_start(n)+1, layout(:,n))
        endif
     end do

     ncontact = 0 ! no update is needed for land and ocean model.
 
     call mpp_define_mosaic(global_indices, layout, Lnd_domain, ntile_lnd, ncontact, dummy, dummy, &
                            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end, name="land")
     deallocate( pe_start, pe_end, global_indices, layout )

     allocate(pe_start(ntile_ocn), pe_end(ntile_ocn) )
     allocate(global_indices(4, ntile_ocn), layout(2, ntile_ocn))
     npes_per_tile = npes/ntile_ocn
     do n = 1, ntile_ocn
        pe_start(n) = (n-1)*npes_per_tile
        pe_end(n)   = n*npes_per_tile - 1
        global_indices(:,n) = (/1, ocn_nx(n), 1, ocn_ny(n)/)
        if(ocn_layout(1)*ocn_layout(2) == npes_per_tile ) then
           layout(:,n) = ocn_layout(:)
        else
           call mpp_define_layout( global_indices(:,n), pe_end(n)-pe_start(n)+1, layout(:,n))
        endif
     end do
 
     call mpp_define_mosaic(global_indices, layout, Ocn_domain, ntile_ocn, ncontact, dummy, dummy, &
                            dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, pe_start, pe_end, name="ocean")
     deallocate( pe_start, pe_end, global_indices, layout )
  else
     call mpp_error(FATAL, 'test_xgrid:both AREA_ATM and atm_mosaic does not exist in '//trim(grid_file))
  end if

  deallocate(atm_nx, atm_ny, lnd_nx, lnd_ny, ocn_nx, ocn_ny)

  call mpp_get_compute_domain(atm_domain, isc_atm, iec_atm, jsc_atm, jec_atm)
  call mpp_get_compute_domain(lnd_domain, isc_lnd, iec_lnd, jsc_lnd, jec_lnd)
  call mpp_get_compute_domain(ocn_domain, isc_ocn, iec_ocn, jsc_ocn, jec_ocn)
  call mpp_get_data_domain(atm_domain, isd_atm, ied_atm, jsd_atm, jed_atm)
  call mpp_get_global_domain(atm_domain, xsize = nxa, ysize = nya)
  call mpp_get_global_domain(lnd_domain, xsize = nxl, ysize = nyl)
  nxc_atm = iec_atm - isc_atm + 1
  nyc_atm = jec_atm - jsc_atm + 1

  ! set up atm_grid for second order conservative interpolation and atm grid is cubic grid.
  if(remap_method == SECOND_ORDER ) then  
     ! check if atmos mosaic is cubic grid or not */
!     call mpp_open(unit,trim(atm_mosaic_file),MPP_RDONLY,MPP_NETCDF,threading=MPP_MULTI,fileset=MPP_SINGLE)
!     call mpp_get_att_value(unit, "mosaic", "grid_descriptor", grid_descriptor)
!     call mpp_close(unit)


!     if(trim(grid_descriptor) == "cubic_grid") then
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
  call mpp_clock_end(id_other)

  !--- conservation check is done in setup_xmap. 
  call mpp_clock_begin(id_setup_xmap)
  call setup_xmap(Xmap, (/ 'ATM', 'OCN', 'LND' /), (/ Atm_domain, Ocn_domain, Lnd_domain /), grid_file, atm_grid)
  call setup_xmap(Xmap_runoff, (/ 'LND', 'OCN'/), (/ Lnd_domain, Ocn_domain/), grid_file )
  call mpp_clock_end(id_setup_xmap)
  call mpp_clock_begin(id_frac_area)
  !--- set frac area if nk_lnd or nk_ocn is greater than 1.
  if(nk_lnd > 0) then
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
  if(nk_ocn > 0) then
    allocate(ocn_frac(isc_ocn:iec_ocn, jsc_ocn:jec_ocn, nk_ocn))
    call random_number(ocn_frac)
    ocn_frac = ocn_frac + 0.5
    do j = jsc_ocn, jec_ocn
       do i = isc_ocn, iec_ocn
          tot = sum(ocn_frac(i,j,:))
          do k = 1, nk_ocn
             ocn_frac(i,j,k)=ocn_frac(i,j,k)/tot
          enddo
       enddo
    enddo
    call set_frac_area(ocn_frac, 'OCN', xmap)
  endif

  call mpp_clock_end(id_frac_area)
  call mpp_clock_begin(id_other)

  call set_domain(atm_domain)
  !--- remap realistic data and write the output file when atmos_input_file does exist
  atm_input_file_exist = file_exist(atm_input_file, domain=atm_domain)
  if( atm_input_file_exist ) then
     if(trim(atm_input_file) == trim(atm_output_file) ) call mpp_error(FATAL, &
          "test_xgrid: atm_input_file should have a different name from atm_output_file")
     call field_size(atm_input_file, atm_field_name, siz )
     if(siz(1) .NE. nxa .OR. siz(2) .NE. nya ) call mpp_error(FATAL,"test_xgrid: x- and y-size of field "//trim(atm_field_name) &
            //" in file "//trim(atm_input_file) //" does not compabile with the grid size" )
     if(siz(3) > 1) call mpp_error(FATAL,"test_xgrid: number of vertical level of field "//trim(atm_field_name) &
            //" in file "//trim(atm_input_file) //" should be no larger than 1")

     allocate(atm_data_in (isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(atm_data_out(isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(lnd_data_out(isc_lnd:iec_lnd, jsc_lnd:jec_lnd, nk_lnd) )
     allocate(ocn_data_out(isc_ocn:iec_ocn, jsc_ocn:jec_ocn, nk_ocn) )
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
     ocn_data_out = 0
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
     call get_from_xgrid(ocn_data_out, 'OCN', x_1, xmap)
     call put_to_xgrid(lnd_data_out, 'LND', x_2, xmap)
     call put_to_xgrid(ocn_data_out, 'OCN', x_2, xmap)
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
     call write_data( ocn_output_file, atm_field_name, ocn_data_out, ocn_domain)
     !--- print out checksum
     write(out_unit,*) "chksum for atm_data_in",  mpp_chksum(atm_data_in)
     write(out_unit,*) "chksum for lnd_data_out", mpp_chksum(lnd_data_out)
     write(out_unit,*) "chksum for ocn_data_out", mpp_chksum(ocn_data_out)
     write(out_unit,*) "chksum for atm_data_out", mpp_chksum(atm_data_out)

     ! conservation check 
     allocate(atm_area(isc_atm:iec_atm, jsc_atm:jec_atm ) )
     allocate(lnd_area(isc_lnd:iec_lnd, jsc_lnd:jec_lnd ) )
     allocate(ocn_area(isc_ocn:iec_ocn, jsc_ocn:jec_ocn ) )
     call get_xmap_grid_area("ATM", Xmap, atm_area)
     call get_xmap_grid_area("LND", Xmap, lnd_area)
     call get_xmap_grid_area("OCN", Xmap, ocn_area)

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
     sum_ocn_out = 0
     do k = 1, nk_ocn
        sum_ocn_out = sum_ocn_out + mpp_global_sum(ocn_domain, ocn_area * ocn_data_out(:,:,k))
     enddo
     sum_atm_out = mpp_global_sum(atm_domain, atm_area * atm_data_out)
     write(out_unit,*) "********************** check conservation *********************** "
     write(out_unit,*) "the global area sum of atmos input data is                    : ", sum_atm_in 
     write(out_unit,*) "the global area sum of atmos output data is                   : ", sum_atm_out
     write(out_unit,*) "the global area sum of land output data + ocean output data is: ", sum_lnd_out+sum_ocn_out
     write(out_unit,*) "The min of atmos input   data is ", min_atm_in
     write(out_unit,*) "The min of xgrid         data is ", min_x
     write(out_unit,*) "The min of atmos output  data is ", min_atm_out
     write(out_unit,*) "The max of atmos input   data is ", max_atm_in
     write(out_unit,*) "The max of xgrid         data is ", max_x
     write(out_unit,*) "The max of atmos output  data is ", max_atm_out


     deallocate(atm_area, lnd_area, ocn_area, atm_data_in, atm_data_out, lnd_data_out, ocn_data_out)
     deallocate(atm_data_out_1, atm_data_out_2, atm_data_out_3)
     deallocate(x_1, x_2)
  else
     write(out_unit,*) "NOTE from test_xgrid ==> file "//trim(atm_input_file)//" does not exist, no check is done for real data sets."
  end if           

  runoff_input_file_exist = file_exist(runoff_input_file)     
  if( runoff_input_file_exist ) then
     if(trim(runoff_input_file) == trim(runoff_output_file) ) call mpp_error(FATAL, &
          "test_xgrid: runoff_input_file should have a different name from runoff_output_file")
     call field_size(runoff_input_file, runoff_field_name, siz )
     if(siz(1) .NE. nxl .OR. siz(2) .NE. nyl ) call mpp_error(FATAL,"test_xgrid: x- and y-size of field "//trim(runoff_field_name) &
            //" in file "//trim(runoff_input_file) //" does not compabile with the grid size" )
     if(siz(3) > 1) call mpp_error(FATAL,"test_xgrid: number of vertical level of field "//trim(runoff_field_name) &
            //" in file "//trim(runoff_input_file) //" should be no larger than 1")

     allocate(runoff_data_in (isc_lnd:iec_lnd, jsc_lnd:jec_lnd   ) ) 
     allocate(runoff_data_out(isc_ocn:iec_ocn, jsc_ocn:jec_ocn, 1) )
     nxgrid = max(xgrid_count(Xmap_runoff), 1)
     allocate(x_1(nxgrid), x_2(nxgrid))

     runoff_data_in  = 0
     runoff_data_out = 0
     ! test one time level should be sufficient
     call read_data(runoff_input_file, runoff_field_name, runoff_data_in, lnd_domain)
     call put_to_xgrid(runoff_data_in, 'LND', x_1, Xmap_runoff)
     call get_from_xgrid(runoff_data_out, 'OCN', x_1, xmap_runoff)
     call write_data( runoff_output_file, runoff_field_name, runoff_data_out, ocn_domain)
     ! conservation check 
     allocate(lnd_area(isc_lnd:iec_lnd, jsc_lnd:jec_lnd ) )
     allocate(ocn_area(isc_ocn:iec_ocn, jsc_ocn:jec_ocn ) )
     call get_xmap_grid_area("LND", Xmap_runoff, lnd_area)
     call get_xmap_grid_area("OCN", Xmap_runoff, ocn_area)

     sum_runoff_in  = mpp_global_sum(lnd_domain, lnd_area * runoff_data_in)
     sum_runoff_out = mpp_global_sum(ocn_domain, ocn_area * runoff_data_out(:,:,1))
     write(out_unit,*) "********************** check conservation *********************** "
     write(out_unit,*) "the global area sum of runoff input data is                    : ", sum_runoff_in 
     write(out_unit,*) "the global area sum of runoff output data is                   : ", sum_runoff_out
  else
     write(out_unit,*) "NOTE from test_xgrid ==> file "//trim(runoff_input_file)//" does not exist, no check is done for real data sets."
  end if           

  call mpp_clock_end(id_other)

  ! when num_iter is greater than 0, create random number as input to test the performance of xgrid_mod.
  if(num_iter > 0) then
     call mpp_clock_begin(id_other)
     allocate(atm_data_in (isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(atm_data_out(isc_atm:iec_atm, jsc_atm:jec_atm   ) )
     allocate(lnd_data_out(isc_lnd:iec_lnd, jsc_lnd:jec_lnd, nk_lnd) )
     allocate(ocn_data_out(isc_ocn:iec_ocn, jsc_ocn:jec_ocn, nk_ocn) )
     nxgrid = max(xgrid_count(Xmap), 1)
     allocate(x_1(nxgrid), x_2(nxgrid))  
     atm_data_in  = 0
     atm_data_out = 0
     lnd_data_out = 0
     ocn_data_out = 0   
     allocate(atm_area(isc_atm:iec_atm, jsc_atm:jec_atm ) )
     allocate(lnd_area(isc_lnd:iec_lnd, jsc_lnd:jec_lnd ) )
     allocate(ocn_area(isc_ocn:iec_ocn, jsc_ocn:jec_ocn ) )
     call get_xmap_grid_area("ATM", Xmap, atm_area)
     call get_xmap_grid_area("LND", Xmap, lnd_area)
     call get_xmap_grid_area("OCN", Xmap, ocn_area)
     call mpp_clock_end(id_other)
     do n = 1, num_iter
        call mpp_clock_begin(id_other)
        call random_number(atm_data_in)
        call mpp_clock_end(id_other)
        call mpp_clock_begin(id_put_side1_to_xgrid)
        call put_to_xgrid(atm_data_in, 'ATM', x_1, Xmap, remap_method=remap_method)
        call mpp_clock_end(id_put_side1_to_xgrid)

        call mpp_clock_begin(id_get_side2_from_xgrid)
        call get_from_xgrid(lnd_data_out, 'LND', x_1, xmap)
        call get_from_xgrid(ocn_data_out, 'OCN', x_1, xmap)
        call mpp_clock_end(id_get_side2_from_xgrid)

        call mpp_clock_begin(id_put_side2_to_xgrid)
        call put_to_xgrid(lnd_data_out, 'LND', x_2, xmap)
        call put_to_xgrid(ocn_data_out, 'OCN', x_2, xmap)
        call mpp_clock_end(id_put_side2_to_xgrid)

        call mpp_clock_begin(id_get_side1_from_xgrid)
        call get_from_xgrid(atm_data_out, 'ATM', x_2, xmap)
        call mpp_clock_end(id_get_side1_from_xgrid)
        call mpp_clock_begin(id_other)
        sum_atm_in  = mpp_global_sum(atm_domain, atm_area * atm_data_in)
        sum_lnd_out = 0
        do k = 1, nk_lnd
           sum_lnd_out = sum_lnd_out + mpp_global_sum(lnd_domain, lnd_area * lnd_data_out(:,:,k))
        enddo
        sum_ocn_out = 0
        do k = 1, nk_ocn
           sum_ocn_out = sum_ocn_out + mpp_global_sum(ocn_domain, ocn_area * ocn_data_out(:,:,k))
        enddo
        sum_atm_out = mpp_global_sum(atm_domain, atm_area * atm_data_out)
        write(out_unit,*) "********************** check conservation *********************** "
        write(out_unit,*) "the global area sum of atmos input data is                    : ", sum_atm_in 
        write(out_unit,*) "the global area sum of atmos output data is                   : ", sum_atm_out
        write(out_unit,*) "the global area sum of land output data + ocean output data is: ", sum_lnd_out+sum_ocn_out
     call mpp_clock_end(id_other)
     enddo
     call mpp_clock_begin(id_other)
     deallocate(atm_area, lnd_area, ocn_area, atm_data_in, atm_data_out, lnd_data_out, ocn_data_out)
     deallocate(x_1, x_2)
     call mpp_clock_end(id_other)
  endif  

  call mpp_clock_begin(id_other)
  write(out_unit,*) "************************************************************************"
  write(out_unit,*) "***********      Finish running program test_xgrid         *************"
  write(out_unit,*) "************************************************************************"

  call mpp_domains_exit
  call fms_io_exit
  call mpp_clock_end(id_other)
  call fms_end

end program xgrid_test

#else
module null_test_xgrid
end module  

#endif /* test_mpp */
