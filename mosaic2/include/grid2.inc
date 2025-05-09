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

!> @file

!> @brief return grid cell area for the specified model component and tile
subroutine GET_GRID_CELL_AREA_SG_(component, tile, cellarea, domain)
  character(len=*), intent(in)    :: component !< Component model (atm, lnd, ocn)
  integer         , intent(in)    :: tile !< Tile number
  real(kind=FMS_MOS_KIND_) , intent(inout) :: cellarea(:,:) !< Cell area
  type(domain2d)  , intent(in), optional :: domain !< Domain

  ! local vars
  integer :: nlon, nlat
  real(kind=r8_kind), allocatable :: glonb(:,:), glatb(:,:)
  real(kind=r8_kind), allocatable :: cellarea8(:,:)

  allocate(cellarea8(size(cellarea,1),size(cellarea,2)))

  select case(grid_version)
  case(VERSION_GEOLON_T,VERSION_X_T)
     if (.not. grid_spec_exists) then
        call mpp_error(FATAL, 'grid2_mod(get_grid_cell_area_SG): grid_spec does not exist')
     end if
     select case(trim(component))
     case('LND')
        call read_data(gridfileobj, 'AREA_LND_CELL', cellarea8)
     case('ATM','OCN')
        call read_data(gridfileobj, 'AREA_'//trim(uppercase(component)),cellarea8)
     case default
        call mpp_error(FATAL, module_name//'/get_grid_cell_area'//&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
     end select
     ! convert area to m2
     cellarea = real( cellarea8*4.0_r8_kind*PI*RADIUS**2, FMS_MOS_KIND_)
  case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES)
     if (present(domain)) then
        call mpp_get_compute_domain(domain,xsize=nlon,ysize=nlat)
     else
        call get_grid_size(component,tile,nlon,nlat)
     endif
     allocate(glonb(nlon+1,nlat+1),glatb(nlon+1,nlat+1))
     call get_grid_cell_vertices(component, tile, glonb, glatb, domain)
     if (great_circle_algorithm) then
        call calc_mosaic_grid_great_circle_area(glonb*PI/180.0_r8_kind, glatb*PI/180_r8_kind, cellarea8)
        cellarea=real(cellarea8,FMS_MOS_KIND_)
     else
        call calc_mosaic_grid_area(glonb*PI/180.0_r8_kind, glatb*PI/180_r8_kind, cellarea8)
        cellarea=real(cellarea8,FMS_MOS_KIND_)
     end if
     deallocate(glonb,glatb)
  end select

  deallocate(cellarea8)

end subroutine GET_GRID_CELL_AREA_SG_

!> @brief get the area of the component per grid cell
subroutine GET_GRID_COMP_AREA_SG_(component,tile,area,domain)
  character(len=*) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real(kind=FMS_MOS_KIND_), intent(inout) :: area(:,:) !< Area of grid cell
  type(domain2d), intent(in), optional :: domain !< Domain
  ! local vars
  integer :: n_xgrid_files ! number of exchange grid files in the mosaic
  integer :: siz(2), nxgrid
  integer :: i,j,m,n
  integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
  real(kind=r8_kind), allocatable :: xgrid_area(:)
  real(kind=r8_kind), allocatable :: rmask(:,:)
  character(len=MAX_NAME) :: &
     xgrid_name, & ! name of the variable holding xgrid names
     tile_name,  & ! name of the tile
     mosaic_name ! name of the mosaic
  character(len=FMS_PATH_LEN) :: &
     tilefile,   & ! name of current tile file
     xgrid_file  ! name of the current xgrid file
  character(len=4096)     :: attvalue
  character(len=MAX_NAME), allocatable :: nest_tile_name(:)
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0, j0 ! offsets for x and y, respectively
  integer :: num_nest_tile, ntiles
  logical :: is_nest
  integer :: found_xgrid_files ! how many xgrid files we actually found in the grid spec
  integer :: ibegin, iend, bsize, l
  type(FmsNetcdfFile_t) :: tilefileobj, xgrid_fileobj

  real(r8_kind),allocatable :: area8(:,:)

  allocate(area8(size(area,1),size(area,2)))

  select case (grid_version   )
  case(VERSION_GEOLON_T,VERSION_X_T)
     if (.not. grid_spec_exists) then
        call mpp_error(FATAL, 'grid2_mod(get_grid_comp_area_SG): grid_spec does not exist')
     end if
     select case(component)
     case('ATM')
        call read_data(gridfileobj,'AREA_ATM',area8)
     case('OCN')
        allocate(rmask(size(area8,1),size(area8,2)))
        call read_data(gridfileobj,'AREA_OCN',area8)
        call read_data(gridfileobj,'wet',     rmask)
        area = real(area8*rmask, FMS_MOS_KIND_)
        deallocate(rmask)
     case('LND')
        call read_data(gridfileobj,'AREA_LND',area8)
     case default
        call mpp_error(FATAL, module_name//'/get_grid_comp_area'//&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
     end select
  case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES) ! mosaic gridspec
     select case (component)
     case ('ATM')
        ! just read the grid cell area and return
        call get_grid_cell_area(component,tile,area8)
        area = real(area8, FMS_MOS_KIND_)
        return
     case ('LND')
        xgrid_name = 'aXl_file'
        if (.not. grid_spec_exists) then
           call mpp_error(FATAL, 'grid2_mod(get_grid_comp_area_SG): grid_spec does not exist')
        end if
        call read_data(gridfileobj, 'lnd_mosaic', mosaic_name)
        tile_name  = trim(mosaic_name)//'_tile'//char(tile+ichar('0'))
     case ('OCN')
        xgrid_name = 'aXo_file'
        if (.not. grid_spec_exists) then
           call mpp_error(FATAL, 'grid2_mod(get_grid_comp_area_SG): grid_spec does not exist')
        end if
        call read_data(gridfileobj, 'ocn_mosaic', mosaic_name)
        tile_name  = trim(mosaic_name)//'_tile'//char(tile+ichar('0'))
     case default
        call mpp_error(FATAL, module_name//'/get_grid_comp_area'//&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
     end select
     ! get the boundaries of the requested domain
     if(present(domain)) then
        call mpp_get_compute_domain(domain,is,ie,js,je)
        i0 = 1-is ; j0=1-js
     else
        call get_grid_size(component,tile,ie,je)
        is = 1 ; i0 = 0
        js = 1 ; j0 = 0
     endif
     if (size(area8,1)/=ie-is+1.or.size(area8,2)/=je-js+1) &
          call mpp_error(FATAL, module_name//'/get_grid_comp_area '//&
          'size of the output argument "area" is not consistent with the domain')

     ! find the nest tile
     if (.not. grid_spec_exists) then
        call mpp_error(FATAL, 'grid2_mod(get_grid_comp_area_SG): grid_spec does not exist')
     end if
     call read_data(gridfileobj, 'atm_mosaic', mosaic_name)
     call get_grid_ntiles('atm', ntiles)
     allocate(nest_tile_name(ntiles))
     num_nest_tile = 0
     do n = 1, ntiles
        tilefile = read_file_name(mosaic_fileobj(1), 'gridfiles', n)
        call open_grid_file(tilefileobj, grid_dir//tilefile)
        if (global_att_exists(tilefileobj, "nest_grid")) then
           call get_global_attribute(tilefileobj, "nest_grid", attvalue)
           if(trim(attvalue) == "TRUE") then
              num_nest_tile = num_nest_tile + 1
              nest_tile_name(num_nest_tile) = trim(mosaic_name)//'_tile'//char(n+ichar('0'))
           else if(trim(attvalue) .NE. "FALSE") then
              call mpp_error(FATAL,module_name//'/get_grid_comp_area value of global attribute nest_grid in file'//&
                   trim(tilefile)//' should be TRUE or FALSE')
           endif
        end if
        call close_file(tilefileobj)
     end do
     area8(:,:) = 0.0_r8_kind
     if (.not. grid_spec_exists) then
        call mpp_error(FATAL, 'grid2_mod(get_grid_comp_area_SG): grid_spec does not exist')
     end if
     if(variable_exists(gridfileobj,xgrid_name)) then
        ! get the number of the exchange-grid files
        call get_variable_size(gridfileobj,xgrid_name,siz)
        n_xgrid_files = siz(2)
        found_xgrid_files = 0
        ! loop through all exchange grid files
        do n = 1, n_xgrid_files
           ! get the name of the current exchange grid file
           xgrid_file = read_file_name(gridfileobj,xgrid_name,n)
           call open_grid_file(xgrid_fileobj, grid_dir//xgrid_file)
           ! skip the rest of the loop if the name of the current tile isn't found
           ! in the file name, but check this only if there is more than 1 tile
           if(n_xgrid_files>1) then
              if(index(xgrid_file,trim(tile_name))==0) cycle
           endif
           found_xgrid_files = found_xgrid_files + 1
           !---make sure the atmosphere grid is not a nested grid
           is_nest = .false.
           do m = 1, num_nest_tile
              if(index(xgrid_file, trim(nest_tile_name(m))) .NE. 0) then
                 is_nest = .true.
                 exit
              end if
           end do
           if(is_nest) cycle

           ! finally read the exchange grid
           nxgrid = get_mosaic_xgrid_size(xgrid_fileobj)
           if(nxgrid < BUFSIZE) then
              allocate(i1(nxgrid), j1(nxgrid), i2(nxgrid), j2(nxgrid), xgrid_area(nxgrid))
           else
              allocate(i1(BUFSIZE), j1(BUFSIZE), i2(BUFSIZE), j2(BUFSIZE), xgrid_area(BUFSIZE))
           endif
           ibegin = 1
           do l = 1,nxgrid,BUFSIZE
              bsize = min(BUFSIZE, nxgrid-l+1)
              iend = ibegin + bsize - 1
              call get_mosaic_xgrid(xgrid_fileobj, i1(1:bsize), j1(1:bsize), i2(1:bsize), j2(1:bsize), &
                   xgrid_area(1:bsize), ibegin, iend)
              ! and sum the exchange grid areas
              do m = 1, bsize
                 i = i2(m); j = j2(m)
                 if (i<is.or.i>ie) cycle
                 if (j<js.or.j>je) cycle
                 area8(i+i0,j+j0) = area8(i+i0,j+j0) + xgrid_area(m)
              end do
              ibegin = iend + 1
           enddo
           deallocate(i1, j1, i2, j2, xgrid_area)
           call close_file(xgrid_fileobj)
        enddo
        if (found_xgrid_files == 0) &
             call mpp_error(FATAL, 'get_grid_comp_area no xgrid files were found for component '&
             //trim(component)//' (mosaic name is '//trim(mosaic_name)//')')

     endif
     deallocate(nest_tile_name)
  end select ! version
  ! convert area to m2
  area = real(area8*4.0_r8_kind*PI*RADIUS**2, FMS_MOS_KIND_)

  deallocate(area8)

end subroutine GET_GRID_COMP_AREA_SG_

!> @brief return grid cell area for the specified model component and tile on an
!! unstructured domain
subroutine GET_GRID_CELL_AREA_UG_(component, tile, cellarea, SG_domain, UG_domain)
  character(len=*),   intent(in)    :: component !< Component model (atm, lnd, ocn)
  integer         ,   intent(in)    :: tile !< Tile number
  real(kind=FMS_MOS_KIND_),   intent(inout) :: cellarea(:) !< Cell area
  type(domain2d)  ,   intent(in)    :: SG_domain !< Structured Domain
  type(domainUG)  ,   intent(in)    :: UG_domain !< Unstructured Domain
  integer :: is, ie, js, je
  real(kind=FMS_MOS_KIND_), allocatable :: SG_area(:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_area(is:ie, js:je))
  call get_grid_cell_area(component, tile, SG_area, SG_domain)
  call mpp_pass_SG_to_UG(UG_domain, SG_area, cellarea)
  deallocate(SG_area)
end subroutine GET_GRID_CELL_AREA_UG_

!> @brief get the area of the component per grid cell for an unstructured domain
subroutine GET_GRID_COMP_AREA_UG_(component, tile, area, SG_domain, UG_domain)
  character(len=*),   intent(in)    :: component !< Component model (atm, lnd, ocn)
  integer         ,   intent(in)    :: tile !< Tile number
  real(kind=FMS_MOS_KIND_), intent(inout) :: area(:) !< Area of the component
  type(domain2d)  ,   intent(in)    :: SG_domain !< Structured domain
  type(domainUG)  ,   intent(in)    :: UG_domain !< Unstructured domain
  integer :: is, ie, js, je
  real(kind=FMS_MOS_KIND_), allocatable :: SG_area(:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_area(is:ie, js:je))
  call get_grid_comp_area(component, tile, SG_area, SG_domain)
  call mpp_pass_SG_to_UG(UG_domain, SG_area, area)
  deallocate(SG_area)

end subroutine GET_GRID_COMP_AREA_UG_

!> @brief returns arrays of global grid cell boundaries for given model component and
!! mosaic tile number.
subroutine GET_GRID_CELL_VERTICES_1D_(component, tile, glonb, glatb)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer,          intent(in) :: tile !< Tile number
  real(kind=FMS_MOS_KIND_),intent(inout) :: glonb(:),glatb(:) !< Grid cell vertices

  integer                      :: nlon, nlat
  integer                      :: start(4), nread(4)
  real(kind=FMS_MOS_KIND_), allocatable :: tmp(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)
  character(len=FMS_PATH_LEN)      :: tilefile
  type(FmsNetcdfFile_t)  :: tilefileobj

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)
  if (size(glonb(:))/=nlon+1) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_vertices_1D '//&
       'Size of argument "glonb" is not consistent with the grid size')
  if (size(glatb(:))/=nlat+1) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_vertices_1D '//&
       'Size of argument "glatb" is not consistent with the grid size')
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call mpp_error(FATAL, module_name//'/get_grid_cell_vertices_1D '//&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
  endif

  select case(grid_version)
  case(VERSION_GEOLON_T)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_cell_vertices_1D): grid_spec does not exist')
     end if
     select case(trim(component))
     case('ATM','LND')
        call read_data(gridfileobj, 'xb'//lowercase(component(1:1)), glonb)
        call read_data(gridfileobj, 'yb'//lowercase(component(1:1)), glatb)
     case('OCN')
        call read_data(gridfileobj, "gridlon_vert_t", glonb)
        call read_data(gridfileobj, "gridlat_vert_t", glatb)
     end select
  case(VERSION_X_T)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_cell_vertices_1D): grid_spec does not exist')
     end if
     select case(trim(component))
     case('ATM','LND')
        call read_data(gridfileobj, 'xb'//lowercase(component(1:1)), glonb)
        call read_data(gridfileobj, 'yb'//lowercase(component(1:1)), glatb)
     case('OCN')
        allocate (x_vert_t(nlon,1,2), y_vert_t(1,nlat,2) )
        start = 1; nread = 1
        nread(1) = nlon; nread(2) = 1; start(3) = 1
        call read_data(gridfileobj, "x_vert_T", x_vert_t(:,:,1), corner=start, edge_lengths=nread)
        nread(1) = nlon; nread(2) = 1; start(3) = 2
        call read_data(gridfileobj, "x_vert_T", x_vert_t(:,:,2), corner=start, edge_lengths=nread)

        nread(1) = 1; nread(2) = nlat; start(3) = 1
        call read_data(gridfileobj, "y_vert_T", y_vert_t(:,:,1), corner=start, edge_lengths=nread)
        nread(1) = 1; nread(2) = nlat; start(3) = 4
        call read_data(gridfileobj, "y_vert_T", y_vert_t(:,:,2), corner=start, edge_lengths=nread)
        glonb(1:nlon) = x_vert_t(1:nlon,1,1)
        glonb(nlon+1) = x_vert_t(nlon,1,2)
        glatb(1:nlat) = y_vert_t(1,1:nlat,1)
        glatb(nlat+1) = y_vert_t(1,nlat,2)
        deallocate(x_vert_t, y_vert_t)
     end select
  case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES)
     ! get the name of the grid file for the component and tile
     tilefile = read_file_name(mosaic_fileobj(get_component_number(trim(component))), 'gridfiles',tile)
     call open_grid_file(tilefileobj, grid_dir//tilefile)

     start = 1; nread = 1
     nread(1) = 2*nlon+1
     allocate( tmp(2*nlon+1,1) )
     call read_data(tilefileobj, "x", tmp, corner=start, edge_lengths=nread)
     glonb(1:nlon+1) = tmp(1:2*nlon+1:2,1)
     deallocate(tmp)
     allocate(tmp(1,2*nlat+1))

     start = 1; nread = 1
     nread(2) = 2*nlat+1
     call read_data(tilefileobj, "y", tmp, corner=start, edge_lengths=nread)
     glatb(1:nlat+1) = tmp(1,1:2*nlat+1:2)
     deallocate(tmp)
     call close_file(tilefileobj)
  end select
end subroutine GET_GRID_CELL_VERTICES_1D_

!> @brief returns cell vertices for the specified model component and mosaic tile number
subroutine GET_GRID_CELL_VERTICES_2D_(component, tile, lonb, latb, domain)
  character(len=*),         intent(in) :: component !< Component model (atm, lnd, ocn)
  integer,                  intent(in) :: tile !< Tile number
  real(kind=FMS_MOS_KIND_), intent(inout) :: lonb(:,:),latb(:,:) !< Cell vertices
  type(domain2d), optional, intent(in) :: domain !< Domain

  ! local vars
  integer :: nlon, nlat
  integer :: i,j
  real(kind=FMS_MOS_KIND_), allocatable :: buffer(:), tmp(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0,j0 ! offsets for coordinates
  integer :: isg, jsg
  integer :: start(4), nread(4)
  character(len=FMS_PATH_LEN)      :: tilefile
  type(FmsNetcdfFile_t)  :: tilefileobj

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)

  if (present(domain)) then
    call mpp_get_compute_domain(domain,is,ie,js,je)
  else
    is = 1 ; ie = nlon
    js = 1 ; je = nlat
    !--- domain normally should be present
    call mpp_error (NOTE, module_name//'/get_grid_cell_vertices '//&
       'domain is not present, global data will be read')
  endif
  i0 = -is+1; j0 = -js+1

  ! verify that lonb and latb sizes are consistent with the size of domain
  if (size(lonb,1)/=ie-is+2.or.size(lonb,2)/=je-js+2) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_vertices '//&
       'Size of argument "lonb" is not consistent with the domain size')
  if (size(latb,1)/=ie-is+2.or.size(latb,2)/=je-js+2) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_vertices '//&
       'Size of argument "latb" is not consistent with the domain size')
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call mpp_error(FATAL, module_name//'/get_grid_cell_vertices '//&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
  endif

     !! use lonb, latb as r4
     select case(grid_version)
     case(VERSION_GEOLON_T)
       if (.not. grid_spec_exists) then
         call mpp_error(FATAL, 'grid2_mod(get_grid_cell_vertices_2D): grid_spec does not exist')
       end if
       select case(component)
       case('ATM','LND')
          allocate(buffer(max(nlon,nlat)+1))
          ! read coordinates of grid cell vertices
          call read_data(gridfileobj, 'xb'//lowercase(component(1:1)), buffer(1:nlon+1))
          do j = js, je+1
             do i = is, ie+1
                lonb(i+i0,j+j0) = buffer(i)
             enddo
          enddo
          call read_data(gridfileobj, 'yb'//lowercase(component(1:1)), buffer(1:nlat+1))
          do j = js, je+1
             do i = is, ie+1
                latb(i+i0,j+j0) = buffer(j)
             enddo
          enddo
          deallocate(buffer)
       case('OCN')
          if (present(domain)) then
             start = 1; nread = 1
             start(1) = is; start(2) = js
             nread(1) = ie-is+2; nread(2) = je-js+2
             call read_data(gridfileobj, "geolon_vert_t", lonb, corner=start, edge_lengths=nread)
             call read_data(gridfileobj, "geolat_vert_t", latb, corner=start, edge_lengths=nread)
           else
             call read_data(gridfileobj, "geolon_vert_t", lonb)
             call read_data(gridfileobj, "geolat_vert_t", latb)
           endif
       end select
     case(VERSION_X_T)
       if (.not. grid_spec_exists) then
         call mpp_error(FATAL, 'grid2_mod(get_grid_cell_vertices_2D): grid_spec does not exist')
       end if
       select case(component)
       case('ATM','LND')
          allocate(buffer(max(nlon,nlat)+1))
          ! read coordinates of grid cell vertices
          call read_data(gridfileobj, 'xb'//lowercase(component(1:1)), buffer(1:nlon+1))
          do j = js, je+1
             do i = is, ie+1
                lonb(i+i0,j+j0) = buffer(i)
             enddo
          enddo
          call read_data(gridfileobj, 'yb'//lowercase(component(1:1)), buffer(1:nlat+1))
          do j = js, je+1
             do i = is, ie+1
                latb(i+i0,j+j0) = buffer(j)
             enddo
          enddo
          deallocate(buffer)
       case('OCN')
          nlon=ie-is+1; nlat=je-js+1
          allocate (x_vert_t(nlon,nlat,4), y_vert_t(nlon,nlat,4) )
          call read_data(gridfileobj, 'x_vert_T', x_vert_t)
          call read_data(gridfileobj, 'y_vert_T', y_vert_t)
          lonb(1:nlon,1:nlat) = x_vert_t(1:nlon,1:nlat,1)
          lonb(nlon+1,1:nlat) = x_vert_t(nlon,1:nlat,2)
          lonb(1:nlon,nlat+1) = x_vert_t(1:nlon,nlat,4)
          lonb(nlon+1,nlat+1) = x_vert_t(nlon,nlat,3)
          latb(1:nlon,1:nlat) = y_vert_t(1:nlon,1:nlat,1)
          latb(nlon+1,1:nlat) = y_vert_t(nlon,1:nlat,2)
          latb(1:nlon,nlat+1) = y_vert_t(1:nlon,nlat,4)
          latb(nlon+1,nlat+1) = y_vert_t(nlon,nlat,3)
          deallocate(x_vert_t, y_vert_t)
       end select
     case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES)
       ! get the name of the grid file for the component and tile
       tilefile = read_file_name(mosaic_fileobj(get_component_number(trim(component))), 'gridfiles',tile)
       call open_grid_file(tilefileobj, grid_dir//tilefile)
       if(PRESENT(domain)) then
          call mpp_get_global_domain(domain, xbegin=isg, ybegin=jsg)
          start = 1; nread = 1
          start(1) = 2*(is-isg+1) - 1; nread(1) = 2*(ie-is)+3
          start(2) = 2*(js-jsg+1) - 1; nread(2) = 2*(je-js)+3
          allocate(tmp(nread(1), nread(2)) )
          call read_data(tilefileobj, "x", tmp, corner=start, edge_lengths=nread)
          do j = 1, je-js+2
             do i = 1, ie-is+2
                lonb(i,j) = tmp(2*i-1,2*j-1)
             enddo
          enddo
          call read_data(tilefileobj, "y", tmp, corner=start, edge_lengths=nread)
          do j = 1, je-js+2
             do i = 1, ie-is+2
                latb(i,j) = tmp(2*i-1,2*j-1)
             enddo
          enddo
       else
          allocate(tmp(2*nlon+1,2*nlat+1))
          call read_data(tilefileobj, "x", tmp)
          do j = js, je+1
             do i = is, ie+1
                lonb(i+i0,j+j0) = tmp(2*i-1,2*j-1)
             end do
          end do
          call read_data(tilefileobj, "y", tmp)
          do j = js, je+1
             do i = is, ie+1
                latb(i+i0,j+j0) = tmp(2*i-1,2*j-1)
             end do
          end do
       endif
       deallocate(tmp)
       call close_file(tilefileobj)
    end select ! end grid_version
  end subroutine GET_GRID_CELL_VERTICES_2D_

!> @brief returns cell vertices for the specified model component and mosaic tile number for
!! an unstructured domain
subroutine GET_GRID_CELL_VERTICES_UG_(component, tile, lonb, latb, SG_domain, UG_domain)
  character(len=*),         intent(in) :: component !< Component model (atm, lnd, ocn)
  integer,                  intent(in) :: tile !< Tile number
  real(kind=FMS_MOS_KIND_),  intent(inout) :: lonb(:,:),latb(:,:) ! The second dimension is 4
  type(domain2d)  ,   intent(in)       :: SG_domain !< Structured domain
  type(domainUG)  ,   intent(in)       :: UG_domain !< Unstructured domain
  integer :: is, ie, js, je, i, j
  real(kind=FMS_MOS_KIND_), allocatable :: SG_lonb(:,:), SG_latb(:,:), tmp(:,:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_lonb(is:ie+1, js:je+1))
  allocate(SG_latb(is:ie+1, js:je+1))
  allocate(tmp(is:ie,js:je,4))
  call get_grid_cell_vertices(component, tile, SG_lonb, SG_latb, SG_domain)
  do j = js, je
     do i = is, ie
        tmp(i,j,1) = SG_lonb(i,j)
        tmp(i,j,2) = SG_lonb(i+1,j)
        tmp(i,j,3) = SG_lonb(i+1,j+1)
        tmp(i,j,4) = SG_lonb(i,j+1)
     enddo
  enddo
  call mpp_pass_SG_to_UG(UG_domain, tmp, lonb)
  do j = js, je
     do i = is, ie
        tmp(i,j,1) = SG_latb(i,j)
        tmp(i,j,2) = SG_latb(i+1,j)
        tmp(i,j,3) = SG_latb(i+1,j+1)
        tmp(i,j,4) = SG_latb(i,j+1)
     enddo
  enddo
  call mpp_pass_SG_to_UG(UG_domain, tmp, latb)


  deallocate(SG_lonb, SG_latb, tmp)
end subroutine GET_GRID_CELL_VERTICES_UG_

!> @brief returns grid cell centers given model component and mosaic tile number
subroutine GET_GRID_CELL_CENTERS_1D_(component, tile, glon, glat)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real(kind=FMS_MOS_KIND_), intent(inout) :: glon(:),glat(:) !< Grid cell centers

  integer                      :: nlon, nlat
  integer                      :: start(4), nread(4)
  real(kind=FMS_MOS_KIND_), allocatable :: tmp(:,:)
  character(len=FMS_PATH_LEN)      :: tilefile
  type(FmsNetcdfFile_t)  :: tilefileobj

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)
  if (size(glon(:))/=nlon) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_centers_1D '//&
       'Size of argument "glon" is not consistent with the grid size')
  if (size(glat(:))/=nlat) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_centers_1D '//&
       'Size of argument "glat" is not consistent with the grid size')
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call mpp_error(FATAL, module_name//'/get_grid_cell_centers_1D '//&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
  endif

  select case(grid_version)
  case(VERSION_GEOLON_T)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_cell_centers_1D): grid_spec does not exist')
     end if
     select case(trim(component))
     case('ATM','LND')
        call read_data(gridfileobj, 'xt'//lowercase(component(1:1)), glon)
        call read_data(gridfileobj, 'yt'//lowercase(component(1:1)), glat)
     case('OCN')
        call read_data(gridfileobj, "gridlon_t", glon)
        call read_data(gridfileobj, "gridlat_t", glat)
     end select
  case(VERSION_X_T)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_cell_centers_1D): grid_spec does not exist')
     end if
     select case(trim(component))
     case('ATM','LND')
        call read_data(gridfileobj, 'xt'//lowercase(component(1:1)), glon)
        call read_data(gridfileobj, 'yt'//lowercase(component(1:1)), glat)
     case('OCN')
        call read_data(gridfileobj, "grid_x_T", glon)
        call read_data(gridfileobj, "grid_y_T", glat)
     end select
  case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES)
     ! get the name of the grid file for the component and tile
     tilefile = read_file_name(mosaic_fileobj(get_component_number(trim(component))), 'gridfiles',tile)
     call open_grid_file(tilefileobj, grid_dir//tilefile)

     start = 1; nread = 1
     nread(1) = 2*nlon+1; start(2) = 2
     allocate( tmp(2*nlon+1,1) )
     call read_data(tilefileobj, "x", tmp, corner=start, edge_lengths=nread)
     glon(1:nlon) = tmp(2:2*nlon:2,1)
     deallocate(tmp)
     allocate(tmp(1, 2*nlat+1))

     start = 1; nread = 1
     nread(2) = 2*nlat+1; start(1) = 2
     call read_data(tilefileobj, "y", tmp, corner=start, edge_lengths=nread)
     glat(1:nlat) = tmp(1,2:2*nlat:2)
     deallocate(tmp)
     call close_file(tilefileobj)
  end select
end subroutine GET_GRID_CELL_CENTERS_1D_

!> @brief returns grid cell centers given model component and mosaic tile number
subroutine GET_GRID_CELL_CENTERS_2D_(component, tile, lon, lat, domain)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real(kind=FMS_MOS_KIND_), intent(inout) :: lon(:,:),lat(:,:) !< Grid cell centers
  type(domain2d), intent(in), optional :: domain !< Domain
  ! local vars
  integer :: nlon, nlat
  integer :: i,j
  real(kind=FMS_MOS_KIND_), allocatable :: buffer(:),tmp(:,:)
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0,j0 ! offsets for coordinates
  integer :: isg, jsg
  integer :: start(4), nread(4)
  character(len=FMS_PATH_LEN)      :: tilefile
  type(FmsNetcdfFile_t)  :: tilefileobj

  call get_grid_size_for_one_tile(component, tile, nlon, nlat)
  if (present(domain)) then
    call mpp_get_compute_domain(domain,is,ie,js,je)
  else
    is = 1 ; ie = nlon
    js = 1 ; je = nlat
    !--- domain normally should be present
    call mpp_error (NOTE, module_name//'/get_grid_cell_centers '//&
       'domain is not present, global data will be read')
  endif
  i0 = -is+1; j0 = -js+1

  ! verify that lon and lat sizes are consistent with the size of domain
  if (size(lon,1)/=ie-is+1.or.size(lon,2)/=je-js+1) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_centers '//&
       'Size of array "lon" is not consistent with the domain size')
  if (size(lat,1)/=ie-is+1.or.size(lat,2)/=je-js+1) &
       call mpp_error (FATAL, module_name//'/get_grid_cell_centers '//&
       'Size of array "lat" is not consistent with the domain size')
  if(trim(component) .NE. 'ATM' .AND. component .NE. 'LND' .AND. component .NE. 'OCN') then
     call mpp_error(FATAL, module_name//'/get_grid_cell_vertices '//&
          'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
  endif

  select case(grid_version)
  case(VERSION_GEOLON_T)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_cell_centers_2D): grid_spec does not exist')
     end if
     select case (trim(component))
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)))
        ! read coordinates of grid cell vertices
        call read_data(gridfileobj, 'xt'//lowercase(component(1:1)), buffer(1:nlon))
        do j = js,je
        do i = is,ie
           lon(i+i0,j+j0) = buffer(i)
        enddo
        enddo
        call read_data(gridfileobj, 'yt'//lowercase(component(1:1)), buffer(1:nlat))
        do j = js,je
        do i = is,ie
           lat(i+i0,j+j0) = buffer(j)
        enddo
        enddo
        deallocate(buffer)
     case('OCN')
        call read_data(gridfileobj, 'geolon_t', lon)
        call read_data(gridfileobj, 'geolat_t', lat)
     end select
  case(VERSION_X_T)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_cell_centers_2D): grid_spec does not exist')
     end if
     select case(trim(component))
     case('ATM','LND')
        allocate(buffer(max(nlon,nlat)))
        ! read coordinates of grid cell vertices
        call read_data(gridfileobj, 'xt'//lowercase(component(1:1)), buffer(1:nlon))
        do j = js,je
        do i = is,ie
           lon(i+i0,j+j0) = buffer(i)
        enddo
        enddo
        call read_data(gridfileobj, 'yt'//lowercase(component(1:1)), buffer(1:nlat))
        do j = js,je
        do i = is,ie
           lat(i+i0,j+j0) = buffer(j)
        enddo
        enddo
        deallocate(buffer)
     case('OCN')
        call read_data(gridfileobj, 'x_T', lon)
        call read_data(gridfileobj, 'y_T', lat)
     end select
  case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES) ! mosaic grid file
     ! get the name of the grid file for the component and tile
     tilefile = read_file_name(mosaic_fileobj(get_component_number(trim(component))), 'gridfiles',tile)
     call open_grid_file(tilefileobj, grid_dir//tilefile)

     if(PRESENT(domain)) then
        call mpp_get_global_domain(domain, xbegin=isg, ybegin=jsg)
        start = 1; nread = 1
        start(1) = 2*(is-isg+1) - 1; nread(1) = 2*(ie-is)+3
        start(2) = 2*(js-jsg+1) - 1; nread(2) = 2*(je-js)+3
        allocate(tmp(nread(1), nread(2)))
        call read_data(tilefileobj, "x", tmp, corner=start, edge_lengths=nread)
        do j = 1, je-js+1
           do i = 1, ie-is+1
              lon(i,j) = tmp(2*i,2*j)
           enddo
        enddo
        call read_data(tilefileobj, "y", tmp, corner=start, edge_lengths=nread)
        do j = 1, je-js+1
           do i = 1, ie-is+1
              lat(i,j) = tmp(2*i,2*j)
           enddo
        enddo
     else
        allocate(tmp(2*nlon+1,2*nlat+1))
        call read_data(tilefileobj, 'x', tmp)
        do j = js,je
           do i = is,ie
              lon(i+i0,j+j0) = tmp(2*i,2*j)
           end do
        end do
        call read_data(tilefileobj, 'y', tmp)
        do j = js,je
           do i = is,ie
              lat(i+i0,j+j0) = tmp(2*i,2*j)
           end do
        end do
        deallocate(tmp)
     endif
     call close_file(tilefileobj)
  end select
end subroutine GET_GRID_CELL_CENTERS_2D_

!> @brief returns grid cell centers given model component and mosaic tile number
!! for unstructured domain
subroutine GET_GRID_CELL_CENTERS_UG_(component, tile, lon, lat, SG_domain, UG_domain)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real(kind=FMS_MOS_KIND_), intent(inout) :: lon(:),lat(:) !< Grid cell centers
  type(domain2d)  ,   intent(in) :: SG_domain !< Structured domain
  type(domainUG)  ,   intent(in) :: UG_domain !< Unstructured domain
  integer :: is, ie, js, je
  real(kind=FMS_MOS_KIND_), allocatable :: SG_lon(:,:), SG_lat(:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_lon(is:ie, js:je))
  allocate(SG_lat(is:ie, js:je))
  call get_grid_cell_centers(component, tile, SG_lon, SG_lat, SG_domain)
  call mpp_pass_SG_to_UG(UG_domain, SG_lon, lon)
  call mpp_pass_SG_to_UG(UG_domain, SG_lat, lat)
  deallocate(SG_lon, SG_lat)
end subroutine GET_GRID_CELL_CENTERS_UG_

!> @}
! close documentation grouping
