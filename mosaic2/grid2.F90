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
!> @defgroup grid2_mod grid2_mod
!> @ingroup mosaic2
!> @brief Routines for grid calculations, using @ref fms2_io

!> @file
!> @brief File for @ref grid2_mod

module grid2_mod

use mpp_mod, only : mpp_root_pe, mpp_error, uppercase, lowercase, FATAL, NOTE
use constants_mod, only : PI, radius
use fms2_io_mod, only : get_global_attribute, read_data, global_att_exists, &
                        variable_exists, file_exists,  open_file, close_file, get_variable_size, &
                        FmsNetcdfFile_t, string => string2
use mosaic2_mod, only : get_mosaic_ntiles, get_mosaic_xgrid_size, get_mosaic_grid_sizes, &
     get_mosaic_xgrid, calc_mosaic_grid_area, calc_mosaic_grid_great_circle_area

! the following two use statement are only needed for define_cube_mosaic
use mpp_domains_mod, only : domain2d, mpp_define_mosaic, mpp_get_compute_domain, &
                            mpp_get_global_domain, domainUG, mpp_pass_SG_to_UG
use mosaic2_mod, only : get_mosaic_ncontacts, get_mosaic_contact

implicit none;private

! ==== public interfaces =====================================================
! grid dimension inquiry subroutines
public :: get_great_circle_algorithm ! returns great_circle_algorithm
public :: get_grid_ntiles ! returns number of tiles
public :: get_grid_size   ! returns horizontal sizes of the grid
! grid geometry inquiry subroutines
public :: get_grid_cell_centers
public :: get_grid_cell_vertices
! grid area inquiry subroutines
public :: get_grid_cell_area
public :: get_grid_comp_area
! decompose cubed sphere domains -- probably does not belong here, but it should
! be in some place available for component models
public :: define_cube_mosaic
public :: grid_init
public :: grid_end
! ==== end of public interfaces ==============================================

!> Gets the size of the grid for one or all tiles
!> @ingroup grid2_mod
interface get_grid_size
   module procedure get_grid_size_for_all_tiles
   module procedure get_grid_size_for_one_tile
end interface

!> Gets arrays of global grid cell boundaries for given model component and
!! mosaic tile number
!> @ingroup grid2_mod
interface get_grid_cell_vertices
   module procedure get_grid_cell_vertices_1D
   module procedure get_grid_cell_vertices_2D
   module procedure get_grid_cell_vertices_UG
end interface

!> Gets grid cell centers
!> @ingroup grid2_mod
interface get_grid_cell_centers
   module procedure get_grid_cell_centers_1D
   module procedure get_grid_cell_centers_2D
   module procedure get_grid_cell_centers_UG
end interface

!> Finds area of a grid cell
!> @ingroup grid2_mod
interface get_grid_cell_area
   module procedure get_grid_cell_area_SG
   module procedure get_grid_cell_area_UG
end interface get_grid_cell_area

!> Gets the area of a given component per grid cell
!> @ingroup grid2_mod
interface get_grid_comp_area
   module procedure get_grid_comp_area_SG
   module procedure get_grid_comp_area_UG
end interface get_grid_comp_area

!> @addtogroup grid2_mod
!> @{
! ==== module constants ======================================================
character(len=*), parameter :: &
     module_name = 'grid2_mod'

! Include variable "version" to be written to log file.
#include<file_version.h>

character(len=*), parameter :: &
     grid_dir  = 'INPUT/',     &      !< root directory for all grid files
     grid_file = 'INPUT/grid_spec.nc' !< name of the grid spec file

integer, parameter :: &
     MAX_NAME = 256,  & !< max length of the variable names
     MAX_FILE = 1024, & !< max length of the file names
     VERSION_0 = 0,   &
     VERSION_1 = 1,   &
     VERSION_2 = 2,   &
     VERSION_3 = 3

integer, parameter :: BUFSIZE = 1048576  !< This is used to control memory usage in get_grid_comp_area
                                         !! We may change this to a namelist variable is needed.

! ==== module variables ======================================================
integer :: grid_version = -1
logical :: great_circle_algorithm = .FALSE.
logical :: module_is_initialized = .FALSE.
logical :: grid_spec_exists = .TRUE.
type(FmsNetcdfFile_t) :: gridfileobj
type(FmsNetcdfFile_t), dimension(3) :: mosaic_fileobj

contains

!> @brief Initialize the grid2 module
subroutine grid_init
   if (module_is_initialized) return
   if (.not. file_exists(grid_file)) then
       module_is_initialized = .TRUE.
       grid_spec_exists = .FALSE.
       return
   endif
   call open_grid_file(gridfileobj, grid_file)
   great_circle_algorithm = get_great_circle_algorithm()
   grid_version = get_grid_version(gridfileobj)
   if (grid_version == VERSION_2) call open_component_mosaics
   if (grid_version == VERSION_3) call assign_component_mosaics
   module_is_initialized = .TRUE.
end subroutine grid_init

!> @brief Shutdown the grid2 module
subroutine grid_end
   if (grid_spec_exists) then
       if (grid_version == VERSION_2) call close_component_mosaics
       call close_file(gridfileobj)
   endif
end subroutine grid_end

!> @brief Determine if we are using the great circle algorithm
!! @return Logical flag describing if we are using the great circlealgorithm
function get_great_circle_algorithm()
   character(len=128)     :: attvalue
   logical :: get_great_circle_algorithm

   get_great_circle_algorithm = .false.
   if (.not. grid_spec_exists) return
   if (global_att_exists(gridfileobj, "great_circle_algorithm")) then
      call get_global_attribute(gridfileobj, "great_circle_algorithm", attvalue)
      if(trim(attvalue) == "TRUE") then
         get_great_circle_algorithm = .true.
      else if(trim(attvalue) .NE. "FALSE") then
         call mpp_error(FATAL, module_name//'/get_great_circle_algorithm value of global attribute "great_circle_algorthm" in file'// &
                   trim(grid_file)//' should be TRUE or FALSE')
      endif
   endif
end function get_great_circle_algorithm

!> @brief Open a grid file
subroutine open_grid_file(myfileobj, myfilename)
  type(FmsNetcdfFile_t), intent(out)  :: myfileobj !< File object of grid file
  character(len=*), intent(in) :: myfilename!< Name of the grid file
   if(.not. open_file(myfileobj, myfilename, 'read')) then
      call mpp_error(FATAL, 'grid2_mod(open_grid_file):Error in opening file '//trim(myfilename))
   endif
end subroutine open_grid_file

!> @brief Open a mosaic file
subroutine open_mosaic_file(mymosaicfileobj, component)
  type(FmsNetcdfFile_t), intent(out)  :: mymosaicfileobj !< File object returned
  character(len=3), intent(in)        :: component !< Component (atm, lnd, etc.)

  character(len=MAX_FILE) :: mosaicfilename
  if (.not. grid_spec_exists) then
    call mpp_error(FATAL, 'grid2_mod(open_mosaic_file): grid_spec does not exist')
  end if
  call read_data(gridfileobj,trim(lowercase(component))//'_mosaic_file', mosaicfilename)
  call open_grid_file(mymosaicfileobj, grid_dir//trim(mosaicfilename))
end subroutine open_mosaic_file

!> @brief Read a tile file name from a netcdf file
!! @return Name of the file as a string
function read_file_name(thisfileobj, filevar, level)
  type(FmsNetcdfFile_t), intent(in)  :: thisfileobj !< File object of file
  character(len=*), intent(in) :: filevar!< Variable containing file names
  integer, intent(in) :: level !< Level of tile file
  integer, dimension(2) :: file_list_size
  character(len=MAX_FILE) :: read_file_name
  character(len=MAX_FILE), dimension(:), allocatable :: file_names

  call get_variable_size(thisfileobj, filevar, file_list_size)
  allocate(file_names(file_list_size(2)))
  call read_data(thisfileobj, filevar, file_names)
  read_file_name = file_names(level)
  deallocate(file_names)
end function read_file_name

!> @brief Get the grid version from a file object
!! @return An integer representation of the grid version
function get_grid_version(fileobj)
  type(FmsNetcdfFile_t), intent(in) :: fileobj !< File object of grid file
  integer :: get_grid_version

  if(grid_version<0) then
    if(variable_exists(fileobj, 'geolon_t')) then
       get_grid_version = VERSION_0
    else if(variable_exists(fileobj, 'x_T')) then
       get_grid_version = VERSION_1
    else if(variable_exists(fileobj, 'ocn_mosaic_file') ) then
       get_grid_version = VERSION_2
    else if(variable_exists(fileobj, 'gridfiles') ) then
       get_grid_version = VERSION_3
    else
       call mpp_error(FATAL, module_name//'/get_grid_version '//&
            'Can''t determine the version of the grid spec: none of "x_T", "geolon_t", or "ocn_mosaic_file" exist in file "'//trim(grid_file)//'"')
    endif
  endif
end function get_grid_version

!> @brief Assign the component mosaic files if grid_spec is Version 3
subroutine assign_component_mosaics
    if (.not. grid_spec_exists) then
      call mpp_error(FATAL, 'grid2_mod(assign_component_mosaics): grid_spec does not exist')
    end if
    mosaic_fileobj(1) = gridfileobj
    mosaic_fileobj(2) = gridfileobj
    mosaic_fileobj(3) = gridfileobj
end subroutine assign_component_mosaics

!> @brief Open the component mosaic files for atm, lnd, and ocn
subroutine open_component_mosaics
    if (.not. grid_spec_exists) then
      call mpp_error(FATAL, 'grid2_mod(open_component_mosaics): grid_spec does not exist')
    end if
    if (variable_exists(gridfileobj, 'atm_mosaic_file')) call open_mosaic_file(mosaic_fileobj(1), 'atm')
    if (variable_exists(gridfileobj, 'ocn_mosaic_file')) call open_mosaic_file(mosaic_fileobj(2), 'ocn')
    if (variable_exists(gridfileobj, 'lnd_mosaic_file')) call open_mosaic_file(mosaic_fileobj(3), 'lnd')
end subroutine open_component_mosaics

!> @brief Close the component mosaic files for atm, lnd, and ocn
subroutine close_component_mosaics
    if (.not. grid_spec_exists) then
      call mpp_error(FATAL, 'grid2_mod(close_component_mosaics): grid_spec does not exist')
    end if
    if (variable_exists(gridfileobj, 'atm_mosaic_file')) call close_file(mosaic_fileobj(1))
    if (variable_exists(gridfileobj, 'ocn_mosaic_file')) call close_file(mosaic_fileobj(2))
    if (variable_exists(gridfileobj, 'lnd_mosaic_file')) call close_file(mosaic_fileobj(3))
end subroutine close_component_mosaics

!> @brief Get the component number of a model component (atm, lnd, ocn)
!! @return Integer component number
function get_component_number(component)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer :: get_component_number
    select case(lowercase(component))
    case('atm')
        get_component_number = 1
    case('ocn')
        get_component_number = 2
    case('lnd')
        get_component_number = 3
    end select
end function get_component_number

!> @brief returns number of tiles for a given component
subroutine get_grid_ntiles(component,ntiles)
  character(len=*)     :: component !< Component model (atm, lnd, ocn)
  integer, intent(out) :: ntiles !< Number of tiles

  select case (grid_version)
  case(VERSION_0,VERSION_1)
     ntiles = 1
  case(VERSION_2, VERSION_3)
     ntiles = get_mosaic_ntiles(mosaic_fileobj(get_component_number(trim(component))))
  end select
end subroutine get_grid_ntiles

!> @brief returns size of the grid for each of the tiles
subroutine get_grid_size_for_all_tiles(component,nx,ny)
  character(len=*)     :: component !< Component model (atm, lnd, ocn)
  integer, intent(inout) :: nx(:),ny(:) !< Grid size in x and y

  ! local vars
  integer :: siz(2) ! for the size of external fields
  character(len=MAX_NAME) :: varname1, varname2

  varname1 = 'AREA_'//trim(uppercase(component))

  select case (grid_version)
  case(VERSION_0,VERSION_1)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_size_for_all_tiles): grid_spec does not exist')
     end if
     call get_variable_size(gridfileobj, varname1, siz)
     nx(1) = siz(1); ny(1)=siz(2)
  case(VERSION_2, VERSION_3) ! mosaic file
     call get_mosaic_grid_sizes(mosaic_fileobj(get_component_number(trim(component))),nx,ny)
  end select
end subroutine get_grid_size_for_all_tiles

!> @brief returns size of the grid for one of the tiles
subroutine get_grid_size_for_one_tile(component,tile,nx,ny)
  character(len=*)       :: component !< Component model (atm, lnd, ocn)
  integer, intent(in)    :: tile !< Tile number
  integer, intent(inout) :: nx,ny !< Grid size in x and y

  ! local vars
  integer, allocatable :: nnx(:), nny(:)
  integer :: ntiles

  call get_grid_ntiles(component, ntiles)
  if(tile>0.and.tile<=ntiles) then
     allocate(nnx(ntiles),nny(ntiles))
     call get_grid_size_for_all_tiles(component,nnx,nny)
     nx = nnx(tile); ny = nny(tile)
     deallocate(nnx,nny)
  else
     call mpp_error(FATAL, 'get_grid_size'//&
          'requested tile index '//trim(string(tile))//' is out of bounds (1:'//trim(string(ntiles))//')')
  endif
end subroutine get_grid_size_for_one_tile

!> @brief return grid cell area for the specified model component and tile
subroutine get_grid_cell_area_SG(component, tile, cellarea, domain)
  character(len=*), intent(in)    :: component !< Component model (atm, lnd, ocn)
  integer         , intent(in)    :: tile !< Tile number
  real            , intent(inout) :: cellarea(:,:) !< Cell area
  type(domain2d)  , intent(in), optional :: domain !< Domain

  ! local vars
  integer :: nlon, nlat
  real, allocatable :: glonb(:,:), glatb(:,:)

  select case(grid_version)
  case(VERSION_0,VERSION_1)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_cell_area_SG): grid_spec does not exist')
     end if
     select case(trim(component))
     case('LND')
        call read_data(gridfileobj, 'AREA_LND_CELL', cellarea)
     case('ATM','OCN')
        call read_data(gridfileobj, 'AREA_'//trim(uppercase(component)),cellarea)
     case default
        call mpp_error(FATAL, module_name//'/get_grid_cell_area'//&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
     end select
     ! convert area to m2
     cellarea = cellarea*4.*PI*radius**2
  case(VERSION_2, VERSION_3)
     if (present(domain)) then
        call mpp_get_compute_domain(domain,xsize=nlon,ysize=nlat)
     else
        call get_grid_size(component,tile,nlon,nlat)
     endif
     allocate(glonb(nlon+1,nlat+1),glatb(nlon+1,nlat+1))
     call get_grid_cell_vertices(component, tile, glonb, glatb, domain)
     if (great_circle_algorithm) then
        call calc_mosaic_grid_great_circle_area(glonb*pi/180.0, glatb*pi/180.0, cellarea)
     else
        call calc_mosaic_grid_area(glonb*pi/180.0, glatb*pi/180.0, cellarea)
     end if
     deallocate(glonb,glatb)
  end select

end subroutine get_grid_cell_area_SG

!> @brief get the area of the component per grid cell
subroutine get_grid_comp_area_SG(component,tile,area,domain)
  character(len=*) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real, intent(inout) :: area(:,:) !< Area of grid cell
  type(domain2d), intent(in), optional :: domain !< Domain
  ! local vars
  integer :: n_xgrid_files ! number of exchange grid files in the mosaic
  integer :: siz(2), nxgrid
  integer :: i,j,m,n
  integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
  real, allocatable :: xgrid_area(:)
  real, allocatable :: rmask(:,:)
  character(len=MAX_NAME) :: &
     xgrid_name, & ! name of the variable holding xgrid names
     tile_name,  & ! name of the tile
     xgrid_file, & ! name of the current xgrid file
     mosaic_name,& ! name of the mosaic
     tilefile
  character(len=4096)     :: attvalue
  character(len=MAX_NAME), allocatable :: nest_tile_name(:)
  character(len=MAX_NAME) :: varname1, varname2
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0, j0 ! offsets for x and y, respectively
  integer :: num_nest_tile, ntiles
  logical :: is_nest
  integer :: found_xgrid_files ! how many xgrid files we actually found in the grid spec
  integer :: ibegin, iend, bsize, l
  type(FmsNetcdfFile_t) :: tilefileobj, xgrid_fileobj

  select case (grid_version)
  case(VERSION_0,VERSION_1)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_comp_area_SG): grid_spec does not exist')
     end if
     select case(component)
     case('ATM')
        call read_data(gridfileobj,'AREA_ATM',area)
     case('OCN')
        allocate(rmask(size(area,1),size(area,2)))
        call read_data(gridfileobj,'AREA_OCN',area)
        call read_data(gridfileobj,'wet',     rmask)
        area = area*rmask
        deallocate(rmask)
     case('LND')
        call read_data(gridfileobj,'AREA_LND',area)
     case default
        call mpp_error(FATAL, module_name//'/get_grid_comp_area'//&
             'Illegal component name "'//trim(component)//'": must be one of ATM, LND, or OCN')
     end select
  case(VERSION_2, VERSION_3) ! mosaic gridspec
     select case (component)
     case ('ATM')
        ! just read the grid cell area and return
        call get_grid_cell_area(component,tile,area)
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
     if (size(area,1)/=ie-is+1.or.size(area,2)/=je-js+1) &
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
              call mpp_error(FATAL, module_name//'/get_grid_comp_area value of global attribute nest_grid in file'// &
                   trim(tilefile)//' should be TRUE or FALSE')
           endif
        end if
        call close_file(tilefileobj)
     end do
     area(:,:) = 0.
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
                 area(i+i0,j+j0) = area(i+i0,j+j0) + xgrid_area(m)
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
  area = area*4.*PI*radius**2
end subroutine get_grid_comp_area_SG

!> @brief return grid cell area for the specified model component and tile on an
!! unstructured domain
subroutine get_grid_cell_area_UG(component, tile, cellarea, SG_domain, UG_domain)
  character(len=*),   intent(in)    :: component !< Component model (atm, lnd, ocn)
  integer         ,   intent(in)    :: tile !< Tile number
  real            ,   intent(inout) :: cellarea(:) !< Cell area
  type(domain2d)  ,   intent(in)    :: SG_domain !< Structured Domain
  type(domainUG)  ,   intent(in)    :: UG_domain !< Unstructured Domain
  integer :: is, ie, js, je
  real, allocatable :: SG_area(:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_area(is:ie, js:je))
  call get_grid_cell_area_SG(component, tile, SG_area, SG_domain)
  call mpp_pass_SG_to_UG(UG_domain, SG_area, cellarea)
  deallocate(SG_area)
end subroutine get_grid_cell_area_UG

!> @brief get the area of the component per grid cell for an unstructured domain
subroutine get_grid_comp_area_UG(component, tile, area, SG_domain, UG_domain)
  character(len=*),   intent(in)    :: component !< Component model (atm, lnd, ocn)
  integer         ,   intent(in)    :: tile !< Tile number
  real            ,   intent(inout) :: area(:) !< Area of the component
  type(domain2d)  ,   intent(in)    :: SG_domain !< Structured domain
  type(domainUG)  ,   intent(in)    :: UG_domain !< Unstructured domain
  integer :: is, ie, js, je
  real, allocatable :: SG_area(:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_area(is:ie, js:je))
  call get_grid_comp_area_SG(component, tile, SG_area, SG_domain)
  call mpp_pass_SG_to_UG(UG_domain, SG_area, area)
  deallocate(SG_area)

end subroutine get_grid_comp_area_UG

!> @brief returns arrays of global grid cell boundaries for given model component and
!! mosaic tile number.
subroutine get_grid_cell_vertices_1D(component, tile, glonb, glatb)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer,          intent(in) :: tile !< Tile number
  real,          intent(inout) :: glonb(:),glatb(:) !< Grid cell vertices

  integer                      :: nlon, nlat
  integer                      :: start(4), nread(4)
  real, allocatable            :: tmp(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)
  character(len=MAX_FILE)      :: tilefile
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
  case(VERSION_0)
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
  case(VERSION_1)
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
  case(VERSION_2, VERSION_3)
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
end subroutine get_grid_cell_vertices_1D

!> @brief returns cell vertices for the specified model component and mosaic tile number
subroutine get_grid_cell_vertices_2D(component, tile, lonb, latb, domain)
  character(len=*),         intent(in) :: component !< Component model (atm, lnd, ocn)
  integer,                  intent(in) :: tile !< Tile number
  real,                  intent(inout) :: lonb(:,:),latb(:,:) !< Cell vertices
  type(domain2d), optional, intent(in) :: domain !< Domain

  ! local vars
  integer :: nlon, nlat
  integer :: i,j
  real, allocatable :: buffer(:), tmp(:,:), x_vert_t(:,:,:), y_vert_t(:,:,:)
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0,j0 ! offsets for coordinates
  integer :: isg, jsg
  integer :: start(4), nread(4)
  character(len=MAX_FILE)      :: tilefile
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

  select case(grid_version)
  case(VERSION_0)
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
  case(VERSION_1)
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
  case(VERSION_2, VERSION_3)
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
  end select
end subroutine get_grid_cell_vertices_2D

!> @brief returns cell vertices for the specified model component and mosaic tile number for
!! an unstructured domain
subroutine get_grid_cell_vertices_UG(component, tile, lonb, latb, SG_domain, UG_domain)
  character(len=*),         intent(in) :: component !< Component model (atm, lnd, ocn)
  integer,                  intent(in) :: tile !< Tile number
  real,                  intent(inout) :: lonb(:,:),latb(:,:) ! The second dimension is 4
  type(domain2d)  ,   intent(in)       :: SG_domain !< Structured domain
  type(domainUG)  ,   intent(in)       :: UG_domain !< Unstructured domain
  integer :: is, ie, js, je, i, j
  real, allocatable :: SG_lonb(:,:), SG_latb(:,:), tmp(:,:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_lonb(is:ie+1, js:je+1))
  allocate(SG_latb(is:ie+1, js:je+1))
  allocate(tmp(is:ie,js:je,4))
  call get_grid_cell_vertices_2D(component, tile, SG_lonb, SG_latb, SG_domain)
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
end subroutine get_grid_cell_vertices_UG

!> @brief returns grid cell centers given model component and mosaic tile number
subroutine get_grid_cell_centers_1D(component, tile, glon, glat)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real, intent(inout) :: glon(:),glat(:) !< Grid cell centers

  integer                      :: nlon, nlat
  integer                      :: start(4), nread(4)
  real, allocatable            :: tmp(:,:)
  character(len=MAX_FILE)      :: tilefile
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
  case(VERSION_0)
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
  case(VERSION_1)
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
  case(VERSION_2, VERSION_3)
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
end subroutine get_grid_cell_centers_1D

!> @brief returns grid cell centers given model component and mosaic tile number
subroutine get_grid_cell_centers_2D(component, tile, lon, lat, domain)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real, intent(inout) :: lon(:,:),lat(:,:) !< Grid cell centers
  type(domain2d), intent(in), optional :: domain !< Domain
  ! local vars
  integer :: nlon, nlat
  integer :: i,j
  real, allocatable :: buffer(:),tmp(:,:)
  integer :: is,ie,js,je ! boundaries of our domain
  integer :: i0,j0 ! offsets for coordinates
  integer :: isg, jsg
  integer :: start(4), nread(4)
  character(len=MAX_FILE)      :: tilefile
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
  case(VERSION_0)
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
  case(VERSION_1)
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
  case(VERSION_2, VERSION_3) ! mosaic grid file
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
end subroutine get_grid_cell_centers_2D

!> @brief returns grid cell centers given model component and mosaic tile number
!! for unstructured domain
subroutine get_grid_cell_centers_UG(component, tile, lon, lat, SG_domain, UG_domain)
  character(len=*), intent(in) :: component !< Component model (atm, lnd, ocn)
  integer, intent(in) :: tile !< Tile number
  real, intent(inout) :: lon(:),lat(:) !< Grid cell centers
  type(domain2d)  ,   intent(in) :: SG_domain !< Structured domain
  type(domainUG)  ,   intent(in) :: UG_domain !< Unstructured domain
  integer :: is, ie, js, je
  real, allocatable :: SG_lon(:,:), SG_lat(:,:)

  call mpp_get_compute_domain(SG_domain, is, ie, js, je)
  allocate(SG_lon(is:ie, js:je))
  allocate(SG_lat(is:ie, js:je))
  call get_grid_cell_centers_2D(component, tile, SG_lon, SG_lat, SG_domain)
  call mpp_pass_SG_to_UG(UG_domain, SG_lon, lon)
  call mpp_pass_SG_to_UG(UG_domain, SG_lat, lat)
  deallocate(SG_lon, SG_lat)
end subroutine get_grid_cell_centers_UG

!> @brief given a model component, a layout, and (optionally) a halo size, returns a
!! domain for current processor
subroutine define_cube_mosaic(component, domain, layout, halo, maskmap)
  character(len=*) , intent(in)    :: component !< Component model (atm, lnd, ocn)
  type(domain2d)   , intent(inout) :: domain !< Domain
  integer          , intent(in)    :: layout(2) !< Layout
  integer, optional, intent(in)    :: halo !< Halo
  logical, optional, intent(in)    :: maskmap(:,:,:) !< Maskmap

  ! ---- local vars
  integer :: ntiles     ! number of tiles
  integer :: ncontacts  ! number of contacts between mosaic tiles
  integer :: n
  integer :: ng, pe_pos, npes         ! halo size
  integer, allocatable :: nlon(:), nlat(:), global_indices(:,:)
  integer, allocatable :: pe_start(:), pe_end(:), layout_2d(:,:)
  integer, allocatable :: tile1(:),tile2(:)
  integer, allocatable :: is1(:),ie1(:),js1(:),je1(:)
  integer, allocatable :: is2(:),ie2(:),js2(:),je2(:)

  call get_grid_ntiles(component,ntiles)
  allocate(nlon(ntiles), nlat(ntiles))
  allocate(global_indices(4,ntiles))
  allocate(pe_start(ntiles),pe_end(ntiles))
  allocate(layout_2d(2,ntiles))
  call get_grid_size(component,nlon,nlat)

  pe_pos = mpp_root_pe()
  do n = 1, ntiles
     global_indices(:,n) = (/ 1, nlon(n), 1, nlat(n) /)
     layout_2d     (:,n) = layout
     if(present(maskmap)) then
        npes = count(maskmap(:,:,n))
     else
        npes = layout(1)*layout(2)
     endif
     pe_start(n) = pe_pos
     pe_end  (n) = pe_pos + npes - 1
     pe_pos      = pe_end(n) + 1
  enddo

  ! get the contact information from mosaic file
  ncontacts = get_mosaic_ncontacts(mosaic_fileobj(get_component_number(trim(component))))
  allocate(tile1(ncontacts),tile2(ncontacts))
  allocate(is1(ncontacts),ie1(ncontacts),js1(ncontacts),je1(ncontacts))
  allocate(is2(ncontacts),ie2(ncontacts),js2(ncontacts),je2(ncontacts))
  call get_mosaic_contact(mosaic_fileobj(get_component_number(trim(component))), tile1, tile2, &
       is1, ie1, js1, je1, is2, ie2, js2, je2)

  ng = 0
  if(present(halo)) ng = halo
  ! create the domain2d variable
  call mpp_define_mosaic ( global_indices, layout_2d, domain, &
       ntiles, ncontacts, tile1, tile2,                  &
       is1, ie1, js1, je1, &
       is2, ie2, js2, je2, &
       pe_start=pe_start, pe_end=pe_end, symmetry=.true.,  &
       shalo = ng, nhalo = ng, whalo = ng, ehalo = ng,     &
       maskmap = maskmap,                                  &
       name = trim(component)//'Cubic-Sphere Grid' )

  deallocate(nlon,nlat,global_indices,pe_start,pe_end,layout_2d)
  deallocate(tile1,tile2)
  deallocate(is1,ie1,js1,je1)
  deallocate(is2,ie2,js2,je2)
end subroutine define_cube_mosaic

end module grid2_mod
!> @}
! close documentation grouping
