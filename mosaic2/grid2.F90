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

!> @addtogroup grid2_mod
!> @{
module grid2_mod

use mpp_mod, only : mpp_root_pe, mpp_error, uppercase, lowercase, FATAL, NOTE
use constants_mod, only : PI, radius
use fms2_io_mod, only : get_global_attribute, read_data, global_att_exists, &
                        variable_exists, file_exists,  open_file, close_file, get_variable_size, &
                        FmsNetcdfFile_t
use fms_string_utils_mod, only: string
use mosaic2_mod, only : get_mosaic_ntiles, get_mosaic_xgrid_size, get_mosaic_grid_sizes, &
     get_mosaic_xgrid, calc_mosaic_grid_area, calc_mosaic_grid_great_circle_area

! the following two use statement are only needed for define_cube_mosaic
use mpp_domains_mod, only : domain2d, mpp_define_mosaic, mpp_get_compute_domain, &
                            mpp_get_global_domain, domainUG, mpp_pass_SG_to_UG
use mosaic2_mod, only : get_mosaic_ncontacts, get_mosaic_contact
use platform_mod

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
   module procedure get_grid_cell_vertices_1D_r4
   module procedure get_grid_cell_vertices_1D_r8
   module procedure get_grid_cell_vertices_2D_r4
   module procedure get_grid_cell_vertices_2D_r8
   module procedure get_grid_cell_vertices_UG_r4
   module procedure get_grid_cell_vertices_UG_r8
end interface

!> Gets grid cell centers
!> @ingroup grid2_mod
interface get_grid_cell_centers
   module procedure get_grid_cell_centers_1D_r4
   module procedure get_grid_cell_centers_1D_r8
   module procedure get_grid_cell_centers_2D_r4
   module procedure get_grid_cell_centers_2D_r8
   module procedure get_grid_cell_centers_UG_r4
   module procedure get_grid_cell_centers_UG_r8
end interface

!> Finds area of a grid cell
!> @ingroup grid2_mod
interface get_grid_cell_area
   module procedure get_grid_cell_area_SG_r4
   module procedure get_grid_cell_area_SG_r8
   module procedure get_grid_cell_area_UG_r4
   module procedure get_grid_cell_area_UG_r8
end interface get_grid_cell_area

!> Gets the area of a given component per grid cell
!> @ingroup grid2_mod
interface get_grid_comp_area
   module procedure get_grid_comp_area_SG_r4
   module procedure get_grid_comp_area_SG_r8
   module procedure get_grid_comp_area_UG_r4
   module procedure get_grid_comp_area_UG_r8
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
     VERSION_GEOLON_T        = 0,   & !< indicates gelon_t variable is present in grid_file
     VERSION_X_T             = 1,   & !< indicates x_t variable is present in grid_file
     VERSION_OCN_MOSAIC_FILE = 2,   & !< indicates ocn_mosaic_file variable is present in grid_file
     VERSION_GRIDFILES       = 3      !< indicates gridfiles variable is present in grid_file

integer, parameter :: BUFSIZE = 1048576  !< This is used to control memory usage in get_grid_comp_area
                                         !! We may change this to a namelist variable is needed.

! ==== module variables ======================================================
integer :: grid_version = -1 !< Value to indicate what type of grid file is being read,
                             !! based on which variables are present
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
   if (grid_version == VERSION_OCN_MOSAIC_FILE) call open_component_mosaics
   if (grid_version == VERSION_GRIDFILES) call assign_component_mosaics
   module_is_initialized = .TRUE.
end subroutine grid_init

!> @brief Shutdown the grid2 module
subroutine grid_end
   if (grid_spec_exists) then
       if (grid_version == VERSION_OCN_MOSAIC_FILE) call close_component_mosaics
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
         call mpp_error(FATAL, module_name//&
                   '/get_great_circle_algorithm value of global attribute "great_circle_algorthm" in file'// &
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

  character(len=FMS_PATH_LEN) :: mosaicfilename
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
  character(len=FMS_PATH_LEN) :: read_file_name
  character(len=FMS_PATH_LEN), dimension(:), allocatable :: file_names

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
       get_grid_version = VERSION_GEOLON_T
    else if(variable_exists(fileobj, 'x_T')) then
       get_grid_version = VERSION_X_T
    else if(variable_exists(fileobj, 'ocn_mosaic_file') ) then
       get_grid_version = VERSION_OCN_MOSAIC_FILE
    else if(variable_exists(fileobj, 'gridfiles') ) then
       get_grid_version = VERSION_GRIDFILES
    else
       call mpp_error(FATAL, module_name//'/get_grid_version Can''t determine the version of the grid spec:'// &
                  & ' none of "x_T", "geolon_t", or "ocn_mosaic_file" exist in file "'//trim(grid_file)//'"')
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
  case(VERSION_GEOLON_T,VERSION_X_T)
     ntiles = 1
  case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES)
     ntiles = get_mosaic_ntiles(mosaic_fileobj(get_component_number(trim(component))))
  end select
end subroutine get_grid_ntiles

!> @brief returns size of the grid for each of the tiles
subroutine get_grid_size_for_all_tiles(component,nx,ny)
  character(len=*)     :: component !< Component model (atm, lnd, ocn)
  integer, intent(inout) :: nx(:),ny(:) !< Grid size in x and y

  ! local vars
  integer :: siz(2) ! for the size of external fields
  character(len=MAX_NAME) :: varname1

  varname1 = 'AREA_'//trim(uppercase(component))

  select case (grid_version)
  case(VERSION_GEOLON_T,VERSION_X_T)
     if (.not. grid_spec_exists) then
       call mpp_error(FATAL, 'grid2_mod(get_grid_size_for_all_tiles): grid_spec does not exist')
     end if
     call get_variable_size(gridfileobj, varname1, siz)
     nx(1) = siz(1); ny(1)=siz(2)
  case(VERSION_OCN_MOSAIC_FILE, VERSION_GRIDFILES) ! mosaic file
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

#include "grid2_r4.fh"
#include "grid2_r8.fh"

end module grid2_mod
!> @}
! close documentation grouping
