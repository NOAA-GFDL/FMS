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

!> @brief  This programs tests the functionality that writes the boundary
!conditions restarts
program test_bc_restart

use   mpp_mod,         only : mpp_init, mpp_exit, mpp_pe, mpp_root_pe, mpp_sync
use   fms2_io_mod,     only : FmsNetcdfFile_t, fms2_io_init, open_file, register_restart_field, &
                              read_restart_bc, write_restart_bc, close_file
use   mpp_domains_mod, only : mpp_get_global_domain, mpp_get_data_domain, mpp_get_compute_domain, &
                              mpp_define_domains, mpp_get_global_domain, domain2d, CORNER

implicit none

type atm_type
   type(FmsNetcdfFile_t)              :: fileobj_sw       !< fms2io netcdf file obj for south west bc
   type(FmsNetcdfFile_t)              :: fileobj_ne       !< fms2io netcdf file obj for north east bc
   logical                            :: BCfile_sw_open   !< flag indicating if the sourth west file is
                                                          !! opened
   logical                            :: BCfile_ne_open   !< flag indicating if the sourth west file is
                                                          !! opened
   type(domain2d)                     :: Domain           !< Domain with halos
   real, allocatable, dimension(:,:)  :: var2d            !< 2d variable data
   real, allocatable, dimension(:,:,:):: var3d            !< 3d variable data
end type

integer, dimension(2)                 :: layout = (/4,4/) !< Domain layout
integer                               :: nlon             !< Number of points in x axis
integer                               :: nlat             !< Number of points in y axis
integer                               :: isd, jsd         !< Starting x/y index (data_domain)
integer                               :: ied, jed         !< Ending x/y index (data_domain)
integer, allocatable, dimension(:)    :: all_pelist       !< List of pelist associated with the test
integer                               :: n                !< No description
type(atm_type)                        :: atm              !< No description

call mpp_init
call fms2_io_init

nlon = 144
nlat = 144

call mpp_define_domains( (/1,nlon,1,nlat/), layout, atm%Domain, xhalo=3, yhalo=3, symmetry=.true., name='test_bc_restart')
call mpp_get_data_domain(atm%domain, isd, ied, jsd, jed )

allocate(atm%var2d(isd:ied,jsd:jed))
allocate(atm%var3d(isd:ied,jsd:jed,5))

atm%var2d = real(mpp_pe()+100)
atm%var3d = real(mpp_pe())

!< Create array of current pelist
allocate  (all_pelist(layout(1)*layout(2)))
do n = 1, layout(1)*layout(2)
   all_pelist(n) = n - 1
end do

atm%BCfile_sw_open = .false.
atm%BCfile_ne_open = .false.

!< Try to write a BC restart file:
atm%BCfile_sw_open = open_file(atm%fileobj_sw, "BCfile_sw.nc", "overwrite", is_restart=.true., pelist=all_pelist)
atm%BCfile_ne_open = open_file(atm%fileobj_ne, "BCfile_ne.nc", "overwrite", is_restart=.true., pelist=all_pelist)

call register_bcs(atm, atm%fileobj_ne, atm%fileobj_sw, "sst", layout)

if (atm%BCfile_sw_open) then
    call write_restart_bc(atm%fileobj_sw)
    call close_file(atm%fileobj_sw)
endif

if (atm%BCfile_ne_open) then
    call write_restart_bc(atm%fileobj_ne)
    call close_file(atm%fileobj_ne)
endif

call mpp_sync()

atm%var2d = real(999.)
atm%var3d = real(999.)

atm%BCfile_sw_open = .false.
atm%BCfile_ne_open = .false.

!< Try to read the boundary condition restart back
atm%BCfile_sw_open = open_file(atm%fileobj_sw, "BCfile_sw.nc", "read", is_restart=.true., pelist=all_pelist)
atm%BCfile_ne_open = open_file(atm%fileobj_ne, "BCfile_ne.nc", "read", is_restart=.true., pelist=all_pelist)

call register_bcs(atm, atm%fileobj_ne, atm%fileobj_sw, "sst", layout)

if (atm%BCfile_sw_open) then
    call read_restart_bc(atm%fileobj_sw)
    call close_file(atm%fileobj_sw)
endif

if (atm%BCfile_ne_open) then
    call read_restart_bc(atm%fileobj_ne)
    call close_file(atm%fileobj_ne)
endif

call mpp_sync()
call mpp_exit()

contains

!> @brief Dummy wrapper for the register_restart_variable calls for boundary
!! condition restarts
subroutine register_bcs(atm, fileobj_ne, fileobj_sw, var_name, layout, istag, jstag)
  type(atm_type), intent(inout) :: atm
  type(FmsNetcdfFile_t), intent(inout) :: fileobj_sw  !< fms2io netcdf file obj for south west bc
  type(FmsNetcdfFile_t), intent(inout) :: fileobj_ne  !< fms2io netcdf file obj for north east bc
  character(len=*),         intent(in) :: var_name    !< Name of the variable
  integer, dimension(2),    intent(in) :: layout      !< Domain layout

  integer,                  intent(in), optional :: istag, jstag

  integer                              :: isd, jsd         !< Starting x/y index (data_domain)
  integer                              :: ied, jed         !< Ending x/y index (data_domain)
  integer                              :: is, js           !< Starting x/y index (compute_domain)
  integer                              :: ie, je           !< Ending x/y index (compute_domain)
  integer                              :: i_stag, j_stag   !< Extra x/y?
  integer                              :: npx, npy         !< Number of points in x/y (global_domain)
  integer                              :: x_halo, y_halo   !< Number of halos in x and y
  integer                              :: x_halo_ns        !< Halo for north and south variables
  integer, allocatable, dimension(:)   :: x1_pelist        !< Pelist corresponding to the bottom region
                                                           !! of the domain
  integer, allocatable, dimension(:)   :: y1_pelist        !< Pelist corresponding to the right region
                                                           !! of the domain
  integer, allocatable, dimension(:)   :: x2_pelist        !< Pelist corresponding to the top region
                                                           !! of the domain
  integer, allocatable, dimension(:)   :: y2_pelist        !< Pelist corresponding to the left region
                                                           !! of the domain
  integer                              :: n                !< No description
  integer, dimension(3)                :: global_size      !< Size of the domain
  integer, dimension(4)                :: indices          !< Starting/Ending indices of the current
  logical                              :: is_root_pe       !< Flag indicating if this is the root pe

  !< .----.----.----.----.
  !< |PE 0|PE 1|PE 2|PE 3| <- x2_pelist
  !< .----.----.----.----.
  !< |PE 4|PE 5|PE 6|PE 7|
  !< .----.----.----.----.
  !< |PE 8|PE 9|PE10|PE11|
  !< .----.----.----.----.
  !< |PE12|PE13|PE14|PE15| <- x1_pelist
  !< .----.----.----.----.
  !<    ^y2_pelist     ^y1_pelist

  i_stag = 0
  j_stag = 0
  if (present(istag)) i_stag = i_stag
  if (present(jstag)) j_stag = j_stag

  call mpp_get_global_domain(atm%domain, xsize = npx, ysize = npy, position=CORNER )
  call mpp_get_data_domain(atm%domain, isd, ied, jsd, jed )
  call mpp_get_compute_domain(atm%domain, is, ie, js, je )

  !< Defaults
  x_halo = is-isd
  y_halo = js-jsd
  i_stag = 0
  j_stag = 0

  allocate (x1_pelist(layout(1)))
  allocate (y1_pelist(layout(2)))
  allocate (x2_pelist(layout(1)))
  allocate (y2_pelist(layout(2)))

  !< Define west and east pelist
  do n = 1,layout(2)
     y1_pelist(n)=mpp_root_pe()+layout(1)*n-1
     y2_pelist(n)=mpp_root_pe()+layout(1)*(n-1)
  enddo

  !< Define south and north pelist
  do n = 1,layout(1)
     x1_pelist(n)=mpp_root_pe()+layout(1)*(layout(2)-1)+(n-1)
     x2_pelist(n)=mpp_root_pe()+(n-1)
  enddo

  !< EAST & WEST
  !< Set defaults for west/east halo regions
  indices(1) = 1
  indices(2) = x_halo
  indices(3) = jsd
  indices(4) = jed + j_stag
  global_size(1) = x_halo
  global_size(2) = npy-1+2*y_halo + j_stag

  global_size(3) = size(atm%var3d, 3)

  !< Define west root_pe
  is_root_pe = .FALSE.
  if (is.eq.1 .and. js.eq.1) is_root_pe = .TRUE.

  if (atm%BCfile_sw_open) call register_restart_field(fileobj_sw, trim(var_name)//'_west', atm%var2d, &
                                                indices, global_size(1:2), y2_pelist, &
                                                is_root_pe, jshift=y_halo)
  if (atm%BCfile_sw_open) call register_restart_field(fileobj_sw, trim(var_name)//'3d_west', atm%var3d, &
                                                indices, global_size(1:3), y2_pelist, &
                                                is_root_pe, jshift=y_halo)

  !< Define east root_pe
  is_root_pe = .FALSE.
  if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.

  !< Reset indices for prognostic variables in the east halo
  indices(1) = ied-x_halo+1+i_stag
  indices(2) = ied+i_stag

  if (atm%BCfile_ne_open) call register_restart_field(fileobj_ne, trim(var_name)//'east', atm%var2d, &
                                                indices, global_size(1:2), y1_pelist, &
                                                is_root_pe, jshift=y_halo, &
                                                x_halo=(size(atm%var2d,1)-x_halo), ishift=-(ie+i_stag))
  if (atm%BCfile_ne_open) call register_restart_field(fileobj_ne, trim(var_name)//'3d_east', atm%var3d, &
                                                indices, global_size(1:3), y1_pelist, &
                                                is_root_pe, jshift=y_halo, &
                                                x_halo=(size(atm%var3d,1)-x_halo), ishift=-(ie+i_stag))
  !< NORTH & SOUTH
  !< set defaults for north/south halo regions
  indices(1) = isd
  indices(2) = ied+i_stag
  indices(3) = 1
  indices(4) = y_halo
  global_size(1) = npx-1+i_stag
  global_size(2) = y_halo

  !< Modify starts and ends for certain pes
  if (is.eq.1)     indices(1) = is
  if (ie.eq.npx-1) indices(2) = ie+i_stag
  x_halo_ns = 0
  if (is.eq.1) x_halo_ns=x_halo

  !define south root_pe
  is_root_pe = .FALSE.
  if (is.eq.1 .and. js.eq.1) is_root_pe = .TRUE.

  if (atm%BCfile_sw_open) call register_restart_field(fileobj_sw, trim(var_name)//'south', atm%var2d, &
                                                indices, global_size(1:2), x2_pelist, &
                                                is_root_pe, x_halo=x_halo_ns)
  if (atm%BCfile_sw_open) call register_restart_field(fileobj_sw, trim(var_name)//'3d_south', atm%var3d, &
                                                indices, global_size(1:3), x2_pelist, &
                                                is_root_pe, x_halo=x_halo_ns)

  !< Define north root_pe
  is_root_pe = .FALSE.
  if (ie.eq.npx-1 .and. je.eq.npy-1) is_root_pe = .TRUE.

  !< Reset indices for prognostic variables in the north halo
  indices(3) = jed-y_halo+1+j_stag
  indices(4) = jed+j_stag

  if (atm%BCfile_ne_open) call register_restart_field(fileobj_ne, trim(var_name)//'north', atm%var2d, &
                                                indices, global_size(1:2), x1_pelist, &
                                                is_root_pe, x_halo=x_halo_ns, &
                                                y_halo=(size(atm%var2d,2)-y_halo), jshift=-(je+j_stag))
  if (atm%BCfile_ne_open) call register_restart_field(fileobj_ne, trim(var_name)//'3d_north', atm%var3d, &
                                                indices, global_size(1:3), x1_pelist, &
                                                is_root_pe, x_halo=x_halo_ns, &
                                                y_halo=(size(atm%var3d,2)-y_halo), jshift=-(je+j_stag))

end subroutine register_bcs

end program test_bc_restart
