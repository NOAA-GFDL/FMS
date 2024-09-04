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
!> @defgroup get_grid_version_mod get_grid_version_mod
!> @ingroup data_override
!> @brief get_grid implementations and helper routines for @ref data_override_mod

!> @addtogroup get_grid_version_mod
!> @{
module get_grid_version_mod
use constants_mod, only: DEG_TO_RAD
use platform_mod, only: r4_kind, r8_kind, FMS_PATH_LEN
use mpp_mod, only : mpp_error,FATAL,NOTE, mpp_min, mpp_max
use mpp_domains_mod, only : domain2d, operator(.NE.),operator(.EQ.)
use mpp_domains_mod, only : mpp_get_global_domain, mpp_get_data_domain
use fms2_io_mod,     only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t, open_file, close_file, &
                            variable_exists, read_data, get_variable_size, get_variable_num_dimensions
use mosaic2_mod,     only : get_mosaic_tile_grid

implicit none

interface get_grid_version_1
  module procedure get_grid_version_1_r4
  module procedure get_grid_version_1_r8
end interface get_grid_version_1

interface get_grid_version_2
  module procedure get_grid_version_2_r4
  module procedure get_grid_version_2_r8
end interface get_grid_version_2

contains

!> Get lon and lat of three model (target) grids from grid_spec.nc
subroutine check_grid_sizes(domain_name, Domain, nlon, nlat)
character(len=12), intent(in) :: domain_name
type (domain2d),   intent(in) :: Domain
integer,           intent(in) :: nlon, nlat

character(len=184) :: error_message
integer            :: xsize, ysize

call mpp_get_global_domain(Domain, xsize=xsize, ysize=ysize)
if(nlon .NE. xsize .OR. nlat .NE. ysize) then
  error_message = 'Error in data_override_init. Size of grid as specified by '// &
                  '             does not conform to that specified by grid_spec.nc.'// &
                  '  From             :     by      From grid_spec.nc:     by    '
  error_message( 59: 70) = domain_name
  error_message(130:141) = domain_name
  write(error_message(143:146),'(i4)') xsize
  write(error_message(150:153),'(i4)') ysize
  write(error_message(174:177),'(i4)') nlon
  write(error_message(181:184),'(i4)') nlat
  call mpp_error(FATAL,error_message)
endif
end subroutine check_grid_sizes

#include "get_grid_version_r4.fh"
#include "get_grid_version_r8.fh"

end module get_grid_version_mod
!> @}
! close documentation grouping
