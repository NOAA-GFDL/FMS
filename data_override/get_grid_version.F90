!***********************************************************************
!*                             Apache License 2.0
!*
!* This file is part of the GFDL Flexible Modeling System (FMS).
!*
!* Licensed under the Apache License, Version 2.0 (the "License");
!* you may not use this file except in compliance with the License.
!* You may obtain a copy of the License at
!*
!*     http://www.apache.org/licenses/LICENSE-2.0
!*
!* FMS is distributed in the hope that it will be useful, but WITHOUT
!* WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied;
!* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!* PARTICULAR PURPOSE. See the License for the specific language
!* governing permissions and limitations under the License.
!***********************************************************************
!> @defgroup get_grid_version_mod get_grid_version_mod
!> @ingroup data_override
!> @brief This module contains helper routines for @ref data_override_mod,
!! and is not intended to be used directly. It provides interfaces for reading the global grid dimensions
!! by reading a grid spec or mosaic file, and also a routine for checking that the read-in sizes
!! are valid. The exact variables read to get this information are controlled by the model domain
!! (ie. atm, ocn, land), which is passed in as a string.

!> @addtogroup get_grid_version_mod
!> @{
module get_grid_version_mod
use constants_mod, only: DEG_TO_RAD
use platform_mod, only: r4_kind, r8_kind, FMS_PATH_LEN
use mpp_mod, only : mpp_error,FATAL,NOTE, mpp_min, mpp_max
use mpp_domains_mod, only : domain2d, operator(.NE.),operator(.EQ.)
use mpp_domains_mod, only : mpp_get_compute_domain, mpp_get_global_domain
use fms2_io_mod,     only : FmsNetcdfDomainFile_t, FmsNetcdfFile_t, open_file, close_file, &
                            variable_exists, read_data, get_variable_size, get_variable_num_dimensions
use mosaic2_mod,     only : get_mosaic_tile_grid

implicit none

!> @brief get_grid_version_1 is intended to read in what is typically named "grid_spec.nc", typically created by
!! make_mosaic/make_coupler_mosaic in fre-nctools
interface get_grid_version_1
  module procedure get_grid_version_1_r4
  module procedure get_grid_version_1_r8
end interface get_grid_version_1

!> @brief get_grid_version_2 reads grid coordinate data from an already-opened grid_spec.nc or mosaic
!! file object, typically created by make_mosaic/make_coupler_mosaic in fre-nctools
interface get_grid_version_2
  module procedure get_grid_version_2_r4
  module procedure get_grid_version_2_r8
end interface get_grid_version_2

contains

!> @brief Validates that the domain size matches the grid dimensions from the grid file, calling FATAL if not.
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
