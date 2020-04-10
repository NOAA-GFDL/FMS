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
module fms2_io_namelist_mod
!> @author Tom Robinson
!> @email gfdl.climate.model.info@noaa.gov
!> @description This module contains the namelist for the fms2_io routines.  The variables are all
!! module variables, so they can be used elsewhere.  The routine fms2_io_init should be called by
!! fms_init

use mpp_mod, only: input_nml_file 
use fms_io_utils_mod, only: error
implicit none


public :: fms2_io_init
public :: fms2_ncblksz

integer :: fms2_ncblksz = 1*1024*1024 !< Used as chunksize argument in netcdf file creation calls.

!< Namelist variables
integer :: ncblksz = 1*1024*1024 !< User defined chunksize argument in netcdf file creation calls.
                                 !! Replaces setting the NC_CHKSZ environment variable.

namelist / fms2_io_nml / &
                      ncblksz

contains

!> @brief Reads the fms2_io_nml
subroutine fms2_io_init ()

integer :: mystat

READ (input_nml_file, NML=fms2_io_nml, IOSTAT=mystat)

  fms2_ncblksz = ncblksz

end subroutine fms2_io_init

end module fms2_io_namelist_mod
