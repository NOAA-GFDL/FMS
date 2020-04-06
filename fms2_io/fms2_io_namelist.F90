module fms2_io_namelist_mod
!> @author Tom Robinson
!> @email thomas.robinson@noaa.gov
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
