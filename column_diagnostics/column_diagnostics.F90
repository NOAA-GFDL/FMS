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
!> @defgroup column_diagnostics_mod column_diagnostics_mod
!> @ingroup column_diagnostics
!! @brief Module to locate and mark desired diagnostic columns

!> @addtogroup column_diagnostics_mod
!> @{
module column_diagnostics_mod

use fms_mod,                only:  fms_init, mpp_pe, mpp_root_pe, &
                                   mpp_npes, check_nml_error, &
                                   error_mesg, FATAL, NOTE, WARNING, &
                                   stdlog, write_version_number
use time_manager_mod,       only:  time_manager_init, month_name, &
                                   get_date, time_type
use constants_mod,          only:  constants_init, PI, RADIAN
use mpp_mod,                only:  input_nml_file
use platform_mod,           only:  r4_kind, r8_kind, FMS_FILE_LEN
!-------------------------------------------------------------------

implicit none
private

!---------------------------------------------------------------------
!      module to locate and mark desired diagnostic columns
!
!
!--------------------------------------------------------------------




!---------------------------------------------------------------------
!----------- ****** VERSION NUMBER ******* ---------------------------


! Include variable "version" to be written to log file.
#include<file_version.h>



!---------------------------------------------------------------------
!-------  interfaces --------

public    column_diagnostics_init,  &
          initialize_diagnostic_columns,  &
          column_diagnostics_header,   &
          close_column_diagnostics_units


interface initialize_diagnostic_columns
  module procedure initialize_diagnostic_columns_r4
  module procedure initialize_diagnostic_columns_r8
end interface initialize_diagnostic_columns

interface column_diagnostics_header
  module procedure column_diagnostics_header_r4
  module procedure column_diagnostics_header_r8
end interface column_diagnostics_header

!private

!--------------------------------------------------------------------
!----    namelist -----

real(kind=r8_kind) :: crit_xdistance = 4.0_r8_kind !< model grid points must be within crit_xdistance in
                                      !! longitude of the requested diagnostics point
                                      !! coordinates in order to be flagged as the desired
                                      !! point
                                      !! [ degrees ]
real(kind=r8_kind) :: crit_ydistance = 4.0_r8_kind !< model grid points must be within crit_ydistance in
                                      !! latitude of the requested diagnostics point
                                      !! coordinates in order to be flagged as the desired
                                      !! point
                                      !! [ degrees ]

namelist / column_diagnostics_nml /                   &
                                      crit_xdistance, &
                                      crit_ydistance

!--------------------------------------------------------------------
!-------- public data  -----


!--------------------------------------------------------------------
!------ private data ------


logical    :: module_is_initialized = .false.

!-------------------------------------------------------------------
!-------------------------------------------------------------------



                        contains



!####################################################################

!> @brief Initialization routine for column_diagnostics_mod.
!!
!> Reads namelist and writes to log.
subroutine column_diagnostics_init

!--------------------------------------------------------------------
!    column_diagnostics_init is the constructor for
!    column_diagnostics_mod.
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    local variables:
!
      integer    :: iunit !< unit number for nml file
      integer    :: ierr !< error return flag
      integer    :: io   !< error return code

!--------------------------------------------------------------------
!   local variables:
!
!       unit       unit number for nml file
!       ierr       error return flag
!       io         error return code
!
!---------------------------------------------------------------------

!--------------------------------------------------------------------
!    if routine has already been executed, return.
!--------------------------------------------------------------------
      if (module_is_initialized) return

!---------------------------------------------------------------------
!    verify that all modules used by this module have been initialized.
!----------------------------------------------------------------------
      call fms_init
      call time_manager_init
      call constants_init

!---------------------------------------------------------------------
!    read namelist.
!---------------------------------------------------------------------
      read (input_nml_file, column_diagnostics_nml, iostat=io)
      ierr = check_nml_error (io, 'column_diagnostics_nml')
!---------------------------------------------------------------------
!    write version number and namelist to logfile.
!---------------------------------------------------------------------
      call write_version_number("COLUMN_DIAGNOSTICS_MOD", version)
      if (mpp_pe() == mpp_root_pe())    then
                    iunit = stdlog()
                    write (iunit, nml=column_diagnostics_nml)
      endif
!--------------------------------------------------------------------
      module_is_initialized = .true.


end subroutine column_diagnostics_init


!######################################################################
!> @brief close_column_diagnostics_units closes any open column_diagnostics
!!    files associated with the calling module.
subroutine close_column_diagnostics_units (diag_units)

!---------------------------------------------------------------------
!    close_column_diagnostics_units closes any open column_diagnostics
!    files associated with the calling module.
!----------------------------------------------------------------------

!----------------------------------------------------------------------
integer, dimension(:), intent(in)  :: diag_units !< array of column diagnostic unit numbers
!----------------------------------------------------------------------

!--------------------------------------------------------------------
!    intent(in) variable:
!
!      diag_units    array of column diagnostic unit numbers
!
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!    local variable

      integer   :: nn    !< do loop index
      integer   :: io
!--------------------------------------------------------------------
!    close the unit associated with each diagnostic column.
!--------------------------------------------------------------------
      do nn=1, size(diag_units(:))
        if (diag_units(nn) /= -1) then
          close(diag_units(nn), iostat=io )
          if(io/=0) call error_mesg('column_diagnostics_mod', 'Error in closing file ', FATAL)
        endif
      end do

!---------------------------------------------------------------------


end subroutine close_column_diagnostics_units


!#####################################################################

#include "column_diagnostics_r4.fh"
#include "column_diagnostics_r8.fh"


               end module column_diagnostics_mod
!@}
! close documentation grouping
