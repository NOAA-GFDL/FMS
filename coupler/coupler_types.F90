module coupler_types_mod  !{
!-----------------------------------------------------------------------
!                   GNU General Public License                        
! This file is a part of MOM.                                                                 
!                                                                      
! MOM is free software; you can redistribute it and/or modify it and  
! are expected to follow the terms of the GNU General Public License  
! as published by the Free Software Foundation; either version 2 of   
! the License, or (at your option) any later version.                 
!                                                                      
! MOM is distributed in the hope that it will be useful, but WITHOUT    
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY  
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public    
! License for more details.                                           
!                                                                      
! For the full text of the GNU General Public License,                
! write to: Free Software Foundation, Inc.,                           
!           675 Mass Ave, Cambridge, MA 02139, USA.                   
! or see:   http://www.gnu.org/licenses/gpl.html                      
!-----------------------------------------------------------------------
!
!<CONTACT EMAIL="Richard.Slater@noaa.gov">
! Richard D. Slater 
!</CONTACT>
!
! <REVIEWER EMAIL="John.Dunne@noaa.gov">
! John Dunne
! </REVIEWER>
!
!<OVERVIEW>
! This module contains type declarations for the coupler.
!</OVERVIEW>
!
!<DESCRIPTION>
! This module contains type declarations for the coupler.
!</DESCRIPTION>
!

!
! /coupler_mod/
!              types/
!                    air_sea_gas_flux_generic/
!                         implementation/
!                                        ocmip2/
!                                                    num_parameters = 2
!                         num_flags = 0
!                         use_atm_pressure = t
!                         use_10m_wind_speed = t
!                         pass_through_ice = f
!                         atm/
!                             name/
!                                  pcair, u10, psurf
!                             long_name/
!                                       'Atmospheric concentration'
!                                       'Wind speed at 10 m'
!                                       'Surface atmospheric pressure'
!                             units/
!                                   'mol/mol', 'm/s', 'Pa'
!                         ice/
!                             name/
!                                  alpha, csurf, sc_no
!                             long_name/
!                                       'Solubility from atmosphere'
!                                       'Surface concentration from ocean'
!                                       'Schmidt number'
!                             units/
!                                   'mol/m^3/atm', 'mol/m^3', 'dimensionless'
!                         flux/
!                              name/
!                                   flux, deltap, kw
!                              long_name/
!                                        'Surface gas flux'
!                                        'ocean-air delta pressure'
!                                        'piston velocity'
!                              units/
!                                    'mol/m^2/s', 'uatm', 'm/s'
!                    air_sea_gas_flux/
!                         implementation/
!                                        ocmip2/
!                                                    num_parameters = 2
!                                        ocmip2_data/
!                                                    num_parameters = 2
!                                        linear/
!                                                    num_parameters = 3
!                         num_flags = 0
!                         use_atm_pressure = t
!                         use_10m_wind_speed = t
!                         pass_through_ice = f
!                         atm/
!                             name/
!                                  pcair, u10, psurf
!                             long_name/
!                                       'Atmospheric concentration'
!                                       'Wind speed at 10 m'
!                                       'Surface atmospheric pressure'
!                             units/
!                                   'mol/mol', 'm/s', 'Pa'
!                         ice/
!                             name/
!                                  alpha, csurf
!                             long_name/
!                                       'Solubility from atmosphere'
!                                       'Surface concentration from ocean'
!                             units/
!                                   'mol/m^3/atm', 'mol/m^3'
!                         flux/
!                              name/
!                                   flux
!                              long_name/
!                                        'Surface gas flux'
!                              units/
!                                    'mol/m^2/s'
!                    air_sea_deposition/
!                         implementation/
!                                        dry/
!                                            num_parameters = 1
!                                        wet/
!                                            num_parameters = 1
!                         num_flags = 0
!                         use_atm_pressure = f
!                         use_10m_wind_speed = f
!                         pass_through_ice = t
!                         atm/
!                             name/
!                                  depostion
!                             long_name/
!                                       'Atmospheric deposition'
!                             units/
!                                   'kg/m^2/s'
!                         ice/
!                             name/
!                             long_name/
!                             units/
!                         flux/
!                              name/
!                                   flux
!                              long_name/
!                                        'Surface deposition'
!                              units/
!                                    'mol/m^2/s'
!                    land_sea_runoff/
!                         implementation/
!                                        river/
!                                              num_parameters = 1
!                         num_flags = 0
!                         use_atm_pressure = f
!                         use_10m_wind_speed = f
!                         pass_through_ice = t
!                         atm/                  ! really land (perhaps should change this?)
!                             name/
!                                  runoff
!                             long_name/
!                                       'Concentration in land runoff'
!                             units/
!                                   'kg/m^3'
!                         ice/
!                             name/
!                             long_name/
!                             units/
!                         flux/
!                              name/
!                                   flux
!                              long_name/
!                                        'Concentration in land runoff'
!                              units/
!                                    'mol/m^3'
!

use fms_mod,           only: write_version_number
use field_manager_mod, only: fm_field_name_len, fm_string_len, fm_dump_list

implicit none
!
!-----------------------------------------------------------------------
! Include variable "version" to be written to log file.
#include<file_version.h>
!-----------------------------------------------------------------------
real, parameter :: bound_tol = 1e-7

private

!
!----------------------------------------------------------------------
!
!       Public routines
!
!----------------------------------------------------------------------
!

public  coupler_types_init
public  coupler_type_copy
public  coupler_type_copy_1d_2d
public  coupler_type_copy_1d_3d

!
!----------------------------------------------------------------------
!
!       Private routines
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Public parameters
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Private parameters
!
!----------------------------------------------------------------------
!

character(len=48), parameter                    :: mod_name = 'coupler_types_mod'

!
!----------------------------------------------------------------------
!
!       Public types
!
!----------------------------------------------------------------------
!

!
!       3-d fields
!

type, public    :: coupler_3d_values_type
  character(len=fm_field_name_len)                      :: name = ' '
  real, pointer, dimension(:,:,:)                       :: values => NULL()
  logical                                               :: mean = .true.
  logical                                               :: override = .false.
  integer                                               :: id_diag = 0
  character(len=fm_string_len)                          :: long_name = ' '
  character(len=fm_string_len)                          :: units = ' '
end type coupler_3d_values_type

type, public    :: coupler_3d_field_type  !{
  character(len=fm_field_name_len)                      :: name = ' '
  integer                                               :: num_fields = 0
  type(coupler_3d_values_type), pointer, dimension(:)   :: field => NULL()
  character(len=fm_string_len)                          :: flux_type = ' '
  character(len=fm_string_len)                          :: implementation = ' '
  real, pointer, dimension(:)                           :: param => NULL()
  logical, pointer, dimension(:)                        :: flag => NULL()
  integer                                               :: atm_tr_index = 0
  character(len=fm_string_len)                          :: ice_restart_file = ' '
  character(len=fm_string_len)                          :: ocean_restart_file = ' '
  logical                                               :: use_atm_pressure
  logical                                               :: use_10m_wind_speed
  logical                                               :: pass_through_ice
  real                                                  :: mol_wt = 0.0
end type coupler_3d_field_type

type, public    :: coupler_3d_bc_type  !{
  integer                                               :: num_bcs = 0
  type(coupler_3d_field_type), pointer, dimension(:)    :: bc => NULL()
end type coupler_3d_bc_type

!
!       2-d fields
!

type, public    :: coupler_2d_values_type
  character(len=fm_field_name_len)                      :: name = ' '
  real, pointer, dimension(:,:)                         :: values => NULL()
  logical                                               :: mean = .true.
  logical                                               :: override = .false.
  integer                                               :: id_diag = 0
  character(len=fm_string_len)                          :: long_name = ' '
  character(len=fm_string_len)                          :: units = ' '
end type coupler_2d_values_type

type, public    :: coupler_2d_field_type  !{
  character(len=fm_field_name_len)                      :: name = ' '
  integer                                               :: num_fields = 0
  type(coupler_2d_values_type), pointer, dimension(:)   :: field => NULL()
  character(len=fm_string_len)                          :: flux_type = ' '
  character(len=fm_string_len)                          :: implementation = ' '
  real, pointer, dimension(:)                           :: param => NULL()
  logical, pointer, dimension(:)                        :: flag => NULL()
  integer                                               :: atm_tr_index = 0
  character(len=fm_string_len)                          :: ice_restart_file = ' '
  character(len=fm_string_len)                          :: ocean_restart_file = ' '
  logical                                               :: use_atm_pressure
  logical                                               :: use_10m_wind_speed
  logical                                               :: pass_through_ice
  real                                                  :: mol_wt = 0.0
end type coupler_2d_field_type

type, public    :: coupler_2d_bc_type  !{
  integer                                               :: num_bcs = 0
  type(coupler_2d_field_type), pointer, dimension(:)    :: bc => NULL()
end type coupler_2d_bc_type

!
!       1-d fields
!

type, public    :: coupler_1d_values_type
  character(len=fm_field_name_len)                      :: name = ' '
  real, pointer, dimension(:)                           :: values => NULL()
  logical                                               :: mean = .true.
  logical                                               :: override = .false.
  integer                                               :: id_diag = 0
  character(len=fm_string_len)                          :: long_name = ' '
  character(len=fm_string_len)                          :: units = ' '
end type coupler_1d_values_type

type, public    :: coupler_1d_field_type  !{
  character(len=fm_field_name_len)                      :: name = ' '
  integer                                               :: num_fields = 0
  type(coupler_1d_values_type), pointer, dimension(:)   :: field => NULL()
  character(len=fm_string_len)                          :: flux_type = ' '
  character(len=fm_string_len)                          :: implementation = ' '
  real, pointer, dimension(:)                           :: param => NULL()
  logical, pointer, dimension(:)                        :: flag => NULL()
  integer                                               :: atm_tr_index = 0
  character(len=fm_string_len)                          :: ice_restart_file = ' '
  character(len=fm_string_len)                          :: ocean_restart_file = ' '
  logical                                               :: use_atm_pressure
  logical                                               :: use_10m_wind_speed
  logical                                               :: pass_through_ice
  real                                                  :: mol_wt = 0.0
end type coupler_1d_field_type

type, public    :: coupler_1d_bc_type  !{
  integer                                               :: num_bcs = 0
  type(coupler_1d_field_type), pointer, dimension(:)    :: bc => NULL()
end type coupler_1d_bc_type

!
!----------------------------------------------------------------------
!
!       Private types
!
!----------------------------------------------------------------------
!

!
!----------------------------------------------------------------------
!
!       Public variables
!
!----------------------------------------------------------------------
!

integer, public :: ind_u10
integer, public :: ind_psurf
integer, public :: ind_pcair
integer, public :: ind_csurf
integer, public :: ind_alpha
integer, public :: ind_sc_no
integer, public :: ind_flux
integer, public :: ind_deltap
integer, public :: ind_kw
integer, public :: ind_deposition
integer, public :: ind_runoff

!
!----------------------------------------------------------------------
!
!       Private variables
!
!----------------------------------------------------------------------
!

logical, save   :: module_is_initialized = .false.

!
!----------------------------------------------------------------------
!
!        Interface definitions for overloaded routines
!
!----------------------------------------------------------------------
!

interface  coupler_type_copy  !{
  module procedure  coupler_type_copy_1d_2d
  module procedure  coupler_type_copy_1d_3d
end interface  coupler_type_copy  !}

!
!-----------------------------------------------------------------------
!
!       Subroutine and function definitions
!
!-----------------------------------------------------------------------
!

contains


!#######################################################################
! <SUBROUTINE NAME="coupler_types_init">
!  <OVERVIEW>
!   Initialize the coupler types
!  </OVERVIEW>
!  <DESCRIPTION>
!   Initialize the coupler types
!  </DESCRIPTION>
!  <TEMPLATE>
!   call coupler_tpyes_init
!  </TEMPLATE>
!

subroutine coupler_types_init

!
!-----------------------------------------------------------------------
!     modules
!-----------------------------------------------------------------------
!

use mpp_mod,           only: stdout, mpp_error, FATAL
use fm_util_mod,       only: fm_util_set_value, fm_util_set_no_overwrite
use fm_util_mod,       only: fm_util_set_caller, fm_util_reset_no_overwrite
use fm_util_mod,       only: fm_util_reset_caller
use field_manager_mod, only: fm_new_list, fm_change_list

implicit none

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'coupler_types_init'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

integer                 :: field_index, outunit
character(len=128)      :: error_msg

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!       Return if already intialized
!

if (module_is_initialized) then  !{
  return
endif  !}

!
!       Write out the version of the file to the log file
!
call write_version_number(trim(mod_name), version)
!
!       Set other defaults for the fm_util_set_value routines
!

call fm_util_set_no_overwrite(.true.)
call fm_util_set_caller(sub_name)

!
!       Be sure that the various lists and fields are defined in the field manager tree
!

if (fm_new_list('/coupler_mod') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "coupler_mod" list')
endif  !}

if (fm_new_list('/coupler_mod/GOOD') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "GOOD" list')
endif  !}
call fm_util_set_value('/coupler_mod/GOOD/good_coupler_mod_list', 'GOOD', append = .true.)

if (fm_new_list('/coupler_mod/fluxes') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "/coupler_mod/fluxes" list')
endif  !}
call fm_util_set_value('/coupler_mod/GOOD/good_coupler_mod_list', 'fluxes', append = .true.)

if (fm_new_list('/coupler_mod/types') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "/coupler_mod/types" list')
endif  !}
call fm_util_set_value('/coupler_mod/GOOD/good_coupler_mod_list', 'types', append = .true.)

!
!       change to the "/coupler_mod/types" list
!

if (.not. fm_change_list('/coupler_mod/types')) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to "/coupler_mod/types"')
endif  !}

!
!       Define the air_sea_gas_flux_generic type
!

!       add the new type

if (fm_new_list('air_sea_gas_flux_generic') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic" list')
endif  !}

!       add the implementation list

if (fm_new_list('air_sea_gas_flux_generic/implementation') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/implementation" list')
endif  !}

!       add the names of the different implementations

if (fm_new_list('air_sea_gas_flux_generic/implementation/ocmip2') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/implementation/ocmip2" list')
endif  !}
call fm_util_set_value('air_sea_gas_flux_generic/implementation/ocmip2/num_parameters', 2)

!       add some scalar quantaties

call fm_util_set_value('air_sea_gas_flux_generic/num_flags', 0)
call fm_util_set_value('air_sea_gas_flux_generic/use_atm_pressure', .true.)
call fm_util_set_value('air_sea_gas_flux_generic/use_10m_wind_speed', .true.)
call fm_util_set_value('air_sea_gas_flux_generic/pass_through_ice', .false.)

!       add required fields that will come from the atmosphere model

if (fm_new_list('air_sea_gas_flux_generic/atm') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/atm" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_pcair = field_index
call fm_util_set_value('air_sea_gas_flux_generic/atm/name',      'pcair',                     index = ind_pcair)
call fm_util_set_value('air_sea_gas_flux_generic/atm/long_name', 'Atmospheric concentration', index = ind_pcair)
call fm_util_set_value('air_sea_gas_flux_generic/atm/units',     'mol/mol',                   index = ind_pcair)

field_index = field_index + 1
ind_u10 = field_index
call fm_util_set_value('air_sea_gas_flux_generic/atm/name',      'u10',                index = ind_u10)
call fm_util_set_value('air_sea_gas_flux_generic/atm/long_name', 'Wind speed at 10 m', index = ind_u10)
call fm_util_set_value('air_sea_gas_flux_generic/atm/units',     'm/s',                index = ind_u10)

field_index = field_index + 1
ind_psurf = field_index
call fm_util_set_value('air_sea_gas_flux_generic/atm/name',      'psurf',                        index = ind_psurf)
call fm_util_set_value('air_sea_gas_flux_generic/atm/long_name', 'Surface atmospheric pressure', index = ind_psurf)
call fm_util_set_value('air_sea_gas_flux_generic/atm/units',     'Pa',                           index = ind_psurf)

!       add required fields that will come from the ice model

if (fm_new_list('air_sea_gas_flux_generic/ice') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/ice" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_alpha = field_index
call fm_util_set_value('air_sea_gas_flux_generic/ice/name',      'alpha',                                                index = ind_alpha)
call fm_util_set_value('air_sea_gas_flux_generic/ice/long_name', 'Solubility w.r.t. atmosphere', index = ind_alpha)
call fm_util_set_value('air_sea_gas_flux_generic/ice/units',     'mol/m^3/atm',                                          index = ind_alpha)

field_index = field_index + 1
ind_csurf = field_index
call fm_util_set_value('air_sea_gas_flux_generic/ice/name',      'csurf',                                         index = ind_csurf)
call fm_util_set_value('air_sea_gas_flux_generic/ice/long_name', 'Ocean concentration', index = ind_csurf)
call fm_util_set_value('air_sea_gas_flux_generic/ice/units',     'mol/m^3',                                       index = ind_csurf)

field_index = field_index + 1
ind_sc_no = field_index
call fm_util_set_value('air_sea_gas_flux_generic/ice/name',      'sc_no',                                         index = ind_sc_no)
call fm_util_set_value('air_sea_gas_flux_generic/ice/long_name', 'Schmidt number', index = ind_sc_no)
call fm_util_set_value('air_sea_gas_flux_generic/ice/units',     'dimensionless',                                       index = ind_sc_no)

!       add the flux output field(s)

if (fm_new_list('air_sea_gas_flux_generic/flux') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux_generic/flux" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_flux = field_index
call fm_util_set_value('air_sea_gas_flux_generic/flux/name',      'flux',         index = ind_flux)
call fm_util_set_value('air_sea_gas_flux_generic/flux/long_name', 'Surface flux', index = ind_flux)
call fm_util_set_value('air_sea_gas_flux_generic/flux/units',     'mol/m^2/s',    index = ind_flux)

field_index = field_index + 1
ind_deltap = field_index
call fm_util_set_value('air_sea_gas_flux_generic/flux/name',      'deltap',         index = ind_deltap)
call fm_util_set_value('air_sea_gas_flux_generic/flux/long_name', 'Ocean-air delta pressure', index = ind_deltap)
call fm_util_set_value('air_sea_gas_flux_generic/flux/units',     'uatm',    index = ind_deltap)

field_index = field_index + 1
ind_kw = field_index
call fm_util_set_value('air_sea_gas_flux_generic/flux/name',      'kw',         index = ind_kw)
call fm_util_set_value('air_sea_gas_flux_generic/flux/long_name', 'Piston velocity', index = ind_kw)
call fm_util_set_value('air_sea_gas_flux_generic/flux/units',     'm/s',    index = ind_kw)

!
!       Define the air_sea_gas_flux type
!

!       add the new type

if (fm_new_list('air_sea_gas_flux') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux" list')
endif  !}

!       add the implementation list

if (fm_new_list('air_sea_gas_flux/implementation') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation" list')
endif  !}

!       add the names of the different implementations

if (fm_new_list('air_sea_gas_flux/implementation/ocmip2') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation/ocmip2" list')
endif  !}
call fm_util_set_value('air_sea_gas_flux/implementation/ocmip2/num_parameters', 2)
if (fm_new_list('air_sea_gas_flux/implementation/ocmip2_data') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation/ocmip2_data" list')
endif  !}
call fm_util_set_value('air_sea_gas_flux/implementation/ocmip2_data/num_parameters', 2)
if (fm_new_list('air_sea_gas_flux/implementation/linear') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/implementation/linear" list')
endif  !}
call fm_util_set_value('air_sea_gas_flux/implementation/linear/num_parameters', 3)

!       add some scalar quantaties

call fm_util_set_value('air_sea_gas_flux/num_flags', 0)
call fm_util_set_value('air_sea_gas_flux/use_atm_pressure', .true.)
call fm_util_set_value('air_sea_gas_flux/use_10m_wind_speed', .true.)
call fm_util_set_value('air_sea_gas_flux/pass_through_ice', .false.)

!       add required fields that will come from the atmosphere model

if (fm_new_list('air_sea_gas_flux/atm') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/atm" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_pcair = field_index
call fm_util_set_value('air_sea_gas_flux/atm/name',      'pcair',                     index = ind_pcair)
call fm_util_set_value('air_sea_gas_flux/atm/long_name', 'Atmospheric concentration', index = ind_pcair)
call fm_util_set_value('air_sea_gas_flux/atm/units',     'mol/mol',                   index = ind_pcair)

field_index = field_index + 1
ind_u10 = field_index
call fm_util_set_value('air_sea_gas_flux/atm/name',      'u10',                index = ind_u10)
call fm_util_set_value('air_sea_gas_flux/atm/long_name', 'Wind speed at 10 m', index = ind_u10)
call fm_util_set_value('air_sea_gas_flux/atm/units',     'm/s',                index = ind_u10)

field_index = field_index + 1
ind_psurf = field_index
call fm_util_set_value('air_sea_gas_flux/atm/name',      'psurf',                        index = ind_psurf)
call fm_util_set_value('air_sea_gas_flux/atm/long_name', 'Surface atmospheric pressure', index = ind_psurf)
call fm_util_set_value('air_sea_gas_flux/atm/units',     'Pa',                           index = ind_psurf)

!       add required fields that will come from the ice model

if (fm_new_list('air_sea_gas_flux/ice') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/ice" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_alpha = field_index
call fm_util_set_value('air_sea_gas_flux/ice/name',      'alpha',                                                index = ind_alpha)
call fm_util_set_value('air_sea_gas_flux/ice/long_name', 'Solubility from atmosphere times Schmidt number term', index = ind_alpha)
call fm_util_set_value('air_sea_gas_flux/ice/units',     'mol/m^3/atm',                                          index = ind_alpha)

field_index = field_index + 1
ind_csurf = field_index
call fm_util_set_value('air_sea_gas_flux/ice/name',      'csurf',                                         index = ind_csurf)
call fm_util_set_value('air_sea_gas_flux/ice/long_name', 'Ocean concentration times Schmidt number term', index = ind_csurf)
call fm_util_set_value('air_sea_gas_flux/ice/units',     'mol/m^3',                                       index = ind_csurf)

!       add the flux output field(s)

if (fm_new_list('air_sea_gas_flux/flux') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_gas_flux/flux" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_flux = field_index
call fm_util_set_value('air_sea_gas_flux/flux/name',      'flux',         index = ind_flux)
call fm_util_set_value('air_sea_gas_flux/flux/long_name', 'Surface flux', index = ind_flux)
call fm_util_set_value('air_sea_gas_flux/flux/units',     'mol/m^2/s',    index = ind_flux)

!
!       Define the air_sea_deposition type
!

!       add the new type

if (fm_new_list('air_sea_deposition') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition" list')
endif  !}

!       add the implementation list

if (fm_new_list('air_sea_deposition/implementation') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/implementation" list')
endif  !}

!       add the names of the different implementations

if (fm_new_list('air_sea_deposition/implementation/dry') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/implementation/dry" list')
endif  !}
call fm_util_set_value('air_sea_deposition/implementation/dry/num_parameters', 1)
if (fm_new_list('air_sea_deposition/implementation/wet') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/implementation/wet" list')
endif  !}
call fm_util_set_value('air_sea_deposition/implementation/wet/num_parameters', 1)

!       add some scalar quantaties

call fm_util_set_value('air_sea_deposition/num_flags', 0)
call fm_util_set_value('air_sea_deposition/use_atm_pressure', .false.)
call fm_util_set_value('air_sea_deposition/use_10m_wind_speed', .false.)
call fm_util_set_value('air_sea_deposition/pass_through_ice', .true.)

!       add required fields that will come from the atmosphere model

if (fm_new_list('air_sea_deposition/atm') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/atm" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_deposition = field_index
call fm_util_set_value('air_sea_deposition/atm/name',      'deposition',             index = ind_deposition)
call fm_util_set_value('air_sea_deposition/atm/long_name', 'Atmospheric deposition', index = ind_deposition)
call fm_util_set_value('air_sea_deposition/atm/units',     'kg/m^2/s',               index = ind_deposition)

!       add required fields that will come from the ice model

if (fm_new_list('air_sea_deposition/ice') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/ice" list')
endif  !}

field_index = 0

call fm_util_set_value('air_sea_deposition/ice/name',      ' ', index = 0)
call fm_util_set_value('air_sea_deposition/ice/long_name', ' ', index = 0)
call fm_util_set_value('air_sea_deposition/ice/units',     ' ', index = 0)

!       add the flux output field(s)

if (fm_new_list('air_sea_deposition/flux') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "air_sea_deposition/flux" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_flux = field_index
call fm_util_set_value('air_sea_deposition/flux/name',      'flux',               index = ind_flux)
call fm_util_set_value('air_sea_deposition/flux/long_name', 'Surface deposition', index = ind_flux)
call fm_util_set_value('air_sea_deposition/flux/units',     'mol/m^2/s',          index = ind_flux)

!
!       Define the land_sea_runoff type
!

!       add the new type

if (fm_new_list('land_sea_runoff') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff" list')
endif  !}

!       add the implementation list

if (fm_new_list('land_sea_runoff/implementation') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/implementation" list')
endif  !}

!       add the names of the different implementations

if (fm_new_list('land_sea_runoff/implementation/river') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/implementation/river" list')
endif  !}
call fm_util_set_value('land_sea_runoff/implementation/river/num_parameters', 1)

!       add some scalar quantaties

call fm_util_set_value('land_sea_runoff/num_flags', 0)
call fm_util_set_value('land_sea_runoff/use_atm_pressure', .false.)
call fm_util_set_value('land_sea_runoff/use_10m_wind_speed', .false.)
call fm_util_set_value('land_sea_runoff/pass_through_ice', .true.)

!       add required fields that will come from the land model (the array name is still called "atm")

if (fm_new_list('land_sea_runoff/atm') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/atm" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_runoff = field_index
call fm_util_set_value('land_sea_runoff/atm/name',      'runoff',                       index = ind_runoff)
call fm_util_set_value('land_sea_runoff/atm/long_name', 'Concentration in land runoff', index = ind_runoff)
call fm_util_set_value('land_sea_runoff/atm/units',     'mol/m^3',                      index = ind_runoff)

!       add required fields that will come from the ice model

if (fm_new_list('land_sea_runoff/ice') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/ice" list')
endif  !}

field_index = 0

call fm_util_set_value('land_sea_runoff/ice/name',      ' ', index = 0)
call fm_util_set_value('land_sea_runoff/ice/long_name', ' ', index = 0)
call fm_util_set_value('land_sea_runoff/ice/units',     ' ', index = 0)

!       add the flux output field(s)

if (fm_new_list('land_sea_runoff/flux') .le. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not set the "land_sea_runoff/flux" list')
endif  !}

field_index = 0

field_index = field_index + 1
ind_flux = field_index
call fm_util_set_value('land_sea_runoff/flux/name',      'flux',                         index = ind_flux)
call fm_util_set_value('land_sea_runoff/flux/long_name', 'Concentration in land runoff', index = ind_flux)
call fm_util_set_value('land_sea_runoff/flux/units',     'mol/m^3',                      index = ind_flux)

!
!       change back to root list
!

if (.not. fm_change_list('/')) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Could not change to "/"')
endif  !}

!
!       Reset the defaults for the fm_util_set_value calls
!

call fm_util_reset_no_overwrite
call fm_util_reset_caller

module_is_initialized = .true.

!
!       Dump the coupler_mod types list
!
outunit = stdout()
write (outunit,*)
write (outunit,*) 'Dumping coupler_mod/types tree'
if (.not. fm_dump_list('/coupler_mod/types', recursive = .true.)) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Problem dumping /coupler_mod/types tree')
endif  !}

return

end subroutine  coupler_types_init  !}
! </SUBROUTINE> NAME="coupler_types_init"



!#######################################################################
! <SUBROUTINE NAME="coupler_type_copy_1d_2d">
!  <OVERVIEW>
!   Copy fields from one coupler type to another.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Copy fields from one coupler type to another.
!   Specific version for generic coupler_type_copy.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call coupler_type_copy(var_in, var_out, is, ie, js, je,     &
!        diag_name, axes, time, suffix = 'something')
!		
!  </TEMPLATE>
!  <IN NAME="var_in" TYPE="coupler_1d_bc_type">
!   variable to copy information from
!  </IN>
!  <IN NAME="var_out" TYPE="coupler_2d_bc_type">
!   variable to copy information to
!  </IN>
!  <IN NAME="is" TYPE="integer">
!   lower bound of first dimension
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   upper bound of first dimension
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   lower bound of second dimension
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   upper bound of second dimension
!  </IN>
!  <IN NAME="diag_name" TYPE="character">
!   name for diagnostic file--if blank, then don't register the fields
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   array of axes identifiers for diagnostic variable registration
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   model time variable for registering diagnostic field
!  </IN>
!  <IN NAME="suffix" TYPE="character">
!   optional suffix to make the name identifier unique
!  </IN>
!

subroutine coupler_type_copy_1d_2d(var_in, var_out, is, ie, js, je,     &
     diag_name, axes, time, suffix)  !{

!
!-----------------------------------------------------------------------
!     modules
!-----------------------------------------------------------------------
!

use time_manager_mod, only: time_type
use diag_manager_mod, only: register_diag_field
use mpp_mod,          only: mpp_error, FATAL

implicit none

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

type(coupler_1d_bc_type), intent(in)    :: var_in
type(coupler_2d_bc_type), intent(inout) :: var_out
integer, intent(in)                     :: is
integer, intent(in)                     :: ie
integer, intent(in)                     :: js
integer, intent(in)                     :: je
character(len=*), intent(in)            :: diag_name
integer, dimension(:), intent(in)       :: axes
type(time_type), intent(in)             :: time
character(len=*), intent(in), optional  :: suffix
 
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'coupler_type_copy_1d_2d'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

character(len=128)      :: error_msg
integer                 :: m
integer                 :: n

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!       Error if output fields is not zero
!

if (var_out%num_bcs .ne. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Number of output fields is non-zero')
endif  !}
  
var_out%num_bcs = var_in%num_bcs

!
!       Return if no input fields
!

if (var_in%num_bcs .ne. 0) then  !{
  if (associated(var_out%bc)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' var_out%bc already associated')
  endif  !}
  allocate ( var_out%bc(var_out%num_bcs) )
  do n = 1, var_out%num_bcs  !{
    var_out%bc(n)%name = var_in%bc(n)%name
    var_out%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
    var_out%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
    var_out%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
    var_out%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
    var_out%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
    var_out%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
    var_out%bc(n)%mol_wt = var_in%bc(n)%mol_wt
    var_out%bc(n)%num_fields = var_in%bc(n)%num_fields
    if (associated(var_out%bc(n)%field)) then  !{
      write (error_msg, *) trim(error_header), ' var_out%bc(', n, ')%field already associated'
      call mpp_error(FATAL, trim(error_msg))
    endif  !}
    allocate ( var_out%bc(n)%field(var_out%bc(n)%num_fields) )
    do m = 1, var_out%bc(n)%num_fields  !{
      if (present(suffix)) then  !{
        var_out%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // trim(suffix)
      else  !}{
        var_out%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
      endif  !}
      var_out%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
      var_out%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
      var_out%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
      if (associated(var_out%bc(n)%field(m)%values)) then  !{
        write (error_msg, *) trim(error_header), ' var_out%bc(', n, ')%field(', m, ')%values already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif  !}
      allocate ( var_out%bc(n)%field(m)%values(is:ie,js:je) )
      var_out%bc(n)%field(m)%values = 0.0
      if (diag_name .ne. ' ') then  !{
        if (size(axes) .lt. 2) then  !{
          call mpp_error(FATAL, trim(error_header) // ' axes less than 2 elements')
        endif  !}
        var_out%bc(n)%field(m)%id_diag = register_diag_field(diag_name,                &
             var_out%bc(n)%field(m)%name, axes(1:2), Time,                             &
             var_out%bc(n)%field(m)%long_name, var_out%bc(n)%field(m)%units )
      endif  !}
    enddo  !} m
  enddo  !} n

endif  !}

return

end subroutine  coupler_type_copy_1d_2d  !}
! </SUBROUTINE> NAME="coupler_type_copy_1d_2d


!#######################################################################
! <SUBROUTINE NAME="coupler_type_copy_1d_3d">
!  <OVERVIEW>
!   Copy fields from one coupler type to another.
!  </OVERVIEW>
!  <DESCRIPTION>
!   Copy fields from one coupler type to another.
!   Specific version for generic coupler_type_copy.
!  </DESCRIPTION>
!  <TEMPLATE>
!   call coupler_type_copy(var_in, var_out, is, ie, js, je, kd, &
!        diag_name, axes, time, suffix = 'something')
!		
!  </TEMPLATE>
!  <IN NAME="var_in" TYPE="coupler_1d_bc_type">
!   variable to copy information from
!  </IN>
!  <IN NAME="var_out" TYPE="coupler_3d_bc_type">
!   variable to copy information to
!  </IN>
!  <IN NAME="is" TYPE="integer">
!   lower bound of first dimension
!  </IN>
!  <IN NAME="ie" TYPE="integer">
!   upper bound of first dimension
!  </IN>
!  <IN NAME="js" TYPE="integer">
!   lower bound of second dimension
!  </IN>
!  <IN NAME="je" TYPE="integer">
!   upper bound of second dimension
!  </IN>
!  <IN NAME="kd" TYPE="integer">
!   third dimension
!  </IN>
!  <IN NAME="diag_name" TYPE="character">
!   name for diagnostic file--if blank, then don't register the fields
!  </IN>
!  <IN NAME="axes" TYPE="integer">
!   array of axes identifiers for diagnostic variable registration
!  </IN>
!  <IN NAME="time" TYPE="time_type">
!   model time variable for registering diagnostic field
!  </IN>
!  <IN NAME="suffix" TYPE="character">
!   optional suffix to make the name identifier unique
!  </IN>
!

subroutine coupler_type_copy_1d_3d(var_in, var_out, is, ie, js, je, kd, &
     diag_name, axes, time, suffix)  !{

!
!-----------------------------------------------------------------------
!     modules
!-----------------------------------------------------------------------
!

use time_manager_mod, only: time_type
use diag_manager_mod, only: register_diag_field
use mpp_mod,          only: mpp_error, FATAL

implicit none

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!

type(coupler_1d_bc_type), intent(in)    :: var_in
type(coupler_3d_bc_type), intent(inout) :: var_out
integer, intent(in)                     :: is
integer, intent(in)                     :: ie
integer, intent(in)                     :: js
integer, intent(in)                     :: je
integer, intent(in)                     :: kd
character(len=*), intent(in)            :: diag_name
integer, dimension(:), intent(in)       :: axes
type(time_type), intent(in)             :: time
character(len=*), intent(in), optional  :: suffix
 
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

character(len=64), parameter    :: sub_name = 'coupler_type_copy_1d_3d'
character(len=256), parameter   :: error_header =                               &
     '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

!
!-----------------------------------------------------------------------
!     local variables
!-----------------------------------------------------------------------
!

character(len=128)      :: error_msg
integer                 :: m
integer                 :: n

!
! =====================================================================
!     begin executable code
! =====================================================================
!

!
!       Error if output fields is not zero
!

if (var_out%num_bcs .ne. 0) then  !{
  call mpp_error(FATAL, trim(error_header) // ' Number of output fields is non-zero')
endif  !}
  
var_out%num_bcs = var_in%num_bcs

!
!       Return if no input fields
!

if (var_in%num_bcs .ne. 0) then  !{
  if (associated(var_out%bc)) then  !{
    call mpp_error(FATAL, trim(error_header) // ' var_out%bc already associated')
  endif  !}
  allocate ( var_out%bc(var_out%num_bcs) )
  do n = 1, var_out%num_bcs  !{
    var_out%bc(n)%name = var_in%bc(n)%name
    var_out%bc(n)%atm_tr_index = var_in%bc(n)%atm_tr_index
    var_out%bc(n)%ice_restart_file = var_in%bc(n)%ice_restart_file
    var_out%bc(n)%ocean_restart_file = var_in%bc(n)%ocean_restart_file
    var_out%bc(n)%use_atm_pressure = var_in%bc(n)%use_atm_pressure
    var_out%bc(n)%use_10m_wind_speed = var_in%bc(n)%use_10m_wind_speed
    var_out%bc(n)%pass_through_ice = var_in%bc(n)%pass_through_ice
    var_out%bc(n)%mol_wt = var_in%bc(n)%mol_wt
    var_out%bc(n)%num_fields = var_in%bc(n)%num_fields
    if (associated(var_out%bc(n)%field)) then  !{
      write (error_msg, *) trim(error_header), ' var_out%bc(', n, ')%field already associated'
      call mpp_error(FATAL, trim(error_msg))
    endif  !}
    allocate ( var_out%bc(n)%field(var_out%bc(n)%num_fields) )
    do m = 1, var_out%bc(n)%num_fields  !{
      if (present(suffix)) then  !{
        var_out%bc(n)%field(m)%name = trim(var_in%bc(n)%field(m)%name) // suffix
      else  !}{
        var_out%bc(n)%field(m)%name = var_in%bc(n)%field(m)%name
      endif  !}
      var_out%bc(n)%field(m)%long_name = var_in%bc(n)%field(m)%long_name
      var_out%bc(n)%field(m)%units = var_in%bc(n)%field(m)%units
      var_out%bc(n)%field(m)%mean = var_in%bc(n)%field(m)%mean
      if (associated(var_out%bc(n)%field(m)%values)) then  !{
        write (error_msg, *) trim(error_header), ' var_out%bc(', n, ')%field(', m, ')%values already associated'
        call mpp_error(FATAL, trim(error_msg))
      endif  !}
      allocate ( var_out%bc(n)%field(m)%values(is:ie,js:je,kd) )
      var_out%bc(n)%field(m)%values = 0.0
      if (diag_name .ne. ' ') then  !{
        if (size(axes) .lt. 3) then  !{
          call mpp_error(FATAL, trim(error_header) // ' axes less than 3 elements')
        endif  !}
        var_out%bc(n)%field(m)%id_diag = register_diag_field(diag_name,                &
             var_out%bc(n)%field(m)%name, axes(1:3), Time,                             &
             var_out%bc(n)%field(m)%long_name, var_out%bc(n)%field(m)%units )
      endif  !}
    enddo  !} m
  enddo  !} n

endif  !}

return

end subroutine  coupler_type_copy_1d_3d  !}
! </SUBROUTINE> NAME="coupler_type_copy_1d_3d

end module coupler_types_mod  !}
