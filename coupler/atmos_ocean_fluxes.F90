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

!> \author Richard Slater <Richard.Slater@noaa.gov>
!! \author John Dunne <John.Dunne@noaa.gov>
!
!! \brief Ocean Carbon Model Intercomparison Study II: Gas exchange coupler.
!! Implementation of routines to solve the gas fluxes at the
!! ocean surface for a coupled model as outlined in the Biotic-HOWTO
!! documentation, revision 1.7, 1999/10/05.
!
!! \link http://www.ipsl.jussieu.fr/OCMIP/phase2/simulations/Biotic/HOWTO-Biotic.html \endlink
!
!! This module will take fields from an atmospheric and an
!! oceanic model and calculate ocean surface fluxes for
!! CO2, O2, CFC-11 or CFC-12 as outlined in the various
!! HOWTO documents at the OCMIP2 website. Multiple instances
!! of a given tracer may be given, resulting in multiple
!! surface fluxes. Additionally, data may be overridden at
!! the individual fields, or fluxes. This could be used in
!! the absence of an atmospheric or oceanic model.
module  atmos_ocean_fluxes_mod

!! Modules Included:
!
!! <table>
!!   <tr>
!!     <th>Module Name</th>
!!     <th>Functions Included</th>
!!   </tr>
!!   <tr>
!!     <td>mpp_mod</td>
!!     <td>stdout, stdlog, mpp_error, FATAL, mpp_sum, mpp_npes</td>
!!   </tr>
!!   <tr>
!!     <td>fms_mod</td>
!!     <td>write_version_number</td>
!!   </tr>
!!   <tr>
!!     <td>coupler_types_mod</td>
!!     <td>coupler_1d_bc_type, ind_alpha, ind_csurf, ind_sc_no, ind_pcair, ind_u10, ind_psurf,
!!         ind_deposition, ind_runoff, ind_flux, ind_deltap, ind_kw</td>
!!   </tr>
!!   <tr>
!!     <td>field_manager_mod</td>
!!     <td>fm_path_name_len, fm_string_len, fm_exists, fm_get_index, fm_new_list, fm_get_current_list,
!!         fm_change_list, fm_field_name_len, fm_type_name_len, fm_dump_list, fm_loop_over_list</td>
!!   </tr>
!!   <tr>
!!     <td>fm_util_mod</td>
!!     <td>fm_util_default_caller, fm_util_get_length, fm_util_set_value, fm_util_set_good_name_list,
!!         fm_util_set_no_overwrite, fm_util_set_caller, fm_util_reset_good_name_list,
!!         fm_util_reset_no_overwrite, fm_util_reset_caller, fm_util_get_string_array,
!!         fm_util_check_for_bad_fields, fm_util_get_string, fm_util_get_real_array, fm_util_get_real,
!!         fm_util_get_integer, fm_util_get_logical, fm_util_get_logical_array</td>
!!   </tr>
!! </table>
use mpp_mod,           only: stdout, stdlog, mpp_error, FATAL, mpp_sum, mpp_npes
use mpp_mod,           only: mpp_pe,mpp_root_pe
use fms_mod,           only: write_version_number

use coupler_types_mod, only: coupler_1d_bc_type
use coupler_types_mod, only: ind_alpha, ind_csurf, ind_sc_no
use coupler_types_mod, only: ind_pcair, ind_u10, ind_psurf
use coupler_types_mod, only: ind_deposition
use coupler_types_mod, only: ind_runoff
use coupler_types_mod, only: ind_flux, ind_deltap, ind_kw!, ind_flux0
use constants_mod,     only: WTMAIR,rdgas,vonkarm

use field_manager_mod, only: fm_path_name_len, fm_string_len, fm_exists, fm_get_index
use field_manager_mod, only: fm_new_list, fm_get_current_list, fm_change_list
use field_manager_mod, only: fm_field_name_len, fm_type_name_len, fm_dump_list
use field_manager_mod, only: fm_loop_over_list

use fm_util_mod,       only: fm_util_default_caller
use fm_util_mod,       only: fm_util_get_length
use fm_util_mod,       only: fm_util_set_value, fm_util_set_good_name_list, fm_util_set_no_overwrite
use fm_util_mod,       only: fm_util_set_caller, fm_util_reset_good_name_list, fm_util_reset_no_overwrite
use fm_util_mod,       only: fm_util_reset_caller, fm_util_get_string_array, fm_util_check_for_bad_fields
use fm_util_mod,       only: fm_util_get_string, fm_util_get_real_array, fm_util_get_real, fm_util_get_integer
use fm_util_mod,       only: fm_util_get_logical, fm_util_get_logical_array

implicit none ; private

public  :: atmos_ocean_dep_fluxes_calc
public  :: atmos_ocean_fluxes_calc
public  :: atmos_ocean_fluxes_init
public  :: aof_set_coupler_flux

character(len=48), parameter    :: mod_name = 'atmos_ocean_fluxes_mod'

! Include variable "version" to be written to log file.
#include<file_version.h>

!-----------------------------------------------------------------------

!       Subroutine and function definitions

!-----------------------------------------------------------------------

contains

!> \brief Set the values for a coupler flux and return its index (0 on error)
!
!! \throw FATAL, "Empty name given"
!!     Name is empty
!! \throw FATAL, "Could not get coupler flux"
!!     coupler_index is less than 1
!! \throw FATAL, "Could not set coupler flux"
!!     coupler_index is less than 1
!! \throw FATAL, "Could not get the current list"
!!     Current list is empty
!! \throw FATAL, "Could not change to the new list"
!!     fm_change_list(coupler_list) returns false
!! \throw FATAL, "Blank flux_type given"
!!     flux_type or implementation is empty
!! \throw FATAL, "Undefined flux_type given from field_table"
!! \throw FATAL, "Undefined flux_type given as argument to the subroutine"
!! \throw FATAL, "Undefined flux_type/implementation (implementation given from field_table)"
!!     flux_type does not equal flux_type_test
!! \throw FATAL, "Undefined flux_type/implementation (flux_type given from field_table)"
!! \throw FATAL, "Undefined flux_type/implementation (both given from field_table)"
!! \throw FATAL, "Undefined flux_type/implementation given as argument to the subroutine"
!! \throw NOTE, "Number of parameters provided for [variable] does not match the number of parameters required"
!!     Mismatch between parameter input and the parameters being replaced
!! \throw FATAL, "Could not change back to [current_list]"
!! \throw FATAL, "Empty [name] list"
function aof_set_coupler_flux(name, flux_type, implementation, atm_tr_index, param, flag, &
     mol_wt, ice_restart_file, ocean_restart_file, units, caller, verbosity) &
         result (coupler_index)

  character(len=*), intent(in)                :: name !< name
  character(len=*), intent(in)                :: flux_type !< flux_type
  character(len=*), intent(in)                :: implementation !< implementation
  integer, intent(in), optional               :: atm_tr_index !< atm_tr_index
  real, intent(in), dimension(:), optional    :: param !< param
  logical, intent(in), dimension(:), optional :: flag !< flag
  real, intent(in), optional                  :: mol_wt !< mol_wt
  character(len=*), intent(in), optional      :: ice_restart_file !< ice_restart_file
  character(len=*), intent(in), optional      :: ocean_restart_file !< ocean_restart_file
  character(len=*), intent(in), optional      :: units !< units
  character(len=*), intent(in), optional      :: caller !< caller
  integer,          intent(in), optional      :: verbosity  !< A 0-9 integer indicating a level of verbosity.

  integer :: coupler_index

  character(len=48), parameter  :: sub_name = 'aof_set_coupler_flux'

  integer                                                 :: n
  integer                                                 :: length
  integer                                                 :: num_parameters
  integer                                                 :: outunit
  character(len=fm_path_name_len)                         :: coupler_list
  character(len=fm_path_name_len)                         :: current_list
  character(len=fm_string_len)                            :: flux_type_test
  character(len=fm_string_len)                            :: implementation_test
  character(len=256)                                      :: error_header
  character(len=256)                                      :: warn_header
  character(len=256)                                      :: note_header
  character(len=128)                                      :: flux_list
  character(len=128)                                      :: caller_str
  character(len=fm_string_len), pointer, dimension(:)     :: good_list => NULL()
  character(len=256)                                      :: long_err_msg
  integer :: verbose ! An integer indicating the level of verbosity.

  verbose = 5 ; if (present(verbosity)) verbose = verbosity

  !>       Set the caller string and headers.
  if (present(caller)) then
    caller_str = '[' // trim(caller) // ']'
  else
    caller_str = fm_util_default_caller
  endif

  error_header = '==>Error from ' // trim(mod_name) // &
                 '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
  warn_header = '==>Warning from ' // trim(mod_name) // &
                '(' // trim(sub_name) // ')' // trim(caller_str) // ':'
  note_header = '==>Note from ' // trim(mod_name) // &
                '(' // trim(sub_name) // ')' // trim(caller_str) // ':'

  !>       Check that a name is given (fatal if not).
  if (name .eq. ' ') then
    call mpp_error(FATAL, trim(error_header) // ' Empty name given')
  endif
  outunit = stdout()
  if (verbose >= 5) then
    write (outunit,*)
    write (outunit,*) trim(note_header), ' Processing coupler fluxes ', trim(name)
  endif

  !>       Define the coupler list name.
  coupler_list = '/coupler_mod/fluxes/' // trim(name)

  !>       Check whether a flux has already been set for this name, and if so, return
  !!       the index for it (this is because the fluxes may be defined in both the atmosphere
  !!       and ocean models) (check whether the good_list list exists, since this will
  !!       indicate that this routine has already been called, and not just that
  !!       the field table input has this list defined)

  if (fm_exists('/coupler_mod/GOOD/fluxes/' // trim(name) // '/good_list')) then
    if (verbose >= 5) then
      write (outunit,*)
      write (outunit,*) trim(note_header), ' Using previously defined coupler flux'
    endif
    coupler_index = fm_get_index(coupler_list)
    if (coupler_index .le. 0) then
      call mpp_error(FATAL, trim(error_header) // ' Could not get coupler flux ')
    endif

  !>       Allow atm_tr_index to be set here, since it will only be set from atmospheric
  !!       PEs, and the atmospheric routines call this routine last, thus overwriting the
  !!       current value is safe (furthermore, this is not a value which could have any meaningful
  !!       value set from the run script.

    if (present(atm_tr_index)) then
      if (verbose >= 5) &
        write (outunit,*) trim(note_header), ' Redefining atm_tr_index to ', atm_tr_index
      call fm_util_set_value(trim(coupler_list) // '/atm_tr_index', atm_tr_index, no_create = .true., &
           no_overwrite = .false., caller = caller_str)
    endif
    return
  endif

  !>       Set a new coupler flux and get its index.
  coupler_index = fm_new_list(coupler_list)
  if (coupler_index .le. 0) then
    call mpp_error(FATAL, trim(error_header) // ' Could not set coupler flux ')
  endif

  !>       Change to the new list, first saving the current list.
  current_list = fm_get_current_list()
  if (current_list .eq. ' ') then
    call mpp_error(FATAL, trim(error_header) // ' Could not get the current list')
  endif

  if (.not. fm_change_list(coupler_list)) then
    call mpp_error(FATAL, trim(error_header) // ' Could not change to the new list')
  endif

  !>       Set the array in which to save the valid names for this list,
  !!       used later for a consistency check. This is used in the fm_util_set_value
  !!       routines to make the list of valid values.
  call fm_util_set_good_name_list('/coupler_mod/GOOD/fluxes/' // trim(name) // '/good_list')

  !>       Set other defaults for the fm_util_set_value routines.
  call fm_util_set_no_overwrite(.true.)
  call fm_util_set_caller(caller_str)

  !>       Set various values to given values, or to defaults if not given.
  if (flux_type .eq. ' ') then
    call mpp_error(FATAL, trim(error_header) // ' Blank flux_type given')
  else
    if (fm_exists('/coupler_mod/types/' // trim(flux_type))) then
      call fm_util_set_value('flux_type', flux_type)

  !>       Check that the flux_type that we will use (possibly given from the field_table)
  !!       is defined.
      flux_type_test = fm_util_get_string('flux_type', scalar = .true.)
      if (.not. fm_exists('/coupler_mod/types/' // trim(flux_type_test))) then
        call mpp_error(FATAL, trim(error_header) // ' Undefined flux_type given from field_table: ' // trim(flux_type_test))
      endif
    else
      call mpp_error(FATAL, trim(error_header) // ' Undefined flux_type given as argument to the subroutine: ' // trim(flux_type))
    endif
  endif

  if (implementation .eq. ' ') then
    call mpp_error(FATAL, trim(error_header) // ' Blank flux_type given')
  else
    if (fm_exists('/coupler_mod/types/' // trim(flux_type) // '/implementation/' // trim(implementation))) then
      call fm_util_set_value('implementation', implementation)

  !>       Check that the flux_type/implementation that we will use
  !!       (both possibly given from the field_table) is defined
      implementation_test = fm_util_get_string('implementation', scalar = .true.)
      if (.not. fm_exists('/coupler_mod/types/' // trim(flux_type_test) //  '/implementation/' // trim(implementation_test))) then
        if (flux_type .eq. flux_type_test) then
          if (implementation .eq. implementation_test) then
            call mpp_error(FATAL, trim(error_header) // ' Should not get here, as it is tested for above')
          else
            call mpp_error(FATAL, trim(error_header) // &
                 ' Undefined flux_type/implementation (implementation given from field_table): ' // &
                 trim(flux_type_test) // '/implementation/' // trim(implementation_test))
          endif
        else
          if (implementation .eq. implementation_test) then
            long_err_msg = 'Undefined flux_type/implementation (flux_type given from field_table): '
            long_err_msg = long_err_msg // trim(flux_type_test) // '/implementation/' // trim(implementation_test)
            call mpp_error(FATAL, trim(error_header) // long_err_msg)
          else
            long_err_msg = ' Undefined flux_type/implementation (both given from field_table): '
            long_err_msg = long_err_msg //  trim(flux_type_test) // '/implementation/' // trim(implementation_test)
            call mpp_error(FATAL, trim(error_header) // long_err_msg)
          endif
        endif
      endif
    else
      call mpp_error(FATAL, trim(error_header) // ' Undefined flux_type/implementation given as argument to the subroutine: ' // &
           trim(flux_type) // '/implementation/' // trim(implementation))
    endif
  endif

  if (present(atm_tr_index)) then
    call fm_util_set_value('atm_tr_index', atm_tr_index)
  else
    call fm_util_set_value('atm_tr_index', 0)
  endif

  if (present(mol_wt)) then
    call fm_util_set_value('mol_wt', mol_wt)
  else
    call fm_util_set_value('mol_wt', 0.0)
  endif

  if (present(ice_restart_file)) then
    call fm_util_set_value('ice_restart_file', ice_restart_file)
  else
    call fm_util_set_value('ice_restart_file', 'ice_coupler_fluxes.res.nc')
  endif

  if (present(ocean_restart_file)) then
    call fm_util_set_value('ocean_restart_file', ocean_restart_file)
  else
    call fm_util_set_value('ocean_restart_file', 'ocean_coupler_fluxes.res.nc')
  endif

  if (present(param)) then
    num_parameters = fm_util_get_integer('/coupler_mod/types/' // &
         trim(fm_util_get_string('flux_type', scalar = .true.)) // '/implementation/' // &
         trim(fm_util_get_string('implementation', scalar = .true.)) // '/num_parameters', scalar = .true.)
    length = min(size(param(:)),num_parameters)
    if ((length .ne. num_parameters) .and. (verbose >= 5)) then
      write (outunit,*) trim(note_header), ' Number of parameters provided for ', trim(name), ' does not match the'
      write (outunit,*) 'number of parameters required (', size(param(:)), ' != ', num_parameters, ').'
      write (outunit,*) 'This could be an error, or more likely is just a result of the implementation being'
      write (outunit,*) 'overridden by the field table input'
    endif
    if (length .gt. 0) then
      call fm_util_set_value('param', param(1:length), length)
    else
      call fm_util_set_value('param', 'null', index = 0)
    endif
  else
    call fm_util_set_value('param', 'null', index = 0)
  endif

  if (present(flag)) then
    call fm_util_set_value('flag', flag, size(flag(:)))
  else
    call fm_util_set_value('flag', .false., index = 0)
  endif

  flux_list = '/coupler_mod/types/' // trim(flux_type) // '/'

  if (present(units)) then
    call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'flux/name', index = ind_flux)) // '-units', units)
  else
    call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'flux/name', index = ind_flux)) // '-units', &
                           fm_util_get_string(trim(flux_list) // 'flux/units', index = ind_flux))
  endif

  do n = 1, fm_util_get_length(trim(flux_list) // 'flux/name')
    if (n .ne. ind_flux) then
      call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'flux/name', index = n)) // '-units', &
                             fm_util_get_string(trim(flux_list) // 'flux/units', index = n))
    endif
    call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'flux/name', index = n)) // '-long_name', &
                           fm_util_get_string(trim(flux_list) // 'flux/long_name', index = n))
  enddo  ! n

  do n = 1, fm_util_get_length(trim(flux_list) // 'atm/name')
    call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'atm/name', index = n)) // '-units', &
                           fm_util_get_string(trim(flux_list) // 'atm/units', index = n))
    call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'atm/name', index = n)) // '-long_name', &
                           fm_util_get_string(trim(flux_list) // 'atm/long_name', index = n))
  enddo  ! n

  do n = 1, fm_util_get_length(trim(flux_list) // 'ice/name')
    call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'ice/name', index = n)) // '-units', &
                           fm_util_get_string(trim(flux_list) // 'ice/units', index = n))
    call fm_util_set_value(trim(fm_util_get_string(trim(flux_list) // 'ice/name', index = n)) // '-long_name', &
                           fm_util_get_string(trim(flux_list) // 'ice/long_name', index = n))
  enddo  ! n

  !>       Reset the defaults for the fm_util_set_value calls.

  call fm_util_reset_good_name_list
  call fm_util_reset_no_overwrite
  call fm_util_reset_caller

  !>       Change back to the saved current list.

  if (.not. fm_change_list(current_list)) then
    call mpp_error(FATAL, trim(error_header) // ' Could not change back to ' // trim(current_list))
  endif

  !>       Check for any errors in the number of fields in this list.

  if (caller_str .eq. ' ') then
    caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'
  endif
  good_list => fm_util_get_string_array('/coupler_mod/GOOD/fluxes/' // trim(name) // '/good_list', &
       caller = caller_str)
  if (associated(good_list)) then
    call fm_util_check_for_bad_fields(trim(coupler_list), good_list, caller = caller_str)
    deallocate(good_list)
  else
    call mpp_error(FATAL, trim(error_header) // ' Empty "' // trim(name) // '" list')
  endif

  return

end function aof_set_coupler_flux

!> \brief Initialize gas flux structures
!
!! \throw FATAL, "Could not get number of fluxes"
!!     Number of gas fluxes is not a valid number
!! \throw NOTE, "No gas fluxes"
!!     No gas fluxes were found
!! \throw NOTE, "Processing [gas_fluxes%num_bcs] gas fluxes"
!!     Gas fluxes were found
!! \throw FATAL, "[name] is not a list"
!!     name needs to be a list, or typ is incorrectly defined
!! \throw FATAL, "Flux index, [ind] does not match array index, [n] for [name]"
!! \throw FATAL, "Problem changing to [name]"
!! \throw FATAL, "Undefined flux_type given for [name]: [gas_fluxes%bc(n)%flux_type]"
!! \throw FATAL, "Undefined implementation given for [name]: [gas_fluxes%bc(n)%flux_type]/implementation/[gas_fluxes%bc(n)%implementation]"
!! \throw FATAL, "No param for [name]: need [num_parameters]"
!! \throw FATAL, "Wrong number of param for [name]: [size(gas_fluxes%bc(n)%param(:))] given, need [num_parameters]"
!! \throw FATAL, "No params needed for [name] but has size of [size(gas_fluxes%bc(n)%param(:))]"
!! \throw FATAL, "Num_parameters is negative for [name]: [num_parameters]"
!! \throw FATAL, "No flag for [name]: need [num_flags]"
!! \throw FATAL, "Wrong number of flag for [name]: [size(gas_fluxes%bc(n)%flag(:))] given, need [num_flags]"
!! \throw FATAL, "No flags needed for [name] but has size of [size(gas_fluxes%bc(n)%flag(:))]"
!! \throw FATAL, "Num_flags is negative for [name]: [num_flags]"
!! \throw FATAL, "Problem dumping fluxes tracer tree"
!! \throw FATAL, "Number of fluxes does not match across the processors: [gas_fluxes%num_bcs] fluxes"
subroutine atmos_ocean_fluxes_init(gas_fluxes, gas_fields_atm, gas_fields_ice, verbosity)
#ifdef __APPLE__
! This directive is needed for compilation with -O2 using ifort 15.0.3 20150408 (Mac OSX)
! because otherwise the model crashes with error "malloc: pointer being freed was not allocated"
! at the end of this subroutine.
!DIR$ OPTIMIZE:0
#endif

  type(coupler_1d_bc_type), intent(inout) :: gas_fluxes !< Structure containing the gas fluxes between
                                                        !! the atmosphere and the ocean and parameters
                                                        !! related to the calculation of these fluxes.
                                                        !! The properties stored in this type are set
                                                        !! here, but the actual value arrays are set later.
  type(coupler_1d_bc_type), intent(inout) :: gas_fields_atm !< Structure containing atmospheric surface
                                                        !! variables that are used in the calculation
                                                        !! of the atmosphere-ocean gas fluxes.
                                                        !! The properties stored in this type are set
                                                        !! here, but the actual value arrays are set later.
  type(coupler_1d_bc_type), intent(inout) :: gas_fields_ice !< Structure containing ice-top and ocean
                                                        !! surface variables that are used in the
                                                        !! calculation of the atmosphere-ocean gas fluxes.
                                                        !! The properties stored in this type are set
                                                        !! here, but the actual value arrays are set later.
  integer,        optional, intent(in)    :: verbosity  !< A 0-9 integer indicating a level of verbosity.

  character(len=64), parameter    :: sub_name = 'atmos_ocean_fluxes_init'
  character(len=256), parameter   :: error_header = &
       '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: warn_header = &
       '==>Warning from ' // trim(mod_name) // '(' // trim(sub_name) // '):'
  character(len=256), parameter   :: note_header = &
       '==>Note from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

  integer                                 :: num_parameters
  integer                                 :: num_flags
  integer                                 :: n
  integer                                 :: m
  character(len=128)                      :: caller_str
  character(len=fm_type_name_len)         :: typ
  character(len=fm_field_name_len)        :: name
  integer                                 :: ind
  integer                                 :: outunit
  integer                                 :: total_fluxes
  character(len=8)                        :: string
  character(len=128)                      :: error_string
  character(len=128)                      :: flux_list
  logical, save                           :: initialized = .false.
  integer :: verbose ! An integer indicating the level of verbosity.

  if (initialized) return

  verbose = 5 ; if (present(verbosity)) verbose = verbosity

  !>       Write out the version of the file to the log file.
  call write_version_number(trim(mod_name), version)

  initialized = .true.
  outunit = stdout()
  if (verbose >= 9) then
    write (outunit,*)
    write (outunit,*) 'Dumping field manager tree'
    if (.not. fm_dump_list('/', recursive = .true.)) &
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping field manager tree')
  endif

  caller_str = trim(mod_name) // '(' // trim(sub_name) // ')'

  !>       Set other defaults for the fm_util_set_value routines.
  call fm_util_set_no_overwrite(.true.)
  call fm_util_set_caller(caller_str)

  !>       Determine the number of flux fields.
  gas_fluxes%num_bcs = fm_util_get_length('/coupler_mod/fluxes/')
  gas_fluxes%set = .true.
  gas_fields_atm%num_bcs = gas_fluxes%num_bcs ; gas_fields_atm%set = .true.
  gas_fields_ice%num_bcs = gas_fluxes%num_bcs ; gas_fields_ice%set = .true.
  if (gas_fluxes%num_bcs .lt. 0) then
    call mpp_error(FATAL, trim(error_header) // ' Could not get number of fluxes')
  elseif (gas_fluxes%num_bcs .eq. 0) then
    if (verbose >= 5) &
      write (outunit,*) trim(note_header), ' No gas fluxes'
    return
  else
    if (verbose >= 5) &
      write (outunit,*) trim(note_header), ' Processing ', gas_fluxes%num_bcs, ' gas fluxes'
  endif

  !       allocate the arrays
  allocate (gas_fluxes%bc(gas_fluxes%num_bcs))
  allocate (gas_fields_atm%bc(gas_fields_atm%num_bcs))
  allocate (gas_fields_ice%bc(gas_fields_ice%num_bcs))

  !>       Loop over the input fields, setting the values in the flux_type.
  n = 0
  do while (fm_loop_over_list('/coupler_mod/fluxes', name, typ, ind))

    if (typ .ne. 'list') then
      call mpp_error(FATAL, trim(error_header) // ' ' // trim(name) // ' is not a list')
    else

      n = n + 1  ! increment the array index

      if (n .ne. ind) then
        if (verbose >= 3) &
          write (outunit,*) trim(warn_header), ' Flux index, ', ind, &
               ' does not match array index, ', n, ' for ', trim(name)
      endif

  !>       Change list to the new flux.
      if (.not. fm_change_list('/coupler_mod/fluxes/' // trim(name))) then
        call mpp_error(FATAL, trim(error_header) // ' Problem changing to ' // trim(name))
      endif

  !>       Save and check the flux_type.
      gas_fluxes%bc(n)%flux_type = fm_util_get_string('flux_type', scalar = .true.)
      if (.not. fm_exists('/coupler_mod/types/' // trim(gas_fluxes%bc(n)%flux_type))) then
        call mpp_error(FATAL, trim(error_header) // ' Undefined flux_type given for ' // &
             trim(name) // ': ' // trim(gas_fluxes%bc(n)%flux_type))
      endif
      gas_fields_atm%bc(n)%flux_type = gas_fluxes%bc(n)%flux_type
      gas_fields_ice%bc(n)%flux_type = gas_fluxes%bc(n)%flux_type

  !>       Save and check the implementation.
      gas_fluxes%bc(n)%implementation = fm_util_get_string('implementation', scalar = .true.)
      if (.not. fm_exists('/coupler_mod/types/' // trim(gas_fluxes%bc(n)%flux_type) // &
           '/implementation/' // trim(gas_fluxes%bc(n)%implementation))) then
        call mpp_error(FATAL, trim(error_header) // ' Undefined implementation given for ' // &
             trim(name) // ': ' // trim(gas_fluxes%bc(n)%flux_type) // '/implementation/' // &
             trim(gas_fluxes%bc(n)%implementation))
      endif
      gas_fields_atm%bc(n)%implementation = gas_fluxes%bc(n)%implementation
      gas_fields_ice%bc(n)%implementation = gas_fluxes%bc(n)%implementation

  !>       Set the flux list name.
      flux_list = '/coupler_mod/types/' // trim(gas_fluxes%bc(n)%flux_type) // '/'

  !       allocate the arrays
      gas_fluxes%bc(n)%num_fields = fm_util_get_length(trim(flux_list) // 'flux/name')
      allocate (gas_fluxes%bc(n)%field(gas_fluxes%bc(n)%num_fields))
      gas_fields_atm%bc(n)%num_fields = fm_util_get_length(trim(flux_list) // 'atm/name')
      allocate (gas_fields_atm%bc(n)%field(gas_fields_atm%bc(n)%num_fields))
      gas_fields_ice%bc(n)%num_fields = fm_util_get_length(trim(flux_list) // 'ice/name')
      allocate (gas_fields_ice%bc(n)%field(gas_fields_ice%bc(n)%num_fields))

  !>       Save the name and generate unique field names for Flux, Ice and Atm.
      gas_fluxes%bc(n)%name = name
      do m = 1, fm_util_get_length(trim(flux_list) // 'flux/name')
        gas_fluxes%bc(n)%field(m)%name = trim(name) // "_" // fm_util_get_string(trim(flux_list) // 'flux/name', index = m)
        gas_fluxes%bc(n)%field(m)%override = .false.
        gas_fluxes%bc(n)%field(m)%mean     = .false.
      enddo

      gas_fields_atm%bc(n)%name = name
      do m = 1, fm_util_get_length(trim(flux_list) // 'atm/name')
        gas_fields_atm%bc(n)%field(m)%name = trim(name) // "_" // fm_util_get_string(trim(flux_list) // 'atm/name', index = m)
        gas_fields_atm%bc(n)%field(m)%override = .false.
        gas_fields_atm%bc(n)%field(m)%mean     = .false.
      enddo

      gas_fields_ice%bc(n)%name = name
      do m = 1, fm_util_get_length(trim(flux_list) // 'ice/name')
        gas_fields_ice%bc(n)%field(m)%name = trim(name) // "_" // fm_util_get_string(trim(flux_list) // 'ice/name', index = m)
        gas_fields_ice%bc(n)%field(m)%override = .false.
        gas_fields_ice%bc(n)%field(m)%mean     = .false.
      enddo

  !>       Save the units.
      do m = 1, fm_util_get_length(trim(flux_list) // 'flux/name')
        gas_fluxes%bc(n)%field(m)%units = &
             fm_util_get_string(trim(fm_util_get_string(trim(flux_list) // 'flux/name', index = m)) // '-units', scalar = .true.)
      enddo
      do m = 1, fm_util_get_length(trim(flux_list) // 'atm/name')
        gas_fields_atm%bc(n)%field(m)%units = &
             fm_util_get_string(trim(fm_util_get_string(trim(flux_list) // 'atm/name', index = m)) // '-units')
      enddo
      do m = 1, fm_util_get_length(trim(flux_list) // 'ice/name')
        gas_fields_ice%bc(n)%field(m)%units = &
             fm_util_get_string(trim(fm_util_get_string(trim(flux_list) // 'ice/name', index = m)) // '-units')
      enddo

  !>       Save the long names.
      do m = 1, fm_util_get_length(trim(flux_list) // 'flux/name')
        gas_fluxes%bc(n)%field(m)%long_name = &
             fm_util_get_string(trim(fm_util_get_string(trim(flux_list) // 'flux/name', index = m)) // '-long_name', scalar = .true.)
        gas_fluxes%bc(n)%field(m)%long_name = trim(gas_fluxes%bc(n)%field(m)%long_name) // ' for ' // name
      enddo
      do m = 1, fm_util_get_length(trim(flux_list) // 'atm/name')
        gas_fields_atm%bc(n)%field(m)%long_name = &
             fm_util_get_string(trim(fm_util_get_string(trim(flux_list) // 'atm/name', index = m)) // '-long_name')
        gas_fields_atm%bc(n)%field(m)%long_name = trim(gas_fields_atm%bc(n)%field(m)%long_name) // ' for ' // name
      enddo
      do m = 1, fm_util_get_length(trim(flux_list) // 'ice/name')
        gas_fields_ice%bc(n)%field(m)%long_name = &
             fm_util_get_string(trim(fm_util_get_string(trim(flux_list) // 'ice/name', index = m)) // '-long_name')
        gas_fields_ice%bc(n)%field(m)%long_name = trim(gas_fields_ice%bc(n)%field(m)%long_name) // ' for ' // name
      enddo

  !>       Save the atm_tr_index.
      gas_fluxes%bc(n)%atm_tr_index = fm_util_get_integer('atm_tr_index', scalar = .true.)

  !>       Save the molecular weight.
      gas_fluxes%bc(n)%mol_wt = fm_util_get_real('mol_wt', scalar = .true.)
      gas_fields_atm%bc(n)%mol_wt = gas_fluxes%bc(n)%mol_wt
      gas_fields_ice%bc(n)%mol_wt = gas_fluxes%bc(n)%mol_wt

  !>       Save the ice_restart_file.
      gas_fluxes%bc(n)%ice_restart_file = fm_util_get_string('ice_restart_file', scalar = .true.)
      gas_fields_atm%bc(n)%ice_restart_file = gas_fluxes%bc(n)%ice_restart_file
      gas_fields_ice%bc(n)%ice_restart_file = gas_fluxes%bc(n)%ice_restart_file

  !>       Save the ocean_restart_file.
      gas_fluxes%bc(n)%ocean_restart_file = fm_util_get_string('ocean_restart_file', scalar = .true.)
      gas_fields_atm%bc(n)%ocean_restart_file = gas_fluxes%bc(n)%ocean_restart_file
      gas_fields_ice%bc(n)%ocean_restart_file = gas_fluxes%bc(n)%ocean_restart_file

  !>       Save the params.
      gas_fluxes%bc(n)%param => fm_util_get_real_array('param')

  !>       Save the flags.
      gas_fluxes%bc(n)%flag => fm_util_get_logical_array('flag')

  !>       Perform some integrity checks.
      num_parameters = fm_util_get_integer(trim(flux_list) // 'implementation/' // &
           trim(gas_fluxes%bc(n)%implementation) // '/num_parameters', scalar = .true.)
      if (num_parameters .gt. 0) then
        if (.not. associated(gas_fluxes%bc(n)%param)) then
          write (error_string,'(a,i2)') ': need ', num_parameters
          call mpp_error(FATAL, trim(error_header) // ' No param for ' // trim(name) // trim(error_string))
        elseif (size(gas_fluxes%bc(n)%param(:)) .ne. num_parameters) then
          write (error_string,'(a,i2,a,i2)') ': ', size(gas_fluxes%bc(n)%param(:)), ' given, need ', num_parameters
          call mpp_error(FATAL, trim(error_header) // ' Wrong number of param for ' // trim(name) // trim(error_string))
        endif
      elseif (num_parameters .eq. 0) then
        if (associated(gas_fluxes%bc(n)%param)) then
          write (error_string,'(a,i3)') ' but has size of ', size(gas_fluxes%bc(n)%param(:))
          call mpp_error(FATAL, trim(error_header) // ' No params needed for ' // trim(name) // trim(error_string))
        endif
      else
        write (error_string,'(a,i2)') ': ', num_parameters
        call mpp_error(FATAL, trim(error_header) // 'Num_parameters is negative for ' // trim(name) // trim(error_string))
      endif
      num_flags = fm_util_get_integer(trim(flux_list) // '/num_flags', scalar = .true.)
      if (num_flags .gt. 0) then
        if (.not. associated(gas_fluxes%bc(n)%flag)) then
          write (error_string,'(a,i2)') ': need ', num_flags
          call mpp_error(FATAL, trim(error_header) // ' No flag for ' // trim(name) // trim(error_string))
        elseif (size(gas_fluxes%bc(n)%flag(:)) .ne. num_flags) then
          write (error_string,'(a,i2,a,i2)') ': ', size(gas_fluxes%bc(n)%flag(:)), ' given, need ', num_flags
          call mpp_error(FATAL, trim(error_header) // ' Wrong number of flag for ' // trim(name) // trim(error_string))
        endif
      elseif (num_flags .eq. 0) then
        if (associated(gas_fluxes%bc(n)%flag)) then
          write (error_string,'(a,i3)') ' but has size of ', size(gas_fluxes%bc(n)%flag(:))
          call mpp_error(FATAL, trim(error_header) // ' No flags needed for ' // trim(name) // trim(error_string))
        endif
      else
        write (error_string,'(a,i2)') ': ', num_flags
        call mpp_error(FATAL, trim(error_header) // 'Num_flags is negative for ' // trim(name) // trim(error_string))
      endif

  !>       Set some flags for this flux_type.

      gas_fluxes%bc(n)%use_atm_pressure = fm_util_get_logical(trim(flux_list) // '/use_atm_pressure')
      gas_fields_atm%bc(n)%use_atm_pressure = gas_fluxes%bc(n)%use_atm_pressure
      gas_fields_ice%bc(n)%use_atm_pressure = gas_fluxes%bc(n)%use_atm_pressure

      gas_fluxes%bc(n)%use_10m_wind_speed = fm_util_get_logical(trim(flux_list) // '/use_10m_wind_speed')
      gas_fields_atm%bc(n)%use_10m_wind_speed = gas_fluxes%bc(n)%use_10m_wind_speed
      gas_fields_ice%bc(n)%use_10m_wind_speed = gas_fluxes%bc(n)%use_10m_wind_speed

      gas_fluxes%bc(n)%pass_through_ice = fm_util_get_logical(trim(flux_list) // '/pass_through_ice')
      gas_fields_atm%bc(n)%pass_through_ice = gas_fluxes%bc(n)%pass_through_ice
      gas_fields_ice%bc(n)%pass_through_ice = gas_fluxes%bc(n)%pass_through_ice

    endif

  enddo ! while loop

  if (verbose >= 5) then
    write (outunit,*)
    write (outunit,*) 'Dumping fluxes tracer tree'
    if (.not. fm_dump_list('/coupler_mod/fluxes', recursive = .true.)) then
      call mpp_error(FATAL, trim(error_header) // ' Problem dumping fluxes tracer tree')
    endif
  endif

  !>       Check that the number of fluxes is the same on all processors
  !!       If they are, then the sum of the number of fluxes across all processors
  !!       should equal to the number of fluxes on each processor times the number of processors
  total_fluxes = gas_fluxes%num_bcs
  call mpp_sum(total_fluxes)
  if (total_fluxes .ne. mpp_npes() * gas_fluxes%num_bcs) then
    write (string, '(i4)') gas_fluxes%num_bcs
    call mpp_error(FATAL, trim(error_header) // &
         ' Number of fluxes does not match across the processors: ' // trim(string) // ' fluxes')
  endif

  !>       Reset the defaults for the fm_util_set_value calls.
  call fm_util_reset_no_overwrite
  call fm_util_reset_caller

end subroutine  atmos_ocean_fluxes_init

!> \brief Calculate the ocean gas fluxes. Units should be mol/m^2/s, upward flux is positive.
!
!! \throw FATAL, "Number of gas fluxes not zero"
!! \throw FATAL, "Lengths of flux fields do not match"
!! \throw FATAL, "Unknown implementation ([implementation]) for [name]"
!! \throw FATAL, "Lengths of flux fields do not match"
!! \throw FATAL, "Bad parameter ([gas_fluxes%bc(n)%param(1)]) for land_sea_runoff for [name]"
!! \throw FATAL, "Unknown flux type ([flux_type]) for [name]"
subroutine atmos_ocean_fluxes_calc(gas_fields_atm, gas_fields_ice, &
                                   gas_fluxes, seawater, tsurf, ustar, cd_m)
  type(coupler_1d_bc_type), intent(in)    :: gas_fields_atm !< Structure containing atmospheric surface
                                                        !! variables that are used in the calculation
                                                        !! of the atmosphere-ocean gas fluxes.
  type(coupler_1d_bc_type), intent(in)    :: gas_fields_ice !< Structure containing ice-top and ocean
                                                        !! surface variables that are used in the
                                                        !! calculation of the atmosphere-ocean gas fluxes.
  type(coupler_1d_bc_type), intent(inout) :: gas_fluxes !< Structure containing the gas fluxes between
                                                        !! the atmosphere and the ocean and parameters
                                                        !! related to the calculation of these fluxes.
  real, dimension(:),       intent(in)    :: seawater   !< 1 for the open water category, 0 if ice or land.
  real, dimension(:),       intent(in)    :: tsurf
  real, dimension(:), intent(in), optional :: ustar, cd_m

  character(len=64), parameter    :: sub_name = 'atmos_ocean_fluxes_calc'
  character(len=256), parameter   :: error_header = &
       '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

  integer                                 :: n
  integer                                 :: i
  integer                                 :: length
  real, dimension(:), allocatable         :: kw
  real, dimension(:), allocatable         :: cair
  character(len=128)                      :: error_string

  real, parameter :: epsln=1.0e-30
  real, parameter :: permeg=1.0e-6
  real :: salt !f1p should be dynamic

  !       Return if no fluxes to be calculated
  if (gas_fluxes%num_bcs .le. 0) return

  if (.not. associated(gas_fluxes%bc)) then
    if (gas_fluxes%num_bcs .ne. 0) then
      call mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
    else
      return
    endif
  endif

  do n = 1, gas_fluxes%num_bcs

  !       only do calculations if the flux has not been overridden
    if ( .not. gas_fluxes%bc(n)%field(ind_flux)%override) then

      if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux_generic') then

        length = size(gas_fluxes%bc(n)%field(1)%values(:))

        if (.not. allocated(kw)) then
          allocate( kw(length) )
          allocate ( cair(length) )
        elseif (size(kw(:)) .ne. length) then
          call mpp_error(FATAL, trim(error_header) // ' Lengths of flux fields do not match')
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
          do i = 1, length
            if (seawater(i) == 1.) then
               gas_fluxes%bc(n)%field(ind_kw)%values(i) = gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(ind_u10)%values(i)**2
              cair(i) = &
                   gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = gas_fluxes%bc(n)%field(ind_kw)%values(i) * &
                   sqrt(660. / (gas_fields_ice%bc(n)%field(ind_sc_no)%values(i) + epsln)) * &
                   (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
              gas_fluxes%bc(n)%field(ind_deltap)%values(i) = (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i)) / &
                   (gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * permeg + epsln)
            else
              gas_fluxes%bc(n)%field(ind_kw)%values(i) = 0.0
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
              gas_fluxes%bc(n)%field(ind_deltap)%values(i) = 0.0
              cair(i) = 0.0
            endif
          enddo  ! i

        elseif (gas_fluxes%bc(n)%implementation .eq. 'duce') then
           
          do i = 1, length
            if (seawater(i) == 1.) then

               gas_fluxes%bc(n)%field(ind_kw)%values(i) = &
                    gas_fields_atm%bc(n)%field(ind_u10)%values(i)/(770.+45.*gas_fluxes%bc(n)%param(1)**(1./3.)) * &
                    101325./(rdgas*wtmair*1e-3*tsurf(i)*gas_fields_ice%bc(n)%field(ind_alpha)%values(i))
!alpha: mol/m3/atm
              cair(i) = &
                   gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * 9.86923e-6
              cair(i) = max(cair(i),0.)
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = gas_fluxes%bc(n)%field(ind_kw)%values(i)  * (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i))
!              gas_fluxes%bc(n)%field(ind_flux0)%values(i) = gas_fluxes%bc(n)%field(ind_kw)%values(i) * max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i) = (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i)) / &
                   (gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * permeg + epsln)
            else
              gas_fluxes%bc(n)%field(ind_kw)%values(i) = 0.0
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
!              gas_fluxes%bc(n)%field(ind_flux0)%values(i) = 0.0
              gas_fluxes%bc(n)%field(ind_deltap)%values(i) = 0.0
              cair(i) = 0.0
            endif
          enddo  ! i


        elseif (gas_fluxes%bc(n)%implementation .eq. 'johnson') then

           !f1p: not sure how to pass salinity. For now, just force at 35.
           salt = 35.
           
          do i = 1, length
            if (seawater(i) == 1.) then

!               calc_kw(tk,p,u10,salt,h,vb,mw,ustar,cd_m) 
               gas_fluxes%bc(n)%field(ind_kw)%values(i) =                &
                    calc_kw(tsurf(i),                                    &
                    gas_fields_atm%bc(n)%field(ind_psurf)%values(i),     &
                    gas_fields_atm%bc(n)%field(ind_u10)%values(i),       &
                    salt,                                                &
                    101325./(rdgas*wtmair*1e-3*tsurf(i)*gas_fields_ice%bc(n)%field(ind_alpha)%values(i)),&
                    gas_fluxes%bc(n)%param(2), &
                    gas_fluxes%bc(n)%param(1), &
                    ustar=ustar(i), cd_m=cd_m(i))

              cair(i) = &
                   gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * 9.86923e-6
              cair(i) = max(cair(i),0.)
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = gas_fluxes%bc(n)%field(ind_kw)%values(i)  * (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i))
!              gas_fluxes%bc(n)%field(ind_flux0)%values(i) = gas_fluxes%bc(n)%field(ind_kw)%values(i) * max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.)
              gas_fluxes%bc(n)%field(ind_deltap)%values(i) = (max(gas_fields_ice%bc(n)%field(ind_csurf)%values(i),0.) - cair(i)) / &
                   (gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * permeg + epsln)
            else
              gas_fluxes%bc(n)%field(ind_kw)%values(i) = 0.0
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
!              gas_fluxes%bc(n)%field(ind_flux0)%values(i) = 0.0
              gas_fluxes%bc(n)%field(ind_deltap)%values(i) = 0.0
              cair(i) = 0.0
            endif
          enddo  ! i

        else
          call mpp_error(FATAL, ' Unknown implementation (' // trim(gas_fluxes%bc(n)%implementation) // &
               ') for ' // trim(gas_fluxes%bc(n)%name))
        endif

      elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_gas_flux') then

        length = size(gas_fluxes%bc(n)%field(1)%values(:))

        if (.not. allocated(kw)) then
          allocate( kw(length) )
          allocate ( cair(length) )
        elseif (size(kw(:)) .ne. length) then
          call mpp_error(FATAL, trim(error_header) // ' Lengths of flux fields do not match')
        endif

        if (gas_fluxes%bc(n)%implementation .eq. 'ocmip2_data') then
          do i = 1, length
            if (seawater(i) == 1.) then
              kw(i) = gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(ind_u10)%values(i)
              cair(i) = &
                   gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = kw(i) * &
                   (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
            else
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
              cair(i) = 0.0
              kw(i) = 0.0
            endif
          enddo  ! i
        elseif (gas_fluxes%bc(n)%implementation .eq. 'ocmip2') then
          do i = 1, length
            if (seawater(i) == 1.) then
              kw(i) = gas_fluxes%bc(n)%param(1) * gas_fields_atm%bc(n)%field(ind_u10)%values(i)**2
              cair(i) = &
                   gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(2)
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = kw(i) * &
                   (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
            else
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
              cair(i) = 0.0
              kw(i) = 0.0
            endif
          enddo  ! i
        elseif (gas_fluxes%bc(n)%implementation .eq. 'linear') then
          do i = 1, length
            if (seawater(i) == 1.) then
              kw(i) = gas_fluxes%bc(n)%param(1) * max(0.0, gas_fields_atm%bc(n)%field(ind_u10)%values(i) - gas_fluxes%bc(n)%param(2))
              cair(i) = &
                   gas_fields_ice%bc(n)%field(ind_alpha)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_pCair)%values(i) * &
                   gas_fields_atm%bc(n)%field(ind_psurf)%values(i) * gas_fluxes%bc(n)%param(3)
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = kw(i) * &
                   (gas_fields_ice%bc(n)%field(ind_csurf)%values(i) - cair(i))
            else
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
              cair(i) = 0.0
              kw(i) = 0.0
            endif
          enddo  ! i
        else
          call mpp_error(FATAL, ' Unknown implementation (' // trim(gas_fluxes%bc(n)%implementation) // &
               ') for ' // trim(gas_fluxes%bc(n)%name))
        endif

      elseif (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then

       cycle !air_sea_deposition is done in another subroutine

      elseif (gas_fluxes%bc(n)%flux_type .eq. 'land_sea_runoff') then

        if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
          write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
          call mpp_error(FATAL, ' Bad parameter (' // trim(error_string) // &
               ') for land_sea_runoff for ' // trim(gas_fluxes%bc(n)%name))
        endif

        length = size(gas_fluxes%bc(n)%field(1)%values(:))

        if (gas_fluxes%bc(n)%implementation .eq. 'river') then
          do i = 1, length
            if (seawater(i) == 1.) then
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = &
                   gas_fields_atm%bc(n)%field(ind_deposition)%values(i) / gas_fluxes%bc(n)%param(1)
            else
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
            endif
          enddo  ! i
        else
          call mpp_error(FATAL, ' Unknown implementation (' // trim(gas_fluxes%bc(n)%implementation) // &
               ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      else
        call mpp_error(FATAL, ' Unknown flux_type (' // trim(gas_fluxes%bc(n)%flux_type) // &
             ') for ' // trim(gas_fluxes%bc(n)%name))
      endif
    endif
  enddo  ! n

  if (allocated(kw)) then
    deallocate(kw)
    deallocate(cair)
  endif
    
end subroutine  atmos_ocean_fluxes_calc

!> \brief atmos_ocean_dep_fluxes_calc
!
!! \throw FATAL, "Number of gas fluxes not zero"
!! \throw FATAL, "atmos_ocean_dep_fluxes_calc: Bad parameter ([gas_fluxes%bc(n)%param(1)]) for air_sea_deposition for [gas_fluxes%bc(n)%name]"
!! \throw FATAL, "atmos_ocean_dep_fluxes_calc: Unknown implementation ([gas_fluxes%bc(n)%implementation] for [gas_fluxes%bc(n)%name]"
subroutine atmos_ocean_dep_fluxes_calc(gas_fields_atm, gas_fields_ice, &
                                       gas_fluxes, seawater)
  type(coupler_1d_bc_type), intent(in)    :: gas_fields_atm !< Structure containing atmospheric surface
                                                        !! variables that are used in the calculation
                                                        !! of the atmosphere-ocean gas fluxes.
  type(coupler_1d_bc_type), intent(in)    :: gas_fields_ice !< Structure containing ice-top and ocean
                                                        !! surface variables that are used in the
                                                        !! calculation of the atmosphere-ocean gas fluxes.
  type(coupler_1d_bc_type), intent(inout) :: gas_fluxes !< Structure containing the gas fluxes between
                                                        !! the atmosphere and the ocean and parameters
                                                        !! related to the calculation of these fluxes.
  real, dimension(:),       intent(in)    :: seawater   !< 1 for the open water category, 0 if ice or land.

  character(len=64), parameter    :: sub_name = 'atmos_ocean_dep_fluxes_calc'
  character(len=256), parameter   :: error_header = &
       '==>Error from ' // trim(mod_name) // '(' // trim(sub_name) // '):'

  integer                                 :: n
  integer                                 :: i
  integer                                 :: length
  real, dimension(:), allocatable         :: kw
  real, dimension(:), allocatable         :: cair
  character(len=128)                      :: error_string

  real, parameter :: epsln=1.0e-30
  real, parameter :: permeg=1.0e-6

  !       Return if no fluxes to be calculated

  if (gas_fluxes%num_bcs .le. 0) return

  if (.not. associated(gas_fluxes%bc)) then
    if (gas_fluxes%num_bcs .ne. 0) then
      call mpp_error(FATAL, trim(error_header) // ' Number of gas fluxes not zero')
    else
      return
    endif
  endif

  do n = 1, gas_fluxes%num_bcs
  !       only do calculations if the flux has not been overridden

    if ( .not. gas_fluxes%bc(n)%field(ind_flux)%override) then

      if (gas_fluxes%bc(n)%flux_type .eq. 'air_sea_deposition') then

        if (gas_fluxes%bc(n)%param(1) .le. 0.0) then
          write (error_string, '(1pe10.3)') gas_fluxes%bc(n)%param(1)
          call mpp_error(FATAL, 'atmos_ocean_dep_fluxes_calc: Bad parameter (' // trim(error_string) // &
               ') for air_sea_deposition for ' // trim(gas_fluxes%bc(n)%name))
        endif

        length = size(gas_fluxes%bc(n)%field(1)%values(:))

        if (gas_fluxes%bc(n)%implementation .eq. 'dry') then
          do i = 1, length
            if (seawater(i) == 1.) then
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = &
                   gas_fields_atm%bc(n)%field(ind_deposition)%values(i) / gas_fluxes%bc(n)%param(1)
            else
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
            endif
          enddo  ! i
        elseif (gas_fluxes%bc(n)%implementation .eq. 'wet') then
          do i = 1, length
            if (seawater(i) == 1.) then
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = &
                   gas_fields_atm%bc(n)%field(ind_deposition)%values(i) / gas_fluxes%bc(n)%param(1)
            else
              gas_fluxes%bc(n)%field(ind_flux)%values(i) = 0.0
            endif
          enddo  ! i
        else
          call mpp_error(FATAL, 'atmos_ocean_dep_fluxes_calc: Unknown implementation ('&
               // trim(gas_fluxes%bc(n)%implementation) // ') for ' // trim(gas_fluxes%bc(n)%name))
        endif
      else
        cycle
      endif
    endif

  enddo  ! n

end subroutine  atmos_ocean_dep_fluxes_calc

!from johnson ocean science 2010
function calc_kw(tk,p,u10,salt,h,vb,mw,ustar,cd_m)  result(kw)
  !from liss (1974)                                                                             
  !f = kg(cg - hcl)                                                                            
  !where cg and cl are the bulk gas and liquid concentrations and h is the henry's constant (h = cgs/cls, unitless)                                                                            
  !1/kg = 1/ka + h/kl    
  !tc: in c (temperature at surface)                                                   
  !p : in pa (pressure at surface)   
  !u10 : in m/s (wind speed at surface)                                                        
  !vb : molar volume                 
  !mw : molecular weight (g/mol)     
  real, intent(in)  :: tk,p,u10,salt,h,vb,mw
  real, intent(in), optional :: ustar,cd_m
  real :: ra,rl,tc,kw
  real, parameter :: epsln=1.0e-30

  tc = tk-273.15
  ra=1./max(h*saltout_correction(h,vb,salt)*calc_ka(tc,p,mw,vb,u10,ustar,cd_m),epsln)
  rl=1./max(calc_kl(tc,u10,salt,vb),epsln)
  kw = 1./max(ra+rl,epsln)
end function calc_kw

function calc_kg(tk,p,u10,salt,h,vb,mw,ustar,cd_m)  result(kg)
  !from liss (1974)                                                                             
  !f = kg(cg - hcl)                                                                            
  !where cg and cl are the bulk gas and liquid concentrations and h is the henry's constant (h = cgs/cls, unitless)                                                                            
  !1/kg = 1/ka + h/kl    
  !tc: in c (temperature at surface)                                                   
  !p : in pa (pressure at surface)   
  !u10 : in m/s (wind speed at surface)                                                        
  !vb : molar volume                 
  !mw : molecular weight (g/mol)     
  real, intent(in)  :: tk,p,u10,salt,h,vb,mw
  real, intent(in), optional :: ustar,cd_m
  real :: ra,rl,tc,kg
  tc = tk-273.15  
  ra=1./calc_ka(tc,p,mw,vb,u10,ustar,cd_m)
  rl=saltout_correction(h,vb,salt)*h/calc_kl(tc,u10,salt,vb)
  kg = 1./(ra+rl)
end function calc_kg

function calc_ka(t,p,mw,vb,u10,ustar,cd_m) result(k)
  real, intent(in) :: t,p, mw, vb, u10 !t in c, p in pa                               
  real, intent(in), optional :: ustar,cd_m
  real             :: sc, k  
  real             :: ustar_t,cd_m_t

  if (.not. present(ustar)) then
     !drag coefficient                                                                         
     cd_m_t     = 0.61e-3 +0.063e-3*u10
     !friction velocity        
     ustar_t  = u10*sqrt(cd_m_t)
  else
     cd_m_t=cd_m
     ustar_t=ustar
  end if
  sc=schmidt_g(t,p,mw,vb)    
  k=1e-3+ustar_t/(13.3*sqrt(sc)+1/sqrt(cd_m_t)-5.+log(sc)/(2.*vonkarm))  
end function calc_ka

function calc_kl(t,v,s,vb)result(k)
  real,intent(in)::t,s,v,vb
  real::k,sc,scco2
  !nightingale2000
  sc=schmidt_w(t,s,vb)
  k = (((0.222*v**2)+0.333*v)*(sc/600.)**(-0.5))/(100*3600)
end function calc_kl

function schmidt_g(t,p,mw,vb) result(sc)
  !schmidt number of the gas in the air                                                     
  !t in c                                                                                    
  real, intent(in) :: t,p,mw,vb
  real             :: sc,d,v
  d=d_air(t,p,mw,vb)
  v=v_air(t)
  sc = v / d
end function schmidt_g

function d_air(t,p,mw,vb) result(d)
  !from fuller 1966
  real, intent(in) :: t  !t in c                                                             
  real, intent(in) :: p  !t in pa                                                             
  real, intent(in) :: mw !mw in g/mol                                                     
  real, intent(in) :: vb !diffusion coefficient (cm3/mol)                                     
  real, parameter :: ma = 28.97d0 !mw air in g/mol                                            
  real, parameter :: va = 20.1d0  !cm3/mol (diffusion volume for air)                         
  real            :: d
  real            :: pa
  !convert p to atm                                                                          
  pa = 9.8692d-6*p
  d = 1d-3 * (t+273.15d0)**(1.75d0)*sqrt(1d0/ma + 1d0/mw)/(pa*(va**(1d0/3d0)+vb**(1d0/3d0))**2d0)
  !d is in cm2/s convert to m2/s                                                             
  d=d*1d-4
end function d_air

function p_air(t) result(p)
  !kinematic viscosity in air                                                                                                                                                                    
  real, intent(in) :: t
  real             :: p,sd_0,sd_1,sd_2,sd_3
  sd_0 = 1.293393662d0
  sd_1 = -5.538444326d-3
  sd_2 = 3.860201577d-5
  sd_3 = -5.2536065d-7
  p = sd_0+(sd_1*t)+(sd_2*t**2)+(sd_3*t**3)
end function p_air

function v_air(t) result(v)
  !kinematic viscosity in air (m2/s)                                                                                                                                                             
  real, intent(in) :: t
  real             :: v
  v = n_air(t)/p_air(t)
end function v_air

function d_hm(t,s,vb) result(d)
  real, intent(in) :: t,s,vb
  real             :: d, epsilonstar
  ! hayduk 1982                                                                              
  epsilonstar = (9.58d0/vb)-1.12d0
  d=1.25d-8*(vb**(-0.19d0)-0.292d0)*((t+273.15d0)**(1.52d0))*((n_sw(t,s))**epsilonstar)
end function d_hm

function schmidt_w(t,s,vb) result(sc)
  !schmidt number of the gas in the water                                                    
  real, intent(in) :: t,s,vb
  real             :: sc
  sc=2d0*v_sw(t,s)/(d_hm(t,s,vb)+d_wc(t,s,vb))
end function schmidt_w

function n_air(t) result(n)
  !dynamic viscosity in air                                                                  
  real, intent(in) :: t
  real             :: sv_0,sv_1,sv_2,sv_3,sv_4,n
  sv_0 = 1.715747771d-5
  sv_1 = 4.722402075d-8
  sv_2 = -3.663027156d-10
  sv_3 = 1.873236686d-12
  sv_4 = -8.050218737d-14
  ! in n.s/m^2 (pa.s)                                                                        
  n = sv_0+(sv_1*t)+(sv_2*t**2)+(sv_3*t**3)+(sv_4*t**4)
end function n_air

function v_sw(t,s) result(v)
  real, intent(in) :: t,s
  real             :: n,p,v
  n=n_sw(t,s)*1d-3
  p=p_sw(t,s)
  v = 1d4*n/p
end function  v_sw

function d_wc(t,s,vb) result(d)
  real, intent(in) :: t,s,vb
  real             :: d
  real, parameter  :: phi = 2.6
  !wilkie and chang 1955
  d = ((t+273.15)*7.4d-8*(phi*18.01)**0.5)/((n_sw(t,s))*(vb**0.6))
end function d_wc

function p_sw(t,s) result(p)
  !density of sea water                                                                           !millero and poisson (1981)                                                               
  real, intent(in) :: t,s
  real             :: p, a, b, c
  a = 0.824493d0-(4.0899d-3*t)+(7.6438d-5*(t**2))-(8.2467d-7*(t**3))+(5.3875d-9*(t**4))
  b = -5.72466d-3+(1.0277d-4*t)-(1.6546d-6*(t**2))
  c = 4.8314d-4
  ! density of pure water
  p = 999.842594d0+(6.793952d-2*t)-(9.09529d-3*(t**2))+(1.001685d-4*(t**3))-(1.120083d-6*(t**4))+(6.536332d-9*(t**5))
  !salinity correction
  p = (p+(a*s)+(b*(s**(1.5d0)))+(c*s))
end function p_sw

function n_sw(t,s) result(n)
  !dynamic viscosity                                   
  !laliberte 2007                                      
  real :: n
  real, intent(in) :: t,s !temperature (c) and salinity
  !salt in the order nacl,kcl,cacl2,mgcl2,mgso4      
  real, parameter :: mass_fraction(5) = (/ 0.798,0.022,0.033,0.047,0.1 /)
  real, parameter :: v1(5) = (/ 16.22 , 6.4883, 32.028, 24.032, 72.269/)
  real, parameter :: v2(5) = (/ 1.3229 , 1.3175, 0.78792, 2.2694, 2.2238/)
  real, parameter :: v3(5) = (/ 1.4849 , -0.7785, -1.1495,  3.7108, 6.6037/)
  real, parameter :: v4(5) = (/ 0.0074691 , 0.09272, 0.0026995,  0.021853, 0.0079004/)
  real, parameter :: v5(5) = (/ 30.78 , -1.3, 780860., -1.1236, 3340.1/)
  real, parameter :: v6(5) = (/ 2.0583 , 2.0811, 5.8442,0.14474, 6.1304/)

  real :: n_0,ln_n_m,w_i_ln_n_i_tot,ni,w_i_tot,w_i
  integer :: i
  w_i_tot=0
  w_i_ln_n_i_tot=0
  do i=1,5
     w_i = mass_fraction(i)*s/1000
     w_i_tot = w_i+w_i_tot
  enddo
  do i=1,5
     w_i = mass_fraction(i)*s/1000
     ni = (exp(((v1(i)*w_i_tot**v2(i))+v3(i))/((v4(i)*t) + 1)))/((v5(i)*(w_i_tot**v6(i)))+1)
     w_i_ln_n_i_tot = w_i_ln_n_i_tot + (w_i*log(ni))
  enddo
  n_0 = (t+246)/(137.37+(5.2842*t)+(0.05594*(t**2)))
  ln_n_m = (1-w_i_tot)*log(n_0)+w_i_ln_n_i_tot
  n = exp(ln_n_m)
end function n_sw

function saltout_correction(kh,vb,salt) result(C)
  real, intent(in) :: Kh,vb,salt
  real*8 :: T,log_kh,theta2    
  real :: theta,C
  log_kh = log(kh)
  theta = (7.3353282561828962e-04 + (3.3961477466551352e-05*log_kh) + (-2.4088830102075734e-06*(log_kh)**2) + (1.5711393120941302e-07*(log_kh)**3))*log(vb)
  C = 10**(theta*salt)    
end function saltout_correction


end module  atmos_ocean_fluxes_mod
