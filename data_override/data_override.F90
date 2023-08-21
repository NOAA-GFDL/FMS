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

#include "data_override_impl_r4.fh"
#include "data_override_impl_r8.fh"

!> @defgroup data_override_mod data_override_mod
!> @ingroup data_override
!! @brief Routines to get data in a file whose path is described in a user-provided data_table
!! and do spatial and temporal interpolation if necessary to convert data to model's grid and time.
!! @author Z. Liang, M.J. Harrison, M. Winton
!!
!! Before using @ref data_override a data_table must be created with the following entries:
!! gridname, fieldname_code, fieldname_file, file_name, ongrid, factor.
!!
!! More explainations about data_table entries can be found in the source code (defining data_type)
!!
!! If user wants to override fieldname_code with a const, set fieldname_file in data_table = ""
!! and factor = const
!!
!! If user wants to override fieldname_code with data from a file, set fieldname_file = name in
!! the netCDF data file, factor then will be for unit conversion (=1 if no conversion required)
!!
!! A field can be overriden globally (by default) or users can specify one or two regions in which
!! data_override will take place, field values outside the region will not be affected.

module data_override_mod
  use data_override_r4
  use data_override_r8
  use platform_mod, only: r4_kind, r8_kind
  use mpp_mod, only: mpp_error, FATAL
  use mpp_domains_mod, only : domain2d, domainUG
  use time_manager_mod, only: time_type

implicit none
private

!> Interface for inserting and interpolating data into a file
!! for a model's grid and time. Data path must be described in
!! a user-provided data_table, see @ref data_override_mod "module description"
!! for more information.
!> @ingroup data_override_mod
interface data_override
     module procedure data_override_0d_r4
     module procedure data_override_0d_r8
     module procedure data_override_2d_r4
     module procedure data_override_2d_r8
     module procedure data_override_3d_r4
     module procedure data_override_3d_r8
end interface

!> Version of @ref data_override for unstructured grids
!> @ingroup data_override_mod
interface data_override_UG
     module procedure data_override_UG_1d_r4
     module procedure data_override_UG_1d_r8
     module procedure data_override_UG_2d_r4
     module procedure data_override_UG_2d_r8
end interface

integer :: atm_mode = 0 !> Atmosphere mode - possible values are 0 (uninitialized), r4_kind, or r8_kind
integer :: ocn_mode = 0 !> Ocean mode - possible values are 0 (uninitialized), r4_kind, or r8_kind
integer :: lnd_mode = 0 !> Land mode - possible values are 0 (uninitialized), r4_kind, or r8_kind
integer :: ice_mode = 0 !> Ice mode - possible values are 0 (uninitialized), r4_kind, or r8_kind

!> @addtogroup data_override_mod
!> @{

public :: data_override_init, data_override, data_override_unset_domains
public :: data_override_UG

contains

!> @brief Initialize either data_override_r4 or data_override_r8
!! Users should call data_override_init before calling data_override
!!
!! This subroutine should be called in coupler_init after
!! (ocean/atmos/land/ice)_model_init have been called.
!!
!! data_override_init can be called more than once. In one call the user can pass
!! up to 4 domains of component models. At least one domain must be present in
!! any call. The real precision of initialized domains can be specified via the
!! optional mode argument. If no mode is specified, r8_kind is assumed. Mixed mode
!! operation can be accomplished via multiple calls to data_override_init with
!! different mode arguments.
!!
!! Data_table is initialized with default values in DATA_OVERRIDE_INIT_IMPL_. Users should
!! provide "real" values that will override the default values. Real values can be
!! specified in either data_table or data_table.yaml. Each line of data_table contains one
!! data_entry. Items of data_entry are comma-separated.
subroutine data_override_init(Atm_domain_in, Ocean_domain_in, Ice_domain_in, Land_domain_in, Land_domainUG_in, mode)
  type (domain2d), intent(in), optional :: Atm_domain_in !< Atmosphere domain
  type (domain2d), intent(in), optional :: Ocean_domain_in !< Ocean domain
  type (domain2d), intent(in), optional :: Ice_domain_in !< Ice domain
  type (domain2d), intent(in), optional :: Land_domain_in !< Land domain
  type(domainUG) , intent(in), optional :: Land_domainUG_in !< Land domain, unstructured grid
  integer, intent(in), optional :: mode !< Real precision of initialized domains. Possible values are r4_kind or r8_kind.
  integer :: mode_selector

  if (present(mode)) then
    mode_selector = mode
  else
    mode_selector = r8_kind
  endif

  select case (mode_selector)
    case (r4_kind)
      call data_override_init_r4(Atm_domain_in, Ocean_domain_in, Ice_domain_in, Land_domain_in, Land_domainUG_in)
    case (r8_kind)
      call data_override_init_r8(Atm_domain_in, Ocean_domain_in, Ice_domain_in, Land_domain_in, Land_domainUG_in)
    case default
      call mpp_error(FATAL, "data_override_init: unsupported mode argument")
  end select

  if (present(Atm_domain_in))   atm_mode = mode_selector
  if (present(Ocean_domain_in)) ocn_mode = mode_selector
  if (present(Ice_domain_in))   ice_mode = mode_selector
  if (present(Land_domain_in))  lnd_mode = mode_selector
end subroutine data_override_init

!> @brief Unset domains that had previously been set for use by data_override.
!!
!! This subroutine deallocates any data override domains that have been set.
subroutine data_override_unset_domains(unset_Atm, unset_Ocean, &
                                      unset_Ice, unset_Land, must_be_set)
  logical, intent(in), optional :: unset_Atm, unset_Ocean, unset_Ice, unset_Land !< Set to true to unset the
                                                                                 !! respective domain
  logical, intent(in), optional :: must_be_set !< Set to false to suppress the error when attempting to unset
                                               !! an uninitialized domain
  logical :: fail_if_not_set

  fail_if_not_set = .true. ; if (present(must_be_set)) fail_if_not_set = must_be_set

  if (present(unset_Atm)) then ; if (unset_Atm) then
    select case (atm_mode)
      case (r4_kind)
        call data_override_unset_atm_r4
      case (r8_kind)
        call data_override_unset_atm_r8
      case default
        if (fail_if_not_set) call mpp_error(FATAL, &
          "data_override_unset_domains: attempted to unset an Atm_domain that has not been set.")
    end select
    atm_mode = 0
  endif ; endif
  if (present(unset_Ocean)) then ; if (unset_Ocean) then
    select case (ocn_mode)
      case (r4_kind)
        call data_override_unset_ocn_r4
      case (r8_kind)
        call data_override_unset_ocn_r8
      case default
        if (fail_if_not_set) call mpp_error(FATAL, &
          "data_override_unset_domains: attempted to unset an Ocn_domain that has not been set.")
    end select
    ocn_mode = 0
  endif ; endif
  if (present(unset_Land)) then ; if (unset_Land) then
    select case (lnd_mode)
      case (r4_kind)
        call data_override_unset_lnd_r4
      case (r8_kind)
        call data_override_unset_lnd_r8
      case default
        if (fail_if_not_set) call mpp_error(FATAL, &
          "data_override_unset_domains: attempted to unset an Land_domain that has not been set.")
    end select
    lnd_mode = 0
  endif ; endif
  if (present(unset_Ice)) then ; if (unset_Ice) then
    select case (ice_mode)
      case (r4_kind)
        call data_override_unset_ice_r4
      case (r8_kind)
        call data_override_unset_ice_r8
      case default
        if (fail_if_not_set) call mpp_error(FATAL, &
          "data_override_unset_domains: attempted to unset an Ice_domain that has not been set.")
    end select
    ice_mode = 0
  endif ; endif
end subroutine data_override_unset_domains

end module data_override_mod
!> @}
! close documentation grouping
