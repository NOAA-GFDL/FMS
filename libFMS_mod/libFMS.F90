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
!> @file
!! @author Ryan Mulhall 2/21
!! @brief A convenience file to use all public FMS routines
!! @description Imports all FMS modules so that when this module is 
!! used, all public FMS routines are avaiable for use

module libFMS_mod
  !> import all FMS modules
  !! for ones we might not need
  use fms_affinity_mod
  use amip_interp_mod
  use astronomy_mod
  use axis_utils_mod !! might be deprecated for 2
  use axis_utils2_mod 
  use block_control_mod
  use column_diagnostics_mod
  use constants_mod
  use coupler_types_mod !! coupler ones needed?
  use ensemble_manager_mod !!
  use atmos_ocean_fluxes_mod !!
  use data_override_mod
  use get_grid_version_mpp_mod !!
  use get_grid_version_fms2io_mod !
  use diag_integral_mod
  use diag_manager_mod
  use diag_table_mod
  use diag_util_mod
  use diag_output_mod
  use diag_axis_mod
  use diag_data_mod
  use drifters_mod
  use drifters_comm_mod
  use drifters_core_mod
  use drifters_input_mod
  use drifters_io_mod
  use cloud_interpolator_mod 
  use xgrid_mod
  use stock_constants_mod
  use field_manager_mod
  use fm_util_mod
  use fms2_io_mod
  use fms_io_utils_mod
  use fms_netcdf_domain_io_mod
  use netcdf_io_mod
  use fms_netcdf_unstructured_domain_io_mod
  use blackboxio
  use fms_mod
  use horiz_interp_type_mod
  use horiz_interp_bicubic_mod
  use horiz_interp_bilinear_mod
  use horiz_interp_conserve_mod
  use horiz_interp_spherical_mod
  use horiz_interp_mod
  use interpolator_mod
  use memutils_mod
  use monin_obukhov_mod
  use monin_obukhov_inter
  use mosaic_mod
  use mosaic2_mod
  use grid_mod
  use gradient_mod
  use mpp_mod
  use mpp_parameter_mod
  use mpp_data_mod
  use mpp_utilities_mod
  use mpp_memutils_mod
  use mpp_efp_mod
  use mpp_domains_mod
  use mpp_io_mod
  use platform_mod
  use mersennetwister_mod
  use random_numbers_mod
  use sat_vapor_pres_k_mod
  use sat_vapor_pres_mod
  use time_interp_mod
  use time_interp_external_mod
  use time_interp_external2_mod
  use time_manager_mod
  use get_cal_time_mod
  use gaussian_topog_mod
  use topography_mod
  use tracer_manager_mod
  use tridiagonal_mod

  implicit none
  !! 
  logical :: usingLibFMS = .true.

end module libFMS_mod
