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
!> @defgroup FMS FMS
!> @ingroup libfms
!> @brief A convenience module to use any FMS routines, functions, values
!> @author Ryan Mulhall
!!
!! @date 2/2021
!!
!! Imports all supported FMS modules so that any public interfaces,
!! variables or routines can be used via this module. Excludes mpp_io modules
!! and routines. Overloaded type operators/assignments cannot be imported individually
!! (ie. `use fms, only: OPERATOR(*)` includes any defined '*' operators within FMS).
!!
!! Renaming scheme:
!!   Routines and variables: fms_<module_name>_routine_name
!!   Types:                  FmsModuleNameTypeName
!!
!! Exceptions (mainly for rep:
!!   - Parameter values are kept their original names
!!   - If module name is already included (like in init routines) only fms prefix will be added.
!!   - Similarly if theres a redundant module name included already included it will not be repeated
!!     (ie. mpp_update_domains => fms_mpp_domains_update_domains)
!!   - Override interfaces for operators and assignment are provided
!!
!! Remappings due to name conflicts:
!!
!!           ZERO from interpolator_mod(mpp_parameter)  => INTERPOLATOR_ZERO
!!
!!           version from fms_mod                       => version_FMS
!!
!! Not in this module:
!!
!!                     axis_utils_mod, fms_io_mod, time_interp_external_mod
!!                     get_grid_version_mpp_mod, mpp_io_mod, mosaic_mod, &
!!                     fms_mod(partial, old io excluded), drifters modules
!!                     constants_mod (FMSconstants should be used externally)
!!                     grid_mod, mosaic_mod
!!
!! A full list of supported interfaces and public types intended for use via
!! this module is provided in the [supported_interfaces.md](../../supported_interfaces.md)
!! file.

!> @file
!> @brief File for @ref FMS

!> @addtogroup FMS
!> @{
module fms

  !> import each FMS module's public routines/functions, interfaces, and variables
  !! done explicitly to avoid including any unwanted/depracated routines/modules

  !> affinity
  use fms_affinity_mod, only: fms_affinity_init, fms_affinity_get, &
                              fms_affinity_set

  !> amip_interp
  use amip_interp_mod, only: fms_amip_interp_init         => amip_interp_init, &
                             fms_amip_interp_get_amip_sst             => get_amip_sst, &
                             fms_amip_interp_get_amip_ice             => get_amip_ice, &
                             fms_amip_interp_new          => amip_interp_new, &
                             fms_amip_interp_del          => amip_interp_del, &
                             FmsAmipInterp_type => amip_interp_type, &
                             assignment(=), &
                             fms_amip_interp_i_sst        => i_sst, &
                             fms_amip_interp_j_sst        => j_sst, &
                             fms_amip_interp_sst_ncep     => sst_ncep, &
                             fms_amip_interp_sst_anom     => sst_anom, &
                             fms_amip_interp_forecast_mode=> forecast_mode, &
                             fms_amip_interp_use_ncep_sst => use_ncep_sst
  !> astronomy
  use astronomy_mod, only: fms_astronomy_init                   => astronomy_init, &
                           fms_astronomy_get_period             => get_period, &
                           fms_astronomy_set_period             => set_period, &
                           fms_astronomy_set_orbital_parameters => set_orbital_parameters, &
                           fms_astronomy_get_orbital_parameters => get_orbital_parameters, &
                           fms_astronomy_set_ref_date_of_ae     => set_ref_date_of_ae, &
                           fms_astronomy_get_ref_date_of_ae     => get_ref_date_of_ae, &
                           fms_astronomy_diurnal_solar          => diurnal_solar, &
                           fms_astronomy_daily_mean_solar       => daily_mean_solar, &
                           fms_astronomy_annual_mean_solar      => annual_mean_solar,  &
                           fms_astronomy_end                    => astronomy_end, &
                           fms_astronomy_universal_time         => universal_time, &
                           fms_astronomy_orbital_time           => orbital_time

  !> axis_utils
  use axis_utils2_mod, only: fms_axis_utils2_get_axis_cart         => get_axis_cart, &
                             fms_axis_utils2_get_axis_modulo       => get_axis_modulo, &
                             fms_axis_utils2_lon_in_range          => lon_in_range, &
                             fms_axis_utils2_tranlon               => tranlon, &
                             fms_axis_utils2_frac_index            => frac_index, &
                             fms_axis_utils2_nearest_index         => nearest_index, &
                             fms_axis_utils2_interp_1d             => interp_1d, &
                             fms_axis_utils2_get_axis_modulo_times => get_axis_modulo_times, &
                             fms_axis_utils2_axis_edges            => axis_edges

  !>block_control
  use block_control_mod, only: FmsBlockControl_type                   => block_control_type, &
                               fms_block_control_define_blocks        => define_blocks, &
                               fms_block_control_define_blocks_packed => define_blocks_packed

  !> column_diagnostics
  use column_diagnostics_mod, only: fms_column_diagnostics_init        => column_diagnostics_init, &
                                    fms_column_diagnostics_initialize_diagnostic_columns => &
                                                                          initialize_diagnostic_columns, &
                                    fms_column_diagnostics_header      => column_diagnostics_header, &
                                    fms_column_diagnostics_close_units => close_column_diagnostics_units

  !> coupler
  use coupler_types_mod, only: fms_coupler_types_init             => coupler_types_init, &
                               fms_coupler_type_copy              => coupler_type_copy, &
                               fms_coupler_type_spawn             => coupler_type_spawn, &
                               fms_coupler_type_set_diags         => coupler_type_set_diags, &
                               fms_coupler_type_write_chksums     => coupler_type_write_chksums, &
                               fms_coupler_type_send_data         => coupler_type_send_data, &
                               fms_coupler_type_data_override     => coupler_type_data_override, &
                               fms_coupler_type_register_restarts => coupler_type_register_restarts, &
                               fms_coupler_type_restore_state     => coupler_type_restore_state, &
                               fms_coupler_type_increment_data    => coupler_type_increment_data, &
                               fms_coupler_type_rescale_data      => coupler_type_rescale_data, &
                               fms_coupler_type_copy_data         => coupler_type_copy_data, &
                               fms_coupler_type_redistribute_data => coupler_type_redistribute_data, &
                               fms_coupler_type_destructor        => coupler_type_destructor, &
                               fms_coupler_type_initialized       => coupler_type_initialized, &
                               fms_coupler_type_extract_data      => coupler_type_extract_data, &
                               fms_coupler_type_set_data          => coupler_type_set_data, &
                               fms_coupler_type_copy_1d_2d        => coupler_type_copy_1d_2d, &
                               fms_coupler_type_copy_1d_3d        => coupler_type_copy_1d_3d, &
                               FmsCoupler3dValues_type         => coupler_3d_values_type, &
                               FmsCoupler3dField_type          => coupler_3d_field_type, &
                               FmsCoupler3dBC_type             => coupler_3d_bc_type, &
                               FmsCoupler2dValues_type         => coupler_2d_values_type, &
                               FmsCoupler2dField_type          => coupler_2d_field_type, &
                               FmsCoupler2dBC_type             => coupler_2d_bc_type, &
                               FmsCoupler1dValues_type         => coupler_1d_values_type, &
                               FmsCoupler1dField_type          => coupler_1d_field_type, &
                               FmsCoupler1dBC_type             => coupler_1d_bc_type, &
                               fms_coupler_ind_pcair                      => ind_pcair, &
                               fms_coupler_ind_u10                        => ind_u10, &
                               fms_coupler_ind_psurf                      => ind_psurf, &
                               fms_coupler_ind_alpha                      => ind_alpha, &
                               fms_coupler_ind_csurf                      => ind_csurf, &
                               fms_coupler_ind_sc_no                      => ind_sc_no, &
                               fms_coupler_ind_flux                       => ind_flux, &
                               fms_coupler_ind_deltap                     => ind_deltap, &
                               fms_coupler_ind_kw                         => ind_kw, &
                               fms_coupler_ind_flux0                      => ind_flux0, &
                               fms_coupler_ind_deposition                 => ind_deposition,&
                               fms_coupler_ind_runoff                     => ind_runoff
  use ensemble_manager_mod, only: fms_ensemble_manager_init     => ensemble_manager_init, &
                               fms_ensemble_manager_get_ensemble_id              => get_ensemble_id, &
                               fms_ensemble_manager_get_ensemble_size            => get_ensemble_size, &
                               fms_ensemble_manager_get_ensemble_pelist          => get_ensemble_pelist, &
                               fms_ensemble_manager_ensemble_pelist_setup        => ensemble_pelist_setup, &
                               fms_ensemble_manager_get_ensemble_filter_pelist   => get_ensemble_filter_pelist
  use atmos_ocean_fluxes_mod, only: fms_atmos_ocean_fluxes_init => atmos_ocean_fluxes_init, &
                               fms_atmos_ocean_type_fluxes_init => atmos_ocean_type_fluxes_init, &
                               fms_atmos_ocean_fluxes_set_coupler_flux         => aof_set_coupler_flux

  !> data_override
  use data_override_mod, only: fms_data_override_init          => data_override_init, &
                               fms_data_override               => data_override, &
                               fms_data_override_unset_domains => data_override_unset_domains, &
                               fms_data_override_UG            => data_override_UG

  !> diag_integral
  use diag_integral_mod, only: fms_diag_integral_init       => diag_integral_init, &
                               fms_diag_integral_field_init => diag_integral_field_init, &
                               fms_sum_diag_integral_field  => sum_diag_integral_field, &
                               fms_diag_integral_output     => diag_integral_output, &
                               fms_diag_integral_end        => diag_integral_end

  !> diag_manager
  !! includes imports from submodules made public
  use diag_manager_mod, only: fms_diag_init => diag_manager_init, &
                              fms_diag_send_data => send_data, &
                              fms_diag_send_tile_averaged_data => send_tile_averaged_data, &
                              fms_diag_end => diag_manager_end, &
                              fms_diag_register_diag_field => register_diag_field, &
                              fms_diag_register_static_field => register_static_field, &
                              fms_diag_axis_init => diag_axis_init, &
                              fms_diag_get_base_time => get_base_time, &
                              fms_diag_get_base_date => get_base_date, &
                              fms_diag_need_data => need_data, &
                              DIAG_ALL, &
                              DIAG_OCEAN, &
                              DIAG_OTHER, &
                              fms_get_date_dif => get_date_dif, &
                              DIAG_SECONDS,&
                              DIAG_MINUTES, &
                              DIAG_HOURS, &
                              DIAG_DAYS, &
                              DIAG_MONTHS, &
                              DIAG_YEARS, &
                              fms_diag_get_global_att => get_diag_global_att, &
                              fms_diag_set_global_att => set_diag_global_att, &
                              fms_diag_field_add_attribute => diag_field_add_attribute, &
                              fms_diag_field_add_cell_measures => diag_field_add_cell_measures, &
                              fms_diag_get_field_id => get_diag_field_id, &
                              fms_diag_axis_add_attribute => diag_axis_add_attribute, &
                              fms_diag_grid_init => diag_grid_init, &
                              fms_diag_grid_end => diag_grid_end, &
                              fms_diag_set_time_end => diag_manager_set_time_end, &
                              fms_diag_send_complete => diag_send_complete, &
                              fms_diag_send_complete_instant => diag_send_complete_instant, &
                              DIAG_FIELD_NOT_FOUND, &
                              CMOR_MISSING_VALUE, &
                              null_axis_id

  !> exchange
  use xgrid_mod, only: FmsXgridXmap_type => xmap_type, &
                       fms_xgrid_setup_xmap => setup_xmap, &
                       fms_xgrid_set_frac_area => set_frac_area, &
                       fms_xgrid_put_to_xgrid => put_to_xgrid, &
                       fms_xgrid_get_from_xgrid => get_from_xgrid, &
                       fms_xgrid_count => xgrid_count, &
                       fms_xgrid_some => some, &
                       fms_xgrid_conservation_check => conservation_check, &
                       fms_xgrid_init => xgrid_init, &
                       AREA_ATM_SPHERE, AREA_OCN_SPHERE, AREA_ATM_MODEL, AREA_OCN_MODEL, &
                       fms_xgrid_get_ocean_model_area_elements => get_ocean_model_area_elements, &
                       FmsXgridGridBox_type => grid_box_type, &
                       fms_xgrid_get_xmap_grid_area => get_xmap_grid_area, &
                       fms_xgrid_put_to_xgrid_ug => put_to_xgrid_ug, &
                       fms_xgrid_get_from_xgrid_ug => get_from_xgrid_ug, &
                       fms_xgrid_set_frac_area_ug => set_frac_area_ug, &
                       FIRST_ORDER, SECOND_ORDER, &
                       fms_xgrid_stock_move_ug => stock_move_ug, &
                       fms_xgrid_stock_move => stock_move, &
                       FmsXgridStock_type => stock_type, &
                       fms_xgrid_stock_print => stock_print, &
                       fms_xgrid_get_index_range => get_index_range, &
                       fms_xgrid_stock_integrate_2d => stock_integrate_2d
  use stock_constants_mod, only: NELEMS, ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT, &
                       ISTOCK_TOP, ISTOCK_BOTTOM, ISTOCK_SIDE, &
                       fms_stock_constants_stocks_file => stocks_file, &
                       fms_stock_constants_stocks_report => stocks_report, &
                       fms_stocks_report_init => stocks_report_init, &
                       fms_stocks_set_init_time => stocks_set_init_time, &
                       fms_stock_constants_atm_stock => atm_stock, &
                       fms_stock_constants_ocn_stock => ocn_stock, &
                       fms_stock_constants_lnd_stock => lnd_stock, &
                       fms_stock_constants_ice_stock => ice_stock

  !> field manager
  use field_manager_mod, only: fms_field_manager_init => field_manager_init, &
                         fms_field_manager_end => field_manager_end, &
                         fms_field_manager_find_field_index => find_field_index, &
                         fms_field_manager_get_field_info => get_field_info, &
                         fms_field_manager_get_field_method => get_field_method, &
                         fms_field_manager_get_field_methods => get_field_methods, &
                         fms_field_manager_parse => parse, &
                         fms_field_manager_fm_change_list => fm_change_list, &
                         fms_field_manager_fm_change_root => fm_change_root, &
                         fms_field_manager_fm_dump_list => fm_dump_list, &
                         fms_field_manager_fm_exists => fm_exists, &
                         fms_field_manager_fm_get_index => fm_get_index, &
                         fms_field_manager_fm_get_current_list => fm_get_current_list, &
                         fms_field_manager_fm_get_length => fm_get_length, &
                         fms_field_manager_fm_get_type => fm_get_type, &
                         fms_field_manager_fm_get_value => fm_get_value, &
                         fms_field_manager_fm_init_loop => fm_init_loop, &
                         fms_field_manager_fm_loop_over_list => fm_loop_over_list, &
                         fms_field_manager_fm_new_list => fm_new_list, &
                         fms_field_manager_fm_new_value => fm_new_value, &
                         fms_field_manager_fm_reset_loop => fm_reset_loop, &
                         fms_field_manager_fm_return_root => fm_return_root, &
                         fms_field_manager_fm_modify_name => fm_modify_name, &
                         fms_field_manager_fm_query_method => fm_query_method, &
                         fms_field_manager_fm_find_methods => fm_find_methods, &
                         fms_field_manager_fm_copy_list => fm_copy_list, &
                         fms_field_manager_fm_field_name_len => fm_field_name_len, &
                         fms_field_manager_fm_path_name_len => fm_path_name_len, &
                         fms_field_manager_fm_string_len => fm_string_len, &
                         fms_field_manager_fm_type_name_len => fm_type_name_len, &
                         NUM_MODELS, NO_FIELD, MODEL_ATMOS, MODEL_OCEAN, MODEL_LAND, MODEL_ICE, MODEL_COUPLER, &
                         FmsFieldManagerMethod_type          => method_type, &
                         FmsFieldManagerMethodShort_type     => method_type_short, &
                         FmsFieldManagerMethodVeryShort_type => method_type_very_short, &
                         FmsFieldManagerListIterator_type => fm_list_iter_type, &
                         fms_field_manager_default_method => default_method
  use fm_util_mod, only: fms_fm_util_start_namelist => fm_util_start_namelist, &
                         fms_fm_util_end_namelist => fm_util_end_namelist, &
                         fms_fm_util_check_for_bad_fields => fm_util_check_for_bad_fields, &
                         fms_fm_util_set_caller => fm_util_set_caller, &
                         fms_fm_util_reset_caller => fm_util_reset_caller, &
                         fms_fm_util_set_no_overwrite => fm_util_set_no_overwrite, &
                         fms_fm_util_reset_no_overwrite => fm_util_reset_no_overwrite, &
                         fms_fm_util_set_good_name_list => fm_util_set_good_name_list, &
                         fms_fm_util_reset_good_name_list => fm_util_reset_good_name_list, &
                         fms_fm_util_get_length => fm_util_get_length, &
                         fms_fm_util_get_integer => fm_util_get_integer, &
                         fms_fm_util_get_logical => fm_util_get_logical, &
                         fms_fm_util_get_real => fm_util_get_real, &
                         fms_fm_util_get_string => fm_util_get_string, &
                         fms_fm_util_get_integer_array => fm_util_get_integer_array, &
                         fms_fm_util_get_logical_array => fm_util_get_logical_array, &
                         fms_fm_util_get_real_array => fm_util_get_real_array, &
                         fms_fm_util_get_string_array => fm_util_get_string_array, &
                         fms_fm_util_set_value => fm_util_set_value, &
                         fms_fm_util_set_value_integer_array => fm_util_set_value_integer_array, &
                         fms_fm_util_set_value_logical_array => fm_util_set_value_logical_array, &
                         fms_fm_util_set_value_real_array => fm_util_set_value_real_array, &
                         fms_fm_util_set_value_string_array => fm_util_set_value_string_array, &
                         fms_fm_util_set_value_integer => fm_util_set_value_integer, &
                         fms_fm_util_set_value_logical => fm_util_set_value_logical, &
                         fms_fm_util_set_value_real => fm_util_set_value_real, &
                         fms_fm_util_set_value_string => fm_util_set_value_string, &
                         fms_fm_util_get_index_list => fm_util_get_index_list, &
                         fms_fm_util_get_index_string => fm_util_get_index_string, &
                         fms_fm_util_default_caller => fm_util_default_caller

  !> fms2_io
  !! TODO need to see opinions on these
  !! not sure if we need fms_ prefix for routines
  !! types do not follow our typical naming convention(no _type and uses camel case instead of _ spacing)
  use fms2_io_mod, only: unlimited, FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                         FmsNetcdfUnstructuredDomainFile_t, &
                         Valid_t, &
                         fms2_io_open_file => open_file, &
                         fms2_io_open_virtual_file => open_virtual_file, &
                         fms2_io_close_file => close_file, &
                         fms2_io_register_axis => register_axis, &
                         fms2_io_register_field => register_field, &
                         fms2_io_register_restart_field => register_restart_field, &
                         fms2_io_write_data => write_data, &
                         fms2_io_read_data => read_data, &
                         fms2_io_write_restart => write_restart, &
                         fms2_io_write_new_restart => write_new_restart, &
                         fms2_io_read_restart => read_restart, &
                         fms2_io_read_new_restart => read_new_restart, &
                         fms2_io_global_att_exists => global_att_exists, &
                         fms2_io_variable_att_exists => variable_att_exists, &
                         fms2_io_register_global_attribute => register_global_attribute, &
                         fms2_io_register_variable_attribute => register_variable_attribute, &
                         fms2_io_get_global_attribute => get_global_attribute, &
                         fms2_io_get_variable_attribute => get_variable_attribute, &
                         fms2_io_get_num_dimensions => get_num_dimensions, &
                         fms2_io_get_dimension_names => get_dimension_names, &
                         fms2_io_dimension_exists => dimension_exists, &
                         fms2_io_is_dimension_unlimited => is_dimension_unlimited, &
                         fms2_io_get_dimension_size => get_dimension_size, &
                         fms2_io_get_num_variables =>   get_num_variables, &
                         fms2_io_get_variable_names =>   get_variable_names, &
                         fms2_io_variable_exists => variable_exists, &
                         fms2_io_get_variable_num_dimensions => get_variable_num_dimensions, &
                         fms2_io_get_variable_dimension_names => get_variable_dimension_names, &
                         fms2_io_get_variable_size =>   get_variable_size, &
                         fms2_io_get_compute_domain_dimension_indices => get_compute_domain_dimension_indices, &
                         fms2_io_get_global_io_domain_indices => get_global_io_domain_indices, &
                         fms2_io_get_valid => get_valid, &
                         fms2_io_is_valid => is_valid, &
                         fms2_io_get_unlimited_dimension_name => get_unlimited_dimension_name, &
                         fms2_io_get_variable_unlimited_dimension_index => get_variable_unlimited_dimension_index, &
                         fms2_io_file_exists => file_exists, &
                         fms2_io_compressed_start_and_count => compressed_start_and_count, &
                         fms2_io_get_variable_sense =>   get_variable_sense, &
                         fms2_io_get_variable_missing =>   get_variable_missing, &
                         fms2_io_get_variable_units =>   get_variable_units, &
                         fms2_io_get_time_calendar =>   get_time_calendar, &
                         fms2_io_open_check =>   open_check, &
                         fms2_io_is_registered_to_restart =>  is_registered_to_restart, &
                         fms2_io_check_if_open =>   check_if_open, &
                         fms2_io_set_fileobj_time_name => set_fileobj_time_name, &
                         fms2_io_is_dimension_registered => is_dimension_registered, &
                         fms2_io_fms2_io_init =>   fms2_io_init, &
                         fms2_io_get_mosaic_tile_grid =>   get_mosaic_tile_grid, &
                         fms2_io_write_restart_bc =>   write_restart_bc, &
                         fms2_io_read_restart_bc =>   read_restart_bc, &
                         fms2_io_get_filename_appendix =>   get_filename_appendix, &
                         fms2_io_set_filename_appendix =>   set_filename_appendix, &
                         fms2_io_get_instance_filename =>   get_instance_filename, &
                         fms2_io_nullify_filename_appendix =>   nullify_filename_appendix, &
                         fms2_io_ascii_read =>   ascii_read, &
                         fms2_io_get_mosaic_tile_file =>   get_mosaic_tile_file, &
                         fms2_io_parse_mask_table => parse_mask_table
  ! used via fms2_io
  ! fms_io_utils_mod, fms_netcdf_domain_io_mod, netcdf_io_mod, &
  ! fms_netcdf_unstructured_domain_io_mod, blackboxio

  !> fms
  !! routines that don't conflict with fms2_io
  use fms_mod, only: fms_init, fms_end, error_mesg, fms_error_handler, &
                     check_nml_error, &
                     fms_monotonic_array => monotonic_array, fms_string_array_index => string_array_index, &
                     fms_clock_flag_default => clock_flag_default, fms_print_memory_usage => print_memory_usage, &
                     fms_write_version_number => write_version_number

  !> horiz_interp
  use horiz_interp_mod, only: fms_horiz_interp => horiz_interp, fms_horiz_interp_new => horiz_interp_new, &
                              fms_horiz_interp_del => horiz_interp_del, fms_horiz_interp_init => horiz_interp_init, &
                              fms_horiz_interp_end => horiz_interp_end
  use horiz_interp_type_mod, only: FmsHorizInterp_type => horiz_interp_type, &
                              assignment(=), CONSERVE, BILINEAR, SPHERICA, BICUBIC, &
                              fms_horiz_interp_type_stats => stats
  !! used via horiz_interp
  ! horiz_interp_bicubic_mod, horiz_interp_bilinear_mod
  ! horiz_interp_conserve_mod, horiz_interp_spherical_mod

  !> interpolator
  use interpolator_mod, only: fms_interpolator_init => interpolator_init, &
                              fms_interpolator => interpolator, &
                              fms_interpolate_type_eq => interpolate_type_eq, &
                              fms_interpolator_obtain_interpolator_time_slices => obtain_interpolator_time_slices, &
                              fms_interpolator_unset_interpolator_time_flag => unset_interpolator_time_flag, &
                              fms_interpolator_end => interpolator_end, &
                              fms_interpolator_init_clim_diag => init_clim_diag, &
                              fms_interpolator_query_interpolator => query_interpolator, &
                              FmsInterpolate_type => interpolate_type, &
                              CONSTANT, INTERP_WEIGHTED_P, INTERP_LINEAR_P, INTERP_LOG_P, &
                              FMS_INTERPOLATOR_ZERO=>ZERO, & !! conflicts with mpp_parameter's ZERO
                              fms_interpolator_read_data=>read_data

  !> memutils
  use memutils_mod, only: fms_memutils_init => memutils_init, &
                          fms_memutils_print_memuse_stats => print_memuse_stats

  !> monin_obukhov
  use monin_obukhov_mod, only: fms_monin_obukhov_init => monin_obukhov_init, &
                               fms_monin_obukhov_end => monin_obukhov_end, &
                               fms_monin_obukhov_mo_drag => mo_drag, &
                               fms_monin_obukhov_mo_profile => mo_profile, &
                               fms_monin_obukhov_mo_diff => mo_diff, &
                               fms_monin_obukhov_stable_mix => stable_mix
  use monin_obukhov_inter, only: fms_monin_obukhov_inter_diff => monin_obukhov_diff, &
                                 fms_monin_obukhov_inter_drag_1d => monin_obukhov_drag_1d, &
                                 fms_monin_obukhov_inter_solve_zeta => monin_obukhov_solve_zeta, &
                                 fms_monin_obukhov_inter_derivative_t => monin_obukhov_derivative_t, &
                                 fms_monin_obukhov_inter_derivative_m => monin_obukhov_derivative_m, &
                                 fms_monin_obukhov_inter_profile_1d => monin_obukhov_profile_1d, &
                                 fms_monin_obukhov_inter_integral_m => monin_obukhov_integral_m, &
                                 fms_monin_obukhov_inter_integral_tq => monin_obukhov_integral_tq, &
                                 fms_monin_obukhov_inter_stable_mix => monin_obukhov_stable_mix

  !> mosaic
  use mosaic2_mod, only: fms_mosaic2_get_mosaic_ntiles => get_mosaic_ntiles, &
                      fms_mosaic2_get_mosaic_ncontacts => get_mosaic_ncontacts, &
                      fms_mosaic2_get_mosaic_grid_sizes => get_mosaic_grid_sizes, &
                      fms_mosaic2_get_mosaic_contact => get_mosaic_contact, &
                      fms_mosaic2_get_mosaic_xgrid_size => get_mosaic_xgrid_size, &
                      fms_mosaic2_get_mosaic_xgrid => get_mosaic_xgrid, &
                      fms_mosaic2_calc_mosaic_grid_area => calc_mosaic_grid_area, &
                      fms_mosaic2_calc_mosaic_grid_great_circle_area => calc_mosaic_grid_great_circle_area, &
                      fms_mosaic2_is_inside_polygon => is_inside_polygon, &
                      fms_mosaic2_get_mosaic_tile_grid => get_mosaic_tile_grid !overloaded in fms2_io
  use grid2_mod, only: fms_grid2_get_grid_ntiles => get_grid_ntiles, &
                       fms_grid2_get_grid_size => get_grid_size, &
                       fms_grid2_get_grid_cell_centers => get_grid_cell_centers, &
                       fms_grid2_get_grid_cell_vertices => get_grid_cell_vertices, &
                       fms_grid2_get_grid_cell_Area => get_grid_cell_Area, &
                       fms_grid2_get_grid_comp_area => get_grid_comp_area, &
                       fms_grid2_define_cube_mosaic => define_cube_mosaic, &
                       fms_grid2_get_great_circle_algorithm => get_great_circle_algorithm, &
                       fms_grid2_grid_init => grid_init, &
                       fms_grid2_end => grid_end
  use gradient_mod, only: fms_gradient_cubic => gradient_cubic, &
                          fms_gradient_calc_cubic_grid_info => calc_cubic_grid_info

  !> mpp
  use mpp_mod, only: fms_mpp_stdin => stdin, &
                     fms_mpp_stdout => stdout, &
                     fms_mpp_stderr => stderr, &
                     fms_mpp_stdlog => stdlog, &
                     fms_mpp_lowercase => lowercase, &
                     fms_mpp_uppercase => uppercase, &
                     fms_mpp_error => mpp_error, &
                     fms_mpp_error_state => mpp_error_state, &
                     fms_mpp_set_warn_level => mpp_set_warn_level, &
                     fms_mpp_sync => mpp_sync, &
                     fms_mpp_sync_self => mpp_sync_self, &
                     fms_mpp_set_stack_size => mpp_set_stack_size, &
                     fms_mpp_pe => mpp_pe, &
                     fms_mpp_npes => mpp_npes, &
                     fms_mpp_root_pe => mpp_root_pe, &
                     fms_mpp_set_root_pe => mpp_set_root_pe, &
                     fms_mpp_declare_pelist => mpp_declare_pelist, &
                     fms_mpp_get_current_pelist => mpp_get_current_pelist, &
                     fms_mpp_set_current_pelist => mpp_set_current_pelist, &
                     fms_mpp_get_current_pelist_name => mpp_get_current_pelist_name, &
                     fms_mpp_clock_id => mpp_clock_id, &
                     fms_mpp_clock_set_grain => mpp_clock_set_grain, &
                     fms_mpp_record_timing_data => mpp_record_timing_data, &
                     fms_mpp_get_unit => get_unit, &
                     fms_mpp_read_ascii_file => read_ascii_file, &
                     fms_mpp_read_input_nml => read_input_nml, &
                     fms_mpp_clock_begin => mpp_clock_begin, &
                     fms_mpp_clock_end => mpp_clock_end, &
                     fms_mpp_get_ascii_file_num_lines => get_ascii_file_num_lines, &
                     fms_mpp_record_time_start => mpp_record_time_start, &
                     fms_mpp_record_time_end => mpp_record_time_end, &
                     fms_mpp_chksum => mpp_chksum, &
                     fms_mpp_max => mpp_max, &
                     fms_mpp_min => mpp_min, &
                     fms_mpp_sum => mpp_sum, &
                     fms_mpp_transmit => mpp_transmit, &
                     fms_mpp_send => mpp_send, &
                     fms_mpp_recv => mpp_recv, &
                     fms_mpp_sum_ad => mpp_sum_ad, &
                     fms_mpp_broadcast => mpp_broadcast, &
                     fms_mpp_init => mpp_init, &
                     fms_mpp_exit => mpp_exit, &
                     fms_mpp_gather => mpp_gather, &
                     fms_mpp_scatter => mpp_scatter, &
                     fms_mpp_alltoall => mpp_alltoall, &
                     FmsMpp_type => mpp_type, &
                     FmsMpp_byte => mpp_byte, &
                     fms_mpp_type_create => mpp_type_create, &
                     fms_mpp_type_free => mpp_type_free, &
                     fms_mpp_input_nml_file => input_nml_file
  use mpp_parameter_mod,only:MAXPES, MPP_VERBOSE, MPP_DEBUG, ALL_PES, ANY_PE, NULL_PE, &
                     NOTE, WARNING, FATAL, MPP_WAIT, MPP_READY, MAX_CLOCKS, &
                     MAX_EVENT_TYPES, MAX_EVENTS, MPP_CLOCK_SYNC, MPP_CLOCK_DETAILED, &
                     CLOCK_COMPONENT, CLOCK_SUBCOMPONENT, CLOCK_MODULE_DRIVER, &
                     CLOCK_MODULE, CLOCK_ROUTINE, CLOCK_LOOP, CLOCK_INFRA, MAX_BINS, &
                     EVENT_ALLREDUCE, EVENT_BROADCAST, EVENT_RECV, EVENT_SEND, EVENT_WAIT, &
                     EVENT_ALLTOALL, EVENT_TYPE_CREATE, EVENT_TYPE_FREE, &
                     DEFAULT_TAG, COMM_TAG_1, COMM_TAG_2, COMM_TAG_3, COMM_TAG_4, &
                     COMM_TAG_5,  COMM_TAG_6,  COMM_TAG_7,  COMM_TAG_8, &
                     COMM_TAG_9,  COMM_TAG_10, COMM_TAG_11, COMM_TAG_12, &
                     COMM_TAG_13, COMM_TAG_14, COMM_TAG_15, COMM_TAG_16, &
                     COMM_TAG_17, COMM_TAG_18, COMM_TAG_19, COMM_TAG_20, &
                     MPP_FILL_INT, MPP_FILL_FLOAT, MPP_FILL_DOUBLE, &
                     GLOBAL_DATA_DOMAIN, CYCLIC_GLOBAL_DOMAIN, BGRID_NE, &
                     BGRID_SW, CGRID_NE, CGRID_SW, DGRID_NE, DGRID_SW, &
                     FOLD_WEST_EDGE, FOLD_EAST_EDGE, FOLD_SOUTH_EDGE, &
                     FOLD_NORTH_EDGE, WUPDATE, EUPDATE, SUPDATE, NUPDATE, &
                     XUPDATE, YUPDATE, BITWISE_EXACT_SUM, NON_BITWISE_EXACT_SUM, &
                     MPP_DOMAIN_TIME, WEST, EAST, SOUTH, NORTH, SCALAR_BIT, SCALAR_PAIR, &
                     BITWISE_EFP_SUM, NORTH_EAST, SOUTH_EAST, SOUTH_WEST, NORTH_WEST, &
                     AGRID, GLOBAL, CYCLIC, DOMAIN_ID_BASE, CENTER, CORNER, &
                     MAX_DOMAIN_FIELDS, MAX_TILES, ZERO, NINETY, MINUS_NINETY, &
                     ONE_HUNDRED_EIGHTY, NONBLOCK_UPDATE_TAG, EDGEUPDATE, EDGEONLY, &
                     NONSYMEDGEUPDATE, NONSYMEDGE
  ! this should really only be used internally
  !use mpp_data_mod, only: stat, mpp_stack, ptr_stack, status, ptr_status, sync, &
  !                   ptr_sync, mpp_from_pe, ptr_from, remote_Data_loc, &
  !                   ptr_remote, mpp_domains_stack, ptr_domains_stack, &
  !                   mpp_domains_stack_nonblock, ptr_domains_stack_nonblock
  use mpp_utilities_mod, only: fms_mpp_utilities_array_global_min_max => mpp_array_global_min_max
  use mpp_memutils_mod, only: fms_mpp_memutils_print_memuse_stats => mpp_print_memuse_stats, &
                              fms_mpp_memutils_mem_dump => mpp_mem_dump, &
                              fms_mpp_memutils_memuse_begin => mpp_memuse_begin, &
                              fms_mpp_memutils_memuse_end => mpp_memuse_end
  use mpp_efp_mod, only: fms_mpp_efp_reproducing_sum => mpp_reproducing_sum, &
                         fms_mpp_efp_list_sum_across_PEs => mpp_efp_list_sum_across_PEs, &
                         fms_mpp_efp_plus => mpp_efp_plus, &
                         fms_mpp_efp_minus => mpp_efp_minus, &
                         fms_mpp_efp_to_real => mpp_efp_to_real, &
                         fms_mpp_efp_real_to_efp => mpp_real_to_efp, &
                         fms_mpp_efp_real_diff => mpp_efp_real_diff, &
                         operator(+), operator(-), assignment(=), &
                         fms_mpp_efp_query_overflow_error => mpp_query_efp_overflow_error, &
                         fms_mpp_efp_reset_overflow_error => mpp_reset_efp_overflow_error, &
                         FmsMppEfp_type => mpp_efp_type
  use mpp_domains_mod, only: FmsMppDomains_axis_spec => domain_axis_spec, &
                             FmsMppDomain1D => domain1D, &
                             FmsMppDomain2D => domain2D, &
                             FmsMppDomainCommunicator2D => DomainCommunicator2D, &
                             FmsMppDomainsNestDomain_type => nest_domain_type, &
                             FmsMppDomainsGroupUpdate_type => mpp_group_update_type, &
                             fms_mpp_domains_domains_set_stack_size => mpp_domains_set_stack_size, &
                             fms_mpp_domains_get_compute_domain => mpp_get_compute_domain, &
                             fms_mpp_domains_get_compute_domains => mpp_get_compute_domains, &
                             fms_mpp_domains_get_data_domain => mpp_get_data_domain, &
                             fms_mpp_domains_get_global_domain => mpp_get_global_domain, &
                             fms_mpp_domains_get_domain_components => mpp_get_domain_components, &
                             fms_mpp_domains_get_layout => mpp_get_layout, &
                             fms_mpp_domains_get_pelist => mpp_get_pelist, &
                             operator(.EQ.), operator(.NE.), &
                             fms_mpp_domains_domain_is_symmetry => mpp_domain_is_symmetry, &
                             fms_mpp_domains_domain_is_initialized => mpp_domain_is_initialized, &
                             fms_mpp_domains_get_neighbor_pe => mpp_get_neighbor_pe, &
                             fms_mpp_domains_nullify_domain_list => mpp_nullify_domain_list, &
                             fms_mpp_domains_set_compute_domain => mpp_set_compute_domain, &
                             fms_mpp_domains_set_data_domain => mpp_set_data_domain, &
                             fms_mpp_domains_set_global_domain => mpp_set_global_domain, &
                             fms_mpp_domains_get_memory_domain => mpp_get_memory_domain, &
                             fms_mpp_domains_get_domain_shift => mpp_get_domain_shift, &
                             fms_mpp_domains_domain_is_tile_root_pe => mpp_domain_is_tile_root_pe, &
                             fms_mpp_domains_get_tile_id => mpp_get_tile_id, &
                             fms_mpp_domains_get_domain_extents => mpp_get_domain_extents, &
                             fms_mpp_domains_get_current_ntile => mpp_get_current_ntile, &
                             fms_mpp_domains_get_ntile_count => mpp_get_ntile_count, &
                             fms_mpp_domains_get_tile_list => mpp_get_tile_list, &
                             fms_mpp_domains_get_tile_npes => mpp_get_tile_npes, &
                             fms_mpp_domains_get_domain_root_pe => mpp_get_domain_root_pe, &
                             fms_mpp_domains_get_tile_pelist => mpp_get_tile_pelist, &
                             fms_mpp_domains_get_tile_compute_domains => mpp_get_tile_compute_domains, &
                             fms_mpp_domains_get_num_overlap => mpp_get_num_overlap, &
                             fms_mpp_domains_get_overlap => mpp_get_overlap, &
                             fms_mpp_domains_get_io_domain => mpp_get_io_domain, &
                             fms_mpp_domains_get_domain_pe => mpp_get_domain_pe, &
                             fms_mpp_domains_get_domain_tile_root_pe => mpp_get_domain_tile_root_pe, &
                             fms_mpp_domains_get_domain_name => mpp_get_domain_name, &
                             fms_mpp_domains_get_io_domain_layout => mpp_get_io_domain_layout, &
                             fms_mpp_domains_copy_domain => mpp_copy_domain, &
                             fms_mpp_domains_set_domain_symmetry => mpp_set_domain_symmetry, &
                             fms_mpp_domains_get_update_pelist => mpp_get_update_pelist, &
                             fms_mpp_domains_get_update_size => mpp_get_update_size, &
                             fms_mpp_domains_get_domain_npes => mpp_get_domain_npes, &
                             fms_mpp_domains_get_domain_pelist => mpp_get_domain_pelist, &
                             fms_mpp_domains_clear_group_update => mpp_clear_group_update, &
                             fms_mpp_domains_group_update_initialized => mpp_group_update_initialized, &
                             fms_mpp_domains_group_update_is_set => mpp_group_update_is_set, &
                             fms_mpp_domains_get_global_domains => mpp_get_global_domains, &
                             fms_mpp_domains_global_field => mpp_global_field, &
                             fms_mpp_domains_global_max => mpp_global_max, &
                             fms_mpp_domains_global_min => mpp_global_min, &
                             fms_mpp_domains_global_sum => mpp_global_sum, &
                             fms_mpp_domains_global_sum_tl => mpp_global_sum_tl, &
                             fms_mpp_domains_global_sum_ad => mpp_global_sum_ad, &
                             fms_mpp_domains_broadcast_domain => mpp_broadcast_domain, &
                             fms_mpp_domains_init => mpp_domains_init, &
                             fms_mpp_domains_exit => mpp_domains_exit, &
                             fms_mpp_domains_redistribute => mpp_redistribute, &
                             fms_mpp_domains_update_domains => mpp_update_domains, &
                             fms_mpp_domains_check_field => mpp_check_field, &
                             fms_mpp_domains_start_update_domains => mpp_start_update_domains, &
                             fms_mpp_domains_complete_update_domains => mpp_complete_update_domains, &
                             fms_mpp_domains_create_group_update => mpp_create_group_update, &
                             fms_mpp_domains_do_group_update => mpp_do_group_update, &
                             fms_mpp_domains_start_group_update => mpp_start_group_update, &
                             fms_mpp_domains_complete_group_update => mpp_complete_group_update, &
                             fms_mpp_domains_reset_group_update_field => mpp_reset_group_update_field, &
                             fms_mpp_domains_update_nest_fine => mpp_update_nest_fine, &
                             fms_mpp_domains_update_nest_coarse => mpp_update_nest_coarse, &
                             fms_mpp_domains_get_boundary => mpp_get_boundary, &
                             fms_mpp_domains_update_domains_ad => mpp_update_domains_ad, &
                             fms_mpp_domains_get_boundary_ad => mpp_get_boundary_ad, &
                             fms_mpp_domains_pass_SG_to_UG => mpp_pass_SG_to_UG, &
                             fms_mpp_domains_pass_UG_to_SG => mpp_pass_UG_to_SG, &
                             fms_mpp_domains_define_layout => mpp_define_layout, &
                             fms_mpp_domains_define_domains => mpp_define_domains, &
                             fms_mpp_domains_modify_domain => mpp_modify_domain, &
                             fms_mpp_domains_define_mosaic => mpp_define_mosaic, &
                             fms_mpp_domains_define_mosaic_pelist => mpp_define_mosaic_pelist, &
                             fms_mpp_domains_define_null_domain => mpp_define_null_domain, &
                             fms_mpp_domains_mosaic_defined => mpp_mosaic_defined, &
                             fms_mpp_domains_define_io_domain => mpp_define_io_domain, &
                             fms_mpp_domains_deallocate_domain => mpp_deallocate_domain, &
                             fms_mpp_domains_compute_extent => mpp_compute_extent, &
                             fms_mpp_domains_compute_block_extent => mpp_compute_block_extent, &
                             fms_mpp_domains_define_unstruct_domain => mpp_define_unstruct_domain, &
                             fmsMppDomainUG => domainUG, &
                             fms_mpp_domains_get_UG_io_domain => mpp_get_UG_io_domain, &
                             fms_mpp_domains_get_UG_domain_npes => mpp_get_UG_domain_npes, &
                             fms_mpp_domains_get_UG_compute_domain => mpp_get_UG_compute_domain, &
                             fms_mpp_domains_get_UG_domain_tile_id => mpp_get_UG_domain_tile_id, &
                             fms_mpp_domains_get_UG_domain_pelist => mpp_get_UG_domain_pelist, &
                             fms_mpp_domains_get_ug_domain_grid_index => mpp_get_ug_domain_grid_index, &
                             fms_mpp_domains_get_UG_domain_ntiles => mpp_get_UG_domain_ntiles, &
                             fms_mpp_domains_get_UG_global_domain => mpp_get_UG_global_domain, &
                             fms_mpp_domains_global_field_ug => mpp_global_field_ug, &
                             fms_mpp_domains_get_ug_domain_tile_list => mpp_get_ug_domain_tile_list, &
                             fms_mpp_domains_get_UG_compute_domains => mpp_get_UG_compute_domains, &
                             fms_mpp_domains_define_null_UG_domain => mpp_define_null_UG_domain, &
                             fms_mpp_domains_NULL_DOMAINUG => NULL_DOMAINUG, &
                             fms_mpp_domains_get_UG_domains_index => mpp_get_UG_domains_index, &
                             fms_mpp_domains_get_UG_SG_domain => mpp_get_UG_SG_domain, &
                             fms_mpp_domains_get_UG_domain_tile_pe_inf => mpp_get_UG_domain_tile_pe_inf, &
                             fms_mpp_domains_define_nest_domains => mpp_define_nest_domains, &
                             fms_mpp_domains_get_C2F_index => mpp_get_C2F_index, &
                             fms_mpp_domains_get_F2C_index => mpp_get_F2C_index, &
                             fms_mpp_domains_get_nest_coarse_domain => mpp_get_nest_coarse_domain, &
                             fms_mpp_domains_get_nest_fine_domain => mpp_get_nest_fine_domain, &
                             fms_mpp_domains_is_nest_coarse => mpp_is_nest_coarse, &
                             fms_mpp_domains_is_nest_fine => mpp_is_nest_fine, &
                             fms_mpp_domains_get_nest_pelist => mpp_get_nest_pelist, &
                             fms_mpp_domains_get_nest_npes => mpp_get_nest_npes, &
                             fms_mpp_domains_get_nest_fine_pelist => mpp_get_nest_fine_pelist, &
                             fms_mpp_domains_get_nest_fine_npes => mpp_get_nest_fine_npes, &
                             fms_mpp_domains_domain_UG_is_tile_root_pe => mpp_domain_UG_is_tile_root_pe, &
                             fms_mpp_domains_deallocate_domainUG => mpp_deallocate_domainUG, &
                             fms_mpp_domains_get_io_domain_UG_layout => mpp_get_io_domain_UG_layout, &
                             NULL_DOMAIN1D, &
                             NULL_DOMAIN2D, &
                             fms_mpp_domains_create_super_grid_domain => mpp_create_super_grid_domain, &
                             fms_mpp_domains_shift_nest_domains => mpp_shift_nest_domains
  !> parser
#ifdef use_yaml
  use yaml_parser_mod, only: fms_yaml_parser_open_and_parse_file => open_and_parse_file, &
                             fms_yaml_parser_get_num_blocks => get_num_blocks, &
                             fms_yaml_parser_get_block_ids => get_block_ids, &
                             fms_yaml_parser_get_value_from_key => get_value_from_key, &
                             fms_yaml_parser_get_nkeys => get_nkeys, &
                             fms_yaml_parser_get_key_ids => get_key_ids, &
                             fms_yaml_parser_get_key_name => get_key_name, &
                             fms_yaml_parser_get_key_value => get_key_value
#endif

  !> platform
  use platform_mod, only: r8_kind, r4_kind, i8_kind, i4_kind, c8_kind, c4_kind, &
                          l8_kind, l4_kind, i2_kind, ptr_kind

  !> random_numbers
  use random_numbers_mod, only: fms_random_numbers_randomNumberStream => randomNumberStream, &
                                fms_random_numbers_initializeRandomNumbersStream => initializeRandomNumberStream, &
                                fms_random_numbers_getRandomNumbers => getRandomNumbers, &
                                fms_random_numbers_constructSeed => constructSeed

  !> sat_vapor_pres
  use sat_vapor_pres_mod, only: fms_sat_vapor_pres_lookup_es => lookup_es, &
                                fms_sat_vapor_pres_lookup_des => lookup_des, &
                                fms_sat_vapor_pres_init => sat_vapor_pres_init, &
                                fms_sat_vapor_pres_lookup_es2 => lookup_es2, &
                                fms_sat_vapor_pres_lookup_des2 => lookup_des2, &
                                fms_sat_vapor_pres_lookup_es2_des2 => lookup_es2_des2, &
                                fms_sat_vapor_pres_lookup_es3 => lookup_es3, &
                                fms_sat_vapor_pres_lookup_des3 => lookup_des3, &
                                fms_sat_vapor_pres_lookup_es3_des3 => lookup_es3_des3, &
                                fms_sat_vapor_pres_lookup_es_des => lookup_es_des, &
                                fms_sat_vapor_pres_compute_qs => compute_qs, &
                                fms_sat_vapor_pres_compute_mrs => compute_mrs, &
                                fms_sat_vapor_pres_escomp => escomp, &
                                fms_sat_vapor_pres_descomp => descomp
  !> string_utils
  use fms_string_utils_mod, only: fms_string_utils_string => string, &
                                  fms_string_utils_array_to_pointer => fms_array_to_pointer, &
                                  fms_string_utils_fms_pointer_to_array => fms_pointer_to_array, &
                                  fms_string_utils_sort_this => fms_sort_this, &
                                  fms_string_utils_find_my_string => fms_find_my_string, &
                                  fms_string_utils_find_unique => fms_find_unique, &
                                  fms_string_utils_c2f_string => fms_c2f_string, &
                                  fms_string_utils_cstring2cpointer => fms_cstring2cpointer, &
                                  fms_string_utils_copy => string_copy

  !> time_interp
  use time_interp_mod, only: fms_time_interp_init => time_interp_init, &
                             fms_time_interp => time_interp, fms_fraction_of_year=> fraction_of_year, &
                             NONE, YEAR, MONTH, DAY
  use time_interp_external2_mod, only: fms_time_interp_external_init_external_field => init_external_field, &
                             fms_time_interp_external => time_interp_external, &
                             fms_time_interp_external_init => time_interp_external_init, &
                             fms_time_interp_external_exit => time_interp_external_exit, &
                             fms_time_interp_external_get_external_field_size => get_external_field_size, &
                             fms_time_interp_external_get_time_axis => get_time_axis, &
                             fms_time_interp_external_get_external_field_missing => get_external_field_missing, &
                             fms_time_interp_external_set_override_region => set_override_region, &
                             fms_time_interp_external_reset_src_data_region => reset_src_data_region, &
                             fms_time_interp_external_get_external_fileobj => get_external_fileobj, &
                             NO_REGION, INSIDE_REGION, OUTSIDE_REGION, &
                             SUCCESS, ERR_FIELD_NOT_FOUND

  !> time_manager
  use time_manager_mod, only: FmsTime_type => time_type, &
                              operator(+), operator(-), operator(*), assignment(=),&
                              operator(/), operator(>), operator(>=), operator(==), &
                              operator(/=), operator(<), operator(<=), operator(//), &
                              fms_time_manager_set_time => set_time, &
                              fms_time_manager_increment_time => increment_time, &
                              fms_time_manager_decrement_time => decrement_time, &
                              fms_time_manager_get_time => get_time, &
                              fms_time_manager_interval_alarm => interval_alarm, &
                              fms_time_manager_repeat_alarm => repeat_alarm, &
                              fms_time_manager_time_type_to_real => time_type_to_real, &
                              fms_time_manager_real_to_time_type => real_to_time_type, &
                              fms_time_manager_time_list_error => time_list_error, &
                              THIRTY_DAY_MONTHS, JULIAN, GREGORIAN, NOLEAP, NO_CALENDAR, INVALID_CALENDAR, &
                              fms_time_manager_set_calendar_type => set_calendar_type, &
                              fms_time_manager_get_calendar_type => get_calendar_type, &
                              fms_time_manager_set_ticks_per_second => set_ticks_per_second, &
                              fms_time_manager_get_ticks_per_second => get_ticks_per_second, &
                              fms_time_manager_set_date => set_date, &
                              fms_time_manager_get_date => get_date, &
                              fms_time_manager_increment_date => increment_date, &
                              fms_time_manager_decrement_date => decrement_date, &
                              fms_time_manager_days_in_month => days_in_month, &
                              fms_time_manager_leap_year => leap_year, &
                              fms_time_manager_length_of_year => length_of_year, &
                              fms_time_manager_days_in_year => days_in_year, &
                              fms_time_manager_day_of_year => day_of_year, &
                              fms_time_manager_month_name => month_name, &
                              fms_time_manager_valid_calendar_types => valid_calendar_types, &
                              fms_time_manager_init => time_manager_init, &
                              fms_time_manager_print_time => print_time, &
                              fms_time_manager_print_date => print_date, &
                              fms_time_manager_set_date_julian => set_date_julian, &
                              fms_time_manager_get_date_julian => get_date_julian, &
                              fms_time_manager_get_date_no_leap => get_date_no_leap, &
                              fms_time_manager_date_to_string => date_to_string
  use get_cal_time_mod, only: fms_get_cal_time => get_cal_time

  !> topography
  use gaussian_topog_mod, only: fms_gaussian_topog_init => gaussian_topog_init, &
                                fms_get_gaussian_topog => get_gaussian_topog
  use topography_mod, only: fms_topography_init => topography_init, &
                            fms_topography_get_topog_mean => get_topog_mean, &
                            fms_topography_get_topog_stdev => get_topog_stdev, &
                            fms_topography_get_ocean_frac => get_ocean_frac, &
                            fms_topography_get_ocean_mask => get_ocean_mask, &
                            fms_topography_get_water_frac => get_water_frac, &
                            fms_topography_get_water_mask => get_water_mask

  !> tracer_manager
  use tracer_manager_mod, only: fms_tracer_manager_init => tracer_manager_init, &
                                fms_tracer_manager_end => tracer_manager_end, &
                                fms_tracer_manager_check_if_prognostic => check_if_prognostic, &
                                fms_tracer_manager_get_tracer_indices => get_tracer_indices, &
                                fms_tracer_manager_get_tracer_index => get_tracer_index, &
                                fms_tracer_manager_get_tracer_names => get_tracer_names, &
                                fms_tracer_manager_get_tracer_name => get_tracer_name, &
                                fms_tracer_manager_query_method => query_method, &
                                fms_tracer_manager_set_tracer_atts => set_tracer_atts, &
                                fms_tracer_manager_set_tracer_profile => set_tracer_profile, &
                                fms_tracer_manager_register_tracers => register_tracers, &
                                fms_tracer_manager_get_number_tracers => get_number_tracers, &
                                fms_tracer_manager_adjust_mass => adjust_mass, &
                                fms_tracer_manager_adjust_positive_def => adjust_positive_def, &
                                NO_TRACER, MAX_TRACER_FIELDS

  !> tridiagonal
  use tridiagonal_mod, only: fms_tridiagonal_tri_invert => tri_invert, &
                             fms_tridiagonal_close_tridiagonal => close_tridiagonal

  implicit none

#include <file_version.h>
  character(len=*), parameter, public :: version_FMS = version
  private :: version

end module fms
!> @}
! close documentation grouping
