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
!! @brief A convenience file to use all public FMS routines, functions, values
!! @description Imports all FMS modules so that any public interfaces,
!! variables and routines are usable via this module
!! Does not include routines from/for mpp_io

module fms

  !> import each FMS module's public routines/functions, interfaces, and variables
  !! done explicitly to avoid including any unwanted/depracated routines
  use fms_affinity_mod, only: fms_affinity_init, fms_affinity_get, &
                              fms_affinity_set
  
  use amip_interp_mod, only: amip_interp_init, get_amip_sst, get_amip_ice, &
                             amip_interp_new,amip_interp_del, amip_interp_type, &
                             assignment(=), i_sst, j_sst, sst_ncep, sst_anom, &
                             forecast_mode, use_ncep_sst

  use astronomy_mod, only: astronomy_init, get_period, set_period, &
                           set_orbital_parameters, get_orbital_parameters, &
                           set_ref_date_of_ae, get_ref_date_of_ae,  &
                           diurnal_solar, daily_mean_solar, annual_mean_solar,  &
                           astronomy_end, universal_time, orbital_time

  !> axis_utils for fms2_io
  use axis_utils2_mod, only: get_axis_cart, get_axis_modulo, lon_in_range, &
                             tranlon, frac_index, nearest_index, interp_1d, &
                             get_axis_modulo_times, axis_edges
!! use axis_utils_mod, only: get_axis_bounds !! this may be needed, not in 2 

  use block_control_mod, only: block_control_type, define_blocks, &
                               define_blocks_packed

  use column_diagnostics_mod, only: column_diagnostics_init, &
                                    initialize_diagnostic_columns, &
                                    column_diagnostics_header, &
                                    close_column_diagnostics_units

  !> no other modules used so only list not neccessary
  use constants_mod
   
  !> coupler modules 
  use coupler_types_mod, only: coupler_types_init, coupler_type_copy, &
                               coupler_type_spawn, coupler_type_set_diags, &
                               coupler_type_write_chksums, coupler_type_send_data, &
                               coupler_type_data_override, coupler_type_register_restarts, &
                               coupler_type_restore_state, coupler_type_increment_data, &
                               coupler_type_rescale_data, coupler_type_copy_data, &
                               coupler_type_redistribute_data, coupler_type_destructor, &
                               coupler_type_initialized, coupler_type_extract_data, &
                               coupler_type_set_data,coupler_type_copy_1d_2d, &
                               coupler_type_copy_1d_3d
  use ensemble_manager_mod, only: ensemble_manager_init, get_ensemble_id, get_ensemble_size, &
                                  get_ensemble_pelist, ensemble_pelist_setup, &
                                  get_ensemble_filter_pelist
  use atmos_ocean_fluxes_mod, only: atmos_ocean_fluxes_init, atmos_ocean_type_fluxes_init, &
                                    aof_set_coupler_flux 

  !> data override modules
  use data_override_mod, only: data_override_init, data_override, &
                               data_override_unset_domains, data_override_UG 
  use get_grid_version_fms2io_mod, only: check_grid_sizes, deg_to_radian,&
                                         get_grid_version_1, get_grid_version_2

  use diag_integral_mod, only: diag_integral_init, diag_integral_field_init, &
                               sum_diag_integral_field, diag_integral_output, &
                               diag_integral_end

  !> diag manager modules
  use diag_manager_mod, only: diag_manager_init, send_data,send_tile_averaged_data, &
                              diag_manager_end, register_diag_field, register_static_field, &
                              diag_axis_init, get_base_time, get_base_date, need_data, &
                              DIAG_ALL, DIAG_OCEAN, DIAG_OTHER, get_date_dif, &
                              DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, &
                              DIAG_YEARS, get_diag_global_att, set_diag_global_att, &
                              diag_field_add_attribute, diag_field_add_cell_measures, &
                              get_diag_field_id, &
                              !> from diag_grid_mod
                              diag_grid_init, diag_grid_end, diag_manager_set_time_end, & 
                              diag_send_complete, diag_send_complete_instant
  use diag_table_mod, only: parse_diag_table
  use diag_util_mod, only: get_subfield_size, log_diag_field_info, update_bounds, &
                           check_out_of_bounds, check_bounds_are_exact_dynamic, &
                           check_bounds_are_exact_static, init_file, diag_time_inc, &
                           find_input_field, init_input_field, init_output_field, &
                           diag_data_out, write_static, check_duplicate_output_fields, &
                           get_date_dif, get_subfield_vert_size, sync_file_times, &
                           prepend_attribute, attribute_init, diag_util_init
  use diag_output_mod, only: diag_output_init, write_axis_meta_data, write_field_meta_data, &
                             done_meta_data, diag_fieldtype, get_diag_global_att, &
                             set_diag_global_att, diag_field_write, diag_write_time
                             !!diag_field_out, done_meta_data_use_mpp_io !< mpp_io
  use diag_axis_mod, only: diag_axis_init, get_diag_axis, get_domain1d, &
                           get_domain2d, get_axis_length, get_axis_global_length, &
                           diag_subaxes_init, get_diag_axis_cart, get_diag_axis_data, &
                           max_axes, get_axis_aux, get_tile_count, get_axes_shift, &
                           get_diag_axis_name, get_axis_num, get_diag_axis_domain_name, &
                           diag_axis_add_attribute, get_domainUG, axis_compatible_check, &
                           axis_is_compressed, get_compressed_axes_ids, get_axis_reqfld, &
                           NORTH, EAST, CENTER
  use diag_data_mod, only: MAX_FIELDS_PER_FILE, DIAG_OTHER, DIAG_OCEAN, DIAG_ALL, &
                           VERY_LARGE_FILE_FREQ, VERY_LARGE_AXIS_LENGTH, EVERY_TIME, &
                           END_OF_RUN, DIAG_SECONDS, DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, &
                           DIAG_MONTHS, DIAG_YEARS, MAX_SUBAXES, GLO_REG_VAL, &
                           GLO_REG_VAL_ALT, CMOR_MISSING_VALUE, DIAG_FIELD_NOT_FOUND, &
                           num_files, num_input_fields, num_output_fields, null_axis_id, &
                           diag_grid, diag_fieldtype, diag_atttype, coord_type, file_type, &
                           input_field_type, output_field_type, diag_axis_type, &
                           diag_global_att_type

  !> drifters modules
  use drifters_mod, only: drifters_type, assignment(=), drifters_push, &
                          drifters_compute_k, drifters_set_field, drifters_new, &
                          drifters_del, drifters_set_domain, drifters_set_pe_neighbors, &
                          drifters_set_v_axes, drifters_set_domain_bounds, & 
                          drifters_positions2lonlat, drifters_print_checksums, &
                          drifters_save, drifters_write_restart, drifters_distribute
  use drifters_comm_mod, only: drifters_comm_Type, drifters_comm_new, &
                               drifters_comm_del, drifters_comm_set_pe_neighbors, &
                               drifters_comm_set_domain, drifters_comm_update, &
                               drifters_comm_gather
  use drifters_core_mod, only: drifters_core_type, drifters_core_new, drifters_core_del, &
                               drifters_core_set_ids, drifters_core_remove_and_add, &
                               drifters_core_set_positions, assignment(=), &
                               drifters_core_print,  drifters_core_resize 
  use drifters_input_mod, only: drifters_input_type, drifters_input_new, &
                                drifters_input_del, drifters_input_save, &
                                assignment(=)
  use drifters_io_mod, only: drifters_io_type, drifters_io_new, drifters_io_del, &
                             drifters_io_set_time_units, drifters_io_set_position_names, &
                             drifters_io_set_position_units, drifters_io_set_field_names, &
                             drifters_io_set_field_units, drifters_io_write
  use cloud_interpolator_mod, only: cld_ntrp_linear_cell_interp, cld_ntrp_locate_cell, &
                                    cld_ntrp_get_cell_values, cld_ntrp_expand_index, &
                                    cld_ntrp_contract_indices

  !> exchange modules
  !! TODO test
  use xgrid_mod, only: xmap_type, setup_xmap, set_frac_area, put_to_xgrid, &
                       get_from_xgrid, xgrid_count, some, conservation_check, &
                       xgrid_init, AREA_ATM_SPHERE, AREA_OCN_SPHERE, AREA_ATM_MODEL, &
                       AREA_OCN_MODEL, get_ocean_model_area_elements, grid_box_type, &
                       get_xmap_grid_area, put_to_xgrid_ug, get_from_xgrid_ug, &
                       set_frac_area_ug, FIRST_ORDER, SECOND_ORDER, stock_move_ug, & 
                       stock_move, stock_type, stock_print, get_index_range, &
                       stock_integrate_2d
  use stock_constants_mod, only: NELEMS, ISTOCK_WATER, ISTOCK_HEAT, ISTOCK_SALT, &
                                 ISTOCK_TOP, ISTOCK_BOTTOM, ISTOCK_SIDE, stocks_file, &
                                 stocks_report, stocks_report_init, stocks_set_init_time

  !> field manager modules
  use field_manager_mod, only: field_manager_init, field_manager_end, find_field_index, &
                               find_field_index_old, find_field_index_new, get_field_info, &
                               get_field_method, get_field_methods, parse, fm_change_list, &
                               fm_change_root, fm_dump_list, fm_exists, fm_get_index, &
                               fm_get_current_list, fm_get_length, fm_get_type, fm_get_value, &
                               fm_get_value_integer, fm_get_value_logical, fm_get_value_real, &
                               fm_get_value_string, fm_intersection, fm_init_loop, &
                               fm_loop_over_list, fm_new_list, fm_new_value, &
                               fm_new_value_integer, fm_new_value_logical,fm_new_value_real, &
                               fm_new_value_string, fm_reset_loop, fm_return_root, &
                               fm_modify_name, fm_query_method, fm_find_methods, fm_copy_list, &
                               fm_set_verbosity, fm_field_name_len, fm_path_name_len, &
                               fm_string_len, fm_type_name_len, NUM_MODELS, NO_FIELD, &
                               MODEL_ATMOS, MODEL_OCEAN, MODEL_LAND, MODEL_ICE, MODEL_COUPLER, &
                               fm_array_list_def, method_type, method_type_short, &
                               method_type_very_short, fm_list_iter_type, default_method
  use fm_util_mod, only: fm_util_start_namelist, fm_util_end_namelist, &
                         fm_util_check_for_bad_fields, fm_util_set_caller, &
                         fm_util_reset_caller, fm_util_set_no_overwrite, &
                         fm_util_reset_no_overwrite, fm_util_set_good_name_list, &
                         fm_util_reset_good_name_list, fm_util_get_length, &
                         fm_util_get_integer, fm_util_get_logical, fm_util_get_real, &
                         fm_util_get_string, fm_util_get_integer_array, &
                         fm_util_get_logical_array, fm_util_get_real_array, &
                         fm_util_get_string_array, fm_util_set_value, &
                         fm_util_set_value_integer_array, fm_util_set_value_logical_array, &
                         fm_util_set_value_real_array, fm_util_set_value_string_array, &
                         fm_util_set_value_integer, fm_util_set_value_logical, &
                         fm_util_set_value_real, fm_util_set_value_string, &
                         fm_util_get_index_list, fm_util_get_index_string, &
                         fm_util_default_caller

  !>fms2_io modules TODO testing
  use fms2_io_mod, only: unlimited, FmsNetcdfFile_t, FmsNetcdfDomainFile_t, &
                         FmsNetcdfUnstructuredDomainFile_t, open_file, open_virtual_file, &
                         close_file, register_axis, register_field, register_restart_field, &
                         write_data, read_data, write_restart, write_new_restart, &
                         read_restart, read_new_restart, global_att_exists, &
                         variable_att_exists, register_global_attribute, &
                         register_variable_attribute, get_global_attribute, &
                         get_variable_attribute, get_num_dimensions, &
                         get_dimension_names, dimension_exists, is_dimension_unlimited, &
                         get_dimension_size, get_num_variables, get_variable_names, &
                         variable_exists, get_variable_num_dimensions, &
                         get_variable_dimension_names, get_variable_size, &
                         get_compute_domain_dimension_indices, &
                         get_global_io_domain_indices, Valid_t, get_valid, is_valid, &
                         get_unlimited_dimension_name, get_variable_unlimited_dimension_index, &
                         file_exists, compressed_start_and_count, get_variable_sense, &
                         get_variable_missing, get_variable_units, get_time_calendar, &
                         open_check, is_registered_to_restart, check_if_open, &
                         set_fileobj_time_name, is_dimension_registered, &
                         fms2_io_init !!, get_mosaic_tile_grid TODO mosaic2 which/one
!! TODO fms2_io has interfaces for most(or all) of these
!  use fms_io_utils_mod, only: char_linked_list, error, file_exists, openmp_thread_trap, &
!                              string_copy, is_in_list, append_to_list, destroy_list, &
!                              domain_tile_filepath_mangle, io_domain_tile_filepath_mangle, &
!                              allocate_array, put_array_section, get_array_section, &
!                              get_data_type_string, get_checksum, open_check, &
!                              string_compare, restart_filepath_mangle
!  use fms_netcdf_domain_io_mod, only: FmsNetcdfDomainFile_t, open_domain_file, &
!                         close_domain_file, register_domain_decomposed_dimension, &
!                         register_domain_variable, register_domain_restart_variable_0d, &
!                         register_domain_restart_variable_1d, &
!                         register_domain_restart_variable_2d, &
!                         register_domain_restart_variable_3d, &
!                         register_domain_restart_variable_4d, &
!                         register_domain_restart_variable_5d, &
!                         domain_read_0d, domain_read_1d, domain_read_2d, domain_read_3d, &
!                         domain_read_4d, domain_read_5d, domain_write_0d, domain_write_1d, &
!                         domain_write_2d, domain_write_3d, domain_write_4d, domain_write_5d, &
!                         save_domain_restart, restore_domain_state, &
!                         get_compute_domain_dimension_indices, get_global_io_domain_indices, &
!                         is_dimension_registered, get_mosaic_tile_grid
  use netcdf_io_mod, only: no_unlimited_dimension, define_mode, data_mode, &
                           max_num_restart_vars, unlimited, max_num_compressed_dims, &
                           FmsNetcdfFile_t, Valid_t, netcdf_io_init, netcdf_file_open, &
                           netcdf_file_close, netcdf_add_dimension, netcdf_add_variable, &
                           netcdf_add_restart_variable, global_att_exists, variable_att_exists, &
                           register_global_attribute, register_variable_attribute, &
                           get_global_attribute, get_variable_attribute, get_num_dimensions, &
                           get_dimension_names, dimension_exists, is_dimension_unlimited, &
                           get_dimension_size, get_num_variables, get_variable_names, &
                           variable_exists, get_variable_num_dimensions, &
                           get_variable_dimension_names, get_variable_size, &
                           get_variable_unlimited_dimension_index, netcdf_read_data, &
                           netcdf_write_data, compressed_write, netcdf_save_restart, &
                           netcdf_restore_state, get_valid, is_valid, &
                           get_unlimited_dimension_name, netcdf_file_open_wrap, &
                           netcdf_file_close_wrap, netcdf_add_variable_wrap, &
                           netcdf_save_restart_wrap, compressed_write_0d_wrap, &
                           compressed_write_1d_wrap, compressed_write_2d_wrap, &
                           compressed_write_3d_wrap, compressed_write_4d_wrap, &
                           compressed_write_5d_wrap, compressed_read_0d, compressed_read_1d, &
                           compressed_read_2d, compressed_read_3d, compressed_read_4d, &
                           compressed_read_5d, register_compressed_dimension, &
                           netcdf_add_restart_variable_0d_wrap, &
                           netcdf_add_restart_variable_1d_wrap, &
                           netcdf_add_restart_variable_2d_wrap, &
                           netcdf_add_restart_variable_3d_wrap, &
                           netcdf_add_restart_variable_4d_wrap, &
                           netcdf_add_restart_variable_5d_wrap, &
                           compressed_start_and_count, get_fill_value, get_variable_sense, &
                           get_variable_missing, get_variable_units, get_time_calendar, &
                           is_registered_to_restart, set_netcdf_mode, check_netcdf_code, &
                           check_if_open, set_fileobj_time_name
  !use fms_netcdf_unstructured_domain_io_mod
  !use blackboxio

  !> fms modules (most unused from old io)
  use fms_mod, only: fms_init, fms_end
                     !open_namelist_file, open_restart_file, &
                     !open_ieee32_file, close_file, open_file, open_direct_file, &
                     !set_domain, read_data

  !> horiz interp modules
  use horiz_interp_mod, only: horiz_interp, horiz_interp_new, horiz_interp_del, &
                              horiz_interp_init, horiz_interp_end
  use horiz_interp_type_mod, only: horiz_interp_type, assignment(=)
!  following used for interfaces in ^
!  use horiz_interp_bicubic_mod
!  use horiz_interp_bilinear_mod
!  use horiz_interp_conserve_mod
!  use horiz_interp_spherical_mod

  use interpolator_mod, only: interpolator_init, interpolator, interpolate_type_eq, &
                              obtain_interpolator_time_slices, unset_interpolator_time_flag, &
                              interpolator_end, init_clim_diag, query_interpolator, &
                              read_data !! this one may cause conflicts

  use memutils_mod, only: memutils_init, print_memuse_stats

  !> monin obukhov modules
  use monin_obukhov_mod, only: monin_obukhov_init, monin_obukhov_end, &
                               mo_drag, mo_profile, mo_diff, stable_mix
  use monin_obukhov_inter, only: monin_obukhov_diff, monin_obukhov_drag_1d, &
                           monin_obukhov_solve_zeta, monin_obukhov_derivative_t, &
                           monin_obukhov_derivative_m, monin_obukhov_profile_1d, &
                           monin_obukhov_integral_m, monin_obukhov_integral_tq, &
                           monin_obukhov_stable_mix

  !> mosaic modules
  use mosaic2_mod, only: get_mosaic_ntiles, get_mosaic_ncontacts, &
                         get_mosaic_grid_sizes, get_mosaic_contact, &
                         get_mosaic_xgrid_size, get_mosaic_xgrid, & !! get_mosaic_tile_grid, 
                         calc_mosaic_grid_area, calc_mosaic_grid_great_circle_area, &
                         is_inside_polygon
  use grid_mod, only: get_grid_ntiles, get_grid_size, get_grid_cell_centers, &
                      get_grid_cell_vertices, get_grid_cell_Area, get_grid_comp_area, &
                      define_cube_mosaic
  use gradient_mod, only: gradient_cubic, calc_cubic_grid_info

  !> mpp modules
  use mpp_mod, only: mpp_init_test_full_init, mpp_init_test_init_true_only, &
                     mpp_init_test_peset_allocated, mpp_init_test_clocks_init, &
                     mpp_init_test_datatype_list_init, mpp_init_test_logfile_init, &
                     mpp_init_test_read_namelist, mpp_init_test_etc_unit, &
                     mpp_init_test_requests_allocated, stdin, stdout, stderr, &
                     stdlog, lowercase, uppercase, mpp_error, mpp_error_state, &
                     mpp_set_warn_level, mpp_sync, mpp_sync_self, mpp_set_stack_size, &
                     mpp_pe, mpp_npes, mpp_root_pe, mpp_set_root_pe, mpp_declare_pelist, &
                     mpp_get_current_pelist, mpp_set_current_pelist, &
                     mpp_get_current_pelist_name, mpp_clock_id, mpp_clock_set_grain, &
                     mpp_record_timing_data, get_unit, read_ascii_file, read_input_nml, &
                     mpp_clock_begin, mpp_clock_end, get_ascii_file_num_lines, &
                     mpp_record_time_start, mpp_record_time_end, mpp_chksum, &
                     mpp_max, mpp_min, mpp_sum, mpp_transmit, mpp_send, mpp_recv, &
                     mpp_sum_ad, mpp_broadcast, mpp_init, mpp_exit, mpp_gather, &
                     mpp_scatter, mpp_alltoall, mpp_type, mpp_byte, mpp_type_create, &
                     mpp_type_free
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
                             MPP_FILL_INT, MPP_FILL_FLOAT, MPP_FILL_DOUBLE, & !! mpp_domains
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
  use mpp_data_mod, only: stat, mpp_stack, ptr_stack, status, ptr_status, sync, &
                          ptr_sync, mpp_from_pe, ptr_from, remote_Data_loc, &
                          ptr_remote, mpp_domains_stack, ptr_domains_stack, &
                          mpp_domains_stack_nonblock, ptr_domains_stack_nonblock
  use mpp_utilities_mod, only: mpp_array_global_min_max
  use mpp_memutils_mod, only: mpp_print_memuse_stats, mpp_mem_dump, &
                              mpp_memuse_begin, mpp_memuse_end
  use mpp_efp_mod, only: mpp_reproducing_sum, mpp_efp_list_sum_across_PEs, &
                         mpp_efp_plus, mpp_efp_minus, mpp_efp_to_real, &
                         mpp_real_to_efp, mpp_efp_real_diff, operator(+), &
                         operator(-), assignment(=), mpp_query_efp_overflow_error, &
                         mpp_reset_efp_overlow_error
  use mpp_domains_mod, only: NULL_DOMAIN1D, NULL_DOMAIN2D 

  use platform_mod, only: r8_kind, r4_kind, i8_kind, i4_kind, c8_kind, c4_kind, &
                          l8_kind, l4_kind, i2_kind, ptr_kind

  !> random_numbers modules
  use mersennetwister_mod, only: randomNumberSequence, new_RandomNumberSequence, &
                                 finalize_RandomNumberSequence, getRandomInt, &
                                 getRandomPositiveInt, getRandomReal
  use random_numbers_mod, only: randomNumberStream, initializeRandomNumberStream, &
                                getRandomNumbers, constructSeed

  !> sat_vapor_pres
  use sat_vapor_pres_k_mod, only: sat_vapor_pres_init_k, lookup_es_k, lookup_des_k, &
                                  lookup_es_des_k, lookup_es2_k, lookup_des2_k, &
                                  lookup_es2_des2_k, lookup_es3_k, lookup_des3_k, &
                                  lookup_es3_des3_k, compute_qs_k, compute_mrs_k
  use sat_vapor_pres_mod, only: lookup_es, lookup_des, sat_vapor_pres_init, &
                                lookup_es2, lookup_des2, lookup_es2_des2, &
                                lookup_es3, lookup_des3, lookup_es3_des3, &
                                lookup_es_des, compute_qs, compute_mrs, &
                                escomp, descomp 
  !> time_interp 
  use time_interp_mod, only: time_interp_init, time_interp, fraction_of_year
  use time_interp_external2_mod, only: init_external_field, time_interp_external, &
                                       time_interp_external_init, time_interp_external_exit, &
                                       get_external_field_size, get_time_axis, &
                                       get_external_field_missing, set_override_region, &
                                       reset_src_data_region, get_external_fileobj, &
                                       NO_REGION, INSIDE_REGION, OUTSIDE_REGION, &
                                       SUCCESS, ERR_FIELD_NOT_FOUND
 
  !> time_manager
  use time_manager_mod, only: time_type, operator(+), operator(-), operator(*), &
                              operator(/), operator(>), operator(>=), operator(==), &
                              operator(/=), operator(<), operator(<=), operator(//), &
                              assignment(=), set_time, increment_time, decrement_time, &
                              get_time, interval_alarm, repeat_alarm, time_type_to_real, &
                              real_to_time_type, time_list_error, THIRTY_DAY_MONTHS, &
                              JULIAN, GREGORIAN, NOLEAP, NO_CALENDAR, INVALID_CALENDAR, &
                              set_calendar_type, get_calendar_type, set_ticks_per_second, &
                              get_ticks_per_second, set_date, get_date, increment_date, &
                              decrement_date, days_in_month, leap_year, length_of_year, &
                              days_in_year, day_of_year, month_name, valid_calendar_types, &
                              time_manager_init, print_time, print_date, set_date_julian, &
                              get_date_julian, get_date_no_leap, date_to_string
  use get_cal_time_mod, only: get_cal_time
  
  !> topography modules
  use gaussian_topog_mod, only: gaussian_topog_init, get_gaussian_topog
  use topography_mod, only: topography_init, get_topog_mean, get_topog_stdev, &
                            get_ocean_frac, get_ocean_mask, get_water_frac, &
                            get_water_mask

  use tracer_manager_mod, only: tracer_manager_init, tracer_manager_end, &
                                check_if_prognostic, get_tracer_indices,  &
                                get_tracer_index, get_tracer_names, &
                                get_tracer_name, query_method, &
                                set_tracer_atts, set_tracer_profile, &
                                register_tracers, get_number_tracers,  &
                                adjust_mass, adjust_positive_def, NO_TRACER, MAX_TRACER_FIELDS

  use tridiagonal_mod, only: tri_invert, close_tridiagonal
  !! these are the depracated overloaded routines from mpp_io
  !! if using any of the following, they must be imported directly
  !use axis_utils_mod, only: interp_1_fms_io => interp_1d, &
  !                          get_axis_cart_fms_io => get_axis_cart, &
  !                          get_axis_modulo_fms_io => get_axis_modulo, & 
  !                          lon_in_range_fms_io => lon_in_range, &
  !                          tranlon_fms_io => tranlon, &
  !                          frac_index_fms_io => frac_index, &
  !                          nearest_index_fms_io => nearest_index, &
  !                          get_axis_bounds
  !use fms_io_mod, only: restart_file_type, save_restart, restore_state, &
  !                      fms_io_init, fms_io_exit, open_file_fms_io=> open_file, &
  !                      register_restart_field_fms_io=> register_restart_field, &
  !                      read_data_fms_io=> read_data,write_data_fms_io=> write_data, &
  !                      close_file_fms_io=> close_file, &
  !                      register_restart_field_fms_io=> register_restart_field, &
  !                      get_mosaic_tile_grid_fms_io=> get_mosaic_tile_grid
  !use time_interp_external_mod, only:init_external_field_fms_io=> init_external_field, &
  !                                   time_interp_external_fms_io=> time_interp_external, &
  !                                   time_interp_external_init_fms_io=> time_interp_external_init, &
  !                                   time_interp_external_exit_fms_io=> time_interp_external_exit, &
  !                                   get_external_field_size_fms_io=> get_external_field_size, &
  !                                   get_time_axis_fms_io=> get_time_axis, &
  !                                   get_external_field_missing_fms_io=> get_external_field_missing, &
  !                                   set_override_region_fms_io=> set_override_region, &
  !                                   reset_src_data_region_fms_io=> reset_src_data_region, &
  !                                   get_external_field_axes
  !use get_grid_version_mpp_mod
  !use mpp_io_mod
  !use mosaic_mod
  implicit none

end module fms
