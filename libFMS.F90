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
!! Remappings due to conflicts:
!!
!!           get_mosaic_tile_grid from mosaic2(fms2_io) => mosaic2_get_mosaic_tile_grid
!!
!!           read_data from interpolator_mod(fms2_io)   => interpolator_read_data
!!
!!           ZERO from interpolator_mod(mpp_parameter)  => INTERPOLATOR_ZERO
!!
!!           version from fms_mod                       => version_FMS
!!
!! Not in this module:
!!
!!                     axis_utils_mod, fms_io_mod, time_interp_external_mod
!!                     get_grid_version_mpp_mod, mpp_io_mod, mosaic_mod,
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
  use amip_interp_mod, only: amip_interp_init, get_amip_sst, get_amip_ice, &
                             amip_interp_new,amip_interp_del, amip_interp_type, &
                             assignment(=), i_sst, j_sst, sst_ncep, sst_anom, &
                             forecast_mode, use_ncep_sst
  !> astronomy
  use astronomy_mod, only: astronomy_init, get_period, set_period, &
                           set_orbital_parameters, get_orbital_parameters, &
                           set_ref_date_of_ae, get_ref_date_of_ae,  &
                           diurnal_solar, daily_mean_solar, annual_mean_solar,  &
                           astronomy_end, universal_time, orbital_time

  !> axis_utils
  use axis_utils2_mod, only: get_axis_cart, get_axis_modulo, lon_in_range, &
                             tranlon, frac_index, nearest_index, interp_1d, &
                             get_axis_modulo_times, axis_edges

  !>block_control
  use block_control_mod, only: block_control_type, define_blocks, &
                               define_blocks_packed

  !> column_diagnostics
  use column_diagnostics_mod, only: column_diagnostics_init, &
                                    initialize_diagnostic_columns, &
                                    column_diagnostics_header, &
                                    close_column_diagnostics_units

  !> coupler
  use coupler_types_mod, only: coupler_types_init, coupler_type_copy, &
                               coupler_type_spawn, coupler_type_set_diags, &
                               coupler_type_write_chksums, coupler_type_send_data, &
                               coupler_type_data_override, coupler_type_register_restarts, &
                               coupler_type_restore_state, coupler_type_increment_data, &
                               coupler_type_rescale_data, coupler_type_copy_data, &
                               coupler_type_redistribute_data, coupler_type_destructor, &
                               coupler_type_initialized, coupler_type_extract_data, &
                               coupler_type_set_data,coupler_type_copy_1d_2d, &
                               coupler_type_copy_1d_3d, coupler_3d_values_type, &
                               coupler_3d_field_type, coupler_3d_bc_type, &
                               coupler_2d_values_type, coupler_2d_field_type, &
                               coupler_2d_bc_type, coupler_1d_values_type, &
                               coupler_1d_field_type, coupler_1d_bc_type, &
                               ind_pcair, ind_u10, ind_psurf, ind_alpha, ind_csurf, &
                               ind_sc_no, ind_flux, ind_deltap, ind_kw, ind_flux0, &
                               ind_deposition, ind_runoff
  use ensemble_manager_mod, only: ensemble_manager_init, get_ensemble_id, get_ensemble_size, &
                               get_ensemble_pelist, ensemble_pelist_setup, &
                               get_ensemble_filter_pelist
  use atmos_ocean_fluxes_mod, only: atmos_ocean_fluxes_init, atmos_ocean_type_fluxes_init, &
                               aof_set_coupler_flux

  !> data_override
  use data_override_mod, only: data_override_init, data_override, &
                               data_override_unset_domains, data_override_UG

  !> diag_integral
  use diag_integral_mod, only: diag_integral_init, diag_integral_field_init, &
                               sum_diag_integral_field, diag_integral_output, &
                               diag_integral_end

  !> diag_manager
  !! includes imports from submodules made public
  use diag_manager_mod, only: diag_manager_init, send_data, send_tile_averaged_data, &
                           diag_manager_end, register_diag_field, register_static_field, &
                           diag_axis_init, get_base_time, get_base_date, need_data, &
                           DIAG_ALL, DIAG_OCEAN, DIAG_OTHER, get_date_dif, DIAG_SECONDS,&
                           DIAG_MINUTES, DIAG_HOURS, DIAG_DAYS, DIAG_MONTHS, DIAG_YEARS, &
                           get_diag_global_att, set_diag_global_att, diag_field_add_attribute, &
                           diag_field_add_cell_measures, get_diag_field_id, &
                           diag_axis_add_attribute, diag_grid_init, diag_grid_end, &
                           diag_manager_set_time_end, diag_send_complete, &
                           diag_send_complete_instant, DIAG_FIELD_NOT_FOUND, &
                           CMOR_MISSING_VALUE, null_axis_id

  !> exchange
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
                       stocks_report, stocks_report_init, stocks_set_init_time, &
                       atm_stock, ocn_stock, lnd_stock, ice_stock

  !> field manager
  use field_manager_mod, only: field_manager_init, field_manager_end, find_field_index, &
                         get_field_info, &
                         get_field_method, get_field_methods, parse, fm_change_list, &
                         fm_change_root, fm_dump_list, fm_exists, fm_get_index, &
                         fm_get_current_list, fm_get_length, fm_get_type, fm_get_value, &
                         fm_intersection, fm_init_loop, &
                         fm_loop_over_list, fm_new_list, fm_new_value, &
                         fm_reset_loop, fm_return_root, &
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

  !> fms2_io
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
                         fms2_io_init, get_mosaic_tile_grid, &
                         write_restart_bc, read_restart_bc, get_filename_appendix, & !> 2021.02-a1
                         set_filename_appendix, get_instance_filename, &
                         nullify_filename_appendix, ascii_read, get_mosaic_tile_file, &
                         parse_mask_table
  ! used via fms2_io
  ! fms_io_utils_mod, fms_netcdf_domain_io_mod, netcdf_io_mod,
  ! fms_netcdf_unstructured_domain_io_mod, blackboxio

  !> fms
  !! routines that don't conflict with fms2_io
  use fms_mod, only: fms_init, fms_end, error_mesg, fms_error_handler, check_nml_error, &
                     monotonic_array, string_array_index, clock_flag_default, &
                     print_memory_usage, write_version_number

  !> horiz_interp
  use horiz_interp_mod, only: horiz_interp, horiz_interp_new, horiz_interp_del, &
                              horiz_interp_init, horiz_interp_end
  use horiz_interp_type_mod, only: horiz_interp_type, assignment(=), CONSERVE, &
                              BILINEAR, SPHERICA, BICUBIC, stats
  !! used via horiz_interp
  ! horiz_interp_bicubic_mod, horiz_interp_bilinear_mod
  ! horiz_interp_conserve_mod, horiz_interp_spherical_mod

  !> interpolator
  use interpolator_mod, only: interpolator_init, interpolator, interpolate_type_eq, &
                              obtain_interpolator_time_slices, unset_interpolator_time_flag, &
                              interpolator_end, init_clim_diag, query_interpolator, &
                              interpolate_type, CONSTANT, &
                              INTERP_WEIGHTED_P, INTERP_LINEAR_P, INTERP_LOG_P, &
                              INTERPOLATOR_ZERO=>ZERO, & !! conflicts with mpp_parameter's ZERO
                              interpolator_read_data=>read_data !! conflicts with fms2_io interface

  !> memutils
  use memutils_mod, only: memutils_init, print_memuse_stats

  !> monin_obukhov
  use monin_obukhov_mod, only: monin_obukhov_init, monin_obukhov_end, &
                               mo_drag, mo_profile, mo_diff, stable_mix
  use monin_obukhov_inter, only: monin_obukhov_diff, monin_obukhov_drag_1d, &
                               monin_obukhov_solve_zeta, monin_obukhov_derivative_t, &
                               monin_obukhov_derivative_m, monin_obukhov_profile_1d, &
                               monin_obukhov_integral_m, monin_obukhov_integral_tq, &
                               monin_obukhov_stable_mix

  !> mosaic
  use mosaic2_mod, only: get_mosaic_ntiles, get_mosaic_ncontacts, &
                      get_mosaic_grid_sizes, get_mosaic_contact, &
                      get_mosaic_xgrid_size, get_mosaic_xgrid, &
                      calc_mosaic_grid_area, calc_mosaic_grid_great_circle_area, &
                      is_inside_polygon, &
                      mosaic2_get_mosaic_tile_grid => get_mosaic_tile_grid !overloaded in fms2_io
  use grid2_mod, only: get_grid_ntiles, get_grid_size, get_grid_cell_centers, &
                      get_grid_cell_vertices, get_grid_cell_Area, get_grid_comp_area, &
                      define_cube_mosaic, get_great_circle_algorithm, grid_init, grid_end
  use gradient_mod, only: gradient_cubic, calc_cubic_grid_info

  !> mpp
  use mpp_mod, only: stdin, stdout, stderr, &
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
                     mpp_type_free, input_nml_file
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
                     mpp_reset_efp_overflow_error, mpp_efp_type
  use mpp_domains_mod, only: domain_axis_spec, domain1D, domain2D, DomainCommunicator2D, &
                     nest_domain_type, mpp_group_update_type, &
                     mpp_domains_set_stack_size, mpp_get_compute_domain, &
                     mpp_get_compute_domains, mpp_get_data_domain, &
                     mpp_get_global_domain, mpp_get_domain_components, &
                     mpp_get_layout, mpp_get_pelist, operator(.EQ.), operator(.NE.), &
                     mpp_domain_is_symmetry, mpp_domain_is_initialized, &
                     mpp_get_neighbor_pe, mpp_nullify_domain_list, &
                     mpp_set_compute_domain, mpp_set_data_domain, mpp_set_global_domain, &
                     mpp_get_memory_domain, mpp_get_domain_shift, &
                     mpp_domain_is_tile_root_pe, mpp_get_tile_id, &
                     mpp_get_domain_extents, mpp_get_current_ntile, &
                     mpp_get_ntile_count, mpp_get_tile_list, mpp_get_tile_npes, &
                     mpp_get_domain_root_pe, mpp_get_tile_pelist, &
                     mpp_get_tile_compute_domains, mpp_get_num_overlap, &
                     mpp_get_overlap, mpp_get_io_domain, mpp_get_domain_pe, &
                     mpp_get_domain_tile_root_pe, mpp_get_domain_name, &
                     mpp_get_io_domain_layout, mpp_copy_domain, mpp_set_domain_symmetry, &
                     mpp_get_update_pelist, mpp_get_update_size, &
                     mpp_get_domain_npes, mpp_get_domain_pelist, &
                     mpp_clear_group_update, mpp_group_update_initialized, &
                     mpp_group_update_is_set, mpp_get_global_domains, &
                     mpp_global_field, mpp_global_max, mpp_global_min, mpp_global_sum, &
                     mpp_global_sum_tl, mpp_global_sum_ad, mpp_broadcast_domain, &
                     mpp_domains_init, mpp_domains_exit, mpp_redistribute, &
                     mpp_update_domains, mpp_check_field, mpp_start_update_domains, &
                     mpp_complete_update_domains, mpp_create_group_update, &
                     mpp_do_group_update, mpp_start_group_update, &
                     mpp_complete_group_update, mpp_reset_group_update_field, &
                     mpp_update_nest_fine, mpp_update_nest_coarse, mpp_get_boundary, &
                     mpp_update_domains_ad, mpp_get_boundary_ad, mpp_pass_SG_to_UG, &
                     mpp_pass_UG_to_SG, mpp_define_layout, mpp_define_domains, &
                     mpp_modify_domain, mpp_define_mosaic, mpp_define_mosaic_pelist, &
                     mpp_define_null_domain, mpp_mosaic_defined, &
                     mpp_define_io_domain, mpp_deallocate_domain, &
                     mpp_compute_extent, mpp_compute_block_extent, &
                     mpp_define_unstruct_domain, domainUG, mpp_get_UG_io_domain, &
                     mpp_get_UG_domain_npes, mpp_get_UG_compute_domain, &
                     mpp_get_UG_domain_tile_id, mpp_get_UG_domain_pelist, &
                     mpp_get_ug_domain_grid_index, mpp_get_UG_domain_ntiles, &
                     mpp_get_UG_global_domain, mpp_global_field_ug, &
                     mpp_get_ug_domain_tile_list, mpp_get_UG_compute_domains, &
                     mpp_define_null_UG_domain, NULL_DOMAINUG, mpp_get_UG_domains_index, &
                     mpp_get_UG_SG_domain, mpp_get_UG_domain_tile_pe_inf, &
                     mpp_define_nest_domains, mpp_get_C2F_index, mpp_get_F2C_index, &
                     mpp_get_nest_coarse_domain, mpp_get_nest_fine_domain, &
                     mpp_is_nest_coarse, mpp_is_nest_fine, &
                     mpp_get_nest_pelist, mpp_get_nest_npes, &
                     mpp_get_nest_fine_pelist, mpp_get_nest_fine_npes, &
                     mpp_domain_UG_is_tile_root_pe, mpp_deallocate_domainUG, &
                     mpp_get_io_domain_UG_layout, NULL_DOMAIN1D, NULL_DOMAIN2D, &
                     mpp_create_super_grid_domain

  !> platform
  use platform_mod, only: r8_kind, r4_kind, i8_kind, i4_kind, c8_kind, c4_kind, &
                          l8_kind, l4_kind, i2_kind, ptr_kind

  !> random_numbers
  use random_numbers_mod, only: randomNumberStream, initializeRandomNumberStream, &
                                getRandomNumbers, constructSeed

  !> sat_vapor_pres
  use sat_vapor_pres_mod, only: lookup_es, lookup_des, sat_vapor_pres_init, &
                                lookup_es2, lookup_des2, lookup_es2_des2, &
                                lookup_es3, lookup_des3, lookup_es3_des3, &
                                lookup_es_des, compute_qs, compute_mrs, &
                                escomp, descomp

  !> time_interp
  use time_interp_mod, only: time_interp_init, time_interp, fraction_of_year, &
                             NONE, YEAR, MONTH, DAY
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

  !> topography
  use gaussian_topog_mod, only: gaussian_topog_init, get_gaussian_topog
  use topography_mod, only: topography_init, get_topog_mean, get_topog_stdev, &
                            get_ocean_frac, get_ocean_mask, get_water_frac, &
                            get_water_mask

  !> tracer_manager
  use tracer_manager_mod, only: tracer_manager_init, tracer_manager_end, &
                                check_if_prognostic, get_tracer_indices,  &
                                get_tracer_index, get_tracer_names, &
                                get_tracer_name, query_method, &
                                set_tracer_atts, set_tracer_profile, &
                                register_tracers, get_number_tracers,  &
                                adjust_mass, adjust_positive_def, NO_TRACER, MAX_TRACER_FIELDS

  !> tridiagonal
  use tridiagonal_mod, only: tri_invert, close_tridiagonal

  implicit none

#include <file_version.h>
  character(len=*), parameter, public :: version_FMS = version
  private :: version

end module fms
!> @}
! close documentation grouping
