# ChangeLog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0),
and this project uses `yyyy.rr[.pp]`, where `yyyy` is the year a patch is released,
`rr` is a sequential release number (starting from `01`), and an optional two-digit
sequential patch number (starting from `01`).

## [2024.03] - 2024-08-22

### Known Issues
- Diag Manager Rewrite: See [below](#20240102---2024-06-14) for known output file differences regarding the new diag manager. The new diag_manager is disabled by default, so this differences will only be present if `use_modern_diag` is set to true in the `diag_manager_nml`.
- BUILD(HDF5): HDF5 version 1.14.3 generates floating point exceptions, and will cause errors if FMS is built with FPE traps enabled. FPE traps are turned on when using the debug target in mkmf.
- GCC: version 14.1.0 is unsupported due to a bug with strings that has come up previously in earlier versions. This will be caught by the configure script, but will cause compilation errors if using other build systems.
- INTEL: The `-check uninit` flag for the Intel Oneapi Fortran compiler (ifx) is unsupported due to a bug causing false positives when using external libraries. If using the `-check all` flag, `-check all,nouninit` should be used instead.

### Added
- DATA_OVERRIDE: Adds a namelist flag `use_center_grid_points` which, if true, enables reading the centroid values from the grid spec for ocean and ice models. This fixes issues with grid files that have longitudes ranges from 0:180 and then -180:0, but will cause answer changes if enabled. (#1566)
- DATA_OVERRIDE: Adds support for reading external weight files. Currently only supported for fregrid generated files while using the bilinear interpolation method. Documentation for this feature can be found in `data_override/README.md`. (#1556)
- PLATFORM: Adds two constants `FMS_PATH_LEN` and `FMS_FILE_LEN` and uses it across the code for any file path or file name character strings. They default to 1024 and 255 but can also be set by the `FMS_MAX_PATH_LEN` and `FMS_MAX_FILE_LEN` CPP macros. (#1567)
- CMAKE: Adds Cmake option to support building shared libraries (#1559)

### Changed
- DIAG_MANAGER: Simplifies the diag_table.yaml format by allowing `module`, `reduction`, and `kind` to be set on a file level, with the ability to override for a specific field. (#1545)
- FIELD_MANAGER: Updated and refactored the `fm_yaml_mod` module for a new table format to remove the `subparams` key. (#1547)
- DATA_OVERRIDE: Updates the yaml format to be improve readablity and to be more consistent in key names. Additional documentation for the new format can be found in `data_override/README.md`. (#1556)

### Fixed
- PARSER: Adds error code checks to the yaml parser to prevent code from hanging on invalid yamls (#1563)
- PARSER: Adds ability to read in generic blocks to support `field_manager` updates. (#1519)
- CRAY COMPILER SUPPORT: Updated any multi-line string literals to be compatible with the cray compiler. (#1554)

### Removed
- DIAG_MANAGER: Removed diag_table schemas from FMS and moved them to the `gfdl_msd_schemas` repository, and updates the `diag_manager` documentation markdowns. (#1543)

### Tag Commit Hashes
- 2024.03-beta1 a5de6a54abeb00be2443db4cf07aa267b7faa724

## [2024.02] - 2024-07-11

### Known Issues
- Diag Manager Rewrite: See [below](#20240102---2024-06-14) for known output file differences regarding the new diag manager. The new diag_manager is disabled by default, so this differences will only be present if `use_modern_diag` is set to true in the `diag_manager_nml`.
- BUILD(HDF5): HDF5 version 1.14.3 generates floating point exceptions, and will cause errors if FMS is built with FPE traps enabled. FPE traps are turned on when using the debug target in mkmf.
- GCC: version 14.1.0 is unsupported due to a bug with strings that has come up previously in earlier versions. This will be caught by the configure script, but will cause compilation errors if using other build systems.

### Added
- TIME_INTERP: Enables use of `verbose` option in `time_interp_external2` calls from `data_override`. The option is enabled in `data_override_nml` by setting `debug_data_override` to true. (#1516)
- COUPLER: Adds optional argument to `coupler_types_send_data` routine that contains the return statuses for any calls made to the diag_manager's `send_data` routine. (#1530)
- MPP: Adds a separate error log file `warnfile.<root pe num>.out` that only holds output from any `mpp_error` calls made during a run (#1544)
### Changed
- DIAG_MANAGER: The `diag_field_log.out` output file of all registered fields will now include the PE number of the root PE at the time of writing (ie. diag_field_log.out.0). This is to prevent overwritting the file in cases where the root PE may change. (#1497)

### Fixed
- CMAKE: Fixes real kind flags being overwritten when using the Debug release type (#1532)
- HORIZ_INTERP: Fixes allocation issues when using method-specific horiz_interp_new routines (such as `horiz_interp_bilinear_new`) by setting `is_allocated` and the `method_type` during initialization for each method. (#1538)


### Tag Commit Hashes
- 2024.02-alpha1 5757c7813f1170efd28f5a4206395534894095b4
- 2024.02-alpha2 5757c7813f1170efd28f5a4206395534894095b4
- 2024.02-beta1  ca592ef8f47c246f4dc56d348d62235bd0ceaa9d
- 2024.02-beta2  ca592ef8f47c246f4dc56d348d62235bd0ceaa9d

## [2024.01.02] - 2024-06-14

### Known Issues
- Diag Manager Rewrite:
	- Expected output file changes:
		- If the model run time is less than the output frequency, old diag_manager would write a specific value (9.96921e+36). The new diag_manager will not, so only fill values will be present.
		- A `scalar_axis` dimension will not be added to scalar variables
		- The `average_*` variables will no longer be added as they are non-standard conventions
		- Attributes added via `diag_field_add_attributes` in the old code were saved as `NF90_FLOAT` regardless of precision, but will now be written as the precision that is passed in
		- Subregional output will have a global attribute `is_subregional = True` set for non-global history files.
		- The `grid_type` and `grid_tile` global attributes will no longer be added for all files, and some differences may be seen in the exact order of the `associated_files` attribute

- DIAG_MANAGER: When using the `do_diag_field_log` nml option, the output log file may be ovewritten if using a multiple root pe's
- BUILD(HDF5): HDF5 version 1.14.3 generates floating point exceptions, and will cause errors if FMS is built with FPE traps enabled.
- GCC: version 14.1.0 is unsupported due to a bug with strings that has come up previously in earlier versions. This will be caught by the configure script, but will cause compilation errors if using other build systems.

### Fixed
- DIAG_MANAGER: Fixes incorrect dates being appended to static file names

## [2024.01.01] - 2024-05-30

### Known Issues
- Diag Manager Rewrite:
	- Expected output file changes:
		- If the model run time is less than the output frequency, old diag_manager would write a specific value (9.96921e+36). The new diag_manager will not, so only fill values will be present.
		- A `scalar_axis` dimension will not be added to scalar variables
		- The `average_*` variables will no longer be added as they are non-standard conventions
		- Attributes added via `diag_field_add_attributes` in the old code were saved as `NF90_FLOAT` regardless of precision, but will now be written as the precision that is passed in
		- Subregional output will have a global attribute `is_subregional = True` set for non-global history files.
		- The `grid_type` and `grid_tile` global attributes will no longer be added for all files, and some differences may be seen in the exact order of the `associated_files` attribute

- DIAG_MANAGER: When using the `do_diag_field_log` nml option, the output log file may be ovewritten if using a multiple root pe's
- BUILD(HDF5): HDF5 version 1.14.3 generates floating point exceptions, and will cause errors if FMS is built with FPE traps enabled.
- GCC: version 14.1.0 is unsupported due to a bug with strings that has come up previously in earlier versions. This will be caught by the configure script, but will cause compilation errors if using other build systems.

### Added
- DIAG_MANAGER: Implements `flush_nc_files` functionality from legacy diag_manager.

### Changed
- FMS2_IO: Changed `register_unlimited_compressed_axis` to use a collective gather rather than send and recieves to improve efficiency when reading in iceberg restarts.

### Fixed
- DIAG_MANAGER: Fixes 0 day output frequencies causing error stating a time_step was skipped. Also adds checks to crash if averaged fields have -1 or 0 day frequencies or if mixing averaged and non-averaged fields in the same file.
- DIAG_MANAGER: Fixes issue with the weight argument not getting passed through to reduction methods.
- DIAG_MANAGER: Allocation errors when using two empty files.
- DIAG_MANAGER: `time` and `time_bnds` being larger than expected when running for 1 day and using daily data.
- DIAG_MANAGER: Allows for mixing static and non-static fields when frequency is 0 days.
- TESTS: Fixes compile failure with ifort 2024.01 from test_mpp_gatscat.F90.

### Removed
- DIAG_MANAGER: The `mix_snapshot_average_fields` option is deprecated for the rewritten diag_manager only.

### Tag Commit Hashes
- 2024.01.01-beta2 c00367fa810960e87610162f0f012c5da724c5a9
- 2024.01.01-beta1 42f8506512e1b5b43982320f5b9d4ca1ca9cbebd

## [2024.01] - 2024-05-03

### Known Issues
- Diag Manager Rewrite:
	- If two empty files are present in the diag_table.yaml file the code will crash with a allocation error (#1506)
	- Setting an output frequency of '0 days' does not work as expected and may cause an error stating a time_step has been skipped (#1502)
	- The `flush_nc_files` and `mix_snapshot_average_fields` nml options are not yet functional. The `mix_snapshot_average_fields` option is planned to be deprecated (for the rewritten diag_manager only).
	- Expected output file changes:
		- If the model run time is less than the output frequency, old diag_manager would write a specific value (9.96921e+36). The new diag_manager will not, so only fill values will be present.
		- A `scalar_axis` dimension will not be added to scalar variables
		- The `average_*` variables will no longer be added as they are non-standard conventions
 		- Attributes added via `diag_field_add_attributes` in the old code were saved as `NF90_FLOAT` regardless of precision, but will now be written as the precision that is passed in
		- Subregional output will have a global attribute `is_subregional = True` set for non-global history files.
		- The `grid_type` and `grid_tile` global attributes will no longer be added for all files, and some differences may be seen in the exact order of the `associated_files` attribute

- DIAG_MANAGER: When using the `do_diag_field_log` nml option, the output log file may be ovewritten if using a multiple root pe's
- TESTS: `test_mpp_gatscat.F90` fails to compile with the Intel Oneapi 2024.01's version of ifort
- BUILD(HDF5): HDF5 version 1.14.3 generates floating point exceptions, and will cause errors if FMS is built with FPE traps enabled.

### Added
- DIAG_MANAGER: The diag manager has been rewritten with a object oriented design. The old diag_manager code has been kept intact and will be used by default. The rewritten diag manager can be enabled via `use_modern_diag = .true.` to your `diag_manager_nml`. New features include:
	- Self-describing YAML formatting for diag_table's
	- Allows 4d variables
  - Support defining subregions with indices
  - More flexibility when adding metadata and defining output frequency
- FMS2_IO: Adds support for collective parallel reads to improve model startup time. The collective reads are disabled by default and enabled via the `use_collective` flag in `netcdf_io_mod`.
- DATA_OVERRIDE: Adds option to use multiple data files for one field within data_override in order to use annual data files in yearly runs without having to append/prepend timesteps from previous and next year. With the legacy data_table, filenames  can be set in order and separated with `:` ie. `prev_year.nc:curr_year.nc:next_year.nc`. With the data_table.yaml format, the key `is_multi_file` enables the functionality and `prev_file_name` and `next_file_name` sets the file paths.

- INTERPOLATOR: Adds support for yearly/annual data
- DATA_OVERRIDE: Adds support for monotonically decreasing arrays for `nearest_index`, `axis_edges`, `horiz_interp`(bilinear), and `data_override` (#1388)
- DOCS: Add documentation for the exchange grid (xgrid_mod) and update the contribution guide to add a section on code reviews
- MPP: MPI sub-communicators for domains are now accessible via `mpp_get_domain_tile_commid` and `mpp_get_domain_commid` in `mpp_domains_mod`

### Changed
- DATA_OVERRIDE: Changes behavior to crash if both data_table and data_table.yaml are present and adds error checking when reading in yaml files
- FIELD_MANAGER: Changes behavior to crash if both field_table and field_table.yaml are present as well as adds a namelist flag (`use_field_table_yaml`) to enable support for the yaml input.

### Fixed
- DATA_OVERRIDE: Fixes allocation error with scalar routine and replaces pointers with allocatables
- INTERPOLATOR: Increase max string size for file paths
- AXIS_UTILS: Improves performance of `nearest_index` routine
- CMAKE: Fixes macOS linking issues with OpenMP

### Tag Commit Hashes
- 2024.01-beta5  d3bab5a84b6a51eddd46ab6fb65eaa532830c6c7
- 2024.01-beta4  ac363ddfd3075637cecae30ddfbae7a78751197b
- 2024.01-alpha6 2ace94564a08aec4d7ab7eca0e57c0289e52d5b1
- 2024.01-alpha5 5ed0bd373cc59a9681052fa837cb83a67169d102
- 2024.01-alpha4 8dd90d72b58f0de3632dc62920f8adfb996b2265
- 2024.01-beta3  f71405a075102aef42f5811dc09e239ddd002637
- 2024.01-beta2  bb6de937f70a08a440f5e63b8553b047c1921509
- 2024.01-beta1  913f8aaecca374d5e10280056de862d5e4a7a668
- 2024.01-alpha3 085c6bfc945a6f1c586b842ca6268fca442884d8
- 2024.01-alpha2 38bfde30e1cb8bf5222410a9c37e71529567bf69
- 2024.01-alpha1 ac0d086296ea8b9196552463655cb9a848db39fe

## [2023.04] - 2023-12-04
### Known Issues
- GCC 9 and below as well as GCC 11.1.0 are unsupported due to compilation issues. See prior releases for more details.
- `NO_QUAD_PRECISION` macro is no longer set by FMS, the `ENABLE_QUAD_PRECISION` macro has replaced prior usage of `NO_QUAD_PRECISION`. `-DENABLE_QUAD_PRECISION` should be set if quad precision is to be used, otherwise FMS will not use quad precision reals where applicable.

### Added
- DATA_OVERRIDE: A new namelist flag `use_data_table_yaml` has been added to enable usage of the yaml format data_override tables. This allows an executable built with yaml support be able to accept either format.

### Changed
- RESERVED KEYWORD CHANGES: Various routines in FMS have been updated to not use fortran keywords for variable names. The names changed were: `data`, `unit`, and `value`. This may affect usage of external code if argument names are explicitly used. Only required arguement names were changed to mitigate any breaking changes.
- TESTS: Changes the testing scripts to allow for the `MPI_LAUNCHER` environment variable override to work with any provided arguments.

### Fixed
- CMAKE: Fixed build issue with CMake where precision default flags were being overwritten when using GNU and MPICH.
- AUTOTOOLS: Fixes issue affecting installs where the global libFMS.F90 module was not being installed correctly and adds post-install message.
- DIAG_MANAGER: Fixes issue with incorrect start_time functionality (from the 2023.02.01 patch)

### Tag Commit Hashes
- 2023.04-beta1 be1856c45accfe2fb15953c5f51e0d58a8816882

## [2023.03] - 2023-10-27
### Known Issues
- GCC 9 and below as well as GCC 11.1.0 are unsupported due to compilation issues. See prior releases for more details.
- `NO_QUAD_PRECISION` macro is no longer set by FMS, the `ENABLE_QUAD_PRECISION` macro has replaced prior usage of `NO_QUAD_PRECISION`. `-DENABLE_QUAD_PRECISION` should be set if quad precision is to be used, otherwise FMS will not use quad precision reals where applicable.

### Added
- UNIT_TESTS: New unit tests have been created or and existing ones expanded on for any modules utilizing mixed precision support.

### Changed
- MIXED PRECISION: Most subroutines and functions in FMS have been updated to simultaneously accept both 4 byte and 8 byte reals as arguments. This deprecates the `--enable-mixed-mode` option, which enabled similar functionality but was limited to certain directories and was not enabled by default. To facilitate easier testing of these code changes, the CMake precision options for default real size were left in (along with an equivalent `--disable-r8-default` flag for autotools). The resulting libraries will support mixed-precision real kinds regardless of default real size. It should also be noted that many routines that accept real arguments have been moved to include files along with headers in order to be compiled with both kinds. Most module level variables were explicitly declared as r8_kind for these updates.
- Some type/module changes were made to facilitate mixed precision support. They are **intended** to have minimal impact to other codebases:
  - COUPLER_TYPES: In coupler_types.F90,  `coupler_nd_field_type` and `coupler_nd_values_type` have been renamed to indicate real kind value: `coupler_nd_real4/8_field_type` and `coupler_nd_real4/8_values_type`. The `bc` field within `coupler_nd_bc_type` was modified to use r8_kind within the value and field types, and an additional field added `bc_r4` to use r4_kind values.
  - TRIDIAGONAL: Module state between r4 and r8 calls are distinct (ie. subsequent calls will only be affected by calls of the same precision). This behaviour can be changed via the `save_both_kinds` optional argument to `tri_invert`.
- CODE_STYLE: has been updated to reflect the formatting used for the mixed precision support updates.

### Fixed
- DIAG_MANAGER: Tile number (ie. tileX) will now be added to filenames for sub-regional diagnostics.
- MPP: Bug affecting non-intel compilers coming from uninitialized pointer in the `nest_domain_type`
- MPP: Bug fix for unallocated field causing seg faults in `mpp_check_field`
- FMS2_IO: Fixed segfault occuring from use of cray pointer remapping along with mpp_scatter/gather
- TEST_FMS: Added various fixes for different compilers within test programs for fms2_io, mpp, diag_manager, parser, and sat_vapor_pres.
- INTERPOLATOR: Deallocates fields in the type that were previously left out in `interpolator_end`

### Removed
- CPP MACROS:
  - `no_4byte_reals` was removed and will not set any additional macros if used. `no_8byte_integers` is still functional.
  - `NO_QUAD_PRECISION` was removed. It was conditionally set if ENABLE_QUAD_PRECISION was undefined. ENABLE_QUAD_PRECISION should be used in model components instead (logic is flipped)
  - `use_netCDF` was set by autotools previously but wasn't consistently used in the code. FMS should always be compiled with netcdf installed so this was removed with the exception of its use in deprecated IO modules.
- DRIFTERS: The drifters subdirectory has been deprecated. It will only be compiled if using the `-Duse_drifters` CPP flag.

### Tag Commit Hashes
- 2023.03-beta1  06b94a7f574e7794684b8584391744ded68e2989
- 2023.03-alpha3 b25a7c52a27dfd52edc10bc0ebe12776af0f03df
- 2023.03-alpha2 9983ce308e62e9f7215b04c227cebd30fd75e784
- 2023.03-alpha1 a46bd94fd8dd1f6f021501e29179003ff28180ec


## [2023.02.01] - 2023-10-13
### Fixed
- DIAG_MANAGER: Fixes issue with incorrect start_time functionality


## [2023.02] - 2023-07-27
### Known Issues
- GCC 11.1.0 is unsupported due to compilation issues with select type. The issue is resolved in later GCC releases.
- When outputting sub-region diagnostics, the current diag_manager does not add "tileX" to the filename when using a cube sphere. This leads to trouble when trying to combine the files and regrid them (if the region is in two different tiles)
- GCC 10 and greater causing io issues when compiled using O2 optimization flags
- GNU compilers prior to the GCC 9.0 release are unsupported for this release due to lack of support for the findloc intrinsic function. This will result in an error saying 'findloc' has no IMPLICIT type and can be resolved by compiling with gcc version 9.0 or greater.

### Added
- MPP/EXCHANGE: Adds association checks before pointer deallocations in mpp includes and xgrid

### Changed
- LIBFMS: The libFMS.F90 file (module name `fms`) meant to provide global access has been updated to include 'fms' and it's module/subdirectory name as prefixes for all names. This will only affect external codes that are already using the global module (via `use fms`) and not individual modules.
- MIXED PRECISION: Updates the axis_utils2, horiz_interp, sat_vapor_pressure, and axis_utils subdirectories to support mixed precision real values.
- FMS2_IO: Added in mpp_scatter and mpp_gather performance changes from the 2023.01.01 patch. See below for more details.
- FMS2_IO: Improved error messages to give more debugging information
- FMS_MOD: Changed fms_init to include a system call to set the stack size to unlimited, removed previously added stack size fixes
- MONIN_OBUKHOV: Restructures the subroutines in `stable_mix` interface so that 1d calls the underlying implementation, and 2 and 3d call it on 1d slices of the data as opposed to passing in mismatched arrays.
- MPP: Updates from JEDI for ajoint version the mpp halo filling (mpp_do_update_ad.fh), adds checkpoint for forward buffer information.

### Fixed
- MPP: mpp_broadcast causing an unintended error message due to checking the wrong pe value
- MPP: Added workaround for GCC 12 issues causing errors with string lengths in fms2_io
- FMS2_IO: Fixed support for 'packed' data when using NF_SHORT variables. Scale_factor and add_offset attributes will now be applied if present.
- DOCS: Improved doxygen comments for tranlon, updated deployment action for site
- TESTS: Workaround added for ICE coming from mpp_alltoall test with intel 2022.3, and fixes for any test scripts missing input.nml creation. Fixes for mpp/test_global_array failures.
- TIME_INTERP: Fixes crashes when calling with a non-existant field
- DIAG_MANAGER: Fixes a module dependency issue causing failures during parallel builds
- AXIS_UTILS2: Fixes an out of bounds memory index

### Removed
- FMS_IO/MPP_IO: The two older io modules, fms_io_mod and mpp_io_mod, have been deprecated and will not be compiled by default. If you wish to compile these modules, you must use the -Duse_deprecated_io CPP flag or the --enable-deprecated-io configure option if building with autotools.

### Tag Commit Hashes
- 2023.02-beta1  2be8aa452ad3e5f43e92c38a64f12d1ae6c43fb8
- 2023.02-alpha3 8c73bd18dc1d580f2ee524c37cf903ff54d40501
- 2023.02-alpha2 783019fdec89a8db2b26247c2f63d4782e1495c0
- 2023.02-alpga1 419c66be31f82ebb13a91ea5e837c707eb54473b


## [2023.01.01] - 2023-06-06
### Changed
- FMS2_IO: Performance changes for domain_reads_2d and domain_reads_3d:
  - Root pe reads the data
  - Uses mpp_scatter to send the data to the other pes
  - Added unit tests to test all of the domain_read/domain_write interfaces

- FMS2_IO: Performance changes for compressed_writes_1d/2d/3d
  - Uses mpp_gather to get data for write
  - Added unit tests to test all of the compressed writes interfaces
  - Compressed_writes_4d/5d were unchanged

- FMS2_IO: Extended mpp_scatter and mpp_gather to work for int8; added a kludge for scatter since the data is assumed to be (x,y,z)


## [2023.01] - 2023-04-03
### Known Issues
- If using GCC 10 or higher as well as MPICH, compilation errors will occur unless `-fallow-argument-mismatch` is included in the Fortran compiler flags(the flag will now be added automatically if building with autotools or CMake).
- GCC 11.1.0 is unsupported due to compilation issues with select type. The issue is resolved in later GCC releases.
- When outputting sub-region diagnostics, the current diag_manager does not add "tileX" to the filename when using a cube sphere. This leads to trouble when trying to combine the files and regrid them (if the region is in two different tiles)

### Added
- DIAG_MANAGER: Added code refactored as part of larger diag_manager rewrite for the send_data routines. The refactored code is disabled by default and enabled by setting  `use_refactored_send` to true in the diag_manager_nml, and should mirror current behaviour.
- FMS2_IO: Added the ability to set deflate_level and shuffle netcdf options in `fms2_io_nml`. Also added functionality for registering dimensions as unlimited compressed.
- YAML_PARSER: Added support for emitting multiple tabbed section keys to allow diag manager yaml output

### Changed
- STRING_UTILS: Extended the `string` interface in fms_string_utils_mod to accept reals of 4 or 8 kind, as well as 1, 2, and 3 dimensional real arrays
- DIAG_MANAGER: Changed the `log_diag_field_info` routine to allow for specifying seperator
- INTERPOLATOR(s): In horiz_interp, amip_interp and interpolator, changed pointers arrays into allocatables

### Fixed
- TRIDIAGONAL: Added OMP directives to prevent race conditions
- DIAG_MANAGER: Added `diag_send_data` routine to fix class(\*) related compiler issues from the refactor update
- SAT_VAPOR_PRES_K: Removed implied saves causing issues with class(\*) type checking
- TIME_INTERP: Fixed naming conflicts between module level and local variables
- YAML_PARSER: Fixed typo in variable name, rename variables to avoid fortran keywords
- DOCS: Fixed incorrect serial build instructions
- COMPILER SUPPORT: Fixed compilation errors with Intel's llvm-based compiler and added support for the CMake build. Also fixed mpp_checksum unit test failures with openmpi and nvhpc compilation issues.
- TIME_MANAGER: Fixed an bug from PR #1169 that was causing answer changes in land models

### Tag Commit Hashes
- 2023.01-beta4		(63626578cb8ed4bed1ce670b88acd6a1ec438e32)
- 2023.01-beta3		(0ff254e409b74d7d17ab234abe5ecd985967256c)
- 2023.01-beta2		(74d8e734bd43b0ce043003da74896e5d747afc2f)
- 2023.01-beta1		(6255971af28381fad22547bdc2c538fc3ea2e8bf)
- 2023.01-alpha4	(4526cc94a3e19fe8fa151f54b0db432e1fb2f7d0)
- 2023.01-alpha3	(f0e8cab3d8e58195f7c2663b84fd0bed12fa8b64)
- 2023.01-alpha2	(91e732473f7cffce070f9ce239f8ffa22c081261)
- 2023.01-alpha1	(203c8bf464ff26fe0fe39b1451caedd026bbce55)


## [2022.04] - 2022-10-13
### Known Issues
- If using GCC 10 or higher as well as MPICH, compilation errors will occur unless `-fallow-argument-mismatch` is included in the Fortran compiler flags(the flag will now be added automatically if building with autotools or CMake).
- GCC 11.1.0 is unsupported due to compilation issues with select type. The issue is resolved in later GCC releases.
- When outputting sub-region diagnostics, the current diag_manager does not add "tileX" to the filename when using a cube sphere. This leads to trouble when trying to combine the files and regrid them (if the region is in two different tiles)

### Added
- FIELD MANAGER: Adds support for reading field tables in the yaml format (field_table.yaml).
 Yaml input only be used if compiled with `-Duse_yaml` preprocessor flag, otherwise it will default
 to previous behaviour. The converter script to convert current field tables to the new format is
 publicly available [here](https://github.com/NOAA-GFDL/fms_yaml_tools/) although the conversions will also be done automatically in FRE.
- FMS2_IO: Adds options to enable data compression and chunking for netcdf output in `register_restart_field`

### Changed
- BUILD: Improves the configuration check for the MPICH/GCC 10+ argument mismatch bug by replacing
it with a m4 compile test

### Fixed
- Compiler Support: allows for compilation via the Cray/HP CCE compilers by fixing string concatenation compilation errors

### Tag Commit Hashes
- 2022.04-beta2   163cb3e434dba05933c3d2151dea5d770758a2f3
- 2022.04-beta1   1099a2890a06d279df0abe1f383b71279643bcdd
- 2022.04-alpha4  8036d8d8448b0da8416a76ee0820314da27d5711
- 2022.04-alpha3  7fafa4f7fb7a89c6f22da7ae19dc4e61a8073451
- 2022.04-alpha2  0c4b3cc98f4bce5b39c6e9f6404ea32f5bf719e5
- 2022.04-alpha1  ec57a48aeefb62b475a11cad7e30ebe460fa0d9f

## [2022.03] - 2022-08-01
### Known Issues
- If using GCC 10 or higher as well as MPICH, compilation errors will occur unless `-fallow-argument-mismatch` is included in the Fortran compiler flags(the flag will now be added automatically if building with autotools or CMake).
- GCC 11.1.0 is unsupported due to compilation issues with select type. The issue is resolved in later GCC releases.
- When outputting sub-region diagnostics, the current diag_manager does not add "tileX" to the filename when using a cube sphere. This leads to trouble when trying to combine the files and regrid them (if the region is in two different tiles)
### Added
- BUILD: Adds checks to autotools and cmake build files to solve compilation issues with GCC 10 and greater. Also adds a debug build type for CMake to allow for overriding compiler flags, and individual override flags for mixed precision routines.
- DOCS: Additional information added for building and testing FMS with the build systems; renamed and moved autotools build document.
- YAML: Adds support for writing yaml files through the `fms_yaml_output_mod` module
### Changed
- MIXED MODE: Expands support for mixed precision reals to the constants files, diag_manager, sat_vapor_pres, time_manager, and tracer_manager
- FMS_IO: Increased the character length for restart file names to allow for longer paths
### Fixed
- COUPLER: Fixes global checksum being written to stdout by every core instead of just the root
- DOCS: Fixed parsing issues with include and header files, adds class diagrams and layout improvements
### Tag Commit Hashes
- 2022.03-alpha1 62588548a5ecbdce7dbf857542ed272f7b2c971f
- 2022.03-beta1  8a4ad847122c7cc597a1f2626290b46af44b143a

## [2022.02] - 2022-04-29
### Known Issues
- If using GCC 10 or higher as well as MPICH, compilation errors will occur unless `-fallow-argument-mismatch` is included in the Fortran compiler flags
- GCC 11.1.0 is unsupported due to compilation issues with select type. The issue is resolved in later GCC releases.
- Current diag_manager does not add "tileX" to the filename when using a cube sphere, which leads to trouble when trying to combine the files and regrid them (if the region is in two different tiles)
### Added
- STRING_UTILS: Adds a module, `fms_string_utils_mod`, for common string operations throughout FMS
- LIBFMS: makes recently added routines available through the global `fms` module
- CMAKE: Adds build option for position independent code
- CONSTANTS: Adds macros to load constants for different modeling systems/uses between GFDL, GEOS and GFS. Can be selected in cmake with `-DCONSTANTS=<GEOS|GFDL|GFS>`
### Changed
- STRING_UTILS: Refactored string routine definitions from fms_mod and fms2_io_mod to be located in fms_string_utils_mod
- CONSTANTS: Makes fmsconstants.F90 contain the constant definitions, with constants_mod refactored to hold the same values
- MOSAIC2: changes grid 'version' names and documentation to be more descriptive
### Removed
- FMS_MOD: Removes fms_c.c and fms_c.h files from the fms directory
### Fixed
- FMS2_IO: Fixed bug casuing non-root pe's to fail during the flush_file routine
### Tag Commit Hashes
- 2022.02-alpha1 270c2a4e1a94229a2ae6b1e431c473589b6e15c3
- 2022.02-alpha2 7768ad1d4941b92ec8f40d34b1b517f5bde3df4e
- 2022.02-beta1  689579eea6bf7a25c64e8b823551ec588be90984

## [2022.01] - 2022-03-25
### Known Issues
- The MPICH MPI implementation is unsupported when used alongside GCC 10 or 11 due to compilation issues with the mixed precision reals. MPICH can still be used to compile FMS with GCC 9 or earlier, or with other compilers.
- GCC 11.1.0 is unsupported due to compilation issues with `select type`. The issue appears to be resolved in later GCC releases
### Added
- FMS2_IO: Added a macro `MAX_NUM_RESTART_VARS_` to allow the max amount of restart variables to be set at compile time
- TESTING: Adds a configure option, `--enable-code-coverage`, to build a code coverage report using intel's codecov
- AFFINITY: Adds an initialization check to `fms_affinity_set`, and updates test program with init/end routines
- FMS2_IO: Adds an optional argument to ignore embedded checksum checks when reading restart files
### Changed
- TESTING: Changes the testing suite scripts for various improvements such improved output, tests with input files, and adding/fixing new tests
- MPP: Change variable names in mpp to use more inclusive language
- DOCS: Updates to correct branch name and doxygen guide for functions, and adds CI information page
### Fixed
- Fixes compilation warnings throughout the code, mainly for uninitialized or unused variables
- Fixed any code not adhering to the projects style guide (mainly line length fixes) so that all future changes can be checked with a linter
- MOSAIC2: Adds `r8_kind` casts to calls to C routines in order match precision of doubles
- TESTS: Fixes crashes in fms2_io tests from namelist read errors
### Tag Commit Hashes
- 2022.01-alpha1 516a5efa681e5ae954c11c0c90677b4444e28ec4
- 2022.01-beta1  12da12884f8dc8bde47b478c997b0e5d49260a1c
- 2022.01-alpha2 28e8e3e751a6d5d81b640fb779304329f3edb82d
- 2022.01-beta2  7b78a73a5ba7acf5d3d932ecfe081e5040e2c778
## [2021.04] - 2021-12-23
### Known Issues
- GCC 11.1.0 is unsupported due to compilation issues with `select type`. The issue appears to be resolved in later GCC releases
### Added
- PARSER: Adds a parser using the libyaml C library to support yaml format input files.
  Currently implemented in data override and can be enabled with the configure option  `--with-yaml` or with CMake option `-DWITH_YAML`
- FMS: Adds an interface, `fms_c2f_string`, to convert C strings and C pointers to Fortran strings
- MPP: Adds a routine `mpp_shift_nest_domains` and a field to `nest_domain_type` to allow for modifying the position of a given nest domain
- FMS2_IO: Reintroduces the option to flush_nc_files with fms2_io
### Changed
- DIAG_MANAGER: Cleans up IO code and replaces any remaining dependencies to mpp_io with fms2_io
- FMS_IO: Changes to allow for custom paths for namelists, field_table, and the INPUT directory
- EXCHANGE: Changes real sizes in xgrid and gradient modules to be explicitly r8_kind to prevent runtime issues with mixed precision
### Deprecated
- MPP: `get_unit` has been deperecated in favor of the Fortran intrinsic `newunit` and will now generate a warning if used
### Removed
- TIME_MANAGER: Removes deprecated array-based gregorian calender calculations that were replaced in 2021.02
### Fixed
- DIAG_MANAGER: Fixes issues with 3D diurnal diagnostic output and removes a redundant write_data call
- TIME_INTERP: Fixes load_record read_data call for 3d variables with fms2_io and eliminates redundant data loading and validity checking for on-grid interpolations.
- MPP: Fixed a bug with non-blocking domain updates failing on GNU compilers from uninitialized values
- MPP: Fixed issues with the `mpp_type_free` function causing errors and memory leaks when freeing the `mpp_byte` type

### Tag Commit Hashes
- 2021.04-alpha1 (e0b998321611f80f2d0c587a13b8c03c173d5520)
- 2021.04-alpha2 (ab1b0a4cb2beac72d889d94a628e0d02092723b2)
- 2021.04-alpha3 (90583aeb369831b01296ab4b0e7e6a1b69ed91b1)
- 2021.04-beta1  (6d179fcdc189070f74d49e0025d072fa304e96d6)

## [2021.03] - 2021-08-16
### Known Issues
- DIAG_MANAGER: 3D diurnal diagnostic variables are not supported in this version of FMS
### Added
- FMS2_IO: Documentation was added for FMS2_io to help users convert from fms_io/mpp_io
### Changed
- FMS2_IO: The error messages in FMS2_io were updated to give more useful information
- TEST_FMS: The unit tests in mosaic, axis_utils, and time_interp_external were updated to use the FMS2_io version of these routines and are no longer skipped
- DOCS: The doxygen generated documentation has been improved with more doxygen comments added and a more cohesive layout
- TEST_FMS: Unit tests for time_manager were updated to use the new get/set_date_gregorian routines
### Removed
- MPP_IO: The namelist variable use_mpp_io was removed from interpolator, amip_interp, diag_manager, topography, xgrid, and data_override
- MPP_IO: Any remaining fms_io/mpp_io calls from the source and test code were removed
- FMS: Removes the hardcoded path for input.nml, path now may be specified in the call to `fms_init`
### Fixed
- MPP: Fixes algorithm used with nested grid updates to properly coalesce x-dir and y-dir pelists for vector quantities
- CMAKE/AUTOTOOLS: Fixes for minor issues with filenames in both the CMake and autotools build systems
- MPP: Restored deleted pset functionality needed by GFDL SCM by reinstating mpp_pset.F90
- MPP: Fixed uninitialized variables for data domains in mpp domains broadcast routines
- MPP: Minor memory leaks from deallocating domains
- AXIS_UTILS: Fix PGI related error with string length sizes
### Tag Commit Hashes
- 2021.03-alpha1 (87d945d8dba6341f1f56631047ae5d3e5b4ab828)
- 2021.03-beta1  (6d6ff9595ede12ea0a342ae014442708a27041d2)

## [2021.02] - 2021-05-20
### Added
- FMS2_IO: Added fms2_io support for boundary condition restarts. `register_restart_region_2d` and `register_restart_region_3d` were added to fms2_io’s `register_restart_field` interface and `read_restart_bc` and `write_restart_bc` subroutines were added to read and write boundary conditions restarts. See [test_fms/fms2_io/test_bc_restart.F90](https://github.com/NOAA-GFDL/FMS/blob/9d55115a331685e4c6e01f2dfb3b770a9f80fa37/test_fms/fms2_io/test_bc_restart.F90) for sample usage.
- FMS2_IO: Added fms2_io’s version of set_filename_appendix. The string sent in will be appended to the filename of a restart/history file before the *.nc or before *.tile if tile is in the filename.
- FMS2_IO: Added workarounds to get fms2_io code working with pgi
- COUPLER_TYPES: Added fms2_io’s version of `register_restarts_2/3d` and `CT_restore_state_2/3d` to the `coupler_type_register_restarts` and `coupler_type_restore_state` interfaces in coupler_types. The fms_io’s versions were renamed as mpp_io_* and both versions may still be used. See [test_fms/coupler/test_coupler_*d.F90](https://github.com/NOAA-GFDL/FMS/tree/9d55115a331685e4c6e01f2dfb3b770a9f80fa37/test_fms/coupler) for sample usage.
- MPP_DOMAINS: Added a subroutine `mpp_create_super_grid_domain`, which sets the indices of the input domain to match the supergrid domain.
- MOSAIC2/GRID2: Added a version of grid.F90 that is suitable for use with fms2_io called grid2.F90, and moved it into the mosaic2 folder along with mosaic2.F90
- FMS_MOD: Added grid_init and grid_end calls for use with grid2.F90
- TIME_MANAGER: New `set_date_gregorian` and `get_date_gregorian` function/subroutine that do not use coded_date and date_to_day arrays are added and set as default.  The original get_date_gregorian and set_date_gregorian have been renamed to `get_date_gregorian_old` and `set_date_gregorian_old` and can be used by adding `old_method=.true.` optional argument to set_date and get_date subroutine calls.  This is the first step out of three to remove the memory-consuming `coded_date` and `date_to_day` arrays from the time_manager_mod.
- FMS GLOBAL MODULE: Adds a new module (libFMS.F90) to be used as a global import for all supported routines/types/variables within FMS. Also adds a separate module (fmsconstants.F90) to be used for constant values from `constants_mod`.
### Changed
- FMS2_IO: Adds logic so that the domain decomposition variable attribute is only added if the io_layout is not 1,1.
- FMS_IO: Moved get_great_circle_algorithm from fms_io to grid2.F90
- FMS_IO: Moved routines and submodules that were not IO related to fms_mod
###Removed
- DIAG_MANAGER: Removed the variable attributes `_FillValue` and `missing_value` from the `time_bounds` variable to be cf compliant
- FV3GFS: Removes the unused fv3gfs directory
- MOSAIC/GRID: Removed references in grid.F90 to fms_mod(replaced with mpp_mod direct calls)
- MOSAIC/GRID: Moved mosaic/mosaic2.F90 to mosaic2 folder
###Fixed
- DIAG_MANAGER: Fixed a bug where the variable type of `Time` and `Time_bounds` were different (float vs double) when compiling with 32 bit reals
- FMS2_IO: Fixed a bug where the code was crashing when you were trying to read/write scalar variables with the domain decomposed fileobj
### Tag Commit Hashes
- 2021.02-alpha1 (1a653fcb86a826251e6c8d0a90db897377acc49e)
- 2021.02-alpha2 (81a5b6ea2559e2c31edbcab32a3230dfc31287be)
- 2021.02-beta1 (d7e564adee07073febeb263b89158f768d665689)
- 2021.02-beta2 (9d55115a331685e4c6e01f2dfb3b770a9f80fa37)

## [2021.01] - 2021-03-08
### Added
- MPP: A counter for timers to report how many times a timer section is run
- MPP: Adds missing interfaces to be consistent with interfaces that use the OVERLOAD and no_8byte_integer macros in order to allow building without MPI
- MPP: Extends interfaces for read and write routines to include 32-bit and 64-bit real data arrays
- MPP: Adds unit tests for mpp and mpp_io for all public routines with mixed-precision interfaces and expands on existing tests for mixed-precision
- Adds an .editorconfig file with the project's preferred editor configuration
- A variable MODDIR in configure.ac for use in Makefiles to find required Fortran module files
- Adds FMS description web page as a markdown file
### Changed
- DOCS: Updates various modules to doxygen style comments and makes adjustments to correctly generate doxygen documentation through the build system
- PLATFORM: changes usage of platform.h to platform_mod and it's associated data types
- Changes all previous uses of flush subroutine calls to function calls
- Changes travis CI to Github actions CI and removed all trailing whitespace
### Removed
### Fixed
- MPP: Fixed a bug causing mpp_get_UG_domain_tile_pe_inf to seg fault from the incorrect assignment of an optional argument
- FMS: Fixes issues with FMS unit tests failing from pointer allocations by reworking deallocate_unstruct_pass_type
- MPP_IO: Fixes unintentional printing of file attributes
- An issue with the automake build system causing unnecessary rebuilds of source files
- Fixes CMake build of the FMS library to install configuration files in the appropriate directories; and for OpenMP dependencies to the private
### Tag Commit Hashes
- 2021.01-alpha1 (dbe8a1060fb33167c2d12239484226b40fb01fd0)
- 2021.01-alpha2 (b94eb18fe8e686c5958cbbacc4cf9130873afc85)
- 2021.01-beta1  (4dcc9a795d9ba0cc959ebd93dda5be5f8473545c)


## [2020.04] - 2020-12-07
### Added
- DIAG_MANAGER: A namelist flag called `use_mpp_io` if set to .true. will use mpp_io. The default is .false. and will use fms2_io.
- DIAG_MANAGER: A check is added before a time axis is registered to check if the time axis is registered as a variable
- DIAG_MANAGER: A unit test was added to test the optional `new_file_freq` functionality
- XGRID: A namelist flag called `use_mpp_io` if set to .true. will use mpp_io.  The default is .false. and will use fms2_io.
- INTERPOLATOR: A namelist flag called `use_mpp_io` if set to .true. will use mpp_io. The default is .false. and will use fms2_io.
- AMIP_INTERP: A namelist flag called `use_mpp_io` if set to .true. will use mpp_io. The default is .false. and will use fms2_io.
- TOPOGRAPHY: A namelist flag called `use_mpp_io` if set to .true. will use mpp_io. The default is .false. and will use fms2_io.
- DATA_OVERRIDE: A namelist flag called `use_mpp_bug` if set to .true. will use mpp_io. The default is .false. and will use fms2_io.
- DATA_OVERRIDE: A namelist flag called `reproduce_null_char_bug_flag` if set to .true. and fms2_io is being used, it will reproduce the mpp_io bug where the axis bounds were calculated instead of read. The default is .false.
A unit test was added to test the functionality of `get_grid_version_1`
- FMS2_IO: A unit test was added to test the functionality of `get_valid` and `is_valid`
### Changed
- The autotools build has been changed to copy each subdirectory module (.mod) files to a common .mod directory located at the top of the source directory.  This change simplifies the include path specifications.
- Use F90 module files for external libraries (MPI and NetCDF) for improved interface checking, thereby removing the reliance on library header include files.
- FMS2_IO: Changed how nest file names are created to be consistent with mpp_io
- CMAKE: Changed visibility of FMS OpenMP libraries to private in order to avoid conflicts with model libraries
### Removed
- LIBFMS: The flag -Duse_mpp_io should not be used and will cause a crash
- LIBFMS: Macros and logic for interfacing to the Flexible File I/O library
- LIBFMS: Macros for SGI MIPSpro compilers, including: mpp_node function and SGI Irix specific high resolution timer
- LIBFMS: Macros for IBM AIX compilers
- LIBFMS: Files in mpp supporting the CRAY SHMEM communications library
- LIBFMS: Files in mpp for the SGI PSET approach for communication via GSM
- MPP_IO: removed left over #ifdefs from backwards compatibility changes
### Fixed
- DATA_OVERRIDE: Fixed a bug in `get_grid_version_1` where the variable_size calls were not correct
- FMS2_IO: Fixed a bug in `get_valid` where the mpp_broadcast calls were done inside `if (root_pe)` blocks
The fms2_io unit tests were modified so they can work with the AOCC compiler
- XGRID: Fixed a bug in `load_xgrid` by checking if a dimension exists before calling `get_dimension_size` to avoid `FATAL: NetCDF: Invalid dimension ID or name` crashes
- DIAG_MANAGER: Fixed a bug where files were getting written with redundant time_bounds data
- CMAKE: Fixed a bug causing the build to fail on latest versions of Ubuntu
- CMAKE: Fixed a bug that caused the build to fail from incorrect source paths
### Tag Commit Hashes
- 2020.04-alpha1 (2428bb182133b8062432ee1b15974739753ca470)
- 2020.04-alpha2 (ad9915d83a2f34610cd748fd85889ba1f0f1fc02)
- 2020.04-alpha3 (e7e48839ab25cc9405abe5d8418e1b6ee6fb2d69)
- 2020.04-beta1 (6d3ac620423e7ede412e4bb1b36c338775d95193)

## [2020.03] - 2020-10-08
### Added
- FMS2_IO: Adds header_buffer_val to the fms2io namelist which sets the netcdf header size in bytes. The default value is 16kb
- FMS2_IO: Adds netcdf_default_format to the fms2io namelist which allows the user to change the netcdf file type. The default value is 64bit.
- FMS2_IO: Adds support to read netcdf string global attributes
- FMS2_IO: Adds an optional argument to open_file, `dont_add_res_to_filename`, which indicates that the filename should not be modified (default adds .res to restart file name)
- FMS2_IO: Modifies the `register_variable_attribute` and `register_global_attribute` interfaces by adding str_len as an argument. This is a workaround to get fms2io to work with PGI because they don't support class (*) with len=*.
- FMS2_IO: Adds unit test that tests `write_data` and `read_data` when using a domain with a mask table
- FMS2_IO: Adds fms2io’s version of get_mosaic_tile_grid
- MPP_IO: Adds `-Duse_mpp_io` compile option for data_override, interpolator, amip_interp, diag_manager, topography, and xgrid to select using mpp_io instead of fms2_io
- MPP_INIT: Adds unit tests for routines/functions that are called in mpp_init
- CMAKE: A cmake build system has been added with a CI build using cmake
### Changed
- FMS2_IO: Improves performance of previous release by gathering the domain decomposed data into one global buffer and doing one write rather than doing multiple reads
- DATA_OVERRIDE: Changes line in time_interp_external2 to enable 3D overrides

### Fixed
- DATA_OVERRIDE[2]: Fixes a crash when doing ongrid data_override calls with a domain with halos
- DIAG_MANAGER[2]: Fixes an issue where time_bnds were written incorrectly for the last time stamp
- DIAG_MANAGER: Regional diagnostics with a mask table now work
- FMS2_IO: Unit test includes fms2io_init call to improve functionality
- MPP: BOZ literals that are used in variable declaration are converted to integers using the int() function.
- MPP_DOMAINS2: Fixed unit test
- FMS_IO: Changes the logic in get_tile_string to fix bug where tile numbers 9 and 99 produce an inappropriate error

### Tag Commit Hashes
- 2020.03-beta4 (4d38679c1e18e920feb03d69f8a9762eb6a047aa)
- 2020.03-beta3 (521a15135a99d1f2da7d82f238353945f82ce1dd)
- 2020.03-beta2 (3dae0dfa405d555ecc09bbd2d60a1be24461f69e)
- 2020.03-beta1 (f7f1c1c73c1f478a53e84caee6aff2fa840ad086)
- 2020.03-alpha1 (2dd30b7ca0ac75a4a38b969e4a6d446eb395b4dd)


## [2020.02] - 2020-05-01
### Added
- FMS2_IO:  An fms2_io_nml namelist has been created.  It includes the variable ncchksz. This is the replacement for the environment variable NC_BLKSZ set in model run scripts and used by mpp_io.  The default value is 64 KB.  Any time a file is opened in fms2_io (nf90_open or nf90_create), the optional argument `chunksize=ncchksz` is passed to the NetCDF library.  NetCDF attempts to use this value to control the blocksize utilized for reads and writes of data from the filesystem.
- FMS2_IO: Adds support to `compute_global_cheksum.inc` for `real32`, assuming the flag `-DOVERLOAD_R4` is used when compiling.

### Changed
- MPP_DOMAINS - nesting:  The logic supporting nested domains for mosaic grids has been overhauled and extended.  FMS now supports multiple nests and telescoping nests (nest embedded within a nest).  The requirement for a nest to lie wholly within a single tile has been relaxed and a first-level nest may cross tile boundaries, but may not contain a tile corner.  Communications for two-way nesting have also been improved.  The 2019 December Public Release and 2020.02 GFDL Release within the  [GFDL_atmos_cubed_sphere] (https://github.com/NOAA-GFDL/GFDL_atmos_cubed_sphere) have been updated and are compatible with this release of FMS.
- FMS2_IO:  The intent of fileobj is changed from (in) to (inout) in netcdf_restore_state_wrap, restore_domain_state, and netcdf_restore_state because the file object type has a pointer that is being reassigned in one of the routines lower in the call stack.

### Deprecated
- MPP_DOMAINS - nesting:  The initial nesting implementation is no longer supported.  Please see the Changed::MPP_DOMAINS sub-entry under.

### Removed
- GENERAL:  References to the macro _ALLOCATABLE have been replaced with “allocatable”, _ALLOCATED has been replaced with “allocated”, and _NULL has been removed.  It is now assumed that all compilers support the Fortran 2003 standard.  The macros still exist in fms_platforms.h for compatibility within other components.
- DIAG_MANAGER:  “fms_platform.h” is no longer included in any of the diag_manager routines.  Instead, fms_platform_mod is now being use-associated where necessary.  This fixes an issue for debuggers not providing correct line numbers.

### Tag Commit Hashes
- 2020.02-beta1 (bbc6f8d33cfb75a411bbcd3f8423fa74b8b7cdfd)
- 2020.02-beta2 (6242941a632f6e261234f3a575e59efa1bfb1b36)
- 2020.02-beta3 (e6fd03b070eb1ba8ffbca740fd69681ee16d6be7)


## [2020.01] - 2020-03-13
### Added
- Adds the modules `axis_utils2`, `mosaic2`, and `time_interp_external2` that use `fms2_io`
- Adds unit tests for thread affinity in test_fms/affinity
- Autotools unit tests are now run with `srun, mpirun`, or `aprun` (whichever is present on your system)
- `fms2_io` provides three new derived types, which target the different I/O paradigms used in GFDL GCMs:
  - type(FmsNetcdfFile_t) - This type provides a thin wrapper over the netCDF4 library, but allows the user to assign a “pelist” to the file.  If a pelist is assigned, only the first rank on the list directly interacts with the netCDF library, and performs broadcasts to relay the information to the rest of the ranks on the list.  This derived type also allows the user to perform “compressed” reads/writes and restart variable “registration” (i.e. the ability to store pointers to user-defined/allocated buffers) and reads/writes in a manner similar to the existing calls in fms_io.
  - type(FmsNetcdfUnstructuredDomainFile_t) - This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user defined mpp_domains unstructured grid.  Users are required to inform the FMS2_io module which dimensions correspond to the unstructured grid using the appropriate `register_axis` call before any domain-decomposed reads/writes can be performed.
  - type(FmsNetcdfDomainFile_t) - This type extends the FmsNetcdfFile_t type to support “domain-decomposed” reads/writes on a user-defined `mpp_domains` two-dimensional lon-lat or cubed-sphere grid.  Users are required to inform the `fms2_io` module which dimensions (‘x’ or ‘y’) correspond to the two-dimensional grid using the appropriate `register_axis` calls before any domain-decomposed reads/writes can be performed.
- `fms2_io` requires the user to manage objects of the types described above.  Calls to `open_file` (`close_file`) act as constructors (destructors) for each of the objects, and must be explicitly made by the user.  Each object must be constructed before any I/O can be performed.  Examples describing how to use these new types are available in [test_fms/fms2_io](https://github.com/NOAA-GFDL/FMS/tree/master/test_fms/fms2_io)
- `fms2_io` treats *_FillValue* attributes as *valid_max* or *valid_min* range specifiers if none of *valid_range*, *valid_min*, and *valid_max* are specified as described in these [netcdf conventions](https://www.unidata.ucar.edu/software/netcdf/docs/attribute_conventions.html)

### Changed
- The `diag_manager` IO is handled by fms2_io instead of mpp_io. Default behavior assumes that the mpp_io namelist variable is set to *cf_compliant = .true.*
- The user must specify the diagnostic attributes that they want to write to the output files.  Example: If there is no *units* attribute, then the variable metadata will not include *units*, and it will not automatically add *units = “none”*.
- Calls to `register_diag_axis` for an X or Y axis that is shifted from the *CENTER* position need to include the optional argument *domain_position* and should be equal to *EAST* or *NORTH* based on the position relative to the domain. EAST and *NORTH* are exposed through `diag_manager_mod`.
- Changed the handling of *average_T* and *time_bnds* variables so that they are set to values that are sent in and are not manipulated as was the case in mpp_io.
- `interpolator`, `xgrid`, `data_override`, and `amip_interp` now call  `fms2_io` routines
- Support for enabling/disabling quad-precision (used in certain calculations) has been changed to a hard on/off switch.  Default behavior is quad-precision disabled. To enable, add the following CPP macro -DENABLE_QUAD_PRECISION.  This change was necessary to remove guessing at the proper setting via a mix of compiler vendor and system-defined environment variables which resulted in different behaviors on machines unbeknownst to the user.

### Deprecated
- fms2_io does **NOT** use the *scale_factor*, *add_offset*, or other attributes to manipulate the data.  The variable/data is returned to the caller as it appears in the file.  All post-read data manipulations should be handled by the caller.
### Removed
- Removes the use of bats when running unit tests

### Tag Commit Hashes
- 2020.01-alpha1 (09dc8e9e0f1c852e9e9190834176d16943cd3729)
- 2020.01-beta1 (e1c0d9d01d844938adc0d18afa09532f336bcdfe)
- 2020.01-beta2 (b68de5382a5ce631ddd6167de8d85f7c9ae54351)

## [2019.01] - 2019-11-26
### Added
- switch from "city" versioning style to `yyyy.<2_digit_version_number>[.<2-digit-patch number>]` style
- main development branch is `master` instead of `dev/master`
- affinity handling moved to the affinity directory
- support for building with autotools
- [fms_platform.h](include/fms_platform.h) contains directives that support building on macOS
- unit and build tests are available in the [test_fms](test_fms) directory
- updated [fv3gfs/makefile](fv3gfs/makefile) for use with current EMC build system
### Fixed
- Fixed `time_interp_missing` parameter in [time_interp/time_interp_external.F90](time_interp/time_interp_external.F90) to be within range when compiled in mixed-mode.
- reverted `QUAD_PRECISION` cpp macro behavior to pre-Xanadu behavior in [include/fms_platform.h](include/fms_platform.h)
- Fixed a GNU compiler issue with the logical check to set the netCDF fill value in `mpp_io_write` by separating the logical `.AND.` into nested `if` statements.
