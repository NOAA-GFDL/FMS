# ChangeLog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0),
and this project uses `yyyy.rr[.pp]`, where `yyyy` is the year a patch is released,
`rr` is a sequential release number (starting from `01`), and an optional two-digit
sequential patch number (starting from `01`).

## [2021.03.01] - 2021-09-07
### Fixed
- TIME_INTERP: Fixes issue in load_record when reading 3d variables with fms2_io and elimates redundant loads and validity checks for on-grid interpolations
### Changed
- Changes configure script to only check for gcc 11.1.0, not any gcc 11 version

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
