# ChangeLog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0),
and this project uses `yyyy.rr[.pp]`, where `yyyy` is the year a patch is released,
`rr` is a sequential release number (starting from `01`), and an optional two-digit
sequential patch number (starting from `01`).

## [Unreleased]
### Added
### Changed
- The diag_manager IO is now handled by fms2_io instead of mpp_io.  This always acts as if the the mpp_io namelist variabe was set to `cf_compliant = .true.`
- Calls to register_diag_axis using an X or Y axis that is shifted from the `CENTER` position need to include the optional argument `domain_position` and should be equal to `EAST` or `NORTH` based on the position relative to the domain.  `EAST` and `NORTH` are exposed through `diag_manager_mod`.
### Deprecated
### Removed
### Fixed

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
