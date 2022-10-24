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
!! This header file is used exclusively for doxygen documentation
!! Defines groups for each subdirectory to add their modules into
!! Additional doxygen pages can be added here as well

!> @defgroup affinity Affinity
!> @brief Modules and associated files in the affinity directory

!> @defgroup astronomy Astronomy
!> @brief Modules and associated files in the astronomy directory

!> @defgroup axis_utils Axis Utilities
!> @brief Modules and associated files in the axis_utils directory

!> @defgroup amip_interp AMIP Interpolator
!> @brief Modules and associated files in the amip_interp directory

!> @defgroup block_control Block Control
!> @brief Modules and associated files in the block_control directory

!> @defgroup column_diagnostics Column Diagnostics
!> @brief Modules and associated files in the column_diagnostics directory

!> @defgroup constants Constants
!> @brief Modules and associated files in the constants directory

!> @defgroup coupler Coupler
!> @brief Modules and associated files in the coupler directory

!> @defgroup data_override Data Override
!> @brief Modules and associated files in the data_override directory

!> @defgroup diag_integral Diag Integral
!> @brief Modules and associated files in the diag_integral directory

!> @defgroup diag_manager Diag Manager
!> @brief Modules and associated files in the diag_manager directory.
!! See below for additional information on diag_tables.

!> @defgroup drifters Drifters
!> @brief Modules and associated files in the drifters directory

!> @defgroup exchange Exchange
!> @brief Modules and associated files in the exchange directory

!> @defgroup field_manager Field Manager
!> @brief Modules and associated files in the field_manager directory

!> @defgroup fms FMS
!> @brief Modules and associated files in the fms directory

!> @defgroup fms2_io FMS2 IO
!> @brief Modules and associated files in the fms2_io directory
!!
!> Updated IO modules for parallel IO via netcdf files. Replaces the functionality of the IO
!! routines in mpp_io. fms2_io_mod is the main module for external usage and provides public
!! interfaces for routines defined throughout this directory, dependent on the
!! type of file.

!> @defgroup horiz_interp Horizontal Interpolator
!> @brief Modules and associated files in the horiz_interp directory

!> @defgroup interpolator Interpolator
!> @brief Modules and associated files in the interpolator directory

!> @defgroup memutils Memory Utilities
!> @brief Modules and associated files in the memutils directory

!> @defgroup monin_obukhov Monin Obukhov
!> @brief Modules and associated files in the monin_obukhov directory

!> @defgroup mosaic Mosaic
!> @brief Modules and associated files in the mosaic directory

!> @defgroup mosaic2 Mosaic2
!> @brief Modules and associated files in the mosaic2 directory
!!
!> Provides a fms2_io equivalent to the mpp_io dependent routines in mosaic

!> @defgroup mpp MPP
!> @brief Modules and associated files in the mpp directory
!!
!> Provides interfaces to facilitate common tasks for parallel computing and
!! modeling. Originally written to provide interfaces across different message-passing libraries,
!! it is now used mainly with the MPI standard, and provides wrappers to many MPI
!! routines. Documentation on serial implementations of routines (usually denoted by _nocommm in
!! the file name) will be excluded on these module pages, as they are currently unused but mo, but their documentation is still
!! available in the Files tab.

!> @defgroup platform Platform
!> @brief Modules and associated files in the platform directory

!> @defgroup random_numbers Random Numbers
!> @brief Modules and associated files in the random_numbers directory

!> @defgroup sat_vapor_pres Saturation Vapor Pressure
!> @brief Modules and associated files in the sat_vapor_pres directory

!> @defgroup time_interp Time Interpolator
!> @brief Modules and associated files in the time_interp directory

!> @defgroup time_manager Time Manager
!> @brief Modules and associated files in the time_manager directory

!> @defgroup topography Topography
!> @brief Modules and associated files in the topography directory

!> @defgroup string_utils String Utils
!> @brief Modules and associated files in the string_utils directory

!> @defgroup tracer_manager Tracer Manager
!> @brief Modules and associated files in the tracer_manager directory

!> @defgroup tridiagonal Tridiagonal
!> @brief Modules and associated files in the tridiagonal directory

!> @defgroup libfms FMS Global Module
!> @brief Modules and associated files in the libfms directory

!> @defgroup parser Parser
!> @brief Modules and associated files for the yaml parser

!> @page build Building and Installation
!> @brief Information about the build systems and how to build and install
!! @subpage install
!!
!! @subpage autotools
!!
!! @subpage cmake
