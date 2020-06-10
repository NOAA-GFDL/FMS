#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

# This is an automake file for the libFMS directory of the FMS
# package. This is the final packaging of the library.

# Ed Hartnett 2/22/19

# This builds the FMS library file.
lib_LTLIBRARIES = libFMS.la

# These linker flags specify libtool version info.
# See http://www.gnu.org/software/libtool/manual/libtool.html#Libtool-versioning
# for information regarding incrementing `-version-info`.
libFMS_la_LDFLAGS = -version-info 4:0:0

# Add the convenience libraries to the FMS library.
libFMS_la_LIBADD = ${top_builddir}/platform/libplatform.la
libFMS_la_LIBADD += ${top_builddir}/constants/libconstants.la
libFMS_la_LIBADD += ${top_builddir}/tridiagonal/libtridiagonal.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_parameter.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_data.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_utilities.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_domains.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_io.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_pset.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_efp.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_c.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_memutils.la
libFMS_la_LIBADD += ${top_builddir}/mpp/libmpp_memuse_c.la
libFMS_la_LIBADD += ${top_builddir}/memutils/libmemutils.la
libFMS_la_LIBADD += ${top_builddir}/fms/libfms_io.la
libFMS_la_LIBADD += ${top_builddir}/fms/libfms.la
libFMS_la_LIBADD += ${top_builddir}/fms2_io/libfms2_io.la
libFMS_la_LIBADD += ${top_builddir}/fms2_io/libfms_io_utils.la
libFMS_la_LIBADD += ${top_builddir}/fms2_io/libnetcdf_io.la
libFMS_la_LIBADD += ${top_builddir}/fms2_io/libfms_netcdf_domain_io.la
libFMS_la_LIBADD += ${top_builddir}/fms2_io/libfms_netcdf_unstructured_domain_io.la
libFMS_la_LIBADD += ${top_builddir}/fms2_io/libblackboxio.la
libFMS_la_LIBADD += ${top_builddir}/affinity/libfms_affinity.la
libFMS_la_LIBADD += ${top_builddir}/affinity/libfms_affinity_c.la
libFMS_la_LIBADD += ${top_builddir}/mosaic/libmosaic.la
libFMS_la_LIBADD += ${top_builddir}/mosaic/libmosaic2.la
libFMS_la_LIBADD += ${top_builddir}/mosaic/libgrid.la
libFMS_la_LIBADD += ${top_builddir}/mosaic/libgradient.la
libFMS_la_LIBADD += ${top_builddir}/mosaic/libmosaic_c.la
libFMS_la_LIBADD += ${top_builddir}/coupler/libcoupler_types.la
libFMS_la_LIBADD += ${top_builddir}/coupler/libensemble_manager.la
libFMS_la_LIBADD += ${top_builddir}/coupler/libatmos_ocean_fluxes.la
libFMS_la_LIBADD += ${top_builddir}/drifters/libdrifters_core.la
libFMS_la_LIBADD += ${top_builddir}/drifters/libdrifters_input.la
libFMS_la_LIBADD += ${top_builddir}/drifters/libdrifters_comm.la
libFMS_la_LIBADD += ${top_builddir}/drifters/libdrifters_io.la
libFMS_la_LIBADD += ${top_builddir}/drifters/libdrifters.la
libFMS_la_LIBADD += ${top_builddir}/drifters/libcloud_interpolator.la
libFMS_la_LIBADD += ${top_builddir}/drifters/libquicksort.la
libFMS_la_LIBADD += ${top_builddir}/axis_utils/libaxis_utils.la
libFMS_la_LIBADD += ${top_builddir}/axis_utils/libaxis_utils2.la
libFMS_la_LIBADD += ${top_builddir}/horiz_interp/libhoriz_interp_type.la
libFMS_la_LIBADD += ${top_builddir}/horiz_interp/libhoriz_interp_bicubic.la
libFMS_la_LIBADD += ${top_builddir}/horiz_interp/libhoriz_interp_bilinear.la
libFMS_la_LIBADD += ${top_builddir}/horiz_interp/libhoriz_interp_conserve.la
libFMS_la_LIBADD += ${top_builddir}/horiz_interp/libhoriz_interp_spherical.la
libFMS_la_LIBADD += ${top_builddir}/horiz_interp/libhoriz_interp.la
libFMS_la_LIBADD += ${top_builddir}/time_manager/libtime_manager.la
libFMS_la_LIBADD += ${top_builddir}/time_manager/libget_cal_time.la
libFMS_la_LIBADD += ${top_builddir}/time_interp/libtime_interp.la
libFMS_la_LIBADD += ${top_builddir}/time_interp/libtime_interp_external.la
libFMS_la_LIBADD += ${top_builddir}/time_interp/libtime_interp_external2.la
libFMS_la_LIBADD += ${top_builddir}/diag_manager/libdiag_data.la
libFMS_la_LIBADD += ${top_builddir}/diag_manager/libdiag_axis.la
libFMS_la_LIBADD += ${top_builddir}/diag_manager/libdiag_grid.la
libFMS_la_LIBADD += ${top_builddir}/diag_manager/libdiag_output.la
libFMS_la_LIBADD += ${top_builddir}/diag_manager/libdiag_util.la
libFMS_la_LIBADD += ${top_builddir}/diag_manager/libdiag_table.la
libFMS_la_LIBADD += ${top_builddir}/diag_manager/libdiag_manager.la
libFMS_la_LIBADD += ${top_builddir}/diag_integral/libdiag_integral.la
libFMS_la_LIBADD += ${top_builddir}/data_override/libdata_override.la
libFMS_la_LIBADD += ${top_builddir}/column_diagnostics/libcolumn_diagnostics.la
libFMS_la_LIBADD += ${top_builddir}/block_control/libblock_control.la
libFMS_la_LIBADD += ${top_builddir}/astronomy/libastronomy.la
libFMS_la_LIBADD += ${top_builddir}/field_manager/libfield_manager.la
libFMS_la_LIBADD += ${top_builddir}/field_manager/libfm_util.la
libFMS_la_LIBADD += ${top_builddir}/monin_obukhov/libmonin_obukhov.la
libFMS_la_LIBADD += ${top_builddir}/monin_obukhov/libmonin_obukhov_kernel.la
libFMS_la_LIBADD += ${top_builddir}/interpolator/libinterpolator.la
libFMS_la_LIBADD += ${top_builddir}/amip_interp/libamip_interp.la
libFMS_la_LIBADD += ${top_builddir}/exchange/libstock_constants.la
libFMS_la_LIBADD += ${top_builddir}/exchange/libxgrid.la
libFMS_la_LIBADD += ${top_builddir}/topography/libtopography.la
libFMS_la_LIBADD += ${top_builddir}/tracer_manager/libtracer_manager.la
libFMS_la_LIBADD += ${top_builddir}/sat_vapor_pres/libsat_vapor_pres_k.la
libFMS_la_LIBADD += ${top_builddir}/sat_vapor_pres/libsat_vapor_pres.la
libFMS_la_LIBADD += ${top_builddir}/random_numbers/libMersenneTwister.la
libFMS_la_LIBADD += ${top_builddir}/random_numbers/librandom_numbers.la

# At least one source file must be included to please Automake.
libFMS_la_SOURCES = ${top_builddir}/include/file_version.h
