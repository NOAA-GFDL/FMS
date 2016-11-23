# Makefile created by mkmf 19.2.0

SRCROOT = ../../

BUILDROOT = ./

CPPDEFS = -DINTERNAL_FILE_NML -g -traceback -Duse_libMPI -Duse_netCDF -D_FILE_VERSION="`git-version-string $<`"

#OTHERFLAGS = -Dtest_diag_manager

MK_TEMPLATE = /ncrc/home2/fms/local/opt/fre-commands/bronx-11/site/ncrc3/intel.mk
include $(MK_TEMPLATE)


.DEFAULT:
	-echo $@ does not exist.
all: libfms.a
MersenneTwister.o: $(SRCROOT)shared/random_numbers/MersenneTwister.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/random_numbers/MersenneTwister.F90
affinity.o: $(SRCROOT)shared/mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/mpp/affinity.c
amip_interp.o: $(SRCROOT)shared/amip_interp/amip_interp.F90 $(SRCROOT)shared/include/file_version.h time_interp.o time_manager.o get_cal_time.o mpp_io.o horiz_interp.o fms.o fms_io.o constants.o platform.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/amip_interp/amip_interp.F90
astronomy.o: $(SRCROOT)shared/astronomy/astronomy.F90 $(SRCROOT)shared/include/file_version.h fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/astronomy/astronomy.F90
atmos_ocean_fluxes.o: $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)shared/include/file_version.h mpp.o fms.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90
axis_utils.o: $(SRCROOT)shared/axis_utils/axis_utils.F90 $(SRCROOT)shared/include/file_version.h mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/axis_utils/axis_utils.F90
block_control.o: $(SRCROOT)shared/block_control/block_control.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/block_control/block_control.F90
cloud_interpolator.o: $(SRCROOT)shared/drifters/cloud_interpolator.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/drifters/cloud_interpolator.F90
column_diagnostics.o: $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 $(SRCROOT)shared/include/file_version.h mpp_io.o fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/column_diagnostics/column_diagnostics.F90
constants.o: $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)shared/include/file_version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/constants/constants.F90
coupler_types.o: $(SRCROOT)shared/coupler/coupler_types.F90 $(SRCROOT)shared/include/file_version.h fms.o field_manager.o mpp.o fm_util.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/coupler/coupler_types.F90
create_xgrid.o: $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)shared/mosaic/mosaic_util.h $(SRCROOT)shared/mosaic/create_xgrid.h $(SRCROOT)shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/create_xgrid.c
data_override.o: $(SRCROOT)shared/data_override/data_override.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h platform.o constants.o mpp_io.o mpp.o horiz_interp.o time_interp_external.o fms_io.o fms.o axis_utils.o mpp_domains.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/data_override/data_override.F90
diag_axis.o: $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_axis.F90
diag_data.o: $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_data.F90
diag_grid.o: $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_grid.F90
diag_manager.o: $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h time_manager.o mpp_io.o mpp.o fms.o fms_io.o diag_axis.o diag_util.o diag_data.o diag_table.o diag_output.o diag_grid.o diag_manifest.o constants.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_manager.F90
diag_manifest.o: $(SRCROOT)shared/diag_manager/diag_manifest.F90 diag_data.o mpp.o fms.o fms_io.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_manifest.F90
diag_output.o: $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_output.F90
diag_table.o: $(SRCROOT)shared/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/diag_manager/diag_table.F90
diag_util.o: $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/diag_manager/diag_util.F90
drifters.o: $(SRCROOT)shared/drifters/drifters.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/drifters/fms_switches.h $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/drifters/drifters_push.h $(SRCROOT)shared/drifters/drifters_set_field.h $(SRCROOT)shared/drifters/drifters_compute_k.h mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/drifters	$(SRCROOT)shared/drifters/drifters.F90
drifters_comm.o: $(SRCROOT)shared/drifters/drifters_comm.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/drifters/fms_switches.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/drifters	$(SRCROOT)shared/drifters/drifters_comm.F90
drifters_core.o: $(SRCROOT)shared/drifters/drifters_core.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/drifters/drifters_core.F90
drifters_input.o: $(SRCROOT)shared/drifters/drifters_input.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/drifters/drifters_input.F90
drifters_io.o: $(SRCROOT)shared/drifters/drifters_io.F90 $(SRCROOT)shared/include/file_version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/drifters/drifters_io.F90
ensemble_manager.o: $(SRCROOT)shared/coupler/ensemble_manager.F90 fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/coupler/ensemble_manager.F90
fft.o: $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)shared/include/file_version.h platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/fft/fft.F90
fft99.o: $(SRCROOT)shared/fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/fft/fft99.F90
field_manager.o: $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/field_manager	$(SRCROOT)shared/field_manager/field_manager.F90
fm_util.o: $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)shared/include/file_version.h field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/field_manager/fm_util.F90
fms.o: $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)shared/include/file_version.h mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/fms/fms.F90
fms_io.o: $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/fms/read_data_2d.inc $(SRCROOT)shared/fms/read_data_3d.inc $(SRCROOT)shared/fms/read_data_4d.inc $(SRCROOT)shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/fms	$(SRCROOT)shared/fms/fms_io.F90
gaussian_topog.o: $(SRCROOT)shared/topography/gaussian_topog.F90 $(SRCROOT)shared/include/file_version.h fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/topography/gaussian_topog.F90
get_cal_time.o: $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)shared/include/file_version.h fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/time_manager/get_cal_time.F90
gradient.o: $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)shared/include/file_version.h mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mosaic/gradient.F90
gradient_c2l.o: $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)shared/mosaic/constant.h $(SRCROOT)shared/mosaic/mosaic_util.h $(SRCROOT)shared/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/gradient_c2l.c
grid.o: $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)shared/include/file_version.h mpp.o constants.o fms.o fms_io.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mosaic/grid.F90
horiz_interp.o: $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)shared/include/file_version.h fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)shared/include/file_version.h mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)shared/include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)shared/include/file_version.h mpp.o fms.o fms_io.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/include/file_version.h mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
interp.o: $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)shared/mosaic/mosaic_util.h $(SRCROOT)shared/mosaic/interp.h $(SRCROOT)shared/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/interp.c
interpolator.o: $(SRCROOT)shared/interpolator/interpolator.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h mpp.o mpp_io.o mpp_domains.o diag_manager.o fms.o horiz_interp.o time_manager.o time_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/interpolator/interpolator.F90
memuse.o: $(SRCROOT)shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/memutils/memuse.c
memutils.o: $(SRCROOT)shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/memutils/memutils.F90
mosaic.o: $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)shared/include/file_version.h fms.o mpp.o mpp_io.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mosaic/mosaic.F90
mosaic_util.o: $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/mosaic/mosaic_util.h $(SRCROOT)shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/mosaic_util.c
mpp.o: $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/mpp/include/system_clock.h $(SRCROOT)shared/mpp/include/mpp_util.inc $(SRCROOT)shared/mpp/include/mpp_util_sma.inc $(SRCROOT)shared/mpp/include/mpp_util_mpi.inc $(SRCROOT)shared/mpp/include/mpp_util_nocomm.inc $(SRCROOT)shared/mpp/include/mpp_error_a_a.h $(SRCROOT)shared/mpp/include/mpp_error_a_s.h $(SRCROOT)shared/mpp/include/mpp_error_s_a.h $(SRCROOT)shared/mpp/include/mpp_error_s_s.h $(SRCROOT)shared/mpp/include/mpp_comm.inc $(SRCROOT)shared/mpp/include/mpp_comm_sma.inc $(SRCROOT)shared/mpp/include/mpp_transmit_sma.h $(SRCROOT)shared/mpp/include/mpp_transmit.inc $(SRCROOT)shared/mpp/include/mpp_reduce_sma.h $(SRCROOT)shared/mpp/include/mpp_sum_sma.h $(SRCROOT)shared/mpp/include/mpp_sum.inc $(SRCROOT)shared/mpp/include/mpp_alltoall_sma.h $(SRCROOT)shared/mpp/include/mpp_comm_mpi.inc $(SRCROOT)shared/mpp/include/mpp_transmit_mpi.h $(SRCROOT)shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)shared/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)shared/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)shared/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)shared/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)shared/mpp/include/mpp_chksum_int.h $(SRCROOT)shared/mpp/include/mpp_chksum_scalar.h $(SRCROOT)shared/mpp/include/mpp_chksum.h $(SRCROOT)shared/mpp/include/mpp_gather.h $(SRCROOT)shared/mpp/include/mpp_scatter.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp.F90
mpp_data.o: $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/mpp/include/mpp_data_sma.inc $(SRCROOT)shared/mpp/include/mpp_data_mpi.inc $(SRCROOT)shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp_data.F90
mpp_domains.o: $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)shared/mpp/include/mpp_domains_util.inc $(SRCROOT)shared/mpp/include/mpp_domains_comm.inc $(SRCROOT)shared/mpp/include/mpp_domains_define.inc $(SRCROOT)shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)shared/mpp/include/mpp_update_domains2D.h $(SRCROOT)shared/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)shared/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)shared/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)shared/mpp/include/mpp_do_update.h $(SRCROOT)shared/mpp/include/mpp_do_updateV.h $(SRCROOT)shared/mpp/include/mpp_do_check.h $(SRCROOT)shared/mpp/include/mpp_do_checkV.h $(SRCROOT)shared/mpp/include/mpp_update_nest_domains.h $(SRCROOT)shared/mpp/include/mpp_do_update_nest.h $(SRCROOT)shared/mpp/include/mpp_do_redistribute.h $(SRCROOT)shared/mpp/include/mpp_get_boundary.h $(SRCROOT)shared/mpp/include/mpp_do_get_boundary.h $(SRCROOT)shared/mpp/include/mpp_group_update.h $(SRCROOT)shared/mpp/include/group_update_pack.inc $(SRCROOT)shared/mpp/include/group_update_unpack.inc $(SRCROOT)shared/mpp/include/mpp_domains_reduce.inc $(SRCROOT)shared/mpp/include/mpp_global_reduce.h $(SRCROOT)shared/mpp/include/mpp_global_sum.h $(SRCROOT)shared/mpp/include/mpp_global_sum_tl.h $(SRCROOT)shared/mpp/include/mpp_global_field.h $(SRCROOT)shared/mpp/include/mpp_do_global_field.h $(SRCROOT)shared/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)shared/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)shared/mpp/include/mpp_global_field_ug.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o mpp_efp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp_domains.F90
mpp_efp.o: $(SRCROOT)shared/mpp/mpp_efp.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/mpp_efp.F90
mpp_io.o: $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/mpp/include/mpp_io_util.inc $(SRCROOT)shared/mpp/include/mpp_io_misc.inc $(SRCROOT)shared/mpp/include/mpp_io_connect.inc $(SRCROOT)shared/mpp/include/mpp_io_read.inc $(SRCROOT)shared/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)shared/mpp/include/mpp_read_compressed.h $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)shared/mpp/include/mpp_io_write.inc $(SRCROOT)shared/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)shared/mpp/include/mpp_write_compressed.h $(SRCROOT)shared/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include -I$(SRCROOT)shared/mpp/include	$(SRCROOT)shared/mpp/mpp_io.F90
mpp_memutils.o: $(SRCROOT)shared/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/mpp/mpp_memutils.F90
mpp_parameter.o: $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/mpp_parameter.F90
mpp_pset.o: $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/mpp_pset.F90
mpp_utilities.o: $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)shared/include/file_version.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/mpp_utilities.F90
nsclock.o: $(SRCROOT)shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/mpp/nsclock.c
oda_core.o: $(SRCROOT)shared/oda_tools/oda_core.F90 fms.o mpp.o mpp_domains.o time_manager.o get_cal_time.o axis_utils.o constants.o oda_types.o write_ocean_data.o mpp_io.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/oda_core.F90
oda_core_ecda.o: $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 fms.o mpp.o mpp_io.o mpp_domains.o mpp_memutils.o time_manager.o get_cal_time.o axis_utils.o horiz_interp_type.o horiz_interp_bilinear.o constants.o oda_types.o xbt_drop_rate_adjust.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/oda_core_ecda.F90
oda_types.o: $(SRCROOT)shared/oda_tools/oda_types.F90 time_manager.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/oda_types.F90
platform.o: $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/platform/platform.F90
quicksort.o: $(SRCROOT)shared/drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/drifters/quicksort.F90
random_numbers.o: $(SRCROOT)shared/random_numbers/random_numbers.F90 MersenneTwister.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/random_numbers/random_numbers.F90
read_mosaic.o: $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/mosaic/read_mosaic.h $(SRCROOT)shared/mosaic/constant.h $(SRCROOT)shared/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c -I$(SRCROOT)shared/mosaic	$(SRCROOT)shared/mosaic/read_mosaic.c
sat_vapor_pres.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)shared/include/file_version.h constants.o fms.o mpp_io.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)shared/include/file_version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
station_data.o: $(SRCROOT)shared/station_data/station_data.F90 $(SRCROOT)shared/include/file_version.h axis_utils.o mpp_io.o fms.o mpp.o mpp_domains.o diag_axis.o diag_output.o diag_manager.o diag_util.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/station_data/station_data.F90
stock_constants.o: $(SRCROOT)shared/exchange/stock_constants.F90 $(SRCROOT)shared/include/file_version.h mpp.o fms.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/exchange/stock_constants.F90
test_fms_io.o: $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/fms/test_fms_io.F90
test_horiz_interp.o: $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 mpp.o mpp_io.o mpp_domains.o horiz_interp.o axis_utils.o fms_io.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/horiz_interp/test_horiz_interp.F90
test_mpp.o: $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp.F90
test_mpp_domains.o: $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp_domains.F90
test_mpp_io.o: $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp_io.F90
test_mpp_pset.o: $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)shared/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/mpp/test_mpp_pset.F90
test_xgrid.o: $(SRCROOT)shared/exchange/test_xgrid.F90 mpp.o mpp_domains.o mpp_io.o fms.o fms_io.o constants.o xgrid.o mosaic.o gradient.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/exchange/test_xgrid.F90
threadloc.o: $(SRCROOT)shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) $(OTHER_CFLAGS) -c	$(SRCROOT)shared/mpp/threadloc.c
time_interp.o: $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)shared/include/file_version.h time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/time_interp/time_interp.F90
time_interp_external.o: $(SRCROOT)shared/time_interp/time_interp_external.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h fms.o mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o time_interp.o axis_utils.o platform.o horiz_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/time_interp/time_interp_external.F90
time_manager.o: $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/time_manager/time_manager.F90
topography.o: $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)shared/include/file_version.h gaussian_topog.o horiz_interp.o fms.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/topography/topography.F90
tracer_manager.o: $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/include/file_version.h mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/tracer_manager/tracer_manager.F90
tridiagonal.o: $(SRCROOT)shared/tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/tridiagonal/tridiagonal.F90
write_ocean_data.o: $(SRCROOT)shared/oda_tools/write_ocean_data.F90 mpp_io.o mpp.o oda_types.o time_manager.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/write_ocean_data.F90
xbt_drop_rate_adjust.o: $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 oda_types.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	$(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90
xgrid.o: $(SRCROOT)shared/exchange/xgrid.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/include/file_version.h fms.o fms_io.o mpp.o mpp_domains.o mpp_io.o constants.o mosaic.o stock_constants.o gradient.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c -I$(SRCROOT)shared/include	$(SRCROOT)shared/exchange/xgrid.F90
./mpp_write.h: $(SRCROOT)shared/mpp/include/mpp_write.h
	cp $(SRCROOT)shared/mpp/include/mpp_write.h .
./read_data_2d.inc: $(SRCROOT)shared/fms/read_data_2d.inc
	cp $(SRCROOT)shared/fms/read_data_2d.inc .
./mpp_io_write.inc: $(SRCROOT)shared/mpp/include/mpp_io_write.inc
	cp $(SRCROOT)shared/mpp/include/mpp_io_write.inc .
./system_clock.h: $(SRCROOT)shared/mpp/include/system_clock.h
	cp $(SRCROOT)shared/mpp/include/system_clock.h .
./mpp_chksum_scalar.h: $(SRCROOT)shared/mpp/include/mpp_chksum_scalar.h
	cp $(SRCROOT)shared/mpp/include/mpp_chksum_scalar.h .
./mpp_do_updateV_nonblock.h: $(SRCROOT)shared/mpp/include/mpp_do_updateV_nonblock.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_updateV_nonblock.h .
./group_update_unpack.inc: $(SRCROOT)shared/mpp/include/group_update_unpack.inc
	cp $(SRCROOT)shared/mpp/include/group_update_unpack.inc .
./mpp_global_sum_tl.h: $(SRCROOT)shared/mpp/include/mpp_global_sum_tl.h
	cp $(SRCROOT)shared/mpp/include/mpp_global_sum_tl.h .
./read_data_3d.inc: $(SRCROOT)shared/fms/read_data_3d.inc
	cp $(SRCROOT)shared/fms/read_data_3d.inc .
./mpp_domains_define.inc: $(SRCROOT)shared/mpp/include/mpp_domains_define.inc
	cp $(SRCROOT)shared/mpp/include/mpp_domains_define.inc .
./test_horiz_interp.F90: $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90
	cp $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 .
./mpp_util.inc: $(SRCROOT)shared/mpp/include/mpp_util.inc
	cp $(SRCROOT)shared/mpp/include/mpp_util.inc .
./mpp_pset.F90: $(SRCROOT)shared/mpp/mpp_pset.F90
	cp $(SRCROOT)shared/mpp/mpp_pset.F90 .
./mpp_error_a_a.h: $(SRCROOT)shared/mpp/include/mpp_error_a_a.h
	cp $(SRCROOT)shared/mpp/include/mpp_error_a_a.h .
./test_mpp_pset.F90: $(SRCROOT)shared/mpp/test_mpp_pset.F90
	cp $(SRCROOT)shared/mpp/test_mpp_pset.F90 .
./drifters_io.F90: $(SRCROOT)shared/drifters/drifters_io.F90
	cp $(SRCROOT)shared/drifters/drifters_io.F90 .
./fms_io.F90: $(SRCROOT)shared/fms/fms_io.F90
	cp $(SRCROOT)shared/fms/fms_io.F90 .
./mpp_global_reduce.h: $(SRCROOT)shared/mpp/include/mpp_global_reduce.h
	cp $(SRCROOT)shared/mpp/include/mpp_global_reduce.h .
./group_update_pack.inc: $(SRCROOT)shared/mpp/include/group_update_pack.inc
	cp $(SRCROOT)shared/mpp/include/group_update_pack.inc .
./interpolator.F90: $(SRCROOT)shared/interpolator/interpolator.F90
	cp $(SRCROOT)shared/interpolator/interpolator.F90 .
./test_xgrid.F90: $(SRCROOT)shared/exchange/test_xgrid.F90
	cp $(SRCROOT)shared/exchange/test_xgrid.F90 .
./mpp_define_nest_domains.inc: $(SRCROOT)shared/mpp/include/mpp_define_nest_domains.inc
	cp $(SRCROOT)shared/mpp/include/mpp_define_nest_domains.inc .
./mpp_alltoall_sma.h: $(SRCROOT)shared/mpp/include/mpp_alltoall_sma.h
	cp $(SRCROOT)shared/mpp/include/mpp_alltoall_sma.h .
./mpp_alltoall_mpi.h: $(SRCROOT)shared/mpp/include/mpp_alltoall_mpi.h
	cp $(SRCROOT)shared/mpp/include/mpp_alltoall_mpi.h .
./platform.F90: $(SRCROOT)shared/platform/platform.F90
	cp $(SRCROOT)shared/platform/platform.F90 .
./mpp_chksum.h: $(SRCROOT)shared/mpp/include/mpp_chksum.h
	cp $(SRCROOT)shared/mpp/include/mpp_chksum.h .
./mpp_sum_nocomm.h: $(SRCROOT)shared/mpp/include/mpp_sum_nocomm.h
	cp $(SRCROOT)shared/mpp/include/mpp_sum_nocomm.h .
./mosaic.F90: $(SRCROOT)shared/mosaic/mosaic.F90
	cp $(SRCROOT)shared/mosaic/mosaic.F90 .
./mpp_read_compressed.h: $(SRCROOT)shared/mpp/include/mpp_read_compressed.h
	cp $(SRCROOT)shared/mpp/include/mpp_read_compressed.h .
./gaussian_topog.F90: $(SRCROOT)shared/topography/gaussian_topog.F90
	cp $(SRCROOT)shared/topography/gaussian_topog.F90 .
./mpp_get_boundary.h: $(SRCROOT)shared/mpp/include/mpp_get_boundary.h
	cp $(SRCROOT)shared/mpp/include/mpp_get_boundary.h .
./mpp_comm_sma.inc: $(SRCROOT)shared/mpp/include/mpp_comm_sma.inc
	cp $(SRCROOT)shared/mpp/include/mpp_comm_sma.inc .
./mpp_io_misc.inc: $(SRCROOT)shared/mpp/include/mpp_io_misc.inc
	cp $(SRCROOT)shared/mpp/include/mpp_io_misc.inc .
./write_data.inc: $(SRCROOT)shared/fms/write_data.inc
	cp $(SRCROOT)shared/fms/write_data.inc .
./fft.F90: $(SRCROOT)shared/fft/fft.F90
	cp $(SRCROOT)shared/fft/fft.F90 .
./cloud_interpolator.F90: $(SRCROOT)shared/drifters/cloud_interpolator.F90
	cp $(SRCROOT)shared/drifters/cloud_interpolator.F90 .
./drifters_compute_k.h: $(SRCROOT)shared/drifters/drifters_compute_k.h
	cp $(SRCROOT)shared/drifters/drifters_compute_k.h .
./mpp_chksum_int.h: $(SRCROOT)shared/mpp/include/mpp_chksum_int.h
	cp $(SRCROOT)shared/mpp/include/mpp_chksum_int.h .
./mpp_util_sma.inc: $(SRCROOT)shared/mpp/include/mpp_util_sma.inc
	cp $(SRCROOT)shared/mpp/include/mpp_util_sma.inc .
./mpp_read_2Ddecomp.h: $(SRCROOT)shared/mpp/include/mpp_read_2Ddecomp.h
	cp $(SRCROOT)shared/mpp/include/mpp_read_2Ddecomp.h .
./oda_types.F90: $(SRCROOT)shared/oda_tools/oda_types.F90
	cp $(SRCROOT)shared/oda_tools/oda_types.F90 .
./mpp_do_update.h: $(SRCROOT)shared/mpp/include/mpp_do_update.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_update.h .
./atmos_ocean_fluxes.F90: $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90
	cp $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 .
./random_numbers.F90: $(SRCROOT)shared/random_numbers/random_numbers.F90
	cp $(SRCROOT)shared/random_numbers/random_numbers.F90 .
./mpp_data.F90: $(SRCROOT)shared/mpp/mpp_data.F90
	cp $(SRCROOT)shared/mpp/mpp_data.F90 .
./amip_interp.F90: $(SRCROOT)shared/amip_interp/amip_interp.F90
	cp $(SRCROOT)shared/amip_interp/amip_interp.F90 .
./mpp_gather.h: $(SRCROOT)shared/mpp/include/mpp_gather.h
	cp $(SRCROOT)shared/mpp/include/mpp_gather.h .
./mpp_global_sum.h: $(SRCROOT)shared/mpp/include/mpp_global_sum.h
	cp $(SRCROOT)shared/mpp/include/mpp_global_sum.h .
./gradient.F90: $(SRCROOT)shared/mosaic/gradient.F90
	cp $(SRCROOT)shared/mosaic/gradient.F90 .
./horiz_interp_type.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 .
./mpp_data_sma.inc: $(SRCROOT)shared/mpp/include/mpp_data_sma.inc
	cp $(SRCROOT)shared/mpp/include/mpp_data_sma.inc .
./mpp_reduce_sma.h: $(SRCROOT)shared/mpp/include/mpp_reduce_sma.h
	cp $(SRCROOT)shared/mpp/include/mpp_reduce_sma.h .
./mpp_update_domains2D.h: $(SRCROOT)shared/mpp/include/mpp_update_domains2D.h
	cp $(SRCROOT)shared/mpp/include/mpp_update_domains2D.h .
./mpp_sum.inc: $(SRCROOT)shared/mpp/include/mpp_sum.inc
	cp $(SRCROOT)shared/mpp/include/mpp_sum.inc .
./drifters_set_field.h: $(SRCROOT)shared/drifters/drifters_set_field.h
	cp $(SRCROOT)shared/drifters/drifters_set_field.h .
./station_data.F90: $(SRCROOT)shared/station_data/station_data.F90
	cp $(SRCROOT)shared/station_data/station_data.F90 .
./fms_platform.h: $(SRCROOT)shared/include/fms_platform.h
	cp $(SRCROOT)shared/include/fms_platform.h .
./mpp_sum_mpi.h: $(SRCROOT)shared/mpp/include/mpp_sum_mpi.h
	cp $(SRCROOT)shared/mpp/include/mpp_sum_mpi.h .
./mpp_utilities.F90: $(SRCROOT)shared/mpp/mpp_utilities.F90
	cp $(SRCROOT)shared/mpp/mpp_utilities.F90 .
./file_version.h: $(SRCROOT)shared/include/file_version.h
	cp $(SRCROOT)shared/include/file_version.h .
./drifters_comm.F90: $(SRCROOT)shared/drifters/drifters_comm.F90
	cp $(SRCROOT)shared/drifters/drifters_comm.F90 .
./mpp_comm_nocomm.inc: $(SRCROOT)shared/mpp/include/mpp_comm_nocomm.inc
	cp $(SRCROOT)shared/mpp/include/mpp_comm_nocomm.inc .
./mpp_write_unlimited_axis.h: $(SRCROOT)shared/mpp/include/mpp_write_unlimited_axis.h
	cp $(SRCROOT)shared/mpp/include/mpp_write_unlimited_axis.h .
./mpp_transmit_sma.h: $(SRCROOT)shared/mpp/include/mpp_transmit_sma.h
	cp $(SRCROOT)shared/mpp/include/mpp_transmit_sma.h .
./axis_utils.F90: $(SRCROOT)shared/axis_utils/axis_utils.F90
	cp $(SRCROOT)shared/axis_utils/axis_utils.F90 .
./read_data_4d.inc: $(SRCROOT)shared/fms/read_data_4d.inc
	cp $(SRCROOT)shared/fms/read_data_4d.inc .
./mpp_do_get_boundary.h: $(SRCROOT)shared/mpp/include/mpp_do_get_boundary.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_get_boundary.h .
./data_override.F90: $(SRCROOT)shared/data_override/data_override.F90
	cp $(SRCROOT)shared/data_override/data_override.F90 .
./time_manager.F90: $(SRCROOT)shared/time_manager/time_manager.F90
	cp $(SRCROOT)shared/time_manager/time_manager.F90 .
./diag_manifest.F90: $(SRCROOT)shared/diag_manager/diag_manifest.F90
	cp $(SRCROOT)shared/diag_manager/diag_manifest.F90 .
./test_mpp_domains.F90: $(SRCROOT)shared/mpp/test_mpp_domains.F90
	cp $(SRCROOT)shared/mpp/test_mpp_domains.F90 .
./quicksort.F90: $(SRCROOT)shared/drifters/quicksort.F90
	cp $(SRCROOT)shared/drifters/quicksort.F90 .
./interp.h: $(SRCROOT)shared/mosaic/interp.h
	cp $(SRCROOT)shared/mosaic/interp.h .
./diag_axis.F90: $(SRCROOT)shared/diag_manager/diag_axis.F90
	cp $(SRCROOT)shared/diag_manager/diag_axis.F90 .
./mpp_scatter.h: $(SRCROOT)shared/mpp/include/mpp_scatter.h
	cp $(SRCROOT)shared/mpp/include/mpp_scatter.h .
./mpp_domains_comm.inc: $(SRCROOT)shared/mpp/include/mpp_domains_comm.inc
	cp $(SRCROOT)shared/mpp/include/mpp_domains_comm.inc .
./horiz_interp_conserve.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 .
./mpp_reduce_nocomm.h: $(SRCROOT)shared/mpp/include/mpp_reduce_nocomm.h
	cp $(SRCROOT)shared/mpp/include/mpp_reduce_nocomm.h .
./fft99.F90: $(SRCROOT)shared/fft/fft99.F90
	cp $(SRCROOT)shared/fft/fft99.F90 .
./mpp_unstruct_pass_data.h: $(SRCROOT)shared/mpp/include/mpp_unstruct_pass_data.h
	cp $(SRCROOT)shared/mpp/include/mpp_unstruct_pass_data.h .
./mpp_do_global_field.h: $(SRCROOT)shared/mpp/include/mpp_do_global_field.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_global_field.h .
./mpp_do_updateV.h: $(SRCROOT)shared/mpp/include/mpp_do_updateV.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_updateV.h .
./oda_core_ecda.F90: $(SRCROOT)shared/oda_tools/oda_core_ecda.F90
	cp $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 .
./mpp_error_s_s.h: $(SRCROOT)shared/mpp/include/mpp_error_s_s.h
	cp $(SRCROOT)shared/mpp/include/mpp_error_s_s.h .
./mpp_domains_util.inc: $(SRCROOT)shared/mpp/include/mpp_domains_util.inc
	cp $(SRCROOT)shared/mpp/include/mpp_domains_util.inc .
./mpp_do_update_nest.h: $(SRCROOT)shared/mpp/include/mpp_do_update_nest.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_update_nest.h .
./diag_util.F90: $(SRCROOT)shared/diag_manager/diag_util.F90
	cp $(SRCROOT)shared/diag_manager/diag_util.F90 .
./diag_manager.F90: $(SRCROOT)shared/diag_manager/diag_manager.F90
	cp $(SRCROOT)shared/diag_manager/diag_manager.F90 .
./mpp_write_2Ddecomp.h: $(SRCROOT)shared/mpp/include/mpp_write_2Ddecomp.h
	cp $(SRCROOT)shared/mpp/include/mpp_write_2Ddecomp.h .
./column_diagnostics.F90: $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90
	cp $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 .
./mpp_io_connect.inc: $(SRCROOT)shared/mpp/include/mpp_io_connect.inc
	cp $(SRCROOT)shared/mpp/include/mpp_io_connect.inc .
./horiz_interp_bicubic.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 .
./fms.F90: $(SRCROOT)shared/fms/fms.F90
	cp $(SRCROOT)shared/fms/fms.F90 .
./astronomy.F90: $(SRCROOT)shared/astronomy/astronomy.F90
	cp $(SRCROOT)shared/astronomy/astronomy.F90 .
./mpp_sum_sma.h: $(SRCROOT)shared/mpp/include/mpp_sum_sma.h
	cp $(SRCROOT)shared/mpp/include/mpp_sum_sma.h .
./block_control.F90: $(SRCROOT)shared/block_control/block_control.F90
	cp $(SRCROOT)shared/block_control/block_control.F90 .
./mpp_util_nocomm.inc: $(SRCROOT)shared/mpp/include/mpp_util_nocomm.inc
	cp $(SRCROOT)shared/mpp/include/mpp_util_nocomm.inc .
./mpp_domains_reduce.inc: $(SRCROOT)shared/mpp/include/mpp_domains_reduce.inc
	cp $(SRCROOT)shared/mpp/include/mpp_domains_reduce.inc .
./field_manager.F90: $(SRCROOT)shared/field_manager/field_manager.F90
	cp $(SRCROOT)shared/field_manager/field_manager.F90 .
./xbt_drop_rate_adjust.f90: $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90
	cp $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 .
./fm_util.F90: $(SRCROOT)shared/field_manager/fm_util.F90
	cp $(SRCROOT)shared/field_manager/fm_util.F90 .
./mpp_error_s_a.h: $(SRCROOT)shared/mpp/include/mpp_error_s_a.h
	cp $(SRCROOT)shared/mpp/include/mpp_error_s_a.h .
./sat_vapor_pres.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 .
./mpp_io.F90: $(SRCROOT)shared/mpp/mpp_io.F90
	cp $(SRCROOT)shared/mpp/mpp_io.F90 .
./mpp_unstruct_domain.inc: $(SRCROOT)shared/mpp/include/mpp_unstruct_domain.inc
	cp $(SRCROOT)shared/mpp/include/mpp_unstruct_domain.inc .
./mpp_util_mpi.inc: $(SRCROOT)shared/mpp/include/mpp_util_mpi.inc
	cp $(SRCROOT)shared/mpp/include/mpp_util_mpi.inc .
./ensemble_manager.F90: $(SRCROOT)shared/coupler/ensemble_manager.F90
	cp $(SRCROOT)shared/coupler/ensemble_manager.F90 .
./mpp.F90: $(SRCROOT)shared/mpp/mpp.F90
	cp $(SRCROOT)shared/mpp/mpp.F90 .
./tridiagonal.F90: $(SRCROOT)shared/tridiagonal/tridiagonal.F90
	cp $(SRCROOT)shared/tridiagonal/tridiagonal.F90 .
./memuse.c: $(SRCROOT)shared/memutils/memuse.c
	cp $(SRCROOT)shared/memutils/memuse.c .
./drifters_core.F90: $(SRCROOT)shared/drifters/drifters_core.F90
	cp $(SRCROOT)shared/drifters/drifters_core.F90 .
./mpp_data_nocomm.inc: $(SRCROOT)shared/mpp/include/mpp_data_nocomm.inc
	cp $(SRCROOT)shared/mpp/include/mpp_data_nocomm.inc .
./diag_table.F90: $(SRCROOT)shared/diag_manager/diag_table.F90
	cp $(SRCROOT)shared/diag_manager/diag_table.F90 .
./mpp_data_mpi.inc: $(SRCROOT)shared/mpp/include/mpp_data_mpi.inc
	cp $(SRCROOT)shared/mpp/include/mpp_data_mpi.inc .
./mpp_domains_misc.inc: $(SRCROOT)shared/mpp/include/mpp_domains_misc.inc
	cp $(SRCROOT)shared/mpp/include/mpp_domains_misc.inc .
./sat_vapor_pres_k.F90: $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90
	cp $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 .
./gradient_c2l.c: $(SRCROOT)shared/mosaic/gradient_c2l.c
	cp $(SRCROOT)shared/mosaic/gradient_c2l.c .
./test_mpp_io.F90: $(SRCROOT)shared/mpp/test_mpp_io.F90
	cp $(SRCROOT)shared/mpp/test_mpp_io.F90 .
./diag_grid.F90: $(SRCROOT)shared/diag_manager/diag_grid.F90
	cp $(SRCROOT)shared/diag_manager/diag_grid.F90 .
./mosaic_util.h: $(SRCROOT)shared/mosaic/mosaic_util.h
	cp $(SRCROOT)shared/mosaic/mosaic_util.h .
./mpp_domains.F90: $(SRCROOT)shared/mpp/mpp_domains.F90
	cp $(SRCROOT)shared/mpp/mpp_domains.F90 .
./mpp_read_distributed_ascii.inc: $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.inc
	cp $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.inc .
./diag_data.F90: $(SRCROOT)shared/diag_manager/diag_data.F90
	cp $(SRCROOT)shared/diag_manager/diag_data.F90 .
./mpp_memutils.F90: $(SRCROOT)shared/mpp/mpp_memutils.F90
	cp $(SRCROOT)shared/mpp/mpp_memutils.F90 .
./mpp_parameter.F90: $(SRCROOT)shared/mpp/mpp_parameter.F90
	cp $(SRCROOT)shared/mpp/mpp_parameter.F90 .
./mpp_update_nest_domains.h: $(SRCROOT)shared/mpp/include/mpp_update_nest_domains.h
	cp $(SRCROOT)shared/mpp/include/mpp_update_nest_domains.h .
./mpp_write_compressed.h: $(SRCROOT)shared/mpp/include/mpp_write_compressed.h
	cp $(SRCROOT)shared/mpp/include/mpp_write_compressed.h .
./mpp_do_checkV.h: $(SRCROOT)shared/mpp/include/mpp_do_checkV.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_checkV.h .
./mpp_group_update.h: $(SRCROOT)shared/mpp/include/mpp_group_update.h
	cp $(SRCROOT)shared/mpp/include/mpp_group_update.h .
./mpp_do_check.h: $(SRCROOT)shared/mpp/include/mpp_do_check.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_check.h .
./test_fms_io.F90: $(SRCROOT)shared/fms/test_fms_io.F90
	cp $(SRCROOT)shared/fms/test_fms_io.F90 .
./gradient_c2l.h: $(SRCROOT)shared/mosaic/gradient_c2l.h
	cp $(SRCROOT)shared/mosaic/gradient_c2l.h .
./horiz_interp_spherical.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 .
./memutils.F90: $(SRCROOT)shared/memutils/memutils.F90
	cp $(SRCROOT)shared/memutils/memutils.F90 .
./MersenneTwister.F90: $(SRCROOT)shared/random_numbers/MersenneTwister.F90
	cp $(SRCROOT)shared/random_numbers/MersenneTwister.F90 .
./mpp_comm_mpi.inc: $(SRCROOT)shared/mpp/include/mpp_comm_mpi.inc
	cp $(SRCROOT)shared/mpp/include/mpp_comm_mpi.inc .
./constants.F90: $(SRCROOT)shared/constants/constants.F90
	cp $(SRCROOT)shared/constants/constants.F90 .
./tracer_manager.F90: $(SRCROOT)shared/tracer_manager/tracer_manager.F90
	cp $(SRCROOT)shared/tracer_manager/tracer_manager.F90 .
./topography.F90: $(SRCROOT)shared/topography/topography.F90
	cp $(SRCROOT)shared/topography/topography.F90 .
./test_mpp.F90: $(SRCROOT)shared/mpp/test_mpp.F90
	cp $(SRCROOT)shared/mpp/test_mpp.F90 .
./mpp_io_read.inc: $(SRCROOT)shared/mpp/include/mpp_io_read.inc
	cp $(SRCROOT)shared/mpp/include/mpp_io_read.inc .
./horiz_interp_bilinear.F90: $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 .
./affinity.c: $(SRCROOT)shared/mpp/affinity.c
	cp $(SRCROOT)shared/mpp/affinity.c .
./mosaic_util.c: $(SRCROOT)shared/mosaic/mosaic_util.c
	cp $(SRCROOT)shared/mosaic/mosaic_util.c .
./xgrid.F90: $(SRCROOT)shared/exchange/xgrid.F90
	cp $(SRCROOT)shared/exchange/xgrid.F90 .
./constant.h: $(SRCROOT)shared/mosaic/constant.h
	cp $(SRCROOT)shared/mosaic/constant.h .
./read_mosaic.c: $(SRCROOT)shared/mosaic/read_mosaic.c
	cp $(SRCROOT)shared/mosaic/read_mosaic.c .
./mpp_alltoall_nocomm.h: $(SRCROOT)shared/mpp/include/mpp_alltoall_nocomm.h
	cp $(SRCROOT)shared/mpp/include/mpp_alltoall_nocomm.h .
./mpp_do_update_nonblock.h: $(SRCROOT)shared/mpp/include/mpp_do_update_nonblock.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_update_nonblock.h .
./interp.c: $(SRCROOT)shared/mosaic/interp.c
	cp $(SRCROOT)shared/mosaic/interp.c .
./mpp_comm.inc: $(SRCROOT)shared/mpp/include/mpp_comm.inc
	cp $(SRCROOT)shared/mpp/include/mpp_comm.inc .
./parse.inc: $(SRCROOT)shared/field_manager/parse.inc
	cp $(SRCROOT)shared/field_manager/parse.inc .
./mpp_global_field.h: $(SRCROOT)shared/mpp/include/mpp_global_field.h
	cp $(SRCROOT)shared/mpp/include/mpp_global_field.h .
./mpp_read_distributed_ascii.h: $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.h
	cp $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.h .
./mpp_error_a_s.h: $(SRCROOT)shared/mpp/include/mpp_error_a_s.h
	cp $(SRCROOT)shared/mpp/include/mpp_error_a_s.h .
./mpp_update_domains2D_nonblock.h: $(SRCROOT)shared/mpp/include/mpp_update_domains2D_nonblock.h
	cp $(SRCROOT)shared/mpp/include/mpp_update_domains2D_nonblock.h .
./drifters_input.F90: $(SRCROOT)shared/drifters/drifters_input.F90
	cp $(SRCROOT)shared/drifters/drifters_input.F90 .
./threadloc.c: $(SRCROOT)shared/mpp/threadloc.c
	cp $(SRCROOT)shared/mpp/threadloc.c .
./fms_switches.h: $(SRCROOT)shared/drifters/fms_switches.h
	cp $(SRCROOT)shared/drifters/fms_switches.h .
./read_mosaic.h: $(SRCROOT)shared/mosaic/read_mosaic.h
	cp $(SRCROOT)shared/mosaic/read_mosaic.h .
./mpp_reduce_mpi.h: $(SRCROOT)shared/mpp/include/mpp_reduce_mpi.h
	cp $(SRCROOT)shared/mpp/include/mpp_reduce_mpi.h .
./mpp_io_util.inc: $(SRCROOT)shared/mpp/include/mpp_io_util.inc
	cp $(SRCROOT)shared/mpp/include/mpp_io_util.inc .
./mpp_do_redistribute.h: $(SRCROOT)shared/mpp/include/mpp_do_redistribute.h
	cp $(SRCROOT)shared/mpp/include/mpp_do_redistribute.h .
./mpp_efp.F90: $(SRCROOT)shared/mpp/mpp_efp.F90
	cp $(SRCROOT)shared/mpp/mpp_efp.F90 .
./create_xgrid.h: $(SRCROOT)shared/mosaic/create_xgrid.h
	cp $(SRCROOT)shared/mosaic/create_xgrid.h .
./oda_core.F90: $(SRCROOT)shared/oda_tools/oda_core.F90
	cp $(SRCROOT)shared/oda_tools/oda_core.F90 .
./get_cal_time.F90: $(SRCROOT)shared/time_manager/get_cal_time.F90
	cp $(SRCROOT)shared/time_manager/get_cal_time.F90 .
./write_ocean_data.F90: $(SRCROOT)shared/oda_tools/write_ocean_data.F90
	cp $(SRCROOT)shared/oda_tools/write_ocean_data.F90 .
./mpp_transmit_mpi.h: $(SRCROOT)shared/mpp/include/mpp_transmit_mpi.h
	cp $(SRCROOT)shared/mpp/include/mpp_transmit_mpi.h .
./time_interp.F90: $(SRCROOT)shared/time_interp/time_interp.F90
	cp $(SRCROOT)shared/time_interp/time_interp.F90 .
./horiz_interp.F90: $(SRCROOT)shared/horiz_interp/horiz_interp.F90
	cp $(SRCROOT)shared/horiz_interp/horiz_interp.F90 .
./diag_output.F90: $(SRCROOT)shared/diag_manager/diag_output.F90
	cp $(SRCROOT)shared/diag_manager/diag_output.F90 .
./nsclock.c: $(SRCROOT)shared/mpp/nsclock.c
	cp $(SRCROOT)shared/mpp/nsclock.c .
./drifters.F90: $(SRCROOT)shared/drifters/drifters.F90
	cp $(SRCROOT)shared/drifters/drifters.F90 .
./mpp_transmit.inc: $(SRCROOT)shared/mpp/include/mpp_transmit.inc
	cp $(SRCROOT)shared/mpp/include/mpp_transmit.inc .
./grid.F90: $(SRCROOT)shared/mosaic/grid.F90
	cp $(SRCROOT)shared/mosaic/grid.F90 .
./mpp_transmit_nocomm.h: $(SRCROOT)shared/mpp/include/mpp_transmit_nocomm.h
	cp $(SRCROOT)shared/mpp/include/mpp_transmit_nocomm.h .
./coupler_types.F90: $(SRCROOT)shared/coupler/coupler_types.F90
	cp $(SRCROOT)shared/coupler/coupler_types.F90 .
./stock_constants.F90: $(SRCROOT)shared/exchange/stock_constants.F90
	cp $(SRCROOT)shared/exchange/stock_constants.F90 .
./mpp_global_field_ug.h: $(SRCROOT)shared/mpp/include/mpp_global_field_ug.h
	cp $(SRCROOT)shared/mpp/include/mpp_global_field_ug.h .
./drifters_push.h: $(SRCROOT)shared/drifters/drifters_push.h
	cp $(SRCROOT)shared/drifters/drifters_push.h .
./create_xgrid.c: $(SRCROOT)shared/mosaic/create_xgrid.c
	cp $(SRCROOT)shared/mosaic/create_xgrid.c .
./time_interp_external.F90: $(SRCROOT)shared/time_interp/time_interp_external.F90
	cp $(SRCROOT)shared/time_interp/time_interp_external.F90 .
SRC = $(SRCROOT)shared/mpp/mpp_efp.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)shared/time_interp/time_interp_external.F90 $(SRCROOT)shared/data_override/data_override.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/coupler/coupler_types.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)shared/drifters/cloud_interpolator.F90 $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/exchange/xgrid.F90 $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)shared/random_numbers/random_numbers.F90 $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)shared/astronomy/astronomy.F90 $(SRCROOT)shared/oda_tools/write_ocean_data.F90 $(SRCROOT)shared/topography/gaussian_topog.F90 $(SRCROOT)shared/oda_tools/oda_core.F90 $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)shared/station_data/station_data.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)shared/drifters/quicksort.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/axis_utils/axis_utils.F90 $(SRCROOT)shared/drifters/drifters_comm.F90 $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/amip_interp/amip_interp.F90 $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/diag_manager/diag_table.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)shared/drifters/drifters.F90 $(SRCROOT)shared/exchange/stock_constants.F90 $(SRCROOT)shared/tridiagonal/tridiagonal.F90 $(SRCROOT)shared/exchange/test_xgrid.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)shared/drifters/drifters_input.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)shared/drifters/drifters_io.F90 $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/oda_tools/oda_types.F90 $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/interpolator/interpolator.F90 $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 $(SRCROOT)shared/block_control/block_control.F90 $(SRCROOT)shared/diag_manager/diag_manifest.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)shared/random_numbers/MersenneTwister.F90 $(SRCROOT)shared/drifters/drifters_core.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)shared/coupler/ensemble_manager.F90 $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)shared/mpp/include/mpp_write.h $(SRCROOT)shared/fms/read_data_2d.inc $(SRCROOT)shared/mpp/include/mpp_io_write.inc $(SRCROOT)shared/mpp/include/mpp_io_connect.inc $(SRCROOT)shared/mpp/include/system_clock.h $(SRCROOT)shared/mpp/include/mpp_chksum_scalar.h $(SRCROOT)shared/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)shared/mpp/include/group_update_unpack.inc $(SRCROOT)shared/mpp/include/mpp_sum_sma.h $(SRCROOT)shared/mpp/include/mpp_util_nocomm.inc $(SRCROOT)shared/mpp/include/mpp_global_sum_tl.h $(SRCROOT)shared/mpp/include/mpp_domains_reduce.inc $(SRCROOT)shared/fms/read_data_3d.inc $(SRCROOT)shared/mpp/include/mpp_domains_define.inc $(SRCROOT)shared/mpp/include/mpp_error_s_a.h $(SRCROOT)shared/mpp/include/mpp_util.inc $(SRCROOT)shared/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)shared/mpp/include/mpp_util_mpi.inc $(SRCROOT)shared/mpp/include/mpp_error_a_a.h $(SRCROOT)shared/mpp/include/mpp_global_reduce.h $(SRCROOT)shared/mpp/include/group_update_pack.inc $(SRCROOT)shared/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)shared/mpp/include/mpp_data_nocomm.inc $(SRCROOT)shared/mpp/include/mpp_alltoall_sma.h $(SRCROOT)shared/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)shared/mpp/include/mpp_data_mpi.inc $(SRCROOT)shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)shared/mpp/include/mpp_chksum.h $(SRCROOT)shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)shared/mosaic/mosaic_util.h $(SRCROOT)shared/mpp/include/mpp_read_compressed.h $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)shared/mpp/include/mpp_get_boundary.h $(SRCROOT)shared/mpp/include/mpp_comm_sma.inc $(SRCROOT)shared/mpp/include/mpp_io_misc.inc $(SRCROOT)shared/fms/write_data.inc $(SRCROOT)shared/mpp/include/mpp_update_nest_domains.h $(SRCROOT)shared/drifters/drifters_compute_k.h $(SRCROOT)shared/mpp/include/mpp_chksum_int.h $(SRCROOT)shared/mpp/include/mpp_write_compressed.h $(SRCROOT)shared/mpp/include/mpp_util_sma.inc $(SRCROOT)shared/mpp/include/mpp_do_checkV.h $(SRCROOT)shared/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)shared/mpp/include/mpp_group_update.h $(SRCROOT)shared/mpp/include/mpp_do_check.h $(SRCROOT)shared/mosaic/gradient_c2l.h $(SRCROOT)shared/mpp/include/mpp_comm_mpi.inc $(SRCROOT)shared/mpp/include/mpp_do_update.h $(SRCROOT)shared/mpp/include/mpp_gather.h $(SRCROOT)shared/mpp/include/mpp_global_sum.h $(SRCROOT)shared/mpp/include/mpp_io_read.inc $(SRCROOT)shared/mpp/include/mpp_data_sma.inc $(SRCROOT)shared/mpp/include/mpp_reduce_sma.h $(SRCROOT)shared/mpp/include/mpp_update_domains2D.h $(SRCROOT)shared/drifters/drifters_set_field.h $(SRCROOT)shared/mpp/include/mpp_sum.inc $(SRCROOT)shared/mosaic/constant.h $(SRCROOT)shared/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)shared/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/mpp/include/mpp_comm.inc $(SRCROOT)shared/field_manager/parse.inc $(SRCROOT)shared/mpp/include/mpp_global_field.h $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/mpp/include/mpp_error_a_s.h $(SRCROOT)shared/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)shared/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)shared/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)shared/drifters/fms_switches.h $(SRCROOT)shared/mosaic/read_mosaic.h $(SRCROOT)shared/fms/read_data_4d.inc $(SRCROOT)shared/mpp/include/mpp_transmit_sma.h $(SRCROOT)shared/mpp/include/mpp_do_get_boundary.h $(SRCROOT)shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)shared/mpp/include/mpp_io_util.inc $(SRCROOT)shared/mpp/include/mpp_do_redistribute.h $(SRCROOT)shared/mosaic/create_xgrid.h $(SRCROOT)shared/mosaic/interp.h $(SRCROOT)shared/mpp/include/mpp_transmit_mpi.h $(SRCROOT)shared/mpp/include/mpp_scatter.h $(SRCROOT)shared/mpp/include/mpp_domains_comm.inc $(SRCROOT)shared/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)shared/mpp/include/mpp_transmit.inc $(SRCROOT)shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)shared/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)shared/mpp/include/mpp_do_global_field.h $(SRCROOT)shared/mpp/include/mpp_do_updateV.h $(SRCROOT)shared/mpp/include/mpp_global_field_ug.h $(SRCROOT)shared/mpp/include/mpp_error_s_s.h $(SRCROOT)shared/mpp/include/mpp_domains_util.inc $(SRCROOT)shared/drifters/drifters_push.h $(SRCROOT)shared/mpp/include/mpp_do_update_nest.h $(SRCROOT)shared/mpp/include/mpp_write_2Ddecomp.h
OBJ = mpp_efp.o mpp_data.o affinity.o gradient_c2l.o xbt_drop_rate_adjust.o fms_io.o time_interp_external.o data_override.o tracer_manager.o coupler_types.o mpp_memutils.o cloud_interpolator.o diag_util.o mosaic_util.o xgrid.o time_interp.o random_numbers.o mosaic.o memuse.o memutils.o threadloc.o horiz_interp_bicubic.o astronomy.o write_ocean_data.o gaussian_topog.o oda_core.o diag_axis.o station_data.o diag_manager.o mpp_utilities.o mpp_pset.o fms.o oda_core_ecda.o test_mpp_io.o create_xgrid.o interp.o diag_output.o horiz_interp_conserve.o quicksort.o mpp_parameter.o axis_utils.o drifters_comm.o atmos_ocean_fluxes.o platform.o amip_interp.o mpp_domains.o diag_grid.o diag_table.o time_manager.o drifters.o stock_constants.o tridiagonal.o test_xgrid.o mpp_io.o sat_vapor_pres_k.o horiz_interp_spherical.o topography.o test_horiz_interp.o horiz_interp.o gradient.o fft.o fm_util.o drifters_input.o horiz_interp_bilinear.o drifters_io.o diag_data.o oda_types.o grid.o get_cal_time.o read_mosaic.o constants.o sat_vapor_pres.o test_mpp_pset.o test_fms_io.o field_manager.o interpolator.o column_diagnostics.o block_control.o diag_manifest.o nsclock.o mpp.o fft99.o MersenneTwister.o drifters_core.o horiz_interp_type.o test_mpp_domains.o ensemble_manager.o test_mpp.o
OFF = $(SRCROOT)shared/mpp/include/mpp_write.h $(SRCROOT)shared/fms/read_data_2d.inc $(SRCROOT)shared/mpp/include/mpp_io_write.inc $(SRCROOT)shared/mpp/include/system_clock.h $(SRCROOT)shared/mpp/include/mpp_chksum_scalar.h $(SRCROOT)shared/mpp/include/mpp_do_updateV_nonblock.h $(SRCROOT)shared/mpp/include/group_update_unpack.inc $(SRCROOT)shared/mpp/include/mpp_global_sum_tl.h $(SRCROOT)shared/fms/read_data_3d.inc $(SRCROOT)shared/mpp/include/mpp_domains_define.inc $(SRCROOT)shared/horiz_interp/test_horiz_interp.F90 $(SRCROOT)shared/mpp/include/mpp_util.inc $(SRCROOT)shared/mpp/mpp_pset.F90 $(SRCROOT)shared/mpp/include/mpp_error_a_a.h $(SRCROOT)shared/mpp/test_mpp_pset.F90 $(SRCROOT)shared/drifters/drifters_io.F90 $(SRCROOT)shared/fms/fms_io.F90 $(SRCROOT)shared/mpp/include/mpp_global_reduce.h $(SRCROOT)shared/mpp/include/group_update_pack.inc $(SRCROOT)shared/interpolator/interpolator.F90 $(SRCROOT)shared/exchange/test_xgrid.F90 $(SRCROOT)shared/mpp/include/mpp_define_nest_domains.inc $(SRCROOT)shared/mpp/include/mpp_alltoall_sma.h $(SRCROOT)shared/mpp/include/mpp_alltoall_mpi.h $(SRCROOT)shared/platform/platform.F90 $(SRCROOT)shared/mpp/include/mpp_chksum.h $(SRCROOT)shared/mpp/include/mpp_sum_nocomm.h $(SRCROOT)shared/mosaic/mosaic.F90 $(SRCROOT)shared/mpp/include/mpp_read_compressed.h $(SRCROOT)shared/topography/gaussian_topog.F90 $(SRCROOT)shared/mpp/include/mpp_get_boundary.h $(SRCROOT)shared/mpp/include/mpp_comm_sma.inc $(SRCROOT)shared/mpp/include/mpp_io_misc.inc $(SRCROOT)shared/fms/write_data.inc $(SRCROOT)shared/fft/fft.F90 $(SRCROOT)shared/drifters/cloud_interpolator.F90 $(SRCROOT)shared/drifters/drifters_compute_k.h $(SRCROOT)shared/mpp/include/mpp_chksum_int.h $(SRCROOT)shared/mpp/include/mpp_util_sma.inc $(SRCROOT)shared/mpp/include/mpp_read_2Ddecomp.h $(SRCROOT)shared/oda_tools/oda_types.F90 $(SRCROOT)shared/mpp/include/mpp_do_update.h $(SRCROOT)shared/coupler/atmos_ocean_fluxes.F90 $(SRCROOT)shared/random_numbers/random_numbers.F90 $(SRCROOT)shared/mpp/mpp_data.F90 $(SRCROOT)shared/amip_interp/amip_interp.F90 $(SRCROOT)shared/mpp/include/mpp_gather.h $(SRCROOT)shared/mpp/include/mpp_global_sum.h $(SRCROOT)shared/mosaic/gradient.F90 $(SRCROOT)shared/horiz_interp/horiz_interp_type.F90 $(SRCROOT)shared/mpp/include/mpp_data_sma.inc $(SRCROOT)shared/mpp/include/mpp_reduce_sma.h $(SRCROOT)shared/mpp/include/mpp_update_domains2D.h $(SRCROOT)shared/mpp/include/mpp_sum.inc $(SRCROOT)shared/drifters/drifters_set_field.h $(SRCROOT)shared/station_data/station_data.F90 $(SRCROOT)shared/include/fms_platform.h $(SRCROOT)shared/mpp/include/mpp_sum_mpi.h $(SRCROOT)shared/mpp/mpp_utilities.F90 $(SRCROOT)shared/include/file_version.h $(SRCROOT)shared/drifters/drifters_comm.F90 $(SRCROOT)shared/mpp/include/mpp_comm_nocomm.inc $(SRCROOT)shared/mpp/include/mpp_write_unlimited_axis.h $(SRCROOT)shared/mpp/include/mpp_transmit_sma.h $(SRCROOT)shared/axis_utils/axis_utils.F90 $(SRCROOT)shared/fms/read_data_4d.inc $(SRCROOT)shared/mpp/include/mpp_do_get_boundary.h $(SRCROOT)shared/data_override/data_override.F90 $(SRCROOT)shared/time_manager/time_manager.F90 $(SRCROOT)shared/diag_manager/diag_manifest.F90 $(SRCROOT)shared/mpp/test_mpp_domains.F90 $(SRCROOT)shared/drifters/quicksort.F90 $(SRCROOT)shared/mosaic/interp.h $(SRCROOT)shared/diag_manager/diag_axis.F90 $(SRCROOT)shared/mpp/include/mpp_scatter.h $(SRCROOT)shared/mpp/include/mpp_domains_comm.inc $(SRCROOT)shared/horiz_interp/horiz_interp_conserve.F90 $(SRCROOT)shared/mpp/include/mpp_reduce_nocomm.h $(SRCROOT)shared/fft/fft99.F90 $(SRCROOT)shared/mpp/include/mpp_unstruct_pass_data.h $(SRCROOT)shared/mpp/include/mpp_do_global_field.h $(SRCROOT)shared/mpp/include/mpp_do_updateV.h $(SRCROOT)shared/oda_tools/oda_core_ecda.F90 $(SRCROOT)shared/mpp/include/mpp_error_s_s.h $(SRCROOT)shared/mpp/include/mpp_domains_util.inc $(SRCROOT)shared/mpp/include/mpp_do_update_nest.h $(SRCROOT)shared/diag_manager/diag_util.F90 $(SRCROOT)shared/diag_manager/diag_manager.F90 $(SRCROOT)shared/mpp/include/mpp_write_2Ddecomp.h $(SRCROOT)shared/column_diagnostics/column_diagnostics.F90 $(SRCROOT)shared/mpp/include/mpp_io_connect.inc $(SRCROOT)shared/horiz_interp/horiz_interp_bicubic.F90 $(SRCROOT)shared/fms/fms.F90 $(SRCROOT)shared/astronomy/astronomy.F90 $(SRCROOT)shared/mpp/include/mpp_sum_sma.h $(SRCROOT)shared/block_control/block_control.F90 $(SRCROOT)shared/mpp/include/mpp_util_nocomm.inc $(SRCROOT)shared/mpp/include/mpp_domains_reduce.inc $(SRCROOT)shared/field_manager/field_manager.F90 $(SRCROOT)shared/oda_tools/xbt_drop_rate_adjust.f90 $(SRCROOT)shared/field_manager/fm_util.F90 $(SRCROOT)shared/mpp/include/mpp_error_s_a.h $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres.F90 $(SRCROOT)shared/mpp/mpp_io.F90 $(SRCROOT)shared/mpp/include/mpp_unstruct_domain.inc $(SRCROOT)shared/mpp/include/mpp_util_mpi.inc $(SRCROOT)shared/coupler/ensemble_manager.F90 $(SRCROOT)shared/mpp/mpp.F90 $(SRCROOT)shared/tridiagonal/tridiagonal.F90 $(SRCROOT)shared/memutils/memuse.c $(SRCROOT)shared/drifters/drifters_core.F90 $(SRCROOT)shared/mpp/include/mpp_data_nocomm.inc $(SRCROOT)shared/diag_manager/diag_table.F90 $(SRCROOT)shared/mpp/include/mpp_data_mpi.inc $(SRCROOT)shared/mpp/include/mpp_domains_misc.inc $(SRCROOT)shared/sat_vapor_pres/sat_vapor_pres_k.F90 $(SRCROOT)shared/mosaic/gradient_c2l.c $(SRCROOT)shared/mpp/test_mpp_io.F90 $(SRCROOT)shared/diag_manager/diag_grid.F90 $(SRCROOT)shared/mosaic/mosaic_util.h $(SRCROOT)shared/mpp/mpp_domains.F90 $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.inc $(SRCROOT)shared/diag_manager/diag_data.F90 $(SRCROOT)shared/mpp/mpp_memutils.F90 $(SRCROOT)shared/mpp/mpp_parameter.F90 $(SRCROOT)shared/mpp/include/mpp_update_nest_domains.h $(SRCROOT)shared/mpp/include/mpp_write_compressed.h $(SRCROOT)shared/mpp/include/mpp_do_checkV.h $(SRCROOT)shared/mpp/include/mpp_group_update.h $(SRCROOT)shared/mpp/include/mpp_do_check.h $(SRCROOT)shared/fms/test_fms_io.F90 $(SRCROOT)shared/mosaic/gradient_c2l.h $(SRCROOT)shared/horiz_interp/horiz_interp_spherical.F90 $(SRCROOT)shared/memutils/memutils.F90 $(SRCROOT)shared/random_numbers/MersenneTwister.F90 $(SRCROOT)shared/mpp/include/mpp_comm_mpi.inc $(SRCROOT)shared/constants/constants.F90 $(SRCROOT)shared/tracer_manager/tracer_manager.F90 $(SRCROOT)shared/topography/topography.F90 $(SRCROOT)shared/mpp/test_mpp.F90 $(SRCROOT)shared/mpp/include/mpp_io_read.inc $(SRCROOT)shared/horiz_interp/horiz_interp_bilinear.F90 $(SRCROOT)shared/mpp/affinity.c $(SRCROOT)shared/mosaic/mosaic_util.c $(SRCROOT)shared/exchange/xgrid.F90 $(SRCROOT)shared/mosaic/constant.h $(SRCROOT)shared/mosaic/read_mosaic.c $(SRCROOT)shared/mpp/include/mpp_alltoall_nocomm.h $(SRCROOT)shared/mpp/include/mpp_do_update_nonblock.h $(SRCROOT)shared/mosaic/interp.c $(SRCROOT)shared/mpp/include/mpp_comm.inc $(SRCROOT)shared/field_manager/parse.inc $(SRCROOT)shared/mpp/include/mpp_global_field.h $(SRCROOT)shared/mpp/include/mpp_read_distributed_ascii.h $(SRCROOT)shared/mpp/include/mpp_error_a_s.h $(SRCROOT)shared/mpp/include/mpp_update_domains2D_nonblock.h $(SRCROOT)shared/drifters/drifters_input.F90 $(SRCROOT)shared/mpp/threadloc.c $(SRCROOT)shared/drifters/fms_switches.h $(SRCROOT)shared/mosaic/read_mosaic.h $(SRCROOT)shared/mpp/include/mpp_reduce_mpi.h $(SRCROOT)shared/mpp/include/mpp_io_util.inc $(SRCROOT)shared/mpp/include/mpp_do_redistribute.h $(SRCROOT)shared/mpp/mpp_efp.F90 $(SRCROOT)shared/mosaic/create_xgrid.h $(SRCROOT)shared/oda_tools/oda_core.F90 $(SRCROOT)shared/time_manager/get_cal_time.F90 $(SRCROOT)shared/oda_tools/write_ocean_data.F90 $(SRCROOT)shared/mpp/include/mpp_transmit_mpi.h $(SRCROOT)shared/time_interp/time_interp.F90 $(SRCROOT)shared/horiz_interp/horiz_interp.F90 $(SRCROOT)shared/diag_manager/diag_output.F90 $(SRCROOT)shared/mpp/nsclock.c $(SRCROOT)shared/drifters/drifters.F90 $(SRCROOT)shared/mpp/include/mpp_transmit.inc $(SRCROOT)shared/mosaic/grid.F90 $(SRCROOT)shared/mpp/include/mpp_transmit_nocomm.h $(SRCROOT)shared/coupler/coupler_types.F90 $(SRCROOT)shared/exchange/stock_constants.F90 $(SRCROOT)shared/mpp/include/mpp_global_field_ug.h $(SRCROOT)shared/drifters/drifters_push.h $(SRCROOT)shared/mosaic/create_xgrid.c $(SRCROOT)shared/time_interp/time_interp_external.F90
clean: neat
	-rm -f .libfms.a.cppdefs $(OBJ) libfms.a zzz.exe

neat:
	-rm -f $(TMPFILES)

localize: $(OFF)
	cp $(OFF) .

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

libfms.a: $(OBJ)
	$(AR) $(ARFLAGS) libfms.a $(OBJ)

prog:
	ftn -i4 -r8 -g -traceback -O0 -I$(SRCROOT)shared/include test_unstructured_io_plus_diag.F90 -B./ -openmp  -L/opt/intel/composer_xe_2015.2.164/15.0/2.164/lib/intel64 -lifcoremt -lnetcdf -L./ -lfms -o zzz.exe
#test_diag.F90 -B./ -openmp  -L/opt/intel/composer_xe_2015.2.164/15.0/2.164/lib/intel64 -lifcoremt -lnetcdf -L./ -lfms -o zzz.exe

diag_test:
	ftn -i4 -r8 -B./ -openmp  -L/opt/intel/composer_xe_2015.2.164/15.0/2.164/lib/intel64 -lifcoremt -lnetcdf -L./ -lfms -o zzz.exe

